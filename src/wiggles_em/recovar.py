"""Reading RECOVAR's conformational-density job — the deconvolved route to occupancy.

RECOVAR is currently the one method in the compendium that ships a
**deconvolved** conformational density: ``recovar estimate_conformational_density``
deconvolves the particle distribution in latent space rather than histogramming
it, and ``recovar estimate_stable_states`` then finds the local maxima of that
density. That is the estimator family Evans *et al.* show recovering a middle
mode which histogramming loses, so it is the first place a trustworthy occupancy
can come from without new code on the method side.

**What is actually on disk**, read off RECOVAR's own source rather than guessed:

===============================================  ============================
``stable_state_all_coords.txt``                  local maxima, one row each
``stable_state_<i>_coords.txt``                  one maximum, one value a line
``density/data/deconv_density_knee.pkl``         the density — **a pickle**
``density/all_densities/deconv_density_<n>.pkl`` per-regularization, pickles
===============================================  ============================

The coordinates are written with ``numpy.savetxt`` and are plain text, so they
are read here directly. The density is not: it is a pickle holding a dict with
``density`` and ``latent_space_bounds`` keys, and this package does not unpickle
— reading a pickle runs whatever code is inside it, and a job directory can come
from a collaborator or a shared scratch volume. Same position as
:func:`~wiggles_em.heterogeneity.read_latent_table` takes, for the same reason.

**Peak height is not occupancy, and this module will not pretend otherwise.**
This is the subtlety that makes an obvious implementation wrong. RECOVAR gives a
density over latent space; a state's occupancy is that density *integrated over
the state's basin*, not its value at the peak. A sharp narrow peak and a broad
low one can have identical heights and very different populations. RECOVAR does
not ship that integration, and neither does this module — it reads what exists
and tells you what the remaining step is. Inventing a basin boundary here would
produce numbers that look measured and are not.

So the honest flow is: RECOVAR gives you the states, you produce weights for
them by whatever integration your analysis justifies, and
:func:`read_deconvolved_weights` reads them back in labelled
:attr:`~wiggles_em.populations.WeightSource.DECONVOLVED`.
"""

from __future__ import annotations

import ast
import re
import struct
from dataclasses import dataclass
from pathlib import Path

from wiggles_em.populations import Populations, WeightSource

#: What ``estimate_stable_states`` writes for the full set of maxima. A 2-D
#: array through ``np.savetxt``, so one state per line.
STABLE_STATES_FILE = "stable_state_all_coords.txt"

#: Per-state files. ``np.savetxt`` on a **1-D** array writes one value per
#: line, so each of these is a column of coordinates, not a row — a difference
#: that silently transposes a latent point if assumed the other way.
PER_STATE_PATTERN = re.compile(r"^stable_state_(\d+)_coords\.txt$")

#: The density files, which are pickles and are deliberately not read.
DENSITY_PATTERN = re.compile(r"^deconv_density_.*\.pkl$")


@dataclass(frozen=True)
class DensityJob:
    """What a RECOVAR density job left on disk that this package can use."""

    directory: Path
    stable_states: tuple[tuple[float, ...], ...] = ()
    """Latent coordinates of the density's local maxima, one per state."""

    density_files: tuple[Path, ...] = ()
    """Deconvolved densities found but **not read**. Their presence is the
    evidence that a deconvolution was run at all."""

    notes: tuple[str, ...] = ()

    @property
    def n_states(self) -> int:
        return len(self.stable_states)

    @property
    def deconvolution_was_run(self) -> bool:
        """Did this job deconvolve, or is it a bare pipeline output?

        The distinction that decides whether weights derived from it may be
        called :attr:`WeightSource.DECONVOLVED` at all.
        """
        return bool(self.density_files) or bool(self.stable_states)

    def report(self) -> str:
        lines = [f"RECOVAR density job at {self.directory}"]
        if self.stable_states:
            lines.append(f"  {self.n_states} stable state(s), from the density's local maxima:")
            for i, point in enumerate(self.stable_states):
                coords = ", ".join(f"{v:.4g}" for v in point)
                lines.append(f"    state {i}: ({coords})")
        else:
            lines.append("  No stable states found. Run `recovar estimate_stable_states`.")

        if self.density_files:
            names = ", ".join(p.name for p in self.density_files[:3])
            more = "" if len(self.density_files) <= 3 else f" (+{len(self.density_files) - 3} more)"
            lines += [
                "",
                f"  {len(self.density_files)} deconvolved density file(s) present: {names}{more}",
                "  NOT READ. These are Python pickles, and unpickling runs whatever code",
                "  the file contains — a job directory is not always your own. Nothing",
                "  here opens one.",
            ]
        lines += ["", *(f"  {n}" for n in self.notes)] if self.notes else []
        lines += [
            "",
            "  NO OCCUPANCY IS REPORTED HERE, and that is not an oversight. A state's",
            "  occupancy is the density integrated over that state's basin, not the",
            "  density at its peak — a sharp narrow peak and a broad low one can be the",
            "  same height and hold very different populations. RECOVAR does not ship",
            "  that integration and this package will not invent a basin boundary to",
            "  fake it. Integrate the density however your analysis justifies, then",
            "  read the weights back with read_deconvolved_weights().",
            "",
            export_instructions(),
        ]
        return "\n".join(lines)


def export_instructions() -> str:
    """The exact commands to get weights out of RECOVAR's pickle.

    Key names are taken from RECOVAR's ``estimate_stable_states.py``, which does
    ``dens_pkl["density"]`` and ``dens_pkl["latent_space_bounds"]``. Guessing
    them would have produced a snippet that fails on the user's first try.
    """
    return (
        "  To export, in an environment that has RECOVAR and numpy:\n"
        "      import pickle, numpy as np\n"
        "      d = pickle.load(open('density/data/deconv_density_knee.pkl', 'rb'))\n"
        "      density, bounds = d['density'], d['latent_space_bounds']\n"
        "      # integrate `density` over each state's basin -> one weight per state\n"
        "      np.save('state_weights.npy', weights)   # or np.savetxt(...)\n"
        "  Then: read_deconvolved_weights('state_weights.npy')"
    )


def _floats(text: str) -> list[list[float]]:
    """Rows of floats from a ``np.savetxt`` file, blank and comment lines dropped."""
    rows: list[list[float]] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            rows.append([float(part) for part in re.split(r"[,\t ]+", stripped)])
        except ValueError:
            continue  # a header line
    return rows


def read_stable_states(directory: str | Path) -> tuple[tuple[float, ...], ...]:
    """Latent coordinates of the density's local maxima.

    Prefers ``stable_state_all_coords.txt``, which holds every state as one row.
    Falls back to the per-state files, which are **column** vectors because
    ``np.savetxt`` writes a 1-D array one value per line — read as rows they
    would transpose each latent point into as many one-dimensional states.

    Returns an empty tuple when the job has no stable-states output, rather
    than raising: a density job without them is a normal intermediate state.
    """
    directory = Path(directory)
    combined = next(iter(sorted(directory.rglob(STABLE_STATES_FILE))), None)
    if combined is not None:
        rows = _floats(combined.read_text(errors="replace"))
        if rows:
            return tuple(tuple(r) for r in rows)

    per_state: dict[int, tuple[float, ...]] = {}
    for path in directory.rglob("stable_state_*_coords.txt"):
        match = PER_STATE_PATTERN.match(path.name)
        if match is None:
            continue
        rows = _floats(path.read_text(errors="replace"))
        if not rows:
            continue
        # One value per line is the 1-D savetxt shape; a single line means the
        # writer used a 2-D array with one row. Both describe one point.
        flat = [v for row in rows for v in row]
        per_state[int(match.group(1))] = tuple(flat)
    return tuple(per_state[k] for k in sorted(per_state))


def read_density_job(directory: str | Path) -> DensityJob:
    """Survey a RECOVAR density job without opening anything unsafe."""
    directory = Path(directory)
    if not directory.is_dir():
        return DensityJob(directory=directory, notes=(f"{directory} is not a directory",))

    density_files = tuple(
        sorted(p for p in directory.rglob("*.pkl") if DENSITY_PATTERN.match(p.name))
    )
    states = read_stable_states(directory)

    notes: list[str] = []
    if density_files and not states:
        notes.append(
            "A deconvolved density is present but no stable states were found. "
            "`recovar estimate_stable_states <density.pkl> -o <dir>` finds the maxima."
        )
    if states and not density_files:
        notes.append(
            "Stable states are present without the density they came from. The "
            "coordinates are usable; the weights cannot be derived without it."
        )
    return DensityJob(
        directory=directory,
        stable_states=states,
        density_files=density_files,
        notes=tuple(notes),
    )


def _read_npy_1d(path: Path) -> list[float]:
    """A 1-D ``.npy`` of floats, without numpy and without unpickling.

    Narrower than :mod:`wiggles_em.heterogeneity`'s 2-D latent reader on
    purpose: a weight vector is one-dimensional, and accepting a 2-D array here
    would quietly flatten a grid into a list of weights.
    """
    raw = path.read_bytes()
    if raw[:6] != b"\x93NUMPY":
        raise ValueError(f"{path} is not a .npy file")
    major = raw[6]
    size = 2 if major == 1 else 4
    start = 8 + size
    (header_len,) = struct.unpack("<H" if major == 1 else "<I", raw[8:start])
    header = ast.literal_eval(raw[start : start + header_len].decode("latin1").strip())

    descr, shape = header.get("descr"), header.get("shape")
    if not isinstance(descr, str) or descr.startswith(("|O", "O")):
        raise ValueError(
            f"{path} holds Python objects, which can only be read by unpickling, and "
            f"unpickling runs whatever code the file contains"
        )
    fmt = {"<f4": "f", "<f8": "d", "<i4": "i", "<i8": "q"}.get(descr)
    if fmt is None or header.get("fortran_order"):
        raise ValueError(f"{path}: unsupported layout {descr!r}")
    if not (isinstance(shape, tuple) and len(shape) == 1):
        raise ValueError(
            f"{path} has shape {shape}, but a weight vector is one-dimensional. If "
            f"this is the density grid, it still needs integrating over each state's "
            f"basin — one weight per state — before it can be read as occupancy."
        )
    body = raw[start + header_len :]
    (n,) = shape
    return [float(v) for v in struct.unpack(f"<{n}{fmt}", body[: n * struct.calcsize(fmt)])]


def read_deconvolved_weights(
    path: str | Path,
    *,
    uncertainty: list[float] | None = None,
    temperature_k: float | None = None,
) -> Populations:
    """Read exported per-state weights and label them as deconvolved.

    Accepts a 1-D ``.npy`` or a text file of one weight per line. The label is
    :attr:`WeightSource.DECONVOLVED` because the numbers came out of RECOVAR's
    deconvolution — **the caller is asserting that**, exactly as
    :mod:`wiggles_em.provenance` has callers assert how a volume was made.
    Nothing here inspects the values to decide; a flat histogram and a
    deconvolved density are the same list of floats.

    Raises:
        ValueError: unreadable, empty, or not one-dimensional.
    """
    path = Path(path)
    if path.suffix == ".npy":
        weights = _read_npy_1d(path)
    else:
        rows = _floats(path.read_text(errors="replace"))
        if any(len(r) > 1 for r in rows):
            raise ValueError(
                f"{path} has more than one value per line. A weight vector is one "
                f"weight per state; a grid needs integrating over each basin first."
            )
        weights = [r[0] for r in rows]
    if not weights:
        raise ValueError(f"{path}: no weights found")
    return Populations.declare(
        weights,
        WeightSource.DECONVOLVED,
        uncertainty=uncertainty,
        temperature_k=temperature_k,
    )
