"""Tests for the RECOVAR density reader.

Fixtures mirror what RECOVAR's own `estimate_stable_states.py` writes —
`np.savetxt` of a 2-D array for the combined file, and of a **1-D** array for
each per-state file, which puts one value per line. Getting that backwards
transposes a latent point into several one-dimensional states, so it is
asserted directly rather than assumed.
"""

from __future__ import annotations

import struct

import pytest

from wiggles_em.populations import WeightSource
from wiggles_em.recovar import (
    export_instructions,
    read_deconvolved_weights,
    read_density_job,
    read_stable_states,
)


def write_npy(path, values, descr="<f8"):
    """A minimal 1-D .npy, matching what numpy would write."""
    fmt = {"<f8": "d", "<f4": "f"}[descr]
    header = f"{{'descr': '{descr}', 'fortran_order': False, 'shape': ({len(values)},), }}"
    header += " " * ((64 - (10 + len(header)) % 64) % 64) + "\n"
    path.write_bytes(
        b"\x93NUMPY\x01\x00"
        + struct.pack("<H", len(header))
        + header.encode("latin1")
        + struct.pack(f"<{len(values)}{fmt}", *values)
    )
    return path


@pytest.fixture
def density_job(tmp_path):
    """A job directory shaped like RECOVAR's, with pickles left unopened."""

    def _make(*, states=((1.0, 2.0, 3.0, 4.0), (-1.0, 0.5, 0.0, 2.0)), densities=True):
        root = tmp_path / "recovar_out"
        (root / "density" / "data").mkdir(parents=True)
        if states:
            (root / "stable_state_all_coords.txt").write_text(
                "\n".join(" ".join(f"{v:.18e}" for v in s) for s in states) + "\n"
            )
            for i, s in enumerate(states):
                # 1-D savetxt: one value per line.
                (root / f"stable_state_{i}_coords.txt").write_text(
                    "\n".join(f"{v:.18e}" for v in s) + "\n"
                )
        if densities:
            (root / "density" / "data" / "deconv_density_knee.pkl").write_bytes(b"\x80\x04junk")
        return root

    return _make


# ── the stable states ────────────────────────────────────────────────────────


def test_states_are_read_from_the_combined_file(density_job):
    states = read_stable_states(density_job())

    assert states == ((1.0, 2.0, 3.0, 4.0), (-1.0, 0.5, 0.0, 2.0))


def test_per_state_files_are_columns_not_rows(density_job):
    """`np.savetxt` on a 1-D array writes one value per line.

    Read as rows, a single 4-D latent point becomes four 1-D states — a silent
    transpose that would put a trajectory endpoint in the wrong place.
    """
    root = density_job()
    (root / "stable_state_all_coords.txt").unlink()

    states = read_stable_states(root)

    assert states == ((1.0, 2.0, 3.0, 4.0), (-1.0, 0.5, 0.0, 2.0))
    assert len(states) == 2, "two states, not eight coordinates"
    assert len(states[0]) == 4, "four dimensions, not four separate states"


def test_per_state_files_are_ordered_numerically(tmp_path):
    """`stable_state_10` must not sort before `stable_state_2`."""
    root = tmp_path / "job"
    root.mkdir()
    for i in (0, 2, 10):
        (root / f"stable_state_{i}_coords.txt").write_text(f"{float(i)}\n")

    assert read_stable_states(root) == ((0.0,), (2.0,), (10.0,))


def test_a_job_with_no_stable_states_is_not_an_error(density_job):
    job = read_density_job(density_job(states=()))

    assert job.stable_states == ()
    assert job.n_states == 0
    assert "estimate_stable_states" in job.report()


def test_a_missing_directory_is_reported_rather_than_raising(tmp_path):
    job = read_density_job(tmp_path / "nope")

    assert job.stable_states == ()
    assert any("not a directory" in n for n in job.notes)


# ── the pickles stay shut ────────────────────────────────────────────────────


def test_density_pickles_are_found_but_never_opened(density_job):
    """The file in the fixture is deliberate junk. Opening it would raise."""
    job = read_density_job(density_job())

    assert [p.name for p in job.density_files] == ["deconv_density_knee.pkl"]
    assert job.deconvolution_was_run is True
    assert "NOT READ" in job.report()
    assert "unpickling runs whatever code" in job.report()


def test_the_report_says_how_to_export_with_the_real_key_names(density_job):
    """Key names are read off RECOVAR's source; a guess would fail on first use."""
    text = read_density_job(density_job()).report()

    assert "d['density']" in text
    assert "d['latent_space_bounds']" in text
    assert "deconv_density_knee.pkl" in text


def test_export_instructions_are_available_without_a_job():
    assert "read_deconvolved_weights" in export_instructions()


# ── peak height is not occupancy ─────────────────────────────────────────────


def test_no_occupancy_is_reported_and_the_report_says_why(density_job):
    """The subtlety that makes the obvious implementation wrong.

    Occupancy is the density integrated over a basin, not its peak value. A
    reader that returned peak heights as weights would look right and be wrong.
    """
    text = read_density_job(density_job()).report()

    assert "NO OCCUPANCY IS REPORTED HERE" in text
    assert "integrated over" in text
    assert "sharp narrow peak" in text


# ── reading exported weights ─────────────────────────────────────────────────


def test_weights_from_npy_are_labelled_deconvolved(tmp_path):
    pops = read_deconvolved_weights(write_npy(tmp_path / "w.npy", [3.0, 1.0]))

    assert pops.source is WeightSource.DECONVOLVED
    assert pops.is_quantitative is True
    assert pops.probabilities == pytest.approx((0.75, 0.25))


def test_weights_from_text_are_read_too(tmp_path):
    path = tmp_path / "w.txt"
    path.write_text("# exported\n3.0\n1.0\n")

    assert read_deconvolved_weights(path).probabilities == pytest.approx((0.75, 0.25))


def test_free_energy_survives_the_round_trip(tmp_path):
    """The whole point of the label: these weights may carry a number."""
    pops = read_deconvolved_weights(
        write_npy(tmp_path / "w.npy", [0.75, 0.25]), temperature_k=298.15
    )

    assert pops.relative_free_energy()[1][0] > 0


def test_a_2d_npy_is_refused_rather_than_flattened(tmp_path):
    """A density grid is not a weight vector, and flattening one would produce
    a plausible-looking distribution over nothing."""
    path = tmp_path / "grid.npy"
    header = "{'descr': '<f8', 'fortran_order': False, 'shape': (2, 2), }"
    header += " " * ((64 - (10 + len(header)) % 64) % 64) + "\n"
    path.write_bytes(
        b"\x93NUMPY\x01\x00"
        + struct.pack("<H", len(header))
        + header.encode("latin1")
        + struct.pack("<4d", 1.0, 2.0, 3.0, 4.0)
    )

    with pytest.raises(ValueError, match="one-dimensional"):
        read_deconvolved_weights(path)


def test_a_multi_column_text_file_is_refused(tmp_path):
    path = tmp_path / "grid.txt"
    path.write_text("1.0 2.0\n3.0 4.0\n")

    with pytest.raises(ValueError, match="one weight per state"):
        read_deconvolved_weights(path)


def test_an_object_array_is_refused_rather_than_unpickled(tmp_path):
    path = tmp_path / "obj.npy"
    header = "{'descr': '|O', 'fortran_order': False, 'shape': (2,), }"
    header += " " * ((64 - (10 + len(header)) % 64) % 64) + "\n"
    path.write_bytes(
        b"\x93NUMPY\x01\x00" + struct.pack("<H", len(header)) + header.encode("latin1")
    )

    with pytest.raises(ValueError, match="unpickling runs whatever code"):
        read_deconvolved_weights(path)


def test_an_empty_file_is_refused(tmp_path):
    path = tmp_path / "empty.txt"
    path.write_text("# nothing but a comment\n")

    with pytest.raises(ValueError, match="no weights found"):
        read_deconvolved_weights(path)
