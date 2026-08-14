# wiggles-em

[![CI](https://github.com/chemrich/wiggles-em/actions/workflows/ci.yml/badge.svg)](https://github.com/chemrich/wiggles-em/actions/workflows/ci.yml)

Cryo-EM views for the things a structure viewer usually throws away — partial
occupancy, alternate conformations, ensemble spread, published Q-scores, map
geometry, the per-voxel resolution field, and the latent spaces of
heterogeneity methods.

Each view exists because a plausible-looking render was hiding something.
They share one rule:

> **Say what is being shown, and refuse when the picture would be
> meaningful-looking and wrong.**

`morph_states` declines to interpolate states that do not share a topology,
because a morph across independently reconstructed volumes animates a
correspondence nobody established. `local_resolution_view` declines to colour
one map by another that does not share its voxel grid, because sampling a
resolution field at the wrong coordinates does not render visibly broken — it
renders smooth and plausible. `load_map` records provenance and never infers
it: a measured reconstruction and a network-enhanced one are the same
isosurface once drawn.

## Viewer-neutral

This package hosts nothing. No MCP server, no viewer dependency, and **no
runtime dependencies at all** — the MRC reader is `struct`, the `.npy` reader
parses the header rather than importing numpy, and the ensemble maths is a
mean and a root-mean-square.

A host supplies a port and wraps the views:

- [MCPymol](https://github.com/chemrich/MCPymol) drives PyMOL.
- [protean](https://github.com/chemrich/protean) drives Mol\* in a browser.

```python
from wiggles_em import occupancy_view
from wiggles_em.port import BridgePort

print(occupancy_view(BridgePort(), "6xyz"))
```

Every view takes a port as its first argument, which is right for testing and
wrong for a tool a model calls — so hosts wrap them. That boundary is the whole
integration surface.

## The views

| Tier | Views |
|---|---|
| Occupancy and ensembles | `occupancy_view`, `altloc_view`, `ensemble_spread_view`, `morph_states`, `qscore_view` |
| Maps | `map_info`, `load_map`, `density_view`, `local_resolution_view` |
| Heterogeneity | `load_ensemble`, `latent_traverse_view`, `deformation_view`, `composition_view` |

Plus `restore_bfactors`, because several views push a scalar through the
B-factor column and saying so is not the same as undoing it.

### Two things named "occupancy"

`occupancy_view` reads per-atom crystallographic `q`. `composition_view` reads
the fraction of imaged *particles* containing a subunit. A model can be
`q = 1.0` everywhere while its subunit is present in half the particles — both
true, different questions, and a render that conflates them looks entirely
normal. They stay separate, neither infers the other, and every legend names
its sense.

## Invariants

- **I1** — provenance defaults to `UNKNOWN` and is never inferred. The one
  deliberate exception is `load_ensemble` declaring `GENERATED`, which asserts
  the conservative direction.
- **I2** — no unlabelled latent plot. `latent_traverse_view` draws no latent
  scatter and estimates no density, because motion is recoverable from a
  heterogeneity method and populations are not.
- **I3** — a gap is not an absence. Missing density means the method did not
  resolve something there; it does not mean nothing is there.

## Install

```bash
pip install wiggles-em
```

## Development

```bash
uv sync
uv run pytest          # the whole suite, with no viewer installed
uv run ruff check
uv run mypy
```

CI runs all three: `pytest` on Python 3.10, 3.11, 3.12 and 3.13, and `ruff` and
`mypy` once each — both are configured to target 3.10 regardless of the
interpreter they run on, so a matrix would repeat one answer four times.

Two suites stay out of CI and have to be run by hand:

```bash
uv run pytest -m live           # drives a running PyMOL; CLEARS ITS SESSION
WIGGLES_LIVE=1 uv run pytest    # downloads real maps from EMDB
```

`FakePort` records commands and replays canned query results, so every view is
testable without a viewer. It raises on an unstubbed *data* query rather than
returning `None`, so a view asking for something the test did not set up fails
loudly instead of silently rendering nothing.

`map_info` also runs as a CLI, which is the quickest way to check a real
deposition:

```bash
uv run wiggles-map-info emd_30913.map.gz
```

## Where this came from

The research it was distilled from — a cited compendium of cryo-EM
heterogeneity methods, and the specification the views reconcile to — lives in
a separate, private repository. Docstrings here name the entry an argument came
from (`occupancy-two-senses`, `local-resolution`) rather than a path that would
not resolve.

The code lived inside MCPymol as `mcpymol.wiggles` between 2026-08-09 and the
split; it moved out when a second viewer turned up.

## Licence

MIT. See [LICENSE](LICENSE).
