# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- The package, extracted from MCPymol's `mcpymol.wiggles`. Fourteen views
  across three tiers: occupancy and ensembles, maps and local resolution, and
  heterogeneity-method output.
- No runtime dependencies. Two hosts depend on this package, so a pin here
  would become a constraint in both of them.
- `wiggles-map-info`, a CLI over `map_info` — `mapinfo.py` was already written
  with a `main()`, and reading a real deposition's header is the quickest check
  that a map is what somebody says it is.

- `wiggles_em.scene` — the viewer-neutral seam. Twelve ops, a five-shape
  selection algebra, and `ScalarField`. Views return `(report, Scene)` and
  never call a viewer.
- `wiggles_em.backends` — `PymolBackend` and `FakeBackend`, the only code in
  the package that knows a viewer exists. A backend refuses an op it cannot
  honour rather than skipping it: a dropped op leaves a picture that looks fine
  and means something else.

### Changed

- **The package no longer hosts its own tools.** `tools.py` held MCPymol's
  `@mcp.tool` wrappers and did not come across; the views are the public
  surface and a host wraps them. `__init__` re-exports each view from the
  module it lives in rather than from the tool layer.

- **Views take data, not a port.** `occupancy_view(atoms, obj)` rather than
  `occupancy_view(port, obj)`. Reading atoms is the host's job — in MCPymol
  that is an `iterate_to_list` round trip, in protean the atoms are already in
  Python because biotite parsed the file.
- **Loaders are not views.** `load_map` and `load_ensemble` issue a load and
  then query back to confirm the object arrived, which a one-shot `Scene`
  cannot express, so they keep taking a port.
- **Unit conversion moved into the backend.** Views state levels in the unit
  the data is in and tag it; `PymolBackend` converts to σ against that map's
  own header. Mol\* takes absolute levels natively and will not convert at all.
- **The B-factor stash moved into `PymolBackend`**, along with the note about
  it. Both viewers route per-atom scalars through the B-factor column, but
  PyMOL has one copy of an object and must stash the originals while protean
  re-sends a display copy — so the "call `restore_bfactors`" advice is PyMOL's,
  not occupancy's.
- `local_resolution_view` lost `ramp_name`. A ramp is a PyMOL object; no other
  viewer has one to name.

### Notes

Between 2026-08-09 and this split the code lived inside MCPymol as
`mcpymol.wiggles` (PR #54, tiers 1–2; PR #55, tier 3). It moved out when a
second viewer — [protean](https://github.com/chemrich/protean), driving Mol\*
in a browser — made "a submodule of one viewer" the wrong shape. No released
MCPymol ever carried it: `v1.5.1` was cut before PR #54 landed, so nothing was
removed from anybody's install.
