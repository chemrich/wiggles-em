"""Wiggles: cryo-EM heterogeneity, occupancy and resolution views.

Fourteen views for the things a structure viewer usually throws away — partial
occupancy, alternate conformations, ensemble spread, published Q-scores, map
geometry, the per-voxel resolution field, and the latent spaces of
heterogeneity methods. Each one exists because a plausible-looking render was
hiding something a user needed to know.

The design rule the whole package follows: **say what is being shown, and
refuse when the picture would be meaningful-looking and wrong.** ``morph_states``
declines to interpolate states that do not share a topology.
``local_resolution_view`` declines to colour one map by another that does not
share its grid. ``load_map`` records provenance and never guesses it, because a
measured reconstruction and a neural-network hallucination are the same
isosurface once they are drawn.

This package is **viewer-neutral and hosts nothing.** It has no MCP tools, no
server and no viewer dependency; a host wraps these functions and supplies a
port. MCPymol drives PyMOL, protean drives Mol\\* in a browser.

Layout:

======================  ====================================================
``port``                the coupling surface to a viewer
``atoms``, ``bfactors`` shared read layer and the B-factor stash
``occupancy``           ``occupancy_view``, ``altloc_view``
``ensembles``           ``ensemble_spread_view``, ``morph_states``
``qscore``              validation-report parsing and ``qscore_view``
``mapinfo``             MRC/CCP4 header reading; no viewer, no network
``maps``, ``density``   loading volumes, and contouring them honestly
``localres``            ``local_resolution_view``
``heterogeneity``       ensembles and method detection; ``load_ensemble``
``latent``              ``latent_traverse_view`` and invariant I2
``deformation``         ``deformation_view``
``composition``         ``composition_view`` — occupancy in the *other* sense
``provenance``          where a volume came from; the default is UNKNOWN
======================  ====================================================

Everything above ``port`` is pure: the whole suite runs against ``FakePort``
with no viewer installed at all.

The research this was distilled from — a cited compendium of cryo-EM
heterogeneity methods — lives in a separate repository. Docstrings here name
the entry they came from (`occupancy-two-senses`, `local-resolution`) rather
than a path that would not resolve.
"""

from wiggles_em.atoms import Atom, fetch_atoms
from wiggles_em.bfactors import clear_stash, restore_bfactors
from wiggles_em.composition import composition_view
from wiggles_em.deformation import deformation_view
from wiggles_em.density import density_view, to_absolute, to_sigma
from wiggles_em.ensembles import ensemble_spread_view, morph_states
from wiggles_em.heterogeneity import Ensemble, Method, load_ensemble, loaded_ensemble
from wiggles_em.latent import contains_absence_claim, latent_traverse_view
from wiggles_em.localres import grid_differences, local_resolution_view
from wiggles_em.mapinfo import MapHeader, map_info, read_map_header
from wiggles_em.maps import load_map
from wiggles_em.occupancy import altloc_view, occupancy_view
from wiggles_em.port import BridgePort, FakePort, PortError, PymolPort, SendRequestPort
from wiggles_em.provenance import Provenance, declare, provenance_of
from wiggles_em.qscore import qscore_view

__all__ = [  # noqa: RUF022
    # the fourteen views, grouped by tier rather than alphabetically
    "occupancy_view",
    "altloc_view",
    "ensemble_spread_view",
    "morph_states",
    "qscore_view",
    "restore_bfactors",
    "map_info",
    "load_map",
    "density_view",
    "local_resolution_view",
    # tier 3 — ensembles
    "load_ensemble",
    "latent_traverse_view",
    "deformation_view",
    "composition_view",
    # the pieces worth reaching for directly
    "MapHeader",
    "read_map_header",
    "to_sigma",
    "to_absolute",
    "grid_differences",
    "Ensemble",
    "Method",
    "loaded_ensemble",
    "contains_absence_claim",
    "Provenance",
    "declare",
    "provenance_of",
    "Atom",
    "fetch_atoms",
    "clear_stash",
    "PymolPort",
    "FakePort",
    "BridgePort",
    "SendRequestPort",
    "PortError",
]
