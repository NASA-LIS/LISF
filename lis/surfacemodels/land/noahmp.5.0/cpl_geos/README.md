# GEOS coupling (`geos_cpl`) run mode — wiring completion notes

This records the files touched to finish wiring the `GEOS coupling` run mode
for NoahMP.5.0, cross-checked against `reference/LIS_FieldCatalog_RAR.F90`
(the reference-only catalog of required import/export variables).

## Files edited

- `lis/plugins/LIS_lsmcpl_pluginMod.F90` — widened the outer `#if ( defined
  RM_WRF_COUPLING )` guard around the whole `LIS_lsmcpl_plugin` body to
  `#if ( defined RM_WRF_COUPLING ) || ( defined RM_GEOS_COUPLING )`. Without
  this, a GEOS-only build (no WRF coupling enabled) would never register
  `NoahMP50_setgeosexport`.
- `lis/surfacemodels/land/noahmp.5.0/cpl_geos/NoahMP50_geosCplMod.F90` —
  made `NoahMP50_geoscpl_init` idempotent (skips allocation if the nest's
  `tile(:)` store already exists), and added two new derived-variable
  members: `canopint`, `tair_t2m`.
- `lis/surfacemodels/land/noahmp.5.0/cpl_geos/noahMP50_convert_to_GEOSexport.F90`
  — added the computation of `canopint` (`canliq + canice`) and `tair_t2m`
  (fveg-weighted `t2mv`/`t2mb`).
- `lis/surfacemodels/land/noahmp.5.0/cpl_geos/noahMP50_setgeosexport.F90` —
  this previously read `NoahMP50_geoscpl(n)%tile(:)` without ever allocating
  or populating it (would have failed at runtime). Now calls
  `NoahMP50_geoscpl_init(n)` + `NoahMP50_convert_to_GEOSexport(n)` at entry.
  Also added population of: `avgsurft_t` (from `tg`), `emiss_t` (from
  `emiss`), `qle_t` (from `lh`) — these export slots existed but were never
  written — plus newly-added exports `tair_t2m_t`, `canopint_t`, `evap_t`
  (from `qfx`), `greenness_t` (from `fveg`), `landcover_t` (raw `vegetype`
  passthrough — see note below), `nee_t` (from `nee`), `avgsurfq_t` (from
  `qsfcmean`). Every addition is traceable to an `adExport=.true.` entry in
  the field catalog.
- `lis/runmodes/geos_cpl_mode/LISGEOSexport_module.F90` — added 2D/1D-tile
  member pairs for the 7 newly-wired fields: `evap`, `greenness`,
  `landcover`, `nee`, `avgsurfq`, `canopint`, `tair_t2m`.
- `lis/runmodes/geos_cpl_mode/LISGEOSGridCompMod.F90` — added matching
  `allocate(...)` (alloc_states) and `= LIS_rc%udef` (reset_states) lines
  for the same 7 fields.
- `lis/core/LIS_surfaceModelMod.F90` — in `surfaceModel_setexport_noesmf`'s
  `#if (defined RM_GEOS_COUPLING)` block, added `LIS_tile2grid` calls for
  `avgsurft`, `emiss`, `qle` (previously allocated but never gridded) and
  for the 7 new fields above.
- `lis/configs/lis.config.adoc` — added a `"GEOS coupling"` row to the
  `Running mode:` options table (was missing; `"WRF coupling"` was already
  documented).

## Files created

- `lis/surfacemodels/land/noahmp.5.0/cpl_geos/README.md` — this file.

## Files moved

- `lis/surfacemodels/land/noahmp.5.0/cpl_geos/LIS_FieldCatalog_RAR.F90` →
  `lis/surfacemodels/land/noahmp.5.0/cpl_geos/reference/LIS_FieldCatalog_RAR.F90`.
  Confirmed (via `make/plugins.py` and `make/Makefile`) that the build's
  source-directory list is non-recursive and exact-path, so this
  subdirectory is never swept into compilation — the file stays
  documentation-only without touching `make/default.cfg`. Added a short
  implementation-status note at the top of the file.

## Files removed

- `lis/runmodes/geos_cpl_mode/.DS_Store` — stray macOS metadata file that
  had been accidentally committed.

## Not changed (flagged for awareness)

- **Submodule pointer**: `surfacemodels/land/noahmp.5.0/phys/noahmp` is
  already modified-but-uncommitted in the working tree, pinned to the
  `feature/add-energy-vars-for-geos` branch that provides the
  `ALBSFCDIRXY`/`ALBSFCDIFXY` fields `NoahMP50_main.F90` depends on.
  Left as-is — committing a submodule pointer bump is a separate,
  cross-repo decision.
- **`MODEL_OUTPUT_LIST.TBL.adoc`**: the four albedo output rows added in an
  earlier commit (`AlbDirVis`, `AlbDirNIR`, `AlbDifVis`, `AlbDifNIR`) all
  reuse numeric token `084`. Pre-existing, unrelated to run-mode wiring;
  not changed here.
- **`lis/runmodes/0Intro_runmodes.txt`**: reviewed but not edited — it's a
  generic description paragraph that doesn't enumerate individual run modes
  by name (not even `WRF coupling` is listed there), so a GEOS-specific
  addition would be inconsistent with its existing style.
- **`Landcover` export**: exported as a raw NoahMP `vegetype` passthrough
  (decided with the requester) since no NoahMP→CLSM land-cover category
  remap table exists yet; flagged with a code comment in
  `noahMP50_setgeosexport.F90`.
