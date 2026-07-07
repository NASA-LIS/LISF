#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.8
#
# Copyright (c) 2026 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
===============================================================================
Dynamic Lapse Rate Calculation Script
===============================================================================

Compute dynamic near-surface air-temperature lapse rates from gridded
temperature, elevation, and land-mask inputs using a local linear regression
within a 3 x 3 neighborhood.

Summary
-------
This script computes lapse rates on the target temperature grid by fitting the
local relationship between temperature differences (dT) and elevation
differences (dZ) around each land grid cell. The method follows the framework
of Rouf et al. (2020) and the configuration used in Whitney et al. (2025).

For each land grid cell and time step, the script:
1. Builds a 3 x 3 neighborhood around the target cell.
2. Screens the neighborhood for valid elevation data and valid land neighbors.
3. Computes a local regression slope dT/dZ.
4. Optionally applies user-specified lapse-rate cutoffs.
5. Converts units if requested.
6. Imprints a fallback lapse rate over shoreline-buffer and inland-water cells.

For NLDAS-3 production workflows, lapse-rate magnitude constraints are typically
applied downstream within the LIS meteorological downscaling workflow rather
than during dynamic lapse-rate generation.

The script supports:
- MERRA-2 daily files containing multiple hourly time steps
- GEOS-IT individual hourly files
- GEOS-IT daily multi-time-step files
- Optional clipping to a reference-grid bounding box
- Optional nearest-neighbor resampling of elevation and land mask onto the
  target temperature grid when grids differ but overlap spatially
- Daily, hourly, or match-input output modes
- Skip-existing, overwrite, and dry-run modes

Scientific background
---------------------
This implementation follows the dynamic lapse-rate framework introduced by:

Rouf, T., Y. Mei, V. Maggioni, P. Houser, and M. Noonan, 2020:
A physically based atmospheric variables downscaling technique.
Journal of Hydrometeorology, 21(1), 93-108.
https://doi.org/10.1175/JHM-D-19-0109.1

and the configuration used in:

Whitney, K. M., S. V. Kumar, D. Mocko, J. Pflug, J. D. Bolten,
F. Z. Maina, M. L. Wrzesien, and C. R. Hain, 2025:
Quantifying the impacts of dynamic lapse regimes on snow simulations.
Journal of Hydrometeorology, 26(10), 1525-1560.
https://doi.org/10.1175/JHM-D-25-0021.1

Whitney et al. (2025) improved robustness and physical realism relative to the
foundational method by:
- requiring a minimum number of valid surrounding land neighbors before a
  dynamic estimate is attempted, and
- applying additional lapse-rate quality-control procedures during the
  NLDAS-3 forcing-generation workflow.

The dynamic lapse-rate generation script implemented here applies the
minimum-neighbor requirement and fallback lapse-rate behavior. These choices
reduce unstable or nonphysical lapse rates in areas with sparse valid
neighbors, weak terrain contrast, coastlines, and other challenging settings.

Citation
--------
If you use this code, the resulting dynamic lapse-rate fields, or NLDAS-3 forcing datasets
generated using these lapse rates, please cite:

Maina, F., S. V. Kumar, C. Hain, D. Mocko, K. Locke, R. Wade,
A. Getirana, J. Case, K. M. Whitney, S. Ahmad, T. Lahmers,
J. Bolten, J. Pflug, R. Junod, M. Dodson, M. Wrzesien,
V. Mishra, W. Nie, and M. G. Z. Hashemi:
Retrospective and Near–Real-Time Surface Meteorology in the North
American Land Data Assimilation System Phase 3 (NLDAS-3).
Journal of Hydrometeorology, in review.

The dynamic lapse-rate methodology and configuration should also be cited as:

Whitney, K. M., S. V. Kumar, D. Mocko, J. Pflug, J. D. Bolten,
F. Z. Maina, M. L. Wrzesien, and C. R. Hain, 2025:
Quantifying the impacts of dynamic lapse regimes on snow simulations.
Journal of Hydrometeorology, 26(10), 1525–1560.
https://doi.org/10.1175/JHM-D-25-0021.1


Method overview
---------------
At each land grid cell, the script computes a local lapse rate from a 3 x 3
window centered on that cell.

Internal computation units:
- Regression slopes are computed in K/m.

Neighborhood behavior:
- A full 3 x 3 neighborhood is required.
- The outermost grid cells are skipped because a full window cannot be centered
  there.
- If any elevation value in the 3 x 3 window is NaN, that target cell is
  skipped.
- Land-neighbor screening uses the 8 surrounding cells only.
- The regression itself uses all valid dT and dZ pairs in the 3 x 3 window.

Fallback behavior:
- If the number of valid surrounding land neighbors is below
  ``--min_land_neighbors``, the script assigns the fallback lapse rate
  ``--LR_forced_val``.
- If fewer than 2 valid regression samples remain, the script assigns the
  fallback lapse rate.
- If the regression slope is undefined (for example, zero effective elevation
  variance), the output remains NaN at that land cell unless later overwritten
  by water imprinting.

Water handling:
- Water cells are not dynamically regressed.
- Optionally, a fallback lapse rate is imprinted over:
  * shoreline-buffer water cells within ``--buffer_width`` source pixels of land
  * inland water bodies not connected to the domain edge

Cutoffs:
- ``--min_lapse_rate_cutoff`` and ``--max_lapse_rate_cutoff`` are interpreted
  in the requested output units:
  * K/km when ``--save_units_as_per_km 1``
  * K/m when ``--save_units_as_per_km 0``

Target grid and clipping
------------------------
The target output grid is the grid of the temperature files.

Elevation and land-mask inputs must spatially cover that grid. If their grid
differs but overlaps spatially, the script resamples elevation and land mask to
the temperature grid using nearest-neighbor indexing.

If ``--clip_to_reference_bbox`` is used, the script clips by the bounding box
of a reference NetCDF file before processing. This is a bounding-box operation,
not a polygon mask. The reference file defines only the clip region; the output
remains on the clipped target temperature grid.

Input preprocessing
-------------------
For the NLDAS-3 workflows, the elevation and land-mask inputs are derived from
the corresponding MERRA-2 and GEOS-IT constant parameter files:

- MERRA-2: MERRA2.const_2d_asm_Nx.00000000.nc4
- GEOS-IT: GEOS.fpit.asm.const_2d_asm_Nx.GEOS5124.00000000_0000.V01.nc4

Elevation is derived from PHIS by dividing by 9.80665 to convert surface
geopotential to metres. The binary LANDMASK field is derived from the
fractional land and land-ice fields as:

    LANDMASK_FRAC = FRLAND + FRLANDICE
    LANDMASK = 1 where LANDMASK_FRAC >= 0.5, otherwise 0

The time dimension is removed after verifying that the fields are constant.

For NLDAS-3 production, the resulting elevation and land-mask files are
typically named:

- merra2.elevation_and_landmask.nc
- geosit.elevation_and_landmask.nc

Input arguments
---------------
Required arguments
~~~~~~~~~~~~~~~~~~
--infile_Z : str
    Path to the elevation NetCDF file.

--infile_landmask : str
    Path to the land-mask NetCDF file.

--inpath_T : str
    Directory containing temperature NetCDF files.

--outpath : str
    Output directory for lapse-rate NetCDF files.

--yyyymmdd_start : str
    Start date in YYYYMMDD format.

--yyyymmdd_end : str
    End date in YYYYMMDD format.

Variable-name arguments
~~~~~~~~~~~~~~~~~~~~~~~
--Z_varin : str, default='ELEVATION'
    Elevation variable name in ``--infile_Z``.

--T_varin : str, default='T2M'
    Temperature variable name in temperature files.

--landmask_varin : str, default='LANDMASK'
    Land-mask variable name in ``--infile_landmask``.

--time_varin : str, default='time'
    Time variable name in temperature files.

Filename and date-format arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--fname_in_pfx_T : str, default=''
    Filename prefix used to identify temperature input files.

--fname_in_sfx_T : str, default='.nc'
    Filename suffix used to identify temperature input files.

--fname_out_pfx : str, default='lapse_rate.'
    Prefix for output filenames.

--fname_out_sfx : str, default='.nc'
    Suffix for output filenames.

--infile_date_format : {'YYYYMMDD', 'YYYY-MM-DD', 'YYYY-MM-DDTHH', 'YYYY-MM-DD-HH'},
default='YYYYMMDD'
    Date token format embedded in input temperature filenames.

--outfile_date_format :
{'YYYYMMDD', 'YYYY-MM-DD', 'YYYY-MM-DDTHH', 'YYYY-MM-DD-HH'},
default='YYYYMMDD'
    Date token format used in output filenames.

--output_mode : {'daily', 'hourly', 'match_input'}, default='daily'
    Output granularity.
    - 'daily': write one file per day with all time steps preserved
    - 'hourly': write one file per time step
    - 'match_input': infer from ``--infile_date_format``

Algorithm arguments
~~~~~~~~~~~~~~~~~~~
--landmask_trueval : float, default=1
    Value in the land-mask field that indicates land.

--min_land_neighbors : int, default=5
    Minimum number of valid surrounding land neighbors required before a
    dynamic estimate is attempted.

--LR_forced_val : float, default=-0.0065
    Static fallback lapse rate in K/m.

--min_lapse_rate_cutoff : float, optional, default=None
    Minimum allowed lapse rate. Interpreted in K/km when
    ``--save_units_as_per_km 1``, otherwise in K/m. If omitted, no minimum
    cutoff is applied.

--max_lapse_rate_cutoff : float, optional, default=None
    Maximum allowed lapse rate. Interpreted in K/km when
    ``--save_units_as_per_km 1``, otherwise in K/m. If omitted, no maximum
    cutoff is applied.

--save_units_as_per_km : {0, 1}, default=1
    Output units flag.
    - 1: write output in K/km
    - 0: write output in K/m

--buffer_width : int, default=1
    Width, in source-grid pixels, of the land-adjacent shoreline water buffer.

--no_fill_water_mask : flag, optional
    If provided, disable shoreline and inland-water imprinting.

--min_elev_diff_std : float, optional, default=None
    Optional elevation-variability threshold (metres). If the standard
    deviation of elevation differences among valid land neighbours is below
    this value, the cell receives the fallback static lapse rate. If omitted,
    no elevation-variability filter is applied.

Pre-existing lapse-rate arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--prefer_preexisting : flag, optional
    If set, the pipeline first looks for a pre-existing lapse-rate file for
    each input date. If found it is clipped, unit-aligned, and imprinted;
    if not found the pipeline falls back to computing from temperature inputs.

--preexisting_lr_dir : str, optional, default=None
    Directory containing pre-existing lapse-rate NetCDF files.

--preexisting_lr_pfx : str, optional, default=None
    Filename prefix for pre-existing lapse-rate files.

--preexisting_lr_sfx : str, optional, default=None
    Filename suffix for pre-existing lapse-rate files.

--preexisting_lr_var : str, default='lapse_rate'
    Variable name for lapse rate inside pre-existing files.

--preexisting_lr_date_format : str, optional, default=None
    Date-token format for pre-existing lapse-rate filenames. Defaults to
    ``--infile_date_format`` if not specified.

--preexisting_time_var : str, default='time'
    Time variable name in pre-existing lapse-rate files.

Reference-bbox clipping arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--clip_to_reference_bbox : flag, optional
    Enable clipping to the bounding box of a reference grid file.

--reference_grid_file : str, optional, default=None
    NetCDF file whose latitude/longitude bounding box defines the clip region.
    Required if ``--clip_to_reference_bbox`` is used.

--reference_lat_var : str, optional, default=None
    Latitude variable name in the reference file. If omitted, common names are
    auto-detected.

--reference_lon_var : str, optional, default=None
    Longitude variable name in the reference file. If omitted, common names are
    auto-detected.

--bbox_pad_pixels : int, default=0
    Expand the clip window by this many source-grid pixels on each side.

Control arguments
~~~~~~~~~~~~~~~~~
--force : flag, optional
    Overwrite existing output files. By default, existing outputs are skipped.

--dryrun : flag, optional
    Report planned outputs without writing files.

--debug : flag, optional
    Enable verbose logging.

--max_files : int, optional, default=None
    Optional cap on the number of temperature input files processed. Useful
    for testing. If omitted, all discovered files are processed.

--start_datetime : str, optional, default=None
    Optional start datetime in 'YYYY-MM-DD HH:MM' format. When provided
    together with ``--end_datetime``, only timesteps whose containing hour
    falls within the window are processed. Minutes are ignored.

--end_datetime : str, optional, default=None
    Optional end datetime in 'YYYY-MM-DD HH:MM' format. See
    ``--start_datetime``.

--time_shift_minutes : int, default=0
    Optional uniform shift in minutes applied to all input timestamps before
    hour filtering and output filename construction. Negative values shift
    earlier (e.g. ``-30`` shifts 06:30 to 06:00).

Optional output metadata arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--output_source_data : str, optional, default=None
    Global attribute describing the source dataset.

--output_source_variable : str, optional, default=None
    Global attribute giving the source variable name.

--output_source_variable_description : str, optional, default=None
    Global attribute describing the source variable.

--output_input_description : str, optional, default=None
    Global attribute describing the inputs used.

--output_method_reference : str, optional, default=None
    Global attribute describing the method reference.

Outputs
-------
Each output NetCDF file contains a variable named ``lapse_rate``.

Output units:
- K/km when ``--save_units_as_per_km 1``
- K/m when ``--save_units_as_per_km 0``

Output shapes:
- Hourly output: 2-D lapse-rate field plus scalar time coordinate
- Daily output: 3-D array with dimensions (time, lat, lon) for rectilinear
  grids, or (time, y, x) for curvilinear grids

Coordinate conventions:
- 1-D rectilinear lat/lon are preserved as dimension coordinates
- 2-D curvilinear lat/lon are preserved as auxiliary coordinate variables

Global attributes written to every output:
- title = 'Dynamic lapse-rate output'
- source = 'Generated by dynamic-lapse-rates'
- history = UTC creation timestamp

Optional metadata attributes are written when corresponding CLI arguments are
provided.

Typical usage patterns
----------------------

1) Minimal run using defaults
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Useful when input files already match the intended target domain.

python compute_dynamic_lapse_rates.py \
  --infile_Z ./data/MERRA2/merra2.elevation_and_landmask.nc \
  --infile_landmask ./data/MERRA2/merra2.elevation_and_landmask.nc \
  --inpath_T ./data/MERRA2/temperature/ \
  --outpath ./outputs/lapse_rates/ \
  --yyyymmdd_start 20250701 \
  --yyyymmdd_end 20250731 \
  --fname_in_pfx_T MERRA2.tavg1_2d_slv_Nx. \
  --fname_in_sfx_T .nc4

2) NLDAS-3 retrospective (2001-2024), MERRA-2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Each MERRA-2 input file contains multiple hourly time steps. One lapse-rate
output file is written per day, while preserving all hourly time steps present
in the source file.

The example below includes options for clipping the MERRA-2 domain to the
NLDAS-3 modeling domain using ``--clip_to_reference_bbox``. The clipping
extent is defined using the LIS input file
``lis_input.nldas3.noahmp401.1km.land_and_inland_water.m2interp.nc``,
which is available from the NLDAS-3 S3 bucket.

For additional information on accessing NLDAS-3 data products, see:
https://github.com/NASAWaterInsight/NLDAS-3

python compute_dynamic_lapse_rates.py \
  --infile_Z ./merra2.elevation_and_landmask.nc \
  --infile_landmask ./merra2.elevation_and_landmask.nc \
  --inpath_T /discover/nobackup/projects/eis_nldas3/kwhitney/RUN/NLDAS3/rouf_dynamic_lapse_rate/global_domain/data/MERRA2/temperature_hourly/ \
  --outpath ./outputs/lapse_rates_merra2/ \
  --yyyymmdd_start 20001231 \
  --yyyymmdd_end 20241231 \
  --infile_date_format YYYYMMDD \
  --outfile_date_format YYYYMMDD \
  --output_mode daily \
  --fname_in_pfx_T MERRA2.tavg1_2d_slv_Nx. \
  --fname_in_sfx_T .nc4 \
  --fname_out_pfx MERRA2.lapse_rate.hourly. \
  --fname_out_sfx .nldas3.nc \
  --landmask_trueval 1 \
  --min_land_neighbors 5 \
  --LR_forced_val -0.0065 \
  --save_units_as_per_km 1 \
  --clip_to_reference_bbox \
  --reference_grid_file ./lis_input.nldas3.noahmp401.1km.land_and_inland_water.m2interp.nc

3) NLDAS-3 near-real-time, GEOS-ITbias
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Computes dynamic lapse rates from GEOS-ITbias 2-m air-temperature inputs.

GEOS-ITbias inputs are generated by bias-correcting GEOS-IT using
MERRA-2 monthly climatologies derived from the 2001–2024 period.

The example below demonstrates processing of 24 hourly GEOS-ITbias files
covering a single simulation day. Because the GEOS-ITbias inputs have
already been clipped to the NLDAS-3 modeling domain, the command does not
include the ``--clip_to_reference_bbox`` or ``--reference_grid_file``
options.

python compute_dynamic_lapse_rates.py \
  --infile_Z ./geosit.elevation_and_landmask.nc \
  --infile_landmask ./geosit.elevation_and_landmask.nc \
  --inpath_T ./data/GEOSIT/temperature/ \
  --outpath ./outputs/lapse_rates_geosit/ \
  --yyyymmdd_start 20230430 \
  --yyyymmdd_end 20230501 \
  --start_datetime "2023-04-30 12:00" \
  --end_datetime "2023-05-01 12:00" \
  --Z_varin "ELEVATION" \
  --T_varin "T2M" \
  --landmask_varin "LANDMASK" \
  --time_varin "time" \
  --infile_date_format YYYY-MM-DDTHH \
  --outfile_date_format YYYY-MM-DDTHH \
  --fname_in_pfx_T "T2M_" \
  --fname_in_sfx_T ".nc" \
  --fname_out_pfx "GEOS-ITbias.lapse_rate.hourly." \
  --fname_out_sfx ".nldas3.nc" \
  --output_mode "hourly" \
  --landmask_trueval 1 \
  --min_land_neighbors 5 \
  --LR_forced_val -0.0065 \
  --save_units_as_per_km 1

4) Overwrite or preview outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Overwrite existing files:
python compute_dynamic_lapse_rates.py ... --force

Preview planned outputs only:
python compute_dynamic_lapse_rates.py ... --dryrun

Notes
-----
- Existing outputs are skipped unless ``--force`` is provided.
- Daily outputs preserve all time steps; no temporal averaging is applied.
- Bounding-box clipping and target-grid resampling are separate steps:
  clipping limits the domain, while resampling aligns elevation and land-mask
  fields to the temperature grid when needed.
- For MERRA-2 and GEOS-IT datasets already on a common grid, no special
  preprocessing is usually required beyond correct variable names and filename
  settings.
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import re
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Iterable
from collections import defaultdict
import gc

import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
import xarray as xr
from scipy.spatial import cKDTree


# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------

def configure_logging(debug: bool = False) -> None:
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        force=True,
    )


LOGGER = logging.getLogger("dynamic_lapse_rates_single")


# -----------------------------------------------------------------------------
# Data containers
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class FileSpec:
    prefix: str
    suffix: str
    date_format: str


@dataclass(frozen=True)
class ClipWindow:
    j0: int
    j1: int
    i0: int
    i1: int


@dataclass(frozen=True)
class Config:
    infile_Z: str
    infile_landmask: str
    inpath_T: str
    outpath: str
    yyyymmdd_start: str
    yyyymmdd_end: str
    start_datetime: str | None
    end_datetime: str | None
    time_shift_minutes: int

    Z_varin: str
    T_varin: str
    landmask_varin: str
    time_varin: str

    fname_in_pfx_T: str
    fname_in_sfx_T: str
    fname_out_pfx: str
    fname_out_sfx: str
    infile_date_format: str
    outfile_date_format: str
    output_mode: str

    landmask_trueval: float
    min_land_neighbors: int
    LR_forced_val: float
    min_lapse_rate_cutoff: float | None
    max_lapse_rate_cutoff: float | None
    save_units_as_per_km: int
    buffer_width: int
    no_fill_water_mask: bool
    min_elev_diff_std: float | None

    prefer_preexisting: bool
    preexisting_lr_dir: str | None
    preexisting_lr_pfx: str | None
    preexisting_lr_sfx: str | None
    preexisting_lr_var: str
    preexisting_lr_date_format: str | None
    preexisting_time_var: str

    max_files: int | None

    clip_to_reference_bbox: bool
    reference_grid_file: str | None
    reference_lat_var: str | None
    reference_lon_var: str | None
    bbox_pad_pixels: int

    force: bool
    dryrun: bool

    output_source_data: str | None
    output_source_variable: str | None
    output_source_variable_description: str | None
    output_input_description: str | None
    output_method_reference: str | None


# -----------------------------------------------------------------------------
# Filename and time helpers
# -----------------------------------------------------------------------------

def parse_date_token(filename: str, fmt: str) -> datetime:
    """Parse a date token from a filename according to the specified format.
    Parameters
    ----------
    filename : str
        The filename to parse.
    fmt : str
        The format string for the date token.

    Returns
    -------
    datetime
        The parsed datetime object.

    Raises
    ------
    ValueError
        If the date token cannot be parsed according to the specified format.
    Notes
    -----
    This function supports the following date formats:
    - "YYYYMMDD": Parses an 8-digit date string (e.g., "20230101").
    - "YYYY-MM-DD": Parses a date string in YYYY-MM-DD format (e.g., "2023-01-01").
    - "YYYY-MM-DDTHH": Parses a date-time string in YYYY-MM-DDTHH format (e.g., "2023-01-01T12").
    - "YYYY-MM-DD-HH": Parses a date-time string in YYYY-MM-DD-HH format (e.g., "2023-01-01-12").
    """
    if fmt == "YYYYMMDD":
        match = re.search(r"(\d{8})", filename)
        if not match:
            raise ValueError(f"Could not parse YYYYMMDD date token from {filename!r}.")
        return datetime.strptime(match.group(1), "%Y%m%d")

    if fmt == "YYYY-MM-DD":
        match = re.search(r"(\d{4}-\d{2}-\d{2})(?!T)", filename)
        if not match:
            raise ValueError(f"Could not parse YYYY-MM-DD date token from {filename!r}.")
        return datetime.strptime(match.group(1), "%Y-%m-%d")

    if fmt == "YYYY-MM-DDTHH":
        match = re.search(r"(\d{4}-\d{2}-\d{2})T(\d{2})(\d{2})?", filename)
        if not match:
            raise ValueError(f"Could not parse YYYY-MM-DDTHH date token from {filename!r}.")
        hh = match.group(2)
        mm = match.group(3) or "00"
        return datetime.strptime(f"{match.group(1)}T{hh}:{mm}", "%Y-%m-%dT%H:%M")
    
    if fmt == "YYYY-MM-DD-HH":
        match = re.search(r"(\d{4}-\d{2}-\d{2})-(\d{2})", filename)
        if not match:
            raise ValueError(f"Could not parse YYYY-MM-DD-HH date token from {filename!r}.")
        return datetime.strptime(f"{match.group(1)}-{match.group(2)}", "%Y-%m-%d-%H")
    
    raise ValueError(f"Unsupported date format: {fmt}")


def build_dated_filename(spec: FileSpec, dt: datetime) -> str:
    """Construct a filename by embedding a date token according to the specified format.
    Parameters
    ----------
    spec : FileSpec
        The file specification containing prefix, suffix, and date format.
        dt : datetime
        The datetime object to embed in the filename.
    
    Returns
        -------
        str
            The constructed filename.
     Raises
     ------
     ValueError
        If the date format in the specification is unsupported.
    
    Notes
        This function supports the following date formats:
        - "YYYYMMDD": Embeds the date as an 8-digit string (e.g., "20230101").
        - "YYYY-MM-DD": Embeds the date in YYYY-MM-DD format (e.g., "2023-01-01").
        - "YYYY-MM-DDTHH": Embeds the date-time in YYYY-MM-DDTHH format (e.g., "2023-01-01T12").
        - "YYYY-MM-DD-HH": Embeds the date-time in YYYY-MM-DD-HH format (e.g., "2023-01-01-12").
    """
    if spec.date_format == "YYYYMMDD":
        token = dt.strftime("%Y%m%d")
    elif spec.date_format == "YYYY-MM-DD":
        token = dt.strftime("%Y-%m-%d")
    elif spec.date_format == "YYYY-MM-DDTHH":
        token = dt.strftime("%Y-%m-%dT%H%M")
    elif spec.date_format == "YYYY-MM-DD-HH":
        token = dt.strftime("%Y-%m-%d-%H")
    else:
        raise ValueError(f"Unsupported date format: {spec.date_format}")
    return f"{spec.prefix}{token}{spec.suffix}"


def discover_temperature_files(
    inpath: str | Path,
    spec: FileSpec,
    start_date: datetime,
    end_date: datetime,
) -> list[Path]:
    """Discover temperature files within a specified date range.
    Parameters
    ----------
    inpath : str | Path
        The directory to search for temperature files.
    spec : FileSpec
        The file specification containing prefix, suffix, and date format.
    start_date : datetime
        The start date of the range to search.
    end_date : datetime
        The end date of the range to search.

    Returns
    -------
    list[Path]
        A list of discovered temperature files.
    """

    base = Path(inpath)
    if not base.exists():
        raise FileNotFoundError(f"Temperature input directory does not exist: {base}")

    matches: list[tuple[datetime, Path]] = []
    for path in sorted(base.iterdir()):
        if not path.is_file():
            continue
        if not path.name.startswith(spec.prefix) or not path.name.endswith(spec.suffix):
            continue
        try:
            stamp = parse_date_token(path.name, spec.date_format)
        except ValueError:
            continue
        if start_date.date() <= stamp.date() <= end_date.date():
            matches.append((stamp, path))

    files = [path for _, path in sorted(matches, key=lambda item: item[0])]
    LOGGER.info("Discovered %s temperature files in requested date range.", len(files))
    return files

def apply_time_shift(time_values: list[datetime], shift_minutes: int = 0) -> list[datetime]:
    """Apply a uniform time shift in minutes to all timestamps.

    Parameters
    ----------
    time_values : list[datetime]
        Input timestamps.
    shift_minutes : int, default=0
        Minutes to shift each timestamp. Negative values shift earlier,
        positive values shift later.

    Returns
    -------
    list[datetime]
        Shifted timestamps.
    """
    if shift_minutes == 0:
        return list(time_values)

    delta = timedelta(minutes=shift_minutes)
    LOGGER.info("Applying explicit time shift of %d minutes to input timestamps.", shift_minutes)
    return [t + delta for t in time_values]

# def _normalize_hourly_times(time_values: list[datetime]) -> list[datetime]:
#     """Detect and correct for common non-hourly timestamp offsets in hourly data.
#     This function checks for a consistent offset of approximately 30 minutes in the
#     timestamps, which can occur when hourly data are timestamped at the center of
#     the hour rather than the start. If such an offset is detected, the function shifts
#     all timestamps back by the median offset to align them with whole hours.
    
#     Parameters
#     ----------
#     time_values : list[datetime]
#         A list of datetime objects representing the time values to check and potentially adjust.
    
#     Returns
#     -------
#     list[datetime]
#         A list of datetime objects with timestamps adjusted to whole hours if a consistent
#         offset was detected, or the original timestamps if no such offset was found.
#     """
#     offsets = np.array([t.minute + t.second / 60.0 for t in time_values], dtype=float)
#     median_offset = float(np.median(offsets))
#     if 20.0 <= median_offset <= 40.0:
#         LOGGER.debug(
#             "Detected ~%.1f min timestamp offset; shifting times by -%.0f min.",
#             median_offset,
#             round(median_offset),
#         )
#         delta = timedelta(minutes=round(median_offset))
#         return [t - delta for t in time_values]
#     return list(time_values)


def get_time_values(
    ds: xr.Dataset,
    time_var: str,
    time_shift_minutes: int = 0,
) -> tuple[list[datetime], dict]:
    """Extract time values from a dataset, handling CF time conventions.
    
    Parameters
    ----------
    ds : xr.Dataset
        The dataset from which to extract time values.
    time_var : str
        The name of the time variable to extract.

    Returns
    -------
    tuple[list[datetime], dict]
        A list of datetime objects representing the time values and a dictionary of their attributes.
    """
    if time_var not in ds:
        raise KeyError(f"Time variable {time_var!r} not present in dataset.")

    raw = ds[time_var].values
    attrs = dict(ds[time_var].attrs)

    if np.issubdtype(raw.dtype, np.datetime64):
        decoded = pd.to_datetime(raw)
    
        if np.ndim(decoded) == 0:
            time_list = [pd.Timestamp(decoded).to_pydatetime()]
        else:
            time_list = [pd.Timestamp(t).to_pydatetime() for t in decoded]
    
        return apply_time_shift(time_list, time_shift_minutes), attrs

    units_str = attrs.get("units", "")
    if "since" in units_str:
        try:
            unit_part, ref_part = units_str.split("since", 1)
            unit_part = unit_part.strip().lower()
            ref_part = ref_part.strip()
            ref_time: datetime | None = None
            for fmt in (
                "%Y-%m-%d %H:%M:%S",
                "%Y-%m-%d %H:%M",
                "%Y-%m-%d",
                "%Y-%m-%dT%H:%M:%S",
                "%Y-%m-%dT%H:%M",
                "%Y-%m-%dT%H",
            ):
                try:
                    ref_time = datetime.strptime(ref_part, fmt)
                    break
                except ValueError:
                    continue
            if ref_time is None:
                raise ValueError(f"Could not parse CF reference time: {ref_part!r}")
        except Exception as exc:
            LOGGER.warning(
                "Failed to parse CF time units %r (%s); defaulting to hours since 2010-01-01.",
                units_str,
                exc,
            )
            unit_part = "hours"
            ref_time = datetime(2010, 1, 1)
    else:
        LOGGER.warning(
            "Time variable %r has no recognisable units attribute; assuming hours since 2010-01-01.",
            time_var,
        )
        unit_part = "hours"
        ref_time = datetime(2010, 1, 1)

    numeric = np.atleast_1d(np.asarray(raw, dtype=float))
    if unit_part in ("minutes", "minute"):
        time_list = [ref_time + timedelta(minutes=float(v)) for v in numeric]
    elif unit_part in ("hours", "hour"):
        time_list = [ref_time + timedelta(hours=float(v)) for v in numeric]
    elif unit_part in ("days", "day"):
        time_list = [ref_time + timedelta(days=float(v)) for v in numeric]
    else:
        raise ValueError(f"Unsupported CF time unit: {unit_part!r}")

    if time_shift_minutes == 0:
        offsets = np.array([t.minute + t.second / 60.0 for t in time_list], dtype=float)
        median_offset = float(np.median(offsets))
        if 20.0 <= median_offset <= 40.0:
            delta = timedelta(minutes=round(median_offset))
            LOGGER.debug("Auto-correcting ~%.0f min centred timestamps.", median_offset)
            time_list = [t - delta for t in time_list]

    return apply_time_shift(time_list, time_shift_minutes), attrs


# -----------------------------------------------------------------------------
# Grid helpers
# -----------------------------------------------------------------------------

def normalize_lon180(lon: np.ndarray) -> np.ndarray:
    """Normalize longitude values to the range (-180, 180], handling edge cases.
    Parameters
    ----------
    lon : np.ndarray
        An array of longitude values to normalize.
    
    Returns
    -------
    np.ndarray
        An array of longitude values normalized to the range (-180, 180], with -180.0 replaced by 180.0.
    
    Notes
    -----
    This function normalizes longitude values to the range (-180, 180] using modulo arithmetic.
    It also handles the edge case where values exactly equal to -180.0 are replaced with 180.0 to maintain a consistent representation.
    """
    values = np.asarray(lon, dtype=float).copy()
    values = ((values + 180.0) % 360.0) - 180.0
    values[values == -180.0] = 180.0
    return values


def get_lat_lon(
    ds_like: xr.Dataset | xr.DataArray,
    lat_var: str | None = None,
    lon_var: str | None = None,
) -> tuple[np.ndarray, np.ndarray, str, str]:
    """Extract latitude and longitude arrays from a dataset, with auto-detection of variable names.
    
    Parameters
    ----------
    ds_like : xr.Dataset | xr.DataArray
        The dataset or data array from which to extract latitude and longitude.
    lat_var : str | None, optional
        The name of the latitude variable. If None, common names will be auto-detected.
    lon_var : str | None, optional
        The name of the longitude variable. If None, common names will be auto-detected.
    
    Returns
    -------
    tuple[np.ndarray, np.ndarray, str, str]
        A tuple containing the latitude array, longitude array, latitude variable name, and longitude variable name.
    
    Raises
    ------
     KeyError
        If latitude or longitude variables cannot be found in the dataset.
    
    Notes
    -----
    - The function first checks for the presence of the specified latitude and longitude variable names.
    - If either is not provided, it attempts to auto-detect common variable names for latitude (e.g., "lat", "latitude", "y")
      and longitude (e.g., "lon", "longitude", "x").
    - If the required variables are not found, a KeyError is raised.
    
    """
    if lat_var is None:
        for candidate in ("lat", "latitude", "y"):
            if candidate in ds_like:
                lat_var = candidate
                break
    if lon_var is None:
        for candidate in ("lon", "longitude", "x"):
            if candidate in ds_like:
                lon_var = candidate
                break
    if lat_var is None or lat_var not in ds_like:
        raise KeyError("Could not locate a latitude variable in the dataset.")
    if lon_var is None or lon_var not in ds_like:
        raise KeyError("Could not locate a longitude variable in the dataset.")
    return ds_like[lat_var].values, ds_like[lon_var].values, lat_var, lon_var


def to_2d_lat_lon(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert 1-D lat/lon arrays to 2-D if necessary, ensuring they are compatible.
    
    Parameters
    ----------
    lat : np.ndarray
        The latitude array, which can be either 1-D or 2-D.
    lon : np.ndarray
        The longitude array, which can be either 1-D or 2-D.
    
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        A tuple containing the latitude and longitude arrays, both in 2-D format.
    
    Raises
    ------
    ValueError
        If the latitude and longitude arrays are not both 1-D or both 2-D, or if their shapes are incompatible for broadcasting.
    
    Notes
    -----
    - If both latitude and longitude are 1-D, they are converted to 2-D using broadcasting, resulting in arrays of shape (lat.size, lon.size).
    - If both latitude and longitude are already 2-D, they are returned as-is.
    - If the arrays do not meet these conditions, a ValueError is raised indicating the issue.
    """
    if lat.ndim == 1 and lon.ndim == 1:
        lat2d = np.broadcast_to(lat[:, None], (lat.size, lon.size))
        lon2d = np.broadcast_to(lon[None, :], (lat.size, lon.size))
        return lat2d, lon2d
    if lat.ndim == 2 and lon.ndim == 2:
        return lat, lon
    raise ValueError("Latitude/longitude must both be 1-D or both be 2-D.")


def compute_clip_window_from_reference_bbox(
    ds_like: xr.Dataset,
    reference_grid_file: str,
    reference_lat_var: str | None = None,
    reference_lon_var: str | None = None,
    pad_pixels: int = 0,
) -> ClipWindow:
    """Compute a clip window for the source grid based on the bounding box of a reference grid file.
    
    Parameters
    ----------
    ds_like : xr.Dataset
        A dataset from which to extract the source grid's latitude and longitude for comparison.
    reference_grid_file : str
        The path to the reference NetCDF file containing latitude and longitude variables that define the bounding box.
    reference_lat_var : str | None, optional
        The name of the latitude variable in the reference file. If None, common names will be auto-detected.
    reference_lon_var : str | None, optional
        The name of the longitude variable in the reference file. If None, common names will be auto-detected.
    pad_pixels : int, default=0
        The number of pixels to expand the clip window on each side beyond the bounding box of the reference grid.
    
    Returns
    -------
    ClipWindow
        A ClipWindow object containing the indices of the source grid that correspond to the bounding box of
        the reference grid, expanded by the specified padding.
    
    Raises
    ------
    ValueError
        If the bounding box of the reference grid does not overlap with the source grid.
    
    Notes
    -----
    - The function first extracts the latitude and longitude from the reference grid file, normalizing longitudes to the range (-180, 180].
    - It then extracts the latitude and longitude from the source dataset, also normalizing longitudes.
    - The function computes a boolean mask of the source grid points that fall within the bounding box defined by the reference grid.
    - If no points in the source grid fall within the reference bounding box, a ValueError is raised.
    - The function then computes the minimum and maximum indices of the source grid that encompass the reference bounding box, applying
      the specified padding while ensuring the indices remain within the bounds of the source grid."""
    with xr.open_dataset(reference_grid_file) as ds_ref:
        ref_lat, ref_lon_raw, _, _ = get_lat_lon(ds_ref, lat_var=reference_lat_var, lon_var=reference_lon_var)
        ref_lon = normalize_lon180(ref_lon_raw)

    src_lat, src_lon, _, _ = get_lat_lon(ds_like)
    src_lon = normalize_lon180(src_lon)
    src_lat2d, src_lon2d = to_2d_lat_lon(src_lat, src_lon)

    lat_min = float(np.nanmin(ref_lat))
    lat_max = float(np.nanmax(ref_lat))
    lon_min = float(np.nanmin(ref_lon))
    lon_max = float(np.nanmax(ref_lon))

    inside = (
        (src_lat2d >= lat_min)
        & (src_lat2d <= lat_max)
        & (src_lon2d >= lon_min)
        & (src_lon2d <= lon_max)
    )
    if not np.any(inside):
        raise ValueError("Reference-file bounding box does not overlap the source grid.")

    jj, ii = np.where(inside)
    j0 = max(0, int(np.min(jj)) - pad_pixels)
    j1 = min(src_lat2d.shape[0] - 1, int(np.max(jj)) + pad_pixels)
    i0 = max(0, int(np.min(ii)) - pad_pixels)
    i1 = min(src_lat2d.shape[1] - 1, int(np.max(ii)) + pad_pixels)
    return ClipWindow(j0=j0, j1=j1, i0=i0, i1=i1)


def expand_clip_window(
    window: ClipWindow,
    expand_pixels: int,
    max_rows: int,
    max_cols: int,
) -> ClipWindow:
    """Return a clip window expanded by ``expand_pixels`` on each side."""
    return ClipWindow(
        j0=max(0, window.j0 - expand_pixels),
        j1=min(max_rows - 1, window.j1 + expand_pixels),
        i0=max(0, window.i0 - expand_pixels),
        i1=min(max_cols - 1, window.i1 + expand_pixels),
    )


def compute_clip_window_from_bbox(
    ds_like: xr.Dataset,
    lat_min: float,
    lat_max: float,
    lon_min: float,
    lon_max: float,
    pad_pixels: int = 0,
) -> ClipWindow:
    """Return a ClipWindow for ``ds_like`` covering the given lat/lon bbox."""
    src_lat, src_lon, _, _ = get_lat_lon(ds_like)
    src_lon = normalize_lon180(src_lon)
    src_lat2d, src_lon2d = to_2d_lat_lon(src_lat, src_lon)
    inside = (
        (src_lat2d >= lat_min) & (src_lat2d <= lat_max)
        & (src_lon2d >= lon_min) & (src_lon2d <= lon_max)
    )
    if not np.any(inside):
        return ClipWindow(
            j0=0,
            j1=src_lat2d.shape[0] - 1,
            i0=0,
            i1=src_lat2d.shape[1] - 1,
        )
    jj, ii = np.where(inside)
    j0 = max(0, int(np.min(jj)) - pad_pixels)
    j1 = min(src_lat2d.shape[0] - 1, int(np.max(jj)) + pad_pixels)
    i0 = max(0, int(np.min(ii)) - pad_pixels)
    i1 = min(src_lat2d.shape[1] - 1, int(np.max(ii)) + pad_pixels)
    return ClipWindow(j0=j0, j1=j1, i0=i0, i1=i1)


def slice_2d(arr: np.ndarray, window: ClipWindow) -> np.ndarray:
    """Slice a 2-D array according to the specified clip window.
    
    Parameters
    ----------
    arr : np.ndarray
        A 2-D array to be sliced.
    window : ClipWindow
        A ClipWindow object containing the indices for slicing the array.
    
    Returns
    -------
    np.ndarray
        A sliced 2-D array corresponding to the specified clip window."""
    return arr[window.j0 : window.j1 + 1, window.i0 : window.i1 + 1]


def sample_src_on_target_grid(
    src_data: np.ndarray,
    src_lat: np.ndarray,
    src_lon: np.ndarray,
    tgt_lat: np.ndarray,
    tgt_lon: np.ndarray,
) -> np.ndarray:
    """Sample source data onto the target grid using nearest-neighbor interpolation.
    
    Parameters
    ----------
    src_data : np.ndarray
        The source data array to be sampled, with shape (M, N).
    src_lat : np.ndarray
        The latitude array corresponding to the source data, which can be either 1-D or 2-D.
    src_lon : np.ndarray
        The longitude array corresponding to the source data, which can be either 1-D or 2-D.
    tgt_lat : np.ndarray
        The latitude array of the target grid, which can be either 1-D or 2-D.
    tgt_lon : np.ndarray
        The longitude array of the target grid, which can be either 1-D or 2-D.
    
    Returns
    -------
    np.ndarray
        An array of the same shape as the target grid containing the source data sampled at the target grid points.
    
    Raises
    ------
    ValueError
        If the target grid extends beyond the coverage of the source grid when both are 1-D,
        or if the latitude and longitude arrays are not compatible for interpolation.
    """
    src_lon = normalize_lon180(src_lon)
    tgt_lon = normalize_lon180(tgt_lon)

    if src_lat.ndim == 1 and src_lon.ndim == 1 and tgt_lat.ndim == 1 and tgt_lon.ndim == 1:
        if np.min(tgt_lat) < np.min(src_lat) or np.max(tgt_lat) > np.max(src_lat):
            raise ValueError("Target latitude outside source coverage")
        if np.min(tgt_lon) < np.min(src_lon) or np.max(tgt_lon) > np.max(src_lon):
            raise ValueError("Target longitude outside source coverage")
        j_idx = np.array([np.argmin(np.abs(src_lat - v)) for v in tgt_lat])
        i_idx = np.array([np.argmin(np.abs(src_lon - v)) for v in tgt_lon])
        return src_data[np.ix_(j_idx, i_idx)]

    if src_lat.ndim == 1 and src_lon.ndim == 1:
        src_lon2d, src_lat2d = np.meshgrid(src_lon, src_lat)
    else:
        src_lat2d, src_lon2d = src_lat, src_lon
    pts_src = np.column_stack([src_lat2d.ravel(), src_lon2d.ravel()])
    vals_src = src_data.ravel()

    if tgt_lat.ndim == 1 and tgt_lon.ndim == 1:
        tgt_lon2d, tgt_lat2d = np.meshgrid(tgt_lon, tgt_lat)
    else:
        tgt_lat2d, tgt_lon2d = tgt_lat, tgt_lon
    pts_tgt = np.column_stack([tgt_lat2d.ravel(), tgt_lon2d.ravel()])

    tree = cKDTree(pts_src)
    _, idx = tree.query(pts_tgt)
    return vals_src[idx].reshape(tgt_lat2d.shape)

def resample_to_target_grid(
    src_data: np.ndarray,
    ds_src: xr.Dataset,
    tgt_lat: np.ndarray,
    tgt_lon: np.ndarray,
) -> np.ndarray:
    """Resample a 2-D source array onto a target lat/lon grid using nearest-neighbour.

    If source and target grids are identical (same shape and values), returns
    ``src_data`` directly without copying.
    """
    src_lat, src_lon, _, _ = get_lat_lon(ds_src)
    src_lon_n = normalize_lon180(src_lon)
    tgt_lon_n = normalize_lon180(tgt_lon)

    if (
        src_lat.shape == tgt_lat.shape
        and src_lon_n.shape == tgt_lon_n.shape
        and np.allclose(src_lat, tgt_lat, atol=1e-6)
        and np.allclose(src_lon_n, tgt_lon_n, atol=1e-6)
    ):
        return src_data

    if src_lat.ndim == 1 and src_lon_n.ndim == 1 and tgt_lat.ndim == 1 and tgt_lon_n.ndim == 1:
        j_idx = np.array([int(np.argmin(np.abs(src_lat - v))) for v in tgt_lat])
        i_idx = np.array([int(np.argmin(np.abs(src_lon_n - v))) for v in tgt_lon_n])
        return src_data[np.ix_(j_idx, i_idx)]

    if src_lat.ndim == 1 and src_lon_n.ndim == 1:
        src_lon2d, src_lat2d = np.meshgrid(src_lon_n, src_lat)
    else:
        src_lat2d, src_lon2d = src_lat, src_lon_n

    if tgt_lat.ndim == 1 and tgt_lon_n.ndim == 1:
        tgt_lon2d, tgt_lat2d = np.meshgrid(tgt_lon_n, tgt_lat)
    else:
        tgt_lat2d, tgt_lon2d = tgt_lat, tgt_lon_n

    tree = cKDTree(np.column_stack([src_lat2d.ravel(), src_lon2d.ravel()]))
    _, idx = tree.query(np.column_stack([tgt_lat2d.ravel(), tgt_lon2d.ravel()]))
    return src_data.ravel()[idx].reshape(tgt_lat2d.shape)


def _subset_array_to_window(
    arr: np.ndarray,
    ds: xr.Dataset,
    window: ClipWindow,
) -> tuple[np.ndarray, xr.Dataset]:
    """Return a spatially subsetted array and a lightweight in-memory dataset
    whose lat/lon coordinates reflect the subset, for use with
    ``resample_to_target_grid``."""
    sub_arr = slice_2d(arr, window)
    lat_src, lon_src, lat_name, lon_name = get_lat_lon(ds)
    lat_sub = (
        lat_src[window.j0 : window.j1 + 1]
        if lat_src.ndim == 1
        else slice_2d(lat_src, window)
    )
    lon_sub = (
        lon_src[window.i0 : window.i1 + 1]
        if lon_src.ndim == 1
        else slice_2d(lon_src, window)
    )
    lat_dims = (lat_name,) if lat_src.ndim == 1 else ("y", "x")
    lon_dims = (lon_name,) if lon_src.ndim == 1 else ("y", "x")
    ds_sub = xr.Dataset(
        {
            lat_name: (lat_dims, lat_sub),
            lon_name: (lon_dims, lon_sub),
        }
    )
    return sub_arr, ds_sub


# -----------------------------------------------------------------------------
# Algorithm
# -----------------------------------------------------------------------------

def calculate_slope(x: np.ndarray, y: np.ndarray) -> float:
    """Calculate the slope of a linear fit to the given x and y data, ignoring NaNs.
    
    Parameters
    ----------
    x : np.ndarray
        The x-coordinates of the data points.
    y : np.ndarray
        The y-coordinates of the data points.
    
    Returns
    -------
    float
        The slope of the linear fit. If the denominator is zero or NaN, returns NaN.
    
    Notes
    -----
    - The function first converts the input arrays to float and computes their means, ignoring NaN values.
    - It then calculates the centered x values and the denominator for the slope calculation.
    - If the denominator is zero or NaN, the function returns NaN to indicate that the slope cannot be computed.
    - Otherwise, it computes and returns the slope as the sum of the product of the centered x values and the centered y values, divided by the denominator.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_mean = np.nanmean(x)
    y_mean = np.nanmean(y)
    xc = x - x_mean
    denom = np.nansum(xc**2)
    if np.isnan(denom) or np.isclose(denom, 0.0):
        return np.nan
    return float(np.nansum(xc * (y - y_mean)) / denom)


def apply_cutoff(arr: np.ndarray, min_value: float | None, max_value: float | None) -> np.ndarray:
    """Apply minimum and/or maximum cutoffs to an array, ignoring NaN values.
    
    Parameters
    ----------
    arr : np.ndarray
        The input array.
    min_value : float | None
        The minimum value to apply.
    max_value : float | None
        The maximum value to apply.
    
    Returns
    -------
    np.ndarray
        The array with cutoffs applied.
    """
    out = np.array(arr, copy=True)
    valid = ~np.isnan(out)
    if min_value is not None:
        out[valid] = np.maximum(out[valid], min_value)
    if max_value is not None:
        out[valid] = np.minimum(out[valid], max_value)
    return out


def get_water_buffer_mask(
    landmask: np.ndarray,
    buffer_width: int,
    landmask_trueval: int | float,
) -> np.ndarray:
    """Get a boolean mask for water pixels that are within a specified buffer width of land pixels.
    
    Parameters
    ----------
    landmask : np.ndarray
        A 2-D array representing the land mask, where land pixels are indicated by a specific
        true value and water pixels are indicated by other values or NaNs.
    buffer_width : int
        The width of the buffer around land pixels, in terms of the number of pixels.
    landmask_trueval : int | float
        The value in the landmask array that indicates land pixels.
    
    Returns
    -------
    np.ndarray
        A boolean array of the same shape as the input landmask, where True values indicate water pixels
        that are within the specified buffer width of land pixels, and False values indicate all other pixels.
    
    Notes
    -----
    - The function first creates boolean arrays to identify land and water pixels based on the provided landmask
      and true value.
    - It then uses binary dilation to expand the land pixels by the specified buffer width, creating a mask of land pixels
      that includes the buffer zone.
    - Finally, the function returns a boolean array where True values correspond to water pixels that are adjacent to the
      dilated land mask, effectively identifying water pixels that are within the buffer zone of land pixels."""
    is_land = np.asarray(landmask) == landmask_trueval
    is_water = ~is_land & ~np.isnan(landmask)
    structure = np.ones((2 * buffer_width + 1, 2 * buffer_width + 1), dtype=bool)
    land_dilated = ndimage.binary_dilation(is_land, structure=structure)
    return is_water & land_dilated


def get_inland_water_mask(landmask: np.ndarray, landmask_trueval: int | float) -> np.ndarray:
    """Get a boolean mask for water pixels that are not connected to the border of the array, indicating inland water bodies.
    
    Parameters
    ----------
    landmask : np.ndarray
        A 2-D array representing the land mask, where land pixels are indicated by a specific
        true value and water pixels are indicated by other values or NaNs.
    landmask_trueval : int | float
        The value in the landmask array that indicates land pixels.
    
    Returns
    -------
    np.ndarray
        A boolean array of the same shape as the input landmask, where True values indicate water pixels
        that are not connected to the border of the array, and False values indicate all other pixels. 
    Notes
    -----
    - The function first creates boolean arrays to identify land and water pixels based on the provided landmask and true value.
    - It then uses connected component labeling to identify contiguous regions of water pixels in the array.
    - The function identifies which of these water regions are connected to the border of the array,
      which are considered to be non-inland water bodies (e.g., oceans, large lakes).
    - Finally, the function returns a boolean array where True values correspond to water pixels that are
      not connected to the border, effectively identifying inland water bodies such as small lakes and rivers.
    
    """
    is_land = np.asarray(landmask) == landmask_trueval
    is_water = ~is_land & ~np.isnan(landmask)
    structure = np.ones((3, 3), dtype=bool)
    labels, n_labels = ndimage.label(is_water, structure=structure)
    if n_labels == 0:
        return np.zeros_like(is_water, dtype=bool)

    border_ids = np.unique(np.concatenate([labels[0, :], labels[-1, :], labels[:, 0], labels[:, -1]]))
    inland = is_water.copy()
    for label_id in border_ids:
        if label_id == 0:
            continue
        inland[labels == label_id] = False
    return inland


def apply_water_imprinting(
    lapse_rate: np.ndarray,
    landmask: np.ndarray,
    fallback_lapse_rate_m: float,
    buffer_width: int,
    landmask_trueval: int | float,
    save_units_as_per_km: bool,
) -> np.ndarray:
    """Apply water imprinting to the lapse rate array by filling in values for water pixels based on proximity to land and a fallback lapse rate.
    
    Parameters
    ----------
    lapse_rate : np.ndarray
        A 2-D array containing the computed lapse rate values, where land pixels have valid values and water pixels may contain NaNs or
        other placeholder values.
    landmask : np.ndarray
        A 2-D array representing the land mask, where land pixels are indicated by a specific true value and water pixels are indicated
        by other values or NaNs.
    fallback_lapse_rate_m : float
        The lapse rate value to use for water pixels that are within the buffer zone of land pixels or are identified as inland water bodies.
        This value is typically a standard lapse rate (e.g., -0.0065 K/m).
    buffer_width : int
        The width of the buffer around land pixels, in terms of the number of pixels, to determine which water pixels should be filled with
        the fallback lapse rate.
    landmask_trueval : int | float
        The value in the landmask array that indicates land pixels.
    save_units_as_per_km : bool
        A boolean flag indicating whether the lapse rate values are saved in units of K/km (True) or K/m (False). If True, the fallback lapse
        rate will be multiplied by 100
        to convert from K/m to K/km before being applied to the water pixels.
    
    Returns
    -------
    np.ndarray
        A 2-D array of the same shape as the input lapse_rate, where water pixels that are within the buffer zone of land pixels or are identified as inland water bodies have been filled with the specified fallback lapse rate, and all other pixels retain their original values.
    
    Notes
    -----
    - The function first identifies water pixels that are within the buffer zone of land pixels using the `get_water_buffer_mask` function.
    - It then identifies inland water bodies using the `get_inland_water_mask` function.
    - The fallback lapse rate is applied to these identified water pixels, converting units if necessary based on the `save_units_as_per_km` flag.
    """
    out = np.array(lapse_rate, copy=True)
    shoreline = get_water_buffer_mask(landmask, buffer_width, landmask_trueval)
    inland = get_inland_water_mask(landmask, landmask_trueval)
    fill_value = fallback_lapse_rate_m * (1000.0 if save_units_as_per_km else 1.0)
    out[shoreline | inland] = fill_value
    return out


def compute_lapse_rate_grid(
    temperature: np.ndarray,
    elevation: np.ndarray,
    landmask: np.ndarray,
    landmask_trueval: int | float = 1,
    min_land_neighbors: int = 5,
    fallback_lapse_rate_m: float = -0.0065,
    min_lapse_rate_cutoff: float | None = None,
    max_lapse_rate_cutoff: float | None = None,
    save_units_as_per_km: bool = True,
    buffer_width: int = 1,
    fill_water_mask: bool = True,
    min_elev_diff_std: float | None = None,
) -> np.ndarray:
    """Compute a dynamic lapse rate grid based on local temperature and elevation differences, applying various criteria for land/water masking and cutoffs.
    
    Parameters
    ----------
    temperature : np.ndarray
        A 2-D array containing temperature values at each grid point, where valid values correspond to land pixels and may contain NaNs for water pixels.
    elevation : np.ndarray
        A 2-D array containing elevation values at each grid point, where valid values correspond to land pixels and may contain NaNs for water pixels.
    landmask : np.ndarray
        A 2-D array representing the land mask, where land pixels are indicated by a specific true value and water pixels are indicated by other
        values or NaNs.
    landmask_trueval : int | float, default=1
        The value in the landmask array that indicates land pixels. This value is used to identify which pixels in the landmask correspond to land and
        which correspond to water.
    min_land_neighbors : int, default=5
        The minimum number of valid land neighbors required to compute a lapse rate for a given pixel.
        If a pixel has fewer than this number of valid land neighbors, the fallback lapse rate will be used instead.
    fallback_lapse_rate_m : float, default=-0.0065
        The lapse rate value to use for pixels that do not meet the criteria for valid land neighbors or have insufficient data for slope calculation.
        This value is typically a standard lapse rate (e.g., -0.0065 K/m).
    min_lapse_rate_cutoff : float | None, optional
        The minimum lapse rate value to apply as a cutoff. If None, no minimum cutoff will be applied.
    max_lapse_rate_cutoff : float | None, optional
        The maximum lapse rate value to apply as a cutoff. If None, no maximum cutoff will be applied.
    save_units_as_per_km : bool, default=True
        A boolean flag indicating whether the computed lapse rate values should be saved in units of K/km (True) or K/m (False).
        If True, the computed lapse rates will be multiplied by 1000.0 to convert from K/m to K/km before being returned.
    buffer_width : int, default=1
        The width of the buffer around land pixels, in terms of the number of pixels, to determine which water pixels should be filled with the fallback
        lapse rate when `fill_water_mask` is True.
    fill_water_mask : bool, default=True
        A boolean flag indicating whether to apply water imprinting to the computed lapse rate grid. If True, water pixels that are within the buffer
        zone of land pixels or are identified as inland water bodies will be filled with the specified fallback lapse rate using the `apply_water_imprinting`
        function.
    
    Returns
    -------
    np.ndarray
        A 2-D array containing the computed lapse rate values for each grid point, with water pixels filled according to the specified criteria and cutoffs.
        The shape of the output array will be the same as the input temperature, elevation, and landmask arrays.
    """
    if temperature.shape != elevation.shape or temperature.shape != landmask.shape:
        raise ValueError("Temperature, elevation, and landmask arrays must share the same shape.")

    output = np.full(elevation.shape, np.nan, dtype=np.float32)
    half_window = 1
    is_land = np.asarray(landmask) == landmask_trueval
    nrows, ncols = elevation.shape

    for row in range(half_window, nrows - half_window):
        for col in range(half_window, ncols - half_window):
            if not is_land[row, col]:
                continue

            elev_window = elevation[row - half_window : row + half_window + 1, col - half_window : col + half_window + 1]
            temp_window = temperature[row - half_window : row + half_window + 1, col - half_window : col + half_window + 1]
            land_window = landmask[row - half_window : row + half_window + 1, col - half_window : col + half_window + 1]

            if np.isnan(elev_window).any() or np.isnan(temp_window).all():
                continue

            center_index = 4
            dz = (elev_window - elev_window[1, 1]).ravel()
            dt = (temp_window - temp_window[1, 1]).ravel()
            land_flat = land_window.ravel()

            neighbor_land = np.delete(land_flat, center_index) == landmask_trueval
            neighbor_dz = np.delete(dz, center_index)
            valid_neighbors = neighbor_land & ~np.isnan(neighbor_dz)

            if int(np.sum(valid_neighbors)) < min_land_neighbors:
                output[row, col] = fallback_lapse_rate_m
                continue
            
            if min_elev_diff_std is not None:
                elev_std = float(np.std(neighbor_dz[valid_neighbors]))
                if elev_std < min_elev_diff_std:
                    output[row, col] = fallback_lapse_rate_m
                    continue
            
            valid_all = ~np.isnan(dz) & ~np.isnan(dt)
            valid_all[center_index] = True

            if int(np.sum(valid_all)) < 2:
                output[row, col] = fallback_lapse_rate_m
                continue

            slope = calculate_slope(dz[valid_all], dt[valid_all])
            output[row, col] = slope

    if save_units_as_per_km:
        output = apply_cutoff(
            output,
            None if min_lapse_rate_cutoff is None else min_lapse_rate_cutoff / 1000.0,
            None if max_lapse_rate_cutoff is None else max_lapse_rate_cutoff / 1000.0,
        )
        output = output * 1000.0
    else:
        output = apply_cutoff(output, min_lapse_rate_cutoff, max_lapse_rate_cutoff)

    if fill_water_mask:
        output = apply_water_imprinting(
            output,
            landmask=landmask,
            fallback_lapse_rate_m=fallback_lapse_rate_m,
            buffer_width=buffer_width,
            landmask_trueval=landmask_trueval,
            save_units_as_per_km=save_units_as_per_km,
        )

    output = output.astype(np.float32, copy=False)
    LOGGER.debug("Finished lapse-rate computation for one time slice.")
    return output


# -----------------------------------------------------------------------------
# IO helpers
# -----------------------------------------------------------------------------

def open_array(file_path: str | Path, var_name: str) -> tuple[np.ndarray, xr.Dataset]:
    """Open a NetCDF file and extract the specified variable as a NumPy array, along with the dataset for further use.
    
    Parameters
    ----------
    file_path : str | Path
        The path to the NetCDF file to be opened.
    var_name : str
        The name of the variable to extract from the dataset.
    
    Returns
    -------
    tuple[np.ndarray, xr.Dataset]
        A tuple containing the extracted variable as a NumPy array and the opened xarray Dataset object.
    """
    ds = xr.open_dataset(file_path, decode_times=True)
    if var_name not in ds:
        ds.close()
        raise KeyError(f"Variable {var_name!r} not found in {file_path!s}")
    return ds[var_name].values, ds


def save_lapse_rate_netcdf(
    output_path: str | Path,
    lapse_rate: np.ndarray,
    time_values: Iterable[datetime] | None,
    lat: np.ndarray,
    lon: np.ndarray,
    lat_name: str,
    lon_name: str,
    save_units_as_per_km: bool,
    time_attrs: dict | None = None,
    extra_attrs: dict | None = None,
    lat_attrs_in: dict | None = None,
    lon_attrs_in: dict | None = None,
) -> None:
    """Save the computed lapse rate grid to a NetCDF file with appropriate metadata and coordinate variables.
    
    Parameters
    ----------
    output_path : str | Path
        The path where the output NetCDF file should be saved.
    lapse_rate : np.ndarray
        A 2-D or 3-D array containing the computed lapse rate values, where the dimensions correspond to time (optional), latitude, and longitude.
    time_values : Iterable[datetime] | None
        An iterable of datetime objects representing the time coordinates for the lapse rate data. If None, no time coordinate will be included in the
        output.
    lat : np.ndarray
        A 1-D or 2-D array containing the latitude values corresponding to the lapse rate grid.
    lon : np.ndarray
        A 1-D or 2-D array containing the longitude values corresponding to the lapse rate grid.
    lat_name : str
        The name to use for the latitude coordinate variable in the output NetCDF file.
    lon_name : str
        The name to use for the longitude coordinate variable in the output NetCDF file.
    save_units_as_per_km : bool
        A boolean flag indicating whether the lapse rate values are saved in units of K/km (True) or K/m (False). This will affect the units attribute
        of the lapse rate variable in the output file.
    time_attrs : dict | None, optional
        A dictionary of attributes to include for the time coordinate variable. If None, no additional attributes will be included beyond the automatically
        generated ones.
    extra_attrs : dict | None, optional
        A dictionary of additional global attributes to include in the output NetCDF file. If None, no additional global attributes will be included beyond
        the automatically generated ones.
    lat_attrs_in : dict | None, optional
        A dictionary of attributes to include for the latitude coordinate variable. If None, default attributes will be used.
    lon_attrs_in : dict | None, optional
        A dictionary of attributes to include for the longitude coordinate variable. If None, default attributes will be used.
    
    Returns
    -------
    None
        The function saves the lapse rate data to a NetCDF file at the specified output path and does not return any value.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    units = "K/km" if save_units_as_per_km else "K/m"

    time_list = list(time_values) if time_values is not None else None
    t_attrs = dict(time_attrs) if time_attrs else {}
    for _k in ("valid_range", "vmin", "vmax"):
        t_attrs.pop(_k, None)

    if time_list is not None and len(time_list) > 0:
        t0 = pd.Timestamp(time_list[0])
        t_attrs["begin_date"] = int(t0.strftime("%Y%m%d"))
        t_attrs["begin_time"] = int(t0.strftime("%H%M%S"))
        t_attrs["time_increment"] = 10000
        t_attrs["long_name"] = "time"

    if time_list is not None:
        t0 = pd.Timestamp(time_list[0])
        units_str = f"minutes since {t0.strftime('%Y-%m-%dT%H:%M:%S')}"
        minutes_vals = np.array(
            [int(round((pd.Timestamp(t) - t0).total_seconds() / 60.0)) for t in time_list],
            dtype=np.int32,
        )
        t_attrs["units"] = units_str
        t_attrs["calendar"] = "proleptic_gregorian"
        time_numeric = minutes_vals
    else:
        time_numeric = None

    lr_attrs = {"units": units}
    lat_attrs = dict(lat_attrs_in) if lat_attrs_in else {}
    lat_attrs.setdefault("units", "degrees_north")
    lat_attrs.setdefault("long_name", "latitude")
    lat_attrs.setdefault("vmin", float(np.nanmin(lat)))
    lat_attrs.setdefault("vmax", float(np.nanmax(lat)))

    lon_attrs = dict(lon_attrs_in) if lon_attrs_in else {}
    lon_attrs.setdefault("units", "degrees_east")
    lon_attrs.setdefault("long_name", "longitude")
    lon_attrs.setdefault("vmin", float(np.nanmin(lon)))
    lon_attrs.setdefault("vmax", float(np.nanmax(lon)))

    if lat.ndim == 1 and lon.ndim == 1:
        if lapse_rate.ndim == 3:
            ds_out = xr.Dataset(
                {"lapse_rate": (("time", lat_name, lon_name), lapse_rate, lr_attrs)},
                coords={
                    "time": ("time", time_numeric, t_attrs),
                    lat_name: (lat_name, lat, lat_attrs),
                    lon_name: (lon_name, lon, lon_attrs),
                },
            )
        else:
            ds_out = xr.Dataset(
                {"lapse_rate": ((lat_name, lon_name), lapse_rate, lr_attrs)},
                coords={
                    lat_name: (lat_name, lat, lat_attrs),
                    lon_name: (lon_name, lon, lon_attrs),
                },
            )
            if time_numeric is not None and len(time_numeric) == 1:
                ds_out["time"] = xr.DataArray(time_numeric[0], attrs=t_attrs)
    else:
        if lapse_rate.ndim == 3:
            ds_out = xr.Dataset(
                {
                    "lapse_rate": (("time", "y", "x"), lapse_rate, lr_attrs),
                    lat_name: (("y", "x"), lat, lat_attrs),
                    lon_name: (("y", "x"), lon, lon_attrs),
                },
                coords={"time": ("time", time_numeric, t_attrs)},
            )
        else:
            ds_out = xr.Dataset(
                {
                    "lapse_rate": (("y", "x"), lapse_rate, lr_attrs),
                    lat_name: (("y", "x"), lat, lat_attrs),
                    lon_name: (("y", "x"), lon, lon_attrs),
                }
            )
            if time_numeric is not None and len(time_numeric) == 1:
                ds_out["time"] = xr.DataArray(time_numeric[0], attrs=t_attrs)

    ds_out.attrs.update(
        {
            "title": "Dynamic lapse rate output",
            "source": "Generated by compute_dynamic_lapse_rates.py",
            "history": f"Created on {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%S.%f')}Z",
        }
    )
    if extra_attrs:
        ds_out.attrs.update(extra_attrs)

    encoding = {"lapse_rate": {"dtype": "float32", "_FillValue": np.nan}}
    ds_out.to_netcdf(output_path, encoding=encoding, format="NETCDF3_64BIT")
    LOGGER.info("Wrote %s", output_path)


# -----------------------------------------------------------------------------
# Pipeline helpers
# -----------------------------------------------------------------------------

def _output_mode_has_hour(date_fmt: str) -> bool:
    """Determine if the output date format includes an hour component based on the presence of 'HH' in the format string."""
    return "HH" in date_fmt


def _validate_output_mode(output_mode: str, outfile_date_format: str) -> None:
    """Validate the compatibility of the output mode with the specified output date format, ensuring that daily output does not include an hour component
    and hourly output does.
    
    Parameters
    ----------
    output_mode : str
        The output mode, which can be 'daily', 'hourly', or 'match_input'.
    outfile_date_format : str
        The date format string used for the output file naming, which should include 'HH' if the output mode is 'hourly' and should not include 'HH' if the
        output mode is 'daily'.
    
    Raises
    ------
    ValueError
        If the output mode is 'daily' but the output date format includes an hour component, or if the output mode is 'hourly' but the output date format
        does not include an hour component, indicating an incompatibility between the output mode and the date format.
    """
    has_hour = _output_mode_has_hour(outfile_date_format)
    if output_mode == "daily" and has_hour:
        raise ValueError(
            f"output_mode='daily' is incompatible with outfile_date_format={outfile_date_format!r}. "
            "Use 'YYYYMMDD' or 'YYYY-MM-DD' for daily output."
        )
    if output_mode == "hourly" and not has_hour:
        raise ValueError(
            f"output_mode='hourly' requires an hour component in outfile_date_format, but got {outfile_date_format!r}. "
            "Use 'YYYY-MM-DDTHH' or 'YYYY-MM-DD-HH'."
        )


def _parse_yyyymmdd(token: str) -> datetime:
    """Parse a date string in 'YYYYMMDD' format into a datetime object."""
    return datetime.strptime(token, "%Y%m%d")


def _parse_optional_datetime(token: str | None) -> datetime | None:
    """Parse an optional datetime string in 'YYYY-MM-DD HH:MM' format."""
    if token is None:
        return None
    return datetime.strptime(token, "%Y-%m-%d %H:%M")


def _floor_to_hour(dt: datetime) -> datetime:
    """Return datetime floored to the start of its hour."""
    return dt.replace(minute=0, second=0, microsecond=0)


def _get_extra_attrs(config: Config) -> dict:
    """Collect extra attributes for the output dataset based on the provided configuration, including source data, variable information, input description,
    and method reference if available."""
    extra_attrs: dict = {}
    if config.output_source_data is not None:
        extra_attrs["source_data"] = config.output_source_data
    if config.output_source_variable is not None:
        extra_attrs["source_variable"] = config.output_source_variable
    if config.output_source_variable_description is not None:
        extra_attrs["source_variable_description"] = config.output_source_variable_description
    if config.output_input_description is not None:
        extra_attrs["input_description"] = config.output_input_description
    if config.output_method_reference is not None:
        extra_attrs["method_reference"] = config.output_method_reference
    return extra_attrs


def _save_hourly(
    lapse: np.ndarray,
    valid_time: datetime,
    lat: np.ndarray,
    lon: np.ndarray,
    lat_name: str,
    lon_name: str,
    outfile_spec: FileSpec,
    outpath: Path,
    save_km: bool,
    force: bool,
    dryrun: bool,
    time_attrs: dict | None = None,
    extra_attrs: dict | None = None,
    lat_attrs_in: dict | None = None,
    lon_attrs_in: dict | None = None,
) -> None:
    """Save the computed lapse rate grid for a specific valid time to a NetCDF file, ensuring that the output file name is generated based on the provided
    date format and that existing files are not overwritten unless the force flag is set.
    
    Parameters
    ----------
    lapse : np.ndarray
        A 2-D array containing the computed lapse rate values for the specific valid time.
    valid_time : datetime
        The valid time corresponding to the lapse rate data, which will be used to generate the output file name based on the provided date format.
    lat : np.ndarray
        A 1-D or 2-D array containing the latitude values corresponding to the lapse rate grid.
    lon : np.ndarray
        A 1-D or 2-D array containing the longitude values corresponding to the lapse rate grid.
    lat_name : str
        The name to use for the latitude coordinate variable in the output NetCDF file.
    lon_name : str
        The name to use for the longitude coordinate variable in the output NetCDF file.
    outfile_spec : FileSpec
        A FileSpec object containing the prefix, suffix, and date format for generating the output file name based on the valid time.
    outpath : Path
        The directory path where the output NetCDF file should be saved.
    save_km : bool
        A boolean flag indicating whether the lapse rate values are saved in units of K/km (True) or K/m (False). This will affect the units attribute of
        the lapse rate variable in the output file.
    force : bool
        A boolean flag indicating whether to overwrite existing files. If False, the function will skip saving if the output file already exists.
        If True, it will overwrite existing files.
    dryrun : bool
        A boolean flag indicating whether to perform a dry run. If True, the function will log the intended output file path without actually saving the
        file. If False, it will proceed to save the file.
    time_attrs : dict | None, optional
        A dictionary of attributes to include for the time coordinate variable. If None, no additional attributes will be included beyond the automatically
        generated ones.
    extra_attrs : dict | None, optional
        A dictionary of additional global attributes to include in the output NetCDF file. If None, no additional global attributes will be included beyond
        the automatically generated ones.
    lat_attrs_in : dict | None, optional
        A dictionary of attributes to include for the latitude coordinate variable. If None, default attributes will be used.
    lon_attrs_in : dict | None, optional
        A dictionary of attributes to include for the longitude coordinate variable. If None, default attributes will be used.
    
    Returns
    -------
    None
        The function saves the lapse rate data for the specific valid time to a NetCDF file at the specified output path and does not return any value.
    """
    output_name = build_dated_filename(outfile_spec, valid_time)
    output_path = outpath / output_name
    if output_path.exists() and not force:
        LOGGER.info("Skipping existing file %s", output_path)
        return
    if dryrun:
        LOGGER.info("[Dry run] Would write %s", output_path)
        return
    save_lapse_rate_netcdf(
        output_path=output_path,
        lapse_rate=lapse,
        time_values=[valid_time],
        lat=lat,
        lon=lon,
        lat_name=lat_name,
        lon_name=lon_name,
        save_units_as_per_km=save_km,
        time_attrs=time_attrs,
        extra_attrs=extra_attrs,
        lat_attrs_in=lat_attrs_in,
        lon_attrs_in=lon_attrs_in,
    )


_UNITS_MAP: dict[str, str] = {
    "k/m": "per_m",
    "k per m": "per_m",
    "k m-1": "per_m",
    "k m^-1": "per_m",
    "k/km": "per_km",
    "k per km": "per_km",
    "k km-1": "per_km",
    "k km^-1": "per_km",
}


def _normalise_units_string(raw: str) -> str:
    return re.sub(r"\s+", " ", raw.strip().lower())


def _read_and_clip_preexisting(
    pre_path,
    pre_var: str,
    pre_time_var: str,
    clip_window,
    landmask: np.ndarray,
    fallback_lapse_rate_m: float,
    save_units_as_per_km: bool,
    min_lapse_rate_cutoff,
    max_lapse_rate_cutoff,
    buffer_width: int,
    fill_water_mask: bool,
    landmask_trueval,
):
    lr_data, ds_pre = open_array(pre_path, pre_var)
    time_vals, _ = get_time_values(ds_pre, pre_time_var, time_shift_minutes=0)
    lat_src, lon_src, lat_name, lon_name = get_lat_lon(ds_pre)

    if lr_data.ndim == 2:
        lr_data = lr_data[np.newaxis, :, :]

    if clip_window is not None:
        lat_out = (
            lat_src[clip_window.j0 : clip_window.j1 + 1]
            if lat_src.ndim == 1
            else slice_2d(lat_src, clip_window)
        )
        lon_out = (
            lon_src[clip_window.i0 : clip_window.i1 + 1]
            if lon_src.ndim == 1
            else slice_2d(lon_src, clip_window)
        )
        lr_clipped = lr_data[
            :,
            clip_window.j0 : clip_window.j1 + 1,
            clip_window.i0 : clip_window.i1 + 1,
        ]
    else:
        lat_out, lon_out = lat_src, lon_src
        lr_clipped = lr_data.copy()

    raw_units = ds_pre[pre_var].attrs.get("units", "") if pre_var in ds_pre else ""
    kind = _UNITS_MAP.get(_normalise_units_string(raw_units))
    want_km = bool(save_units_as_per_km)

    if kind == "per_m" and want_km:
        lr_clipped = lr_clipped * 1000.0
    elif kind == "per_km" and not want_km:
        lr_clipped = lr_clipped / 1000.0
    elif kind is None:
        LOGGER.warning(
            "Unrecognised or missing units %r in pre-existing LR file; assuming output units are already %s.",
            raw_units,
            "K/km" if want_km else "K/m",
        )

    for t in range(lr_clipped.shape[0]):
        lr_clipped[t] = apply_cutoff(lr_clipped[t], min_lapse_rate_cutoff, max_lapse_rate_cutoff)

    if fill_water_mask:
        fill_val = fallback_lapse_rate_m * (1000.0 if want_km else 1.0)
        water_buf = get_water_buffer_mask(landmask, buffer_width, landmask_trueval)
        inland_w = get_inland_water_mask(landmask, landmask_trueval)
        combined = water_buf | inland_w
        for t in range(lr_clipped.shape[0]):
            lr_clipped[t][combined] = fill_val

    lat_out_attrs = dict(ds_pre[lat_name].attrs)
    lon_out_attrs = dict(ds_pre[lon_name].attrs)
    ds_pre.close()
    return lr_clipped, time_vals, lat_out, lon_out, lat_name, lon_name, lat_out_attrs, lon_out_attrs


def _flush_daily_group(
    daily_groups: dict,
    day_key: pd.Timestamp,
    outfile_spec: FileSpec,
    outpath: Path,
    lat: np.ndarray,
    lon: np.ndarray,
    lat_name: str,
    lon_name: str,
    save_km: bool,
    force: bool,
    dryrun: bool,
    time_attrs: dict,
    extra_attrs: dict,
    lat_attrs_in: dict,
    lon_attrs_in: dict,
) -> None:
    entries = daily_groups.get(day_key)
    if not entries:
        return
    output_name = build_dated_filename(outfile_spec, day_key.to_pydatetime())
    output_path = outpath / output_name
    if output_path.exists() and not force:
        LOGGER.info("Skipping existing file %s", output_path)
        return
    if dryrun:
        LOGGER.info("[Dry run] Would write %s", output_path)
        return
    lapse_stack = np.stack([e[1] for e in entries], axis=0)
    time_list = [e[0] for e in entries]
    save_lapse_rate_netcdf(
        output_path=output_path,
        lapse_rate=lapse_stack,
        time_values=time_list,
        lat=lat,
        lon=lon,
        lat_name=lat_name,
        lon_name=lon_name,
        save_units_as_per_km=save_km,
        time_attrs=time_attrs,
        extra_attrs=extra_attrs,
        lat_attrs_in=lat_attrs_in,
        lon_attrs_in=lon_attrs_in,
    )


def run_pipeline(config: Config) -> None:
    """Run the main pipeline for computing dynamic lapse rates based on the provided configuration, including input validation, data loading, processing,
    and output saving.
    
    Parameters
    ----------
    config : Config
        A Config object containing all necessary parameters and settings for running the pipeline, including input file paths, variable names, date range,
        output settings, and other options.
    
    Returns
    -------
    None
        The function executes the pipeline to compute dynamic lapse rates and saves the output to NetCDF files according to the specified configuration,
        but does not return any value."""


    # Validate date range
    start_date = _parse_yyyymmdd(config.yyyymmdd_start)
    end_date = _parse_yyyymmdd(config.yyyymmdd_end)
    if end_date < start_date:
        raise ValueError("yyyymmdd_end must be greater than or equal to yyyymmdd_start.")

    # Optional hour-aware filtering. Minutes are intentionally ignored:
    # 06:15 and 06:30 both map to hour 06.
    start_dt_filter = _parse_optional_datetime(config.start_datetime)
    end_dt_filter = _parse_optional_datetime(config.end_datetime)

    if (start_dt_filter is None) != (end_dt_filter is None):
        raise ValueError("start_datetime and end_datetime must either both be provided or both be omitted.")

    if start_dt_filter is not None and end_dt_filter is not None:
        start_hour_filter = _floor_to_hour(start_dt_filter)
        end_hour_filter = _floor_to_hour(end_dt_filter)
        if end_hour_filter < start_hour_filter:
            raise ValueError("end_datetime must be on or after start_datetime.")
        LOGGER.info(
            "Applying hour-window filter from %s through %s (minutes ignored; compared by hour).",
            start_hour_filter.strftime("%Y-%m-%d %H:%M"),
            end_hour_filter.strftime("%Y-%m-%d %H:%M"),
        )
    else:
        start_hour_filter = None
        end_hour_filter = None

    # Determine output mode and validate compatibility with output date format
    output_mode = config.output_mode
    if output_mode == "match_input":
        output_mode = "hourly" if "THH" in config.infile_date_format else "daily"
        LOGGER.info(
            "output_mode='match_input' resolved to '%s' based on infile_date_format='%s'.",
            output_mode,
            config.infile_date_format,
        )

    _validate_output_mode(output_mode, config.outfile_date_format)

    # Set up file specifications for input and output files based on the provided configuration
    infile_spec = FileSpec(
        prefix=config.fname_in_pfx_T,
        suffix=config.fname_in_sfx_T,
        date_format=config.infile_date_format,
    )
    outfile_spec = FileSpec(
        prefix=config.fname_out_pfx,
        suffix=config.fname_out_sfx,
        date_format=config.outfile_date_format,
    )

    # Discover temperature files that match the specified date range and naming pattern, and raise an error if no files are found
    temp_files = discover_temperature_files(config.inpath_T, infile_spec, start_date, end_date)
    if config.max_files is not None:
        temp_files = temp_files[: config.max_files]
    if not temp_files:
        raise FileNotFoundError("No temperature files matched the requested date range and naming pattern.")

    # Load static fields
    elevation_full, ds_elev = open_array(config.infile_Z, config.Z_varin)
    landmask_full, ds_land = open_array(config.infile_landmask, config.landmask_varin)

    # Temperature defines the target grid. Read one file to get lat/lon.
    _first_temp, ds_first_temp = open_array(temp_files[0], config.T_varin)
    output_lat_full, output_lon_full, output_lat_name, output_lon_name = get_lat_lon(ds_first_temp)
    output_lat_attrs = dict(ds_first_temp[output_lat_name].attrs)
    output_lon_attrs = dict(ds_first_temp[output_lon_name].attrs)
    ds_first_temp.close()
    del _first_temp

    # Compute final output clip window on the temperature grid, then expand by 1
    # cell to obtain the internal compute window needed by the 3x3 stencil.
    clip_window = None
    compute_window = None
    if config.clip_to_reference_bbox:
        if not config.reference_grid_file:
            raise ValueError("reference_grid_file is required when clip_to_reference_bbox is enabled.")
        lat_dims = (output_lat_name,) if output_lat_full.ndim == 1 else ("y", "x")
        lon_dims = (output_lon_name,) if output_lon_full.ndim == 1 else ("y", "x")
        ds_temp_grid = xr.Dataset(
            {
                output_lat_name: (lat_dims, output_lat_full),
                output_lon_name: (lon_dims, output_lon_full),
            }
        )
        clip_window = compute_clip_window_from_reference_bbox(
            ds_like=ds_temp_grid,
            reference_grid_file=config.reference_grid_file,
            reference_lat_var=config.reference_lat_var,
            reference_lon_var=config.reference_lon_var,
            pad_pixels=config.bbox_pad_pixels,
        )
        del ds_temp_grid
        LOGGER.info("Reference-bbox clip window (temp grid): %s", clip_window)

        if output_lat_full.ndim == 2:
            nrows_full, ncols_full = output_lat_full.shape
        else:
            nrows_full, ncols_full = output_lat_full.size, output_lon_full.size

        if config.bbox_pad_pixels > 0:
            compute_window = expand_clip_window(
                clip_window,
                expand_pixels=1,
                max_rows=nrows_full,
                max_cols=ncols_full,
            )
            LOGGER.info("Computation window (expanded by 1 for stencil halo): %s", compute_window)
        else:
            compute_window = clip_window
            LOGGER.info("No bbox padding requested; computing directly on clip window.")

    # Final output coordinates
    if clip_window is not None:
        output_lat = (
            output_lat_full[clip_window.j0 : clip_window.j1 + 1]
            if output_lat_full.ndim == 1
            else slice_2d(output_lat_full, clip_window)
        )
        output_lon = (
            output_lon_full[clip_window.i0 : clip_window.i1 + 1]
            if output_lon_full.ndim == 1
            else slice_2d(output_lon_full, clip_window)
        )
    else:
        output_lat = output_lat_full
        output_lon = output_lon_full

    # Internal compute-grid coordinates
    if compute_window is not None:
        tgt_lat = (
            output_lat_full[compute_window.j0 : compute_window.j1 + 1]
            if output_lat_full.ndim == 1
            else slice_2d(output_lat_full, compute_window)
        )
        tgt_lon = (
            output_lon_full[compute_window.i0 : compute_window.i1 + 1]
            if output_lon_full.ndim == 1
            else slice_2d(output_lon_full, compute_window)
        )
    else:
        tgt_lat = output_lat_full
        tgt_lon = output_lon_full

    # Subset static fields near the compute target, then resample onto the
    # compute grid.
    tgt_lat_min = float(np.nanmin(tgt_lat))
    tgt_lat_max = float(np.nanmax(tgt_lat))
    tgt_lon_n = normalize_lon180(tgt_lon)
    tgt_lon_min = float(np.nanmin(tgt_lon_n))
    tgt_lon_max = float(np.nanmax(tgt_lon_n))

    elev_pre_window = compute_clip_window_from_bbox(
        ds_elev, tgt_lat_min, tgt_lat_max, tgt_lon_min, tgt_lon_max, pad_pixels=2
    )
    land_pre_window = compute_clip_window_from_bbox(
        ds_land, tgt_lat_min, tgt_lat_max, tgt_lon_min, tgt_lon_max, pad_pixels=2
    )

    elevation_sub, ds_elev_sub = _subset_array_to_window(elevation_full, ds_elev, elev_pre_window)
    landmask_sub, ds_land_sub = _subset_array_to_window(landmask_full, ds_land, land_pre_window)

    elevation = resample_to_target_grid(elevation_sub, ds_elev_sub, tgt_lat, tgt_lon)
    landmask = resample_to_target_grid(landmask_sub, ds_land_sub, tgt_lat, tgt_lon)
    ds_elev_sub.close()
    ds_land_sub.close()

    outpath = Path(config.outpath)
    outpath.mkdir(parents=True, exist_ok=True)
    save_km = bool(config.save_units_as_per_km)
    fill_water = not config.no_fill_water_mask
    extra_attrs = _get_extra_attrs(config)
    daily_groups: dict[pd.Timestamp, list[tuple[datetime, np.ndarray]]] = defaultdict(list)
    current_day_key: pd.Timestamp | None = None
    time_attrs: dict = {}

    pre_dfmt = config.preexisting_lr_date_format or config.infile_date_format
    pre_spec = (
        FileSpec(
            prefix=config.preexisting_lr_pfx or "",
            suffix=config.preexisting_lr_sfx or "",
            date_format=pre_dfmt,
        )
        if config.prefer_preexisting
        and config.preexisting_lr_dir
        and config.preexisting_lr_pfx
        and config.preexisting_lr_sfx
        else None
    )

    # landmask on the final output domain for pre-existing LR imprinting
    if clip_window is not None and compute_window is not None:
        row_off = clip_window.j0 - compute_window.j0
        col_off = clip_window.i0 - compute_window.i0
        nrows_out = clip_window.j1 - clip_window.j0 + 1
        ncols_out = clip_window.i1 - clip_window.i0 + 1
        landmask_out = landmask[row_off : row_off + nrows_out, col_off : col_off + ncols_out]
    else:
        landmask_out = landmask

    # Process each temperature file, extracting the temperature data and corresponding coordinates, applying clipping if specified, and computing the lapse rate grid for each time slice,
    # saving the output according to the specified output mode and date format.
    for file_path in temp_files:
        LOGGER.info("Processing %s", file_path)

        # --------------------------------------------------
        # Pre-check skip logic BEFORE reading file
        # --------------------------------------------------
        if pre_spec is not None:
            try:
                file_date = parse_date_token(file_path.name, config.infile_date_format)
            except ValueError:
                file_date = None

            if file_date is not None:
                pre_name = build_dated_filename(pre_spec, file_date)
                pre_path = Path(config.preexisting_lr_dir) / pre_name  # type: ignore[arg-type]
                if pre_path.exists():
                    LOGGER.info("Using pre-existing LR file: %s", pre_path)
                    lr_pre, time_vals_pre, lat_pre, lon_pre, lat_name_pre, lon_name_pre, lat_attrs_pre, lon_attrs_pre = _read_and_clip_preexisting(
                        pre_path=pre_path,
                        pre_var=config.preexisting_lr_var,
                        pre_time_var=config.preexisting_time_var,
                        clip_window=clip_window,
                        landmask=landmask_out,
                        fallback_lapse_rate_m=config.LR_forced_val,
                        save_units_as_per_km=save_km,
                        min_lapse_rate_cutoff=config.min_lapse_rate_cutoff,
                        max_lapse_rate_cutoff=config.max_lapse_rate_cutoff,
                        buffer_width=config.buffer_width,
                        fill_water_mask=fill_water,
                        landmask_trueval=config.landmask_trueval,
                    )
                    for t_idx, t_dt in enumerate(time_vals_pre):
                        lapse_2d = lr_pre[t_idx]
                        if output_mode == "hourly":
                            _save_hourly(
                                lapse=lapse_2d,
                                valid_time=t_dt,
                                lat=lat_pre,
                                lon=lon_pre,
                                lat_name=lat_name_pre,
                                lon_name=lon_name_pre,
                                outfile_spec=outfile_spec,
                                outpath=outpath,
                                save_km=save_km,
                                force=config.force,
                                dryrun=config.dryrun,
                                extra_attrs=extra_attrs,
                                lat_attrs_in=lat_attrs_pre,
                                lon_attrs_in=lon_attrs_pre,
                            )
                        else:
                            day_key = pd.Timestamp(t_dt).normalize()
                            daily_groups[day_key].append((t_dt, lapse_2d))
                    continue
                else:
                    LOGGER.info("Pre-existing LR file not found (%s); falling back to compute.", pre_path)

        if output_mode == "daily":
            try:
                file_date = parse_date_token(file_path.name, config.infile_date_format)
                day_key = pd.Timestamp(file_date).normalize()
                output_name = build_dated_filename(outfile_spec, day_key.to_pydatetime())
                output_path = outpath / output_name

                if output_path.exists() and not config.force:
                    LOGGER.info("Skipping existing daily file %s", output_path)
                    continue

                if config.dryrun:
                    LOGGER.info("[Dry run] Would write %s", output_path)
                    continue
            except ValueError:
                pass


        # Extract temperature data and corresponding coordinates from the current file,
        # applying clipping to a reference bounding box if specified in the configuration,
        # and check for expected units
        temp_values, ds_temp = open_array(file_path, config.T_varin)

        # Hourly merged MERRA2bias files store T2M as a 2-D scalar-time field.
        # Normalize to (time, lat, lon) so the downstream loop can index by time.
        if temp_values.ndim == 2:
            temp_values = temp_values[np.newaxis, :, :]
        elif temp_values.ndim != 3:
            ds_temp.close()
            raise ValueError(
                f"Expected 2-D or 3-D temperature array, got shape {temp_values.shape}"
            )

        raw_t_units = ds_temp[config.T_varin].attrs.get("units", "").strip().lower()
        if raw_t_units not in {"", "k", "kelvin", "kelvins", "degree_kelvin", "degk"}:
            LOGGER.warning(
                "Unexpected units '%s' for '%s' in %s; expected Kelvin.",
                raw_t_units,
                config.T_varin,
                file_path,
            )

        # Extract time values and attributes.
        time_values, time_attrs = get_time_values(
            ds_temp,
            config.time_varin,
            time_shift_minutes=config.time_shift_minutes,
        )
        ds_temp.close()

        mismatch = np.isnan(elevation) != np.isnan(landmask)
        if np.any(mismatch):
            LOGGER.warning(
                "%d cells have mismatched NaNs between elevation and landmask.",
                int(np.sum(mismatch)),
            )

        for index, valid_time in enumerate(time_values[: temp_values.shape[0]]):
            valid_time_bucket = _floor_to_hour(valid_time)

            # If an hour filter was provided, keep only timesteps whose containing
            # hour falls within the requested window. Exact minutes are ignored.
            if start_hour_filter is not None and end_hour_filter is not None:
                if valid_time_bucket < start_hour_filter or valid_time_bucket > end_hour_filter:
                    LOGGER.debug(
                        "Skipping timestep %s because its hour bucket %s is outside requested window.",
                        valid_time,
                        valid_time_bucket,
                    )
                    continue

            temp_slice = temp_values[index]
            if compute_window is not None:
                temp_slice = slice_2d(temp_slice, compute_window)

            if output_mode == "hourly":
                output_name = build_dated_filename(outfile_spec, valid_time_bucket)
                output_path = outpath / output_name
                if output_path.exists() and not config.force:
                    LOGGER.info("Skipping existing file %s", output_path)
                    continue
                if config.dryrun:
                    LOGGER.info("[Dry run] Would write %s", output_path)
                    continue
            
            # Compute the lapse rate grid for the current time slice using the sampled elevation and landmask data on the target grid, applying the specified parameters for land neighbor
            # requirements, fallback lapse rate, cutoffs, and water mask filling.
            lapse = compute_lapse_rate_grid(
                temperature=temp_slice,
                elevation=elevation,
                landmask=landmask,
                landmask_trueval=config.landmask_trueval,
                min_land_neighbors=config.min_land_neighbors,
                fallback_lapse_rate_m=config.LR_forced_val,
                min_lapse_rate_cutoff=config.min_lapse_rate_cutoff,
                max_lapse_rate_cutoff=config.max_lapse_rate_cutoff,
                save_units_as_per_km=save_km,
                buffer_width=config.buffer_width,
                fill_water_mask=fill_water,
                min_elev_diff_std=config.min_elev_diff_std,
            )

            # If clipping is enabled, apply the clip window to the computed lapse rate grid to obtain the
            # final output grid that matches the desired output domain, which may be smaller than the
            if clip_window is not None and compute_window is not None:
                row_off = clip_window.j0 - compute_window.j0
                col_off = clip_window.i0 - compute_window.i0
                nrows_out = clip_window.j1 - clip_window.j0 + 1
                ncols_out = clip_window.i1 - clip_window.i0 + 1
                lapse = lapse[row_off : row_off + nrows_out, col_off : col_off + ncols_out]
            
            # Depending on the output mode, either save the lapse rate grid for the current time slice immediately as an hourly file, or append it to the list of daily slices for later stacking
            # and saving as a single daily file.
            if output_mode == "hourly":
                _save_hourly(
                    lapse=lapse,
                    valid_time=valid_time_bucket,
                    lat=output_lat,
                    lon=output_lon,
                    lat_name=output_lat_name,
                    lon_name=output_lon_name,
                    outfile_spec=outfile_spec,
                    outpath=outpath,
                    save_km=save_km,
                    force=config.force,
                    dryrun=config.dryrun,
                    time_attrs=time_attrs,
                    extra_attrs=extra_attrs,
                    lat_attrs_in=output_lat_attrs,
                    lon_attrs_in=output_lon_attrs,
                )
            else:
                day_key = pd.Timestamp(valid_time).normalize()
                if current_day_key is not None and day_key != current_day_key:
                    _flush_daily_group(
                        daily_groups, current_day_key, outfile_spec, outpath,
                        output_lat, output_lon, output_lat_name, output_lon_name,
                        save_km, config.force, config.dryrun,
                        time_attrs, extra_attrs, output_lat_attrs, output_lon_attrs,
                    )
                    daily_groups.pop(current_day_key, None)
                    gc.collect()
                current_day_key = day_key
                daily_groups[day_key].append((valid_time, lapse))

    # After processing all files and time slices, if the output mode is daily,
    # iterate over the grouped daily slices, stack them into a 3-D array, and save each daily file
    # with the appropriate metadata and coordinate variables.
    if output_mode == "daily" and current_day_key is not None:
        _flush_daily_group(
            daily_groups, current_day_key, outfile_spec, outpath,
            output_lat, output_lon, output_lat_name, output_lon_name,
            save_km, config.force, config.dryrun,
            time_attrs, extra_attrs, output_lat_attrs, output_lon_attrs,
        )

    ds_elev.close()
    ds_land.close()


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """ Build and return an ArgumentParser for the command-line interface, defining all necessary arguments and options for running the dynamic lapse rate computation pipeline,
    including input file paths, variable names, date range, output settings, and other parameters.
    
    Returns
    -------
    argparse.ArgumentParser
        An ArgumentParser object configured with all necessary arguments and options for running the dynamic lapse rate computation pipeline from the command line.
    """
    parser = argparse.ArgumentParser(
        description="Compute dynamic temperature lapse rates from gridded meteorological inputs."
    )

    parser.add_argument("--infile_Z", required=True, help="Elevation NetCDF input file.")
    parser.add_argument("--infile_landmask", required=True, help="Land-mask NetCDF input file.")
    parser.add_argument("--inpath_T", required=True, help="Directory containing temperature NetCDF files.")
    parser.add_argument("--outpath", required=True, help="Output directory.")
    parser.add_argument("--yyyymmdd_start", required=True, help="Start date in YYYYMMDD format.")
    parser.add_argument("--yyyymmdd_end", required=True, help="End date in YYYYMMDD format.")
    parser.add_argument(
        "--start_datetime",
        default=None,
        help="Optional start datetime in 'YYYY-MM-DD HH:MM'. Hour filtering ignores minutes and uses the containing hour.",
    )
    parser.add_argument(
        "--end_datetime",
        default=None,
        help="Optional end datetime in 'YYYY-MM-DD HH:MM'. Hour filtering ignores minutes and uses the containing hour.",
    )
    parser.add_argument(
        "--time_shift_minutes",
        type=int,
        default=0,
        help=(
            "Optional uniform shift, in minutes, applied to input timestamps before "
            "hour filtering and output writing. Default is 0 (no shift). "
            "Example: --time_shift_minutes -30 shifts 06:30 to 06:00."
        ),
    )
    
    parser.add_argument("--Z_varin", default="ELEVATION", help="Elevation variable name.")
    parser.add_argument("--T_varin", default="T2M", help="Temperature variable name.")
    parser.add_argument("--landmask_varin", default="LANDMASK", help="Land-mask variable name.")
    parser.add_argument("--time_varin", default="time", help="Time variable name in temperature files.")

    parser.add_argument("--fname_in_pfx_T", default="", help="Temperature input filename prefix.")
    parser.add_argument("--fname_in_sfx_T", default=".nc", help="Temperature input filename suffix.")
    parser.add_argument("--fname_out_pfx", default="lapse_rate.", help="Output filename prefix.")
    parser.add_argument("--fname_out_sfx", default=".nc", help="Output filename suffix.")
    parser.add_argument(
        "--infile_date_format",
        default="YYYYMMDD",
        choices=["YYYYMMDD", "YYYY-MM-DD", "YYYY-MM-DDTHH", "YYYY-MM-DD-HH"],
        help="Date-token format used in temperature input filenames.",
    )
    parser.add_argument(
        "--outfile_date_format",
        default="YYYYMMDD",
        choices=["YYYYMMDD", "YYYY-MM-DD", "YYYY-MM-DDTHH", "YYYY-MM-DD-HH"],
        help="Date-token format used in output filenames.",
    )
    parser.add_argument(
        "--output_mode",
        default="daily",
        choices=["daily", "hourly", "match_input"],
        help=(
            "Controls output granularity. 'daily' saves one file per day containing all time steps; "
            "'hourly' saves one file per hour containing a single 2-D slice; "
            "'match_input' infers the mode from --infile_date_format."
        ),
    )

    parser.add_argument("--landmask_trueval", type=float, default=1, help="Value that indicates land in the land mask.")
    parser.add_argument("--min_land_neighbors", type=int, default=5, help="Minimum number of valid neighbouring land cells required in the 3 x 3 window.")
    parser.add_argument("--LR_forced_val", type=float, default=-0.0065, help="Fallback static lapse rate in K/m.")
    parser.add_argument(
        "--min_lapse_rate_cutoff",
        type=float,
        default=None,
        help="Optional minimum allowed lapse rate. Interpreted in K/km when --save_units_as_per_km=1, otherwise in K/m.",
    )
    parser.add_argument(
        "--max_lapse_rate_cutoff",
        type=float,
        default=None,
        help="Optional maximum allowed lapse rate. Interpreted in K/km when --save_units_as_per_km=1, otherwise in K/m.",
    )
    parser.add_argument("--save_units_as_per_km", type=int, default=1, choices=[0, 1], help="Save output in K/km if 1, K/m if 0.")
    parser.add_argument("--buffer_width", type=int, default=1, help="Width (in source-grid pixels) of the land-adjacent water buffer.")
    parser.add_argument("--no_fill_water_mask", action="store_true", help="Disable shoreline and inland-water imprinting.")

    parser.add_argument("--clip_to_reference_bbox", action="store_true", help="Clip outputs to the bounding box of a reference lat/lon file.")
    parser.add_argument("--reference_grid_file", default=None, help="Reference NetCDF file whose latitude/longitude bounding box is used to clip all inputs.")
    parser.add_argument("--reference_lat_var", default=None, help="Latitude variable name in the reference file.")
    parser.add_argument("--reference_lon_var", default=None, help="Longitude variable name in the reference file.")
    parser.add_argument("--bbox_pad_pixels", type=int, default=0, help="Expand the reference bounding-box clip window by this many source pixels.")
    parser.add_argument("--min_elev_diff_std", type=float, default=None,
        help="Optional elevation-variability threshold (metres).")
    parser.add_argument("--prefer_preexisting", action="store_true")
    parser.add_argument("--preexisting_lr_dir", default=None)
    parser.add_argument("--preexisting_lr_pfx", default=None)
    parser.add_argument("--preexisting_lr_sfx", default=None)
    parser.add_argument("--preexisting_lr_var", default="lapse_rate")
    parser.add_argument("--preexisting_lr_date_format", default=None)
    parser.add_argument("--preexisting_time_var", default="time")
    parser.add_argument("--max_files", type=int, default=None)
    parser.add_argument("--force", action="store_true", help="Overwrite existing output files.")
    parser.add_argument("--dryrun", action="store_true", help="Log planned outputs without writing any files.")
    parser.add_argument("--debug", action="store_true", help="Enable verbose logging.")

    parser.add_argument("--output_source_data", default=None, help="Optional global attribute: name of the source dataset.")
    parser.add_argument("--output_source_variable", default=None, help="Optional global attribute: source variable name.")
    parser.add_argument("--output_source_variable_description", default=None, help="Optional global attribute: description of the source variable.")
    parser.add_argument("--output_input_description", default=None, help="Optional global attribute: description of inputs used to generate lapse rates.")
    parser.add_argument("--output_method_reference", default=None, help="Optional global attribute: method reference.")

    return parser


def main() -> None:
    """Main entry point for the dynamic lapse rate computation script, handling argument parsing, logging configuration, and pipeline execution."""
    import warnings

    # Suppress specific warnings that may arise from underlying libraries (e.g., PROJ data directory not found, gini engine loading failure) to avoid
    # cluttering the output.
    warnings.filterwarnings("ignore", message="Valid PROJ data directory not found")
    warnings.filterwarnings("ignore", message=".*Engine 'gini' loading failed.*")

    # Parse command-line arguments and configure logging based on the debug flag. Then, create a Config object from the parsed arguments and run the main
    # pipeline to compute dynamic lapse rates and save the outputs.
    parser = build_parser()
    args = parser.parse_args()
    configure_logging(debug=args.debug)

    config = Config(
        infile_Z=args.infile_Z,
        infile_landmask=args.infile_landmask,
        inpath_T=args.inpath_T,
        outpath=args.outpath,
        yyyymmdd_start=args.yyyymmdd_start,
        yyyymmdd_end=args.yyyymmdd_end,
        start_datetime=args.start_datetime,
        end_datetime=args.end_datetime,
        time_shift_minutes=args.time_shift_minutes,
        Z_varin=args.Z_varin,
        T_varin=args.T_varin,
        landmask_varin=args.landmask_varin,
        time_varin=args.time_varin,
        fname_in_pfx_T=args.fname_in_pfx_T,
        fname_in_sfx_T=args.fname_in_sfx_T,
        fname_out_pfx=args.fname_out_pfx,
        fname_out_sfx=args.fname_out_sfx,
        infile_date_format=args.infile_date_format,
        outfile_date_format=args.outfile_date_format,
        output_mode=args.output_mode,
        landmask_trueval=args.landmask_trueval,
        min_land_neighbors=args.min_land_neighbors,
        LR_forced_val=args.LR_forced_val,
        min_lapse_rate_cutoff=args.min_lapse_rate_cutoff,
        max_lapse_rate_cutoff=args.max_lapse_rate_cutoff,
        save_units_as_per_km=args.save_units_as_per_km,
        buffer_width=args.buffer_width,
        no_fill_water_mask=args.no_fill_water_mask,
        min_elev_diff_std=args.min_elev_diff_std,
        prefer_preexisting=args.prefer_preexisting,
        preexisting_lr_dir=args.preexisting_lr_dir,
        preexisting_lr_pfx=args.preexisting_lr_pfx,
        preexisting_lr_sfx=args.preexisting_lr_sfx,
        preexisting_lr_var=args.preexisting_lr_var,
        preexisting_lr_date_format=args.preexisting_lr_date_format,
        preexisting_time_var=args.preexisting_time_var,
        max_files=args.max_files,
        clip_to_reference_bbox=args.clip_to_reference_bbox,
        reference_grid_file=args.reference_grid_file,
        reference_lat_var=args.reference_lat_var,
        reference_lon_var=args.reference_lon_var,
        bbox_pad_pixels=args.bbox_pad_pixels,
        force=args.force,
        dryrun=args.dryrun,
        output_source_data=args.output_source_data,
        output_source_variable=args.output_source_variable,
        output_source_variable_description=args.output_source_variable_description,
        output_input_description=args.output_input_description,
        output_method_reference=args.output_method_reference,
    )

    LOGGER.debug("===== Parameter Summary =====")
    for key, value in vars(args).items():
        LOGGER.debug("  %s: %s", key, value)
    LOGGER.debug("=============================")

    run_pipeline(config)


if __name__ == "__main__":
    main()
