# NLDAS-3 Retrospective Precipitation Processing Codes

This directory contains Fortran utilities used to generate, downscale, rescale, assimilate, and temporally disaggregate precipitation forcing fields for the NLDAS-3 retrospective processing workflow.

## Codes

### `cloud_frequencies.f90`

Computes monthly MODIS cloud-frequency fields at 1 km resolution using the MODIS `state_1km_1` quality flag. The code extracts the MODIS internal cloud flag, masks ocean and ocean-adjacent coastline pixels, averages valid daily cloud flags over each month, and writes monthly NetCDF cloud-frequency files.

### `downscaled_merra2_cloud_4km.f90`

Downscales daily MERRA-2/LIS precipitation fields using monthly MODIS Terra and Aqua cloud-frequency information. The code combines MODIS and MYD cloud-frequency fields, aggregates precipitation and cloud frequency to a coarse grid, and redistributes daily precipitation back to the 4 km grid.

### `downscaled_nldas3_cloud_1km.f90`

Downscales daily NLDAS-3 precipitation from 4 km to 1 km using MODIS Terra and Aqua cloud-frequency fields. Missing 4 km NLDAS-3 precipitation values are filled using MERRA-2, and the final output is written as daily 1 km NetCDF precipitation files.

### `rescale_imergv07.f90`

Applies a GPCP-based bias correction to daily IMERG V07 precipitation. The code multiplies each daily IMERG grid cell by a monthly GPCP-derived scale factor and writes daily corrected precipitation NetCDF files.

### `semivar.f90`

Computes empirical semivariograms of daily precipitation errors between gridded precipitation fields and station observations. The output text files provide distance bins, semivariogram values, and pair counts for use in precipitation assimilation.

### `code_bratseth_nldas3.f90`

Performs MPI-parallel precipitation assimilation using a BRASH analysis framework, which is based on the Bratseth algorithm and follows the precipitation merging approach of Kemp et al. (2022). The code assimilates IMERG and CAPA precipitation estimates into a MERRA-2/LIS background field and writes daily analyzed precipitation fields on the NLDAS-3 grid.

### `temporal_disaggregation_1km.f90`

Converts daily 1 km NLDAS-3 precipitation into hourly precipitation fields. The code uses hourly MERRA-2 and IMERG precipitation fractions to distribute each daily 1 km precipitation total across 24 hourly fields.

## General Compilation

Most serial codes can be compiled on NASA Discover using the Intel Fortran compiler and the LISF NetCDF/HDF5 libraries:

```bash
ifort -g -check all -traceback -names lowercase -convert big_endian \
  -assume byterecl \
  -I/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/include \
  <code_name>.f90 -o <executable_name> \
  -L/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/lib \
  -L/discover/nobackup/projects/lis/libs/sles-12.3/hdf5/1.12.1_intel-2021.4.0/lib \
  -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl \
  -Wl,--no-relax -shared-intel
```

For MPI-based codes such as `code_bratseth_nldas3.f90`, use `mpif90` instead of `ifort`:

```bash
mpif90 -g -check all -traceback -names lowercase -convert big_endian \
  -assume byterecl \
  -I/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/include \
  code_bratseth_nldas3.f90 -o code_bratseth_nldas3 \
  -L/discover/nobackup/projects/lis/libs/sles-12.3/netcdf/4.8.1_intel-2021.4.0/lib \
  -L/discover/nobackup/projects/lis/libs/sles-12.3/hdf5/1.12.1_intel-2021.4.0/lib \
  -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl \
  -Wl,--no-relax -shared-intel
```

## Notes

Input and output directories are defined in the `USER SETTINGS` block near the top of each source file. Users should modify only that block unless changes to the algorithm or grid configuration are required.

The processing workflow assumes that input files follow the expected LIS-style naming convention, such as:

```text
LIS_HIST_YYYYMMDDHH00.d01.nc
```
