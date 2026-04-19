"""
SCRIPT: traditional_predictor.py

Traditional physical algorithms for snow depth retrieval from AMSR2 data.
  - Foster (Yong et al. 2022)
  - Kelly, 2009

REVISION HISTORY:
"""
# pylint: disable=import-error
# pylint: disable=invalid-name, too-many-locals
# Standard modules
import logging
import os

# Third party modules
import numpy as np
import xarray as xr
from rasterio.enums import Resampling

from .run_prediction import SnowDepthPredictor

logger = logging.getLogger('SnowDepthPredictor')

# Encoding used for all NetCDF outputs
_ENCODING = {
    'snow_depth': {
        'zlib': True,
        'complevel': 4,
        'shuffle': True,
        '_FillValue': -9999.0,
        'dtype': 'float32'
    },
    'y': {'dtype': 'float32'},
    'x': {'dtype': 'float32'}
}


class SnowDepthPredictorTraditional(SnowDepthPredictor):
    """
    Extends SnowDepthPredictor with traditional physical algorithms
    (Foster 2022, Kelly) instead of an ML model.
    """

    def __init__(self, config=None):
        super().__init__(config=config)

    # ──────────────────────────────────────────────────────────────
    # Algorithm implementations — pure calculation, no I/O
    # ──────────────────────────────────────────────────────────────
    @staticmethod
    def _extract_yx(da: xr.DataArray):
        """
        Extract 1-D y and x arrays from a DataArray.
        Handles both 'lat'/'lon' and 'y'/'x' dimension names.

        Returns:
            y     : 1-D numpy array of y (latitude) values
            x     : 1-D numpy array of x (longitude) values
            y_dim : dimension name used for y ('lat' or 'y')
            x_dim : dimension name used for x ('lon' or 'x')
        """
        if 'y' in da.dims and 'x' in da.dims:
            y = da['y'].values
            x = da['x'].values
            y_dim = 'y'
            x_dim = 'x'
        elif 'lat' in da.dims and 'lon' in da.dims:
            y = da['lat'].values
            x = da['lon'].values
            y_dim = 'lat'
            x_dim = 'lon'
        else:
            raise ValueError(
                f'Cannot find y/x or lat/lon dims in DataArray. '
                f'Available dims: {da.dims}'
            )
        return y, x, y_dim, x_dim

    def cal_foster_SD(self, ds_pmw):
        """
        Foster algorithm (Yong et al. 2022):
            SD = 1.59 * (Tb18H - Tb36H) / (1 - ff)
        ff: fractional forest cover from MCD12Q1 (band_data)
        """

        data = ds_pmw
        ds_ff = xr.open_dataset(self.config.fraction_forest_cover).squeeze()

        # reproject forest variables to data
        data = data.rio.write_crs(self.config.proj)
        ds_ff = ds_ff.rio.write_crs(self.config.proj)
        ds_ff_repo = ds_ff.rio.reproject_match(
            data,
            resampling=Resampling.nearest
        )

        # rename spatial dims to match data if needed
        if 'y' in ds_ff_repo.dims and 'lat' in data.dims:
            data = data.rename({'lat': 'y', 'lon': 'x'})

        ff = ds_ff_repo['band_data'].squeeze()

        ff = xr.where(ff < 0, np.nan, ff)
        ff = xr.where(ff > 1, 1.0, ff)

        tb_diff = data['tb_18h'] - data['tb_36h']
        tb_diff = xr.where(tb_diff < 0, np.nan, tb_diff)

        denominator = 1 - ff
        denominator = xr.where(denominator <= 0, np.nan, denominator)

        sd = 1.59 * tb_diff / denominator / 100  # cm to meter
        sd = xr.where(sd < 0, 0.0, sd)  # No negative snow depth
        sd = xr.where(sd > 5, np.nan, sd)  # Max realistic snow depth 5 cm

        # assign coords
        sd = sd.squeeze('time') if 'time' in sd.dims else sd
        y, x, y_dim, x_dim = self._extract_yx(data['tb_18h'])
        sd = sd.rename({y_dim: 'y', x_dim: 'x'})
        sd = sd.assign_coords(
            y=('y', y),
            x=('x', x)
        )
        sd.coords['y'].attrs = {
            'units': 'degrees_north',
            'long_name': 'latitude',
            'standard_name': 'latitude',
            'axis': 'Y'
        }
        sd.coords['x'].attrs = {
            'units': 'degrees_east',
            'long_name': 'longitude',
            'standard_name': 'longitude',
            'axis': 'X'
        }
        logger.info('Foster SD calculated.')
        return sd

    def cal_kelly_SD(self, ds_pmw):
        """
        Kelly algorithm:
            SDf = 1/log10(pol36) * (tb18v - tb36v) / (1 - fd*0.6)
            SD0 = 1/log10(pol36) * (tb10v - tb36v)
                + 1/log10(pol18) * (tb10v - tb18v)
            SD  = ff * SDf + (1 - ff) * SD0

        ff: fractional forest cover from MCD12Q1
        fd: forest density from MOD44B
        """
        data = ds_pmw
        ds_ff = xr.open_dataset(self.config.fraction_forest_cover).squeeze()
        ds_fd = xr.open_dataset(self.config.forest_density).squeeze()

        # reproject forest variables to data
        data = data.rio.write_crs(self.config.proj)
        ds_ff = ds_ff.rio.write_crs(self.config.proj)
        ds_fd = ds_fd.rio.write_crs(self.config.proj)
        ds_ff_repo = ds_ff.rio.reproject_match(
            data,
            resampling=Resampling.nearest
        )
        ds_fd_repo = ds_fd.rio.reproject_match(
            data,
            resampling=Resampling.nearest
        )

        if 'y' in ds_ff_repo.dims and 'lat' in data.dims:
            data = data.rename({'lat': 'y', 'lon': 'x'})

        ff = ds_ff_repo['band_data'].squeeze()
        fd = ds_fd_repo['band_data'].squeeze()

        ff = xr.where(ff < 0, np.nan, ff)
        ff = xr.where(ff > 1, 1.0, ff)

        fd = xr.where(fd < 0, np.nan, fd)
        fd = xr.where(fd > 1, 1.0, fd)

        pol36 = data['tb_36v'] - data['tb_36h']
        pol18 = data['tb_18v'] - data['tb_18h']
        pol36 = xr.where(pol36 > 0, pol36, np.nan)
        pol18 = xr.where(pol18 > 0, pol18, np.nan)

        denominator = 1 - fd * 0.6
        denominator = xr.where(denominator <= 0, np.nan, denominator)

        SDf = (1 / np.log10(pol36) *
               (data['tb_18v'] - data['tb_36v']) /
               denominator)

        SD0 = ((1 / np.log10(pol36) *
                (data['tb_10v'] - data['tb_36v'])) +
               (1 / np.log10(pol18) *
                (data['tb_10v'] - data['tb_18v'])))

        sd = (ff * SDf + (1 - ff) * SD0) / 100  # cm to meter

        sd = xr.where(sd < 0, 0.0, sd)
        sd = xr.where(sd > 5, np.nan, sd)

        # assign coords
        sd = sd.squeeze('time') if 'time' in sd.dims else sd
        y, x, y_dim, x_dim = self._extract_yx(data['tb_18h'])
        sd = sd.rename({y_dim: 'y', x_dim: 'x'})
        sd = sd.assign_coords(
            y=('y', y),
            x=('x', x)
        )
        sd.coords['y'].attrs = {
            'units': 'degrees_north',
            'long_name': 'latitude',
            'standard_name': 'latitude',
            'axis': 'Y'
        }
        sd.coords['x'].attrs = {
            'units': 'degrees_east',
            'long_name': 'longitude',
            'standard_name': 'longitude',
            'axis': 'X'
        }

        logger.info('Kelly SD calculated.')
        return sd

    # ──────────────────────────────────────────────────────────────
    # Save helper — shared by both algorithms
    # ──────────────────────────────────────────────────────────────

    def _save_sd(self, sd: xr.DataArray, method: str) -> str:
        """
        Save a snow depth DataArray to NetCDF.

        Args:
            sd     : snow depth DataArray returned by a cal_*_SD method
            method : 'foster' or 'kelly' — used in the output filename

        Returns:
            output_file path
        """

        target_datetime = self.target_datetime.strftime("%Y%m%d%H")
        dir_out = self.config.project_path / self.config.output_dir
        os.makedirs(dir_out, exist_ok=True)

        output_file = os.path.join(
            dir_out,
            f'amsr2_snip_0p1deg_{target_datetime}_{method}.nc'
        )

        # apply mask from ML to these two methods
        mask_file = os.path.join(
            dir_out,
            f'amsr2_snip_0p1deg_{target_datetime}.nc'
        )
        if os.path.exists(mask_file):
            ds_mask = xr.open_dataset(mask_file)
            mask_sd = ds_mask['snow_depth'].squeeze(
                [d for d in ds_mask['snow_depth'].dims
                 if ds_mask['snow_depth'].sizes[d] == 1],
                drop=True
            )
            mask_sd = mask_sd.assign_coords(
                y=sd['y'],
                x=sd['x']
            )
            sd_coords = sd.coords
            sd = xr.where(mask_sd >= 0, sd, np.nan, keep_attrs=True)
            sd = sd.assign_coords(sd_coords)
            logger.info('ML mask applied from %s', mask_file)
        else:
            logger.warning(
                'Mask file not found: %s — saving without mask', mask_file
            )

        # Wrap in a Dataset so encoding keys match variable names
        ds_out = sd.to_dataset(name='snow_depth')
        ds_out.to_netcdf(output_file, encoding=_ENCODING, format='NETCDF4')
        logger.info('%s output saved to %s',
                    method.capitalize(), output_file)
        return output_file

    # ──────────────────────────────────────────────────────────────
    # Public entry point
    # ──────────────────────────────────────────────────────────────

    def predict(self, method: str = 'foster', save: bool = True):
        """
        Run the selected traditional algorithm and optionally save to NetCDF.

        Args:
            method : 'foster' or 'kelly'
            save   : if True (default) write output to NetCDF

        Returns:
            xarray.DataArray of predicted snow depth
        """
        if self.pmw_file is None:
            raise RuntimeError(
                'pmw_file is not set. Assign the merged PMW file '
                'to predictor.ds_result before calling predict().'
            )

        pmw_file = self.pmw_file
        with xr.open_dataset(pmw_file) as ds_pmw:
            logger.info('%s file opened', pmw_file)

        if method == 'foster':
            logger.info('Running Foster algorithm ...')
            sd = self.cal_foster_SD(ds_pmw)

        elif method == 'kelly':
            logger.info('Running Kelly algorithm ...')
            sd = self.cal_kelly_SD(ds_pmw)
        else:
            raise ValueError(
                f"Unknown method '{method}'. Choose 'foster' or 'kelly'."
            )

        if save:
            self._save_sd(sd, method)
            logger.info("Save %s SD results", method)

        return sd
