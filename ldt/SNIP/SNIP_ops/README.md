# USAF Snow & Ice Product (SNIP) sub-System

An AI/ML-enhanced snow and ice data production subsystem that processes AMSR2 L1R satellite data to generate next-generation snow and ice products for U.S. Air Force (USAF) operational use. This repository contains the Phase 0 implementation of the Snow and Ice Products (SNIP) subsystem. It  produces snow depth estimates for the Phase 1: LDT.

## Overview

The AMSR2 SNIP system addresses the need for accurate, automated snow depth mapping from JAXA's GCOM-W/AMSR2 microwave radiometer data. This system:
- Developed as part of NASA's next-generation Earth observation capabilities, leveraging AMSR2 L1R brightness temperature data for enhanced snow depth estimation
- Provides reliable, near real-time snow depth estimation with improved accuracy over traditional methods
- Python-based processing pipeline with satellite data reading, reprojection, and machine learning predictions
- Snow depth quality control and gap filling in LIS LDT module
- LIS LVT post processing

## Features

- ✅ Read and resample AMSR2 L1R data to geographic gridcells
- ✅ Save data quality flag including flags for cold desert, frozen ground, glaciers, rain, and C bands RFI 
- ✅ AI/ML-based snow depth prediction
- ✅ Snow depth prediction quality control and gap filling in LDT
- ✅ LVT post processing to write outputs to legacy GRIB-1 and GRIB-2 formats

## Python module: AI/ML snow depth retrieval sub-system
### Configuration Files

| File | Description |
|------|-------------|
| `config/SNIP_config.json` | Template configuration file containing project settings and parameters |
| `config/load_config.py` | Utility script to load and parse configuration file for current time step |

#### Basic configuration in SNIP_config.json
- **target_datetime**: Project running time with the format of YYYYMMDDHHMM. For example, "202501200600".
- **project_path**: Project root directory containing source code and configuration files
- **amsr2_path**: Primary AMSR2 L1R data storage location
- **amsr2_merge_path**: Output directory for merged AMSR2 data products (serves as input for ML prediction pipeline)
- **output_dir**: Final snow depth data output location and storage path for reprojected data products
- **model_path**: Pre-trained machine learning model directory path
- **template_path**:  Reference template data path used for spatial reprojection operations
- **AMSR2_template_path**: Reference template file used for snow depth reprojection retrieved using traditional approaches
- **fraction_forest_cover**: Fractional forest cover dataset ("data/template/MCD12Q1_LC_frac_2019_reprojected.tif")
- **forest_density**: Forest density ("data/template/MOD44B_FF_2019_reprojected.tif")
- **viirs_path**: Directory to the VIIRS snow cover files ("../data/input/viirs")
- **flag_output_kelly**: Boolean flag to control calculating snow depth using Kelly, (2009) approach (default: false)
- **flag_output_foster**: Boolean flag to control calculating snow depth using Foster et al., (2005) approach (default: false). The same approach is used by USAFSI.
- **apply_viirs_mask**: Boolean flag to control applying VIIRS snow mask (default: true)
- **reproject_USAF**: Boolean flag to control saving of reprojected data products (default: true)
- **source**: AMSR2 L1R data source specification - accepts "NOAA" or "JAXA"

## Data Processing 

| File | Description |
|------|-------------|
| `data_processing/amsr2_reader.py` | Core data processing module that:<br>• Reads AMSR2 L1R data files<br>• Reprojects data to common grids<br>• Processes data quality flags<br>• Merges multiple data extents<br>• Saves processed data to NetCDF format |

## Machine Learning Prediction 

| File                                        | Description                                                                                                                                                                                                                                                                                       |
|---------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `ml_prediction/run_prediction.py`           | ML inference pipeline that:<br>• Loads pretrained machine learning models<br>• Runs ML predictions on processed data<br>• Applies VIIRS Snow Covered Area (SCA) masking<br>• Handles data quality flagging<br>• Reprojects results to Analysis Framework (AF) grids<br>• Saves prediction results |
| `ml_prediction/run_traditional_approach.py` | Traditional snow depth retrival pipeline that:<br>• Calculate snow depth  using Kelly, (2009) approach and Foster et al, (2005) approach. <br>• Options provided in config file to output these two snow depth results. <br>• Saves prediction results at PMW native resolution        |


## Main Execution

| File     | Description                                                          |
|----------|----------------------------------------------------------------------|
| `main.py` | Main execution script that orchestrates the entire SNIP_ops workflow |
| `amsr2_reader_main.py`    | Main execution script to read, merge and reproject AMSR2 data        |


## Usage
### Activate the virtual environment 
```shell
conda activate afw-py311-202501
```
### Process single datetime
```shell
# cd to the script directory
python main.py ./config/SNIP_config.json
```

## LIS LDT: snow depth prediction quality control and gap filling
This use case runs LDT to generate SNIP snow and ice analyses for the 557WW ~10-km global domain.

### Configuration Files

| File         | Description |
|--------------|-------------|
| `ldt.config` | Template configuration file containing project settings and parameters |

## Usage
* Generate the LDT executable following the build instructions with the
  software.  Be sure to select "little endian".
* Run the LDT executable using the lis_input.nrt_streamflow.noah39.nc file
  and the input datasets.
* Output netCDF files will be in the ../data/output/snip directory with the log file saved in ../data/output/logs

## Operational usage
Under the parent folder, we provide a SLURM job submission example to run the workflow for a single time step.

| File                   | Description                                       |
|------------------------|---------------------------------------------------|
| 'job_submit.sh`        | Bash script to submit a SLURM job to run the workflow |
| 'SNIP_LDT_template.sh` | Template bash script for workflow execution                           |
 
For example, to run model at 06:00 UTC on 2025-01-20, execute the following command in the terminal: 
```shell
./job_submit.sh 202501200600
```
## Support
- Technical Contact: kehan.yang@nasa.gov

## Reference
- Yoon, Y., Kemp, E.M., Kumar, S.V., Wegiel, J.W., Vuyovich, C.M., Peters-Lidard, C., 2022. Development of a global operational snow analysis: The US Air Force Snow and Ice Analysis. Remote Sensing of Environment 278, 113080. https://doi.org/10.1016/j.rse.2022.113080
- Foster, J.L., Sun, C., Walker, J.P., Kelly, R., Chang, A., Dong, J., Powell, H., 2005. Quantifying the uncertainty in passive microwave snow water equivalent observations. Remote Sensing of Environment 94, 187–203. https://doi.org/10.1016/j.rse.2004.09.012
- KELLY, R., 2009. AMSR-E積雪深アルゴリズム : 解説と初期成果〔英文〕. 日本リモートセンシング学会誌 29, 307–317. https://doi.org/10.11440/rssj.29.307
