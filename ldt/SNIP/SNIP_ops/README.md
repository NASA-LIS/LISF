# USAF Snow & Ice Product (SNIP) sub-System

An AI/ML-enhanced snow and ice data production subsystem that processes AMSR2 L1R and WSF-M SDR data to generate next-generation snow and ice products for U.S. Air Force (USAF) operational use. This repository contains the Phase 0 implementation of the Snow and Ice Products (SNIP) subsystem. It produces snow depth estimates for the Phase 1: LDT.

## Overview

The SNIP system addresses the need for accurate, automated snow depth mapping from JAXA's GCOM-W/AMSR2 microwave radiometer data, and the WSF-M data. This system:
- Developed as part of NASA's next-generation Earth observation capabilities, leveraging brightness temperature data for enhanced snow depth estimation
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


The SNIP pipeline is controlled via a central JSON configuration file. Below is a breakdown of the available parameters:
#### Basic configuration 
* **`target_datetime`**: Target processing time in `YYYYMMDDHHMM` format (e.g., `"202501200600"`).
* **`project_path`**: Root directory of the project containing source code and configuration files.
* **`input_SD`**: Passive microwave data source (`"WSF"` or `"AMSR2"`). *Note: This value can be overridden at runtime using the `--input` flag in `workflow.py`.*
* **`mpirun`**: Command used to execute MPI jobs (typically `"mpirun"`).

#### Model & Output Paths
* **`model_path`**: Path to the pre-trained machine learning model directory/file.
* **`output_dir`**: Directory for saving the final ML-based snow depth retrievals and reprojected data products.
**`resampled_base`**: Directory for storing the resampled WSF SDR output. Save directory as it in the WSF read and resampling config file.

#### AMSR2 Settings
* **`AMSR2_source`**: AMSR2 L1R data source specification (`"NOAA"` or `"JAXA"`).
* **`amsr2_path`**: Primary storage directory for AMSR2 L1R input data.
* **`amsr2_merge_path`**: Output directory for merged AMSR2 data products (serves as the input for the ML prediction pipeline).
* **`AMSR2_template_path`**: Path to the reference template file used for reprojecting snow depth retrieved via traditional baseline approaches.

#### Ancillary Data & Templates
* **`template_path`**: Path to the reference template data used for spatial reprojection operations.
* **`viirs_path`**: Directory containing the VIIRS snow cover input files.
* **`fraction_forest_cover`**: Path to the fractional forest cover dataset (e.g., MCD12Q1).
* **`forest_density`**: Path to the forest density dataset (e.g., MOD44B).

#### Processing Flags
* **`flag_output_kelly`**: Enable snow depth calculation using the Kelly (2009) baseline approach (default: `false`).
* **`flag_output_foster`**: Enable snow depth calculation using the Foster et al. (2005) baseline approach, which is the method used by USAFSI (default: `false`).
* **`apply_viirs_mask`**: Enable the application of the VIIRS snow mask for AMSR2 processing (default: `false`).


## Data Processing 

| File | Description |
|------|-------------|
| `data_processing/amsr2_reader.py` | Core data processing module that:<br>• Reads AMSR2 L1R data files<br>• Reprojects data to common grids<br>• Processes data quality flags<br>• Merges multiple data extents<br>• Saves processed data to NetCDF format |
| `data_processing/WSF_reader.py` | Core data processing module to call LDT for WSF-M data reading and resampling. Data processing includes<br>• Search and reads WSF SDR data files<br>• Resample and reprojects data to AF grids<br>• Processes data quality flags<br>• Saves data to NetCDF format |

## Machine Learning Prediction 

| File                                        | Description                                                                                                                                                                                                                                                                                       |
|---------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `ml_prediction/run_prediction.py`           | ML inference pipeline that:<br>• Loads pretrained machine learning models<br>• Runs ML predictions on processed data<br>• Applies VIIRS Snow Covered Area (SCA) masking (default config false) <br>• Handles data quality flagging<br>• Reprojects results to Analysis Framework (AF) grids<br>• Saves prediction results |
| `ml_prediction/run_traditional_approach.py` | Traditional snow depth retrival pipeline that:<br>• Calculate snow depth  using Kelly, (2009) approach and Foster et al, (2005) approach. <br>• Options provided in config file to output these two snow depth results. <br>• Saves prediction results at PMW native resolution        |
| `ml_prediction/run_prediction_WSF.py` |   ML inference pipeline that:<br>• Loads pretrained machine learning models<br>• Runs ML predictions on hourly resampled WSF data<br>• Merges predicted snow depth every six hours for LDT post-processing<br>• Applies VIIRS Snow Covered Area (SCA) masking (default config false) <br>• Handles data quality flagging|


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
| `03_ldt_integration.config` | Template configuration file containing project settings and parameters |

## Usage
* Generate the LDT executable following the build instructions with the
  software. Be sure to select "little endian".
* Run the LDT executable using the `lis_input.nrt_streamflow.noah39.nc` file
  and the input datasets.
* Output netCDF files will be in the `../data/output/snip` directory with the log file saved in `../data/output/logs`.

## Operational SLURM Workflow

The operational SNIP subsystem is organized as a three-step SLURM workflow. This generates SNIP snow and ice analyses for the 557WW ~10-km global domain.

1. **Step 1: WSF-M data reading and resampling** — Uses LDT to read and resample incoming WSF brightness temperature inputs to the common snow analysis grid.
2. **Step 2: ML snow depth prediction** — Runs the AI/ML Python pipeline to generate snow depth estimates from the resampled PMW input data.
3. **Step 3: LDT data integration** — Performs LDT quality control, gap filling, and conversion of the snow depth output for integration with the LIS workflow. The final product includes snow depth, snow age, ice concentration, ice mask, and ice age in NetCDF format. 

### Job Submission Files

Under the parent folder, we provide SLURM job submission script samples to run the workflow for a single time step.

| File | Description |
|------|-------------|
| `01_wsf_resample.sh` | Step 1 SLURM script that runs WSF resampling via `SNIP_ops/data_processing/WSF_reader.py`
| `01_wsf_resample.ldt.config` | Config file for WSF resampling. Key parameters include:<br>• **`ldt_config_template`**: Path to the LDT configuration template file<br>• **`v522_sdr_base`**: Base directory for WSF-M v5.2.2 SDR input data<br>• **`raw_sdr_base`**: Base directory for WSF-M SDR input data for newer versions<br>• **`resampled_base`**: Directory for storing the resampled WSF SDR output |
| `02_snip_python.sh` | Step 2 SLURM script that runs the Python module for the ML-based snow depth prediction pipeline. |
| `./SNIP_ops/config/SNIP_config.json` | Configuration file for the Python prediction step with the details introduced in the previous table. Make sure the `project_path`is provided correctly.|
| `03_ldt_integration.sh` | Step 3 SLURM script that runs LDT and outputs the final SNIP product. |
| `03_ldt_integration.config` | LDT config file. Similar to the USAFSI config file, updated to use WSF-M and AMSR2 as inputs. |

### Command Line Execution Examples

**To run the full workflow for snow depth retrieval at 00:00 UTC on 2025-01-21 on HPC11:** 
```shell
sbatch 01_wsf_resample.sh 202501210000
sbatch 02_snip_python.sh 202501210000
sbatch 03_ldt_integration.sh 202501210000
```

## Support
- Technical Contact: kehan.yang@nasa.gov

## Reference
- Yoon, Y., Kemp, E.M., Kumar, S.V., Wegiel, J.W., Vuyovich, C.M., Peters-Lidard, C., 2022. Development of a global operational snow analysis: The US Air Force Snow and Ice Analysis. Remote Sensing of Environment 278, 113080. https://doi.org/10.1016/j.rse.2022.113080
- Foster, J.L., Sun, C., Walker, J.P., Kelly, R., Chang, A., Dong, J., Powell, H., 2005. Quantifying the uncertainty in passive microwave snow water equivalent observations. Remote Sensing of Environment 94, 187–203. https://doi.org/10.1016/j.rse.2004.09.012
- KELLY, R., 2009. AMSR-E積雪深アルゴリズム : 解説と初期成果〔英文〕. 日本リモートセンシング学会誌 29, 307–317. https://doi.org/10.11440/rssj.29.307
