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
- ✅ VIIRS Snow Cover Area (SCA) masking integration
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
- **viirs_path**: VIIRS snow mask data directory path
- **template_path**:  Reference template data path used for spatial reprojection operations
- **apply_viirs_mask**: Boolean flag to enable/disable VIIRS snow extent masking (default: True)
- **reproject_USAF**: Boolean flag to control saving of reprojected data products (default: True)
- **source**: AMSR2 L1R data source specification - accepts "NOAA" or "JAXA"

## Data Processing 

| File | Description |
|------|-------------|
| `data_processing/amsr2_reader.py` | Core data processing module that:<br>• Reads AMSR2 L1R data files<br>• Reprojects data to common grids<br>• Processes data quality flags<br>• Merges multiple data extents<br>• Saves processed data to NetCDF format |

## Machine Learning Prediction 

| File | Description |
|------|-------------|
| `ml_prediction/run_prediction.py` | ML inference pipeline that:<br>• Loads pretrained machine learning models<br>• Runs ML predictions on processed data<br>• Applies VIIRS Snow Covered Area (SCA) masking<br>• Handles data quality flagging<br>• Reprojects results to Analysis Framework (AF) grids<br>• Saves prediction results |


## Main Execution

| File | Description |
|------|-------------|
| `main.py` | Main execution script that orchestrates the entire SNIP_ops workflow |


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

This directory contains:
* This README file
* An "input" tarball file that should be un-tarred to put the required
  input parameter datasets on disk where LDT will be able to find it.
  To un-tar this file, please run:
      tar xvfz ldt_usafsi_use-case_input_v76.tar.gz
* An "output" tarball file that should be un-tarred to put the target
  output on disk.  Running this LDT use case should produce output that
  matches data contained in this tarball.  To un-tar this file, run:
      tar xvfz ldt_usafsi_use-case_output_v76.tar.gz

## Usage
* Generate the LDT executable following the build instructions with the
  software.  Be sure to select "little endian".
* Run the LDT executable using the ldt.config.foc.nrt.noah39.param.76 file
  and the input datasets.
* Output netCDF files will be in the current directory.

## Operational usage
Under the job folder, we provide a few slurm job submission examples to run the workflow for either a specific period or one time step.

| File                       | Description         |
|----------------------------|---------------------|
| 'run_SNIP_LDT_workflow.sh` | Template bash file  |

## Support
- Technical Contact: kehan.yang@nasa.gov