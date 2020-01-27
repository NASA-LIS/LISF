#!/usr/bin/env bash

### Step 1: delete ./src/util/logging_mod.F90 from JULES code
### logging_mod produces a lot of screen outputs. A LIS version 
### of logging_mod.F90 has been developed as part of the LIS-JULES
### container. Therefore, the original JULES version has to be removed. 
rm ./src/util/logging_mod.F90 

### Step 2: disable 4 JULES initialization subroutine calls 
sed -i 's/CALL\ init_output(nml_dir)/!\ CALL\ init_output(nml_dir)/g' \
  ./src/initialisation/standalone/init.F90
sed -i 's/CALL\ init_params(nml_dir)/!\ CALL\ init_params(nml_dir)/g' \
  ./src/initialisation/standalone/init.F90
sed -i 's/CALL\ init_parms()/!\ CALL\ init_parms()/g' \
  ./src/initialisation/standalone/init.F90
sed -i 's/CALL\ init_ic(nml_dir)/!\ CALL\ init_ic(nml_dir)/g' \
  ./src/initialisation/standalone/init.F90
sed -i 's/CALL\ init_dump(nml_dir)/!\ CALL\ init_dump(nml_dir)/g' \
  ./src/initialisation/standalone/init.F90

### Step 3: disable MPI in JULES 
### rename mpi mode into lis_mpi_dummy
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/control/standalone/parallel/decompose_domain.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/control/standalone/parallel/gather_land_field.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/control/standalone/parallel/is_master_task.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/control/standalone/parallel/scatter_land_field.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/control/standalone/spinup/spinup_check.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/initialisation/standalone/ancillaries/init_overbank.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/initialisation/standalone/ancillaries/init_rivers_props.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/initialisation/standalone/initial_conditions/init_ic.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/io/output/internal_open_output_file.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/io/dump/read_dump.inc

### Step 4: delete src/params/shared/cable_maths_constants_mod.F90, which has a bug and not necessary
rm ./src/params/shared/cable_maths_constants_mod.F90
### Step 5: delete jules.F9o0 
rm ./src/control/standalone/jules.F90
