#!/usr/bin/env bash

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

### Step 1: disable write_to_log in ./src/util/logging_mod.F90 from JULES code.
### logging_mod produces a lot of screen output. Replace it with a LIS version. 
vim -N -u NONE -i NONE -e -s \
    -c 'silent /\c\vsubroutine\s+write_to_log/|normal O#if 0' \
    -c 'silent /\c\vend\s*subroutine\s+write_to_log/|normal o#endif' \
    -c 'silent /\c\vsubroutine\s+write_to_log/-1|normal O
! Disable JULES logging
subroutine write_to_log(log_level, proc_name, message)
   integer, intent(in) :: log_level
   character(len=*), intent(in) :: proc_name
   character(len=*), intent(in) :: message
   return
end subroutine write_to_log
' \
    -c 'wq' ./src/util/logging_mod.F90

### Step 2: disable 4 JULES initialization subroutine calls 
sed -i 's/CALL\ init_output(nml_dir)/!\ CALL\ init_output(nml_dir)/g' \
  ./src/initialisation/standalone/init.F90
sed -i 's/CALL\ init_params(nml_dir)/!\ CALL\ init_params(nml_dir)/g' \
  ./src/initialisation/standalone/init.F90
sed -i 's/CALL\ init_parms()/!\ CALL\ init_parms()/g' \
  ./src/initialisation/standalone/init.F90
sed -i 's/CALL\ init_ic(nml_dir)/!\ CALL\ init_ic(nml_dir)/g' \
  ./src/initialisation/standalone/init.F90
sed -i 's/CALL\ init_dump()/!\ CALL\ init_dump()/g' \
  ./src/initialisation/standalone/init.F90

### Step 3: disable MPI in JULES 
### and rename mpi module to lis_mpi_dummy
vim -N -u NONE -i NONE -e -s \
    -c '/\c\v^\s*subroutine\s+mpi_bcast/,/\c\v^\s*end\s*subroutine\s+mpi_bcast/delete' \
    -c '/\c\v^\s*subroutine\s+mpi_type_create_resized/,/\c\v^\s*end\s*subroutine\s+mpi_type_create_resized/delete' \
    -c '/\c\v^\s*subroutine\s+mpi_type_get_extent/,/\c\v^\s*end\s*subroutine\s+mpi_type_get_extent/delete' \
    -c 'silent $ | normal O
subroutine mpi_bcast1(buffer, count, datatype, root, comm, error)
  implicit none
  integer, intent(in) :: count
  real, intent(inout) :: buffer(count)
  integer, intent(in) :: datatype
  integer, intent(in) :: root
  integer, intent(in) :: comm
  integer, intent(out) :: error
  ! since there is only one task, we don'\''t need to do anything to broadcast a value
  error = 0
  return
end subroutine mpi_bcast1

subroutine mpi_bcast2(buffer, count, datatype, root, comm, error)
  implicit none
  integer, intent(in) :: count
  real, intent(inout) :: buffer
  integer, intent(in) :: datatype
  integer, intent(in) :: root
  integer, intent(in) :: comm
  integer, intent(out) :: error
! since there is only one task, we don'\''t need to do anything to broadcast a value
  error = 0
  return
end subroutine mpi_bcast2

subroutine mpi_bcast3(buffer, count, datatype, root, comm, error)
  implicit none
  integer, intent(in) :: count
  logical, intent(inout) :: buffer
  integer, intent(in) :: datatype
  integer, intent(in) :: root
  integer, intent(in) :: comm
  integer, intent(out) :: error
! since there is only one task, we don'\''t need to do anything to broadcast a value
  error = 0
  return
end subroutine mpi_bcast3

subroutine mpi_bcast4(buffer, count, datatype, root, comm, error)
  implicit none
  integer, intent(in) :: count
  integer, intent(in) :: buffer(count)
  integer, intent(in) :: datatype
  integer, intent(in) :: root
  integer, intent(in) :: comm
  integer, intent(out) :: error
! since there is only one task, we don'\''t need to do anything to broadcast a value
  error = 0
  return
end subroutine mpi_bcast4

subroutine mpi_bcast5(buffer, count, datatype, root, comm, error)
  implicit none
  integer, intent(in) :: count
  integer, intent(in) :: buffer
  integer, intent(in) :: datatype
  integer, intent(in) :: root
  integer, intent(in) :: comm
  integer, intent(out) :: error
! since there is only one task, we don'\''t need to do anything to broadcast a value
  error = 0
  return
end subroutine mpi_bcast5

subroutine mpi_type_create_resized1(old_type, lb, extent, new_type, error)
  implicit none
  integer, intent(in) :: old_type
  integer, intent(in) :: lb
  integer, intent(in) :: extent
  integer, intent(out) :: new_type
  integer, intent(out) :: error
! all information about mpi types is ignored by the dummy mpi routines, so
! it doesn'\''t matter what we return
  new_type = 1
  error    = 0
  return
end subroutine mpi_type_create_resized1

subroutine mpi_type_create_resized2(old_type, lb, extent, new_type, error)
  implicit none
  integer, intent(in) :: old_type
  integer(kind=mpi_address_kind) :: lb
  integer(kind=mpi_address_kind),intent(in) ::  extent
  integer, intent(out) :: new_type
  integer, intent(out) :: error
! all information about mpi types is ignored by the dummy mpi routines, so
! it doesn'\''t matter what we return
  new_type = 1
  error    = 0
  return
end subroutine mpi_type_create_resized2

subroutine mpi_type_get_extent(type, lb, extent, error)
  implicit none
  integer, intent(in) :: type
  integer(kind=mpi_address_kind), intent(out) :: lb
  integer(kind=mpi_address_kind), intent(out) :: extent
  integer, intent(out) :: error
! all information about mpi types is ignored by the dummy mpi routines, so
! it doesn'\''t matter what we return
  lb     = 0
  extent = 4
  error  = 0
  return
end subroutine mpi_type_get_extent
' \
    -c 'wq' ./utils/mpi_dummy/mpi_routines.F90

vim -N -u NONE -i NONE -e -s \
    -c '%s/\v(module\s+)(mpi)/\1lis_mpi_dummy/i' \
    -c '/\c\vend\s*module\s+lis_mpi_dummy/normal O
interface mpi_bcast
  procedure mpi_bcast1, mpi_bcast2, mpi_bcast3, mpi_bcast4, mpi_bcast5
end interface
interface mpi_type_create_resized
  procedure mpi_type_create_resized1, mpi_type_create_resized2
end interface

contains
' \
    -c '/\c\vend\s*module\s+lis_mpi_dummy/-1 read ./utils/mpi_dummy/mpi_routines.F90' \
    -c '%s/\V#if !defined(UM_JULES) && defined(MPI_DUMMY)/#if 1/' \
    -c 'wq' ./utils/mpi_dummy/mpi_mod.F90

rm ./utils/mpi_dummy/mpi_routines.F90

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
sed -i -e 's/USE\ mpi/USE lis_mpi_dummy/g' \
       -e '/lis_mpi_dummy,/ s/$/,mpi_comm_size/g' \
  ./src/initialisation/standalone/ancillaries/init_rivers_props.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/initialisation/standalone/initial_conditions/init_ic.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/io/output/internal_open_output_file.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/io/dump/read_dump.inc
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/util/logging_mod.F90
sed -i 's/USE\ mpi/USE lis_mpi_dummy/g' \
  ./src/initialisation/standalone/init_rivers.F90
vim -N -u NONE -i NONE -e -s \
    -c '/\c\v^\s*function\s+file_ascii_open' \
    -c '/\c\v^\s*implicit\s+none' \
    -c 'normal Ouse lis_mpi_dummy' \
    -c 'wq' ./src/io/file_handling/core/drivers/ascii/file_ascii_open.inc

### Step 4: delete src/params/shared/cable_maths_constants_mod.F90, which has a bug and not necessary
# rm ./src/params/shared/cable_maths_constants_mod.F90

### Step 5: delete jules.F90
rm ./src/control/standalone/jules.F90
