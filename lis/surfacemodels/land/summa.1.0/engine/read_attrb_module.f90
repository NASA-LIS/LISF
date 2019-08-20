! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module read_attrb_module
USE nrtype
implicit none
private
public::read_dimension
public::read_dimension_LIS
public::read_attrb
public::read_attrb_LIS
contains


 ! ************************************************************************************************
 ! public subroutine read_dimension: read HRU and GRU dimension information on local attributes
 ! ************************************************************************************************
 subroutine read_dimension_LIS(attrFile,fileGRU,fileHRU,nGRU,nHRU,err,message,startGRU,checkHRU)
 USE netcdf
 USE netcdf_util_module,only:nc_file_open                   ! open netcdf file
 USE netcdf_util_module,only:nc_file_close                  ! close netcdf file
 USE nr_utility_module ,only:arth
 ! provide access to global data
 USE globalData,only:gru_struc                              ! gru->hru mapping structure
 USE globalData,only:index_map                              ! hru->gru mapping structure
 implicit none 

 character(*),intent(in)              :: attrFile           ! name of attributed file
 integer(i4b),intent(out)             :: fileGRU            ! number of GRUs in the input file
 integer(i4b),intent(out)             :: fileHRU            ! number of HRUs in the input file
 integer(i4b),intent(inout)           :: nGRU               ! number of GRUs in the run domain
 integer(i4b),intent(inout)           :: nHRU               ! number of HRUs in the run domain
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message
 integer(i4b),intent(in),optional     :: startGRU           ! index of the starting GRU for parallelization run
 integer(i4b),intent(in),optional     :: checkHRU           ! index of the HRU for a single HRU run

 ! locals
 integer(i4b)                         :: sGRU               ! starting GRU
 integer(i4b)                         :: iHRU               ! HRU couinting index
 integer(i4b)                         :: iGRU               ! GRU loop index
 integer(i4b),allocatable             :: gru_id(:),hru_id(:)! read gru/hru IDs in from attributes file
 integer(i4b),allocatable             :: hru2gru_id(:)      ! read hru->gru mapping in from attributes file
 integer(i4b),allocatable             :: hru_ix(:)          ! hru index for search

 ! define variables for NetCDF file operation
 integer(i4b)                         :: ncID               ! NetCDF file ID
 integer(i4b)                         :: varID              ! NetCDF variable ID
 integer(i4b)                         :: gruDimId           ! variable id of GRU dimension from netcdf file
 integer(i4b)                         :: hruDimId           ! variable id of HRU dimension from netcdf file
 character(len=256)                   :: cmessage           ! error message for downwind routine

#if 0 
 ! Start procedure here
 err=0; message="read_dimension/" 

 ! check that we do not have conflicting flags
 if(present(startGRU).and.present(checkHRU))then; message=trim(message)//'startGRU and checkHRU both exist'; return; end if

 ! open nc file
 call nc_file_open(trim(attrFile),nf90_noWrite,ncID,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! *********************************************************************************************
 ! read and set GRU dimensions 
 ! **********************************************************************************************
 ! get gru dimension of whole file
 err = nf90_inq_dimid(ncID,"gru",gruDimId);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding gru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID, gruDimId, len = fileGRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading gru dimension/'//trim(nf90_strerror(err)); return; end if

 ! get hru dimension of whole file
 err = nf90_inq_dimid(ncID,"hru",hruDimId);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID, hruDimId, len = fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if 
#endif

 ! get runtime GRU dimensions
 if     (present(startGRU)) then
  if (nGRU < 1) then; err=20; message=trim(message)//'nGRU < 1 for a startGRU run'; return; end if
  sGRU = startGRU
 elseif (present(checkHRU)) then
  nGRU = 1
 else
  sGRU = 1
  nGRU = fileGRU
 endif

 ! check dimensions
 if ((present(startGRU)).and.(startGRU + nGRU > fileGRU)) then; err=20; message=trim(message)//'startGRU + nGRU is larger than then the GRU dimension'; return; end if
 if ((present(checkHRU)).and.(checkHRU        > fileHRU)) then; err=20; message=trim(message)//'checkHRU is larger than then the HRU dimension'       ; return; end if

 ! *********************************************************************************************
 ! read mapping vectors and populate mapping structures
 ! **********************************************************************************************
 ! allocate space for GRU indices
 allocate(gru_id(fileGRU))
 allocate(hru_ix(fileHRU),hru_id(fileHRU),hru2gru_id(fileHRU))

#if 0 
 ! read gru_id from netcdf file
 err = nf90_inq_varid(ncID,"gruId",varID);     if (err/=0) then; message=trim(message)//'problem finding gruId'; return; end if
 err = nf90_get_var(ncID,varID,gru_id);        if (err/=0) then; message=trim(message)//'problem reading gruId'; return; end if

 ! read hruIndex from netcdf file
 err = nf90_inq_varid(ncID,"hruId",varID);     if (err/=0) then; message=trim(message)//'problem finding hruId'; return; end if
 err = nf90_get_var(ncID,varID,hru_id);        if (err/=0) then; message=trim(message)//'problem reading hruId'; return; end if

 ! read hru2gru_id from netcdf file
 err = nf90_inq_varid(ncID,"hru2gruId",varID); if (err/=0) then; message=trim(message)//'problem finding hru2gruId'; return; end if
 err = nf90_get_var(ncID,varID,hru2gru_id);    if (err/=0) then; message=trim(message)//'problem reading hru2gruId'; return; end if

 ! close netcdf file
 call nc_file_close(ncID,err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
#endif

 ! array from 1 to total # of HRUs in attributes file
 hru_ix=arth(1,1,fileHRU)

! check that the mappings are not alreaday allocated
if (allocated(gru_struc)) then; message=trim(message)//'gru_struc is unexpectedly allocated'; return; end if
if (allocated(index_map)) then; message=trim(message)//'index_map is unexpectedly allocated'; return; end if

! allocate first level of gru to hru mapping
allocate(gru_struc(nGRU))  

! set gru to hru mapping  
if (present(checkHRU)) then                                                                    ! allocate space for single-HRU run

 ! gru to hru maping
 iGRU = 1
 gru_struc(iGRU)%hruCount             = 1                                                      ! number of HRUs in each GRU
 gru_struc(iGRU)%gruId                = hru2gru_id(checkHRU)                                   ! set gru id
 allocate(gru_struc(iGRU)%hruInfo(gru_struc(iGRU)%hruCount))                                   ! allocate second level of gru to hru map
 gru_struc(iGRU)%hruInfo(iGRU)%hru_nc = checkHRU                                               ! set hru id in attributes netcdf file
 gru_struc(iGRU)%hruInfo(iGRU)%hru_ix = 1                                                      ! set index of hru in run domain
 gru_struc(iGRU)%hruInfo(iGRU)%hru_id = hru_id(checkHRU)                                       ! set id of hru

else ! allocate space for anything except a single HRU run
 
 iHRU = 1
 do iGRU = 1,nGRU 

!  if (count(hru2gru_Id == gru_id(iGRU+sGRU-1)) < 1) then; err=20; message=trim(message)//'problem finding HRUs belonging to GRU'; return; end if
!  gru_struc(iGRU)%hruCount          = count(hru2gru_Id == gru_id(iGRU+sGRU-1))                 ! number of HRUs in each GRU
  gru_struc(iGRU)%hruCount          = 1
!  gru_struc(iGRU)%gruId             = gru_id(iGRU+sGRU-1)                                      ! set gru id
  gru_struc(iGRU)%gruId             = iGRU
  allocate(gru_struc(iGRU)%hruInfo(gru_struc(iGRU)%hruCount))                                  ! allocate second level of gru to hru map
!  gru_struc(iGRU)%hruInfo(:)%hru_nc = pack(hru_ix,hru2gru_id == gru_struc(iGRU)%gruId)         ! set hru id in attributes netcdf file
  gru_struc(iGRU)%hruInfo(:)%hru_nc = iGRU
!  gru_struc(iGRU)%hruInfo(:)%hru_ix = arth(iHRU,1,gru_struc(iGRU)%hruCount)                    ! set index of hru in run domain
  gru_struc(iGRU)%hruInfo(:)%hru_ix = iGRU
!  gru_struc(iGRU)%hruInfo(:)%hru_id = hru_id(gru_struc(iGRU)%hruInfo(:)%hru_nc)                ! set id of hru
  gru_struc(iGRU)%hruInfo(:)%hru_id = iGRU
  iHRU = iHRU + gru_struc(iGRU)%hruCount
 enddo ! iGRU = 1,nGRU

end if ! not checkHRU

! set hru to gru mapping 
nHRU = sum(gru_struc%hruCount)                                                                 ! total number of HRUs
allocate(index_map(nHRU))                                                                      ! allocate first level of hru to gru mapping

if (present(checkHRU)) then                                                                    ! allocate space for single-HRU run
 if (nHRU/=1) then; err=-20; message=trim(message)//'wrong # of HRUs for checkHRU run'; return; end if
 iGRU = 1; 
 index_map(1)%gru_ix   = iGRU                                                                  ! index of gru in run domain to which the hru belongs
 index_map(1)%localHRU = hru_ix(1)                                                             ! index of hru within the gru

else ! anything other than a single HRU run
 do iGRU = 1,nGRU 
  index_map(gru_struc(iGRU)%hruInfo(:)%hru_ix)%gru_ix   = iGRU                                 ! index of gru in run domain to which the hru belongs
  index_map(gru_struc(iGRU)%hruInfo(:)%hru_ix)%localHRU = hru_ix(1:gru_struc(iGRU)%hruCount)   ! index of hru within the gru 
 enddo ! iGRU = 1,nGRU

end if ! not checkHRU

end subroutine read_dimension_LIS


 ! ************************************************************************************************
 ! public subroutine read_dimension: read HRU and GRU dimension information on local attributes
 ! ************************************************************************************************
 subroutine read_dimension(attrFile,fileGRU,fileHRU,nGRU,nHRU,err,message,startGRU,checkHRU)
 USE netcdf
 USE netcdf_util_module,only:nc_file_open                   ! open netcdf file
 USE netcdf_util_module,only:nc_file_close                  ! close netcdf file
 USE nr_utility_module ,only:arth
 ! provide access to global data
 USE globalData,only:gru_struc                              ! gru->hru mapping structure
 USE globalData,only:index_map                              ! hru->gru mapping structure
 implicit none 

 character(*),intent(in)              :: attrFile           ! name of attributed file
 integer(i4b),intent(out)             :: fileGRU            ! number of GRUs in the input file
 integer(i4b),intent(out)             :: fileHRU            ! number of HRUs in the input file
 integer(i4b),intent(inout)           :: nGRU               ! number of GRUs in the run domain
 integer(i4b),intent(inout)           :: nHRU               ! number of HRUs in the run domain
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message
 integer(i4b),intent(in),optional     :: startGRU           ! index of the starting GRU for parallelization run
 integer(i4b),intent(in),optional     :: checkHRU           ! index of the HRU for a single HRU run

 ! locals
 integer(i4b)                         :: sGRU               ! starting GRU
 integer(i4b)                         :: iHRU               ! HRU couinting index
 integer(i4b)                         :: iGRU               ! GRU loop index
 integer(i4b),allocatable             :: gru_id(:),hru_id(:)! read gru/hru IDs in from attributes file
 integer(i4b),allocatable             :: hru2gru_id(:)      ! read hru->gru mapping in from attributes file
 integer(i4b),allocatable             :: hru_ix(:)          ! hru index for search

 ! define variables for NetCDF file operation
 integer(i4b)                         :: ncID               ! NetCDF file ID
 integer(i4b)                         :: varID              ! NetCDF variable ID
 integer(i4b)                         :: gruDimId           ! variable id of GRU dimension from netcdf file
 integer(i4b)                         :: hruDimId           ! variable id of HRU dimension from netcdf file
 character(len=256)                   :: cmessage           ! error message for downwind routine

 ! Start procedure here
 err=0; message="read_dimension/" 

 ! check that we do not have conflicting flags
 if(present(startGRU).and.present(checkHRU))then; message=trim(message)//'startGRU and checkHRU both exist'; return; end if

 ! open nc file
 call nc_file_open(trim(attrFile),nf90_noWrite,ncID,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! *********************************************************************************************
 ! read and set GRU dimensions 
 ! **********************************************************************************************
 ! get gru dimension of whole file
 err = nf90_inq_dimid(ncID,"gru",gruDimId);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding gru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID, gruDimId, len = fileGRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading gru dimension/'//trim(nf90_strerror(err)); return; end if

 ! get hru dimension of whole file
 err = nf90_inq_dimid(ncID,"hru",hruDimId);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID, hruDimId, len = fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if 

 ! get runtime GRU dimensions
 if     (present(startGRU)) then
  if (nGRU < 1) then; err=20; message=trim(message)//'nGRU < 1 for a startGRU run'; return; end if
  sGRU = startGRU
 elseif (present(checkHRU)) then
  nGRU = 1
 else
  sGRU = 1
  nGRU = fileGRU
 endif

 ! check dimensions
 if ((present(startGRU)).and.(startGRU + nGRU > fileGRU)) then; err=20; message=trim(message)//'startGRU + nGRU is larger than then the GRU dimension'; return; end if
 if ((present(checkHRU)).and.(checkHRU        > fileHRU)) then; err=20; message=trim(message)//'checkHRU is larger than then the HRU dimension'       ; return; end if

 ! *********************************************************************************************
 ! read mapping vectors and populate mapping structures
 ! **********************************************************************************************
 ! allocate space for GRU indices
 allocate(gru_id(fileGRU))
 allocate(hru_ix(fileHRU),hru_id(fileHRU),hru2gru_id(fileHRU))

 ! read gru_id from netcdf file
 err = nf90_inq_varid(ncID,"gruId",varID);     if (err/=0) then; message=trim(message)//'problem finding gruId'; return; end if
 err = nf90_get_var(ncID,varID,gru_id);        if (err/=0) then; message=trim(message)//'problem reading gruId'; return; end if

 ! read hruIndex from netcdf file
 err = nf90_inq_varid(ncID,"hruId",varID);     if (err/=0) then; message=trim(message)//'problem finding hruId'; return; end if
 err = nf90_get_var(ncID,varID,hru_id);        if (err/=0) then; message=trim(message)//'problem reading hruId'; return; end if

 ! read hru2gru_id from netcdf file
 err = nf90_inq_varid(ncID,"hru2gruId",varID); if (err/=0) then; message=trim(message)//'problem finding hru2gruId'; return; end if
 err = nf90_get_var(ncID,varID,hru2gru_id);    if (err/=0) then; message=trim(message)//'problem reading hru2gruId'; return; end if

 ! close netcdf file
 call nc_file_close(ncID,err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! array from 1 to total # of HRUs in attributes file
 hru_ix=arth(1,1,fileHRU)

! check that the mappings are not alreaday allocated
if (allocated(gru_struc)) then; message=trim(message)//'gru_struc is unexpectedly allocated'; return; end if
if (allocated(index_map)) then; message=trim(message)//'index_map is unexpectedly allocated'; return; end if

! allocate first level of gru to hru mapping
allocate(gru_struc(nGRU))  

! set gru to hru mapping  
if (present(checkHRU)) then                                                                    ! allocate space for single-HRU run

 ! gru to hru maping
 iGRU = 1
 gru_struc(iGRU)%hruCount             = 1                                                      ! number of HRUs in each GRU
 gru_struc(iGRU)%gruId                = hru2gru_id(checkHRU)                                   ! set gru id
 allocate(gru_struc(iGRU)%hruInfo(gru_struc(iGRU)%hruCount))                                   ! allocate second level of gru to hru map
 gru_struc(iGRU)%hruInfo(iGRU)%hru_nc = checkHRU                                               ! set hru id in attributes netcdf file
 gru_struc(iGRU)%hruInfo(iGRU)%hru_ix = 1                                                      ! set index of hru in run domain
 gru_struc(iGRU)%hruInfo(iGRU)%hru_id = hru_id(checkHRU)                                       ! set id of hru

else ! allocate space for anything except a single HRU run
 
 iHRU = 1
 do iGRU = 1,nGRU 

  if (count(hru2gru_Id == gru_id(iGRU+sGRU-1)) < 1) then; err=20; message=trim(message)//'problem finding HRUs belonging to GRU'; return; end if
  gru_struc(iGRU)%hruCount          = count(hru2gru_Id == gru_id(iGRU+sGRU-1))                 ! number of HRUs in each GRU
  gru_struc(iGRU)%gruId             = gru_id(iGRU+sGRU-1)                                      ! set gru id
  allocate(gru_struc(iGRU)%hruInfo(gru_struc(iGRU)%hruCount))                                  ! allocate second level of gru to hru map
  gru_struc(iGRU)%hruInfo(:)%hru_nc = pack(hru_ix,hru2gru_id == gru_struc(iGRU)%gruId)         ! set hru id in attributes netcdf file
  gru_struc(iGRU)%hruInfo(:)%hru_ix = arth(iHRU,1,gru_struc(iGRU)%hruCount)                    ! set index of hru in run domain
  gru_struc(iGRU)%hruInfo(:)%hru_id = hru_id(gru_struc(iGRU)%hruInfo(:)%hru_nc)                ! set id of hru
  iHRU = iHRU + gru_struc(iGRU)%hruCount
 enddo ! iGRU = 1,nGRU

end if ! not checkHRU

! set hru to gru mapping 
nHRU = sum(gru_struc%hruCount)                                                                 ! total number of HRUs
allocate(index_map(nHRU))                                                                      ! allocate first level of hru to gru mapping

if (present(checkHRU)) then                                                                    ! allocate space for single-HRU run
 if (nHRU/=1) then; err=-20; message=trim(message)//'wrong # of HRUs for checkHRU run'; return; end if
 iGRU = 1; 
 index_map(1)%gru_ix   = iGRU                                                                  ! index of gru in run domain to which the hru belongs
 index_map(1)%localHRU = hru_ix(1)                                                             ! index of hru within the gru

else ! anything other than a single HRU run
 do iGRU = 1,nGRU 
  index_map(gru_struc(iGRU)%hruInfo(:)%hru_ix)%gru_ix   = iGRU                                 ! index of gru in run domain to which the hru belongs
  index_map(gru_struc(iGRU)%hruInfo(:)%hru_ix)%localHRU = hru_ix(1:gru_struc(iGRU)%hruCount)   ! index of hru within the gru 
 enddo ! iGRU = 1,nGRU

end if ! not checkHRU

end subroutine read_dimension

 ! ************************************************************************************************
 ! public subroutine read_attrb: read information on local attributes
 ! ************************************************************************************************
 subroutine read_attrb(attrFile,nGRU,attrStruct,typeStruct,err,message)
 ! provide access to subroutines
 USE netcdf
 USE netcdf_util_module,only:nc_file_open                   ! open netcdf file
 USE netcdf_util_module,only:nc_file_close                  ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                     ! netcdf error handling function
 ! provide access to derived data types
 USE data_types,only:gru_hru_int                            ! x%gru(:)%hru(:)%var(:)     (i4b)
 USE data_types,only:gru_hru_double                         ! x%gru(:)%hru(:)%var(:)     (dp)
 ! provide access to global data
 USE globalData,only:gru_struc                              ! gru-hru mapping structure
 USE globalData,only:attr_meta,type_meta                    ! metadata structures
 USE get_ixname_module,only:get_ixAttr,get_ixType           ! access function to find index of elements in structure
 implicit none

 ! io vars
 character(*)                         :: attrFile           ! input filename
 integer(i4b),intent(in)              :: nGRU               ! number of grouped response units
 type(gru_hru_double),intent(inout)   :: attrStruct         ! local attributes for each HRU
 type(gru_hru_int),intent(inout)      :: typeStruct         ! local classification of soil veg etc. for each HRU
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message

 ! define local variables
 character(len=256)                   :: cmessage           ! error message for downwind routine
 integer(i4b)                         :: iVar               ! loop through varibles in the netcdf file
 integer(i4b)                         :: iHRU               ! index of an HRU within a GRU
 integer(i4b)                         :: iGRU               ! index of an GRU
 integer(i4b)                         :: varType            ! type of variable (categorical or numerical)
 integer(i4b)                         :: varIndx            ! index of variable within its data structure

 ! check structures
 integer(i4b)                         :: iCheck             ! index of an attribute name
 logical(lgt),allocatable             :: checkType(:)       ! vector to check if we have all desired categorical values
 logical(lgt),allocatable             :: checkAttr(:)       ! vector to check if we have all desired local attributes
 
 ! netcdf variables
 integer(i4b)                         :: ncID               ! netcdf file id
 character(LEN=nf90_max_name)         :: varName            ! character array of netcdf variable name
 integer(i4b)                         :: nVar               ! number of variables in netcdf local attribute file
 integer(i4b),parameter               :: categorical=101    ! named variable to denote categorical data
 integer(i4b),parameter               :: numerical=102      ! named variable to denote numerical data
 integer(i4b)                         :: categorical_var(1) ! temporary categorical variable from local attributes netcdf file
 real(dp)                             :: numeric_var(1)     ! temporary numeric variable from local attributes netcdf file

 ! define mapping variables

 ! Start procedure here
 err=0; message="read_attrb/"

 ! **********************************************************************************************
 ! (1) prepare check vectors 
 ! **********************************************************************************************
 allocate(checkType(size(type_meta)),checkAttr(size(attr_meta)),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for variable check vectors'; return; endif
 checkType(:) = .false.
 checkAttr(:) = .false.

 ! **********************************************************************************************
 ! (2) open netcdf file 
 ! **********************************************************************************************
 ! open file
 call nc_file_open(trim(attrFile),nf90_noWrite,ncID,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get number of variables total in netcdf file
 err = nf90_inquire(ncID,nvariables=nVar)
 call netcdf_err(err,message); if (err/=0) return

 ! **********************************************************************************************
 ! (3) read local attributes
 ! **********************************************************************************************
 ! loop through variables in netcdf file and pull out local attributes
 iCheck = 1
 do iVar = 1,nVar

  ! inqure about current variable name, type, number of dimensions
  err = nf90_inquire_variable(ncID,iVar,name=varName)
  if(err/=nf90_noerr)then; message=trim(message)//'problem inquiring variable: '//trim(varName)//'/'//trim(nf90_strerror(err)); return; endif

  ! find attribute name
  select case(trim(varName))

   ! ** categorical data
   case('hruId','vegTypeIndex','soilTypeIndex','slopeTypeIndex','downHRUindex')

    ! get the index of the variable
    varType = categorical
    if(trim(varName) == "hruId") then
     varIndx = get_ixType("hruIndex")
    else
     varIndx = get_ixType(varName)
    end if
    checkType(varIndx) = .true.

    ! check that the variable could be identified in the data structure
    if(varIndx < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(varName)//'] in data structure'; return; endif

    ! get data from netcdf file and store in vector
    do iGRU=1,nGRU
     do iHRU = 1,gru_struc(iGRU)%hruCount
      err = nf90_get_var(ncID,iVar,categorical_var,start=(/gru_struc(iGRU)%hruInfo(iHRU)%hru_nc/),count=(/1/)) 
      if(err/=nf90_noerr)then; message=trim(message)//'problem reading: '//trim(varName); return; end if
      typeStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = categorical_var(1)
     end do
    end do

   ! ** numerical data
   case('latitude','longitude','elevation','tan_slope','contourLength','HRUarea','mHeight')

    ! get the index of the variable
    varType = numerical
    varIndx = get_ixAttr(varName)
    checkAttr(varIndx) = .true.

    ! check that the variable could be identified in the data structure
    if(varIndx < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(varName)//'] in data structure'; return; endif

    ! get data from netcdf file and store in vector
    do iGRU=1,nGRU
     do iHRU = 1, gru_struc(iGRU)%hruCount      
      err = nf90_get_var(ncID,iVar,numeric_var,start=(/gru_struc(iGRU)%hruInfo(iHRU)%hru_nc/),count=(/1/))
      if(err/=nf90_noerr)then; message=trim(message)//'problem reading: '//trim(varName); return; end if
      attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = numeric_var(1)
     end do
    end do

   ! for mapping varibles, do nothing (information read above)
   case('hru2gruId','gruId'); cycle

   ! check that variables are what we expect
   case default; message=trim(message)//'unknown variable ['//trim(varName)//'] in local attributes file'; err=20; return

  end select ! select variable

 end do ! (looping through netcdf local attribute file)

 ! **********************************************************************************************
 ! (4) check that we have all the desired varaibles 
 ! **********************************************************************************************
 ! check that we have all desired categorical variables
 if(any(.not.checkType))then
  do iCheck = 1,size(type_meta)
   if(.not.checkType(iCheck))then; err=20; message=trim(message)//'missing variable ['//trim(type_meta(iCheck)%varname)//'] in local attributes file'; return; endif
  end do
 endif

 ! check that we have all desired local attributes
 if(any(.not.checkAttr))then
  do iCheck = 1,size(attr_meta)
   if(.not.checkAttr(iCheck))then; err=20; message=trim(message)//'missing variable ['//trim(attr_meta(iCheck)%varname)//'] in local attributes file'; return; endif
  end do
 endif

 ! **********************************************************************************************
 ! (5) close netcdf file 
 ! **********************************************************************************************
 call nc_file_close(ncID,err,cmessage)
 if (err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! free memory
 deallocate(checkType)
 deallocate(checkAttr)

 end subroutine read_attrb

 ! ************************************************************************************************
 ! public subroutine read_attrb: read information on local attributes
 ! ************************************************************************************************
 subroutine read_attrb_LIS(n, attrFile,nGRU,attrStruct,typeStruct,err,message)
 ! provide access to subroutines
 USE LIS_coreMod
 USE LIS_fileIOMod
 USE netcdf
 USE netcdf_util_module,only:nc_file_open                   ! open netcdf file
 USE netcdf_util_module,only:nc_file_close                  ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                     ! netcdf error handling function
 ! provide access to derived data types
 USE data_types,only:gru_hru_int                            ! x%gru(:)%hru(:)%var(:)     (i4b)
 USE data_types,only:gru_hru_double                         ! x%gru(:)%hru(:)%var(:)     (dp)
 ! provide access to global data
 USE globalData,only:gru_struc                              ! gru-hru mapping structure
 USE globalData,only:attr_meta,type_meta                    ! metadata structures
 USE get_ixname_module,only:get_ixAttr,get_ixType           ! access function to find index of elements in structure
 implicit none

 ! io vars
 integer                              :: n
 character(*)                         :: attrFile           ! input filename
 integer(i4b),intent(in)              :: nGRU               ! number of grouped response units
 type(gru_hru_double),intent(inout)   :: attrStruct         ! local attributes for each HRU
 type(gru_hru_int),intent(inout)      :: typeStruct         ! local classification of soil veg etc. for each HRU
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message

 ! define local variables
 character(len=256)                   :: cmessage           ! error message for downwind routine
 integer(i4b)                         :: iVar               ! loop through varibles in the netcdf file
 integer(i4b)                         :: iHRU               ! index of an HRU within a GRU
 integer(i4b)                         :: iGRU               ! index of an GRU
 integer(i4b)                         :: varType            ! type of variable (categorical or numerical)
 integer(i4b)                         :: varIndx            ! index of variable within its data structure

 ! check structures
 integer(i4b)                         :: iCheck             ! index of an attribute name
 logical(lgt),allocatable             :: checkType(:)       ! vector to check if we have all desired categorical values
 logical(lgt),allocatable             :: checkAttr(:)       ! vector to check if we have all desired local attributes
 
 ! netcdf variables
 integer(i4b)                         :: ncID               ! netcdf file id
 character(LEN=nf90_max_name)         :: varName            ! character array of netcdf variable name
 integer(i4b)                         :: nVar               ! number of variables in netcdf local attribute file
 integer(i4b),parameter               :: categorical=101    ! named variable to denote categorical data
 integer(i4b),parameter               :: numerical=102      ! named variable to denote numerical data
 integer(i4b)                         :: categorical_var(1) ! temporary categorical variable from local attributes netcdf file
 real(dp)                             :: numeric_var(1)     ! temporary numeric variable from local attributes netcdf file
 integer(i4b)                         :: col,row,t
 real(dp)                             :: var2d(LIS_rc%lnc(n),LIS_rc%lnr(n))
 real(dp)                             :: var1d(LIS_rc%ntiles(n))
 ! define mapping variables

 ! Start procedure here
 err=0; message="read_attrb/"

 varIndx = get_ixType("hruIndex")
 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1,gru_struc(iGRU)%hruCount
       typeStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = iGRU
    end do
 end do

 varIndx = get_ixType("vegTypeIndex")
 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1,gru_struc(iGRU)%hruCount
       typeStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = &
            LIS_domain(n)%tile(iGRU)%vegt
    end do
 end do

 varIndx = get_ixType("soilTypeIndex")
 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1,gru_struc(iGRU)%hruCount
       typeStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = &
            LIS_domain(n)%tile(iGRU)%soilt
    end do
 end do

 varIndx = get_ixType("slopeTypeIndex")
 call LIS_read_param(n,"SLOPETYPE",var2d)

 do t=1,LIS_rc%ntiles(n)
    col = LIS_domain(n)%tile(t)%col
    row = LIS_domain(n)%tile(t)%row
    if(LIS_domain(n)%gindex(col,row).ne.-1) then 
       var1d(t) = var2d(col,row)
    endif
 enddo

 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1,gru_struc(iGRU)%hruCount
       typeStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = &
            var1d(iGRU)
    end do
 end do

 varIndx = get_ixType("downHRUindex")

 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1,gru_struc(iGRU)%hruCount
       typeStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = 0            
    end do
 end do

 varIndx = get_ixAttr('latitude')
 call LIS_read_param(n,"lat",var2d)

 do t=1,LIS_rc%ntiles(n)
    col = LIS_domain(n)%tile(t)%col
    row = LIS_domain(n)%tile(t)%row
    if(LIS_domain(n)%gindex(col,row).ne.-1) then 
       var1d(t) = var2d(col,row)
    endif
 enddo
 
 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1, gru_struc(iGRU)%hruCount      
       attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) =&
            var1d(iGRU)
    end do
 end do

 varIndx = get_ixAttr('longitude')
 call LIS_read_param(n,"lon",var2d)

 do t=1,LIS_rc%ntiles(n)
    col = LIS_domain(n)%tile(t)%col
    row = LIS_domain(n)%tile(t)%row
    if(LIS_domain(n)%gindex(col,row).ne.-1) then 
       var1d(t) = var2d(col,row)
    endif
 enddo
 
 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1, gru_struc(iGRU)%hruCount      
       attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) =&
            var1d(iGRU)
    end do
 end do

 varIndx = get_ixAttr('elevation')
 call LIS_read_param(n,"ELEVATION",var2d)

 do t=1,LIS_rc%ntiles(n)
    col = LIS_domain(n)%tile(t)%col
    row = LIS_domain(n)%tile(t)%row
    if(LIS_domain(n)%gindex(col,row).ne.-1) then 
       var1d(t) = var2d(col,row)
    endif
 enddo
 
 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1, gru_struc(iGRU)%hruCount      
       attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) =&
            var1d(iGRU)
    end do
 end do

 varIndx = get_ixAttr('tan_slope')
 
 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1, gru_struc(iGRU)%hruCount      
       attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) =0.1
    end do
 end do

 varIndx = get_ixAttr('contourLength')
 call LIS_read_param(n,"HYMAP_river_length",var2d)

 do t=1,LIS_rc%ntiles(n)
    col = LIS_domain(n)%tile(t)%col
    row = LIS_domain(n)%tile(t)%row
    if(LIS_domain(n)%gindex(col,row).ne.-1) then 
       var1d(t) = var2d(col,row)
    endif
 enddo
 
 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1, gru_struc(iGRU)%hruCount      
       attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) =&
            var1d(iGRU)
    end do
 end do

 varIndx = get_ixAttr('HRUarea')
 call LIS_read_param(n,"HYMAP_grid_area",var2d)

 do t=1,LIS_rc%ntiles(n)
    col = LIS_domain(n)%tile(t)%col
    row = LIS_domain(n)%tile(t)%row
    if(LIS_domain(n)%gindex(col,row).ne.-1) then 
       var1d(t) = var2d(col,row)
    endif
 enddo
 
 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1, gru_struc(iGRU)%hruCount      
       attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) =&
            var1d(iGRU)
    end do
 end do

 varIndx = get_ixAttr('mHeight')
 
 ! get data from netcdf file and store in vector
 do iGRU=1,nGRU
    do iHRU = 1, gru_struc(iGRU)%hruCount      
       attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = 10
    end do
 end do

 end subroutine read_attrb_LIS

end module read_attrb_module
