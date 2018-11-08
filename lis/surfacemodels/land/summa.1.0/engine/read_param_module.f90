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

module read_param_module

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! runtime options
USE globalData,only:iRunModeFull,iRunModeGRU,iRunModeHRU ! run modes

! common modules
USE nrtype
USE netcdf
USE netcdf_util_module,only:nc_file_close  ! close netcdf file
USE netcdf_util_module,only:nc_file_open   ! open netcdf file
USE netcdf_util_module,only:netcdf_err     ! netcdf error handling function

! data types
USE data_types,only:gru_double             ! spatial double data type:  x%gru(:)%var(:)
USE data_types,only:gru_hru_int            ! spatial integer data type: x%gru(:)%hru(:)%var(:)
USE data_types,only:gru_hru_doubleVec      ! spatial double data type:  x%gru(:)%hru(:)%var(:)%dat(:)

implicit none
private
public::read_param
contains


 ! ************************************************************************************************
 ! public subroutine read_param: read trial model parameter values
 ! ************************************************************************************************
 subroutine read_param(iRunMode,checkHRU,startGRU,nHRU,nGRU,typeStruct,mparStruct,bparStruct,err,message)
 ! used to read model initial conditions
 USE summaFileManager,only:SETNGS_PATH               ! path for metadata files
 USE summaFileManager,only:PARAMETER_TRIAL           ! file with parameter trial values
 USE get_ixname_module,only:get_ixparam,get_ixbpar   ! access function to find index of elements in structure
 USE globalData,only:index_map,gru_struc             ! mapping from global HRUs to the elements in the data structures
 USE var_lookup,only:iLookPARAM,iLookTYPE            ! named variables to index elements of the data vectors
 implicit none
 ! define input
 integer(i4b),        intent(in)       :: iRunMode         ! run mode
 integer(i4b),        intent(in)       :: checkHRU         ! index of single HRU if runMode = checkHRU
 integer(i4b),        intent(in)       :: startGRU         ! index of single GRU if runMode = startGRU
 integer(i4b),        intent(in)       :: nHRU             ! number of global HRUs
 integer(i4b),        intent(in)       :: nGRU             ! number of global GRUs
 type(gru_hru_int),   intent(in)       :: typeStruct       ! local classification of soil veg etc. for each HRU
 ! define output
 type(gru_hru_doubleVec),intent(inout) :: mparStruct       ! model parameters
 type(gru_double)    ,intent(inout)    :: bparStruct       ! basin parameters
 integer(i4b),        intent(out)      :: err              ! error code
 character(*),        intent(out)      :: message          ! error message
 ! define local variables
 character(len=1024)                   :: cmessage         ! error message for downwind routine
 character(LEN=1024)                   :: infile           ! input filename
 integer(i4b)                          :: iHRU             ! index of HRU within data vector
 integer(i4b)                          :: localHRU,iGRU    ! index of HRU and GRU within data structure
 integer(i4b)                          :: ixParam          ! index of the model parameter in the data structure
 ! indices/metadata in the NetCDF file
 integer(i4b)                          :: ncid             ! netcdf id
 integer(i4b)                          :: nDims            ! number of dimensions
 integer(i4b)                          :: nVars            ! number of variables
 integer(i4b)                          :: idimid           ! dimension index
 integer(i4b)                          :: ivarid           ! variable index
 character(LEN=64)                     :: dimName          ! dimension name
 character(LEN=64)                     :: parName          ! parameter name
 integer(i4b)                          :: dimLength        ! dimension length
 integer(i4b)                          :: nHRU_file        ! number of HRUs in the parafile
 integer(i4b)                          :: nGRU_file        ! number of GRUs in the parafile
 integer(i4b)                          :: nSoil_file       ! number of soil layers in the file
 integer(i4b)                          :: idim_list(2)     ! list of dimension ids
 ! data in the netcdf file
 integer(i4b)                          :: parLength        ! length of the parameter data
 integer(i4b),allocatable              :: hruId(:)      ! HRU identifier in the file
 real(dp),allocatable                  :: parVector(:)     ! model parameter vector
 logical                               :: fexist           ! inquire whether the parmTrial file exists 
 integer(i4b)                          :: fHRU             ! index of HRU in input file

 ! Start procedure here
 err=0; message="read_param/"

 ! **********************************************************************************************
 ! * open files, etc.
 ! **********************************************************************************************

 ! build filename
 infile = trim(SETNGS_PATH)//trim(PARAMETER_TRIAL)

 ! do we need the file?
 inquire(file=trim(infile),exist=fexist)
 if (.not.fexist) return

 ! open file
 call nc_file_open(trim(infile),nf90_nowrite,ncid,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! get the number of variables in the parameter file
 err=nf90_inquire(ncid, nDimensions=nDims, nVariables=nVars)
 call netcdf_err(err,message); if (err/=0) then; err=20; return; end if

 ! initialize the number of HRUs
 nHRU_file=integerMissing
 nGRU_file=integerMissing

 ! get the length of the dimensions
 do idimid=1,nDims
  ! get the dimension name and length
  err=nf90_inquire_dimension(ncid, idimid, name=dimName, len=dimLength)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  ! get the number of HRUs
  if(trim(dimName)=='hru') nHRU_file=dimLength
  if(trim(dimName)=='gru') nGRU_file=dimLength
 end do

 ! allocate hruID vector
 allocate(hruId(nHRU_file))

 ! check HRU dimension exists
 if(nHRU_file==integerMissing)then
  message=trim(message)//'unable to identify HRU dimension in file '//trim(infile)
  err=20; return
 endif

 ! check have the correct number of HRUs
 if ((irunMode==irunModeFull).and.(nHRU_file/=nHRU)) then
  message=trim(message)//'incorrect number of HRUs in file '//trim(infile)
  err=20; return
 endif
 if ((irunMode==irunModeHRU).and.(nHRU_file<checkHRU)) then
  message=trim(message)//'not enough HRUs in file '//trim(infile)
  err=20; return
 endif
 
 ! check have the correct number of GRUs
 if ((irunMode==irunModeGRU).and.(nGRU_file<startGRU).and.(nGRU_file/=integerMissing)) then
  message=trim(message)//'not enough GRUs in file '//trim(infile)
  err=20; return
 endif
 if ((irunMode==irunModeFull).and.(nGRU_file/=nGRU).and.(nGRU_file/=integerMissing)) then
  message=trim(message)//'incorrect number of GRUs in file '//trim(infile)
  err=20; return
 endif
 
 ! loop through the parameters in the NetCDF file
 do ivarid=1,nVars

  ! get the parameter name
  err=nf90_inquire_variable(ncid, ivarid, name=parName)
  call netcdf_err(err,message); if (err/=0) then; err=20; return; end if

  ! **********************************************************************************************
  ! * read the HRU index
  ! **********************************************************************************************

  ! special case of the HRU id
  if(trim(parName)=='hruIndex' .or. trim(parName)=='hruId')then

   ! read HRUs
   err=nf90_get_var(ncid, ivarid, hruId)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

   ! check HRUs  -- expect HRUs to be in the same order as the local attributes
   if (iRunMode==iRunModeFull) then
    do iHRU=1,nHRU
     iGRU=index_map(iHRU)%gru_ix
     localHRU=index_map(iHRU)%localHRU
     if((hruId(iHRU)>0).and.(hruId(iHRU)/=typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%hruIndex)))then
      write(message,'(a,i0,a,i0,a)') trim(message)//'mismatch for HRU ', typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%hruIndex), '(param HRU = ', hruId(iHRU), ')'
      err=20; return
     endif
    end do  ! looping through HRUs

   else if (iRunMode==iRunModeGRU) then
    do iHRU=1,nHRU
     iGRU=index_map(iHRU)%gru_ix
     localHRU=index_map(iHRU)%localHRU
     fHRU = gru_struc(iGRU)%hruInfo(localHRU)%hru_nc
     if(hruId(fHRU)/=typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%hruIndex))then
     write(message,'(a,i0,a,i0,a)') trim(message)//'mismatch for HRU ', typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%hruIndex), '(param HRU = ', hruId(iHRU), ')'
     err=20; return
    endif
   enddo

   else if (iRunMode==iRunModeHRU) then
    iGRU=index_map(1)%gru_ix
    localHRU=index_map(1)%localHRU
    if(hruId(checkHRU)/=typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%hruIndex))then
     write(message,'(a,i0,a,i0,a)') trim(message)//'mismatch for HRU ', typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%hruIndex), '(param HRU = ', hruId(iHRU), ')'
     err=20; return
    endif

   else 
    err = 20; message = 'run mode not recognized'; return;
   end if

  ! all other variables
  else

   ! **********************************************************************************************
   ! * read the local parameters
   ! **********************************************************************************************

   ! get the local parameters
   ixParam = get_ixparam( trim(parName) )
   if(ixParam/=integerMissing)then

    ! get the variable shape
    err=nf90_inquire_variable(ncid, ivarid, nDims=nDims, dimids=idim_list)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
    ! get the length of the depth dimension (if it exists)
    if(nDims==2)then
  
     ! get the information on the 2nd dimension for 2-d variables
     err=nf90_inquire_dimension(ncid, idim_list(2), dimName, nSoil_file)
     if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     
     ! check that it is the depth dimension
     if(trim(dimName)/='depth')then
      message=trim(message)//'expect 2nd dimension of 2-d variable to be depth (dimension name = '//trim(dimName)//')'
      err=20; return
     endif
  
     ! check that the dimension length is correct
     if(size(mparStruct%gru(iGRU)%hru(localHRU)%var(ixParam)%dat) /= nSoil_file)then
      message=trim(message)//'unexpected number of soil layers in parameter file'
      err=20; return
     endif
 
     ! define parameter length
     parLength = nSoil_file
 
    else
     parLength = 1
    endif  ! if two dimensions
  
    ! allocate space for model parameters
    allocate(parVector(parLength),stat=err)
    if(err/=0)then
     message=trim(message)//'problem allocating space for parameter vector'
     err=20; return
    endif

    ! loop through HRUs
    do iHRU=1,nHRU
  
     ! map to the GRUs and HRUs    
     iGRU=index_map(iHRU)%gru_ix
     localHRU=index_map(iHRU)%localHRU
     fHRU = gru_struc(iGRU)%hruInfo(localHRU)%hru_nc

     ! read parameter data
     select case(nDims)
      case(1); err=nf90_get_var(ncid, ivarid, parVector, start=(/fHRU/), count=(/1/) )
      case(2); err=nf90_get_var(ncid, ivarid, parVector, start=(/fHRU,1/), count=(/1,nSoil_file/) )
      case default; err=20; message=trim(message)//'unexpected number of dimensions for parameter '//trim(parName)
     end select

     ! error check for the parameter read
     if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

     ! populate parameter structures
     select case(nDims)
      case(1); mparStruct%gru(iGRU)%hru(localHRU)%var(ixParam)%dat(:) = parVector(1)  ! also distributes scalar across depth dimension 
      case(2); mparStruct%gru(iGRU)%hru(localHRU)%var(ixParam)%dat(:) = parVector(:)
      case default; err=20; message=trim(message)//'unexpected number of dimensions for parameter '//trim(parName)
     end select

    end do  ! looping through HRUs

    ! deallocate space for model parameters
    deallocate(parVector,stat=err)
    if(err/=0)then
     message=trim(message)//'problem deallocating space for parameter vector'
     err=20; return
    endif

   ! **********************************************************************************************
   ! * read the basin parameters
   ! **********************************************************************************************

   ! get the basin parameters
   else

    ! get the parameter index
    ixParam = get_ixbpar( trim(parName) )

    ! check that we found it
    if(ixParam==integerMissing)then
     message=trim(message)//'parameter '//trim(parName)//' does not exist in the local or basin parameter structure'
     err=20; return
    endif

    ! allocate space for model parameters
    allocate(parVector(nGRU_file),stat=err)
    if(err/=0)then
     message=trim(message)//'problem allocating space for parameter vector'
     err=20; return
    endif

    ! read parameter data
    err=nf90_get_var(ncid, ivarid, parVector )
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

    ! populate parameter structures
    if (iRunMode==iRunModeGRU) then
     do iGRU=1,nGRU
      bparStruct%gru(iGRU)%var(ixParam) = parVector(iGRU+startGRU-1) 
     end do  ! looping through GRUs
    else if (iRunMode==iRunModeFull) then
     do iGRU=1,nGRU
      bparStruct%gru(iGRU)%var(ixParam) = parVector(iGRU) 
     end do  ! looping through GRUs
    else if (iRunMode==iRunModeHRU) then
     err = 20; message='checkHRU run mode not working'; return; 
    endif

    ! deallocate space for model parameters
    deallocate(parVector,stat=err)
    if(err/=0)then
     message=trim(message)//'problem deallocating space for parameter vector'
     err=20; return
    endif

   endif  ! reading the basin parameters

  endif  ! if a "regular" parameter (i.e., not the HRU index)

 end do ! (looping through the parameters in the NetCDF file)

 end subroutine read_param

end module read_param_module
