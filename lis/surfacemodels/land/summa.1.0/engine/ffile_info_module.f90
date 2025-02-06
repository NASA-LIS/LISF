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

module ffile_info_module
USE nrtype
USE netcdf
USE globalData,only:integerMissing
implicit none
private
public::ffile_info
contains


 ! ************************************************************************************************
 ! public subroutine ffile_info: read information on model forcing files
 ! ************************************************************************************************
 subroutine ffile_info(nGRU,err,message)
 ! used to read metadata on the forcing data file
 USE ascii_util_module,only:file_open
 USE netcdf_util_module,only:nc_file_open    ! open netCDF file
 USE netcdf_util_module,only:netcdf_err      ! netcdf error handling function
 USE summaFileManager,only:SETNGS_PATH       ! path for metadata files
 USE summaFileManager,only:INPUT_PATH        ! path for forcing files
 USE summaFileManager,only:FORCING_FILELIST  ! list of model forcing files
 USE globalData,only:forcFileInfo,data_step  ! info on model forcing file
 USE globalData,only:forc_meta               ! forcing metadata
 USE get_ixname_module,only:get_ixtime,get_ixforce  ! identify index of named variable
 USE ascii_util_module,only:get_vlines       ! get a vector of non-comment lines
 USE ascii_util_module,only:split_line       ! split a line into words
 USE globalData,only:gru_struc               ! gru-hru mapping structure
 implicit none
 ! define input & output
 integer(i4b),intent(in)              :: nGRU           ! number of grouped response units
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 ! netcdf file i/o related
 integer(i4b)                         :: ncid           ! netcdf file id
 integer(i4b)                         :: mode           ! netCDF file open mode
 integer(i4b)                         :: varid          ! netcdf variable id
 integer(i4b)                         :: dimId          ! netcdf dimension id
 character(LEN=nf90_max_name)         :: varName        ! character array of netcdf variable name
 integer(i4b)                         :: iNC            ! index of a variable in netcdf file
 integer(i4b)                         :: nvar           ! number of variables in netcdf local attribute file
 ! the rest
 character(LEN=1024),allocatable      :: dataLines(:)   ! vector of lines of information (non-comment lines)
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=256)                   :: infile         ! input filename
 integer(i4b)                         :: unt            ! file unit (free unit output from file_open)
 character(LEN=256)                   :: filenameData   ! name of forcing datafile
 integer(i4b)                         :: ivar           ! index of model variable
 integer(i4b)                         :: iFile          ! counter for forcing files
 integer(i4b)                         :: nFile          ! number of forcing files in forcing file list
 integer(i4b)                         :: file_nHRU      ! number of HRUs in current forcing file
 integer(i4b)                         :: nForcing       ! number of forcing variables
 integer(i4b)                         :: iGRU,localHRU  ! index of GRU and HRU
 integer(i4b)                         :: ncHruId(1)     ! hruID from the forcing files
 real(dp)                             :: dataStep_iFile ! data step for a given forcing data file
 logical(lgt)                         :: xist           ! .TRUE. if the file exists

 ! Start procedure here
 err=0; message="ffile_info/"
 ! ------------------------------------------------------------------------------------------------------------------
 ! (1) read from the list of forcing files
 ! ------------------------------------------------------------------------------------------------------------------
 ! build filename for forcing file list
 infile = trim(SETNGS_PATH)//trim(FORCING_FILELIST)

 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! get a list of character strings from non-comment lines
 call get_vlines(unt,dataLines,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
 nFile = size(dataLines)

 ! allocate space for forcing information
 if(allocated(forcFileInfo)) deallocate(forcFileInfo)
 allocate(forcFileInfo(nFile), stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for forcFileInfo'; return; end if

 ! poputate the forcingInfo structure with filenames 
 do iFile=1,nFile
  ! split the line into "words" (expect one word: the file describing forcing data for that index)
  read(dataLines(iFile),*,iostat=err) filenameData
  if(err/=0)then; message=trim(message)//'problem reading a line of data from file ['//trim(infile)//']'; return; end if
  ! set forcing file name attribute
  forcFileInfo(iFile)%filenmData = trim(filenameData)
 end do  ! (looping through files)

 ! close ascii file
 close(unit=unt,iostat=err); if(err/=0)then;message=trim(message)//'problem closing forcing file list'; return; end if

 ! ------------------------------------------------------------------------------------------------------------------
 ! (2) pull descriptive information from netcdf forcing file and check number of HRUs in each forcing file matches nHRU
 ! ------------------------------------------------------------------------------------------------------------------

 ! get the number of forcing variables
 nForcing = size(forc_meta)

 ! loop through files, and read descriptive information from each file
 do iFile=1,nFile

  ! ensure allocatable structure components are deallocated
  if(allocated(forcFileInfo(iFile)%data_id)) deallocate(forcFileInfo(iFile)%data_id)
  if(allocated(forcFileInfo(iFile)%varName)) deallocate(forcFileInfo(iFile)%varName)

  ! allocate space for structure components
  allocate(forcFileInfo(iFile)%data_id(nForcing), forcFileInfo(iFile)%varName(nForcing), stat=err)
  if(err/=0)then; err=41; message=trim(message)//"problemAllocateStructureElement"; return; end if

  ! initialize variable ids to missing
  forcFileInfo(iFile)%data_id(:) = integerMissing

  ! build filename for actual forcing file
  infile = trim(INPUT_PATH)//trim(forcFileInfo(iFile)%filenmData)
  ! check if file exists
  inquire(file=trim(infile),exist=xist)
  if(.not.xist)then
   message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
   err=10; return
  end if
  
  ! open file
  mode=nf90_NoWrite
  call nc_file_open(trim(infile), mode, ncid, err, cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! how many variables are there?
  err = nf90_inquire(ncid, nvariables=nVar)
  call netcdf_err(err,message); if (err/=0) return

  ! set nVar attribute
  forcFileInfo(iFile)%nVars = nVar

  ! inquire nhru dimension size
  err = nf90_inq_dimid(ncid,'hru',dimId);                 if(err/=0)then; message=trim(message)//'cannot find dimension hru'; return; endif
  err = nf90_inquire_dimension(ncid,dimId,len=file_nHRU); if(err/=0)then; message=trim(message)//'cannot read dimension hru'; return; endif 

  ! inquire time dimension size
  err = nf90_inq_dimid(ncid,'time',dimId);                                     if(err/=0)then; message=trim(message)//'cannot find dimension time'; return; end if
  err = nf90_inquire_dimension(ncid,dimId,len=forcFileInfo(iFile)%nTimeSteps); if(err/=0)then; message=trim(message)//'cannot read dimension time'; return; end if

  ! loop through all variables in netcdf file, check to see if everything needed to run the model exists and data_step is correct
  do iNC=1,nVar

   ! inqure about current variable name, type, number of dimensions
   err = nf90_inquire_variable(ncid,iNC,name=varName)
   if(err/=0)then; message=trim(message)//'problem inquiring variable: '//trim(varName); return; end if

   ! process variable
   select case(trim(varName))

    ! if variable is in the forcing vector
    case('time','pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres','spechum')

     ! get variable index
     ivar = get_ixforce(trim(varname))
     if(ivar < 0)then;                               err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return; end if
     if(ivar>size(forcFileInfo(iFile)%data_id))then; err=35; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return; end if

     ! put netcdf file variable index in the forcing file metadata structure
     err = nf90_inq_varid(ncid, trim(varName), forcFileInfo(iFile)%data_id(ivar))
     if(err/=0)then; message=trim(message)//"problem inquiring forcing variable[var="//trim(varName)//"]"; return; end if

     ! put variable name in forcing file metadata structure
     forcFileInfo(iFile)%varName(ivar) = trim(varName)

     ! get first time from file, place into forcFileInfo
     if(trim(varname)=='time')then
      err = nf90_get_var(ncid,forcFileInfo(iFile)%data_id(ivar),forcFileInfo(iFile)%firstJulDay,start=(/1/))
      if(err/=0)then; message=trim(message)//'problem reading first Julian day'; return; end if
     end if  ! if the variable name is time

    ! data step
    case('data_step' )

     ! read data_step from netcdf file
     err = nf90_inq_varid(ncid, "data_step", varId); if(err/=0)then; message=trim(message)//'cannot find data_step'; return; end if
     err = nf90_get_var(ncid,varid,dataStep_iFile);  if(err/=0)then; message=trim(message)//'cannot read data_step'; return; end if

     ! check data_step is the same for all forcing files
     if(iFile == 1)then
      data_step = dataStep_iFile
     else
      if(abs(dataStep_iFile - data_step) > epsilon(dataStep_iFile))then
       write(message,'(a,i0,a)') trim(message)//'data step for forcing file ',iFile,'differs from the datastep of the first forcing file'
       err=20; return
      end if
     end if

    ! HRU id -- required
    case('hruId')

     ! check to see if hruId exists as a variable, this is a required variable
     err = nf90_inq_varid(ncid,trim(varname),varId)
     if(err/=0)then; message=trim(message)//'hruID variable not present'; return; endif

     ! check that the hruId is what we expect
     ! NOTE: we enforce that the HRU order in the forcing files is the same as in the zLocalAttributes files (too slow otherwise)
     do iGRU=1,nGRU
      do localHRU=1,gru_struc(iGRU)%hruCount
       err = nf90_get_var(ncid,varId,ncHruId,start=(/gru_struc(iGRU)%hruInfo(localHRU)%hru_nc/),count=(/1/))
       if(gru_struc(iGRU)%hruInfo(localHRU)%hru_id /= ncHruId(1))then
        write(message,'(a,i0,a,i0,a,i0,a,a)') trim(message)//'hruId for global HRU: ',gru_struc(iGRU)%hruInfo(localHRU)%hru_nc,' - ',  &
            ncHruId(1), ' differs from the expected: ',gru_struc(iGRU)%hruInfo(localHRU)%hru_id, ' in file ', trim(infile)
        write(message,'(a)') trim(message)//' order of hruId in forcing file needs to match order in zLocalAttributes.nc'
        err=40; return
       endif
      end do
     end do
  
    ! OK to have additional variables in the forcing file that are not used
    case default; cycle
   end select  ! select variable name
  end do ! (end of netcdf file variable loop)

  ! check to see if any forcing variables are missed
  if(any(forcFileInfo(iFile)%data_id(:)==integerMissing))then
   do iVar=1,size(forcFileInfo(iFile)%data_id)
    if(forcFileInfo(iFile)%data_id(iVar)==integerMissing)then; err=40; message=trim(message)//"variable missing [var='"//trim(forcFileInfo(iFile)%varname(iVar))//"']"; return; end if
   end do
  end if

 end do ! (loop through files)

 end subroutine ffile_info

end module ffile_info_module
