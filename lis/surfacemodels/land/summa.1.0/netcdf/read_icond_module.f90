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

module read_icond_module
USE nrtype
USE netcdf
implicit none
private
public::read_icond
public::read_icond_nlayers
contains

 ! ************************************************************************************************
 ! public subroutine read_icond_nlayers: read model initial conditions file for number of snow/soil layers
 ! ************************************************************************************************
 subroutine read_icond_nlayers(iconFile,nGRU,indx_meta,err,message)
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookIndex                        ! variable lookup structure
 USE globalData,only:gru_struc                         ! gru-hru mapping structures
 USE netcdf_util_module,only:nc_file_close             ! close netcdf file
 USE netcdf_util_module,only:nc_file_open              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                ! netcdf error handling
 USE data_types,only:gru_hru_intVec                    ! actual data
 USE data_types,only:var_info                          ! metadata 
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)        ,intent(in)     :: iconFile       ! name of input (restart) file 
 integer(i4b)        ,intent(in)     :: nGRU           ! total # of GRUs in run domain 
 type(var_info)      ,intent(in)     :: indx_meta(:)   ! metadata 
 integer(i4b)        ,intent(out)    :: err            ! error code
 character(*)        ,intent(out)    :: message        ! returned error message

 ! locals
 integer(i4b)             :: ncID                       ! netcdf file id
 integer(i4b)             :: dimID                      ! netcdf file dimension id
 integer(i4b)             :: fileHRU                    ! number of HRUs in netcdf file
 integer(i4b)             :: snowID, soilID             ! netcdf variable ids
 integer(i4b)             :: iGRU, iHRU                 ! loop indexes
 integer(i4b),allocatable :: snowData(:)                ! number of snow layers in all HRUs
 integer(i4b),allocatable :: soilData(:)                ! number of soil layers in all HRUs
 character(len=256)       :: cmessage                   ! downstream error message 

 ! --------------------------------------------------------------------------------------------------------
 ! initialize error message
 err=0
 message = 'read_icond_nlayers/'

 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncid,err,cmessage);
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncID,"hru",dimId);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID,dimId,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! allocate sotrage for reading from file
 allocate(snowData(fileHRU))
 allocate(soilData(fileHRU))
 snowData = 0
 soilData = 0

 ! get variable ids
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookIndex%nSnow)%varName),snowid); call netcdf_err(err,message)
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookIndex%nSoil)%varName),soilid); call netcdf_err(err,message)

 ! get data
 err = nf90_get_var(ncid,snowid,snowData); call netcdf_err(err,message) 
 err = nf90_get_var(ncid,soilid,soilData); call netcdf_err(err,message) 

 ! assign to index structure - gru by hru
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   gru_struc(iGRU)%hruInfo(iHRU)%nSnow = snowData(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc)
   gru_struc(iGRU)%hruInfo(iHRU)%nSoil = soilData(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc)
  end do
 end do

 ! close file
 call nc_file_close(ncid,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(snowData,soilData)

 end subroutine read_icond_nlayers


 ! ************************************************************************************************
 ! public subroutine read_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine read_icond(iconFile,                      & ! name of initial conditions file
                       nGRU,                          & ! number of GRUs 
                       prog_meta,                     & ! metadata
                       progData,                      & ! model prognostic (state) variables
                       indxData,                      & ! layer index data
                       err,message)                     ! error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookVarType                       ! variable lookup structure
 USE var_lookup,only:iLookProg                          ! variable lookup structure
 USE var_lookup,only:iLookIndex                         ! variable lookup structure
 USE globalData,only:gru_struc                          ! gru-hru mapping structures
 USE globaldata,only:iname_soil,iname_snow              ! named variables to describe the type of layer
 USE netcdf_util_module,only:nc_file_close              ! close netcdf file
 USE netcdf_util_module,only:nc_file_open               ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                 ! netcdf error handling
 USE data_types,only:gru_hru_doubleVec                  ! actual data
 USE data_types,only:gru_hru_intVec                     ! actual data
 USE data_types,only:var_info                           ! metadata 
 USE get_ixName_module,only:get_varTypeName             ! to access type strings for error messages
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)           ,intent(in)     :: iconFile     ! name of netcdf file containing the initial conditions
 integer(i4b)           ,intent(in)     :: nGRU         ! number of grouped response units in simulation domain
 type(var_info)         ,intent(in)     :: prog_meta(:) ! prognostic metadata
 type(gru_hru_doubleVec),intent(inout)  :: progData     ! prognostic vars 
 type(gru_hru_intVec)   ,intent(inout)  :: indxData     ! layer indexes 
 integer(i4b)           ,intent(out)    :: err          ! error code
 character(*)           ,intent(out)    :: message      ! returned error message

 ! locals
 character(len=256)                     :: cmessage     ! downstream error message
 integer(i4b)                           :: fileHRU      ! number of HRUs in file
 integer(i4b)                           :: iVar         ! loop index 
 integer(i4b)                           :: iGRU         ! loop index 
 integer(i4b)                           :: iHRU         ! loop index 
 integer(i4b)                           :: dimID        ! varible dimension ids
 integer(i4b)                           :: ncVarID      ! variable ID in netcdf file
 character(256)                         :: dimName      ! not used except as a placeholder in call to inq_dim function
 integer(i4b)                           :: dimLen       ! data dimensions
 integer(i4b)                           :: ncID         ! netcdf file ID
 real(dp),allocatable                   :: varData(:,:) ! variable data storage        
 integer(i4b)                           :: nSoil, nSnow, nToto ! # layers
 integer(i4b),parameter                 :: nBand=2      ! number of spectral bands

 character(len=32),parameter            :: scalDimName   ='scalarv'  ! dimension name for scalar data
 character(len=32),parameter            :: midSoilDimName='midSoil'  ! dimension name for soil-only layers
 character(len=32),parameter            :: midTotoDimName='midToto'  ! dimension name for layered varaiables
 character(len=32),parameter            :: ifcTotoDimName='ifcToto'  ! dimension name for layered varaiables

 ! --------------------------------------------------------------------------------------------------------

 ! Start procedure here
 err=0; message="read_icond/"

 ! --------------------------------------------------------------------------------------------------------
 ! (1) read the file
 ! --------------------------------------------------------------------------------------------------------
 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncID,err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncID,"hru",dimID);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID,dimID,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! loop through prognostic variables
 do iVar = 1,size(prog_meta)

  ! skip variables that are computed later
  if(prog_meta(iVar)%varName=='scalarCanopyWat'           .or. &
     prog_meta(iVar)%varName=='spectralSnowAlbedoDiffuse' .or. &
     prog_meta(iVar)%varName=='scalarSurfaceTemp'         .or. &
     prog_meta(iVar)%varName=='mLayerVolFracWat'          .or. &
     prog_meta(iVar)%varName=='mLayerHeight'                   ) cycle

  ! get variable id
  err = nf90_inq_varid(ncID,trim(prog_meta(iVar)%varName),ncVarID); call netcdf_err(err,message)
  if(err/=0)then
   message=trim(message)//': problem with getting variable id, var='//trim(prog_meta(iVar)%varName)
   return
  endif

  ! get variable dimension IDs
  select case (prog_meta(iVar)%varType)
   case (iLookVarType%scalarv); err = nf90_inq_dimid(ncID,trim(scalDimName)   ,dimID); call netcdf_err(err,message)
   case (iLookVarType%midSoil); err = nf90_inq_dimid(ncID,trim(midSoilDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%midToto); err = nf90_inq_dimid(ncID,trim(midTotoDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%ifcToto); err = nf90_inq_dimid(ncID,trim(ifcTotoDimName),dimID); call netcdf_err(err,message)
   case default
    message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
    err=20; return
  end select
  
  ! check errors
  if(err/=0)then
   message=trim(message)//': problem with dimension ids, var='//trim(prog_meta(iVar)%varName)
   return
  endif

  ! get the dimension length
  err = nf90_inquire_dimension(ncID,dimID,dimName,dimLen); call netcdf_err(err,message)
  if(err/=0)then; message=trim(message)//': problem getting the dimension length'; return; endif

  ! iniitialize the variable data
  allocate(varData(fileHRU,dimLen),stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating variable data'; return; endif

  ! get data
  err = nf90_get_var(ncID,ncVarID,varData); call netcdf_err(err,message) 
  if(err/=0)then; message=trim(message)//': problem getting the data'; return; endif

  ! check data are not set to the fill value
  if( any( abs(varData - nf90_fill_double) < epsilon(varData) ) )then   
   message=trim(message)//"data set to the fill value (name='"//trim(prog_meta(iVar)%varName)//"')"
   err=20; return
  endif

  ! store data in prognostics structure 
  ! loop through GRUs
  do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount

    ! get ther number of layers
    nSnow = gru_struc(iGRU)%hruInfo(iHRU)%nSnow 
    nSoil = gru_struc(iGRU)%hruInfo(iHRU)%nSoil 
    nToto = nSnow + nSoil 
   
    ! put the data into data structures
    select case (prog_meta(iVar)%varType)
     case (iLookVarType%scalarv); progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1        ) = varData(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc,1        ) 
     case (iLookVarType%midSoil); progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1:nSoil  ) = varData(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc,1:nSoil  ) 
     case (iLookVarType%midToto); progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1:nToto  ) = varData(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc,1:nToto  ) 
     case (iLookVarType%ifcToto); progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(0:nToto  ) = varData(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc,1:nToto+1) 
     case default
      message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
      err=20; return
    end select

    ! initialize the spectral albedo
    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nBand) = progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1)

   end do ! iHRU
  end do ! iGRU

  ! deallocate storage vector for next variable
  deallocate(varData, stat=err)
  if(err/=0)then; message=trim(message)//'problem deallocating variable data'; return; endif

 end do ! iVar 

 ! --------------------------------------------------------------------------------------------------------
 ! (2) set number of layers 
 ! --------------------------------------------------------------------------------------------------------
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount

   ! save the number of layers
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)   = gru_struc(iGRU)%hruInfo(iHRU)%nSnow 
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)   = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nLayers)%dat(1) = gru_struc(iGRU)%hruInfo(iHRU)%nSnow + gru_struc(iGRU)%hruInfo(iHRU)%nSoil

   ! initalize the indices for midSnow, midSoil, midToto, and ifcToto
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midSnowStartIndex)%dat(1) = 1
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midSoilStartIndex)%dat(1) = 1
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midTotoStartIndex)%dat(1) = 1
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcSnowStartIndex)%dat(1) = 1
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcSoilStartIndex)%dat(1) = 1
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcTotoStartIndex)%dat(1) = 1

   ! set layer type
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat(1:gru_struc(iGRU)%hruInfo(iHRU)%nSnow) = iname_snow
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat((gru_struc(iGRU)%hruInfo(iHRU)%nSnow+1):(gru_struc(iGRU)%hruInfo(iHRU)%nSnow+gru_struc(iGRU)%hruInfo(iHRU)%nSoil)) = iname_soil

  end do
 end do

 end subroutine read_icond

end module read_icond_module
