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

module modelwrite_module
USE netcdf
USE netcdf_util_module,only:netcdf_err                    ! netcdf error handling function
USE nrtype, integerMissing=>nr_integerMissing             ! top-level data types
implicit none
private
public::writeParm
public::writeData
public::writeBasin
public::writeTime
public::writeRestart
! define dimension lengths
integer(i4b),parameter      :: maxSpectral=2              ! maximum number of spectral bands
contains

 ! **********************************************************************************************************
 ! public subroutine writeParm: write model parameters
 ! **********************************************************************************************************
 subroutine writeParm(iHRU,struct,meta,err,message)
 USE globalData,only:ncid                        ! netcdf file ids
 USE globalData,only:integerMissing              ! missing value
 USE data_types,only:var_info                    ! metadata info
 USE data_types,only:var_i,var_d,var_dlength     ! derived data types
 USE var_lookup,only:iLookStat                   ! to index into write flag
 implicit none

 ! declare input variables
 integer(i4b)  ,intent(in)   :: iHRU             ! hydrologic response unit
 class(*)      ,intent(in)   :: struct           ! data structure
 type(var_info),intent(in)   :: meta(:)          ! metadata structure
 integer(i4b)  ,intent(out)  :: err              ! error code
 character(*)  ,intent(out)  :: message          ! error message
 ! local variables
 integer(i4b)                :: iVar             ! loop through variables
 integer(i4b)  ,parameter    :: modelTime=1      ! these particular data are only output in the timestep file

 ! initialize error control
 err=0;message="f-writeParm/"

 ! loop through local column model parameters
 do iVar = 1,size(meta)

  ! check that the variable is desired
  if (.not.meta(iVar)%statFlag(iLookStat%inst)) cycle

  ! initialize message
  message=trim(message)//trim(meta(iVar)%varName)//'/'

  ! write data
  if (iHRU.ne.integerMissing) then
   select type (struct)
    type is (var_i)
     err = nf90_put_var(ncid(modelTime),meta(iVar)%ncVarID(iLookStat%inst),(/struct%var(iVar)/),start=(/iHRU/),count=(/1/))
    type is (var_d)
     err = nf90_put_var(ncid(modelTime),meta(iVar)%ncVarID(iLookStat%inst),(/struct%var(iVar)/),start=(/iHRU/),count=(/1/))
    type is (var_dlength)
     err = nf90_put_var(ncid(modelTime),meta(iVar)%ncVarID(iLookStat%inst),(/struct%var(iVar)%dat/),start=(/iHRU,1/),count=(/1,size(struct%var(iVar)%dat)/))
    class default; err=20; message=trim(message)//'unkonwn variable type (with HRU)'; return
   end select
  else
   select type (struct)
    type is (var_d)
     err = nf90_put_var(ncid(modelTime),meta(iVar)%ncVarID(iLookStat%inst),(/struct%var(iVar)/),start=(/1/),count=(/1/))
    class default; err=20; message=trim(message)//'unkonwn variable type (no HRU)'; return
   end select
  end if
  call netcdf_err(err,message); if (err/=0) return

  ! re-initialize message
  message="f-writeParm/"
 end do  ! looping through local column model parameters

 end subroutine writeParm

 ! **************************************************************************************
 ! public subroutine writeData: write model time-dependent data
 ! **************************************************************************************
 subroutine writeData(modelTimestep,outputTimestep,meta,stat,dat,map,indx,iHRU,err,message)
 USE data_types,only:var_info,dlength,ilength       ! type structures for passing
 USE var_lookup,only:maxVarStat                     ! index into stats structure
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE var_lookup,only:iLookIndex                     ! index into index structure
 USE var_lookup,only:iLookStat                      ! index into stat structure
 USE globalData,only:outFreq,nFreq,ncid             ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none

 ! declare dummy variables
 type(var_info),intent(in)     :: meta(:)           ! meta data
 class(*)      ,intent(in)     :: stat(:)           ! stats data
 class(*)      ,intent(in)     :: dat(:)            ! timestep data
 type(ilength) ,intent(in)     :: indx(:)           ! index data
 integer(i4b)  ,intent(in)     :: map(:)            ! map into stats child struct
 integer(i4b)  ,intent(in)     :: iHRU              ! hydrologic response unit
 integer(i4b)  ,intent(in)     :: modelTimestep     ! model time step
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iStat             ! statistics index
 integer(i4b)                  :: iFreq             ! frequency index
 integer(i4b)                  :: ncVarID           ! used only for time
 integer(i4b)                  :: nSnow             ! number of snow layers
 integer(i4b)                  :: nSoil             ! number of soil layers
 integer(i4b)                  :: nLayers           ! total number of layers
 integer(i4b)                  :: midSnowStartIndex ! start index of the midSnow vector for a given timestep
 integer(i4b)                  :: midSoilStartIndex ! start index of the midSoil vector for a given timestep
 integer(i4b)                  :: midTotoStartIndex ! start index of the midToto vector for a given timestep
 integer(i4b)                  :: ifcSnowStartIndex ! start index of the ifcSnow vector for a given timestep
 integer(i4b)                  :: ifcSoilStartIndex ! start index of the ifcSoil vector for a given timestep
 integer(i4b)                  :: ifcTotoStartIndex ! start index of the ifcToto vector for a given timestep

 ! initialize error control
 err=0;message="writeData/"

 ! model layers
 nSoil             = indx(iLookIndex%nSoil)%dat(1)
 nSnow             = indx(iLookIndex%nSnow)%dat(1)
 nLayers           = indx(iLookIndex%nLayers)%dat(1)
 ! model indices
 midSnowStartIndex = indx(iLookIndex%midSnowStartIndex)%dat(1)
 midSoilStartIndex = indx(iLookIndex%midSoilStartIndex)%dat(1)
 midTotoStartIndex = indx(iLookIndex%midTotoStartIndex)%dat(1)
 ifcSnowStartIndex = indx(iLookIndex%ifcSnowStartIndex)%dat(1)
 ifcSoilStartIndex = indx(iLookIndex%ifcSoilStartIndex)%dat(1)
 ifcTotoStartIndex = indx(iLookIndex%ifcTotoStartIndex)%dat(1)

 ! loop through output frequencies
 do iFreq = 1,nFreq

  ! check that the timestep is desired
  if (mod(modelTimestep,outFreq(iFreq)).ne.0) cycle

   ! loop through model variables
   do iVar = 1,size(meta)

    ! handle time first
    if (meta(iVar)%varName=='time') then    
     select type(stat)
      type is (dlength)
       err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID) 
       call netcdf_err(err,message); if (err/=0) return
       err = nf90_put_var(ncid(iFreq),ncVarID,(/stat(iVar)%dat(iLookStat%inst)/),start=(/outputTimestep(iFreq)/),count=(/1,1/))
       call netcdf_err(err,message); if (err/=0) return
       cycle
     class default; err=20; message=trim(message)//'time variable must be of type dlength'; return; 
     end select
    end if

    ! check that the variable is desired
    if (meta(iVar)%outFreq.ne.iFreq) cycle

    ! loop through output stats
    do iStat = 1,maxVarStat

     ! check that the variable is desired
     if ((.not.meta(iVar)%statFlag(iStat)).or.(trim(meta(iVar)%varName)=='unknown')) cycle

     ! stats/data output - select data type
     if (meta(iVar)%varType==iLookVarType%scalarv) then
       select type(stat)
        type is (ilength)
         err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/stat(map(iVar))%dat(iStat)/),start=(/iHRU,outputTimestep(iFreq)/),count=(/1,1/))
        type is (dlength)
         err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/stat(map(iVar))%dat(iStat)/),start=(/iHRU,outputTimestep(iFreq)/),count=(/1,1/))
        class default; err=20; message=trim(message)//'stats must be scalarv and either ilength of dlength'; return
       end select  ! stat

     ! non-scalar variables
     else
      select type (dat)
       type is (dlength)
        select case (meta(iVar)%varType)
         case(iLookVarType%wLength); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,1,outputTimestep(iFreq)/),count=(/1,maxSpectral,1/))
         case(iLookVarType%midToto); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midTotoStartIndex/),count=(/1,nLayers/))
         case(iLookVarType%midSnow); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midSnowStartIndex/),count=(/1,nSnow/))
         case(iLookVarType%midSoil); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midSoilStartIndex/),count=(/1,nSoil/))
         case(iLookVarType%ifcToto); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcTotoStartIndex/),count=(/1,nLayers+1/))
         case(iLookVarType%ifcSnow); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcSnowStartIndex/),count=(/1,nSnow+1/))
         case(iLookVarType%ifcSoil); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcSoilStartIndex/),count=(/1,nSoil+1/))
        end select ! vartype
       type is (ilength)
        select case (meta(iVar)%varType)
         case(iLookVarType%wLength); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,1,outputTimestep(iFreq)/),count=(/1,maxSpectral,1/))
         case(iLookVarType%midToto); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midTotoStartIndex/),count=(/1,nLayers/))
         case(iLookVarType%midSnow); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midSnowStartIndex/),count=(/1,nSnow/))
         case(iLookVarType%midSoil); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,midSoilStartIndex/),count=(/1,nSoil/))
         case(iLookVarType%ifcToto); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcTotoStartIndex/),count=(/1,nLayers+1/))
         case(iLookVarType%ifcSnow); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcSnowStartIndex/),count=(/1,nSnow+1/))
         case(iLookVarType%ifcSoil); err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/iHRU,ifcSoilStartIndex/),count=(/1,nSoil+1/))
        end select ! vartype
      end select ! dat
     end if ! sacalarv

     ! process error code
     if (err.ne.0) message=trim(message)//trim(meta(iVar)%varName)//'_'//trim(get_statName(iStat))
     call netcdf_err(err,message); if (err/=0) return

    end do ! iStat
   end do ! iVar
  end do ! iFreq

 end subroutine writeData

 ! **************************************************************************************
 ! public subroutine writeBasin: write basin-average variables
 ! **************************************************************************************
 subroutine writeBasin(modelTimestep,outputTimestep,meta,stat,dat,map,err,message)
 USE data_types,only:var_info,dlength,ilength       ! type structures for passing
 USE var_lookup,only:maxVarStat                     ! index into stats structure
 USE var_lookup,only:iLookVarType                   ! index into type structure
 USE globalData,only:outFreq,nFreq,ncid             ! output file information
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixName_module,only:get_statName            ! to access type strings for error messages
 implicit none

 ! declare dummy variables
 type(var_info),intent(in)     :: meta(:)           ! meta data
 type(dlength) ,intent(in)     :: stat(:)           ! stats data
 type(dlength) ,intent(in)     :: dat(:)            ! timestep data
 integer(i4b)  ,intent(in)     :: map(:)            ! map into stats child struct
 integer(i4b)  ,intent(in)     :: modelTimestep     ! model time step
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iStat             ! statistics index
 integer(i4b)                  :: iFreq             ! frequency index
 ! initialize error control
 err=0;message="f-writeBasin/"

 do iFreq = 1,nFreq
  ! check that the timestep is desired
  if (mod(modelTimestep,outFreq(iFreq)).ne.0) cycle

   ! loop through model variables
   do iVar = 1,size(meta)

    ! check that the variable is desired
    if (meta(iVar)%outFreq.ne.iFreq) cycle

    ! loop through output stats
    do iStat = 1,maxVarStat
     ! check that the variable is desired
     if ((.not.meta(iVar)%statFlag(iStat)).or.(trim(meta(iVar)%varName)=='unknown')) cycle

     ! stats/dats output - select data type
     select case (meta(iVar)%varType)

      case (iLookVarType%scalarv)
       err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/stat(map(iVar))%dat(iStat)/),start=(/outputTimestep(iFreq)/),count=(/1/))

      case (iLookVarType%routing)
       if (modelTimestep==1) then
        err = nf90_put_var(ncid(iFreq),meta(iVar)%ncVarID(iStat),(/dat(iVar)%dat/),start=(/1/),count=(/1000/))
       end if

      case default
       err=40; message=trim(message)//"unknownVariableType[name='"//trim(meta(iVar)%varName)//"';type='"//trim(get_varTypeName(meta(iVar)%varType))//    "']"; return
      end select ! variable type

     ! process error code
     if (err.ne.0) message=trim(message)//trim(meta(iVar)%varName)//'_'//trim(get_statName    (iStat))
     call netcdf_err(err,message); if (err/=0) return

    end do ! iStat
   end do ! iVar
  end do ! iFreq

 end subroutine writeBasin

 ! **************************************************************************************
 ! public subroutine writeTime: write current time to all files 
 ! **************************************************************************************
 subroutine writeTime(modelTimestep,outputTimestep,meta,dat,err,message)
 USE data_types,only:var_info,dlength,ilength       ! type structures for passing
 USE globalData,only:outFreq,nFreq,ncid             ! output file information
 USE var_lookup,only:iLookStat                      ! index into stat structure
 implicit none

 ! declare dummy variables
 type(var_info),intent(in)     :: meta(:)           ! meta data
 integer       ,intent(in)     :: dat(:)            ! timestep data
 integer(i4b)  ,intent(in)     :: modelTimestep     ! model time step
 integer(i4b)  ,intent(in)     :: outputTimestep(:) ! output time step
 integer(i4b)  ,intent(out)    :: err               ! error code
 character(*)  ,intent(out)    :: message           ! error message
 ! local variables
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iFreq             ! frequency index
 integer(i4b)                  :: ncVarID           ! used only for time
 ! initialize error control
 err=0;message="f-writeTime/"

 do iFreq = 1,nFreq
  ! check that the timestep is desired
  if (mod(modelTimestep,outFreq(iFreq)).ne.0) cycle

   ! loop through model variables
   do iVar = 1,size(meta)

    ! if variable is desired
    if (.not.meta(iVar)%statFlag(iLookStat%inst)) cycle

    ! get variable id in file
    err = nf90_inq_varid(ncid(iFreq),trim(meta(iVar)%varName),ncVarID) 
    if (err/=0) message=trim(message)//trim(meta(iVar)%varName)
    call netcdf_err(err,message)
    if (err/=0) then; err=20; return; end if

    ! add to file
    err = nf90_put_var(ncid(iFreq),ncVarID,(/dat(iVar)/),start=(/outputTimestep(iFreq)/),count=(/1/))
    if (err/=0) message=trim(message)//trim(meta(iVar)%varName)
    call netcdf_err(err,message)
    if (err/=0) then; err=20; return; end if

   end do ! iVar
  end do ! iFreq

 end subroutine writeTime 

 ! *********************************************************************************************************
 ! public subroutine printRestartFile: print a re-start file
 ! *********************************************************************************************************
 subroutine writeRestart(filename,         & ! intent(in): name of restart file
                         nGRU,             & ! intent(in): number of GRUs
                         nHRU,             & ! intent(in): number of HRUs
                         prog_meta,        & ! intent(in): prognostics metadata 
                         prog_data,        & ! intent(in): prognostics data
                         indx_meta,        & ! intent(in): index metadata 
                         indx_data,        & ! intent(in): index data
                         err,message)        ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 ! access the derived types to define the data structures
! USE data_types,only:hru_d                  ! length of substep
 USE data_types,only:gru_hru_doubleVec      ! actual data
 USE data_types,only:gru_hru_intVec         ! actual data
 USE data_types,only:var_info               ! metadata 
 ! access named variables defining elements in the data structures
 USE var_lookup,only:iLookINDEX             ! named variables for structure elements
 USE var_lookup,only:iLookVarType           ! named variables for structure elements
 ! constants
 USE globalData,only:gru_struc              ! gru-hru mapping structures
 ! external routines
 USE netcdf_util_module,only:nc_file_close  ! close netcdf file
 USE netcdf_util_module,only:nc_file_open   ! open netcdf file
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input
 character(len=256),intent(in)      :: filename      ! name of the restart file
 integer(i4b),intent(in)            :: nGRU          ! number of GRUs
 integer(i4b),intent(in)            :: nHRU          ! number of HRUs
 type(var_info),intent(in)          :: prog_meta(:)  ! metadata 
 type(gru_hru_doubleVec),intent(in) :: prog_data     ! prognostic vars 
 type(var_info),intent(in)          :: indx_meta(:)  ! metadata 
 type(gru_hru_intVec),intent(in)    :: indx_data     ! indexing vars 
 ! output: error control
 integer(i4b),intent(out)           :: err           ! error code
 character(*),intent(out)           :: message       ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                       :: ncid          ! netcdf file id
 integer(i4b),allocatable           :: ncVarID(:)    ! netcdf variable id
 integer(i4b)                       :: ncSnowID      ! index variable id
 integer(i4b)                       :: ncSoilID      ! index variable id

 integer(i4b)                       :: nSoil         ! number of soil layers
 integer(i4b)                       :: nSnow         ! number of snow layers
 integer(i4b)                       :: maxSnow       ! maximum number of snow layers
 integer(i4b)                       :: maxSoil       ! maximum number of soil layers
 integer(i4b)                       :: nLayers       ! number of total layers
 integer(i4b)                       :: maxLayers     ! maximum number of total layers
 integer(i4b),parameter             :: nSpectral=2   ! number of spectal bands
 integer(i4b),parameter             :: nScalar=1     ! size of a scalar

 integer(i4b)                       :: hruDimID      ! variable dimension ID
 integer(i4b)                       :: scalDimID     ! variable dimension ID
 integer(i4b)                       :: specDimID     ! variable dimension ID
 integer(i4b)                       :: midSnowDimID  ! variable dimension ID
 integer(i4b)                       :: midSoilDimID  ! variable dimension ID
 integer(i4b)                       :: midTotoDimID  ! variable dimension ID
 integer(i4b)                       :: ifcSnowDimID  ! variable dimension ID
 integer(i4b)                       :: ifcSoilDimID  ! variable dimension ID
 integer(i4b)                       :: ifcTotoDimID  ! variable dimension ID

 character(len=32),parameter        :: hruDimName    ='hru'      ! dimension name for HRUs
 character(len=32),parameter        :: scalDimName   ='scalarv'  ! dimension name for scalar data
 character(len=32),parameter        :: specDimName   ='spectral' ! dimension name for spectral bands
 character(len=32),parameter        :: midSnowDimName='midSnow'  ! dimension name for snow-only layers
 character(len=32),parameter        :: midSoilDimName='midSoil'  ! dimension name for soil-only layers
 character(len=32),parameter        :: midTotoDimName='midToto'  ! dimension name for layered varaiables
 character(len=32),parameter        :: ifcSnowDimName='ifcSnow'  ! dimension name for snow-only layers
 character(len=32),parameter        :: ifcSoilDimName='ifcSoil'  ! dimension name for soil-only layers
 character(len=32),parameter        :: ifcTotoDimName='ifcToto'  ! dimension name for layered varaiables

 integer(i4b)                       :: cHRU          ! count of HRUs
 integer(i4b)                       :: iHRU          ! index of HRUs
 integer(i4b)                       :: iGRU          ! index of GRUs
 integer(i4b)                       :: iVar          ! variable index
 logical(lgt)                       :: okLength      ! flag to check if the vector length is OK
 character(len=256)                 :: cmessage      ! downstream error message
 ! --------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='writeRestart/'

 ! size of prog vector
 allocate(ncVarID(size(prog_meta)))

 ! maximum number of soil layers
 maxSoil = 0
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   maxSoil = max(maxSoil,gru_struc(iGRU)%hruInfo(iHRU)%nSoil)
  end do
 end do

 ! maximum number of snow layers
 maxSnow = 0
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   maxSnow = max(maxSnow,gru_struc(iGRU)%hruInfo(iHRU)%nSnow)
  end do
 end do
 
 ! total number of layers
 maxLayers = maxSnow+maxSoil

 ! create file 
 err = nf90_create(trim(filename),nf90_classic_model,ncid)
 message='iCreate[create]'; call netcdf_err(err,message); if(err/=0)return

 ! define dimensions
                err = nf90_def_dim(ncid,trim(hruDimName)    ,nHRU       ,   hruDimID) ; message='iCreate[hru]'     ;call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(scalDimName)   ,nScalar    ,   scalDimID); message='iCreate[scalar]'  ;call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(specDimName)   ,nSpectral  ,   specDimID); message='iCreate[spectral]';call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(midSoilDimName),maxSoil    ,midSoilDimID); message='iCreate[ifcSoil]' ;call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(midTotoDimName),maxLayers  ,midTotoDimID); message='iCreate[midToto]' ;call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(ifcSoilDimName),maxSoil+1  ,ifcSoilDimID); message='iCreate[ifcSoil]' ;call netcdf_err(err,message); if(err/=0)return
                err = nf90_def_dim(ncid,trim(ifcTotoDimName),maxLayers+1,ifcTotoDimID); message='iCreate[ifcToto]' ;call netcdf_err(err,message); if(err/=0)return
 if (maxSnow>0) err = nf90_def_dim(ncid,trim(midSnowDimName),maxSnow    ,midSnowDimID); message='iCreate[ifcSnow]' ;call netcdf_err(err,message); if(err/=0)return
 if (maxSnow>0) err = nf90_def_dim(ncid,trim(ifcSnowDimName),maxSnow+1  ,ifcSnowDimID); message='iCreate[ifcSnow]' ;call netcdf_err(err,message); if(err/=0)return
 ! re-initialize error control
 err=0; message='writeRestart/'

 ! define prognostic variables
 do iVar = 1,size(prog_meta)
  if (prog_meta(iVar)%varType==iLookvarType%unknown) cycle

  ! define variable
  select case(prog_meta(iVar)%varType)
   case(iLookvarType%scalarv);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,  scalDimID /),ncVarID(iVar)) 
   case(iLookvarType%wLength);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,  specDimID /),ncVarID(iVar)) 
   case(iLookvarType%midSoil);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,midSoilDimID/),ncVarID(iVar)) 
   case(iLookvarType%midToto);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,midTotoDimID/),ncVarID(iVar)) 
   case(iLookvarType%ifcSoil);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,ifcSoilDimID/),ncVarID(iVar)) 
   case(iLookvarType%ifcToto);                err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,ifcTotoDimID/),ncVarID(iVar)) 
   case(iLookvarType%midSnow); if (maxSnow>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,midSnowDimID/),ncVarID(iVar)) 
   case(iLookvarType%ifcSnow); if (maxSnow>0) err = nf90_def_var(ncid,trim(prog_meta(iVar)%varname),nf90_double,(/hruDimID,ifcSnowDimID/),ncVarID(iVar)) 
  end select
 
  ! check errors
  if(err/=0)then
   message=trim(message)//trim(cmessage)//' [variable '//trim(prog_meta(iVar)%varName)//']'
   return
  end if

  ! add parameter description
  err = nf90_put_att(ncid,ncVarID(iVar),'long_name',trim(prog_meta(iVar)%vardesc))
  call netcdf_err(err,message)

  ! add parameter units
  err = nf90_put_att(ncid,ncVarID(iVar),'units',trim(prog_meta(iVar)%varunit))
  call netcdf_err(err,message)

 end do ! iVar 

 ! define index variables - snow
 err = nf90_def_var(ncid,trim(indx_meta(iLookIndex%nSnow)%varName),nf90_int,(/hruDimID/),ncSnowID); call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSnowID,'long_name',trim(indx_meta(iLookIndex%nSnow)%vardesc));           call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSnowID,'units'    ,trim(indx_meta(iLookIndex%nSnow)%varunit));           call netcdf_err(err,message)

 ! define index variables - soil
 err = nf90_def_var(ncid,trim(indx_meta(iLookIndex%nSoil)%varName),nf90_int,(/hruDimID/),ncSoilID); call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSoilID,'long_name',trim(indx_meta(iLookIndex%nSoil)%vardesc));           call netcdf_err(err,message)
 err = nf90_put_att(ncid,ncSoilID,'units'    ,trim(indx_meta(iLookIndex%nSoil)%varunit));           call netcdf_err(err,message)

 ! end definition phase
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return

 ! write variables
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   cHRU = gru_struc(iGRU)%hruInfo(iHRU)%hru_ix 
   do iVar = 1,size(prog_meta)

    ! excape if this variable is not used
    if (prog_meta(iVar)%varType==iLookvarType%unknown) cycle

    ! actual number of layers
    nSnow = gru_struc(iGRU)%hruInfo(iHRU)%nSnow
    nSoil = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
    nLayers = nSoil + nSnow

    ! check size
    ! NOTE: this may take time that we do not wish to use
    okLength=.true.
    select case (prog_meta(iVar)%varType)
     case(iLookVarType%scalarv);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nScalar  )
     case(iLookVarType%wlength);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSpectral)
     case(iLookVarType%midSoil);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSoil    )
     case(iLookVarType%midToto);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nLayers  )
     case(iLookVarType%ifcSoil);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSoil+1  )
     case(iLookVarType%ifcToto);              okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nLayers+1)
     case(iLookVarType%midSnow); if (nSnow>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSnow    )
     case(iLookVarType%ifcSnow); if (nSnow>0) okLength = (size(prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat) == nSnow+1  )
     case default; err=20; message=trim(message)//'unknown var type'; return
    end select 

    ! error check
    if(.not.okLength)then
     message=trim(message)//'bad vector length for variable '//trim(prog_meta(iVar)%varname)
     err=20; return
    endif

    ! write data 
    select case (prog_meta(iVar)%varType)
     case(iLookVarType%scalarv);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nScalar  /))
     case(iLookVarType%wlength);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSpectral/))
     case(iLookVarType%midSoil);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSoil    /))
     case(iLookVarType%midToto);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nLayers  /))
     case(iLookVarType%ifcSoil);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSoil+1  /))
     case(iLookVarType%ifcToto);              err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nLayers+1/))
     case(iLookVarType%midSnow); if (nSnow>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSnow    /))
     case(iLookVarType%ifcSnow); if (nSnow>0) err=nf90_put_var(ncid,ncVarID(iVar),(/prog_data%gru(iGRU)%hru(iHRU)%var(iVar)%dat/),start=(/cHRU,1/),count=(/1,nSnow+1  /))
     case default; err=20; message=trim(message)//'unknown var type'; return
    end select 

    ! error check
    if (err.ne.0) message=trim(message)//'writing variable:'//trim(prog_meta(iVar)%varName)
    call netcdf_err(err,message); if (err/=0) return
    err=0; message='writeRestart/'

   end do ! iVar 

   ! write index variables 
   err=nf90_put_var(ncid,ncSnowID,(/indx_data%gru(iGRU)%hru(iHRU)%var(iLookIndex%nSnow)%dat/),start=(/cHRU/),count=(/1/))
   err=nf90_put_var(ncid,ncSoilID,(/indx_data%gru(iGRU)%hru(iHRU)%var(iLookIndex%nSoil)%dat/),start=(/cHRU/),count=(/1/))
 
  end do ! iGRU
 end do ! iHRU

 ! close file 
 call nc_file_close(ncid,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(ncVarID)

 end subroutine writeRestart

end module modelwrite_module
