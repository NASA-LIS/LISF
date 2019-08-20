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

module def_output_module
USE netcdf
USE netcdf_util_module,only:netcdf_err        ! netcdf error handling function
USE nrtype, integerMissing=>nr_integerMissing ! top-level data types
USE f2008funcs_module,only:cloneStruc         ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
implicit none
private
public :: def_output

! define dimension names
character(len=32),parameter :: hru_DimName            = 'hru'              ! dimension name for the HRUs
character(len=32),parameter :: depth_DimName          = 'depth'            ! dimension name for soil depth
character(len=32),parameter :: scalar_DimName         = 'scalar'           ! dimension name for scalar variables
character(len=32),parameter :: wLength_dimName        = 'spectral_bands'   ! dimension name for the number of spectral bands
character(len=32),parameter :: timestep_DimName       = 'time'             ! dimension name for the time step
character(len=32),parameter :: routing_DimName        = 'timeDelayRouting' ! dimension name for thetime delay routing vectors
character(len=32),parameter :: midSnowAndTime_DimName = 'midSnowAndTime'   ! dimension name for midSnow-time
character(len=32),parameter :: midSoilAndTime_DimName = 'midSoilAndTime'   ! dimension name for midSoil-time
character(len=32),parameter :: midTotoAndTime_DimName = 'midTotoAndTime'   ! dimension name for midToto-time
character(len=32),parameter :: ifcSnowAndTime_DimName = 'ifcSnowAndTime'   ! dimension name for ifcSnow-time
character(len=32),parameter :: ifcSoilAndTime_DimName = 'ifcSoilAndTime'   ! dimension name for ifcSoil-time
character(len=32),parameter :: ifcTotoAndTime_DimName = 'ifcTotoAndTime'   ! dimension name for ifcToto-time

! define the dimension IDs
integer(i4b)                :: hru_DimID                               ! dimension name for the HRUs
integer(i4b)                :: depth_DimID                             ! dimension name for the soil depth
integer(i4b)                :: scalar_DimID                            ! dimension name for scalar variables
integer(i4b)                :: wLength_dimID                           ! dimension name for the number of spectral bands
integer(i4b)                :: timestep_DimID                          ! dimension name for the time step
integer(i4b)                :: routing_DimID                           ! dimension name for thetime delay routing vectors
integer(i4b)                :: midSnowAndTime_DimID                    ! dimension name for midSnow-time
integer(i4b)                :: midSoilAndTime_DimID                    ! dimension name for midSoil-time
integer(i4b)                :: midTotoAndTime_DimID                    ! dimension name for midToto-time
integer(i4b)                :: ifcSnowAndTime_DimID                    ! dimension name for ifcSnow-time
integer(i4b)                :: ifcSoilAndTime_DimID                    ! dimension name for ifcSoil-time
integer(i4b)                :: ifcTotoAndTime_DimID                    ! dimension name for ifcToto-time

! define named variables to specify dimensions
integer(i4b),parameter  :: needHRU=1,noHRU=2    ! define if there is an HRU dimension
integer(i4b),parameter  :: needTime=1,noTime=2  ! define if there is a time dimension

contains

 ! **********************************************************************************************************
 ! public subroutine def_output: define model output file
 ! **********************************************************************************************************
 subroutine def_output(nHRU,nSoil,infile,err,message)
 USE globalData,only:structInfo                               ! information on the data structures
 USE globalData,only:forc_meta,attr_meta,type_meta            ! metaData structures
 USE globalData,only:prog_meta,diag_meta,flux_meta,deriv_meta ! metaData structures
 USE globalData,only:mpar_meta,indx_meta                      ! metaData structures
 USE globalData,only:bpar_meta,bvar_meta,time_meta            ! metaData structures
 USE globalData,only:ncid
 USE globalData,only:nFreq,outFreq                            ! output frequencies
 ! declare dummy variables
 integer(i4b),intent(in)     :: nHRU                          ! number of HRUs
 integer(i4b),intent(in)     :: nSoil                         ! number of soil layers in the first HRU (used to define fixed length dimensions)
 character(*),intent(in)     :: infile                        ! file suffix
 integer(i4b),intent(out)    :: err                           ! error code
 character(*),intent(out)    :: message                       ! error message
 ! local variables
 integer(i4b)                :: iFreq                         ! loop through output frequencies
 integer(i4b)                :: iStruct                       ! loop through structure types 
 integer(i4b),parameter      :: modelTime=1                   ! model timestep output frequency
 character(len=5)            :: fstring                       ! string to hold model output freuqnecy
 character(len=1000)         :: fname                         ! temporary filename
 character(len=256)          :: cmessage                      ! temporary error message

 ! initialize errors
 err=0; message="def_output/"

 ! create initial file
 ! each file will have a master name with a frequency appended at the end:
 ! e.g., xxxxxxxxx_1.nc  (for output at every model timestep)
 ! e.g., xxxxxxxxx_24.nc (for daily output with hourly model timestep)
 do iFreq = 1,nFreq

  ! create file
  write(fstring,'(i5)') outFreq(iFreq)
  fstring = adjustl(fstring)
  fname = trim(infile)//'_'//trim(fstring)//'.nc'
  call ini_create(nHRU,nSoil,trim(fname),ncid(iFreq),err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  print "(A,A)",'Created output file:',trim(fname)

  ! ensure that all time variables are written to all files
  time_meta(:)%outFreq = iFreq

  do iStruct = 1,size(structInfo)
   select case (trim(structInfo(iStruct)%structName))
    case('attr' ); call def_variab(ncid(iFreq),iFreq,needHRU,  noTime,attr_meta, nf90_double,err,cmessage)  ! local attributes HRU
    case('type' ); call def_variab(ncid(iFreq),iFreq,needHRU,  noTime,type_meta, nf90_int,   err,cmessage)  ! local classification
    case('mpar' ); call def_variab(ncid(iFreq),iFreq,needHRU,  noTime,mpar_meta, nf90_double,err,cmessage)  ! model parameters
    case('bpar' ); call def_variab(ncid(iFreq),iFreq,  noHRU,  noTime,bpar_meta, nf90_double,err,cmessage)  ! basin-average param
    case('indx' ); call def_variab(ncid(iFreq),iFreq,needHRU,needTime,indx_meta, nf90_int,   err,cmessage)  ! model variables
    case('deriv'); call def_variab(ncid(iFreq),iFreq,needHRU,needTime,deriv_meta,nf90_double,err,cmessage)  ! model derivatives
    case('time' ); call def_variab(ncid(iFreq),iFreq,  noHRU,needTime,time_meta,nf90_int,    err,cmessage)  ! model derivatives
    case('forc' ); call def_variab(ncid(iFreq),iFreq,needHRU,needTime,forc_meta, nf90_double,err,cmessage)  ! model forcing data
    case('prog' ); call def_variab(ncid(iFreq),iFreq,needHRU,needTime,prog_meta, nf90_double,err,cmessage)  ! model prognostics
    case('diag' ); call def_variab(ncid(iFreq),iFreq,needHRU,needTime,diag_meta, nf90_double,err,cmessage)  ! model diagnostic variables
    case('flux' ); call def_variab(ncid(iFreq),iFreq,needHRU,needTime,flux_meta, nf90_double,err,cmessage)  ! model fluxes
    case('bvar' ); call def_variab(ncid(iFreq),iFreq,  noHRU,needTime,bvar_meta, nf90_double,err,cmessage)  ! basin-average variables
    case default; err=20; message=trim(message)//'unable to identify lookup structure';
   end select
   ! error handling
   if(err/=0)then;err=20;message=trim(message)//trim(cmessage)//'[structure =  '//trim(structInfo(iStruct)%structName);return;end if
  end do ! iStruct 

  ! write HRU dimension for each output file
  call write_hru_dim(ncid(iFreq), err, cmessage); if(err/=0) then; message=trim(message)//trim(cmessage); return; end if
 end do ! iFreq 

 end subroutine def_output

 ! **********************************************************************************************************
 ! private subroutine ini_create: initial create
 ! **********************************************************************************************************
 subroutine ini_create(nHRU,nSoil,infile,ncid,err,message)
 ! variables to define number of steps per file (total number of time steps, step length, etc.)
 USE multiconst,only:secprday           ! number of seconds per day
 USE globalData,only:data_step          ! time step of model forcing data (s)
 USE globalData,only:numtim             ! number of time steps
 ! model decisions
 USE globalData,only:model_decisions    ! model decision structure
 USE var_lookup,only:iLookDECISIONS     ! named variables for elements of the decision structure
 USE mDecisions_module,only:&
  sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
  rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
 implicit none
 ! declare dummy variables
 integer(i4b),intent(in)     :: nHRU                       ! number of HRUs
 integer(i4b),intent(in)     :: nSoil                      ! number of soil layers in the first HRU (used to define fixed length dimensions)
 character(*),intent(in)     :: infile                     ! filename
 integer(i4b),intent(out)    :: ncid                       ! netcdf file id
 integer(i4b),intent(out)    :: err                        ! error code
 character(*),intent(out)    :: message                    ! error message
 ! define local variables
 integer(i4b)                :: maxRouting=1000            ! maximum length of routing vector
 integer(i4b),parameter      :: maxSpectral=2              ! maximum number of spectral bands
 integer(i4b),parameter      :: scalarLength=1             ! length of scalar variable
 integer(i4b)                :: meanSnowLayersPerStep      ! mean number of snow layers per time step
 integer(i4b)                :: maxStepsPerFile            ! maximum number of time steps to be stored in each file
 integer(i4b)                :: maxLength                  ! maximum length of the variable vector
 ! initialize error control
 err=0;message="f-iniCreate/"
 ! identify length of the variable vector
 maxStepsPerFile = min(numtim,nint(366._dp * secprday/data_step))
 if(maxStepsPerFile < numtim) maxStepsPerFile=numtim
 select case(model_decisions(iLookDECISIONS%snowLayers)%iDecision)
  case(sameRulesAllLayers);    meanSnowLayersPerStep = 100
  case(rulesDependLayerIndex); meanSnowLayersPerStep = 5
  case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
 end select ! (option to combine/sub-divide snow layers)
 maxLength = maxStepsPerFile*(nSoil+1 + meanSnowLayersPerStep)
 print*, 'maxStepsPerFile, maxLength = ', maxStepsPerFile, maxLength

 ! create output file
 err = nf90_create(trim(infile),nf90_classic_model,ncid)
 message='iCreate[create]'; call netcdf_err(err,message); if (err/=0) return

 ! create dimensions
 err = nf90_def_dim(ncid, trim(           hru_DimName), nHRU,                  hru_DimID); message='iCreate[HRU]';      call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(      timestep_DimName), nf90_unlimited,   timestep_DimID); message='iCreate[time]';     call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(         depth_DimName), nSoil,               depth_DimID); message='iCreate[depth]';    call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(        scalar_DimName), scalarLength,       scalar_DimID); message='iCreate[scalar]';   call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(       wLength_DimName), maxSpectral,       wLength_DimID); message='iCreate[spectral]'; call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(       routing_DimName), maxRouting,        routing_DimID); message='iCreate[routing]';  call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(midSnowAndTime_DimName), maxLength,  midSnowAndTime_DimID); message='iCreate[midSnow]';  call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(midSoilAndTime_DimName), maxLength,  midSoilAndTime_DimID); message='iCreate[midSoil]';  call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(midTotoAndTime_DimName), maxLength,  midTotoAndTime_DimID); message='iCreate[midToto]';  call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(ifcSnowAndTime_DimName), maxLength,  ifcSnowAndTime_DimID); message='iCreate[ifcSnow]';  call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(ifcSoilAndTime_DimName), maxLength,  ifcSoilAndTime_DimID); message='iCreate[ifcSoil]';  call netcdf_err(err,message); if (err/=0) return
 err = nf90_def_dim(ncid, trim(ifcTotoAndTime_DimName), maxLength,  ifcTotoAndTime_DimID); message='iCreate[ifcToto]';  call netcdf_err(err,message); if (err/=0) return

 ! Leave define mode of NetCDF files
 err = nf90_enddef(ncid);  message='nf90_enddef'; call netcdf_err(err,message); if (err/=0) return
 
 end subroutine ini_create

 ! **********************************************************************************************************
 ! private subroutine put_attrib: put global attributes as character string
 ! **********************************************************************************************************
 subroutine put_attrib(ncid,attname,attvalue,err,message)
 USE data_types,only:var_info              ! derived type for metaData
 implicit none
 ! declare dummy variables
 integer(i4b), intent(in)   :: ncid        ! netcdf file ID
 character(*), intent(in)   :: attname     ! attribute name
 character(*), intent(in)   :: attvalue    ! attribute vaue
 integer(i4b),intent(out)   :: err         ! error code
 character(*),intent(out)   :: message     ! error message
 ! initialize error control
 err=0;message="put_attrib/"//trim(attname)//"/"//trim(attvalue)//"/"
 ! allow re-definition of variables
 err = nf90_redef(ncid); call netcdf_err(err,message); if (err/=0) return
 ! put the attribute
 err = nf90_put_att(ncid,nf90_global,trim(attname),trim(attvalue))
 call netcdf_err(err,message); if (err/=0) return
 ! close output file
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine put_attrib

 ! **********************************************************************************************************
 ! private subroutine def_variab: define variables
 ! **********************************************************************************************************
 subroutine def_variab(ncid,iFreq,hruDesire,timeDesire,metaData,ivtype,err,message)
 USE var_lookup,only:iLookvarType                   ! look up structure for variable typed
 USE data_types,only:var_info                       ! derived type for metaData
 USE var_lookup,only:iLookStat                      ! index into stats structure
 USE var_lookup,only:maxVarStat                     ! # of available stats
 USE get_ixName_module,only:get_varTypeName         ! to access type strings for error messages
 USE get_ixname_module,only:get_statName            ! statistics names for variable defs in output file
 implicit none
 ! input
 integer(i4b)  ,intent(in)     :: ncid              ! netcdf file id
 integer(i4b)  ,intent(in)     :: iFreq             ! frequency of current file
 integer(i4b)  ,intent(in)     :: hruDesire         ! variable to define if we desire the HRU dimension
 integer(i4b)  ,intent(in)     :: timeDesire        ! variable to define if we desire the time dimension
 type(var_info),intent(inout)  :: metaData(:)       ! metaData structure for a given variable
 integer(i4b)  ,intent(in)     :: ivtype            ! variable type
 ! output
 integer(i4b),intent(out)      :: err               ! error code
 character(*),intent(out)      :: message           ! error message
 ! local
 integer(i4b)                  :: iVar              ! variable index
 integer(i4b)                  :: iStat             ! stat index
 integer(i4b),allocatable      :: dimensionIDs(:)   ! vector of dimension IDs
 integer(i4b)                  :: iVarId            ! variable ID
! integer                       :: index             ! intrinsic function to find substring index
 integer(i4b)                  :: timePosition      ! extrinsic variable to hold substring index
 character(LEN=256)            :: cmessage          ! error message of downwind routine
 character(LEN=256)            :: catName           ! full variable name
 ! initialize error control
 err=0; message='def_variab/'

 ! allow re-definition of variables
 err = nf90_redef(ncid); call netcdf_err(err,message); if (err/=0) return

 ! loop through metaData
 do iVar = 1,size(metaData)

  ! check that the variable is desired
  if (metaData(iVar)%varType==iLookvarType%unknown) cycle
  if ((iFreq.ne.metaData(iVar)%outFreq).and.(metaData(iVar)%varName.ne.'time')) cycle

  ! special case of the time variable
  if(metaData(iVar)%varName == 'time')then
   call cloneStruc(dimensionIDs, lowerBound=1, source=(/Timestep_DimID/),err=err,message=cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage)//' [variable '//trim(metaData(iVar)%varName)//']'; return; end if

  ! standard case
  else
   select case(metaData(iVar)%varType)
    ! (scalar variable -- many different types)
    case(iLookvarType%scalarv)
     if(hruDesire==needHRU .and. timeDesire==needTime) call cloneStruc(dimensionIDs, lowerBound=1, source=(/     hru_DimID,Timestep_DimID/), err=err, message=cmessage)
     if(hruDesire==needHRU .and. timeDesire==  noTime) call cloneStruc(dimensionIDs, lowerBound=1, source=(/     hru_DimID/)               , err=err, message=cmessage)
     if(hruDesire==  noHRU .and. timeDesire==needTime) call cloneStruc(dimensionIDs, lowerBound=1, source=(/Timestep_DimID/)               , err=err, message=cmessage)
     if(hruDesire==  noHRU .and. timeDesire==  noTime) call cloneStruc(dimensionIDs, lowerBound=1, source=(/  scalar_DimID/)               , err=err, message=cmessage)

    ! (other variables)
    case(iLookvarType%wLength); call cloneStruc(dimensionIDs, lowerBound=1, source=(/hru_DimID, wLength_DimID,       Timestep_DimID/), err=err, message=cmessage)
    case(iLookvarType%midSnow); call cloneStruc(dimensionIDs, lowerBound=1, source=(/hru_DimID, midSnowAndTime_DimID               /), err=err, message=cmessage)
    case(iLookvarType%midSoil); call cloneStruc(dimensionIDs, lowerBound=1, source=(/hru_DimID, midSoilAndTime_DimID               /), err=err, message=cmessage)
    case(iLookvarType%midToto); call cloneStruc(dimensionIDs, lowerBound=1, source=(/hru_DimID, midTotoAndTime_DimID               /), err=err, message=cmessage)
    case(iLookvarType%ifcSnow); call cloneStruc(dimensionIDs, lowerBound=1, source=(/hru_DimID, ifcSnowAndTime_DimID               /), err=err, message=cmessage)
    case(iLookvarType%ifcSoil); call cloneStruc(dimensionIDs, lowerBound=1, source=(/hru_DimID, ifcSoilAndTime_DimID               /), err=err, message=cmessage)
    case(iLookvarType%ifcToto); call cloneStruc(dimensionIDs, lowerBound=1, source=(/hru_DimID, ifcTotoAndTime_DimID               /), err=err, message=cmessage)
    case(iLookvarType%parSoil); call cloneStruc(dimensionIDs, lowerBound=1, source=(/hru_DimID, depth_DimID                        /), err=err, message=cmessage)
    case(iLookvarType%routing); call cloneStruc(dimensionIDs, lowerBound=1, source=(/routing_DimID                                 /), err=err, message=cmessage)
   end select
   ! check errors
   if(err/=0)then
    message=trim(message)//trim(cmessage)//' [variable '//trim(metaData(iVar)%varName)//']'
    return
   end if
  end if  ! check if we are processing the time variable

  ! check that we got the shape
  if(.not.allocated(dimensionIDs))then
   message=trim(message)//'problem defining dimensions for variable '//trim(metaData(iVar)%varName)
   err=20; return
  end if

  ! loop through statistics
  do iStat = 1,maxvarStat

   ! if requested
   if ((.not.metaData(iVar)%statFlag(iStat)).and.(metaData(iVar)%varName.ne.'time'))  cycle
   if ((metaData(iVar)%varName=='time').and.(iStat.ne.iLookStat%inst)) cycle

   ! create full variable name
   catName = trim(metaData(iVar)%varName)
   if (iStat.ne.iLookStat%inst) catName = trim(metaData(iVar)%varName)//'_'//trim(get_statName(iStat))

   ! define variable
   err = nf90_def_var(ncid,trim(catName),ivtype,dimensionIDs,iVarId)
   call netcdf_err(err,message); if (err/=0) return

   ! add parameter description
   catName = trim(metaData(iVar)%vardesc)//' ('//trim(get_statName(iStat))
   catName = trim(catName)//')' 
   err = nf90_put_att(ncid,iVarId,'long_name',trim(catName))
   call netcdf_err(err,message); if (err/=0) return

   ! add parameter units
   catName = trim(metaData(iVar)%varunit) 
   if (iStat==iLookStat%totl) then
    ! make sure that the units of this varaible allow for integration
    if ((index(catName,'s-1')<=0).and.(index(catName,'s-2')<=0).and.(index(catName,'W m-2')<=0)) then 
     err=20
     message=trim(message)//'trying to integrate a non-time variable: '//trim(metaData(iVar)%varName)//' - units: '//trim(catName)
     return
    endif
    ! change to integrated units
    if (index(catName,'s-1')>0)       then 
     timePosition = index(catName,'s-1')
     catName(timePosition:(timePosition+3)) = '   '
    elseif (index(catName,'s-2')>0)   then 
     timePosition = index(catName,'s-2')
     catName(timePosition:(timePosition+3)) = 's-1'
    elseif (index(catName,'W m-2')>0) then 
     timePosition = index(catName,'W') 
     catName(timePosition:(timePosition+1)) = 'J'
    end if
   end if
   err = nf90_put_att(ncid,iVarId,'units',trim(catName))
   call netcdf_err(err,message); if (err/=0) return

   ! add file info to metadata structure
   metaData(iVar)%ncVarID(iStat) = iVarID

  end do ! looping through statistics
 end do  ! looping through variables
  
 ! close output file
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return

 end subroutine def_variab

 ! **********************************************************************************************************
 ! internal subroutine write_hru_dim: write HRU dimension
 ! **********************************************************************************************************
 subroutine write_hru_dim(ncid, err, message)  
 use globalData,only:gru_struc                    ! gru-hru mapping structures  
 ! input
 integer(i4b)  ,intent(in)   :: ncid              ! netcdf file id
 ! output
 integer(i4b),intent(out)    :: err               ! error code
 character(*),intent(out)    :: message           ! error message
 ! define local variables
 integer(i4b)                :: iHRU              ! local HRU index 
 integer(i4b)                :: iGRU              ! GRU index 
 integer(i4b)                :: hruVarID          ! HRU varID in netcdf 
 
 ! initialize error control
 err=0; message='write_hru_dim/'
 
 ! allow re-definition of variables
 err = nf90_redef(ncid); call netcdf_err(err, message); if (err/=nf90_NoErr) return
 
 ! define HRU var
 err = nf90_def_var(ncid, trim(hru_DimName), nf90_int, hru_DimID, hruVarID);     if (err/=nf90_NoErr) then; message=trim(message)//'nf90_define_hruVar'  ;  call netcdf_err(err,message); return; end if 
 err = nf90_put_att(ncid, hruVarID, 'long_name', 'hru index in the input file'); if (err/=nf90_NoErr) then; message=trim(message)//'write_hruVar_longname'; call netcdf_err(err,message); return; end if 
 err = nf90_put_att(ncid, hruVarID, 'units',     '-'                          ); if (err/=nf90_NoErr) then; message=trim(message)//'write_hruVar_unit';     call netcdf_err(err,message); return; end if 
 
 ! Leave define mode of NetCDF files
 err = nf90_enddef(ncid);  message=trim(message)//'nf90_enddef'; call netcdf_err(err,message); if (err/=nf90_NoErr) return
 
 ! write the HRU dimension to record position in the input netcdf file for concatenation of outputs of a parallelized run.
 do iGRU = 1, size(gru_struc)
  do iHRU = 1, gru_struc(iGRU)%hruCount
   err = nf90_put_var(ncid, hruVarID, gru_struc(iGRU)%hruInfo(iHRU)%hru_nc, start=(/gru_struc(iGRU)%hruInfo(iHRU)%hru_ix/))
   if (err/=nf90_NoErr) then; message=trim(message)//'nf90_write_hruVar'; call netcdf_err(err,message); return; end if
  end do
 end do 
 
 end subroutine
 
end module def_output_module
