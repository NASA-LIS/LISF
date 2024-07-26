!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: PMW_snow_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to handle PMW-based 
!   swe or snow depth retrievals (e.g., AMSR-E, ANSA, SSMI, SMMR).
!   Currently supported data products include:
!
!   1) LATLON grid in HDF5 file format (e.g., ANSA, the original dataset 
!       or bias corrected dataset; see references below) \newline
!   2) EASE grid (NL \& SL) in HDF-EOS file format (e.g., AMSR-E Leve 3 SWE) \newline
!   3) EASE grid (NL \& SL) in HDF4 file format (e.g., SMMR and SSMI)
!  
!   NOTE: As of 4/23/2014, only option 1) is rigorously tested and confirmed 
!   to work as expected. Options 2) and 3) compile successfully but have not 
!   been tested.
! 
!   Note the data file names must have one of the following two conventions: \newline
!   - *YYYYMMDD* \newline
!   - *YYYYDOY*  \newline
!   In addition, data files should be stored in corresponding year directory
!   as follows: datadir/YYYY/YYYMMDD
!
!   ANSA snow depth data is primarily based on PMW, with observation
!   gaps filled in using previous overpasses. For more details, please see 
!   the ANSA reference: 
!
!   Foster et al. ``A blended global snow product using visible, 
!   passive microwave and scatterometer satellite data'', International
!   Journal of Remote Sensing, 2010. 
!  
!  For PMW snow depth bias correction, please refer to:
!  
!  Liu, Y., C.D. Peters-Lidard, S. Kumar, J. Foster, M. Shaw, Y. Tian,
!  and G. Fall, 2013: 
!  Assimilating satellite-based snow depth and snow cover products for
!  improving snow predictions in Alaska. Advances in Water Resources,
!  54, 208-227.
!
! !REVISION HISTORY: 
!  1 Jun 09   Sujay Kumar;  Initial Specification
!  12 Apr 2012 Yuqiong Liu; adapted to assimilate PMW-based snow depth or
!                           SWE data

 
module PMW_snow_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: PMW_snow_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: PMW_snow_struc

  type, public ::  PMW_snow_dec
     logical             :: startMode
     real                :: ssdev
     integer             :: mi, mo
     real                :: cornerlat1,cornerlat2
     real                :: cornerlon1,cornerlon2
     real                :: minlat,minlon
     real                :: maxlat,maxlon
     integer             :: nc,nr
     integer             :: offset1, offset2 !for HDF5 data subsetting
     real                :: gridDesc(6)
     real,    allocatable    :: rlat2(:,:),rlon2(:,:)
     integer, allocatable    :: n11(:), n12(:), n21(:), n22(:), n112(:,:), n113(:)
     real,    allocatable    :: w11(:), w12(:), w21(:), w22(:)
     real,    allocatable    :: snow(:) !PMW snow data
     character*100       :: snow_field_name, snow_flag_field_name !HDF5 snow data and flag fields
     integer             :: data_nl_sdsid, data_sl_sdsid, flag_nl_sdsid, flag_sl_sdsid !HDF4 snow data and flag fields
     character*100       :: data_nl_grid, data_sl_grid, data_nl_sds, data_sl_sds, flag_nl_sds, flag_sl_sds !HDF-ESO
     integer             :: flag_n1_invalid_value, data_n1_invalid_value !invlaid values in flag field
     real, allocatable       :: flag_invalid_value(:), data_invalid_value(:) !additional invalid snow data values
     real                :: data_min, data_max, data_scale !snow data min, max amd scale factor
     integer             :: gvf_mask, lctype_mask, temp_mask, use_flag, use_minmax !whether to use these masks/flags
     character*50        :: dataset, data_var, dataunit, data_format, data_coordsys !infomration on the dataset
     character*100       :: data_fn_conv  !data file convention
     real                :: assim_lhour  !local time to conduct assimilation
     integer             :: ihemi, nhemi !number of hemisphere
     

  end type PMW_snow_dec

  type(PMW_snow_dec), allocatable :: PMW_snow_struc(:)

contains
!BOP
! 
! !ROUTINE: PMW_snow_setup
! \label{PMW_snow_setup}
! 
! !INTERFACE: 
  subroutine PMW_snow_setup(k, OBS_State, OBS_Pert_State)
! !USES: 

    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
    use LIS_logMod
   
    implicit none 

! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for PMW snow data. 
!
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n
    integer                ::  ftn
    integer                ::  i
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  PMW_snowobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real,         allocatable  ::  varmin(:)
    real,         allocatable  ::  varmax(:)
    type(pert_dec_type)    :: obs_pert
    real, pointer          :: obs_temp(:,:)
    real, allocatable          :: ssdev(:)
    real                   :: dx, dy

    allocate(PMW_snow_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snowobsdir,&
            rc=status)
       call LIS_verify(status,'PMW snow data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            PMW_snowobsdir, rc=status)
       call LIS_verify(status)
    enddo

    do n=1,LIS_rc%nnest
       call ESMF_AttributeSet(OBS_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Update Time",&
            -99.0, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Assimilate Status",&
            .false., rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeSet(OBS_State(n),"Number Of Observations",&
            LIS_rc%ngrid(n),rc=status)
       call LIS_verify(status)
       
    enddo

    write(LIS_logunit,*) 'read PMW snow data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. For this synthetic case, it is assumed that the 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%ngrid(n). 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(grid=LIS_obsVecGrid(n,k),&
            arrayspec=realarrspec,&
            name="Observation"//vid(1)//vid(2), rc=status)
       call LIS_verify(status)
       
!Perturbations State
       write(LIS_logunit,*) 'Opening attributes for observations ',&
            trim(LIS_rc%obsattribfile(k))
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(LIS_rc%obsattribfile(k)),status='old')
       read(ftn,*)
       read(ftn,*) LIS_rc%nobtypes(k)
       read(ftn,*)
    
       allocate(vname(LIS_rc%nobtypes(k)))
       allocate(varmax(LIS_rc%nobtypes(k)))
       allocate(varmin(LIS_rc%nobtypes(k)))
       
       do i=1,LIS_rc%nobtypes(k)
          read(ftn,fmt='(a40)') vname(i)
          read(ftn,*) varmin(i),varmax(i)
          write(LIS_logunit,*) vname(i),varmin(i),varmax(i)
       enddo
       call LIS_releaseUnitNumber(ftn)  
       
       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

       allocate(ssdev(LIS_rc%ngrid(n)))

       if(trim(LIS_rc%perturb_obs(k)).ne."none") then 

          allocate(obs_pert%vname(1))
          allocate(obs_pert%perttype(1))
          allocate(obs_pert%ssdev(1))
          allocate(obs_pert%stdmax(1))
          allocate(obs_pert%zeromean(1))
          allocate(obs_pert%tcorr(1))
          allocate(obs_pert%xcorr(1))
          allocate(obs_pert%ycorr(1))
          allocate(obs_pert%ccorr(1,1))

          call LIS_readPertAttributes(1,LIS_rc%obspertAttribfile(k),&
               obs_pert)

          ssdev = obs_pert%ssdev(1)
          PMW_snow_struc(n)%ssdev =obs_pert%ssdev(1) 

          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_obsEnsOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status)
          
! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
          call LIS_verify(status)
          obs_temp(:,:) = 0 

          call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
               obs_pert%perttype(1), rc=status)
          call LIS_verify(status)
          
          if(LIS_rc%ngrid(n).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
             call LIS_verify(status)
          endif

          call ESMF_AttributeSet(pertField(n),"Std Normal Max",&
               obs_pert%stdmax(1), rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean",&
               obs_pert%zeromean(1),rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Temporal Correlation Scale",&
               obs_pert%tcorr(1),rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"X Correlation Scale",&
               obs_pert%xcorr(1),rc=status)
          
          call ESMF_AttributeSet(pertField(n),"Y Correlation Scale",&
               obs_pert%ycorr(1),rc=status)

          call ESMF_AttributeSet(pertField(n),"Cross Correlation Strength",&
               obs_pert%ccorr(1,:),itemCount=1,rc=status)

          call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
          call LIS_verify(status)
       endif
          
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)
    enddo

! read PMW snow data configuration 
    call readPMWSnowDataConfig()
    

! compute interpolation wights
    do n=1, LIS_rc%nnest

       PMW_snow_struc(n)%mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

       if (PMW_snow_struc(n)%data_coordsys .eq. "LATLON") then          
          call computeInterpWeights_latlon()
       else if (PMW_snow_struc(n)%data_coordsys .eq. "EASE") then
          call computeInterpWeights_ease()
       end if

       allocate(PMW_snow_struc(n)%snow(PMW_snow_struc(n)%mo))
       PMW_snow_struc(n)%snow = LIS_rc%udef
       
       call LIS_registerAlarm("PMW snow data read alarm",&
            86400.0, 86400.0)

       PMW_snow_struc(n)%startMode = .true. 
    enddo

    write(LIS_logunit,*) 'Created ESMF States to hold PMW observations data'

  end subroutine PMW_snow_setup


!BOP
! !ROUTINE: readPMWSnowDataConfig
! \label{PMW_snow_readPMWSnowDataConfig}
! 
! !INTERFACE: 
  subroutine readPMWSnowDataConfig()
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config
    use LIS_logMod, only : LIS_logunit, LIS_verify

    implicit none

    integer       :: n, i
    integer       :: status

! 
! !DESCRIPTION: 
!   This subroutine reads the configuration for PMW snow assimilation 
!   from the lis.config file 
! 
!EOP
!-------------------------------------------------------------
! read in dataset parameters
!-------------------------------------------------------------
    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data file format (HDF4, HDF-EOS, HDF5):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_format,&
            rc=status)
       call LIS_verify(status,'PMW snow data file format (HDF4, HDF-EOS, HDF5): not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data coordinate system (EASE, LATLON):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_coordsys,&
            rc=status)
       call LIS_verify(status,'PMW snow data coordinate system (EASE, LATLON): not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data use flag (1=yes, 0=no):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%use_flag,&
            rc=status)
       call LIS_verify(status,'PMW snow data use flag (1=yes, 0=no): not defined')
    enddo
      
    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data variable (SWE, snow depth):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_var,&
            rc=status)
       call LIS_verify(status,'PMW snow data variable (SWE, snow depth): not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data unit (m, cm, mm, inch):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%dataunit,&
            rc=status)
       call LIS_verify(status,'PMW snow data unit (m, cm, mm, inch): not defined')
    enddo

!!! HDF5, latitude-longitude coordinate system (e.g., ANSA) 
!!! currently all nests must have the same PMW snow data format and coordinate system
    if (trim(PMW_snow_struc(1)%data_format) .eq. "HDF5" .and. trim(PMW_snow_struc(1)%data_coordsys) .eq. "LATLON") then

       call ESMF_ConfigFindLabel(LIS_config,"PMW snow data lower left lat:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%gridDesc(1),&
               rc=status)
          call LIS_verify(status,'PMW snow data lower left lat: not defined')
       enddo
   
       call ESMF_ConfigFindLabel(LIS_config,"PMW snow data lower left lon:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%gridDesc(2),&
               rc=status)
          call LIS_verify(status,'PMW snow data lower left lon: not defined')
       enddo
       
       call ESMF_ConfigFindLabel(LIS_config,"PMW snow data upper right lat:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%gridDesc(3),&
               rc=status)
          call LIS_verify(status,'PMW snow data upper right lat: not defined')
       enddo
   
       call ESMF_ConfigFindLabel(LIS_config,"PMW snow data upper right lon:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%gridDesc(4),&
               rc=status)
          call LIS_verify(status,'PMW snow data upper right lon: not defined')
       enddo
   
       call ESMF_ConfigFindLabel(LIS_config,"PMW snow data resolution (dx):",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%gridDesc(5),&
               rc=status)
          call LIS_verify(status,'PMW or ANSA snow data resolution (dx): not defined')
       enddo
   
       call ESMF_ConfigFindLabel(LIS_config,"PMW snow data resolution (dy):",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%gridDesc(6),&
               rc=status)
          call LIS_verify(status,'PMW snow data resolution (dy): not defined')
       enddo
   
       call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF5) snow data field name:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%snow_field_name,&
               rc=status)
          call LIS_verify(status,'PMW (HDF5) snow data field name: not defined')
       enddo
   
       if (PMW_snow_struc(1)%use_flag .eq. 1) then
          call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF5) snow data flag field name:",&
               rc=status)
          do n=1,LIS_rc%nnest
             call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%snow_flag_field_name,&
                  rc=status)
             call LIS_verify(status,'PMW (HDF5) snow data flag field name: not defined')
          enddo
       end if

!!! HDF4, EASE coordinate system (e.g., SMMR, SSMI)
    else if (trim(PMW_snow_struc(1)%data_format) .eq. "HDF4" .and. trim(PMW_snow_struc(1)%data_coordsys) .eq. "EASE") then
 
       call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF4) snow data NL SDS index (-1, 0, 1, 2, ...):",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_nl_sdsid,&
               rc=status)
          call LIS_verify(status,'PMW (HDF4) snow data NL SDS index (-1, 0, 1, 2, ...): not defined')
       enddo

       call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF4) snow data SL SDS index (-1, 0, 1, 2, ...):",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_sl_sdsid,&
               rc=status)
          call LIS_verify(status,'PMW (HDF4) snow data SL SDS index (-1, 0, 1, 2, ...): not defined')
       enddo
   
       if (PMW_snow_struc(1)%use_flag .eq. 1) then
   
          if (PMW_snow_struc(1)%data_nl_sdsid .ge. 0) then
             call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF4) snow data flag NL SDS index (-1, 0, 1, 2, ...):",&
                 rc=status)
             do n=1,LIS_rc%nnest
                call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%flag_nl_sdsid,&
                    rc=status)
                call LIS_verify(status,'PMW (HDF4) snow data flag NL SDS index (-1, 0, 1, 2, ...): not defined')
             enddo
          end if
   
          if (PMW_snow_struc(1)%data_sl_sdsid .ge. 0) then
             call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF4) snow data flag SL SDS index (-1, 0, 1, 2, ...):",&
                  rc=status)
             do n=1,LIS_rc%nnest
                call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%flag_sl_sdsid,&
                     rc=status)
                call LIS_verify(status,'PMW (HDF4) snow data flag SL SDS index (-1, 0, 1, 2, ...): not defined')
             enddo
          end if
   
       end if
   
!!! HDF-EOS, EASE coordinate system (e.g., AMSR-E)
    else if (trim(PMW_snow_struc(1)%data_format) .eq. "HDF-EOS" .and. trim(PMW_snow_struc(1)%data_coordsys) .eq. "EASE") then

       call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF-EOS) NL grid name:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_nl_grid,&
               rc=status)
          call LIS_verify(status,'PMW (HDF-EOS) NL grid name: not defined')
       enddo

       call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF-EOS) SL grid name:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_sl_grid,&
               rc=status)
          call LIS_verify(status,'PMW (HDF-EOS) SL grid name: not defined')
       enddo
   
       if (PMW_snow_struc(1)%data_nl_grid .ne. "none") then
          call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF-EOS) NL SDS name:",&
               rc=status)
          do n=1,LIS_rc%nnest
             call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_nl_sds,&
                  rc=status)
             call LIS_verify(status,'PMW (HDF-EOS) NL SDS name: not defined')
          enddo
       end if 
   
       if (PMW_snow_struc(1)%data_sl_grid .ne. "none") then
          call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF-EOS) SL SDS name:",&
               rc=status)
          do n=1,LIS_rc%nnest
             call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_sl_sds,&
                  rc=status)
             call LIS_verify(status,'PWm (HDF-EOS) SL SDS name: not defined')
          enddo
       end if 
   
       if (PMW_snow_struc(1)%use_flag .eq. 1) then
          if (PMW_snow_struc(1)%data_nl_sds .ne. "none") then
             call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF-EOS) NL snow data flag SDS name:",&
                  rc=status)
             do n=1,LIS_rc%nnest
                call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%flag_nl_sds,&
                     rc=status)
                call LIS_verify(status,'PMW (HDF-EOS) NL snow data flag SDS name: not defined')
             enddo
          end if 
   
          if (PMW_snow_struc(1)%data_sl_sds .ne. "none") then
             call ESMF_ConfigFindLabel(LIS_config,"PMW (HDF-EOS) SL snow data flag SDS name:",&
                  rc=status)
             do n=1,LIS_rc%nnest
                call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%flag_sl_sds,&
                     rc=status)
                call LIS_verify(status,'PMW (HDF-EOS) SL snow data flag SDS name: not defined')
             enddo
          end if 
   
       end if 
!!! otherwise, un-supported dataset
    else
       write(LIS_logunit,*) 'Dataset in format ', PMW_snow_struc(1)%data_format, ' and ', &
           PMW_snow_struc(1)%data_coordsys, ' is currently not supported'
       return
    end if  
!!! end of checking dataset

    if (PMW_snow_struc(1)%use_flag .eq. 1) then
       call ESMF_ConfigFindLabel(LIS_config,"PMW snow data flag - number of invalid values:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%flag_n1_invalid_value,&
               rc=status)
          call LIS_verify(status,'PMW snow data flag - number of invalid values: not defined')
       enddo

       if (PMW_snow_struc(1)%flag_n1_invalid_value .gt. 0) then          
          call ESMF_ConfigFindLabel(LIS_config,"PMW snow data flag - invalid values:",rc=status)
          do n=1,LIS_rc%nnest
             allocate(PMW_snow_struc(n)%flag_invalid_value(PMW_snow_struc(n)%flag_n1_invalid_value))
             do i=1,PMW_snow_struc(n)%flag_n1_invalid_value
                call ESMF_ConfigGetAttribute(LIS_config,PMW_snow_struc(n)%flag_invalid_value(i),rc=status)
                call LIS_verify(status,"PMW snow data flag - invalid values: not defined")
             enddo 
          enddo  
       end if 
    end if

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data - number of additional invalid values:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_n1_invalid_value,&
            rc=status)
       call LIS_verify(status,'PMW snow data - number of additional invalid values: not defined')
    enddo

    if (PMW_snow_struc(1)%data_n1_invalid_value .gt. 0) then       
       call ESMF_ConfigFindLabel(LIS_config,"PMW snow data - additional invalid values:",rc=status)
       do n=1,LIS_rc%nnest
          allocate(PMW_snow_struc(n)%data_invalid_value(PMW_snow_struc(n)%data_n1_invalid_value))
          do i=1,PMW_snow_struc(n)%data_n1_invalid_value
             call ESMF_ConfigGetAttribute(LIS_config,PMW_snow_struc(n)%data_invalid_value(i),rc=status)
             call LIS_verify(status, "PMW snow data - additional invalid values: not defined")
          enddo
       enddo   
    end if  

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data - apply min/max mask:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%use_minmax,&
            rc=status)
       call LIS_verify(status,'PMW snow data - apply min/max mask: not defined')
    enddo


    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data minimum valid value:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_min,&
            rc=status)
       call LIS_verify(status,'PMW snow data minimum valid value: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data maximum valid value:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_max,&
            rc=status)
       call LIS_verify(status,'PMW snow data maximum valid value: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data scale factor:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_scale,&
            rc=status)
       call LIS_verify(status,'PMW snow data scale factor: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data file name convention:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%data_fn_conv,&
            rc=status)
       call LIS_verify(status,'PMW snow data file name convention: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data assimilation local time:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%assim_lhour,&
            rc=status)
       call LIS_verify(status,'PMW snow data assimilation local time: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data - apply mask with GVF (1=yes, 0=no):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%gvf_mask,&
            rc=status)
       call LIS_verify(status,'PMW snow data - apply mask with GVF (1=yes, 0=no): not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data - apply mask with landcover type (1=yes, 0=no):",&
         rc=status)
    do n=1,LIS_rc%nnest       
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%lctype_mask,&
            rc=status)
       call LIS_verify(status,'PMW snow data - apply mask with landcover type (1=yes, 0=no): not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"PMW snow data - apply mask with LSM temperature (1=yes, 0=no):",&
         rc=status)
    do n=1,LIS_rc%nnest       
       call ESMF_ConfigGetattribute(LIS_config,PMW_snow_struc(n)%temp_mask,&
            rc=status)
       call LIS_verify(status,'PMW snow data - apply mask with LSM temperature (1=yes, 0=no): not defined')
    enddo

  end subroutine readPMWSnowDataConfig


!BOP
! !ROUTINE: computeInterpWeights_latlon
! \label{PMW_snow_computeInterpWeights_latlon}
! 
! !INTERFACE: 
  subroutine computeInterpWeights_latlon()
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_domain
    use LIS_logMod, only : LIS_logunit

    implicit none

    integer    :: n, i
    real       :: dx, dy
! 
! !DESCRIPTION: 
!   This subroutine sets up the interpolation weights to transform the 
!   PMW snow data in latlon grid to the grid used in the LIS simulation. 
!   The code employs a bilinear-based interpolation
! 
!EOP  

    real                   :: gridDesci(LIS_rc%nnest,50)

    do n=1,LIS_rc%nnest

       PMW_snow_struc(n)%minlat = PMW_snow_struc(n)%gridDesc(1)
       PMW_snow_struc(n)%minlon = PMW_snow_struc(n)%gridDesc(2)
       PMW_snow_struc(n)%maxlat = PMW_snow_struc(n)%gridDesc(3)
       PMW_snow_struc(n)%maxlon = PMW_snow_struc(n)%gridDesc(4)
       dx = PMW_snow_struc(n)%gridDesc(5)
       dy = PMW_snow_struc(n)%gridDesc(6)
    enddo
    
    do n=1, LIS_rc%nnest
!sets the local domain corner points with additional buffer 
       PMW_snow_struc(n)%cornerlat1 = max(PMW_snow_struc(n)%minlat, &
            nint((LIS_domain(n)%minlat-PMW_snow_struc(n)%minlat)/dx)*dx+PMW_snow_struc(n)%minlat-2*dx)
       PMW_snow_struc(n)%cornerlon1 = max(PMW_snow_struc(n)%minlon, &
            nint((LIS_domain(n)%minlon-PMW_snow_struc(n)%minlon)/dy)*dy+PMW_snow_struc(n)%minlon-2*dy)
       PMW_snow_struc(n)%cornerlat2 = min(PMW_snow_struc(n)%maxlat, &
            nint((LIS_domain(n)%maxlat-PMW_snow_struc(n)%minlat)/dx)*dx+PMW_snow_struc(n)%minlat+2*dx)
       PMW_snow_struc(n)%cornerlon2 = min(PMW_snow_struc(n)%maxlon, &
            nint((LIS_domain(n)%maxlon-PMW_snow_struc(n)%minlon)/dy)*dy+PMW_snow_struc(n)%minlon+2*dy)


       PMW_snow_struc(n)%offset1 = &
            nint((PMW_snow_struc(n)%cornerlon1-&
            PMW_snow_struc(n)%minlon)/dy)
       PMW_snow_struc(n)%offset2 = &
            nint((PMW_snow_struc(n)%cornerlat1-&
            PMW_snow_struc(n)%minlat)/dx)

       PMW_snow_struc(n)%nr = &
            nint((PMW_snow_struc(n)%cornerlat2-&
            PMW_snow_struc(n)%cornerlat1)/dx)+1
       PMW_snow_struc(n)%nc = &
            nint((PMW_snow_struc(n)%cornerlon2-&
            PMW_snow_struc(n)%cornerlon1)/dy)+1

       gridDesci(n,1) = 0 
       gridDesci(n,2) = PMW_snow_struc(n)%nc
       gridDesci(n,3) = PMW_snow_struc(n)%nr
       gridDesci(n,4) = PMW_snow_struc(n)%cornerlat1
       gridDesci(n,5) = PMW_snow_struc(n)%cornerlon1
       gridDesci(n,6) = 128
       gridDesci(n,7) = PMW_snow_struc(n)%cornerlat2
       gridDesci(n,8) = PMW_snow_struc(n)%cornerlon2
       gridDesci(n,9) = PMW_snow_struc(n)%gridDesc(5)
       gridDesci(n,10) = PMW_snow_struc(n)%gridDesc(6)
       gridDesci(n,20) = 64.0

       PMW_snow_struc(n)%mi = PMW_snow_struc(n)%nc*PMW_snow_struc(n)%nr

       allocate(PMW_snow_struc(n)%n11(PMW_snow_struc(n)%mo))
       allocate(PMW_snow_struc(n)%n12(PMW_snow_struc(n)%mo))
       allocate(PMW_snow_struc(n)%n21(PMW_snow_struc(n)%mo))
       allocate(PMW_snow_struc(n)%n22(PMW_snow_struc(n)%mo))

       allocate(PMW_snow_struc(n)%w11(PMW_snow_struc(n)%mo))
       allocate(PMW_snow_struc(n)%w12(PMW_snow_struc(n)%mo))
       allocate(PMW_snow_struc(n)%w21(PMW_snow_struc(n)%mo))
       allocate(PMW_snow_struc(n)%w22(PMW_snow_struc(n)%mo))

       call bilinear_interp_input(n, gridDesci(n,:),&
            PMW_snow_struc(n)%n11,PMW_snow_struc(n)%n12,&
            PMW_snow_struc(n)%n21,PMW_snow_struc(n)%n22,&
            PMW_snow_struc(n)%w11,PMW_snow_struc(n)%w12,&
            PMW_snow_struc(n)%w21,PMW_snow_struc(n)%w22)

       allocate(PMW_snow_struc(n)%n113(PMW_snow_struc(n)%mo))

#if 0 
       call neighbor_interp_input(n, gridDesci(n,:),&
            PMW_snow_struc(n)%n113)

#endif
     enddo
   end subroutine computeInterpWeights_latlon


!BOP
! !ROUTINE: computeInterpWeights_ease
! \label{PMW_snow_computeInterpWeights_ease}
! 
! !INTERFACE: 
  subroutine computeInterpWeights_ease()
! !USES: 
    use LIS_coreMod
    use LIS_logMod, only : LIS_logunit
    use map_utils

    implicit none

! 
! !DESCRIPTION: 
!   This subroutine sets up the interpolation weights to transform the 
!   AMSRE data in EASE grid to the grid used in the LIS simulation. 
!   The code employs a neighbor-based interpolation
! 
!EOP
    integer,parameter  :: ease_nr=721
    integer,parameter  :: ease_nc=721
    integer         :: npts,n
    integer         :: hemi
    real            :: gridDesci(50)
    real            :: max_lat, min_lat, rlat, rlon
    integer         :: c,r


!--------------------------------------------------------------------------------
! figure out the grid span and whether we need to process both southern 
! and northern hemispheres
!--------------------------------------------------------------------------------
    do n=1, LIS_rc%nnest
       max_lat = -10000
       min_lat = 10000 
       
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
             call ij_to_latlon(LIS_domain(n)%lisproj,float(c),float(r),&
                  rlat,rlon)
             
             if(rlat.gt.max_lat) max_lat = rlat
             if(rlat.lt.min_lat) min_lat = rlat
          enddo
       enddo

       if(max_lat.gt.0.and.min_lat.lt.0) then ! domain split in 2 hemispheres
          PMW_snow_struc(n)%ihemi = 1
          PMW_snow_struc(n)%nhemi = 2
       elseif(max_lat.ge.0.and.min_lat.ge.0) then !all northern hemisphere
          PMW_snow_struc(n)%ihemi = 1
          PMW_snow_struc(n)%nhemi = 1
       else !all southern hemisphere
          PMW_snow_struc(n)%ihemi = 2
          PMW_snow_struc(n)%nhemi = 2
       endif

    enddo

!------------------------------------------

    !read the domain specs from LIS_rc% struct

    do n =1, LIS_rc%nnest

       npts= LIS_rc%lnc(n)*LIS_rc%lnr(n)
       PMW_snow_struc(n)%mo=npts
       allocate(PMW_snow_struc(n)%rlat2(npts,2))
       allocate(PMW_snow_struc(n)%rlon2(npts,2))
       allocate(PMW_snow_struc(n)%n112(npts,2))

       do hemi=PMW_snow_struc(n)%ihemi, PMW_snow_struc(n)%nhemi 
          !initialize the entire array
          gridDesci =0.0 
          
      !filling the items needed by the interpolation library
          gridDesci(1) = 9  !input is EASE grid
          gridDesci(2) = ease_nc  !nx
          gridDesci(3) = ease_nr  !ny

          !these  corner coordinates were calculated based on ezlh_convert
          gridDesci(4) = -90.0  !lat
          gridDesci(5) = -179.6096 !lon
          gridDesci(7) = 83.33788  !lat
          gridDesci(8) = 180.1301  !lon

          if(hemi.eq.1) then 
             gridDesci(9)  = 2  !Northern hemi
          else
             gridDesci(9)  = 3  !Southern hemi
          endif
          
          PMW_snow_struc(n)%rlat2(:,hemi)=0.0
          PMW_snow_struc(n)%rlon2(:,hemi)=0.0
          PMW_snow_struc(n)%n112(:,hemi)=0.0
          call neighbor_interp_input_withgrid(gridDesci,&
               LIS_rc%gridDesc(n,:),&
               npts,PMW_snow_struc(n)%rlat2(:,hemi),&
               PMW_snow_struc(n)%rlon2(:,hemi),PMW_snow_struc(n)%n112(:,hemi))
       enddo
       
    end do
      
  end subroutine computeInterpWeights_ease

end module PMW_snow_Mod

