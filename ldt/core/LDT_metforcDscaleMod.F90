!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_metforcDscaleMod
!BOP
!
! !MODULE: LDT_metforcDscaleMod
! 
! !DESCRIPTION:
! 
! !REVISION HISTORY: 
!  16 Oct 2014: KR Aresnault; initial specification
!

  use ESMF
  use LDT_logMod
  use LDT_constantsMod
  use LDT_timeMgrMod
  use LDT_coreMod
  use LDT_metforcingMod
  use LDT_FORC_AttributesMod
  use LDT_fileIOMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LDT_metforcDscaleInit
  PUBLIC :: LDT_diagnoseForMetForcDscale
  PUBLIC :: LDT_computeMetForcDscaleParams
  PUBLIC :: LDT_applyDscaleCorrection
  PUBLIC :: LDT_setupMetDscaleTimeWindow
  PUBLIC :: LDT_endofMetDscaleTimeWindow
  PUBLIC :: LDT_resetMetDscaleTimeWindow
  PUBLIC :: LDT_resetDscaleVars

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!----------------------------------------------------------------------------- 
  PUBLIC :: LDT_metDscale
  PUBLIC :: LDT_FORC_Dscale_State
 
  type, public :: metDscale_dec_type
     real    :: fine_ts
     real    :: coarse_ts
     integer :: numfine2coarse
     real    :: last_validhr
     character(20) :: fine_metname
     character(20) :: coarse_metname
  end type metDscale_dec_type

  type(metDscale_dec_type) :: LDT_metDscale
  type(ESMF_Time)               :: LDT_metProcTWstart
  type(ESMF_Time)               :: LDT_metProcTWend
  type(ESMF_State), allocatable :: LDT_FORC_Dscale_State(:)

!EOP

contains

!BOP
! !ROUTINE: LDT_metforcDscaleInit
! \label{LDT_metforcDscaleInit}
!
! !REVISION HISTORY:
!  14 Oct 2014: KR Arsenault; Initial interface for Metforc downscaling
!
! !INTERFACE:
  subroutine LDT_metforcDscaleInit()

! !USES:

   implicit none
! !ARGUMENTS: 

! !DESCRIPTION:
!
!  Initialize the meteorological forcing downscaling 
!   components to either temporally or spatially 
!   interpolate to a finer target grid.
!
!EOP

   integer          :: rc
   integer          :: n, m
   character*10     :: time
   character*100    :: temp
   character*1      :: nestid(2)
   integer          :: status
   integer          :: swyr, swmo, swda, swhr
   integer          :: ewyr, ewmo, ewda, ewhr
! ___________________________________________

   write(LDT_logunit,*)" "
   write(LDT_logunit,*)" - - - - - - Metforcing Temporal Downscaling - - - - - -"
   
   ! TIME WINDOW PARAMETERS - READ IN FROM LDT.CONGIF FILE ...
   call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%timeDscaleType,&
             label="Temporal downscaling method:",&
             default="Simple weighting", rc=rc)
   call LDT_verify(rc,'Temporal downscaling method: not defined')

   write(LDT_logunit,*) "[INFO] Running Metforcing Downscaling Method: ",&
         trim(LDT_rc%timeDscaleType)//" -- "

   allocate(LDT_FORC_Dscale_State(LDT_rc%nnest))

   ! Read in met forcing processing interval (set for time window):
   call ESMF_ConfigGetAttribute(LDT_config,time,&
             label="Metforcing processing interval:",rc=rc)
   call LDT_verify(rc,'Metforcing processing interval: not defined')
   call LDT_parseTimeString(time,LDT_rc%metForcProcInterval)

#if 0
 ! Read in met forcing processing start and ending hours (set for time window):
 ! (KRA)
   LDT_rc%metForcTWstarthour = 0
   LDT_rc%metForcTWendhour = 0

   call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%metForcTWstarthour,&
             label="Metforcing time window start hour:",rc=rc)
   call LDT_verify(rc,'Metforcing time window start hour: not defined')

   call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%metForcTWendhour,&
             label="Metforcing time window ending hour:",rc=rc)
   call LDT_verify(rc,'Metforcing time window ending hour: not defined')
!
!  - Test version that allows user to specify the above start / end window times:
!   call LDT_setupMetDscaleTimeWindow(LDT_rc%metForcProcInterval, &
!            LDT_metProcTWstart, LDT_metProcTWend)
!  (KRA)
#endif

   ! Set initial time window ("TW") start and end times:
   call LDT_setupTimeWindow(LDT_rc%metForcProcInterval, &
            LDT_metProcTWstart, LDT_metProcTWend)

   call ESMF_TimeGet(LDT_metProcTWstart, yy=swyr, mm=swmo, dd=swda, h=swhr)
   call ESMF_TimeGet(LDT_metProcTWend, yy=ewyr, mm=ewmo, dd=ewda, h=ewhr)

   ! Create Downscaled Forcing/State structures for all forcing fields:
   do n=1,LDT_rc%nnest       
      write(unit=temp,fmt='(i2.2)') n
      read(unit=temp,fmt='(2a1)') nestid
      
      LDT_FORC_Dscale_State(n) = ESMF_StateCreate(name=&
           "Forcing State"//nestid(1)//nestid(2),&
           rc=status)
      call LDT_verify(status, &
           'error in ESMF_StateCreate:LDT_FORC_Dscale_State in LDT_metforcDscaleInit')       
      
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Tair)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Qair)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_SWdown)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_SWdirect)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_SWdiffuse)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_LWdown)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Wind_E)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Wind_N)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Psurf)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Rainf)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Snowf)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_CRainf)
      
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Forc_Hgt)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Ch)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Cm)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Q2sat)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Emiss)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Cosz)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Alb)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Pardr)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Pardf)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_SWnet)
      
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Xice)
      
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_QSFC)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_CHS2)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_CQS2)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_T2)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_Q2)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_TH2)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_TMN)
      
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_lpressure)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_forc_o3)
      
      !<for vic>
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
            LDT_FORC_SNOWFLAG)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_FORC_DENSITY)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_FORC_VAPORPRESS)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_FORC_VAPORPRESSDEFICIT)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_FORC_WIND)
      !</for vic>
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_FORC_PET)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_FORC_RefET)
      call add_forcing_fields(n,LDT_FORC_Dscale_State(n),&
           LDT_FORC_CAPE)
      
   enddo

! ---------

! Determine timesteps based on selected metforcings 
!  and downscaling method:

   LDT_metDscale%fine_ts   = LDT_rc%udef
   LDT_metDscale%coarse_ts = LDT_rc%udef
   LDT_metDscale%numfine2coarse = LDT_rc%udef

   select case( LDT_rc%timeDscaleType ) 

    ! Simple weighting option (based on Alo Clement method):
    case( "Simple weighting" )

      if( LDT_rc%nmetforc > 2 .or. LDT_rc%nmetforc < 2 ) then
        write(LDT_logunit,*) "[ERR] Simple weighting only takes 2 forcing types"
        write(LDT_logunit,*) " at this time.  Program stopping ..."
        call LDT_endrun
      endif

      ! Determine fine and coarse timescales of two forcings selected:
      if( LDT_rc%met_ts(1) < LDT_rc%met_ts(2) ) then
         LDT_metDscale%fine_ts   = LDT_rc%met_ts(1)
         LDT_metDscale%coarse_ts = LDT_rc%met_ts(2)
         LDT_metDscale%fine_metname = LDT_rc%metforc(1)
         LDT_metDscale%coarse_metname = LDT_rc%metforc(2)

         write(LDT_logunit,*) "[INFO] Fine ts:", LDT_metDscale%fine_ts,&
                  "... for Forcing dataset: ",trim(LDT_metDscale%fine_metname)
         write(LDT_logunit,*) "[INFO] Coarse ts:", LDT_metDscale%coarse_ts,&
                  "... for Forcing dataset: ",trim(LDT_metDscale%coarse_metname)

         ! Determine how many timesteps occur for fine_ts to coarse_ts:
         LDT_metDscale%numfine2coarse = &
             nint( LDT_metDscale%coarse_ts/LDT_metDscale%fine_ts )

         write(LDT_logunit,*) "[INFO] Number of fine timesteps to coarse timesteps:",&
               LDT_metDscale%numfine2coarse
      else
        write(LDT_logunit,*) "[ERR]  The fine timescale forcing dataset "
        write(LDT_logunit,*) " needs to be listed first in the ldt.config file."
        write(LDT_logunit,*) " Program stopping ..."
        call LDT_endrun
      endif

    ! Option not available ...
    case default
       write(LDT_logunit,*) "[ERR] This temporal downscaling method, ",&
             trim(LDT_rc%timeDscaleType)
       write(LDT_logunit,*) " is not recognized. Program stopping ..."
       call LDT_endrun

   end select

  end subroutine LDT_metforcDscaleInit


!BOP
! 
! !ROUTINE: add_forcing_fields
! \label{add_forcing_fields}
! 
! !INTERFACE: 
  subroutine add_forcing_fields(n, FORC_State, forc_attrib)
! !USES:     

    implicit none
! !ARGUMENTS: 
    integer                   :: n 
    type(ESMF_State)          :: FORC_State
    type(LDT_forcDataEntry)   :: forc_attrib
! 
! !DESCRIPTION: 
!  This subroutine creates fields for the specified downscaling 
!  forcing variables and initializes them to be undefined values. 
!
!EOP

    integer                :: k,m
    type(ESMF_Field)       :: varField1, varField2
    type(ESMF_Field)       :: sumVarField1,   sumVarField2
    type(ESMF_Field)       :: countVarField1, countVarField2
    type(ESMF_ArraySpec)   :: arrspec1
    real,        pointer   :: var(:)
    integer                :: status

    if(forc_attrib%selectOpt.eq.1) then 

       call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LDT_verify(status,&
            'error in ESMF_ArraySpecSet in add_forcing_fields')       

       do k=1,forc_attrib%vlevels

     ! -- Create ESMF Fields:

        ! Variable Field 1:
          varField1 = ESMF_FieldCreate(grid=LDT_vecTile(n), &
               arrayspec=arrspec1, name = &
               trim(forc_attrib%varname(k))//"_1",&
               rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldCreate in add_forcing_fields')

          sumVarField1 = ESMF_FieldCreate(grid=LDT_vecTile(n), &
               arrayspec=arrspec1, name = &
               "Sum_"//trim(forc_attrib%varname(k))//"_1",&
               rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldCreate in add_forcing_fields')

          countVarField1 = ESMF_FieldCreate(grid=LDT_vecTile(n), &
               arrayspec=arrspec1, name = &
               "Count_"//trim(forc_attrib%varname(k))//"_1",&
               rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldCreate in add_forcing_fields')

        ! Variable Field 2:
          varField2 = ESMF_FieldCreate(grid=LDT_vecTile(n), &
               arrayspec=arrspec1, name = &
               trim(forc_attrib%varname(k))//"_2",&
               rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldCreate in add_forcing_fields')

          sumVarField2 = ESMF_FieldCreate(grid=LDT_vecTile(n), &
               arrayspec=arrspec1, name = &
               "Sum_"//trim(forc_attrib%varname(k))//"_2",&
               rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldCreate in add_forcing_fields')

          countVarField2 = ESMF_FieldCreate(grid=LDT_vecTile(n), &
               arrayspec=arrspec1, name = &
               "Count_"//trim(forc_attrib%varname(k))//"_2",&
               rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldCreate in add_forcing_fields')

     ! -- Initialize ESMF Fields and Associated Pointer Arrays:

          call ESMF_FieldGet(varField1,localDE=0,farrayPtr=var,rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldGet in add_forcing_fields')
          var = 0.0

          call ESMF_FieldGet(varField2,localDE=0,farrayPtr=var,rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldGet in add_forcing_fields')
          var = 0.0

          call ESMF_FieldGet(sumVarField1,localDE=0,farrayPtr=var,rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldGet in add_forcing_fields')
          var = 0.0

          call ESMF_FieldGet(sumVarField2,localDE=0,farrayPtr=var,rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldGet in add_forcing_fields')
          var = 0.0

          call ESMF_FieldGet(countVarField1,localDE=0,farrayPtr=var,rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldGet in add_forcing_fields')
          var = 0.0

          call ESMF_FieldGet(countVarField2,localDE=0,farrayPtr=var,rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldGet in add_forcing_fields')
          var = 0.0

     ! -- Add ESMF State for use later:

          call ESMF_StateAdd(FORC_State,(/varField1/),rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateAdd in add_forcing_fields')

          call ESMF_StateAdd(FORC_State,(/varField2/),rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateAdd in add_forcing_fields')

          call ESMF_StateAdd(FORC_State,(/sumVarField1/),rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateAdd in add_forcing_fields')

          call ESMF_StateAdd(FORC_State,(/sumVarField2/),rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateAdd in add_forcing_fields')

          call ESMF_StateAdd(FORC_State,(/countVarField1/),rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateAdd in add_forcing_fields')

          call ESMF_StateAdd(FORC_State,(/countVarField2/),rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateAdd in add_forcing_fields')

       enddo
    endif

  end subroutine add_forcing_fields

!BOP
! 
! !ROUTINE: LDT_diagnoseForMetForcDscale(n)
! \label{LDT_diagnoseForMetForcDscale}
! 
! !INTERFACE: 
  subroutine LDT_diagnoseForMetForcDscale(n)
! !USES:   

! !ARGUMENTS: 
    integer, intent(in)        :: n 
! 
! !DESCRIPTION: 
!  This subroutine sums and counts each forcing's value
!   at each timestep through the time window.
!
!EOP
    integer                    :: i,t,m,fobjcount
    character*100, allocatable :: forcobjs(:)       ! The saved array of forcing fields
    type(ESMF_Field)           :: base1Field, base2Field
    type(ESMF_Field)           :: dscale1Field, dscale2Field
    type(ESMF_Field)           :: sumDscale1Field,   sumDscale2Field
    type(ESMF_Field)           :: countDscale1Field, countDscale2Field
    real,          pointer     :: forcdata_base1(:), forcdata_base2(:)
    real,          pointer     :: forcdata_dscale1(:), forcdata_dscale2(:)
    real,          pointer     :: forcdata_sum_dscale1(:), forcdata_sum_dscale2(:)
    real,          pointer     :: forcdata_count_dscale1(:), forcdata_count_dscale2(:)
    integer                    :: status, status1, status2
    real                       :: factor1, factor2
    logical                    :: attrib_exists
! ____________________________________________________________________
    
    call ESMF_StateGet(LDT_FORC_State(n),itemCount=fobjcount,rc=status)
    call LDT_verify(status,'ESMF_StateGet failed for fobjcount in LDT_diagnoseForMetForcDscale')

    allocate(forcobjs(fobjcount))

    call ESMF_StateGet(LDT_FORC_State(n),itemNameList=forcobjs,rc=status)
    call LDT_verify(status,'ESMF_StateGet failed for forcobjs in LDT_diagnoseForMetForcDscale')

    ! Loop over the number of forcing types:
    do i=1,fobjcount

       call ESMF_StateGet(LDT_FORC_Base_State(n,1),trim(forcobjs(i)), &
            base1Field, rc=status1 )

       call ESMF_StateGet(LDT_FORC_Base_State(n,2),trim(forcobjs(i)), &
            base2Field, rc=status1 )

       call ESMF_StateGet(LDT_FORC_Dscale_State(n),trim(forcobjs(i))//"_1",&
            dscale1Field, rc=status1 )

       call ESMF_StateGet(LDT_FORC_Dscale_State(n),trim(forcobjs(i))//"_2",&
            dscale2Field, rc=status1 )

       call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Sum_"//trim(forcobjs(i))//"_1",&
            sumDscale1Field, rc=status1 )

       call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Sum_"//trim(forcobjs(i))//"_2",&
            sumDscale2Field, rc=status1 )

       call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Count_"//trim(forcobjs(i))//"_1",&
            countDscale1Field, rc=status1 )

       call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Count_"//trim(forcobjs(i))//"_2",&
            countDscale2Field, rc=status1 )

       if( status1.eq.0 ) then 
 
          call ESMF_FieldGet(base1Field,localDE=0,farrayPtr=forcdata_base1, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet: base1Field failed in LDT_diagnoseForMetForcDscale')

          call ESMF_FieldGet(base2Field,localDE=0,farrayPtr=forcdata_base2, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet: base2Field failed in LDT_diagnoseForMetForcDscale')

          call ESMF_FieldGet(dscale1Field,localDE=0,farrayPtr=forcdata_dscale1, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet: dscale1Field failed in LDT_diagnoseForMetForcDscale')

          call ESMF_FieldGet(dscale2Field,localDE=0,farrayPtr=forcdata_dscale2, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet: dscale2Field failed in LDT_diagnoseForMetForcDscale')

          call ESMF_FieldGet(sumdscale1Field,localDE=0,farrayPtr=forcdata_sum_dscale1, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet: sumdscale1Field failed in LDT_diagnoseForMetForcDscale')

          call ESMF_FieldGet(sumdscale2Field,localDE=0,farrayPtr=forcdata_sum_dscale2, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet: sumdscale2Field failed in LDT_diagnoseForMetForcDscale')

          call ESMF_FieldGet(countdscale1Field,localDE=0,farrayPtr=forcdata_count_dscale1, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet: countdscale1Field failed in LDT_diagnoseForMetForcDscale')

          call ESMF_FieldGet(countdscale2Field,localDE=0,farrayPtr=forcdata_count_dscale2, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet: countdscale2Field failed in LDT_diagnoseForMetForcDscale')

          factor1 = 1.0
          factor2 = 1.0
          ! Assign forcing timescale as factor when accounting for rainfall:
          if( index(forcobjs(i),"Rainfall") > 0. ) then
            factor1 = LDT_metDscale%fine_ts   
            factor2 = LDT_metDscale%coarse_ts 
          endif

          do t=1,LDT_rc%ntiles(n)
           ! Field 1 (F1) -- finer timescale:
             if(forcdata_base1(t).ne.LDT_rc%udef) then 
              ! Sum to estimate weights:
                forcdata_dscale1(t) = forcdata_base1(t) * factor1
                forcdata_sum_dscale1(t) = forcdata_sum_dscale1(t) &
                        + ( forcdata_base1(t)*factor1 )                    ! Sum F1
                forcdata_count_dscale1(t) = forcdata_count_dscale1(t) + 1  ! Count F1
             endif
           ! Field 2 (F2) -- coarser timescale:
             if(forcdata_base2(t).ne.LDT_rc%udef) then 
              ! Sum to estimate weights:
                forcdata_dscale2(t) = forcdata_base2(t) * factor2
                forcdata_sum_dscale2(t) = forcdata_sum_dscale2(t) &
                        + ( forcdata_base2(t)*factor2 )                    ! Sum F2
                forcdata_count_dscale2(t) = forcdata_count_dscale2(t) + 1  ! Count F2
             endif
          enddo
       endif
    enddo

    deallocate(forcobjs)
    
  end subroutine LDT_diagnoseForMetForcDscale


!BOP
! 
! !ROUTINE: LDT_computeMetForcDscaleParams(pass)
! \label{LDT_computeMetForcDscaleParams}
! 
! !INTERFACE: 
  subroutine LDT_computeMetForcDscaleParams(pass)
! !USES:   
!
! !ARGUMENTS: 
    integer, intent(in)  :: pass
! 
! !DESCRIPTION: 
!  This subroutine sums and counts each forcing's value
!   at each timestep through the time window.
!
!EOP
    integer                    :: n 
    integer                    :: i,t,m,fobjcount
    character*100, allocatable :: forcobjs(:)
    type(ESMF_Field)           :: dscale1Field, dscale2Field
    type(ESMF_Field)           :: sumDscale1Field, sumDscale2Field
    type(ESMF_Field)           :: countDscale1Field, countDscale2Field
    real,          pointer     :: forcdata_dscale1(:), forcdata_dscale2(:)
    real,          pointer     :: forcdata_sum_dscale1(:), forcdata_sum_dscale2(:)
    real,          pointer     :: forcdata_count_dscale1(:), forcdata_count_dscale2(:)
    integer                    :: status, status1, status2

    if(pass.eq.1) then    ! Do only for pass 1

       do n=1,LDT_rc%nnest
          call ESMF_StateGet(LDT_FORC_State(n),itemCount=fobjcount,rc=status)
          call LDT_verify(status,&
              'ESMF_StateGet failed for fobjcount in LDT_computeMetForcDscaleParams')
          
          allocate(forcobjs(fobjcount))
          
          call ESMF_StateGet(LDT_FORC_State(n),itemNameList=forcobjs,rc=status)
          call LDT_verify(status,&
              'ESMF_StateGet failed for forcobjs in LDT_computeMetForcDscaleParams')
          
          do i=1,fobjcount  ! Loop over each forcing field type
             
           ! Read-in/"get" forcing fields:
             call ESMF_StateGet(LDT_FORC_Dscale_State(n),trim(forcobjs(i))//"_1",dscale1Field,&
                  rc=status1)
             
             call ESMF_StateGet(LDT_FORC_Dscale_State(n),trim(forcobjs(i))//"_2",dscale2Field,&
                  rc=status1)

             call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Sum_"//trim(forcobjs(i))//"_1",&
                  sumDscale1Field,rc=status1)

             call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Sum_"//trim(forcobjs(i))//"_2",&
                  sumDscale2Field,rc=status1)

             call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Count_"//trim(forcobjs(i))//"_1",&
                  countDscale1Field,rc=status1)

             call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Count_"//trim(forcobjs(i))//"_2",&
                  countDscale2Field,rc=status1)
             
             if( status1.eq.0 ) then 
                call ESMF_FieldGet(dscale1Field,localDE=0,farrayPtr=forcdata_dscale1, &
                     rc=status2)
                call LDT_verify(status2,'ESMF_FieldGet failed in LDT_computeMetForcDscaleParams')
                
                call ESMF_FieldGet(dscale2Field,localDE=0,farrayPtr=forcdata_dscale2, &
                     rc=status2)
                call LDT_verify(status2,'ESMF_FieldGet failed in LDT_computeMetForcDscaleParams')

                call ESMF_FieldGet(sumdscale1Field,localDE=0,farrayPtr=forcdata_sum_dscale1, &
                     rc=status2)
                call LDT_verify(status2,'ESMF_FieldGet failed in LDT_computeMetForcDscaleParams')

                call ESMF_FieldGet(sumdscale2Field,localDE=0,farrayPtr=forcdata_sum_dscale2, &
                     rc=status2)
                call LDT_verify(status2,'ESMF_FieldGet failed in LDT_computeMetForcDscaleParams')

                call ESMF_FieldGet(countdscale1Field,localDE=0,farrayPtr=forcdata_count_dscale1, &
                     rc=status2)
                call LDT_verify(status2,'ESMF_FieldGet failed in LDT_computeMetForcDscaleParams')

                call ESMF_FieldGet(countdscale2Field,localDE=0,farrayPtr=forcdata_count_dscale2, &
                     rc=status2)
                call LDT_verify(status2,'ESMF_FieldGet failed in LDT_computeMetForcDscaleParams')
                
!                do t=1,LDT_rc%ntiles(n)
!                   if(forcdata_dscale2(t).ne.0) then 
!                    forcdata_sum_dscale1(t) = forcdata_dscale1(t)  ! Assign Sum of F1 
!                   endif
!                enddo
                
                forcdata_dscale1 = 0 
                forcdata_dscale2 = 0 
!                forcdata_count_dscale1 = 0 
!                forcdata_count_dscale2 = 0 
             endif
          enddo   ! End forcing type loop
          
          deallocate(forcobjs)
       enddo   ! End nest loop
    endif
    
  end subroutine LDT_computeMetForcDscaleParams


!BOP
! 
! !ROUTINE:  LDT_applyDscaleCorrection(n)
! \label{LDT_applyDscaleCorrection}
! 
! !INTERFACE: 
  subroutine LDT_applyDscaleCorrection(n)
! !USES:   
!
! !ARGUMENTS: 
    integer, intent(in)    :: n 
! 
! !DESCRIPTION: 
!   This routine applies the temporal downscaling method
!   selected.  Then the final basefields are "overlayed"
!   here for writing out in subsequent steps.
!
!EOP
  integer                :: fobjcount
  integer                :: i,t,m
  type(ESMF_Field)       :: mrgField, baseField
  integer                :: status, status1, status2
  real, pointer          :: forcdata_base(:), forcdata_mrg(:)
  character*100, allocatable :: forcobjs(:)
! ______________________________________________________________


  ! Select/apply temporal downscaling method via 
  !  C function-table (in LDT_metforcdscale_FTable.c):
  call applytimedscale( trim(LDT_rc%timeDscaleType)//char(0), n )
  !

  ! Overlay forcings:
  call ESMF_StateGet(LDT_FORC_State(n),itemCount=fobjcount,rc=status)
  call LDT_verify(status,'ESMF_StateGet failed for fobjcount in LDT_applyDscaleCorrection')

  allocate(forcobjs(fobjcount))

  call ESMF_StateGet(LDT_FORC_State(n),itemNameList=forcobjs,rc=status)
  call LDT_verify(status,'ESMF_StateGet failed for forcobjs in LDT_applyDscaleCorrection')

  ! Loop over forcing type (e.g., air temp, precip, etc.):
  do i=1,fobjcount

     call ESMF_StateGet(LDT_FORC_State(n),forcobjs(i),mrgField,&
          rc=status)
     call LDT_verify(status, 'ESMF_StateGet failed for '//trim(forcobjs(i)))
     call ESMF_FieldGet(mrgField,localDE=0,farrayPtr=forcdata_mrg, &
            rc=status)
     call LDT_verify(status,'ESMF_FieldGet failed for forcdata_mrg in LDT_applyDscaleCorrection')

     ! Loop over stored metforcing dataset types:
     do m=1,LDT_rc%nmetforc
        call ESMF_StateGet(LDT_FORC_Base_State(n,m),forcobjs(i),baseField,&
             rc=status1)

        if( status1.eq.0 ) then
           call ESMF_FieldGet(baseField,localDE=0,farrayPtr=forcdata_base, &
                rc=status2)
           call LDT_verify(status2,'ESMF_FieldGet (basefield) failed in LDT_applyDscaleCorrection')

           do t=1,LDT_rc%ntiles(n)
              if(forcdata_base(t).ne.-9999.0) then
                 forcdata_mrg(t) = forcdata_base(t)
              endif
           enddo
        endif

     enddo   ! End metforcing dataset loop
  enddo      ! End forcing type loop

  deallocate(forcobjs)
    
 end subroutine LDT_applyDscaleCorrection


!BOP
! !ROUTINE: LDT_setupMetDscaleTimeWindow(interval, twstart, twend)
! \label{LDT_setupMetDscaleTimeWindow}
!      
! !INTERFACE:            
  subroutine LDT_setupMetDscaleTimeWindow(interval, twstart, twend)

! !USES:
   use LDT_coreMod, only : LDT_rc

! !ARGUMENTS:
   real, intent(in)    :: interval
   type(ESMF_Time)     :: twstart
   type(ESMF_Time)     :: twend
!
! !DESCRIPTION:
!  This routine calculates the initial time-window ("TW") 
!   start and end dates and times, based on the initial 
!   start time specified in the config file.  The initial
!   TWend date and time are established by searching for
!   when the processing time interval (e.g., 86400 sec)
!   reaches 00Z, for example. 
!
!  Another example, would be if the processing time interval 
!   is at 6-hours * 3600 secs and looping over 24-hour period, 
!   where the time window would be valid at 00Z, 6Z, 12Z, 18Z. 
!   
!   The arguments are:
!   \begin{description}
!   \item [interval]
!     Time window interval for meteorological forcing temporal downscaling.
!   \item [twstart]
!     Time window start date/time (ESMF_Time).
!   \item [twend]
!     Time window ending date/time (ESMF_Time).
!   \end{description}
!    
!EOP
    type(ESMF_Time)         :: tTime
    type(ESMF_TimeInterval) :: timestep
    logical                 :: check_flag
    integer                 :: status
    integer                 :: yr, mo, da, hr, mn, ss
    integer                 :: tstep_n_secs
    integer                 :: interval_int
! ___________________________________________________________

    write(LDT_logunit,*) "[INFO] Initializing LDT Time Window Book-ends "

    ! Initialize the time-window (TW) start date/time with
    !  current start time:
    call ESMF_ClockGet(LDT_clock, currTime=twStart, &
         rc=status)

    call ESMF_TimeGet(twStart, yy=YR, mm=MO, dd=DA, h=HR)
    print *, "The clock's start time:", YR, MO, DA, HR

#if 0
! (KRA)
    if( interval == 86400 ) then   ! Day in Seconds 
       print *, " Interval == Day (in seconds) "
       print *, "(init)TW bookends, interval: ",interval
       print *, "(init)TW bookends, twstarthr,twendhr: ", &
             LDT_rc%metForcTWstarthour, LDT_rc%metForcTWendhour
       print *, "(config) Start month, day, hour: ", LDT_rc%mo, LDT_rc%da, LDT_rc%hr

       if( LDT_rc%hr .ne. LDT_rc%metForcTWstarthour ) then
         print *, " LDT_config starting hour and Time Window starting hours don't match "
         print *, "  To ensure proper time-window start and end setting and updating, it "
         print *, "  recommended that the user specify the same starting hours."
         print *, " LDT run stopping ... "
         call LDT_endrun()
       endif

       ! Check if TW start and end hours are the same; 
       !  If so, advance day for end hour book-end.
       interval_int = 0
       if( LDT_rc%metForcTWstarthour == LDT_rc%metForcTWendhour ) then
          interval_int = int(interval)
       endif

       ! Set ending time-window bookend date/time (ESMF):
       call ESMF_TimeSet(twend, yy = LDT_rc%yr, &
            mm = LDT_rc%mo, dd = LDT_rc%da, &
            h = LDT_rc%metForcTWendhour,&
            m = 0, s = interval_int, &
            calendar = LDT_calendar,&
            rc=status)
       call LDT_verify(status,'Error in ESMF_TimeSet(TWend):in LDT_setupMetDscaleTimeWindow')

       call ESMF_TimeGet(twend, yy = yr, &
            mm = mo, dd = da, h = hr, &
            m = mn, s = ss, &
            calendar = LDT_calendar, &
            rc=status)
       call LDT_verify(status,&
           'error in ESMF_TimeGet: LDT_setupMetDscaleTimeWindow)')
       print *, " TWend (LDT_setupMetDscaleTimeWindow): ",mo, da, hr, interval

    else
! (KRA)
#endif

    ! The LDT clock timestep hasn't been set yet. So using a
    ! small timestep (60 seconds) to initialize the calculation:
    call ESMF_TimeIntervalSet(timeStep,s=60,rc=status)

    ! Initialize parameters for loop:
    tTime = twStart
    check_flag = .true.
    do while(check_flag)

       ! Search for time-window ending book-end date/time:
       tTime = tTime + timeStep

       call ESMF_TimeGet(tTime, yy = yr, &
            mm = mo, dd = da, h = hr, &
            m = mn, s = ss, &
            calendar = LDT_calendar, &
            rc=status)
       call LDT_verify(status,&
           'error in ESMF_TimeGet: LDT_timeMgrMod(LDT_setupTimeWindow)')

       ! Set date for the end of the time-window when the 
       !  processing interval (e.g., 86400 secs) reaches
       !  either: 1) hr=00, min=00, sec=00; or
       !          2) processing interval divides nicely into ...?
       ! 
       if(mod(real(hr)*3600+60*real(mn)+    &
            real(ss), interval).eq.0.0) then

          ! Set end date/time for TW:
          twend = tTime
          check_flag = .false.   ! Exit loop
          exit
       endif
    enddo

!   endif ! (KRA)

 end subroutine LDT_setupMetDscaleTimeWindow


!BOP
! !ROUTINE: LDT_endofMetDscaleTimeWindow
! \label{LDT_endofMetDscaleTimeWindow}
!      
! !INTERFACE:            
  function LDT_endofMetDscaleTimeWindow() result(finish)
! !USES:
    use LDT_coreMod, only : LDT_rc

! !ARGUMENTS:
    logical :: finish
!
! !DESCRIPTION:
!  This function checks to see if the runtime clock has reached the
!  specified stop time of the simulation. 
! 
!   The arguments are:
!   \begin{description}
!   \item [finish]
!     boolean value indicating if the end of simulation is reached. 
!   \end{description}
!    
!  The calling sequence is: 
!  \begin{description}
!   \item[ESMF\_ClockGet] (\ref{ESMF_ClockGet}) \newline
!    Retrieves current LDT clock time and timestep.
!   \end{description}
!
!EOP
    type(ESMF_Time)         :: currTime
    type(ESMF_TimeInterval) :: ts
    integer  :: rc
    integer  :: YY, MM, DD, H, M, S
    integer  :: yr, mo, da, hr
    integer  :: ewyr, ewmo, ewda, ewhr
! ___________________________________________________________

    call ESMF_ClockGet(LDT_clock, currTime=currTime, &
         timestep = ts, rc=rc)

!    call ESMF_TimeGet(LDT_metProcTWend, yy=ewyr, mm=ewmo, dd=ewda, h=ewhr)

    ! Check if LDT current time >= ending TW bookend:
    if((currTime).ge.LDT_metProcTWend) then 
!     call ESMF_TimeGet(currTime, yy=yr, mm=mo, dd=da, h=hr, rc=rc) 
!      print *, " Current time    :", yr, mo, da, hr
       finish = .true. 
    else
       finish = .false. 
    endif
 
  end function LDT_endofMetDscaleTimeWindow

!BOP
! !ROUTINE: LDT_resetMetDscaleTimeWindow(pass)
! \label{LDT_resetMetDscaleTimeWindow}
!      
! !INTERFACE:            
  subroutine LDT_resetMetDscaleTimeWindow(pass)

! !USES:

! !ARGUMENTS:
   integer, intent(in)  :: pass
!
! !DESCRIPTION:
!  This routine resets and updates the bookend dates and times
!   of the time window.
! 
!   The arguments are:
!   \begin{description}
!   \item [pass]
!     Time window pass on temporal downscaling steps.
!   \end{description}
!    
!EOP
    integer                 :: rc,status
    integer                 :: yr,mo,da,hr,mn,ss
    integer                 :: interval_int
    logical                 :: check_flag
    type(ESMF_Time)         :: tTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: tstep_n_secs
! ____________________________________________________

    ! Set the clock at the start Time:
    if(pass.eq.1) then 

       call ESMF_TimeGet(LDT_metProcTWstart, &
            yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, calendar=LDT_calendar,&
            rc=status)

       call LDT_timemgr_set(LDT_rc,yr,mo,da,hr,mn,ss,0,0.0)

! ------------------------------------------------

    ! Second "round" or pass through the time window:
    elseif(pass.eq.2) then 


#if 0 
      ! Treat differently for daily, (or monthly) files:
      ! (KRA)
      if( LDT_rc%metForcProcInterval == 86400 ) then
        interval_int = int(LDT_rc%metForcProcInterval)
        ! TWstart ::
        ! Retrieve date/time info for TW start term:
        call ESMF_TimeGet(LDT_metProcTWstart,  yy=yr, &
                mm = mo, dd = da, h = hr, &
                m = mn, s = ss,&
                calendar = LDT_calendar, &
                rc=status)

        ! Update TW start date/time:
        call ESMF_TimeSet(LDT_metProcTWstart, yy=yr, &
                mm = mo, dd = da, h = hr, &
                m = mn, s = interval_int, &
                calendar = LDT_calendar, &
                rc=status)
        call LDT_verify(status, 'error in ESMF_TimeSet: LDT_resetMetDscaleTimeWindow')

        ! Print updated TWstart:
        call ESMF_TimeGet(LDT_metProcTWstart,  yy=yr, &
                mm = mo, dd = da, h = hr, &
                m = mn, s = ss,&
                calendar = LDT_calendar, &
                rc=status)
!         print *, " -- pass 2 updated TWstart:", YR, MO, DA, HR

        ! TWend ::
        ! Retrieve date/time info for TW ending term:
        call ESMF_TimeGet(LDT_metProcTWend,  yy=yr, &
                mm = mo, dd = da, h = hr, &
                m = mn, s = ss,&
                calendar = LDT_calendar, &
                rc=status)
!        print *, " -- pass 2 Orig TWend:", YR, MO, DA, HR

        ! Update TW end date/time:
        call ESMF_TimeSet(LDT_metProcTWend, yy=yr, &
                mm = mo, dd = da, h = hr, &
                m = mn, s = interval_int, &
                calendar = LDT_calendar, &
                rc=status)
        call LDT_verify(status, 'error in ESMF_TimeSet: LDT_resetMetDscaleTimeWindow')

    ! Subdaily processing window code (default for now)
    else
#endif

       ! Main code:
       ! Update TW start date/time using the LDT current "clock" date/time:
       call ESMF_ClockGet(LDT_clock, currTime=LDT_metProcTWstart, &
            timeStep = timeStep, rc=rc)

       tTime = LDT_metProcTWstart

       check_flag = .true. 
       do while(check_flag) 

          ! Search for next valid TW end date/time using LDT run timestep (e.g., 3600s):
          tTime = tTime + timeStep
          call ESMF_TimeGet(tTime,  yy=yr, &
               mm = mo, dd = da, h = hr, &
               m = mn, s = ss,&
               calendar = LDT_calendar, &
               rc=status)
          call LDT_verify(status, 'error in ESMF_TimeGet: LDT_resetMetDscaleTimeWindow')

          ! When the accumulated time is the same as the "processing" interval,
          !  set the final accumulated time to the final "time-window" bookend:
          if(mod(real(hr)*3600+60*real(mn)+    &
               real(ss), LDT_rc%metForcProcInterval).eq.0.0) then 
             LDT_metProcTWend = tTime
             check_flag = .false. 
             exit
          endif
       enddo

    endif   ! End pass conditional check

  end subroutine LDT_resetMetDscaleTimeWindow


!BOP
! !ROUTINE: LDT_resetDscaleVars
! \label{LDT_resetDscaleVars}
! 
! !INTERFACE: 
  subroutine LDT_resetDscaleVars(pass)

! !ARGUMENTS:
  integer, intent(in)  :: pass
! 
! !DESCRIPTION: 
!   This routine resets the specified downscaling variables
!    for the next history output step. 
!
!   The arguments are: 
!   \begin{description}
!   \end{description}
!EOP

    integer                    :: i,n,fobjcount
    character*100, allocatable :: forcobjs(:)
    type(ESMF_Field)           :: sumDscale1Field,   sumDscale2Field
    type(ESMF_Field)           :: countDscale1Field, countDscale2Field
    real,          pointer     :: forcdata_sum_dscale1(:), forcdata_sum_dscale2(:)
    real,          pointer     :: forcdata_count_dscale1(:), forcdata_count_dscale2(:)
    integer                    :: status, status1, status2

    if( pass.eq.2 ) then    ! Do only for pass 2

       do n=1,LDT_rc%nnest
          call ESMF_StateGet(LDT_FORC_State(n),itemCount=fobjcount,rc=status)
          call LDT_verify(status,&
              'ESMF_StateGet failed for fobjcount in LDT_metforcDscaleMod(LDT_resetDscaleVars)')

          allocate(forcobjs(fobjcount))

          call ESMF_StateGet(LDT_FORC_State(n),itemNameList=forcobjs,rc=status)
          call LDT_verify(status,&
              'ESMF_StateGet failed for forcobjs in LDT_metforcDscaleMod(LDT_resetDscaleVars)')

          do i=1,fobjcount  ! Loop over each forcing field type

             call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Sum_"//trim(forcobjs(i))//"_1",&
                  sumDscale1Field,rc=status1)

             call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Sum_"//trim(forcobjs(i))//"_2",&
                  sumDscale2Field,rc=status1)

             call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Count_"//trim(forcobjs(i))//"_1",&
                  countDscale1Field,rc=status1)

             call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Count_"//trim(forcobjs(i))//"_2",&
                  countDscale2Field,rc=status1)
      
             if( status1.eq.0 ) then

                call ESMF_FieldGet(sumdscale1Field,localDE=0,farrayPtr=forcdata_sum_dscale1, &
                     rc=status2)
                call LDT_verify(status2,'ESMF_FieldGet failed in LDT_resetDscaleVars')

                call ESMF_FieldGet(sumdscale2Field,localDE=0,farrayPtr=forcdata_sum_dscale2, &
                     rc=status2)
                call LDT_verify(status2,'ESMF_FieldGet failed in LDT_resetDscaleVars')

                call ESMF_FieldGet(countdscale1Field,localDE=0,farrayPtr=forcdata_count_dscale1, &
                     rc=status2)
                call LDT_verify(status2,'ESMF_FieldGet failed in LDT_resetDscaleVars')

                call ESMF_FieldGet(countdscale2Field,localDE=0,farrayPtr=forcdata_count_dscale2, &
                     rc=status2)
                call LDT_verify(status2,'ESMF_FieldGet failed in LDT_resetDscaleVars')

              ! Reset sums and counts for two fields:
                forcdata_sum_dscale1 = 0
                forcdata_sum_dscale2 = 0
                forcdata_count_dscale1 = 0 
                forcdata_count_dscale2 = 0 

             endif

          enddo   ! End forcing type loop
       enddo
    endif

  end subroutine LDT_resetDscaleVars


end module LDT_metforcDscaleMod

