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
module SACHTET_parmsMod
!BOP
!
! !MODULE: SACHTET_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read SAC-HTET parameter
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  SAC-HTET parameter file data.
!
! !REVISION HISTORY:
!
!  04 Nov 2013: K. Arsenault: Added layers for SAC-HTET model
!
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_soilsMod
  use LDT_paramMaskCheckMod
  use LDT_xmrg_reader

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: SACHTETparms_init    !allocates memory for required structures
  public :: SACHTETparms_writeHeader
  public :: SACHTETparms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: SACHTET_struc

  type, public :: sachtet_type_dec

     real          :: pet_gridDesc(20)
     character*50  :: pet_proj
     character*50  :: pet_gridtransform

     character(len=LDT_CONST_PATH_LEN) :: petdir
     character(len=LDT_CONST_PATH_LEN) :: petfile
     character(len=LDT_CONST_PATH_LEN) :: petadjdir
     character(len=LDT_CONST_PATH_LEN) :: petadjfile
     character*20  :: petInterval

     real          :: sachtetparms_gridDesc(20)
     character*50  :: sachtetparms_proj
     character*50  :: sachtetparms_gridtransform

     character*50  :: tbot_topocorr

! - RDHM 3.5.6
     character*140 :: rdhmconsts_table
     real          :: rdhm_undef

!  - SAC-HTET parameters:
     character*140 :: csoilparm_method
!     character*140 :: cosbysoils_table
     character*140 :: sacsoilparms_table
     character*140 :: sacvegparms_table
     character*140 :: lzfpm_file      ! Lower zone 
     character*140 :: lzfsm_file      !
     character*140 :: lzpk_file       !
     character*140 :: lzsk_file       !
     character*140 :: lztwm_file      !
     character*140 :: uzfwm_file      !
     character*140 :: uztwm_file      !
     character*140 :: uzk_file        !
     character*140 :: pfree_file      !
     character*140 :: rexp_file       !
     character*140 :: zperc_file      !

     character*140 :: sacmask_file    ! SAC HTET mask file
     character*140 :: pctim_file      !
     character*140 :: efc_file        ! Forest fraction cover
     character*140 :: soilalb_file    !
     character*140 :: timeOffset_file !

     character*140 :: stxt_file       !
     character*140 :: tbot_file       !
     character*140 :: zbot_file       !
     character*140 :: adimp_file      !
     character*140 :: side_file       !
     character*140 :: riva_file       !
     character*140 :: rserv_file      !
     character*140 :: rsmax_file      !
     character*140 :: cksl_file       !

! -  SAC-HTET model-specific:
     type(LDT_paramEntry) :: sachtet356  ! SAC-HTET v.3.5.6 model version parameters (collective)
     type(LDT_paramEntry) :: sacmask     ! SAC-HTET landmask (e.g.,mask_land_newGRN-9)

     type(LDT_paramEntry) :: lzfpm       ! Lower zone primary free water (slow) max. storage (mm)
     type(LDT_paramEntry) :: lzfsm       ! Lower zone supplemental free water (fast) max. storage (mm)
     type(LDT_paramEntry) :: lzpk        ! Lower zone primary free water depletion rate (day^-1)
     type(LDT_paramEntry) :: lzsk        ! Lower zone supplemental free water depletion rate (day^-1)
     type(LDT_paramEntry) :: lztwm       ! Lower zone tension water max. storage (mm)
     type(LDT_paramEntry) :: pfree       ! Fraction percolation from upper to lower free water storage (day^-1)
     type(LDT_paramEntry) :: rexp        ! Expon. of the percolation eqtn.
     type(LDT_paramEntry) :: uzfwm       ! Upper zone free water max. storage (mm)
     type(LDT_paramEntry) :: uzk         ! Upper zone free water latent depletion rate (day^-1)
     type(LDT_paramEntry) :: uztwm       ! Upper zone tension water max. storage (mm)
     type(LDT_paramEntry) :: zperc       ! Max. percolation rate (mm)

     type(LDT_paramEntry) :: txtlw       ! Dominant texture in lower zone
     type(LDT_paramEntry) :: txtot       ! Dominant texture in total column
     type(LDT_paramEntry) :: dzup        ! Thickness of depth in upper zone
     type(LDT_paramEntry) :: dzlw        ! Thickness of depth in lower zone

     type(LDT_paramEntry) :: forestfrac  ! forest fraction (EFC)
     type(LDT_paramEntry) :: pctim       ! Impervious fraction of the watershed area (-)
     type(LDT_paramEntry) :: adimp       ! Additional impervious area (-)
     type(LDT_paramEntry) :: riva        ! Riparian vegetation area (-) 
     type(LDT_paramEntry) :: soilalb     ! Soil (snow-free) albedo (-)
     type(LDT_paramEntry) :: rserv       ! LZ free water fraction not transferable to LZ tension water (-)
     type(LDT_paramEntry) :: side        ! Ratio of deep recharge to channel base flow (-)
     type(LDT_paramEntry) :: offsetTime  ! Time offset (to account for UTC->local time conversion) 
     type(LDT_paramEntry) :: frz_soiltext ! Frozen ground - Soil texture
     type(LDT_paramEntry) :: frz_cksl    ! Frozen ground - ratio of frozen to non-frozen surface
     type(LDT_paramEntry) :: frz_rsmax   ! Frozen ground - residual porosity
     type(LDT_paramEntry) :: frz_zbot    ! Frozen ground - lower boundary depth (neg. value)

     type(LDT_paramEntry) :: fxexp       ! bare soil evaporation exponent (-)
     type(LDT_paramEntry) :: rcmax       ! Maximum stomatal resistance (s/m)
     type(LDT_paramEntry) :: topt        ! Optimum air temperature for transpiration (K)
     type(LDT_paramEntry) :: plantcoef   ! Plant coefficient (-)
     type(LDT_paramEntry) :: penpt       ! Potential evaporation constant option (-)
     type(LDT_paramEntry) :: bareadj     ! Bare soil evaporation switch (Ek/Chen) (-)
     type(LDT_paramEntry) :: rdst        ! Constant allows tension water redist. option (-)
     type(LDT_paramEntry) :: czil        ! Zilitinkevich parameter controls the ratio of 
     type(LDT_paramEntry) :: rcmin       ! Minimal stomatal resistance (s/m)
     type(LDT_paramEntry) :: rcminclim   ! Climate dependent minimal stomatal resistance (s/m)
     type(LDT_paramEntry) :: rgl         ! Solar radiation threshold ( )
     type(LDT_paramEntry) :: hs          ! Vapor pressure resistance factor ( )
     type(LDT_paramEntry) :: maxlai      ! Max Leaf Area Index (-)
     type(LDT_paramEntry) :: z0          ! Roughness length (m)
     type(LDT_paramEntry) :: d50         ! The depth at which 50% roots are allocated (cm)
     type(LDT_paramEntry) :: sachtet_croot  ! SAC-HTET root depth distribution (m?)
     type(LDT_paramEntry) :: tbot        ! Bottom temperature (K)

   ! PET parameters
     type(LDT_paramEntry) :: pet         ! Climatology-based PET or PE (mm/day)
     type(LDT_paramEntry) :: petadj      ! Adjusted PET or PE (-)
     type(LDT_paramEntry) :: ksat_table  ! Saturated hydraulic conductivity read in 
                                         !  from table (eg, SAC-HTET)
     type(LDT_paramEntry) :: sand_table  ! Sand fraction read in from table (eg, SAC-HTET)
     type(LDT_paramEntry) :: clay_table  ! Clay fraction read in from table (eg, SAC-HTET)


  end type sachtet_type_dec

  type(sachtet_type_dec), allocatable :: SACHTET_struc(:)

contains

!BOP
! 
! !ROUTINE: SACHTETparms_init
! \label{SACHTETparms_init}
! 
! !INTERFACE:
  subroutine SACHTETparms_init(flag)

! !USES:
   use LDT_logMod,    only : LDT_verify, LDT_endrun, &
             LDT_getNextUnitNumber, LDT_releaseUnitNumber
!   use LDT_paramOptCheckMod, only: LDT_sachtetparmsOptChecks, &
!                      LDT_gridOptChecks
   use LDT_xmrg_reader
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the Sacramento Soil Moisture Accounting (SAC-SMA) 
! Heat transfer and evapotranspiration (HT-ET) parameter 
! datasets.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[sachtetParmssetup](\ref{sachtetParmssetup}) \newline
!    calls the registry to invoke the sachtetParms setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer  :: flag
   integer  :: n
   integer  :: c,r,m,k
   integer  :: numveg
   integer  :: ftn
   integer  :: rc
   real     :: temp
   real     :: xfactor
   integer  :: file_status
   logical  :: file_exists
   logical  :: sachtet_select 
   real, allocatable :: data2d(:,:)
   type(ESMF_Config)  :: rdhmconsts_table

   real :: clay_default(12), sand_default(12), satdk_default(12)
   data clay_default / 0.03,0.06,0.10,0.13,0.05,0.18,0.27,0.34,0.34,0.42,0.47,0.58 /
   data sand_default / 0.92,0.82,0.58,0.17,0.09,0.43,0.58,0.10,0.32,0.52,0.06,0.22 /
   data satdk_default / 0.0176, 0.0156, 0.347, 7.19, 5.56, 6.94, &
        6.31, 1.69, 2.44, 2.17, 1.03, 1.28 /
   ! SATDK values are actually 10^-6.

   real :: rcmin_default(14), rcminclim_default(14), rgl_default(14), lai_default(14)
   real :: hs_default(14), z0_default(14), d50_default(14), croot_default(14)
   data rcmin_default / 70.0, 50.0, 60.0, 70.0, 50.0, 40.0, 40.0, &   
        90.0, 150.0, 150.0, 150.0, 40.0, 100.0, 100.0 /
   data rcminclim_default / 70.0, 50.0, 60.0, 60.0, 50.0, 30.0, 30.0, &   
        20.0, 40.0, 40.0, 30.0, 40.0, 100.0, 100.0 /
   data rgl_default / 30.0, 30.0, 30.0, 30.0, 30.0, 65.0, 100.0, &   
        100.0,100.0,100.0,100.0,100.0,100.0, 30.0 /
   data hs_default / 41.69, 54.53, 51.93, 47.35, 47.35, 54.53, 36.35, &   
        42.0, 42.0, 42.0, 42.0, 36.35, 42.0, 51.75 /
   data lai_default / 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, &   
        5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 4.0 /
   data d50_default / 12.0, 21.0, 12.0, 23.0, 23.0, 23.0, 28.0, &   
        28.0, 27.0, 7.0, 16.0, 16.0, 5.0, 0.0 /
   data croot_default / -1.88, -1.84, -1.88, -1.76, -1.76, -1.76, -1.91, &
        -1.91, -2.05, -1.18, -1.45, -1.45, -1.45, -1.00 /

   character*3 :: months(12)
   character*3 :: sacmonths(12)
   data months /'jan','feb','mar','apr','may','jun','jul','aug',&
                'sep','oct','nov','dec'/
   data sacmonths /'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG',&
                   'SEP','OCT','NOV','DEC'/

!   character*50       :: sacsoilparms_option
   character*140      :: cosbysoils_table
   character*50       :: sachtetparms_proj
   integer            :: type_value
   character(20)      :: type_name 
   type(LDT_fillopts) :: sachtet
   type(LDT_fillopts) :: tbot
   type(LDT_fillopts) :: pet
   character*50       :: pet_proj

   real, allocatable  :: force_elev(:,:)

! _____________________________________________________________________

   allocate(SACHTET_struc(LDT_rc%nnest))

   sachtet_select = .false.
   do n=1,LDT_rc%nnest
      ! - SAC-HTET (v3.5.6) parameters:
      call set_param_attribs(SACHTET_struc(n)%sachtet356,"SACHTET356")
      
      call set_param_attribs(SACHTET_struc(n)%pet,"PET",&
            units="mm", &
            full_name="SAC-HTET Potential ET climatology")
      call set_param_attribs(SACHTET_struc(n)%petadj,"PETADJ",&
            units="-", &
            full_name="SAC-HTET Potential ET clim adjustment factor")
      call set_param_attribs(SACHTET_struc(n)%tbot,"TBOT",&
            units="K", &
            full_name="SAC-HTET LSM bottom temperature")

      if( SACHTET_struc(n)%sachtet356%selectOpt == 1 ) then
         sachtet_select = .true.
      endif
   enddo

   if( sachtet_select ) then
     write(LDT_logunit,*)" - - - - - - - - - - SAC-HTET LSM Parameters - - - - - - - - - - - - -"

   !- Load RDHM CONSTANTS input table file (filepath read-in from ldt.config file)
      call ESMF_ConfigFindLabel(LDT_config,"RDHM356 constants table:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%rdhmconsts_table,rc=rc)
         call LDT_verify(rc,'RDHM356 constants table: not specified')

         inquire(file=trim(SACHTET_struc(n)%rdhmconsts_table), exist=file_exists)
         if( .not. file_exists ) then
            write(LDT_logunit,*) "[ERR] RDHM Parameter Constants Table ",&
                 trim(SACHTET_struc(n)%rdhmconsts_table)," does not exist."
            call LDT_endrun
         endif
         write(LDT_logunit,*) "Reading in RDHM Parameter Constants Table Entries"
         rdhmconsts_table = ESMF_ConfigCreate(rc=rc)
         call ESMF_ConfigLoadFile(rdhmconsts_table, &
              trim(SACHTET_struc(n)%rdhmconsts_table), rc=rc)
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"RDHM356 universal undefined value:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%rdhm_undef,rc=rc)
         call LDT_verify(rc,'RDHM356 universal undefined value: not specified')
      enddo

      ! -----
      
      do n = 1, LDT_rc%nnest

      !- Soil parameters - 

         ! Fill in derived parameter entries:
         ! ( input_parmattribs -> output_parmattribs ) 
         call populate_param_attribs( "RDHM356_LZFPM", &
              "Lower zone primary free water max. storage","mm", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%lzfpm )

         allocate(SACHTET_struc(n)%lzfpm%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%lzfpm%vlevels))       

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_LZFPM=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_LZFPM =", temp
         SACHTET_struc(n)%lzfpm%value = temp

         call populate_param_attribs( "RDHM356_LZFSM", &
              "Lower zone supplemental free water max. storage","mm", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%lzfsm )

         allocate(SACHTET_struc(n)%lzfsm%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%lzfsm%vlevels))

         write(LDT_logunit,'(A20,1X,f12.5)') " sac_LZFSM =", temp
         SACHTET_struc(n)%lzfsm%value = temp

         call populate_param_attribs( "RDHM356_LZPK", &
              "Lower zone primary free water depletion rate","day-1", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%lzpk )

         allocate(SACHTET_struc(n)%lzpk%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%lzpk%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_LZPK=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_LZPK =", temp
         SACHTET_struc(n)%lzpk%value = temp

         call populate_param_attribs( "RDHM356_LZSK", &
              "Lower zone supplemental free water depletion rate","day-1", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%lzsk )

         allocate(SACHTET_struc(n)%lzsk%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%lzsk%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_LZSK=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_LZSK =", temp
         SACHTET_struc(n)%lzsk%value = temp

         call populate_param_attribs( "RDHM356_LZTWM", &
              "Lower zone tension water max storage","mm", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%lztwm )

         allocate(SACHTET_struc(n)%lztwm%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%lztwm%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_LZTWM=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_LZTWM =", temp
         SACHTET_struc(n)%lztwm%value = temp

         call populate_param_attribs( "RDHM356_PFREE", &
              "Fraction percolation free water storage","day-1", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%pfree )

         allocate(SACHTET_struc(n)%pfree%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%pfree%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_PFREE=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_PFREE =", temp
         SACHTET_struc(n)%pfree%value = temp

         call populate_param_attribs( "RDHM356_REXP", &
              "Exponent of percolation equation","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%rexp )

         allocate(SACHTET_struc(n)%rexp%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%rexp%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_REXP=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_REXP =", temp
         SACHTET_struc(n)%rexp%value = temp

         call populate_param_attribs( "RDHM356_UZFWM", &
              "Upper zone free water max. storage","mm", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%uzfwm )

         allocate(SACHTET_struc(n)%uzfwm%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%uzfwm%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_UZFWM=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_UZFWM =", temp
         SACHTET_struc(n)%uzfwm%value = temp

         call populate_param_attribs( "RDHM356_UZTWM", &
              "Upper zone tension water max storage","mm", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%uztwm )

         allocate(SACHTET_struc(n)%uztwm%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%uztwm%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_UZTWM=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_UZTWM =", temp
         SACHTET_struc(n)%uztwm%value = temp

         call populate_param_attribs( "RDHM356_UZK", &
              "Upper zone free water latent depletion rate","day-1", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%uzk )

         allocate(SACHTET_struc(n)%uzk%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%uzk%vlevels))

         write(LDT_logunit,'(A20,1X,f12.5)') " sac_UZK =", temp
         SACHTET_struc(n)%uzk%value = temp

         call populate_param_attribs( "RDHM356_ZPERC", &
              "Max percolation rate","mm", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%zperc )

         allocate(SACHTET_struc(n)%zperc%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%zperc%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_ZPERC=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_ZPERC =", temp
         SACHTET_struc(n)%zperc%value = temp

         call populate_param_attribs( "RDHM356_DZUP", &
              "Depth of upper zone","mm", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%dzup )

         allocate(SACHTET_struc(n)%dzup%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%dzup%vlevels))

         call populate_param_attribs( "RDHM356_DZLW", &
              "Depth of lower zone","mm", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%dzlw )

         allocate(SACHTET_struc(n)%dzlw%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%dzlw%vlevels))

         call populate_param_attribs( "RDHM356_TXTLW", &
              "Dominant texture of lower zone","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%txtlw )

         allocate(SACHTET_struc(n)%txtlw%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%txtlw%vlevels))

         call populate_param_attribs( "RDHM356_TXTOT", &
              "Dominant texture of total column","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%txtot )

         allocate(SACHTET_struc(n)%txtot%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%txtot%vlevels))





! - End of soil parameters

         call populate_param_attribs( "RDHM356_EFC", &
              "Forest fraction","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%forestfrac )

         allocate(SACHTET_struc(n)%forestfrac%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%forestfrac%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_EFC=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_EFC =", temp
         SACHTET_struc(n)%forestfrac%value = temp

         call populate_param_attribs( "RDHM356_PCTIM", &
              "Impervious fraction of the watershed area","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%pctim )

         allocate(SACHTET_struc(n)%pctim%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%pctim%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_PCTIM=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_PCTIM =", temp
         SACHTET_struc(n)%pctim%value = temp

         call populate_param_attribs( "RDHM356_ADIMP", &
              "Additional impervious area","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%adimp )

         allocate(SACHTET_struc(n)%adimp%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%adimp%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_ADIMP=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_ADIMP =", temp
         SACHTET_struc(n)%adimp%value = temp

         call populate_param_attribs( "RDHM356_RIVA", &
              "Riparian vegetation area","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%riva )

         allocate(SACHTET_struc(n)%riva%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%riva%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_RIVA=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_RIVA =", temp
         SACHTET_struc(n)%riva%value = temp

         call populate_param_attribs( "RDHM356_RSERV", &
              "LZ free water fraction not transferable to LZ tension water","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%rserv )

         allocate(SACHTET_struc(n)%rserv%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%rserv%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_RSERV=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_RSERV =", temp
         SACHTET_struc(n)%rserv%value = temp

         call populate_param_attribs( "RDHM356_SIDE", &
              "Ratio of deep recharge to channel base flow","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%side )

         allocate(SACHTET_struc(n)%side%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%side%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_SIDE=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_SIDE =", temp
         SACHTET_struc(n)%side%value = temp

         call populate_param_attribs( "RDHM356_SOILALB", &
              "Soil (snow-free) albedo","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%soilalb )

         allocate(SACHTET_struc(n)%soilalb%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%soilalb%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_SOILALB=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_SOILALB =", temp
         SACHTET_struc(n)%soilalb%value = temp

      !- Frozen soil physics:
         call populate_param_attribs( "RDHM356_TIMEOFFSET", &
              "Time offset (UTC->local time conversion)","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%offsetTime )

         allocate(SACHTET_struc(n)%offsetTime%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%offsetTime%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='sac_TIMEOFFSET=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " sac_TIMEOFFSET =", temp
         SACHTET_struc(n)%offsetTime%value = temp

         call populate_param_attribs( "RDHM356_SOILTEXT", &
              "Soil texture (frz)","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%frz_soiltext )

         allocate(SACHTET_struc(n)%frz_soiltext%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%frz_soiltext%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='frz_STXT=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " frz_STXT =", temp
         SACHTET_struc(n)%frz_soiltext%value = temp

         call populate_param_attribs( "RDHM356_CKSL", &
              "Ratio of frozen to non-frozen surface (frz)","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%frz_cksl )

         allocate(SACHTET_struc(n)%frz_cksl%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%frz_cksl%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='frz_CKSL=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " frz_CKSL =", temp
         SACHTET_struc(n)%frz_cksl%value = temp

         call populate_param_attribs( "RDHM356_RSMAX", &
              "Residual porosity (frz)","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%frz_rsmax )

         allocate(SACHTET_struc(n)%frz_rsmax%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%frz_rsmax%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='frz_RSMAX=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " frz_RSMAX =", temp
         SACHTET_struc(n)%frz_rsmax%value = temp

         call populate_param_attribs( "RDHM356_ZBOT", &
              "Lower boundary depth (frz)","m", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%frz_zbot )

         allocate(SACHTET_struc(n)%frz_zbot%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%frz_zbot%vlevels))

         call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='frz_ZBOT=',rc=rc)
         if( rc .ne. 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " frz_ZBOT =", temp
         SACHTET_struc(n)%frz_zbot%value = temp

         !- Soil parameter table entries
         call populate_param_attribs( "RDHM356_SATDK", &
              "Saturated hydraulic conductivity","m s-1", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%ksat_table )

         allocate(SACHTET_struc(n)%ksat_table%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%ksat_table%vlevels))

         call populate_param_attribs( "RDHM356_SAND", &
              "Sand fraction","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%sand_table )

         allocate(SACHTET_struc(n)%sand_table%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%sand_table%vlevels))

         call populate_param_attribs( "RDHM356_CLAY", &
              "Clay fraction","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%clay_table )

         allocate(SACHTET_struc(n)%clay_table%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%clay_table%vlevels))

      !- Vegetation parameter table entries

         call populate_param_attribs( "RDHM356_RCMIN", &
              "Minimal stomatal resistance","s/m", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%rcmin )

         allocate(SACHTET_struc(n)%rcmin%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
                                !           LDT_LSMparam_struc(n)%landcover%vlevels))
              SACHTET_struc(n)%rcmin%vlevels))

         call populate_param_attribs( "RDHM356_RCMINCLIM", &
              "Climate dependent minimal stomatal resistance","s/m", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%rcminclim )

         allocate(SACHTET_struc(n)%rcminclim%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
                                !           LDT_LSMparam_struc(n)%landcover%vlevels))
              SACHTET_struc(n)%rcminclim%vlevels))

         call populate_param_attribs( "RDHM356_RGL", &
              "Solar radiation threshold","??", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%rgl )

         allocate(SACHTET_struc(n)%rgl%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
                                !           LDT_LSMparam_struc(n)%landcover%vlevels))
              SACHTET_struc(n)%rgl%vlevels))

         call populate_param_attribs( "RDHM356_HS", &
              "Vapor pressure resistance factor","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%hs )

         allocate(SACHTET_struc(n)%hs%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
                                !           LDT_LSMparam_struc(n)%landcover%vlevels))
              SACHTET_struc(n)%hs%vlevels))

         call populate_param_attribs( "RDHM356_MAXLAI", &
              "Maximum LAI","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%maxlai )

         allocate(SACHTET_struc(n)%maxlai%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
                                !           LDT_LSMparam_struc(n)%landcover%vlevels))
              SACHTET_struc(n)%maxlai%vlevels))

         call populate_param_attribs( "RDHM356_Z0", &
              "Roughness length","m", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%z0 )

         allocate(SACHTET_struc(n)%z0%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
                                !           LDT_LSMparam_struc(n)%landcover%vlevels))
              SACHTET_struc(n)%z0%vlevels))

         call populate_param_attribs( "RDHM356_D50", &
              "The depth at which 50% roots are allocated","cm", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%d50 )

         allocate(SACHTET_struc(n)%d50%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),& 
                                !           LDT_LSMparam_struc(n)%landcover%vlevels))
              SACHTET_struc(n)%d50%vlevels))

         call populate_param_attribs( "RDHM356_CROOT", &
              "SAC-HTET root depth distribution","m", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%sachtet_croot )

         allocate(SACHTET_struc(n)%sachtet_croot%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
                                !           LDT_LSMparam_struc(n)%landcover%vlevels))
              SACHTET_struc(n)%sachtet_croot%vlevels))

         call populate_param_attribs( "RDHM356_CZIL", &
              "Zilitinkevich parameter control","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%czil )

         allocate(SACHTET_struc(n)%czil%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%czil%vlevels))

         call populate_param_attribs( "RDHM356_FXEXP", &
              "Bare soil evaporation exponent","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%fxexp )

         allocate(SACHTET_struc(n)%fxexp%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),& 
              SACHTET_struc(n)%fxexp%vlevels))

         call populate_param_attribs( "RDHM356_RCMAX", &
              "Maximum stomatal resistance","s/m", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%rcmax )

         allocate(SACHTET_struc(n)%rcmax%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%rcmax%vlevels))

         call populate_param_attribs( "RDHM356_TOPT", &
              "Optimum air temperature for transpiration","K", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%topt )

         allocate(SACHTET_struc(n)%topt%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%topt%vlevels))

         call populate_param_attribs( "RDHM356_PLANTCOEF", &
              "Plant coefficient","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%plantcoef )

         allocate(SACHTET_struc(n)%plantcoef%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%plantcoef%vlevels))

         call populate_param_attribs( "RDHM356_PENPT", &
              "Potential evaporation constant option","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%penpt )

         allocate(SACHTET_struc(n)%penpt%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%penpt%vlevels))

         call populate_param_attribs( "RDHM356_BAREADJ", &
              "Bare soil evaporation switch (Ek/Chen)","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%bareadj )

         allocate(SACHTET_struc(n)%bareadj%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%bareadj%vlevels))

         call populate_param_attribs( "RDHM356_RDST", &
              "Constant allows tension water redist. option","-", &
              SACHTET_struc(n)%sachtet356, &
              SACHTET_struc(n)%rdst )

         allocate(SACHTET_struc(n)%rdst%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SACHTET_struc(n)%rdst%vlevels))

      end do


   !- Generate SAC-HTET soil parameters within LDT:
      if( LDT_rc%create_soilparms_option == "create" ) then

        write(LDT_logunit,*) " -- Calculating the SAC model soil pararmeters -- "

     !- Read in ldt.config file entries:
        SACHTET_struc(:)%csoilparm_method = "none"
        call ESMF_ConfigFindLabel(LDT_config,"SACHTET soil parameter method:",rc=rc)
        do n=1,LDT_rc%nnest
           call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%csoilparm_method,rc=rc)
           call LDT_verify(rc,'SACHTET soil parameter method: not specified')
        enddo

        call ESMF_ConfigGetAttribute(LDT_config, cosbysoils_table, &
             label="SACHTET Cosby soil parameter table:",rc=rc)
        call LDT_verify(rc,"SACHTET Cosby soil parameter table: option not specified in the config file")

      ! Check if HSG or Soil texture fields are present:
        do n=1,LDT_rc%nnest
          if( LDT_soils_struc(n)%hsg%selectOpt .ne. 1 ) then
            print *, " No HSG dataset has been selected ... needed for SAC soil parameter"
            call LDT_endrun
          endif
          if( LDT_LSMparam_struc(n)%texture%selectOpt .ne. 1 ) then
            print *, " No Soil texture dataset has been selected ... needed for SAC soil parameter"
            call LDT_endrun
          endif
 
        ! Select SAC soil parameter derivation method: 
          select case( SACHTET_struc(n)%csoilparm_method )

            case( "Koren_v1" )

              call create_SACSoilparms_Korenv1(n, &
                      LDT_soils_struc(n)%texture_nlyrs, &
                      LDT_soils_struc(n)%hsg, &
                      cosbysoils_table, &
                      SACHTET_struc(n)%uzfwm, SACHTET_struc(n)%uztwm, &
                      SACHTET_struc(n)%lztwm, SACHTET_struc(n)%lzfsm, &
                      SACHTET_struc(n)%lzfpm, SACHTET_struc(n)%uzk,   &
                      SACHTET_struc(n)%lzsk,  SACHTET_struc(n)%lzpk,  &
                      SACHTET_struc(n)%zperc, SACHTET_struc(n)%rexp,  &
                      SACHTET_struc(n)%pfree, SACHTET_struc(n)%frz_soiltext, &
                      SACHTET_struc(n)%dzup,  SACHTET_struc(n)%dzlw,  &
                      SACHTET_struc(n)%txtlw, SACHTET_struc(n)%txtot ) 

            case default
              write(LDT_logunit,*) "[ERR] No other SAC-HTET soil parameter "
              write(LDT_logunit,*) " derivation method is available at this time."
              write(LDT_logunit,*) " Either select: 'Koren_v1' or read in parameters."
              call LDT_endrun
          end select

        ! 1) Lower zone primary free water (slow) max. storage (mm)
          print *, "LZFPM: ", SACHTET_struc(n)%lzfpm%value(1000,1000,1)

        ! 2) Lower zone supplemental free water (fast) max. storage (mm)
          print *, "LZFSM: ", SACHTET_struc(n)%lzfsm%value(1000,1000,1)
 
        ! 3) Lower zone primary free water depletion rate (day^-1)
          print *, "LZPK: ", SACHTET_struc(n)%lzpk%value(1000,1000,1)

        ! 4) Lower zone supplemental free water depletion rate (day^-1)
          print *, "LZSK: ", SACHTET_struc(n)%lzsk%value(1000,1000,1)

        ! 5) Lower zone tension water max. storage (mm)
          print *, "LZTWM: ", SACHTET_struc(n)%lztwm%value(1000,1000,1)
 
        ! 6) Fraction percolation from upper to lower free water storage (day^-1)
          print *, "PFREE: ", SACHTET_struc(n)%pfree%value(1000,1000,1)

        ! 7) Expon. of the percolation eqtn.
          print *, "REXP: ", SACHTET_struc(n)%rexp%value(1000,1000,1)

        ! 8) Upper zone free water max. storage (mm)
          print *, "UZFWM: ", SACHTET_struc(n)%uzfwm%value(1000,1000,1)

        ! 9) Upper zone tension water max. storage (mm)
          print *, "UZTWM: ", SACHTET_struc(n)%uztwm%value(1000,1000,1)

        ! 10) Upper zone free water latent depletion rate (day^-1)
          print *, "UZK: ", SACHTET_struc(n)%uzk%value(1000,1000,1)

        ! 11) Max. percolation rate (mm)
          print *, "ZPERC: ", SACHTET_struc(n)%zperc%value(1000,1000,1)

        end do


   !- Read in SAC-HTET soil parameters from a-priori gridded maps:
      elseif( LDT_rc%create_soilparms_option == "readin" ) then

         write(LDT_logunit,*) " == Reading in the SAC model soil pararmeters == "

         SACHTET_struc(:)%lzfpm_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET LZFPM map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%lzfpm_file,rc=rc)
            call LDT_verify(rc,'SACHTET LZFPM map: not specified')
         enddo

         SACHTET_struc(:)%lzfsm_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET LZFSM map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%lzfsm_file,rc=rc)
            call LDT_verify(rc,'SACHTET LZFSM map: not specified')
         enddo

         SACHTET_struc(:)%lzpk_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET LZPK map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%lzpk_file,rc=rc)
            call LDT_verify(rc,'SACHTET LZPK map: not specified')
         enddo

         SACHTET_struc(:)%lzsk_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET LZSK map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%lzsk_file,rc=rc)
            call LDT_verify(rc,'SACHTET LZSK map: not specified')
         enddo

         SACHTET_struc(:)%lztwm_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET LZTWM map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%lztwm_file,rc=rc)
            call LDT_verify(rc,'SACHTET LZTWM map: not specified')
         enddo

         SACHTET_struc(:)%uzfwm_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET UZFWM map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%uzfwm_file,rc=rc)
            call LDT_verify(rc,'SACHTET UZFWM map: not specified')
         enddo

         SACHTET_struc(:)%uztwm_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET UZTWM map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%uztwm_file,rc=rc)
            call LDT_verify(rc,'SACHTET UZTWM map: not specified')
         enddo

         SACHTET_struc(:)%uzk_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET UZK map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%uzk_file,rc=rc)
            call LDT_verify(rc,'SACHTET UZK map: not specified')
         enddo

         SACHTET_struc(:)%pfree_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET PFREE map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%pfree_file,rc=rc)
            call LDT_verify(rc,'SACHTET PFREE map: not specified')
         enddo

         SACHTET_struc(:)%rexp_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET REXP map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%rexp_file,rc=rc)
            call LDT_verify(rc,'SACHTET REXP map: not specified')
         enddo

         SACHTET_struc(:)%zperc_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET ZPERC map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%zperc_file,rc=rc)
            call LDT_verify(rc,'SACHTET ZPERC map: not specified')
         enddo

         SACHTET_struc(:)%stxt_file = "none"
         call ESMF_ConfigFindLabel(LDT_config,"SACHTET STXT map:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%stxt_file,rc=rc)
            call LDT_verify(rc,'SACHTET STXT map: not specified')
         enddo

      else
         write(*,*) " That is not a valid SACHTET soil parameter option input ... "
         stop
      endif

      SACHTET_struc(:)%efc_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET EFC map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%efc_file,rc=rc)
         call LDT_verify(rc,'SACHTET EFC map: not specified')
      enddo

      SACHTET_struc(:)%pctim_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET PCTIM map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%pctim_file,rc=rc)
         call LDT_verify(rc,'SACHTET PCTIM map: not specified')
      enddo

      SACHTET_struc(:)%adimp_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET ADIMP map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%adimp_file,rc=rc)
         call LDT_verify(rc,'SACHTET ADIMP map: not specified')
      enddo

      SACHTET_struc(:)%side_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET SIDE map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%side_file,rc=rc)
         call LDT_verify(rc,'SACHTET SIDE map: not specified')
      enddo

      SACHTET_struc(:)%riva_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET RIVA map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%riva_file,rc=rc)
         call LDT_verify(rc,'SACHTET RIVA map: not specified')
      enddo

      SACHTET_struc(:)%rserv_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET RSERV map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%rserv_file,rc=rc)
         call LDT_verify(rc,'SACHTET RSERV map: not specified')
      enddo

      SACHTET_struc(:)%soilalb_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET soil albedo map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%soilalb_file,rc=rc)
         call LDT_verify(rc,'SACHTET soil albedo map: not specified')
      enddo

      SACHTET_struc(:)%timeOffset_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET offset time map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%timeOffset_file,rc=rc)
         call LDT_verify(rc,'SACHTET offset time map: not specified')
      enddo

    ! Frozen soil physics parameters ("frz"):
!      SACHTET_struc(:)%stxt_file = "none"
!      call ESMF_ConfigFindLabel(LDT_config,"SACHTET STXT map:",rc=rc)
!      do n=1,LDT_rc%nnest
!         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%stxt_file,rc=rc)
!         call LDT_verify(rc,'SACHTET STXT map: not specified')
!      enddo

      SACHTET_struc(:)%tbot_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET TBOT map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%tbot_file,rc=rc)
         call LDT_verify(rc,'SACHTET TBOT map: not specified')
      enddo

    ! Read in Tbot topographic correction entries:
      SACHTET_struc(:)%tbot_topocorr = "none"
      call ESMF_ConfigFindLabel(LDT_config,"Bottom temperature topographic downscaling:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%tbot_topocorr,rc=rc)
         call LDT_verify(rc,'Bottom temperature topographic downscaling: not specified')
       ! Allow for mis-entered lapse-rate option entry:
         if( SACHTET_struc(n)%tbot_topocorr == "lapse_rate" ) SACHTET_struc(n)%tbot_topocorr="lapse-rate"
         if( SACHTET_struc(n)%tbot_topocorr == "lapse rate" ) SACHTET_struc(n)%tbot_topocorr="lapse-rate"
         if( SACHTET_struc(n)%tbot_topocorr == "Lapse-rate" ) SACHTET_struc(n)%tbot_topocorr="lapse-rate"
      enddo

      SACHTET_struc(:)%cksl_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET CKSL map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%cksl_file,rc=rc)
         call LDT_verify(rc,'SACHTET CKSL map: not specified')
      enddo

      SACHTET_struc(:)%rsmax_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET RSMAX map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%rsmax_file,rc=rc)
         call LDT_verify(rc,'SACHTET RSMAX map: not specified')
      enddo

      SACHTET_struc(:)%zbot_file = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET ZBOT map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%zbot_file,rc=rc)
         call LDT_verify(rc,'SACHTET ZBOT map: not specified')
      enddo

   !- Parameter lookup tables:
      SACHTET_struc(:)%sacsoilparms_table = "none" 
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET soiltype parameter table:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%sacsoilparms_table,rc=rc)
         call LDT_verify(rc,'SACHTET soiltype parameter table: not specified')
      enddo

      SACHTET_struc(:)%sacvegparms_table = "none"
      call ESMF_ConfigFindLabel(LDT_config,"SACHTET vegetation parameter table:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%sacvegparms_table,rc=rc)
         call LDT_verify(rc,'SACHTET vegetation parameter table: not specified')
      enddo


    ! --  Grid info and inputs:

      call ESMF_ConfigFindLabel(LDT_config,"SACHTET parameter spatial transform:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%sachtetparms_gridtransform,&
              rc=rc)
         call LDT_verify(rc,'SACHTET parameter spatial transform: option not specified in the config file')
      enddo

      sachtet%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, sachtet%filltype, &
           label="SACHTET parameter fill option:",rc=rc)
      call LDT_verify(rc,"SACHTET parameter fill option: option not specified in the config file")

      if( sachtet%filltype == "neighbor" .or. sachtet%filltype == "average" ) then
         call ESMF_ConfigGetAttribute(LDT_config, sachtet%fillvalue, &
              label="SACHTET parameter fill value:",rc=rc)
         call LDT_verify(rc,"SACHTET parameter fill value: option not specified in the config file")

         call ESMF_ConfigGetAttribute(LDT_config, sachtet%fillradius, &
              label="SACHTET parameter fill radius:",rc=rc)
         call LDT_verify(rc,"SACHTET parameter fill radius: option not specified in the config file")
      elseif( sachtet%filltype == "none" ) then
         write(LDT_logunit,*) " -- 'NONE' Parameter-Mask Agreement Option Selected for SAC-HTET parameters"
      else
         write(LDT_logunit,*) "[ERR] Fill option for SAC-HTET is not valid: ",trim(sachtet%filltype)
         write(LDT_logunit,*) " Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) " Programming stopping ..."
         call LDT_endrun
      end if

      call ESMF_ConfigGetAttribute(LDT_config,sachtetparms_proj,&
           label="SACHTET map projection:",rc=rc)
      call LDT_verify(rc,'SACHTET map projection: option not specified in the config file')
      SACHTET_struc(:)%sachtetparms_proj = sachtetparms_proj

      !- Loop over all domain nests:
      do n = 1, LDT_rc%nnest

      !- Read-in Soil Parameters:
         if( LDT_rc%create_soilparms_option == "readin" ) then
            
            if( SACHTET_struc(n)%lzfpm%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%lzfpm%value(1,1,1))
               write(LDT_logunit,*) "Reading LZFPM file: "//trim(SACHTET_struc(n)%lzfpm_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%lzfpm_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%lzfpm%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%lzfpm_file)
               SACHTET_struc(n)%lzfpm%value(:,:,1) = &
                    SACHTET_struc(n)%lzfpm%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%lzfsm%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%lzfsm%value(1,1,1))
               write(LDT_logunit,*) "Reading LZFSM file: "//trim(SACHTET_struc(n)%lzfsm_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%lzfsm_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%lzfsm%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%lzfsm_file)
               SACHTET_struc(n)%lzfsm%value(:,:,1) = &
                    SACHTET_struc(n)%lzfsm%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%lztwm%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%lztwm%value(1,1,1))
               write(LDT_logunit,*) "Reading LZTWM file: "//trim(SACHTET_struc(n)%lztwm_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%lztwm_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%lztwm%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%lztwm_file)
               SACHTET_struc(n)%lztwm%value(:,:,1) = &
                    SACHTET_struc(n)%lztwm%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%lzpk%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%lzpk%value(1,1,1))
               write(LDT_logunit,*) "Reading LZPK file: "//trim(SACHTET_struc(n)%lzpk_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),    &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%lzpk_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%lzpk%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%lzpk_file)
               SACHTET_struc(n)%lzpk%value(:,:,1) = &
                    SACHTET_struc(n)%lzpk%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%lzsk%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%lzsk%value(1,1,1))
               write(LDT_logunit,*) "Reading LZSK file: "//trim(SACHTET_struc(n)%lzsk_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),    &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%lzsk_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%lzsk%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%lzsk_file)
               SACHTET_struc(n)%lzsk%value(:,:,1) = &
                    SACHTET_struc(n)%lzsk%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%uzfwm%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%uzfwm%value(1,1,1))
               write(LDT_logunit,*) "Reading UZFWM file: "//trim(SACHTET_struc(n)%uzfwm_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%uzfwm_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%uzfwm%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%uzfwm_file)
               SACHTET_struc(n)%uzfwm%value(:,:,1) = &
                    SACHTET_struc(n)%uzfwm%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%uztwm%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%uztwm%value(1,1,1))
               write(LDT_logunit,*) "Reading UZTWM file: "//trim(SACHTET_struc(n)%uztwm_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%uztwm_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%uztwm%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%uztwm_file)
               SACHTET_struc(n)%uztwm%value(:,:,1) = &
                    SACHTET_struc(n)%uztwm%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%uzk%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%uzk%value(1,1,1))
               write(LDT_logunit,*) "Reading UZK file: "//trim(SACHTET_struc(n)%uzk_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),   &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%uzk_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%uzk%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%uzk_file)
               SACHTET_struc(n)%uzk%value(:,:,1) = &
                    SACHTET_struc(n)%uzk%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%pfree%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%pfree%value(1,1,1))
               write(LDT_logunit,*) "Reading PFREE file: "//trim(SACHTET_struc(n)%pfree_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%pfree_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%pfree%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%pfree_file)
               SACHTET_struc(n)%pfree%value(:,:,1) = &
                    SACHTET_struc(n)%pfree%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%rexp%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%rexp%value(1,1,1))
               write(LDT_logunit,*) "Reading REXP file: "//trim(SACHTET_struc(n)%rexp_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%rexp_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%rexp%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%rexp_file)
               SACHTET_struc(n)%rexp%value(:,:,1) = &
                    SACHTET_struc(n)%rexp%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%zperc%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%zperc%value(1,1,1))
               write(LDT_logunit,*) "Reading ZPERC file: "//trim(SACHTET_struc(n)%zperc_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%zperc_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%zperc%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%zperc_file)
               SACHTET_struc(n)%zperc%value(:,:,1) = &
                    SACHTET_struc(n)%zperc%value(:,:,1) * xfactor
            endif

            if( SACHTET_struc(n)%frz_soiltext%value(1,1,1) < 0. ) then
               xfactor = abs(SACHTET_struc(n)%frz_soiltext%value(1,1,1))
               write(LDT_logunit,*) "Reading soil texture file: "//trim(SACHTET_struc(n)%stxt_file)
               call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),    &
                    LDT_rc%gridDesc(n,:), SACHTET_struc(n)%stxt_file, LDT_rc%udef,  &
                    SACHTET_struc(n)%frz_soiltext%value(:,:,1) )
               write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%stxt_file)
               SACHTET_struc(n)%frz_soiltext%value(:,:,1) = &
                       SACHTET_struc(n)%frz_soiltext%value(:,:,1) * xfactor
            endif

         end if   ! Done reading in soil parameters

      !- Other SAC-HTET parameters:

         if( SACHTET_struc(n)%forestfrac%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%forestfrac%value(1,1,1))
            write(LDT_logunit,*) "Reading Forest frac file: "//trim(SACHTET_struc(n)%efc_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),   &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%efc_file, LDT_rc%udef,  &
                 SACHTET_struc(n)%forestfrac%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%efc_file)
            SACHTET_struc(n)%forestfrac%value(:,:,1) = &
                 SACHTET_struc(n)%forestfrac%value(:,:,1) * xfactor
         endif

         if( SACHTET_struc(n)%soilalb%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%soilalb%value(1,1,1))
            write(LDT_logunit,*) "Reading soil albedo file: "//trim(SACHTET_struc(n)%soilalb_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),      &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%soilalb_file, LDT_rc%udef, &
                 SACHTET_struc(n)%soilalb%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%soilalb_file)
            SACHTET_struc(n)%soilalb%value(:,:,1) = &
                 SACHTET_struc(n)%soilalb%value(:,:,1) * xfactor
         endif

         if( SACHTET_struc(n)%pctim%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%pctim%value(1,1,1))
            write(LDT_logunit,*) "Reading PCTIM file: "//trim(SACHTET_struc(n)%pctim_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%pctim_file, LDT_rc%udef,  &
                 SACHTET_struc(n)%pctim%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%pctim_file)
            SACHTET_struc(n)%pctim%value(:,:,1) = &
                 SACHTET_struc(n)%pctim%value(:,:,1) * xfactor
         endif

         if( SACHTET_struc(n)%adimp%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%adimp%value(1,1,1))
            write(LDT_logunit,*) "Reading ADIMP file: "//trim(SACHTET_struc(n)%adimp_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%adimp_file, LDT_rc%udef,  &
                 SACHTET_struc(n)%adimp%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%adimp_file)
            SACHTET_struc(n)%adimp%value(:,:,1) = &
                 SACHTET_struc(n)%adimp%value(:,:,1) * xfactor
         endif

         if( SACHTET_struc(n)%riva%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%riva%value(1,1,1))
            write(LDT_logunit,*) "Reading RIVA file: "//trim(SACHTET_struc(n)%riva_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),    &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%riva_file, LDT_rc%udef,  &
                 SACHTET_struc(n)%riva%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%riva_file)
            SACHTET_struc(n)%riva%value(:,:,1) = &
                 SACHTET_struc(n)%riva%value(:,:,1) * xfactor
         endif

         if( SACHTET_struc(n)%rserv%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%rserv%value(1,1,1))
            write(LDT_logunit,*) "Reading RSERV file: "//trim(SACHTET_struc(n)%rserv_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%rserv_file, LDT_rc%udef,  &
                 SACHTET_struc(n)%rserv%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%rserv_file)
            SACHTET_struc(n)%rserv%value(:,:,1) = &
                 SACHTET_struc(n)%rserv%value(:,:,1) * xfactor
         endif

         if( SACHTET_struc(n)%side%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%side%value(1,1,1))
            write(LDT_logunit,*) "Reading SIDE file: "//trim(SACHTET_struc(n)%side_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),    &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%side_file, LDT_rc%udef,  &
                 SACHTET_struc(n)%side%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%side_file)
            SACHTET_struc(n)%side%value(:,:,1) = &
                 SACHTET_struc(n)%side%value(:,:,1) * xfactor
         endif

         write(LDT_logunit,*) " -- Reading in SAC-HTET FRZ Soil Parameters -- "

         if( SACHTET_struc(n)%frz_zbot%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%frz_zbot%value(1,1,1))
            write(LDT_logunit,*) "Reading ZBOT file: "//trim(SACHTET_struc(n)%zbot_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),    &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%zbot_file, LDT_rc%udef,  &
                 SACHTET_struc(n)%frz_zbot%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%zbot_file)
            SACHTET_struc(n)%frz_zbot%value(:,:,1) = &
                 SACHTET_struc(n)%frz_zbot%value(:,:,1) * xfactor
         endif

         if( SACHTET_struc(n)%offsetTime%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%offsetTime%value(1,1,1))
            write(LDT_logunit,*) "Reading Time offset file: "//trim(SACHTET_struc(n)%timeoffset_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),      &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%timeOffset_file, LDT_rc%udef, &
                 SACHTET_struc(n)%offsetTime%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%timeoffset_file)
            SACHTET_struc(n)%offsetTime%value(:,:,1) = &
                 SACHTET_struc(n)%offsetTime%value(:,:,1) * xfactor
         endif

!         if( SACHTET_struc(n)%frz_soiltext%value(1,1,1) < 0. ) then
!            xfactor = abs(SACHTET_struc(n)%frz_soiltext%value(1,1,1))
!            write(LDT_logunit,*) "Reading soil texture file: "//trim(SACHTET_struc(n)%stxt_file)
!            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),    &
!                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%stxt_file, LDT_rc%udef,  &
!                 SACHTET_struc(n)%frz_soiltext%value(:,:,1) )
!            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%stxt_file)
!            SACHTET_struc(n)%frz_soiltext%value(:,:,1) = &
!                 SACHTET_struc(n)%frz_soiltext%value(:,:,1) * xfactor
!         endif

         if( SACHTET_struc(n)%frz_cksl%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%frz_cksl%value(1,1,1))
            write(LDT_logunit,*) "Reading CKSL file: "//trim(SACHTET_struc(n)%cksl_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),    &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%cksl_file, LDT_rc%udef,  &
                 SACHTET_struc(n)%frz_cksl%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%cksl_file)
            SACHTET_struc(n)%frz_cksl%value(:,:,1) = &
                 SACHTET_struc(n)%frz_cksl%value(:,:,1) * xfactor
         endif

         if( SACHTET_struc(n)%frz_rsmax%value(1,1,1) < 0. ) then
            xfactor = abs(SACHTET_struc(n)%frz_rsmax%value(1,1,1))
            write(LDT_logunit,*) "Reading RSMAX file: "//trim(SACHTET_struc(n)%rsmax_file)
            call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                 LDT_rc%gridDesc(n,:), SACHTET_struc(n)%rsmax_file, LDT_rc%udef,  &
                 SACHTET_struc(n)%frz_rsmax%value(:,:,1) )
            write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%rsmax_file)
            SACHTET_struc(n)%frz_rsmax%value(:,:,1) = &
                 SACHTET_struc(n)%frz_rsmax%value(:,:,1) * xfactor
         endif

         !- Read in soil parameter table values:
         if( SACHTET_struc(n)%sacsoilparms_table == "none" .or. &
              SACHTET_struc(n)%stxt_file == "none" ) then
            write(LDT_logunit,*) "[INFO] SAC-HTET -- No soil parameters read in "
            write(LDT_logunit,*) "[INFO] SAC-HTET -- from look-up table."
            write(LDT_logunit,*) "[INFO] If you want SAC/FRZ soil texture and"
            write(LDT_logunit,*) "[INFO] read-in look-up soil parameters, you need"
            write(LDT_logunit,*) "[INFO] to enter full path and filenames. "
         endif
         if( SACHTET_struc(n)%sacsoilparms_table .ne. "none" .and. &
              SACHTET_struc(n)%stxt_file .ne. "none" ) then
            write(LDT_logunit,*) "[INFO] SAC-HTET -- Reading in soils parameter from table."
            ftn = LDT_getNextUnitNumber()
            open( ftn, file=SACHTET_struc(n)%sacsoilparms_table, form="formatted" )
            read( ftn, * ); read( ftn, * )
            do k = 1, 12
               read(ftn,*) type_value, type_name, satdk_default(k), &
                    sand_default(k), clay_default(k)
            enddo
            call LDT_releaseUnitNumber(ftn)

         !- Map parameters to soil texture map:
            do r = 1, LDT_rc%lnr(n)
               do c = 1, LDT_rc%lnc(n)
                  if( SACHTET_struc(n)%frz_soiltext%value(c,r,1) .ge. 0. ) then
                     SACHTET_struc(n)%clay_table%value(c,r,1) = &
                          clay_default(nint(SACHTET_struc(n)%frz_soiltext%value(c,r,1)))
                     SACHTET_struc(n)%sand_table%value(c,r,1) = &
                          sand_default(nint(SACHTET_struc(n)%frz_soiltext%value(c,r,1)))
                     SACHTET_struc(n)%ksat_table%value(c,r,1) = &
                          satdk_default(nint(SACHTET_struc(n)%frz_soiltext%value(c,r,1)))
                  else
                     SACHTET_struc(n)%clay_table%value(c,r,1) = LDT_rc%udef
                     SACHTET_struc(n)%sand_table%value(c,r,1) = LDT_rc%udef
                     SACHTET_struc(n)%ksat_table%value(c,r,1) = LDT_rc%udef
                  endif
               end do
            end do
         endif

         !- Read in vegetation parameter table values:
         if( SACHTET_struc(n)%sacvegparms_table == "none" ) then
            write(LDT_logunit,*) "[INFO] SAC-HTET -- No vegetation parameters read in "
            write(LDT_logunit,*) "[INFO] SAC-HTET -- from look-up table."
            write(LDT_logunit,*) "[INFO] If you want to read in look-up vegetation"
            write(LDT_logunit,*) "[INFO] parameters, you need to enter the full path/filename"
         endif
         if( SACHTET_struc(n)%sacvegparms_table .ne. "none" ) then
            write(LDT_logunit,*) "[INFO] SAC-HTET -- Reading in vegetation parameters from table."
            ftn = LDT_getNextUnitNumber()
            open( ftn, file=SACHTET_struc(n)%sacvegparms_table, form="formatted" )
            do k = 1, 3; read( ftn, * ); enddo
               do k = 1, LDT_LSMparam_struc(n)%landcover%vlevels
                  read(ftn,*) type_value, type_name, &
                       rcmin_default(k), rcminclim_default(k), rgl_default(k), hs_default(k), &
                       lai_default(k), z0_default(k), d50_default(k), croot_default(k)
               enddo
               call LDT_releaseUnitNumber(ftn)

            !- Map parameters to vegetation class map:
               SACHTET_struc(n)%rcmin%value = LDT_rc%udef
               SACHTET_struc(n)%rcminclim%value = LDT_rc%udef
               SACHTET_struc(n)%rgl%value = LDT_rc%udef
               SACHTET_struc(n)%hs%value = LDT_rc%udef
               SACHTET_struc(n)%maxlai%value = LDT_rc%udef
               SACHTET_struc(n)%z0%value = LDT_rc%udef
               SACHTET_struc(n)%d50%value = LDT_rc%udef
               SACHTET_struc(n)%sachtet_croot%value = LDT_rc%udef

               numveg = LDT_LSMparam_struc(n)%landcover%vlevels-1
               do r = 1, LDT_rc%lnr(n)
                  do c = 1, LDT_rc%lnc(n)
                     !         if( sum(LDT_LSMparam_struc(n)%landcover%value(c,r,(1:numveg)),&
                     !                 mask=vegcnt(c,r,1:numveg).ne.LDT_rc%udef ) .gt. 0. ) then
                     do k = 1, LDT_LSMparam_struc(n)%landcover%vlevels
                        if( LDT_LSMparam_struc(n)%landcover%value(c,r,k) > 0. ) then
                           SACHTET_struc(n)%rcmin%value(c,r,1) = &
                                rcmin_default(k)
                           SACHTET_struc(n)%rcminclim%value(c,r,1) = &
                                rcminclim_default(k)
                           SACHTET_struc(n)%rgl%value(c,r,1) = &
                                rgl_default(k)
                           SACHTET_struc(n)%hs%value(c,r,1) = &
                                hs_default(k)
                           SACHTET_struc(n)%maxlai%value(c,r,1) = &
                                lai_default(k)
                           SACHTET_struc(n)%z0%value(c,r,1) = &
                                z0_default(k)
                           SACHTET_struc(n)%d50%value(c,r,1) = &
                                d50_default(k)
                           SACHTET_struc(n)%sachtet_croot%value(c,r,1) = &
                                croot_default(k)
                        endif
                     enddo
                  enddo
               enddo
            endif  ! end vegetation parameter table assignments

         end do    ! End nest loop
      end if       ! End SAC-HTET/SNOW-17 check


! - PET/PET-Adjustment Parameter Section:

   write(LDT_logunit,*)" - - - - - - - - -  PET Parameters  - - - - - - - - -"

   call ESMF_ConfigFindLabel(LDT_config,"PET directory:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%petdir,rc=rc)
      call LDT_verify(rc,'PET directory: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"PET adjustment factor directory:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%petadjdir,rc=rc)
      call LDT_verify(rc,'PET adjustment factor directory: not specified')
   enddo

   call ESMF_ConfigFindLabel(LDT_config,"PET climatology interval:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,SACHTET_struc(n)%petInterval,rc=rc)
      call LDT_verify(rc,'PET climatology interval: not specified')

      if( SACHTET_struc(n)%petInterval == "monthly" ) then
         LDT_rc%monthlyData(n) = .true.
         SACHTET_struc(n)%pet%num_times = 12
         SACHTET_struc(n)%petadj%num_times = 12
         write(LDT_logunit,*) "[INFO] Reading in 'monthly' PET/PETAdj parameters"
      else
         write(LDT_logunit,*) "[ERR] No other time dimension option for PET/PETAdj"
         write(LDT_logunit,*) "  exists at this time, except 'monthly'."
      endif
   enddo

   do n=1,LDT_rc%nnest
      SACHTET_struc(n)%pet%vlevels = SACHTET_struc(n)%pet%num_times
      SACHTET_struc(n)%petadj%vlevels = SACHTET_struc(n)%petadj%num_times

      allocate(SACHTET_struc(n)%PET%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               SACHTET_struc(n)%PET%vlevels))       
      allocate(SACHTET_struc(n)%PETADJ%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               SACHTET_struc(n)%PETADJ%vlevels))
   enddo

   pet%filltype = "none"
   call ESMF_ConfigGetAttribute(LDT_config, pet%filltype, &
             label="PET fill option:",rc=rc)
   call LDT_verify(rc,"PET fill option: option not specified in the config file")

   if( pet%filltype == "average" .or. pet%filltype == "neighbor" ) then
      call ESMF_ConfigGetAttribute(LDT_config, pet%fillradius, &
                label="PET fill radius:",rc=rc)
      call LDT_verify(rc,"PET fill radius: option not specified in the config file")

      call ESMF_ConfigGetAttribute(LDT_config, pet%fillvalue, &
                label="PET fill value:",rc=rc)
      call LDT_verify(rc,"PET fill value: option not specified in the config file")
   elseif( pet%filltype == "none" ) then
      write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for PET"
   else
      write(LDT_logunit,*) "[ERR] Fill option for PET is not valid: ",trim(pet%filltype)
      write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
      write(LDT_logunit,*) "  Programming stopping ..."
      call LDT_endrun
   end if

   call ESMF_ConfigGetAttribute(LDT_config,pet_proj,&
             label="PET map projection:",rc=rc)
   call LDT_verify(rc,'PET map projection: option not specified in the config file')
   SACHTET_struc(:)%pet_proj = pet_proj

   call ESMF_ConfigFindLabel(LDT_config,"PET spatial transform:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,&
                SACHTET_struc(n)%pet_gridtransform, rc=rc)
      call LDT_verify(rc,'PET spatial transform: option not specified in the config file')
   enddo

!- Read in files:
   do n=1,LDT_rc%nnest        

    ! Monthly time-interval:
      if( SACHTET_struc(n)%petInterval.eq."monthly" ) then 

         do k = 1, SACHTET_struc(n)%pet%vlevels   ! months

            call ESMF_ConfigGetAttribute(rdhmconsts_table, temp, &
                 label='pe_'//sacmonths(k)//'=',rc=rc)
            if( rc /= 0 ) temp = SACHTET_struc(n)%rdhm_undef
            write(LDT_logunit,'(A20,1X,f12.5)') " pe_"//sacmonths(k)//" =", temp
            SACHTET_struc(n)%pet%value(:,:,k) = temp
            if( temp >= 0. ) cycle    ! Assign constant value and not grid               

            SACHTET_struc(n)%petfile = trim(SACHTET_struc(n)%petdir)//'_'//&
                                       trim(sacmonths(k))//'.gz'

            write(LDT_logunit,*) "Reading PET file: "//trim(SACHTET_struc(n)%petfile)
            call read_SACHTET356_pet(n, SACHTET_struc(n)%pet%value(:,:,k))
            write(LDT_logunit,*) "Done reading "//trim(SACHTET_struc(n)%petfile)

         enddo     ! End month loop

         do k = 1, SACHTET_struc(n)%petadj%vlevels   ! months

            call ESMF_ConfigGetAttribute(rdhmconsts_table, temp, &
                            label='peadj_'//sacmonths(k)//'=',rc=rc)
            if( rc /= 0 ) temp = SACHTET_struc(n)%rdhm_undef
            write(LDT_logunit,'(A20,1X,f12.5)') " peadj_"//sacmonths(k)//" =", temp
            SACHTET_struc(n)%petadj%value(:,:,k) = temp
            if( temp >= 0. ) cycle    ! Assign constant value and not grid               

            SACHTET_struc(n)%petadjfile = trim(SACHTET_struc(n)%petadjdir)//'_'//&
                                          trim(sacmonths(k))//'.gz'

            write(LDT_logunit,*) "Reading PET-ADJ file: "//trim(SACHTET_struc(n)%petadjfile)
            call read_SACHTET356_petadj(&
                      n,SACHTET_struc(n)%petadj%value(:,:,k))
            write(LDT_logunit,*) "Done reading "//trim(SACHTET_struc(n)%petadjfile)

         enddo  ! End month loop
      endif     ! End interval check

   !- Parameter-Mask agreement section:
      if( pet%filltype == "average" .or. pet%filltype == "neighbor" ) then
          write(LDT_logunit,*) "Checking/filling mask values for: ", &
               trim(SACHTET_struc(n)%pet%short_name)
          write(fill_logunit,*) "Checking/filling mask values for: ", &
               trim(SACHTET_struc(n)%pet%short_name)
          pet%watervalue = LDT_rc%udef
            !           fill_option = "average"
            !           fill_value = 10.0
            !           fill_rad = 2.
          call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                 SACHTET_struc(n)%pet_gridtransform, SACHTET_struc(n)%pet%num_times, &
                 SACHTET_struc(n)%pet%value, pet%watervalue,          &
                 LDT_LSMparam_struc(n)%landmask2%value, pet%filltype,       &
                 pet%fillvalue, pet%fillradius )
      endif   ! end PET read

   enddo      ! End nest loop

!- Bottom soil temperature (K) field:
   do n = 1, LDT_rc%nnest
      if( SACHTET_struc(n)%tbot%selectOpt == 1 ) then 
        write(LDT_logunit,*)" - - - - - - - - - Bottom Temperature Parameter - - - - - - - - - - - -"

         allocate(SACHTET_struc(n)%tbot%value(&
                          LDT_rc%lnc(n),LDT_rc%lnr(n),&
                          SACHTET_struc(n)%tbot%vlevels))       

         allocate(force_elev(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       ! Determine if bottom temperature should be set/read in:
         call ESMF_ConfigGetAttribute(rdhmconsts_table, temp, &
              label='frz_TBOT=',rc=rc)
         if( rc /= 0 ) temp = SACHTET_struc(n)%rdhm_undef
         write(LDT_logunit,'(A20,1X,f12.5)') " frz_TBOT =", temp
         SACHTET_struc(n)%tbot%value(:,:,:) = temp
         if( temp >= 0. ) exit    ! Assign constant value and not grid          

       ! Read in Tbot File:
         write(LDT_logunit,*) "Reading: "//trim(SACHTET_struc(n)%tbot_file)
         call read_SACHTET356_tbot(n, SACHTET_struc(n)%tbot%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(SACHTET_struc(n)%tbot_file)

#if 0 
       ! Fill where parameter values are missing compared to land/water mask:
         if( tbot%filltype == "neighbor" .or. &
              tbot%filltype == "average" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                 trim(SACHTET_struc(n)%tbot%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                 trim(SACHTET_struc(n)%tbot%short_name)
            tbot%watervalue = LDT_rc%udef
            call LDT_contIndivParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
!                 SACHTET_struc(n)%tbot_gridtransform,           &
!                 SACHTET_struc(n)%sachtetparms_gridtransform,   &
                 "none",           &
                 SACHTET_struc(n)%tbot%num_bins,                &
                 SACHTET_struc(n)%tbot%value, tbot%watervalue,  &
                 LDT_LSMparam_struc(n)%landmask2%value,         &
                 tbot%filltype, tbot%fillvalue, tbot%fillradius )
         endif
#endif

         !- Modify final Tbot output with elevation correction:
         if( SACHTET_struc(n)%tbot_topocorr == "lapse-rate" ) then
            if( LDT_LSMparam_struc(n)%elevation%selectOpt == 1 ) then
               write(LDT_logunit,*) "Performing lapse-rate correction to Tbot output."
               force_elev = 0.
               do r = 1, LDT_rc%lnr(n)
                  do c = 1, LDT_rc%lnc(n)
                     if( SACHTET_struc(n)%tbot%value(c,r,1)/=LDT_rc%udef ) &
                          SACHTET_struc(n)%tbot%value(c,r,1) =       &
                          SACHTET_struc(n)%tbot%value(c,r,1)     &
                          + (-0.0065)*(LDT_LSMparam_struc(n)%elevation%value(c,r,1)  &
                          - force_elev(c,r))
                  end do
               end do
            elseif( LDT_LSMparam_struc(n)%elevation%selectOpt == 0 ) then
               write(LDT_logunit,*) "Cannot perform lapse-rate correction to Tbot output,"
               write(LDT_logunit,*) " since no elevation/terrain map option was selected. "
               write(LDT_logunit,*) " Stopping ... "
               call LDT_endrun
            endif
         endif

         deallocate(force_elev)
      end if  ! Tbot block 

      !== Other Noah LSM related parameters ==

   enddo

 end subroutine SACHTETparms_init


 subroutine SACHTETparms_writeHeader(n,ftn,dimID,monthID)
   
   integer   :: n 
   integer   :: ftn
   integer   :: dimID(3)
   integer   :: monthID
   
   integer   :: t_dimID(3)
   
   if( SACHTET_struc(n)%sachtet356%selectOpt == 1 ) then
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%lzfpm)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%lzfsm)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%lztwm)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%lzpk)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%lzsk)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%uzfwm)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%uztwm)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%uzk)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%pfree)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%rexp)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%zperc)
      
      if( SACHTET_struc(n)%sacsoilparms_table .ne. "none" .and. &
           SACHTET_struc(n)%stxt_file .ne. "none" ) then
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%clay_table)
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%sand_table)
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%ksat_table)
      endif
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%forestfrac)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%soilalb)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%pctim)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%adimp)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%riva)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%rserv)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%side)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%offsetTime)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%frz_soiltext)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%frz_cksl)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%frz_rsmax)
      
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%frz_zbot)
      
      if( SACHTET_struc(n)%sacvegparms_table .ne. "none" ) then
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%rcmin)
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%rcminclim)
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%rgl)
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%hs)
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%maxlai)
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%z0)
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%d50)
         call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
              SACHTET_struc(n)%sachtet_croot)
      endif

      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           SACHTET_struc(n)%tbot)
      
      if(SACHTET_struc(n)%petInterval.eq."monthly") then !monthly
        t_dimID(1) = dimID(1)
        t_dimID(2) = dimID(2)
        t_dimID(3) = monthID       
        call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
                                       SACHTET_struc(n)%pet)
        call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
                                       SACHTET_struc(n)%petadj)
   
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"PET_DATA_INTERVAL", &
                        SACHTET_struc(n)%petInterval))
      endif

   endif
   
 end subroutine SACHTETparms_writeHeader
 
  subroutine SACHTETparms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    if( SACHTET_struc(n)%sachtet356%selectOpt == 1 ) then

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%lzfpm)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%lzfsm)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%lztwm)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%lzpk)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%lzsk)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%uzfwm)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%uztwm)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%uzk)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%pfree)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%rexp)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%zperc)

       if( SACHTET_struc(n)%sacsoilparms_table .ne. "none" .and. &
            SACHTET_struc(n)%stxt_file .ne. "none" ) then
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%clay_table)
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%sand_table)
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%ksat_table)
       endif

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%forestfrac)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%soilalb)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%pctim)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%adimp)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%riva)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%rserv)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%side)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%offsetTime)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%frz_soiltext)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%frz_cksl)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%frz_rsmax)

       call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%frz_zbot)

       if( SACHTET_struc(n)%sacvegparms_table .ne. "none" ) then
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%rcmin)
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%rcminclim)
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%rgl)
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%hs)
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%maxlai)
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%z0)
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%d50)
          call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%sachtet_croot)
       endif

    endif

    call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%tbot)

    call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%pet)

    call LDT_writeNETCDFdata(n,ftn,SACHTET_struc(n)%petadj)

  end subroutine SACHTETparms_writeData

!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name, &
                               units, full_name )

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name
   character(len=*),     optional     :: units
   character(len=*),     optional     :: full_name

   character(20) :: unit_temp
   character(100):: name_temp
! ____________________________________________________
   
   if(present(units)) then
      unit_temp = units
   else
      unit_temp = "none"
   endif
   if(present(full_name)) then
      name_temp = full_name
   else
      name_temp = trim(short_name)
   endif

   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = 1
   paramEntry%selectOpt = 1
   paramEntry%source = "SACHTET.3.5.6"
   paramEntry%units = trim(unit_temp)
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

end module SACHTET_parmsMod

