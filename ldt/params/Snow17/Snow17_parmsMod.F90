!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module Snow17_parmsMod
!BOP
!
! !MODULE: Snow17_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read Snow17 parameter
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  Snow17 parameter file data.
!
! !REVISION HISTORY:
!
!  04 Nov 2013: K. Arsenault: Added layers for SNOW-17 model
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: Snow17Parms_init    !allocates memory for required structures
  public :: Snow17Parms_writeHeader
  public :: Snow17Parms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: Snow17_struc

  type, public :: snow17_type_dec

     character*140 :: rdhmconsts_table
     real          :: rdhm_undef

     real          :: snow17parms_gridDesc(20)
     character*50  :: snow17parms_proj
     character*50  :: snow17parms_gridtransform

     character*140 :: mfmax_file      !
     character*140 :: mfmin_file      !
     character*140 :: uadj_file       !
     character*140 :: alat_file       !
     character*140 :: elev_file       !
     character*140 :: scf_file        !
     character*140 :: nmf_file        !
     character*140 :: si_file         ! 
     character*140 :: mbase_file      !
     character*140 :: pxtemp_file     !
     character*140 :: plwhc_file      !
     character*140 :: tipm_file       !
     character*140 :: pgm_file        !
     character*140 :: laec_file       !
     character*140 :: adc_dir         ! 
     character*140 :: adc_file        ! 

! -  Snow-17 model-specific:
     type(LDT_paramEntry) :: snow17      ! SNOW-17 model parameters (collective)
     type(LDT_paramEntry) :: mfmax       ! Max. melt factor (mm/6hr degC)
     type(LDT_paramEntry) :: mfmin       ! Min. melt factor (mm/6hr degC)
     type(LDT_paramEntry) :: uadj        ! Ave. wind function during rain-on-snow periods (mm/mb)
     type(LDT_paramEntry) :: alat        ! Snow-17 latitude parameter  (-)
     type(LDT_paramEntry) :: snow17_elev ! Snow-17 elevation (a priori; m)
     type(LDT_paramEntry) :: scf         ! Snowfall gage catch definiency adjustment factor (-)
     type(LDT_paramEntry) :: nmf         ! Maximum negative melt factor (mm/6hr degC)
     type(LDT_paramEntry) :: si          ! Snow water-equivalent above which 100% cover exists (mm)
     type(LDT_paramEntry) :: mbase       ! Base temperature for non-rain melt factor (C) 
     type(LDT_paramEntry) :: pxtemp      ! Rain/snow delineation temperature (C)
     type(LDT_paramEntry) :: plwhc       ! Decimal percent liquid water holding capacity (-)
     type(LDT_paramEntry) :: pgm         ! Ground melt rate (mm/day)
     type(LDT_paramEntry) :: tipm        ! Antecedent temperature index parameter (-)
     type(LDT_paramEntry) :: laec        ! Rain-snow split temperature, 0=temp threshold, not 0=rain-snow elev (-)
     type(LDT_paramEntry) :: adc         ! 11 points on areal snow depletion curve (-)

  end type snow17_type_dec

  type(snow17_type_dec), allocatable :: Snow17_struc(:)

contains

!BOP
! 
! !ROUTINE: Snow17Parms_init
! \label{Snow17Parms_init}
! 
! !INTERFACE:
  subroutine Snow17Parms_init

! !USES:
!   use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
   use LDT_logMod,    only : LDT_verify, LDT_endrun
!   use LDT_paramOptCheckMod, only: LDT_snow17parmsOptChecks, &
!                      LDT_gridOptChecks
   use LDT_xmrg_reader
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the Research Distributed Hydrology Model (RDHM)
! SNOW-17 model parameter datasets.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[snow17Parmssetup](\ref{snow17Parmssetup}) \newline
!    calls the registry to invoke the snow17Parms setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer  :: n
   integer  :: c,r,m,k
   integer  :: rc
   real     :: temp 
   real     :: xfactor
   logical  :: file_exists
   logical  :: snow17_select
   character(1) :: numadc1
   character(2) :: numadc2
   type(ESMF_Config)  :: rdhmconsts_table
   type(LDT_fillopts) :: snow17
   character*50       :: snow17parms_proj

! _____________________________________________________________________

   allocate(Snow17_struc(LDT_rc%nnest))

   snow17_select = .false.
   do n=1,LDT_rc%nnest

      ! - Snow-17 parameters:
      call set_param_attribs(Snow17_struc(n)%snow17,"SNOW17")
      ! --
      call set_param_attribs(Snow17_struc(n)%adc,"ADC",&
            units="-", &
            full_name="SNOW17 Points on Areal (snow) Depletion Curve")

      if( Snow17_struc(n)%snow17%selectOpt.gt.0 ) then
         snow17_select = .true.
      endif
   enddo

 
  if( snow17_select ) then
    write(LDT_logunit,*)" - - - - - - - - - - SNOW-17 LSM Parameters - - - - - - - - - - - - -"

!- Load RDHM CONSTANTS input table file (filepath read-in from ldt.config file)
   call ESMF_ConfigFindLabel(LDT_config,"RDHM356 constants table:",rc=rc)
   do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%rdhmconsts_table,rc=rc)
     call LDT_verify(rc,'RDHM356 constants table: not specified')

     inquire(file=trim(Snow17_struc(n)%rdhmconsts_table), exist=file_exists)
     if( .not. file_exists ) then
        write(LDT_logunit,*) "[ERR] RDHM Parameter Constants Table ",&
              trim(Snow17_struc(n)%rdhmconsts_table)," does not exist."
        call LDT_endrun
     endif
     write(LDT_logunit,*) "Reading in RDHM Parameter Constants Table Entries"
     rdhmconsts_table = ESMF_ConfigCreate(rc=rc)
     call ESMF_ConfigLoadFile(rdhmconsts_table, &
                     trim(Snow17_struc(n)%rdhmconsts_table), rc=rc)
   enddo

   call ESMF_ConfigFindLabel(LDT_config,"RDHM356 universal undefined value:",rc=rc)
   do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%rdhm_undef,rc=rc)
     call LDT_verify(rc,'RDHM356 universal undefined value: not specified')
   enddo

 ! Enter number of ADC points
   Snow17_struc(:)%adc%num_bins = 11
   call ESMF_ConfigFindLabel(LDT_config,"SNOW17 ADC number of points:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%adc%num_bins,&
           rc=rc)  ! ,default=11,
      call LDT_verify(rc,'SNOW17 ADC number of points: not specified')
      if( Snow17_struc(n)%adc%num_bins > 11 ) then
         write(LDT_logunit,*) " Currently, the number of ADC points should not exceed"
         write(LDT_logunit,*) "  the standard 11 values.  This may be expanded in the future."
         write(LDT_logunit,*) " Calling end run ..."
         call LDT_endrun
      endif
   enddo

! -----

!   allocate(LDT_rc%snow17parms_gridDesc(LDT_rc%nnest,20))

   do n = 1, LDT_rc%nnest
     if( Snow17_struc(n)%snow17%selectOpt == 1 ) then 

     ! Fill in derived parameter entries:
     ! ( input_parmattribs -> output_parmattribs ) 
       call populate_param_attribs( "RDHM356_MFMAX", &
             "Max. melt factor","mm/6hr degC", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%mfmax )

       allocate(Snow17_struc(n)%mfmax%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%mfmax%vlevels))       

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_MFMAX=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_MFMAX =", temp
       Snow17_struc(n)%mfmax%value = temp

       call populate_param_attribs( "RDHM356_MFMIN", &
             "Min. melt factor","mm/6hr degC", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%mfmin )

       allocate(Snow17_struc(n)%mfmin%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%mfmin%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_MFMIN=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_MFMIN =", temp
       Snow17_struc(n)%mfmin%value = temp

       call populate_param_attribs( "RDHM356_UADJ", &
             "Ave. wind function during rain-on-snow periods","mm/mb", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%uadj )

       allocate(Snow17_struc(n)%uadj%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%uadj%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_UADJ=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_UADJ =", temp
       Snow17_struc(n)%uadj%value = temp

       call populate_param_attribs( "RDHM356_ALAT", &
             "Snow-17 latitude parameter","deg", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%alat )

       allocate(Snow17_struc(n)%alat%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%alat%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_ALAT=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_ALAT =", temp
       Snow17_struc(n)%alat%value = temp

       call populate_param_attribs( "RDHM356_ELEV", &
             "Snow-17 elevation (a priori)","m", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%snow17_elev )

       allocate(Snow17_struc(n)%snow17_elev%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%snow17_elev%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_ELEV=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_ELEV =", temp
       Snow17_struc(n)%snow17_elev%value = temp

       call populate_param_attribs( "RDHM356_SCF", &
             "Snowfall gage catch definiency adjustment factor","-", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%scf )

       allocate(Snow17_struc(n)%scf%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%scf%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_SCF=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_SCF =", temp
       Snow17_struc(n)%scf%value = temp

       call populate_param_attribs( "RDHM356_NMF", &
             "Maximum negative melt factor","mm/6hr degC", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%nmf )

       allocate(Snow17_struc(n)%nmf%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%nmf%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_NMF=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_NMF =", temp
       Snow17_struc(n)%nmf%value = temp

       call populate_param_attribs( "RDHM356_SI", &
             "Snow water-equivalent above which 100% cover exists","mm", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%si )

       allocate(Snow17_struc(n)%si%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%si%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_SI=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_SI =", temp
       Snow17_struc(n)%si%value = temp

       call populate_param_attribs( "RDHM356_MBASE", &
             "Base temperature for non-rain melt factor","degC", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%mbase )

       allocate(Snow17_struc(n)%mbase%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%mbase%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_MBASE=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_MBASE =", temp
       Snow17_struc(n)%mbase%value = temp

       call populate_param_attribs( "RDHM356_PXTEMP", &
             "Rain/snow delineation temperature","degC", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%pxtemp )

       allocate(Snow17_struc(n)%pxtemp%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%pxtemp%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_PXTMP=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_PXTMP =", temp
       Snow17_struc(n)%pxtemp%value = temp

       call populate_param_attribs( "RDHM356_PLWHC", &
             "Decimal percent liquid water holding capacity","-", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%plwhc )

       allocate(Snow17_struc(n)%plwhc%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%plwhc%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_PLWHC=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_PLWHC =", temp
       Snow17_struc(n)%plwhc%value = temp

       call populate_param_attribs( "RDHM356_PGM", &
             "Ground melt rate","mm/day", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%pgm )

       allocate(Snow17_struc(n)%pgm%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%pgm%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_PGM=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_PGM =", temp
       Snow17_struc(n)%pgm%value = temp

       call populate_param_attribs( "RDHM356_TIPM", &
             "Antecedent temperature index parameter","-", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%tipm )

       allocate(Snow17_struc(n)%tipm%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%tipm%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_TIPM=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_TIPM =", temp
       Snow17_struc(n)%tipm%value = temp

       call populate_param_attribs( "RDHM356_LAEC", &
             "Rain-snow split temperature: 0=temp threshold, !0=rain-snow elev","-", &
              Snow17_struc(n)%snow17, &
              Snow17_struc(n)%laec )

       allocate(Snow17_struc(n)%laec%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%laec%vlevels))

       call ESMF_ConfigGetAttribute(rdhmconsts_table,temp,label='snow_LAEC=',rc=rc)
       if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
       write(LDT_logunit,'(A20,1X,f12.5)') " snow_LAEC =", temp
       Snow17_struc(n)%laec%value = temp

       Snow17_struc(n)%adc%vlevels = Snow17_struc(n)%adc%num_bins
       allocate(Snow17_struc(n)%adc%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           Snow17_struc(n)%adc%vlevels))

       do k = 1, Snow17_struc(n)%adc%vlevels
          if(k < 10) then
             write( numadc1, '(i1)' ) k
             call ESMF_ConfigGetAttribute(rdhmconsts_table, temp, &
                  label='snow_ADC'//numadc1//'=',rc=rc)
             if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
             write(LDT_logunit,'(A20,1X,f12.5)') " snow_ADC"//numadc1//" =",temp
          elseif( k > 9 ) then
             write( numadc2, '(i2)' ) k
             call ESMF_ConfigGetAttribute(rdhmconsts_table, temp, &
                  label='snow_ADC'//numadc2//'=',rc=rc)
             if( rc .ne. 0 ) temp = Snow17_struc(n)%rdhm_undef
             write(LDT_logunit,'(A20,1X,f12.5)') " snow_ADC"//numadc2//" =",temp
          endif
          Snow17_struc(n)%adc%value(:,:,k) = temp
       end do

     end if
   end do

!- Read in ldt.config file entries:
   if( snow17_select ) then 

   !- Read in SNOW-17 soil parameters from a-priori gridded maps:

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 MFMAX map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%mfmax_file,rc=rc)
         call LDT_verify(rc,'SNOW17 MFMAX map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 MFMIN map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%mfmin_file,rc=rc)
         call LDT_verify(rc,'SNOW17 MFMIN map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 UADJ map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%uadj_file,rc=rc)
         call LDT_verify(rc,'SNOW17 UADJ map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 ALAT map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%alat_file,rc=rc)
         call LDT_verify(rc,'SNOW17 ALAT map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 ELEV map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%elev_file,rc=rc)
         call LDT_verify(rc,'SNOW17 ELEV map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 SCF map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%scf_file,rc=rc)
         call LDT_verify(rc,'SNOW17 SCF map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 NMF map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%nmf_file,rc=rc)
         call LDT_verify(rc,'SNOW17 NMF map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 SI map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%si_file,rc=rc)
         call LDT_verify(rc,'SNOW17 SI map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 MBASE map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%mbase_file,rc=rc)
         call LDT_verify(rc,'SNOW17 MBASE map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 PXTEMP map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%pxtemp_file,rc=rc)
         call LDT_verify(rc,'SNOW17 PXTEMP map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 PLWHC map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%plwhc_file,rc=rc)
         call LDT_verify(rc,'SNOW17 PLWHC map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 PGM map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%pgm_file,rc=rc)
         call LDT_verify(rc,'SNOW17 PGM map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 TIPM map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%tipm_file,rc=rc)
         call LDT_verify(rc,'SNOW17 TIPM map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 LAEC map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%laec_file,rc=rc)
         call LDT_verify(rc,'SNOW17 LAEC map: not specified')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 ADC directory:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%adc_dir,rc=rc)
         call LDT_verify(rc,'SNOW17 ADC directory: not specified')
      enddo

! --  Grid info and inputs:
      call ESMF_ConfigFindLabel(LDT_config,"SNOW17 parameter spatial transform:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Snow17_struc(n)%snow17parms_gridtransform,&
              rc=rc)
         call LDT_verify(rc,'SNOW17 parameter spatial transform: option not specified in the config file')
      enddo

      snow17%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, snow17%filltype, &
           label="SNOW17 parameter fill option:",rc=rc)
      call LDT_verify(rc,"SNOW17 parameter fill option: option not specified in the config file")

      if( snow17%filltype == "neighbor" .or. snow17%filltype == "average" ) then
        call ESMF_ConfigGetAttribute(LDT_config, snow17%fillvalue, &
             label="SNOW17 parameter fill value:",rc=rc)
        call LDT_verify(rc,"SNOW17 parameter fill value: option not specified in the config file")

        call ESMF_ConfigGetAttribute(LDT_config, snow17%fillradius, &
             label="SNOW17 parameter fill radius:",rc=rc)
        call LDT_verify(rc,"SNOW17 parameter fill radius: option not specified in the config file")
      elseif( snow17%filltype == "none" ) then
         write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for SNOW17 parameters"
      else
         write(LDT_logunit,*) "[ERR] Fill option for SNOW17 is not valid: ",trim(snow17%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
      end if

      call ESMF_ConfigGetAttribute(LDT_config,snow17parms_proj,&
           label="SNOW17 map projection:",rc=rc)
      call LDT_verify(rc,'SNOW17 map projection: option not specified in the config file')
      Snow17_struc(:)%snow17parms_proj = snow17parms_proj

#if 0
    ! LDT_rc%snow17parms_gridDesc
!      call LDT_snow17parmsOptChecks( n, "SNOW17", LDT_rc%snow17parms_proj, &
!                      Snow17_struc(n)%snow17parms_gridtransform )
#endif

   end if   ! End SNOW-17 parameter array definitions


!- Loop over all domain nests:
   do n = 1, LDT_rc%nnest
     if( Snow17_struc(n)%snow17%selectOpt == 1 ) then 

!        call LDT_gridOptChecks( n, "SNOW17", Snow17_struc(n)%snow17parms_gridtransform, &
!                     LDT_rc%snow17parms_proj, LDT_rc%snow17parms_gridDesc(n,9) )
      
       if( Snow17_struc(n)%mfmax%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%mfmax%value(1,1,1))
         write(LDT_logunit,*) "Reading MFMAX file: "//trim(Snow17_struc(n)%mfmax_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%mfmax_file, LDT_rc%udef,  &
                  Snow17_struc(n)%mfmax%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%mfmax_file)
         Snow17_struc(n)%mfmax%value(:,:,1) = &
                Snow17_struc(n)%mfmax%value(:,:,1) * xfactor
       end if

       if( Snow17_struc(n)%mfmin%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%mfmin%value(1,1,1))
         write(LDT_logunit,*) "Reading MFMIN file: "//trim(Snow17_struc(n)%mfmin_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%mfmin_file, LDT_rc%udef,  &
                  Snow17_struc(n)%mfmin%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%mfmin_file)
         Snow17_struc(n)%mfmin%value(:,:,1) = &
                Snow17_struc(n)%mfmin%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%uadj%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%uadj%value(1,1,1))
         write(LDT_logunit,*) "Reading UADJ file: "//trim(Snow17_struc(n)%uadj_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%uadj_file, LDT_rc%udef,  &
                  Snow17_struc(n)%uadj%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%uadj_file)
         Snow17_struc(n)%uadj%value(:,:,1) = &
                Snow17_struc(n)%uadj%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%alat%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%alat%value(1,1,1))
         write(LDT_logunit,*) "Reading ALAT file: "//trim(Snow17_struc(n)%alat_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%alat_file, LDT_rc%udef,  &
                  Snow17_struc(n)%alat%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%alat_file)
         Snow17_struc(n)%alat%value(:,:,1) = &
                Snow17_struc(n)%alat%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%snow17_elev%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%snow17_elev%value(1,1,1))
         write(LDT_logunit,*) "Reading SNOW17 Elev file: "//trim(Snow17_struc(n)%elev_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%elev_file, LDT_rc%udef,  &
                  Snow17_struc(n)%snow17_elev%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%elev_file)
         Snow17_struc(n)%snow17_elev%value(:,:,1) = &
                Snow17_struc(n)%snow17_elev%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%scf%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%scf%value(1,1,1))
         write(LDT_logunit,*) "Reading SNOW17 SCF file: "//trim(Snow17_struc(n)%scf_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%scf_file, LDT_rc%udef,  &
                  Snow17_struc(n)%scf%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%scf_file)
         Snow17_struc(n)%scf%value(:,:,1) = &
                Snow17_struc(n)%scf%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%nmf%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%nmf%value(1,1,1))
         write(LDT_logunit,*) "Reading SNOW17 NMF file: "//trim(Snow17_struc(n)%nmf_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%nmf_file, LDT_rc%udef,  &
                  Snow17_struc(n)%nmf%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%nmf_file)
         Snow17_struc(n)%nmf%value(:,:,1) = &
                Snow17_struc(n)%nmf%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%si%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%si%value(1,1,1))
         write(LDT_logunit,*) "Reading SNOW17 SI file: "//trim(Snow17_struc(n)%si_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%si_file, LDT_rc%udef,  &
                  Snow17_struc(n)%si%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%si_file)
         Snow17_struc(n)%si%value(:,:,1) = &
                Snow17_struc(n)%si%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%mbase%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%mbase%value(1,1,1))
         write(LDT_logunit,*) "Reading SNOW17 MBASE file: "//trim(Snow17_struc(n)%mbase_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%mbase_file, LDT_rc%udef,  &
                  Snow17_struc(n)%mbase%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%mbase_file)
         Snow17_struc(n)%mbase%value(:,:,1) =  &
                Snow17_struc(n)%mbase%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%pxtemp%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%pxtemp%value(1,1,1))
         write(LDT_logunit,*) "Reading SNOW17 PXTEMP file: "//trim(Snow17_struc(n)%pxtemp_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%pxtemp_file, LDT_rc%udef,  &
                  Snow17_struc(n)%pxtemp%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%pxtemp_file)
         Snow17_struc(n)%pxtemp%value(:,:,1) = &
                Snow17_struc(n)%pxtemp%value(:,:,1) * xfactor
       endif
   
       if( Snow17_struc(n)%plwhc%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%plwhc%value(1,1,1))
         write(LDT_logunit,*) "Reading SNOW17 PLWHC file: "//trim(Snow17_struc(n)%plwhc_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%plwhc_file, LDT_rc%udef,  &
                  Snow17_struc(n)%plwhc%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%plwhc_file)
         Snow17_struc(n)%plwhc%value(:,:,1) = &
                Snow17_struc(n)%plwhc%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%pgm%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%pgm%value(1,1,1))
         write(LDT_logunit,*) "Reading SNOW17 PGM file: "//trim(Snow17_struc(n)%pgm_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%pgm_file, LDT_rc%udef,  &
                  Snow17_struc(n)%pgm%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%pgm_file)
         Snow17_struc(n)%pgm%value(:,:,1) = Snow17_struc(n)%pgm%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%tipm%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%tipm%value(1,1,1))
         write(LDT_logunit,*) "Reading SNOW17 TIPM file: "//trim(Snow17_struc(n)%tipm_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%tipm_file, LDT_rc%udef,  &
                  Snow17_struc(n)%tipm%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%tipm_file)
         Snow17_struc(n)%tipm%value(:,:,1) = &
                Snow17_struc(n)%tipm%value(:,:,1) * xfactor
       endif

       if( Snow17_struc(n)%laec%value(1,1,1) < 0. ) then
         xfactor = abs(Snow17_struc(n)%laec%value(1,1,1))
         write(LDT_logunit,*) "Reading SNOW17 LAEC file: "//trim(Snow17_struc(n)%laec_file)
         call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                  LDT_rc%gridDesc(n,:), Snow17_struc(n)%laec_file, LDT_rc%udef,  &
                  Snow17_struc(n)%laec%value(:,:,1) )
         write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%laec_file)
         Snow17_struc(n)%laec%value(:,:,1) = &
                Snow17_struc(n)%laec%value(:,:,1) * xfactor
       endif

       do k = 1, Snow17_struc(n)%adc%vlevels
         if( Snow17_struc(n)%adc%value(1,1,k) < 0. ) then
           xfactor = abs(Snow17_struc(n)%adc%value(1,1,k))
           if(k < 10) then 
              write( numadc1, '(i1)' ) k
              Snow17_struc(n)%adc_file = trim(Snow17_struc(n)%adc_dir)//"/ADC"//numadc1//".gz"
           elseif( k > 9 ) then 
              write( numadc2, '(i2)' ) k
              Snow17_struc(n)%adc_file = trim(Snow17_struc(n)%adc_dir)//"/ADC"//numadc2//".gz"
           endif
           write(LDT_logunit,*) "Reading SNOW17 ADC file: "//trim(Snow17_struc(n)%adc_file)
           call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                    LDT_rc%gridDesc(n,:), Snow17_struc(n)%adc_file, LDT_rc%udef,  &
                    Snow17_struc(n)%adc%value(:,:,1) )
           write(LDT_logunit,*) "Done reading: "//trim(Snow17_struc(n)%adc_file)
           Snow17_struc(n)%adc%value(:,:,k) = &
               Snow17_struc(n)%adc%value(:,:,k) * xfactor
         endif
       end do

     end if
   enddo

  endif  ! End SAC-HTET/SNOW-17 check

  end subroutine Snow17Parms_init

  subroutine Snow17Parms_writeHeader(n,ftn,dimID)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
   integer   :: n 
   integer   :: ftn
   integer   :: dimID(3)
   integer   :: adcdimID(3)

   adcdimID(1) = dimID(1)
   adcdimID(2) = dimID(2)

   if( Snow17_struc(1)%snow17%selectOpt == 1 ) then

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%mfmax)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%mfmin)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%uadj)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%alat)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%snow17_elev)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%scf)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%nmf)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%si)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%mbase)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%pxtemp)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%plwhc)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%pgm)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%tipm)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Snow17_struc(n)%laec)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    Snow17_struc(n)%adc%vlevels = &
        Snow17_struc(n)%adc%num_bins

    call LDT_verify(nf90_def_dim(ftn,'adcpoints',&
         Snow17_struc(n)%adc%num_bins,adcdimID(3)))

    call LDT_writeNETCDFdataHeader(n,ftn,adcdimID,&
         Snow17_struc(n)%adc)
#endif

   endif
    
  end subroutine Snow17Parms_writeHeader

  subroutine Snow17Parms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

   if( Snow17_struc(1)%snow17%selectOpt == 1 ) then

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%mfmax)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%mfmin)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%uadj)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%alat)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%snow17_elev)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%scf)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%nmf)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%si)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%mbase)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%pxtemp)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%plwhc)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%pgm)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%tipm)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%laec)

    call LDT_writeNETCDFdata(n,ftn,Snow17_struc(n)%adc)

   endif

  end subroutine Snow17Parms_writeData


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
   paramEntry%source = "Snow17"
   paramEntry%units = trim(unit_temp)
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

end module Snow17_parmsMod

