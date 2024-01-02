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
! !MODULE: ShuffledComplexEvolution
!
! !REVISION HISTORY:
!  09 Jun 2009;  Soni Yatheendradas, Initial version
!
! !DESCRIPTION: 
!   This module contains the different functions required for implementing the
!   SCEUA algorithm into LIS. For details about the SCEUA introduction and
!   LIS implementation, refer to the description below. 
! 
!      \begin{verbatim}
!      This Shuffled Complex Evolution - University of Arizona (SCEUA) code 
!      is adapted into the LIS from the original SAHRA hydroarchive code at: 
!      http://www.sahra.arizona.edu/cgi-bin/hydroware/hydroware.pl?mode=displayAbstract&ID=9
!      The description there is exactly reproduced in following paragraph in quotes:
!      "The SCE-UA method is a general purpose global optimization
!       program. It was originally developed by Dr. Qingyun Duan as 
!       part of his doctoral dissertation work at the Department of 
!       Hydrology and Water Resources, University of Arizona, Tucson, 
!       AZ 85721, USA. The dissertation is entitled "A Global 
!       Optimization Strategy for Efficient and Effective Calibration 
!       of Hydrologic Models". The program has since been modified to 
!       make it easier for use on problems of users' interests. The 
!       algorithm has been described in detail in an article entitled 
!       "Effective and Efficient Global Optimization for Conceptual 
!       Rainfall-Runoff Models", Water Resources Research, Vol 28(4), 
!       pp.1015-1031, 1992; and in an article entitled "A Shuffled 
!       Complex Evolution Approach for Effective and Efficient Global 
!       Minimization" by Q. Duan, V.K. Gupta and S. Sorooshian, Journal 
!       of Optimization Theory and its Applications, Vol 76(3), pp 501-521, 
!       1993. A paper entitled "Optimal Use of the SCE-UA Global 
!       Optimization Method for Calibrating Watershed Models", by Q. Duan, 
!       S. Sorooshian, & V.K. Gupta, Journal of Hydrology, Vol.158, 265-284, 
!       1994, discussed how to use the SCE-UA Method in an efficient and 
!       effective manner."
!
!       Variable/subroutine name mapping between original code and LIS SCEUA is:
!       (Original code variable/function name, then LIS name after colon below)
!        nopt: sceua_ctl%nparam
!        NGS: sceua_ctl%NComplexes
!        NPG: sceua_ctl%NPointsInComplex
!        NPT: LIS_rc%nensem(n) [LIS_rc%nensem(1) here] OR ga_ctl%npopsize
!        NPS: sceua_ctl%NPointsInSubcomplex
!        NSPL: sceua_ctl%NStepsBeforeShuffle
!        MINGS: sceua_ctl%NMinComplexes
!        MAXN: sceua_ctl%MaxNFuncEvals
!        KSTOP: sceua_ctl%NShufflesForMinChange
!        PCENTO: sceua_ctl%MinChangeInNShuffles
!        ideflt: sceua_ctl%iDefault
!        iniflg: sceua_ctl%InitialPointFlag
!        bound: sceua_ctl%pardel
!        bu: sceua_ctl%parmax
!        bl: sceua_ctl%parmin 
!        a: sceuastruc%parsetinit
!        x: sceuastruc%parsets
!        xx: sceuastruc%parset
!        cx: complexParsets
!        s: subcomplexParsets
!        bestx: sceuastruc%bestparsetinloop
!        worstx: sceuastruc%worstparsetinloop
!        xf: sceuastruc%objfuncs
!        cf: complexObjFuncs
!        sf: subcomplexObjFuncs
!        bestf: sceuastruc%bestobjfuncinloop
!        worstf: sceuastruc%worstobjfuncinloop
!        xnstd: sceuastruc%parstd
!        gnrng: sceuastruc%gmparrng
!        criter: sceuastruc%bestlastobjfuncs
!        icall: sceua_ctl%%SumNewMemberFuncEvals
!        ipcnvg: sceuastruc%isConvergedFlag
!        igs: SCEUAOpt_run :: iComplex
!        nloop: sceua_ctl%iEvolutionLoop
!        loop: iStepsBeforeShuffle
!        criter: sceuastruc%BestObjFuncHistory
!        SUBROUTINE cce_NewPointCreationOnly: SCEUA_SubcomplexCreateChildren
!        SUBROUTINE cce_NewPointAcceptanceOnly: SCEUA_SubcomplexNewGeneration
!
!      \end{verbatim}
!
module ShuffledComplexEvolution
! !USES: 
  USE ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  USE SCEUA_varctl 

  IMPLICIT NONE

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SCEUAOpt_init ! Initialization for SCEUA
  public :: SCEUAOpt_setup ! Setup part of Initialization actually runs LSM
  public :: SCEUAOpt_run ! Non-initialization SCEUA iteration
  public :: SCEUAOpt_checkOptStopCriterion ! Stopping Criterion in SCEUA
  public :: SCEUA_getdecSpaceValues ! Copy all algorithm parameter sets into single array 
  public :: SCEUA_setdecSpaceValues ! Copy above array into algorithm parameter sets 
  public :: SCEUA_getNparam ! Get number of parameters in optimized parameter set
!EOP

  type sceuastruc
     ! SY: This is for single spatial tile location  
     real,    allocatable :: parsets(:,:) ! Co-ordinates of points in the population (parents + potential children)  
     real*8,  allocatable :: parset(:) ! Co-ordinates of a single point in parsets  
     real,    allocatable :: bestparsetinloop(:) ! best point at current shuffling loop  
     real,    allocatable :: worstparsetinloop(:) ! worst point at current shuffling loop  
     real,    allocatable :: objfuncs(:) ! objective function values of parsets?
     real             :: bestobjfuncinloop ! objective function value at best point at current shuffling loop  
     real             :: worstobjfuncinloop ! objective function value at worst point at current shuffling loop  
     real*8,  allocatable :: parstd(:) ! standard deviation of parameters in the spatial tile
     real*8           :: gmparrng ! normalized geometric mean of parameter ranges
     real             :: bestlastobjfuncs(20) ! vector containing the best criterion values of the last 20 shuffling loops
     integer          :: isDecSpaceConvergedFlag ! Did this spatial tile decision space converge yet?
     integer          :: isObjFuncConvergedFlag ! Did this spatial tile best objective function converge yet?
     integer, allocatable :: lcses(:,:) ! Locations of subcomplex in complex, further arrayed over all complexes
     real, allocatable    :: complexesParsets(:,:,:) ! Co-ordinates of complex points, further arrayed over all complexes 
     real, allocatable    :: subcomplexesParsets(:,:,:) ! Co-ordinates of subcomplex points, further arrayed over all complexes
     real, allocatable    :: complexesObjFuncs(:,:) ! objective function values of complex points, further arrayed over all complexes
     real, allocatable    :: subcomplexesObjFuncs(:,:) ! objective function values of subcomplex points, further arrayed over all complexes
     real*8           :: BestObjFuncHistory(20) ! Best obj. Func. history in time
     REAL*8,  allocatable :: parsetinit(:) ! SY: Apriori default par. set
     INTEGER, allocatable :: ibound(:) ! SY: If reflectn. pt. in each complex is within bounds
  end type sceuastruc 

  type(sceuactl)            :: sceua_ctl 
  type(sceuastruc), allocatable :: sceua_struc(:) 


  CONTAINS
   
!BOP
! !ROUTINE: SCEUAOpt_init
! \label{SCEUAOpt_init}
! 
! !INTERFACE: 
    subroutine SCEUAOpt_init()
! !USES: 
      use LIS_coreMod,         only : LIS_rc, LIS_config, LIS_vecTile
      use LIS_optUEMod,        only : LIS_decisionSpace 
      use LIS_logMod,          only : LIS_logunit, LIS_getNextUnitNumber, &
           LIS_releaseUnitNumber, LIS_endrun, LIS_verify
!
! !DESCRIPTION: 
!   This routine performs the initialization steps for the SCEUA.  
!
! !REVISION HISTORY: 
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP      
      implicit none 

      character*100,  allocatable     :: vname(:)
      integer                     :: i,j,t
      integer                     :: ftn
      integer                     :: n 
      integer                     :: status
      type(ESMF_Field)            :: varField
      type(ESMF_ArraySpec)        :: arrspec1
      REAL*8,POINTER :: unit(:)
      INTEGER                     :: jseed(10),nthRun
! SY: Begin parsetinit section, complete this later!!
# if 0
      REAL, allocatable       :: vardata(:) ! SY: For parsetinit
      REAL, allocatable    :: decvalues(:,:) ! SY: For parsetinit
# endif
! SY: End parsetinit section, complete this later!!
     
      data jseed/2,3,5,7,11,13,17,19,23,29/

! currently limited to one nest
      n = 1

      CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%decspaceAttribsFile,&
           label="SCEUA Decision Space Attributes File:",rc=status)
      call LIS_verify(status, 'SCEUA decision space attributes file: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%restart,&
           label="SCEUA start mode:",default=2) !SY: Check default value later
      CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%rfile,&
           label="SCEUA restart file:",rc=status)
      call LIS_verify(status, 'SCEUA restart file: not defined')

!     READ sceua_ctl%iDefault-independent CONTROL PARAMETERS and sceua_ctl%iDefault
      CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%MaxNFuncEvals,&
           label="SCEUA Max. Num. of Func. Evals. before Optimization Terminates:",&
           default=10000 )
      CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%NShufflesForMinChange,&
           label="SCEUA Num. of Shuffles to End Opt. if Crit. less than Min.:",&
           default=5)

      CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%MinChangeInNShuffles,&
           label="SCEUA Min. Frac. Crit. Change in Specified Shuffles to Cont. Opt.:",&
           default=0.0001)
      CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%NComplexes,&
           label="SCEUA Number of Optimization Complexes:",default=10)
      CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%iseed,&
           label="SCEUA Seed Value:",default=1969) 
      IF (sceua_ctl%iseed .EQ. 0) sceua_ctl%iseed = 1969
      IF (sceua_ctl%iseed .GT. 0) THEN
        nthRun = min(sceua_ctl%iseed,10)
      else
        nthRun = 1
      end if
!      if (nthRun .ne. 1) sceua_ctl%iseed = jseed(nthRun) ! SY: This line from original SCEUA code now corrected below to have prime number 2 as seed instead of 1 for nthRun == 1
      sceua_ctl%iseed = jseed(nthRun)
      write(LIS_logunit,*) 'Random initial seed value is  ',sceua_ctl%iseed

      CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%iDefault,&
           label="SCEUA Whether to User-specify the Control Parameters:",default=0)

      ftn = LIS_getNextUnitNumber()
      write(LIS_logunit,*) 'Reading decision space attributes ...', &
           sceua_ctl%decspaceAttribsFile
      open(ftn,file=trim(sceua_ctl%decspaceAttribsFile),status='old')
      read(ftn,*)
      read(ftn,*) sceua_ctl%nparam

      allocate(sceua_ctl%parmax(sceua_ctl%nparam)) !SY:Check later if deallocated!!
      allocate(sceua_ctl%parmin(sceua_ctl%nparam)) !SY:Check later if deallocated!!
      allocate(sceua_ctl%pardel(sceua_ctl%nparam)) !SY:Check later if deallocated!!
      allocate(vname(sceua_ctl%nparam))

!     READ THE PARAMETER BOUNDS
      read(ftn,*) 
      do i=1,sceua_ctl%nparam
         read(ftn,fmt='(a100)') vname(i)
         write(LIS_logunit,*) 'vname ',vname(i)
      enddo
      read(ftn,*) 
      read(ftn,*) (sceua_ctl%parmax(i),i=1,sceua_ctl%nparam)
      read(ftn,*) 
      read(ftn,*) (sceua_ctl%parmin(i),i=1,sceua_ctl%nparam)

      call LIS_releaseUnitNumber(ftn)
      write(LIS_logunit,*) 'Finished reading decision space attributes ..'

!     IF sceua_ctl%iDefault IS EQUAL TO 1, READ THE REMAINING SCE CONTROL PARAMETERS
      IF (sceua_ctl%iDefault .EQ. 1) THEN
!        CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%NPointsInComplex,&
!             "SCEUA Number of Points in a Complex:")
        IF (MOD(LIS_rc%nensem(n),sceua_ctl%NComplexes) .EQ. 0) THEN
          sceua_ctl%NPointsInComplexPlus4 = LIS_rc%nensem(n)/sceua_ctl%NComplexes
          sceua_ctl%NPointsInComplex = sceua_ctl%NPointsInComplexPlus4 - 4
        ELSE
         write(LIS_logunit,*) 'No. of ensembles is not a whole multiple of No. of complexes?'
         write(LIS_logunit,*) 'Stopping program ....'
         call LIS_endrun()
        END IF
        CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%NPointsInSubcomplex,&
             label="SCEUA Number of Points in a Subcomplex:",default=sceua_ctl%nparam+1)
        CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%NStepsBeforeShuffle,&
             label="SCEUA Num. of Evolution Steps before Shuffle for a Complex:",&
             default=sceua_ctl%NPointsInComplex)
!        CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%NMinComplexes,&
!             "Minimum to which Number of Complexes possibly Reduces:")
        CALL ESMF_ConfigGetAttribute(LIS_config,sceua_ctl%InitialPointFlag,&
             label="SCEUA Whether Include Initial Point in Population:",default=0)

! IF sceua_ctl%iDefault IS EQUAL TO 0, SET THE REMAINING SCE CONTROL PARAMETERS TO THE DEFAULT VALUES
      ELSE IF (sceua_ctl%iDefault .EQ. 0) THEN
        sceua_ctl%NPointsInComplex = 2*sceua_ctl%nparam + 1
        sceua_ctl%NPointsInComplexPlus4 = sceua_ctl%NPointsInComplex + 4 
!        IF (LIS_rc%nensem(n) .NE. sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex ) THEN ! SY
        IF (LIS_rc%nensem(n) .NE. sceua_ctl%NComplexes*sceua_ctl%NPointsInComplexPlus4 ) THEN ! SY
         write(LIS_logunit,*) 'LIS_rc%nensem(n) .NE. sceua_ctl%NComplexes*sceua_ctl%NPointsInComplexPlus4!'
         write(LIS_logunit,*) 'NOTE: sceua_ctl%NPointsInComplexPlus4=(2*sceua_ctl%nparam+1)+4 if sceua_ctl%IDefault is 0'
         write(LIS_logunit,*) 'Stopping program ....'
         call LIS_endrun()
        END IF
        sceua_ctl%NPointsInSubcomplex = sceua_ctl%nparam + 1
        sceua_ctl%NStepsBeforeShuffle = sceua_ctl%NPointsInComplex
!        sceua_ctl%NMinComplexes = sceua_ctl%NComplexes
        sceua_ctl%InitialPointFlag = 0

      ELSE
        write(LIS_logunit,*) 'sceua_ctl%iDefault should be either 1 or 0!'
        write(LIS_logunit,*) 'Stopping program ....'
        call LIS_endrun()
      END IF

!      sceua_ctl%npopsize = LIS_rc%nensem(n) ! SY

!  START OF CHECK IF THE SCE CONTROL PARAMETERS ARE VALID
      IF (sceua_ctl%NComplexes .LT. 1 .OR. sceua_ctl%NComplexes .GE. 1320) THEN
         write(LIS_logunit,*) 'ERROR: NUMBER OF COMPLEXES IN INITIAL POPULATION:'
         write(LIS_logunit,*) sceua_ctl%NComplexes,' IS NOT A VALID CHOICE'
         write(LIS_logunit,*) 'Should be between 1 and 1320'
         write(LIS_logunit,*) 'Stopping program ....'
         call LIS_endrun()
      END IF

      IF ((sceua_ctl%NShufflesForMinChange .LT. 0) .OR. &
              (sceua_ctl%NShufflesForMinChange .GE. 20)) THEN
         write(LIS_logunit,*) 'WARNING: THE NUMBER OF SHUFFLING LOOPS IN WHICH'
         write(LIS_logunit,*) 'CRITERION VALUE MUST CHANGE SHOULD BE BETWEEN 1 AND 19'
         write(LIS_logunit,*) 'SPECIFIED sceua_ctl%NShufflesForMinChange WAS', &
                              sceua_ctl%NShufflesForMinChange
         write(LIS_logunit,*) 'sceua_ctl%NShufflesForMinChange = 5 WILL BE USED INSTEAD.'
         sceua_ctl%NShufflesForMinChange=5
      END IF

!      IF ((sceua_ctl%NMinComplexes.LT. 1) .OR. (sceua_ctl%NMinComplexes .GT. ngs)) THEN
!         write(LIS_logunit,*) 'WARNING: THE MINIMUM NUMBER OF COMPLEXES SPECIFIED:'
!         write(LIS_logunit,*) sceua_ctl%NMinComplexes,'IS NOT A VALID CHOICE.'
!         write(LIS_logunit,*) 'SETTING IT TO DEFAULT' 
!         sceua_ctl%NMinComplexes = sceua_ctl%NComplexes
!      END IF

      IF ((sceua_ctl%NPointsInComplex .LT. 2) .OR. &
           (sceua_ctl%NPointsInComplex .GT. 1320/MAX(sceua_ctl%NComplexes,1))) THEN
! SY: NOTE: This was only a warning in original SCEUA code.
! SY:       The sceua_ctl%NPointsInComplex adjustment here would affect LIS_rc%nensem(n) 
         write(LIS_logunit,*) 'ERROR: THE NUMBER OF POINTS IN A COMPLEX:'
         write(LIS_logunit,*) sceua_ctl%NPointsInComplex,'IS NOT A VALID CHOICE.'
         write(LIS_logunit,*) 'SETTING IT TO DEFAULT NOW WOULD CHANGE LIS_rc%nensem(n)' 
         write(LIS_logunit,*) 'Stopping program ....'
         call LIS_endrun()
!         sceua_ctl%NPointsInComplex = 2*sceua_ctl%nparam+1
      END IF

      IF ((sceua_ctl%NPointsInSubcomplex .lt. 2) .or. &
          (sceua_ctl%NPointsInSubcomplex .gt. sceua_ctl%NPointsInComplex) .or. &
          (sceua_ctl%NPointsInSubcomplex .gt. 50)) THEN
         write(LIS_logunit,*) 'WARNING: THE NUMBER OF POINTS IN A SUB-COMPLEX:'
         write(LIS_logunit,*) sceua_ctl%NPointsInSubcomplex,'IS NOT A VALID CHOICE.'
         write(LIS_logunit,*) 'SETTING IT TO DEFAULT' 
         sceua_ctl%NPointsInSubcomplex = sceua_ctl%nparam + 1
      END IF

      if (sceua_ctl%NStepsBeforeShuffle .lt. 1) then
         write(LIS_logunit,*) 'WARNING: THE NUMBER OF EVOLUTION STEPS'
         write(LIS_logunit,*) 'TAKEN IN EACH COMPLEX BEFORE SHUFFLING:'
         write(LIS_logunit,*) sceua_ctl%NStepsBeforeShuffle,'IS NOT A VALID CHOICE.'
         write(LIS_logunit,*) 'SETTING IT TO DEFAULT' 
         sceua_ctl%NStepsBeforeShuffle = sceua_ctl%NPointsInComplex
      end if

      if (sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex .gt. 1320) then
! SY: NOTE: This was only a warning in original SCEUA code.
! SY:       The adjustments here of sceua_ctl%NComplexes and sceua_ctl%NPointsInComplex would affect LIS_rc%nensem(n) 
         write(LIS_logunit,*) 'ERROR: THE NUMBER OF POINTS IN INITIAL POPULATION:'
         write(LIS_logunit,*) sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex, &
                             ' EXCEED THE POPULATION LIMIT;'
         write(LIS_logunit,*) 'SETTING NGS TO 2, AND NPG, NPS AND NSPL TO DEFAULTS'
         write(LIS_logunit,*) 'NOW WOULD CHANGE LIS_rc%nensem(n)'
         write(LIS_logunit,*) 'Stopping program ....'
         call LIS_endrun()
!         sceua_ctl%NComplexes = 2
!         sceua_ctl%NPointsInComplex = 2*sceua_ctl%nparam + 1
!         sceua_ctl%NPointsInSubcomplex = sceua_ctl%nparam + 1
!         sceua_ctl%NStepsBeforeShuffle = sceua_ctl%NPointsInComplex
      end if
!  END OF CHECK IF THE SCE CONTROL PARAMETERS ARE VALID

      allocate(sceua_struc(LIS_rc%ntiles(n)/LIS_rc%nensem(n))) !SY:Check later if deallocated!

      call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
           rc=status)
      call LIS_verify(status)
      do i=1,sceua_ctl%nparam
         varField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
              name=trim(vname(i)),rc=status)
         call LIS_verify(status)
         call ESMF_StateAdd(LIS_decisionSpace, (/varField/),rc=status)
         call LIS_verify(status)
      enddo
      deallocate(vname)
 
!      call SCEUA_setup()

      ALLOCATE(unit(sceua_ctl%nparam))

!  INITIALIZE VARIABLES
      sceua_ctl%iEvolutionLoop = 0  

!  INITIALIZE RANDOM SEED TO A NEGATIVE INTEGER, and ran1_iff of ran1 
      sceua_ctl%iseed1 = -abs(sceua_ctl%iseed)
      sceua_ctl%ran1_iff = 0 ! SY: Can change later so that is read from restart file, not required though 
     
!  COMPUTE THE BOUND FOR PARAMETERS BEING OPTIMIZED
      do i = 1, sceua_ctl%nparam
        sceua_ctl%pardel(i) = sceua_ctl%parmax(i)-sceua_ctl%parmin(i)
        unit(i) = 1.0
      end do 

      allocate(sceua_ctl%isDecSpaceConvergedFlagArray(LIS_rc%ntiles(n)/LIS_rc%nensem(n)))
      allocate(sceua_ctl%isObjFuncConvergedFlagArray(LIS_rc%ntiles(n)/LIS_rc%nensem(n)))
      allocate(sceua_ctl%ParentParsets( sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex, sceua_ctl%nparam )) 
      allocate(sceua_ctl%ParentObjFuncs( sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex )) 
      allocate(sceua_ctl%complexParsets( sceua_ctl%NPointsInComplex, sceua_ctl%nparam )) 
      allocate(sceua_ctl%complexObjFuncs( sceua_ctl%NPointsInComplex )) 
!SY:  Check later if above arrays deallocated!

      sceua_ctl%maxObjFunc = 1.0E15

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) ! SY: Begin Looping for each spatial tile 

! BEGIN allocation of required algorithmic variable arrays (for each spatial tile)
         allocate(sceua_struc(t)%parsets( LIS_rc%nensem(n), sceua_ctl%nparam ))
         allocate(sceua_struc(t)%parset(sceua_ctl%nparam))
         allocate(sceua_struc(t)%bestparsetinloop(sceua_ctl%nparam))
         allocate(sceua_struc(t)%worstparsetinloop(sceua_ctl%nparam))
         allocate(sceua_struc(t)%objfuncs(LIS_rc%nensem(n)))
         allocate(sceua_struc(t)%parstd(sceua_ctl%nparam))
         allocate(sceua_struc(t)%lcses(sceua_ctl%NComplexes,sceua_ctl%NPointsInSubcomplex))
         allocate(sceua_struc(t)%complexesParsets(sceua_ctl%NComplexes,sceua_ctl%NPointsInComplex,sceua_ctl%nparam))
         allocate(sceua_struc(t)%subcomplexesParsets(sceua_ctl%NComplexes,sceua_ctl%NPointsInSubcomplex,sceua_ctl%nparam))
         allocate(sceua_struc(t)%complexesObjFuncs(sceua_ctl%NComplexes,sceua_ctl%NPointsInComplex))
         allocate(sceua_struc(t)%subcomplexesObjFuncs(sceua_ctl%NComplexes,sceua_ctl%NPointsInSubcomplex))
         allocate(sceua_struc(t)%parsetinit(sceua_ctl%nparam))
         allocate(sceua_struc(t)%ibound(sceua_ctl%NComplexes))
! END allocation of required algorithmic variable arrays (for each spatial tile)
!SY:Check later if ALL ABOVE deallocated!

         sceua_struc(t)%BestObjFuncHistory = 0.0 

      enddo ! SY: End Looping for each spatial tile 

!     Begin section assigning initial point if InitialPointFlag EQUALS 1, complete this section later!!    
# if 0
      allocate( decvalues(sceua_ctl%nparam, LIS_rc%ntiles(n)) )
#endif
!     End section assigning initial point if InitialPointFlag EQUALS 1, complete this section later!!    

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) ! SY: Begin Looping for each spatial tile 

!  GENERATE 1ST POINT OF INITIAL SET OF LIS_rc%nensem(n) POINTS IN PARAM. SPACE
!  IF sceua_ctl%InitialPointFlag EQUALS 1, 
!  SET sceuastruc%parsets(.,1) TO INITIAL POINT. SY: ALREADY @ INITIAL PT.! 
         if (sceua_ctl%InitialPointFlag .eq. 1) then
           do j = 1, sceua_ctl%nparam
              sceua_struc(t)%parsetinit(j) = sceua_struc(t)%parsets(1,j) 
           end do
!  ELSE, GENERATE A POINT RANDOMLY AND GIVE IT'S VALUE TO sceuastruc%parsets(1,:)
         else
! SY: Generate a point here
           call getpnt(sceua_ctl%nparam,1,sceua_ctl%iseed1, &
                       sceua_struc(t)%parset,sceua_ctl%parmin, &
                       sceua_ctl%parmax,unit,sceua_ctl%parmin)
           do j=1, sceua_ctl%nparam
              sceua_struc(t)%parsets(1,j) = sceua_struc(t)%parset(j)
           end do
         end if
     
!  GENERATE (sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex - 1) RAND. POINTS DISTRIB. UNIFORMLY IN PARAM. SPACE
         do i = 2, sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex
            call getpnt(sceua_ctl%nparam,1,sceua_ctl%iseed1, &
                        sceua_struc(t)%parset,sceua_ctl%parmin, &
                        sceua_ctl%parmax,unit,sceua_ctl%parmin)
            do j = 1, sceua_ctl%nparam
               sceua_struc(t)%parsets(i,j) = sceua_struc(t)%parset(j)
            end do
         end do

! SY: ASSIGN AN ARBITRATY PARENT PAR. SET VALUE TO THE CHILDREN PAR. SETS AS FILLER
         do i = sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex+1,LIS_rc%nensem(n)
            do j = 1, sceua_ctl%nparam
               sceua_struc(t)%parsets(i,j) = sceua_struc(t)%parset(j)
            end do
         end do

      enddo ! SY: End Looping for each spatial tile 

      DEALLOCATE(unit)

    end subroutine SCEUAOpt_init

!BOP
! !ROUTINE: SCEUAOpt_setup
! \label{SCEUAOpt_setup}
!
! !INTERFACE: SCEUAOpt_setup
    subroutine SCEUAOpt_setup()
! !USES: 
      use LIS_coreMod,         only : LIS_rc
!      use LIS_optUEMod, only : LIS_decisionSpace ! SY
      use LIS_optUEMod,        only : LIS_decisionSpace, LIS_ObjectiveFunc 
      use LIS_logMod,          only : LIS_logunit !, LIS_endrun ! SY
! 
! !DESCRIPTION: 
!   This subroutine initializes the required memory structures,  
!   creates an initial random population, and calculates the 
!   objective functions values, all for each grid point. 
!   Note that this setup actually runs the LSM
!   because objective function value is evaluated. 
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP
      integer        :: i,j,l !,k ! SY: Not really used
      integer        :: n
      integer        :: t
      real,  allocatable :: InitialObjFunc(:) ! SY: Initial/default objfuncs value, i.e., objfuncs(1) if sceua_ctl%InitialPointFlag is 1, across all spatial tiles
      type(ESMF_Field)   :: objfuncField
      real, pointer      :: objfuncValue(:)
      REAL*8             :: denomi
      REAL*8             :: timeou

#if 0 
! Currently supporting only one nest. 
      n = 1

      if(sceua_ctl%restart.eq.1) then !restart run
         call readSCEUArestart()
      endif

! translate the updated decision space to the generic objects, required for SCEUA_evaluateObjFunc below
      call setoptuetypedecspace(LIS_rc%optuetype) 
      write(LIS_logunit,*) 'Algorithmic decision space updated & assigned to generic objects'

! SY: Function(LSM) evaluations comes here after readSCEUArestart()!!
      call SCEUA_evaluateObjFunc() 

      if(sceua_ctl%restart .NE. 1) then ! SY: NOT a restart run
         sceua_ctl%SumReplacedFuncEvals = sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex
         sceua_ctl%SumNewMemberFuncEvals = sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex
         sceua_ctl%SumAllFuncEvals = LIS_rc%nensem(n)
      endif
!         if (sceua_ctl%SumNewMemberFuncEvals .ge. sceua_ctl%MaxNFuncEvals) go to 9000 
            !if (sceua_ctl%SumNewMemberFuncEvals .ge. sceua_ctl%MaxNFuncEvals) then 
              !go to 45 
            !end if

      write(LIS_logunit,*) 'Running SCEUA Evolution loop: ', &
                           sceua_ctl%iEvolutionLoop, &
                           ', sceua_ctl%SumReplacedFuncEvals is ', &
                           sceua_ctl%SumReplacedFuncEvals, &
                           ', sceua_ctl%SumNewMemberFuncEvals is ', &
                           sceua_ctl%SumNewMemberFuncEvals, &
                           ', sceua_ctl%SumAllFuncEvals is ', &
                           sceua_ctl%SumAllFuncEvals

      allocate(InitialObjFunc( LIS_rc%ntiles(n)/LIS_rc%nensem(n)  ))
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) ! SY: Begin Looping for each spatial tile 

         do i = 1, sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex
            do j = 1, sceua_ctl%nparam
               sceua_ctl%ParentParsets(i,j) = sceua_struc(t)%parsets(i,j) ! SY
            end do
            sceua_ctl%ParentObjFuncs(i) = sceua_struc(t)%objfuncs(i)
         end do

      if(sceua_ctl%restart .NE. 1) then ! SY: NOT a restart run; restart run files are written after following sorting

!  ARRANGE THE POINTS IN ORDER OF INCREASING FUNCTION VALUE 
         call sort(sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex, &
                   sceua_ctl%nparam, sceua_ctl%ParentParsets, sceua_ctl%ParentObjFuncs) ! SY: This statement is labeled 45 in original SCE-UA code

         do i = 1, sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex
            do j = 1, sceua_ctl%nparam
               sceua_struc(t)%parsets(i,j) = sceua_ctl%ParentParsets(i,j) 
            end do
            sceua_struc(t)%objfuncs(i) = sceua_ctl%ParentObjFuncs(i) 
         end do

      endif


!  RECORD THE BEST AND WORST POINTS
         do j = 1, sceua_ctl%nparam
           sceua_struc(t)%bestparsetinloop(j) = sceua_struc(t)%parsets(1,j)
           sceua_struc(t)%worstparsetinloop(j) = &
                sceua_struc(t)%parsets(sceua_ctl%NComplexes* &
                                       sceua_ctl%NPointsInComplex,j)
         end do
         sceua_struc(t)%bestobjfuncinloop = sceua_struc(t)%objfuncs(1)
         sceua_struc(t)%worstobjfuncinloop = &
                         sceua_struc(t)%objfuncs(sceua_ctl%NComplexes* &
                                                 sceua_ctl%NPointsInComplex)

!  COMPUTE THE PARAMETER RANGE FOR THE INITIAL POPULATION
         call parstt(sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex, &
                   sceua_ctl%nparam, sceua_ctl%ParentParsets, sceua_struc(t)%parstd, &
                   sceua_ctl%pardel,sceua_struc(t)%gmparrng, &
                   sceua_struc(t)%isDecSpaceConvergedFlag)
         sceua_ctl%isDecSpaceConvergedFlagArray(t) = &
                       sceua_struc(t)%isDecSpaceConvergedFlag

         IF (sceua_ctl%InitialPointFlag .EQ. 1) THEN 
           InitialObjFunc(t) = sceua_ctl%ParentObjFuncs(1)  
         END IF 

! SY: INTRODUCING HERE THE Obj. Func. CONVERGENCE CRITERION INITIALIZATION ALSO
         sceua_struc(t)%isObjFuncConvergedFlag = 0

!  COMPUTE THE COUNT ON SUCCESSIVE LOOPS W/O FUNCTION IMPROVEMENT
         sceua_struc(t)%BestObjFuncHistory(20) = &
                         sceua_struc(t)%bestobjfuncinloop 
         write(LIS_logunit,*) 'BestObjFuncHistory is ', &
                              sceua_struc(t)%BestObjFuncHistory
         IF ( sceua_ctl%iEvolutionLoop .GE. &
                       (sceua_ctl%NShufflesForMinChange+1)) THEN 
           denomi = dabs(sceua_struc(t)%BestObjFuncHistory(20- &
                           sceua_ctl%NShufflesForMinChange) + &
                         sceua_struc(t)%BestObjFuncHistory(20)) / 2.
           write(LIS_logunit,*) 'denomi is ', denomi
           timeou = dabs(sceua_struc(t)%BestObjFuncHistory(20- &
                           sceua_ctl%NShufflesForMinChange) - &
                         sceua_struc(t)%BestObjFuncHistory(20)) / denomi
           if (timeou .lt. sceua_ctl%MinChangeInNShuffles) then
             sceua_struc(t)%isObjFuncConvergedFlag = 1
           end if
         END IF
         sceua_ctl%isObjFuncConvergedFlagArray(t) = &
                       sceua_struc(t)%isObjFuncConvergedFlag
         if (sceua_struc(t)%isObjFuncConvergedFlag .EQ. 0 ) then
           do l = 1, 19
             sceua_struc(t)%BestObjFuncHistory(l) = &
                   sceua_struc(t)%BestObjFuncHistory(l+1)
           end do 
         end if

      enddo ! SY: End Looping for each spatial tile 

      IF (sceua_ctl%InitialPointFlag .EQ. 1) THEN 
!  PRINT THE INITIAL POINT'S CRITERION VALUE
         write(LIS_logunit,*) "Array of Initial Point's criterion value "
         write(LIS_logunit,*) 'across all spatial tiles is: '
         write(LIS_logunit,*) InitialObjFunc
      END IF 

      deallocate(InitialObjFunc) ! SY

      IF ( sum(sceua_ctl%isDecSpaceConvergedFlagArray) .EQ. &
                    LIS_rc%ntiles(n)/LIS_rc%nensem(n) ) THEN 
         sceua_ctl%isDecSpaceConvergAllSpatialTiles = 1 ! SY
      ELSE
         sceua_ctl%isDecSpaceConvergAllSpatialTiles = 0 ! SY
      END IF

      IF ( sum(sceua_ctl%isObjFuncConvergedFlagArray) .EQ. &
                    LIS_rc%ntiles(n)/LIS_rc%nensem(n) ) THEN 
         sceua_ctl%isObjFuncConvergAllSpatialTiles = 1 ! SY
      ELSE
         sceua_ctl%isObjFuncConvergAllSpatialTiles = 0 ! SY
      END IF
#endif
    end subroutine SCEUAOpt_setup

!BOP
! 
! !ROUTINE: SCEUAOpt_checkOptStopCriterion
! \label{SCEUAOpt_checkOptStopCriterion}
! 
! !INTERFACE: 
    subroutine SCEUAOpt_checkOptStopCriterion(check)
! !USES: 
      use LIS_logMod,          only : LIS_logunit
! !INPUT/OUTPUT PARAMETERS:: 
      logical, intent(INOUT) :: check
! 
! !DESCRIPTION: 
!  This routine checks to see if the stopping criterion for SCEUA 
!  is reached. In this case, the routine simply checks to see if 
!  [1] specified number of function (LSM) evaluations is reached, or
!  [2] decision space converged for all spatial tiles, or  
!  [3] objective function converged for all spatial tiles
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP

      if( (sceua_ctl%SumNewMemberFuncEvals .ge. sceua_ctl%MaxNFuncEvals) .OR. &
          (sceua_ctl%isObjFuncConvergAllSpatialTiles .EQ. 1) .OR. & ! SY: Uncomment later!!!
          (sceua_ctl%isDecSpaceConvergAllSpatialTiles .EQ. 1)  ) then 
         write(LIS_logunit,*) 'SCEUAOpt_checkOptStopCriterion is true!'
         write(LIS_logunit,*) 'isDecSpaceConvergAllSpatialTiles is', &
                              sceua_ctl%isDecSpaceConvergAllSpatialTiles
         write(LIS_logunit,*) 'isObjFuncConvergAllSpatialTiles is', &
                              sceua_ctl%isObjFuncConvergAllSpatialTiles
         check = .true.
      else
         check = .false.
      endif

    end subroutine SCEUAOpt_checkOptStopCriterion

!BOP
! !ROUTINE: SCEUAOpt_run
! \label{SCEUAOpt_run}
! 
! !INTERFACE: 
    subroutine SCEUAOpt_run()
! !USES: 
      use LIS_coreMod,   only : LIS_rc
      use LIS_logMod,    only : LIS_logunit

! !DESCRIPTION: 
!   This routine performs the run steps for the SCEUA:  
!   assigning points into complexes, random selection of subcomplexes,
!   using subcomplexes to generate new points, evaluating objective 
!   function values, replacing new subcomplexes into complexes, 
!   sorting complexes, and replacing complexes back into the population.
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP  
      integer                  :: t
      integer                  :: n
      integer                  :: iComplex, k1, k2, iStepsBeforeShuffle
      integer                  :: k, j, ii, l
      REAL*8                   :: denomi
      REAL*8                   :: timeou

#if 0 
      n = 1
      sceua_ctl%iEvolutionLoop = sceua_ctl%iEvolutionLoop + 1
      write(LIS_logunit,*) 'Running SCEUA Evolution loop: ', &
                           sceua_ctl%iEvolutionLoop
    
!  BEGIN LOOP ON COMPLEXES
      do iComplex = 1, sceua_ctl%NComplexes

!  BEGIN LOOP ON SPATIAL TILES
        do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)

!  ASSIGN POINTS INTO COMPLEXES
          do k1 = 1, sceua_ctl%NPointsInComplex
            k2 = (k1-1) * sceua_ctl%NComplexes + iComplex
            do j = 1, sceua_ctl%nparam
              sceua_struc(t)%complexesParsets(iComplex,k1,j) = &
                                sceua_struc(t)%parsets(k2,j)
            end do
            sceua_struc(t)%complexesObjFuncs(iComplex,k1) = &
                              sceua_struc(t)%objfuncs(k2)
          end do

        enddo
!  END LOOP ON SPATIAL TILES

      end do
!  END LOOP ON COMPLEXES

!  BEGIN LOOP - RANDOM SELECTION OF SUB-COMPLEXES ---------------
!  SY: NOTE: This was an INNER LOOP in the original SCEUA code!!
      do iStepsBeforeShuffle = 1, sceua_ctl%NStepsBeforeShuffle 

!  BEGIN LOOP ON COMPLEXES
        do iComplex = 1, sceua_ctl%NComplexes

!  BEGIN LOOP ON SPATIAL TILES
          do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)

            call SCEUA_SubcomplexSelection(iComplex,t)

!  USE THE SUB-COMPLEX TO GENERATE NEW POINT(S)
!           call cce(nopt,nps,s,sf,bl,bu,xnstd,icall,maxn,iseed1) ! SY
            call SCEUA_SubcomplexCreateChildren(iComplex, t) ! SY

          enddo
!  END LOOP ON SPATIAL TILES

        end do
!  END LOOP ON COMPLEXES

! translate the updated decision space to the generic objects, required for SCEUA_evaluateObjFunc below 
        call setoptuetypedecspace(LIS_rc%optuetype) 

        call SCEUA_evaluateObjFunc() 

        sceua_ctl%SumReplacedFuncEvals = sceua_ctl%SumReplacedFuncEvals + &
                         sceua_ctl%NComplexes
        sceua_ctl%SumNewMemberFuncEvals = sceua_ctl%SumNewMemberFuncEvals + &
                         sceua_ctl%NComplexes*4
        sceua_ctl%SumAllFuncEvals = sceua_ctl%SumAllFuncEvals + &
                         LIS_rc%nensem(n)

        write(LIS_logunit,*) 'In iStepsBeforeShuffle # ', &
                             iStepsBeforeShuffle, &
                             ',sceua_ctl%SumReplacedFuncEvals is ', &
                             sceua_ctl%SumReplacedFuncEvals, &
                             ',sceua_ctl%SumNewMemberFuncEvals is ', &
                             sceua_ctl%SumNewMemberFuncEvals, &
                             ', sceua_ctl%SumAllFuncEvals is ', &
                             sceua_ctl%SumAllFuncEvals

!  BEGIN RE-LOOP ON COMPLEXES
        do iComplex = 1, sceua_ctl%NComplexes

!  BEGIN RE-LOOP ON SPATIAL TILES
          do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)

            call SCEUA_SubcomplexNewGeneration(iComplex, t) ! SY

!  REPLACE THE NEW SUB-COMPLEX INTO THE COMPLEX
            do k = 1, sceua_ctl%NPointsInSubcomplex
              do j = 1, sceua_ctl%nparam
                sceua_struc(t)%complexesParsets(iComplex, &
                                    sceua_struc(t)%lcses(iComplex,k),j) = &
                           sceua_struc(t)%subcomplexesParsets(iComplex,k,j)
              end do
              sceua_struc(t)%complexesObjFuncs(iComplex, &
                                    sceua_struc(t)%lcses(iComplex,k)) = &
                         sceua_struc(t)%subcomplexesObjFuncs(iComplex,k)
            end do

            do ii = 1, sceua_ctl%NPointsInComplex
              do j = 1, sceua_ctl%nparam
                sceua_ctl%complexParsets(ii,j) = &
                           sceua_struc(t)%complexesParsets(iComplex, ii, j)
              end do
              sceua_ctl%complexObjFuncs(ii) = &
                         sceua_struc(t)%complexesObjFuncs(iComplex, ii)
            end do    

!  SORT THE COMPLEX POINTS
            call sort(sceua_ctl%NPointsInComplex,sceua_ctl%nparam, &
                      sceua_ctl%complexParsets,sceua_ctl%complexObjFuncs)

            do ii = 1, sceua_ctl%NPointsInComplex
              do j = 1, sceua_ctl%nparam
                sceua_struc(t)%complexesParsets(iComplex, ii, j) = &
                           sceua_ctl%complexParsets(ii,j) 
              end do
              sceua_struc(t)%complexesObjFuncs(iComplex, ii) = &
                         sceua_ctl%complexObjFuncs(ii) 
            end do    

!  REPLACE THE NEW COMPLEX INTO ORIGINAL ARRAY x(.,.)
!  SY: NOTE THAT THIS REPLACES NEW PARENTS HERE NOW, 
!      CHILDREN ALREADY REPLACED IN SCEUA_SubcomplexNewGeneration 
            do k1 = 1, sceua_ctl%NPointsInComplex
              k2 = (k1-1) * sceua_ctl%NComplexes + iComplex
              do j = 1, sceua_ctl%nparam
                sceua_struc(t)%parsets(k2,j) = &
                        sceua_struc(t)%complexesParsets(iComplex,k1,j) 
              end do
              sceua_struc(t)%objfuncs(k2) = &
                      sceua_struc(t)%complexesObjFuncs(iComplex,k1)
            end do 

          enddo
!  END RE-LOOP ON SPATIAL TILES

        end do
!  END RE-LOOP ON COMPLEXES

      end do
!  END LOOP - RANDOM SELECTION OF SUB-COMPLEXES ---------------

!  BEGIN RE-RE-LOOP ON SPATIAL TILES
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)

        do ii = 1, sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex
          do j = 1, sceua_ctl%nparam
            sceua_ctl%ParentParsets(ii,j) = sceua_struc(t)%parsets(ii,j) ! SY
          end do
          sceua_ctl%ParentObjFuncs(ii) = sceua_struc(t)%objfuncs(ii)
        end do

!  SORT THE PARENT MEMBERS OF ENSEMBLE
        call sort(sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex, &
                  sceua_ctl%nparam,sceua_ctl%ParentParsets,sceua_ctl%ParentObjFuncs)

        do ii = 1, sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex
          do j = 1, sceua_ctl%nparam
            sceua_struc(t)%parsets(ii,j) = sceua_ctl%ParentParsets(ii,j)
          end do
          sceua_struc(t)%objfuncs(ii) = sceua_ctl%ParentObjFuncs(ii)
        end do


!  RECORD THE BEST AND WORST POINTS
        do j = 1, sceua_ctl%nparam
          sceua_struc(t)%bestparsetinloop(j) = sceua_struc(t)%parsets(1,j)
          sceua_struc(t)%worstparsetinloop(j) = &
               sceua_struc(t)%parsets(sceua_ctl%NComplexes* &
                                      sceua_ctl%NPointsInComplex,j)
        end do
        sceua_struc(t)%bestobjfuncinloop = sceua_struc(t)%objfuncs(1)
        sceua_struc(t)%worstobjfuncinloop = &
             sceua_struc(t)%objfuncs(sceua_ctl%NComplexes* &
                                     sceua_ctl%NPointsInComplex)

!  TEST THE POPULATION FOR PARAMETER CONVERGENCE 
        call parstt(sceua_ctl%NComplexes*sceua_ctl%NPointsInComplex, &
                  sceua_ctl%nparam, sceua_ctl%ParentParsets, sceua_struc(t)%parstd, &
                  sceua_ctl%pardel,sceua_struc(t)%gmparrng, &
                  sceua_struc(t)%isDecSpaceConvergedFlag)
        sceua_ctl%isDecSpaceConvergedFlagArray(t) = &
                      sceua_struc(t)%isDecSpaceConvergedFlag

!  COMPUTE THE COUNT ON SUCCESSIVE LOOPS W/O FUNCTION IMPROVEMENT
        sceua_struc(t)%BestObjFuncHistory(20) = &
                        sceua_struc(t)%bestobjfuncinloop 
        write(LIS_logunit,*) 'BestObjFuncHistory is ', &
                              sceua_struc(t)%BestObjFuncHistory
        IF ( sceua_ctl%iEvolutionLoop .GE. &
                      (sceua_ctl%NShufflesForMinChange+1)) THEN 
          denomi = dabs(sceua_struc(t)%BestObjFuncHistory(20- &
                          sceua_ctl%NShufflesForMinChange) + &
                        sceua_struc(t)%BestObjFuncHistory(20)) / 2.
          write(LIS_logunit,*) 'denomi is ', denomi
          timeou = dabs(sceua_struc(t)%BestObjFuncHistory(20- &
                          sceua_ctl%NShufflesForMinChange) - &
                        sceua_struc(t)%BestObjFuncHistory(20)) / denomi
          if (timeou .lt. sceua_ctl%MinChangeInNShuffles) then
            sceua_struc(t)%isObjFuncConvergedFlag = 1
          end if
        END IF
        sceua_ctl%isObjFuncConvergedFlagArray(t) = &
                      sceua_struc(t)%isObjFuncConvergedFlag
        if (sceua_struc(t)%isObjFuncConvergedFlag .EQ. 0 ) then
          do l = 1, 19 
            sceua_struc(t)%BestObjFuncHistory(l) = &
                  sceua_struc(t)%BestObjFuncHistory(l+1)
          end do
        end if

      enddo
!  END RE-RE-LOOP ON SPATIAL TILES

      IF ( sum(sceua_ctl%isDecSpaceConvergedFlagArray) .EQ. &
                    LIS_rc%ntiles(n)/LIS_rc%nensem(n) ) THEN  
         sceua_ctl%isDecSpaceConvergAllSpatialTiles = 1 ! SY
      ELSE    
         sceua_ctl%isDecSpaceConvergAllSpatialTiles = 0 ! SY
      END IF  

      IF ( sum(sceua_ctl%isObjFuncConvergedFlagArray) .EQ. &
                    LIS_rc%ntiles(n)/LIS_rc%nensem(n) ) THEN  
         sceua_ctl%isObjFuncConvergAllSpatialTiles = 1 ! SY
      ELSE    
         sceua_ctl%isObjFuncConvergAllSpatialTiles = 0 ! SY
      END IF  

! translate the updated decision space to the generic objects, required for LIS_decisionSpace in writeSCEUArestart part of SCEUA_printPopulation below 
      call setoptuetypedecspace(LIS_rc%optuetype)

      call SCEUA_printPopulation() ! SY: NOTE: This does not come at the beginning right after evaulating fitness as in GA. It not comes at the end after evaluating Obj. Func. and ALSO AFTER sort to.
#endif
    end subroutine SCEUAOpt_run

!BOP
! 
! !ROUTINE: SCEUA_evaluateObjFunc
!  \label{SCEUA_evaluateObjFunc}
! 
! !INTERFACE: 
    subroutine SCEUA_evaluateObjFunc()
! !USES:   
    use LIS_coreMod,         only : LIS_rc
    use LIS_optUEMod,        only : LIS_ObjectiveFunc, LIS_feasibleSpace
    use LIS_logMod,          only : LIS_logunit, LIS_verify
! 
! !DESCRIPTION: 
!   This method computes the objective function values for each parameter set,
!   Note that elitism is always in effect since only the worst parameter set 
!   (and not the best) in the subcomplex is always replaced.
! 
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP
    integer            :: n 
    type(ESMF_Field)   :: objfuncField 
    type(ESMF_Field)   :: feasField 
    integer            :: status
    real, pointer      :: objfuncValue(:)
    integer, pointer   :: mod_flag(:)
    integer            :: t
    integer            :: m

    call evaluateobjfunction(LIS_rc%optuetype)

    call ESMF_StateGet(LIS_ObjectiveFunc,"Objective Function Value",objfuncField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(objfuncField,localDE=0, farrayPtr=objfuncValue,rc=status)
    call LIS_verify(status)

    call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
    call LIS_verify(status)

    n = 1
    do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
       do m=1,LIS_rc%nensem(n)
          sceua_struc(t)%objfuncs(m) = objfuncValue((t-1)*LIS_rc%nensem(n)+m)
!penalize infeasible solutions
          if(mod_flag((t-1)*LIS_rc%nensem(n)+m).eq.1) then
             sceua_struc(t)%objfuncs(m) = sceua_ctl%maxObjFunc
          endif
       end do
    enddo
!hack
    call resetobjectivefunctype(LIS_rc%objfuncmethod)

  end subroutine SCEUA_evaluateObjFunc


!BOP
! 
! !ROUTINE: SCEUA_SubcomplexSelection
!  \label{SCEUA_SubcomplexSelection}
! 
! !INTERFACE: 
    subroutine SCEUA_SubcomplexSelection(iComplex,t)
! !USES:   
      use LIS_logMod,          only : LIS_logunit ! SY: del. later
!
      implicit none
! !INPUT PARAMETERS: 
!
      integer, intent(IN)      :: iComplex
      integer, intent(IN)      :: t
!
! !DESCRIPTION: 
!    This routine selects the subcomplex from each complex
!    according to a linear probability distribution
! 
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP
      integer                  :: n
      integer                  :: k, lpos, k1, j, lposcounter ! SY
      integer, allocatable         :: lcs(:) ! Locations of subcomplex in complex
      real*8                   :: rand

      ALLOCATE(lcs(sceua_ctl%NPointsInSubcomplex)) 

!  CHOOSE A SUB-COMPLEX (sceua_ctl%NPointsInSubcomplex points) 
!  ACCORDING TO A LINEAR PROBABILITY DISTRIBUTION
      if (sceua_ctl%NPointsInSubcomplex .eq. sceua_ctl%NPointsInComplex) then
        do k = 1, sceua_ctl%NPointsInSubcomplex
          lcs(k) = k
          sceua_struc(t)%lcses(iComplex,k) = lcs(k)
        end do  
      else
        rand = ran1(sceua_ctl%iseed1)
        lcs(1) = 1 + dint(sceua_ctl%NPointsInComplex + 0.5 - &
                          dsqrt( (sceua_ctl%NPointsInComplex+.5)**2 - &
                                 sceua_ctl%NPointsInComplex * &
                                 (sceua_ctl%NPointsInComplex+1) * rand )) 
        do k = 2, sceua_ctl%NPointsInSubcomplex
          DO
!   60       rand = ran1(iseed1) ! SY
            rand = ran1(sceua_ctl%iseed1) ! SY
            lpos = 1 + dint(sceua_ctl%NPointsInComplex + 0.5 - &
                            dsqrt((sceua_ctl%NPointsInComplex+.5)**2 - &
                                  sceua_ctl%NPointsInComplex * &
                                  (sceua_ctl%NPointsInComplex+1) * rand )) 
            lposcounter = 0 ! SY
            do k1 = 1, k-1
!                if (lpos .eq. lcs(k1)) go to 60 ! SY
!             SY: Begin
              if (lpos .eq. lcs(k1)) then
                lposcounter=1
                EXIT
              end if
!             SY: End
            end do  
            if (lposcounter .eq. 0) EXIT ! SY
          END DO
          lcs(k) = lpos
        end do  

!  ARRANGE THE SUB-COMPLEX IN ORDER OF INCREASING FUNCTION VALUE
        call sort1(sceua_ctl%NPointsInSubcomplex,lcs) 
        do k = 1, sceua_ctl%NPointsInSubcomplex
          sceua_struc(t)%lcses(iComplex,k) = lcs(k)
        end do  
      end if  

!   CREATE THE SUB-COMPLEX ARRAYS
      do k = 1, sceua_ctl%NPointsInSubcomplex
        do j = 1, sceua_ctl%nparam
           sceua_struc(t)%subcomplexesParsets(iComplex,k,j) = &
                  sceua_struc(t)%complexesParsets(iComplex,lcs(k),j)
        end do
        sceua_struc(t)%subcomplexesObjFuncs(iComplex,k) = &
               sceua_struc(t)%complexesObjFuncs(iComplex,lcs(k))
      end do

      DEALLOCATE(lcs)
    end subroutine SCEUA_SubcomplexSelection


!BOP
!
! !ROUTINE: SCEUA_SubcomplexCreateChildren
! \label{SCEUA_SubcomplexCreateChildren}
!
! !INTERFACE:
    subroutine SCEUA_SubcomplexCreateChildren(iComplex, t)
! !USES:   
      use LIS_logMod,          only : LIS_logunit ! SY: del. later
!
      implicit none
!
! !INPUT PARAMETERS: 
      integer, intent(IN)           :: iComplex
      integer, intent(IN)           :: t
!
! !DESCRIPTION:
!   SY: This method creates 4 potential children for the considered 
!   subcomplex: 1 using reflection, 1 using contraction, and 2 using 
!   random mutations in the desired parameter space region. Only one of 
!   four would finally replace the worst parent later.
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP
      INTEGER           :: m,n,i,j
      REAL*8            :: alpha, beta
      REAL*8            :: sb(sceua_ctl%nparam) ! SY: Best sub-complex parent 
      REAL*8            :: sw(sceua_ctl%nparam) ! SY: Worst sub-complex parent
      REAL*8            :: ce(sceua_ctl%nparam) ! SY: Sub-complex parents' centroid
      INTEGER           :: ibound

!  EQUIVALENCE OF VARIABLES FOR READABILTY OF CODE
      n = sceua_ctl%NPointsInSubcomplex
      m = sceua_ctl%nparam
      alpha = 1.0
      beta = 0.5

!  IDENTIFY THE WORST POINT wo OF THE SUB-COMPLEX s
!  COMPUTE THE CENTROID ce OF THE REMAINING POINTS
!  COMPUTE step, THE VECTOR BETWEEN wo AND ce
!  IDENTIFY THE WORST FUNCTION VALUE fw
      do j = 1, m
        sb(j) = sceua_struc(t)%subcomplexesParsets(iComplex,1,j) ! SY
        sw(j) = sceua_struc(t)%subcomplexesParsets(iComplex,n,j) ! SY
        ce(j) = 0.0
        do i = 1, n-1
          ce(j) = ce(j) + sceua_struc(t)%subcomplexesParsets(iComplex,i,j) ! SY
        end do
        ce(j) = ce(j)/dble(n-1)
      end do

!  COMPUTE FIRST NEW POINT BY TRYING A REFLECTION STEP
      do j = 1, m
        sceua_struc(t)%parset(j) = ce(j) + alpha * (ce(j) - sw(j)) 
        sceua_struc(t)%parsets( sceua_ctl%NComplexes * &
                                sceua_ctl%NPointsInComplex + &
                                (iComplex-1)*4+1, j ) = &
                                      sceua_struc(t)%parset(j)
      end do

!  CHECK IF 1st CHILD SATISFIES ALL CONSTRAINTS
      call chkcst(sceua_ctl%nparam,sceua_struc(t)%parset, &
                  sceua_ctl%parmin,sceua_ctl%parmax, ibound)
      sceua_struc(t)%ibound(iComplex) = ibound

!  CHOOSE THE 2nd NEW POINT AT RANDOM WITHIN FEASIBLE REGION ACCORDING TO
!  A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
!  AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
      call getpnt(sceua_ctl%nparam,2,sceua_ctl%iseed1, &
                  sceua_struc(t)%parset,sceua_ctl%parmin, &
                  sceua_ctl%parmax, sceua_struc(t)%parstd, sb) ! SY
      do j = 1, m
        sceua_struc(t)%parsets( sceua_ctl%NComplexes * &
                                sceua_ctl%NPointsInComplex + &
                                (iComplex-1)*4+2, j ) = &
                                      sceua_struc(t)%parset(j)
      end do
!  IF 1st CHILD DOES NOT SATISFY ALL CONSTRAINTS, ASSIGN 2nd POINT VALUE
!  TO IT: IT WOULD NOT BE USED FURTHER IN OPTIMIZATION, ONLY ibound WOULD BE.
      if (sceua_struc(t)%ibound(iComplex) .ge. 1) then
        do j = 1, m
          sceua_struc(t)%parsets( sceua_ctl%NComplexes * &
                                  sceua_ctl%NPointsInComplex + &
                                  (iComplex-1)*4+1, j ) = &
                                      sceua_struc(t)%parset(j)
        end do
      end if

!  FOR 3rd NEW POINT, TRY A CONTRACTION STEP
      do j = 1, m
!        sceua_struc(t)%parset(j) = ce(j) - beta * (ce(j) - sw(j)) ! SY: see below
        sceua_struc(t)%parset(j) = ( ( ce(j)*dble(n-1) )+sw(j) )/dble(n) ! SY: Soni's dimensionality-apt contraction
        sceua_struc(t)%parsets( sceua_ctl%NComplexes * &
                                sceua_ctl%NPointsInComplex + &
                                (iComplex-1)*4+3, j ) = &
                                      sceua_struc(t)%parset(j)
      end do

!  CHOOSE THE 4th NEW POINT AT RANDOM WITHIN FEASIBLE REGION ACCORDING TO
!  A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
!  AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
      call getpnt(sceua_ctl%nparam,2,sceua_ctl%iseed1, &
                  sceua_struc(t)%parset,sceua_ctl%parmin, &
                  sceua_ctl%parmax, sceua_struc(t)%parstd, sb) ! SY
      do j = 1, m
        sceua_struc(t)%parsets( sceua_ctl%NComplexes * &
                                sceua_ctl%NPointsInComplex + &
                                (iComplex-1)*4+4, j ) = &
                                      sceua_struc(t)%parset(j)
      end do
 
    end subroutine SCEUA_SubcomplexCreateChildren


!BOP
! !ROUTINE: SCEUA_SubcomplexNewGeneration
! \label{SCEUA_SubcomplexNewGeneration}
!
! !INTERFACE:
    subroutine SCEUA_SubcomplexNewGeneration(iComplex, t)
!
! !USES:
      use LIS_coreMod,  only : LIS_rc
      use LIS_optUEMod, only : LIS_decisionSpace
      use LIS_logMod,          only : LIS_logunit ! SY: del. later
!
      implicit none
!
! !INPUT PARAMETERS: 
      integer, intent(IN)           :: iComplex
      integer, intent(IN)           :: t
!
! !DESCRIPTION:
!
!   This method performs the computations to update the subcomplex
!   worst parent with the relevant best child out of 4 created earlier 
!   using the subcomplex
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP
      integer      :: n
      integer      :: j,g,k,i
      real         :: rand
      integer      :: irand

      integer      :: ibound
      REAL*8       :: sw(sceua_ctl%nparam) ! SY: Worst sub-complex parent
      REAL*8       :: swObjFunc ! SY: Obj. func. value of sw 
      REAL*8       :: ObjFunc 
      REAL*8       :: swTemp(sceua_ctl%nparam) ! SY: Swapping buffer point
      REAL*8       :: swObjFuncTemp ! SY: Swapping buffer obj. func.

!  EXTRACT WORST PARENT AND 1st CHILD (FROM REFLECTION) OF SUBCOMPLEX, 
!  WITH THEIR CORRESPONDING OBJ. FUNC. VALUES
      do j = 1, sceua_ctl%nparam
        sw(j) = sceua_struc(t)%subcomplexesParsets(iComplex, &
                                    sceua_ctl%NPointsInSubcomplex,j) ! SY
        sceua_struc(t)%parset(j) = &
                         sceua_struc(t)%parsets(sceua_ctl%NComplexes * &
                                                sceua_ctl%NPointsInComplex + &
                                                (iComplex-1)*4+1, j ) 
      end do
      swObjFunc = sceua_struc(t)%subcomplexesObjFuncs(iComplex, &
                                            sceua_ctl%NPointsInSubcomplex )
      ObjFunc = sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                        sceua_ctl%NPointsInComplex + &
                                        (iComplex-1)*4+1 ) 

!  IF 1st CHILD DOES SATISFY CONSTRAINTS, AND ITS Obj. Func. VALUE IS LESS
!  THAN THAT OF THE WORST SUBCOMPLEX PARENT, THEN ACCEPT AND SWAP WITH
!  THAT WORST SUBCOMPLEX PARENT 
      if ((sceua_struc(t)%ibound(iComplex) .lt. 1) .AND. &
          (ObjFunc .LE. swObjFunc) .AND. (ObjFunc .LE. sceua_ctl%maxObjFunc/2.0)) then ! SY
        do j = 1, sceua_ctl%nparam
          swTemp(j) = sw(j)
! SY: Assigned to sw-equivalent sceua_struc(t)%subcomplexesParsets
! now below instead of to sw or sceua_struc(t)%parsets parent, will be 
! assigned through that to sceua_struc(t)%parsets parent later by 
! intermediately passing to sceua_struc(t)%complexesParsets 
          sceua_struc(t)%subcomplexesParsets(iComplex, & 
                              sceua_ctl%NPointsInSubcomplex,j ) = &
                   sceua_struc(t)%parsets(sceua_ctl%NComplexes * & 
                                          sceua_ctl%NPointsInComplex + &
                                          (iComplex-1)*4+1, j ) 
! SY: Assigned directly to sceua_struc(t)%parsets child instead of to 
! the equivalent sceua_struc(t)%parset
          sceua_struc(t)%parsets(sceua_ctl%NComplexes * &
                                 sceua_ctl%NPointsInComplex + &
                                 (iComplex-1)*4+1, j ) = swTemp(j)
        end do
        swObjFuncTemp = swObjFunc
! SY: Assigned to swObjFunc-equivalent sceua_struc(t)%subcomplexesObjFuncs
! now below instead of to swObjFuncTemp or sceua_struc(t)%objfuncs parent, 
! will be assigned through that to sceua_struc(t)%objfuncs parent later by
! intermediately passing to sceua_struc(t)%complexesObjFuncs 
        sceua_struc(t)%subcomplexesObjFuncs(iComplex, &
                                            sceua_ctl%NPointsInSubcomplex) = &
                   sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                           sceua_ctl%NPointsInComplex + &
                                           (iComplex-1)*4+1 )
! SY: Assigned directly to sceua_struc(t)%objfuncs child instead of to 
! the equivalent ObjFunc
        sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                sceua_ctl%NPointsInComplex + &
                                (iComplex-1)*4+1 ) = swObjFuncTemp


!        write(LIS_logunit,*) 'In SCEUA_SubcomplexNewGeneration, point 1'
        RETURN
      end if

!  EXTRACT 2nd CHILD (FROM MUTATION) OF SUBCOMPLEX, 
!  WITH CORRESPONDING OBJ. FUNC. VALUE
      do j = 1, sceua_ctl%nparam
        sceua_struc(t)%parset(j) = &
                         sceua_struc(t)%parsets(sceua_ctl%NComplexes * &
                                                sceua_ctl%NPointsInComplex + &
                                                (iComplex-1)*4+2, j ) 
      end do
      ObjFunc = sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                        sceua_ctl%NPointsInComplex + &
                                        (iComplex-1)*4+2 ) 

!  IF 2nd CHILD Obj. Func. VALUE IS LESS
!  THAN THAT OF THE WORST SUBCOMPLEX PARENT, THEN ACCEPT AND SWAP WITH
!  THAT WORST SUBCOMPLEX PARENT 
      if ((ObjFunc .LE. swObjFunc) .AND. (ObjFunc .LE. sceua_ctl%maxObjFunc/2.0)) then ! SY
        do j = 1, sceua_ctl%nparam
          swTemp(j) = sw(j)
! SY: Assigned to sw-equivalent sceua_struc(t)%subcomplexesParsets
! now below instead of to sw or sceua_struc(t)%parsets parent, will be 
! assigned through that to sceua_struc(t)%parsets parent later by 
! intermediately passing to sceua_struc(t)%complexesParsets 
          sceua_struc(t)%subcomplexesParsets(iComplex, & 
                              sceua_ctl%NPointsInSubcomplex,j ) = &
                   sceua_struc(t)%parsets(sceua_ctl%NComplexes * & 
                                          sceua_ctl%NPointsInComplex + &
                                          (iComplex-1)*4+2, j ) 
! SY: Assigned directly to sceua_struc(t)%parsets child instead of to 
! the equivalent sceua_struc(t)%parset
          sceua_struc(t)%parsets(sceua_ctl%NComplexes * &
                                 sceua_ctl%NPointsInComplex + &
                                 (iComplex-1)*4+2, j ) = swTemp(j)
        end do
        swObjFuncTemp = swObjFunc
! SY: Assigned to swObjFunc-equivalent sceua_struc(t)%subcomplexesObjFuncs
! now below instead of to swObjFuncTemp or sceua_struc(t)%objfuncs parent,
! will be assigned through that to sceua_struc(t)%objfuncs parent later by
! intermediately passing to sceua_struc(t)%complexesObjFuncs 
        sceua_struc(t)%subcomplexesObjFuncs(iComplex, &
                                            sceua_ctl%NPointsInSubcomplex) = &
                   sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                           sceua_ctl%NPointsInComplex + &
                                           (iComplex-1)*4+2 )
! SY: Assigned directly to sceua_struc(t)%objfuncs child instead of to 
! the equivalent ObjFunc
        sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                sceua_ctl%NPointsInComplex + &
                                (iComplex-1)*4+2 ) = swObjFuncTemp

!        write(LIS_logunit,*) 'In SCEUA_SubcomplexNewGeneration, point 2'
        RETURN
      end if

!  EXTRACT 3rd CHILD (FROM CONTRACTION) OF SUBCOMPLEX, 
!  WITH CORRESPONDING OBJ. FUNC. VALUE
      do j = 1, sceua_ctl%nparam
        sceua_struc(t)%parset(j) = &
                         sceua_struc(t)%parsets(sceua_ctl%NComplexes * &
                                                sceua_ctl%NPointsInComplex + &
                                                (iComplex-1)*4+3, j ) 
      end do
      ObjFunc = sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                        sceua_ctl%NPointsInComplex + &
                                        (iComplex-1)*4+3 ) 

!  IF 3rd CHILD Obj. Func. VALUE IS LESS
!  THAN THAT OF THE WORST SUBCOMPLEX PARENT, THEN ACCEPT AND SWAP WITH
!  THAT WORST SUBCOMPLEX PARENT 
      if ((ObjFunc .LE. swObjFunc) .AND. (ObjFunc .LE. sceua_ctl%maxObjFunc/2.0)) then ! SY
        do j = 1, sceua_ctl%nparam
          swTemp(j) = sw(j)
! SY: Assigned to sw-equivalent sceua_struc(t)%subcomplexesParsets
! now below instead of to sw or sceua_struc(t)%parsets parent, will be 
! assigned through that to sceua_struc(t)%parsets parent later by 
! intermediately passing to sceua_struc(t)%complexesParsets 
          sceua_struc(t)%subcomplexesParsets(iComplex, & 
                              sceua_ctl%NPointsInSubcomplex,j ) = &
                   sceua_struc(t)%parsets(sceua_ctl%NComplexes * & 
                                          sceua_ctl%NPointsInComplex + &
                                          (iComplex-1)*4+3, j ) 
! SY: Assigned directly to sceua_struc(t)%parsets child instead of to 
! the equivalent sceua_struc(t)%parset
          sceua_struc(t)%parsets(sceua_ctl%NComplexes * &
                                 sceua_ctl%NPointsInComplex + &
                                 (iComplex-1)*4+3, j ) = swTemp(j)
        end do
        swObjFuncTemp = swObjFunc
! SY: Assigned to swObjFunc-equivalent sceua_struc(t)%subcomplexesObjFuncs
! now below instead of to swObjFuncTemp or sceua_struc(t)%objfuncs parent,
! will be assigned through that to sceua_struc(t)%objfuncs parent later by
! intermediately passing to sceua_struc(t)%complexesObjFuncs 
        sceua_struc(t)%subcomplexesObjFuncs(iComplex, &
                                            sceua_ctl%NPointsInSubcomplex) = &
                   sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                           sceua_ctl%NPointsInComplex + &
                                           (iComplex-1)*4+3 )
! SY: Assigned directly to sceua_struc(t)%objfuncs child instead of to 
! the equivalent ObjFunc
        sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                sceua_ctl%NPointsInComplex + &
                                (iComplex-1)*4+3 ) = swObjFuncTemp

!        write(LIS_logunit,*) 'In SCEUA_SubcomplexNewGeneration, point 3'
        RETURN
      end if

!  EXTRACT 4th CHILD (FROM MUTATION) OF SUBCOMPLEX, 
!  WITH CORRESPONDING OBJ. FUNC. VALUE
      do j = 1, sceua_ctl%nparam
        sceua_struc(t)%parset(j) = &
                         sceua_struc(t)%parsets(sceua_ctl%NComplexes * &
                                                sceua_ctl%NPointsInComplex + &
                                                (iComplex-1)*4+4, j ) 
      end do
      ObjFunc = sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                        sceua_ctl%NPointsInComplex + &
                                        (iComplex-1)*4+4 ) 

!  SINCE ABOVE THREE CHILDREN ARE REJECTED, THIS 4th CHILD IS
!  COMPULSORILY ACCEPTED AND SWAPPED WITH THE WORST SUBCOMPLEX PARENT 
! SY: This now changed to make sure that an infeasible space child is
!     not picked up, especially after discarding a feasible space parent!!
!     So, at this stage, all 4 children can be potentially rejected 
      if (ObjFunc .LE. sceua_ctl%maxObjFunc/2.0) then ! SY
        do j = 1, sceua_ctl%nparam
          swTemp(j) = sw(j)
! SY: Assigned to sw-equivalent sceua_struc(t)%subcomplexesParsets
! now below instead of to sw or sceua_struc(t)%parsets parent, will be 
! assigned through that to sceua_struc(t)%parsets parent later by 
! intermediately passing to sceua_struc(t)%complexesParsets 
          sceua_struc(t)%subcomplexesParsets(iComplex, & 
                              sceua_ctl%NPointsInSubcomplex,j ) = &
                   sceua_struc(t)%parsets(sceua_ctl%NComplexes * & 
                                          sceua_ctl%NPointsInComplex + &
                                          (iComplex-1)*4+4, j ) 
! SY: Assigned directly to sceua_struc(t)%parsets child instead of to 
! the equivalent sceua_struc(t)%parset
          sceua_struc(t)%parsets(sceua_ctl%NComplexes * &
                                 sceua_ctl%NPointsInComplex + &
                                 (iComplex-1)*4+4, j ) = swTemp(j)
        end do
        swObjFuncTemp = swObjFunc
! SY: Assigned to swObjFunc-equivalent sceua_struc(t)%subcomplexesObjFuncs
! now below instead of to swObjFuncTemp or sceua_struc(t)%objfuncs parent,
! will be assigned through that to sceua_struc(t)%objfuncs parent later by
! intermediately passing to sceua_struc(t)%complexesObjFuncs 
        sceua_struc(t)%subcomplexesObjFuncs(iComplex, &
                                            sceua_ctl%NPointsInSubcomplex) = &
                   sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                           sceua_ctl%NPointsInComplex + &
                                           (iComplex-1)*4+4 )
! SY: Assigned directly to sceua_struc(t)%objfuncs child instead of to 
! the equivalent ObjFunc
        sceua_struc(t)%objfuncs(sceua_ctl%NComplexes * &
                                sceua_ctl%NPointsInComplex + &
                                (iComplex-1)*4+4 ) = swObjFuncTemp

!        write(LIS_logunit,*) 'In SCEUA_SubcomplexNewGeneration, point 4'
        RETURN
      end if

! SY: NO translating the updated decision space to generic objects done here
!     because they are done later anyway in SCEUAOpt_run before 
!     SCEUA_evaluateObjFunc where it is really required. 
!      call setoptuetypedecspace(LIS_rc%optuetype)

    end subroutine SCEUA_SubcomplexNewGeneration


!BOP
! !ROUTINE: SCEUA_printPopulation
! \label{SCEUA_printPopulation}
! 
! !INTERFACE:
  subroutine SCEUA_printPopulation()
!
! !DESCRIPTION: 
!  This routine prints out the population with the associated
!  optimization iteration information. 
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP

    integer                 :: n

    call writeSCEUArestart()
    call writeSCEUAoutput()

  end subroutine SCEUA_printPopulation


!BOP
! !ROUTINE: SCEUA_getdecSpaceValues
! \label{SCEUA_getdecSpaceValues}
! 
! !INTERFACE:
    subroutine SCEUA_getdecSpaceValues(n, decvals)
! !USES:
      use LIS_coreMod,    only : LIS_rc
!      
      implicit none
!
! !INPUT PARAMETERS: 
      integer, intent(IN)            :: n
! !INPUT/OUTPUT PARAMETERS: 
      real               :: decvals(sceua_ctl%nparam, LIS_rc%ntiles(n))
! !DESCRIPTION: 
!  This routine copies all algorithm parameter sets into single array
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP

      integer            :: i, t, m 

      do i=1, sceua_ctl%nparam
         do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            do m=1,LIS_rc%nensem(n)
               decvals(i,(t-1)*LIS_rc%nensem(n)+m) = sceua_struc(t)%parsets(m,i) 
            enddo
         enddo
      enddo

    end subroutine SCEUA_getdecSpaceValues

!BOP
! !ROUTINE: SCEUA_setdecSpaceValues
! \label{SCEUA_setdecSpaceValues}
! 
! !INTERFACE:
    subroutine SCEUA_setdecSpaceValues(n, decvals)
! !USES:
      use LIS_coreMod,    only : LIS_rc
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(IN)            :: n
! !INPUT/OUTPUT PARAMETERS:
      real               :: decvals(sceua_ctl%nparam, LIS_rc%ntiles(n))
! !DESCRIPTION:
!  This routine acquires updated algorithm parameter sets from the single array
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
! 
!EOP 

      integer            :: i, t, m

      do i=1, sceua_ctl%nparam
         do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            do m=1,LIS_rc%nensem(n)
               sceua_struc(t)%parsets(m,i) = decvals(i,(t-1)*LIS_rc%nensem(n)+m) 
            enddo
         enddo
      enddo

    end subroutine SCEUA_setdecSpaceValues

!BOP
! !ROUTINE: SCEUA_getNparam
! \label{SCEUA_getNparam}
! 
! !INTERFACE:
    subroutine SCEUA_getNparam(nparam)
!
! !OUTPUT PARAMETERS:
      integer, intent(OUT)   :: nparam
!
! !DESCRIPTION:
! Get number of parameters in an optimized parameter set
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
! 
!EOP 

      nparam = sceua_ctl%nparam

    end subroutine SCEUA_getNparam


!BOP
! !ROUTINE: getpnt
! \label{getpnt}
!
! !INTERFACE:
    subroutine getpnt(nopt,idist,iseed,x,bl,bu,std,xi)
! !USES: 
      use LIS_logMod,          only : LIS_logunit ! SY
!
      implicit real*8 (a-h,o-z)
!
! !INPUT PARAMETERS: 
!      dimension x(16),bl(16),bu(16),std(16),xi(16) ! SY
      INTEGER, intent(IN)     :: nopt,idist
      integer                 :: iseed
      REAL*8, intent(IN)      :: bl(nopt),bu(nopt),std(nopt),xi(nopt) !SY
! !OUTPUT PARAMETERS: 
      REAL*8, intent(OUT)     :: x(nopt) ! SY
!
! !DESCRIPTION: 
!     This subroutine generates a new point within feasible region
!
!     x(.) = new point
!     xi(.) = focal point
!     bl(.) = lower bound
!     bu(.) = upper bound
!     std(.) = standard deviation of probability distribution
!     idist = probability flag
!           = 1 - uniform distribution
!           = 2 - Gaussian distribution
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 
      INTEGER     :: j, ibound  
      REAL*8      :: rand ! SY  
 
      DO ! SY
        do j=1, nopt ! SY
!    1 do j=1, nopt ! SY

          DO ! SY

            if (idist .eq. 1) rand = ran1(iseed) ! SY
!    2   if (idist .eq. 1) rand = ran1(iseed) ! SY
!            if (idist .eq. 2) rand = gasdev(iseed) ! SY: Replace back later with below!
            if (idist .eq. 2) then
              rand = gasdev(iseed)
!              write(LIS_logunit,*) 'rand is ', rand
            end if
            x(j) = xi(j) + std(j) * rand * (bu(j) - bl(j))

!     Check explicit constraints

            call chkcst(1,x(j),bl(j),bu(j),ibound)
!            if (ibound .ge. 1) go to 2 ! SY
            if (ibound .LT. 1) EXIT ! SY
          END DO ! SY
        end do

!     Check implicit constraints
      
        call chkcst(nopt,x,bl,bu,ibound)
!        if (ibound .ge. 1) go to 1 ! SY
        if (ibound .LT. 1) EXIT ! SY
      
      END DO ! SY

      return

    end subroutine getpnt

!BOP
! !ROUTINE: ran1
! \label{ran1}
!
! !INTERFACE:
    real*8 function ran1(idum)
! !USES: 
      use LIS_logMod,          only : LIS_logunit ! SY
!
      implicit real*8 (a-h,o-z)
!
! !INPUT/OUTPUT PARAMETERS: 
      INTEGER, intent(INOUT)     :: idum 
!
! !DESCRIPTION: 
!  THIS FUNCTION IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
!  This function generates random number from a seed
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 

!      INTEGER            :: idum,iff,ix1,ix2,ix3,j ! SY: See below
      INTEGER            :: ix1,ix2,ix3,j ! SY: iff moved to sceua_ctl
!      dimension r(97)
      REAL*8             :: r(97)
      INTEGER, parameter :: m1 = 259200, ia1 = 7141, ic1 = 54773 ! SY: Check later about machine-limiting this!!
      REAL*8, parameter  :: rm1 = 3.8580247e-6 ! SY: Check later about machine-limiting this!!
      INTEGER, parameter :: m2 = 134456, ia2 = 8121, ic2 = 28411 ! SY: Check later about machine-limiting this!!
      REAL*8, parameter  ::  rm2 = 7.4373773e-6 ! SY: Check later about machine-limiting this!!
      INTEGER, parameter :: m3 = 243000, ia3 = 4561, ic3 = 51349
      save
!      data iff / 0 / ! SY: Now taken care of in SCEUAOpt_init

      if ((idum .lt. 0) .or. (sceua_ctl%ran1_iff .eq. 0)) then
         sceua_ctl%ran1_iff = 1
         ix1 = mod(ic1 - idum,m1)
         ix1 = mod((ia1 * ix1) + ic1,m1)
         ix2 = mod(ix1,m2)
         ix1 = mod((ia1 * ix1) + ic1,m1)
         ix3 = mod(ix1,m3)
!         do 11 j = 1, 97 ! SY
         do j = 1, 97 ! SY
            ix1 = mod((ia1 * ix1) + ic1,m1)
            ix2 = mod((ia2 * ix2) + ic2,m2)
            r(j) = (dble(ix1) + (dble(ix2) * rm2)) * rm1
         end do ! SY
!         11 continue ! SY
         idum = 1
      end if
      ix1 = mod((ia1 * ix1) + ic1,m1)
      ix2 = mod((ia2 * ix2) + ic2,m2)
      ix3 = mod((ia3 * ix3) + ic3,m3)
      j = 1 + ((97 * ix3) / m3)
!      if ((j .gt. 97) .or. (j .lt. 1)) pause ! SY
      ! SY: Begin
      if ((j .gt. 97) .or. (j .lt. 1)) then 
        write(LIS_logunit,*) 'WARNING: Supposed to PAUSE here in ran1, '
        write(LIS_logunit,*) 'in module ShuffledComplexEvolution: '
        write(LIS_logunit,*) '(j .gt. 97) .or. (j .lt. 1) '
      end if
      ! SY: End
      ran1 = r(j)
      r(j) = (dble(ix1) + (dble(ix2) * rm2)) * rm1 ! SY: This seems not used: Fishy!!

!  END OF FUNCTION RAN1
      return

    end function ran1


!BOP
! !ROUTINE: gasdev
! \label{gasdev}
!
! !INTERFACE:
    real*8 function gasdev(idum)
! !USES: 
      use LIS_logMod,          only : LIS_logunit ! SY
!
      implicit real*8 (a-h,o-z)
!
! !INPUT/OUTPUT PARAMETERS: 
      INTEGER, intent(INOUT)     :: idum 
!
! !DESCRIPTION: 
!  THIS FUNCTION IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
!  This function generates random number from a seed
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 
!      common /gasblk/ iset ! SY
      INTEGER :: iset
      data iset / 0 /
      REAL*8  :: v1,v2,r,fac,gset
      data v1,v2,r,fac,gset / 0.0,0.0,0.0,0.0,0.0 /

      if (iset .eq. 0) then

        DO ! SY

!    1 v1 = (2. * ran1(idum)) - 1. ! SY
          v1 = (2. * ran1(idum)) - 1. ! SY
          v2 = (2. * ran1(idum)) - 1.
          r = (v1 ** 2) + (v2 ** 2)
!          if (r .ge. 1.) goto 1 ! SY
          if (r .LT. 1.) EXIT ! SY

        END DO ! SY

        fac = sqrt(- ((2. * log(r)) / r))
        gset = v1 * fac
        gasdev = v2 * fac
        iset = 1
      else
        gasdev = gset
        iset = 0
      end if

! SY: Note that v2,fac both going near machine-precision zero means gasdev becomes Nan!! Not printing out gasdev value seems ok, investigate why later.

!  END OF FUNCTION GASDEV
      return

    end function gasdev


!BOP
! !ROUTINE: chkcst
! \label{chkcst}
!
! !INTERFACE:
    subroutine chkcst(nopt,x,bl,bu,ibound)
!
      implicit real*8 (a-h,o-z)
! !INPUT PARAMETERS: 
      INTEGER, intent(IN)     :: nopt 
      REAL*8, intent(IN)      :: x(nopt),bl(nopt),bu(nopt)
! !INPUT/OUTPUT PARAMETERS: 
      INTEGER, intent(INOUT)     :: ibound 
!
! !DESCRIPTION: 
!     This subroutine check if the trial point satisfies all
!     constraints.
!
!     ibound - violation indicator
!            = -1 initial value ! NOTE: This ALSO means no violation \& nopt gt 1 !!
!            = 0  no violation
!            = 1  violation
!     nopt = number of optimizing variables
!     ii = the ii'th variable of the arrays x, bl, and bu
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 
!      common /gasblk/ iset ! SY
      INTEGER :: ii

      ibound = -1

!     Check if explicit constraints are violated

      do ii=1, nopt
!        if (x(ii) .lt. bl(ii) .or. x(ii) .gt. bu(ii)) go to 10 ! SY
!       SY: Begin
        if (x(ii) .lt. bl(ii) .or. x(ii) .gt. bu(ii)) then
          ibound = 1
          return
        end if
!       SY: End
      end do
!      if (nopt .eq. 1) go to 9 ! SY
!     SY: Begin
      if (nopt .eq. 1) then
        ibound = 0
        return
      end if
!     SY: End

!     Check if implicit constraints are violated
!     (no implicit constraints for this function)

!     No constraints are violated

!    9 ibound = 0 ! SY
!      return ! SY

!     At least one of the constraints are violated

!   10 ibound = 1 ! SY
!      return ! SY

      return ! SY

    end subroutine chkcst


!BOP
! !ROUTINE: sort
! \label{sort}
!
! !INTERFACE:
      subroutine sort(n,m,rb,ra)
!
        implicit real*8 (a-h,o-z)
! !INPUT PARAMETERS: 
      INTEGER, intent(IN)     :: n,m 
! !INPUT/OUTPUT PARAMETERS: 
!        dimension ra(2000),rb(2000,16),wk(2000,16),iwk(2000) ! SY
      REAL*8, intent(INOUT)      :: ra(n),rb(n,m) ! SY
!
! !DESCRIPTION: 
!  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
!  BY W.H. PRESS ET AL., pp. 233-234
!
!  LIST OF VARIABLES
!     ra(.) = array to be sorted
!     rb(.,.) = arrays ordered corresponding to rearrangement of ra(.)
!     wk(.,.), iwk(.) = local varibles
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 

        INTEGER :: iwk(n),i,j ! SY
        REAL*8 :: wk(n,m) ! SY

        call indexx(n, ra, iwk)
!      do 11 i = 1, n ! SY
        do i = 1, n ! SY
          wk(i,1) = ra(i)
        enddo ! SY
!   11 continue ! SY
!      do 12 i = 1, n ! SY
        do i = 1, n ! SY
          ra(i) = wk(iwk(i),1)
        enddo ! SY
!   12 continue ! SY
!      do 14 j = 1, m ! SY 
!      do 13 i = 1, n ! SY
        do j = 1, m ! SY
          do i = 1, n ! SY
            wk(i,j) = rb(i,j)
          enddo ! SY
        enddo ! SY
!   13 continue ! SY
!   14 continue ! SY
!      do 16 j = 1, m ! SY
!      do 15 i = 1, n ! SY
        do j = 1, m ! SY
          do i = 1, n ! SY
            rb(i,j) = wk(iwk(i),j)
          enddo ! SY
        enddo ! SY
!   15 continue ! SY
!   16 continue  ! SY

!  END OF SUBROUTINE SORT
        return
      end subroutine sort


!BOP
! !ROUTINE: indexx
! \label{indexx}
!
! !INTERFACE:
      subroutine indexx(n, arrin, indx)
!
        implicit real*8 (a-h,o-z)
!
! !INPUT PARAMETERS: 
        INTEGER, intent(IN)     :: n
        REAL*8, intent(IN)      :: arrin(n)
! !OUTPUT PARAMETERS:
        INTEGER, intent(OUT)     :: indx(n)
!
! !DESCRIPTION: 
!   THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
!   Calculates the sorting index order for the array 
!   to be sorted in the 'sort' routine 
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 
        INTEGER :: j, l, ir, indxt, i
        REAL*8 :: q

!      do 11 j = 1, n ! SY
        do j = 1, n ! SY
          indx(j) = j
        enddo ! SY
!   11 continue ! SY
        l = (n / 2) + 1
        ir = n
!   10 continue ! SY
        DO ! SY 
          if (l .gt. 1) then
            l = l - 1
            indxt = indx(l)
            q = arrin(indxt)
          else
            indxt = indx(ir)
            q = arrin(indxt)
            indx(ir) = indx(1)
            ir = ir - 1
            if (ir .eq. 1) then
              indx(1) = indxt
!      return ! SY
              EXIT ! SY
            end if
          end if
          i = l
          j = l + l
          DO ! SY 
            if (j .le. ir) then ! SY
!   20 if (j .le. ir) then
              if (j .lt. ir) then
                if (arrin(indx(j)) .lt. arrin(indx(j + 1))) j = j + 1
              end if
              if (q .lt. arrin(indx(j))) then
                indx(i) = indx(j)
                i = j
                j = j + j
              else
                j = ir + 1
              end if
!      goto 20 ! SY
            else ! SY
              EXIT ! SY
            end if
          END DO ! SY 
          indx(i) = indxt
!      goto 10 ! SY
        END DO ! SY 

!  END OF SUBROUTINE INDEXX
        return ! SY

      end subroutine indexx


!BOP
! !ROUTINE: parstt
! \label{parstt}
!
! !INTERFACE:
      subroutine parstt(npt,nopt,x,xnstd,bound,gnrng,ipcnvg)
!
        implicit real*8 (a-h,o-z)
!
! !INPUT PARAMETERS: 
        INTEGER, intent(IN)     :: npt, nopt
!        dimension x(2000,16),xmax(16),xmin(16) ! SY
!        dimension xmean(16),xnstd(16),bound(16) ! SY
        REAL*8, intent(IN)     :: x(npt,nopt),bound(nopt) ! SY
! !OUTPUT PARAMETERS: 
        REAL*8, intent(OUT)     :: xnstd(nopt),gnrng ! SY
        INTEGER, intent(OUT)     :: ipcnvg
!
! !DESCRIPTION: 
!   SUBROUTINE CHECKING FOR PARAMETER CONVERGENCE
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 
        INTEGER :: k, i 
        REAL*8 :: xmax(nopt),xmin(nopt) ! SY
        REAL*8 :: xmean(nopt) ! SY
        REAL*8 :: gsum,xsum1,xsum2
        REAL*8, parameter :: delta = 1.0d-20,peps=1.0d-3 ! SY: Check later about machine-limiting these!!

!  COMPUTE MAXIMUM, MINIMUM AND STANDARD DEVIATION OF PARAMETER VALUES
        gsum = 0.d0
        do k = 1, nopt 
          xmax(k) = -1.0d+20 ! SY: Check later about machine-limiting this!!
          xmin(k) = 1.0d+20 ! SY: Check later about machine-limiting this!!
          xsum1 = 0.d0
          xsum2 = 0.d0
          do i = 1, npt
            xmax(k) = dmax1(x(i,k), xmax(k))
            xmin(k) = dmin1(x(i,k), xmin(k))
            xsum1 = xsum1 + x(i,k)
            xsum2 = xsum2 + x(i,k)*x(i,k)
          end do
          xmean(k) = xsum1 / dble(npt)
          xnstd(k) = (xsum2 / dble(npt) - xmean(k)*xmean(k))
          if (xnstd(k) .le. delta) xnstd(k) = delta 
          xnstd(k) = dsqrt(xnstd(k))
          xnstd(k) = xnstd(k) / bound(k)
          gsum = gsum + dlog( delta + (xmax(k)-xmin(k))/bound(k) )
        end do  
        gnrng = dexp(gsum/dble(nopt))

!  CHECK IF NORMALIZED STANDARD DEVIATION OF PARAMETER IS <= eps
        ipcnvg = 0
        if (gnrng .le. peps) then
          ipcnvg = 1
        end if  

!  END OF SUBROUTINE PARSTT
        return
  
      end subroutine parstt    


!BOP
! !ROUTINE: sort1
! \label{sort1}
!
! !INTERFACE:
      subroutine sort1(n,ra)
!     
        implicit real*8 (a-h,o-z)
!
! !INPUT PARAMETERS: 
        INTEGER, intent(IN)     :: n
! !INPUT/OUTPUT PARAMETERS: 
!        dimension ra(n) ! SY
!        integer ra, rra ! SY
        INTEGER, intent(INOUT)     :: ra(n) ! SY
!
! !DESCRIPTION: 
!  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
!  BY W.H. PRESS ET AL., pp. 231
!  
!  LIST OF VARIABLES
!     ra(.) = integer array to be sorted
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 
     
        INTEGER :: rra,l,ir,j,i  ! SY
  
        l = (n / 2) + 1
        ir = n
        DO ! SY
!   10 continue
          if (l .gt. 1) then 
            l = l - 1
            rra = ra(l)
          else
            rra = ra(ir)
            ra(ir) = ra(1)
            ir = ir - 1
            if (ir .eq. 1) then
              ra(1) = rra
              RETURN
            end if
          end if
          i = l
          j = l + l
          DO ! SY
!   20 if (j .le. ir) then ! SY
            if (j .le. ir) then
              if (j .lt. ir) then
                if (ra(j) .lt. ra(j + 1)) j = j + 1
              end if
              if (rra .lt. ra(j)) then
                ra(i) = ra(j)
                i = j
                j = j + j
              else
                j = ir + 1
              end if
 !             goto 20 ! SY
            else
              EXIT ! SY
            end if
          END DO ! SY
          ra(i) = rra
!      goto 10 ! SY
        END DO ! SY

!  END OF SUBROUTINE SORT1
      end subroutine sort1


!BOP
!
! !ROUTINE: readSCEUArestart
! \label{readSCEUArestart}
!
! !INTERFACE:
  subroutine readSCEUArestart
! !USES:
    use LIS_optUEMod,        only : LIS_decisionSpace 
    use LIS_coreMod,         only : LIS_rc
    use LIS_logMod,          only : LIS_logunit, LIS_verify
    use LIS_historyMod,      only : LIS_readvar_restart
!   
! !DESCRIPTION:
! 
!   This routine reads the checkpoint data for a SCEUA restart
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 
    integer             :: n
    integer             :: i, t
    integer             :: status
    character*100       :: vnames(sceua_ctl%nparam)
    type(ESMF_Field)    :: varField(sceua_ctl%nparam)
    real, pointer       :: vardata(:)
    real,    allocatable    :: decvalues(:,:) 
   
    allocate( decvalues(sceua_ctl%nparam, LIS_rc%ntiles(n)) ) 
 
    write(LIS_logunit,*) 'Reading the SCEUA restart file ..'
    n = 1
    call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
         rc=status)
    call LIS_verify(status)

    open(40,file=sceua_ctl%rfile,form='unformatted')

    read(40) sceua_ctl%SumReplacedFuncEvals
    write(LIS_logunit,*) 'Num. of REPLACED member function (LSM) evaluations finished: ',sceua_ctl%SumReplacedFuncEvals

    read(40) sceua_ctl%SumNewMemberFuncEvals
    write(LIS_logunit,*) 'Num. of NEW member function (LSM) evaluations finished: ',sceua_ctl%SumNewMemberFuncEvals

    read(40) sceua_ctl%SumAllFuncEvals
    write(LIS_logunit,*) 'Num. of OVERALL member function (LSM) evaluations finished: ',sceua_ctl%SumAllFuncEvals

    do i=1,sceua_ctl%nparam
       call ESMF_StateGet(LIS_decisionSpace, trim(vnames(i)), &
            varField(i), rc=status)
       call LIS_verify(status)

       call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata,rc=status)
       call LIS_verify(status)

       call LIS_readvar_restart(40,n,vardata)

!      SY: following mimics updateDecSpaceObject
!       call updateDecSpaceObject(sceua_ctl%nparam, LIS_rc%ntiles(n), decvalues, .false.)
       do t=1,LIS_rc%ntiles(n)
         decvalues(i,t) = vardata(t)
       enddo 

    enddo

    close(40)
    write(LIS_logunit,*) 'Finished reading the SCEUA restart file ..'

!   set them to the algorithm data structures
    call setOptUEAlgDecSpace(LIS_rc%optUEAlg, n, decvalues) 

    deallocate(decvalues) 

  end subroutine readSCEUArestart


!BOP
!
! !ROUTINE: writeSCEUArestart
! \label{writeSCEUArestart}
!
! !INTERFACE:
  subroutine writeSCEUArestart
! !USES:
    use LIS_optUEMod, only : LIS_decisionSpace 
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_fileIOMod, only : LIS_create_output_directory
    use LIS_logMod,    only : LIS_logunit, LIS_verify
    use LIS_historyMod, only : LIS_writevar_restart
!
! !DESCRIPTION:
!
! This routine writes the checkpoint data for a SCEUA restart
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 
    integer             :: n
    integer             :: i
    integer             :: status
    character(len=LIS_CONST_PATH_LEN) :: filen
!    character (len=3)   :: fTrial 
    character (len=6)   :: fTrial 
    character*100       :: vnames(sceua_ctl%nparam)
    type(ESMF_Field)    :: varField(sceua_ctl%nparam)
    real, pointer       :: vardata(:)

    n = 1

    call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
         rc=status)
    call LIS_verify(status)

    if(LIS_masterproc) then
       call LIS_create_output_directory('SCEUA')
!       write(unit=fTrial, fmt='(i3.3)') sceua_ctl%genNo 
       write(unit=fTrial, fmt='(i6.6)') sceua_ctl%SumNewMemberFuncEvals 
       filen = trim(LIS_rc%odir)//'/SCEUA/SCEUA.'&
            //trim(fTrial)//'.SCEUArst'
       open(40,file=filen,status='unknown',form='unformatted')
       write(40) sceua_ctl%SumReplacedFuncEvals
       write(40) sceua_ctl%SumNewMemberFuncEvals
       write(40) sceua_ctl%SumAllFuncEvals
    endif

    do i=1,sceua_ctl%nparam
       call ESMF_StateGet(LIS_decisionSpace, trim(vnames(i)), &
            varField(i), rc=status)
       call LIS_verify(status)

       call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata, rc=status)
       call LIS_verify(status)

       call LIS_writevar_restart(40,n,vardata)
    enddo

    if(LIS_masterproc) then
       close(40)
       write(LIS_logunit,*) 'SCEUA checkpoint file written ',trim(filen)
    endif
  end subroutine writeSCEUArestart


!BOP
!
! !ROUTINE: writeSCEUAoutput
! \label{writeSCEUAoutput}
!
! !INTERFACE:
  subroutine writeSCEUAoutput
! !USES:
    use LIS_optUEMod,        only : LIS_decisionSpace 
    use LIS_coreMod,         only : LIS_rc, LIS_masterproc
    use LIS_fileIOMod,       only : LIS_create_output_directory
    use LIS_logMod,          only : LIS_logunit
    use LIS_historyMod,      only : LIS_writevar_gridded
!
! !DESCRIPTION:
!
! This routine writes the checkpoint data for a SCEUA restart
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
!EOP 
    integer             :: n
    character(len=LIS_CONST_PATH_LEN) :: filen
!    character (len=3)   :: fTrial 
    character (len=6)   :: fTrial 
    real, allocatable       :: GridBestObjFunc(:) ! SY
    real, allocatable       :: GridEnsembleAvgObjFunc(:) ! SY
!    real, allocatable       :: objs(:,:) ! SY
    real, allocatable       :: GridBestParSet(:,:) 
    integer             :: k,iparam, t,m,j,l
! Write the best parameter set from each population (i.e., grid ensemble?). 
! true if subgrid tiling is not turned on. The idea
! here is to write a 'gridded objective function' field. 

    n = 1
    allocate(GridBestObjFunc(LIS_rc%ngrid(n))) ! SY
    allocate(GridEnsembleAvgObjFunc(LIS_rc%ngrid(n))) ! SY
    allocate(GridBestParSet(LIS_rc%ngrid(n),sceua_ctl%nparam)) ! SY
    GridEnsembleAvgObjFunc = 0.0

    if(LIS_masterproc) then
       call LIS_create_output_directory('SCEUA')
!       write(unit=fTrial, fmt='(i3.3)') sceua_ctl%SumNewMemberFuncEvals 
       write(unit=fTrial, fmt='(i6.6)') sceua_ctl%SumNewMemberFuncEvals
       filen = trim(LIS_rc%odir)//'/SCEUA/SCEUA.'&
            //trim(fTrial)//'.1gd4r'
       open(40,file=filen,status='unknown',form='unformatted')
    endif

    do t=1,LIS_rc%ngrid(n)
       GridBestObjFunc(t) = sceua_struc(t)%objfuncs(1)
       do m=1,LIS_rc%nensem(n)
          GridEnsembleAvgObjFunc(t) = GridEnsembleAvgObjFunc(t) + sceua_struc(t)%objfuncs(m)
       enddo
       GridEnsembleAvgObjFunc(t) = GridEnsembleAvgObjFunc(t)/float(LIS_rc%nensem(n))
       do k=1, sceua_ctl%nparam
          GridBestParSet(t,k) = sceua_struc(t)%parsets(1,k) ! SY
       enddo
    enddo
    call LIS_writevar_gridded(40,n,GridBestObjFunc)
    call LIS_writevar_gridded(40,n,GridEnsembleAvgObjFunc)

!write the best parameter set 
    do k=1,sceua_ctl%nparam
       call LIS_writevar_gridded(40,n,GridBestParSet(:,k))
    enddo
    if(LIS_masterproc) then
       close(40)
    endif

    write(LIS_logunit,*) 'GridBestObjFunc is ', GridBestObjFunc
    write(LIS_logunit,*) 'GridEnsembleAvgObjFunc is ', &
                         GridEnsembleAvgObjFunc
    write(LIS_logunit,*) 'GridBestParSet is ', GridBestParSet

    deallocate(GridBestObjFunc)
    deallocate(GridEnsembleAvgObjFunc)
    deallocate(GridBestParSet) 

  end subroutine writeSCEUAoutput


 end module ShuffledComplexEvolution

  
