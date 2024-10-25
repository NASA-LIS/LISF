!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module LevenbergMarquardt
!BOP
!
! !MODULE: LevenbergMarquardt
!
! !DESCRIPTION:  Implementation of the Levenberg-Marquardt (LM)
!  algorithm for minimizing the sum of squared residuals.
!  The algorithm interpolates bewteen two numerical methods:
!  the Gauss-Newton method and the method of steepest descent.  The LM
!  metric, or objective function, must be selected in LIS  as 
!  LM requires at each iteration the vector of residuals
!  (as opposed to the sum of squared residuals). This implementation
!  is based on the MINPACK code (single precision code).  MINPACK does not handle bounds or
!  constraints.  Here, bounds are handled by transforming
!  the values in the bounded space to an unbounded space (logit function); 
!  values are passed to the LSM by transforming back to the bounded space 
!  (sigmoid function, the inverse of the logit).  Constraints are handled 
!  by artificially inflating the residuals when an infeasible 
!  solution is proposed.
!  
!  MINPACK's  lmdif.f and fdjac2.f code are rewritten here to conform to
!  LIS's architecture.
!  The code draws on the following MINPACK and FORTRAN routines:
!       MINPACK-SUPPLIED ... SPMPAR,LMPAR,QRFAC,SNRM2,QRSOLV
!
!       FORTRAN-SUPPLIED ... ABS,AMAX1,AMIN1,SQRT,MOD, MIN0
!  
! !REVISION HISTORY:
!  05 Jul 2008; Ken Harrison; Initial Specification
!
!EOP
  use ESMF
  use LM_varctl
  use LIS_coreMod
  use LIS_optUEMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!MINPACK-DEFINED PARAMETERS
  real, parameter :: ONE = 1.0E0
  real, parameter :: P1 =  1.0E-1
  real, parameter :: P5 = 5.0E-1
  real, parameter :: P25 = 2.5E-1
  real, parameter :: P75 = 7.5E-1
  real, parameter :: P0001 = 1.0E-4
  real, parameter :: ZERO = 0.0E0

!Large penalty for constraint violation
  real, parameter :: LargeNum = 100000

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LMOpt_init
  public :: LMOpt_setup
  public :: LMOpt_run
  public :: LMOpt_checkConvergence
  public :: LMOpt_getdecSpaceValues
!  public :: LM_setdecSpaceValues  !pre-LIS7 way
  public :: LMOpt_getNparam
  public :: LMOpt_readrestart
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
  type lmstruc
 
!LMDIF arg list
      integer           ::    M
      integer           ::    INFO
      integer, allocatable  ::    IPVT(:)
      real, allocatable     ::    X(:)
      real, allocatable     ::    FVEC(:)
      real, allocatable     ::    DIAG(:)
      real, allocatable     ::    FJAC(:,:)
      real, allocatable     ::    QTF(:)
      real, allocatable     ::    WA1(:)
      real, allocatable     ::    WA2(:)
      real, allocatable     ::    WA3(:)
      real, allocatable     ::    WA4(:)

!Step sized used for Jacobian evaluation
      real, allocatable  ::    H(:)

!Variables local to LMDIF
      integer :: I   
      integer :: ITER
      integer :: J   
      integer :: L   
      real ::    ACTRED
      real ::    DELTA
      real ::    DIRDER
      real ::    FNORM
      real ::    FNORM1
      real ::    GNORM
      real ::    PAR
      real ::    PNORM 
      real ::    PRERED
      real ::    RATIO
      real ::    SUM
      real ::    TEMP
      real ::    TEMP1 
      real ::    TEMP2
      real ::    XNORM
      real ::    SPMPAR
!      real ::    ENORM

!Needed to replicate the MINPACK program flow of MINPACK's lmdif.f
      logical :: BeginInOuterLoop
      logical :: IsTerminated
      integer :: NumLISRuns

  end type lmstruc

  type(lmctl)            :: lm_ctl
  type(lmstruc), allocatable :: lm_struc(:)

  contains
   
!BOP
! !ROUTINE: LMOpt_init
! \label{LMOpt_init}
! 
! !INTERFACE: 
    subroutine LMOpt_init()
! !USES: 
     use LIS_logMod,          only : LIS_logunit, LIS_verify, LIS_endrun
! !DESCRIPTION: 
!   This routine performs the initialization steps for the LM.  
!   Arrays not dependent on the number of observations are initialized
!   Those dependent on the number of observations are handled in LMOpt\_setup.
!EOP      
     implicit none 

     character*100,  allocatable     :: vname(:)
     integer                     :: i, t
     integer                     :: ftn
     integer                     :: n 
     integer                     :: status
     type(ESMF_Field)            :: varField  !this is the variable name field
     type(ESMF_ArraySpec)        :: arrspec1

#if (defined USE_MINPACK)
  !nothing
#else
     !call LIS_endrun
     write(LIS_logunit,*) 'No link to MINPACK library'
     write(LIS_logunit,*) 'program stopping ..'
     call LIS_endrun()     
#endif


! currently limited to one nest
      n = 1

!      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%decspaceAttribsFile,&
!           label="LM decision space attributes file:")
!      call LIS_verify(status, 'LM decision space attributes file: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%restart, &
           label="LM start mode:",default=2, rc=status)      
      call LIS_verify(status, 'LM start mode: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%rfile, &
           label="LM restart file:",rc=status)
      call LIS_verify(status, 'LM restart file: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%MMAX, &
           label="LM maximum number of observations:",rc=status)      
      call LIS_verify(status, 'LM maximum number of observations: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%MAXITER, &
           label="LM maximum iterations:",rc=status)      
      call LIS_verify(status, 'LM maximum iterations: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%MODE, &
           label="LM mode:",default=1,rc=status)      
      call LIS_verify(status, 'LM mode: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%FTOL, &
           label="LM objective function tolerance:",default=0.000001,rc=status)      
      call LIS_verify(status, 'LM objective function tolerance: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%XTOL, &
           label="LM decision space tolerance:",default=0.000001,rc=status)      
      call LIS_verify(status, 'LM decision space tolerance: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%GTOL, &
           label="LM orthogonality tolerance:",default=0.0,rc=status)      
      call LIS_verify(status, 'LM orthogonality tolerance: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%FACTOR, &
           label="LM step bound factor:",default=100.0,rc=status)
      call LIS_verify(status, 'LM step bound factor: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,lm_ctl%EPSFCN, &
           label="LM forward difference step length:",default=0.001,rc=status)      
      call LIS_verify(status, 'LM forward difference step length: not defined')

      call LMOpt_setup()
    end subroutine LMOpt_init
!BOP
!
! !ROUTINE: LMOpt_checkConvergence
! \label{LMOpt_checkConvergence}
!
! !INTERFACE:
    subroutine LMOpt_checkConvergence(check)
! !ARGUMENTS:
      logical, intent(INOUT) :: check
!
! !DESCRIPTION:
!  This routine checks to see if the convergence criteria for LM
!  is reached.

!EOP
      check = lm_ctl%check
      !last check may be redundant

    end subroutine LMOpt_checkConvergence


!BOP
! !ROUTINE: LMOpt_setup
! \label{LMOpt_setup}
!
! !INTERFACE: LMOpt_setup
    subroutine LMOpt_setup()
! !USES: 
! 
! !DESCRIPTION: 
! 
!   This subroutine initializes the required memory structures
!   dependent on the number of observations.  This requires
!   an evaluation of the starting solution. 
!EOP
      type(ESMF_Field)           :: varField
      real,   pointer            :: vardata(:)
      integer                    :: k,j,n,t,m
      integer                    :: status
     real                        :: SPMPAR

! Currently supporting only one nest. 
      n = 1

      call ESMF_StateGet(LIS_decisionSpace, itemCount=lm_ctl%nparam, &
           rc=status)
      call LIS_verify(status)

      allocate(lm_ctl%parmax(lm_ctl%nparam))
      allocate(lm_ctl%parmin(lm_ctl%nparam))
      allocate(lm_ctl%vname(lm_ctl%nparam))
!      allocate(lm_ctl%parstart(lm_ctl%nparam))  used to hold params

      allocate(lm_struc(LIS_rc%ntiles(n)/LIS_rc%nensem(n))) 

      if(LIS_rc%nensem(n).ne.lm_ctl%nparam + 1) then 
         write(LIS_logunit,*)'The number of ensembles should be exactly'
         write(LIS_logunit,*) 'N + 1 to run LM.'
         write(LIS_logunit,*) 'Stopping program'
         call LIS_endrun()
      endif

      lm_ctl%check = .false.

! Size arrays dependent only on num parameters
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) !num tiles
         allocate(lm_struc(t)%X(lm_ctl%nparam))
         allocate(lm_struc(t)%DIAG(lm_ctl%nparam))
         allocate(lm_struc(t)%IPVT(lm_ctl%nparam))
         allocate(lm_struc(t)%QTF(lm_ctl%nparam))
         allocate(lm_struc(t)%WA1(lm_ctl%nparam))
         allocate(lm_struc(t)%WA2(lm_ctl%nparam))
         allocate(lm_struc(t)%WA3(lm_ctl%nparam))
         allocate(lm_struc(t)%IPVT(lm_ctl%nparam))
         allocate(lm_struc(t)%H(lm_ctl%nparam))

         call ESMF_StateGet(LIS_decisionSpace, itemNameList=lm_ctl%vname,&
              rc=status)
         call LIS_verify(status)
         
         do k=1,lm_ctl%nparam
            call ESMF_StateGet(LIS_decisionSpace, trim(lm_ctl%vname(k)), &
                 varField, rc=status)
            call LIS_verify(status)
            
            call ESMF_AttributeGet(varField, 'MinRange',lm_ctl%parmin(k),&
                 rc=status)
            call LIS_verify(status, 'setting minrange to decspace obj in LM')
            call ESMF_AttributeGet(varField, 'MaxRange',lm_ctl%parmax(k),&
                 rc=status)
            call LIS_verify(status, 'setting maxrange to decspace obj in LM')
            
            call ESMF_FieldGet(varField,localDE=0, farrayPtr=vardata,rc=status)
            call LIS_verify(status)
            
            m=1  !  Ensemble only to eval jacobian -> only look at 1st ensemble
            lm_struc(t)%X(k) = vardata((t-1)*LIS_rc%nensem(n)+m)
         enddo
         
         ! Transform bounded to unbounded
         call bounded2unbounded(lm_struc(t)%X,&
              lm_ctl%parmin,lm_ctl%parmax,lm_ctl%nparam)
      enddo
   
      !
      !     EPSMCH IS THE MACHINE PRECISION.
      !
#if (defined USE_MINPACK)
      lm_ctl%EPSMCH = SPMPAR(1)
#endif
      !
      !     
      !     CHECK THE INPUT PARAMETERS FOR ERRORS.
      !     
         if (lm_ctl%nparam .LE. 0  &
             .OR. lm_ctl%FTOL .LT. ZERO .OR. lm_ctl%XTOL .LT. ZERO &
             .OR. lm_ctl%GTOL .LT. ZERO .OR. lm_ctl%MAXITER .LE. 0 &
             .OR. lm_ctl%FACTOR .LE. ZERO) then
                  write(LIS_logunit,*)'Error in LM input parameters'
                  write(LIS_logunit,*) 'program stopping ..'
                  call LIS_endrun()
         endif

         do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) !num tiles
            lm_struc(t)%INFO = 0

            !  **Revisit if want to allow scaling of parameters--related to MODE; 
            !  **  if so, then have scaling in dec space file
            !  (The transformation may well  negate scaling issues).
            !         IF (MODE .NE. 2) GO TO 20
            !         DO 10 J = 1, N
            !            IF (DIAG(J) .LE. ZERO) GO TO 300  !terminate as above with check input
            
            !Flow control
            lm_struc(t)%BeginInOuterLoop=.true.
            lm_struc(t)%ITER = 0
         enddo
         lm_ctl%ITER=0
end subroutine LMOpt_setup

!BOP
! !ROUTINE: firstIteration
! \label{firstIteration}
!
! !INTERFACE: firstIteration
    subroutine firstIteration()
! !USES: 
! 
! !DESCRIPTION: 
! 
!   Evaluate starting solution.  With count of observations (by tile), 
!   allocate those arrays dependent on number of observations.
!EOP

      integer            :: i,j,n,t,m
      type(ESMF_Field)   :: errField, numobsField, feasField
      integer, pointer   :: numobs(:), mod_flag(:)
      real, pointer      :: err(:,:)

      real               :: SNRM2
      integer            :: status
! Currently supporting only one nest. 
    n = 1
    call ESMF_StateGet(LIS_ObjectiveFunc,"Number of observations",numobsField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(numobsField, localDE=0, farrayPtr=numobs, rc=status)
    call LIS_verify(status)

    call ESMF_StateGet(LIS_ObjectiveFunc,"Error field",errField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(errField, localDE=0, farrayPtr=err, rc=status)
    call LIS_verify(status)

    call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
    call LIS_verify(status)

    lm_ctl%ITER=1
    do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) !num tiles
       lm_struc(t)%ITER = 1
    enddo   

! Allocate those arrays dependent on num obs
    do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) !num tiles
       lm_struc(t)%M=numobs(t)
       if (lm_ctl%nparam .gt. lm_struc(t)%M) then
                write(LIS_logunit,*)'Error in LM input parameters:'
                write(LIS_logunit,*) 'Num parameters exceeds num obs'
                write(LIS_logunit,*) 'program stopping ..'
                call LIS_endrun()
       endif
       allocate(lm_struc(t)%FVEC(lm_struc(t)%M))
       allocate(lm_struc(t)%WA4(lm_struc(t)%M))
       allocate(lm_struc(t)%FJAC(lm_struc(t)%M,lm_ctl%nparam))
    end do

!Retrieve FVEC from first ensemble for the tile
    do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
       if(mod_flag((t-1)*LIS_rc%nensem(n)+1).eq.1) then !infeasible
          write(LIS_logunit,*) 'Error in LM:'
          write(LIS_logunit,*) 'Starting solution is not feasible'
          write(LIS_logunit,*) 'program stopping ..'
          call LIS_endrun()
       endif

       do i=1,lm_struc(t)%M
          lm_struc(t)%FVEC(i) = err((t-1)*LIS_rc%nensem(n)+1,i)
       enddo
#if (defined USE_MINPACK) 
       lm_struc(t)%FNORM = SNRM2(lm_struc(t)%M,lm_struc(t)%FVEC,1)
#endif
!
!      INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNTER.
!
       lm_struc(t)%PAR = ZERO

       call printCurrent()
    enddo
end subroutine firstIteration


!BOP
! !ROUTINE: jacobian
! \label{jacobian}
!
! !INTERFACE: jacobian
    subroutine jacobian(t)
! !USES: 
! 
! !DESCRIPTION: 
! 
!   This subroutine computes the Jacobian from the N Jacobian-related 
!   ensembles (one for each parameter requiring fitting) associated
!   with the most recently evaluated solution. 
!EOP

      integer            :: i,j,n,t,m, status
      type(ESMF_Field)   :: errField
      real, pointer      :: err(:,:)
      real               :: fprime, f, change

! Currently supporting only one nest. 
    n = 1

    call ESMF_StateGet(LIS_ObjectiveFunc,"Error field",errField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(errField, localDE=0, farrayPtr=err, rc=status)
    call LIS_verify(status)

    !Compute Jacobian from ensemble 2 to nparam+1
    do i=1,lm_struc(t)%M
       f=err((t-1)*LIS_rc%nensem(n)+1,i)
       do j=1,lm_ctl%nparam
          fprime=err((t-1)*LIS_rc%nensem(n)+j+1,i)
          change=fprime-f
          lm_struc(t)%FJAC(i,j) = change/(lm_struc(t)%H(j))
          write(*,*) change
       enddo
    enddo

!    write(*,*) lm_struc(t)%H
!    write(*,*) lm_struc(t)%FJAC
end subroutine jacobian

!BOP
! !ROUTINE: fvecAsW4
! \label{fvecAsW4}
!
! !INTERFACE: fvecAsW4
    subroutine fvecAsW4
! !USES: 
! 
! !DESCRIPTION: 
! 
!   This subroutine computes the fills the residuals array (FVEC) 
!   from the most recent run and computes the norm.
!EOP

      integer            :: i,j,n,t,m, status
      type(ESMF_Field)   :: errField, feasField
      real, pointer      :: err(:,:)
      integer, pointer   :: mod_flag(:)
      real               :: f
      real               :: SNRM2
! Currently supporting only one nest. 
    n = 1

    call ESMF_StateGet(LIS_ObjectiveFunc,"Error field",errField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(errField, localDE=0, farrayPtr=err, rc=status)
    call LIS_verify(status)

    call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
    call LIS_verify(status)

    do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
       do i=1,lm_struc(t)%M
          if(mod_flag((t-1)*LIS_rc%nensem(n)+1).eq.1) then !infeasible
             lm_struc(t)%WA4(i) = LargeNum
          else !feasible
             f=err((t-1)*LIS_rc%nensem(n)+1,i)
             lm_struc(t)%WA4(i)  =f
          endif
       enddo
#if(defined USE_MINPACK)
       lm_struc(t)%FNORM1=SNRM2(lm_struc(t)%M,lm_struc(t)%WA4,1)    
#endif
    enddo
end subroutine fvecAsW4

!BOP
! !ROUTINE: LMOpt_getdecSpaceValues
! \label{LMOpt_getdecSpaceValues}
!
! !INTERFACE: LMOpt_getdecSpaceValues
    subroutine LMOpt_getdecSpaceValues(n)
! !USES: 
! 
! !DESCRIPTION: 
! 
!   This subroutine supplies the proposed solution to the LSM
!EOP
      
      implicit none

      real               :: EPS
      integer            :: n
      real               :: Z(lm_ctl%nparam)
      real               :: temp(lm_ctl%nparam)
      real               :: Y(lm_ctl%nparam)
      real               :: decvals(lm_ctl%nparam, LIS_rc%ntiles(n))
      integer            :: i, t, m

      type(ESMF_Field)   :: varField(lm_ctl%nparam)
      real, pointer      :: vardata(:)
      integer            :: status

      do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         !Copy 'X' array into temp; 1st iter starting X; else proposed new X which is WA2
         do i=1, lm_ctl%nparam
            if (lm_ctl%ITER .le. 1) then
                  temp(i)=lm_struc(t)%X(i)
            else
                  temp(i)=lm_struc(t)%WA2(i)
            endif
         enddo     

         !Transform to original bounded values; paste across all ensembles
         call unbounded2bounded(temp,lm_ctl%parmin,lm_ctl%parmax,lm_ctl%nparam,Y)
         do i=1, lm_ctl%nparam
            do m=1, lm_ctl%nparam+1
               call ESMF_StateGet(LIS_decisionSpace, trim(lm_ctl%vname(i)), &
                    varField(i), rc=status)
               call LIS_verify(status)
               
               call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata,rc=status)
               call LIS_verify(status)
               
               !decvals(i,(t-1)*LIS_rc%nensem(n)+m) = Y(i)
               vardata((t-1)*LIS_rc%nensem(n)+m)= Y(i)
            enddo
         enddo
      
!Jacobian-related ensembles
         ! Determine step size H for each param; add to make Z
         do i=1, lm_ctl%nparam
!            EPS = SQRT(AMAX1(lm_ctl%EPSFCN,lm_ctl%EPSMCH))
!            lm_struc(t)%H(i)=EPS*ABS(temp(i))
!            if (lm_struc(t)%H(i) .eq. ZERO) lm_struc(t)%H(i) = EPS
            lm_struc(t)%H(i)=lm_ctl%EPSFCN
            Z(i)=temp(i)+lm_struc(t)%H(i)
         enddo

         !Transform perturbation vector to bounded values
         call unbounded2bounded(Z,lm_ctl%parmin,lm_ctl%parmax,lm_ctl%nparam,Y)

         !Stick in perturbation vector
         do i=1, lm_ctl%nparam
            m=i+1
            call ESMF_StateGet(LIS_decisionSpace, trim(lm_ctl%vname(i)), &
                 varField(i), rc=status)
            call LIS_verify(status)
            
            call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata,rc=status)
            call LIS_verify(status)

            !decvals(i,(t-1)*LIS_rc%nensem(n)+m) = Y(i)
            vardata((t-1)*LIS_rc%nensem(n)+m)= Y(i)
         enddo
      enddo
    end subroutine LMOpt_getdecSpaceValues


!BOP
! !ROUTINE: LMOpt_run
! \label{LMOpt_run}
! 
! !INTERFACE: 
    subroutine LMOpt_run()
! !USES: 
! !DESCRIPTION: 
!   This routine performs the run steps for the LM.  
!   LIS has been run once, in LMOpt\_setup, prior to this point.
!EOP  
      integer            :: t
      integer            :: n
      integer            :: m
      integer            :: status
      integer            :: i

      n = 1


      write(LIS_logunit,*) 'Running LM, generation: ',lm_ctl%ITER

      if((lm_ctl%ITER.eq.0).or.(lm_ctl%restart.eq.1)) then 
         call firstIteration()
         do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            lm_struc(t)%IsTerminated=.false.
         enddo
      else
         !         lm_ctl%AllTerminated=.false.
         !         do while (.not.lm_ctl%AllTerminated)
         if(lm_ctl%ITER.gt.1) then 
            !Function evaluation
            !         call LM_reset_err
            !            call setoptuetypedecspace(LIS_rc%optuetype)
            !            call evaluateobjfunction(LIS_rc%optuetype)
            call fvecAsW4

            !After function evaluation
            do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
               if (.not.lm_struc(t)%IsTerminated) then
                  call FCNEVAL2END(t)
                  if (.not.lm_struc(t)%IsTerminated) then
                     lm_ctl%AllTerminated=.false.
                  endif
               else
                  !nothing
               endif
               if (t .eq. 1) call printCurrent()
            enddo

            !call write restart and write X as it stands
            !will be starting X on restart

         endif

         !Write restart
         call writeLMrestart

         !Set termination flag
         if (lm_ctl%AllTerminated.or.(lm_ctl%ITER.eq.lm_ctl%MAXITER)) then
            lm_ctl%check=.true.
         endif

         !reset AllTerminated so can be flipped on next run
         lm_ctl%AllTerminated=.true.
         lm_ctl%ITER = lm_ctl%ITER + 1

         !Before function evaluation
         do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            if (.not.lm_struc(t)%IsTerminated) then
               if (lm_struc(t)%BeginInOuterLoop) then
                  call OUTER2INNER(t)
               else
                  !Begin in inner loop
               endif
               call INNER2FCNEVAL(t)
            else
               !simply run as left off
            endif
         enddo
      endif
    end subroutine LMOpt_run


!BOP
! !ROUTINE: OUTER2INNER
! \label{OUTER2INNER}
! 
! !INTERFACE: 
    subroutine OUTER2INNER(t)
! !USES: 
! !DESCRIPTION: 
!   This routine carries out the MINPACK code
!   from the beginning of the MINPACK 'outer loop' to the
!   beginning of the 'inner loop'.
!EOP  
       integer    :: I,J,t
       real       :: SNRM2
!
!        CALCULATE THE JACOBIAN MATRIX.
!
         call jacobian(t)

!
!        IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
!
         !Minpack prints at this time--presumably to avoid printing unsuccessful solutions
!
!        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
!
#if (defined USE_MINPACK)
         CALL QRFAC(lm_struc(t)%M,lm_ctl%nparam, lm_struc(t)%FJAC,&
                      lm_struc(t)%M, .TRUE.,lm_struc(t)%IPVT,&
                      lm_ctl%nparam, lm_struc(t)%WA1,&
                      lm_struc(t)%WA2, lm_struc(t)%WA3)
#endif
!
!      ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
!      TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
!
       if (lm_struc(t)%ITER .eq. 1) then
          if (lm_ctl%MODE .eq. 1) then
             do j=1,lm_ctl%nparam
                lm_struc(t)%DIAG(j)=lm_struc(t)%WA2(j)
                if (lm_struc(t)%WA2(j) .eq. ZERO) lm_struc(t)%DIAG(j)=ONE
             enddo
          endif
   !
   !      ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
   !      AND INITIALIZE THE STEP BOUND DELTA.
   !
          do j=1,lm_ctl%nparam
             lm_struc(t)%WA3(j) = lm_struc(t)%DIAG(j)*lm_struc(t)%X(j)
          enddo
#if(defined USE_MINPACK)
          lm_struc(t)%XNORM = SNRM2(lm_ctl%nparam,lm_struc(t)%WA3,1)
#endif
          lm_struc(t)%DELTA = lm_ctl%FACTOR*lm_struc(t)%XNORM
          if (lm_struc(t)%DELTA .eq. ZERO) lm_struc(t)%DELTA = lm_ctl%FACTOR
       endif
!
!        FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN
!        QTF.
!
         do I = 1, lm_struc(t)%M
            lm_struc(t)%WA4(I) = lm_struc(t)%FVEC(I)
         enddo
         do J = 1, lm_ctl%nparam
            IF (lm_struc(t)%FJAC(J,J) .EQ. ZERO) then
                !nothing
            else
                lm_struc(t)%SUM = ZERO
                do I = J, lm_struc(t)%M
                   lm_struc(t)%SUM = lm_struc(t)%SUM + lm_struc(t)%FJAC(I,J)*lm_struc(t)%WA4(I)
                enddo
                lm_struc(t)%TEMP = -lm_struc(t)%SUM/lm_struc(t)%FJAC(J,J)
                do I = J, lm_struc(t)%M
                   lm_struc(t)%WA4(I) = lm_struc(t)%WA4(I) + lm_struc(t)%FJAC(I,J)*lm_struc(t)%TEMP
                enddo
            endif
            lm_struc(t)%FJAC(J,J) = lm_struc(t)%WA1(J)
            lm_struc(t)%QTF(J) = lm_struc(t)%WA4(J)
         enddo
!
!        COMPUTE THE NORM OF THE SCALED GRADIENT.
!
         lm_struc(t)%GNORM = ZERO
         if (lm_struc(t)%FNORM .EQ. ZERO) then
            !nothing
         else
            do J = 1, lm_ctl%nparam
               lm_struc(t)%L = lm_struc(t)%IPVT(J)
               IF (lm_struc(t)%WA2(lm_struc(t)%L) .EQ. ZERO) then
                  !nothing
               else
                  lm_struc(t)%SUM = ZERO
                  do I = 1, J
                     lm_struc(t)%SUM = lm_struc(t)%SUM + lm_struc(t)%FJAC(I,J)*(lm_struc(t)%QTF(I)/lm_struc(t)%FNORM)
                  enddo
                  lm_struc(t)%GNORM = AMAX1(lm_struc(t)%GNORM,ABS(lm_struc(t)%SUM/lm_struc(t)%WA2(lm_struc(t)%L)))
               endif
            enddo
         endif
!
!        TEST FOR CONVERGENCE OF THE GRADIENT NORM.
!
         IF (lm_struc(t)%GNORM .LE. lm_ctl%GTOL) lm_struc(t)%INFO = 4
         IF (lm_struc(t)%INFO .NE. 0) then
            !do nothing
         else ! info=0
!
!        RESCALE IF NECESSARY.
!
            IF (lm_ctl%MODE .EQ. 2) then
               !do nothing
            else
               DO J = 1, lm_ctl%nparam
                  lm_struc(t)%DIAG(J) = AMAX1(lm_struc(t)%DIAG(J),lm_struc(t)%WA2(J))
               enddo
            endif
         endif

    end subroutine OUTER2INNER

!BOP
! !ROUTINE: INNER2FCNEVAL
! \label{INNER2FCNEVAL}
! 
! !INTERFACE: 
    subroutine INNER2FCNEVAL(t)
! !USES: 
! !DESCRIPTION: 
!   This routine carries out the MINPACK code
!   from the beginning of the MINPACK 'inner loop' to the
!   proposed solution evaluation step of the 'inner loop'.
!EOP  
       integer    :: J,t
       real       :: SNRM2
!
!           DETERMINE THE LEVENBERG-MARQUARDT PARAMETER.
!
#if(defined USE_MINPACK)
            CALL LMPAR(lm_ctl%nparam,lm_struc(t)%FJAC, &
                         lm_struc(t)%M,lm_struc(t)%IPVT, &
                         lm_struc(t)%DIAG, lm_struc(t)%QTF, &
                         lm_struc(t)%DELTA,lm_struc(t)%PAR, &
                         lm_struc(t)%WA1,lm_struc(t)%WA2, &
                         lm_struc(t)%WA3,lm_struc(t)%WA4)
#endif
!
!           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
!
            DO J = 1, lm_ctl%nparam
               lm_struc(t)%WA1(J) = -lm_struc(t)%WA1(J)
               lm_struc(t)%WA2(J) = lm_struc(t)%X(J) + lm_struc(t)%WA1(J)
               lm_struc(t)%WA3(J) = lm_struc(t)%DIAG(J)*lm_struc(t)%WA1(J)
            enddo
#if(defined USE_MINPACK)
            lm_struc(t)%PNORM = SNRM2(lm_ctl%nparam,lm_struc(t)%WA3,1)
#endif
!
!           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
!
            IF (lm_struc(t)%ITER .EQ. 1) then
               lm_struc(t)%DELTA = AMIN1(lm_struc(t)%DELTA,lm_struc(t)%PNORM)
            endif
    end subroutine INNER2FCNEVAL

!BOP
! !ROUTINE: FCNEVAL2END
! \label{FCNEVAL2END}
! 
! !INTERFACE: 
    subroutine FCNEVAL2END(t)
! !USES: 
! !DESCRIPTION: 
!   This routine carries out the MINPACK code
!   from the proposed solution evaluation of the MINPACK 'inner loop' 
!   to the end of the 'inner loop'.
!EOP  
       integer    :: I,J,t
       real       :: SNRM2
!
!           COMPUTE THE SCALED ACTUAL REDUCTION.
!
            lm_struc(t)%ACTRED = -ONE
            IF (P1*lm_struc(t)%FNORM1 .LT. lm_struc(t)%FNORM) then
               lm_struc(t)%ACTRED = ONE - (lm_struc(t)%FNORM1/lm_struc(t)%FNORM)**2
            end if
!
!           COMPUTE THE SCALED PREDICTED REDUCTION AND
!           THE SCALED DIRECTIONAL DERIVATIVE.
!
            DO J = 1, lm_ctl%nparam
               lm_struc(t)%WA3(J) = ZERO
               lm_struc(t)%L = lm_struc(t)%IPVT(J)
               lm_struc(t)%TEMP = lm_struc(t)%WA1(lm_struc(t)%L)
               DO I = 1, J
                  lm_struc(t)%WA3(I) = lm_struc(t)%WA3(I) &
                     + lm_struc(t)%FJAC(I,J)*lm_struc(t)%TEMP
               enddo
            enddo
#if(defined USE_MINPACK)
            lm_struc(t)%TEMP1 = SNRM2(lm_ctl%nparam,lm_struc(t)%WA3,1)/lm_struc(t)%FNORM
#endif
            lm_struc(t)%TEMP2 = (SQRT(lm_struc(t)%PAR)*lm_struc(t)%PNORM)/lm_struc(t)%FNORM
            lm_struc(t)%PRERED = lm_struc(t)%TEMP1**2 + lm_struc(t)%TEMP2**2/P5
            lm_struc(t)%DIRDER = -(lm_struc(t)%TEMP1**2 + lm_struc(t)%TEMP2**2)
!
!           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
!           REDUCTION.
!
            lm_struc(t)%RATIO = ZERO
            IF (lm_struc(t)%PRERED .NE. ZERO) then
               lm_struc(t)%RATIO = lm_struc(t)%ACTRED/lm_struc(t)%PRERED
            endif
!
!           UPDATE THE STEP BOUND.
!
            IF (lm_struc(t)%RATIO .GT. P25) then
               IF (lm_struc(t)%PAR .NE. ZERO .AND. lm_struc(t)%RATIO .LT. P75) then
                  !nothing
               else
                  lm_struc(t)%DELTA = lm_struc(t)%PNORM/P5
                  lm_struc(t)%PAR = P5*lm_struc(t)%PAR
               endif
            else
               IF (lm_struc(t)%ACTRED .GE. ZERO) then
                  lm_struc(t)%TEMP = P5
               endif
               IF (lm_struc(t)%ACTRED .LT. ZERO) then
                  lm_struc(t)%TEMP = P5*lm_struc(t)%DIRDER/(lm_struc(t)%DIRDER + P5*lm_struc(t)%ACTRED)
               endif
               IF (P1*lm_struc(t)%FNORM1 .GE. lm_struc(t)%FNORM .OR. lm_struc(t)%TEMP .LT. P1) then
                  lm_struc(t)%TEMP = P1
               endif
               lm_struc(t)%DELTA = lm_struc(t)%TEMP*AMIN1(lm_struc(t)%DELTA,lm_struc(t)%PNORM/P1)
               lm_struc(t)%PAR = lm_struc(t)%PAR/lm_struc(t)%TEMP
            endif
!
!           TEST FOR SUCCESSFUL ITERATION.
!
            IF (lm_struc(t)%RATIO .LT. P0001) then
               ! Is unsuccessful iteration
               lm_struc(t)%BeginInOuterLoop=.false.
            else
               lm_struc(t)%BeginInOuterLoop=.true.               
   !           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
   !
               DO J = 1, lm_ctl%nparam
                  lm_struc(t)%X(J) = lm_struc(t)%WA2(J)
                  lm_struc(t)%WA2(J) = lm_struc(t)%DIAG(J)*lm_struc(t)%X(J)
               enddo
               DO I = 1, lm_struc(t)%M
                  lm_struc(t)%FVEC(I) = lm_struc(t)%WA4(I)
               enddo
#if(defined USE_MINPACK)
               lm_struc(t)%XNORM = SNRM2(lm_ctl%nparam,lm_struc(t)%WA2,1)
#endif
               lm_struc(t)%FNORM = lm_struc(t)%FNORM1
               lm_struc(t)%ITER = lm_struc(t)%ITER + 1
            endif
!
!           TESTS FOR CONVERGENCE.
!
            IF (ABS(lm_struc(t)%ACTRED) .LE. lm_ctl%FTOL .AND. lm_struc(t)%PRERED .LE. lm_ctl%FTOL &
               .AND. P5*lm_struc(t)%RATIO .LE. ONE) then
               lm_struc(t)%INFO = 1
            endif
            IF (lm_struc(t)%DELTA .LE. lm_ctl%XTOL*lm_struc(t)%XNORM) then
               lm_struc(t)%INFO = 2
            endif
            IF (ABS(lm_struc(t)%ACTRED) .LE. lm_ctl%FTOL .AND. lm_struc(t)%PRERED .LE. lm_ctl%FTOL &
                .AND. P5*lm_struc(t)%RATIO .LE. ONE .AND. lm_struc(t)%INFO .EQ. 2) then
               lm_struc(t)%INFO = 3
            endif
            IF (lm_struc(t)%INFO .NE. 0) lm_struc(t)%IsTerminated = .true.
            if (lm_struc(t)%IsTerminated) then
               !nothing
            else
   !
   !           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
   !
   !            IF (NFEV .GE. MAXFEV) INFO = 5
               IF (ABS(lm_struc(t)%ACTRED) .LE. lm_ctl%EPSMCH .AND. lm_struc(t)%PRERED .LE. lm_ctl%EPSMCH &
                   .AND. P5*lm_struc(t)%RATIO .LE. ONE) lm_struc(t)%INFO = 6
               IF (lm_struc(t)%DELTA .LE. lm_ctl%EPSMCH*lm_struc(t)%XNORM) lm_struc(t)%INFO = 7
               IF (lm_struc(t)%GNORM .LE. lm_ctl%EPSMCH) lm_struc(t)%INFO = 8
            endif
            IF (lm_struc(t)%INFO .NE. 0) lm_struc(t)%IsTerminated = .true.
    end subroutine FCNEVAL2END

!BOP
! 
! !ROUTINE: LMOpt_getNparam
! \label{LMOpt_getNparam}
! 
! !INTERFACE: 
    subroutine LMOpt_getNparam(nparam)
! !USES: 
! !DESCRIPTION: 
!   This routine returns the number of parameters
!   requiring fitting.
!EOP  
      
      integer   :: nparam
      
      nparam = lm_ctl%nparam

    end subroutine LMOpt_getNparam


!BOP
! 
! !ROUTINE: writeLMrestart
! \label{writeLMrestart}
! 
! !INTERFACE: 
  subroutine writeLMrestart
! !USES: 
    use LIS_fileIOMod, only : LIS_create_output_directory
    use LIS_historyMod, only : LIS_writevar_restart
! 
! !DESCRIPTION: 
! 
! This routine writes the checkpoint data for an LM restart.
! 
!EOP
    integer             :: n 
    integer             :: i,t,m
    integer             :: status
    character(len=LIS_CONST_PATH_LEN) :: filen
    character (len=4)   :: fgen
    character*100       :: vnames(lm_ctl%nparam)
    type(ESMF_Field)    :: varField(lm_ctl%nparam)
    real, pointer       :: vardata(:)
    real, allocatable       :: vardata2(:)
    real                :: Y(lm_ctl%nparam)

    n = 1

    call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
         rc=status)
    call LIS_verify(status)
    
    if(LIS_masterproc) then 
       call LIS_create_output_directory('LM')
       write(unit=fgen, fmt='(i4.4)') lm_ctl%ITER
       filen = trim(LIS_rc%odir)//'/LM/LM.'&
            //trim(fgen)//'.LMrst'
       open(40,file=filen,status='unknown',form='unformatted')
       write(40) lm_ctl%ITER
    endif

!    allocate(vardata(LIS_rc%ntiles(n)/LIS_rc%nensem(n))

    do i=1,lm_ctl%nparam
       call ESMF_StateGet(LIS_decisionSpace, trim(vnames(i)), &
            varField(i), rc=status)
       call LIS_verify(status)
       
       call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata, rc=status)
       call LIS_verify(status)
       
       !vardata will contain unsuccessful trial solutions; just using tile layout
       !therefore copy vardata into vardata2
       !and swap in the current solution lm_struc(t)%X
       !into all slots and then write; note: completely wasting rest of ensemble
       vardata2=vardata
       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
          call unbounded2bounded(lm_struc(t)%X,lm_ctl%parmin, lm_ctl%parmax,lm_ctl%nparam,Y)
          do m=1,LIS_rc%nensem(n)
             vardata2((t-1)*LIS_rc%nensem(n)+m) = Y(i)
          enddo
       enddo

       call LIS_writevar_restart(40,n,vardata2)
    enddo
    
    if(LIS_masterproc) then 
       close(40)
       write(LIS_logunit,*) 'LM checkpoint file written ',trim(filen)
    endif
  end subroutine writeLMrestart

!BOP
! 
! !ROUTINE: LMOpt_readrestart
! \label{LMOpt_readrestart}
! 
! !INTERFACE: 
  subroutine LMOpt_readrestart
! !USES: 
    use LIS_historyMod,      only : LIS_readvar_restart
! 
! !DESCRIPTION: 
! 
!   This routine reads the checkpoint data for a LM restart
!EOP
    integer             :: n, m
    integer             :: k, t 
    integer             :: status
    character*100       :: vnames(lm_ctl%nparam)
    type(ESMF_Field)    :: varField(lm_ctl%nparam)
    real, pointer       :: vardata(:)
    real, allocatable       :: vardataX(:)
    real                :: SPMPAR

    if(lm_ctl%restart.eq.1) then !restart run
       
       write(LIS_logunit,*) 'Reading the LM restart file ..'
       n = 1
       call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
            rc=status)
       call LIS_verify(status)
       
       open(40,file=lm_ctl%rfile,form='unformatted')        
       
       read(40) lm_ctl%ITER
       write(LIS_logunit,*) 'Iteration number ',lm_ctl%ITER

       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
          do k=1,lm_ctl%nparam
             call ESMF_StateGet(LIS_decisionSpace, trim(vnames(k)), &
                  varField(k), rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField(k),localDE=0, farrayPtr=vardata,rc=status)
             call LIS_verify(status)
             
             call LIS_readvar_restart(40,n,vardataX)       
             
             m=1 !  Ensemble only to eval jacobian -> only look at 1st ensemble
             lm_struc(t)%X=vardataX((t-1)*LIS_rc%nensem(n)+m)
             vardata((t-1)*LIS_rc%nensem(n)+m) = lm_struc(t)%X(k)
          enddo
          
          !  REMAINDER OF SUBROUTINE SHOULD MATCH THE CODE IN OPT_SETUP (THAT ASSUMED STARTING SOL'N NOT RESTART) 
          ! Transform bounded to unbounded
          call bounded2unbounded(lm_struc(t)%X,&
               lm_ctl%parmin,lm_ctl%parmax,lm_ctl%nparam)
       enddo
       !
       !     EPSMCH IS THE MACHINE PRECISION.
       !
#if (defined USE_MINPACK)
       lm_ctl%EPSMCH = SPMPAR(1)
#endif
    !
       !     
       !     CHECK THE INPUT PARAMETERS FOR ERRORS.
       !     
       if (lm_ctl%nparam .LE. 0  &
            .OR. lm_ctl%FTOL .LT. ZERO .OR. lm_ctl%XTOL .LT. ZERO &
            .OR. lm_ctl%GTOL .LT. ZERO .OR. lm_ctl%MAXITER .LE. 0 &
            .OR. lm_ctl%FACTOR .LE. ZERO) then
          write(LIS_logunit,*)'Error in LM input parameters'
          write(LIS_logunit,*) 'program stopping ..'
          call LIS_endrun()
       endif
       
       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) !num tiles
          lm_struc(t)%INFO = 0
          
          !  **Revisit if want to allow scaling of parameters--related to MODE; 
          !  **  if so, then have scaling in dec space file
          !  (The transformation may well  negate scaling issues).
          !         IF (MODE .NE. 2) GO TO 20
          !         DO 10 J = 1, N
          !            IF (DIAG(J) .LE. ZERO) GO TO 300  !terminate as above with check input
          
          !Flow control
          lm_struc(t)%BeginInOuterLoop=.true.
          lm_struc(t)%ITER = 0
       enddo
       lm_ctl%ITER=0
       
       write(LIS_logunit,*) 'Finished reading the LM restart file ..'
       close(40)
    endif
  end subroutine LMOpt_readrestart


!BOP
! 
! !ROUTINE: writeLMoutput
! \label{writeLMoutput}
! 
! !INTERFACE: 
  subroutine writeLMoutput
! !USES: 
    use LIS_fileIOMod,       only : LIS_create_output_directory
    use LIS_historyMod,      only : LIS_writevar_gridded
! 
! !DESCRIPTION: 
! 
! This routine writes the RMSE value in gridded format
! 
!EOP
    integer             :: n 
    character(len=LIS_CONST_PATH_LEN) :: filen
    character (len=3)   :: fgen
    real, allocatable       :: rmse(:)
    real, allocatable       :: objs(:,:)
    integer             :: k,iparam, t,m,j,l

    n = 1
    allocate(rmse(LIS_rc%ngrid(n)))
    allocate(objs(LIS_rc%ngrid(n),lm_ctl%nparam))

    if(LIS_masterproc) then 
       call LIS_create_output_directory('LM')
       write(unit=fgen, fmt='(i4.4)') lm_ctl%ITER
       filen = trim(LIS_rc%odir)//'/LM/LM.'&
            //trim(fgen)//'.1gd4r'
       open(40,file=filen,status='unknown',form='unformatted')
    endif

    do t=1,LIS_rc%ngrid(n)
       rmse(t) = lm_struc(t)%FNORM/sqrt(real(lm_ctl%nparam))
    enddo
    call LIS_writevar_gridded(40,n,rmse)

    if(LIS_masterproc) then 
       close(40)
    endif
    deallocate(rmse)
  end subroutine writeLMoutput

!BOP
! 
! !ROUTINE: bounded2unbounded
! \label{bounded2unbounded}
!
! !INTERFACE: 
  subroutine bounded2unbounded(X, a, b, n)
! !ARGUMENTS: 
    integer  :: n
    real     :: a(n),b(n)
    real     :: X(n)
! 
! !DESCRIPTION: 
!
!  This routine converts an array X of values
!  bounded by a(i) and b(i) to array Y bounded by -inf to +inf.  
!  The logit function is used for this transformation.
! 
!EOP
    integer  :: i

!   Transform to range 0 to 1
    do i=1,n
       X(i) = (X(i)-a(i))/(b(i)-a(i))
    enddo

!   Apply logit function
    do i=1,n
       X(i) = log(X(i))-log(1-X(i))
    enddo
  end subroutine bounded2unbounded

!BOP
! 
! !ROUTINE: unbounded2bounded
! \label{unbounded2bounded}
!
! !INTERFACE: 
  subroutine unbounded2bounded(X, a, b, n, Y)
! !ARGUMENTS: 
    integer  :: n
    real     :: a(n),b(n)
    real     :: X(n), Y(n)
! 
! !DESCRIPTION: 
!
!  This routine takes unbounded array of real values and converts
!  them to a bounded range from a(i) to b(i). 
!  The logistic function, the inverse of the logit, is used.
!
!EOP
    integer  :: i
    real, parameter     :: AA=14.0
!   Apply logistic function
    do i=1,n
       if (-X(i) .gt. AA) then
          Y(i)=0.0
       else
          Y(i) = 1/(1+exp(-X(i)))
       endif
    enddo

!   Scale to range a to b
    do i=1,n
       Y(i) = a(i) + Y(i)*(b(i)-a(i))
    enddo
  end subroutine unbounded2bounded

!BOP
! 
! !ROUTINE: printCurrent
! \label{printCurrent}
!
! !INTERFACE: 
  subroutine printCurrent()
! !ARGUMENTS: 
    integer  :: t

! 
! !DESCRIPTION: 
!
!  This temporary debugging routine prints
!  X, FVEC, and FNORM
!
!EOP

    real  :: Y(lm_ctl%nparam)
    t=1
    call unbounded2bounded(lm_struc(t)%X,lm_ctl%parmin, lm_ctl%parmax,lm_ctl%nparam,Y)
    write (*,*) Y
    write (*,*) lm_struc(t)%FVEC
    write (*,*) lm_struc(t)%FNORM/sqrt(real(lm_struc(t)%M))
    write (*,*) lm_struc(t)%INFO
  end subroutine printCurrent

!!$! This needs to be eliminated in the 'merge'
!!$    subroutine LM_setdecSpaceValues(n, decvals)
!!$    integer :: n
!!$    real    :: decvals(lm_ctl%nparam, LIS_rc%ntiles(n))
!!$    end subroutine LM_setdecSpaceValues

end module LevenbergMarquardt
