!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module PobjFunc_Mod
!BOP
! !MODULE: PobjFunc_Mod
! 
! !DESCRIPTION: 
!  This module provides the objects and methods to compute a "(P)osterior probability" 
!  metric useful for Bayesian algorithms (and uncertainty propagation).  As described below,while 
!  related to posterior probability, it is not equated to posterior probability.  
!  Posterior probability is defined as \verb"p(theta|y)", where theta is the vector of parameters for which uncertainty
!  estimation is sought (think parameters requiring calibration) and y is data (e.g., remote sensing data).  
!  Its evaluation is dictated by Bayes' rule of probability:
!
! \begin{verbatim}
!        p(theta|y)=p(y|theta)*p(theta)
!                   -------------------
!                        p(y)
! \end{verbatim}
!
!  where \verb"p(y|theta)" is the likelihood of data y given theta; p(theta) is the initial, or prior, probability of theta;
!  and p(y) is the probability of observing y.
!
!  The metric used is the natural log of the numerator of Bayes' rule, and thus ignores the term,  
!  p(y), in the denominator.  p(y) can be thought of as a 
!  normalizing constant to ensure that the 
!  posterior probability integrated over all theta is equal to one.  The difficulty in carrying out the integration
!  to compute p(y) (=int over theta of \verb"p(y|theta)*p(theta)") 
!  is what motivates the development of the Bayesian algorithms that use this metric.  All such algorithms involve
!  computing the ratio of the posterior probability (density) of two values of theta. 
!  In this ratio-ing, the problematic denominator conveniently cancels out.  Therefore
!  this metric only concerns the numerator of Bayes' rule.
!
!  We work in log space as the 
!  the evaluation of the likelihood term --\verb"p(y|theta)"-- is subject to truncation error (involving 
!  much product-ing of very low values, any one of which may be truncated to zero). Moreover, the operations,
!  specifically, the aformentioned ratio-ing can be carried out conveniently in log space.
! 
! !REVISION HISTORY: 
!  
! 7 Sep 2010: Ken Harrison
! 
! !USES:
  use ESMF

  implicit none

  PRIVATE
  integer, parameter :: dist_uniformID     = 0
  integer, parameter :: dist_normalID      = -1
  integer, parameter :: dist_lognormalID   = -2
  integer, parameter :: dist_betaID        = -3
  integer, parameter :: dist_triangularID  = -4
  integer, parameter :: dist_discreteID    = -5
  integer, parameter :: dist_mvnormalID    = -6

  real, parameter :: mindensity=-1e+20
  real, parameter :: PI=3.14159265359
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public ::  initializePObjFunc
  public ::  dist_sample
  public ::  dist_lnprob

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
  public :: distctl
  public :: dist_ctl
  public :: diststruc
  public :: dist_struc

  type diststruc
     integer :: dist_ID
     character*100, allocatable :: param_name(:)
     integer :: numparam
     integer, allocatable :: param(:)
     integer :: numstatparam
     real, allocatable :: statparam(:)  
  end type diststruc

  type distctl
     integer :: ndists
  end type distctl

  type(diststruc), allocatable :: dist_struc(:) !for now, assume same across tilespace
  type(distctl) :: dist_ctl 
  integer :: Pobjfunc_minobs  !unused for now; less problem for P than LS

contains

!BOP
! !ROUTINE: intializePObjFunc
! \label{initializePObjFunc}
! 
! !INTERFACE: 
  subroutine initializePObjFunc()
! !USES: 
    use LIS_coreMod,         only : LIS_vecTile, LIS_config
    use LIS_optUEMod,        only : LIS_ObjectiveFunc, LIS_decisionSpace
    use LIS_logMod,          only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify
    use LIS_constantsMod,    only : LIS_CONST_PATH_LEN
! 
! !DESCRIPTION:
!  This method initializes the objects to be used in the posterior probability
!  metric for use by parameter and uncertainty estimation routines.
!EOP    
    implicit none

    type(ESMF_ArraySpec)   :: arrspec1
    type(ESMF_Field)       :: objfuncField
    type(ESMF_Field)       :: minField
    type(ESMF_Field)       :: maxField
    real,   pointer        :: objfunc(:)
    integer                :: n 
    integer                :: status
    integer                :: ftn
    integer                :: i,j,k,m,kk,t,p !k is #params whereas p is # values to read in
    real                   :: rand
    real                   :: seed
    character*100, allocatable :: vname(:)
    character(len=LIS_CONST_PATH_LEN) :: priorAttribsFile  !link to lis.config entry for optUE
    character*200          :: line
    integer                :: ios
    integer                :: nvars  ! from LIS_decisionSpace as check
    integer                :: var_count !count of variables covered by prior; should equal

    n = 1

    seed = -1000

    ! Read prior distribution attributes file; assume variables
    ! are defined in a decision space file; here referred to solely by index

!!!!    call ESMF_ConfigGetAttribute(LIS_config,Pobjfunc_minobs, &
!!!!         label="Probability objective function minimum number of obs:",rc=status)
!!!!    call LIS_verify(status, &
!!!!         'Probability objective function minimum number of obs:')
!!!!
    call ESMF_ConfigGetAttribute(LIS_config,priorAttribsFile,&
         label="Prior distribution attributes file:",rc=status)
    call LIS_verify(status, 'Prior distribution attributes file: not defined')

    ftn = LIS_getNextUnitNumber()
    write(LIS_logunit,*) 'Reading prior distribution attributes ...', &
         trim(priorAttribsFile)
    open(ftn,file=trim(priorAttribsFile),status='old')
    read(ftn,*)  !demarcate section
    read(ftn,*) dist_ctl%ndists
    read(ftn,*) line!demarcate section


    allocate(dist_struc(dist_ctl%ndists))

    ! Tailor reading on distribution type
    ! Multivariate and discrete can be accommodated
    var_count=0
    do j=1,dist_ctl%ndists
       read(ftn,'(a)') line
!This fails with gfortran compilation ... need to be fixed.
       print*, 'the code fails with gfortran'
       print*, 'stopping purposely'
       stop
!       read(line,'(i)',iostat=ios) dist_struc(j)%dist_ID

!Determine # of parameters over which distribution applies, allocate and read
       select case (dist_struc(j)%dist_ID)
       case (dist_mvnormalID)
!          read(ftn,*) dist_struc(j)%numparam
!This fails with gfortran compilation ... need to be fixed.
       print*, 'the code fails with gfortran'
       print*, 'stopping purposely'
       stop

!          read(line,'(i)',iostat=ios) dist_struc(j)%dist_ID, dist_struc(j)%numparam
       case default
          dist_struc(j)%numparam=1
       end select
       allocate(dist_struc(j)%param(dist_struc(j)%numparam))
       allocate(dist_struc(j)%param_name(dist_struc(j)%numparam))
!       read(ftn,'(a)',advance='no') (dist_struc(j)%param_name(p),p=1,dist_struc(j)%numparam)

!Allocate and Read in the statistical parameters (e.g., mean, sigma for normal) describing the distribution
       select case (dist_struc(j)%dist_ID)
       case (dist_triangularID)
          dist_struc(j)%numstatparam=3
       case (dist_betaID)
          dist_struc(j)%numstatparam=4
       case default
          dist_struc(j)%numstatparam=2
       end select
       allocate(dist_struc(j)%statparam(dist_struc(j)%numstatparam))
!       read(ftn,*) (dist_struc(j)%statparam(m),m=1,dist_struc(j)%numstatparam)

       read(line,*) dist_struc(j)%dist_ID, &
            (dist_struc(j)%param_name(p),p=1,dist_struc(j)%numparam), &
            (dist_struc(j)%statparam(m),m=1,dist_struc(j)%numstatparam)

       var_count=var_count+dist_struc(j)%numparam
    enddo


    call LIS_releaseUnitNumber(ftn)
    write(LIS_logunit,*) 'Finished reading prior probability distribution attributes file ..'

    ! Initialize objective function data structures
    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    objfuncField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Objective Function Value",rc=status)
    call LIS_verify(status)

    minField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Min Criteria Value",rc=status)
    call LIS_verify(status)

    maxField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Max Criteria Value",rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(objfuncField, localDE=0, farrayPtr=objfunc, rc=status)
    call LIS_verify(status)
    objfunc = 0.0

    call ESMF_StateAdd(LIS_ObjectiveFunc, (/objfuncField/), rc=status)
    call LIS_verify(status)  

    call ESMF_StateAdd(LIS_ObjectiveFunc, (/minField/), rc=status)
    call LIS_verify(status)  

    call ESMF_StateAdd(LIS_ObjectiveFunc, (/maxField/), rc=status)
    call LIS_verify(status)  

    ! Could add test here (or elsewhere) to warn if bounds set in dec space file are 
    ! truncating distributions set here.  Recall that the bounds set
    ! in the decision space files are passed along to the lsm (in use 
    ! to date); lsm, in turn, sets a flag (modflag) to true if a bound
    ! is violated.  Use of truncated distributions can be a useful feature 
    ! but we want to avoid  inadvertent use of truncated distributions.
    ! 
    !Test
   call ESMF_StateGet(LIS_decisionSpace, itemCount = nvars, rc=status)
   call LIS_verify(status)

   if (var_count.ne.nvars) then
      write(LIS_logunit,*) 'The number of variables covered by the prior (',var_count,') &
           not equal to number of decision space variable (',nvars,')'
      write(LIS_logunit,*) 'Program stopping....'
      call LIS_endrun
   endif

  end subroutine initializePObjFunc

!Sample from distribution
  subroutine dist_sample(distrib,sample)
    type(diststruc) :: distrib
    real  :: sample(distrib%numparam)
    select case (distrib%dist_ID)
    case (dist_uniformID)
       call sampleUniform(distrib%statparam(1), &
            distrib%statparam(2), sample(1))
    case (dist_normalID)
       call sampleNormal(distrib%statparam(1), &
            distrib%statparam(2),sample(1))
    case (dist_lognormalID)
       call sampleLognormal(distrib%statparam(1), &
            distrib%statparam(2),sample(1))
    case (dist_triangularID)
       call sampleTriangular(distrib%statparam(1), &
            distrib%statparam(2),distrib%statparam(3), sample(1))
    case (dist_betaID)
       call sampleBeta(distrib%statparam(1), &
            distrib%statparam(2),distrib%statparam(3), distrib%statparam(4), sample(1))
    end select
  end subroutine dist_sample
  
  subroutine sampleNormal(mean, stddev, sample)
    real :: mean
    real :: stddev
    real :: sample(1)
    real :: rand
    real :: randsum
    integer :: i
    randsum=0
    do i=1,12
       call ran3(1,rand)
       randsum = randsum + rand
    enddo
    sample(1) = (randsum-6)*stddev+mean
  end subroutine sampleNormal
  
  subroutine sampleLognormal(median, GSD, sample)
    real :: median
    real :: GSD
    real :: mean
    real :: stddev
    real :: temp
    real :: sample(1)
    real :: rand
    real :: randsum
    integer :: i
    mean = log(median)
    stddev = log(GSD)
    randsum=0
    do i=1,12
       call ran3(1,rand)
       randsum = randsum + rand
    enddo
    temp = (randsum-6)*stddev+mean
    sample(1) = exp(temp)
  end subroutine sampleLognormal
  
  subroutine sampleUniform(lb, ub, sample)
    real :: lb
    real :: ub
    real :: sample(1)
    real :: rand
    call ran3(1,rand)
    sample(1) = rand*(ub-lb)+lb
  end subroutine sampleUniform

  subroutine sampleTriangular(lb, ub, mode, sample)
    real :: lb
    real :: ub
    real :: mode
    real :: sample(1)
    real :: rand
    real :: F_mode
    real :: X
    call ran3(1,rand)
    F_mode=(mode-lb)/(ub-lb)
    if (rand<=F_mode) then
       X=lb+sqrt(rand*(ub-lb)*(mode-lb))
    else
       X=ub-sqrt((1-rand)*(ub-lb)*(ub-mode))
    endif
    sample(1) = X
  end subroutine sampleTriangular

 subroutine sampleBeta(alpha, beta, lb, ub, sample)
    real :: lb
    real :: ub
    real :: alpha
    real :: beta
    real :: sample(1)
    real :: rand1
    real :: rand2
    real :: X
    logical :: done
    real :: temp

    done = .false.
    do while (.not.done)
       call ran3(1,rand1)
       call ran3(1,rand2)
       temp=rand1**(1/alpha)+rand2**(1/beta)
       if (temp<=1.0) then
          done=.true.
       end if
    end do
    X=rand1**(1/alpha)/temp
    sample(1) = X*(ub-lb)+lb
  end subroutine sampleBeta

  !  Evaluate probability (density)
  subroutine dist_prob(distrib,X, prob)
    type(diststruc) :: distrib
    real, intent(in) :: X(distrib%numparam)
    real, intent(out) :: prob
    select case (distrib%dist_ID)
    case (dist_uniformID)
       call densityUniform(distrib%statparam(1), &
            distrib%statparam(2), X(1), prob)
    case (dist_NormalID)
       call densityNormal(distrib%statparam(1), &
            distrib%statparam(2),X(1), prob)
    case (dist_lognormalID)
       call densityLognormal(distrib%statparam(1), &
            distrib%statparam(2),X(1), prob)
    case (dist_triangularID)
       call densityTriangular(distrib%statparam(1), &
            distrib%statparam(2),distrib%statparam(3), X(1), prob)
    case (dist_betaID)
       call densityBeta(distrib%statparam(1), &
            distrib%statparam(2),distrib%statparam(3),distrib%statparam(4), X(1), prob)
    end select
  end subroutine dist_prob

  subroutine densityNormal(mean, stddev, X, prob)
    real, intent(in) :: mean
    real, intent(in) :: stddev
    real, intent(in) :: X
    real, intent(out) :: prob
    prob=(1/(stddev*sqrt(2.0*PI)))*exp(-0.5*((X-mean)**2/(stddev**2)))
  end subroutine densityNormal

  subroutine densityLognormal(median, GSD, X, prob)
    real :: median
    real :: GSD
    real  :: mean
    real  :: stddev
    real, intent(in) :: X
    real, intent(out) :: prob
    mean = log(median)
    stddev = log(GSD)
    prob=(1/(X*stddev*sqrt(2.0*PI)))*exp(-0.5*((log(X)-mean)**2/(stddev**2)))
  end subroutine densityLognormal
  
  subroutine densityUniform(lb, ub, X, prob)
    real, intent(in) :: lb
    real, intent(in) :: ub
    real, intent(in) :: X
    real, intent(out) :: prob
    if ((X<lb) .or. (X>ub)) then
       prob = 0
    else
       prob = 1/(ub-lb)
    endif
  end subroutine densityUniform

  subroutine densityTriangular(lb, ub, mode, X, prob)
    real, intent(in) :: lb
    real, intent(in) :: ub
    real, intent(in) :: mode
    real, intent(in) :: X
    real, intent(out) :: prob
    if ((X<lb) .or. (X>ub)) then
       prob = 0
    elseif (X<=mode) then
       prob = 2.0*(X-lb)/((ub-lb)*(mode-lb))
    else ! X>mode
       prob = 2.0*(ub-X)/((ub-lb)*(ub-mode))
    endif
  end subroutine densityTriangular

subroutine densityBeta(alpha, beta, lb, ub, X, prob)
    real :: lb
    real :: ub
    real :: alpha
    real :: beta
    real :: prob
    real :: X
    real :: X_tfm
    real :: lnprob
integer :: ifault
    real :: B
    X_tfm=(X-lb)/(ub-lb)
    B=exp(alogam(dble(alpha),ifault)+alogam(dble(beta),ifault)-alogam(dble(alpha+beta),ifault))
    !The divide by (ub-lb) is so that prob sums to unity
    prob=((1/B)*X_tfm**(alpha-1)*(1-X_tfm)**(beta-1))/(ub-lb)
  end subroutine densityBeta

!  Evaluate nat log of probability (density)
  subroutine dist_lnprob(distrib,X, lnprob)
    type(diststruc) :: distrib
    real, intent(in) :: X(distrib%numparam)
    real, intent(out) :: lnprob
    select case (distrib%dist_ID)
    case (dist_uniformID)
       call lndensityUniform(distrib%statparam(1), &
            distrib%statparam(2), X(1), lnprob)
    case (dist_normalID)
       call lndensityNormal(distrib%statparam(1), &
            distrib%statparam(2),X(1), lnprob)
    case (dist_lognormalID)
       call lndensityLognormal(distrib%statparam(1), &
            distrib%statparam(2),X(1), lnprob)
    case (dist_triangularID)
       call lndensityTriangular(distrib%statparam(1), &
            distrib%statparam(2),distrib%statparam(3),X(1), lnprob)
    case (dist_betaID)
       call lndensityBeta(distrib%statparam(1), &
            distrib%statparam(2),distrib%statparam(3),distrib%statparam(4),X(1), lnprob)
    end select
  end subroutine dist_lnprob

  subroutine lndensityNormal(mean, stddev, X, lnprob)
    real, intent(in) :: mean
    real, intent(in) :: stddev
    real, intent(in) :: X
    real, intent(out) :: lnprob
    lnprob= 0 - 0.5* log(2.0*PI) -log(stddev) + &
          (-1/2)*((X - mean)**2)/(stddev**2)
  end subroutine lndensityNormal
  
  subroutine lndensityLognormal(median, GSD, X, lnprob)
    real :: median
    real :: GSD
    real :: mean
    real :: stddev
    real, intent(in) :: X
    real, intent(out) :: lnprob
    mean = log(median)
    stddev = log(GSD)
    if (X<=0) then
       lnprob = mindensity
    else
       lnprob= 0 -log(X) -log(stddev) -0.5* log(2.0*PI)  &
            -0.5*((log(X)-mean)**2/(stddev**2))
    endif
  end subroutine lndensityLognormal

  subroutine lndensityUniform(lb, ub, X, lnprob)
    real, intent(in) :: lb
    real, intent(in) :: ub
    real, intent(in) :: X
    real, intent(out) :: lnprob

    if ((X<lb) .or. (X>ub)) then
       lnprob = mindensity
    else
       lnprob = 0 - log(ub-lb)
    endif
  end subroutine lndensityUniform

   subroutine lndensityTriangular(lb, ub, mode, X, lnprob)
    real, intent(in) :: lb
    real, intent(in) :: ub
    real, intent(in) :: mode
    real, intent(in) :: X
    real, intent(out) :: lnprob
    if ((X<lb) .or. (X>ub)) then
       lnprob = mindensity
    elseif (X<=mode) then
       lnprob = log(2.0)+log(X-lb)-log(ub-lb)-log(mode-lb)
    else ! X>mode
       lnprob = log(2.0)+log(ub-X)-log(ub-lb)-log(ub-mode)
    endif
  end subroutine lndensityTriangular

subroutine lndensityBeta(alpha, beta, lb, ub, X, prob)
    real :: lb
    real :: ub
    real :: alpha
    real :: beta
    real :: prob
    real :: X
    real :: X_tfm
    integer :: ifault
    real :: lnprob
    real :: B
    X_tfm=(X-lb)/(ub-lb)
    B=alogam(dble(alpha),ifault)+alogam(dble(beta),ifault)-alogam(dble(alpha+beta),ifault)
    lnprob=0-B+(alpha-1)*log(X_tfm)+(beta-1)*log(1-X_tfm)
    prob=((1/B)*X_tfm**(alpha-1)*(1-X_tfm)**(beta-1))/(ub-lb)
  end subroutine lndensityBeta

!  Retrieve lower fractiles
  subroutine dist_lowbound(distrib,lowbound)
    type(diststruc) :: distrib
    real, intent(out) :: lowbound
    select case (distrib%dist_ID)
    case (dist_uniformID)
       call lowboundUniform(distrib%statparam(1), &
            distrib%statparam(2), lowbound)
    case (dist_normalID)
       call lowboundNormal(distrib%statparam(1), &
            distrib%statparam(2),lowbound)
    case (dist_lognormalID)
       call lowboundLognormal(distrib%statparam(1), &
            distrib%statparam(2),lowbound)
    case (dist_triangularID)
       call lowboundTriangular(distrib%statparam(1), &
            distrib%statparam(2),distrib%statparam(3), lowbound)
    case (dist_betaID)
       call lowboundBeta(distrib%statparam(1), &
            distrib%statparam(2),distrib%statparam(3),distrib%statparam(4), lowbound)
    end select
  end subroutine dist_lowbound

  subroutine lowboundNormal(mean, stddev, lowbound)
    real, intent(in) :: mean
    real, intent(in) :: stddev
    real, intent(out) :: lowbound
    lowbound= mean - 3*stddev
  end subroutine lowboundNormal
  
  subroutine lowboundLognormal(median, GSD, lowbound)
    real :: median
    real :: GSD
    real :: mean
    real :: stddev
    real :: temp
    real, intent(out) :: lowbound
    mean = log(median)
    stddev = log(GSD)
    temp=mean - 3*stddev
    lowbound= exp(temp)
  end subroutine lowboundLognormal
  
  subroutine lowboundUniform(lb, ub, lowbound)
    real, intent(in) :: lb
    real, intent(in) :: ub
    real, intent(out) :: lowbound
    lowbound = lb
  end subroutine lowboundUniform

  subroutine lowboundTriangular(lb, ub, mode, lowbound)
    real, intent(in) :: lb
    real, intent(in) :: ub
    real, intent(in) :: mode
    real, intent(out) :: lowbound
    lowbound = lb
  end subroutine lowboundTriangular

  subroutine lowboundBeta(alpha, beta, lb, ub, lowbound)
    real, intent(in) :: lb
    real, intent(in) :: ub
    real, intent(in) :: alpha
    real, intent(in) :: beta
    real, intent(out) :: lowbound
    lowbound = lb
  end subroutine lowboundBeta

 !  Retrieve lower fractiles
  subroutine dist_upbound(distrib,upbound)
    type(diststruc) :: distrib
    real, intent(out) :: upbound
    select case (distrib%dist_ID)
    case (dist_uniformID)
       call upboundUniform(distrib%statparam(1), &
            distrib%statparam(2), upbound)
    case (dist_normalID)
       call upboundNormal(distrib%statparam(1), &
            distrib%statparam(2),upbound)
    case (dist_lognormalID)
       call upboundLognormal(distrib%statparam(1), &
            distrib%statparam(2),upbound)
    case (dist_triangularID)
       call upboundTriangular(distrib%statparam(1), &
            distrib%statparam(2),distrib%statparam(3), upbound)
    case (dist_betaID)
       call upboundBeta(distrib%statparam(1), &
            distrib%statparam(2),distrib%statparam(3),distrib%statparam(4), upbound)
    end select
  end subroutine dist_upbound

  subroutine upboundNormal(mean, stddev, upbound)
    real, intent(in) :: mean
    real, intent(in) :: stddev
    real, intent(out) :: upbound
    upbound= mean + 3*stddev
  end subroutine upboundNormal
  
  subroutine upboundLognormal(median, GSD, upbound)
    real :: median
    real :: GSD
    real :: mean
    real :: stddev
    real :: temp
    real, intent(out) :: upbound
    mean = log(median)
    stddev = log(GSD)
    temp=mean + 3*stddev
    upbound= exp(temp)
  end subroutine upboundLognormal
  
  subroutine upboundUniform(lb, ub, upbound)
    real, intent(in) :: lb
    real, intent(in) :: ub
    real, intent(out) :: upbound
    upbound = ub
  end subroutine upboundUniform

  subroutine upboundTriangular(lb, ub, mode, upbound)
    real, intent(in) :: lb
    real, intent(in) :: ub
    real, intent(in) :: mode
    real, intent(out) :: upbound
    upbound = ub
  end subroutine upboundTriangular

  subroutine upboundBeta(alpha, beta, lb, ub, upbound)
    real, intent(in) :: lb
    real, intent(in) :: ub
    real, intent(in) :: alpha
    real, intent(in) :: beta
    real, intent(out) :: upbound
    upbound = ub
  end subroutine upboundBeta

subroutine ran3(idum,rand)

!  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to
!  any negative value to initialize or reinitialize the sequence.
!  This function is taken from W.H. Press', "Numerical Recipes" p. 199.

    implicit none
    integer :: seed
    integer :: count, count_rate, count_max
    integer :: idum
    real   :: rand
      
    call system_clock(count,count_rate,count_max)
    seed = -count-count_max
    call random_seed(seed)
    call random_number(rand)
#if 0 
    real, parameter :: mbig=4000000.,mseed=1618033.,mz=0.
    real, parameter :: fac=1./4000000.
    
!  According to Knuth, any large mbig, and any smaller (but still large)
!  mseed can be substituted for the above values.
    integer         ::  ma(55)
    integer         ::  iff
    integer         ::  idum
    integer         ::  i,j,k,ii,jj
    real            ::  mk
    integer         ::  inext
    integer         ::  inextp
    real            ::  mj
    real            ::  rand
    
    iff = 0 
    ma = 0 
    
    if (idum.lt.0 .or. iff.eq.0) then
       iff=1
       mj=mseed-float(abs(idum))
       mj=mod(mj,mbig)      
       ma(55)=mj
       mk=1
       do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz) mk=mk+mbig
          mj=ma(ii)
       enddo
       do k=1,4
          do i=1,55
             ma(i)=ma(i)-ma(1+mod(i+30,55))
             if(ma(i).lt.mz) ma(i)=ma(i)+mbig
          enddo
       enddo
       inext=0
       inextp=31
       idum=1
    endif
    inext=inext+1
    if(inext.eq.56) inext=1
    inextp=inextp+1
    if(inextp.eq.56) inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.mz) mj=mj+mbig
    ma(inext)=mj
    rand=mj*fac
    return
#endif
  end subroutine ran3

function alogam ( x, ifault )

!*****************************************************************************80
!
!! ALOGAM computes the logarithm of the Gamma function.
!
!  Modified:
!
!    28 March 1999
!
!  Author:
!
!    Malcolm Pike, David Hill
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Malcolm Pike, David Hill,
!    Algorithm 291: 
!    Logarithm of Gamma Function,
!    Communications of the ACM,
!    Volume 9, Number 9, September 1966, page 684.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!    X should be greater than 0.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    Output, real ( kind = 8 ) ALOGAM, the logarithm of the Gamma 
!    function of X.
!
  implicit none

  real    ( kind = 8 ) alogam
  real    ( kind = 8 ) f
  integer ( kind = 4 ) ifault
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z

  if ( x <= 0.0D+00 ) then
    ifault = 1
    alogam = 0.0D+00
    return
  end if

  ifault = 0
  y = x

  if ( x < 7.0D+00 ) then

    f = 1.0D+00
    z = y

    do while ( z < 7.0D+00 )
      f = f * z
      z = z + 1.0D+00
    end do

    y = z
    f = - log ( f )

  else

    f = 0.0D+00

  end if

  z = 1.0D+00 / y / y
    
  alogam = f + ( y - 0.5D+00 ) * log ( y ) - y &
    + 0.918938533204673D+00 + &
    ((( &
    - 0.000595238095238D+00   * z &
    + 0.000793650793651D+00 ) * z &
    - 0.002777777777778D+00 ) * z &
    + 0.083333333333333D+00 ) / y

  return
end function alogam
end module PobjFunc_Mod
