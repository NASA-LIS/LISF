!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module DEMCAlgorithm
!BOP
!
! !MODULE: DEMCAlgorithm
!
! !DESCRIPTION: 
! This module contains routines that define the operations
! of the DEMC algorithm.  It is a parallel implementation of 
! the ter Braak (2006)
! algorithm
!  
! !REVISION HISTORY:
!  08 July 2010; Sujay Kumar, Ken Harrison; Initial Specification
!
!
! TODO: 
!  1. Initialization capability from a GA file
!  2. Modify the output so that it can be used to initialize the LSM
!     
  use ESMF
  use DEMC_varctl
  use LIS_coreMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_optUEMod
  use LIS_logMod

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: DEMC_init
  public :: DEMC_setup
  public :: DEMC_run
  public :: DEMC_checkConvergence
  public :: DEMC_getdecSpaceValues
  public :: DEMC_readrestart
  public :: DEMC_getNparam
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
  type demcstruc
     real,    allocatable :: pre_fitness(:)
     real,    allocatable :: temp_fitness(:)
     real,    allocatable :: curr_fitness(:)
     real,    allocatable :: pre_soln(:,:)
     real,    allocatable :: temp_soln(:,:)
     real,    allocatable :: curr_soln(:,:)
     integer, allocatable :: iparent1(:)
     integer, allocatable :: iparent2(:)
     real             :: acceptcount
     logical, allocatable :: SlotFilled(:)
     logical, allocatable :: Provisional(:)
     logical, allocatable :: CanBeParent(:)
  end type demcstruc

  type(demcctl)            :: demc_ctl
  type(demcstruc), allocatable :: demc_struc(:)

  contains
   
!BOP
! !ROUTINE: DEMC_init
! \label{DEMC_init}
! 
! !INTERFACE: 
    subroutine DEMC_init()
! !USES: 

! !DESCRIPTION: This routine performs the initialization steps for DEMC.  
!   It initializes the required memory structures, and 
!   creates an initial random population, for each grid point. 

!   
!
!EOP      
      implicit none 
      
      integer                :: status
      integer                :: ftn
      integer                :: i,t,j,m
      real                   :: rand
      type(ESMF_Field)       :: varField
      type(ESMF_ArraySpec)   :: arrspec1
      integer                :: iparent1, iparent2
      integer                :: dummy
      real,          allocatable :: var(:)

      demc_ctl%iterNo = 0 
      demc_ctl%zerothRun = .true. 
      demc_ctl%seed = -1000
      demc_ctl%minfitness = -1E20
      demc_ctl%overall_check = .false. 

      call ESMF_ConfigGetAttribute(LIS_config,demc_ctl%decspaceAttribsFile,&
           label="DEMC decision space attributes file:",rc=status)
      call LIS_verify(status, 'DEMC decision space attributes file: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,demc_ctl%restart,&
           label="DEMC start mode:",rc=status)
      call LIS_verify(status, 'DEMC start mode: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,demc_ctl%rfile,&
           label="DEMC restart file:",rc=status)

      call ESMF_ConfigGetAttribute(LIS_config,demc_ctl%maxiter,&
           label="DEMC number of iterations:",rc=status)
      call LIS_verify(status, 'DEMC number of generations: not defined')

! The parameter, 'pert_spread', is used to specify the sigma for all parameters 
! instead of specifying sigma for each parameter.  It does this by basing sigma, used in 'perturb_step',
! on this parameter and the width (ie, max-min) of the parameter range.

      call ESMF_ConfigGetAttribute(LIS_config,demc_ctl%pert_spread,&
           label="DEMC perturbation factor:",default=0.001,rc=status)
      call LIS_verify(status, 'DEMC perturbation factor: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,demc_ctl%modehopfreq,&
           label="DEMC mode hopping frequency:",default=0.1,rc=status)
      call LIS_verify(status, 'DEMC mode hopping frequency: not defined')

      call DEMC_setup()
      
   end subroutine DEMC_init

!BOP
! !ROUTINE: DEMC_setup
! \label{DEMC_setup}
!
! !INTERFACE: DEMC_setup
    subroutine DEMC_setup()
! !USES: 

! 
! !DESCRIPTION: 
!   This subroutine performs the second part of the DEMC initialization, 
!   The routine obtains the decision space object (from the models used
!   in optimization instance) and assigns the decision space variables
!   to the DEMC data structures
!  
!EOP
      implicit none
      
      type(ESMF_Field)           :: varField
      real,   pointer            :: vardata(:)
      integer                    :: t, k, n, m, j
      integer                    :: status

      n = 1

      call ESMF_StateGet(LIS_decisionSpace, itemCount=demc_ctl%nparam, &
           rc=status)
      call LIS_verify(status)

      allocate(demc_ctl%vname(demc_ctl%nparam))
      allocate(demc_ctl%parmax(demc_ctl%nparam))
      allocate(demc_ctl%parmin(demc_ctl%nparam))

      allocate(demc_struc(LIS_rc%ntiles(n)/LIS_rc%nensem(n)))

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
! population (for each tile point), each ensemble member represents a Markov chain. 

         allocate(demc_struc(t)%pre_soln(demc_ctl%nparam, LIS_rc%nensem(n)))
         allocate(demc_struc(t)%temp_soln(demc_ctl%nparam, LIS_rc%nensem(n)))
         allocate(demc_struc(t)%curr_soln(demc_ctl%nparam, LIS_rc%nensem(n)))
         allocate(demc_struc(t)%pre_fitness(LIS_rc%nensem(n)))
         allocate(demc_struc(t)%temp_fitness(LIS_rc%nensem(n)))
         allocate(demc_struc(t)%curr_fitness(LIS_rc%nensem(n)))

         allocate(demc_struc(t)%iparent1(LIS_rc%nensem(n)))
         allocate(demc_struc(t)%iparent2(LIS_rc%nensem(n)))

         allocate(demc_struc(t)%SlotFilled(LIS_rc%nensem(n))) !fitness is evaluated or not
         allocate(demc_struc(t)%Provisional(LIS_rc%nensem(n))) !fitness is evaluated or not
         allocate(demc_struc(t)%CanBeParent(LIS_rc%nensem(n))) !fitness is evaluated or not

         demc_struc(t)%acceptcount=0
      enddo

      call ESMF_StateGet(LIS_decisionSpace, itemNameList=demc_ctl%vname,&
           rc=status)
      call LIS_verify(status)

      do k=1,demc_ctl%nparam
         call ESMF_StateGet(LIS_decisionSpace, trim(demc_ctl%vname(k)), &
              varField, rc=status)
         call LIS_verify(status)
         
         call ESMF_AttributeGet(varField, 'MinRange',demc_ctl%parmin(k),&
              rc=status)
         call LIS_verify(status, 'setting minrange to decspace obj in DEMC')
         call ESMF_AttributeGet(varField, 'MaxRange',demc_ctl%parmax(k),&
              rc=status)
         call LIS_verify(status, 'setting maxrange to decspace obj in DEMC')

         call ESMF_FieldGet(varField,localDE=0, farrayPtr=vardata,rc=status)
         call LIS_verify(status)
         
         do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            do m=1,LIS_rc%nensem(n)
               demc_struc(t)%pre_soln(k,m) = vardata((t-1)*LIS_rc%nensem(n)+m)
            enddo
         enddo
      enddo      

    end subroutine DEMC_setup

!BOP
! 
! !ROUTINE: DEMC_checkConvergence
! \label{DEMC_checkConvergence}
! 
! !INTERFACE: 
    subroutine DEMC_checkConvergence(check)
! !ARGUMENTS: 
      logical, intent(INOUT) :: check
! 
! !DESCRIPTION: 
!  This routine checks to see if the convergence criteria for DEMC 
!  is reached. In this case, the routine simply checks to see if the
!  specified number of generations is reached. 
!EOP
      if(demc_ctl%iterNo.ge.demc_ctl%maxiter) then 
         check = .true.
      else
         check = .false.
      endif

    end subroutine DEMC_checkConvergence

!BOP
! !ROUTINE: DEMC_run
! \label{DEMC_run}
! 
! !INTERFACE: 
    subroutine DEMC_run()
! !USES: 

! !DESCRIPTION: 
!
!
!EOP  
      implicit none

      integer                  :: n

      n = 1
      write(LIS_logunit, *)  'DEMC iteration ',demc_ctl%iterNo 

      if(demc_ctl%iterNo.eq.0.or.demc_ctl%restart.eq.1) then 
         call DEMC_setPre()
         call DEMC_initNewIteration()
      else                                 !the crux of the DEMC processing
         call DEMC_evaluatePopulation()
         if(demc_ctl%overall_check .eqv. .true.) then 
            call DEMC_writeoutput()
            call DEMC_setCurrent()
            call DEMC_initNewIteration()
         endif
      endif
      call DEMC_perturbPopulation()
    end subroutine DEMC_run

!BOP
! 
! !ROUTINE: DEMC_setPre
! \label{DEMC_setPre}
!
! !INTERFACE: 
    subroutine  DEMC_setPre()
!
! !DESCRIPTION: 
!
!EOP
      implicit none

      integer            :: status
      type(ESMF_Field)   :: fitField
      type(ESMF_Field)   :: feasField
      integer,pointer    :: mod_flag(:)
      real, pointer      :: fitValue(:)
      integer            :: t,j,m,n

      n = 1
      call ESMF_StateGet(LIS_ObjectiveFunc,"Max Criteria Value",fitField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(fitField,localDE=0, farrayPtr=fitValue,rc=status)
      call LIS_verify(status)
      
      call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
      call LIS_verify(status)


      ! This iteration is setting the fitness values for the pre_solns
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            demc_struc(t)%pre_fitness(m) = fitValue((t-1)*LIS_rc%nensem(n)+m)
            !penalize infeasible solutions
            if(mod_flag((t-1)*LIS_rc%nensem(n)+m).eq.1) then 
               demc_struc(t)%pre_fitness(m) = demc_ctl%minfitness
            endif
         end do
      enddo

! Double check parameter bounds; Actually needed for the error variable (sigma) that is
! not currently thought of as a parameter of the LSM
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            do j=1,demc_ctl%nparam
               if((demc_struc(t)%pre_soln(j,m) .lt. demc_ctl%parmin(j))&
                    .or.(demc_struc(t)%pre_soln(j,m) .gt. demc_ctl%parmax(j))) then
                  demc_struc(t)%pre_fitness(m) = demc_ctl%minfitness
               endif
            enddo
         enddo
         
      enddo
    end subroutine DEMC_setPre

!BOP
! 
! !ROUTINE: DEMC_initNewIteration
! \label{DEMC_initNewIteration}
!
! !INTERFACE: 
    subroutine DEMC_initNewIteration()
!
! !DESCRIPTION: 
!
!EOP
      implicit none 

      integer                  :: t,n,m
      integer                  :: iparent1, iparent2
      real                     :: rand      

      n = 1

      demc_ctl%overall_check = .false.

      demc_ctl%iterNo = demc_ctl%iterNo + 1
      write(LIS_logunit, *)  'DEMC iteration ',demc_ctl%iterNo 

      !Initialize flags
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            demc_struc(t)%SlotFilled(m) = .false.
         enddo
      enddo
      
      !Find parent indices randomly
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            iparent1 = m
            
            !iparent1
            do while(iparent1.eq.m) 
               call ran3(1,rand)
               iparent1= int(rand*LIS_rc%nensem(n))+1
            enddo
            
            !iparent2
            iparent2 = m
            do while((iparent2.eq.m).or.(iparent2.eq.iparent1)) 
               call ran3(1,rand)
               iparent2 = int(rand*LIS_rc%nensem(n))+1
            enddo
            
            demc_struc(t)%iparent1(m) = iparent1
            demc_struc(t)%iparent2(m) = iparent2
         enddo
      enddo
    end subroutine DEMC_initNewIteration

!BOP
! 
! !ROUTINE: DEMC_perturbPopulation
! \label{DEMC_perturbPopulation}
!
! !INTERFACE: 
    subroutine DEMC_perturbPopulation()
!
! !DESCRIPTION: 
! 
!EOP      
      implicit none
      integer                  :: n, t, m
      logical                  :: parent1valid
      logical                  :: parent2valid

      n = 1

      !if(demc_ctl%overall_check.eqv..true) then
      !do while(.not.demc_ctl%overall_check)
      !Initialize flags
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            demc_struc(t)%Provisional(m) = .false.
            demc_struc(t)%CanBeParent(m) = .false.
         enddo
      enddo
      
      !CanBeParent is different than SlotFilled
      !In the 'change' while-loop, if Slot becomes filled
      !CanBeParent is set to true only if rejected
      !Provisionals are only accepted with this condition
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            if(demc_struc(t)%SlotFilled(m)) then
               demc_struc(t)%CanBeParent(m) = .true.
            endif
         enddo
      enddo
      
      !In 'change' while-loop below, need to know whether a slot is provisional or not
      !If provisional, there is the additional check to see if CanBeParent is true
      !for parents of provisional
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            if(.not.demc_struc(t)%SlotFilled(m)) then
               parent1valid= (demc_struc(t)%CanBeParent(demc_struc(t)%iparent1(m))) .or. &
                    (demc_struc(t)%iparent1(m).gt. m)
               parent2valid= (demc_struc(t)%CanBeParent(demc_struc(t)%iparent2(m))) .or. &
                    (demc_struc(t)%iparent2(m).gt. m)
               
               if (parent1valid .and. parent2valid) then
                  demc_struc(t)%Provisional(m) = .false.
               else
                  demc_struc(t)%Provisional(m) = .true.
               endif
            endif
         enddo
      enddo
      
      !Find proposed candidates for all unfilled slots,
      !For provisionals, this assumes that unfilled parents will be rejected
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            if (.not.demc_struc(t)%SlotFilled(m)) then
               call perturb_step(t,m)
            endif
         enddo
      enddo
      
    end subroutine DEMC_perturbPopulation
    
!BOP
! 
! !ROUTINE: DEMC_evaluatePopulation
! \label{DEMC_evaluatePopulation}
!
! !INTERFACE: 
    subroutine DEMC_evaluatePopulation()
!
! !DESCRIPTION: 
!
!EOP

      implicit none

      integer                  :: status
      type(ESMF_Field)         :: fitField
      type(ESMF_Field)         :: feasField
      integer,pointer          :: mod_flag(:)
      real, pointer            :: fitValue(:)
      logical                  :: parent1valid
      logical                  :: parent2valid
      logical                  :: accept
      logical                  :: change
      logical                  :: proceed
      integer                  :: t,m,n,j
      integer                  :: count
      logical                  :: temp_boolean

      n = 1

      call ESMF_StateGet(LIS_ObjectiveFunc,"Max Criteria Value",fitField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(fitField,localDE=0, farrayPtr=fitValue,rc=status)
      call LIS_verify(status)

      call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
      call LIS_verify(status)


      ! Set the fitness values for the solutions in the current LIS run (includes candidates and unchanged)
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            demc_struc(t)%temp_fitness(m) = fitValue((t-1)*LIS_rc%nensem(n)+m)
            !penalize infeasible solutions
            if(mod_flag((t-1)*LIS_rc%nensem(n)+m).eq.1) then 
               demc_struc(t)%temp_fitness(m) = demc_ctl%minfitness
            endif
         end do
      enddo

      ! Double check parameter bounds; Actually needed for the error variable (sigma) that is
      ! not currently thought of as a parameter of the LSM
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            do j=1,demc_ctl%nparam
               if((demc_struc(t)%temp_soln(j,m) .lt. demc_ctl%parmin(j)) &
                    .or.(demc_struc(t)%temp_soln(j,m) .gt. demc_ctl%parmax(j))) then
                  demc_struc(t)%temp_fitness(m) = demc_ctl%minfitness
               endif
            enddo
         enddo
      enddo


      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
      !'change' while loop
         change = .true.
         do while (change)
            change=.false.
            do m=1,LIS_rc%nensem(n)
               proceed=.false.
               if(.not.demc_struc(t)%SlotFilled(m)) then
                  if (demc_struc(t)%Provisional(m)) then
                     parent1valid= (demc_struc(t)%CanBeParent(demc_struc(t)%iparent1(m))) .or. &
                          (demc_struc(t)%iparent1(m).gt. m)
                     parent2valid= (demc_struc(t)%CanBeParent(demc_struc(t)%iparent2(m))) .or. &
                          (demc_struc(t)%iparent2(m).gt. m)
                     if(parent1valid.and.parent2valid) then
                        proceed=.true.
                     else
                        proceed=.false.
                     endif
                  else  !not provisional
                     proceed=.true.
                  endif
                  if (proceed) then
                     call accept_reject_test(t,m, accept)
                     if (accept) then 
                        demc_struc(t)%curr_soln(:,m)=demc_struc(t)%temp_soln(:,m)                     
                        demc_struc(t)%curr_fitness(m)=demc_struc(t)%temp_fitness(m)                     
                        demc_struc(t)%SlotFilled(m) = .true.
                        demc_struc(t)%acceptcount=demc_struc(t)%acceptcount+1
                     else
                        !curr_soln=prev_soln; but already was initialized as such
                        demc_struc(t)%SlotFilled(m) = .true.
                        demc_struc(t)%CanBeParent(m)=.true.
                     endif
                     change=.true.
                  else
                     !pass
                  endif
               else
                  !Skip; slot already filled
               endif
            enddo
         enddo
      enddo !ends change while loop

      ! Current iteration complete?
      temp_boolean = .true.  !will flip to false if any of the slots are unfilled 
      count = 0 
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            temp_boolean = temp_boolean.and.demc_struc(t)%SlotFilled(m)
            if(demc_struc(t)%SlotFilled(m)) then 
               count = count + 1
            endif
         enddo
      enddo
      !print*, 'Percent completed ', real(count*100.0/LIS_rc%nensem(n))
      !print*, 'Accept rate (%) ', real(demc_struc(t)%acceptcount*100.0/demc_ctl%iterNo*50)
      write(LIS_logunit,*) 'Percent completed ', real(count*100.0/LIS_rc%nensem(n))

      demc_ctl%overall_check = temp_boolean
    end subroutine DEMC_evaluatePopulation

!BOP
! 
! !ROUTINE: DEMC_setCurrent
! \label{DEMC_setCurrent}
!
! !INTERFACE: 
    subroutine DEMC_setCurrent
!
! !DESCRIPTION: 
!
!EOP   

      implicit none
      
      integer           :: n, t, m

      n = 1
      !1) If parent index r<i, must refer to current iteration; others (r>i) to pre_solution
      !2) Provisionals assume proposed candidate to be rejected
      !3) Hence, need to initialize curr_soln to pre_soln values
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         !do m=1,LIS_rc%nensem(n)
         demc_struc(t)%curr_soln=demc_struc(t)%pre_soln
         demc_struc(t)%curr_fitness=demc_struc(t)%pre_fitness
         !enddo
      enddo
      !endif
      
    end subroutine DEMC_setCurrent
    
!BOP
! 
! !ROUTINE: perturb_setup
! \label{perturb_setup}
!
! !INTERFACE: 
    subroutine perturb_step(t,i)      
!
! !DESCRIPTION: 
!
!EOP      
      implicit none

      integer :: n 
      integer :: t, i, j,k
      real    :: R1_j, R2_j
      integer :: iparent1, iparent2
      real    :: gamma
      real    :: epsilon, b_param
      real    :: rand

      n = 1

      !Set gamma = gamma*, the optimal for multivariate post. dist. as per ter Braak (2006)
      if (demc_ctl%modehopfreq.ne.0) then
         if (mod(demc_ctl%iterNo,int(1.0/demc_ctl%modehopfreq)).eq.0) then
            gamma = 0.98  ! see ter Braak (2006) for explanation why preferable to 1.0; 
         else
            gamma = 2.38/sqrt(2.0*dble(demc_ctl%nparam))
         endif
      else  !is equal to zero
         gamma = 2.38/sqrt(2.0*dble(demc_ctl%nparam))
      endif

      do j=1,demc_ctl%nparam

         iparent1 = demc_struc(t)%iparent1(i)
         iparent2 = demc_struc(t)%iparent2(i)


         b_param = (demc_ctl%parmax(j) - demc_ctl%parmin(j))&
              *demc_ctl%pert_spread
         call ran3(1,rand)
         epsilon = ((rand-0.5)*2)*b_param

         if (iparent1>i) then
            R1_j=demc_struc(t)%pre_soln(j,iparent1)
         else
            R1_j=demc_struc(t)%curr_soln(j,iparent1)
         endif

         if (iparent2>i) then
            R2_j=demc_struc(t)%pre_soln(j,iparent2)
         else
            R2_j=demc_struc(t)%curr_soln(j,iparent2)
         endif

         demc_struc(t)%temp_soln(j,i) = demc_struc(t)%pre_soln(j,i) + &
              gamma*(R2_j - R1_j) +  epsilon                  
      enddo
    end subroutine perturb_step

!BOP
! 
! !ROUTINE: accept_reject_test
! \label{accept_reject_test}
!
! !INTERFACE:       
    subroutine accept_reject_test(t,i,accept)
!
! !DESCRIPTION: 
!
!EOP
      implicit none

      integer :: n 
      integer :: t, i
      logical, intent(out) :: accept
      real    :: rand, rand_t
      real    :: ln_p_ratio
      real    :: p_temp_soln
      real    :: p_pre_soln

      n = 1
      
!compare fitnesses of pre_soln and candidates. 
      call ran3(1,rand)
      
      p_temp_soln = demc_struc(t)%temp_fitness(i)
      p_pre_soln = demc_struc(t)%pre_fitness(i)

!assume that ll values in log form
      ln_p_ratio = p_temp_soln - p_pre_soln
      rand_t = log(rand)
      
      if(rand_t.lt.ln_p_ratio) then 
         accept = .true.
      else
         accept = .false.
      endif
      
  end subroutine accept_reject_test

!BOP
! 
! !ROUTINE: DEMC_writeOutput
! \label{DEMC_writeOutput}
!
! !INTERFACE: 
  subroutine DEMC_writeOutput()
!
! !DESCRIPTION: 
! 
!EOP

!svk: do we need a seperate output file?    
    call writeDEMCrestart()
!    call writeDEMCoutput()

  end subroutine DEMC_writeOutput

!BOP
! 
! !ROUTINE: writeDEMCrestart
! \label{writeDEMCrestart}
! 
! !INTERFACE: 
  subroutine writeDEMCrestart
! !USES: 
    use LIS_fileIOMod,  only : LIS_create_output_directory
    use LIS_historyMod, only : LIS_writevar_restart
! 
! !DESCRIPTION: 
! 
! This routine writes the checkpoint data for a DEMC restart
! 
!EOP

    integer             :: n 
    integer             :: i,m,t
    integer             :: status
    character(len=LIS_CONST_PATH_LEN) :: filen
    character (len=4)   :: fgen
    character*100       :: vnames(demc_ctl%nparam)
    real, allocatable       :: vardata(:)

    n = 1

    allocate(vardata(LIS_rc%ntiles(n)))

    if(LIS_masterproc) then 
       call LIS_create_output_directory('DEMC')
       write(unit=fgen, fmt='(i4.4)') demc_ctl%iterNo
       filen = trim(LIS_rc%odir)//'/DEMC/DEMC.'&
            //trim(fgen)//'.DEMCrst'
       open(40,file=filen,status='unknown',form='unformatted')
       write(40) demc_ctl%iterNo
    endif

    do i=1,demc_ctl%nparam
       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
          do m=1,LIS_rc%nensem(n)
             vardata((t-1)*LIS_rc%nensem(n)+m) = demc_struc(t)%pre_soln(i,m)
          enddo
       enddo
       if(LIS_masterproc) write(40) demc_ctl%vname(i)
       call LIS_writevar_restart(40,n,vardata)
    enddo

    do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
       do m=1,LIS_rc%nensem(n)
          vardata((t-1)*LIS_rc%nensem(n)+m) = demc_struc(t)%pre_fitness(m)
       enddo
    enddo
    call LIS_writevar_restart(40,n,vardata)
       
    if(LIS_masterproc) then 
       close(40)
       write(LIS_logunit,*) 'DEMC checkpoint file written ',trim(filen)
    endif
    deallocate(vardata)

  end subroutine writeDEMCrestart


!BOP
! 
! !ROUTINE: DEMC_readrestart
! \label{DEMC_readrestart}
! 
! !INTERFACE: 
  subroutine DEMC_readrestart
! !USES: 
    use LIS_historyMod,      only : LIS_readvar_restart
! 
! !DESCRIPTION: 
! 
!   This routine reads the checkpoint data for a DEMC restart
!EOP
    integer             :: n 
    integer             :: i, t, m
    integer             :: status
    real, pointer       :: vardata(:) 
    character*100       :: vname_in_file
    character*100       :: vnames(demc_ctl%nparam)
    type(ESMF_Field)    :: varField(demc_ctl%nparam)

    if(demc_ctl%restart.eq.1) then 
       n = 1    
       allocate(vardata(LIS_rc%ntiles(n)))
       
       call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
            rc=status)
       call LIS_verify(status)
       
       write(LIS_logunit,*) 'Reading the DEMC restart file ..', trim(demc_ctl%rfile)
       
       open(40,file=demc_ctl%rfile,form='unformatted')        
       
       read(40) demc_ctl%iterNo
       write(LIS_logunit,*) 'Iteration Number ',demc_ctl%iterNo
       
       do i=1,demc_ctl%nparam
          
          read(40) vname_in_file

!!!! should be picked up by set decision space routine
!!!!          call ESMF_StateGet(LIS_decisionSpace, trim(vname_in_file), &
!!!!               varField(i), rc=status)
!!!!          call LIS_verify(status)
!!!!          
!!!!          call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata1,rc=status)
!!!!          call LIS_verify(status)
          
          call LIS_readvar_restart(40,n,vardata)       
          
!!!!          vardata1 = vardata 
          
          do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
             do m=1,LIS_rc%nensem(n)
                demc_struc(t)%pre_soln(i,m) = vardata((t-1)*LIS_rc%nensem(n)+m) 
             enddo
          enddo
       enddo
!fitness
       call LIS_readvar_restart(40,n,vardata)            
       
       close(40)
       write(LIS_logunit,*) 'Finished reading the DEMC restart file ..'
       
! update the DEMC objects. 

       deallocate(vardata)
    endif

  end subroutine DEMC_readrestart

!BOP
! 
! !ROUTINE: DEMC_getNparam
! \label{DEMC_getNparam}
! 
! !INTERFACE: 
    subroutine DEMC_getNparam(nparam)
! 
! !DESCRIPTION: 
!  This method returns the number of decision space variables 
!
!EOP      
      integer   :: nparam
      
      nparam = demc_ctl%nparam

    end subroutine DEMC_getNparam

!BOP
! !ROUTINE: DEMC_getdecSpaceValues
! \label{DEMC_getdecSpaceValues}
! 
! !INTERFACE: 
    subroutine DEMC_getdecSpaceValues(n)
! !USES: 

! 
! !DESCRIPTION: 
!  This routine returns the array of decision space variables from 
!  the GA data structures
!EOP
      implicit none

      integer            :: n
      type(ESMF_Field)   :: varField(demc_ctl%nparam)
      real, pointer      :: vardata(:)
      integer            :: i, t, m
      integer            :: status
      
      do i=1,demc_ctl%nparam
         call ESMF_StateGet(LIS_decisionSpace, trim(demc_ctl%vname(i)), &
              varField(i), rc=status)
         call LIS_verify(status)
         
         call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata,rc=status)
         call LIS_verify(status)
         
         if(demc_ctl%restart.eq.1) then 
            do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
               do m=1,LIS_rc%nensem(n)
                  vardata((t-1)*LIS_rc%nensem(n)+m)= demc_struc(t)%pre_soln(i,m)
               enddo
            enddo
         else         
            if(demc_ctl%zerothRun) then 
               
               do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                  do m=1,LIS_rc%nensem(n)
                     vardata((t-1)*LIS_rc%nensem(n)+m) = demc_struc(t)%pre_soln(i,m)
                  enddo
               enddo
            else
               do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                  do m=1,LIS_rc%nensem(n)
                     vardata((t-1)*LIS_rc%nensem(n)+m) = demc_struc(t)%temp_soln(i,m)
                  enddo
               enddo
            endif
         endif
      enddo

      if(demc_ctl%restart.eq.1) then 
         demc_ctl%restart = 2
         demc_ctl%zerothRun = .false. 
      endif

    end subroutine DEMC_getdecSpaceValues


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

end module DEMCAlgorithm
