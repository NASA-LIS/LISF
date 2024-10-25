!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module RWMCMCAlgorithm
!BOP
!
! !MODULE: RWMCMCAlgorithm
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!  04 Feb 2008; Sujay Kumar; Initial Specification
!
  use ESMF
  use RWMCMC_varctl
  use LIS_coreMod
  use LIS_optUEMod
  use LIS_logMod

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: RWMCMC_init
  public :: RWMCMC_setup
  public :: RWMCMC_run
  public :: RWMCMC_checkConvergence
  public :: RWMCMC_getdecSpaceValues
  public :: RWMCMC_readrestart
  public :: RWMCMC_getNparam
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
  type mcmcstruc
     real,    allocatable :: pre_fitness(:)
     real,    allocatable :: curr_fitness(:)
     real,    allocatable :: pre_soln(:,:)
     real,    allocatable :: candidate(:,:)
     integer          :: acceptcount
  end type mcmcstruc

  type(mcmcctl)            :: mcmc_ctl
  type(mcmcstruc), allocatable :: mcmc_struc(:)

  contains
   
!BOP
! !ROUTINE: RWMCMC_init
! \label{RWMCMC_init}
! 
! !INTERFACE: 
    subroutine RWMCMC_init()
! !USES: 
!!$
!!$
! !DESCRIPTION: 
!   
!
!EOP      
      implicit none 
      
      integer                :: n 
      integer                :: status
      integer                :: ftn
      integer                :: i,t,j,m
      real                   :: rand



      mcmc_ctl%iterNo = 0 
      mcmc_ctl%zerothRun = .true. 
      mcmc_ctl%seed = -1000
      mcmc_ctl%minFitness = -1E20

      call ESMF_ConfigGetAttribute(LIS_config,mcmc_ctl%decspaceAttribsFile,&
           label="RWMCMC decision space attributes file:",rc=status)
      call LIS_verify(status, 'RWMCMC decision space attributes file: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,mcmc_ctl%restart,&
           label="RWMCMC start mode:",rc=status)
      call LIS_verify(status, 'RWMCMC start mode: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,mcmc_ctl%rfile,&
           label="RWMCMC restart file:",rc=status)
      call LIS_verify(status, 'RWMCMC restart file: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,mcmc_ctl%maxiter,&
           label="RWMCMC number of iterations:",rc=status)
      call LIS_verify(status, 'RWMCMC number of generations: not defined')

! sigma (used in the perturbations) is based on this value and the 
! width of the parameter range. This single parameter is used to specify 
! the sigma for all parameters instead of specifying sigma for each 
! parameter. 

      call ESMF_ConfigGetAttribute(LIS_config,mcmc_ctl%pert_spread,&
           label="RWMCMC perturbation factor:",default=0.10,rc=status)
      call LIS_verify(status, 'RWMCMC perturbation factor: not defined')

      call RWMCMC_setup()

   end subroutine RWMCMC_init

!BOP
! !ROUTINE: RWMCMC_setup
! \label{RWMCMC_setup}
!
! !INTERFACE: RWMCMC_setup
    subroutine RWMCMC_setup()
! !USES: 
! 
! !DESCRIPTION: 
!   This subroutine performs the second part of the GA initialization, 
!   This routine is currently empty. 
!  
!EOP
      implicit none
      
      type(ESMF_Field)       :: varField
      real,          pointer :: vardata(:)
      integer                    :: t, k, n, m, j
      integer                    :: status

      n = 1

      call ESMF_StateGet(LIS_decisionSpace, itemCount=mcmc_ctl%nparam, &
           rc=status)
      call LIS_verify(status)

      allocate(mcmc_ctl%vname(mcmc_ctl%nparam))
      allocate(mcmc_ctl%parmax(mcmc_ctl%nparam))
      allocate(mcmc_ctl%parmin(mcmc_ctl%nparam))

      allocate(mcmc_struc(LIS_rc%ntiles(n)/LIS_rc%nensem(n)))

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
! population (for each tile point), each ensemble member represents a Markov chain. 

         allocate(mcmc_struc(t)%pre_soln(mcmc_ctl%nparam, LIS_rc%nensem(n)))
         allocate(mcmc_struc(t)%candidate(mcmc_ctl%nparam, LIS_rc%nensem(n)))
         allocate(mcmc_struc(t)%pre_fitness(LIS_rc%nensem(n)))
         allocate(mcmc_struc(t)%curr_fitness(LIS_rc%nensem(n)))
      enddo


      call ESMF_StateGet(LIS_decisionSpace, itemNameList=mcmc_ctl%vname,&
           rc=status)
      call LIS_verify(status)

      do k=1,mcmc_ctl%nparam
         call ESMF_StateGet(LIS_decisionSpace, trim(mcmc_ctl%vname(k)), &
              varField, rc=status)
         call LIS_verify(status)
         
         call ESMF_AttributeGet(varField, 'MinRange',mcmc_ctl%parmin(k),&
              rc=status)
         call LIS_verify(status, 'setting minrange to decspace obj in RWMCMC')
         call ESMF_AttributeGet(varField, 'MaxRange',mcmc_ctl%parmax(k),&
              rc=status)
         call LIS_verify(status, 'setting maxrange to decspace obj in RWMCMC')

         call ESMF_FieldGet(varField,localDE=0, farrayPtr=vardata,rc=status)
         call LIS_verify(status)
         
         do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            do m=1,LIS_rc%nensem(n)
               mcmc_struc(t)%pre_soln(k,m) = vardata((t-1)*LIS_rc%nensem(n)+m)
            enddo
         enddo
      enddo      

    end subroutine RWMCMC_setup

!BOP
! 
! !ROUTINE: RWMCMC_checkConvergence
! \label{RWMCMC_checkConvergence}
! 
! !INTERFACE: 
    subroutine RWMCMC_checkConvergence(check)
! !ARGUMENTS: 
      logical, intent(INOUT) :: check
! 
! !DESCRIPTION: 
!  This routine checks to see if the convergence criteria for GA 
!  is reached. In this case, the routine simply checks to see if the
!  specified number of generations is reached. 
!EOP
      if(mcmc_ctl%iterNo.ge.mcmc_ctl%maxiter) then 
         check = .true.
      else
         check = .false.
      endif

    end subroutine RWMCMC_checkConvergence

!BOP
! !ROUTINE: RWMCMC_run
! \label{RWMCMC_run}
! 
! !INTERFACE: 
    subroutine RWMCMC_run()
! !USES: 
! !DESCRIPTION: 
!   This routine performs the run steps for the GA, which includes the 
!   following steps
!   \begin{enumerate}
!   \item[Evaluate Fitness] : Evaluates the fitness of potential solutions  
!    \item[Selection] : Selects better solutions for "mating" 
!    \item[Crossover] : Performs recombination between parent solutions to 
!                      produce new potential solutions
!    \item[Mutation]  : Performs the mutation step 
!   \end{enumerate}
!EOP  
      implicit none

      integer                  :: n
!!$      integer                  :: t, m, j
!!$      integer                  :: status
!!$      type(ESMF_Field)         :: fitField
!!$      type(ESMF_Field)         :: feasField
!!$      real,   allocatable          :: fitValue(:)
!!$      integer, allocatable         :: mod_flag(:)

      n = 1
      write(LIS_logunit, *)  'RWMCMC iteration ',mcmc_ctl%iterNo 

      if(mcmc_ctl%iterNo.eq.0) then 
         call RWMCMC_setPre()
      else                                 !the crux of the RWMCMC processing
         call RWMCMC_evaluatePopulation()
      endif
      call RWMCMC_initNewIteration()
      call RWMCMC_perturbPopulation()

    end subroutine RWMCMC_run

    subroutine  RWMCMC_setPre()

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
            mcmc_struc(t)%pre_fitness(m) = fitValue((t-1)*LIS_rc%nensem(n)+m)
            !penalize infeasible solutions
            if(mod_flag((t-1)*LIS_rc%nensem(n)+m).eq.1) then 
               mcmc_struc(t)%pre_fitness(m) = mcmc_ctl%minfitness
            endif
         end do
      enddo

! Double check parameter bounds; Actually needed for the error variable (sigma) that is
! not currently thought of as a parameter of the LSM
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            do j=1,mcmc_ctl%nparam
               if((mcmc_struc(t)%pre_soln(j,m) .lt. mcmc_ctl%parmin(j))&
                    .or.(mcmc_struc(t)%pre_soln(j,m) .gt. mcmc_ctl%parmax(j))) then
                  mcmc_struc(t)%pre_fitness(m) = mcmc_ctl%minfitness
               endif
            enddo
         enddo
         
      enddo
    end subroutine RWMCMC_setPre

    subroutine RWMCMC_initNewIteration()

      implicit none 

      integer                  :: t,n,m
      integer                  :: iparent1, iparent2
      real                     :: rand      

      n = 1

      mcmc_ctl%iterNo = mcmc_ctl%iterNo + 1
      write(LIS_logunit, *)  'RWMCMC iteration ',mcmc_ctl%iterNo 

    end subroutine RWMCMC_initNewIteration

    subroutine RWMCMC_evaluatePopulation()


      implicit none

      integer                  :: status
      type(ESMF_Field)         :: fitField
      type(ESMF_Field)         :: feasField
      integer,pointer          :: mod_flag(:)
      real, pointer            :: fitValue(:)
      logical                  :: accept
      integer                  :: t,m,n,j
      integer                  :: count

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
            mcmc_struc(t)%curr_fitness(m) = fitValue((t-1)*LIS_rc%nensem(n)+m)
            !penalize infeasible solutions
            if(mod_flag((t-1)*LIS_rc%nensem(n)+m).eq.1) then 
               mcmc_struc(t)%curr_fitness(m) = mcmc_ctl%minfitness
            endif
         end do
      enddo

      ! Double check parameter bounds; Actually needed for the error variable (sigma) that is
      ! not currently thought of as a parameter of the LSM
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            do j=1,mcmc_ctl%nparam
               if((mcmc_struc(t)%candidate(j,m) .lt. mcmc_ctl%parmin(j)) &
                    .or.(mcmc_struc(t)%candidate(j,m) .gt. mcmc_ctl%parmax(j))) then
                  mcmc_struc(t)%curr_fitness(m) = mcmc_ctl%minfitness
               endif
            enddo
         enddo
      enddo

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            do m=1,LIS_rc%nensem(n)
               call accept_reject_test(t,m, accept)
               if (accept) then 
                  mcmc_struc(t)%pre_soln(:,m)=mcmc_struc(t)%candidate(:,m)                     
                  mcmc_struc(t)%pre_fitness(m)=mcmc_struc(t)%curr_fitness(m)                     
                  mcmc_struc(t)%acceptcount=mcmc_struc(t)%acceptcount+1
               else
                  !pre_soln=pre_soln; but already was initialized as such
               endif
            enddo
         enddo

    end subroutine RWMCMC_evaluatePopulation

    subroutine RWMCMC_perturbPopulation()      

      integer :: n 
      integer :: t, i, j,k
      real    :: rand, rand_sum, sigma

      n = 1

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)

         do i=1,LIS_rc%nensem(n)
            do j=1,mcmc_ctl%nparam

               rand_sum = 0 
               do k=1,12
                  call ran3(1,rand)
                  rand_sum = rand_sum + rand
               enddo

               sigma = (mcmc_ctl%parmax(j) - mcmc_ctl%parmin(j))*mcmc_ctl%pert_spread

               mcmc_struc(t)%candidate(j,i) = mcmc_struc(t)%pre_soln(j,i) + &
                    ((rand_sum - 6)*sigma)

            enddo
         enddo
              
      enddo
      
    end subroutine RWMCMC_perturbPopulation


    subroutine accept_reject_test(t,i,accept)

      integer :: n 
      integer :: t, i
      logical, intent(out) :: accept
      real    :: rand, rand_t
      real    :: ln_p_ratio
      real    :: p_candidate
      real    :: p_pre_soln

      n = 1
      
!compare fitnesses of pre_soln and candidates. 
      call ran3(1,rand)
            
      p_candidate = mcmc_struc(t)%curr_fitness(i)
      p_pre_soln = mcmc_struc(t)%pre_fitness(i)

      ln_p_ratio = p_candidate - p_pre_soln
      rand_t = log(rand)
      if(rand_t.lt.ln_p_ratio) then 
         accept = .true.
!         mcmc_struc(t)%pre_soln(:,i) = mcmc_struc(t)%candidate(:,i)
!         mcmc_struc(t)%pre_fitness(i) = mcmc_struc(t)%curr_fitness(i)
      else
         accept = .false.
       endif

    end subroutine accept_reject_test


  subroutine writeRWMCMCdata()

!svk: do we need a seperate output file?    
    call writeRWMCMCrestart()
!    call writeRWMCMCoutput()

  end subroutine writeRWMCMCdata

!BOP
! 
! !ROUTINE: writeRWMCMCrestart
! \label{RWMCMC_writeRWMCMCrestart}
! 
! !INTERFACE: 
  subroutine writeRWMCMCrestart
! !USES: 
    use LIS_fileIOMod,  only : LIS_create_output_directory
    use LIS_historyMod, only : LIS_writevar_restart
! 
! !DESCRIPTION: 
! 
! This routine writes the checkpoint data for a GA restart
! 
!EOP
!svk: This routine only writes the "candidate" soln. We need both pre_soln & candidate?
    integer             :: n 
    integer             :: i,m,t
    integer             :: status
    character*100       :: filen
    character (len=4)   :: fgen
    character*100       :: vnames(mcmc_ctl%nparam)
    real, allocatable       :: vardata(:)

    n = 1

    allocate(vardata(LIS_rc%ntiles(n)))

    if(LIS_masterproc) then 
       call LIS_create_output_directory('RWMCMC')
       write(unit=fgen, fmt='(i4.4)') mcmc_ctl%iterNo
       filen = trim(LIS_rc%odir)//'/RWMCMC/RWMCMC.'&
            //trim(fgen)//'.RWMCMCrst'
       open(40,file=filen,status='unknown',form='unformatted')
       write(40) mcmc_ctl%iterNo
    endif

    do i=1,mcmc_ctl%nparam
       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
          do m=1,LIS_rc%nensem(n)
             vardata((t-1)*LIS_rc%nensem(n)+m) = mcmc_struc(t)%pre_soln(i,m)
          enddo
       enddo
       if(LIS_masterproc) write(40) mcmc_ctl%vname(i)
       call LIS_writevar_restart(40,n,vardata)
    enddo

    do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
       do m=1,LIS_rc%nensem(n)
          vardata((t-1)*LIS_rc%nensem(n)+m) = mcmc_struc(t)%pre_fitness(m)
       enddo
    enddo
    call LIS_writevar_restart(40,n,vardata)

       
    if(LIS_masterproc) then 
       close(40)
       write(LIS_logunit,*) 'RWMCMC checkpoint file written ',trim(filen)
    endif
    deallocate(vardata)

  end subroutine writeRWMCMCrestart


!BOP
! 
! !ROUTINE: readRWMCMCrestart
! \label{readRWMCMCrestart}
! 
! !INTERFACE: 
  subroutine RWMCMC_readrestart
! !USES: 
    use LIS_historyMod,      only : LIS_readvar_restart
! 
! !DESCRIPTION: 
! 
!   This routine reads the checkpoint data for a RWMCMC restart
!EOP
    integer             :: n 
    integer             :: i, t, m
    integer             :: status
    real, allocatable       :: vardata(:) 
    real, pointer       :: vardata1(:)
    character*100       :: vnames(mcmc_ctl%nparam)
    type(ESMF_Field)    :: varField(mcmc_ctl%nparam)
    character*100       :: vname_in_file

    n = 1    
    allocate(vardata(LIS_rc%ntiles(n)))

    call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
         rc=status)
    call LIS_verify(status)
    
    write(LIS_logunit,*) 'Reading the RWMCMC restart file ..'

    open(40,file=mcmc_ctl%rfile,form='unformatted')        
    
    read(40) mcmc_ctl%iterNo
    write(LIS_logunit,*) 'Iteration Number ',mcmc_ctl%iterNo

    do i=1,mcmc_ctl%nparam

       read(40) vname_in_file

       call ESMF_StateGet(LIS_decisionSpace, trim(vname_in_file), &
            varField(i), rc=status)
       call LIS_verify(status)
       
       call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata1,rc=status)
       call LIS_verify(status)
       
       call LIS_readvar_restart(40,n,vardata)       
      
       vardata1 = vardata 

       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
          do m=1,LIS_rc%nensem(n)
              mcmc_struc(t)%pre_soln(i,m) = vardata((t-1)*LIS_rc%nensem(n)+m) 
           enddo
        enddo
     enddo
!fitness     
    call LIS_readvar_restart(40,n,vardata)       

    close(40)
    write(LIS_logunit,*) 'Finished reading the RWMCMC restart file ..'

! update the RWMCMC objects. 

    deallocate(vardata)

  end subroutine RWMCMC_readrestart

!BOP
! 
! !ROUTINE: RWMCMC_getNparam
! \label{RWMCMC_getNparam}
! 
! !INTERFACE: 
    subroutine RWMCMC_getNparam(nparam)
! 
! !DESCRIPTION: 
!  This method returns the number of decision space variables 
!
!EOP      
      integer   :: nparam
      
      nparam = mcmc_ctl%nparam

    end subroutine RWMCMC_getNparam

!BOP
! 
! !ROUTINE: RWMCMC_evaluateFitness
!  \label{RWMCMC_evaluateFitness}
! 
! !INTERFACE: 
    subroutine RWMCMC_evaluateFitness()
! !USES:   
! 
! !DESCRIPTION: 
!   This method computes the fitness values for each organism, 
!   selects the best solution and ensures its selection in the 
!   subsequent generation (if elitism is enabled). 
! 
!EOP
    integer            :: n 
    type(ESMF_Field)   :: fitField
    type(ESMF_Field)   :: feasField
    integer            :: status
    real, pointer      :: fitValue(:)
    integer, pointer   :: mod_flag(:)
    integer            :: t
    integer            :: m
    logical            :: best_status

    call evaluateobjfunction(LIS_rc%optuetype)

    call ESMF_StateGet(LIS_ObjectiveFunc,"Max Criteria Value",fitField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(fitField,localDE=0, farrayPtr=fitValue,rc=status)
    call LIS_verify(status)

    call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
    call LIS_verify(status)

  end subroutine RWMCMC_evaluateFitness

!BOP
! !ROUTINE: RWMCMC_getdecSpaceValues
! \label{RWMCMC_getdecSpaceValues}
! 
! !INTERFACE: 
    subroutine RWMCMC_getdecSpaceValues(n)
! !USES: 
! 
! !DESCRIPTION: 
!  This routine returns the array of decision space variables from 
!  the GA data structures
!EOP
      implicit none

      integer            :: n
      type(ESMF_Field)   :: varField(mcmc_ctl%nparam)
      real, pointer      :: vardata(:)
      integer            :: i, t, m
      integer            :: status
      
      do i=1,mcmc_ctl%nparam
         call ESMF_StateGet(LIS_decisionSpace, trim(mcmc_ctl%vname(i)), &
              varField(i), rc=status)
         call LIS_verify(status)
         
         call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata,rc=status)
         call LIS_verify(status)
         
         if(mcmc_ctl%restart.eq.1) then 
            do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
               do m=1,LIS_rc%nensem(n)
                  vardata((t-1)*LIS_rc%nensem(n)+m)= mcmc_struc(t)%pre_soln(i,m)
               enddo
            enddo
         else         
            if(mcmc_ctl%zerothRun) then 
               
               do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                  do m=1,LIS_rc%nensem(n)
                     vardata((t-1)*LIS_rc%nensem(n)+m) = mcmc_struc(t)%pre_soln(i,m)
                  enddo
               enddo
            else
               do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                  do m=1,LIS_rc%nensem(n)
                     vardata((t-1)*LIS_rc%nensem(n)+m) = mcmc_struc(t)%candidate(i,m)
                  enddo
               enddo
            endif
         endif
      enddo

      if(mcmc_ctl%restart.eq.1) then 
         mcmc_ctl%restart = 2
         mcmc_ctl%zerothRun = .false. 
      endif
    end subroutine RWMCMC_getdecSpaceValues

!!$!BOP
!!$! 
!!$! !ROUTINE: RWMCMC_setdecSpaceValues
!!$! 
!!$! !INTERFACE: 
!!$    subroutine RWMCMC_setdecSpaceValues(n, decvals)
!!$! !USES: 
!!$! 
!!$! !DESCRIPTION: 
!!$!   This routine sets the RWMCMC decision space variables based on the 
!!$!   input array. 
!!$!EOP      
!!$      implicit none
!!$
!!$      integer            :: n
!!$      real               :: decvals(mcmc_ctl%nparam, LIS_rc%ntiles(n))
!!$
!!$      integer            :: i, t, k, m
!!$      
!!$      do i=1, mcmc_ctl%nparam
!!$         do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
!!$            do m=1,LIS_rc%nensem(n)
!!$               mcmc_struc(t)%candidate(i,m) = decvals(i,(t-1)*LIS_rc%nensem(n)+m) 
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$    end subroutine RWMCMC_setdecSpaceValues


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

end module RWMCMCAlgorithm
