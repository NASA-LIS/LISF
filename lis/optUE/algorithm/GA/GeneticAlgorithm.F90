!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module GeneticAlgorithm
!BOP
!
! !MODULE: GeneticAlgorithm
!
! !DESCRIPTION: This module contains routines that define the operations
!  of a Genetic Algorithm (GA)
!  
! !REVISION HISTORY:
!  04 Feb 2008; Sujay Kumar; Initial Specification
!
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use GA_varctl

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: GAOpt_init
  public :: GAOpt_setup
  public :: GAOpt_run
  public :: GAOpt_checkConvergence
  public :: GAOpt_getdecSpaceValues
  public :: GAOpt_getNparam
  public :: GAOpt_readrestart
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
  type gastruc
     real,    allocatable :: fitness(:)
     integer, allocatable :: iparent(:,:)
     integer, allocatable :: ichild(:,:)
     real, allocatable    :: child(:,:)
     real,    allocatable :: parent(:,:)
     integer, allocatable :: ibest(:)
     real             :: best
     logical          :: fileopen
  end type gastruc

  type(gactl)            :: ga_ctl
  type(gastruc), allocatable :: ga_struc(:)

  contains
   
!BOP
! !ROUTINE: GAOpt_init
! \label{GAOpt_init}
! 
! !INTERFACE: 
    subroutine GAOpt_init()
! !USES: 
      use LIS_coreMod,         only : LIS_rc, LIS_config, LIS_vecTile
      use LIS_optUEMod,        only : LIS_decisionSpace
      use LIS_logMod,          only : LIS_logunit, LIS_getNextUnitNumber, &
           LIS_releaseUnitNumber, LIS_endrun, LIS_verify

! !DESCRIPTION: 
!   This routine performs the initialization steps for the GA.  
!   It initializes the required memory structures, and 
!   creates an initial random population, for each grid point. 
! 
!
!EOP      
      implicit none 

      integer                     :: i
      integer                     :: ftn
      integer                     :: n 
      integer                     :: status

! currently limited to one nest
      n = 1
      if(LIS_rc%nensem(n).le.1) then 
         write(LIS_logunit,*) '[ERR] The number of ensembles should be greater than 1'
         write(LIS_logunit,*) '[ERR] to run GA. A desired number is = 50'
         write(LIS_logunit,*) '[ERR] Stopping program ....'
         !call LIS_endrun()
      endif

      call ESMF_ConfigGetAttribute(LIS_config,ga_ctl%rfile,&
           label="GA restart file:",rc=status)
      call LIS_verify(status, 'GA restart file: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,ga_ctl%ngens,&
           label="GA number of generations:",default=100,rc=status)
      call LIS_verify(status, 'GA number of generations: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,ga_ctl%nchild,&
           label="GA number of children per parent:",default=1,rc=status)
      call LIS_verify(status, 'GA number of children per parent: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,ga_ctl%xoverscheme, &
           label="GA crossover scheme:",rc=status)
      call LIS_verify(status, 'GA crossover scheme: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,ga_ctl%pcross,&
           label="GA crossover probability:",default=0.5, rc=status)
      call LIS_verify(status, 'GA crossover probability: not defined')

      call ESMF_ConfigGetAttribute(LIS_config, ga_ctl%icreep, &
           label="GA mutation scheme:",default=0,rc=status)
      call LIS_verify(status, 'GA mutation scheme: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,ga_ctl%pcreep,&
           label="GA creep mutation probability:",default=0.04,rc=status)
      call LIS_verify(status, 'GA creep mutation probability: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,ga_ctl%pmutate,&
           label="GA jump mutation probability:",default=0.02,rc=status)
      call LIS_verify(status, 'GA jump mutation probability: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,ga_ctl%ielite,&
           label="GA use elitism:",default=1,rc=status)
      call LIS_verify(status, 'GA use elitism: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,ga_ctl%restart,&
           label="GA start mode:",rc=status)
      call LIS_verify(status, 'GA start mode: not defined')

      ga_ctl%npopsize = LIS_rc%nensem(n)

      call GAOpt_setup()
    end subroutine GAOpt_init

!BOP
! !ROUTINE: GAOpt_setup
! \label{GAOpt_setup}
!
! !INTERFACE: GAOpt_setup
    subroutine GAOpt_setup()
! !USES: 
      use LIS_coreMod,         only : LIS_rc, LIS_masterproc

      use LIS_optUEMod,        only : LIS_decisionSpace
      use LIS_logMod,          only : LIS_logunit, LIS_verify, LIS_endrun

      implicit none
! 
! !DESCRIPTION: 
!   This subroutine performs the second part of the GA initialization, 
!   The routine obtains the decision space object (from the models used
!   in optimization instance) and assigns the decision space variables
!   to the GA data structures
!  
!EOP
      integer                     :: n 
      integer                     :: status
      real                        :: rand      
      integer                     :: i
      integer                     :: j, nn, t,k,m
      type(ESMF_Field)            :: varField
      real,   pointer             :: vardata(:)

      n = 1
      
      call ESMF_StateGet(LIS_decisionSpace, itemCount=ga_ctl%nparam, &
           rc=status)
      call LIS_verify(status)
      
      allocate(ga_ctl%parmax(ga_ctl%nparam))
      allocate(ga_ctl%parmin(ga_ctl%nparam))
!      allocate(ga_ctl%pardel(ga_ctl%nparam))
!      allocate(ga_ctl%inp(ga_ctl%nparam))
      allocate(ga_ctl%nbp(ga_ctl%nparam))
      allocate(ga_ctl%npossbl(ga_ctl%nparam))
      allocate(ga_ctl%useSingleParamSet(ga_ctl%nparam))
      allocate(ga_ctl%useIntegerValues(ga_ctl%nparam))
      allocate(ga_ctl%vname(ga_ctl%nparam))

      ga_ctl%npossbl = 32768  ! 2^15 --> initialize all parameters to 15 bits/parameter

!      ga_ctl%seed = -1000  !kwh: unused

      allocate(ga_struc(LIS_rc%ntiles(n)/LIS_rc%nensem(n)))

      call ESMF_StateGet(LIS_decisionSpace, itemNameList=ga_ctl%vname,&
           rc=status)
      call LIS_verify(status)

      do i=1,ga_ctl%nparam
         call ESMF_StateGet(LIS_decisionSpace, trim(ga_ctl%vname(i)), &
              varField, rc=status)
         call LIS_verify(status)
         
         call ESMF_AttributeGet(varField,'MinRange',ga_ctl%parmin(i),&
              rc=status)
         call LIS_verify(status)

         call ESMF_AttributeGet(varField,'MaxRange',ga_ctl%parmax(i),&
              rc=status)
         call LIS_verify(status)
      enddo

!     do i=1,ga_ctl%nparam
!         ga_ctl%pardel(i) = ga_ctl%parmax(i)-ga_ctl%parmin(i)
!         ga_ctl%inp(i) = ga_ctl%pardel(i)/dble(ga_ctl%npossbl(i)-1)
!      enddo

      do i=1,ga_ctl%nparam
         !kwh: replace with: j=ga_ctl%nbp(i)=int(log(ga_ctl%npossbl(i)/log(2))) ?
         do j=1,30
            nn = 2**j
            if(nn.ge.ga_ctl%npossbl(i)) then 
               ga_ctl%nbp(i) = j
               exit
            endif
            if(j.ge.30) then 
               write(LIS_logunit,*) '[ERR] The chromosome size is set to 30.'
               write(LIS_logunit,*) '[ERR] program stopping ..'
               call LIS_endrun()
            endif
         enddo
      enddo
      
      !Check for inconsistencies
      ga_ctl%ngenes = 0 
      ga_ctl%sum_nposs_param = 0 
      ga_ctl%sum_nbp = 0 
      
      do i=1,ga_ctl%nparam
         ga_ctl%ngenes = ga_ctl%ngenes + ga_ctl%nbp(i)
         ga_ctl%sum_nposs_param = &
              ga_ctl%sum_nposs_param + ga_ctl%npossbl(i)
         ga_ctl%sum_nbp = ga_ctl%sum_nbp + (2**ga_ctl%nbp(i))
      enddo
      
      if(ga_ctl%sum_nposs_param.lt.ga_ctl%sum_nbp) then 
         write(LIS_logunit,*) '[ERR] Mismatch in sum_nposs_param and sum_nbp'
         write(LIS_logunit,*) '[ERR] Program stopping... '
         call LIS_endrun()
      endif
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         ! population (for each tile point)
         allocate(ga_struc(t)%iparent(ga_ctl%ngenes, LIS_rc%nensem(n)))
         allocate(ga_struc(t)%ichild(ga_ctl%ngenes, LIS_rc%nensem(n)))
         allocate(ga_struc(t)%child(ga_ctl%ngenes, LIS_rc%nensem(n)))
         allocate(ga_struc(t)%parent(ga_ctl%nparam, LIS_rc%nensem(n)))
         allocate(ga_struc(t)%fitness(LIS_rc%nensem(n)))
      enddo
      
      do k=1,ga_ctl%nparam
         call ESMF_StateGet(LIS_decisionSpace, trim(ga_ctl%vname(k)), &
              varField, rc=status)
         call LIS_verify(status)
         call ESMF_FieldGet(varField,localDE=0, farrayPtr=vardata,rc=status)
         call LIS_verify(status)
         do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            do m=1,LIS_rc%nensem(n)
               ga_struc(t)%parent(k,m) = vardata((t-1)*LIS_rc%nensem(n)+m)
            enddo
         enddo
         do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            do m=1,LIS_rc%nensem(n)            
               call encode(m,k,ga_struc(t)%parent, ga_struc(t)%iparent)
            enddo
         enddo
      enddo

      ga_ctl%genNo = 0
      
      ga_ctl%minfitness = -1.0E+20

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         allocate(ga_struc(t)%ibest(ga_ctl%ngenes))
         ga_struc(t)%best = ga_ctl%minfitness 
         ga_struc(t)%fitness = ga_ctl%minfitness 
      enddo
      
    end subroutine GAOpt_setup

!BOP
! 
! !ROUTINE: GAOpt_checkConvergence
! \label{GAOpt_checkConvergence}
! 
! !INTERFACE: 
    subroutine GAOpt_checkConvergence(check)
! !ARGUMENTS: 
      logical, intent(INOUT) :: check
! 
! !DESCRIPTION: 
!  This routine checks to see if the convergence criteria for GA 
!  is reached. In this case, the routine simply checks to see if the
!  specified number of generations is reached. 
!EOP
      if(ga_ctl%genNo.ge.ga_ctl%ngens) then 
         check = .true.
      else
         check = .false.
      endif

    end subroutine GAOpt_checkConvergence

!BOP
! !ROUTINE: GAOpt_run
! \label{GAOpt_run}
! 
! !INTERFACE: 
    subroutine GAOpt_run()
! !USES: 
      use LIS_coreMod,   only : LIS_rc
      use LIS_logMod,    only : LIS_logunit

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
      integer                  :: t
      integer                  :: n
      integer                  :: g,j

      n = 1

      write(LIS_logunit,*) '[INFO] Running GA, generation: ',ga_ctl%genNo

! kwh: calling writeGArestart *before* GA_recordFitness
! so that if time-consuming evaluation is not completed
! will start where left off; 
! otherwise already evaluated solution is re-evaluated
! also note cannot put it after GA_newGeneration() as
! writeGArestart() currently writes variables as LIS_dec_space knows
! them and not as the GA knows them; note: will overwrite initial 
! restart file but not big cost (note: could be avoided)

      call GA_recordFitness() 
      call writeGAoutput()
      write(LIS_logunit,*) '[INFO] wrote GA output'

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         call GA_selection(t)
         call GA_crossover(t)
         call GA_mutate(t)
      enddo
      write(LIS_logunit,*) '[INFO] finished GA operations'
      call GA_newGeneration()
      write(LIS_logunit,*) '[INFO] generated new GA generation'
      ga_ctl%genNo = ga_ctl%genNo + 1
      call writeGArestart()  
    end subroutine GAOpt_run

!BOP
! !ROUTINE: GAOpt_getdecSpaceValues
! \label{GAOpt_getdecSpaceValues}
! 
! !INTERFACE: 
    subroutine GAOpt_getdecSpaceValues(n)
! !USES: 
      use LIS_coreMod,    only : LIS_rc
      use LIS_optUEMod,   only : LIS_decisionSpace
      use LIS_logMod,     only : LIS_verify
! 
! !DESCRIPTION: 
!  This routine returns the array of decision space variables from 
!  the GA data structures
!EOP
      implicit none

      integer            :: n
      type(ESMF_Field)   :: varField(ga_ctl%nparam)
      real, pointer      :: vardata(:)
      integer            :: i, t, m
      integer            :: status
      do i=1,ga_ctl%nparam
         call ESMF_StateGet(LIS_decisionSpace, trim(ga_ctl%vname(i)), &
              varField(i), rc=status)
         call LIS_verify(status)
         
         call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata,rc=status)
         call LIS_verify(status)
         
         do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            do m=1,LIS_rc%nensem(n)
               vardata((t-1)*LIS_rc%nensem(n)+m)= ga_struc(t)%parent(i,m)
            enddo
         enddo
      enddo

    end subroutine GAOpt_getdecSpaceValues

!BOP
! 
! !ROUTINE: GAOpt_getNparam
! \label{GAOpt_getNparam}
! 
! !INTERFACE: 
    subroutine GAOpt_getNparam(nparam)
! 
! !DESCRIPTION: 
!  This method returns the number of decision space variables 
!
!EOP      
      integer   :: nparam
      
      nparam = ga_ctl%nparam

    end subroutine GAOpt_getNparam

!BOP
! 
! !ROUTINE: GA_recordFitness
!  \label{GA_recordFitness}
! 
! !INTERFACE: 
    subroutine GA_recordFitness()
! !USES:   
      use LIS_coreMod
      use LIS_optUEMod,        only : LIS_ObjectiveFunc, LIS_feasibleSpace
      use LIS_logMod,          only : LIS_logunit, LIS_verify, LIS_endrun

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

    call ESMF_StateGet(LIS_ObjectiveFunc,"Max Criteria Value",fitField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(fitField,localDE=0, farrayPtr=fitValue,rc=status)
    call LIS_verify(status)

    call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
    call LIS_verify(status)

    n = 1

    do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
       do m=1,LIS_rc%nensem(n)
          ga_struc(t)%fitness(m) = fitValue((t-1)*LIS_rc%nensem(n)+m)
          !penalize infeasible solutions
          if(mod_flag((t-1)*LIS_rc%nensem(n)+m).eq.1) then
             ga_struc(t)%fitness(m) = ga_ctl%minfitness
             !kwh: this should not occur as GA solutions are normalized to bounded range
             !but (for now) set feas flag to 1 in Least Squares obj if simply too few observations on which to opt
          endif
       end do
    enddo
    
    if (ga_ctl%genNo.eq.0) then !output gen 0 (default) values/fitness
       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
          ga_struc(t)%ibest = ga_struc(t)%iparent(:,1) !slot 1 is location of defaults
          ga_struc(t)%best = ga_struc(t)%fitness(1) 
       enddo
       call writeGAoutput()  !write out gen 0 'default' fitness and params
       ga_ctl%genNo = 1
       call writeGArestart()  !for testing the restarts, nothing more
    endif
    
!update best parameter set and associated fitness
    do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
       do m=1,LIS_rc%nensem(n)
          if(.not.(mod_flag((t-1)*LIS_rc%nensem(n)+m).eq.1)) then 
             if(ga_struc(t)%fitness(m).gt.ga_struc(t)%best) then 
                ga_struc(t)%ibest = ga_struc(t)%iparent(:,m)
                ga_struc(t)%best = ga_struc(t)%fitness(m)
             endif
          endif
       enddo
    enddo
  end subroutine GA_recordFitness
  
!BOP
! 
! 
! !ROUTINE: GA_selection
! \label{GA_selection}
!
! !INTERFACE: 
    subroutine GA_selection(t)
!
! !DESCRIPTION: 
!  This routine chooses pool of prospective parents via binary selection  
! 
!EOP
      implicit none
      integer       :: t
      integer       :: parent
      integer       :: ipick
      integer       :: j 
      integer       :: g

      ipick = ga_ctl%npopsize
      do j=1, ga_ctl%npopsize  
         call select(t, parent, ipick)
         do g=1,ga_ctl%ngenes
            ! this initializes ichild j to parent; in crossover will get some from mate2
            ga_struc(t)%ichild(g,j) = ga_struc(t)%iparent(g,parent)
!!!            if(ga_ctl%nchild.eq.2) then 
!!!               ga_struc(t)%ichild(g,j+1) = ga_struc(t)%iparent(g,mate2)
!!!            endif
         enddo
      enddo

      !kwh: why assigned to iparent prior to crossover step?
      !kwh: answer: appears ichild was just temp storage
      do j=1,ga_ctl%npopsize
         do g=1,ga_ctl%ngenes
            ga_struc(t)%iparent(g,j) = ga_struc(t)%ichild(g,j)
         enddo
      enddo
    end subroutine GA_selection
    
!BOP
! 
! !ROUTINE: GA_crossover
! \label{GA_crossover}
! 
! !INTERFACE: 
    subroutine GA_crossover(t)
! 
! !DESCRIPTION: 
!   This method performs the crossover or recombination operation, using 
!   either a two point crossover or a uniform crossover method. 
!   
!EOP      
    
      implicit none
      
      integer           :: t
      integer           :: j
      integer           :: mate1
      integer           :: mate2  
      logical           :: xovertemp(ga_ctl%ngenes) !xover template
      integer           :: k  !index variable 
      integer           :: K_cap  !number of times to xover and back again  (a multixover scheme)
      real     :: rand
      integer  :: g,i
      integer  :: icross1, icross2 ! two point crossover

      do j=1,ga_ctl%npopsize !, ga_ctl%nchild
         !kwh: why re-drawing mates after selection?
         call ran3(1,rand)
         mate1 = 1 + dint(dble(ga_ctl%npopsize-1)*rand)

         call ran3(1,rand)
         mate2 = 1 + dint(dble(ga_ctl%npopsize-1)*rand)

! two point crossover
         if(ga_ctl%xoverscheme.eq.1) then 
            xovertemp=.false.
            call ran3(1,rand)
            if(rand.le.ga_ctl%pcross) then
               call ran3(1,rand)
               K_cap=dint(dble(rand*4.0))+1  !up to four cross-crossovers
               do k=1,K_cap
                  call ran3(1,rand)
                  icross1 = dint(dble(ga_ctl%ngenes-1)*rand)+1
                  
                  call ran3(1,rand)
                  icross2 = dint(dble(ga_ctl%ngenes-1)*rand)+1
                  
                  i=icross1-1  !minus 1 to work in zero-based array space for modulo ops
                  g=icross1
                  do while (g.ne.icross2)
                     xovertemp(g)=.not.xovertemp(g)
!                     ga_struc(t)%ichild(g,j) = ga_struc(t)%iparent(g,mate2)
                     i=i+1
                     g=mod(i,ga_ctl%ngenes)+1
                  enddo
               enddo
               do g=1, ga_ctl%ngenes
                  if(xovertemp(g)) then
                     ga_struc(t)%ichild(g,j) = ga_struc(t)%iparent(g,mate2)
                  endif
               enddo
            endif
! uniform xover 
         elseif(ga_ctl%xoverscheme.eq.2) then 
            do g=1,ga_ctl%ngenes
               call ran3(1,rand)
               if(rand.le.ga_ctl%pcross) then 
                  ga_struc(t)%ichild(g,j) = ga_struc(t)%iparent(g,mate2)
!!!                  if(ga_ctl%nchild.eq.2) then 
!!!                     ga_struc(t)%ichild(g,j+1) = ga_struc(t)%iparent(g,mate1)
!!!                  endif
               endif
            enddo
         endif
      enddo

    end subroutine GA_crossover

!BOP
! 
! !ROUTINE: GA_mutate
! \label{GA_mutate}
! 
! !INTERFACE:
    subroutine GA_mutate(t)
      
      implicit none
      integer         :: t
! 
! !DESCRIPTION: 
! 
!  This subroutine performs mutations on the children generation. 
!  Random jump mutation or random creep mutation is performed based
!  on the mutation probability.
!EOP
      
      integer      :: nmutate
      integer      :: ncreep
      integer      :: j,k
      real         :: rand
      real         :: creep

      nmutate = 0 
      ncreep  = 0 
      
      do j=1,ga_ctl%npopsize  
         do k=1,ga_ctl%ngenes
! jump mutation
            call ran3(1,rand)
            if(rand.le.ga_ctl%pmutate) then 
               nmutate = nmutate + 1
               if(ga_struc(t)%ichild(k,j).eq.0) then 
                  ga_struc(t)%ichild(k,j) = 1
               else
                  ga_struc(t)%ichild(k,j) = 0 
               endif
            endif
         enddo
! creep mutation
         if(ga_ctl%icreep.ne.0) then 
            do k=1,ga_ctl%nparam
               call ran3(1,rand)
               if(rand.le.ga_ctl%pcreep) then 
                  call decode(j,ga_struc(t)%child,ga_struc(t)%ichild)
                  ncreep = ncreep + 1
                  creep = 1.0
                  call ran3(1,rand)
                  if(rand.lt.0.5) creep = -1.0
                  ga_struc(t)%child(k,j) = ga_struc(t)%child(k,j)+&
                       (ga_ctl%parmax(k)-ga_ctl%parmin(k))/dble(ga_ctl%npossbl(k)-1)*creep
                  if(ga_struc(t)%child(k,j).gt.ga_ctl%parmax(k)) then 
                     ga_struc(t)%child(k,j) = ga_ctl%parmax(k)-1.0*(ga_ctl%parmax(k)-ga_ctl%parmin(k))/dble(ga_ctl%npossbl(k)-1)
                  elseif(ga_struc(t)%child(k,j).lt.ga_ctl%parmin(k)) then 
                     ga_struc(t)%child(k,j) = ga_ctl%parmin(k)+1.0*(ga_ctl%parmax(k)-ga_ctl%parmin(k))/dble(ga_ctl%npossbl(k)-1)
                  endif
                  call encode(j,k,ga_struc(t)%child, ga_struc(t)%ichild)
               endif
            enddo
         endif
      end do

    end subroutine GA_mutate

!BOP
! !ROUTINE: GA_newGeneration
! \label{GA_newGeneration}
!
! !INTERFACE: 
    subroutine GA_newGeneration 
! !USES:
      use LIS_coreMod,  only : LIS_rc
      use LIS_optUEMod, only : LIS_decisionSpace
      use LIS_PE_HandlerMod, only : LIS_setPEDecisionSpace
      

      implicit none
! 
! !DESCRIPTION: 
! 
!   This method performs the computations to update the decision 
!   space variables based on the new generation organisms. 
!EOP
      integer      :: n 
      integer      :: t,j,g,k,i
      real         :: rand
      integer      :: irand
      
      n = 1
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do j=1,ga_ctl%npopsize 
            do g=1,ga_ctl%ngenes
               ga_struc(t)%iparent(g,j) = ga_struc(t)%ichild(g,j)
            enddo
         enddo
      enddo
     
!!!   BEGIN TEMP STRATEGY UNTIL LDT CAN PROCESS GA/UNC OUTPUT FILES
!!!      if(ga_ctl%ielite.eq.1) then !preserve best in population 
!!!         do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
!!!            call ran3(1,rand)
!!!            irand = 1 + dint(dble(ga_ctl%npopsize-1)*rand)
!!!            do g=1,ga_ctl%ngenes
!!!               ga_struc(t)%iparent(g,irand) = ga_struc(t)%ibest(g)  
!!!            enddo
!!!         enddo
!!!      endif
!!!
      if(ga_ctl%ielite.eq.1) then !preserve best in population 
         do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            call ran3(1,rand)
            irand = 1 + dint(dble(ga_ctl%npopsize-1)*rand)
            !COPY SLOT 1 TO IRAND...
            do g=1,ga_ctl%ngenes
               ga_struc(t)%iparent(g,irand)=ga_struc(t)%iparent(g,1)
            enddo
            !...THEN COPY BEST TO SLOT 1; EFFECTIVELY OVERWRITES IRAND SLOT
            do g=1,ga_ctl%ngenes
               ga_struc(t)%iparent(g,1) = ga_struc(t)%ibest(g)  
            enddo
         enddo
      endif
!!!
!!!
!!!   END TEMP STRATEGY UNTIL LDT CAN PROCESS GA/UNC OUTPUT FILES
!!!

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do j=1,ga_ctl%npopsize
            call decode(j,ga_struc(t)%parent, ga_struc(t)%iparent)
         enddo
      enddo


!integer parameter check
!      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) 
!         do k=1,ga_ctl%nparam
!            if(ga_ctl%useIntegerValues(k).eq.1) then !integer values
!               ga_struc(t)%parent(k,:) = nint(ga_struc(t)%parent(k,:))
!            endif
!         enddo
!      enddo

!!!    kwh:  NOT INITIALIZED ANYWHERE !!!!  JUST ALLOCATED.  FOR NOW, COMMENTING OUT
!!!      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) 
!!!         do k=1,ga_ctl%nparam
!!!            if(ga_ctl%useSingleParamSet(k).eq.1) then !domain-wide param
!!!               ga_struc(t)%parent(k,:) = ga_struc(1)%parent(k,:)
!!!            endif
!!!         enddo
!!!      enddo
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n) 
         !make sure that the binary code is consistent 
         do k=1,ga_ctl%nparam
            do i=1,LIS_rc%nensem(n) 
               call encode(i,k,ga_struc(t)%parent, ga_struc(t)%iparent)
            enddo
         enddo
      enddo

    end subroutine GA_newGeneration



!BOP
! 
! !ROUTINE: select
! \label{select}
! 
! !INTERFACE: 
    subroutine select(t, mate, ipick)
! 
! !DESCRIPTION: 
!  This subroutine selects the better of two possible parents for mating. 
! 
!EOP
      implicit none
      integer    :: t
      integer    :: mate
      integer    :: ipick
      integer    :: ifirst
      integer    :: isecond
      real       :: rand
      
      call ran3(1,rand)
      ifirst = 1 + dint(dble(ga_ctl%npopsize-1)*rand)
      call ran3(1,rand)
      isecond = 1 + dint(dble(ga_ctl%npopsize-1)*rand)
      if(ga_struc(t)%fitness(ifirst).ge.ga_struc(t)%fitness(isecond)) then 
         mate = ifirst
      else
         mate = isecond
      endif      
#if 0 
      if(ipick+1.gt.ga_ctl%npopsize) call shuffle(t, ipick)
      ifirst = ipick
      isecond = ipick+1
      ipick = ipick+2
      if(ga_struc(t)%fitness(ifirst).gt.ga_struc(t)%fitness(isecond)) then 
         mate = ifirst
      else
         mate = isecond
      endif      
#endif
      
    end subroutine select

!!!!BOP
!!!! !ROUTINE: shuffle
!!!! \label{shuffle}
!!!! 
!!!! !INTERFACE: 
!!!    subroutine shuffle(t, ipick)
!!!
!!!      implicit none
!!!
!!!      integer :: t
!!!      integer :: ipick
!!!! 
!!!! !DESCRIPTION: 
!!!!  This routine shuffles the parent array and its corresponding fitness
!!!! 
!!!!EOP
!!!      integer :: j,g
!!!      real    :: rand
!!!      integer :: iother
!!!      integer :: itemp
!!!      real    :: temp
!!!      
!!!      print*, 'in shuffle ',t
!!!      ipick = 1
!!!      do j=1,ga_ctl%npopsize-1
!!!         call ran3(1,rand)
!!!         iother = j + 1 + dint(dble(ga_ctl%npopsize-j)*rand)
!!!         do g = 1, ga_ctl%ngenes
!!!            itemp = ga_struc(t)%iparent(g,iother)
!!!            ga_struc(t)%iparent(g,iother) = ga_struc(t)%iparent(g,j)
!!!            ga_struc(t)%iparent(g,j) = itemp
!!!         enddo
!!!         temp = ga_struc(t)%fitness(iother)
!!!         ga_struc(t)%fitness(iother) = ga_struc(t)%fitness(j)
!!!         ga_struc(t)%fitness(j) = temp
!!!      enddo
!!!    end subroutine shuffle

!BOP
! 
! !ROUTINE: encode
! \label{encode}
!
! !INTERFACE: 
  subroutine encode(j,k,array,iarray)
! !ARGUMENTS: 
    integer  :: j
    integer  :: k
    real     :: array(ga_ctl%nparam, ga_ctl%npopsize)
    integer  :: iarray(ga_ctl%ngenes, ga_ctl%npopsize)
! 
! !DESCRIPTION: 
!
!  This routine converts a real number to a binary string
!
!EOP
    integer  :: m,i,iparam,istart
    real     :: f !fraction
    istart = 1
    do i=1,k-1
       istart = istart + ga_ctl%nbp(k)
    enddo

    m = ga_ctl%nbp(k)-1
!!!!!    if((ga_ctl%parmax(k)-ga_ctl%parmin(k))/dble(ga_ctl%npossbl(k)-1).eq.0) return
!!!!!    iparam= nint((array(k,j)-ga_ctl%parmin(k))/(ga_ctl%parmax(k)-ga_ctl%parmin(k))/dble(ga_ctl%npossbl(k)-1))
    f=(array(k,j)-ga_ctl%parmin(k))/(ga_ctl%parmax(k)-ga_ctl%parmin(k))
    iparam= nint(f*dble(ga_ctl%npossbl(k)-1))
    do i=istart, istart + ga_ctl%nbp(k) - 1
       iarray(i,j) = 0 
       if((iparam+1).gt.(2**m)) then 
          iarray(i,j) = 1
          iparam = iparam - 2**m 
       endif
       m = m-1
    enddo
  end subroutine encode
  
!BOP
! 
! !ROUTINE: decode
! \label{decode}
!
! !INTERFACE: 
  subroutine decode(i,array,iarray)
! !ARGUMENTS: 
    integer  :: i 
    real     :: array(ga_ctl%nparam, ga_ctl%npopsize)
    integer  :: iarray(ga_ctl%ngenes, ga_ctl%npopsize)
! 
! !DESCRIPTION: 
!
!  This routine converts a binary string to a real number.  
!
!EOP
    integer  :: k,iparam, l,m,j
    real     :: f  !fraction
    l = 1
    do k=1, ga_ctl%nparam
       iparam = 0 
       m = l 
       do j=m, m + ga_ctl%nbp(k) - 1
          l = l + 1
          iparam = iparam + iarray(j,i)*(2**(m+ga_ctl%nbp(k)-1-j))
       enddo
       f=float(iparam)/dble(ga_ctl%npossbl(k)-1)
       array(k,i) = ga_ctl%parmin(k) + (ga_ctl%parmax(k)-ga_ctl%parmin(k))*f
    enddo
  end subroutine decode

!BOP
! 
! !ROUTINE: writeGArestart
! \label{writeGArestart}
! 
! !INTERFACE: 
  subroutine writeGArestart
! !USES: 
    use LIS_optUEMod,   only : LIS_decisionSpace
    use LIS_coreMod,    only : LIS_rc, LIS_masterproc
    use LIS_fileIOMod,  only : LIS_create_output_directory
    use LIS_logMod,     only : LIS_logunit, LIS_verify
    use LIS_historyMod, only : LIS_writevar_restart
! 
! !DESCRIPTION: 
! 
! This routine writes the checkpoint data for a GA restart
! 
!EOP
    integer             :: n 
    integer             :: i,m,t
    integer             :: status
    character(len=LIS_CONST_PATH_LEN) :: filen
    character (len=4)   :: fgen
    character*100       :: vnames(ga_ctl%nparam)
    type(ESMF_Field)    :: varField(ga_ctl%nparam)
    real, allocatable       :: vardata(:)

    n = 1
    allocate(vardata(LIS_rc%ntiles(n)))
    call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
         rc=status)
    call LIS_verify(status)
    
    if(LIS_masterproc) then 
       write(LIS_logunit,*) '[INFO] writing GA restart.. ',trim(filen)
       call LIS_create_output_directory('GA')
       write(unit=fgen, fmt='(i4.4)') ga_ctl%genNo
       filen = trim(LIS_rc%odir)//'/GA/GA.'&
            //trim(fgen)//'.GArst'
       open(40,file=filen,status='unknown',form='unformatted')
       write(40) ga_ctl%genNo
    endif

!write out in order of ga_ctl parameters list
    do i=1,ga_ctl%nparam
       do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
          do m=1,LIS_rc%nensem(n)
             vardata((t-1)*LIS_rc%nensem(n)+m)= ga_struc(t)%parent(i,m)
          enddo
       enddo
       
      call LIS_writevar_restart(40,n,vardata)
    enddo
       
    if(LIS_masterproc) then 
       close(40)
       write(LIS_logunit,*) '[INFO] GA checkpoint file written ',trim(filen)
    endif
  end subroutine writeGArestart


!BOP
! 
! !ROUTINE: GAopt_readrestart
! \label{GAopt_readrestart}
! 
! !INTERFACE: 
  subroutine GAopt_readrestart
! !USES: 
    use LIS_optUEMod,        only : LIS_decisionSpace
    use LIS_coreMod,         only : LIS_rc
    use LIS_logMod,          only : LIS_logunit, LIS_verify
    use LIS_historyMod,      only : LIS_readvar_restart
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
! 
!   This routine reads the checkpoint data for a GA restart
!EOP
    integer             :: n 
    integer             :: i, t, m
    integer             :: status
    character*100       :: vnames(ga_ctl%nparam)
    type(ESMF_Field)    :: varField(ga_ctl%nparam)
    real, allocatable       :: vardata(:)
!    real, allocatable       :: tempdata(:,:)
    if(ga_ctl%restart.eq."restart") then !restart run
       n = 1    
       write(LIS_logunit,*) '[INFO] Reading the GA restart file ..'
       allocate(vardata(LIS_rc%ntiles(n)))
!       allocate(tempdata(ga_ctl%nparam,LIS_rc%nensem(n)))

       call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
            rc=status)
       call LIS_verify(status)
       
       open(40,file=ga_ctl%rfile,form='unformatted')        
       
       read(40) ga_ctl%genNo
       write(LIS_logunit,*) '[INFO] Generation Number ',ga_ctl%genNo

!read in, as with write_restart, in order of ga_ctl parameters list
!which should be same as ESMF decision space state       
       do i=1,ga_ctl%nparam
          call LIS_readvar_restart(40,n,vardata)       
          
          do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
             do m=1,LIS_rc%nensem(n)
                ga_struc(t)%parent(i,m) = vardata((t-1)*LIS_rc%nensem(n)+m)
                ! make sure that the binary representation is consistent
                call encode(m,i,ga_struc(t)%parent, ga_struc(t)%iparent)
!                call decode(m,tempdata, ga_struc(t)%iparent) !******temp******kwh
             enddo
          enddo
      enddo



      
      close(40)
      write(LIS_logunit,*) '[INFO] Finished reading the GA restart file ..'
   endif
 end subroutine GAopt_readrestart
 

!BOP
! 
! !ROUTINE: writeGAoutput
! \label{writeGAoutput}
! 
! !INTERFACE: 
  subroutine writeGAoutput
! !USES: 
    use LIS_optUEMod
    use LIS_coreMod
    use LIS_fileIOMod
    use LIS_logMod
    use LIS_historyMod
! 
! !DESCRIPTION: 
! 
! This routine writes a gridded output file, including GA fitness and GA best set of parameters
! 
!EOP
    integer             :: n 
    character(len=LIS_CONST_PATH_LEN) :: filen
    character (len=4)   :: fgen
    real, allocatable       :: fitness(:)
    real, allocatable       :: avgfitness(:)
    integer             :: fitcount
    real, allocatable       :: objs(:,:)
    integer             :: ftn 
    integer             :: k,iparam, t,m,j,l
    real                :: f  !fraction
! Write the best organism from each population. 
! true if subgrid tiling is not turned on. The idea 
! here is to write a 'gridded fitness' field. 

    n = 1
    allocate(fitness(LIS_rc%ngrid(n)))
    allocate(avgfitness(LIS_rc%ngrid(n)))
    allocate(objs(LIS_rc%ngrid(n),ga_ctl%nparam))
    avgfitness = 0.0

    if(LIS_masterproc) then 
       ftn  = LIS_getNextUnitNumber()
       call LIS_create_output_directory('GA')
       write(unit=fgen, fmt='(i4.4)') ga_ctl%genNo
       filen = trim(LIS_rc%odir)//'/GA/GA.'&
            //trim(fgen)//'.1gd4r'
       open(ftn,file=filen,status='unknown',form='unformatted')
       write(ftn) ga_ctl%nparam
    endif

!decoding the ibest objectives
    do t=1,LIS_rc%ngrid(n)
       l = 1
       do k=1, ga_ctl%nparam
          iparam = 0 
          m = l 
          do j=m, m + ga_ctl%nbp(k) - 1
             l = l + 1
             iparam = iparam + ga_struc(t)%ibest(j)*(2**(m+ga_ctl%nbp(k)-1-j))
          enddo
          f=float(iparam)/dble(ga_ctl%npossbl(k)-1)
          objs(t,k) = ga_ctl%parmin(k) + (ga_ctl%parmax(k)-ga_ctl%parmin(k))*f
       enddo
    enddo
!write the ibest objectives
    do k=1,ga_ctl%nparam
       if(LIS_masterproc) then
          write(ftn) ga_ctl%vname(k)
       endif
       call LIS_writevar_gridded(ftn,n,objs(:,k),wopt="2d gridspace")
    enddo
    do t=1,LIS_rc%ngrid(n)
       fitness(t) = ga_struc(t)%best
       fitcount = 0 
       do m=1,LIS_rc%nensem(n)
          if(ga_struc(t)%fitness(m).ne.ga_ctl%minfitness) then 
             avgfitness(t) = avgfitness(t) + ga_struc(t)%fitness(m)
             fitcount = fitcount + 1
          endif
       enddo

       if(fitcount.ne.0) avgfitness(t) = avgfitness(t)/float(fitcount)
    enddo
    call LIS_writevar_gridded(ftn,n,fitness, wopt="2d gridspace")
    call LIS_writevar_gridded(ftn,n,avgfitness, wopt="2d gridspace")

    if(LIS_masterproc) then 
       call LIS_releaseUnitNumber(ftn)
    endif
    deallocate(fitness)
    deallocate(avgfitness)
    deallocate(objs)
  end subroutine writeGAoutput

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


end module GeneticAlgorithm
