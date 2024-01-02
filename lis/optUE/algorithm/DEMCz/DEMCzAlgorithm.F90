!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module DEMCzAlgorithm
!BOP
!
! !MODULE: DEMCzAlgorithm
!
! !DESCRIPTION: 
! This module contains routines that define the operations
! of the DEMCz algorithm.  It is a parallel implementation of 
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
  use DEMCz_varctl
  use LIS_coreMod
  use LIS_optUEMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: DEMCz_init
  public :: DEMCz_setup
  public :: DEMCz_run
  public :: DEMCz_checkConvergence
  public :: DEMCz_getdecSpaceValues
  public :: DEMCz_readrestart
  public :: DEMCz_getNparam
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
  type demczstruc
     real,    allocatable :: X(:,:)
     real,    allocatable :: X_cand(:,:)
     real,    allocatable :: LnZ(:)
!     real,    allocatable :: LnZ_cand(:)
     integer, allocatable :: iparent1(:)
     integer, allocatable :: iparent2(:)
     real             :: acceptcount
  end type demczstruc

  type(demczctl)        :: demcz_ctl
  type(demczstruc), allocatable :: demcz_struc(:)

  contains
   
!BOP
! !ROUTINE: DEMCz_init
! \label{DEMCz_init}
! 
! !INTERFACE: 
    subroutine DEMCz_init()
! !USES: 

! !DESCRIPTION: This routine performs the initialization steps for DEMCz.  
!   It initializes the required memory structures, and 
!   creates an initial random population, for each grid point. 

!   
!
!EOP      
      implicit none 
      
      integer                :: status

      demcz_ctl%iterNo = 0 
      demcz_ctl%zerothRun = .true. 
      demcz_ctl%seed = -1000
      demcz_ctl%minfitness = -1E20
      demcz_ctl%optrfile = ''

!      call ESMF_ConfigGetAttribute(LIS_config,demcz_ctl%decspaceAttribsFile,&
!           label="DEMCz decision space attributes file:",rc=status)
!      call LIS_verify(status, 'DEMCz decision space attributes file: not defined')
!
      call ESMF_ConfigGetAttribute(LIS_config,demcz_ctl%restart,&
           label="DEMCz start mode:",rc=status)
      call LIS_verify(status, 'DEMCz start mode: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,demcz_ctl%optrfile,&
           label="DEMCz GA restart file:",rc=status)

      call ESMF_ConfigGetAttribute(LIS_config,demcz_ctl%rfile,&
           label="DEMCz restart file:",rc=status)

      call ESMF_ConfigGetAttribute(LIS_config,demcz_ctl%maxiter,&
           label="DEMCz number of iterations:",rc=status)
      call LIS_verify(status, 'DEMCz number of generations: not defined')

! The parameter, 'pert_spread', is used to specify the sigma for all parameters 
! instead of specifying sigma for each parameter.  It does this by basing sigma, used in 'perturb_step',
! on this parameter and the width (ie, max-min) of the parameter range.

      call ESMF_ConfigGetAttribute(LIS_config,demcz_ctl%pert_spread,&
           label="DEMCz perturbation factor:",default=0.001,rc=status)
      call LIS_verify(status, 'DEMCz perturbation factor: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,demcz_ctl%modehopfreq,&
           label="DEMCz mode hopping frequency:",default=0.1,rc=status)
      call LIS_verify(status, 'DEMCz mode hopping frequency: not defined')

      call DEMCz_setup()
      
   end subroutine DEMCz_init

!BOP
! !ROUTINE: DEMCz_setup
! \label{DEMCz_setup}
!
! !INTERFACE: DEMCz_setup
    subroutine DEMCz_setup()
! !USES: 
! 
! !DESCRIPTION: 
!   This subroutine performs the second part of the DEMCz initialization, 
!   The routine obtains the decision space object (from the models used
!   in optimization instance) and assigns the decision space variables
!   to the DEMCz data structures
!  
!EOP
      implicit none
      
      type(ESMF_Field)           :: varField
      real,   pointer            :: vardata(:)
      integer                    :: t, k, n, m, j
      integer                    :: status
      integer                    :: tiles_per_ens

      n = 1

      tiles_per_ens = LIS_rc%ntiles(n) / LIS_rc%nensem(n)

      call ESMF_StateGet(LIS_decisionSpace, itemCount=demcz_ctl%nparam, &
           rc=status)
      call LIS_verify(status)

      allocate(demcz_ctl%vname(demcz_ctl%nparam))
      allocate(demcz_ctl%parmax(demcz_ctl%nparam))
      allocate(demcz_ctl%parmin(demcz_ctl%nparam))

      allocate(demcz_struc(tiles_per_ens))

      do t=1,tiles_per_ens
! population (for each tile point), each ensemble member represents a Markov chain. 

         allocate(demcz_struc(t)%X(demcz_ctl%nparam, LIS_rc%nensem(n)))
         allocate(demcz_struc(t)%X_cand(demcz_ctl%nparam, LIS_rc%nensem(n)))
         allocate(demcz_struc(t)%LnZ(LIS_rc%nensem(n)))

         demcz_struc(t)%acceptcount=0
      enddo

      call ESMF_StateGet(LIS_decisionSpace, itemNameList=demcz_ctl%vname,&
           rc=status)
      call LIS_verify(status)

      do k=1,demcz_ctl%nparam
         !Assumed that the ESMF ordering remains unchanged through run
         call ESMF_StateGet(LIS_decisionSpace, trim(demcz_ctl%vname(k)), &
              varField, rc=status)
         call LIS_verify(status)
         
         call ESMF_AttributeGet(varField, 'MinRange',demcz_ctl%parmin(k),&
              rc=status)
         call LIS_verify(status, 'setting minrange to decspace obj in DEMCz')
         call ESMF_AttributeGet(varField, 'MaxRange',demcz_ctl%parmax(k),&
              rc=status)
         call LIS_verify(status, 'setting maxrange to decspace obj in DEMCz')

         call ESMF_FieldGet(varField,localDE=0, farrayPtr=vardata,rc=status)
         call LIS_verify(status)
         
         do t=1,tiles_per_ens
            do m=1,LIS_rc%nensem(n)
               demcz_struc(t)%X(k,m) = vardata((t-1)*LIS_rc%nensem(n)+m)
            enddo
         enddo
      enddo      

    end subroutine DEMCz_setup

!BOP
! 
! !ROUTINE: DEMCz_checkConvergence
! \label{DEMCz_checkConvergence}
! 
! !INTERFACE: 
    subroutine DEMCz_checkConvergence(check)
! !ARGUMENTS: 
      logical, intent(INOUT) :: check
! 
! !DESCRIPTION: 
!  This routine checks to see if the convergence criteria for DEMCz 
!  is reached. In this case, the routine simply checks to see if the
!  specified number of generations is reached.
!EOP
      if(demcz_ctl%iterNo.ge.demcz_ctl%maxiter) then 
         check = .true.
      else
         check = .false.
      endif

    end subroutine DEMCz_checkConvergence

!BOP
! !ROUTINE: DEMCz_run
! \label{DEMCz_run}
! 
! !INTERFACE: 
    subroutine DEMCz_run()
! !USES: 

! !DESCRIPTION: 
!
!
!EOP  
      implicit none

      integer                  :: n

      n = 1
      write(LIS_logunit, *)  'DEMCz iteration ',demcz_ctl%iterNo 

      if(demcz_ctl%iterNo.eq.0) then 
         call DEMCz_recordLnZ_X0()
      else                                
         call DEMCz_acceptReject()
      endif

      call writeDEMCzrestart()  !really an output file, it overwrites restart file with accepted values

      call DEMCz_drawCandidates()
      demcz_ctl%iterNo = demcz_ctl%iterNo + 1

      call writeDEMCzrestart()    !will write OLD LnZ and CANDIDATE X, the needed info for restart 
    end subroutine DEMCz_run

!BOP
! 
! !ROUTINE: DEMCz_setPreLnZ
! \label{DEMCz_setPreLnZ}
!
! !INTERFACE: 
    subroutine  DEMCz_recordLnZ_X0()
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
      integer            :: m_best
      real               :: p_best

      n = 1

      call ESMF_StateGet(LIS_ObjectiveFunc,"Max Criteria Value",fitField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(fitField,localDE=0, farrayPtr=fitValue,rc=status)
      call LIS_verify(status)
      
      call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
      call LIS_verify(status)


      ! This iteration is setting the fitness values for the Xs
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            demcz_struc(t)%LnZ(m) = fitValue((t-1)*LIS_rc%nensem(n)+m)
            !penalize infeasible solutions
            if(mod_flag((t-1)*LIS_rc%nensem(n)+m).eq.1) then 
               demcz_struc(t)%LnZ(m) = demcz_ctl%minfitness
            endif
         end do
      enddo

! Double check parameter bounds; Actually needed for the error variable (sigma) that is
! not currently thought of as a parameter of the LSM
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            do j=1,demcz_ctl%nparam
               if((demcz_struc(t)%X(j,m) .lt. demcz_ctl%parmin(j))&
                    .or.(demcz_struc(t)%X(j,m) .gt. demcz_ctl%parmax(j))) then
                  demcz_struc(t)%LnZ(m) = demcz_ctl%minfitness
               endif
            enddo
         enddo
         
      enddo

!TEMP IMPLEMENTATION OF OPTIMIZATION RESTART; LATER MOVE TO LDT
      if(demcz_ctl%optrfile.ne.'') then !clean up lisconfig entry
         do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            p_best=demcz_ctl%minfitness  !initialize to worst      
            do m=1,LIS_rc%nensem(n)
               if(demcz_struc(t)%LnZ(m).gt.p_best) then
                  p_best=demcz_struc(t)%LnZ(m)
                  m_best=m
               endif
            end do
            !cycle through to set all at best for seeding
            do m=1,LIS_rc%nensem(n)
               demcz_struc(t)%LnZ(m)=p_best
               demcz_struc(t)%X(:,m)=demcz_struc(t)%X(:,m_best)
            end do
         enddo
      endif
    end subroutine DEMCz_recordLnZ_X0

!BOP
! 
! !ROUTINE: DEMCz_drawCandidates
! \label{DEMCz_drawCandidates}
!
! !INTERFACE: 
    subroutine DEMCz_drawCandidates()
! !USES: 
    use LIS_historyMod,      only : LIS_readvar_restart
! !DESCRIPTION: 
! 
!EOP      
      implicit none
      integer                  :: n, t, m, i
      integer                  :: iparentGen  !To draw random parent, draw generation...
      integer                  :: m_p  !then draw member from within Generation
      real                     :: rand      
      real, pointer            :: vardata(:) 
      character*100            :: vnames(demcz_ctl%nparam)
      character(len=LIS_CONST_PATH_LEN) :: filen
      integer :: iteratNo
      real    :: gamma
      real    :: epsilon, b_param
      integer            :: status
      integer :: q
      integer :: SIGN
      character (len=4)   :: fgen
      character*100       :: vname_in_file

      n = 1

      !Initialize candidate solution to current; later add difference vector
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
         do m=1,LIS_rc%nensem(n)
            demcz_struc(t)%X_cand(:,m)=demcz_struc(t)%X(:,m)
         enddo
      enddo

      !Set gamma = gamma*, the optimal for multivariate post. dist. as per ter Braak (2006)
      if (demcz_ctl%modehopfreq.ne.0) then
         if (mod(demcz_ctl%iterNo,int(1.0/demcz_ctl%modehopfreq)).eq.0) then
            gamma = 0.98  ! see ter Braak (2006) for explanation why preferable to 1.0; 
         else
            gamma = 2.38/sqrt(2.0*dble(demcz_ctl%nparam))
         endif
      else  !is equal to zero
         gamma = 2.38/sqrt(2.0*dble(demcz_ctl%nparam))
      endif

      call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
           rc=status)
      call LIS_verify(status)
      
      !main transition formula: cand =cand + gamma*(parent1 - parent2) + epsilon
      !Here done in stages to limit memory requirements
      !cand=cand+gamma*parent2-gamma*parent1 +epsilon
      allocate(vardata(LIS_rc%ntiles(n)))
      do q=1,2  !iterate over parents
         if (q.eq.1) then 
            SIGN =  1
         elseif (q.eq.2) then
            SIGN = -1
         endif
         
         !Find parent indices randomly
         call ran3(1,rand)
         iparentGen= int(rand*demcz_ctl%iterNo) !+1; because also have gen = 0; which is needed initially
         
         write(unit=fgen, fmt='(i4.4)') iparentGen
         filen = trim(LIS_rc%odir)//'/DEMCz/DEMCz.'&
              //trim(fgen)//'.DEMCzrst'
         
         write(LIS_logunit,*) 'Reading DEMCz parent ', q
         
         open(40,file=filen,form='unformatted')        
         
         read(40) iteratNo
         write(LIS_logunit,*) 'Parent drawn from: ',iteratNo
         
         !Read random parent1
         do i=1,demcz_ctl%nparam
            
            read(40) vname_in_file
            
            call LIS_readvar_restart(40,n,vardata)       
            
            do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
               do m=1,LIS_rc%nensem(n)
                  !Find random member of ensemble
                  call ran3(1,rand)
                  m_p= int(rand*LIS_rc%nensem(n))+1
                  demcz_struc(t)%X_cand(i,m) = demcz_struc(t)%X_cand(i,m) &
                       + SIGN*gamma*vardata((t-1)*LIS_rc%nensem(n)+m_p) 
               enddo
            enddo
         enddo
         
         close(40)
         write(LIS_logunit,*) 'Finished reading DEMCz parent ', q
      enddo
      deallocate(vardata)

      !Add epsilon
      do i=1,demcz_ctl%nparam
         b_param = (demcz_ctl%parmax(i) - demcz_ctl%parmin(i))&
              *demcz_ctl%pert_spread
!         call ran3(1,rand)
!         epsilon = ((rand-0.5)*2)*b_param
         do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
            do m=1,LIS_rc%nensem(n)
               call ran3(1,rand)
               epsilon = ((rand-0.5)*2)*b_param
               demcz_struc(t)%X_cand(i,m) = demcz_struc(t)%X_cand(i,m) + epsilon
            enddo
         enddo
      enddo
      
    end subroutine DEMCz_drawCandidates
    
!BOP
! 
! !ROUTINE: DEMCz_acceptReject
! \label{DEMCz_acceptReject}
!
! !INTERFACE: 
    subroutine DEMCz_acceptReject()
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
      logical                  :: accept
      integer                  :: t,m,n,j
      integer                  :: count
      logical                  :: temp_boolean
      real                     :: LnZ_cand
      real                     :: rand, rand_t
      real                     :: ln_p_ratio

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
            LnZ_cand = fitValue((t-1)*LIS_rc%nensem(n)+m)
            !penalize infeasible solutions
            if(mod_flag((t-1)*LIS_rc%nensem(n)+m).eq.1) then 
               LnZ_cand = demcz_ctl%minfitness
            endif
            do j=1,demcz_ctl%nparam
               if((demcz_struc(t)%X_cand(j,m) .lt. demcz_ctl%parmin(j)) &
                    .or.(demcz_struc(t)%X_cand(j,m) .gt. demcz_ctl%parmax(j))) then
                  LnZ_cand = demcz_ctl%minfitness
               endif
            end do

            !Accept/reject test
            call ran3(1,rand)
            
            !following assumes log-transformed probability (Z)
            ln_p_ratio = LnZ_cand - demcz_struc(t)%LnZ(m)
            rand_t = log(rand)
            
            if(rand_t.lt.ln_p_ratio) then 
               accept = .true.
            else
               accept = .false.
            endif
            
            if (accept) then 
               demcz_struc(t)%X(:,m)=demcz_struc(t)%X_cand(:,m)                     
               demcz_struc(t)%LnZ(m)=LnZ_cand                     
               demcz_struc(t)%acceptcount=demcz_struc(t)%acceptcount+1
            else
               !Leave unchanged, ie, X=X 
            endif
         enddo
      enddo
            
 end subroutine DEMCz_acceptReject

!!!!!BOP
!!!!! 
!!!!! !ROUTINE: DEMCz_accept_reject
!!!!! \label{DEMCz_accept_reject}
!!!!!
!!!!! !INTERFACE: 
!!!!    subroutine DEMCz_accept_reject()
!!!!!
!!!!! !DESCRIPTION: 
!!!!!
!!!!!EOP
!!!!
!!!!      implicit none
!!!!      
!!!!      integer                  :: status
!!!!      logical                  :: accept
!!!!      integer                  :: t,m,n,j
!!!!
!!!!      n = 1
!!!!
!!!!      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
!!!!         do m=1,LIS_rc%nensem(n)
!!!!            call accept_reject_test(t,m, accept)
!!!!            if (accept) then 
!!!!               demcz_struc(t)%X(:,m)=demcz_struc(t)%X_cand(:,m)                     
!!!!               demcz_struc(t)%LnZ(m)=demcz_struc(t)%LnZ_cand(m)                     
!!!!               demcz_struc(t)%acceptcount=demcz_struc(t)%acceptcount+1
!!!!            else
!!!!               !Leave unchanged, ie, X=X 
!!!!            endif
!!!!         enddo
!!!!      enddo
!!!!      
!!!! end subroutine DEMCz_accept_reject


!!!!!BOP
!!!!! 
!!!!! !ROUTINE: accept_reject_test
!!!!! \label{accept_reject_test}
!!!!!
!!!!! !INTERFACE:       
!!!!    subroutine accept_reject_test(t,i,accept)
!!!!!
!!!!! !DESCRIPTION: 
!!!!!
!!!!!EOP
!!!!      implicit none
!!!!
!!!!      integer :: n 
!!!!      integer :: t, i
!!!!      logical, intent(out) :: accept
!!!!      real    :: rand, rand_t
!!!!      real    :: ln_p_ratio
!!!!      real    :: LnZ_cand
!!!!      real    :: LnZ
!!!!
!!!!      n = 1
!!!!      
!!!!!compare fitnesses of X and candidates. 
!!!!      call ran3(1,rand)
!!!!      
!!!!      LnZ_cand = demcz_struc(t)%LnZ_cand(i)
!!!!      LnZ = demcz_struc(t)%LnZ(i)
!!!!
!!!!!assumes that ll values in log form
!!!!      ln_p_ratio = LnZ_cand - LnZ
!!!!      rand_t = log(rand)
!!!!      
!!!!      if(rand_t.lt.ln_p_ratio) then 
!!!!         accept = .true.
!!!!      else
!!!!         accept = .false.
!!!!      endif
!!!!      
!!!!  end subroutine accept_reject_test


!BOP
! 
! !ROUTINE: writeDEMCzrestart
! \label{writeDEMCzrestart}
! 
! !INTERFACE: 
  subroutine writeDEMCzrestart
! !USES: 
    use LIS_fileIOMod,  only : LIS_create_output_directory
    use LIS_historyMod, only : LIS_writevar_restart
! 
! !DESCRIPTION: 
! 
! This routine writes the checkpoint data for a DEMCz restart
! 
!EOP

    integer             :: n 
    integer             :: i,m,t
    integer             :: status
    character*100       :: filen
    character (len=4)   :: fgen
    character*100       :: vnames(demcz_ctl%nparam)
    real, allocatable       :: vardata(:)

    n = 1

    allocate(vardata(LIS_rc%ntiles(n)))

    if(LIS_masterproc) then 
       call LIS_create_output_directory('DEMCz')
       write(unit=fgen, fmt='(i4.4)') demcz_ctl%iterNo
       filen = trim(LIS_rc%odir)//'/DEMCz/DEMCz.'&
            //trim(fgen)//'.DEMCzrst'
       open(40,file=filen,status='unknown',form='unformatted')
       write(40) demcz_ctl%iterNo
    endif

    do i=1,demcz_ctl%nparam
       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
          do m=1,LIS_rc%nensem(n)
             vardata((t-1)*LIS_rc%nensem(n)+m) = demcz_struc(t)%X(i,m)
          enddo
       enddo
       if(LIS_masterproc) then
          write(40) demcz_ctl%vname(i)
       endif
       call LIS_writevar_restart(40,n,vardata)
    enddo

    do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
       do m=1,LIS_rc%nensem(n)
          vardata((t-1)*LIS_rc%nensem(n)+m) = demcz_struc(t)%LnZ(m)
       enddo
    enddo
    call LIS_writevar_restart(40,n,vardata)
       
    if(LIS_masterproc) then 
       close(40)
       write(LIS_logunit,*) 'DEMCz checkpoint file written ',trim(filen)
    endif
    deallocate(vardata)

  end subroutine writeDEMCzrestart


!BOP
! 
! !ROUTINE: DEMCz_readrestart
! \label{DEMCz_readrestart}
! 
! !INTERFACE: 
  subroutine DEMCz_readrestart
! !USES: 
    use LIS_historyMod,      only : LIS_readvar_restart
! 
! !DESCRIPTION: 
! 
!   This routine reads the checkpoint data for a DEMCz restart
!EOP
    integer             :: n 
    integer             :: i, t, m
    integer             :: status
    real, pointer       :: vardata(:) 
    real, pointer       :: vardata1(:)
    character*100       :: vname_in_file
    character*100       :: vnames(demcz_ctl%nparam)
    type(ESMF_Field)    :: varField(demcz_ctl%nparam)
    integer             :: junk

    if(demcz_ctl%restart.eq."restart") then 
       n = 1    
       allocate(vardata(LIS_rc%ntiles(n)))
       
       call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
            rc=status)
       call LIS_verify(status)
       
       write(LIS_logunit,*) 'Reading the DEMCz restart file ..', trim(demcz_ctl%rfile)
       
       open(40,file=demcz_ctl%rfile,form='unformatted')        
       
       read(40) demcz_ctl%iterNo
       write(LIS_logunit,*) 'Iteration Number ',demcz_ctl%iterNo
       
       do i=1,demcz_ctl%nparam
          
          read(40) vname_in_file

          call LIS_readvar_restart(40,n,vardata)       
          
          do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
             do m=1,LIS_rc%nensem(n)
                demcz_struc(t)%X(i,m) = vardata((t-1)*LIS_rc%nensem(n)+m) 
             enddo
          enddo
       enddo

       !fitness
       call LIS_readvar_restart(40,n,vardata)            
       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
          do m=1,LIS_rc%nensem(n)
             demcz_struc(t)%LnZ(m) = vardata((t-1)*LIS_rc%nensem(n)+m) 
          enddo
       enddo

       close(40)
       write(LIS_logunit,*) 'Finished reading the DEMCz restart file ..'

       deallocate(vardata)
    elseif(demcz_ctl%restart.eq."coldstart") then
       if(demcz_ctl%optrfile.ne.'') then !start from optimization restart
          n = 1    
          write(LIS_logunit,*) 'Reading the DEMCz GA restart file ..'
          allocate(vardata(LIS_rc%ntiles(n)))

          call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
               rc=status)
          call LIS_verify(status)
          
          open(40,file=demcz_ctl%optrfile,form='unformatted')        
          
          read(40) junk !ga_ctl%genNo; KEEP iterNo=0
          write(LIS_logunit,*) 'Generation Number ',demcz_ctl%iterNo
          
          !read in, as with write_restart, in order of demcz_ctl (and ga_ctl) parameters list
          !which should be same as ESMF decision space state       
          do i=1,demcz_ctl%nparam
             call LIS_readvar_restart(40,n,vardata)       
             
             do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                do m=1,LIS_rc%nensem(n)
                   demcz_struc(t)%X(i,m) = vardata((t-1)*LIS_rc%nensem(n)+m) 
                   !ga_struc(t)%parent(i,m) = vardata((t-1)*LIS_rc%nensem(n)+m)
                   ! make sure that the binary representation is consistent
                   !call encode(m,i,ga_struc(t)%parent, ga_struc(t)%iparent)
                   !                call decode(m,tempdata, ga_struc(t)%iparent) !******temp******kwh
                enddo
             enddo
          enddo
          close(40)
          write(LIS_logunit,*) 'Finished reading the DEMCz GA restart file ..'
          
       endif
    endif

  end subroutine DEMCz_readrestart

!BOP
! 
! !ROUTINE: DEMCz_getNparam
! \label{DEMCz_getNparam}
! 
! !INTERFACE: 
    subroutine DEMCz_getNparam(nparam)
! 
! !DESCRIPTION: 
!  This method returns the number of decision space variables 
!
!EOP      
      integer   :: nparam
      
      nparam = demcz_ctl%nparam

    end subroutine DEMCz_getNparam

!BOP
! !ROUTINE: DEMCz_getdecSpaceValues
! \label{DEMCz_getdecSpaceValues}
! 
! !INTERFACE: 
    subroutine DEMCz_getdecSpaceValues(n)
! !USES: 

! 
! !DESCRIPTION: 
!  This routine returns the array of decision space variables from 
!  the GA data structures
!EOP
      implicit none

      integer            :: n
      type(ESMF_Field)   :: varField(demcz_ctl%nparam)
      real, pointer      :: vardata(:)
      integer            :: i, t, m
      integer            :: status
      
      do i=1,demcz_ctl%nparam
         call ESMF_StateGet(LIS_decisionSpace, trim(demcz_ctl%vname(i)), &
              varField(i), rc=status)
         call LIS_verify(status)
         
         call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata,rc=status)
         call LIS_verify(status)
         
         if(demcz_ctl%restart.eq."restart") then 
            do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
               do m=1,LIS_rc%nensem(n)
                  vardata((t-1)*LIS_rc%nensem(n)+m)= demcz_struc(t)%X(i,m)
               enddo
            enddo
         elseif(demcz_ctl%restart.eq."coldstart") then
            if(demcz_ctl%zerothRun) then 
               do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                  do m=1,LIS_rc%nensem(n)
                     vardata((t-1)*LIS_rc%nensem(n)+m) = demcz_struc(t)%X(i,m)
                  enddo
               enddo
            else
               do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                  do m=1,LIS_rc%nensem(n)
                     vardata((t-1)*LIS_rc%nensem(n)+m) = demcz_struc(t)%X_cand(i,m)
                  enddo
               enddo
            endif
         else
            write(LIS_logunit,*) 'Restart : ',demcz_ctl%restart, ' not valid'
            write(LIS_logunit,*) 'Must be either restart or coldstart'
            write(LIS_logunit,*) 'Stopping program ....'
            call LIS_endrun()
         endif
      enddo
      
      if(demcz_ctl%restart.eq."restart") then 
         demcz_ctl%restart = "coldstart"  !shouldn't be needed with zeroth run now
         demcz_ctl%zerothRun = .false. 
      endif

    end subroutine DEMCz_getdecSpaceValues


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

end module DEMCzAlgorithm
