!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module MCSIMAlgorithm
!BOP
!
! !MODULE: MCSIMAlgorithm
!
! !DESCRIPTION: This module contains routines that define the operations
! of the MCSIM algorithm.
!  
! !REVISION HISTORY:
!  08 July 2010; Sujay Kumar, Ken Harrison; Initial Specification
!
! !USES:
  use ESMF
  use MCSIM_varctl
  use PobjFunc_Mod 
  use LIS_coreMod
  use LIS_optUEMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: MCSIM_init
  public :: MCSIM_setup
  public :: MCSIM_run
  public :: MCSIM_checkConvergence
  public :: MCSIM_getdecSpaceValues
  public :: MCSIM_getNparam
  public :: MCSIM_readrestart
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
  type mcsimstruc
     real,    allocatable :: pre_fitness(:)
     real,    allocatable :: pre_soln(:,:)
  end type mcsimstruc

  type(mcsimctl)            :: mcsim_ctl
  type(mcsimstruc), allocatable :: mcsim_struc(:)

  contains
   
!BOP
! !ROUTINE: MCSIM_init
! \label{MCSIM_init}
! 
! !INTERFACE: 
    subroutine MCSIM_init()
! !USES: 

! !DESCRIPTION: This routine performs the initialization steps for MCSIM.  
!   It initializes the required memory structures, and 
!   creates an initial random population, for each grid point. 

!   
!
!EOP      
      implicit none 
      
      integer                :: n 
      integer                :: status

!!$   type(diststruc) :: d
!!$   real, allocatable :: X(:)

      n = 1

      mcsim_ctl%iterNo = 0 

!      call ESMF_ConfigGetAttribute(LIS_config,mcsim_ctl%decspaceAttribsFile,&
!           label="MCSIM decision space attributes file:",rc=status)
!      call LIS_verify(status, 'MCSIM decision space attributes file: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,mcsim_ctl%maxiter,&
           label="MCSIM number of iterations:",rc=status)
      call LIS_verify(status, 'MCSIM number of iterations: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,mcsim_ctl%rfile,&
           label="MCSIM restart file:",rc=status)
      call LIS_verify(status, 'MCSIM restart file: not defined')

      call ESMF_ConfigGetAttribute(LIS_config,mcsim_ctl%restart,&
           label="MCSIM start mode:",rc=status)
      call LIS_verify(status, 'MCSIM start mode: not defined')

!      mcsim_ctl%maxiter = 1  !removing this so as to have ability 
!      to run without enormous ensemble sizes, and to conform to output of other algorithms
      call MCSIM_setup()

   end subroutine MCSIM_init

!BOP
! !ROUTINE: MCSIM_setup
! \label{MCSIM_setup}
!
! !INTERFACE: MCSIM_setup
    subroutine MCSIM_setup()
! !USES: 
! 
! !DESCRIPTION: 
!   This subroutine performs the second part of the MCMC initialization, 
!   This routine is currently empty. 
!  
!EOP

      type(ESMF_Field)           :: varField
      real,   pointer            :: vardata(:)
      integer                    :: t, k, n, m, j,kk
      integer                    :: status
      type(diststruc) :: d
      character*100            :: vname
      real, allocatable :: X(:)
      logical :: found

      n = 1
      
      call ESMF_StateGet(LIS_decisionSpace, itemCount=mcsim_ctl%nparam, &
           rc=status)
      call LIS_verify(status)

      allocate(mcsim_ctl%parmax(mcsim_ctl%nparam))
      allocate(mcsim_ctl%parmin(mcsim_ctl%nparam))
      allocate(mcsim_ctl%vname(mcsim_ctl%nparam))

      allocate(mcsim_struc(LIS_rc%ntiles(n)/LIS_rc%nensem(n)))

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         allocate(mcsim_struc(t)%pre_soln(mcsim_ctl%nparam, LIS_rc%nensem(n)))
         allocate(mcsim_struc(t)%pre_fitness(LIS_rc%nensem(n)))
      enddo

      call ESMF_StateGet(LIS_decisionSpace, itemNameList=mcsim_ctl%vname,&
           rc=status)
      call LIS_verify(status)
      
      do k=1,mcsim_ctl%nparam
         call ESMF_StateGet(LIS_decisionSpace, trim(mcsim_ctl%vname(k)), &
              varField, rc=status)
         call LIS_verify(status)
         
         call ESMF_AttributeGet(varField, 'MinRange',mcsim_ctl%parmin(k),&
              rc=status)
         call LIS_verify(status, 'setting minrange to decspace obj in MCSIM')
         call ESMF_AttributeGet(varField, 'MaxRange',mcsim_ctl%parmax(k),&
              rc=status)
         call LIS_verify(status, 'setting maxrange to decspace obj in MCSIM')

         call ESMF_FieldGet(varField,localDE=0, farrayPtr=vardata,rc=status)
         call LIS_verify(status)
         do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            do m=1,LIS_rc%nensem(n)
               mcsim_struc(t)%pre_soln(k,m) = vardata((t-1)*LIS_rc%nensem(n)+m)
            enddo
         enddo
      enddo

!sample and initialize

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            do j=1,dist_ctl%ndists
               d=dist_struc(j)
               allocate(X(d%numparam))
               call dist_sample(d,X)
               do k=1,d%numparam
!                  i=d%param(k)
                  vname = d%param_name(k)
                  found=.false.
                  do kk=1,mcsim_ctl%nparam
                     if(trim(vname).eq.trim(mcsim_ctl%vname(kk))) then
                        found=.true.
                        mcsim_struc(t)%pre_soln(kk,m) = X(k)
                        exit;
                     endif
                  enddo
                  if (.not.found) then
                     write(LIS_logunit,*) 'Prior for parameter ',trim(vname),' not found'
                     write(LIS_logunit,*) 'Program stopping ...'
                     call LIS_endrun
                  endif
               enddo
               deallocate(X)
            enddo
         enddo
      enddo
      call writeMCSIMdata()
     
    end subroutine MCSIM_setup

!BOP
! 
! !ROUTINE: MCSIM_checkConvergence
! \label{MCSIM_checkConvergence}
! 
! !INTERFACE: 
    subroutine MCSIM_checkConvergence(check)
! !ARGUMENTS: 
      logical, intent(INOUT) :: check
! 
! !DESCRIPTION: 
!  This routine checks to see if the convergence criteria for GA 
!  is reached. In this case, the routine simply checks to see if the
!  specified number of generations is reached. 
!EOP
      if(mcsim_ctl%iterNo.ge.mcsim_ctl%maxiter) then 
         check = .true.
      else
         check = .false.
      endif

    end subroutine MCSIM_checkConvergence

!BOP
! !ROUTINE: MCSIM_run
! \label{MCSIM_run}
! 
! !INTERFACE: 
    subroutine MCSIM_run()
! !USES: 
      use LIS_coreMod,   only : LIS_rc
      use LIS_optUEMod,  only : LIS_ObjectiveFunc, LIS_feasibleSpace
      use LIS_logMod,    only : LIS_logunit, LIS_verify

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
      integer                  :: t, m, j, k, i,kk
      integer                  :: status
      type(ESMF_Field)         :: fitField
      type(ESMF_Field)         :: feasField
      logical                  :: bounds_violation
      real,   pointer          :: fitValue(:)
      integer, pointer         :: mod_flag(:)
      type(diststruc) :: d
      character*100            :: vname
      real, allocatable :: X(:)

      n = 1
      write(LIS_logunit, *)  'MCSIM iteration attempt # ',mcsim_ctl%iterNo+1

! Evaluate those solutions that can be evaluated
!      call setoptuetypedecspace(LIS_rc%optuetype)
!      call evaluateobjfunction(LIS_rc%optuetype)
      
      call ESMF_StateGet(LIS_ObjectiveFunc,"Max Criteria Value",fitField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(fitField,localDE=0, farrayPtr=fitValue,rc=status)
      call LIS_verify(status)
      
      call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
      call LIS_verify(status)

      bounds_violation=.false.

! Set the fitness values for the solutions in the current LIS run (includes candidates and unchanged)
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            mcsim_struc(t)%pre_fitness(m) = fitValue((t-1)*LIS_rc%nensem(n)+m)
            !penalize infeasible solutions
            if(mod_flag((t-1)*LIS_rc%nensem(n)+m).eq.1) then 
               bounds_violation = .true.
            endif
         end do
      enddo

! Double check parameter bounds; Actually needed for the error variable (sigma) that is
! not currently thought of as a parameter of the LSM
      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            do j=1,mcsim_ctl%nparam
               if((mcsim_struc(t)%pre_soln(j,m) .lt. mcsim_ctl%parmin(j)) &
                    .or.(mcsim_struc(t)%pre_soln(j,m) .gt. &
                    mcsim_ctl%parmax(j))) then
                  bounds_violation=.true.
               endif
            enddo
         enddo
      enddo

      do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            do j=1,dist_ctl%ndists
               d=dist_struc(j)
               allocate(X(d%numparam))
               call dist_sample(d,X)
               do k=1,d%numparam
!                  i=d%param(k)
                  vname = d%param_name(k)
                  do kk=1,mcsim_ctl%nparam
                     if(trim(vname).eq.trim(mcsim_ctl%vname(kk))) then 
                        mcsim_struc(t)%pre_soln(kk,m) = X(k)
                        exit;
                     endif
                  enddo
               enddo
               deallocate(X)
            enddo
         enddo
      enddo

      if (.not.bounds_violation) then
         !advance iteration
         mcsim_ctl%iterNo = mcsim_ctl%iterNo +1
!         call writeMCSIMdata()
      else
         !redo full iteration
         mcsim_ctl%iterNo = mcsim_ctl%iterNo +0
      end if
      call writeMCSIMdata()

    end subroutine MCSIM_run

  subroutine writeMCSIMdata()
    call writeMCSIMrestart()
  end subroutine writeMCSIMdata

!BOP
! 
! !ROUTINE: writeMCSIMrestart
! \label{writeMCSIMrestart}
! 
! !INTERFACE: 
  subroutine writeMCSIMrestart
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_masterproc
    use LIS_fileIOMod,  only : LIS_create_output_directory
    use LIS_logMod,     only : LIS_logunit, LIS_verify
    use LIS_historyMod, only : LIS_writevar_restart
! 
! !DESCRIPTION: 
! 
! This routine writes the checkpoint data for a MCMC restart.
!Actually, the restart file is not necessary but instead serves
!as an output capability in the same format as the other optUE algorithms
! 
!EOP

    integer             :: n 
    integer             :: i,m,t
    integer             :: status
    character(len=LIS_CONST_PATH_LEN) :: filen
    character (len=4)   :: fgen
    character*100       :: vnames(mcsim_ctl%nparam)
    real, allocatable       :: vardata(:)

    n = 1
    allocate(vardata(LIS_rc%ntiles(n)))

    if(LIS_masterproc) then 
       call LIS_create_output_directory('MCSIM')
       write(unit=fgen, fmt='(i4.4)') mcsim_ctl%iterNo
       filen = trim(LIS_rc%odir)//'/MCSIM/MCSIM.'&
            //trim(fgen)//'.MCSIMrst'
       open(40,file=filen,status='unknown',form='unformatted')
       write(40) mcsim_ctl%iterNo
    endif

    do i=1,mcsim_ctl%nparam
       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
          do m=1,LIS_rc%nensem(n)
             vardata((t-1)*LIS_rc%nensem(n)+m) = mcsim_struc(t)%pre_soln(i,m)
          enddo
       enddo
       if(LIS_masterproc) then 
          write(40) mcsim_ctl%vname(i)
       endif
       call LIS_writevar_restart(40,n,vardata)
    enddo

   do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)       
       do m=1,LIS_rc%nensem(n)
          vardata((t-1)*LIS_rc%nensem(n)+m) = mcsim_struc(t)%pre_fitness(m)
       enddo
    enddo
    call LIS_writevar_restart(40,n,vardata)
       
    if(LIS_masterproc) then 
       close(40)
       write(LIS_logunit,*) 'MCSIM checkpoint file written ',trim(filen)
    endif
    deallocate(vardata)

  end subroutine writeMCSIMrestart


!BOP
! 
! !ROUTINE: MCSIM_readrestart
! \label{MCSIM_readrestart}
! 
! !INTERFACE: 
  subroutine MCSIM_readrestart
! !USES: 
    use LIS_historyMod,      only : LIS_readvar_restart

! 
! !DESCRIPTION: 
! 
!   This routine reads the checkpoint data for a MCSIM restart
!EOP
    integer             :: n 
    integer             :: i, t, m
    integer             :: status
    real, allocatable       :: vardata(:) 
    real, allocatable       :: vardata1(:)
    character*100       :: vnames(mcsim_ctl%nparam)
    type(ESMF_Field)    :: varField(mcsim_ctl%nparam)

    if(mcsim_ctl%restart.eq."restart") then 
       n = 1    
       write(LIS_logunit,*) 'Reading the MCSIM restart file ..'
       allocate(vardata(LIS_rc%ntiles(n)))
       
       call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
            rc=status)
       call LIS_verify(status)
       
       open(40,file=mcsim_ctl%rfile,form='unformatted')        
       
       read(40) mcsim_ctl%iterNo
       write(LIS_logunit,*) 'Iteration Number ',mcsim_ctl%iterNo
       
       !If surface model run not complete, will need, in addition to iteration number, the data.
       !read in, as with write_restart, in order of mcsim_ctl parameters list
       !which should be same as ESMF decision space state       

       do i=1,mcsim_ctl%nparam
          read(40) mcsim_ctl%vname(i)

          call LIS_readvar_restart(40,n,vardata)       

          do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
             do m=1,LIS_rc%nensem(n)
                mcsim_struc(t)%pre_soln(i,m) = vardata((t-1)*LIS_rc%nensem(n)+m)
             enddo
          enddo
      enddo

      call LIS_readvar_restart(40,n,vardata)       
      do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
         do m=1,LIS_rc%nensem(n)
            mcsim_struc(t)%pre_fitness(m) = vardata((t-1)*LIS_rc%nensem(n)+m)
         enddo
      enddo

      close(40)
      write(LIS_logunit,*) 'Finished reading the MCSIM restart file ..'
       
! update the MCSIM objects. 

       deallocate(vardata)
    endif

  end subroutine MCSIM_readrestart

!BOP
! 
! !ROUTINE: MCSIM_getNparam
! \label{MCSIM_getNparam}
! 
! !INTERFACE: 
    subroutine MCSIM_getNparam(nparam)
! 
! !DESCRIPTION: 
!  This method returns the number of decision space variables 
!
!EOP      
      integer   :: nparam
      
      nparam = mcsim_ctl%nparam

    end subroutine MCSIM_getNparam

!BOP
! !ROUTINE: MCSIM_getdecSpaceValues
! \label{MCSIM_getdecSpaceValues}
! 
! !INTERFACE: 
    subroutine MCSIM_getdecSpaceValues(n)
! !USES: 
      use LIS_coreMod,    only : LIS_rc
! 
! !DESCRIPTION: 
!  This routine returns the array of decision space variables from 
!  the GA data structures
!EOP
      implicit none

      integer            :: n
      type(ESMF_Field)   :: varField(mcsim_ctl%nparam)
      real, pointer      :: vardata(:)
      integer            :: i, t, m
      integer            :: status
      
      do i=1,mcsim_ctl%nparam
         call ESMF_StateGet(LIS_decisionSpace, trim(mcsim_ctl%vname(i)), &
              varField(i), rc=status)
         call LIS_verify(status)
         
         call ESMF_FieldGet(varField(i),localDE=0, farrayPtr=vardata,rc=status)
         call LIS_verify(status)
         
         do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
            do m=1,LIS_rc%nensem(n)
               vardata((t-1)*LIS_rc%nensem(n)+m)= &
                    mcsim_struc(t)%pre_soln(i,m)
            enddo
         enddo         
      enddo

      if(mcsim_ctl%restart.eq."restart") then 
         mcsim_ctl%restart = "coldstart"
      endif
    end subroutine MCSIM_getdecSpaceValues


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

end module MCSIMAlgorithm
