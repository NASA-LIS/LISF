!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
module LIS_perturbMod

!BOP
!
! !MODULE: LIS_perturbMod
!
! !DESCRIPTION:
!  The code in this file provides abstractions to manage different 
!  perturbation implementations. This includes methods to perturb
!  forcing, observations, and land surface states. 
!  
! !REVISION HISTORY:
!
!  25 Jun 2006: Sujay Kumar; Initial implementation
!  24 Nov 2010: Sujay Kumar; Added support for time varying perturbations
!   9 Jun 2020: David Mocko; Write restart file at end of the LIS run
!
!EOP
  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_perturb_init    ! initialization for perturbation routines
  public :: LIS_perturb_readrestart ! read perturbation restart files
  public :: LIS_perturb_writerestart ! write perturbation restart files
  public :: LIS_readPertAttributes !read perturbation attributes
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------

!EOP

  type, public :: pert_dec_type
     character*40,allocatable :: vname(:)     !name of the observation variables
     integer     ,allocatable :: perttype(:)  !type of perturbation (0-additive, 1-multiplicative)
     real        ,allocatable :: ssdev(:) !standard deviation of perturbation
     real        ,allocatable :: stdmax(:) !standard normal max for the perturbation distribution
     integer     ,allocatable :: zeromean(:) !flag to specify if zero mean is to be ensured
     real        ,allocatable :: tcorr(:) !temporal correlation scale (seconds)
     real        ,allocatable :: xcorr(:) !horizontal x-correlation scale
     real        ,allocatable :: ycorr(:) !horizontal y-correlation scale
     real        ,allocatable :: ccorr(:,:) !cross correlations (specified as a matrix)
  end type pert_dec_type 

contains
 
!BOP
! 
! !ROUTINE: LIS_perturb_init
! \label{LIS_perturb_init}
! 
! !INTERFACE:  
  subroutine LIS_perturb_init
! 
! !DESCRIPTION: 
! 
! The routines invoked are: 
! \begin{description}
!  \end{description}
! 
! !USES: 
    use LIS_logMod,             only : LIS_verify
    use LIS_timeMgrMod,         only : LIS_registerAlarm
!EOP
    integer                   :: rc

    TRACE_ENTER("perturb_init")
!    if(LIS_rc%nperts.gt.0) then 
!       call LIS_registerAlarm("LIS pert restart alarm",&
!            900.0 ,&  ! 1800 must be modified with the min perturbation timestep
!            LIS_rc%pertrestartInterval)
!    endif
    TRACE_EXIT("perturb_init")

  end subroutine LIS_perturb_init
!BOP
! 
! !ROUTINE: LIS_perturb_readrestart
! \label{LIS_perturb_readrestart}
! 
! !INTERFACE: 
  subroutine LIS_perturb_readrestart
! !USES: 
! !ARGUMENTS:     

! 
! !DESCRIPTION: 
! 
!EOP
    integer      :: k
    logical      :: rstFlag
    character*50 :: pertname

    TRACE_ENTER("perturb_readrst")
    rstFlag = .false. 
    pertname = "none"
    if(LIS_rc%pertrestart.eq."restart".and.LIS_rc%nperts.gt.0) then 
       rstFlag = .true. 
    endif
!
! In the smoother mode, the run may start with coldstart, but the perturbation
! restart files should be read after the initial iteration. 
!
! Currently we assume that all the nests are in sync, so only checking for
! the first nest iteration id
! 
    if(LIS_rc%runmode.eq."ensemble smoother".and.&
         LIS_rc%pertrestart.eq."coldstart") then 
       if(LIS_rc%iterationId(1).gt.1) then 
          rstFlag = .true. 
       endif
    endif
    
! The following is to take care of instances when only certain set of 
! perturbations are enabled. It is implictly assumed that the same
! perturbation algorithm is used across forcing, state and observation
! perturbations

    if(LIS_rc%perturb_forcing.ne."none") then 
       pertname = trim(LIS_rc%perturb_forcing)
    else
       do k=1,LIS_rc%nperts
          if(LIS_rc%perturb_state(k).ne."none") then 
             pertname = trim(LIS_rc%perturb_state(k))
          endif
       enddo
       
    endif
    if(pertname.eq."none") then 
       do k=1,LIS_rc%nperts
          if(LIS_rc%perturb_obs(k).ne."none") then 
             pertname = trim(LIS_rc%perturb_obs(k))
          endif
       enddo
    endif

    if (rstFlag) then 
       call readpertrestart(trim(pertname)//char(0))
    endif
    TRACE_EXIT("perturb_readrst")

  end subroutine LIS_perturb_readrestart
!BOP
! 
! !ROUTINE: LIS_perturb_writerestart
! \label{LIS_perturb_writerestart}
! 
! !INTERFACE: 
  subroutine LIS_perturb_writerestart(n)
! !USES: 
! !ARGUMENTS:     
    integer, intent(IN) :: n 
! 
! !DESCRIPTION: 
! 
!EOP
    integer           :: rc
    logical           :: alarmCheck
    character*50      :: pertname
    integer           :: k

    TRACE_ENTER("perturb_writerst")
    if(LIS_rc%nperts.gt.0) then 
       
       pertname = "none"

       alarmCheck = LIS_isAlarmRinging(LIS_rc, "LIS pert restart alarm")

       if (alarmCheck.or.(LIS_rc%endtime==1)) then
!  NOTE: 
!  A single restart file will be written for forcing, state and observation
!  perturbations. Here we assume that the same perturbation algorithm will be 
!  used for all three. 
          pertname = "none"
          if(LIS_rc%perturb_forcing.ne."none") then 
             pertname = trim(LIS_rc%perturb_forcing)
          else
             do k=1,LIS_rc%nperts
                if(LIS_rc%perturb_state(k).ne."none") then 
                   pertname = trim(LIS_rc%perturb_state(k))
                endif
             enddo
             
          endif
          if(pertname.eq."none") then 
             do k=1,LIS_rc%nperts
                if(LIS_rc%perturb_obs(k).ne."none") then 
                   pertname = trim(LIS_rc%perturb_obs(k))
                endif
             enddo
          endif
          
          if(pertname.ne."none") then 
             call writepertrestart(trim(pertname)//char(0), n)
          endif
       endif
    endif
    TRACE_EXIT("perturb_writerst")

  end subroutine LIS_perturb_writerestart

!BOP
! !ROUTINE: LIS_readPertAttributes
! \label{LIS_readPertAttributes}
!
! !REVISION HISTORY:
!  14 Oct 2006: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine LIS_readPertAttributes(nvars, fname,pert_struc)
! !USES: 
    use LIS_logMod,    only : LIS_logunit, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber

    implicit none
! !ARGUMENTS: 
    integer             :: nvars
    character(len=*)    :: fname
    type(pert_dec_type) :: pert_struc
! 
! !DESCRIPTION: 
!  This routine reads the observation perturbation attributes for each
!  observation type. The format of the perturbation attributes file
!  is: 
!
!\begin{verbatim}
!  <Variable name>  
!  <perturbation type> <perturbation standard deviation type>
!  <standard deviation> <coeff standard deviation (optional, only if 
!  perturbation standard deviation type is > 0) <standard normal max> (contd.)
!  <ensure zeromean?> <temporal correlation scale> (contd.)
!  <x correlation scale>  <y correlation scale> <cross correlations> 
!\end{verbatim}
!
!  The arguments are: 
!  \begin{description}
!   \item[n]           index of nest
!   \item[vname]       name of the observation variables
!   \item[perttype]    type of perturbation (0-additive, 1-multiplicative)
!   \item[ssdev]       standard deviation of perturbation
!   \item[stdmax]      standard normal max for the perturbation distribution
!   \item[zeromean]    flag to specify if zero mean is to be ensured
!   \item[tcorr]       temporal correlation scale (seconds)
!   \item[xcorr]       horizontal x-correlation scale
!   \item[ycorr]       horizontal y-correlation scale
!   \item[ccorr]       cross correlations
!  \end{description}
!EOP
    integer            :: i,j
    integer            :: ftn 

    TRACE_ENTER("perturb_readAtt")
    write(LIS_logunit,*) '[INFO] Opening perturbation attributes for observations ',&
         trim(fname)
    ftn = LIS_getNextUnitNumber()
    open(ftn,file=fname,status='old')
    read(ftn,*)
    
    do i=1,nvars

       read(ftn,fmt='(a40)') pert_struc%vname(i)
       read(ftn,*) pert_struc%perttype(i) , & 
            pert_struc%ssdev(i), pert_struc%stdmax(i), &
            pert_struc%zeromean(i), &
            pert_struc%tcorr(i), pert_struc%xcorr(i), pert_struc%ycorr(i), &
            (pert_struc%ccorr(i,j),j=1,nvars)
       
       write(LIS_logunit,*) '[INFO] ',pert_struc%vname(i), pert_struc%perttype(i), &
            pert_struc%ssdev(i),pert_struc%stdmax(i),&
            pert_struc%zeromean(i), &
            pert_struc%tcorr(i), pert_struc%xcorr(i),pert_struc%ycorr(i), &
            (pert_struc%ccorr(i,j),j=1,nvars)
              
    enddo

    write(LIS_logunit,*) &
         '[INFO] Successfully read perturbation attributes file ', trim(fname)
    call LIS_releaseUnitNumber(ftn)
    TRACE_EXIT("perturb_readAtt")

  end subroutine LIS_readPertAttributes

end module LIS_perturbMod
