!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_routingMod
!BOP
! 
! !MODULE: LIS_routingMod
! 
! !DESCRIPTION: 
!  The code in this file provides the top level calls to manage the 
!  operation of different runoff-routing algorithms and models. 
! 
! !REVISION HISTORY: 
!  6 May 2011: Sujay Kumar, Initial implementation
! 
! !USES: 
  use ESMF

  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_routing_init         ! initialize the routing model
  public :: LIS_routing_readrestart  ! read the routing model restart file
  public :: LIS_routing_run           ! execute the routing model 
  public :: LIS_routing_writeoutput   ! write the output file
  public :: LIS_routing_writerestart   ! write the restart file
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LIS_runoff_state
!EOP
  
  type(ESMF_State), allocatable :: LIS_runoff_state(:)

  contains
!BOP
! !ROUTINE: LIS_routing_init
! \label{LIS_routing_init}
! 
! !INTERFACE: 
  subroutine LIS_routing_init
! !USES: 
    use LIS_coreMod,  only : LIS_rc, LIS_Config
    use LIS_logMod,   only : LIS_verify
! 
! !DESCRIPTION: 
!EOP
    integer           :: rc
    integer           :: n, status
    character*1       :: nestid(2)
    character*100     :: temp

    call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%routingmodel, &
         label="Routing model:",default="none", rc=rc)

    if(LIS_rc%routingmodel.ne."none") then 

       allocate(LIS_runoff_state(LIS_rc%nnest))

       do n=1, LIS_rc%nnest

          write(unit=temp,fmt='(i2.2)') n
          read(unit=temp,fmt='(2a1)') nestid

          LIS_runoff_state(n) = ESMF_StateCreate(name="LIS Runoff State"//&
               nestid(1)//nestid(2), rc=status)
          call LIS_verify(status, "ESMF_StateCreate failed in LIS_routing_init")
               
       enddo

       if(LIS_rc%lsm.eq."none") then 
          call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%runoffdatasource, &
               label="External runoff data source:", rc=rc)
          call LIS_verify(rc,"External runoff data source: not defined")
       endif

       call routinginit(trim(LIS_rc%routingmodel)//char(0))
    endif



  end subroutine LIS_routing_init


!BOP
! !ROUTINE: LIS_routing_readrestart
! \label{LIS_routing_readrestart}
! 
! !INTERFACE: 
  subroutine LIS_routing_readrestart
! !USES: 
    use LIS_coreMod,  only : LIS_rc
! 
! !DESCRIPTION: 
!EOP
    if(LIS_rc%routingmodel.ne."none") then 
       call routingreadrestart(trim(LIS_rc%routingmodel)//char(0))
    endif
  
  end subroutine LIS_routing_readrestart

!BOP
! !ROUTINE: LIS_routing_run
! \label{LIS_routing_run}
! 
! !INTERFACE: 
  subroutine LIS_routing_run(n)
! !USES: 
    use LIS_coreMod
! 
! !DESCRIPTION: 
!EOP
    integer, intent(in) :: n

    if(LIS_rc%routingmodel.ne."none") then 
       call lsmroutinggetrunoff(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%routingmodel)//char(0),n)
       if(LIS_rc%glaciermodel.ne."none") then
          call glacierroutinggetrunoff(trim(LIS_rc%glaciermodel)//"+"//&
               trim(LIS_rc%routingmodel)//char(0),n)
       endif

       call routingrun(trim(LIS_rc%routingmodel)//char(0),n)
    endif

  end subroutine LIS_routing_run

!BOP
! !ROUTINE: LIS_routing_writeoutput
! \label{LIS_routing_writeoutput}
! 
! !INTERFACE: 
  subroutine LIS_routing_writeoutput(n)
! !USES: 
    use LIS_coreMod,  only : LIS_rc
! 
! !DESCRIPTION: 
!EOP

    integer, intent(in) :: n
  
    if(LIS_rc%routingmodel.ne."none") then 
       call routingoutput(trim(LIS_rc%routingmodel)//char(0),n)
    endif

  end subroutine LIS_routing_writeoutput


!BOP
! !ROUTINE: LIS_routing_writerestart
! \label{LIS_routing_writerestart}
! 
! !INTERFACE: 
  subroutine LIS_routing_writerestart(n)
! !USES: 
    use LIS_coreMod,  only : LIS_rc
! 
! !DESCRIPTION: 
!EOP

    integer, intent(in) :: n
  
    if(LIS_rc%wopt_rst.ne.0) then 
       if(LIS_rc%routingmodel.ne."none") then 
          call routingwriterestart(trim(LIS_rc%routingmodel)//char(0),n)
       endif
    endif

  end subroutine LIS_routing_writerestart



end module LIS_routingMod
