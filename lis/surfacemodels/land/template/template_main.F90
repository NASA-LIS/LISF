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
#undef TBOT_TESTING
!BOP
!
! !ROUTINE: template_main
! \label{template_main}
!
! !ROUTINE: template_main.F90
!    Apr 2003; Sujay Kumar, Initial Code
! 23 Oct 2007; Kristi Arsenault, Updated code for LISv50
! 
! !INTERFACE:
subroutine template_main(n)

#if ( defined TBOT_TESTING )
  use LIS_coreMod,     only : LIS_rc, LIS_domain
  use LIS_histDataMod, only : LIS_diagnoseSurfaceOutputVar, & 
                              LIS_MOC_TEMPBOT
#endif
  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
! 
!  Calls the run routines for the forcing-only option (template)
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

#if ( defined TBOT_TESTING )
! Since tbot is not yet available as a parameter output option,
! make this routine read tbot at each time-step to enable writing tbot
! to output file for debugging.
  integer :: t
  real, allocatable, dimension(:,:) :: tbot
  real, allocatable, dimension(:)   :: tbot1

  allocate(tbot(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  allocate(tbot1(LIS_rc%npatch(n)))

  tbot  = LIS_rc%udef
  tbot1 = LIS_rc%udef

  call readtbot(n, LIS_rc%tbotsrc(n), tbot)

  do t = 1, LIS_rc%npatch(n)
     if ( tbot(LIS_domain(n)%tile(t)%col, &
               LIS_domain(n)%tile(t)%row) /= LIS_rc%udef ) then
        tbot1(t) = tbot(LIS_domain(n)%tile(t)%col,LIS_domain(n)%tile(t)%row)
     endif
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TEMPBOT, vlevel=1, &
                                value=tbot1(t),unit="K",direction="-",&
                                surface_type=LIS_rc%lsm_index)
  enddo

  deallocate(tbot)
  deallocate(tbot1)
#endif

end subroutine template_main
