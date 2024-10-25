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
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
!BOP
!
! !ROUTINE: lisdrv
!  \label{lisdrv}
!   Main program for LIS in an offline simulation
!  
! !DESCRIPTION: 
!  Main driver program for LIS. It performs four main functions
!
!  \begin{description}
!  \item[LIS\_config\_init]  
!      calls the routines to read the runtime configurations
!  \item[LIS\_Init]  
!      calls the initializations based on the runmode
!  \item[LIS\_run]
!      calls the run routines based on the runmode
!  \item[LIS\_finalize] 
!     calls the cleanup routines 
!  \end{description}
!
!  The routines invoked are :
!  \begin{description}
!   \item[LIS\_config\_init](\ref{LIS_config_init}) \newline
!    call to initialize configuration tool and read model independent
!    options
!   \item[lisinit](\ref{lisinit}) \newline
!    call to initialize lis based on the runmode
!   \item[lisrun](\ref{lisrun}) \newline
!    call to run lis based on the runmode
!   \item[lisfinalize](\ref{lisfinalize}) \newline
!    call to cleanup lis structures based on the runmode
!  \end{description}
! !REVISION HISTORY: 
!  14Nov02    Sujay Kumar  Initial Specification
!  21Oct05    Sujay Kumar  Modified to include the runmodes. Switched
!                          to a init,run,finalize mode
program lisdrv
! !USES:       
  use LIS_coreMod, only : LIS_config_init, LIS_rc
#ifdef ESMF_TRACE
  use ESMF
#endif
!EOP
  implicit none

!BOC  
  call LIS_config_init(cmd_args=.true.)
  TRACE_ENTER("LIS_init")
  call lisinit(trim(LIS_rc%runmode)//char(0))
  TRACE_EXIT("LIS_init")
  TRACE_ENTER("LIS_run")
  call lisrun(trim(LIS_rc%runmode)//char(0))
  TRACE_EXIT("LIS_run")
  call lisfinalize(trim(LIS_rc%runmode)//char(0))
!EOC
end program lisdrv


