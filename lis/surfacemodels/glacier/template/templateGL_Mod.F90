!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
module templateGL_Mod
!BOP
!
! !MODULE: templateGL_Mod
!
! !DESCRIPTION:
!
! !USES:
  
  implicit none
  
  PRIVATE
  !-------------------------------------------------------------------------
  ! PUBLIC MEMBER FUNCTIONS
  !-------------------------------------------------------------------------
  public :: templateGL_ini
  !-------------------------------------------------------------------------
  ! PUBLIC TYPES
  !-------------------------------------------------------------------------
    public :: templateGL_struc
!EOP
    type, public :: templateGL_type_dec
       real               :: ts
    end type templateGL_type_dec
    
    type(templateGL_type_dec), pointer :: templateGL_struc(:)
 
contains 

!BOP
!
! !ROUTINE: templateGL_ini
! \label{templateGL_ini}
!
! !INTERFACE:
  subroutine templateGL_ini()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use LIS_logMod, only : LIS_verify
    use LIS_timeMgrMod, only : LIS_clock,  LIS_calendar, &
         LIS_update_timestep, LIS_registerAlarm
    use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for templateGL-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for templateGL from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[templateGL\_readcrd](\ref{templateGL_readcrd}) \newline
!    reads the runtime options for templateGL model
!  \end{description}
!EOP
    implicit none        
    integer  :: n, t     
    integer  :: status   

    allocate(templateGL_struc(LIS_rc%nnest))

    call templateGL_readcrd()
    
    do n=1,LIS_rc%nnest
       call LIS_update_timestep(LIS_rc, n, templateGL_struc(n)%ts)
       
       call LIS_registerAlarm("templateGL model alarm",&
            templateGL_struc(n)%ts, &
            templateGL_struc(n)%ts)
       
       LIS_sfmodel_struc(n)%ts = templateGL_struc(n)%ts
    enddo

  end subroutine templateGL_ini
end module templateGL_Mod
