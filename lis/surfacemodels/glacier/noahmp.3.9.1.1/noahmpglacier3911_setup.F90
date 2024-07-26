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
!BOP
!
! !ROUTINE: noahmpglacier3911_setup
! \label{noahmpglacier3911_setup}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine noahmpglacier3911_setup()
! !USES:
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    !use NOAHMP_VEG_PARAMETERS_36, only: read_mp_veg_parameters
    use MODULE_SF_NOAHMPLSM_36, only: read_mp_veg_parameters
    use noahmpglacier3911_Mod

!
! !DESCRIPTION:
!
!  This routine is the entry point to set up the parameters
!  required for noahmpglacier3911.  These include: 
!    vegetype     - land cover type index [-]
!    soiltype     - soil type index [-]
!    slopetype    - slope type for Noah baseflow [-]
!    tbot         - deep-layer soil temperature [K]
!    pblh         - planetary boundary layer height [m]
! 
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
!    retrieves LIS parameter data from NetCDF file
!  \item[NOAHMP36\_read\_MULTILEVEL\_param](\ref{NOAHMP36_read_MULTILEVEL_param}) \newline
!    retrieves MULTILEVEL spatial parameter from NetCDF file
!  \end{description}
!EOP
    implicit none
    integer           :: mtype
    integer           :: t, k, n
    integer           :: col, row
    real, allocatable :: placeholder(:,:)

    mtype = LIS_rc%glacier_index
    
    do n=1, LIS_rc%nnest

       allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
       ! read: tbot
       write(LIS_logunit,*) "[INFO] noahmpglacier3911: reading parameter TBOT from ", trim(LIS_rc%paramfile(n))
       call LIS_read_param(n, trim(Noahmpgl3911_struc(n)%LDT_ncvar_tbot), placeholder)
       do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          noahmpgl3911_struc(n)%noahmpgl(t)%tbot = placeholder(col, row)
       enddo
       deallocate(placeholder)
    enddo

  end subroutine noahmpglacier3911_setup
