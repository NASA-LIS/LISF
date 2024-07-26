!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module metForcTemplate_forcingMod
!BOP
! !MODULE: metForcTemplate_forcingMod
! 
! !USES:
! none

  implicit none

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_metForcTemplate      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: metForcTemplate_struc
!EOP

  type metForcTemplate_type_dec
     integer      :: ncold, nrold   !AWIPS 212 dimensions

  end type metForcTemplate_type_dec
  
  type(metForcTemplate_type_dec), allocatable :: metForcTemplate_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_metForcTemplate
!  \label{init_metForcTemplate}
! !INTERFACE:
  subroutine init_metForcTemplate(findex)
! !USES: 
    use LIS_coreMod, only : LIS_rc
    use LIS_logMod,  only : LIS_logunit, LIS_endrun
!EOP

    implicit none
    integer,   intent(in)  :: findex

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the Template forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    LIS_rc%met_nf(findex) = 0

  end subroutine init_metForcTemplate
end module metForcTemplate_forcingMod


