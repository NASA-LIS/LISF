!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: template_obsMod
! \label(template_obsMod)
!
! !INTERFACE:
module template_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This is an empty template for defining new observational plugins into
!  LVT
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 Apr 2009   Sujay Kumar  Initial Specification
! 
!EOP

  PUBLIC :: template_obsinit
  PUBLIC :: templateobs

  type, public :: templateobsdec
     character*100 :: odir
  end type templateobsdec

  type(templateobsdec) :: templateobs

contains
  
!BOP
! 
! !ROUTINE: template_obsInit
! \label{template_obsInit}
!
! !INTERFACE:
  subroutine template_obsinit(i)
! 
! !USES:   
    use LVT_coreMod, only : LVT_rc, LVT_config
    use LVT_logMod,  only : LVT_verify

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine should contain code to initialize and set up the 
!  data structures required for reading the specific data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer :: status
!------------------------------------------------------------------------------------
! Read any runtime specifications from the lvt.config file. 
!------------------------------------------------------------------------------------

!    call ESMF_ConfigGetAttribute(LVT_config, templateobs%odir, &
!         'Template observation directory:',rc=status)
!    call LVT_verify(status, 'Template observation directory: not defined')
  
!------------------------------------------------------------------------------------
! Initialize any other variables.   
!------------------------------------------------------------------------------------

  end subroutine Template_obsinit


end module Template_obsMod
