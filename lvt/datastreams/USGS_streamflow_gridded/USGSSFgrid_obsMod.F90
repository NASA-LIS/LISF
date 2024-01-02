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
! !MODULE: USGSSFgrid_obsMod
! \label(USGSSFgrid_obsMod)
!
! !INTERFACE:
module USGSSFgrid_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 

  PUBLIC :: USGSSFgrid_obsinit
  PUBLIC :: USGSSFgridobs

  type, public :: USGSSFgridobsdec
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  13 May 2011   Sujay Kumar  Initial Specification
! 
!EOP

     character*100           :: odir
     integer                 :: nc,nr
     integer                 :: mo
     integer                 :: nts
     integer                 :: startFlag
     real,      allocatable  :: q(:,:,:)

  end type USGSSFgridobsdec

  type(USGSSFgridobsdec), allocatable :: USGSSFgridobs(:)

contains
  
!BOP
! 
! !ROUTINE: USGSSFgrid_obsInit
! \label{USGSSFgrid_obsInit}
!
! !INTERFACE: 
 subroutine USGSSFgrid_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    character*100 :: stnlist_file
    integer       :: ftn, k, status

    if(.not.allocated(USGSSFgridobs)) then 
       allocate(USGSSFgridobs(LVT_rc%nDataStreams))
    endif
!------------------------------------------------------------------------------
! Read any runtime specifications from the lvt.config file. 
!------------------------------------------------------------------------------

    write(LVT_logunit,*) '[INFO] Initializing USGS gridded streamflow data reader....'
    call ESMF_ConfigGetAttribute(LVT_config, USGSSFgridobs(i)%odir, &
         label='USGS gridded streamflow observation directory:',rc=status)
    call LVT_verify(status, 'USGS gridded streamflow observation directory: not defined')

    USGSSFgridobs(i)%nc = 464
    USGSSFgridobs(i)%nr = 224
    USGSSFgridobs(i)%nts = 417

    allocate(USGSSFgridobs(i)%q(&
         USGSSFgridobs(i)%nc,&
         USGSSFgridobs(i)%nr,&
         USGSSFgridobs(i)%nts))
    USGSSFgridobs(i)%q = -9999.0

    USGSSFgridobs(i)%startFlag = -1
    USGSSFgridobs(i)%mo        = LVT_rc%mo

  end subroutine USGSSFgrid_obsinit


end module USGSSFgrid_obsMod
