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
! !ROUTINE: readpetusgscrd
! \label{readpetusgscrd}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 10 Mar 2012; K. Arsenault, modified for PET USGS
! 25 Oct 2013; K. Arsenault, added PET USGS to LIS7
! 02 May 2014; K. Arsenault, added climatological PET file option
!
! !INTERFACE:    
subroutine readpetusgscrd()

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_config, LIS_rc
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
  use petusgs_forcingMod, only : petusgs_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to USGS PET forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"USGS PET forcing directory:",rc=rc)
  call LIS_verify(rc, 'USGS PET forcing directory: not defined')
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,petusgs_struc(n)%petdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"USGS PET forcing type:",rc=rc)
  call LIS_verify(rc, 'USGS PET forcing type: not defined')
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,petusgs_struc(n)%pettype,rc=rc)
     if( trim(petusgs_struc(n)%pettype) .ne. "current" .and. &
         trim(petusgs_struc(n)%pettype) .ne. "climatology" ) then
        write(LIS_logunit,*) "ERR: This USGS PET forcing type is not supported: ",&
                              petusgs_struc(n)%pettype
        write(LIS_logunit,*) " ... Please select: 'current' or 'climatology'."
        write(LIS_logunit,*) " LIS end run being called ... "
        call LIS_endrun
     endif
  enddo

  do n=1, LIS_rc%nnest
     
    write(LIS_logunit,*) "Using USGS PET forcing"
    write(LIS_logunit,*) "USGS PET forcing directory :: ",trim(petusgs_struc(n)%petdir)
    write(LIS_logunit,*) "USGS PET forcing file type :: ",petusgs_struc(n)%pettype
    write(LIS_logunit,*) " "

!------------------------------------------------------------------------
! Setting global observed PET times to zero to ensure 
! data is read in during first time step
!------------------------------------------------------------------------
    petusgs_struc(n)%pettime = 0.0

  enddo

end subroutine readpetusgscrd

