!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_ecmwfreanal
!  \label{readcrd_ecmwfreanal}
!
! !REVISION HISTORY:
! 26Jan2004; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_ecmwfreanal()
! !USES:
  use ESMF
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc, LIS_config
  use ecmwfreanal_forcingMod, only : ecmwfreanal_struc
! !DESCRIPTION:
!
!  This routine reads the options specific to ECMWF reanalysis forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"ECMWF Reanalysis forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwfreanal_struc(n)%ecmwfreanaldir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"ECMWF Reanalysis maskfile:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwfreanal_struc(n)%emaskfile,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF Reanalysis elevation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwfreanal_struc(n)%elevfile,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF Reanalysis domain x-dimension size:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwfreanal_struc(n)%ncold,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF Reanalysis domain y-dimension size:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwfreanal_struc(n)%nrold,rc=rc)
     write(LIS_logunit,*)'Using ECMWF Reanalysis forcing'
     ecmwfreanal_struc(n)%fmodeltime1 = 3000.0
     ecmwfreanal_struc(n)%fmodeltime2 = 0.0
  enddo

end subroutine readcrd_ecmwfreanal
