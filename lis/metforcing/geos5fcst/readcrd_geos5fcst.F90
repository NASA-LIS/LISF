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
! !ROUTINE: readcrd_geos5fcst
!  \label{readcrd_geos5fcst}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
!
!
! !INTERFACE:    
subroutine readcrd_geos5fcst()
! !USES:
  use ESMF
  use LIS_coreMod,          only : LIS_rc, LIS_config
  use LIS_logMod,           only : LIS_logunit, LIS_verify
  use geos5fcst_forcingMod, only : geos5fcst_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to GEOS5 forecast forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n,rc

  
  call ESMF_ConfigFindLabel(LIS_config,"GEOS5 forecast forcing directory:",rc=rc)
  call LIS_verify(rc, 'GEOS5 forecast forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,geos5fcst_struc(n)%geos5fcstdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEOS5 forecast forcing number of ensemble members:",rc=rc)
  call LIS_verify(rc, 'GEOS5 forecast forcing number of ensemble members: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,geos5fcst_struc(n)%max_ens_members,rc=rc)
  enddo

  write(unit=LIS_logunit,fmt=*)'[INFO] Using GEOS5 forecast forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] GEOS5 forecast forcing directory :',&
          trim(geos5fcst_struc(n)%geos5fcstdir)

     geos5fcst_struc(n)%fcsttime1 = 3000.0
     geos5fcst_struc(n)%fcsttime2 = 0.0
  enddo

end subroutine readcrd_geos5fcst
