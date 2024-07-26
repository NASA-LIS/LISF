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
! !ROUTINE: readcrd_era5
! \label{readcrd_era5}
!
! !REVISION HISTORY:
! 23 Dec 2019: Sujay Kumar, initial code 
!
! !INTERFACE:    
subroutine readcrd_era5()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod
  use era5_forcingMod, only : era5_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to ERA5 forcing
!  from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n,t,rc

  call ESMF_ConfigFindLabel(LIS_config,"ERA5 forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,era5_struc(n)%era5dir,&
          rc=rc)
     call LIS_verify(rc,&
          'ERA5 forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"ERA5 forcing tile to grid mapping file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,era5_struc(n)%mapfile,&
          rc=rc)
     call LIS_verify(rc,&
          'ERA5 forcing tile to grid mapping file: not defined')
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using ERA5 forcing'
     write(LIS_logunit,*) '[INFO] ERA5 forcing directory: ',&
          trim(era5_struc(n)%era5DIR)

     era5_struc(n)%era5time1 = 3000.0
     era5_struc(n)%era5time2 = 0.0

  enddo
end subroutine readcrd_era5
