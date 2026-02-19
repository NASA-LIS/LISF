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
! !ROUTINE: readcrd_era5cds
! \label{readcrd_era5cds}
!
! !REVISION HISTORY:
! 23 Dec 2019: Sujay Kumar, initial code 
! 17 Apr 2025: Hiroko Beudoing, adopted ERA5 routines for the public CDS
!                               data format
!
! !INTERFACE:    
subroutine readcrd_era5cds()
! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod
  use era5cds_forcingMod, only : era5cds_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to ERA5 forcing 
!  (downloaded from the Climate Data) from the LDT configuration file. 
!  
!EOP
  implicit none

  integer :: n,t,rc

  call ESMF_ConfigFindLabel(LDT_config,"ERA5CDS forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,era5cds_struc(n)%era5cdsdir,&
          rc=rc)
     call LDT_verify(rc,&
          'ERA5CDS forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"ERA5CDS forcing terrain height file:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,era5cds_struc(n)%era5cdshgt_file,&
          rc=rc)
     call LDT_verify(rc,&
          'ERA5CDS forcing terrain height file: not defined')
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) '[INFO] Using ERA5CDS forcing'
     write(LDT_logunit,*) '[INFO] ERA5CDS forcing directory: ',&
          trim(era5cds_struc(n)%era5cdsDIR)

     era5cds_struc(n)%era5cdstime1 = 3000.0
     era5cds_struc(n)%era5cdstime2 = 0.0

  enddo
end subroutine readcrd_era5cds
