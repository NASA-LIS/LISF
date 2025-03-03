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
! !ROUTINE: readcrd_gfs
! \label{readcrd_gfs}
!
!
! !REVISION HISTORY:
!  16 Mar 2008: Sujay Kumar; Initial specification
!
! !INTERFACE:    
subroutine readcrd_gfs()

! !USES:
  use ESMF
  use LDT_coreMod,   only : LDT_rc, LDT_config
  use LDT_logMod,    only : LDT_logunit
  use gfs_forcingMod, only: gfs_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GFS forcing from 
!  the LDT configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LDT_config,"GFS forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,gfs_struc(n)%gfsdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"GFS T126 elevation map:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,gfs_struc(n)%t126elevfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"GFS T170 elevation map:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,gfs_struc(n)%t170elevfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"GFS T254 elevation map:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,gfs_struc(n)%t254elevfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"GFS T382 elevation map:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,gfs_struc(n)%t382elevfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"GFS T574 elevation map:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,gfs_struc(n)%t574elevfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"GFS number of forcing variables:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,gfs_struc(n)%nmif,rc=rc)

     write(LDT_logunit,*) 'Using GFS forcing'
     write(LDT_logunit,*) 'GFS forcing directory :',gfs_struc(n)%GFSDIR
     gfs_struc(n)%GFSTIME1  = 3000.0
     gfs_struc(n)%GFSTIME2  = 0.0
  enddo

end subroutine readcrd_gfs
