!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_gdas3d
! \label{readcrd_gdas3d}
!
!
! !REVISION HISTORY:
! 18 Mar 2009; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_gdas3d()
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_config
  use LIS_logMod,    only : LIS_logunit, LIS_verify
  use gdas3d_forcingMod, only: gdas3d_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GDAS3D forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"GDAS3D forcing directory:",rc=rc)
  call LIS_verify(rc,'GDAS3D forcing directory: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdas3d_struc(n)%gdasdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GDAS3D domain x-dimension size:",rc=rc)
  call LIS_verify(rc,'GDAS3D domain x-dimension size: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdas3d_struc(n)%nc,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GDAS3D domain y-dimension size:",rc=rc)
  call LIS_verify(rc,'GDAS3D domain y-dimension size: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdas3d_struc(n)%nr,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GDAS3D domain z-dimension size:",rc=rc)
  call LIS_verify(rc,'GDAS3D domain z-dimension size: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdas3d_struc(n)%nlayer,rc=rc)
  enddo

  do n=1,LIS_rc%nnest

     write(LIS_logunit,*) 'Using GDAS3D forcing'
     write(LIS_logunit,*) 'GDAS3D forcing directory :',gdas3d_struc(n)%GDASDIR
     gdas3d_struc(n)%GDASTIME1  = 3000.0
     gdas3d_struc(n)%GDASTIME2  = 0.0
  enddo

end subroutine readcrd_gdas3d
