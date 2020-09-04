!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_gefs
!  \label{readcrd_gefs}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
! 1 Jul 2019: K. Arsenault, expand support for GEFS forecasts
!
! !INTERFACE:    
subroutine readcrd_gefs()
! !USES:
  use ESMF
  use LIS_coreMod,     only : LIS_rc, LIS_config
  use LIS_logMod,      only : LIS_logunit, LIS_verify
  use gefs_forcingMod, only : gefs_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to GEFS forecast forcing 
!   from the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n,rc

  write(LIS_logunit,*) ' --- '
  write(unit=LIS_logunit,fmt=*)'[INFO] Using GEFS forecast forcing'

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast directory:",rc=rc)
  call LIS_verify(rc, 'GEFS forecast directory: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_dir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast type:",rc=rc)
  call LIS_verify(rc, 'GEFS forecast type: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_fcsttype,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast run mode:",rc=rc)
  call LIS_verify(rc, 'GEFS forecast run mode: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_runmode,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast grid projection:",rc=rc)
  call LIS_verify(rc, 'GEFS forecast grid projection: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_proj,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast number of ensemble members:",rc=rc)
  call LIS_verify(rc, 'GEFS forecast number of ensemble members: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%max_ens_members,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEFS pressure level field:",rc=rc)
  call LIS_verify(rc, 'GEFS pressure level field: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_preslevel,rc=rc)
  enddo

  do n=1,LIS_rc%nnest

     write(LIS_logunit,*) '[INFO] GEFS forecast directory:  ',&
          trim(gefs_struc(n)%gefs_dir)
     write(LIS_logunit,*) '[INFO] GEFS forecast type:  ',&
          gefs_struc(n)%gefs_fcsttype
     write(LIS_logunit,*) '[INFO] GEFS forecast run mode:  ',&
          gefs_struc(n)%gefs_runmode
     write(LIS_logunit,*) '[INFO] GEFS forecast forecast projection:  ',&
          gefs_struc(n)%gefs_proj
     write(LIS_logunit,*) '[INFO] GEFS forecast number of ensemble members:',&
          gefs_struc(n)%max_ens_members
     write(LIS_logunit,*) '[INFO] GEFS pressure level field:  ',&
          gefs_struc(n)%gefs_preslevel

     gefs_struc(n)%fcsttime1 = 3000.0
     gefs_struc(n)%fcsttime2 = 0.0
  enddo

end subroutine readcrd_gefs


