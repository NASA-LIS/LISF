!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_mogrepsg
! \label{readcrd_mogrepsg}
!
! !REVISION HISTORY:
! 26 Jan 2023; Yeosang Yoon, Initial Code
! 01 Jan 2024; Yeosang Yoon; update codes for precpi. bias-correction
!
! !INTERFACE:
subroutine readcrd_mogrepsg()
! !USES:
  use LIS_logMod
  use LIS_coreMod
  use mogrepsg_forcingMod, only : mogrepsg_struc
  use ESMF
!
! !DESCRIPTION:
!
!  This routine reads the options specific to MOGREPS-G forecast forcing from
!  the LIS configuration file.
!
!EOP

  implicit none

  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config, &
       "MOGREPS-G forecast forcing directory:", rc=rc)
  do n = 1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, mogrepsg_struc(n)%odir, &
          rc=rc)
     call LIS_verify(rc, &
          'MOGREPS-G forecast forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "MOGREPS-G forecast run mode:", &
       rc=rc)
  call LIS_verify(rc, 'MOGREPS-G forecast run mode: not defined ')
  do n = 1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          mogrepsg_struc(n)%runmode, rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, &
       "MOGREPS-G forecast number of ensemble members:", rc=rc)
  call LIS_verify(rc, &
       'MOGREPS-G forecast number of ensemble members: not defined')
  do n = 1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          mogrepsg_struc(n)%max_ens_members, rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, &
       "Apply MOGREPS-G precipitation bias correction:", rc=rc)
  call LIS_verify(rc, &
       'Apply MOGREPS-G precipitation bias correction: not defined')
  do n = 1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, mogrepsg_struc(n)%bc, rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "MOGREPS-G model CDF file:", rc=rc)
  call LIS_verify(rc, 'MOGREPS-G model CDF file: not defined')
  do n = 1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          mogrepsg_struc(n)%cdf_fname, rc=rc)
  enddo

  do n = 1, LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using MOGREPS-G forecast forcing'
     write(LIS_logunit,*) &
          '[INFO] MOGREPS-G forecast forcing directory: ', &
          trim(mogrepsg_struc(n)%odir)
     write(LIS_logunit,*) '[INFO] MOGREPS-G forecast run mode: ', &
          mogrepsg_struc(n)%runmode
     write(LIS_logunit,*) &
          '[INFO] MOGREPS-G forecast number of ensemble members:', &
           mogrepsg_struc(n)%max_ens_members
     write(LIS_logunit,*) &
          '[INFO] Using MOGREPS-G precipitation bias correction:',&
           mogrepsg_struc(n)%bc
     if (mogrepsg_struc(n)%bc == 1) then
        write(LIS_logunit,*) '[INFO] MOGREPS-G model CDF file: ', &
             trim(mogrepsg_struc(n)%cdf_fname)
     endif
  enddo
end subroutine readcrd_mogrepsg
