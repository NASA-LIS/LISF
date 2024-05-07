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
! !ROUTINE: readcrd_plumber2
! \label{readcrd_plumber2}
!
! !REVISION HISTORY:
! 15 Sep 2021: Mark Beauharnois, based on Bondville and GSWP2 readers
!
! !INTERFACE:
      subroutine readcrd_plumber2()
! !USES:
      use ESMF
      use LIS_logMod, only : LIS_logunit
      use LIS_coreMod, only : LIS_rc, LIS_config
      use plumber2_forcingMod, only : plumber2_struc

! !DESCRIPTION:
!  This routine reads the options specific to a PLUMBER2 test
!  case station data forcing from the LIS configuration file.
!
!EOP
      implicit none
      integer :: n, rc

      call ESMF_ConfigFindLabel(LIS_config,                            &
                                "PLUMBER2 forcing file:",rc=rc)
      do n = 1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,                      &
                            plumber2_struc(n)%plumber2file,rc=rc)
      enddo

      call ESMF_ConfigFindLabel(LIS_config,                           &
                                "PLUMBER2 Station ID:",rc=rc)
      do n = 1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,                      &
                            plumber2_struc(n)%stnid,rc=rc)
      enddo

      call ESMF_ConfigFindLabel(LIS_config,                            &
                                "PLUMBER2 Time Delta:",rc=rc)

      do n = 1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,plumber2_struc(n)%ts,&
                              rc=rc)
      enddo

      do n = 1,LIS_rc%nnest
         write(LIS_logunit,*) "[INFO] Using PLUMBER2 forcing"
         write(LIS_logunit,*) "[INFO] PLUMBER2 forcing file: ",       &
                               trim(plumber2_struc(n)%plumber2file)

         plumber2_struc(n)%plumber2time2 = 3000.0
         plumber2_struc(n)%plumber2time1 = 0.0

      enddo

end subroutine readcrd_plumber2
