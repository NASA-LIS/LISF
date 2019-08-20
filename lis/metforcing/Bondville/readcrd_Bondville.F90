!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_Bondville
! \label{readcrd_Bondville}
!
! !REVISION HISTORY:
! 13 Apr 2007: Bailing Li, Initial Code
! 05 Oct 2010: David Mocko, Updated for Bondville test case
! 26 Oct 2018: David Mocko, Updated for Noah-MP-4.0.1 HRLDAS test case
!
! !INTERFACE:
      subroutine readcrd_Bondville()
! !USES:
      use ESMF
      use LIS_logMod, only : LIS_logunit
      use LIS_coreMod, only : LIS_rc, LIS_config
      use Bondville_forcingMod, only : Bondville_struc

! !DESCRIPTION:
!  This routine reads the options specific to the Bondville
!  test case station data forcing from the LIS configuration file.
!
!EOP
      implicit none
      integer :: n, rc

      call ESMF_ConfigFindLabel(LIS_config,                            &
                                "Bondville forcing file:",rc=rc)
      do n = 1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,                      &
                            Bondville_struc(n)%Bondvillefile,rc=rc)
      enddo

      do n = 1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,Bondville_struc(n)%MP,&
                              label="Bondville Noah-MP-4.0.1 forcing:",&
                              default=0,rc=rc)
      enddo

      do n = 1,LIS_rc%nnest
         write(LIS_logunit,*) "Using Bondville forcing"
         write(LIS_logunit,*) "Bondville forcing file: ",             &
                               trim(Bondville_struc(n)%Bondvillefile)
         Bondville_struc(n)%Bondvilletime2 = 3000.0
         Bondville_struc(n)%Bondvilletime1 = 0.0
         Bondville_struc(n)%startRead = .false.
      enddo

      end subroutine readcrd_Bondville

