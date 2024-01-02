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
! \label{cable_coldstart}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar, Initial version
!   8 May 2009: Sujay Kumar, additions for Noah3.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!  27 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!  14 Dec 2011: Claire Carouge (ccc), CABLE LSM improvements
!  23 May 2013: David Mocko, latest CABLE v1.4b version for LIS6.2
!
! !INTERFACE:
subroutine cable_coldstart
! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit
  use LIS_timeMgrMod,     only : LIS_date2time
  use cable_dimensions,   only : ms,msn,ncp,ncs
  use cable_lsmMod
!
! !DESCRIPTION:
!
!  This subroutine initializes the CABLE state variables with some
!  predefined values uniformly for the entire domain, except for
!  cplant and csoil, which are initialized by vegetation type.
!
!EOP
  implicit none
  
  integer :: t,l,n,i
  real :: csice, rhowat, cswat, tfrz, xx

!ccc Move initialisation from cable_soilsnow to here. Dirty: constants are hard-wired,
!    see where to put them to get a nice implementation.
  csice  = 2.100e3
  rhowat = 1000.0
  cswat  = 4.218e3
  tfrz   = 273.16

  do n = 1,LIS_rc%nnest
     if (LIS_rc%startcode.eq."coldstart") then
        write(LIS_logunit,*)                                       &
             'MSG: cable_coldstart -- cold-starting CABLE'

        ! Initialise time step counter
        cable_struc(n)%ktau = 1

        do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           cable_struc(n)%cable(t)%cansto = 0.0
           cable_struc(n)%cable(t)%rtsoil = 100.0
           cable_struc(n)%cable(t)%ssdnn = 140.0
           cable_struc(n)%cable(t)%snowd = 0.0
           cable_struc(n)%cable(t)%osnowd = 0.0
           cable_struc(n)%cable(t)%snage = 0.0
           cable_struc(n)%cable(t)%isflag = 0
           do l = 1,ms
              cable_struc(n)%cable(t)%tgg(l) = 300.0
              cable_struc(n)%cable(t)%wb(l) = 0.3
              if (cable_struc(n)%cable(t)%tgg(l).lt.tfrz) then
                 cable_struc(n)%cable(t)%wbice(l) = 0.8 *          &
                      cable_struc(n)%cable(t)%wb(l)
              else
                 cable_struc(n)%cable(t)%wbice(l) = 0.0
              endif
           enddo
           do l = 1,msn
              cable_struc(n)%cable(t)%tggsn(l) = 273.1
              cable_struc(n)%cable(t)%ssdn(l) = 140.0
              cable_struc(n)%cable(t)%smass(l) = 0.0
           enddo
           cable_struc(n)%cable(t)%albsoilsn(1) = 0.1 ! soil+snow albedo
           cable_struc(n)%cable(t)%albsoilsn(2) = 0.3
           cable_struc(n)%cable(t)%albsoilsn(3) = 0.05
           do l = 1,ncp
              cable_struc(n)%cable(t)%cplant(l) =                  &
                   cable_struc(n)%cable(t)%cplant_init(l)
           enddo
           do l = 1,ncs
              cable_struc(n)%cable(t)%csoil(l) =                   &
                   cable_struc(n)%cable(t)%csoil_init(l)
           enddo

           !ccc Real initialization of pwb_min in cable_soilsnow.f90
           !ccc               cable_struc(n)%cable(t)%pwb_min = 0.0
           cable_struc(n)%cable(t)%gammzz  = 0.0
           cable_struc(n)%cable(t)%wbtot   = 0.0

           DO l = 1, ms
              cable_struc(n)%cable(t)%wbtot = cable_struc(n)%cable(t)%wbtot + &
                   cable_struc(n)%cable(t)%wb(l)*1000.0*   &
                   cable_struc(n)%cable(t)%zse(l)
              ! Update soil ice:
              IF (cable_struc(n)%cable(t)%tgg(l) <= tfrz .AND. &
                   cable_struc(n)%cable(t)%wbice(l) <= 0.01) then
                 cable_struc(n)%cable(t)%wbice(l) = 0.1 *   &
                      cable_struc(n)%cable(t)%wb(l)
              END IF
              IF (cable_struc(n)%cable(t)%tgg(l) < tfrz) THEN
                 cable_struc(n)%cable(t)%wbice(l) = MIN(0.98 * cable_struc(n)%cable(t)%wb(l), &
                      cable_struc(n)%cable(t)%ssat)
              END IF
           END DO
           xx=cable_struc(n)%cable(t)%css * cable_struc(n)%cable(t)%rhosoil
           cable_struc(n)%cable(t)%gammzz(1) = MAX(      &
                (1.0 - cable_struc(n)%cable(t)%ssat) *    &
                cable_struc(n)%cable(t)%css * cable_struc(n)%cable(t)%rhosoil &
                + (cable_struc(n)%cable(t)%wb(1) - cable_struc(n)%cable(t)%wbice(1) )  &
                * cswat * rhowat &
                + cable_struc(n)%cable(t)%wbice(1) * csice * rhowat * .9, xx )   &
                * cable_struc(n)%cable(t)%zse(1)

           !ccc Initialisation moved from soilsnow
           cable_struc(n)%cable(t)%ssdn = 120.0
           cable_struc(n)%cable(t)%ssdnn = 120.0
           cable_struc(n)%cable(t)%tggsn = tfrz
           cable_struc(n)%cable(t)%isflag = 0

        enddo
     endif

     LIS_rc%yr = LIS_rc%syr
     LIS_rc%mo = LIS_rc%smo
     LIS_rc%da = LIS_rc%sda
     LIS_rc%hr = LIS_rc%shr
     LIS_rc%mn = LIS_rc%smn
     LIS_rc%ss = LIS_rc%sss

     call LIS_date2time(LIS_rc%time,LIS_rc%doy,LIS_rc%gmt,         &
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da,             &
          LIS_rc%hr,LIS_rc%mn,LIS_rc%ss)
     write(LIS_logunit,*) 'MSG: cable_coldstart -- ',              &
          'Using the specified start time ',LIS_rc%time

  enddo

end subroutine cable_coldstart
