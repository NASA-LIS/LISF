!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: noahmpglacier3911_coldstart
! \label{noahmpglacier3911_coldstart}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine noahmpglacier3911_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use noahmpglacier3911_Mod
    use module_sf_noahmpglacier3911, only : snow_init
   
!
! !DESCRIPTION:
!
!  This routine initializes the noahmpglacier3911 state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
 
    implicit none
    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    
    ! added by Shugong Wang
    integer :: isnow(1) 
    real, allocatable, dimension(:) :: zsnso 
    real, allocatable, dimension(:) :: tsno
    real, allocatable, dimension(:) :: snice
    real, allocatable, dimension(:) :: snliq 
    real, allocatable, dimension(:) :: zsoil 
    ! end add

    do n=1, LIS_rc%nnest
       allocate(zsnso(-noahmpgl3911_struc(n)%nsnow+1:noahmpgl3911_struc(n)%nsoil))
       allocate(tsno(-noahmpgl3911_struc(n)%nsnow+1:0))
       allocate(snice(-noahmpgl3911_struc(n)%nsnow+1:0))
       allocate(snliq(-noahmpgl3911_struc(n)%nsnow+1:0))
       allocate(zsoil(noahmpgl3911_struc(n)%nsoil))
       zsoil(1) = -noahmpgl3911_struc(n)%sldpth(1)
       do l=2, noahmpgl3911_struc(n)%nsoil
          zsoil(l) = zsoil(l-1) - noahmpgl3911_struc(n)%sldpth(l) 
       enddo
       ! end add 
       
       if (trim(LIS_rc%startcode) .eq. "coldstart") then
          write(LIS_logunit,*) "[INFO] noahmpglacier3911_coldstart -- cold-starting noahmpglacier3911"
          do t=1, LIS_rc%npatch(n,mtype)
             noahmpgl3911_struc(n)%noahmpgl(t)%zlvl = &
                  noahmpgl3911_struc(n)%init_zlvl 
             noahmpgl3911_struc(n)%noahmpgl(t)%qsnow = &
                  noahmpgl3911_struc(n)%init_qsnow

             noahmpgl3911_struc(n)%noahmpgl(t)%sneqvo = &
                  noahmpgl3911_struc(n)%init_sneqvo
             noahmpgl3911_struc(n)%noahmpgl(t)%tg = &
                  noahmpgl3911_struc(n)%init_tg(1)
             noahmpgl3911_struc(n)%noahmpgl(t)%sneqv =&
                  noahmpgl3911_struc(n)%init_sneqv(1)
             zsnso = 0.0 
             call snow_init(1, 1, 1, 1, 1, 1, 1, 1,           & !input 
                  noahmpgl3911_struc(n)%nsnow,          & !input 
                  noahmpgl3911_struc(n)%nsoil,          & !input 
                  zsoil,                            & !input
                  noahmpgl3911_struc(n)%init_sneqv,     & !input
                  noahmpgl3911_struc(n)%init_tg,        & !input
                  noahmpgl3911_struc(n)%init_snowh,     & !input
                  zsnso, tsno, snice, snliq, isnow) ! output 

             noahmpgl3911_struc(n)%noahmpgl(t)%snowice(1:noahmpgl3911_struc(n)%nsnow) = snice(-noahmpgl3911_struc(n)%nsnow+1:0)
             noahmpgl3911_struc(n)%noahmpgl(t)%snowliq(1:noahmpgl3911_struc(n)%nsnow) = snliq(-noahmpgl3911_struc(n)%nsnow+1:0)
             noahmpgl3911_struc(n)%noahmpgl(t)%zss(1:noahmpgl3911_struc(n)%nsnow+noahmpgl3911_struc(n)%nsoil) = zsnso(-noahmpgl3911_struc(n)%nsnow+1:noahmpgl3911_struc(n)%nsoil) 
             noahmpgl3911_struc(n)%noahmpgl(t)%sstc(noahmpgl3911_struc(n)%nsnow+isnow(1)+1:noahmpgl3911_struc(n)%nsnow) = tsno(isnow(1)+1:0) 
             noahmpgl3911_struc(n)%noahmpgl(t)%isnow = isnow(1)

             noahmpgl3911_struc(n)%noahmpgl(t)%albold = &
                  noahmpgl3911_struc(n)%init_albold
             noahmpgl3911_struc(n)%noahmpgl(t)%cm = &
                  noahmpgl3911_struc(n)%init_cm
             noahmpgl3911_struc(n)%noahmpgl(t)%ch = &
                  noahmpgl3911_struc(n)%init_ch

             noahmpgl3911_struc(n)%noahmpgl(t)%snowh = &
                  noahmpgl3911_struc(n)%init_snowh(1)

             noahmpgl3911_struc(n)%noahmpgl(t)%tauss = &
                  noahmpgl3911_struc(n)%init_tauss

             do l=1, noahmpgl3911_struc(n)%nsoil
                noahmpgl3911_struc(n)%noahmpgl(t)%sstc(noahmpgl3911_struc(n)%nsnow+l) = &
                     noahmpgl3911_struc(n)%init_stc(l)
             enddo
             do l=1, noahmpgl3911_struc(n)%nsoil
                noahmpgl3911_struc(n)%noahmpgl(t)%sh2o(l) = &
                     noahmpgl3911_struc(n)%init_sh2o(l)
             enddo
             do l=1, noahmpgl3911_struc(n)%nsoil
                noahmpgl3911_struc(n)%noahmpgl(t)%smc(l) = &
                     noahmpgl3911_struc(n)%init_smc(l)
             enddo
          enddo
       endif
    enddo
  end subroutine noahmpglacier3911_coldstart
