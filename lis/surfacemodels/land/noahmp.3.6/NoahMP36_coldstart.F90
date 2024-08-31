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
! !ROUTINE: NoahMP36_coldstart
! \label{NoahMP36_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   9/4/14: Shugong Wang; initial implementation for LIS 7 and NoahMP36
!
! !INTERFACE:
subroutine NoahMP36_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use NoahMP36_lsmMod
   
!
! !DESCRIPTION:
!
!  This routine initializes the NoahMP36 state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
 
    implicit none
    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    
    ! added by Shugong Wang
    integer :: isnow
    real, allocatable, dimension(:) :: zsnso 
    real, allocatable, dimension(:) :: tsno
    real, allocatable, dimension(:) :: snice
    real, allocatable, dimension(:) :: snliq 
    real, allocatable, dimension(:) :: zsoil 
    ! end add

    !EMK...Temporary arrays.
    real :: tmp_swe(1,1),tmp_tgxy(1,1),tmp_snodep(1,1)
    integer :: tmp_isnowxy(1,1)

    do n=1, LIS_rc%nnest
       isnow = -NOAHMP36_struc(n)%nsnow
        ! added by shugong
        allocate(zsnso(-NOAHMP36_struc(n)%nsnow+1:NOAHMP36_struc(n)%nsoil))
        allocate(tsno(-NOAHMP36_struc(n)%nsnow+1:0))
        allocate(snice(-NOAHMP36_struc(n)%nsnow+1:0))
        allocate(snliq(-NOAHMP36_struc(n)%nsnow+1:0))
        allocate(zsoil(NOAHMP36_struc(n)%nsoil))
        zsoil(1) = -NOAHMP36_struc(n)%sldpth(1)
        do l=2, NOAHMP36_struc(n)%nsoil
          zsoil(l) = zsoil(l-1) - NOAHMP36_struc(n)%sldpth(l) 
        enddo
        ! end add 

        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "MSG: NoahMP36_coldstart -- cold-starting NoahMP36"
            do t=1, LIS_rc%npatch(n,mtype)
                NOAHMP36_struc(n)%noahmp36(t)%albold = NOAHMP36_struc(n)%init_albold
                NOAHMP36_struc(n)%noahmp36(t)%sneqvo = NOAHMP36_struc(n)%init_sneqvo
                ! only soil temperature is intialized, snow temperature is calculated by snow_init 
                do l=1, NOAHMP36_struc(n)%nsoil
                    NOAHMP36_struc(n)%noahmp36(t)%sstc(NOAHMP36_struc(n)%nsnow+l) = NOAHMP36_struc(n)%init_stc(l)
                enddo
                do l=1, NOAHMP36_struc(n)%nsoil
                    NOAHMP36_struc(n)%noahmp36(t)%sh2o(l) = NOAHMP36_struc(n)%init_sh2o(l)
                enddo
                do l=1, NOAHMP36_struc(n)%nsoil
                    NOAHMP36_struc(n)%noahmp36(t)%smc(l) = NOAHMP36_struc(n)%init_smc(l)
                enddo
                NOAHMP36_struc(n)%noahmp36(t)%tah = NOAHMP36_struc(n)%init_tah
                NOAHMP36_struc(n)%noahmp36(t)%eah = NOAHMP36_struc(n)%init_eah
                NOAHMP36_struc(n)%noahmp36(t)%fwet = NOAHMP36_struc(n)%init_fwet
                NOAHMP36_struc(n)%noahmp36(t)%canliq = NOAHMP36_struc(n)%init_canliq
                NOAHMP36_struc(n)%noahmp36(t)%canice = NOAHMP36_struc(n)%init_canice
                NOAHMP36_struc(n)%noahmp36(t)%tv = NOAHMP36_struc(n)%init_tv
                NOAHMP36_struc(n)%noahmp36(t)%tg = NOAHMP36_struc(n)%init_tg
                NOAHMP36_struc(n)%noahmp36(t)%qsnow = NOAHMP36_struc(n)%init_qsnow
                !do l=1, NOAHMP36_struc(n)%nsoil + NOAHMP36_struc(n)%nsnow
                !    NOAHMP36_struc(n)%noahmp36(t)%zss(l) = NOAHMP36_struc(n)%init_zss(l)
                !enddo
                NOAHMP36_struc(n)%noahmp36(t)%snowh = NOAHMP36_struc(n)%init_snowh
                NOAHMP36_struc(n)%noahmp36(t)%sneqv = NOAHMP36_struc(n)%init_sneqv
                !do l=1, NOAHMP36_struc(n)%nsnow
                !    NOAHMP36_struc(n)%noahmp36(t)%snowice(l) = NOAHMP36_struc(n)%init_snowice(l)
                !enddo
                !do l=1, NOAHMP36_struc(n)%nsnow
                !    NOAHMP36_struc(n)%noahmp36(t)%snowliq(l) = NOAHMP36_struc(n)%init_snowliq(l)
                !enddo
                NOAHMP36_struc(n)%noahmp36(t)%zwt = NOAHMP36_struc(n)%init_zwt
                NOAHMP36_struc(n)%noahmp36(t)%wa = NOAHMP36_struc(n)%init_wa
                NOAHMP36_struc(n)%noahmp36(t)%wt = NOAHMP36_struc(n)%init_wt
                NOAHMP36_struc(n)%noahmp36(t)%wslake = NOAHMP36_struc(n)%init_wslake
                NOAHMP36_struc(n)%noahmp36(t)%lfmass = NOAHMP36_struc(n)%init_lfmass
                NOAHMP36_struc(n)%noahmp36(t)%rtmass = NOAHMP36_struc(n)%init_rtmass
                NOAHMP36_struc(n)%noahmp36(t)%stmass = NOAHMP36_struc(n)%init_stmass
                NOAHMP36_struc(n)%noahmp36(t)%wood = NOAHMP36_struc(n)%init_wood
                NOAHMP36_struc(n)%noahmp36(t)%stblcp = NOAHMP36_struc(n)%init_stblcp
                NOAHMP36_struc(n)%noahmp36(t)%fastcp = NOAHMP36_struc(n)%init_fastcp
                NOAHMP36_struc(n)%noahmp36(t)%lai = NOAHMP36_struc(n)%init_lai
                NOAHMP36_struc(n)%noahmp36(t)%sai = NOAHMP36_struc(n)%init_sai
                NOAHMP36_struc(n)%noahmp36(t)%cm = NOAHMP36_struc(n)%init_cm
                NOAHMP36_struc(n)%noahmp36(t)%ch = NOAHMP36_struc(n)%init_ch
                NOAHMP36_struc(n)%noahmp36(t)%tauss = NOAHMP36_struc(n)%init_tauss
                NOAHMP36_struc(n)%noahmp36(t)%smcwtd = NOAHMP36_struc(n)%init_smcwtd
                NOAHMP36_struc(n)%noahmp36(t)%deeprech = NOAHMP36_struc(n)%init_deeprech
                NOAHMP36_struc(n)%noahmp36(t)%rech = NOAHMP36_struc(n)%init_rech
                NOAHMP36_struc(n)%noahmp36(t)%zlvl = NOAHMP36_struc(n)%init_zlvl 

                ! added by shugong 
                zsnso = 0.0 
                !EMK...snow_init_36 is expecting several arrays which
                !are being passed as scalars.  Although no memory corruption
                !occurs here because of the declared array dimensions (all 1),
                !this is still technically a syntax error.  So, we will
                !copy the required fields to temporary arrays with the
                !correct declarations and pass those instead.
                !call snow_init_36(1, 1, 1, 1, 1, 1, 1, 1,           & !input 
                !                  NOAHMP36_struc(n)%nsnow,          & !input 
                !                  NOAHMP36_struc(n)%nsoil,          & !input 
                !                  zsoil,                            & !input
                !                  NOAHMP36_struc(n)%init_sneqv,     & !input
                !                  NOAHMP36_struc(n)%init_tg,        & !input
                !                  NOAHMP36_struc(n)%init_snowh,     & !input
                !                  zsnso, tsno, snice, snliq, isnow) ! output 
                tmp_swe(1,1) = NOAHMP36_struc(n)%init_sneqv
                tmp_tgxy(1,1) = NOAHMP36_struc(n)%init_tg
                tmp_snodep(1,1) = NOAHMP36_struc(n)%init_snowh
                call snow_init_36(1, 1, 1, 1, 1, 1, 1, 1,           & !input 
                                  NOAHMP36_struc(n)%nsnow,          & !input 
                                  NOAHMP36_struc(n)%nsoil,          & !input 
                                  zsoil,                            & !input
                                  tmp_swe,         & !input
                                  tmp_tgxy,        & !input
                                  tmp_snodep,      & !input
                                  zsnso, tsno, snice, snliq, &
                                  tmp_isnowxy) ! output
                isnow = tmp_isnowxy(1,1)
                NOAHMP36_struc(n)%noahmp36(t)%snowice(1:NOAHMP36_struc(n)%nsnow) = snice(-NOAHMP36_struc(n)%nsnow+1:0)
                NOAHMP36_struc(n)%noahmp36(t)%snowliq(1:NOAHMP36_struc(n)%nsnow) = snliq(-NOAHMP36_struc(n)%nsnow+1:0)
                NOAHMP36_struc(n)%noahmp36(t)%zss(1:NOAHMP36_struc(n)%nsnow+NOAHMP36_struc(n)%nsoil) = zsnso(-NOAHMP36_struc(n)%nsnow+1:NOAHMP36_struc(n)%nsoil) 
                NOAHMP36_struc(n)%noahmp36(t)%sstc(NOAHMP36_struc(n)%nsnow+isnow+1:NOAHMP36_struc(n)%nsnow) = tsno(isnow+1:0) 
                NOAHMP36_struc(n)%noahmp36(t)%isnow = isnow                

                
                ! end add
            enddo
        endif
    
        LIS_rc%yr = LIS_rc%syr
        LIS_rc%mo = LIS_rc%smo
        LIS_rc%da = LIS_rc%sda
        LIS_rc%hr = LIS_rc%shr
        LIS_rc%mn = LIS_rc%smn
        LIS_rc%ss = LIS_rc%sss
        
        call LIS_date2time(LIS_rc%time, LIS_rc%doy, LIS_rc%gmt, LIS_rc%yr,      &
                           LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
        write(LIS_logunit,*) "MSG: NoahMP36_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
        deallocate(zsnso)
        deallocate(tsno)
        deallocate(snice)
        deallocate(snliq)
        deallocate(zsoil) 
    enddo
end subroutine NoahMP36_coldstart
