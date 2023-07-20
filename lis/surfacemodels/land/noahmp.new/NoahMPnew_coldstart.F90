!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LIS_misc.h"
!BOP
!
! !ROUTINE: NoahMPnew_coldstart
! \label{NoahMPnew_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  10/25/18: Shugong Wang, Zhuo Wang; initial implementation for LIS 7 and NoahMP401
!  May 2023: Cenlin He; update to work with refactored NoahMP (v5.0 and newer)
!
! !INTERFACE:

subroutine NoahMPnew_coldstart(mtype)
! !USES:
!   use LIS_coreMod, only: LIS_rc
!   use LIS_logMod, only: LIS_logunit
!   use LIS_timeMgrMod, only: LIS_date2time

    use LIS_coreMod
    use LIS_logMod
    use LIS_histDataMod

    use LIS_timeMgrMod, only: LIS_date2time
    use NoahMPnew_lsmMod
    use NoahmpIOVarType
    use NoahmpInitMainMod, only : NoahmpInitMain
!
! !DESCRIPTION:
!
!  This routine initializes the NoahMP state variables with
!  some predefined values constantly for the entire domain.
!
!EOP
 
    implicit none

    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    integer ::     row, col

!---------------------------------------------------------------------

    !-------- initialize NoahmpIO dimension (1-D)
    NoahmpIO%xstart = 1
    NoahmpIO%xend   = 1
    NoahmpIO%ystart = 1
    NoahmpIO%yend   = 1
    NoahmpIO%ids    = NoahmpIO%xstart
    NoahmpIO%ide    = NoahmpIO%xend
    NoahmpIO%jds    = NoahmpIO%ystart
    NoahmpIO%jde    = NoahmpIO%yend
    NoahmpIO%kds    = 1
    NoahmpIO%kde    = 2
    NoahmpIO%its    = NoahmpIO%xstart
    NoahmpIO%ite    = NoahmpIO%xend
    NoahmpIO%jts    = NoahmpIO%ystart
    NoahmpIO%jte    = NoahmpIO%yend
    NoahmpIO%kts    = 1
    NoahmpIO%kte    = 2
    NoahmpIO%ims    = NoahmpIO%xstart
    NoahmpIO%ime    = NoahmpIO%xend
    NoahmpIO%jms    = NoahmpIO%ystart
    NoahmpIO%jme    = NoahmpIO%yend
    NoahmpIO%kms    = 1
    NoahmpIO%kme    = 2

    ! start initialization of NoahmpIO and NoahmpNew_struc
    do n=1, LIS_rc%nnest

        ! initialize NoahmpIO variables
        NoahmpIO%nsoil        = NoahMPnew_struc(n)%nsoil
        NoahmpIO%nsnow        = NoahMPnew_struc(n)%nsnow
        NoahmpIO%restart_flag = .false.
        NoahmpIO%FNDSNOWH     = .true.

        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) &
             "[INFO] NoahMPnew_coldstart -- cold-starting Noah-MP.New"

           do t=1, LIS_rc%npatch(n,mtype)

             NoahmpIO%ISNOWXY(1,1) = -NoahmpNew_struc(n)%nsnow
             NoahmpIO%ZSOIL(1) = -NoahmpNew_struc(n)%sldpth(1)
             do l=2, NoahmpNew_struc(n)%nsoil
               NoahmpIO%ZSOIL(l) = NoahmpIO%ZSOIL(l-1) - NoahmpNew_struc(n)%sldpth(l)
             enddo

             do l=1, NoahmpNew_struc(n)%nsoil
               NoahmpIO%TSLB(1,l,1)  = NoahmpNew_struc(n)%init_tslb(l)
               NoahmpIO%SMOIS(1,l,1) = NoahmpNew_struc(n)%init_smc(l)
               NoahmpIO%SH2O(1,l,1)  = 0.0
               NoahmpIO%DZS(l)       = NoahmpNew_struc(n)%sldpth(l)
             enddo

             NoahmpIO%croptype(1,:,1) = 0.0
             row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
             col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
             NoahmpIO%XLAT(1,1) = LIS_domain(n)%grid(LIS_domain(n)%gindex(col,row))%lat

             NoahmpIO%CROPCAT(1,1)     = LIS_rc%cropclass
             NoahmpIO%CANWAT(1,1)      = NoahmpNew_struc(n)%init_canwat
             NoahmpIO%IVGTYP(1,1)      = NoahmpNew_struc(n)%noahmpnew(t)%vegetype
             NoahmpIO%ISLTYP(1,1)      = NoahmpNew_struc(n)%noahmpnew(t)%soiltype
             NoahmpIO%TSK(1,1)         = NoahmpNew_struc(n)%init_tskin 
             NoahmpIO%tvxy(1,1)        = NoahmpNew_struc(n)%init_tskin
             NoahmpIO%tgxy(1,1)        = NoahmpNew_struc(n)%init_tskin
             NoahmpIO%canicexy(1,1)    = 0.0 
             NoahmpIO%canliqxy(1,1)    = 0.0
             NoahmpIO%TMN(1,1)         = NoahmpNew_struc(n)%noahmpnew(t)%tbot
             NoahmpIO%XICE(1,1)        = 0.0
             NoahmpIO%eahxy(1,1)       = 0.0
             NoahmpIO%tahxy(1,1)       = 0.0
             NoahmpIO%cmxy(1,1)        = 0.0   
             NoahmpIO%chxy(1,1)        = 0.0   
             NoahmpIO%fwetxy(1,1)      = 0.0
             NoahmpIO%sneqvoxy(1,1)    = 0.0 
             NoahmpIO%alboldxy(1,1)    = 0.0 
             NoahmpIO%qsnowxy(1,1)     = 0.0
             NoahmpIO%wslakexy(1,1)    = 0.0 
             NoahmpIO%zwtxy(1,1)       = NoahmpNew_struc(n)%init_zwt 
             NoahmpIO%waxy(1,1)        = NoahmpNew_struc(n)%init_wa 
             NoahmpIO%wtxy(1,1)        = NoahmpNew_struc(n)%init_wt
             NoahmpIO%lfmassxy(1,1)    = 0.0
             NoahmpIO%rtmassxy(1,1)    = 0.0
             NoahmpIO%stmassxy(1,1)    = 0.0 
             NoahmpIO%woodxy(1,1)      = 0.0 
             NoahmpIO%stblcpxy(1,1)    = 0.0 
             NoahmpIO%fastcpxy(1,1)    = 0.0 
             NoahmpIO%xsaixy(1,1)      = 0.0 
             NoahmpIO%lai(1,1)         = NoahmpNew_struc(n)%init_lai
             NoahmpIO%grainxy(1,1)     = 0.0
             NoahmpIO%gddxy(1,1)       = 0.0
             NoahmpIO%cropcat(1,1)     = 0.0
             NoahmpIO%t2mvxy(1,1)      = 0.0 
             NoahmpIO%t2mbxy(1,1)      = 0.0 
             NoahmpIO%IOPT_RUNSUB      = NoahmpNew_struc(n)%runsub_opt
             NoahmpIO%IOPT_CROP        = NoahmpNew_struc(n)%crop_opt
             NoahmpIO%sf_urban_physics = NoahmpNew_struc(n)%urban_opt
             NoahmpIO%IOPT_IRR         = NoahmpNew_struc(n)%irr_opt
             NoahmpIO%IOPT_IRRM        = NoahmpNew_struc(n)%irrm_opt

             ! The following variables are optional for groundwater dynamics iopt_run=5, 
             !  so some random values are set temporarily.
             do l=1, NoahmpNew_struc(n)%nsoil
               NoahmpIO%smoiseq(1,l,1) = 0.0
             end do
             NoahmpIO%smcwtdxy(1,1)    = 0.0
             NoahmpIO%rechxy(1,1)      = 0.0
             NoahmpIO%deeprechxy(1,1)  = 0.0
             NoahmpIO%areaxy(1,1)      = 100.0
             NoahmpIO%dx               = 10.0
             NoahmpIO%dy               = 10.0
             NoahmpIO%msftx(1,1)       = 0.0
             NoahmpIO%msfty(1,1)       = 0.0
             NoahmpIO%wtddt            = 30.0
             NoahmpIO%stepwtd          = 0
             NoahmpIO%dtbl             = NoahMPnew_struc(n)%ts
             NoahmpIO%qrfsxy(1,1)      = 0.0
             NoahmpIO%qslatxy(1,1)     = 0.0
             NoahmpIO%fdepthxy(1,1)    = 0.0
             NoahmpIO%riverbedxy(1,1)  = 0.0
             NoahmpIO%eqzwt(1,1)       = 0.0
             NoahmpIO%rivercondxy(1,1) = 0.0
             NoahmpIO%pexpxy(1,1)      = 0.0
             NoahmpIO%rechclim(1,1)    = 0.0
             NoahmpIO%sfcrunoff(1,1)   = 0.0
             NoahmpIO%udrunoff(1,1)    = 0.0
             NoahmpIO%acsnom(1,1)      = 0.0
             NoahmpIO%acsnow(1,1)      = 0.0
             NoahmpIO%taussxy(1,1)     = 0.0
             NoahmpIO%pgsxy(1,1)       = 0
             NoahmpIO%zsnsoxy(1,:,1)   = 0.0
             NoahmpIO%SNOW(1,1)        = NoahMPnew_struc(n)%init_sneqv
             NoahmpIO%SNOWH(1,1)       = NoahMPnew_struc(n)%init_snowh
             NoahmpIO%snicexy(1,:,1)   = 0.0
             NoahmpIO%snliqxy(1,:,1)   = 0.0
             NoahmpIO%tsnoxy(1,:,1)    = 0.0

             ! main Noah-MP initialization module
             call NoahmpInitMain(NoahmpIO)

             ! update NoahmpNew_struc initial values
             NoahmpNew_struc(n)%noahmpnew(t)%sfcrunoff  = NoahmpIO%sfcrunoff(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%udrrunoff  = NoahmpIO%udrunoff(1,1)
             do l=1, NoahmpNew_struc(n)%nsoil
                NoahmpNew_struc(n)%noahmpnew(t)%smc(l)  = NoahmpIO%SMOIS(1,l,1) 
             enddo
             do l=1, NoahmpNew_struc(n)%nsoil
                NoahmpNew_struc(n)%noahmpnew(t)%sh2o(l) = NoahmpIO%SH2O(1,l,1)
             enddo
             do l=1, NoahmpNew_struc(n)%nsoil
                NoahmpNew_struc(n)%noahmpnew(t)%tslb(l) = NoahmpIO%TSLB(1,l,1)
             enddo
             NoahmpNew_struc(n)%noahmpnew(t)%sneqv    = NoahmpIO%SNOW(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%snowh    = NoahmpIO%SNOWH(1,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%canwat   = NoahmpIO%canwat(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%acsnom   = NoahmpIO%acsnom(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%acsnow   = NoahmpIO%acsnow(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%isnow    = NoahmpIO%isnowxy(1,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%tv       = NoahmpIO%tvxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%tg       = NoahmpIO%tgxy(1,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%canice   = NoahmpIO%canicexy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%canliq   = NoahmpIO%canliqxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%fwet     = NoahmpIO%fwetxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%sneqvo   = NoahmpIO%sneqvoxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%albold   = NoahmpIO%alboldxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%qsnow    = NoahmpIO%qsnowxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%wslake   = NoahmpIO%wslakexy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%zwt      = NoahmpIO%zwtxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%wa       = NoahmpIO%waxy(1,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%wt       = NoahmpIO%wtxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%lfmass   = NoahmpIO%lfmassxy(1,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%rtmass   = NoahmpIO%rtmassxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%stmass   = NoahmpIO%stmassxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%wood     = NoahmpIO%woodxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%stblcp   = NoahmpIO%stblcpxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%fastcp   = NoahmpIO%fastcpxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%lai      = NoahmpIO%lai(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%sai      = NoahmpIO%xsaixy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%tauss    = NoahmpIO%taussxy(1,1)
             do l=1, NoahmpNew_struc(n)%nsoil
                NoahmpNew_struc(n)%noahmpnew(t)%smoiseq(l) = NoahmpIO%smoiseq(1,l,1) 
             enddo
             NoahmpNew_struc(n)%noahmpnew(t)%smcwtd   = NoahmpIO%smcwtdxy(1,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%deeprech = NoahmpIO%deeprechxy(1,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%rech     = NoahmpIO%rechxy(1,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%grain    = NoahmpIO%grainxy(1,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%gdd      = NoahmpIO%gddxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%pgs      = NoahmpIO%pgsxy(1,1)
             NoahmpNew_struc(n)%noahmpnew(t)%snowice(1:NoahmpNew_struc(n)%nsnow) = &
                                                        NoahmpIO%snicexy(1,-NoahmpNew_struc(n)%nsnow+1:0,1)
             NoahmpNew_struc(n)%noahmpnew(t)%snowliq(1:NoahmpNew_struc(n)%nsnow) = &
                                                        NoahmpIO%snliqxy(1,-NoahmpNew_struc(n)%nsnow+1:0,1)
             NoahmpNew_struc(n)%noahmpnew(t)%zss(1:NoahmpNew_struc(n)%nsnow+NoahmpNew_struc(n)%nsoil) = &
                                                        NoahmpIO%zsnsoxy(1,-NoahmpNew_struc(n)%nsnow+1:NoahmpNew_struc(n)%nsoil,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%isnow = NoahmpIO%isnowxy(1,1) 
             NoahmpNew_struc(n)%noahmpnew(t)%tsno(1:NoahmpNew_struc(n)%nsnow) = &
                                                        NoahmpIO%tsnoxy(1,-NoahmpNew_struc(n)%nsnow+1:0,1)
!-----------------------------------------------------------------------
            enddo   ! t=1,1
        endif       ! coldstart

        LIS_rc%yr = LIS_rc%syr
        LIS_rc%mo = LIS_rc%smo
        LIS_rc%da = LIS_rc%sda
        LIS_rc%hr = LIS_rc%shr
        LIS_rc%mn = LIS_rc%smn
        LIS_rc%ss = LIS_rc%sss
        
        call LIS_date2time(LIS_rc%time, LIS_rc%doy, LIS_rc%gmt, LIS_rc%yr,      &
                           LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
        write(LIS_logunit,*) "[INFO] NoahMPnew_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
       enddo        ! nnest

end subroutine NoahMPnew_coldstart
