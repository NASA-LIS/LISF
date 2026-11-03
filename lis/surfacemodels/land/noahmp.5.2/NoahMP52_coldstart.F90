!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!BOP
!
! !ROUTINE: NoahMP52_coldstart
! \label{NoahMP52_coldstart}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; update to work with refactored NoahMP (v5.0 and newer)
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:

subroutine NoahMP52_coldstart(mtype)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_histDataMod
  use LIS_timeMgrMod, only: LIS_date2time
  use NoahMP52_lsmMod
  use NoahmpInitMainMod, only : NoahmpInitMain
  use NoahmpIOVarType

!
! !DESCRIPTION:
!
!  This routine initializes the NoahMP state variables with
!  some predefined values constantly for the entire domain.
!
!EOP

  implicit none

  integer :: mtype
  integer :: t, l, n
  integer :: row, col

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

  ! start initialization of NoahmpIO and Noahmp52_struc
  do n=1, LIS_rc%nnest

     ! initialize NoahmpIO variables
     NoahmpIO%nsoil        = NoahMP52_struc(n)%nsoil
     NoahmpIO%nsnow        = NoahMP52_struc(n)%nsnow
     NoahmpIO%restart_flag = .false.
     NoahmpIO%FNDSNOWH     = .true.

     if (trim(LIS_rc%startcode) .eq. "coldstart") then
        write(LIS_logunit,*) &
             "[INFO] NoahMP52_coldstart -- cold-starting Noah-MP.5.2"

        do t=1, LIS_rc%npatch(n,mtype)

           NoahmpIO%ISNOWXY(1,1) = -Noahmp52_struc(n)%nsnow
           NoahmpIO%ZSOIL(1) = -Noahmp52_struc(n)%sldpth(1)
           do l=2, Noahmp52_struc(n)%nsoil
              NoahmpIO%ZSOIL(l) = NoahmpIO%ZSOIL(l-1) - &
                   Noahmp52_struc(n)%sldpth(l)
           enddo

           do l=1, Noahmp52_struc(n)%nsoil
              NoahmpIO%TSLB(1,l,1)  = Noahmp52_struc(n)%init_tslb(l)
              NoahmpIO%SMOIS(1,l,1) = Noahmp52_struc(n)%init_smc(l)
              NoahmpIO%SH2O(1,l,1)  = Noahmp52_struc(n)%init_smc(l)
              NoahmpIO%DZS(l)       = Noahmp52_struc(n)%sldpth(l)
           enddo

           NoahmpIO%croptype(1,:,1) = 0.0
           row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
           col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
           NoahmpIO%XLAT(1,1) = &
                LIS_domain(n)%grid(LIS_domain(n)%gindex(col,row))%lat

           NoahmpIO%CROPCAT(1,1)     = LIS_rc%cropclass
           NoahmpIO%CANWAT(1,1)      = Noahmp52_struc(n)%init_canwat
           NoahmpIO%IVGTYP(1,1)      = &
                Noahmp52_struc(n)%noahmp52(t)%vegetype
           NoahmpIO%ISLTYP(1,1)      = &
                Noahmp52_struc(n)%noahmp52(t)%soiltype
           NoahmpIO%TSK(1,1)         = Noahmp52_struc(n)%init_tskin
           NoahmpIO%tvxy(1,1)        = Noahmp52_struc(n)%init_tskin
           NoahmpIO%tgxy(1,1)        = Noahmp52_struc(n)%init_tskin
           NoahmpIO%canicexy(1,1)    = 0.0
           NoahmpIO%canliqxy(1,1)    = 0.0
           NoahmpIO%TMN(1,1)         = Noahmp52_struc(n)%noahmp52(t)%tbot
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
           NoahmpIO%zwtxy(1,1)       = Noahmp52_struc(n)%init_zwt
           NoahmpIO%waxy(1,1)        = Noahmp52_struc(n)%init_wa
           NoahmpIO%wtxy(1,1)        = Noahmp52_struc(n)%init_wt
           NoahmpIO%lfmassxy(1,1)    = 0.0
           NoahmpIO%rtmassxy(1,1)    = 0.0
           NoahmpIO%stmassxy(1,1)    = 0.0
           NoahmpIO%woodxy(1,1)      = 0.0
           NoahmpIO%stblcpxy(1,1)    = 0.0
           NoahmpIO%fastcpxy(1,1)    = 0.0
           NoahmpIO%xsaixy(1,1)      = 0.0
           NoahmpIO%lai(1,1)         = Noahmp52_struc(n)%init_lai
           NoahmpIO%grainxy(1,1)     = 0.0
           NoahmpIO%gddxy(1,1)       = 0.0
           NoahmpIO%cropcat(1,1)     = 0.0
           NoahmpIO%t2mvxy(1,1)      = 0.0
           NoahmpIO%t2mbxy(1,1)      = 0.0
           NoahmpIO%IOPT_RUNSUB      = Noahmp52_struc(n)%runsub_opt
           NoahmpIO%IOPT_CROP        = Noahmp52_struc(n)%crop_opt
           NoahmpIO%sf_urban_physics = Noahmp52_struc(n)%urban_opt
           NoahmpIO%IOPT_IRR         = Noahmp52_struc(n)%irr_opt
           NoahmpIO%IOPT_IRRM        = Noahmp52_struc(n)%irrm_opt

           ! The following variables are optional for groundwater dynamics
           ! iopt_run=5, so some random values are set temporarily.
           do l=1, Noahmp52_struc(n)%nsoil
              NoahmpIO%smoiseq(1,l,1) = 0.0
           end do
           NoahmpIO%smcwtdxy(1,1)    = 0.0
           NoahmpIO%rechxy(1,1)      = 0.0
           NoahmpIO%deeprechxy(1,1)  = 0.0
           NoahmpIO%areaxy(1,1)      = 100.0
           NoahmpIO%wtddt            = 30.0
           NoahmpIO%stepwtd          = 0
           NoahmpIO%dtbl             = NoahMP52_struc(n)%ts
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
           NoahmpIO%SNOW(1,1)        = NoahMP52_struc(n)%init_sneqv
           NoahmpIO%SNOWH(1,1)       = NoahMP52_struc(n)%init_snowh
           NoahmpIO%snicexy(1,:,1)   = 0.0
           NoahmpIO%snliqxy(1,:,1)   = 0.0
           NoahmpIO%tsnoxy(1,:,1)    = 0.0

           ! main Noah-MP initialization module
           call NoahmpInitMain(NoahmpIO)

           ! update Noahmp52_struc initial values
           Noahmp52_struc(n)%noahmp52(t)%sfcrunoff  = &
                NoahmpIO%sfcrunoff(1,1)
           Noahmp52_struc(n)%noahmp52(t)%udrrunoff  = &
                NoahmpIO%udrunoff(1,1)
           do l=1, Noahmp52_struc(n)%nsoil
              Noahmp52_struc(n)%noahmp52(t)%smc(l)  = &
                   NoahmpIO%SMOIS(1,l,1)
           enddo
           do l=1, Noahmp52_struc(n)%nsoil
              Noahmp52_struc(n)%noahmp52(t)%sh2o(l) = NoahmpIO%SH2O(1,l,1)
           enddo
           do l=1, Noahmp52_struc(n)%nsoil
              Noahmp52_struc(n)%noahmp52(t)%tslb(l) = NoahmpIO%TSLB(1,l,1)
           enddo
           Noahmp52_struc(n)%noahmp52(t)%sneqv    = NoahmpIO%SNOW(1,1)
           Noahmp52_struc(n)%noahmp52(t)%snowh    = NoahmpIO%SNOWH(1,1)
           Noahmp52_struc(n)%noahmp52(t)%canwat   = NoahmpIO%canwat(1,1)
           Noahmp52_struc(n)%noahmp52(t)%acsnom   = NoahmpIO%acsnom(1,1)
           Noahmp52_struc(n)%noahmp52(t)%acsnow   = NoahmpIO%acsnow(1,1)
           Noahmp52_struc(n)%noahmp52(t)%isnow    = NoahmpIO%isnowxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%tv       = NoahmpIO%tvxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%tg       = NoahmpIO%tgxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%canice   = NoahmpIO%canicexy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%canliq   = NoahmpIO%canliqxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%fwet     = NoahmpIO%fwetxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%sneqvo   = NoahmpIO%sneqvoxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%albold   = NoahmpIO%alboldxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%qsnow    = NoahmpIO%qsnowxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%wslake   = NoahmpIO%wslakexy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%zwt      = NoahmpIO%zwtxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%wa       = NoahmpIO%waxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%wt       = NoahmpIO%wtxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%lfmass   = NoahmpIO%lfmassxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%rtmass   = NoahmpIO%rtmassxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%stmass   = NoahmpIO%stmassxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%wood     = NoahmpIO%woodxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%stblcp   = NoahmpIO%stblcpxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%fastcp   = NoahmpIO%fastcpxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%lai      = NoahmpIO%lai(1,1)
           Noahmp52_struc(n)%noahmp52(t)%sai      = NoahmpIO%xsaixy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%tauss    = NoahmpIO%taussxy(1,1)
           do l=1, Noahmp52_struc(n)%nsoil
              Noahmp52_struc(n)%noahmp52(t)%smoiseq(l) = &
                   NoahmpIO%smoiseq(1,l,1)
           enddo
           Noahmp52_struc(n)%noahmp52(t)%smcwtd   = NoahmpIO%smcwtdxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%deeprech = &
                NoahmpIO%deeprechxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%rech     = NoahmpIO%rechxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%grain    = NoahmpIO%grainxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%gdd      = NoahmpIO%gddxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%pgs      = NoahmpIO%pgsxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%snowice(1:Noahmp52_struc(n)%nsnow) = &
                NoahmpIO%snicexy(1,-Noahmp52_struc(n)%nsnow+1:0,1)
           Noahmp52_struc(n)%noahmp52(t)%snowliq(1:Noahmp52_struc(n)%nsnow) = &
                NoahmpIO%snliqxy(1,-Noahmp52_struc(n)%nsnow+1:0,1)
           Noahmp52_struc(n)%noahmp52(t)%zss(1:Noahmp52_struc(n)%nsnow+ &
                Noahmp52_struc(n)%nsoil) = &
                NoahmpIO%zsnsoxy(1,-Noahmp52_struc(n)%nsnow+1: &
                Noahmp52_struc(n)%nsoil,1)
           Noahmp52_struc(n)%noahmp52(t)%isnow = NoahmpIO%isnowxy(1,1)
           Noahmp52_struc(n)%noahmp52(t)%tsno(1:Noahmp52_struc(n)%nsnow) = &
                NoahmpIO%tsnoxy(1,-Noahmp52_struc(n)%nsnow+1:0,1)
           !-------------------------------------------------------------------
        enddo   ! t=1,1
     endif      ! coldstart

     LIS_rc%yr = LIS_rc%syr
     LIS_rc%mo = LIS_rc%smo
     LIS_rc%da = LIS_rc%sda
     LIS_rc%hr = LIS_rc%shr
     LIS_rc%mn = LIS_rc%smn
     LIS_rc%ss = LIS_rc%sss

     call LIS_date2time(LIS_rc%time, LIS_rc%doy, LIS_rc%gmt, LIS_rc%yr,  &
          LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
     write(LIS_logunit,*) "[INFO] NoahMP52_coldstart -- ",     &
          "Using the specified start time ", LIS_rc%time
  enddo        ! nnest

end subroutine NoahMP52_coldstart
