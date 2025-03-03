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
! !ROUTINE: clm2_main
! \label{clm2_main}
! 
! !INTERFACE: 
subroutine clm2_main (n)
! !USES: 
  use LIS_precisionMod
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use clm2_lsmMod
  use clm2_varpar, only : nlevsoi
  use clm2_varcon     , only : doalb, eccen, obliqr, lambm0, mvelpp, & 
       denh2o, denice, hvap, hsub, hfus, istwet
  use LIS_coreMod, only : LIS_rc
  use LIS_timeMgrMod  , only : LIS_get_curr_calday, LIS_get_nstep
#if (defined RTM)
  use RtmMod        , only : Rtmriverflux
#endif

#if (defined COUP_CSM)
  use clm_csmMod    , only : csm_dosndrcv, csm_recv, csm_send, csm_flxave, &
                             dorecv, dosend, csmstop_now
#endif
  use clm2_shr_sys_mod   , only : clm2_shr_sys_flush
  use LIS_histDataMod
  use LIS_constantsMod,  only : LIS_CONST_RHOFW

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
!
! !DESCRIPTION: 
! 
! This is the entry point for calling the CLM physics. This routine
! calls various physics routines for CLM that performs the land surface
! computations, to solve for water and energy equations. For a more
! detailed documentation of CLM, please contact the NCAR repository
!
! Calling sequence:
! \begin{verbatim}
! -> loop over patch points calling for each patch point:
!    -> Hydrology1          canopy interception and precip on ground
!    -> Biogeophysics1      leaf temperature and surface fluxes
!    -> Biogeophysics_Lake  lake temperature and surface fluxes
!    -> Biogeophysics2      soil/snow and ground temp and update surface fluxes
!    -> Hydrology2          surface and soil hydrology
!    -> Hydrology_Lake      lake hydrology
!    -> Biogeochemistry     surface biogeochemical fluxes (LSM)
!    -> EcosystemDyn:       ecosystem dynamics: phenology, vegetation,
!                           soil carbon
!    -> SurfaceAlbedo:      albedos for next time step
!      -> SnowAlbedo:       snow albedos: direct beam
!      -> SnowAlbedo:       snow albedos: diffuse
!      -> SoilAlbedo:       soil/lake albedos
!      -> TwoStream:        absorbed, reflected, transmitted 
!                           solar fluxes (vis dir)
!      -> TwoStream:        absorbed, reflected, transmitted 
!                           solar fluxes (vis dif)
!      -> TwoStream:        absorbed, reflected, transmitted 
!                           solar fluxes (nir dir)
!      -> TwoStream:        absorbed, reflected, transmitted 
!                           solar fluxes (nir dif)
!    -> BalanceCheck        check for errors in energy and water balances
! \end{verbatim}
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

  real :: soilm(LIS_rc%npatch(n,LIS_rc%lsm_index),1:nlevsoi)
  real :: soilmtc(LIS_rc%npatch(n,LIS_rc%lsm_index))
  

! ------------------- arguments -----------------------------------
!  logical , intent(in) :: doalb   !true if time for surface albedo calculation
!  real(r8), intent(in) :: eccen   !Earth's orbital eccentricity
!  real(r8), intent(in) :: obliqr  !Earth's obliquity in radians
!  real(r8), intent(in) :: lambm0  !Mean longitude of perihelion at the vernal equinox (radians)
!  real(r8), intent(in) :: mvelpp  !Earth's moving vernal equinox long. of perihelion + pi (radians)
! -----------------------------------------------------------------

! ---------------------- local variables --------------------------
  integer  :: j,k,t,m           !loop/array indices
!  real :: tvegb(144,76)
  real(r8) :: caldayp1            !calendar day for nstep+1
  integer  :: dtime               !timestep size [seconds]
  integer  :: nstep               !timestep index
!  real(r8) :: buf1d(numpatch)     !temporary buffer 
!  real(r8) :: tsxyav              !average ts for diagnostic output
#if (defined SPMD)
!  integer :: numrecvv(0:npes-1)   !vector of items to be received  
!  integer :: displsv(0:npes-1)    !displacement vector
!  integer :: numsend              !number of items to be sent
#endif
  real(r8) :: cgrnd
  real(r8) :: cgrndl
  real(r8) :: cgrnds
  real(r8) :: tg
  real(r8) :: emg
  real(r8) :: htvp
  real(r8) :: dlrad
  real(r8) :: ulrad
  real(r8) :: tssbef(-5:10)
  
  REAL     :: soilmr
  REAL     :: roottemp
  real :: vol_ice(1:nlevsoi)
  real :: vol_liq(1:nlevsoi)
  real :: eff_porosity(1:nlevsoi)
  real :: snowtemp, snowt, totaldepth, asurf, tempvar, soimtc
  integer :: i
  logical             :: alarmCheck
  integer             :: iret
  character*3   :: fnest
  
  write(fnest,'(i3.3)') n
! -----------------------------------------------------------------
!  call t_startf('clm_clm2_main')
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "CLM2 model alarm "//trim(fnest))
  if(alarmCheck) then 

     nstep = LIS_get_nstep(LIS_rc,n)
     
     dtime = LIS_rc%nts(n)
     caldayp1 = LIS_get_curr_calday(LIS_rc,offset=dtime)  
#if (!defined DGVM)
! ----------------------------------------------------------------------
! Determine weights for time interpolation of monthly vegetation data.
! This also determines whether it is time to read new monthly vegetation and
! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
! weights obtained here are used in subroutine ecosystemdyn to obtain time
! interpolated values.
! ----------------------------------------------------------------------
!===LDAS modification: Commented out because LAI etc data read in elsewhere
!!  if (doalb) call interpMonthlyVeg (fsurdat, monp1, dayp1)
#endif

! ----------------------------------------------------------------------
! LOOP 1
! ----------------------------------------------------------------------
  do k = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)    ! begin 1st loop over patches
     clm2_struc(n)%clm(k)%forc_t = clm2_struc(n)%clm(k)%forc_t/&
          clm2_struc(n)%forc_count
     clm2_struc(n)%clm(k)%forc_q  = clm2_struc(n)%clm(k)%forc_q/&
          clm2_struc(n)%forc_count
     clm2_struc(n)%clm(k)%forc_solad = clm2_struc(n)%clm(k)%forc_solad/&
          clm2_struc(n)%forc_count
     clm2_struc(n)%clm(k)%forc_solai = clm2_struc(n)%clm(k)%forc_solai/&
          clm2_struc(n)%forc_count
     clm2_struc(n)%clm(k)%forc_lwrad = clm2_struc(n)%clm(k)%forc_lwrad/&
          clm2_struc(n)%forc_count
     clm2_struc(n)%clm(k)%forc_u     = clm2_struc(n)%clm(k)%forc_u/&
          clm2_struc(n)%forc_count
     clm2_struc(n)%clm(k)%forc_v     = clm2_struc(n)%clm(k)%forc_v/&
          clm2_struc(n)%forc_count
     clm2_struc(n)%clm(k)%forc_pbot  = clm2_struc(n)%clm(k)%forc_pbot/&
          clm2_struc(n)%forc_count
     clm2_struc(n)%clm(k)%forc_rain  = clm2_struc(n)%clm(k)%forc_rain/&
          clm2_struc(n)%forc_count
     clm2_struc(n)%clm(k)%forc_snow  = clm2_struc(n)%clm(k)%forc_snow/&
          clm2_struc(n)%forc_count

     clm2_struc(n)%clm(k)%nstep = nstep
     clm2_struc(n)%clm(k)%h2osno_old = clm2_struc(n)%clm(k)%h2osno  ! snow mass at previous time step
     clm2_struc(n)%clm(k)%frac_veg_nosno = clm2_struc(n)%clm(k)%frac_veg_nosno_alb
     if (clm2_struc(n)%clm(k)%h2osno > 1000.) then
        clm2_struc(n)%clm(k)%do_capsnow = .true.
     else
        clm2_struc(n)%clm(k)%do_capsnow = .false.
     endif

     if (.not. clm2_struc(n)%clm(k)%lakpoi) then

        do j = clm2_struc(n)%clm(k)%snl+1, 0       ! ice fraction of snow at previous time step
           clm2_struc(n)%clm(k)%frac_iceold(j) =  &
           clm2_struc(n)%clm(k)%h2osoi_ice(j)/    &
           (clm2_struc(n)%clm(k)%h2osoi_liq(j)+clm2_struc(n)%clm(k)%h2osoi_ice(j))
        enddo
!
! Determine beginning water balance (water balance at previous time step)
!
        clm2_struc(n)%clm(k)%begwb = clm2_struc(n)%clm(k)%h2ocan + &
                                    clm2_struc(n)%clm(k)%h2osno
        do j = 1, nlevsoi
           clm2_struc(n)%clm(k)%begwb = clm2_struc(n)%clm(k)%begwb +         &
                                       clm2_struc(n)%clm(k)%h2osoi_ice(j) + &
                                       clm2_struc(n)%clm(k)%h2osoi_liq(j)

        enddo
! Determine canopy interception and precipitation onto ground surface.
! Determine the fraction of foliage covered by water and the fraction
! of foliage that is dry and transpiring. Initialize snow layer if the
! snow accumulation exceeds 10 mm.
        call Hydrology1(clm2_struc(n)%clm(k))
! Determine leaf temperature and surface fluxes based on ground
! temperature from previous time step.
!

        call Biogeophysics1(clm2_struc(n)%clm(k),cgrnd,cgrndl,cgrnds,tg,emg,htvp, dlrad,ulrad,tssbef)
     else if (clm2_struc(n)%clm(k)%lakpoi) then
!
! Determine lake temperature and surface fluxes
!
        call Biogeophysics_Lake (clm2_struc(n)%clm(k))

     endif
 
 44   format(i6,1x,i4,1x,i4,1x,i2,1x,f6.2,1x,4(f8.4,1x),f12.8,1x,f12.8)
     if (.not. clm2_struc(n)%clm(k)%lakpoi) then
!
! Surface biogeochemical fluxes: co2 respiration and plant production
!
#if (defined BGC)
        call Biogeochemistry (clm2_struc(n)%clm(k))
#endif
!
! Ecosystem dynamics: phenology, vegetation, soil carbon.
! Also updates snow fraction
!
        call EcosystemDyn (clm2_struc(n)%clm(k), doalb, .false.)
     else if (clm2_struc(n)%clm(k)%lakpoi) then


     endif
!
! Albedos for next time step
!

     if (doalb) then
        call SurfaceAlbedo (clm2_struc(n)%clm(k), caldayp1, eccen, obliqr, lambm0, mvelpp)
     endif
!
! THIS WILL EVENTUALLY MARK THE END OF THE PATCH LOOP AND
! THE BEGINNING OF THE SINGLE COLUMN SOIL LOOP(S)
!
! Determine soil/snow temperatures including ground temperature and
! update surface fluxes for new ground temperature.
!
     if (.not. clm2_struc(n)%clm(k)%lakpoi) then
        call Biogeophysics2(clm2_struc(n)%clm(k),cgrnd,cgrndl,cgrnds,tg,emg,htvp, dlrad,ulrad,tssbef)
     endif
  end do     

! ----------------------------------------------------------------------
! LOOP 2
! ----------------------------------------------------------------------

  do k = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)   ! begin 2nd loop over patches
! Vertical (column) soil and surface hydrology
!
     if (.not. clm2_struc(n)%clm(k)%lakpoi) call Hydrology2 (clm2_struc(n)%clm(k))
     
! Lake hydrology
!
     if (clm2_struc(n)%clm(k)%lakpoi) call Hydrology_Lake (clm2_struc(n)%clm(k))
     
! Update Snow Age (needed for surface albedo calculation - but is
! really a column type property
!
     call SnowAge (clm2_struc(n)%clm(k))
     
! Fraction of soil covered by snow - really a column property
!
     clm2_struc(n)%clm(k)%frac_sno = clm2_struc(n)%clm(k)%snowdp/(0.1 + clm2_struc(n)%clm(k)%snowdp)
     
!
! Check the energy and water balance
!
     call BalanceCheck (clm2_struc(n)%clm(k))
  end do    

! ----------------------------------------------------------------------
! Write global average diagnostics to standard output
! ----------------------------------------------------------------------
  do k=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_SWNET,&
          value=clm2_struc(n)%clm(k)%fsa,&
          vlevel=1,unit="W m-2", direction="DN",&
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_LWNET,value=&
          ((-1.0)*clm2_struc(n)%clm(k)%eflx_lwrad_net),vlevel=1,unit="W m-2",&
          direction="DN",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_QLE,&
          value=clm2_struc(n)%clm(k)%eflx_lh_tot,vlevel=1,unit="W m-2",&
          direction="UP",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_QH,&
          value=clm2_struc(n)%clm(k)%eflx_sh_tot,&
          vlevel=1,unit="W m-2", direction="UP",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_QG,&
          value=clm2_struc(n)%clm(k)%eflx_soil_grnd,vlevel=1,unit="W m-2",&
          direction="DN",surface_type=LIS_rc%lsm_index)

!Bowen Ratio - sensible/latent
     if(clm2_struc(n)%clm(k)%eflx_lh_tot.gt.0) then 
        call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_BR,&
          value=clm2_struc(n)%clm(k)%eflx_sh_tot/clm2_struc(n)%clm(k)%eflx_lh_tot, &
          vlevel=1,unit="-", direction="-",surface_type=LIS_rc%lsm_index)
     else
        call LIS_diagnoseSurfaceOutputVar(n,k,LIS_MOC_BR,value=0.0,vlevel=1,unit="-",&
             direction="-",surface_type=LIS_rc%lsm_index)
     endif

!Evaporative Fraction
     if( (clm2_struc(n)%clm(k)%eflx_lh_tot + &
          clm2_struc(n)%clm(k)%eflx_sh_tot ) .ne. 0 ) then 
        call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_EF,&
          value=clm2_struc(n)%clm(k)%eflx_lh_tot/(clm2_struc(n)%clm(k)%eflx_lh_tot + &
          clm2_struc(n)%clm(k)%eflx_sh_tot), vlevel=1,unit="-",&
          direction="-",surface_type=LIS_rc%lsm_index)
     else
!double check
        call LIS_diagnoseSurfaceOutputVar(n,k,LIS_MOC_EF,value=1.0,vlevel=1,unit="-",&
             direction="-",surface_type=LIS_rc%lsm_index)
     endif

     call LIS_diagnoseSurfaceOutputVar(n,k,LIS_MOC_TOTALPRECIP,&
          value=clm2_struc(n)%clm(k)%forc_snow + clm2_struc(n)%clm(k)%forc_rain,&
          vlevel=1,unit="kg m-2 s-1", direction="DN",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,k,LIS_MOC_TOTALPRECIP,&
          value=(clm2_struc(n)%clm(k)%forc_snow + &
          clm2_struc(n)%clm(k)%forc_rain)*LIS_rc%nts(n),&
          vlevel=1,unit="kg m-2", direction="DN") ! EMK
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_SNOWF,&
          value=clm2_struc(n)%clm(k)%forc_snow,&
          vlevel=1,unit="kg m-2 s-1", direction="DN",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_RAINF,&
          value=clm2_struc(n)%clm(k)%forc_rain,&
          vlevel=1,unit="kg m-2 s-1", direction="DN",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_EVAP,&
          value=clm2_struc(n)%clm(k)%qflx_evap_tot,&
          vlevel=1,unit="kg m-2 s-1", direction="UP",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_QS,&
          value=(clm2_struc(n)%clm(k)%qflx_surf+&
          clm2_struc(n)%clm(k)%qflx_qrgwl), vlevel=1,unit="kg m-2 s-1", direction="OUT",surface_type=LIS_rc%lsm_index)

     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_QS,&
          value=(clm2_struc(n)%clm(k)%qflx_surf+&
          clm2_struc(n)%clm(k)%qflx_qrgwl)*LIS_rc%nts(n), &
          vlevel=1,unit="kg m-2", direction="OUT",surface_type=LIS_rc%lsm_index) ! EMK
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_QSB,&
          value=clm2_struc(n)%clm(k)%qflx_drain,&
          vlevel=1,unit="kg m-2 s-1", direction="OUT",surface_type=LIS_rc%lsm_index)

     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_QSB,&
          value=clm2_struc(n)%clm(k)%qflx_drain*LIS_rc%nts(n),&
          vlevel=1,unit="kg m-2", direction="OUT",surface_type=LIS_rc%lsm_index) ! EMK
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_QSM,&
          value=clm2_struc(n)%clm(k)%qflx_snomelt,&
          vlevel=1,unit="kg m-2 s-1", direction="S2L",surface_type=LIS_rc%lsm_index)

     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_QSM,&
          value=clm2_struc(n)%clm(k)%qflx_snomelt*LIS_rc%nts(n),&
          vlevel=1,unit="kg m-2", direction="S2L",surface_type=LIS_rc%lsm_index) ! EMK
     soimtc = 0.0
     soilmr = 0.0
     roottemp = 0.0
     do m=1,nlevsoi
        vol_ice(m) = min(clm2_struc(n)%clm(k)%watsat(m),&
             clm2_struc(n)%clm(k)%h2osoi_ice(m)/&
             (clm2_struc(n)%clm(k)%dz(m)*denice))
        eff_porosity(m) = clm2_struc(n)%clm(k)%watsat(m)-vol_ice(m)
        vol_liq(m) = min(eff_porosity(m), &
             clm2_struc(n)%clm(k)%h2osoi_liq(m)/&
             (clm2_struc(n)%clm(k)%dz(m)*denh2o))       
        call LIS_diagnoseSurfaceOutputVar(n,k,LIS_MOC_SOILMOIST,&
             value=(vol_liq(m)+vol_ice(m)),&
             vlevel = m, unit="m^3 m-3", direction="-",surface_type=LIS_rc%lsm_index)
        soimtc = soimtc + vol_liq(m)+vol_ice(m)
        if(clm2_struc(n)%clm(k)%lakpoi) then 
           tempvar = clm2_struc(n)%clm(k)%t_lake(m)
        else
           tempvar = clm2_struc(n)%clm(k)%t_soisno(m)
        endif
        call LIS_diagnoseSurfaceOutputVar(n,k,LIS_MOC_SOILTEMP, value=tempvar, &
             vlevel=m, unit="K",direction="-",surface_type=LIS_rc%lsm_index)
        
        if(m.le.8) then 
           if(m.eq.8) then 
              soilmr = soilmr +(vol_liq(m)+vol_ice(m))*0.17110706
           else
              soilmr=soilmr + (vol_liq(m)+vol_ice(m))*clm2_struc(n)%clm(k)%dz(m)
           endif
        endif
        
        if(m.le.8) then 
           if(m.eq.8) then 
              roottemp = roottemp + tempvar*0.17110706
           else
              roottemp=roottemp + tempvar*clm2_struc(n)%clm(k)%dz(m)
           endif
        endif
        
     enddo
   
! rootmoist = 1m soil moisture
     call LIS_diagnoseSurfaceOutputVar(n,k,LIS_MOC_ROOTMOIST, &
          value=1000.0*soilmr, vlevel=1,&
          unit="kg m-2",direction="-",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,k,LIS_MOC_ROOTTEMP, &
          value=roottemp, vlevel=1,&
          unit="K",direction="-",surface_type=LIS_rc%lsm_index)

     call LIS_diagnoseSurfaceOutputVar(n,k,LIS_MOC_DELSOILMOIST, value=&
          clm2_struc(n)%clm(k)%soilmtc_prev-soimtc, vlevel=1, unit="kg m-2",&
          direction="INC",surface_type=LIS_rc%lsm_index)

     snowtemp = 0 
     if (clm2_struc(n)%clm(k)%itypwat/=istwet)then
        if(clm2_struc(n)%clm(k)%snl < 0)then
           totaldepth=0.
           do i=clm2_struc(n)%clm(k)%snl+1,0    ! Compute total depth of snow layers
              totaldepth=totaldepth+clm2_struc(n)%clm(k)%dz(i)
           enddo
           
           do i=clm2_struc(n)%clm(k)%snl+1,0    ! Compute snow temperature
              snowtemp=snowtemp+&
                   (clm2_struc(n)%clm(k)%t_soisno(i)*clm2_struc(n)%clm(k)%dz(i))
           enddo
           snowtemp=snowtemp/totaldepth
        endif
        if(snowtemp.eq.0)snowtemp=LIS_rc%udef
     endif
     call LIS_diagnoseSurfaceOutputVar(n, k, LIS_MOC_SNOWT, &
          value=snowtemp, vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k, LIS_MOC_VEGT, &
          value=clm2_struc(n)%clm(k)%t_veg,&
          vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k, LIS_MOC_BARESOILT, &
          value=clm2_struc(n)%clm(k)%t_grnd,&
          vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
!----------------------------------------------------------------------  
! AvgSurfT is the average surface temperature which depends on
! the snow temperature, bare soil temperature and canopy temperature
!----------------------------------------------------------------------
     snowt=0.
     if (clm2_struc(n)%clm(k)%itypwat/=istwet)then 
        if(clm2_struc(n)%clm(k)%snl < 0)then
           snowt =clm2_struc(n)%clm(k)%t_soisno(clm2_struc(n)%clm(k)%snl+1)
        endif
     endif
     if(snowt ==0.)snowt =LIS_rc%udef  
     if(snowt.ne.LIS_rc%udef)then
        asurf=clm2_struc(n)%clm(k)%frac_sno*snowt+ & 
             clm2_struc(n)%clm(k)%frac_veg_nosno*clm2_struc(n)%clm(k)%t_veg+  & 
             (1-(clm2_struc(n)%clm(k)%frac_sno+clm2_struc(n)%clm(k)%frac_veg_nosno))* & 
             clm2_struc(n)%clm(k)%t_grnd
     else
        asurf=clm2_struc(n)%clm(k)%frac_veg_nosno*clm2_struc(n)%clm(k)%t_veg+ & 
             (1-clm2_struc(n)%clm(k)%frac_veg_nosno)*clm2_struc(n)%clm(k)%t_grnd
     endif

     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_AVGSURFT, &
          value=asurf,vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_RADT, &
          value=clm2_struc(n)%clm(k)%t_rad,&
          vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
     
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_ALBEDO, &
          value=clm2_struc(n)%clm(k)%surfalb,&
          vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)

     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_ECANOP, value=&
          (clm2_struc(n)%clm(k)%qflx_evap_veg-clm2_struc(n)%clm(k)%qflx_tran_veg),&
          vlevel=1,unit="kg m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_TVEG, value=&
          clm2_struc(n)%clm(k)%qflx_tran_veg,vlevel=1,unit="kg m-2 s-1",&
          direction="UP",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_ESOIL, value=&
          clm2_struc(n)%clm(k)%qflx_evap_veg,vlevel=1,unit="kg m-2 s-1",&
          direction="UP",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_SUBSNOW, value=&
          clm2_struc(n)%clm(k)%qflx_sub_snow, vlevel=1,unit="kg m-2 s-1",&
          direction="-",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_SWE, &
          value=clm2_struc(n)%clm(k)%h2osno,&
          vlevel=1,unit="kg m-2",direction="-",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_SWE, &
          value=clm2_struc(n)%clm(k)%h2osno/LIS_CONST_RHOFW,&
          vlevel=1,unit="m",direction="-",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_SNOWDEPTH, &
          value=clm2_struc(n)%clm(k)%snowdp,&
          vlevel=1,unit="m",direction="-",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_SNOWCOVER, vlevel=1,unit="-",&
          value=clm2_struc(n)%clm(k)%frac_sno,direction="-",&
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_CANOPINT, value=&
          clm2_struc(n)%clm(k)%h2ocan,vlevel=1,unit="kg m-2",&
          direction="-",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k,LIS_MOC_DELSWE,vlevel=1,unit="kg m-2",value=&
          clm2_struc(n)%clm(k)%h2osno-clm2_struc(n)%clm(k)%h2osno_prev,&
          direction="INC",surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n, k, LIS_MOC_ACOND,&
          value=clm2_struc(n)%clm(k)%acond,&
          vlevel=1,unit="m s-1",direction="-",surface_type=LIS_rc%lsm_index)
  enddo
  clm2_struc(n)%count=clm2_struc(n)%count+1
  soilmtc = 0.0

  if(LIS_rc%tscount(n)==0 .or. LIS_rc%tscount(n) ==1 &
       .or. LIS_rc%rstflag(n).eq.1 ) then 
     do m=1,nlevsoi 
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
           soilm(t,m)=clm2_struc(n)%clm(t)%h2osoi_liq(m)+&
                clm2_struc(n)%clm(t)%h2osoi_ice(m)
        enddo
     enddo
     do m=1,nlevsoi
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
           soilmtc(t)=soilmtc(t)+soilm(t,m)
        enddo
     enddo
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        clm2_struc(n)%clm(t)%soilmtc_prev = soilmtc(t)
        clm2_struc(n)%clm(t)%h2osno_prev = clm2_struc(n)%clm(t)%h2osno
     enddo
  endif

  do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     clm2_struc(n)%clm(k)%forc_t = 0
     clm2_struc(n)%clm(k)%forc_q  = 0 
     clm2_struc(n)%clm(k)%forc_solad = 0
     clm2_struc(n)%clm(k)%forc_solai = 0
     clm2_struc(n)%clm(k)%forc_lwrad = 0
     clm2_struc(n)%clm(k)%forc_u     = 0
     clm2_struc(n)%clm(k)%forc_v     = 0
     clm2_struc(n)%clm(k)%forc_pbot  = 0
     clm2_struc(n)%clm(k)%forc_rain  = 0
     clm2_struc(n)%clm(k)%forc_snow  = 0
  enddo
  clm2_struc(n)%forc_count = 0 
endif
  return
end subroutine clm2_main

