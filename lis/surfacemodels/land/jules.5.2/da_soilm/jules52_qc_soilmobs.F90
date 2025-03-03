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
! !ROUTINE: jules52_qc_soilmobs
! \label{jules52_qc_soilmobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
! 20 Dec 2018: Mahdi Navari; Modified for JULES 5.2
!
! !INTERFACE:
subroutine jules52_qc_soilmobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify, LIS_logunit, LIS_endrun
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_DAobservationsMod
  use jules52_lsmMod
  use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
  use jules_surface_mod,      only: l_aggregate

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the soil moisture observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  soil is frozen or (3) ground is fully or partially covered 
!  with snow. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_sm_field

  real, pointer            :: smobs(:)
  integer                  :: t,j,pft,l
  integer                  :: gid
  integer                  :: status
  real                     :: lat,lon
  real                     :: smc1(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! total soil moisture
  real                     :: smc2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc4(LIS_rc%npatch(n,LIS_rc%lsm_index))

  real                     :: sthu1(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! unfrozen soil moisture
  real                     :: sthu2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sthu3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sthu4(LIS_rc%npatch(n,LIS_rc%lsm_index))

  real                     :: sthf1(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! frozen soil moisture
  real                     :: sthf2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sthf3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sthf4(LIS_rc%npatch(n,LIS_rc%lsm_index))

  real                     :: stc1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc4(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: vegt(LIS_rc%npatch(n,LIS_rc%lsm_index))

! real 		 	   :: frac_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN  fveg
 real 		 	   :: lai_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN  
! real 		 	   :: nsnow_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN 

  real			   :: t_skin(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN
  real 		 	   :: fveg(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN  
  real 		 	   :: sneqv(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN
  real 		 	   :: fsno(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN 
  real			   :: BDSNO(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN
  real			   :: FMELT(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN
  real 		 	   :: p_s_smvcwt_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN
  real 		 	   :: p_s_smvcst_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN

  real                     :: rainf_obs(LIS_rc%obs_ngrid(k))
  real                     :: sneqv_obs(LIS_rc%obs_ngrid(k))
!  real                     :: nsnow_obs(LIS_rc%obs_ngrid(k))!sca_obs
  real                     :: fsno_obs(LIS_rc%obs_ngrid(k))!sca_obs
  real                     :: shdfac_obs(LIS_rc%obs_ngrid(k))
  real                     :: t1_obs(LIS_rc%obs_ngrid(k))
  real                     :: smcwlt_obs(LIS_rc%obs_ngrid(k))
  real                     :: smcmax_obs(LIS_rc%obs_ngrid(k))
!  real                     :: smc1_obs(LIS_rc%obs_ngrid(k))
!  real                     :: smc2_obs(LIS_rc%obs_ngrid(k))
!  real                     :: smc3_obs(LIS_rc%obs_ngrid(k))
!  real                     :: smc4_obs(LIS_rc%obs_ngrid(k))
!  real                     :: sthu1_obs(LIS_rc%obs_ngrid(k))
!  real                     :: sthu2_obs(LIS_rc%obs_ngrid(k))
!  real                     :: sthu3_obs(LIS_rc%obs_ngrid(k))
!  real                     :: sthu4_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthf1_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthf2_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthf3_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthf4_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc1_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc2_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc3_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc4_obs(LIS_rc%obs_ngrid(k))
  real                     :: vegt_obs(LIS_rc%obs_ngrid(k))

  call ESMF_StateGet(OBS_State,"Observation01",obs_sm_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in jules52_qc_soilmobs")
  call ESMF_FieldGet(obs_sm_field,localDE=0,farrayPtr=smobs,rc=status) ! [m3/m3] ?
  call LIS_verify(status,& 
       "ESMF_FieldGet failed in jules52_qc_soilmobs")

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
!     In Jules unfrozn soil moisture content are stored as an independent variable
!     so we can use that variable instead of computing frozn soil moisture using smc and sthu?
#if 0 
     smc1(t) = jules52_struc(n)%jules52(t)%smcl_soilt(1) ! [kg/m2]
     smc2(t) = jules52_struc(n)%jules52(t)%smcl_soilt(2)
     smc3(t) = jules52_struc(n)%jules52(t)%smcl_soilt(3)
     smc4(t) = jules52_struc(n)%jules52(t)%smcl_soilt(4)
     ! Unfrozen soil moisture
     sthu1(t) = jules52_struc(n)%jules52(t)%p_s_sthu(1) ![-]    Noah -->sh2o
     sthu2(t) = jules52_struc(n)%jules52(t)%p_s_sthu(2)
     sthu3(t) = jules52_struc(n)%jules52(t)%p_s_sthu(3)
     sthu4(t) = jules52_struc(n)%jules52(t)%p_s_sthu(4)
#endif

     ! Frozen soil moisture
     sthf1(t) = jules52_struc(n)%jules52(t)%p_s_sthf(1) ![-]    
     sthf2(t) = jules52_struc(n)%jules52(t)%p_s_sthf(2)
     sthf3(t) = jules52_struc(n)%jules52(t)%p_s_sthf(3)
     sthf4(t) = jules52_struc(n)%jules52(t)%p_s_sthf(4)

     stc1(t) = jules52_struc(n)%jules52(t)%t_soil(1) ! stc
     stc2(t) = jules52_struc(n)%jules52(t)%t_soil(2)
     stc3(t) = jules52_struc(n)%jules52(t)%t_soil(3)
     stc4(t) = jules52_struc(n)%jules52(t)%t_soil(4)

  enddo



fsno = 0 
fveg = 0
sneqv = 0
  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
      !pft= int(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt)   
      pft= jules52_struc(n)%jules52(t)%pft
      vegt(t) = pft  
      sneqv(t) = jules52_struc(n)%jules52(t)%snowdepth(pft) * &
                     jules52_struc(n)%jules52(t)%rho_snow_grnd(pft) ![kg m-2]

      !! MN : Note : compute snow cover fraction 
      !! module_sf_noahmplsm_36.F90 
      !! ground snow cover fraction [Niu and Yang, 2007, JGR]
      !! REAL, PARAMETER :: M  = 2.50   ! melting factor (-)
      !! REAL, PARAMETER :: Z0 = 0.01   ! Bare-soil roughness length (m) (i.e., under the canopy)
      !! FSNO = 0.
      !! IF(SNOWH.GT.0.)  THEN
      !!   BDSNO    = SNEQV / SNOWH ! SNOWH[mm], SNEQV[mm]
      !!   FMELT    = (BDSNO/100.)**M
      !!   FSNO     = TANH( SNOWH /(2.5* Z0 * FMELT))
      !! ENDIF
      !if (jules52_struc(n)%jules52(t)%snowdepth(pft).GT.0.) then ! m
      !   BDSNO(t)    = sneqv(t) / jules52_struc(n)%jules52(t)%snowdepth(pft) * 100 !kg/m2(=mm)/(m*100) 
      !   FMELT(t)    = (BDSNO(t)/100.)**2.50 !M = 2.50
      !   fsno(t) = TANH( jules52_struc(n)%jules52(t)%snowdepth(pft) * 100 /(2.5* 0.01 * FMELT(t))) ! Z0 = 0.01 , convert the snow depth from m to mm 
      !endif

      !shugong (5/17/2018) suggest to use this to be consistant with the main routine.      
      !snow frac according to the algorithm in JULES output interface 
      !snow_frac = 0.0

       if(l_aggregate) then
         if(jules52_struc(n)%jules52(t)%snow_tile(1)+jules52_struc(n)%jules52(t)%snow_grnd(1)>1.0) then
            fsno(t) = 1.0
         endif
       else
         if(jules52_struc(n)%jules52(t)%snow_tile(pft)+jules52_struc(n)%jules52(t)%snow_grnd(pft)>1.0) then
            fsno(t) = jules52_struc(n)%jules52(t)%frac(pft)
         endif
       endif
      t_skin(t) = jules52_struc(n)%jules52(t)%tstar_tile(pft)


      !MN:fveg is used for 12-month green vegetation fraction (i.e., noah33_struc(n)%noah(:)%shdfac)  
      !print*,'l_aggregate', l_aggregate
      if(.NOT. l_aggregate) then   
         write(LIS_logunit,*) 'Please set the l_aggregate to .true. in the jules_surface.nml ' 
         call LIS_endrun
      else
       do l=1,5 !Broadleaf trees, Needleleaf trees, C3 (temperate) grass, C4 (tropical) grass, Shrubs
         fveg(t) = fveg(t) + jules52_struc(n)%jules52(t)%surft_frac(l)   
         !print*,'t, l',t, l,jules52_struc(n)%jules52(t)%surft_frac(l),fveg(t)
       enddo 
      endif          
      !print*,'t, l',t, l,jules52_struc(n)%jules52(t)%surft_frac(l),fveg(t)  
      !print*,''

#if 0 
!       !print*, jules52_struc(n)%ntype
!       !print*, LIS_rc%npatch(n,LIS_rc%lsm_index)
!       !print*, 'l_aggregate', l_aggregate
!       do l=1,5 !Broadleaf trees, Needleleaf trees, C3 (temperate) grass, C4 (tropical) grass, Shrubs
!         !print*, 't,l', t,l
!         !fveg(t) = fveg(t) + jules52_struc(n)%jules52(t)%frac(l)  
!         !print*,'frac', jules52_struc(n)%jules52(t)%frac(l) 
!         fveg(t) = fveg(t) + jules52_struc(n)%jules52(t)%surft_frac(l)   
!         print*,'t, l',t, l,jules52_struc(n)%jules52(t)%surft_frac(l),fveg(t)
!       enddo         
#endif 

      !!frac_tmp(t)= jules52_struc(n)%jules52(t)%frac(pft)
      !frac_tmp(t)= 1-EXP(-0.5*jules52_struc(n)%jules52(t)%lai(pft))
      lai_tmp(t)=jules52_struc(n)%jules52(t)%lai(pft)
      !sliq_tmp(t)= sum(jules52_struc(n)%jules52(t)%sliq(pft,:))
      !nsnow_tmp(t)= jules52_struc(n)%jules52(t)%nsnow(pft)
      p_s_smvcwt_tmp(t)= jules52_struc(n)%jules52(t)%p_s_smvcwt(1) ! m3/m3
      p_s_smvcst_tmp(t)= jules52_struc(n)%jules52(t)%p_s_smvcst(1) ! m3/m3
 enddo


  call LIS_convertPatchSpaceToObsSpace(n,k,&       
       LIS_rc%lsm_index, &
       jules52_struc(n)%jules52(:)%rainf,&
       rainf_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sneqv,&
       sneqv_obs)  ! sneqv
!  call LIS_convertPatchSpaceToObsSpace(n,k,&
!       LIS_rc%lsm_index, &
!       nsnow_tmp ,&
!       nsnow_obs) ! noah33_struc(n)%noah(:)%sca 
	! MN: jules does not store the snow cover fraction. We use this 
	! variable to determine the presence of snow on the ground. I think  
	! we can use any variable that shows the presence of snow on ground
	! here number of snow layer was replaced with "sca"
!#if 0
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       fsno,&
       fsno_obs)! noah33_struc(n)%noah(:)%sca 
!#endif 
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       fveg,&
       shdfac_obs) ! noah33_struc(n)%noah(:)%shdfac
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       t_skin,&
       t1_obs) ! noah33_struc(n)%noah(:)%t1
! tstar_tile is a function of pft (plant functional types)   
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       p_s_smvcst_tmp,&
       smcmax_obs) !smcmax Here I assumed p_s_smvcst for all layers are the same
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       p_s_smvcwt_tmp,&
       smcwlt_obs) !smcwlt  Here I assumed p_s_smvcwt for all layers are the same
#if 0 
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc1,&
       smc1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc2,&
       smc2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc3,&
       smc3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc4,&
       smc4_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthu1,&
       sthu1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthu2,&
       sthu2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthu3,&
       sthu3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthu4,&
       sthu4_obs)
#endif 

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthf1,&
       sthf1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthf2,&
       sthf2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthf3,&
       sthf3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthf4,&
       sthf4_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc1,&
       stc1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc2,&
       stc2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc3,&
       stc3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc4,&
       stc4_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       vegt,&
       vegt_obs)


  do t = 1,LIS_rc%obs_ngrid(k)
     if(smobs(t).ne.LIS_rc%udef) then 
        if(rainf_obs(t).gt.3E-6) then 
           smobs(t) = LIS_rc%udef
! MN: use sthf(1, :) Frozen soil moisture content of layers as a 
! fraction of saturation (instead of smc4_obs(t)- sthu4_obs(t))
!        elseif(abs(smc4_obs(t)- &
!             sthu4_obs(t)).gt.0.0001) then
!           smobs(t) = LIS_rc%udef

 
        elseif(sthf1_obs(t).gt.0.0001) then 
           smobs(t) = LIS_rc%udef 
        elseif(sthf2_obs(t).gt.0.0001) then
           smobs(t) = LIS_rc%udef
        elseif(sthf3_obs(t).gt.0.0001) then
           smobs(t) = LIS_rc%udef
        elseif(sthf4_obs(t).gt.0.0001) then
           smobs(t) = LIS_rc%udef
! layer temperature --> jules52_struc(n)%jules52(t)%t_soil(i) 
        elseif(stc1_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(stc2_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(stc3_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(stc4_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
! MN: skin temperature 
! jules52_struc(n)%jules52(t)%tstar_tile(pft)
        elseif(t1_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
 


!MN
!I got this from namelist (jules_surface_types.nml)
!ice=9,
!lake=7,
!nnvg=4,
!npft=5,
!soil=8,
!urban=6,)
! I got this from website (comparing the number with the following surface type makes sense. 
!Five Plant Functional Types (PFTs)
!Broadleaf trees
!Needle leaf trees
!C3 (temperate) grass
!C4 (tropical) grass
!Shrubs
!Four non-vegetation types
!Urban
!Inland water
!Bare soil
!Land-ice
!for Noah we just compare against forest hear we have to consider the pft of 1,2,6,7,9 
        elseif(vegt_obs(t).eq.6) then !Urban
           smobs(t) = LIS_rc%udef
        elseif(vegt_obs(t).eq.7) then !Inland water
           smobs(t) = LIS_rc%udef
        elseif(vegt_obs(t).eq.9) then !Land-ice
           smobs(t) = LIS_rc%udef
        elseif(sneqv_obs(t).gt.0.001) then 
           smobs(t) = LIS_rc%udef
        elseif(fsno_obs(t).gt.0) then   
           smobs(t) = LIS_rc%udef
!too close to the tails, could be due to scaling, so reject. 
!MN Note: By the time this routine is called, the obs soil moisture has already been rescaled into the volumetric units.
        elseif(smcmax_obs(t)-smobs(t).lt.0.02) then 
           smobs(t) = LIS_rc%udef
!#if 0
        elseif(shdfac_obs(t).gt.0.7) then ! vegetation fraction 
           smobs(t) = LIS_rc%udef    
!#endif  
!In some soil types wilting point is very high e.g. 0.237 m3/m3
        elseif(smobs(t) - smcwlt_obs(t).lt.0.02) then  ! changed from 0.02 to ... 
            smobs(t) = LIS_rc%udef

        endif
     endif
  enddo

end subroutine jules52_qc_soilmobs

