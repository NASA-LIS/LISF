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
! !ROUTINE: noahmp36_settws
!  \label{noahmp36_settws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp36_settws(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp36_lsmMod
  use module_sf_noahlsm_36

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space.
!
!EOP
  real, parameter        :: MIN_GWS_THRESHOLD = 0.00
  real, parameter        :: MAX_GWS_THRESHOLD = 7000.0
  real, parameter        :: MAX_WA = 7000.0
  real, parameter        :: ZSOIL = 2 !mm
  real, parameter        :: ROUS = 0.2 ! specific yield
  !Bailing changed this to be WLTSMC
  real, parameter        :: MIN_THRESHOLD = 0.02
!  real                   :: MIN_THRESHOLD
  real                   :: MAX_THRESHOLD
  real                   :: sm_threshold
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  type(ESMF_Field)       :: gwField
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField

  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real, pointer          :: gws(:)
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  integer                :: t, j,i, gid, m, t_unpert
  integer                :: status
  real                   :: delta(5)
  real                   :: delta1,delta2,delta3,delta4,delta5
  real                   :: tmpval
  logical                :: bounds_violation
  integer                :: nIter
  logical                :: update_flag(LIS_rc%ngrid(n))
  logical                :: ens_flag(LIS_rc%nensem(n))
! mn
  integer                :: SOILTYP           ! soil type index [-]
  real                   :: SMCMAX , SMCWLT
  real                   :: tmp(LIS_rc%nensem(n)), tmp0(LIS_rc%nensem(n))
  real                   :: tmp1(LIS_rc%nensem(n)),tmp2(LIS_rc%nensem(n)), &
       tmp3(LIS_rc%nensem(n)),tmp4(LIS_rc%nensem(n))
  real                   :: tmp5(LIS_rc%nensem(n))
  logical                :: update_flag_tile(LIS_rc%npatch(n,LIS_rc%lsm_index))
  logical                :: flag_ens(LIS_rc%ngrid(n))
  logical                :: flag_tmp(LIS_rc%nensem(n))
  logical                :: update_flag_ens(LIS_rc%ngrid(n))
  logical                :: update_flag_new(LIS_rc%ngrid(n))
  integer                :: RESULT, icount
  real                   :: MaxEnsSM1 ,MaxEnsSM2 ,MaxEnsSM3 ,MaxEnsSM4, &
       MaxEnsGWS !Wanshu
  real                   :: MinEnsSM1 ,MinEnsSM2 ,MinEnsSM3 ,MinEnsSM4, &
       MinEnsGWS !Wanshu
  real                   :: MaxEns_sh2o1, MaxEns_sh2o2, MaxEns_sh2o3, &
       MaxEns_sh2o4
  real                   :: smc_rnd, smc_tmp, gws_tmp !Wanshu
  real                   :: sh2o_tmp, sh2o_rnd
  real                   :: dsneqv,dsnowh,swe_old, snowh_old

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in noahmp36_settws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 2 failed in noahmp36_settws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 3 failed in noahmp36_settws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 4 failed in noahmp36_settws")
  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp36_settws")
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: SWE failed in noahmp36_settws")
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Snowdepth failed in noahmp36_settws")


  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmp36_settws")
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in noahmp36_settws")
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in noahmp36_settws")
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in noahmp36_settws")
  call ESMF_FieldGet(gwField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Groundwater Storage failed in noahmp36_settws")
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: SWE failed in noahmp36_settws")
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Snowdepth failed in noahmp36_settws")


  update_flag = .true.
  update_flag_tile= .true.

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     ! MN: NOTE: SMCMAX and SMCWLT are not stored in the data structure but we
     !       can get module variables MAXSMC and WLTSMC from the
     !       module_sf_noahlsm_36
     SOILTYP = NOAHMP36_struc(n)%noahmp36(t)%soiltype
     MAX_THRESHOLD = MAXSMC (SOILTYP)
     !MIN_THRESHOLD = WLTSMC (SOILTYP)
     sm_threshold = MAXSMC (SOILTYP) - 0.02
     !sm_threshold = NOAHMP36_struc(n)%noahmp36(t)%smcmax - 0.02
     !MAX_THRESHOLD = NOAHMP36_struc(n)%noahmp36(t)%smcmax

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     !MN: delta = X(+) - X(-)
     !NOTE: "noahmp36_updatesoilm.F90" updates the soilm_(t)
     delta1 = soilm1(t)-NOAHMP36_struc(n)%noahmp36(t)%smc(1)
     delta2 = soilm2(t)-NOAHMP36_struc(n)%noahmp36(t)%smc(2)
     delta3 = soilm3(t)-NOAHMP36_struc(n)%noahmp36(t)%smc(3)
     delta4 = soilm4(t)-NOAHMP36_struc(n)%noahmp36(t)%smc(4)

     !Wanshu
     delta5 = gws(t)-NOAHMP36_struc(n)%noahmp36(t)%wa


     ! MN: check  MIN_THRESHOLD < volumetric liquid soil moisture < threshold
     if(soilm1(t).ge.MIN_THRESHOLD .and. soilm1(t).lt. sm_threshold) then
        update_flag(gid) = update_flag(gid).and.(.true.)
        ! MN save the flag for each tile (col*row*ens)   (64*44)*20
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        write(LIS_logunit,*) 'problem updating sm1 ', soilm1(t), &
             NOAHMP36_struc(n)%noahmp36(t)%smc(1), &
             NOAHMP36_struc(n)%noahmp36(t)%sh2o(1)
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     if(soilm2(t).ge.MIN_THRESHOLD .and. soilm2(t).lt.sm_threshold) then
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
        write(LIS_logunit,*) 'problem updating sm2 ', soilm2(t), &
             NOAHMP36_struc(n)%noahmp36(t)%smc(2), &
             NOAHMP36_struc(n)%noahmp36(t)%sh2o(2),WLTSMC(SOILTYP)
     endif
     if(soilm3(t).ge.MIN_THRESHOLD .and. soilm3(t).lt.sm_threshold) then
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
        write(LIS_logunit,*) 'problem updating sm3 ', soilm3(t), &
             NOAHMP36_struc(n)%noahmp36(t)%smc(3), &
             NOAHMP36_struc(n)%noahmp36(t)%sh2o(3),WLTSMC(SOILTYP)
     endif
     if(soilm4(t).ge.MIN_THRESHOLD .and. soilm4(t).lt.sm_threshold) then
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
        write(LIS_logunit,*) 'problem updating sm4 ', soilm4(t), &
             NOAHMP36_struc(n)%noahmp36(t)%smc(4), &
             NOAHMP36_struc(n)%noahmp36(t)%sh2o(4),WLTSMC(SOILTYP)
     endif
     !     noahmp36_struc(n)%noahmp36(t)%wa = gws(t)
     !Wanshu
     if(gws(t).gt. MIN_GWS_THRESHOLD .and. gws(t).lt.MAX_GWS_THRESHOLD) then
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
        !-------
     endif
  enddo

!-----------------------------------------------------------------------------
! MN create new falg: if update falg for 50% of the ensemble members is true
! then update the state variabels
!-----------------------------------------------------------------------------
  update_flag_ens = .True.
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
     gid =i
     flag_tmp=update_flag_tile((i-1)*LIS_rc%nensem(n)+1:(i)*LIS_rc%nensem(n))
     RESULT = COUNT(flag_tmp)
     if (RESULT.lt.LIS_rc%nensem(n)*0.5) then   !>10/20*100 = 50%
        update_flag_ens(gid)= .False.
     endif
     update_flag_new(gid)= update_flag(gid).or.update_flag_ens(gid)  ! new flag
  enddo

  ! update step
  ! loop over grid points
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)

     gid =i

     !if(update_flag(gid)) then
     if(update_flag_new(gid)) then
        !----------------------------------------------------------------------
        ! 1- update the states
        ! 1-1- if update flag for tile is TRUE --> apply the DA update
        ! 1-2- if update flag for tile is FALSE --> resample the states
        ! 2- adjust the states
        !---------------------------------------------------------------------
        ! store update value for  cases that update_flag_tile &
        ! update_flag_new are TRUE
        ! update_flag_tile = TRUE --> means met the both min and max threshold

        tmp1 = LIS_rc%udef
        tmp2 = LIS_rc%udef
        tmp3 = LIS_rc%udef
        tmp4 = LIS_rc%udef
        !Wanshu
        tmp5 = LIS_rc%udef

        do m=1,LIS_rc%nensem(n)
           t = (i-1)*LIS_rc%nensem(n)+m

           if(update_flag_tile(t)) then

              tmp1(m) = soilm1(t)
              tmp2(m) = soilm2(t)
              tmp3(m) = soilm3(t)
              tmp4(m) = soilm4(t)
              tmp5(m) = gws(t)

           endif
        enddo

        MaxEnsSM1 = -10000
        MaxEnsSM2 = -10000
        MaxEnsSM3 = -10000
        MaxEnsSM4 = -10000
        MaxEnsGWS = -10000  !Wanshu

        MinEnsSM1 = 10000
        MinEnsSM2 = 10000
        MinEnsSM3 = 10000
        MinEnsSM4 = 10000
        MinEnsGWS = 10000  !Wanshu


        do m=1,LIS_rc%nensem(n)
           if(tmp1(m).ne.LIS_rc%udef) then
              MaxEnsSM1 = max(MaxEnsSM1, tmp1(m))
              MaxEnsSM2 = max(MaxEnsSM2, tmp2(m))
              MaxEnsSM3 = max(MaxEnsSM3, tmp3(m))
              MaxEnsSM4 = max(MaxEnsSM4, tmp4(m))
              MaxEnsGWS = max(MaxEnsGWS, tmp5(m)) !Wanshu

              MinEnsSM1 = min(MinEnsSM1, tmp1(m))
              MinEnsSM2 = min(MinEnsSM2, tmp2(m))
              MinEnsSM3 = min(MinEnsSM3, tmp3(m))
              MinEnsSM4 = min(MinEnsSM4, tmp4(m))
              MinEnsGWS = min(MinEnsGWS, tmp5(m)) !Wanshu

           endif
        enddo

        ! loop over tile
        do m=1,LIS_rc%nensem(n)
           t = (i-1)*LIS_rc%nensem(n)+m
           SOILTYP = NOAHMP36_struc(n)%noahmp36(t)%soiltype
           MAX_THRESHOLD = MAXSMC (SOILTYP)
           !MIN_THRESHOLD = WLTSMC (SOILTYP)
           sm_threshold = MAXSMC (SOILTYP) - 0.02

           ! MN check update status for each tile
           if(update_flag_tile(t)) then

              delta1 = soilm1(t)-NOAHMP36_struc(n)%noahmp36(t)%smc(1)
              delta2 = soilm2(t)-NOAHMP36_struc(n)%noahmp36(t)%smc(2)
              delta3 = soilm3(t)-NOAHMP36_struc(n)%noahmp36(t)%smc(3)
              delta4 = soilm4(t)-NOAHMP36_struc(n)%noahmp36(t)%smc(4)

              !Wanshu
              delta5=gws(t)-NOAHMP36_struc(n)%noahmp36(t)%wa

              NOAHMP36_struc(n)%noahmp36(t)%sh2o(1) =&
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(1)+&
                   delta1
              NOAHMP36_struc(n)%noahmp36(t)%smc(1) = soilm1(t)
              !Wanshu
              !              NOAHMP36_struc(n)%noahmp36(t)%wa = gws(t)

              NOAHMP36_struc(n)%noahmp36(t)%sh2o(2) = &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(2) + &
                   delta2
              NOAHMP36_struc(n)%noahmp36(t)%smc(2) = soilm2(t)
              ! MN: Test shutdown the update for 4th layer
              NOAHMP36_struc(n)%noahmp36(t)%sh2o(3) = &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(3) + &
                   delta3
              NOAHMP36_struc(n)%noahmp36(t)%smc(3) = soilm3(t)

              NOAHMP36_struc(n)%noahmp36(t)%sh2o(4) = &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(4) + &
                   delta4
              NOAHMP36_struc(n)%noahmp36(t)%smc(4) = soilm4(t)

              ! Bailing: add this check for liquid smc
              if(NOAHMP36_struc(n)%noahmp36(t)%sh2o(1)<MIN_THRESHOLD ) &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(1) = MIN_THRESHOLD
              if(NOAHMP36_struc(n)%noahmp36(t)%sh2o(1)>MAX_THRESHOLD ) &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(1) = MAX_THRESHOLD
              if(NOAHMP36_struc(n)%noahmp36(t)%sh2o(2)<MIN_THRESHOLD ) &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(2) = MIN_THRESHOLD
              if(NOAHMP36_struc(n)%noahmp36(t)%sh2o(2)>MAX_THRESHOLD ) &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(2) = MAX_THRESHOLD
              if(NOAHMP36_struc(n)%noahmp36(t)%sh2o(3)<MIN_THRESHOLD ) &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(3) = MIN_THRESHOLD
              if(NOAHMP36_struc(n)%noahmp36(t)%sh2o(3)>MAX_THRESHOLD ) &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(3) = MAX_THRESHOLD
              if(NOAHMP36_struc(n)%noahmp36(t)%sh2o(4)<MIN_THRESHOLD ) &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(4) = MIN_THRESHOLD
              if(NOAHMP36_struc(n)%noahmp36(t)%sh2o(4)>MAX_THRESHOLD ) &
                   NOAHMP36_struc(n)%noahmp36(t)%sh2o(4) = MAX_THRESHOLD
              !Wanshu
              if(NOAHMP36_struc(n)%noahmp36(t)%wa+delta5.gt. &
                   MIN_GWS_THRESHOLD.and.&
                   NOAHMP36_struc(n)%noahmp36(t)%wa+delta5.lt. &
                   MAX_GWS_THRESHOLD) then
                 NOAHMP36_struc(n)%noahmp36(t)%wa=gws(t)
              elseif(gws(t).lt.MIN_GWS_THRESHOLD .or. &
                   gws(t).gt.MAX_GWS_THRESHOLD) then
                 write(LIS_logunit,*) 'setgws',t,gws(t)
              endif
           else
              ! use mean value
              ! Assume sh2o = smc (i.e. ice content=0)
              smc_tmp = (MaxEnsSM1 - MinEnsSM1)/2 + MinEnsSM1
              NOAHMP36_struc(n)%noahmp36(t)%sh2o(1) = smc_tmp
              NOAHMP36_struc(n)%noahmp36(t)%smc(1) = smc_tmp

              smc_tmp = (MaxEnsSM2 - MinEnsSM2)/2 + MinEnsSM2
              NOAHMP36_struc(n)%noahmp36(t)%sh2o(2) = smc_tmp
              NOAHMP36_struc(n)%noahmp36(t)%smc(2) = smc_tmp

              smc_tmp = (MaxEnsSM3 - MinEnsSM3)/2 + MinEnsSM3
              NOAHMP36_struc(n)%noahmp36(t)%sh2o(3) = smc_tmp
              NOAHMP36_struc(n)%noahmp36(t)%smc(3) = smc_tmp

              smc_tmp = (MaxEnsSM4 - MinEnsSM4)/2 + MinEnsSM4
              NOAHMP36_struc(n)%noahmp36(t)%sh2o(4) = smc_tmp
              NOAHMP36_struc(n)%noahmp36(t)%smc(4) = smc_tmp

              gws_tmp = (MaxEnsGWS - MinEnsGWS)/2 + MinEnsGWS
              NOAHMP36_struc(n)%noahmp36(t)%wa = gws_tmp

           endif ! flag for each tile

        enddo ! loop over tile

     else ! if update_flag_new(gid) is FALSE
        if(LIS_rc%pert_bias_corr.eq.1) then
           !-------------------------------------------------------------------
           ! if no update is made, then we need to readjust the ensemble if
           ! pert bias correction is turned on because the forcing
           ! perturbations may cause biases to persist.
           !-------------------------------------------------------------------
           write(LIS_logunit,*) &
                '[WARNING] more than half of the ens violate the bounds',gid

           bounds_violation = .true.
           nIter = 0
           ens_flag = .true.

           do while(bounds_violation)
              niter = niter + 1
              t_unpert = i*LIS_rc%nensem(n)
              do j=1,4
                 delta(j) = 0.0
                 do m=1,LIS_rc%nensem(n)-1
                    t = (i-1)*LIS_rc%nensem(n)+m

                    if(m.ne.LIS_rc%nensem(n)) then
                       delta(j) = delta(j) + &
                            (NOAHMP36_struc(n)%noahmp36(t)%smc(j) - &
                            NOAHMP36_struc(n)%noahmp36(t_unpert)%smc(j))
                    endif

                 enddo
              enddo

              delta(5) = 0
              do m=1,LIS_rc%nensem(n)-1
                 t = (i-1)*LIS_rc%nensem(n)+m

                 if(m.ne.LIS_rc%nensem(n)) then
                    delta(5) = delta(5) + &
                         (NOAHMP36_struc(n)%noahmp36(t)%wa - &
                         NOAHMP36_struc(n)%noahmp36(t_unpert)%wa)
                 endif

              enddo

              do j=1,4
                 delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                    t = (i-1)*LIS_rc%nensem(n)+m
                    SOILTYP = NOAHMP36_struc(n)%noahmp36(t)%soiltype
                    MAX_THRESHOLD = MAXSMC (SOILTYP)
                    !MIN_THRESHOLD = WLTSMC (SOILTYP)
                    sm_threshold = MAXSMC (SOILTYP) - 0.02

                    tmpval = NOAHMP36_struc(n)%noahmp36(t)%smc(j) - &
                         delta(j)
                    if(tmpval.le.MIN_THRESHOLD) then
                       NOAHMP36_struc(n)%noahmp36(t)%sh2o(j) = &
                            max(NOAHMP36_struc(n)%noahmp36(t_unpert)%sh2o(j),&
                            MIN_THRESHOLD)
                       NOAHMP36_struc(n)%noahmp36(t)%smc(j) = &
                            max(NOAHMP36_struc(n)%noahmp36(t_unpert)%smc(j),&
                            MIN_THRESHOLD)
                       ens_flag(m) = .false.
                    elseif(tmpval.ge.sm_threshold) then
                       NOAHMP36_struc(n)%noahmp36(t)%sh2o(j) = &
                            min(NOAHMP36_struc(n)%noahmp36(t_unpert)%sh2o(j),&
                            sm_threshold)
                       NOAHMP36_struc(n)%noahmp36(t)%smc(j) = &
                            min(NOAHMP36_struc(n)%noahmp36(t_unpert)%smc(j),&
                            sm_threshold)
                       ens_flag(m) = .false.
                    endif
                 enddo
              enddo

              delta(5) =delta(5)/(LIS_rc%nensem(n)-1)
              do m=1,LIS_rc%nensem(n)-1
                 t = (i-1)*LIS_rc%nensem(n)+m

                 tmpval = NOAHMP36_struc(n)%noahmp36(t)%wa - &
                      delta(5)
                 if(tmpval.le.MIN_GWS_THRESHOLD) then
                    NOAHMP36_struc(n)%noahmp36(t)%wa = &
                         max(NOAHMP36_struc(n)%noahmp36(t_unpert)%wa,&
                         MIN_GWS_THRESHOLD)
                    ens_flag(m) = .false.
                 elseif(tmpval.ge.MAX_GWS_THRESHOLD) then
                    NOAHMP36_struc(n)%noahmp36(t)%wa = &
                         min(NOAHMP36_struc(n)%noahmp36(t_unpert)%wa,&
                         MAX_GWS_THRESHOLD)
                    ens_flag(m) = .false.
                 endif
              enddo

              !----------------------------------------------------------------
              ! Recalculate the deltas and adjust the ensemble
              !----------------------------------------------------------------
              do j=1,4
                 delta(j) = 0.0
                 do m=1,LIS_rc%nensem(n)-1
                    t = (i-1)*LIS_rc%nensem(n)+m
                    if(m.ne.LIS_rc%nensem(n)) then
                       delta(j) = delta(j) + &
                            (NOAHMP36_struc(n)%noahmp36(t)%smc(j) - &
                            NOAHMP36_struc(n)%noahmp36(t_unpert)%smc(j))
                    endif
                 enddo
              enddo

              delta(5) = 0.0
              do m=1,LIS_rc%nensem(n)-1
                 t = (i-1)*LIS_rc%nensem(n)+m
                 if(m.ne.LIS_rc%nensem(n)) then
                    delta(5) = delta(5) + &
                         (NOAHMP36_struc(n)%noahmp36(t)%wa - &
                         NOAHMP36_struc(n)%noahmp36(t_unpert)%wa)
                 endif
              enddo

              do j=1,4
                 delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                    t = (i-1)*LIS_rc%nensem(n)+m

                    SOILTYP = NOAHMP36_struc(n)%noahmp36(t)%soiltype
                    MAX_THRESHOLD = MAXSMC (SOILTYP)
                    !MIN_THRESHOLD = WLTSMC (SOILTYP)
                    if(ens_flag(m)) then
                       tmpval = NOAHMP36_struc(n)%noahmp36(t)%smc(j) - &
                            delta(j)

                       if(.not.(tmpval.le.MIN_THRESHOLD .or.&
                            tmpval.gt.(MAX_THRESHOLD))) then

                          NOAHMP36_struc(n)%noahmp36(t)%smc(j) = &
                               NOAHMP36_struc(n)%noahmp36(t)%smc(j) - delta(j)
                          NOAHMP36_struc(n)%noahmp36(t)%sh2o(j) = &
                               NOAHMP36_struc(n)%noahmp36(t)%sh2o(j) - delta(j)
                          bounds_violation = .false.
                       endif
                    endif

                    tmpval = NOAHMP36_struc(n)%noahmp36(t)%smc(j)

                    if(tmpval.le.MIN_THRESHOLD .or.&
                         tmpval.gt.(MAX_THRESHOLD)) then
                       bounds_violation = .true.
                    else
                       bounds_violation = .false.
                    endif
                 enddo
              enddo

              delta(5) =delta(5)/(LIS_rc%nensem(n)-1)
              do m=1,LIS_rc%nensem(n)-1
                 t = (i-1)*LIS_rc%nensem(n)+m

                 if(ens_flag(m)) then
                    tmpval = NOAHMP36_struc(n)%noahmp36(t)%wa - &
                         delta(5)

                    if(.not.(tmpval.le.MIN_GWS_THRESHOLD .or.&
                         tmpval.gt.MAX_GWS_THRESHOLD)) then

                       NOAHMP36_struc(n)%noahmp36(t)%wa = &
                            NOAHMP36_struc(n)%noahmp36(t)%wa - delta(5)
                       bounds_violation = .false.
                    endif
                 endif

                 tmpval = NOAHMP36_struc(n)%noahmp36(t)%wa

                 if(tmpval.le.MIN_GWS_THRESHOLD .or.&
                      tmpval.gt.MAX_GWS_THRESHOLD) then
                    bounds_violation = .true.
                 else
                    bounds_violation = .false.
                 endif
              enddo

              if(nIter.gt.10.and.bounds_violation) then
                 !-------------------------------------------------------------
                 ! All else fails, set to the bounds
                 !-------------------------------------------------------------

                 write(LIS_logunit,*) &
                      '[ERR] Ensemble structure violates physical bounds '
                 write(LIS_logunit,*) &
                      '[ERR] Please adjust the perturbation settings ..'

                 do j=1,4
                    do m=1,LIS_rc%nensem(n)
                       t = (i-1)*LIS_rc%nensem(n)+m

                       SOILTYP = NOAHMP36_struc(n)%noahmp36(t)%soiltype
                       MAX_THRESHOLD = MAXSMC (SOILTYP)
                       !MIN_THRESHOLD = WLTSMC (SOILTYP)

                       if(NOAHMP36_struc(n)%noahmp36(t)%sh2o(j).gt. &
                            MAX_THRESHOLD.or.&
                            NOAHMP36_struc(n)%noahmp36(t)%smc(j).gt. &
                            MAX_THRESHOLD) then
                          NOAHMP36_struc(n)%noahmp36(t)%sh2o(j) = MAX_THRESHOLD
                          NOAHMP36_struc(n)%noahmp36(t)%smc(j) = MAX_THRESHOLD
                       endif

                       if(NOAHMP36_struc(n)%noahmp36(t)%sh2o(j).lt. &
                            MIN_THRESHOLD.or.&
                            NOAHMP36_struc(n)%noahmp36(t)%smc(j).lt. &
                            MIN_THRESHOLD) then
                          NOAHMP36_struc(n)%noahmp36(t)%sh2o(j) = MIN_THRESHOLD
                          NOAHMP36_struc(n)%noahmp36(t)%smc(j) = MIN_THRESHOLD
                       endif
                    enddo
                    call LIS_endrun()
                 enddo

#if 0
                 do m=1,LIS_rc%nensem(n)
                    t = (i-1)*LIS_rc%nensem(n)+m

                    if(NOAHMP36_struc(n)%noahmp36(t)%wa.gt. &
                         MAX_GWS_THRESHOLD) then
                       NOAHMP36_struc(n)%noahmp36(t)%wa = MAX_GWS_THRESHOLD
                    endif

                    if(NOAHMP36_struc(n)%noahmp36(t)%wa.lt. &
                         MIN_GWS_THRESHOLD) then
                       NOAHMP36_struc(n)%noahmp36(t)%wa = MIN_GWS_THRESHOLD
                    endif
                 enddo
#endif
                 do m=1,LIS_rc%nensem(n)-1
                    t = (i-1)*LIS_rc%nensem(n)+m
                    t_unpert = i*LIS_rc%nensem(n)

                    if(.not.(NOAHMP36_struc(n)%noahmp36(t_unpert)%wa.le. &
                         MIN_GWS_THRESHOLD .or.&
                         NOAHMP36_struc(n)%noahmp36(t_unpert)%wa.gt. &
                         (MAX_GWS_THRESHOLD))) then

                       NOAHMP36_struc(n)%noahmp36(t)%wa = &
                            NOAHMP36_struc(n)%noahmp36(t_unpert)%wa

                    endif
                 enddo

                 !               call LIS_endrun()
              endif

           end do
        endif

     endif
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     swe_old = noahmp36_struc(n)%noahmp36(t)%sneqv
     snowh_old = noahmp36_struc(n)%noahmp36(t)%snowh

     dsneqv = swe(t)*1000.0 - noahmp36_struc(n)%noahmp36(t)%sneqv !in mm
     dsnowh = snod(t) - noahmp36_struc(n)%noahmp36(t)%snowh  !in m

!alternate option
#if 0
!only update dshowh, update SWE based on density
     dsnowh = snod(t) - noahmp36_struc(n)%noahmp36(t)%snowh  !in m
     snow_dens = noahmp36_struc(n)%noahmp36(t)%sneqv/&
          noahmp36_struc(n)%noahmp36(t)%snowh
     swe_new = snow_dens*snod(t)

     dsneqv = swe_new - noahmp36_struc(n)%noahmp36(t)%sneqv !in mm

#endif
     call noahmp36_snow_update(n, t, dsneqv, dsnowh)

     if(noahmp36_struc(n)%noahmp36(t)%sneqv.lt.0.or.&
          noahmp36_struc(n)%noahmp36(t)%snowh.lt.0) then
        write(LIS_logunit,*)'[ERR] State update leads to negative snow!'
        write(LIS_logunit,*)'[ERR] dsneqv, dsnowh = ', dsneqv, dsnowh
        write(LIS_logunit,*)'[ERR] swe(t), snod(t) = ', swe(t), snod(t)
        write(LIS_logunit,*)'[ERR] Stopping...'
        call LIS_endrun()
        !print*, dsneqv, dsnowh
        !print*, swe(t), snod(t)
        !stop
     endif

     if(noahmp36_struc(n)%noahmp36(t)%sneqv.gt.0.and.&
          noahmp36_struc(n)%noahmp36(t)%snowh.le.0) then
        write(LIS_logunit,*) &
             '[WARN] state update leads to positive SWE but zero snow depth'
        write(LIS_logunit,*) &
             '[WARN] ',noahmp36_struc(n)%noahmp36(t)%sneqv,&
             noahmp36_struc(n)%noahmp36(t)%snowh
        write(LIS_logunit,*) '[WARN] rejecting the update ...'
        noahmp36_struc(n)%noahmp36(t)%sneqv = swe_old
        noahmp36_struc(n)%noahmp36(t)%snowh = snowh_old
     endif

     if(noahmp36_struc(n)%noahmp36(t)%sneqv.le.0.and.&
          noahmp36_struc(n)%noahmp36(t)%snowh.gt.0) then
        write(LIS_logunit,*) &
             '[WARN] state update leads to positive snow depth but zero SWE'
        write(LIS_logunit,*) '[WARN] ',noahmp36_struc(n)%noahmp36(t)%sneqv,&
             noahmp36_struc(n)%noahmp36(t)%snowh
        write(LIS_logunit,*) '[WARN] rejecting the update ...'

        noahmp36_struc(n)%noahmp36(t)%sneqv = swe_old
        noahmp36_struc(n)%noahmp36(t)%snowh = snowh_old
     endif

  enddo

end subroutine noahmp36_settws

