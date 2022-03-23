!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_settws
!  \label{noahmp401_settws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li; created for Noah-MP4.0.1 
!
! !INTERFACE:
subroutine noahmp401_settws(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahMP401_lsmMod
  use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE

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
!  real, parameter        :: MIN_THRESHOLD = 0.02 
  real                   :: MIN_THRESHOLD 
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

  real                   :: delta
  logical                :: diffCheck(LIS_rc%ngrid(n))
  logical                :: ensCheck(LIS_rc%ngrid(n))
  logical                :: largeSM(LIS_rc%ngrid(n))
  integer                :: i, c,r,t,m
  integer                :: SOILTYP           ! soil type index [-]
  real                   :: sh2o_tmp, sh2o_rnd 
  real                   :: dsneqv,dsnowh,swe_old, snowh_old
  integer                :: status
  logical                :: rc1,rc2,rc3,rc4,rc5
  
  
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in noahmp401_settws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 2 failed in noahmp401_settws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 3 failed in noahmp401_settws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 4 failed in noahmp401_settws")
  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp401_settws")
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: SWE failed in noahmp401_settws")
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Snowdepth failed in noahmp401_settws")


  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmp401_settws")
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in noahmp401_settws")
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in noahmp401_settws")
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in noahmp401_settws")
  call ESMF_FieldGet(gwField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Groundwater Storage failed in noahmp401_settws")
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: SWE failed in noahmp401_settws")
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Snowdepth failed in noahmp401_settws")


  ensCheck = .true.
  diffCheck = .false.
  largeSM  = .false.
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     
     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     i = LIS_domain(n)%gindex(c,r)

     SOILTYP = NOAHMP401_struc(n)%noahmp401(t)%soiltype        
     MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP) 

     !locations with large soil moisture values are ice points.
     !we turn off the increments in such locations.
     if(noahmp401_struc(n)%noahmp401(t)%smc(1).gt.MAX_THRESHOLD.or.&
         noahmp401_struc(n)%noahmp401(t)%smc(1).gt.0.50) then 
        largeSM(i) = .true.
     endif
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     
     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     i = LIS_domain(n)%gindex(c,r)
     if(largeSM(i)) then 
        soilm1(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
        soilm2(t) = noahmp401_struc(n)%noahmp401(t)%smc(2)
        soilm3(t) = noahmp401_struc(n)%noahmp401(t)%smc(3)
        soilm4(t) = noahmp401_struc(n)%noahmp401(t)%smc(4)
        gws(t) = NOAHMP401_struc(n)%noahmp401(t)%wa
        snod(t) = noahmp401_struc(n)%noahmp401(t)%snowh 
        swe(t) = noahmp401_struc(n)%noahmp401(t)%sneqv
     endif
  enddo
  
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     i = LIS_domain(n)%gindex(c,r)
     
     SOILTYP = NOAHMP401_struc(n)%noahmp401(t)%soiltype        
     MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP) 
     MIN_THRESHOLD = 0.02 !SMCWLT_TABLE(SOILTYP)
     sm_threshold = MAX_THRESHOLD - 0.02
     
     if((soilm1(t).lt.MIN_THRESHOLD.or.&
          soilm1(t).gt.MAX_THRESHOLD).or.&
          (soilm2(t).lt.MIN_THRESHOLD.or.&
          soilm2(t).gt.MAX_THRESHOLD).or.&
          (soilm3(t).lt.MIN_THRESHOLD.or.&
          soilm3(t).gt.MAX_THRESHOLD).or.&
          (soilm4(t).lt.MIN_THRESHOLD.or.&
          soilm4(t).gt.MAX_THRESHOLD).or.&
          (gws(t).lt.MIN_GWS_THRESHOLD.or.&
          gws(t).gt.MAX_GWS_THRESHOLD)) then
        ensCheck(i) = .false.
     endif
     if((soilm1(t).ne.soilm1(i*LIS_rc%nensem(n))).and.&
          (soilm2(t).ne.soilm2(i*LIS_rc%nensem(n))).and.&
          (soilm3(t).ne.soilm3(i*LIS_rc%nensem(n))).and.&
          (soilm4(t).ne.soilm4(i*LIS_rc%nensem(n))).and.&
          (gws(t).ne.gws(i*LIS_rc%nensem(n)))) then 
        diffCheck(i) = .true.
     endif
  enddo

  do i=1,LIS_rc%ngrid(n)
     rc1 = .true.
     rc2 = .true.
     rc3 = .true.
     rc4 = .true.
     rc5 = .true. 
     if(.not.ensCheck(i).and.diffCheck(i).and.(.not.largeSM(i))) then
        call noahmp401_tws_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             soilm1((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_THRESHOLD, MAX_THRESHOLD,rc1)
        call noahmp401_tws_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             soilm2((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_THRESHOLD, MAX_THRESHOLD,rc2)
        call noahmp401_tws_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             soilm3((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_THRESHOLD, MAX_THRESHOLD,rc3)
        call noahmp401_tws_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             soilm4((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_THRESHOLD, MAX_THRESHOLD,rc4)
        call noahmp401_tws_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             gws((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_GWS_THRESHOLD, MAX_GWS_THRESHOLD,rc5)
     endif
     if(.not.rc1.or.&
          .not.rc2.or.&
          .not.rc3.or.&
          .not.rc4.or.&
          .not.rc5) then

        do m=1,LIS_rc%nensem(n)
           t = (i-1)*LIS_rc%nensem(n)+m
           soilm1(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
           soilm2(t) = noahmp401_struc(n)%noahmp401(t)%smc(2)
           soilm3(t) = noahmp401_struc(n)%noahmp401(t)%smc(3)
           soilm4(t) = noahmp401_struc(n)%noahmp401(t)%smc(4)
           gws(t) = NOAHMP401_struc(n)%noahmp401(t)%wa
           snod(t) = noahmp401_struc(n)%noahmp401(t)%snowh 
           swe(t) = noahmp401_struc(n)%noahmp401(t)%sneqv
        enddo
     endif
  enddo
        
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     delta = soilm1(t) - NOAHMP401_struc(n)%noahmp401(t)%smc(1)    
     NOAHMP401_struc(n)%noahmp401(t)%smc(1) = soilm1(t)
     NOAHMP401_struc(n)%noahmp401(t)%sh2o(1) = &
          NOAHMP401_struc(n)%noahmp401(t)%sh2o(1) + delta

     delta = soilm2(t) - NOAHMP401_struc(n)%noahmp401(t)%smc(2)    
     NOAHMP401_struc(n)%noahmp401(t)%smc(2) = soilm2(t)
     NOAHMP401_struc(n)%noahmp401(t)%sh2o(2) = &
          NOAHMP401_struc(n)%noahmp401(t)%sh2o(2) + delta

     delta = soilm3(t) - NOAHMP401_struc(n)%noahmp401(t)%smc(3)    
     NOAHMP401_struc(n)%noahmp401(t)%smc(3) = soilm3(t)
     NOAHMP401_struc(n)%noahmp401(t)%sh2o(3) = &
          NOAHMP401_struc(n)%noahmp401(t)%sh2o(3) + delta

     delta = soilm4(t) - NOAHMP401_struc(n)%noahmp401(t)%smc(4)    
     NOAHMP401_struc(n)%noahmp401(t)%smc(4) = soilm4(t)
     NOAHMP401_struc(n)%noahmp401(t)%sh2o(4) = &
          NOAHMP401_struc(n)%noahmp401(t)%sh2o(4) + delta
     
     NOAHMP401_struc(n)%noahmp401(t)%wa=gws(t)
  enddo
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     swe_old = noahmp401_struc(n)%noahmp401(t)%sneqv
     snowh_old = noahmp401_struc(n)%noahmp401(t)%snowh

     dsneqv = swe(t) - noahmp401_struc(n)%noahmp401(t)%sneqv !in mm
     dsnowh = snod(t) - noahmp401_struc(n)%noahmp401(t)%snowh  !in m

     !alternate option
     call noahmp401_snow_update(n, t, dsneqv, dsnowh)

     if(noahmp401_struc(n)%noahmp401(t)%sneqv.eq.0.or.&
          noahmp401_struc(n)%noahmp401(t)%snowh.eq.0) then
        noahmp401_struc(n)%noahmp401(t)%sneqv = 0
        noahmp401_struc(n)%noahmp401(t)%snowh = 0
     endif
  enddo

end subroutine noahmp401_settws


subroutine noahmp401_tws_reorderEnsForOutliers(i,nensem, statevec, &
     minvalue,maxvalue, status)
  
  use LIS_coreMod
  
  implicit none
  integer              :: i
  integer              :: nensem
  real                 :: statevec(nensem)
  real                 :: minvalue,maxvalue
  logical              :: status
  
  real                 :: minvT, maxvT, minvG, maxvG
  integer              :: k
  real                 :: spread_total, spread_good, spread_ratio
  
  !Ensemble spread (total and with 'good' ensemble members
  minvT = 1E10
  maxvT = -1E10
  minvG = 1E10
  maxvG = -1E10
  status = .true. 
  
  do k=1,nensem

     if(statevec(k).lt.minvT) then 
        minvT = statevec(k)
     endif
     if(statevec(k).gt.maxvT) then 
        maxvT = statevec(k)
     endif

     if(statevec(k).gt.minvalue.and.statevec(k).lt.maxvalue) then 
        if(statevec(k).lt.minvG) then 
           minvG = statevec(k)
        endif
        if(statevec(k).gt.maxvG) then 
           maxvG = statevec(k)
        endif
     endif
  enddo
  
  if(minvG.eq.1E10.and.maxvG.eq.-1E10) then
     !all members are unphysical.

     statevec = minvalue
     status = .false.
     
  else
     spread_total = (maxvT - minvT)
     spread_good  = (maxvG - minvG)
     
     spread_ratio = spread_good/spread_total
     
     !rescale the ensemble 
     
     do k=1,nensem-1
        statevec(k) = statevec(nensem) + &
             (statevec(k) - statevec(nensem))*spread_ratio 
     enddo
  endif

end subroutine noahmp401_tws_reorderEnsForOutliers
