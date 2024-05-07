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
! !ROUTINE: jules50_setsoilm
!  \label{jules50_setsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 23Apr2018: Mahdi Navari: Modified for JULES 5.0 
! 25Oct2019: Yonghwan Kwon: Modified to fix a bug for glacier grid
! 08Apr2020: Yonghwan Kwon: Modified to consider frozen soil moisture (p_s_sthf)
!
! There are 3 cases 
! 1- If all the ensemble members met the update conditions --> apply the update
! 2- If more than 50% of the ensemble members met the update condition --> 
!    apply the update for that members and set the other member to the mean 
!    value of the ensemble (i.e. mean of the members that met the conditions)
! 3- If less then 50% of the ensemble members met the update conditions --> 
!    adjust the states 
! 8 May 2023: Mahdi Navari; Soil temperature bias bug fix
!                           (add check for frozen soil)
 
! !INTERFACE:
subroutine jules50_setsoilm(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use jules50_lsmMod
 !MN: added to use LIS_sfmodel_struc
  !use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
  use jules_soil_mod, only:  dzsoil !
  use jules_surface_mod,      only: l_aggregate   
  use LIS_constantsMod, only : LIS_CONST_TKFRZ  

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to JULES's
!  model space. 
! 
!EOP
  real                   :: MIN_THRESHOLD(jules50_struc(n)%sm_levels)   !Yonghwan Kwon
  real                   :: MAX_THRESHOLD(jules50_struc(n)%sm_levels)
  real                   :: sm_threshold(jules50_struc(n)%sm_levels)
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  integer                :: t, j,i, gid, m, t_unpert, row, col
  integer                :: status
  real                   :: delta(4)
  real                   :: delta1,delta2,delta3,delta4 , lat, lon
  real                   :: tmpval
  real                   :: tmpval_u, tmpval_f  !Yonghwan Kwon
  real                   :: timenow ! for print 
  logical                :: bounds_violation
  integer                :: nIter, N_ens
  logical                :: update_flag(LIS_rc%ngrid(n))
  logical                :: ens_flag(LIS_rc%nensem(n))
! mn
  real                   :: tmp(LIS_rc%nensem(n)), tmp0(LIS_rc%nensem(n))
  real                   :: tmp1(LIS_rc%nensem(n)),tmp2(LIS_rc%nensem(n)),tmp3(LIS_rc%nensem(n)),tmp4(LIS_rc%nensem(n)),tmp5(LIS_rc%nensem(n)) ,tmp6(LIS_rc%nensem(n)) ,tmp7(LIS_rc%nensem(n)) ,tmp8(LIS_rc%nensem(n))  
  logical                :: update_flag_tile(LIS_rc%npatch(n,LIS_rc%lsm_index))
  logical                :: flag_ens(LIS_rc%ngrid(n))
  logical                :: flag_tmp(LIS_rc%nensem(n))
  logical                :: update_flag_ens(LIS_rc%ngrid(n))
  logical                :: update_flag_new(LIS_rc%ngrid(n))
  integer                :: pcount, icount
  real                   :: MaxEnsSM1 ,MaxEnsSM2 ,MaxEnsSM3 ,MaxEnsSM4
  real                   :: MinEnsSM1 ,MinEnsSM2 ,MinEnsSM3 ,MinEnsSM4 
!  real                   :: MaxEns_p_s_sthu1, MaxEns_p_s_sthu2, MaxEns_p_s_sthu3, MaxEns_p_s_sthu4
  real                   :: smc_rnd, smc_tmp 
  real                   :: p_s_sthu_tmp, p_s_sthu_rnd 
  INTEGER, DIMENSION (1) :: seed
  !-----------------------------------------------------Yonghwan Kwon
  !real, dimension(4)     :: frac_sthu_l, frac_sthf_l
  real, dimension(jules50_struc(n)%sm_levels)  :: sat_p, p_s_sth, frac_sthu, frac_sthf
  real, dimension(4)     :: delta_u, delta_f
  real                   :: delta1_sthu, delta2_sthu, delta3_sthu, delta4_sthu
  real                   :: delta1_sthf, delta2_sthf, delta3_sthf, delta4_sthf
  integer                :: ilat, ilon
  logical                :: violation_new 
  !-----------------------------------------------------

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in jules50_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 2 failed in jules50_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 3 failed in jules50_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 4 failed in jules50_setsoilm")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in jules50_setsoilm")
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in jules50_setsoilm")
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in jules50_setsoilm")
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in jules50_setsoilm")

  update_flag = .true. 
  update_flag_tile = .true. 

! MN: NOTE: be careful with the unit. The unit of soil moisture diagnostic variable (smcl) [kg/m2]
! differs from that of the prognostic variable (p_s_sthu, p_s_smvcst) [m3/m3].

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     !-------------------------------------Yonghwan Kwon
     do j=1,jules50_struc(n)%sm_levels  
        MAX_THRESHOLD(j) = jules50_struc(n)%jules50(t)%p_s_smvcst(j) ! Volumetric saturation point (m^3 m-3 of soil)
        MIN_THRESHOLD(j) = jules50_struc(n)%jules50(t)%p_s_smvcwt(j) ! Volumetric wilting point (m^3 m-3 of soil) !Yonghwan Kwon 
        !sm_threshold(j)  = MAX_THRESHOLD(j) - MIN_THRESHOLD(j)      
        sm_threshold(j)  = MAX_THRESHOLD(j)

        sat_p(j) = jules50_struc(n)%jules50(t)%p_s_smvcst(j) ! Volumetric saturation point (m^3 m-3 of soil)
        p_s_sth(j) = jules50_struc(n)%jules50(t)%p_s_sthu(j) + jules50_struc(n)%jules50(t)%p_s_sthf(j)  !saturated fraction
     enddo     
     !-------------------------------------
     gid = LIS_domain(n)%gindex(&
           LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
           LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     !MN: delta = X(+) - X(-)
     !NOTE: "jules50_updatesoilm.F90" updates the soilmx(t)   
     delta1 = (soilm1(t)-jules50_struc(n)%jules50(t)%smcl_soilt(1))* 1/dzsoil(1)*1/1000  ! [kg/m2]*1/m*1/(kg/m3) --> [m3/m3]
     delta2 = (soilm2(t)-jules50_struc(n)%jules50(t)%smcl_soilt(2))* 1/dzsoil(2)*1/1000 
     delta3 = (soilm3(t)-jules50_struc(n)%jules50(t)%smcl_soilt(3))* 1/dzsoil(3)*1/1000
     delta4 = (soilm4(t)-jules50_struc(n)%jules50(t)%smcl_soilt(4))* 1/dzsoil(4)*1/1000
     
     ! MN: check MIN_THRESHOLD < volumetric liquid soil moisture < threshold 
     ! unit conversion -->   ..%p_s_sthu(1)  *  ..%p_s_smvcst(1)-->   [-] * [m3/m3] 

     if (jules50_struc(n)%jules50(t)%p_s_smvcst(1) == 0) then     !Yonghwan Kwon: set update_flag to false for glacier land 
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     else      
        ! frozen soil moisture content has been added (Yonghwan Kwon)
        if (p_s_sth(1) * sat_p(1) + delta1.gt.MIN_THRESHOLD(1) .and.&
            p_s_sth(1) * sat_p(1) + delta1.lt.sm_threshold(1) .and.&
            jules50_struc(n)%jules50(t)%t_soil(1) .gt. LIS_CONST_TKFRZ) then  
           update_flag(gid) = update_flag(gid).and.(.true.)
           ! MN save the flag for each tile (col*row*ens) 
           update_flag_tile(t) = update_flag_tile(t).and.(.true.)
        else
           update_flag(gid) = update_flag(gid).and.(.false.)
           update_flag_tile(t) = update_flag_tile(t).and.(.false.)
        endif

        if (p_s_sth(2) * sat_p(2) + delta2.gt.MIN_THRESHOLD(2) .and.&
            p_s_sth(2) * sat_p(2) + delta2.lt.sm_threshold(2) .and.&
            jules50_struc(n)%jules50(t)%t_soil(2) .gt. LIS_CONST_TKFRZ) then
           update_flag(gid) = update_flag(gid).and.(.true.)
           update_flag_tile(t) = update_flag_tile(t).and.(.true.)
        else
           update_flag(gid) = update_flag(gid).and.(.false.)
           update_flag_tile(t) = update_flag_tile(t).and.(.false.)
        endif

        if (p_s_sth(3) * sat_p(3) + delta3.gt.MIN_THRESHOLD(3) .and.&
            p_s_sth(3) * sat_p(3) + delta3.lt.sm_threshold(3) .and.&
            jules50_struc(n)%jules50(t)%t_soil(3) .gt. LIS_CONST_TKFRZ) then
           update_flag(gid) = update_flag(gid).and.(.true.)
           update_flag_tile(t) = update_flag_tile(t).and.(.true.)
        else
           update_flag(gid) = update_flag(gid).and.(.false.)
           update_flag_tile(t) = update_flag_tile(t).and.(.false.)
        endif

        if (p_s_sth(4) * sat_p(4) + delta4.gt.MIN_THRESHOLD(4) .and.&
            p_s_sth(4) * sat_p(4) + delta4.lt.sm_threshold(4) .and.&
            jules50_struc(n)%jules50(t)%t_soil(4) .gt. LIS_CONST_TKFRZ) then
           update_flag(gid) = update_flag(gid).and.(.true.)
           update_flag_tile(t) = update_flag_tile(t).and.(.true.)
        else
           update_flag(gid) = update_flag(gid).and.(.false.)
           update_flag_tile(t) = update_flag_tile(t).and.(.false.)
        endif
     endif
  enddo

!-----------------------------------------------------------------------------------------
! MN create new flag: if update flag for 50% of the ensemble members is true 
! then update the stats 
!-----------------------------------------------------------------------------------------
  update_flag_ens = .true.
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     gid = LIS_domain(n)%gindex(&
           LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
           LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row) 
     flag_tmp=update_flag_tile(i:i+LIS_rc%nensem(n)-1)
     !flag_tmp=update_flag_tile((i-1)*LIS_rc%nensem(n)+1:(i)*LIS_rc%nensem(n))
     pcount = COUNT(flag_tmp) ! Counts the number of .TRUE. elements
     if (pcount.lt.LIS_rc%nensem(n)*0.5) then   ! 50%
        update_flag_ens(gid)= .False.
     endif
     update_flag_new(gid)= update_flag(gid).or.update_flag_ens(gid)  ! new flag
  enddo
   
  ! update step
  ! loop over grid points 
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)     

     gid =LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)  
     
     !if(update_flag(gid)) then
     if(update_flag_new(gid)) then 
!-----------------------------------------------------------------------------------------
! Update the states
! case 1-1- if the update flag for a given grid and tile are TRUE --> apply the DA update    
! case 1-2- if the update flag for a given grid is true but for the tile is FALSE --> 
! set the update value for that ensemble memeber to the ens. mean  
! case 2- if the update flag for a given grid is FALSE readjust the states
!-----------------------------------------------------------------------------------------
! store update value for cases that flag_tile & update_flag_new are TRUE
! flag_tile = TRUE --> means both the min and max threshold has been satisfied

! compute the ensemble min and max for case 1-2        
        tmp1 = LIS_rc%udef
        tmp2 = LIS_rc%udef
        tmp3 = LIS_rc%udef
        tmp4 = LIS_rc%udef

        do m=1,LIS_rc%nensem(n)
           t = i+m-1
           !t = (i-1)*LIS_rc%nensem(n)+m
           
           if(update_flag_tile(t)) then
              
              tmp1(m) = soilm1(t) 
              tmp2(m) = soilm2(t) 
              tmp3(m) = soilm3(t) 
              tmp4(m) = soilm4(t) 
           endif
        enddo

        MaxEnsSM1 = -10000
        MaxEnsSM2 = -10000
        MaxEnsSM3 = -10000
        MaxEnsSM4 = -10000

        MinEnsSM1 = 10000
        MinEnsSM2 = 10000
        MinEnsSM3 = 10000
        MinEnsSM4 = 10000

        do m=1,LIS_rc%nensem(n)
           if(tmp1(m).ne.LIS_rc%udef) then 
              MaxEnsSM1 = max(MaxEnsSM1, tmp1(m)) ![kg/m2]
              MaxEnsSM2 = max(MaxEnsSM2, tmp2(m))
              MaxEnsSM3 = max(MaxEnsSM3, tmp3(m))
              MaxEnsSM4 = max(MaxEnsSM4, tmp4(m))

              MinEnsSM1 = min(MinEnsSM1, tmp1(m))
              MinEnsSM2 = min(MinEnsSM2, tmp2(m))
              MinEnsSM3 = min(MinEnsSM3, tmp3(m))
              MinEnsSM4 = min(MinEnsSM4, tmp4(m))
              
           endif
        enddo

        ! loop over tile       
        do m=1,LIS_rc%nensem(n)
            t = i+m-1
           !t = (i-1)*LIS_rc%nensem(n)+m

           !-------------------------------------Yonghwan Kwon
           do j=1,jules50_struc(n)%sm_levels
              MAX_THRESHOLD(j) = jules50_struc(n)%jules50(t)%p_s_smvcst(j) ! Volumetric saturation point (m^3 m-3 of soil)
              MIN_THRESHOLD(j) = jules50_struc(n)%jules50(t)%p_s_smvcwt(j) ! Volumetric wilting point (m^3 m-3 of soil) !Yonghwan Kwon 
              !sm_threshold(j)  = MAX_THRESHOLD(j) - MIN_THRESHOLD(j)
              sm_threshold(j)  = MAX_THRESHOLD(j)

              sat_p(j) = jules50_struc(n)%jules50(t)%p_s_smvcst(j) ! Volumetric saturation point (m^3 m-3 of soil)
              p_s_sth(j) = jules50_struc(n)%jules50(t)%p_s_sthu(j) + jules50_struc(n)%jules50(t)%p_s_sthf(j)  !saturated fraction
           enddo
           !-------------------------------------
                      
           ! MN check update status for each tile  
           if(update_flag_tile(t)) then

	      delta1 = (soilm1(t)-jules50_struc(n)%jules50(t)%smcl_soilt(1))* 1/dzsoil(1)*1/1000  ! [kg/m2]*1/m*1/(kg/m3) --> [m3/m3]
	      delta2 = (soilm2(t)-jules50_struc(n)%jules50(t)%smcl_soilt(2))* 1/dzsoil(2)*1/1000 
	      delta3 = (soilm3(t)-jules50_struc(n)%jules50(t)%smcl_soilt(3))* 1/dzsoil(3)*1/1000
	      delta4 = (soilm4(t)-jules50_struc(n)%jules50(t)%smcl_soilt(4))* 1/dzsoil(4)*1/1000   
 
              !------------------------------------- Yonghwan Kwon
              ! Compute delta for unfrozen and frozen soil moiture (Yonghwan Kwon)
              delta1_sthu = delta1 * (jules50_struc(n)%jules50(t)%p_s_sthu(1) / p_s_sth(1))  ![m3/m3]
              delta1_sthf = delta1 * (jules50_struc(n)%jules50(t)%p_s_sthf(1) / p_s_sth(1))  ![m3/m3]
              delta2_sthu = delta2 * (jules50_struc(n)%jules50(t)%p_s_sthu(2) / p_s_sth(2))  ![m3/m3]
              delta2_sthf = delta2 * (jules50_struc(n)%jules50(t)%p_s_sthf(2) / p_s_sth(2))  ![m3/m3]
              delta3_sthu = delta3 * (jules50_struc(n)%jules50(t)%p_s_sthu(3) / p_s_sth(3))  ![m3/m3]
              delta3_sthf = delta3 * (jules50_struc(n)%jules50(t)%p_s_sthf(3) / p_s_sth(3))  ![m3/m3]
              delta4_sthu = delta4 * (jules50_struc(n)%jules50(t)%p_s_sthu(4) / p_s_sth(4))  ![m3/m3]
              delta4_sthf = delta4 * (jules50_struc(n)%jules50(t)%p_s_sthf(4) / p_s_sth(4))  ![m3/m3]

              jules50_struc(n)%jules50(t)%p_s_sthu(1) = (jules50_struc(n)%jules50(t)%p_s_sthu(1) * sat_p(1)&
                                                        + delta1_sthu)/sat_p(1) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction
              jules50_struc(n)%jules50(t)%p_s_sthf(1) = (jules50_struc(n)%jules50(t)%p_s_sthf(1) * sat_p(1)&
                                                        + delta1_sthf)/sat_p(1) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction
              jules50_struc(n)%jules50(t)%smcl_soilt(1) = soilm1(t) ![kg/m2]

              if(soilm1(t).lt.0) then 
                 write(LIS_logunit, *) 'setsoilm1 ',t,soilm1(t)
                 call LIS_endrun
              endif
	      ! I think the following sentence is redundant. 
	      !Because both “update_flag_new” and “update_flag_tile” are TRUE 

              !if(jules50_struc(n)%jules50(t)%p_s_sthu(2)+delta2.gt.MIN_THRESHOLD .and.&   
              !     jules50_struc(n)%jules50(t)%p_s_sthu(2)+delta2.lt.sm_threshold) then 
              !   jules50_struc(n)%jules50(t)%p_s_sthu(2) = jules50_struc(n)%jules50(t)%p_s_sthu(2)+&
              !        soilm2(t)-jules50_struc(n)%jules50(t)%smcl_soilt(2)
              jules50_struc(n)%jules50(t)%p_s_sthu(2) = (jules50_struc(n)%jules50(t)%p_s_sthu(2) * sat_p(2)&
                                                        + delta2_sthu)/sat_p(2) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction
              jules50_struc(n)%jules50(t)%p_s_sthf(2) = (jules50_struc(n)%jules50(t)%p_s_sthf(2) * sat_p(2)&
                                                        + delta2_sthf)/sat_p(2) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction
              jules50_struc(n)%jules50(t)%smcl_soilt(2) = soilm2(t)

              if(soilm2(t).lt.0) then 
                 write(LIS_logunit, *) 'setsoilm2 ',t,soilm2(t)
                 call LIS_endrun
              endif

              !endif
              !if(jules50_struc(n)%jules50(t)%p_s_sthu(3)+delta3.gt.MIN_THRESHOLD .and.&
              !     jules50_struc(n)%jules50(t)%p_s_sthu(3)+delta3.lt.sm_threshold) then 
              !   jules50_struc(n)%jules50(t)%p_s_sthu(3) = jules50_struc(n)%jules50(t)%p_s_sthu(3)+&
              !        soilm3(t)-jules50_struc(n)%jules50(t)%smcl_soilt(3)
              jules50_struc(n)%jules50(t)%p_s_sthu(3) = (jules50_struc(n)%jules50(t)%p_s_sthu(3) * sat_p(3)&
                                                        + delta3_sthu)/sat_p(3) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction
              jules50_struc(n)%jules50(t)%p_s_sthf(3) = (jules50_struc(n)%jules50(t)%p_s_sthf(3) * sat_p(3)&
                                                        + delta3_sthf)/sat_p(3) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction
              jules50_struc(n)%jules50(t)%smcl_soilt(3) = soilm3(t)

              if(soilm3(t).lt.0) then 
                 write(LIS_logunit, *) 'setsoilm3 ',t,soilm3(t)
                 call LIS_endrun
              endif

              !endif
              ! surface layer
              !if(jules50_struc(n)%jules50(t)%p_s_sthu(4)+delta4.gt.MIN_THRESHOLD .and.&
              !     jules50_struc(n)%jules50(t)%p_s_sthu(4)+delta4.lt.sm_threshold) then 
              !   jules50_struc(n)%jules50(t)%p_s_sthu(4) = jules50_struc(n)%jules50(t)%p_s_sthu(4)+&
              !        soilm4(t)-jules50_struc(n)%jules50(t)%smcl_soilt(4)
              jules50_struc(n)%jules50(t)%p_s_sthu(4) = (jules50_struc(n)%jules50(t)%p_s_sthu(4) * sat_p(4)&
                                                        + delta4_sthu)/sat_p(4) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction
              jules50_struc(n)%jules50(t)%p_s_sthf(4) = (jules50_struc(n)%jules50(t)%p_s_sthf(4) * sat_p(4)&
                                                        + delta4_sthf)/sat_p(4) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction
              jules50_struc(n)%jules50(t)%smcl_soilt(4) = soilm4(t)

              if(soilm4(t).lt.0) then 
                 write(LIS_logunit, *) 'setsoilm4 ',t,soilm4(t)
                 call LIS_endrun
              endif
              !endif
           else 
 
!-----------------------------------------------------------------------------------------  
! USE ENSEMBLE MEAN VALUE 
!-----------------------------------------------------------------------------------------
              
              ! Assume p_s_sthu = smc (i.e. ice content=0) 

              !-----------------------------------------------------
              ! Frozen soil moisture has been added (Yonghwan Kwon)
              ! Compute fraction
              do j=1,jules50_struc(n)%sm_levels
                 frac_sthu(j) = jules50_struc(n)%jules50(t)%p_s_sthu(j) / p_s_sth(j)
                 frac_sthf(j) = jules50_struc(n)%jules50(t)%p_s_sthf(j) / p_s_sth(j)
              enddo
              !-----------------------------------------------------

              smc_tmp = (MaxEnsSM1 - MinEnsSM1)/2 + MinEnsSM1 ! [kg/m2]
              jules50_struc(n)%jules50(t)%p_s_sthu(1) = (smc_tmp * 1/dzsoil(1)*1/1000 * frac_sthu(1)) / sat_p(1) 
                                                                                ! ([kg/m2]*(1/m*1/(kg/m3)))/[m3/m3]--> fraction
              jules50_struc(n)%jules50(t)%p_s_sthf(1) = (smc_tmp * 1/dzsoil(1)*1/1000 * frac_sthf(1)) / sat_p(1) 
                                                                                ! ([kg/m2]*(1/m*1/(kg/m3)))/[m3/m3]--> fraction  !Yonghwan Kwon
              jules50_struc(n)%jules50(t)%smcl_soilt(1) = smc_tmp

              smc_tmp = (MaxEnsSM2 - MinEnsSM2)/2 + MinEnsSM2            
              jules50_struc(n)%jules50(t)%p_s_sthu(2) = (smc_tmp * 1/dzsoil(2)*1/1000 * frac_sthu(2)) / sat_p(2)
              jules50_struc(n)%jules50(t)%p_s_sthf(2) = (smc_tmp * 1/dzsoil(2)*1/1000 * frac_sthf(2)) / sat_p(2)
              jules50_struc(n)%jules50(t)%smcl_soilt(2) = smc_tmp

              smc_tmp = (MaxEnsSM3 - MinEnsSM3)/2 + MinEnsSM3
              jules50_struc(n)%jules50(t)%p_s_sthu(3) = (smc_tmp * 1/dzsoil(3)*1/1000 * frac_sthu(3)) / sat_p(3)
              jules50_struc(n)%jules50(t)%p_s_sthf(3) = (smc_tmp * 1/dzsoil(3)*1/1000 * frac_sthf(3)) / sat_p(3)
              jules50_struc(n)%jules50(t)%smcl_soilt(3) = smc_tmp

              smc_tmp = (MaxEnsSM4 - MinEnsSM4)/2 + MinEnsSM4
              jules50_struc(n)%jules50(t)%p_s_sthu(4) = (smc_tmp * 1/dzsoil(4)*1/1000 * frac_sthu(4)) / sat_p(4)
              jules50_struc(n)%jules50(t)%p_s_sthf(4) = (smc_tmp * 1/dzsoil(4)*1/1000 * frac_sthf(4)) / sat_p(4)
              jules50_struc(n)%jules50(t)%smcl_soilt(4) = smc_tmp
           endif ! flag for each tile           
        enddo ! loop over tile
        
     else ! if update_flag_new(gid) is FALSE   
        if(LIS_rc%pert_bias_corr.eq.1) then           

!--------------------------------------------------------------------------
! if no update is made, then we need to readjust the ensemble if pert bias
! correction is turned on because the forcing perturbations may cause 
! biases to persist. 
!--------------------------------------------------------------------------
           bounds_violation = .true. 
           violation_new = .false.
           nIter = 0
           ens_flag = .true. 
           
           do while(bounds_violation) 
              niter = niter + 1

              !t_unpert = i*LIS_rc%nensem(n)
	      t_unpert = i+LIS_rc%nensem(n)-1
              !do j=1,4
              do j=1,jules50_struc(n)%sm_levels  !Yonghwan Kwon
                 delta_u(j) = 0.0
                 delta_f(j) = 0.0
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
   
                    if(m.ne.LIS_rc%nensem(n)) then 
                    ! NOTE: be careful with the unit. do not change the unit of delta in the loop
                    ! [m3w/m3s]  + ([-] - [-])*[m3/m3] 
                       delta_u(j) = delta_u(j)+ & !  * 1/dzsoil(j)*1/1000 
                                  (jules50_struc(n)%jules50(t)%p_s_sthu(j) - jules50_struc(n)%jules50(t_unpert)%p_s_sthu(j)) *&
			          jules50_struc(n)%jules50(t)%p_s_smvcst(j)  !) &
			          ! / (1/dzsoil(j)*1/1000) 
                       delta_f(j) = delta_f(j)+ & !  
                                  (jules50_struc(n)%jules50(t)%p_s_sthf(j) - jules50_struc(n)%jules50(t_unpert)%p_s_sthf(j)) *&
                                  jules50_struc(n)%jules50(t)%p_s_smvcst(j) 
                    endif
                    
                 enddo
              enddo
              
              !do j=1,4
              do j=1,jules50_struc(n)%sm_levels  !Yonghwan Kwon
                 delta_u(j) = delta_u(j)/(LIS_rc%nensem(n)-1)
                 delta_f(j) = delta_f(j)/(LIS_rc%nensem(n)-1)
                 delta(j) = delta_u(j) + delta_f(j)
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m

		    MAX_THRESHOLD(j) = jules50_struc(n)%jules50(t)%p_s_smvcst(j) ! assume it is the same for all layers [m3/m3]
                    MIN_THRESHOLD(j) = jules50_struc(n)%jules50(t)%p_s_smvcwt(j) ! Volumetric wilting point (m^3 m-3 of soil) !Yonghwan Kwon 
                    !sm_threshold(j)  = MAX_THRESHOLD(j) - MIN_THRESHOLD(j)
                    sm_threshold(j)  = MAX_THRESHOLD(j)

                    sat_p(j) = jules50_struc(n)%jules50(t)%p_s_smvcst(j) ! Volumetric saturation point (m^3 m-3 of soil) !Yonghwan Kwon
                    p_s_sth(j) = jules50_struc(n)%jules50(t)%p_s_sthu(j) + jules50_struc(n)%jules50(t)%p_s_sthf(j)  !saturated fraction (Yonghwan Kwon)

                    !-----------------------------------------------------
                    ! Frozen soil moisture has been added (Yonghwan Kwon)
                    ! Compute fraction
     
                    if (p_s_sth(j) > 0) then
                       frac_sthu(j) = jules50_struc(n)%jules50(t)%p_s_sthu(j) / p_s_sth(j)
                       frac_sthf(j) = jules50_struc(n)%jules50(t)%p_s_sthf(j) / p_s_sth(j)
                    else
                       frac_sthu(j) = 0
                       frac_sthf(j) = 0
                    endif
                    !-----------------------------------------------------
                    
                    tmpval = p_s_sth(j) * sat_p(j) - delta(j) !* 1/dzsoil(j)*1/1000 ! [-][m3/m3]-[m3/m3] --> [m3/m3]                    
                    if(tmpval.le.MIN_THRESHOLD(j)) then 

                       if (sat_p(j) == 0) then
                          jules50_struc(n)%jules50(t)%p_s_sthu(j) = 0
                       else
                          jules50_struc(n)%jules50(t)%p_s_sthu(j) = &
                               (max(jules50_struc(n)%jules50(t_unpert)%p_s_sthu(j)*sat_p(j),&
                                    0.1*MIN_THRESHOLD(j)*frac_sthu(j))) / sat_p(j) ! max( [-][m3/m3] , [m3/m3] ) / [m3/m3] --> fraction
                          jules50_struc(n)%jules50(t)%p_s_sthf(j) = &
                               (max(jules50_struc(n)%jules50(t_unpert)%p_s_sthf(j)*sat_p(j),&
                                    0.1*MIN_THRESHOLD(j)*frac_sthf(j))) / sat_p(j) ! max( [-][m3/m3] , [m3/m3] ) / [m3/m3] --> fraction 
                       endif
                       jules50_struc(n)%jules50(t)%smcl_soilt(j) = &
                              (jules50_struc(n)%jules50(t)%p_s_sthu(j) + jules50_struc(n)%jules50(t)%p_s_sthf(j)) * sat_p(j)&
                              /(1/dzsoil(j)*1/1000)


                       !jules50_struc(n)%jules50(t)%smcl_soilt(j) = & 
                       !     (max(jules50_struc(n)%jules50(t_unpert)%smcl_soilt(j)&
                       !	        *1/dzsoil(j)*1/1000,&
                       !          0.1*MIN_THRESHOLD(j))) / (1/dzsoil(j)*1/1000) ! max( [kg/m2]*[1/m]*[1/kg/m3] , [m3/m3] ) / ([1/m][kg/m3]) --> [kg/m2]

                       ens_flag(m) = .false. 

                    elseif(tmpval.gt.sm_threshold(j)) then  ! MN ge --> gt

                       if (sat_p(j) == 0) then
                          jules50_struc(n)%jules50(t)%p_s_sthu(j) = 0
                       else
                          jules50_struc(n)%jules50(t)%p_s_sthu(j) = &
                               (min(jules50_struc(n)%jules50(t_unpert)%p_s_sthu(j)*sat_p(j),&
                                   sm_threshold(j)*frac_sthu(j))) / sat_p(j) ! min( [-][m3/m3] , [m3/m3] ) / [m3/m3] --> fraction 
                          jules50_struc(n)%jules50(t)%p_s_sthf(j) = &
                               (min(jules50_struc(n)%jules50(t_unpert)%p_s_sthf(j)*sat_p(j),&
                                   sm_threshold(j)*frac_sthf(j))) / sat_p(j) ! min( [-][m3/m3] , [m3/m3] ) / [m3/m3] --> fraction  
                       endif
                       jules50_struc(n)%jules50(t)%smcl_soilt(j) = &
                              (jules50_struc(n)%jules50(t)%p_s_sthu(j) + jules50_struc(n)%jules50(t)%p_s_sthf(j)) * sat_p(j)&
                              /(1/dzsoil(j)*1/1000)

                       !jules50_struc(n)%jules50(t)%smcl_soilt(j) = &
                       !     (min(jules50_struc(n)%jules50(t_unpert)%smcl_soilt(j)&
                       !	        *1/dzsoil(j)*1/1000,&
                       !         sm_threshold(j))) / (1/dzsoil(j)*1/1000) ! min( [kg/m2]*[1/m]*[1/kg/m3] , [m3/m3] ) / ([1/m][kg/m3]) --> [kg/m2]

                       ens_flag(m) = .false. 
                    endif
                 enddo
              enddo
              



!--------------------------------------------------------------------------
! Recalculate the deltas and adjust the ensemble
!--------------------------------------------------------------------------
              !do j=1,4
              do j=1,jules50_struc(n)%sm_levels  !Yonghwan Kwon
                 delta_u(j) = 0.0
                 delta_f(j) = 0.0
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1

                    !t = (i-1)*LIS_rc%nensem(n)+m
                    if(m.ne.LIS_rc%nensem(n)) then 
		       ! NOTE: be careful with the unit. do not change the unit of delta in the loop
                       ! [m3w/m3s]  + ([-] - [-])*[m3/m3] --> [m3/m3]
                       delta_u(j) = delta_u(j)+ & !  * 1/dzsoil(j)*1/1000 
                            (jules50_struc(n)%jules50(t)%p_s_sthu(j) - jules50_struc(n)%jules50(t_unpert)%p_s_sthu(j)) * &
                                jules50_struc(n)%jules50(t)%p_s_smvcst(j)  !) &
                                ! / (1/dzsoil(j)*1/1000)
                       delta_f(j) = delta_f(j)+ & !  * 1/dzsoil(j)*1/1000 
                            (jules50_struc(n)%jules50(t)%p_s_sthf(j) - jules50_struc(n)%jules50(t_unpert)%p_s_sthf(j)) * &
                                jules50_struc(n)%jules50(t)%p_s_smvcst(j)  !) &
                                ! / (1/dzsoil(j)*1/1000)
                    endif
                 enddo
              enddo
              
              !do j=1,4
              do j=1,jules50_struc(n)%sm_levels  !Yonghwan Kwon
                 delta_u(j) =delta_u(j)/(LIS_rc%nensem(n)-1)  !Yonghwan Kwon
                 delta_f(j) =delta_f(j)/(LIS_rc%nensem(n)-1)  !Yonghwan Kwon
                 delta(j) = delta_u(j) + delta_f(j)  !Yonghwan Kwon
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    
                    MAX_THRESHOLD(j) = jules50_struc(n)%jules50(t)%p_s_smvcst(j) ! assume it is the same for all layers [m3/m3]
                    MIN_THRESHOLD(j) = jules50_struc(n)%jules50(t)%p_s_smvcwt(j) ! Volumetric wilting point (m^3 m-3 of soil) !Yonghwan Kwon 
                    !sm_threshold(j)  = MAX_THRESHOLD(j) - MIN_THRESHOLD(j)
                    sm_threshold(j)  = MAX_THRESHOLD(j)

                    sat_p(j) = jules50_struc(n)%jules50(t)%p_s_smvcst(j) ! Volumetric saturation point (m^3 m-3 of soil) !Yonghwan Kwon
                    p_s_sth(j) = jules50_struc(n)%jules50(t)%p_s_sthu(j) + jules50_struc(n)%jules50(t)%p_s_sthf(j)  !saturated fraction (Yonghwan Kwon)

                    !-----------------------------------------------------
                    ! Frozen soil moisture has been added (Yonghwan Kwon)
                    ! Compute fraction

                    if (p_s_sth(j) > 0) then
                       frac_sthu(j) = jules50_struc(n)%jules50(t)%p_s_sthu(j) / p_s_sth(j)
                       frac_sthf(j) = jules50_struc(n)%jules50(t)%p_s_sthf(j) / p_s_sth(j)
                    else
                       frac_sthu(j) = 0
                       frac_sthf(j) = 0
                    endif
                    !-----------------------------------------------------

                    if(ens_flag(m)) then 
                       tmpval = p_s_sth(j) * sat_p(j) - delta(j) !* 1/dzsoil(j)*1/1000 ! [-][m3/m3]-[m3/m3] --> [m3/m3] 
                       tmpval_u = jules50_struc(n)%jules50(t)%p_s_sthu(j) * sat_p(j) - delta_u(j)
                       tmpval_f = jules50_struc(n)%jules50(t)%p_s_sthf(j) * sat_p(j) - delta_f(j)

!if (t .gt.600 .and. t.lt.613 .and. j==1) then
!WRITE (*, '(A55 , 1x, 6(F10.6,1x) )')'2 tmpval, tmpval_u, tmpval_f, smc, sthu, sthf, MIN, MAX',tmpval,tmpval_u,tmpval_f, &
!jules50_struc(n)%jules50(t)%smcl_soilt(1) *1/dzsoil(1)*1/1000 , & ! MN
!jules50_struc(n)%jules50(t)%p_s_sthu(1) * sat_p(1), &! MN
!jules50_struc(n)%jules50(t)%p_s_sthf(1) * sat_p(1), &! MN
!MIN_THRESHOLD(j) , MAX_THRESHOLD(j)
!endif

                       if(.not.(tmpval.le.0.0 .or.tmpval.gt.(MAX_THRESHOLD(j)))) then 
                          if(.not.(tmpval_u.lt.0.0 .or.tmpval_u.gt.(MAX_THRESHOLD(j)))) then 
                             if(.not.(tmpval_f.lt.0.0 .or.tmpval_f.gt.(MAX_THRESHOLD(j)))) then
                                jules50_struc(n)%jules50(t)%p_s_sthu(j) = &
                                   (jules50_struc(n)%jules50(t)%p_s_sthu(j) * sat_p(j)&
                                    - delta_u(j)) / sat_p(j) !([-][m3/m3]-[m3/m3]) / [m3/m3] --> [-]
                                jules50_struc(n)%jules50(t)%p_s_sthf(j) = &
                                   (jules50_struc(n)%jules50(t)%p_s_sthf(j) * sat_p(j)&
                                    - delta_f(j)) / sat_p(j) !([-][m3/m3]-[m3/m3]) / [m3/m3] --> [-]

                                jules50_struc(n)%jules50(t)%smcl_soilt(j) = &
                                   jules50_struc(n)%jules50(t)%smcl_soilt(j) - delta(j)/(1/dzsoil(j)*1/1000) ! [kg/m2] - [m3/m3] / ([1/m]*[1/kg/m3]) --> [kg/m2]

                                !jules50_struc(n)%jules50(t)%smcl_soilt(j) = &
                                !   (jules50_struc(n)%jules50(t)%p_s_sthu(j) + jules50_struc(n)%jules50(t)%p_s_sthf(j)) * sat_p(j)&
                                !   /(1/dzsoil(j)*1/1000)

                                bounds_violation = .false.
                             endif
                          endif
                       endif           
                    endif

                    p_s_sth(j) = jules50_struc(n)%jules50(t)%p_s_sthu(j) + jules50_struc(n)%jules50(t)%p_s_sthf(j)
                    tmpval = p_s_sth(j) * sat_p(j) ! [-][m3/m3]
                    
                    if(tmpval.le.0.0 .or.&
                         tmpval.gt.(MAX_THRESHOLD(j))) then 
                       bounds_violation = .true.
                    else
                       bounds_violation = .false.
                    endif

                    if (sat_p(j) == 0) then
                       bounds_violation = .false.
                       !glacier grid -> does not conduct assimilation
                    endif

                 enddo
              enddo
              
              if(nIter.gt.10.and.bounds_violation) then

!--------------------------------------------------------------------------
! All else fails, set to the bounds
!--------------------------------------------------------------------------
                 
!                 write(LIS_logunit,*) '[ERR] Ensemble structure violates physical bounds '
!                 write(LIS_logunit,*) '[ERR] Please adjust the perturbation settings ..'
                 !do j=1,4
                 do j=1,jules50_struc(n)%sm_levels !Yonghwan Kwon
                    do m=1,LIS_rc%nensem(n)
                       t = i+m-1
                       !t = (i-1)*LIS_rc%nensem(n)+m
                      
                       MAX_THRESHOLD(j) = jules50_struc(n)%jules50(t)%p_s_smvcst(j) ! assume it is the same for all layers [m3/m3]
                       MIN_THRESHOLD(j) = jules50_struc(n)%jules50(t)%p_s_smvcwt(j) ! Volumetric wilting point (m^3 m-3 of soil) !Yonghwan Kwon 

                       sat_p(j) = jules50_struc(n)%jules50(t)%p_s_smvcst(j) ! Volumetric saturation point (m^3 m-3 of soil) !Yonghwan Kwon
                       p_s_sth(j) = jules50_struc(n)%jules50(t)%p_s_sthu(j) + jules50_struc(n)%jules50(t)%p_s_sthf(j)  !saturated fraction (Yonghwan Kwon)
 
                       !-----------------------------------------------------
                       ! Frozen soil moisture has been added (Yonghwan Kwon)
                       ! Compute fraction

                       if (p_s_sth(j) > 0) then
                          frac_sthu(j) = jules50_struc(n)%jules50(t)%p_s_sthu(j) / p_s_sth(j)
                          frac_sthf(j) = jules50_struc(n)%jules50(t)%p_s_sthf(j) / p_s_sth(j)
                       else
                          frac_sthu(j) = 0
                          frac_sthf(j) = 0
                       endif
                       !-----------------------------------------------------

                       if (p_s_sth(j) * sat_p(j).gt.MAX_THRESHOLD(j).or.&  
                          jules50_struc(n)%jules50(t)%smcl_soilt(j) * &
                             1/dzsoil(j)*1/1000.gt.MAX_THRESHOLD(j)) then
                       
                          jules50_struc(n)%jules50(t)%p_s_sthu(j) = MAX_THRESHOLD(j)*frac_sthu(j) / sat_p(j) ! [m3/m3]/[m3/m3] --> [-] 
                          jules50_struc(n)%jules50(t)%p_s_sthf(j) = MAX_THRESHOLD(j)*frac_sthf(j) / sat_p(j) ! [m3/m3]/[m3/m3] --> [-] 
                          jules50_struc(n)%jules50(t)%smcl_soilt(j) = MAX_THRESHOLD(j) / (1/dzsoil(j)*1/1000)   ! [m3w/m3s] / ([1/m1s][m3w/kg]) --> kg/m2s
                       endif

                       if (p_s_sth(j) * sat_p(j).lt.MIN_THRESHOLD(j).or.&                       
                          jules50_struc(n)%jules50(t)%smcl_soilt(j)* &
                             1/dzsoil(j)*1/1000.lt.MIN_THRESHOLD(j)) then
                     
                          jules50_struc(n)%jules50(t)%p_s_sthu(j) = 0.1*MIN_THRESHOLD(j)*frac_sthu(j) / sat_p(j) ! [m3/m3]/[m3/m3] --> [-]
                          jules50_struc(n)%jules50(t)%p_s_sthf(j) = 0.1*MIN_THRESHOLD(j)*frac_sthf(j) / sat_p(j) ! [m3/m3]/[m3/m3] --> [-]
                          jules50_struc(n)%jules50(t)%smcl_soilt(j) = 0.1*MIN_THRESHOLD(j) / (1/dzsoil(j)*1/1000)   ! [m3w/m3s] / ([1/m1s][m3w/kg]) --> kg/m2s
                       endif
                    enddo
!                 call LIS_endrun()
                 enddo

                 bounds_violation = .false.
              endif
              
           end do

        endif
     endif
  enddo
end subroutine jules50_setsoilm

