!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP401_setsoilm
!  \label{NoahMP401_setsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari: Modified for NoahMP401 
! 
! Apply the update if it met the update conditions
! Update conditions: 
!                  1- Prior SM(sh2o) + increment > MIN_THRESHOLD 
!                  2- Prior SM(sh2o) + increment < sm_threshold
! There are 3 cases 
! 1- If all the ensemble members met the update conditions --> apply the update
! 2- If more than 50% of the ensemble members met the update condition --> 
!    apply the update for that members and set the other member to the mean 
!    value of the ensemble (i.e. mean of the members that met the conditions)
! 3- If less then 50% of the ensemble members met the update conditions --> 
!    adjust the states    


! !INTERFACE:
subroutine NoahMP401_setsoilm(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use NoahMP401_lsmMod
  !use module_sf_noahlsm_36  !MN
  !use module_sf_noahmpdrv_401, only: parameters
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
  real, parameter        :: MIN_THRESHOLD = 0.02 
  real                   :: MAX_threshold
  real                   :: sm_threshold
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  integer                :: t, j,i, gid, m, t_unpert
  integer                :: status
  real                   :: delta(4)
  real                   :: delta1,delta2,delta3,delta4
  real                   :: tmpval
  logical                :: bounds_violation
  integer                :: nIter
  logical                :: update_flag(LIS_rc%ngrid(n))
  logical                :: ens_flag(LIS_rc%nensem(n))
! mn
  integer                :: SOILTYP           ! soil type index [-]
  real                   :: SMCMAX , SMCWLT
  real                   :: tmp(LIS_rc%nensem(n)), tmp0(LIS_rc%nensem(n))
  real                   :: tmp1(LIS_rc%nensem(n)),tmp2(LIS_rc%nensem(n)),tmp3(LIS_rc%nensem(n)),tmp4(LIS_rc%nensem(n)) 
  logical                :: update_flag_tile(LIS_rc%npatch(n,LIS_rc%lsm_index))
  logical                :: flag_ens(LIS_rc%ngrid(n))
  logical                :: flag_tmp(LIS_rc%nensem(n))
  logical                :: update_flag_ens(LIS_rc%ngrid(n))
  logical                :: update_flag_new(LIS_rc%ngrid(n))
  integer                :: RESULT, pcount, icount
  real                   :: MaxEnsSM1 ,MaxEnsSM2 ,MaxEnsSM3 ,MaxEnsSM4
  real                   :: MinEnsSM1 ,MinEnsSM2 ,MinEnsSM3 ,MinEnsSM4 
  real                   :: MaxEns_sh2o1, MaxEns_sh2o2, MaxEns_sh2o3, MaxEns_sh2o4
  real                   :: smc_rnd, smc_tmp 
  real                   :: sh2o_tmp, sh2o_rnd 
  INTEGER, DIMENSION (1) :: seed 

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in NoahMP401_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 2 failed in NoahMP401_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 3 failed in NoahMP401_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 4 failed in NoahMP401_setsoilm")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in NoahMP401_setsoilm")
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in NoahMP401_setsoilm")
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in NoahMP401_setsoilm")
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in NoahMP401_setsoilm")

  update_flag = .true. 
  update_flag_tile= .true. 

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
  
 ! MN: NOTE: SMCMAX and SMCWLT are not stored in the data structure but we 
 !       can get module variables MAXSMC and WLTSMC from the module_sf_noahlsm_36    
 !       In NoahMP401 module variables MAXSMC and WLTSMC from module_sf_noahmpdrv_401      
     SOILTYP = noahmp401_struc(n)%noahmp401(t)%soiltype        
     MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP)   ! MAXSMC (SOILTYP)
     sm_threshold = SMCMAX_TABLE(SOILTYP) - 0.02  ! MAXSMC (SOILTYP) - 0.02
     
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 
     
     !MN: delta = X(+) - X(-)
     !NOTE: "noahmp401_updatesoilm.F90" updates the soilm_(t)     
     delta1 = soilm1(t)-noahmp401_struc(n)%noahmp401(t)%smc(1)
     delta2 = soilm2(t)-noahmp401_struc(n)%noahmp401(t)%smc(2)
     delta3 = soilm3(t)-noahmp401_struc(n)%noahmp401(t)%smc(3)
     delta4 = soilm4(t)-noahmp401_struc(n)%noahmp401(t)%smc(4)

     ! MN: check    MIN_THRESHOLD < volumetric liquid soil moisture < threshold 
     if(noahmp401_struc(n)%noahmp401(t)%sh2o(1)+delta1.gt.MIN_THRESHOLD .and.&
          noahmp401_struc(n)%noahmp401(t)%sh2o(1)+delta1.lt.&
          sm_threshold) then 
        update_flag(gid) = update_flag(gid).and.(.true.)
        ! MN save the flag for each tile (col*row*ens)   (64*44)*20
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     if(noahmp401_struc(n)%noahmp401(t)%sh2o(2)+delta2.gt.MIN_THRESHOLD .and.&
          noahmp401_struc(n)%noahmp401(t)%sh2o(2)+delta2.lt.sm_threshold) then 
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     if(noahmp401_struc(n)%noahmp401(t)%sh2o(3)+delta3.gt.MIN_THRESHOLD .and.&
          noahmp401_struc(n)%noahmp401(t)%sh2o(3)+delta3.lt.sm_threshold) then 
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     if(noahmp401_struc(n)%noahmp401(t)%sh2o(4)+delta4.gt.MIN_THRESHOLD .and.&
          noahmp401_struc(n)%noahmp401(t)%sh2o(4)+delta4.lt.sm_threshold) then 
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
  enddo

!-----------------------------------------------------------------------------------------
! MN create new falg: if update falg for 50% of the ensemble members is true 
! then update the state variabels 
!-----------------------------------------------------------------------------------------
  update_flag_ens = .True.
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
  
  ! MN print 
#if 0   
if(i.eq.66) then !i=66  ! --> domain's center  1376
  if(LIS_rc%hr.eq.12) then
     write(2001,'(I4, 2x, 3(I2,x), 2x, 23(L1,2x))'),&
          i, LIS_rc%mo, LIS_rc%da, LIS_rc%hr,update_flag_tile&
          ((i-1)*LIS_rc%nensem(n)+1:(i)*LIS_rc%nensem(n)),&
          update_flag_ens(i), update_flag_new(i), update_flag(i) 
  endif !mn
  endif
#endif 
  
  ! update step
  ! loop over grid points 
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row) 
     
     !if(update_flag(gid)) then
     if(update_flag_new(gid)) then 
!-----------------------------------------------------------------------------------------
       ! 1- update the states
       ! 1-1- if update flag for tile is TRUE --> apply the DA update    
       ! 1-2- if update flag for tile is FALSE --> resample the states  
       ! 2- adjust the states
!-----------------------------------------------------------------------------------------
        ! store update value for  cases that update_flag_tile & update_flag_new are TRUE        ! update_flag_tile = TRUE --> means met the both min and max threshold 
      
        tmp1 = LIS_rc%udef
        tmp2 = LIS_rc%udef
        tmp3 = LIS_rc%udef
        tmp4 = LIS_rc%udef
	!icount = 1
        do m=1,LIS_rc%nensem(n)
           t = i+m-1
           !t = (i-1)*LIS_rc%nensem(n)+m
 
	  if(update_flag_tile(t)) then
          
           tmp1(m) = soilm1(t) 
           tmp2(m) = soilm2(t) 
           tmp3(m) = soilm3(t) 
           tmp4(m) = soilm4(t) 
           !icount = icount + 1 
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
              MaxEnsSM1 = max(MaxEnsSM1, tmp1(m))
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
           
           ! MN check update status for each tile  
           if(update_flag_tile(t)) then
              
              delta1 = soilm1(t)-noahmp401_struc(n)%noahmp401(t)%smc(1)
              delta2 = soilm2(t)-noahmp401_struc(n)%noahmp401(t)%smc(2)
              delta3 = soilm3(t)-noahmp401_struc(n)%noahmp401(t)%smc(3)
              delta4 = soilm4(t)-noahmp401_struc(n)%noahmp401(t)%smc(4)

              noahmp401_struc(n)%noahmp401(t)%sh2o(1) = noahmp401_struc(n)%noahmp401(t)%sh2o(1)+&
                   delta1
              noahmp401_struc(n)%noahmp401(t)%smc(1) = soilm1(t)
              if(soilm1(t).lt.0) then 
                 print*, 'setsoilm1 ',t,soilm1(t)
                 stop
              endif
              if(noahmp401_struc(n)%noahmp401(t)%sh2o(2)+delta2.gt.MIN_THRESHOLD .and.&
                   noahmp401_struc(n)%noahmp401(t)%sh2o(2)+delta2.lt.sm_threshold) then 
                 noahmp401_struc(n)%noahmp401(t)%sh2o(2) = noahmp401_struc(n)%noahmp401(t)%sh2o(2)+&
                      soilm2(t)-noahmp401_struc(n)%noahmp401(t)%smc(2)
                 noahmp401_struc(n)%noahmp401(t)%smc(2) = soilm2(t)
                 if(soilm2(t).lt.0) then 
                    print*, 'setsoilm2 ',t,soilm2(t)
                    stop
                 endif
              endif
  
              if(noahmp401_struc(n)%noahmp401(t)%sh2o(3)+delta3.gt.MIN_THRESHOLD .and.&
                   noahmp401_struc(n)%noahmp401(t)%sh2o(3)+delta3.lt.sm_threshold) then 
                 noahmp401_struc(n)%noahmp401(t)%sh2o(3) = noahmp401_struc(n)%noahmp401(t)%sh2o(3)+&
                      soilm3(t)-noahmp401_struc(n)%noahmp401(t)%smc(3)
                 noahmp401_struc(n)%noahmp401(t)%smc(3) = soilm3(t)
                 if(soilm3(t).lt.0) then 
                    print*, 'setsoilm3 ',t,soilm3(t)
                    stop
                 endif
              endif

              if(noahmp401_struc(n)%noahmp401(t)%sh2o(4)+delta4.gt.MIN_THRESHOLD .and.&
                   noahmp401_struc(n)%noahmp401(t)%sh2o(4)+delta4.lt.sm_threshold) then 
                 noahmp401_struc(n)%noahmp401(t)%sh2o(4) = noahmp401_struc(n)%noahmp401(t)%sh2o(4)+&
                      soilm4(t)-noahmp401_struc(n)%noahmp401(t)%smc(4)
                 noahmp401_struc(n)%noahmp401(t)%smc(4) = soilm4(t)

                 if(soilm4(t).lt.0) then 
                    print*, 'setsoilm4 ',t,soilm4(t)
                    stop
                 endif
              endif
              
              
!-----------------------------------------------------------------------------------------              
              ! randomly resample the smc from [MIN_THRESHOLD,  Max value from DA @ that tiem step]
!-----------------------------------------------------------------------------------------
           else 
          
!-----------------------------------------------------------------------------------------  
! set the soil moisture to the ensemble mean  
!-----------------------------------------------------------------------------------------
              
              ! use mean value
              ! Assume sh2o = smc (i.e. ice content=0) 
              smc_tmp = (MaxEnsSM1 - MinEnsSM1)/2 + MinEnsSM1
              noahmp401_struc(n)%noahmp401(t)%sh2o(1) = smc_tmp 
              noahmp401_struc(n)%noahmp401(t)%smc(1) = smc_tmp
                          
              smc_tmp = (MaxEnsSM2 - MinEnsSM2)/2 + MinEnsSM2            
              noahmp401_struc(n)%noahmp401(t)%sh2o(2) = smc_tmp
              noahmp401_struc(n)%noahmp401(t)%smc(2) = smc_tmp
               
              smc_tmp = (MaxEnsSM3 - MinEnsSM3)/2 + MinEnsSM3
              noahmp401_struc(n)%noahmp401(t)%sh2o(3) = smc_tmp
              noahmp401_struc(n)%noahmp401(t)%smc(3) = smc_tmp
              
              smc_tmp = (MaxEnsSM4 - MinEnsSM4)/2 + MinEnsSM4
              noahmp401_struc(n)%noahmp401(t)%sh2o(4) = smc_tmp
              noahmp401_struc(n)%noahmp401(t)%smc(4) = smc_tmp
              
 
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
           nIter = 0
           ens_flag = .true. 
           
           do while(bounds_violation) 
              niter = niter + 1
              !t_unpert = i*LIS_rc%nensem(n)
	      t_unpert = i+LIS_rc%nensem(n)-1
              do j=1,4
                 delta(j) = 0.0
                 do m=1,LIS_rc%nensem(n)-1
                     t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    
                    if(m.ne.LIS_rc%nensem(n)) then 
                       delta(j) = delta(j) + &
                            (noahmp401_struc(n)%noahmp401(t)%sh2o(j) - &
                            noahmp401_struc(n)%noahmp401(t_unpert)%sh2o(j))
                    endif
                    
                 enddo
              enddo
              
              do j=1,4
                 delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                     t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    SOILTYP = noahmp401_struc(n)%noahmp401(t)%soiltype  
                    MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP) 
                    sm_threshold = SMCMAX_TABLE(SOILTYP)  - 0.02
                    
                    tmpval = noahmp401_struc(n)%noahmp401(t)%sh2o(j) - &
                         delta(j)
                    if(tmpval.le.MIN_THRESHOLD) then 
                       noahmp401_struc(n)%noahmp401(t)%sh2o(j) = &
                            max(noahmp401_struc(n)%noahmp401(t_unpert)%sh2o(j),&
                            MIN_THRESHOLD)
                       noahmp401_struc(n)%noahmp401(t)%smc(j) = &
                            max(noahmp401_struc(n)%noahmp401(t_unpert)%smc(j),&
                            MIN_THRESHOLD)
                       ens_flag(m) = .false. 
                    elseif(tmpval.ge.sm_threshold) then
                       noahmp401_struc(n)%noahmp401(t)%sh2o(j) = &
                            min(noahmp401_struc(n)%noahmp401(t_unpert)%sh2o(j),&
                            sm_threshold)
                       noahmp401_struc(n)%noahmp401(t)%smc(j) = &
                            min(noahmp401_struc(n)%noahmp401(t_unpert)%smc(j),&
                            sm_threshold)
                       ens_flag(m) = .false. 
                    endif
                 enddo
              enddo
              
              !--------------------------------------------------------------------------
              ! Recalculate the deltas and adjust the ensemble
              !--------------------------------------------------------------------------
              do j=1,4
                 delta(j) = 0.0
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    if(m.ne.LIS_rc%nensem(n)) then 
                       delta(j) = delta(j) + &
                            (noahmp401_struc(n)%noahmp401(t)%sh2o(j) - &
                            noahmp401_struc(n)%noahmp401(t_unpert)%sh2o(j))
                    endif
                 enddo
              enddo
              
              do j=1,4
                 delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    
                    if(ens_flag(m)) then 
                       tmpval = noahmp401_struc(n)%noahmp401(t)%sh2o(j) - &
                            delta(j)
                       SOILTYP = noahmp401_struc(n)%noahmp401(t)%soiltype  
                       MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP)  
                       
                       if(.not.(tmpval.le.0.0 .or.&
                            tmpval.gt.(MAX_THRESHOLD))) then 
                          
                          noahmp401_struc(n)%noahmp401(t)%smc(j) = &
                               noahmp401_struc(n)%noahmp401(t)%smc(j) - delta(j)
                          noahmp401_struc(n)%noahmp401(t)%sh2o(j) = &
                               noahmp401_struc(n)%noahmp401(t)%sh2o(j) - delta(j)
                          bounds_violation = .false.
                       endif
                    endif
                    
                    tmpval = noahmp401_struc(n)%noahmp401(t)%sh2o(j)
                    SOILTYP = noahmp401_struc(n)%noahmp401(t)%soiltype  
                    MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP)  
                    
                    if(tmpval.le.0.0 .or.&
                         tmpval.gt.(MAX_THRESHOLD)) then 
                       bounds_violation = .true. 
                    else
                       bounds_violation = .false.
                    endif
                 enddo
              enddo
              
              if(nIter.gt.10.and.bounds_violation) then 
                 !--------------------------------------------------------------------------
                 ! All else fails, set to the bounds
                 !--------------------------------------------------------------------------
                 
                 write(LIS_logunit,*) '[ERR] Ensemble structure violates physical bounds '
                 write(LIS_logunit,*) '[ERR] Please adjust the perturbation settings ..'

                 do j=1,4
                    do m=1,LIS_rc%nensem(n)
                       t = i+m-1
                       !t = (i-1)*LIS_rc%nensem(n)+m
                       
                       SOILTYP = noahmp401_struc(n)%noahmp401(t)%soiltype  
                       MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP)           
                       
                       if(noahmp401_struc(n)%noahmp401(t)%sh2o(j).gt.MAX_THRESHOLD.or.&
                            noahmp401_struc(n)%noahmp401(t)%smc(j).gt.MAX_THRESHOLD) then                        
                          noahmp401_struc(n)%noahmp401(t)%sh2o(j) = MAX_THRESHOLD
                          noahmp401_struc(n)%noahmp401(t)%smc(j) = MAX_THRESHOLD
                       endif
                       
                       if(noahmp401_struc(n)%noahmp401(t)%sh2o(j).lt.MIN_THRESHOLD.or.&
                            noahmp401_struc(n)%noahmp401(t)%smc(j).lt.MIN_THRESHOLD) then                        
                          noahmp401_struc(n)%noahmp401(t)%sh2o(j) = MIN_THRESHOLD
                          noahmp401_struc(n)%noahmp401(t)%smc(j) = MIN_THRESHOLD
                       endif
!                    print*, i, m
!                    print*, 'smc',t, noahmp401_struc(n)%noahmp401(t)%smc(:)
!                    print*, 'sh2o ',t,noahmp401_struc(n)%noahmp401(t)%sh2o(:)
!                    print*, 'max ',t,MAX_THRESHOLD !noahmp401_struc(n)%noahmp401(t)%smcmax
                    enddo
!                  call LIS_endrun()
                 enddo
              endif          
           end do
        endif        
     endif
  enddo


end subroutine NoahMP401_setsoilm

