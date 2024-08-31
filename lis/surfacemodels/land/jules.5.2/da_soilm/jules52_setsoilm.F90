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
! !ROUTINE: jules52_setsoilm
!  \label{jules52_setsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 20 Dec 2018: Mahdi Navari; Modified for JULES 5.2
!
! There are 3 cases 
! 1- If all the ensemble members met the update conditions --> apply the update
! 2- If more than 50% of the ensemble members met the update condition --> 
!    apply the update for that members and set the other member to the mean 
!    value of the ensemble (i.e. mean of the members that met the conditions)
! 3- If less then 50% of the ensemble members met the update conditions --> 
!    adjust the states  
! !INTERFACE:
subroutine jules52_setsoilm(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use jules52_lsmMod
 !MN: added to use LIS_sfmodel_struc
  !use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
  use jules_soil_mod, only:  dzsoil !
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
  integer                :: t, j,i, gid, m, t_unpert, row, col
  integer                :: status
  real                   :: delta(4)
  real                   :: delta1,delta2,delta3,delta4 , lat, lon
  real                   :: tmpval
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


  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in jules52_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 2 failed in jules52_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 3 failed in jules52_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 4 failed in jules52_setsoilm")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in jules52_setsoilm")
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in jules52_setsoilm")
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in jules52_setsoilm")
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in jules52_setsoilm")

  update_flag = .true. 
  update_flag_tile = .true. 

! MN: NOTE: be careful with the unit. The unit of soil moisture diagnostic variable (smcl) [kg/m2]
! differs from that of the prognostic variable (p_s_sthu, p_s_smvcst) [m3/m3].

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
  
     MAX_THRESHOLD = jules52_struc(n)%jules52(t)%p_s_smvcst(1) ! Volumetric saturation point (m^3 m-3 of soil)
     sm_threshold  = MAX_THRESHOLD-MIN_THRESHOLD       
     
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 
     
     !MN: delta = X(+) - X(-)
     !NOTE: "jules52_updatesoilm.F90" updates the soilmx(t)   
     delta1 = (soilm1(t)-jules52_struc(n)%jules52(t)%smcl_soilt(1))* 1/dzsoil(1)*1/1000  ! [kg/m2]*1/m*1/(kg/m3) --> [m3/m3]
     delta2 = (soilm2(t)-jules52_struc(n)%jules52(t)%smcl_soilt(2))* 1/dzsoil(2)*1/1000 
     delta3 = (soilm3(t)-jules52_struc(n)%jules52(t)%smcl_soilt(3))* 1/dzsoil(3)*1/1000
     delta4 = (soilm4(t)-jules52_struc(n)%jules52(t)%smcl_soilt(4))* 1/dzsoil(4)*1/1000
     
     ! MN: check MIN_THRESHOLD < volumetric liquid soil moisture < threshold 
     ! unit conversion -->   ..%p_s_sthu(1)  *  ..%p_s_smvcst(1)-->   [-] * [m3/m3] 
     if(jules52_struc(n)%jules52(t)%p_s_sthu(1)*jules52_struc(n)%jules52(t)%p_s_smvcst(1)&
	+ delta1.gt.MIN_THRESHOLD .and.&
        jules52_struc(n)%jules52(t)%p_s_sthu(1)*jules52_struc(n)%jules52(t)%p_s_smvcst(1)&
	+delta1.lt.sm_threshold) then 
        update_flag(gid) = update_flag(gid).and.(.true.)
        ! MN save the flag for each tile (col*row*ens)   PILDAS -->(64*44)*20
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     if(jules52_struc(n)%jules52(t)%p_s_sthu(2)*jules52_struc(n)%jules52(t)%p_s_smvcst(2)&
	+delta2.gt.MIN_THRESHOLD .and.&
        jules52_struc(n)%jules52(t)%p_s_sthu(2)*jules52_struc(n)%jules52(t)%p_s_smvcst(2)&
	+delta2.lt.sm_threshold) then 
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     if(jules52_struc(n)%jules52(t)%p_s_sthu(3)*jules52_struc(n)%jules52(t)%p_s_smvcst(3)&
	+delta3.gt.MIN_THRESHOLD .and.&
        jules52_struc(n)%jules52(t)%p_s_sthu(3)*jules52_struc(n)%jules52(t)%p_s_smvcst(3)&
	+delta3.lt.sm_threshold) then 
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     if(jules52_struc(n)%jules52(t)%p_s_sthu(4)*jules52_struc(n)%jules52(t)%p_s_smvcst(4)&
	+delta4.gt.MIN_THRESHOLD .and.&
        jules52_struc(n)%jules52(t)%p_s_sthu(4)*jules52_struc(n)%jules52(t)%p_s_smvcst(4)&
	+delta4.lt.sm_threshold) then 
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
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

           MAX_THRESHOLD = jules52_struc(n)%jules52(t)%p_s_smvcst(1) ! [m3/m3] assumed value for all layers are the same
           sm_threshold  = MAX_THRESHOLD - MIN_THRESHOLD       
                      
           ! MN check update status for each tile  
           if(update_flag_tile(t)) then
              
	      delta1 = (soilm1(t)-jules52_struc(n)%jules52(t)%smcl_soilt(1))* 1/dzsoil(1)*1/1000  ! [kg/m2]*1/m*1/(kg/m3) --> [m3/m3]
	      delta2 = (soilm2(t)-jules52_struc(n)%jules52(t)%smcl_soilt(2))* 1/dzsoil(2)*1/1000 
	      delta3 = (soilm3(t)-jules52_struc(n)%jules52(t)%smcl_soilt(3))* 1/dzsoil(3)*1/1000
	      delta4 = (soilm4(t)-jules52_struc(n)%jules52(t)%smcl_soilt(4))* 1/dzsoil(4)*1/1000   

   
              jules52_struc(n)%jules52(t)%p_s_sthu(1) = &
		(jules52_struc(n)%jules52(t)%p_s_sthu(1)*jules52_struc(n)%jules52(t)%p_s_smvcst(1)&
                 +delta1)/jules52_struc(n)%jules52(t)%p_s_smvcst(1) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction  
              jules52_struc(n)%jules52(t)%smcl_soilt(1) = soilm1(t) ![kg/m2]
              if(soilm1(t).lt.0) then 
                 print*, 'setsoilm1 ',t,soilm1(t)
                 stop
              endif
	      ! I think the following sentence is redundant. 
	      !Because both “update_flag_new” and “update_flag_tile” are TRUE 

              !if(jules52_struc(n)%jules52(t)%p_s_sthu(2)+delta2.gt.MIN_THRESHOLD .and.&   
              !     jules52_struc(n)%jules52(t)%p_s_sthu(2)+delta2.lt.sm_threshold) then 
              !   jules52_struc(n)%jules52(t)%p_s_sthu(2) = jules52_struc(n)%jules52(t)%p_s_sthu(2)+&
              !        soilm2(t)-jules52_struc(n)%jules52(t)%smcl_soilt(2)
              jules52_struc(n)%jules52(t)%p_s_sthu(2) = &
		(jules52_struc(n)%jules52(t)%p_s_sthu(2)*jules52_struc(n)%jules52(t)%p_s_smvcst(2)&
                 +delta2)/jules52_struc(n)%jules52(t)%p_s_smvcst(2) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction  
              jules52_struc(n)%jules52(t)%smcl_soilt(2) = soilm2(t)
              if(soilm2(t).lt.0) then 
                 print*, 'setsoilm2 ',t,soilm2(t)
                   stop
              endif
              !endif
              !if(jules52_struc(n)%jules52(t)%p_s_sthu(3)+delta3.gt.MIN_THRESHOLD .and.&
              !     jules52_struc(n)%jules52(t)%p_s_sthu(3)+delta3.lt.sm_threshold) then 
              !   jules52_struc(n)%jules52(t)%p_s_sthu(3) = jules52_struc(n)%jules52(t)%p_s_sthu(3)+&
              !        soilm3(t)-jules52_struc(n)%jules52(t)%smcl_soilt(3)
              jules52_struc(n)%jules52(t)%p_s_sthu(3) = &
		(jules52_struc(n)%jules52(t)%p_s_sthu(3)*jules52_struc(n)%jules52(t)%p_s_smvcst(3)&
                 +delta3)/jules52_struc(n)%jules52(t)%p_s_smvcst(3) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction  
              jules52_struc(n)%jules52(t)%smcl_soilt(3) = soilm3(t)
              if(soilm3(t).lt.0) then 
                 print*, 'setsoilm3 ',t,soilm3(t)
                   stop
              endif
              !endif
              ! surface layer
              !if(jules52_struc(n)%jules52(t)%p_s_sthu(4)+delta4.gt.MIN_THRESHOLD .and.&
              !     jules52_struc(n)%jules52(t)%p_s_sthu(4)+delta4.lt.sm_threshold) then 
              !   jules52_struc(n)%jules52(t)%p_s_sthu(4) = jules52_struc(n)%jules52(t)%p_s_sthu(4)+&
              !        soilm4(t)-jules52_struc(n)%jules52(t)%smcl_soilt(4)
              jules52_struc(n)%jules52(t)%p_s_sthu(4) = &
		(jules52_struc(n)%jules52(t)%p_s_sthu(4)*jules52_struc(n)%jules52(t)%p_s_smvcst(4)&
                 +delta4)/jules52_struc(n)%jules52(t)%p_s_smvcst(4) ! ([-]*[m3/m3] + [m3/m3])/[m3/m3] --> fraction  
              jules52_struc(n)%jules52(t)%smcl_soilt(4) = soilm4(t)
              if(soilm4(t).lt.0) then 
                 print*, 'setsoilm4 ',t,soilm4(t)
                  stop
              endif
              !endif
           else 
              
!-----------------------------------------------------------------------------------------  
! USE ENSEMBLE MEAN VALUE 
!-----------------------------------------------------------------------------------------
              
              ! Assume p_s_sthu = smc (i.e. ice content=0) 
              smc_tmp = (MaxEnsSM1 - MinEnsSM1)/2 + MinEnsSM1 ! [kg/m2]
              jules52_struc(n)%jules52(t)%p_s_sthu(1) = &
			(smc_tmp * 1/dzsoil(1)*1/1000)/ &
			jules52_struc(n)%jules52(t)%p_s_smvcst(1) ! ([kg/m2]*(1/m*1/(kg/m3)))/[m3/m3]--> fraction
              jules52_struc(n)%jules52(t)%smcl_soilt(1) = smc_tmp
              
              

              smc_tmp = (MaxEnsSM2 - MinEnsSM2)/2 + MinEnsSM2            
              jules52_struc(n)%jules52(t)%p_s_sthu(2) = &
			(smc_tmp * 1/dzsoil(2)*1/1000)/ &
			jules52_struc(n)%jules52(t)%p_s_smvcst(2)
              jules52_struc(n)%jules52(t)%smcl_soilt(2) = smc_tmp
              
              

              smc_tmp = (MaxEnsSM3 - MinEnsSM3)/2 + MinEnsSM3
              jules52_struc(n)%jules52(t)%p_s_sthu(3) = &
			(smc_tmp * 1/dzsoil(3)*1/1000)/ &
			jules52_struc(n)%jules52(t)%p_s_smvcst(3)
              jules52_struc(n)%jules52(t)%smcl_soilt(3) = smc_tmp
              
              

              smc_tmp = (MaxEnsSM4 - MinEnsSM4)/2 + MinEnsSM4
              jules52_struc(n)%jules52(t)%p_s_sthu(4) = &
			(smc_tmp * 1/dzsoil(4)*1/1000)/ &
			jules52_struc(n)%jules52(t)%p_s_smvcst(4)
              jules52_struc(n)%jules52(t)%smcl_soilt(4) = smc_tmp

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
		! NOTE: be careful with the unit. do not change the unit of delta in the loop
              ! [m3w/m3s]  + ([-] - [-])*[m3/m3] 
                       delta(j) = delta(j)+ & !  * 1/dzsoil(j)*1/1000 
                            (jules52_struc(n)%jules52(t)%p_s_sthu(j) - &
                             jules52_struc(n)%jules52(t_unpert)%p_s_sthu(j)) * &
			        jules52_struc(n)%jules52(t)%p_s_smvcst(j)  !) &
			        ! / (1/dzsoil(j)*1/1000)  
                    endif
                    
                 enddo
              enddo
              
              do j=1,4
                 delta(j) = delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                  
		      MAX_THRESHOLD = jules52_struc(n)%jules52(t)%p_s_smvcst(1) ! assume it is the same for all layers [m3/m3]
                    sm_threshold  = MAX_THRESHOLD-MIN_THRESHOLD
                    
                    tmpval = jules52_struc(n)%jules52(t)%p_s_sthu(j) * &
			        jules52_struc(n)%jules52(t)%p_s_smvcst(j) - &
                             delta(j) !* 1/dzsoil(j)*1/1000 ! [-][m3/m3]-[m3/m3] --> [m3/m3]  
                    if(tmpval.le.MIN_THRESHOLD) then 

                       jules52_struc(n)%jules52(t)%p_s_sthu(j) = &
                            (max(jules52_struc(n)%jules52(t_unpert)%p_s_sthu(j)& 
				*jules52_struc(n)%jules52(t)%p_s_smvcst(j),&
                            	 MIN_THRESHOLD)) / jules52_struc(n)%jules52(t)%p_s_smvcst(j)! max( [-][m3/m3] , [m3/m3] ) / [m3/m3] --> fraction

                       jules52_struc(n)%jules52(t)%smcl_soilt(j) = & 
                            (max(jules52_struc(n)%jules52(t_unpert)%smcl_soilt(j)&
			        *1/dzsoil(j)*1/1000,&
                                 MIN_THRESHOLD)) / (1/dzsoil(j)*1/1000) ! max( [kg/m2]*[1/m]*[1/kg/m3] , [m3/m3] ) / ([1/m][kg/m3]) --> [kg/m2]

                       ens_flag(m) = .false. 

                    elseif(tmpval.ge.sm_threshold) then

                       jules52_struc(n)%jules52(t)%p_s_sthu(j) = &
                            (min(jules52_struc(n)%jules52(t_unpert)%p_s_sthu(j)&
				*jules52_struc(n)%jules52(t)%p_s_smvcst(j),&
                            	 sm_threshold)) / jules52_struc(n)%jules52(t)%p_s_smvcst(j) ! min( [-][m3/m3] , [m3/m3] ) / [m3/m3] --> fraction

                       jules52_struc(n)%jules52(t)%smcl_soilt(j) = &
                            (min(jules52_struc(n)%jules52(t_unpert)%smcl_soilt(j)&
			        *1/dzsoil(j)*1/1000,&
                                sm_threshold)) / (1/dzsoil(j)*1/1000) ! min( [kg/m2]*[1/m]*[1/kg/m3] , [m3/m3] ) / ([1/m][kg/m3]) --> [kg/m2]

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
		! NOTE: be careful with the unit. do not change the unit of delta in the loop
              ! [m3w/m3s]  + ([-] - [-])*[m3/m3] --> [m3/m3]
                       delta(j) = delta(j)  + & ! * 1/dzsoil(j)*1/1000
                            (jules52_struc(n)%jules52(t)%p_s_sthu(j) - &
                             jules52_struc(n)%jules52(t_unpert)%p_s_sthu(j)) * &
			        jules52_struc(n)%jules52(t)%p_s_smvcst(j)   !) &
			         !/ (1/dzsoil(j)*1/1000)  
                    endif
                 enddo
              enddo
              
              do j=1,4
                 delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    
                    if(ens_flag(m)) then 
                       tmpval = jules52_struc(n)%jules52(t)%p_s_sthu(j) * &
			           jules52_struc(n)%jules52(t)%p_s_smvcst(j) - &
                                delta(j) !* 1/dzsoil(j)*1/1000 !! [-][m3/m3]-[m3/m3] --> [m3/m3]

		         MAX_THRESHOLD = jules52_struc(n)%jules52(t)%p_s_smvcst(1)
                       if(.not.(tmpval.le.0.0 .or.&
                            tmpval.gt.(MAX_THRESHOLD))) then 
                          
                          jules52_struc(n)%jules52(t)%smcl_soilt(j) = &
                               jules52_struc(n)%jules52(t)%smcl_soilt(j) - delta(j)/(1/dzsoil(j)*1/1000) ! [kg/m2] - [m3/m3] / ([1/m]*[1/kg/m3]) --> [kg/m2]

                          jules52_struc(n)%jules52(t)%p_s_sthu(j) = &
                          (jules52_struc(n)%jules52(t)%p_s_sthu(j)*&
			      jules52_struc(n)%jules52(t)%p_s_smvcst(j)-&
			      delta(j)) /& ! * 1/dzsoil(j)*1/1000)/&
			      jules52_struc(n)%jules52(t)%p_s_smvcst(j) !([-][m3/m3]-[m3/m3]) / [m3/m3] --> [-]

                          bounds_violation = .false.
                       endif
                    endif

           
                    tmpval = jules52_struc(n)%jules52(t)%p_s_sthu(j) * &
			        jules52_struc(n)%jules52(t)%p_s_smvcst(j) ! [-][m3/m3]
                    
		      MAX_THRESHOLD = jules52_struc(n)%jules52(t)%p_s_smvcst(1)                   
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
                 
!                 write(LIS_logunit,*) '[ERR] Ensemble structure violates physical bounds '
!                 write(LIS_logunit,*) '[ERR] Please adjust the perturbation settings ..'
                 do j=1,4
                    do m=1,LIS_rc%nensem(n)
                     t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                       
			  MAX_THRESHOLD = jules52_struc(n)%jules52(t)%p_s_smvcst(1)                       
                       if(jules52_struc(n)%jules52(t)%p_s_sthu(j)* &
			     jules52_struc(n)%jules52(t)%p_s_smvcst(j).gt.MAX_THRESHOLD.or.&
                          jules52_struc(n)%jules52(t)%smcl_soilt(j) * &
			     1/dzsoil(j)*1/1000.gt.MAX_THRESHOLD) then 
                       
                          jules52_struc(n)%jules52(t)%p_s_sthu(j) = MAX_THRESHOLD / &
			     jules52_struc(n)%jules52(t)%p_s_smvcst(j) ! [m3/m3]/[m3/m3] --> [-]
                          jules52_struc(n)%jules52(t)%smcl_soilt(j) = MAX_THRESHOLD / & 
			     (1/dzsoil(j)*1/1000)   ! [m3w/m3s] / ([1/m1s][m3w/kg]) --> kg/m2s
                       endif
                       
                       if(jules52_struc(n)%jules52(t)%p_s_sthu(j)* &
			     jules52_struc(n)%jules52(t)%p_s_smvcst(j).lt.MIN_THRESHOLD.or.&
                          jules52_struc(n)%jules52(t)%smcl_soilt(j)* & 
			     1/dzsoil(j)*1/1000.lt.MIN_THRESHOLD) then   
                     
                          jules52_struc(n)%jules52(t)%p_s_sthu(j) = MIN_THRESHOLD / &
			     jules52_struc(n)%jules52(t)%p_s_smvcst(j) ! [m3/m3]/[m3/m3] --> [-] 
                          jules52_struc(n)%jules52(t)%smcl_soilt(j) = MIN_THRESHOLD / & 
			     (1/dzsoil(j)*1/1000)   ! [m3w/m3s] / ([1/m1s][m3w/kg]) --> kg/m2s
                       endif
                    print*, i, m
                    print*, '2smc',t, jules52_struc(n)%jules52(t)%smcl_soilt(:)
                    print*, '2p_s_sthu ',t,jules52_struc(n)%jules52(t)%p_s_sthu(:)
                    print*, '2max ',t,MAX_THRESHOLD !jules52_struc(n)%jules52(t)%smcmax
                    enddo
!                 call LIS_endrun()
                 enddo
              endif
              
           end do
        endif
     endif
  enddo
 
end subroutine jules52_setsoilm

