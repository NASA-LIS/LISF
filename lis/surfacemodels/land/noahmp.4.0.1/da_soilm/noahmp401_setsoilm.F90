!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_setsoilm
!  \label{noahmp401_setsoilm}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li; created for Noah-MP4.0.1 
!
! !INTERFACE:
subroutine noahmp401_setsoilm(n, LSM_State)
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

  real                   :: MIN_THRESHOLD 
  real                   :: MAX_THRESHOLD
  real                   :: sm_threshold
  type(ESMF_Field)       :: sm1Field

  real, pointer          :: soilm1(:)

  real                   :: delta,delta1
  logical                :: diffCheck(LIS_rc%ngrid(n))
  logical                :: ensCheck(LIS_rc%ngrid(n))
  logical                :: largeSM(LIS_rc%ngrid(n))
  real                   :: snodens(LIS_rc%npatch(n,LIS_rc%lsm_index))
  integer                :: i, c,r,t,m,gid
  integer                :: SOILTYP           ! soil type index [-]
  real                   :: sh2o_tmp, sh2o_rnd 
  real                   :: dsneqv,dsnowh, snowh_new
  integer                :: status
  logical                :: update_flag(LIS_rc%ngrid(n))
  logical                :: rc1,rc2,rc3,rc4,rc5

  
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in noahmp401_setsoilm")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmp401_setsoilm")

  update_flag = .true. 

!  if(LIS_localPet.eq.387) then
!     gid = LIS_domain(n)%gindex(&
!          LIS_surface(n,LIS_rc%lsm_index)%tile(16068)%col,&
!          LIS_surface(n,LIS_rc%lsm_index)%tile(16068)%row) 
!     print*, 'set1 ',NOAHMP401_struc(n)%noahmp401(16068)%smc,&
!          NOAHMP401_struc(n)%noahmp401(16068)%sh2o,&
!          NOAHMP401_struc(n)%noahmp401(16068)%sneqv,&
!          update_flag(gid)
!  endif
  

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
     endif
  enddo
  
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     i = LIS_domain(n)%gindex(c,r)
     
     SOILTYP = NOAHMP401_struc(n)%noahmp401(t)%soiltype        
     MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP) 
     MIN_THRESHOLD = SMCWLT_TABLE(SOILTYP)
     sm_threshold = MAX_THRESHOLD - 0.02
     
     if((soilm1(t).lt.MIN_THRESHOLD.or.&
          soilm1(t).gt.MAX_THRESHOLD)) then
        ensCheck(i) = .false.
     endif
     if((soilm1(t).ne.soilm1(i*LIS_rc%nensem(n)))) then
        diffCheck(i) = .true.
     endif
  enddo

  do i=1,LIS_rc%ngrid(n)
     rc1 = .true.
     if(.not.ensCheck(i).and.diffCheck(i).and.(.not.largeSM(i))) then
        call noahmp401_sm_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             soilm1((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_THRESHOLD, MAX_THRESHOLD,rc1)
     endif
     if(.not.rc1) then 

        do m=1,LIS_rc%nensem(n)
           t = (i-1)*LIS_rc%nensem(n)+m
           soilm1(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
        enddo
     endif
  enddo
        
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
       
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 
     
     delta1 = noahmp401_struc(n)%noahmp401(t)%smc(1) - &
          noahmp401_struc(n)%noahmp401(t)%sh2o(1)

     !if ice layer is present, do not update
     if(abs(delta1).ne.0) then 
        update_flag(gid) = update_flag(gid).and.(.false.)
     else
        update_flag(gid) = update_flag(gid).and.(.true.)
     endif
  enddo  

  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 
     
     if(update_flag(gid)) then 
        delta = soilm1(t) - NOAHMP401_struc(n)%noahmp401(t)%smc(1)    
        NOAHMP401_struc(n)%noahmp401(t)%smc(1) = soilm1(t)
        NOAHMP401_struc(n)%noahmp401(t)%sh2o(1) = &
             NOAHMP401_struc(n)%noahmp401(t)%sh2o(1) + delta
     endif

  enddo

!  if(LIS_localPet.eq.387) then
!     gid = LIS_domain(n)%gindex(&
!          LIS_surface(n,LIS_rc%lsm_index)%tile(16068)%col,&
!          LIS_surface(n,LIS_rc%lsm_index)%tile(16068)%row) 
!     print*, 'set2 ',NOAHMP401_struc(n)%noahmp401(16068)%smc,&
!          NOAHMP401_struc(n)%noahmp401(16068)%sh2o,&
!          NOAHMP401_struc(n)%noahmp401(16068)%sneqv,&
!          update_flag(gid)
!  endif
  
     
!  write(101,fmt='(I4.4, 1x, I2.2, 1x, I2.2, 1x, I2.2, 1x, I2.2,1x,10E14.6)') &
!       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr,LIS_rc%mn,&
!       NOAHMP401_struc(n)%noahmp401(991:1000)%sneqv     
end subroutine noahmp401_setsoilm


subroutine noahmp401_sm_reorderEnsForOutliers(i,nensem, statevec, &
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

end subroutine noahmp401_sm_reorderEnsForOutliers
