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
! !ROUTINE: clsmf25_setsnowvars
!  \label{clsmf25_setsnowvars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 29Sep2011: Ben Zaitchik: Applied to GRACE
!
! !INTERFACE:
subroutine clsmf25_setsnowvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only  : LIS_logunit,LIS_verify
  use clsmf25_lsmMod
  use clsmf25_constants
  use clsmf25_model

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the TWS prognostic variables to CLSM's
!  model space. 
! 
!EOP
!  real, parameter        :: WEMIN  = 13.0   ! [KG/M2]
  type(ESMF_Field)       :: wesn1Field
  type(ESMF_Field)       :: wesn2Field
  type(ESMF_Field)       :: wesn3Field
  integer                :: t
  integer                :: status
  real, pointer          :: wesn1(:)
  real, pointer          :: wesn2(:)
  real, pointer          :: wesn3(:)
  real                   :: delta1,delta2,delta3,delta4,delta5,dz
  integer                :: i,j,k
  real, dimension(3,2)   :: ds
  real, dimension(4)     :: sdold, sdnew
  real, dimension(4,2)   :: h,s
  real                   :: snow_density,areasc,areasc1  
!  integer, parameter ::  N_snow = 3
!  real, parameter :: dz1max = 0.05   ! [m]
  real, parameter :: small  = 1.e-20
  real            :: swediff
  real                   :: check_val 
  logical                :: upd_check
  real                   :: wesn_minus(LIS_rc%npatch(n,LIS_rc%lsm_index),3)
  real                   :: sndz_minus(LIS_rc%npatch(n,LIS_rc%lsm_index),3)
  real                   :: htsn_minus(LIS_rc%npatch(n,LIS_rc%lsm_index),3)
  real                   :: catdef_minus(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                   :: srfexc_minus(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                   :: rzexc_minus(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real, dimension(LIS_rc%npatch(n,LIS_rc%lsm_index),3) :: fices
  real                   :: swe_pchange, sndz_pchange

  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 1",wesn1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 2",wesn2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 3",wesn3Field,rc=status)
  call LIS_verify(status) 

  call ESMF_FieldGet(wesn1Field,localDE=0,farrayPtr=wesn1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn2Field,localDE=0,farrayPtr=wesn2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn3Field,localDE=0,farrayPtr=wesn3,rc=status)
  call LIS_verify(status) 

  swe_pchange = 0.0
  sndz_pchange = 0.0
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     wesn_minus(t,:) = clsmf25_struc(n)%cat_progn(t)%wesn(:)
     sndz_minus(t,:) = clsmf25_struc(n)%cat_progn(t)%sndz(:)
     htsn_minus(t,:) = clsmf25_struc(n)%cat_progn(t)%htsn(:)
     catdef_minus(t)  = clsmf25_struc(n)%cat_progn(t)%catdef
     rzexc_minus(t)  = clsmf25_struc(n)%cat_progn(t)%rzexc
     srfexc_minus(t)  = clsmf25_struc(n)%cat_progn(t)%srfexc

    
!     if(LIS_localPet.eq.147.and.t.eq.1882) then 
!        print*, 'in set b ', clsmf25_struc(n)%cat_progn(t)%catdef,&
!             clsmf25_struc(n)%cat_progn(t)%rzexc, &
!             clsmf25_struc(n)%cat_progn(t)%srfexc, &
!             clsmf25_struc(n)%cat_progn(t)%sndz, &
!             clsmf25_struc(n)%cat_progn(t)%wesn
!     endif
   
     areasc = min(1.0,(clsmf25_struc(n)%cat_progn(t)%wesn(1)+&
          clsmf25_struc(n)%cat_progn(t)%wesn(2)+&
          clsmf25_struc(n)%cat_progn(t)%wesn(3))/wemin)
     areasc1 = min(1.0,((wesn1(t)+wesn2(t)+wesn3(t))/wemin))
     
     swediff = (wesn1(t)+wesn2(t)+wesn3(t)-&
          clsmf25_struc(n)%cat_progn(t)%wesn(1)-&
          clsmf25_struc(n)%cat_progn(t)%wesn(2)-&
          clsmf25_struc(n)%cat_progn(t)%wesn(3))
!     if(t.ge.1.and.t.le.10) then 
!        write(LIS_logunit,*) 'set ',catdef(t), srfexc(t), rzexc(t)
!     endif
!#if 0
     if((wesn1(t) + wesn2(t) + wesn3(t)) > 20.0 .and.& 
          areasc.gt.0 .and. & 
          wesn1(t).gt.0.and.&
          wesn2(t).gt.0.and.&
          wesn3(t).gt.0.and.&
          clsmf25_struc(n)%cat_progn(t)%sndz(1).gt.0.and.&
          clsmf25_struc(n)%cat_progn(t)%sndz(2).gt.0.and.&
          clsmf25_struc(n)%cat_progn(t)%sndz(3).gt.0.and.&
          abs(swediff) < 500.0) then 

        if(wesn1(t).gt.2.0) then

           snow_density = clsmf25_struc(n)%cat_progn(t)%wesn(1)/&
                (clsmf25_struc(n)%cat_progn(t)%sndz(1)*areasc)
           clsmf25_struc(n)%cat_progn(t)%htsn(1) = &
                clsmf25_struc(n)%cat_progn(t)%htsn(1) &
                * wesn1(t)/clsmf25_struc(n)%cat_progn(t)%wesn(1)
           clsmf25_struc(n)%cat_progn(t)%sndz(1) =  &
                clsmf25_struc(n)%cat_progn(t)%wesn(1)/&
                (snow_density*areasc1)

           clsmf25_struc(n)%cat_progn(t)%wesn(1) = wesn1(t)

        else
           clsmf25_struc(n)%cat_progn(t)%catdef &
                = clsmf25_struc(n)%cat_progn(t)%catdef -&
                (wesn1(t) - clsmf25_struc(n)%cat_progn(t)%wesn(1))
        end if
        
        if(wesn2(t).gt.2.0) then
!           print*, t, wesn2(t), clsmf25_struc(n)%cat_progn(t)%wesn(2),&
!                clsmf25_struc(n)%cat_progn(t)%sndz(2),areasc
           snow_density = clsmf25_struc(n)%cat_progn(t)%wesn(2)/&
                (clsmf25_struc(n)%cat_progn(t)%sndz(2)*areasc)
           clsmf25_struc(n)%cat_progn(t)%htsn(2) = &
                clsmf25_struc(n)%cat_progn(t)%htsn(2) &
                * wesn2(t)/clsmf25_struc(n)%cat_progn(t)%wesn(2)
           clsmf25_struc(n)%cat_progn(t)%sndz(2) =  &
                clsmf25_struc(n)%cat_progn(t)%wesn(2)/&
                (snow_density*areasc1)
           clsmf25_struc(n)%cat_progn(t)%wesn(2) = wesn2(t)

!           clsmf25_struc(n)%cat_progn(t)%sndz(2) = clsmf25_struc(n)%cat_progn(t)%sndz(2) &
!                * wesn2(t)/clsmf25_struc(n)%cat_progn(t)%wesn(2)       
        else
           ! write(LIS_logunit,*)'dried to wesn2 low limit, tile',t,'wesn2 is',wesn2(t)
           clsmf25_struc(n)%cat_progn(t)%catdef &
                = clsmf25_struc(n)%cat_progn(t)%catdef - &
                (wesn2(t)-clsmf25_struc(n)%cat_progn(t)%wesn(2))
        end if
        
        if(wesn3(t).gt.2.0) then
           snow_density = clsmf25_struc(n)%cat_progn(t)%wesn(3)/&
                (clsmf25_struc(n)%cat_progn(t)%sndz(3)*areasc)
           clsmf25_struc(n)%cat_progn(t)%htsn(3) = &
                clsmf25_struc(n)%cat_progn(t)%htsn(3) &
                * wesn3(t)/clsmf25_struc(n)%cat_progn(t)%wesn(3)
           clsmf25_struc(n)%cat_progn(t)%sndz(3) = &
                clsmf25_struc(n)%cat_progn(t)%wesn(3)/&
                (snow_density*areasc1)
           clsmf25_struc(n)%cat_progn(t)%wesn(3) = wesn3(t)

!           clsmf25_struc(n)%cat_progn(t)%sndz(3) = clsmf25_struc(n)%cat_progn(t)%sndz(3) &
!                * wesn3(t)/clsmf25_struc(n)%cat_progn(t)%wesn(3)        
           ! write(LIS_logunit,*)'updated snow3 in tile',t
        else
           ! write(LIS_logunit,*)'dried to wesn3 low limit, tile',t,'wesn3 is',wesn3(t)
           clsmf25_struc(n)%cat_progn(t)%catdef &
                = clsmf25_struc(n)%cat_progn(t)%catdef - &
                (wesn3(t) - clsmf25_struc(n)%cat_progn(t)%wesn(3))
        end if
        
     ! relayer in order to avoid exeeding max sndz in top snow layer

     !**** Initialize some variables.
        
        h  = 0.
        s  = 0.
        ds = 0.
        dz = 0.
        
        !**** Compute specific heat & water contents of old layers.
        
        do i=1,N_snow
           if (clsmf25_struc(n)%cat_progn(t)%sndz(i) > 0.) then
              h(i,1) = &
                   clsmf25_struc(n)%cat_progn(t)%htsn(i)/clsmf25_struc(n)%cat_progn(t)%sndz(i)
              h(i,2) =  &
                   clsmf25_struc(n)%cat_progn(t)%wesn(i)/clsmf25_struc(n)%cat_progn(t)%sndz(i)
           endif
        enddo

        !**** Obtain old & new layer thicknesses & boundaries.
        
        sdold = 0.
        sdnew = 0.
        
        do i=N_snow,1,-1
           sdold(i) = sdold(i+1) + clsmf25_struc(n)%cat_progn(t)%sndz(i)
        enddo
        
        clsmf25_struc(n)%cat_progn(t)%sndz = sdold(1)/float(N_snow)
        if(clsmf25_struc(n)%cat_progn(t)%sndz(1) > dz1max) then
           clsmf25_struc(n)%cat_progn(t)%sndz(2:) = (sdold(1)-dz1max)/float(N_snow-1)
           clsmf25_struc(n)%cat_progn(t)%sndz(1)  = dz1max
        endif

        do i=N_snow,1,-1
           sdnew(i) = sdnew(i+1) + clsmf25_struc(n)%cat_progn(t)%sndz(i)
        enddo
        
        !**** Since the snow boundary has moved, redistribute heat
        !     contents & water equivalents of old to new snow layers.
        
        do i=1,N_snow
           
           j = i
           dz=sdnew(i+1)-sdold(i+1)
           if(dz < 0.) j = i + 1
           s(i+1,:) = h(j,:)*dz
           ds(i,:)  = s(i,:) - s(i+1,:)
        enddo
        
        clsmf25_struc(n)%cat_progn(t)%htsn = clsmf25_struc(n)%cat_progn(t)%htsn + ds(:,1)
        clsmf25_struc(n)%cat_progn(t)%wesn  = clsmf25_struc(n)%cat_progn(t)%wesn  + ds(:,2)
        
     !**** End relayer
   
     else
        
        clsmf25_struc(n)%cat_progn(t)%catdef = clsmf25_struc(n)%cat_progn(t)%catdef - &
             ((wesn1(t)- clsmf25_struc(n)%cat_progn(t)%wesn(1)) &
             + (wesn2(t) - clsmf25_struc(n)%cat_progn(t)%wesn(2)) &
             + (wesn3(t) - clsmf25_struc(n)%cat_progn(t)%wesn(3)))
     end if
!#endif
   ! add check for negative catdef
     if(clsmf25_struc(n)%cat_progn(t)%catdef<1.) then
!        clsmf25_struc(n)%cat_progn(t)%catdef=1. !max(1.,min(cdcr2(t),catdef(t)))
!  ignore the update
        clsmf25_struc(n)%cat_progn(t)%catdef = catdef_minus(t)
        clsmf25_struc(n)%cat_progn(t)%rzexc = rzexc_minus(t)
        clsmf25_struc(n)%cat_progn(t)%srfexc = srfexc_minus(t)
        clsmf25_struc(n)%cat_progn(t)%wesn(:) = wesn_minus(t,:)
        clsmf25_struc(n)%cat_progn(t)%sndz(:) = sndz_minus(t,:)
        clsmf25_struc(n)%cat_progn(t)%htsn(:) = htsn_minus(t,:)
     end if

!     if(LIS_localPet.eq.24.and.t.eq.7311) then 
!        print*, 'set: ',clsmf25_struc(n)%cat_progn(t)%htsn(1), &
!             clsmf25_struc(n)%cat_progn(t)%wesn(1), &
!             clsmf25_struc(n)%cat_diagn(t)%tpsn(1)
!     endif

  enddo
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     upd_check = .true. 
     if(clsmf25_struc(n)%cat_progn(t)%sndz(1).lt.0.or.&
          clsmf25_struc(n)%cat_progn(t)%sndz(2).lt.0.or.&
          clsmf25_struc(n)%cat_progn(t)%sndz(3).lt.0.or.&
          clsmf25_struc(n)%cat_progn(t)%wesn(1).lt.0.or.&
          clsmf25_struc(n)%cat_progn(t)%wesn(2).lt.0.or.&
          clsmf25_struc(n)%cat_progn(t)%wesn(3).lt.0) then 
        upd_check = .false. 
     endif
! update leads to unphysical snow values. roll back the update
     if(.not. upd_check) then 
        clsmf25_struc(n)%cat_progn(t)%catdef = catdef_minus(t)
        clsmf25_struc(n)%cat_progn(t)%rzexc = rzexc_minus(t)
        clsmf25_struc(n)%cat_progn(t)%srfexc = srfexc_minus(t)
        clsmf25_struc(n)%cat_progn(t)%wesn(:) = wesn_minus(t,:)
        clsmf25_struc(n)%cat_progn(t)%sndz(:) = sndz_minus(t,:)
        clsmf25_struc(n)%cat_progn(t)%htsn(:) = htsn_minus(t,:)
     endif
  enddo
  
  do k=1,N_snow
     call get_tf_nd( LIS_rc%npatch(n,LIS_rc%lsm_index), &
          clsmf25_struc(n)%cat_progn%htsn(k),&
          clsmf25_struc(n)%cat_progn%wesn(k), &
          clsmf25_struc(n)%cat_diagn%tpsn(k), fices(:,k) )
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     upd_check = .true. 
     do k=1,N_snow
        if((clsmf25_struc(n)%cat_diagn(t)%tpsn(k)+TF).lt.240.0) then 
           upd_check = .false. 
        endif
     enddo
!if temperatures are unphysical, roll back the update
     if(.not. upd_check) then 
        clsmf25_struc(n)%cat_progn(t)%catdef = catdef_minus(t)
        clsmf25_struc(n)%cat_progn(t)%rzexc = rzexc_minus(t)
        clsmf25_struc(n)%cat_progn(t)%srfexc = srfexc_minus(t)
        clsmf25_struc(n)%cat_progn(t)%wesn(:) = wesn_minus(t,:)
        clsmf25_struc(n)%cat_progn(t)%sndz(:) = sndz_minus(t,:)
        clsmf25_struc(n)%cat_progn(t)%htsn(:) = htsn_minus(t,:)
     endif

  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     swe_pchange = 0.0
     sndz_pchange = 0.0
     
     if(sum(wesn_minus(t,:)).gt.0) then 
        swe_pchange = abs((sum(clsmf25_struc(n)%cat_progn(t)%wesn(:)) - & 
             sum(wesn_minus(t,:)))*100.0/sum(wesn_minus(t,:)))
     endif
     if(sum(sndz_minus(t,:)).gt.0) then 
        sndz_pchange = abs((sum(clsmf25_struc(n)%cat_progn(t)%sndz(:)) - & 
             sum(sndz_minus(t,:)))*100.0/sum(sndz_minus(t,:)))
     endif
     
!     if(LIS_localPet.eq.28.and.t.eq.1303) then 
!        write(LIS_logunit,*) 'change ',swe_pchange, sndz_pchange
!        write(LIS_logunit,*) 'swe_ ',wesn_minus(t,:)
!        write(LIS_logunit,*) 'swe+ ',clsmf25_struc(n)%cat_progn(t)%wesn
!        write(LIS_logunit,*) 'sndz_ ',sndz_minus(t,:)
!        write(LIS_logunit,*) 'sndz+ ',clsmf25_struc(n)%cat_progn(t)%sndz
!
!     endif
     if((swe_pchange.gt.60.or.sndz_pchange.gt.60).and.&
          (sum(clsmf25_struc(n)%cat_progn(t)%wesn).gt.4*wemin)) then 
        if(sndz_pchange.eq.0) then
           clsmf25_struc(n)%cat_progn(t)%catdef = catdef_minus(t)
           clsmf25_struc(n)%cat_progn(t)%rzexc = rzexc_minus(t)
           clsmf25_struc(n)%cat_progn(t)%srfexc = srfexc_minus(t) 
           clsmf25_struc(n)%cat_progn(t)%wesn(:) = wesn_minus(t,:)
           clsmf25_struc(n)%cat_progn(t)%sndz(:) = sndz_minus(t,:)
           clsmf25_struc(n)%cat_progn(t)%htsn(:) = htsn_minus(t,:)
        elseif((swe_pchange/sndz_pchange.lt.0.8) .or. &
             (swe_pchange/sndz_pchange.gt.1.2) ) then 
           clsmf25_struc(n)%cat_progn(t)%catdef = catdef_minus(t)
           clsmf25_struc(n)%cat_progn(t)%rzexc = rzexc_minus(t)
           clsmf25_struc(n)%cat_progn(t)%srfexc = srfexc_minus(t)    
           clsmf25_struc(n)%cat_progn(t)%wesn(:) = wesn_minus(t,:)
           clsmf25_struc(n)%cat_progn(t)%sndz(:) = sndz_minus(t,:)
           clsmf25_struc(n)%cat_progn(t)%htsn(:) = htsn_minus(t,:)
           
        endif
     endif
     if(swe_pchange.gt.60.or.sndz_pchange.gt.60) then 
        clsmf25_struc(n)%cat_progn(t)%catdef = catdef_minus(t)
        clsmf25_struc(n)%cat_progn(t)%rzexc = rzexc_minus(t)
        clsmf25_struc(n)%cat_progn(t)%srfexc = srfexc_minus(t)
        clsmf25_struc(n)%cat_progn(t)%wesn(:) = wesn_minus(t,:)
        clsmf25_struc(n)%cat_progn(t)%sndz(:) = sndz_minus(t,:)
        clsmf25_struc(n)%cat_progn(t)%htsn(:) = htsn_minus(t,:)
     endif
  enddo
  CALL CALC_SOIL_MOIST (                       &
       LIS_rc%npatch(n,LIS_rc%lsm_index),      &
       clsmf25_struc(n)%cat_param%vegcls,  &
       clsmf25_struc(n)%cat_param%dzsf,    &
       clsmf25_struc(n)%cat_param%vgwmax,  &
       clsmf25_struc(n)%cat_param%cdcr1,   &
       clsmf25_struc(n)%cat_param%cdcr2,   &
       clsmf25_struc(n)%cat_param%wpwet,   &
       clsmf25_struc(n)%cat_param%poros,   &
       clsmf25_struc(n)%cat_param%psis,    &
       clsmf25_struc(n)%cat_param%bee,     &
       clsmf25_struc(n)%cat_param%ars1,    & 
       clsmf25_struc(n)%cat_param%ars2,    & 
       clsmf25_struc(n)%cat_param%ars3,    & 
       clsmf25_struc(n)%cat_param%ara1,    & 
       clsmf25_struc(n)%cat_param%ara2,    & 
       clsmf25_struc(n)%cat_param%ara3,    & 
       clsmf25_struc(n)%cat_param%ara4,    & 
       clsmf25_struc(n)%cat_param%arw1,    & 
       clsmf25_struc(n)%cat_param%arw2,    & 
       clsmf25_struc(n)%cat_param%arw3,    & 
       clsmf25_struc(n)%cat_param%arw4,    &
       clsmf25_struc(n)%cat_progn%srfexc,  &
       clsmf25_struc(n)%cat_progn%rzexc,   &
       clsmf25_struc(n)%cat_progn%catdef,  &
       clsmf25_struc(n)%cat_diagn%sfmc,    & 
       clsmf25_struc(n)%cat_diagn%rzmc,    & 
       clsmf25_struc(n)%cat_diagn%prmc)

!  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     if(LIS_localPet.eq.147.and.t.eq.1882) then 
!        print*, 'in set a ', clsmf25_struc(n)%cat_progn(t)%catdef,&
!             clsmf25_struc(n)%cat_progn(t)%rzexc, &
!             clsmf25_struc(n)%cat_progn(t)%srfexc, &
!             clsmf25_struc(n)%cat_progn(t)%sndz, &
!             clsmf25_struc(n)%cat_progn(t)%wesn
!     endif
!  enddo

!  if(LIS_localPet.eq.87) then 
!     write(LIS_logunit,*) 'snow ',clsmf25_struc(n)%cat_progn(22527)%sndz(:),&
!          clsmf25_struc(n)%cat_progn(22527)%wesn(:)
!  endif
  
end subroutine clsmf25_setsnowvars

