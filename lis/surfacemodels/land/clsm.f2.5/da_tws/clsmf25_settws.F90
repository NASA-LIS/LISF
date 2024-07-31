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
! !ROUTINE: clsmf25_settws
!  \label{clsmf25_settws}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 29Sep2011: Ben Zaitchik: Applied to GRACE
!
! !INTERFACE:
subroutine clsmf25_settws(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only  : LIS_logunit,LIS_verify
  use clsmf25_lsmMod
  
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
  real, parameter        :: WEMIN  = 13.0   ! [KG/M2]
  type(ESMF_Field)       :: catdefField
  type(ESMF_Field)       :: rzexcField
  type(ESMF_Field)       :: srfexcField
  type(ESMF_Field)       :: wesn1Field
  type(ESMF_Field)       :: wesn2Field
  type(ESMF_Field)       :: wesn3Field
  integer                :: t
  integer                :: status
  real, pointer          :: catdef(:), cdcr2(:)
  real, pointer          :: srfexc(:), rzexc(:)
  real, pointer          :: wesn1(:)
  real, pointer          :: wesn2(:)
  real, pointer          :: wesn3(:)
  real                   :: delta1,delta2,delta3,delta4,delta5,dz
  integer                :: i,j
  real, dimension(3,2)   :: ds
  real, dimension(4)     :: sdold, sdnew
  real, dimension(4,2)   :: h,s
  real                   :: snow_density,areasc,areasc1  
  integer, parameter ::  N_snow = 3
  real, parameter :: dz1max = 0.05   ! [m]
  real, parameter :: small  = 1.e-20
  real            :: swediff
  real                   :: check_val 

  call ESMF_StateGet(LSM_State,"Catchment Deficit",catdefField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Root Zone Excess",rzexcField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Surface Excess",srfexcField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 1",wesn1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 2",wesn2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 3",wesn3Field,rc=status)
  call LIS_verify(status) 

  call ESMF_FieldGet(catdefField,localDE=0,farrayPtr=catdef,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(rzexcField,localDE=0,farrayPtr=rzexc,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(srfexcField,localDE=0,farrayPtr=srfexc,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn1Field,localDE=0,farrayPtr=wesn1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn2Field,localDE=0,farrayPtr=wesn2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn3Field,localDE=0,farrayPtr=wesn3,rc=status)
  call LIS_verify(status)  

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     
     clsmf25_struc(n)%cat_progn(t)%catdef = catdef(t)
     clsmf25_struc(n)%cat_progn(t)%srfexc = srfexc(t)
     if(rzexc(t) > 0.001) then
        clsmf25_struc(n)%cat_progn(t)%rzexc = rzexc(t)
     else
        clsmf25_struc(n)%cat_progn(t)%catdef = catdef(t) - &
             (rzexc(t)-clsmf25_struc(n)%cat_progn(t)%rzexc)
     end if

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
          if (areasc .eq. 1.0) then
             clsmf25_struc(n)%cat_progn(t)%sndz(1) =  &
                clsmf25_struc(n)%cat_progn(t)%wesn(1)/&
                (snow_density*areasc1)
           elseif (areasc .lt. 1.0 .and. areasc1.eq.1.0) then
             clsmf25_struc(n)%cat_progn(t)%sndz(1) =  &
                clsmf25_struc(n)%cat_progn(t)%sndz(1)*wesn1(t)/wemin
           elseif (areasc .eq. 1.0 .and. areasc1.lt.1.0 ) then
             clsmf25_struc(n)%cat_progn(t)%sndz(1) =  &
                clsmf25_struc(n)%cat_progn(t)%sndz(1)*wemin/&
                clsmf25_struc(n)%cat_progn(t)%wesn(1)
           end if

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
           if (areasc .eq. 1.0) then
             clsmf25_struc(n)%cat_progn(t)%sndz(2) =  &
                clsmf25_struc(n)%cat_progn(t)%wesn(2)/&
                (snow_density*areasc1)
           elseif (areasc .lt. 1.0 .and. areasc1.eq.1.0) then
             clsmf25_struc(n)%cat_progn(t)%sndz(2) =  &
                clsmf25_struc(n)%cat_progn(t)%sndz(2)*wesn2(t)/wemin
           elseif (areasc .eq. 1.0 .and. areasc1.lt.1.0 ) then
             clsmf25_struc(n)%cat_progn(t)%sndz(2) =  &
                clsmf25_struc(n)%cat_progn(t)%sndz(2)*wemin/&
                clsmf25_struc(n)%cat_progn(t)%wesn(2)
           end if

           clsmf25_struc(n)%cat_progn(t)%wesn(2) = wesn2(t)

!           clsmf25_struc(n)%cat_progn(t)%sndz(2) = clsmf25_struc(n)%cat_progn(t)%sndz(2) &
!                * wesn2(t)/clsmf25_struc(n)%cat_progn(t)%wesn(2)       
        else
           ! write(LIS_logunit,*)'dried to wesn2 low limit, tile',t,'wesn2 is',wesn2(t)
           clsmf25_struc(n)%cat_progn(t)%catdef &
                = clsmf25_struc(n)%cat_progn(t)%catdef - (wesn2(t)-clsmf25_struc(n)%cat_progn(t)%wesn(2))
        end if
        
        if(wesn3(t).gt.2.0) then
           snow_density = clsmf25_struc(n)%cat_progn(t)%wesn(3)/&
                (clsmf25_struc(n)%cat_progn(t)%sndz(3)*areasc)
           clsmf25_struc(n)%cat_progn(t)%htsn(3) = &
                clsmf25_struc(n)%cat_progn(t)%htsn(3) &
                * wesn3(t)/clsmf25_struc(n)%cat_progn(t)%wesn(3)
           if (areasc .eq. 1.0) then
             clsmf25_struc(n)%cat_progn(t)%sndz(3) =  &
                clsmf25_struc(n)%cat_progn(t)%wesn(3)/&
                (snow_density*areasc1)
           elseif (areasc .lt. 1.0 .and. areasc1.eq.1.0) then
             clsmf25_struc(n)%cat_progn(t)%sndz(3) =  &
                clsmf25_struc(n)%cat_progn(t)%sndz(3)*wesn3(t)/wemin
           elseif (areasc .eq. 1.0 .and. areasc1.lt.1.0 ) then
             clsmf25_struc(n)%cat_progn(t)%sndz(3) =  &
                clsmf25_struc(n)%cat_progn(t)%sndz(3)*wemin/&
                clsmf25_struc(n)%cat_progn(t)%wesn(3)
           end if

           clsmf25_struc(n)%cat_progn(t)%wesn(3) = wesn3(t)

!           clsmf25_struc(n)%cat_progn(t)%sndz(3) = clsmf25_struc(n)%cat_progn(t)%sndz(3) &
!                * wesn3(t)/clsmf25_struc(n)%cat_progn(t)%wesn(3)        
           ! write(LIS_logunit,*)'updated snow3 in tile',t
        else
           ! write(LIS_logunit,*)'dried to wesn3 low limit, tile',t,'wesn3 is',wesn3(t)
           clsmf25_struc(n)%cat_progn(t)%catdef &
                = clsmf25_struc(n)%cat_progn(t)%catdef - (wesn3(t) - clsmf25_struc(n)%cat_progn(t)%wesn(3))
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
              h(i,1) = clsmf25_struc(n)%cat_progn(t)%htsn(i)/clsmf25_struc(n)%cat_progn(t)%sndz(i)
              h(i,2) =  clsmf25_struc(n)%cat_progn(t)%wesn(i)/clsmf25_struc(n)%cat_progn(t)%sndz(i)
           endif
        enddo

        !**** Obtain old & new layer thicknesses & boundaries.
        
        sdold = 0.
        sdnew = 0.
        
        do i=N_snow,1,-1
           sdold(i) = sdold(i+1) + clsmf25_struc(n)%cat_progn(t)%sndz(i)
        enddo
        
!        clsmf25_struc(n)%cat_progn(t)%sndz = sdold(1)/float(N_snow)
!        if(clsmf25_struc(n)%cat_progn(t)%sndz(1) > dz1max) then
!           clsmf25_struc(n)%cat_progn(t)%sndz(2:) = (sdold(1)-dz1max)/float(N_snow-1)
!           clsmf25_struc(n)%cat_progn(t)%sndz(1)  = dz1max
!        endif

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
        write(LIS_logunit,*) 't',t,'catdef(t)',catdef(t)
        clsmf25_struc(n)%cat_progn(t)%catdef=1. !max(1.,min(cdcr2(t),catdef(t)))
     end if

!     if(LIS_localPet.eq.24.and.t.eq.7311) then 
!        print*, 'set: ',clsmf25_struc(n)%cat_progn(t)%htsn(1), &
!             clsmf25_struc(n)%cat_progn(t)%wesn(1), &
!             clsmf25_struc(n)%cat_diagn(t)%tpsn(1)
!     endif
  enddo
  
end subroutine clsmf25_settws

