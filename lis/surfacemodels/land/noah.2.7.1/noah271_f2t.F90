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
! !ROUTINE: noah271_f2t
! \label{noah271_f2t}
!
! !REVISION HISTORY:
!  15 Oct 2002: Sujay Kumar; Initial Code
!  27 Oct 2010: David Mocko, changes for Noah2.7.1 in LIS6.1
!
! !INTERFACE:
subroutine noah271_f2t(n)
! !USES:      
  use ESMF
  use LIS_coreMod
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_logMod,         only : LIS_verify
  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the Noah2.7.1
!  model tiles. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

  integer            :: t,tid,v,status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField,fhgtField,acondField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer       :: fheight1(:,:),fheight(:),acond(:)
  integer, allocatable   :: layer_h(:), layer_m(:)

  if ( LIS_FORC_Forc_Hgt%selectOpt.eq.1 ) then
     allocate(fheight1(LIS_FORC_Forc_Hgt%vlevels,LIS_rc%npatch(n,LIS_rc%lsm_index)))
! forcing heights are specified. Find the layers 
! corresponding to the reference heights.
! if not, use the lowest model layer as the height.
     allocate(layer_h(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(layer_m(LIS_rc%npatch(n,LIS_rc%lsm_index)))

     do v=1,LIS_FORC_Forc_Hgt%vlevels
        call ESMF_StateGet(LIS_FORC_State(n),&
             trim(LIS_FORC_Forc_Hgt%varname(v)),&
             fhgtField, rc=status)
        call LIS_verify(status,'noah271_f2t: error getting Forc_Hgt')
        call ESMF_FieldGet(fhgtField,localDE=0,farrayPtr=fheight,rc=status)
        fheight1(v,:) = fheight(:)
        call LIS_verify(status,'noah271_f2t: error retrieving fheight')
     enddo

     call noah271_find_forcing_heights(LIS_FORC_Forc_Hgt%vlevels, &
          LIS_rc%npatch(n,LIS_rc%lsm_index),&
          fheight1, noah271_struc(n)%zh, noah271_struc(n)%zm,&
          layer_h, layer_m)

     if (LIS_rc%ntiles(n).ne.0) then
        noah271_struc(n)%zh = fheight1(layer_h(1),1)
        noah271_struc(n)%zm = fheight1(layer_m(1),1)
     endif

     deallocate(fheight1)
     deallocate(layer_h)
     deallocate(layer_m)
  endif
! if forcing heights are not specified, then LIS will assume that forcing data
! corresponds to the reference heights specified in the config file (zh, zm)

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status,'noah271_f2t: error getting Tair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status,'noah271_f2t: error getting Qair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),uField,&
       rc=status)
  call LIS_verify(status,'noah271_f2t: error getting Wind_E')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),vField,&
       rc=status)
  call LIS_verify(status,'noah271_f2t: error getting Wind_N')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status,'noah271_f2t: error getting SWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status,'noah271_f2t: error getting LWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status, 'noah271_f2t: error getting PSurf')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status,'noah271_f2t: error getting Rainf')

  if(LIS_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
          rc=status)
     call LIS_verify(status,'noah271_f2t: error getting CRainf')
  endif

  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Forc_Hgt%varname(1)),fhgtField,&
          rc=status)
     call LIS_verify(status,'noah271_f2t: error getting Forc_Hgt')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),snowfField,&
          rc=status)
     call LIS_verify(status,'noah271_f2t: error getting Snowf')
  endif

  if ( LIS_FORC_Ch%selectOpt.eq.1 ) then
     noah271_struc(n)%forcing_ch = 1
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Ch%varname(1)),&
          acondField,rc=status)
     call LIS_verify(status,'noah271_f2t: error getting ch')
  else
     noah271_struc(n)%forcing_ch = 0
  endif

  call ESMF_FieldGet(tmpField, localDE=0, farrayPtr= tmp,rc=status)
  call LIS_verify(status, 'noah271_f2t: error retrieving tmp')

  call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving q2')

  call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving swd')

  call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving lwd')

  call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving u')

  call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving v')

  call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving psurf')

  call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving pcp')

  if(LIS_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_FieldGet(cpcpField,localDE=0, farrayPtr=cpcp,rc=status)
     call LIS_verify(status, 'noah271_f2t: error retrieving cpcp')
  endif

  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
     call ESMF_FieldGet(fhgtField,localDE=0, farrayPtr=fheight,rc=status)
     call LIS_verify(status,'noah271_f2t: error retrieving forc_hgt')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_FieldGet(snowfField,localDE=0, farrayPtr=snowf,rc=status)
     call LIS_verify(status,'noah271_f2t: error retrieving snowf')
  endif

  if ( LIS_FORC_Ch%selectOpt.eq.1 ) then
     call ESMF_FieldGet(acondField,localDE=0, farrayPtr=acond,rc=status)
     call LIS_verify(status,'noah271_f2t: error retrieving acond')
  endif

  noah271_struc(n)%forc_count = noah271_struc(n)%forc_count + 1
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
     noah271_struc(n)%noah(t)%tair=noah271_struc(n)%noah(t)%tair + & 
          tmp(tid)
     noah271_struc(n)%noah(t)%qair=noah271_struc(n)%noah(t)%qair + &
          q2(tid)
     noah271_struc(n)%noah(t)%swdown=noah271_struc(n)%noah(t)%swdown + & 
          swd(tid)
     noah271_struc(n)%noah(t)%lwdown=noah271_struc(n)%noah(t)%lwdown + &
          lwd(tid)
     noah271_struc(n)%noah(t)%uwind=noah271_struc(n)%noah(t)%uwind + & 
          uwind(t)
     noah271_struc(n)%noah(t)%vwind=noah271_struc(n)%noah(t)%vwind + & 
          vwind(t)
     noah271_struc(n)%noah(t)%psurf=noah271_struc(n)%noah(t)%psurf + & 
          psurf(tid)
!no averaging (in time) of forcing height is done
     if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
        noah271_struc(n)%noah(t)%fheight=  fheight(tid)
     endif
     noah271_struc(n)%noah(t)%emiss=1.0

     ! Handle cases when input data does
     ! not have enough precision for this
     ! field
     if ( LIS_FORC_Ch%selectOpt.eq.1) then
        if(acond(t).le.0) then
           noah271_struc(n)%noah(t)%ch=noah271_struc(n)%noah(t)%ch + &
                0.00001
        else
           noah271_struc(n)%noah(t)%ch=noah271_struc(n)%noah(t)%ch + &
                acond(tid)
        endif
     !else
     !   noah271_struc(n)%noah(t)%ch=LIS_rc%udef
     endif

     if(pcp(t).ne.LIS_rc%udef) then
        noah271_struc(n)%noah(t)%rainf=noah271_struc(n)%noah(t)%rainf + pcp(tid)
     else
        noah271_struc(n)%noah(t)%rainf=noah271_struc(n)%noah(t)%rainf + 0.0
     endif
     ! Not all forcings supply convective rainfall, but it is used in
     ! Noah main for diagnosing total precipitation output.
     ! Therefore, initialize rainf_c to 0.

     if ( LIS_FORC_CRainf%selectOpt.eq.1) then
        if(cpcp(t).ne.LIS_rc%udef) then
           noah271_struc(n)%noah(t)%rainf_c=noah271_struc(n)%noah(t)%rainf_c +&
                cpcp(tid)
        else
           noah271_struc(n)%noah(t)%rainf_c=noah271_struc(n)%noah(t)%rainf_c + &
                0.0
        endif
     endif
! If there is snowf add it to precipitation.  Noah2.7.1 does not use
! separate rainf and snowf.  It determines what to do with precipitation.
     if ( LIS_FORC_Snowf%selectOpt.eq.1) then
        if(snowf(t).ne.LIS_rc%udef) then
           noah271_struc(n)%noah(t)%rainf=noah271_struc(n)%noah(t)%rainf + &
                snowf(tid)
        endif
     endif
     noah271_struc(n)%noah(t)%snowf=noah271_struc(n)%noah(t)%snowf+ 0.0
  enddo
end subroutine noah271_f2t

subroutine noah271_find_forcing_heights(vlevels, ntiles, fheight, zh, &
     zm, layer_h, layer_m)

  implicit none

  integer      :: vlevels, ntiles
  real         :: fheight(vlevels, ntiles)
  real         :: zh, zm
  integer      :: layer_h(ntiles)
  integer      :: layer_m(ntiles)
  integer      :: t,v

!assume lowest model layer to start with.
  layer_h = 1
  layer_m = 1

  do t=1,ntiles
     do v=1,vlevels
        if(fheight(v,t).ge.zh) layer_h(t) = v
        if(fheight(v,t).ge.zm) layer_m(t) = v
     enddo
  enddo

end subroutine noah271_find_forcing_heights
