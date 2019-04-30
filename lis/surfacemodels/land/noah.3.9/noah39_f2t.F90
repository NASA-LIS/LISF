!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noah39_f2t
! \label{noah39_f2t}
!
! !REVISION HISTORY:
!  15 Oct 2002: Sujay Kumar; Initial Code
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
!
! !INTERFACE:
subroutine noah39_f2t(n)
! !USES:      
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_surface
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_logMod,         only : LIS_verify
  use noah39_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the Noah-3.9
!  model tiles. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

  integer            :: t,v,status
  integer            :: tid
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField,fhgtField,acondField
  type(ESMF_Field)   :: gvfField, albField, z0Field
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer       :: fheight1(:,:),fheight(:),acond(:)
  real,pointer       :: gvf(:),alb(:),z0(:)
  integer, pointer   :: layer_h(:), layer_m(:)

  if ( LIS_FORC_Forc_Hgt%selectOpt.eq.1 ) then
     allocate(fheight1(LIS_FORC_Forc_Hgt%vlevels,LIS_rc%ntiles(n)))
! forcing heights are specified. Find the layers corresponding to the reference heights.
! if not, use the lowest model layer as the height.
     allocate(layer_h(LIS_rc%ntiles(n)))
     allocate(layer_m(LIS_rc%ntiles(n)))

     do v=1,LIS_FORC_Forc_Hgt%vlevels
        call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Forc_Hgt%varname(v)),&
             fhgtField, rc=status)
        call LIS_verify(status,'noah39_f2t: error getting Forc_Hgt')
        call ESMF_FieldGet(fhgtField,localDE=0,farrayPtr=fheight,rc=status)
        fheight1(v,:) = fheight(:)
        call LIS_verify(status,'noah39_f2t: error retrieving fheight')
     enddo

     call noah39_find_forcing_heights(LIS_FORC_Forc_Hgt%vlevels, LIS_rc%ntiles(n),&
          fheight1, noah39_struc(n)%zh, noah39_struc(n)%zm, layer_h, layer_m)

     if (LIS_rc%ntiles(n).ne.0) then
        noah39_struc(n)%zh = fheight1(layer_h(1),1)
        noah39_struc(n)%zm = fheight1(layer_m(1),1)
     endif

     deallocate(fheight1)
     deallocate(layer_h)
     deallocate(layer_m)
  endif
! if forcing heights are not specified, then LIS will assume that forcing data
! corresponds to the reference heights specified in the config file (zh, zm)

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status,'noah39_f2t: error getting Tair')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status,'noah39_f2t: error getting Qair')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_E%varname(1)),uField,&
       rc=status)
  call LIS_verify(status,'noah39_f2t: error getting Wind_E')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_N%varname(1)),vField,&
       rc=status)
  call LIS_verify(status,'noah39_f2t: error getting Wind_N')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status,'noah39_f2t: error getting SWdown')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status,'noah39_f2t: error getting LWdown')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status, 'noah39_f2t: error getting PSurf')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status,'noah39_f2t: error getting Rainf')

  if(LIS_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_CRainf%varname(1)),cpcpField,&
          rc=status)
     call LIS_verify(status,'noah39_f2t: error getting CRainf')
  endif

  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Forc_Hgt%varname(1)),fhgtField,&
          rc=status)
     call LIS_verify(status,'noah39_f2t: error getting Forc_Hgt')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Snowf%varname(1)),snowfField,&
          rc=status)
     call LIS_verify(status,'noah39_f2t: error getting Snowf')
  endif

  if ( LIS_FORC_Ch%selectOpt.eq.1 ) then
     noah39_struc(n)%forcing_ch = 1
     call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Ch%varname(1)),&
          acondField,rc=status)
     call LIS_verify(status,'noah39_f2t: error getting ch')
  else
     noah39_struc(n)%forcing_ch = 0
  endif

  if ( LIS_FORC_Cm%selectOpt.eq.1 ) then
     noah39_struc(n)%forcing_cm = 1
  else
     noah39_struc(n)%forcing_cm = 0
  endif

  if ( LIS_FORC_GVF%selectOpt.eq.1 ) then
     call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_GVF%varname(1)),&
          gvfField,rc=status)
     call LIS_verify(status,'noah39_f2t: error getting GVF')
  endif
  if ( LIS_FORC_Alb%selectOpt.eq.1 ) then
     call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Alb%varname(1)),&
          albField,rc=status)
     call LIS_verify(status,'noah39_f2t: error getting Alb')
  endif
  if ( LIS_FORC_Z0%selectOpt.eq.1 ) then
     call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Z0%varname(1)),&
          z0Field,rc=status)
     call LIS_verify(status,'noah39_f2t: error getting Z0')
  endif

  call ESMF_FieldGet(tmpField, localDE=0, farrayPtr= tmp,rc=status)
  call LIS_verify(status, 'noah39_f2t: error retrieving tmp')

  call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
  call LIS_verify(status,'noah39_f2t: error retrieving q2')

  call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
  call LIS_verify(status,'noah39_f2t: error retrieving swd')

  call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
  call LIS_verify(status,'noah39_f2t: error retrieving lwd')

  call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
  call LIS_verify(status,'noah39_f2t: error retrieving u')

  call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
  call LIS_verify(status,'noah39_f2t: error retrieving v')

  call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
  call LIS_verify(status,'noah39_f2t: error retrieving psurf')

  call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
  call LIS_verify(status,'noah39_f2t: error retrieving pcp')

  if(LIS_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_FieldGet(cpcpField,localDE=0, farrayPtr=cpcp,rc=status)
     call LIS_verify(status, 'noah39_f2t: error retrieving cpcp')
  endif

  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
     call ESMF_FieldGet(fhgtField,localDE=0, farrayPtr=fheight,rc=status)
     call LIS_verify(status,'noah39_f2t: error retrieving forc_hgt')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_FieldGet(snowfField,localDE=0, farrayPtr=snowf,rc=status)
     call LIS_verify(status,'noah39_f2t: error retrieving snowf')
  endif

  if ( LIS_FORC_Ch%selectOpt.eq.1 ) then
     call ESMF_FieldGet(acondField,localDE=0, farrayPtr=acond,rc=status)
     call LIS_verify(status,'noah39_f2t: error retrieving acond')
  endif

 if ( LIS_FORC_GVF%selectOpt.eq.1 ) then
     call ESMF_FieldGet(gvfField,localDE=0, farrayPtr=gvf,rc=status)
     call LIS_verify(status,'noah39_f2t: error retrieving gvf')
  endif
  if ( LIS_FORC_Alb%selectOpt.eq.1 ) then
     call ESMF_FieldGet(albField,localDE=0, farrayPtr=alb,rc=status)
     call LIS_verify(status,'noah39_f2t: error retrieving alb')
  endif
  if ( LIS_FORC_Z0%selectOpt.eq.1 ) then
     call ESMF_FieldGet(z0Field,localDE=0, farrayPtr=z0,rc=status)
     call LIS_verify(status,'noah39_f2t: error retrieving z0')
  endif

  noah39_struc(n)%forc_count  =   noah39_struc(n)%forc_count  + 1

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     !transform t to the patch
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
     noah39_struc(n)%noah(t)%tair=noah39_struc(n)%noah(t)%tair+ tmp(tid)
     noah39_struc(n)%noah(t)%qair=noah39_struc(n)%noah(t)%qair + q2(tid)
     noah39_struc(n)%noah(t)%swdown=noah39_struc(n)%noah(t)%swdown + swd(tid)
     noah39_struc(n)%noah(t)%lwdown=noah39_struc(n)%noah(t)%lwdown + lwd(tid)
     noah39_struc(n)%noah(t)%uwind=noah39_struc(n)%noah(t)%uwind + uwind(tid)
     noah39_struc(n)%noah(t)%vwind=noah39_struc(n)%noah(t)%vwind + vwind(tid)
     noah39_struc(n)%noah(t)%psurf=noah39_struc(n)%noah(t)%psurf + psurf(tid)
     if ( LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
        noah39_struc(n)%noah(t)%fheight=fheight(tid)
     endif

     ! Handle cases when input data does
     ! not have enough precision for this
     ! field
     if ( LIS_FORC_Ch%selectOpt.eq.1) then
        if(acond(tid).le.0) then
           noah39_struc(n)%noah(t)%ch=noah39_struc(n)%noah(t)%ch + 0.00001
        else
           noah39_struc(n)%noah(t)%ch=noah39_struc(n)%noah(t)%ch + acond(tid)
        endif
     endif

     if(pcp(tid).ne.LIS_rc%udef) then
        noah39_struc(n)%noah(t)%rainf=noah39_struc(n)%noah(t)%rainf + pcp(tid)
     else
        noah39_struc(n)%noah(t)%rainf=noah39_struc(n)%noah(t)%rainf + 0.0
     endif
     ! Not all forcings supply convective rainfall, but it is used in
     ! Noah main for diagnosing total precipitation output.
     ! Therefore, initialize rainf_c to 0.

     if ( LIS_FORC_CRainf%selectOpt.eq.1) then
        if(cpcp(tid).ne.LIS_rc%udef) then
           noah39_struc(n)%noah(t)%rainf_c=noah39_struc(n)%noah(t)%rainf_c+ cpcp(tid)
        else
           noah39_struc(n)%noah(t)%rainf_c=noah39_struc(n)%noah(t)%rainf_c + 0.0
        endif
     endif
! If there is snowf add it to precipitation.  Noah-3.6 does not use
! separate rainf and snowf.  It determines what to do with precipitation.
     if ( LIS_FORC_Snowf%selectOpt.eq.1) then
        if(snowf(tid).ne.LIS_rc%udef) then
           noah39_struc(n)%noah(t)%rainf=noah39_struc(n)%noah(t)%rainf + & 
                snowf(tid)
        endif
     endif
     noah39_struc(n)%noah(t)%snowf = noah39_struc(n)%noah(t)%snowf+ 0.0
     if(LIS_FORC_GVF%selectOpt.eq.1) then 
        noah39_struc(n)%noah(t)%shdfac=noah39_struc(n)%noah(t)%shdfac + & 
             gvf(tid)
     endif

     if(LIS_FORC_ALB%selectOpt.eq.1) then 
        noah39_struc(n)%noah(t)%alb=noah39_struc(n)%noah(t)%alb + & 
             alb(tid)
     endif
     
     if(LIS_FORC_Z0%selectOpt.eq.1) then 
        noah39_struc(n)%noah(t)%z0=noah39_struc(n)%noah(t)%z0 + & 
             z0(tid)
     endif
  enddo
end subroutine noah39_f2t

subroutine noah39_find_forcing_heights(vlevels, ntiles, fheight, zh, &
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

end subroutine noah39_find_forcing_heights
