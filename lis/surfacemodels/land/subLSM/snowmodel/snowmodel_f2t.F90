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
! !ROUTINE: snowmodel_f2t
! \label{snowmodel_f2t}
!
! !REVISION HISTORY:
!  15 Oct 2002: Sujay Kumar; Initial Code
!  14 Apr 2020: Kristi Arsenault; Added G. Liston's SnowModel
!
! !INTERFACE:
subroutine snowmodel_f2t(n)
! !USES:      
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_surface
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_logMod,        only : LIS_logunit, LIS_verify, &
                                LIS_endrun
  use snowmodel_lsmMod
  use snowmodel_vars 

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the SnowModel
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
  type(ESMF_Field)   :: psurfField,pcpField,snowfField,fhgtField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:)
  real,pointer       :: fheight1(:,:),fheight(:)
  integer, pointer   :: layer_windht(:), layer_relhumht(:)
! __________________

  write(LIS_logunit,*) '[INFO] Call to the SnowModel Forcing 2 tile routine ...'

  if ( LIS_FORC_Forc_Hgt%selectOpt.eq.1 ) then

     allocate(fheight1(LIS_FORC_Forc_Hgt%vlevels,LIS_rc%ntiles(n)))
     allocate(layer_windht(LIS_rc%ntiles(n)))
     allocate(layer_relhumht(LIS_rc%ntiles(n)))

   ! Forcing heights are specified. Find the layers corresponding to the reference heights.
   ! If not, use the lowest model layer as the height.

     do v=1,LIS_FORC_Forc_Hgt%vlevels
        call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Forc_Hgt%varname(v)),&
             fhgtField, rc=status)
        call LIS_verify(status,'snowmodel_f2t: error getting Forc_Hgt')

        call ESMF_FieldGet(fhgtField,localDE=0,farrayPtr=fheight,rc=status)
        fheight1(v,:) = fheight(:)
        call LIS_verify(status,'snowmodel_f2t: error retrieving fheight')
     enddo

     call snowmodel_find_forcing_heights(LIS_FORC_Forc_Hgt%vlevels, LIS_rc%ntiles(n),&
              fheight1, snowmodel_struc(n)%ht_rhobs, snowmodel_struc(n)%ht_windobs, &
              layer_windht, layer_relhumht)

     if (LIS_rc%ntiles(n).ne.0) then
        snowmodel_struc(n)%ht_rhobs = fheight1(layer_relhumht(1),1)
        snowmodel_struc(n)%ht_windobs = fheight1(layer_windht(1),1)
     endif

     deallocate(fheight1)
     deallocate(layer_windht)
     deallocate(layer_relhumht)
  endif
  ! If forcing heights are not specified, then LIS will assume that forcing data
  ! corresponds to the reference heights snowmodel.par (tmprh_ht, wind_ht).

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status,'snowmodel_f2t: error getting Tair')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status,'snowmodel_f2t: error getting Qair')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_E%varname(1)),uField,&
       rc=status)
  call LIS_verify(status,'snowmodel_f2t: error getting Wind_E')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_N%varname(1)),vField,&
       rc=status)
  call LIS_verify(status,'snowmodel_f2t: error getting Wind_N')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status,'snowmodel_f2t: error getting SWdown')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status,'snowmodel_f2t: error getting LWdown')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status, 'snowmodel_f2t: error getting PSurf')

  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status,'snowmodel_f2t: error getting Rainf')

  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Forc_Hgt%varname(1)),fhgtField,&
          rc=status)
     call LIS_verify(status,'snowmodel_f2t: error getting Forc_Hgt')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Snowf%varname(1)),snowfField,&
          rc=status)
     call LIS_verify(status,'snowmodel_f2t: error getting Snowf')
  endif

  call ESMF_FieldGet(tmpField, localDE=0, farrayPtr=tmp,rc=status)
  call LIS_verify(status, 'snowmodel_f2t: error retrieving Tair')

  call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
  call LIS_verify(status,'snowmodel_f2t: error retrieving q2')

  call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
  call LIS_verify(status,'snowmodel_f2t: error retrieving swd')

  call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
  call LIS_verify(status,'snowmodel_f2t: error retrieving lwd')

  call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
  call LIS_verify(status,'snowmodel_f2t: error retrieving u')

  call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
  call LIS_verify(status,'snowmodel_f2t: error retrieving v')

  call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
  call LIS_verify(status,'snowmodel_f2t: error retrieving psurf')

  call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
  call LIS_verify(status,'snowmodel_f2t: error retrieving pcp')

  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
     call ESMF_FieldGet(fhgtField,localDE=0, farrayPtr=fheight,rc=status)
     call LIS_verify(status,'snowmodel_f2t: error retrieving forc_hgt')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_FieldGet(snowfField,localDE=0, farrayPtr=snowf,rc=status)
     call LIS_verify(status,'snowmodel_f2t: error retrieving snowf')
  endif

  snowmodel_struc(n)%forc_count = snowmodel_struc(n)%forc_count + 1

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ! Transform tile to the patch
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
     snowmodel_struc(n)%sm(t)%tair=snowmodel_struc(n)%sm(t)%tair + tmp(tid)
     snowmodel_struc(n)%sm(t)%qair=snowmodel_struc(n)%sm(t)%qair + q2(tid)
     snowmodel_struc(n)%sm(t)%swdown=snowmodel_struc(n)%sm(t)%swdown + swd(tid)
     snowmodel_struc(n)%sm(t)%lwdown=snowmodel_struc(n)%sm(t)%lwdown + lwd(tid)
     snowmodel_struc(n)%sm(t)%uwind=snowmodel_struc(n)%sm(t)%uwind + uwind(tid)
     snowmodel_struc(n)%sm(t)%vwind=snowmodel_struc(n)%sm(t)%vwind + vwind(tid)
     snowmodel_struc(n)%sm(t)%psurf=snowmodel_struc(n)%sm(t)%psurf + psurf(tid)

     if ( LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
        snowmodel_struc(n)%sm(t)%fheight=fheight(tid)
     endif

     if(pcp(tid).ne.LIS_rc%udef) then
        snowmodel_struc(n)%sm(t)%rainf=snowmodel_struc(n)%sm(t)%rainf + pcp(tid)
     else
        snowmodel_struc(n)%sm(t)%rainf=snowmodel_struc(n)%sm(t)%rainf + 0.0
     endif

     ! If there is snowf add it to precipitation. 
     ! NOTE: SnowModel has options/code to discriminate between rainf/snowf
     !        within MicroMet.
     !   Will need to address this issue ...
     !   Place options / code in new "sliced out" MicroMet routines ...
     if ( LIS_FORC_Snowf%selectOpt.eq.1) then
        if(snowf(tid).ne.LIS_rc%udef) then
           snowmodel_struc(n)%sm(t)%rainf=snowmodel_struc(n)%sm(t)%rainf + &
                snowf(tid)
        endif
     endif
     snowmodel_struc(n)%sm(t)%snowf = snowmodel_struc(n)%sm(t)%snowf+ 0.0

  enddo

!     if(pcp(100) < 0) then
!   print *, "pcp :: ",snowmodel_struc(n)%forc_count, pcp(100)
!     endif
!  print *, tmp(100), snowmodel_struc(n)%sm(100)%tair
!  print *, uwind(100), snowmodel_struc(n)%sm(100)%uwind
!  print *, vwind(100), snowmodel_struc(n)%sm(100)%vwind
!  print *, psurf(100), snowmodel_struc(n)%sm(100)%psurf
!  print *, pcp(100), snowmodel_struc(n)%sm(100)%rainf

end subroutine snowmodel_f2t


subroutine snowmodel_find_forcing_heights(vlevels, ntiles, fheight, &
     tmprh_ht, wind_ht, layer_tmprh, layer_wind)

  implicit none

  integer      :: vlevels, ntiles
  real         :: fheight(vlevels, ntiles)
  real         :: tmprh_ht, wind_ht
  integer      :: layer_tmprh(ntiles)
  integer      :: layer_wind(ntiles)
  integer      :: t,v

!assume lowest model layer to start with.
  layer_tmprh = 1
  layer_wind = 1

  do t=1,ntiles
     do v=1,vlevels
        if(fheight(v,t).ge.tmprh_ht) layer_tmprh(t) = v
        if(fheight(v,t).ge.wind_ht) layer_wind(t) = v
     enddo
  enddo

end subroutine snowmodel_find_forcing_heights
