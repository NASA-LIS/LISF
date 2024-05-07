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
! !ROUTINE: clsmf25_f2t
!  \label{clsmf25_f2t}
!
!
! !REVISION HISTORY:
! 16 Dec 2005: Sujay Kumar; Initial Specification
! 11 Oct 2006: Sujay Kumar; modified to work with ESMF_State design
! 23 Nov 2012: David Mocko, Added forcing height option
!
! !INTERFACE:
subroutine clsmf25_f2t(n)
! !USES:      
  use ESMF
  use LIS_coreMod , only : LIS_rc,LIS_surface
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_constantsMod,       only : LIS_CONST_TKFRZ
  use LIS_logMod,             only : LIS_logunit,LIS_verify,LIS_endrun
  use clsmf25_lsmMod, only : clsmf25_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 

! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the Catchment
!  model tiles. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP  
  integer            :: t,tid,status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField,RefHField
  type(ESMF_Field)   :: swnetField,pardrField,pardfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:),RefH(:)
  real,pointer       :: swnet(:),pardr(:),pardf(:)

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status,'clsmf25_f2t: error getting Tair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status,'clsmf25_f2t: error getting Qair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status,'clsmf25_f2t: error getting SWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status,'clsmf25_f2t: error getting LWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),uField,&
       rc=status)
  call LIS_verify(status,'clsmf25_f2t: error getting Wind_E')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),vField,&
       rc=status)
  call LIS_verify(status,'clsmf25_f2t: error getting Wind_N')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status, 'clsmf25_f2t: error getting PSurf')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status,'clsmf25_f2t: error getting Rainf')

  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
          rc=status)
     call LIS_verify(status,'clsmf25_f2t: error getting CRainf')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),snowfField,&
          rc=status)
     call LIS_verify(status,'clsmf25_f2t: error getting Snowf')
  endif

  if(LIS_FORC_SWnet%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWnet%varname(1)),&
          swnetField,rc=status)
     call LIS_verify(status,'clsmf25_f2t: error getting SWnet')
  endif

  if(LIS_FORC_Pardr%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Pardr%varname(1)),&
          pardrField,rc=status)
     call LIS_verify(status,'clsmf25_f2t: error getting Pardr')
  endif

  if(LIS_FORC_Pardf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Pardf%varname(1)),&
          pardfField,rc=status)
     call LIS_verify(status,'clsmf25_f2t: error getting Pardf')
  endif

  call ESMF_FieldGet(tmpField, localDE=0, farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2Field, localDE=0, farrayPtr=q2,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swdField, localDE=0, farrayPtr=swd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField, localDE=0, farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField, localDE=0, farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField, localDE=0, farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField, localDE=0, farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(pcpField, localDE=0, farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_FieldGet(cpcpField, localDE=0, farrayPtr=cpcp,rc=status)
     call LIS_verify(status)
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_FieldGet(snowfField, localDE=0, farrayPtr=snowf,rc=status)
     call LIS_verify(status,'clsmf25_f2t: error retrieving snowf')
  endif

  if(LIS_FORC_SWnet%selectOpt.eq.1) then 
     call ESMF_FieldGet(swnetField, localDE=0, farrayPtr=swnet,rc=status)
     call LIS_verify(status,'clsmf25_f2t: error retrieving swnet')
  endif

  if(LIS_FORC_Pardr%selectOpt.eq.1) then 
     call ESMF_FieldGet(pardrField, localDE=0, farrayPtr=pardr,rc=status)
     call LIS_verify(status,'clsmf25_f2t: error retrieving pardr')
  endif

  if(LIS_FORC_Pardf%selectOpt.eq.1) then 
     call ESMF_FieldGet(pardfField, localDE=0, farrayPtr=pardf,rc=status)
     call LIS_verify(status,'clsmf25_f2t: error retrieving pardf')
  endif

  clsmf25_struc(n)%forc_count = clsmf25_struc(n)%forc_count +1
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
     clsmf25_struc(n)%met_force(t)%Tair = &
          clsmf25_struc(n)%met_force(t)%Tair + tmp(tid)
     clsmf25_struc(n)%met_force(t)%Qair = &
          clsmf25_struc(n)%met_force(t)%Qair + q2(tid)
     clsmf25_struc(n)%met_force(t)%Psurf = &
          clsmf25_struc(n)%met_force(t)%Psurf + psurf(tid)
     if(LIS_FORC_CRainf%selectOpt.eq.1) then 
        clsmf25_struc(n)%met_force(t)%Rainf_C = &
             clsmf25_struc(n)%met_force(t)%Rainf_C + cpcp(tid)
     endif
     if(LIS_FORC_SWnet%selectOpt.eq.1) then 
        clsmf25_struc(n)%met_force(t)%swnet = & 
             clsmf25_struc(n)%met_force(t)%swnet + swnet(tid)
     endif
     if(LIS_FORC_Pardr%selectOpt.eq.1) then 
        clsmf25_struc(n)%met_force(t)%pardr = & 
             clsmf25_struc(n)%met_force(t)%pardr + pardr(tid)
     endif
     if(LIS_FORC_Pardf%selectOpt.eq.1) then 
        clsmf25_struc(n)%met_force(t)%pardf = & 
             clsmf25_struc(n)%met_force(t)%pardf + pardf(tid)
     endif

     if ( LIS_FORC_Snowf%selectOpt.eq.1) then
        if(snowf(tid).eq.LIS_rc%udef) then 
           if(tmp(tid).ge.LIS_CONST_TKFRZ) then 
              clsmf25_struc(n)%met_force(t)%Rainf = &
                   clsmf25_struc(n)%met_force(t)%Rainf + pcp(tid)
              clsmf25_struc(n)%met_force(t)%Snowf = &
                   clsmf25_struc(n)%met_force(t)%Snowf + 0.0
           else
              clsmf25_struc(n)%met_force(t)%Rainf = &
                   clsmf25_struc(n)%met_force(t)%Rainf + 0.0
              clsmf25_struc(n)%met_force(t)%Rainf_C = & 
                   clsmf25_struc(n)%met_force(t)%Rainf_C + 0.0
              clsmf25_struc(n)%met_force(t)%Snowf = &
                   clsmf25_struc(n)%met_force(t)%Snowf + pcp(tid)
           endif
        else
           clsmf25_struc(n)%met_force(t)%Rainf = &
                clsmf25_struc(n)%met_force(t)%Rainf + pcp(tid)
           clsmf25_struc(n)%met_force(t)%Snowf = &
                clsmf25_struc(n)%met_force(t)%Snowf + snowf(tid)
        endif
     else
        if(tmp(tid).ge.LIS_CONST_TKFRZ) then 
           clsmf25_struc(n)%met_force(t)%Rainf = &
                clsmf25_struc(n)%met_force(t)%Rainf + pcp(tid)
           clsmf25_struc(n)%met_force(t)%Snowf = &
                clsmf25_struc(n)%met_force(t)%Snowf + 0.0
        else
           clsmf25_struc(n)%met_force(t)%Rainf = & 
                clsmf25_struc(n)%met_force(t)%Rainf + 0.0 
           clsmf25_struc(n)%met_force(t)%Rainf_C = &
                clsmf25_struc(n)%met_force(t)%Rainf_C + 0.0
           clsmf25_struc(n)%met_force(t)%Snowf = &
                clsmf25_struc(n)%met_force(t)%Snowf + pcp(tid)
        endif
     endif

     clsmf25_struc(n)%met_force(t)%LWdown = &
          clsmf25_struc(n)%met_force(t)%LWdown + lwd(tid)
     clsmf25_struc(n)%met_force(t)%SWdown = & 
          clsmf25_struc(n)%met_force(t)%SWdown +swd(tid)
     clsmf25_struc(n)%met_force(t)%Wind = & 
          clsmf25_struc(n)%met_force(t)%Wind + & 
          sqrt(uwind(tid)*uwind(tid)+&
          vwind(tid)*vwind(tid))
     
  ! Ensure that there are no negative precipitation values.
     if ( clsmf25_struc(n)%met_force(t)%Rainf /= LIS_rc%udef ) then
        clsmf25_struc(n)%met_force(t)%Rainf = &
             max(0.0,clsmf25_struc(n)%met_force(t)%Rainf)
     endif
     if ( clsmf25_struc(n)%met_force(t)%Rainf_C /= LIS_rc%udef ) then
        clsmf25_struc(n)%met_force(t)%Rainf_C = &
             max(0.0,clsmf25_struc(n)%met_force(t)%Rainf_C)
        clsmf25_struc(n)%met_force(t)%Rainf_C = &
             min(clsmf25_struc(n)%met_force(t)%Rainf_C,   &
             clsmf25_struc(n)%met_force(t)%Rainf)
     endif
     if ( clsmf25_struc(n)%met_force(t)%Snowf /= LIS_rc%udef ) then
        clsmf25_struc(n)%met_force(t)%Snowf = &
             max(0.0,clsmf25_struc(n)%met_force(t)%Snowf)
     endif
     
     ! zero or negative causes problems in the tubulence parameterization
     if ( clsmf25_struc(n)%met_force(t)%Wind /= LIS_rc%udef ) then
        clsmf25_struc(n)%met_force(t)%Wind =&
             max(0.0001,clsmf25_struc(n)%met_force(t)%Wind)
     endif
     
     if ( LIS_FORC_Forc_Hgt%selectOpt.eq.1 ) then
        call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Forc_Hgt%varname(1)),  &
             RefHField,rc=status)
        call LIS_verify(status,'clsmf25_f2t: error getting Forc_Hgt')
        call ESMF_FieldGet(RefHField,localDE=0,farrayPtr=RefH,rc=status)
        call LIS_verify(status)
        clsmf25_struc(n)%met_force(t)%RefH =  &
             clsmf25_struc(n)%met_force(t)%RefH + RefH(tid)
     else
        if (clsmf25_struc(n)%RefHfixed.lt.0.00001) then
           write(LIS_logunit,*) &
                 'Fixed reference height must not be equal to zero'
           call LIS_endrun
        endif
        clsmf25_struc(n)%met_force(t)%RefH = &
              clsmf25_struc(n)%met_force(t)%RefH + &
              clsmf25_struc(n)%RefHfixed
     endif
  enddo

end subroutine clsmf25_f2t
   
