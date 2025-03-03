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
! !ROUTINE: noah271_gce_f2t
! \label{noah271_gce_f2t}
!
! !REVISION HISTORY:
!  15 Oct 2002: Sujay Kumar; Initial Code
!
! !INTERFACE:
subroutine noah271_gce_f2t(n)
! !USES:      
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_localPet
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_logMod,         only : LIS_verify
  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the Noah
!  model tiles. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

  integer            :: t,status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField,fhgtField,acondField
  real,allocatable       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,allocatable       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,allocatable       :: fheight(:),acond(:)
  real               :: esat,e

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error getting Tair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error getting Qair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error getting SWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error getting LWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),uField,&
       rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error getting Wind_E')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),vField,&
       rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error getting Wind_N')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status, 'noah271_gce_f2t: error getting PSurf')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error getting Rainf')

  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
          rc=status)
     call LIS_verify(status,'noah271_gce_f2t: error getting CRainf')
  endif
  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),snowfField,&
          rc=status)
     call LIS_verify(status,'noah271_gce_f2t: error getting Snowf')
  endif
#if 0 

  if ( noah271_struc(n)%forcing_z > 0 ) then
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Forc_Hgt%varname(1)),&
          fhgtField, rc=status)
     call LIS_verify(status,'noah271_gce_f2t: error getting Forc_Hgt')
  endif

  if ( noah271_struc(n)%forcing_ch > 0 ) then
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Ch%varname(1)),&
          acondField,rc=status)
     call LIS_verify(status,'noah271_gce_f2t: error getting ch')
  endif
#endif

  call ESMF_FieldGet(tmpField, localDE=0, farrayPtr= tmp,rc=status)
  call LIS_verify(status, 'noah271_gce_f2t: error retrieving tmp')

  call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error retrieving q2')
  
  call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error retrieving swd')

  call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error retrieving lwd')

  call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error retrieving u')

  call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error retrieving v')
        
  call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error retrieving psurf')

  call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
  call LIS_verify(status,'noah271_gce_f2t: error retrieving pcp')

  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_FieldGet(cpcpField,localDE=0, farrayPtr=cpcp,rc=status)
     call LIS_verify(status, 'noah271_gce_f2t: error retrieving cpcp')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_FieldGet(snowfField,localDE=0, farrayPtr=snowf,rc=status)
     call LIS_verify(status,'noah271_gce_f2t: error retrieving snowf')
  endif
#if 0 
  if ( noah271_struc(n)%forcing_z > 0 ) then
     call ESMF_FieldGet(fhgtField,localDE=0, farrayPtr=fheight,rc=status)
     call LIS_verify(status,'noah271_gce_f2t: error retrieving fheight')
  endif

  if ( noah271_struc(n)%forcing_ch > 0 ) then
     call ESMF_FieldGet(acondField,localDE=0, farrayPtr=acond,rc=status)
     call LIS_verify(status,'noah271_gce_f2t: error retrieving acond')
  endif
#endif
  do t=1,LIS_rc%ntiles(n)
     noah271_struc(n)%noah(t)%tair=tmp(t)
     noah271_struc(n)%noah(t)%qair=q2(t)
     noah271_struc(n)%noah(t)%swdown=swd(t)
     noah271_struc(n)%noah(t)%lwdown=lwd(t)
     noah271_struc(n)%noah(t)%uwind=uwind(t)
     noah271_struc(n)%noah(t)%vwind=vwind(t)
     noah271_struc(n)%noah(t)%psurf=psurf(t)
     noah271_struc(n)%noah(t)%emiss=1.0

#if 0 
     ! Handle cases when input data does
     ! not have enough precision for this
     ! field
     if ( noah271_struc(n)%forcing_ch > 0 ) then
        if(acond(t).le.0) then 
           noah271_struc(n)%noah(t)%ch=0.00001
        else
           noah271_struc(n)%noah(t)%ch=acond(t)
        endif
     else
        noah271_struc(n)%noah(t)%ch=LIS_rc%udef
     endif

     !Note z is normally read from lis.config
     !In Noah section of card file, overwrite 
     !with forcing data if user so desires  
     if ( noah271_struc(n)%forcing_z > 0 ) then
        noah271_struc(n)%noah(t)%z = fheight(t)
     else
        noah271_struc(n)%noah(t)%z = noah271_struc(n)%obsz   
     endif
#endif
     if(pcp(t).ne.LIS_rc%udef) then 
        noah271_struc(n)%noah(t)%rainf=pcp(t)        
     else
        noah271_struc(n)%noah(t)%rainf=0.0
     endif
     if ( LIS_FORC_CRainf%selectOpt.eq.1) then
        if(cpcp(t).ne.LIS_rc%udef) then 
           noah271_struc(n)%noah(t)%rainf_c=cpcp(t)
        else
           noah271_struc(n)%noah(t)%rainf_c=0.0
        endif
     endif
     ! If there is snowf add it to precipitation.  Noah does not use
     ! separate rainf and snowf.  It determines what to do with precipitation.
     if ( LIS_FORC_Snowf%selectOpt.eq.1) then
        if(snowf(t).ne.LIS_rc%udef) then 
           noah271_struc(n)%noah(t)%rainf=noah271_struc(n)%noah(t)%rainf+snowf(t)
        endif
     endif
     noah271_struc(n)%noah(t)%snowf=0.0

     esat = e(noah271_struc(n)%noah(t)%tair)
     noah271_struc(n)%noah(t)%Q2SAT = 0.622 * ESAT /&
          (noah271_struc(n)%noah(t)%psurf- (1.-0.622)*ESAT)

     noah271_struc(n)%noah(t)%CH=  1.0E-4
     noah271_struc(n)%noah(t)%CM=  1.0E-4
  enddo
end subroutine noah271_gce_f2t
