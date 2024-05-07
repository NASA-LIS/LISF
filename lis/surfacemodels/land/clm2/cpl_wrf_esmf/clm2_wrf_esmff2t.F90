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
! 
! !ROUTINE: clm2_wrf_esmff2t
! \label{clm2_wrf_esmff2t}
!
! !INTERFACE: 
  subroutine clm2_wrf_esmff2t(n)
! !USES: 
    use ESMF
    use LIS_precisionMod
    use clm2_lsmMod
    use clm2_varcon  , only : rair, cpair, po2, pco2, tcrit, tfrz    
    use LIS_coreMod, only : LIS_rc
    use LIS_metforcingMod, only : LIS_FORC_State
    use LIS_logMod,         only : LIS_verify

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n 

!
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the CLM
!  model tiles, when used in the mode coupled to WRF. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP
    real(r8), allocatable :: solar(:)
    real(r8), allocatable :: prcp(:)
    real(r8) :: forc_vp
    integer            :: t,status
    type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
    type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
    type(ESMF_Field)   :: zField, q2satField, chField,coszField
    real,allocatable       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
    real,allocatable       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
    real,allocatable       :: zval(:),q2sat(:),ch(:),cosz(:)
    
    allocate(solar(LIS_rc%ntiles(n)))
    allocate(prcp(LIS_rc%ntiles(n)))    

    call ESMF_StateGet(LIS_FORC_State(n),"Near Surface Air Temperature",tmpField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Near Surface Specific Humidity",q2Field,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Incident Shortwave Radiation",swdField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Incident Longwave Radiation",lwdField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Eastward Wind",uField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Northward Wind",vField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Surface Pressure",psurfField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Rainfall Rate",pcpField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Convective Rainfall Rate",cpcpField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Snowfall Rate",snowfField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_StateGet(LIS_FORC_State(n),"Surface Exchange Coefficient for Heat",chField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Saturated Mixing Ratio",q2satField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_StateGet(LIS_FORC_State(n),"Height of Atmospheric Forcing",zField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_StateGet(LIS_FORC_State(n),"Cosine of Zenith Angle",coszField,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(chField,localDE=0,farrayPtr=ch,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(q2satField,localDE=0,farrayPtr=q2sat,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(zField,localDE=0,farrayPtr=zval,rc=status)
    call LIS_verify(status)   

    call ESMF_FieldGet(coszField,localDE=0,farrayPtr=cosz,rc=status)
    call LIS_verify(status)   
 
    do t=1,LIS_rc%ntiles(n)
       clm2_struc(n)%clm(t)%forc_t        = tmp(t)
       clm2_struc(n)%clm(t)%forc_q        = q2(t)
       solar(t)             = swd(t)
       clm2_struc(n)%clm(t)%forc_solad(1) = solar(t)*35.0/100.0
       clm2_struc(n)%clm(t)%forc_solad(2) = solar(t)*35.0/100.0
       clm2_struc(n)%clm(t)%forc_solai(1) = solar(t)*15.0/100.0
       clm2_struc(n)%clm(t)%forc_solai(2) = solar(t)*15.0/100.0
       clm2_struc(n)%clm(t)%forc_lwrad    = lwd(t)
       clm2_struc(n)%clm(t)%forc_u        = uwind(t)
       clm2_struc(n)%clm(t)%forc_v        = vwind(t)
       clm2_struc(n)%clm(t)%forc_pbot     = psurf(t)
       prcp(t)              = pcp(t)
       
       clm2_struc(n)%clm(t)%forc_hgt      = 0.5*zval(t)
       clm2_struc(n)%clm(t)%forc_ch       = ch(t)
       clm2_struc(n)%clm(t)%forc_q2sat    = q2sat(t)
       clm2_struc(n)%clm(t)%forc_cosz     = cosz(t)
   
       clm2_struc(n)%clm(t)%forc_hgt_u    = clm2_struc(n)%clm(t)%forc_hgt+&
           clm2_struc(n)%clm(t)%displa+clm2_struc(n)%clm(t)%z0m
       clm2_struc(n)%clm(t)%forc_hgt_t    = clm2_struc(n)%clm(t)%forc_hgt+&
           clm2_struc(n)%clm(t)%displa+clm2_struc(n)%clm(t)%z0m
       clm2_struc(n)%clm(t)%forc_hgt_q    = clm2_struc(n)%clm(t)%forc_hgt+&
           clm2_struc(n)%clm(t)%displa+clm2_struc(n)%clm(t)%z0m
   
       if (prcp(t) > 0.) then
          if (clm2_struc(n)%clm(t)%forc_t > (tfrz + tcrit)) then
             clm2_struc(n)%clm(t)%itypprc   = 1
             clm2_struc(n)%clm(t)%forc_rain = prcp(t)
             clm2_struc(n)%clm(t)%forc_snow = 0.
          else
             clm2_struc(n)%clm(t)%itypprc   = 2
             clm2_struc(n)%clm(t)%forc_rain = 0.
             clm2_struc(n)%clm(t)%forc_snow = prcp(t)
          endif
       else
          clm2_struc(n)%clm(t)%itypprc   = 0
          clm2_struc(n)%clm(t)%forc_rain = 0.
          clm2_struc(n)%clm(t)%forc_snow = 0
       endif

! Derive new fields (potential temperature, vapor pressure, 
! air density, CO2, and O2) and copy solar radiations

!=== Cuurently potential temperature set to 2 m temperature since 
!=== we only get surface pressure in our forcing and elevation differences 
!=== are accounted for in for.f

!===LDAS modification: Slight change to be consistent with our forcing dataset
!          clm2_struc(n)%clm(t)%forc_th  = clm2_struc(n)%clm(t)%forc_t * (clm2_struc(n)%clm(t)%forc_psrf/clm2_struc(n)%clm(t)%forc_pbot)**(rair/cpair)
       clm2_struc(n)%clm(t)%forc_th  = clm2_struc(n)%clm(t)%forc_t * &
            (clm2_struc(n)%clm(t)%forc_pbot/clm2_struc(n)%clm(t)%forc_pbot)**(rair/cpair)
       forc_vp  = clm2_struc(n)%clm(t)%forc_q*clm2_struc(n)%clm(t)%forc_pbot / &
            (0.622+0.378*clm2_struc(n)%clm(t)%forc_q)
       clm2_struc(n)%clm(t)%forc_rho = &
            (clm2_struc(n)%clm(t)%forc_pbot-0.378*forc_vp) / &
            (rair*clm2_struc(n)%clm(t)%forc_t)
    enddo
    
    deallocate(solar)
    deallocate(prcp)
  end subroutine clm2_wrf_f2t

