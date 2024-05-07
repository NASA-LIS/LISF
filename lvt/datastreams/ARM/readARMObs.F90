!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readARMObs
! \label{readARMObs}
!
! !INTERFACE: 
subroutine readARMObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_logMod,       only : LVT_logunit, LVT_verify
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_calendar, LVT_tick
  use ARM_obsMod,      only : armobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for ARM station data. 
! The plugin processes three types of data (1) flux and soil mositure 
! estimates from the energy balance bowen ratio (ebbr) stations, (2) 
! flux estimates from the eddy correlaton (ecor) flux measurement
! system, (3) Surface Meteorological Observing System (SMOS) instruments
! and  (4) soil moisture and temperature profiles from 
! the soil water and temperature system (SWATS). When both 
! ebbr and ecor-based fluxes are available, ebbr fluxes are 
! chosen. 
!
!  NOTES: 
!   The use of SWATS data is not recommended
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 Nov 2010: Sujay Kumar, Initial Specification
! 
!EOP
!BOP
! !ARGUMENTS: 

!EOP
  type(ESMF_Time)  :: armtime1
  integer          :: t,c,r
  real             :: gmt
  real                :: time
  type(ESMF_TimeInterval) :: dayInterval
  type(ESMF_Time)         :: initTime
  integer             :: yr, mo, da, hr, mn, ss, doy
  integer             :: status
  real                :: qle(LVT_rc%lnc, LVT_rc%lnr)
  real                :: qh(LVT_rc%lnc, LVT_rc%lnr)
  real                :: qg(LVT_rc%lnc, LVT_rc%lnr)
  real                :: netrad(LVT_rc%lnc, LVT_rc%lnr)
  real                :: solarrad(LVT_rc%lnc, LVT_rc%lnr)
  real                :: br(LVT_rc%lnc, LVT_rc%lnr)
  real                :: sfsm(LVT_rc%lnc, LVT_rc%lnr)
  real                :: rzsm(LVT_rc%lnc, LVT_rc%lnr)
  real                :: sfst(LVT_rc%lnc, LVT_rc%lnr)
  real                :: rzst(LVT_rc%lnc, LVT_rc%lnr)
  
  real                :: pcp(LVT_rc%lnc,LVT_rc%lnr)
  real                :: sh(LVT_rc%lnc,LVT_rc%lnr)
  real                :: prs(LVT_rc%lnc,LVT_rc%lnr)
  real                :: wspd(LVT_rc%lnc,LVT_rc%lnr)
  real                :: tair(LVT_rc%lnc,LVT_rc%lnr)
  real                :: snow(LVT_rc%lnc,LVT_rc%lnr)

  qle  = LVT_rc%udef 
  qh   = LVT_rc%udef 
  qg   = LVT_rc%udef 

  netrad   = LVT_rc%udef 
  solarrad = LVT_rc%udef 

  br   = LVT_rc%udef 
  sfsm  = LVT_rc%udef 
  rzsm  = LVT_rc%udef 
  sfst  = LVT_rc%udef 
  rzst  = LVT_rc%udef 
  pcp = LVT_rc%udef 
  prs = LVT_rc%udef 
  wspd = LVT_rc%udef 
  tair = LVT_rc%udef 
  snow = LVT_rc%udef 

  time = LVT_rc%dhr(source)*3600+LVT_rc%dmn(source)*60+LVT_rc%dss(source)
  if((mod(time,86400.0).eq.0.0).or.(LVT_rc%dda(source).ne.armobs(source)%da)) then 

     call ESMF_TimeSet(armobs(source)%starttime, yy=LVT_rc%dyr(source), &
          mm = LVT_rc%dmo(source), &
          dd = LVT_rc%dda(source), &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status,'error in timeset: readArmobs(Source)')
     armobs(source)%da = LVT_rc%dda(source)

     call process_baebbr_flux_data(source,LVT_rc%dyr(source),&
          LVT_rc%dmo(source),LVT_rc%dda(source))
     call process_ecor_flux_data(source,LVT_rc%dyr(source),&
          LVT_rc%dmo(source),LVT_rc%dda(source))
     call process_swats_sm_data(source,LVT_rc%dyr(source),&
          LVT_rc%dmo(source),LVT_rc%dda(source))
     call process_smos_data(source,LVT_rc%dyr(source),&
          LVT_rc%dmo(source),LVT_rc%dda(source))
     call process_ebbr_data(source,LVT_rc%dyr(source),&
          LVT_rc%dmo(source),LVT_rc%dda(source))

  endif

  call ESMF_TimeSet(armtime1,yy=LVT_rc%dyr(source), &
       mm = LVT_rc%dmo(source), &
       dd = LVT_rc%dda(source), &
       h = LVT_rc%dhr(source), &
       m = LVT_rc%dmn(source), &
       calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'error in timeset: readArmobs(Source)')
  
  t = nint((armtime1 - armobs(source)%starttime)/armobs(source)%ecor_ts) +1

  if(mod(time,1800.0).eq.0.0) then 
     call retrieve_ecor_flux_data(source,t,qle,qh)
  endif

  t = nint((armtime1 - armobs(source)%starttime)/armobs(source)%baebbr_ts) +1
!process the flux data only at 1/2 hour intervals
  if(mod(time,1800.0).eq.0.0) then 
     call retrieve_baebbr_flux_data(source,t,qle,qh,qg,&
          netrad,solarrad)
  endif

  t = nint((armtime1 - armobs(source)%starttime)/armobs(source)%swats_ts) +1
!process the flux data only at 1 hour intervals
  if(mod(time,3600.0).eq.0.0) then 
     call retrieve_swats_smst_data(source,t,sfsm,rzsm,&
          sfst,rzst)
  endif

  t = nint((armtime1 - armobs(source)%starttime)/armobs(source)%smos_ts) +1
  if(mod(time,1800.0).eq.0.0) then 
     call retrieve_smos_data(source,t,pcp,wspd,&
          tair,prs,snow,sh)
  endif

  t = nint((armtime1 - armobs(source)%starttime)/armobs(source)%ebbr_ts) +1
  if(mod(time,1800.0).eq.0.0) then 
     call retrieve_ebbr_data(source,t,sfsm,sfst)
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_qle,source,qle,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_qh,source,qh,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_qg,source,qg,vlevel=1,units="W/m2",dir="DN")
  call LVT_logSingleDataStreamVar(LVT_MOC_rnet,source,netrad,vlevel=1,units="W/m2",dir="DN")
  call LVT_logSingleDataStreamVar(LVT_MOC_swdownforc,source,solarrad,vlevel=1,units="W/m2",dir="DN")
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source,sfsm, vlevel=1,units="m3/m3")
  do r=1, LVT_rc%lnr
     do c=1, LVT_rc%lnc
        if(sfsm(c,r).ne.-9999.0) then 
           sfsm(c,r) = sfsm(c,r)*LVT_rc%lis_sf_d*1000.0
        endif
     enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source,sfsm, vlevel=1,units="kg/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_soiltemp, source,sfst, vlevel=1,units="K")
  do r=1, LVT_rc%lnr
     do c=1, LVT_rc%lnc
        if(rzsm(c,r).ne.-9999.0) then 
           rzsm(c,r) = rzsm(c,r)*LVT_rc%lis_rz_d*1000.0
        endif
     enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_rootmoist, source,rzsm,vlevel=1,units="m3/m3")
  call LVT_logSingleDataStreamVar(LVT_MOC_roottemp, source,rzst,vlevel=1,units="K")

  call LVT_logSingleDataStreamVar(LVT_MOC_snowdepth, source,snow,vlevel=1,units="m")

  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,pcp,vlevel=1,units="kg/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_windforc,source,wspd,vlevel=1,units="m/s")
  call LVT_logSingleDataStreamVar(LVT_MOC_psurfforc,source,prs,vlevel=1,units="Pa")
  call LVT_logSingleDataStreamVar(LVT_MOC_tairforc,source,tair,vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_qairforc,source,sh,vlevel=1,units="kg/kg")

  call LVT_logSingleDataStreamVar(LVT_MOC_t2diag,source,tair,vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_q2diag,source,sh,vlevel=1,units="kg/kg")
  
end subroutine readArmobs


!BOP
! 
! !ROUTINE: process_baebbr_flux_data
!  \label{process_baebbr_flux_data}
!
! !INTERFACE:
subroutine process_baebbr_flux_data(source,yr,mo,da)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod,      only : LVT_logunit
  use ARM_obsMod,      only : armobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the bulk aerodynamic estimates of 
!  sensible,latent and ground heat flux
!  measurements from the Energy Balance Bowen Ratio (BAEBBR) 
!  system. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[yr]      year for which the data is being processed
!   \item[mo]      month for which the data is being processed
!   \item[da]      day for which the data is being processed
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer          :: source
  integer          :: yr
  integer          :: mo
  integer          :: da
!
!EOP
  integer          :: i 
  character*100    :: filename
  integer          :: status

  if(armobs(source)%baebbr_select.eq.1) then 
  ! Read the baebbr data
     do i=1,armobs(source)%n_stns
        call create_arm_baebbr_flux_filename(armobs(source)%odir, &
             armobs(source)%site_id, armobs(source)%stn_name(i), &
             yr, mo, da, filename,status)
        if(status.eq.0) then 
           write(LVT_logunit,*) '[INFO] Reading baebbr file   ',trim(filename)
           call read_baebbr_flux_file(source, i,yr, mo, da, &
                filename)
        else
           armobs(source)%baebbr_qle(i,:) = LVT_rc%udef
           armobs(source)%baebbr_qh(i,:) = LVT_rc%udef
           armobs(source)%baebbr_qg(i,:) = LVT_rc%udef
           armobs(source)%baebbr_netrad(i,:) = LVT_rc%udef
           armobs(source)%baebbr_solarrad(i,:) = LVT_rc%udef
        endif
     enddo
  endif

end subroutine process_baebbr_flux_data

!BOP
! 
! !ROUTINE: process_ecor_flux_data
! \label{process_ecor_flux_data}
!
! !INTERFACE:
subroutine process_ecor_flux_data(source,yr,mo,da)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod,      only : LVT_logunit
  use ARM_obsMod,      only : armobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the sensible and latent heat flux
!  measurements from the Eddy Correlation (ECOR) flux 
!  measurement system.
!
!  The arguments are: 
!  \begin{description}
!   \item[yr]      year for which the data is being processed
!   \item[mo]      month for which the data is being processed
!   \item[da]      day for which the data is being processed
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
  integer          :: source
  integer          :: yr
  integer          :: mo
  integer          :: da
!EOP 
  integer          :: i 
  character*100    :: filename
  integer          :: status

  if(armobs(source)%ecor_select.eq.1) then 
  ! read the ecor data
     do i=1,armobs(source)%n_stns
        call create_arm_ecor_flux_filename(armobs(source)%odir, &
             armobs(source)%site_id, armobs(source)%stn_name(i), &
             yr, mo, da, filename,status)
        if(status.eq.0) then 
           write(LVT_logunit,*) '[INFO] Reading ecor file  ',trim(filename)
           !exclude E14 ECOR data. 
           if(trim(armobs(source)%stn_name(i)).ne."E14") then 
              call read_ecor_flux_file(source,i, yr, mo, da, &
                   filename)
           endif
        else
           armobs(source)%ecor_qh(i,:) = LVT_rc%udef
           armobs(source)%ecor_qle(i,:) = LVT_rc%udef
        endif
     enddo
  endif
end subroutine process_ecor_flux_data

!BOP
! 
! !ROUTINE: process_swats_sm_data
! \label{process_swats_sm_data}
!
! !INTERFACE:
subroutine process_swats_sm_data(source,yr,mo,da)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc
  use ARM_obsMod,      only : armobs
  use LVT_logMod,      only : LVT_logunit
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the soil moisture and temperature
!  measurements from the Soil Water and Temperature (SWATS)
!  system.
!
!  The arguments are: 
!  \begin{description}
!   \item[yr]      year for which the data is being processed
!   \item[mo]      month for which the data is being processed
!   \item[da]      day for which the data is being processed
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
!
! !ARGUMENTS: 
  integer          :: source
  integer          :: yr
  integer          :: mo
  integer          :: da
!
!EOP
  integer          :: i 
  character*100    :: filename
  integer          :: status

  if(armobs(source)%swats_select.eq.1) then 
  ! read the swats data
     do i=1,armobs(source)%n_stns
        call create_arm_swats_filename(armobs(source)%odir, &
             armobs(source)%site_id, armobs(source)%stn_name(i), &
             yr, mo, da, filename,status)
        if(status.eq.0) then 
           write(LVT_logunit,*) '[INFO] Reading swats file ',trim(filename)
           call read_swats_file(source, i, yr, mo, da, &
                filename)
        else
           armobs(source)%swats_sfsm(i,:) = LVT_rc%udef
           armobs(source)%swats_rzsm(i,:) = LVT_rc%udef
           armobs(source)%swats_sfst(i,:) = LVT_rc%udef
           armobs(source)%swats_rzst(i,:) = LVT_rc%udef
        endif
     enddo
  endif
end subroutine process_swats_sm_data

!BOP
! 
! !ROUTINE: process_smos_data
! \label{process_smos_data}
!
! !INTERFACE:
subroutine process_smos_data(source,yr,mo,da)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod,      only : LVT_logunit
  use ARM_obsMod,      only : armobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the surface meteorology 
!  measurements from the Surface Meteorological Observation
!  System (SMOS)
!
!  The arguments are: 
!  \begin{description}
!   \item[yr]      year for which the data is being processed
!   \item[mo]      month for which the data is being processed
!   \item[da]      day for which the data is being processed
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
  integer          :: source
  integer          :: yr
  integer          :: mo
  integer          :: da
!EOP
  integer          :: i 
  character*100    :: filename
  integer          :: status

  if(armobs(source)%smos_select.eq.1) then 
  ! read the ecor data
     do i=1,armobs(source)%n_stns
        call create_arm_smos_filename(armobs(source)%odir, &
             armobs(source)%site_id, armobs(source)%stn_name(i), &
             yr, mo, da, filename,status)
        if(status.eq.0) then 
           write(LVT_logunit,*) '[INFO] Reading SMOS file ',trim(filename)
!exclude E14 ECOR data. 
           if(trim(armobs(source)%stn_name(i)).ne."E14") then 
              call read_smos_file(source, i, yr, mo, da, &
                   filename)
           endif
        else
           armobs(source)%smos_pcp(i,:) = LVT_rc%udef
           armobs(source)%smos_prs(i,:) = LVT_rc%udef
           armobs(source)%smos_wspd(i,:) = LVT_rc%udef
           armobs(source)%smos_temp(i,:) = LVT_rc%udef
           armobs(source)%smos_snow(i,:) = LVT_rc%udef
           armobs(source)%smos_sh(i,:) = LVT_rc%udef
        endif
     enddo
  endif
end subroutine process_smos_data

!BOP
! 
! !ROUTINE: process_ebbr_data
!  \label{process_ebbr_data}
!
! !INTERFACE:
subroutine process_ebbr_data(source,yr,mo,da)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod,      only : LVT_logunit
  use ARM_obsMod,      only : armobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the bulk aerodynamic estimates of 
!  sensible,latent and ground heat flux
!  measurements from the Energy Balance Bowen Ratio (BAEBBR) 
!  system. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[yr]      year for which the data is being processed
!   \item[mo]      month for which the data is being processed
!   \item[da]      day for which the data is being processed
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer          :: source
  integer          :: yr
  integer          :: mo
  integer          :: da
!EOP
  integer          :: i 
  character*100    :: filename
  integer          :: status

  if(armobs(source)%ebbr_select.eq.1) then 
  ! Read the baebbr data
     do i=1,armobs(source)%n_stns
        call create_arm_ebbr_filename(armobs(source)%odir, &
             armobs(source)%site_id, armobs(source)%stn_name(i), &
             yr, mo, da, filename,status)
        if(status.eq.0) then 
           write(LVT_logunit,*) '[INFO] Reading ebbr file ',trim(filename)
           call read_ebbr_file(source, i,yr, mo, da, &
                filename)
        else
           armobs(source)%ebbr_sfsm(i,:) = LVT_rc%udef
           armobs(source)%ebbr_sfst(i,:) = LVT_rc%udef
        endif
     enddo
  endif

end subroutine process_ebbr_data

!BOP
! 
! !ROUTINE: retrieve_baebbr_flux_data
! \label{retrieve_baebbr_flux_data}
!
! !INTERFACE: 
subroutine retrieve_baebbr_flux_data(source,t,qle,qh,qg, &
     netrad,  solarrad)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc, LVT_domain
  use ARM_obsMod,      only : armobs
  use LVT_logMod,      only : LVT_logunit
  use map_utils
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine temporally aggregates the BAEBBR flux measurement
!  data up to the LIS output timestep 
!
!  The arguments are: 
!  \begin{description}
!   \item[t]  time index                
!   \item[qle]  array containting aggregated latent heat flux  
!   \item[qh]   array containting aggregated sensible heat flux    
!   \item[qg]   array containting aggregated ground heat flux  
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
  integer          :: source
  integer          :: t
  real             :: qle(LVT_rc%lnc,LVT_rc%lnr)
  real             :: qh(LVT_rc%lnc,LVT_rc%lnr)
  real             :: qg(LVT_rc%lnc,LVT_rc%lnr)
  real             :: netrad(LVT_rc%lnc,LVT_rc%lnr)
  real             :: solarrad(LVT_rc%lnc,LVT_rc%lnr)
!
!EOP

  integer          :: c,stn_col,stn_row
  real             :: col,row
  
  if(armobs(source)%baebbr_select.eq.1) then 
     do c=1,armobs(source)%n_stns
        call latlon_to_ij(LVT_domain%lvtproj, armobs(source)%stnlat(c), &
             armobs(source)%stnlon(c), col, row)
        stn_col = nint(col)
        stn_row = nint(row)
        
        if(armobs(source)%baebbr_qle(c,t).ne.LVT_rc%udef) then 
           
           qle(stn_col, stn_row) = &
                (-1)*armobs(source)%baebbr_qle(c,t) 
        endif
        
        if(armobs(source)%baebbr_qh(c,t).ne.LVT_rc%udef) then 
           qh(stn_col, stn_row) = &
                (-1)*armobs(source)%baebbr_qh(c,t) 
        endif
        
        if(armobs(source)%baebbr_qg(c,t).ne.LVT_rc%udef) then 
           qg(stn_col, stn_row) = qg(stn_col, stn_row) + & 
                armobs(source)%baebbr_qg(c,t) 
        endif
        if(armobs(source)%baebbr_netrad(c,t).ne.LVT_rc%udef) then 
           netrad(stn_col, stn_row) = & 
                armobs(source)%baebbr_netrad(c,t) 
        endif
        if(armobs(source)%baebbr_solarrad(c,t).ne.LVT_rc%udef) then 
           solarrad(stn_col, stn_row) = &
                armobs(source)%baebbr_solarrad(c,t) 
        endif
     enddo
  endif
end subroutine retrieve_baebbr_flux_data
   
!BOP
! 
! !ROUTINE: retrieve_ecor_flux_data
! \label{retrieve_ecor_flux_data}
!
! !INTERFACE: 
subroutine retrieve_ecor_flux_data(source,t,qle,qh)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc, LVT_domain
  use ARM_obsMod,      only : armobs
  use LVT_logMod,      only : LVT_logunit
  use map_utils
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine temporally aggregates the ECOR flux measurement
!  data up to the LIS output timestep 
!
!  The arguments are: 
!  \begin{description}
!   \item[t]  time index                
!   \item[qle]  array containting aggregated latent heat flux  
!   \item[qh]   array containting aggregated sensible heat flux    
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer          :: source
  integer          :: st,et
  real             :: qle(LVT_rc%lnc,LVT_rc%lnr)
  real             :: qh(LVT_rc%lnc,LVT_rc%lnr)
!EOP

  integer          :: c,t,stn_col,stn_row
  real             :: col,row

  if(armobs(source)%ecor_select.eq.1) then 
     do c=1,armobs(source)%n_stns
        call latlon_to_ij(LVT_domain%lvtproj, armobs(source)%stnlat(c), &
             armobs(source)%stnlon(c), col, row)
        stn_col = nint(col)
        stn_row = nint(row)
        
        if(stn_row.ge.1.and.stn_row.le.LVT_rc%lnr.and.&
             stn_col.ge.1.and.stn_col.le.LVT_rc%lnc) then 
           if(armobs(source)%ecor_qle(c,t).ne.LVT_rc%udef) then 
              qle(stn_col, stn_row) = & 
                   armobs(source)%ecor_qle(c,t) 
           endif
           
           if(armobs(source)%ecor_qh(c,t).ne.LVT_rc%udef) then 
              qh(stn_col, stn_row) = & 
                   armobs(source)%ecor_qh(c,t) 
           endif
        endif
     enddo
  endif
end subroutine retrieve_ecor_flux_data

!BOP
! 
! !ROUTINE: retrieve_swats_smst_data
! \label{retrieve_swats_smst_data}
!
! !INTERFACE: 
subroutine retrieve_swats_smst_data(source,t,sfsm,rzsm,&
     sfst,rzst)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc, LVT_domain
  use ARM_obsMod,      only : armobs
  use LVT_logMod,      only : LVT_logunit
  use map_utils
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine temporally aggregates the SWATS
!  data up to the LIS output timestep 
!
!  The arguments are: 
!  \begin{description}
!   \item[t]  time index                
!   \item[sfsm]  array containting aggregated surface soil moisture
!   \item[rzsm]   array containting aggregated root zone soil moisture
!   \item[sfst]  array containting aggregated surface soil temperature
!   \item[rzst]   array containting aggregated root zone soil temperature
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer          :: source
  integer          :: t
  real             :: sfsm(LVT_rc%lnc,LVT_rc%lnr)
  real             :: rzsm(LVT_rc%lnc,LVT_rc%lnr)
  real             :: sfst(LVT_rc%lnc,LVT_rc%lnr)
  real             :: rzst(LVT_rc%lnc,LVT_rc%lnr)
!EOP

  integer          :: c,stn_col,stn_row
  real             :: col,row

  if(armobs(source)%swats_select.eq.1) then 
     do c=1,armobs(source)%n_stns
        call latlon_to_ij(LVT_domain%lvtproj, armobs(source)%stnlat(c), &
             armobs(source)%stnlon(c), col, row)
        stn_col = nint(col)
        stn_row = nint(row)
        
        if(armobs(source)%swats_sfsm(c,t).ne.LVT_rc%udef) then 
           sfsm(stn_col, stn_row) = &
                armobs(source)%swats_sfsm(c,t) 
        endif
        
        if(armobs(source)%swats_rzsm(c,t).ne.LVT_rc%udef) then 
           rzsm(stn_col, stn_row) = &
                armobs(source)%swats_rzsm(c,t) 
        endif
        
        if(armobs(source)%swats_sfst(c,t).ne.LVT_rc%udef) then 
           sfst(stn_col, stn_row) =& 
                armobs(source)%swats_sfst(c,t) 
        endif
        if(armobs(source)%swats_rzst(c,t).ne.LVT_rc%udef) then 
           rzst(stn_col, stn_row) = &
                armobs(source)%swats_rzst(c,t) 
        endif
     enddo
  endif
end subroutine retrieve_swats_smst_data

!BOP
! 
! !ROUTINE: retrieve_smos_data
! \label{retrieve_smos_data}
!
! !INTERFACE: 
subroutine retrieve_smos_data(source,t,pcp,wspd,&
             tair,prs,snow,sh)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc, LVT_domain
  use ARM_obsMod,      only : armobs
  use LVT_logMod,      only : LVT_logunit
  use map_utils
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine temporally aggregates the SWATS
!  data up to the LIS output timestep 
!
!  The arguments are: 
!  \begin{description}
!   \item[t]  time index                
!   \item[pcp]  array containting aggregated precipitation
!   \item[wspd]   array containting aggregated wind speed
!   \item[tair]  array containting aggregated temperature
!   \item[prs]   array containting aggregated pressure
!   \item[snow]   array containting aggregated snow depth
!   \item[sh]   array containting aggregated specific humidity
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
  integer          :: source
  integer          :: t
  real             :: pcp(LVT_rc%lnc,LVT_rc%lnr)
  real             :: wspd(LVT_rc%lnc,LVT_rc%lnr)
  real             :: tair(LVT_rc%lnc,LVT_rc%lnr)
  real             :: prs(LVT_rc%lnc,LVT_rc%lnr)
  real             :: snow(LVT_rc%lnc,LVT_rc%lnr)
  real             :: sh(LVT_rc%lnc,LVT_rc%lnr)
!
!EOP

  integer          :: c,stn_col,stn_row
  real             :: col,row

  if(armobs(source)%smos_select.eq.1) then 
     do c=1,armobs(source)%n_stns
        call latlon_to_ij(LVT_domain%lvtproj, armobs(source)%stnlat(c), &
             armobs(source)%stnlon(c), col, row)
        stn_col = nint(col)
        stn_row = nint(row)
        
        if(armobs(source)%smos_pcp(c,t).ne.LVT_rc%udef) then 
           pcp(stn_col, stn_row) =&
                armobs(source)%smos_pcp(c,t) 
        endif
        
        if(armobs(source)%smos_wspd(c,t).ne.LVT_rc%udef) then 
           wspd(stn_col, stn_row) = &
                armobs(source)%smos_wspd(c,t) 
        endif
        
        if(armobs(source)%smos_temp(c,t).ne.LVT_rc%udef) then 
           tair(stn_col, stn_row) = &
                armobs(source)%smos_temp(c,t) 
        endif
        if(armobs(source)%smos_prs(c,t).ne.LVT_rc%udef) then 
           prs(stn_col, stn_row) = &
                armobs(source)%smos_prs(c,t) 
        endif
        if(armobs(source)%smos_snow(c,t).ne.LVT_rc%udef) then 
           snow(stn_col, stn_row) = &
                armobs(source)%smos_snow(c,t) 
        endif
        if(armobs(source)%smos_sh(c,t).ne.LVT_rc%udef) then 
           sh(stn_col, stn_row) = &
                armobs(source)%smos_sh(c,t) 
        endif
     enddo
  endif
end subroutine retrieve_smos_data

!BOP
! 
! !ROUTINE: retrieve_ebbr_data
! \label{retrieve_ebbr_data}
!
! !INTERFACE: 
subroutine retrieve_ebbr_data(source,t,sfsm,sfst)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc, LVT_domain
  use ARM_obsMod,      only : armobs
  use LVT_logMod,      only : LVT_logunit
  use map_utils
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine temporally aggregates the BAEBBR flux measurement
!  data up to the LIS output timestep 
!
!  The arguments are: 
!  \begin{description}
!   \item[t]  time index                
!   \item[qle]  array containting aggregated latent heat flux  
!   \item[qh]   array containting aggregated sensible heat flux    
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
  integer          :: source
  integer          :: t
  real             :: sfsm(LVT_rc%lnc,LVT_rc%lnr)
  real             :: sfst(LVT_rc%lnc,LVT_rc%lnr)

!
!EOP

  integer          :: c,stn_col,stn_row
  real             :: col,row
  
  if(armobs(source)%ebbr_select.eq.1) then 
     do c=1,armobs(source)%n_stns
        call latlon_to_ij(LVT_domain%lvtproj, armobs(source)%stnlat(c), &
             armobs(source)%stnlon(c), col, row)
        stn_col = nint(col)
        stn_row = nint(row)
        
        if(armobs(source)%ebbr_sfsm(c,t).ne.LVT_rc%udef) then 
           sfsm(stn_col, stn_row) =&
                armobs(source)%ebbr_sfsm(c,t) 
        endif
        if(armobs(source)%ebbr_sfst(c,t).ne.LVT_rc%udef) then 
           sfst(stn_col, stn_row) =&
                armobs(source)%ebbr_sfst(c,t) 
        endif
     enddo
  endif
end subroutine retrieve_ebbr_data
   

!BOP
! 
! !ROUTINE: read_swats_file
! \label{read_swats_file}
!
! !INTERFACE: 
subroutine read_swats_file(source, k, yr, mo, da, filename)
! 
! !USES: 
  use ESMF     
  use LVT_coreMod,      only : LVT_rc
  use LVT_logMod,       only : LVT_logunit, LVT_verify
  use LVT_timeMgrMod,   only : LVT_calendar
  use ARM_obsMod,      only : armobs
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!DESCRIPTION: 
!  This routine reads the soil moisture and soil temperature
!  profiles from the Soil Water and Temperature System (SWATS) 
!  for a particular station. The data is read from the 
!  native NETCDF files, and is processed by applying the 
!  quality control flags. Only data points with 'good'
!  classification (qc value=0) is chosen. This routine 
!  also computes the temporal offset of the data relative
!  to the 0z of a particular day. 
!
!  The arguments are: 
!  \begin{description}
!   \item[k]         Index of the ARM station
!   \item[yr]        year of ARM data 
!   \item[mo]        month of ARM data
!   \item[da]        day of ARM data
!   \item[filename]  Name of the SWATS file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer,  intent(in) :: source
  integer,  intent(in) :: k 
  integer              :: yr
  integer              :: mo
  integer              :: da
  character(len=*)     :: filename 
  integer              :: bend
!EOP
  integer       :: nid, timeid,btimeid, latid, lonid, dim1ID, dim2ID
  integer       :: depthid,tsoilEid, tsoilWid, qc_tsoilEid, qc_tsoilWid
  integer       :: watcontEid,watcontWid,qc_watcontEid,qc_watcontWid
  integer       :: tdims, ddims
  real          :: lat, lon 
  integer       :: base_time
  integer       :: ios
  real, allocatable :: time(:)
  real, allocatable :: depth(:)
  real, allocatable :: smE(:,:)
  real, allocatable :: smW(:,:)
  real, allocatable :: stE(:,:)
  real, allocatable :: stW(:,:)
  real, allocatable :: qc_smE(:,:)
  real, allocatable :: qc_smW(:,:)
  real, allocatable :: qc_stE(:,:)
  real, allocatable :: qc_stW(:,:)
  integer           :: kk,t
  type(ESMF_Time)   :: reftime, datatime, currtime
  type(ESMF_TimeInterval) :: dt
  integer           :: yr1, mo1, da1, hr1, mn1, ss1
  integer           :: status
  integer           :: data_index
  real              :: sfsm, rzsm, sfst, rzst
  logical           :: qc_sm, qc_st
  real, allocatable     :: sf_wt(:), rz_wt(:)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
  call LVT_verify(ios, 'Error opening file'//trim(filename))

!variable ids
  ios = nf90_inq_varid(nid, 'base_time',btimeid)
  call LVT_verify(ios, 'Error nf90_inq_varid: base_time')

  ios = nf90_inq_varid(nid, 'lat',latid)
  call LVT_verify(ios, 'Error nf90_inq_varid: lat')

  ios = nf90_inq_varid(nid, 'lon',lonid)
  call LVT_verify(ios, 'Error nf90_inq_varid: lon')

  ios = nf90_inq_varid(nid, 'time_offset',timeid)
  call LVT_verify(ios, 'Error nf90_inq_varid: time')

  ios = nf90_inq_varid(nid, 'depth',depthid)
  call LVT_verify(ios, 'Error nf90_inq_varid: depth')

  ios = nf90_inq_varid(nid, 'tsoil_E',tsoilEid)
  call LVT_verify(ios, 'Error nf90_inq_varid: tsoilE')

  ios = nf90_inq_varid(nid, 'qc_tsoil_E',qc_tsoilEid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_tsoilE')

  ios = nf90_inq_varid(nid, 'tsoil_W',tsoilWid)
  call LVT_verify(ios, 'Error nf90_inq_varid: tsoilW')

  ios = nf90_inq_varid(nid, 'qc_tsoil_W',qc_tsoilWid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_tsoilW')

  ios = nf90_inq_varid(nid, 'watcont_E',watcontEid)
  call LVT_verify(ios, 'Error nf90_inq_varid: watcont_E')

  ios = nf90_inq_varid(nid, 'qc_watcont_E',qc_watcontEid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_watcont_E')

  ios = nf90_inq_varid(nid, 'watcont_W',watcontWid)
  call LVT_verify(ios, 'Error nf90_inq_varid: watcont_W')

  ios = nf90_inq_varid(nid, 'qc_watcont_W',qc_watcontWid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_watcont_W')

!dimensions
  ios = nf90_inq_dimid(nid, 'time',dim1Id)
  call LVT_verify(ios, 'Error nf90_inq_dimid: time')

  ios = nf90_inquire_dimension(nid, dim1Id, len=tdims)
  call LVT_verify(ios, 'Error nf90_inquire_dimension:')

  ios = nf90_inq_dimid(nid, 'depth',dim2Id)
  call LVT_verify(ios, 'Error nf90_inq_dimid: depth')

  ios = nf90_inquire_dimension(nid, dim2Id, len=ddims)
  call LVT_verify(ios, 'Error nf90_inquire_dimension:')

!values
  ios = nf90_get_var(nid,latid, lat)
  call LVT_verify(ios, 'Error nf90_get_var: lat')

  ios = nf90_get_var(nid,lonid, lon)
  call LVT_verify(ios, 'Error nf90_get_var: lon')

  ios = nf90_get_var(nid,btimeid, base_time)
  call LVT_verify(ios, 'Error nf90_get_var: base_time')

  call ESMF_TimeSet(refTime, yy=1970,mm=1,dd=1,&
       h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
  call LVT_verify(status, 'error in timeset: readArmobs(Source)')

!  armobs(source)%stnlat(k) = lat
!  armobs(source)%stnlon(k) = lon

  allocate(time(tdims))
  allocate(depth(ddims))
  allocate(smE(ddims,tdims))
  allocate(stE(ddims,tdims))

  allocate(qc_smE(ddims,tdims))
  allocate(qc_stE(ddims,tdims))
  allocate(qc_smW(ddims,tdims))
  allocate(qc_stW(ddims,tdims))

  ios = nf90_get_var(nid,timeid, time)
  call LVT_verify(ios, 'Error nf90_get_var: time')

  ios = nf90_get_var(nid,depthid, depth)
  call LVT_verify(ios, 'Error nf90_get_var: depth')

  ios = nf90_get_var(nid,watcontEid, smE)
  call LVT_verify(ios, 'Error nf90_get_var: smE')

  ios = nf90_get_var(nid,qc_watcontEid, qc_smE)
  call LVT_verify(ios, 'Error nf90_get_var: qc_smE')

  ios = nf90_get_var(nid,watcontWid, smW)
  call LVT_verify(ios, 'Error nf90_get_var: smW')

  ios = nf90_get_var(nid,qc_watcontWid, qc_smW)
  call LVT_verify(ios, 'Error nf90_get_var: qc_smW')

  ios = nf90_get_var(nid,tsoilWid, stW)
  call LVT_verify(ios, 'Error nf90_get_var: stW')

  ios = nf90_get_var(nid,qc_tsoilWid, qc_stW)
  call LVT_verify(ios, 'Error nf90_get_var: qc_stW')

  ios = nf90_get_var(nid,tsoilEid, stE)
  call LVT_verify(ios, 'Error nf90_get_var: stE')

  ios = nf90_get_var(nid,qc_tsoilEid, qc_stE)
  call LVT_verify(ios, 'Error nf90_get_var: qc_stE')

  ios = nf90_close(nid)
  call LVT_verify(ios, 'Error in nf90_close')
  
  call ESMF_TimeIntervalSet(dt,s=base_time,rc=status)
  call LVT_verify(status, 'Error in timeintervalset: readARMobs')

  allocate(sf_wt(ddims))
  allocate(rz_wt(ddims)) 
  
  sf_wt = 0 
  rz_wt = 0 

  call compute_vinterp_weights(ddims,LVT_rc%lis_sf_d, &
       LVT_rc%lis_rz_d, depth(:)/100.0, sf_wt, rz_wt)

  reftime = reftime + dt

  do kk=1,tdims
     call ESMF_TimeIntervalSet(dt,s=nint(time(kk)),rc=status)
     call LVT_verify(status, 'Error in timeintervalset: readARMobs')

     datatime = reftime + dt
     
     call ESMF_TimeGet(datatime, yy=yr1,mm=mo1,dd=da1,&
          h=hr1,m=mn1,s=ss1,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readARMObs')

! data index from 0z with 30mn timestep. All this extensive processing is
! required because some of the ARM files report data at times different 
! from 0z. 
     
     call ESMF_TimeSet(currTime, yy=yr, mm=mo, &
          dd=da, h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readARMObs')

     data_index = (datatime - currtime)/armobs(source)%swats_ts + 1
     armobs(source)%swats_tindex(k,data_index) = data_index

     armobs(source)%swats_sfsm(k,data_index) = LVT_rc%udef
     armobs(source)%swats_rzsm(k,data_index) = LVT_rc%udef
     
     qc_sm = .true. 
     do t=1,ddims              
        qc_sm = qc_sm .and. (qc_smE(t,kk).eq.0)
     enddo
     
     if(qc_sm) then
        sfsm = 0 
        rzsm = 0 
        do t=1,ddims
           sfsm = sfsm + sf_wt(t)*smE(t,kk)
           rzsm = rzsm + rz_wt(t)*smE(t,kk)
        enddo
        
        armobs(source)%swats_sfsm(k,data_index)     = sfsm
        armobs(source)%swats_rzsm(k,data_index)     = rzsm
     else
        armobs(source)%swats_sfsm(k,data_index) = LVT_rc%udef
        armobs(source)%swats_rzsm(k,data_index) = LVT_rc%udef
     endif
     
     qc_st = .true.         
     do t=1,ddims              
        qc_st = qc_st .and. (qc_stE(t,kk).eq.0)
     enddo
     
     if(qc_st) then
        sfst = 0 
        rzst = 0 
        do t=1,ddims
           sfst = sfst + sf_wt(t)*stE(t,kk) 
           rzst = rzst + rz_wt(t)*stE(t,kk) 
        enddo
        armobs(source)%swats_sfst(k,data_index)     = sfst + 273.15
        armobs(source)%swats_rzst(k,data_index)     = rzst + 273.15
     else
        armobs(source)%swats_sfst(k,data_index) = LVT_rc%udef
        armobs(source)%swats_rzst(k,data_index) = LVT_rc%udef
     endif
  end do

  deallocate(sf_wt)
  deallocate(rz_wt)

  deallocate(time)
  deallocate(smE)
  deallocate(stE)

  deallocate(qc_smE)
  deallocate(qc_stE)
  deallocate(qc_smW)
  deallocate(qc_stW)
#endif

end subroutine read_swats_file


!BOP
! 
! !ROUTINE: read_baebbr_flux_file
! \label{read_baebbr_flux_file}
!
! !INTERFACE: 
subroutine read_baebbr_flux_file(source, k, yr, mo, da, filename)
! 
! !USES: 
  use ESMF     
  use LVT_coreMod,      only : LVT_rc
  use LVT_logMod,       only : LVT_logunit, LVT_verify
  use LVT_timeMgrMod,   only : LVT_calendar
  use ARM_obsMod,      only : armobs
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!DESCRIPTION: 
!  This routine reads the bulk aerodynamic estimates of 
!  latent, sensible and ground heat fluxes from the BAEBBR file, 
!  for a particular station. The data is read from the 
!  native NETCDF files, and is processed by applying the 
!  quality control flags. Only data points with 'good'
!  classification (qc value=0) is chosen. This routine 
!  also computes the temporal offset of the data relative
!  to the 0z of a particular day. 
!
!  The arguments are: 
!  \begin{description}
!   \item[k]         Index of the ARM station
!   \item[yr]        year of ARM data 
!   \item[mo]        month of ARM data
!   \item[da]        day of ARM data
!   \item[filename]  Name of the BAEBBR file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer,  intent(in) :: k 
  integer,  intent(in) :: source
  integer              :: yr
  integer              :: mo
  integer              :: da
  character(len=*)     :: filename 
  integer              :: bend
!EOP
  integer       :: nid, btimeid, latid, lonid, qleid
  integer       :: qhid,qgid, timeid, dimId
  integer       :: qc_qhid, qc_qleid, qc_qgid
  integer       :: netradid, solarradid, qc_netradid
  integer       :: ndims
  real          :: lat, lon 
  integer       :: base_time
  integer       :: ios
  real, allocatable :: time(:)
  real, allocatable :: baebbr_qh(:)
  real, allocatable :: baebbr_qle(:)
  real, allocatable :: baebbr_qg(:)
  real, allocatable :: baebbr_netrad(:)
  real, allocatable :: baebbr_solarrad(:)
  integer, allocatable :: baebbr_qc_qh(:)
  integer, allocatable :: baebbr_qc_qle(:)
  integer, allocatable :: baebbr_qc_qg(:)
  integer, allocatable :: baebbr_qc_netrad(:)
  integer           :: kk
  type(ESMF_Time)   :: reftime, datatime, currtime
  type(ESMF_TimeInterval) :: dt
  integer           :: yr1, mo1, da1, hr1, mn1, ss1
  integer           :: status
  integer           :: data_index

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
  call LVT_verify(ios, 'Error opening file'//trim(filename))

!variable ids
  ios = nf90_inq_varid(nid, 'base_time',btimeid)
  call LVT_verify(ios, 'Error nf90_inq_varid: base_time')

  ios = nf90_inq_varid(nid, 'lat',latid)
  call LVT_verify(ios, 'Error nf90_inq_varid: lat')

  ios = nf90_inq_varid(nid, 'lon',lonid)
  call LVT_verify(ios, 'Error nf90_inq_varid: lon')

  ios = nf90_inq_varid(nid, 'time_offset',timeid)
  call LVT_verify(ios, 'Error nf90_inq_varid: time')

  ios = nf90_inq_varid(nid, 'be_sensible_heat_flux',qhid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qh')

  ios = nf90_inq_varid(nid, 'qc_be_sensible_heat_flux',qc_qhid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_qh')

  ios = nf90_inq_varid(nid, 'be_latent_heat_flux',qleid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qle')

  ios = nf90_inq_varid(nid, 'qc_be_latent_heat_flux',qc_qleid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_qle')

  ios = nf90_inq_varid(nid, 'surface_soil_heat_flux_avg',qgid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qg')

  ios = nf90_inq_varid(nid, 'qc_surface_soil_heat_flux_avg',qc_qgid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_qg')

  ios = nf90_inq_varid(nid, 'net_radiation',netradid)
  call LVT_verify(ios,'Error in nf90_inq_varid: net_radiation')

  ios = nf90_inq_varid(nid,'qc_net_radiation',qc_netradid)
  call LVT_verify(ios,'Error in nf90_inq_varid: qc_net_radiation')

  ios = nf90_inq_varid(nid, 'solar_radiation',solarradid)
  call LVT_verify(ios,'Error in nf90_inq_varid: solar_radiation')

!dimensions
  ios = nf90_inq_dimid(nid, 'time',dimId)
  call LVT_verify(ios, 'Error nf90_inq_dimid: time')

  ios = nf90_inquire_dimension(nid, dimId, len=ndims)
  call LVT_verify(ios, 'Error nf90_inquire_dimension:')

!values
  ios = nf90_get_var(nid,latid, lat)
  call LVT_verify(ios, 'Error nf90_get_var: lat')

  ios = nf90_get_var(nid,lonid, lon)
  call LVT_verify(ios, 'Error nf90_get_var: lon')

  ios = nf90_get_var(nid,btimeid, base_time)
  call LVT_verify(ios, 'Error nf90_get_var: base_time')

  call ESMF_TimeSet(refTime, yy=1970,mm=1,dd=1,&
       h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
  call LVT_verify(status, 'error in timeset: readARMObs')

!  armobs%stnlat(k) = lat
!  armobs%stnlon(k) = lon

  allocate(time(ndims))
  allocate(baebbr_qle(ndims))
  allocate(baebbr_qh(ndims))
  allocate(baebbr_qg(ndims))
  allocate(baebbr_netrad(ndims))
  allocate(baebbr_solarrad(ndims))

  allocate(baebbr_qc_qle(ndims))
  allocate(baebbr_qc_qh(ndims))
  allocate(baebbr_qc_qg(ndims))
  allocate(baebbr_qc_netrad(ndims))

  ios = nf90_get_var(nid,timeid, time)
  call LVT_verify(ios, 'Error nf90_get_var: time')

  ios = nf90_get_var(nid,qhid, baebbr_qh)
  call LVT_verify(ios, 'Error nf90_get_var: qh')

  ios = nf90_get_var(nid,qc_qhid, baebbr_qc_qh)
  call LVT_verify(ios, 'Error nf90_get_var: qc_qh')

  ios = nf90_get_var(nid,qleid, baebbr_qle)
  call LVT_verify(ios, 'Error nf90_get_var: qle')

  ios = nf90_get_var(nid,qc_qleid, baebbr_qc_qle)
  call LVT_verify(ios, 'Error nf90_get_var: qc_qle')

  ios = nf90_get_var(nid,qgid, baebbr_qg)
  call LVT_verify(ios, 'Error nf90_get_var: qg')

  ios = nf90_get_var(nid,qc_qgid, baebbr_qc_qg)
  call LVT_verify(ios, 'Error nf90_get_var: qc_qg')

  ios = nf90_get_var(nid,netradid, baebbr_netrad)
  call LVT_verify(ios, 'Error nf90_get_var: net_radiation')

  ios = nf90_get_var(nid,qc_netradid, baebbr_qc_netrad)
  call LVT_verify(ios, 'Error nf90_get_var: qc_net_radiation')

  ios = nf90_get_var(nid,solarradid, baebbr_solarrad)
  call LVT_verify(ios, 'Error nf90_get_var: solar_radiation')

  ios = nf90_close(nid)
  call LVT_verify(ios, 'Error in nf90_close')
  
  call ESMF_TimeIntervalSet(dt,s=base_time,rc=status)
  call LVT_verify(status, 'Error in timeintervalset: readARMobs')

  reftime = reftime + dt

  do kk=1,ndims
     call ESMF_TimeIntervalSet(dt,s=nint(time(kk)),rc=status)
     call LVT_verify(status, 'Error in timeintervalset: readARMobs')

     datatime = reftime + dt
     
     call ESMF_TimeGet(datatime, yy=yr1,mm=mo1,dd=da1,&
          h=hr1,m=mn1,s=ss1,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readARMObs')

! data index from 0z with 30mn timestep. All this extensive processing is
! required because some of the ARM files report data at times different 
! from 0z. 
     
     call ESMF_TimeSet(currTime, yy=yr, mm=mo, &
          dd=da, h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readARMObs')

     data_index = (datatime - currtime)/armobs(source)%baebbr_ts + 1
     armobs(source)%baebbr_tindex(k,data_index) = data_index

     if(baebbr_qc_qh(kk).eq.0) then  
        armobs(source)%baebbr_qh(k,data_index)     = baebbr_qh(kk)
     else
        armobs(source)%baebbr_qh(k,data_index)     = LVT_rc%udef
     endif
     if(baebbr_qc_qle(kk).eq.0) then 
        armobs(source)%baebbr_qle(k,data_index)    = baebbr_qle(kk)
     else
        armobs(source)%baebbr_qle(k,data_index)    = LVT_rc%udef
     endif
     
     if(baebbr_qc_qg(kk).eq.0) then 
        armobs(source)%baebbr_qg(k,data_index)     = baebbr_qg(kk)
     else
        armobs(source)%baebbr_qg(k,data_index)     = LVT_rc%udef
     endif
     
     if(baebbr_qc_netrad(kk).eq.0) then 
        armobs(source)%baebbr_netrad(k,data_index)     = baebbr_netrad(kk)
     else
        armobs(source)%baebbr_netrad(k,data_index)     = LVT_rc%udef
     endif
     armobs(source)%baebbr_solarrad(k,data_index)     = baebbr_solarrad(kk)
  enddo

  deallocate(time)
  deallocate(baebbr_qle)
  deallocate(baebbr_qh)
  deallocate(baebbr_qg)
  deallocate(baebbr_netrad)
  deallocate(baebbr_solarrad)

  deallocate(baebbr_qc_qle)
  deallocate(baebbr_qc_qh)
  deallocate(baebbr_qc_qg)
  deallocate(baebbr_qc_netrad)
#endif

end subroutine read_baebbr_flux_file

!BOP
! 
! !ROUTINE: read_ebbr_file
! \label{read_ebbr_file}
!
! !INTERFACE: 
subroutine read_ebbr_file(source, k, yr, mo, da, filename)
! 
! !USES: 
  use ESMF     
  use LVT_coreMod,      only : LVT_rc
  use LVT_logMod,       only : LVT_logunit, LVT_verify
  use LVT_timeMgrMod,   only : LVT_calendar
  use ARM_obsMod,      only : armobs
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!DESCRIPTION: 
! 
!  The arguments are: 
!  \begin{description}
!   \item[k]         Index of the ARM station
!   \item[yr]        year of ARM data 
!   \item[mo]        month of ARM data
!   \item[da]        day of ARM data
!   \item[filename]  Name of the BAEBBR file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer,  intent(in) :: k 
  integer,  intent(in) :: source
  integer              :: yr
  integer              :: mo
  integer              :: da
  character(len=*)     :: filename 

!EOP
  integer       :: nid, btimeid, latid, lonid
  integer       :: timeid, dimId
  integer       :: sm1id,sm2id,sm3id,sm4id,sm5id
  integer       :: qc_sm1id,qc_sm2id,qc_sm3id,qc_sm4id,qc_sm5id
  integer       :: st1id,st2id,st3id,st4id,st5id
  integer       :: qc_st1id,qc_st2id,qc_st3id,qc_st4id,qc_st5id
  integer       :: ndims
  real          :: lat, lon 
  integer       :: base_time
  integer       :: ios
  real, allocatable :: time(:)
  real, allocatable :: ebbr_sm1(:),ebbr_st1(:)
  real, allocatable :: ebbr_sm2(:),ebbr_st2(:)
  real, allocatable :: ebbr_sm3(:),ebbr_st3(:)
  real, allocatable :: ebbr_sm4(:),ebbr_st4(:)
  real, allocatable :: ebbr_sm5(:),ebbr_st5(:)
  
  integer, allocatable :: ebbr_qc_sm1(:),ebbr_qc_st1(:)
  integer, allocatable :: ebbr_qc_sm2(:),ebbr_qc_st2(:)
  integer, allocatable :: ebbr_qc_sm3(:),ebbr_qc_st3(:)
  integer, allocatable :: ebbr_qc_sm4(:),ebbr_qc_st4(:)
  integer, allocatable :: ebbr_qc_sm5(:),ebbr_qc_st5(:)
  logical           :: qc_sm,qc_st
  integer           :: kk
  type(ESMF_Time)   :: reftime, datatime, currtime
  type(ESMF_TimeInterval) :: dt
  real              :: sfsm,sfst
  integer           :: yr1, mo1, da1, hr1, mn1, ss1
  integer           :: status
  integer           :: data_index

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
  call LVT_verify(ios, 'Error opening file'//trim(filename))

!variable ids
  ios = nf90_inq_varid(nid, 'base_time',btimeid)
  call LVT_verify(ios, 'Error nf90_inq_varid: base_time')

  ios = nf90_inq_varid(nid, 'lat',latid)
  call LVT_verify(ios, 'Error nf90_inq_varid: lat')

  ios = nf90_inq_varid(nid, 'lon',lonid)
  call LVT_verify(ios, 'Error nf90_inq_varid: lon')

  ios = nf90_inq_varid(nid, 'time_offset',timeid)
  call LVT_verify(ios, 'Error nf90_inq_varid: time')

  ios = nf90_inq_varid(nid, 'soil_moisture_1',sm1id)
  call LVT_verify(ios, 'Error nf90_inq_varid: sm1')

  ios = nf90_inq_varid(nid, 'qc_soil_moisture_1',qc_sm1id)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_sm1')

  ios = nf90_inq_varid(nid, 'soil_moisture_2',sm2id)
  call LVT_verify(ios, 'Error nf90_inq_varid: sm2')

  ios = nf90_inq_varid(nid, 'qc_soil_moisture_2',qc_sm2id)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_sm2')

  ios = nf90_inq_varid(nid, 'soil_moisture_3',sm3id)
  call LVT_verify(ios, 'Error nf90_inq_varid: sm3')

  ios = nf90_inq_varid(nid, 'qc_soil_moisture_3',qc_sm3id)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_sm3')

  ios = nf90_inq_varid(nid, 'soil_moisture_4',sm4id)
  call LVT_verify(ios, 'Error nf90_inq_varid: sm4')

  ios = nf90_inq_varid(nid, 'qc_soil_moisture_4',qc_sm4id)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_sm4')

  ios = nf90_inq_varid(nid, 'soil_moisture_5',sm5id)
  call LVT_verify(ios, 'Error nf90_inq_varid: sm5')

  ios = nf90_inq_varid(nid, 'qc_soil_moisture_5',qc_sm5id)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_sm5')

  ios = nf90_inq_varid(nid, 'soil_temp_1',st1id)
  call LVT_verify(ios, 'Error nf90_inq_varid: st1')

  ios = nf90_inq_varid(nid, 'qc_soil_temp_1',qc_st1id)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_st1')

  ios = nf90_inq_varid(nid, 'soil_temp_2',st2id)
  call LVT_verify(ios, 'Error nf90_inq_varid: st2')

  ios = nf90_inq_varid(nid, 'qc_soil_temp_2',qc_st2id)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_st2')

  ios = nf90_inq_varid(nid, 'soil_temp_3',st3id)
  call LVT_verify(ios, 'Error nf90_inq_varid: st3')

  ios = nf90_inq_varid(nid, 'qc_soil_temp_3',qc_st3id)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_st3')

  ios = nf90_inq_varid(nid, 'soil_temp_4',st4id)
  call LVT_verify(ios, 'Error nf90_inq_varid: st4')

  ios = nf90_inq_varid(nid, 'qc_soil_temp_4',qc_st4id)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_st4')

  ios = nf90_inq_varid(nid, 'soil_temp_5',st5id)
  call LVT_verify(ios, 'Error nf90_inq_varid: st5')

  ios = nf90_inq_varid(nid, 'qc_soil_temp_5',qc_st5id)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_st5')
!dimensions
  ios = nf90_inq_dimid(nid, 'time',dimId)
  call LVT_verify(ios, 'Error nf90_inq_dimid: time')

  ios = nf90_inquire_dimension(nid, dimId, len=ndims)
  call LVT_verify(ios, 'Error nf90_inquire_dimension:')

!values
  ios = nf90_get_var(nid,latid, lat)
  call LVT_verify(ios, 'Error nf90_get_var: lat')

  ios = nf90_get_var(nid,lonid, lon)
  call LVT_verify(ios, 'Error nf90_get_var: lon')

  ios = nf90_get_var(nid,btimeid, base_time)
  call LVT_verify(ios, 'Error nf90_get_var: base_time')

  call ESMF_TimeSet(refTime, yy=1970,mm=1,dd=1,&
       h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
  call LVT_verify(status, 'error in timeset: readARMObs')

!  armobs(source)%stnlat(k) = lat
!  armobs(source)%stnlon(k) = lon

  allocate(time(ndims))
  allocate(ebbr_sm1(ndims))
  allocate(ebbr_sm2(ndims))
  allocate(ebbr_sm3(ndims))
  allocate(ebbr_sm4(ndims))
  allocate(ebbr_sm5(ndims))
  allocate(ebbr_st1(ndims))
  allocate(ebbr_st2(ndims))
  allocate(ebbr_st3(ndims))
  allocate(ebbr_st4(ndims))
  allocate(ebbr_st5(ndims))

  allocate(ebbr_qc_sm1(ndims))
  allocate(ebbr_qc_sm2(ndims))
  allocate(ebbr_qc_sm3(ndims))
  allocate(ebbr_qc_sm4(ndims))
  allocate(ebbr_qc_sm5(ndims))
  allocate(ebbr_qc_st1(ndims))
  allocate(ebbr_qc_st2(ndims))
  allocate(ebbr_qc_st3(ndims))
  allocate(ebbr_qc_st4(ndims))
  allocate(ebbr_qc_st5(ndims))

  ios = nf90_get_var(nid,timeid, time)
  call LVT_verify(ios, 'Error nf90_get_var: time')

  ios = nf90_get_var(nid,sm1id, ebbr_sm1)
  call LVT_verify(ios, 'Error nf90_get_var: sm1')

  ios = nf90_get_var(nid,qc_sm1id, ebbr_qc_sm1)
  call LVT_verify(ios, 'Error nf90_get_var: qc_sm1')

  ios = nf90_get_var(nid,sm2id, ebbr_sm2)
  call LVT_verify(ios, 'Error nf90_get_var: sm2')

  ios = nf90_get_var(nid,qc_sm2id, ebbr_qc_sm2)
  call LVT_verify(ios, 'Error nf90_get_var: qc_sm2')

  ios = nf90_get_var(nid,sm3id, ebbr_sm3)
  call LVT_verify(ios, 'Error nf90_get_var: sm3')

  ios = nf90_get_var(nid,qc_sm3id, ebbr_qc_sm3)
  call LVT_verify(ios, 'Error nf90_get_var: qc_sm3')

  ios = nf90_get_var(nid,sm4id, ebbr_sm4)
  call LVT_verify(ios, 'Error nf90_get_var: sm4')

  ios = nf90_get_var(nid,qc_sm4id, ebbr_qc_sm4)
  call LVT_verify(ios, 'Error nf90_get_var: qc_sm4')

  ios = nf90_get_var(nid,sm5id, ebbr_sm5)
  call LVT_verify(ios, 'Error nf90_get_var: sm5')

  ios = nf90_get_var(nid,qc_sm5id, ebbr_qc_sm5)
  call LVT_verify(ios, 'Error nf90_get_var: qc_sm5')

  ios = nf90_get_var(nid,st1id, ebbr_st1)
  call LVT_verify(ios, 'Error nf90_get_var: st1')

  ios = nf90_get_var(nid,qc_st1id, ebbr_qc_st1)
  call LVT_verify(ios, 'Error nf90_get_var: qc_st1')

  ios = nf90_get_var(nid,st2id, ebbr_st2)
  call LVT_verify(ios, 'Error nf90_get_var: st2')

  ios = nf90_get_var(nid,qc_st2id, ebbr_qc_st2)
  call LVT_verify(ios, 'Error nf90_get_var: qc_st2')

  ios = nf90_get_var(nid,st3id, ebbr_st3)
  call LVT_verify(ios, 'Error nf90_get_var: st3')

  ios = nf90_get_var(nid,qc_st3id, ebbr_qc_st3)
  call LVT_verify(ios, 'Error nf90_get_var: qc_st3')

  ios = nf90_get_var(nid,st4id, ebbr_st4)
  call LVT_verify(ios, 'Error nf90_get_var: st4')

  ios = nf90_get_var(nid,qc_st4id, ebbr_qc_st4)
  call LVT_verify(ios, 'Error nf90_get_var: qc_st4')

  ios = nf90_get_var(nid,st5id, ebbr_st5)
  call LVT_verify(ios, 'Error nf90_get_var: st5')

  ios = nf90_get_var(nid,qc_st5id, ebbr_qc_st5)
  call LVT_verify(ios, 'Error nf90_get_var: qc_st5')

  ios = nf90_close(nid)
  call LVT_verify(ios, 'Error in nf90_close')
  
  call ESMF_TimeIntervalSet(dt,s=base_time,rc=status)
  call LVT_verify(status, 'Error in timeintervalset: readARMobs')

  reftime = reftime + dt

  do kk=1,ndims
     call ESMF_TimeIntervalSet(dt,s=nint(time(kk)),rc=status)
     call LVT_verify(status, 'Error in timeintervalset: readARMobs')

     datatime = reftime + dt
     
     call ESMF_TimeGet(datatime, yy=yr1,mm=mo1,dd=da1,&
          h=hr1,m=mn1,s=ss1,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readARMObs')

! data index from 0z with 30mn timestep. All this extensive processing is
! required because some of the ARM files report data at times different 
! from 0z. 
     
     call ESMF_TimeSet(currTime, yy=yr, mm=mo, &
          dd=da, h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readARMObs')

     data_index = (datatime - currtime)/armobs(source)%ebbr_ts + 1
     armobs(source)%ebbr_tindex(k,data_index) = data_index
        
     sfsm = LVT_rc%udef
     sfst = LVT_rc%udef
     qc_sm = .true. 
     qc_sm = qc_sm .and. &
             (ebbr_qc_sm1(kk).eq.0).and.&
             (ebbr_qc_sm2(kk).eq.0).and.&
             (ebbr_qc_sm3(kk).eq.0).and.&
             (ebbr_qc_sm4(kk).eq.0).and.&
             (ebbr_qc_sm5(kk).eq.0)
     
     if(qc_sm) then
        sfsm = 0 
        if(armobs(source)%stnbd(k).ne.-9999.0) then 
           sfsm = sfsm + &
                (((ebbr_sm1(kk) + &
                ebbr_sm2(kk) + &
                ebbr_sm3(kk) + &
                ebbr_sm4(kk) + &
                ebbr_sm5(kk))/5.0)*armobs(source)%stnbd(k))/100.0
        else
           sfsm = LVT_rc%udef
        endif
     endif
     
     qc_st = .true. 
     qc_st = qc_st .and. &
          (ebbr_qc_st1(kk).eq.0).and.&
          (ebbr_qc_st2(kk).eq.0).and.&
          (ebbr_qc_st3(kk).eq.0).and.&
          (ebbr_qc_st4(kk).eq.0).and.&
          (ebbr_qc_st5(kk).eq.0).and.&
          (ebbr_st1(kk).gt.-100).and.&
          (ebbr_st2(kk).gt.-100).and.&
          (ebbr_st3(kk).gt.-100).and.&
          (ebbr_st4(kk).gt.-100).and.&
          (ebbr_st5(kk).gt.-100).and.&
          (ebbr_st1(kk).lt.100).and.&
          (ebbr_st2(kk).lt.100).and.&
          (ebbr_st3(kk).lt.100).and.&
          (ebbr_st4(kk).lt.100).and.&
          (ebbr_st5(kk).lt.100)
     if(qc_st) then
        sfst = 0 
        sfst = sfst + &
             ((ebbr_st1(kk) + & 
             ebbr_st2(kk) + & 
             ebbr_st3(kk) + & 
             ebbr_st4(kk) + & 
             ebbr_st5(kk))/5)+273.15
        if(sfst.lt.0) then 
           print*, ebbr_st1(kk),ebbr_qc_st1(kk),ebbr_st2(kk),&
                ebbr_st3(kk),ebbr_st4(kk),&
                ebbr_st5(kk)
           stop
        endif
     else
        sfst = LVT_rc%udef
     endif
     
     armobs(source)%ebbr_sfsm(k,data_index)     = sfsm
     armobs(source)%ebbr_sfst(k,data_index)     = sfst
     
  enddo

  deallocate(time)
  deallocate(ebbr_sm1)
  deallocate(ebbr_sm2)
  deallocate(ebbr_sm3)
  deallocate(ebbr_sm4)
  deallocate(ebbr_sm5)
  deallocate(ebbr_st1)
  deallocate(ebbr_st2)
  deallocate(ebbr_st3)
  deallocate(ebbr_st4)
  deallocate(ebbr_st5)

  deallocate(ebbr_qc_sm1)
  deallocate(ebbr_qc_sm2)
  deallocate(ebbr_qc_sm3)
  deallocate(ebbr_qc_sm4)
  deallocate(ebbr_qc_sm5)
  deallocate(ebbr_qc_st1)
  deallocate(ebbr_qc_st2)
  deallocate(ebbr_qc_st3)
  deallocate(ebbr_qc_st4)
  deallocate(ebbr_qc_st5)

#endif

end subroutine read_ebbr_file


!BOP
! 
! !ROUTINE: read_ecor_flux_file
! \label{read_ecor_flux_file}
!
! !INTERFACE: 
subroutine read_ecor_flux_file(source, k, yr, mo, da, filename)
! 
! !USES: 
  use ESMF     
  use LVT_coreMod,      only : LVT_rc
  use LVT_logMod,       only : LVT_logunit, LVT_verify
  use LVT_timeMgrMod,   only : LVT_calendar
  use ARM_obsMod,      only : armobs
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!DESCRIPTION: 
!  This routine reads the eddy correlation estimates of 
!  latent and sensible fluxes from the ECOR file, 
!  for a particular station. The data is read from the 
!  native NETCDF files, and is processed by applying the 
!  quality control flags. Only data points with 'good'
!  classification (qc value=0) is chosen. This routine 
!  also computes the temporal offset of the data relative
!  to the 0z of a particular day. 
!
!  The arguments are: 
!  \begin{description}
!   \item[k]         Index of the ARM station
!   \item[yr]        year of ARM data 
!   \item[mo]        month of ARM data
!   \item[da]        day of ARM data
!   \item[filename]  Name of the ECOR file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer,  intent(in) :: k 
  integer,  intent(in) :: source
  character(len=*)  :: filename 
  integer           :: yr, mo, da
  integer          :: bend
!EOP
  integer       :: nid, btimeid, latid, lonid, qleid
  integer       :: qhid,timeid, dimId
  integer       :: qc_qhid, qc_qleid
  integer       :: ndims
  real          :: lat, lon 
  integer       :: base_time
  integer       :: ios
  real, allocatable :: time(:)
  real, allocatable :: ecor_qh(:)
  real, allocatable :: ecor_qle(:)
  integer, allocatable :: ecor_qc_qh(:)
  integer, allocatable :: ecor_qc_qle(:)
  integer           :: kk
  type(ESMF_Time)   :: reftime, datatime, currtime
  type(ESMF_TimeInterval) :: dt
  integer           :: yr1, mo1, da1, hr1, mn1, ss1
  integer           :: status
  integer           :: data_index

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
  call LVT_verify(ios, 'Error opening file'//trim(filename))

!variable ids
  ios = nf90_inq_varid(nid, 'base_time',btimeid)
  call LVT_verify(ios, 'Error nf90_inq_varid: base_time')

  ios = nf90_inq_varid(nid, 'lat',latid)
  call LVT_verify(ios, 'Error nf90_inq_varid: lat')

  ios = nf90_inq_varid(nid, 'lon',lonid)
  call LVT_verify(ios, 'Error nf90_inq_varid: lon')

  ios = nf90_inq_varid(nid, 'time_offset',timeid)
  call LVT_verify(ios, 'Error nf90_inq_varid: time')

  ios = nf90_inq_varid(nid, 'h',qhid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qh')

  ios = nf90_inq_varid(nid, 'qc_h',qc_qhid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_qh')

  ios = nf90_inq_varid(nid, 'lv_e',qleid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qle')

  ios = nf90_inq_varid(nid, 'qc_lv_e',qc_qleid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_qle')

!dimensions
  ios = nf90_inq_dimid(nid, 'time',dimId)
  call LVT_verify(ios, 'Error nf90_inq_dimid: time')

  ios = nf90_inquire_dimension(nid, dimId, len=ndims)
  call LVT_verify(ios, 'Error nf90_inquire_dimension:')

!values
  ios = nf90_get_var(nid,latid, lat)
  call LVT_verify(ios, 'Error nf90_get_var: lat')

  ios = nf90_get_var(nid,lonid, lon)
  call LVT_verify(ios, 'Error nf90_get_var: lon')

  ios = nf90_get_var(nid,btimeid, base_time)
  call LVT_verify(ios, 'Error nf90_get_var: base_time')

  call ESMF_TimeSet(refTime, yy=1970,mm=1,dd=1,&
       h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
  call LVT_verify(status, 'error in timeset: readARMObs')

!  armobs(source)%stnlat(k) = lat
!  armobs(source)%stnlon(k) = lon

  allocate(time(ndims))
  allocate(ecor_qle(ndims))
  allocate(ecor_qh(ndims))

  allocate(ecor_qc_qle(ndims))
  allocate(ecor_qc_qh(ndims))

  ios = nf90_get_var(nid,timeid, time)
  call LVT_verify(ios, 'Error nf90_get_var: time')

  ios = nf90_get_var(nid,qhid, ecor_qh)
  call LVT_verify(ios, 'Error nf90_get_var: qh')

  ios = nf90_get_var(nid,qc_qhid, ecor_qc_qh)
  call LVT_verify(ios, 'Error nf90_get_var: qc_qh')

  ios = nf90_get_var(nid,qleid, ecor_qle)
  call LVT_verify(ios, 'Error nf90_get_var: qle')

  ios = nf90_get_var(nid,qc_qleid, ecor_qc_qle)
  call LVT_verify(ios, 'Error nf90_get_var: qc_qle')

  ios = nf90_close(nid)
  call LVT_verify(ios, 'Error in nf90_close')
  
  call ESMF_TimeIntervalSet(dt,s=base_time,rc=status)
  call LVT_verify(status, 'Error in timeintervalset: readARMobs')

  reftime = reftime + dt

  do kk=1,ndims
     call ESMF_TimeIntervalSet(dt,s=nint(time(kk)),rc=status)
     call LVT_verify(status, 'Error in timeintervalset: readARMobs')

     datatime = reftime + dt
     
     call ESMF_TimeGet(datatime, yy=yr1,mm=mo1,dd=da1,&
          h=hr1,m=mn1,s=ss1,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readARMObs')

! data index from 0z with 30mn timestep. All this extensive processing is
! required because some of the ARM files report data at times different 
! from 0z. 
     
     call ESMF_TimeSet(currTime, yy=yr, mm=mo, &
          dd=da, h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readARMObs')

     data_index = (datatime - currtime)/armobs(source)%ecor_ts + 1
     armobs(source)%ecor_tindex(k,data_index) = data_index

     if(ecor_qc_qh(kk).eq.0) then 
        armobs(source)%ecor_qh(k,data_index)     = ecor_qh(kk)
     else
        armobs(source)%ecor_qh(k,data_index)     = LVT_rc%udef
     endif
     if(ecor_qc_qle(kk).eq.0) then 
        armobs(source)%ecor_qle(k,data_index)    = ecor_qle(kk)
     else
        armobs(source)%ecor_qle(k,data_index)    = LVT_rc%udef
     endif
  enddo

  deallocate(time)
  deallocate(ecor_qle)
  deallocate(ecor_qh)
  deallocate(ecor_qc_qle)
  deallocate(ecor_qc_qh)

#endif

end subroutine read_ecor_flux_file

!BOP
! 
! !ROUTINE: read_smos_file
! \label{read_smos_file}
!
! !INTERFACE: 
subroutine read_smos_file(source, k, yr, mo, da, filename)
! 
! !USES: 
  use ESMF     
  use LVT_coreMod,      only : LVT_rc
  use LVT_logMod,       only : LVT_logunit, LVT_verify
  use LVT_timeMgrMod,   only : LVT_calendar
  use ARM_obsMod,      only : armobs
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!DESCRIPTION: 
!  This routine reads the surface meteorology data
!  from the SMOS file, for a particular station. The data is read from the 
!  native NETCDF files, and is processed by applying the 
!  quality control flags. Only data points with 'good'
!  classification (qc value=0) is chosen. This routine 
!  also computes the temporal offset of the data relative
!  to the 0z of a particular day. 
!
!  The arguments are: 
!  \begin{description}
!   \item[k]         Index of the ARM station
!   \item[yr]        year of ARM data 
!   \item[mo]        month of ARM data
!   \item[da]        day of ARM data
!   \item[filename]  Name of the ECOR file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer,  intent(in) :: k 
  integer,  intent(in) :: source
  character(len=*)  :: filename 
  integer           :: yr, mo, da
!EOP
  integer       :: nid, btimeid, latid, lonid, wspdid
  integer       :: pcpid,timeid, dimId,tempid,prsid
  integer       :: snowid,rhid
  integer       :: qc_pcpid, qc_wspdid,qc_tempid
  integer       :: qc_prsid, qc_snowid,qc_rhid
  integer       :: ndims
  real          :: lat, lon 
  integer       :: base_time
  integer       :: ios
  real, allocatable :: time(:)
  real, allocatable :: smos_pcp(:)
  real, allocatable :: smos_wspd(:)
  real, allocatable :: smos_rh(:)
  real, allocatable :: smos_temp(:)
  real, allocatable :: smos_prs(:)
  real, allocatable :: smos_snow(:)
  integer, allocatable :: smos_qc_pcp(:)
  integer, allocatable :: smos_qc_rh(:)
  integer, allocatable :: smos_qc_wspd(:)
  integer, allocatable :: smos_qc_temp(:)
  integer, allocatable :: smos_qc_prs(:)
  integer, allocatable :: smos_qc_snow(:)
  integer           :: kk
  type(ESMF_Time)   :: reftime, datatime, currtime
  type(ESMF_TimeInterval) :: dt
  integer           :: yr1, mo1, da1, hr1, mn1, ss1
  integer           :: status
  integer           :: data_index

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
  call LVT_verify(ios, 'Error opening file'//trim(filename))

!variable ids
  ios = nf90_inq_varid(nid, 'base_time',btimeid)
  call LVT_verify(ios, 'Error nf90_inq_varid: base_time')

  ios = nf90_inq_varid(nid, 'lat',latid)
  call LVT_verify(ios, 'Error nf90_inq_varid: lat')

  ios = nf90_inq_varid(nid, 'lon',lonid)
  call LVT_verify(ios, 'Error nf90_inq_varid: lon')

  ios = nf90_inq_varid(nid, 'time_offset',timeid)
  call LVT_verify(ios, 'Error nf90_inq_varid: time')

  ios = nf90_inq_varid(nid, 'precip',pcpid)
  call LVT_verify(ios, 'Error nf90_inq_varid: pcp')

  ios = nf90_inq_varid(nid, 'qc_precip',qc_pcpid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_pcp')

  ios = nf90_inq_varid(nid, 'wspd',wspdid)
  call LVT_verify(ios, 'Error nf90_inq_varid: wspd')

  ios = nf90_inq_varid(nid, 'qc_wspd',qc_wspdid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_wspd')

  ios = nf90_inq_varid(nid, 'temp',tempid)
  call LVT_verify(ios, 'Error nf90_inq_varid: temp')

  ios = nf90_inq_varid(nid, 'qc_temp',qc_tempid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_temp')

  ios = nf90_inq_varid(nid, 'rh',rhid)
  call LVT_verify(ios, 'Error nf90_inq_varid: rh')

  ios = nf90_inq_varid(nid, 'qc_rh',qc_rhid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_rh')

  ios = nf90_inq_varid(nid, 'bar_pres',prsid)
  call LVT_verify(ios, 'Error nf90_inq_varid: bar_pres')

  ios = nf90_inq_varid(nid, 'qc_bar_pres',qc_prsid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qc_prs')

  ios = nf90_inq_varid(nid, 'snow',snowid)
  call LVT_verify(ios, 'Error nf90_inq_varid: snow')

  ios = nf90_inq_varid(nid, 'qc_snow',qc_snowid)
  call LVT_verify(ios, 'Error nf90_inq_varid: snow')

!dimensions
  ios = nf90_inq_dimid(nid, 'time',dimId)
  call LVT_verify(ios, 'Error nf90_inq_dimid: time')

  ios = nf90_inquire_dimension(nid, dimId, len=ndims)
  call LVT_verify(ios, 'Error nf90_inquire_dimension:')

!values
  ios = nf90_get_var(nid,latid, lat)
  call LVT_verify(ios, 'Error nf90_get_var: lat')

  ios = nf90_get_var(nid,lonid, lon)
  call LVT_verify(ios, 'Error nf90_get_var: lon')

  ios = nf90_get_var(nid,btimeid, base_time)
  call LVT_verify(ios, 'Error nf90_get_var: base_time')

  call ESMF_TimeSet(refTime, yy=1970,mm=1,dd=1,&
       h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
  call LVT_verify(status, 'error in timeset: readARMObs')

!  armobs(source)%stnlat(k) = lat
!  armobs(source)%stnlon(k) = lon

  allocate(time(ndims))
  allocate(smos_wspd(ndims))
  allocate(smos_temp(ndims))
  allocate(smos_prs(ndims))
  allocate(smos_pcp(ndims))
  allocate(smos_rh(ndims))
  allocate(smos_snow(ndims))

  allocate(smos_qc_wspd(ndims))
  allocate(smos_qc_temp(ndims))
  allocate(smos_qc_prs(ndims))
  allocate(smos_qc_pcp(ndims))
  allocate(smos_qc_rh(ndims))
  allocate(smos_qc_snow(ndims))

  ios = nf90_get_var(nid,timeid, time)
  call LVT_verify(ios, 'Error nf90_get_var: time')

  ios = nf90_get_var(nid,pcpid, smos_pcp)
  call LVT_verify(ios, 'Error nf90_get_var: precip')

  ios = nf90_get_var(nid,qc_pcpid, smos_qc_pcp)
  call LVT_verify(ios, 'Error nf90_get_var: qc_precip')

  ios = nf90_get_var(nid,wspdid, smos_wspd)
  call LVT_verify(ios, 'Error nf90_get_var: wspd')

  ios = nf90_get_var(nid,qc_wspdid, smos_qc_wspd)
  call LVT_verify(ios, 'Error nf90_get_var: qc_wspd')

  ios = nf90_get_var(nid,tempid, smos_temp)
  call LVT_verify(ios, 'Error nf90_get_var: temp')

  ios = nf90_get_var(nid,qc_tempid, smos_qc_temp)
  call LVT_verify(ios, 'Error nf90_get_var: qc_temp')

  ios = nf90_get_var(nid,rhid, smos_rh)
  call LVT_verify(ios, 'Error nf90_get_var: rh')

  ios = nf90_get_var(nid,qc_rhid, smos_qc_rh)
  call LVT_verify(ios, 'Error nf90_get_var: qc_rh')

  ios = nf90_get_var(nid,prsid, smos_prs)
  call LVT_verify(ios, 'Error nf90_get_var: prs')

  ios = nf90_get_var(nid,qc_prsid, smos_qc_prs)
  call LVT_verify(ios, 'Error nf90_get_var: qc_prs')

  ios = nf90_get_var(nid,qc_snowid, smos_snow)
  call LVT_verify(ios, 'Error nf90_get_var: qc_snow')

  ios = nf90_close(nid)
  call LVT_verify(ios, 'Error in nf90_close')
  
  call ESMF_TimeIntervalSet(dt,s=base_time,rc=status)
  call LVT_verify(status, 'Error in timeintervalset: readARMobs')

  reftime = reftime + dt

  do kk=1,ndims
     call ESMF_TimeIntervalSet(dt,s=nint(time(kk)),rc=status)
     call LVT_verify(status, 'Error in timeintervalset: readARMobs')

     datatime = reftime + dt
     
     call ESMF_TimeGet(datatime, yy=yr1,mm=mo1,dd=da1,&
          h=hr1,m=mn1,s=ss1,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readARMObs')

! data index from 0z with 30mn timestep. All this extensive processing is
! required because some of the ARM files report data at times different 
! from 0z. 
     
     call ESMF_TimeSet(currTime, yy=yr, mm=mo, &
          dd=da, h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readARMObs')

     data_index = (datatime - currtime)/armobs(source)%smos_ts + 1
     armobs(source)%smos_tindex(k,data_index) = data_index

     if(smos_qc_pcp(kk).eq.0) then 
        armobs(source)%smos_pcp(k,data_index)     = smos_pcp(kk)
     else
        armobs(source)%smos_pcp(k,data_index)     = LVT_rc%udef
     endif
     if(smos_qc_wspd(kk).eq.0) then 
        armobs(source)%smos_wspd(k,data_index)    = smos_wspd(kk)
     else
        armobs(source)%smos_wspd(k,data_index)    = LVT_rc%udef
     endif
     if(smos_qc_temp(kk).eq.0) then 
        armobs(source)%smos_temp(k,data_index)    = smos_temp(kk) + 273.15
     else
        armobs(source)%smos_temp(k,data_index)    = LVT_rc%udef
     endif
     if(smos_qc_rh(kk).eq.0) then 
        armobs(source)%smos_sh(k,data_index)    = &
             (smos_rh(kk)*0.622*10*0.611*exp(17.27*(&
             (273.15+smos_temp(kk)-273.15)/&
             (273.15+smos_temp(kk)-35.86)))/(smos_prs(kk)))/1000.0
     else
        armobs(source)%smos_sh(k,data_index)    = LVT_rc%udef
     endif
     if(smos_qc_prs(kk).eq.0) then 
        armobs(source)%smos_prs(k,data_index)    = smos_prs(kk)*1000.0
     else
        armobs(source)%smos_prs(k,data_index)    = LVT_rc%udef
     endif
     if(smos_qc_snow(kk).eq.0) then 
        armobs(source)%smos_snow(k,data_index)    = smos_snow(kk)
     else
        armobs(source)%smos_snow(k,data_index)    = LVT_rc%udef
     endif
  enddo

  deallocate(time)
  deallocate(smos_wspd)
  deallocate(smos_temp)
  deallocate(smos_pcp)
  deallocate(smos_rh)
  deallocate(smos_prs)
  deallocate(smos_snow)
  deallocate(smos_qc_wspd)
  deallocate(smos_qc_temp)
  deallocate(smos_qc_pcp)
  deallocate(smos_qc_rh)
  deallocate(smos_qc_prs)
  deallocate(smos_qc_snow)
#endif

end subroutine read_smos_file

!BOP
! 
! !ROUTINE: create_arm_swats_filename
! \label{create_arm_swats_filename}
!
! !INTERFACE: 
subroutine create_arm_swats_filename(odir, site_id, stnid, &
     yr, mo, da, filename, rc)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This routine creates a filename for ARM in-situ data files 
! from the Soil Water and Temperature System (SWATS). The 
! system produces estimates of soil moisture and temperature
! at different soil depths. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      ARM base directory
!   \item[site_id]   ARM site identifier
!   \item[stnid]     Station ID 
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the SWATS file
!   \item[rc]        return code (0-success, 1-failed to generate)
!  \end{description}!
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMRENTS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: site_id
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
  character(len=*), intent(out) :: filename
  integer,          intent(out) :: rc 
!EOP

  character*4       :: fyr
  character*2       :: fmo
  character*2       :: fda
  
  integer           :: fsize
  character*100     :: ls_comm, cmd2

  rc = 1 !fail to find the file

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//&
       trim(site_id)//'swats'//&
       trim(adjustl(stnid))//'.'//'b1.'//&
       trim(fyr)//trim(fmo)//trim(fda)//&
       '*cdf 2>&1 2>/dev/null > swats_file'
  
  cmd2 = 'wc -w swats_file > swats_file.wc'

  call system(ls_comm)
  call system(cmd2)

  open(110,file='swats_file.wc',form='formatted',action='read') 
  read(110,*) fsize
  close(110)

  if(fsize.gt.0) then 
     open(110,file='swats_file',form='formatted',action='read')
     read(110,'(a)') filename
     close(110)
     rc=0
  endif
  
end subroutine create_arm_swats_filename

!BOP
! 
! !ROUTINE: create_arm_ecor_flux_filename
! \label{create_arm_ecor_flux_filename}
!
! !INTERFACE: 
subroutine create_arm_ecor_flux_filename(odir, site_id, stnid, &
     yr, mo, da, filename, rc)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This routine creates a filename for ARM in-situ data files 
! from the eddy correlation (ECOR) flux measurement system. The 
! system produces 30mn estimates of vertical fluxes of 
! sensible and latent heat fluxes. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      ARM base directory
!   \item[site_id]   ARM site identifier
!   \item[stnid]     Station ID 
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the ECOR file
!   \item[rc]        return code (0-success, 1-failed to generate)
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMRENTS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: site_id
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
  character(len=*), intent(out) :: filename
  integer,          intent(out) :: rc 

!EOP

  character*4       :: fyr
  character*2       :: fmo
  character*2       :: fda
  
  integer           :: fsize
  character*100     :: ls_comm, cmd2

  rc = 1 !fail to find the file

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//&
       trim(site_id)//'30ecor'//&
       trim(adjustl(stnid))//'.'//'b1.'//&
       trim(fyr)//trim(fmo)//trim(fda)//&
       '*cdf 2>&1 2>/dev/null > ecor_file'
  
  cmd2 = 'wc -w ecor_file > ecor_file.wc'

  call system(ls_comm)
  call system(cmd2)

  open(110,file='ecor_file.wc',form='formatted',action='read') 
  read(110,*) fsize
  close(110)

  if(fsize.gt.0) then 
     open(110,file='ecor_file',form='formatted',action='read')
     read(110,'(a)') filename
     close(110)
     rc=0
  endif
  
end subroutine create_arm_ecor_flux_filename


!BOP
! 
! !ROUTINE: create_arm_baebbr_flux_filename
! \label{create_arm_baebbr_flux_filename}
!
! !INTERFACE: 
subroutine create_arm_baebbr_flux_filename(odir, site_id, stnid, &
     yr, mo, da, filename, rc)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This routine creates a filename for ARM in-situ data files 
! from the Energy Balance Bowen Ratio (BAEBBR) system. The 
! system produces 30mn estimates of vertical fluxes of 
! sensible and latent heat fluxes. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      ARM base directory
!   \item[site_id]   ARM site identifier
!   \item[stnid]     Station ID 
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the BAEBBR file
!   \item[rc]        return code (0-success, 1-failed to generate)
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMRENTS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: site_id
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
  character(len=*), intent(out) :: filename
  integer,          intent(out) :: rc 
!EOP

  character*4       :: fyr
  character*2       :: fmo
  character*2       :: fda
  
  integer           :: fsize
  character*100     :: ls_comm, cmd2

  rc = 1 !fail to find the file

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//&
       trim(site_id)//'30baebbr'//&
       trim(adjustl(stnid))//'.'//'s1.'//&
       trim(fyr)//trim(fmo)//trim(fda)//&
       '*cdf 2>&1 2>/dev/null > baebbr_file'
  
  cmd2 = 'wc -w baebbr_file > baebbr_file.wc'

  call system(ls_comm)
  call system(cmd2)

  open(110,file='baebbr_file.wc',form='formatted',action='read') 
  read(110,*) fsize
  close(110)

  if(fsize.gt.0) then 
     open(110,file='baebbr_file',form='formatted',action='read')
     read(110,'(a)') filename
     close(110)
     rc=0
  endif
  
end subroutine create_arm_baebbr_flux_filename


!BOP
! 
! !ROUTINE: create_arm_smos_filename
! \label{create_arm_smos_filename}
!
! !INTERFACE: 
subroutine create_arm_smos_filename(odir, site_id, stnid, &
     yr, mo, da, filename, rc)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This routine creates a filename for ARM in-situ data files 
! from the Surface Meteorological Observing System (SMOS). The 
! system produces 30mn estimates of surface meteorology variables.
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      ARM base directory
!   \item[site_id]   ARM site identifier
!   \item[stnid]     Station ID 
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the SMOS file
!   \item[rc]        return code (0-success, 1-failed to generate)
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMRENTS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: site_id
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
  character(len=*), intent(out) :: filename
  integer,          intent(out) :: rc 
!EOP

  character*4       :: fyr
  character*2       :: fmo
  character*2       :: fda
  
  integer           :: fsize
  character*100     :: ls_comm, cmd2

  rc = 1 !fail to find the file

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//&
       trim(site_id)//'30smos'//&
       trim(adjustl(stnid))//'.'//'b1.'//&
       trim(fyr)//trim(fmo)//trim(fda)//&
       '*cdf 2>&1 2>/dev/null > smos_file'
  
  cmd2 = 'wc -w smos_file > smos_file.wc'

  call system(ls_comm)
  call system(cmd2)

  open(110,file='smos_file.wc',form='formatted',action='read') 
  read(110,*) fsize
  close(110)

  if(fsize.gt.0) then 
     open(110,file='smos_file',form='formatted',action='read')
     read(110,'(a)') filename
     close(110)
     rc=0
  endif
  
end subroutine create_arm_smos_filename


!BOP
! 
! !ROUTINE: create_arm_ebbr_filename
! \label{create_arm_ebbr_filename}
!
! !INTERFACE: 
subroutine create_arm_ebbr_filename(odir, site_id, stnid, &
     yr, mo, da, filename, rc)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This routine creates a filename for ARM in-situ data files 
! from the Energy Balance Bowen Ratio (EBBR) system.
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      ARM base directory
!   \item[site_id]   ARM site identifier
!   \item[stnid]     Station ID 
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the EBBR file
!   \item[rc]        return code (0-success, 1-failed to generate)
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMRENTS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: site_id
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
  character(len=*), intent(out) :: filename
  integer,          intent(out) :: rc 
!EOP

  character*4       :: fyr
  character*2       :: fmo
  character*2       :: fda
  
  integer           :: fsize
  character*100     :: ls_comm, cmd2

  rc = 1 !fail to find the file

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//&
       trim(site_id)//'30ebbr'//&
       trim(adjustl(stnid))//'.'//'b1.'//&
       trim(fyr)//trim(fmo)//trim(fda)//&
       '*cdf 2>&1 2>/dev/null > ebbr_file'
  
  cmd2 = 'wc -w ebbr_file > ebbr_file.wc'

  call system(ls_comm)
  call system(cmd2)

  open(110,file='ebbr_file.wc',form='formatted',action='read') 
  read(110,*) fsize
  close(110)

  if(fsize.gt.0) then 
     open(110,file='ebbr_file',form='formatted',action='read')
     read(110,'(a)') filename
     close(110)
     rc=0
  endif
  
end subroutine create_arm_ebbr_filename
