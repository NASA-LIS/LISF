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
! !ROUTINE: readLPRM_AMSREsmObs
! \label{readLPRM_AMSREsmObs}
!
! !INTERFACE: 
subroutine readLPRM_AMSREsmObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit, LVT_getNextUnitNumber, & 
       LVT_releaseUnitNumber
  use LVT_timeMgrMod,   only : LVT_get_julss
  use LPRM_AMSREsm_obsMod, only : LPRM_AMSREsmobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)      :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This subroutine provides the data reader for the standard 
! LPRM soil moisture retrieval product. 
! 
! !NOTES: 
!  The mismatches between the LVT time and LPRM data time is not
!  handled. This is not an issue as long as LVT analysis is not
!  done for sub-daily intervals. 
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
!  25 May 2012: Sujay Kumar, Updated for LPRM version 5. 
! 
!EOP


  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character*100     :: fname_A, fname_D,fname
  real              :: smobs_A(LVT_rc%lnc*LVT_rc%lnr)
  real              :: smobs_D(LVT_rc%lnc*LVT_rc%lnr)
  real              :: smobs(LVT_rc%lnc*LVT_rc%lnr)
  real              :: lat,lon
!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
       LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  smobs  = LVT_rc%udef
  if(LPRM_AMSREsmobs(source)%version.eq."V05") then 
     LPRM_AMSREsmobs(source)%smobs = LVT_rc%udef
     smobs_A= LVT_rc%udef
     smobs_D= LVT_rc%udef
     
     if(LPRM_AMSREsmobs(source)%startmode.or.alarmCheck.or.&
          LVT_rc%resetFlag(source)) then 
        
        LVT_rc%resetFlag(source) = .false. 
        LPRM_AMSREsmobs(source)%startmode = .false. 
        
        call create_LPRM_AMSREsm_filename(LPRM_AMSREsmobs(source)%odir, &
             'A',LVT_rc%dyr(source), LVT_rc%dmo(source), &
             LVT_rc%dda(source), fname_A)
        
        inquire(file=trim(fname_A),exist=file_exists)
        
        if(file_exists) then
           
           write(LVT_logunit,*) '[INFO] Reading ..',trim(fname_A)
           call read_LPRM_data(source, fname_A,smobs_A)
        endif
        
        call create_LPRM_AMSREsm_filename(LPRM_AMSREsmobs(source)%odir, &
             'D',LVT_rc%dyr(source), LVT_rc%dmo(source), &
             LVT_rc%dda(source), fname_D)
        
        inquire(file=trim(fname_D),exist=file_exists)
        if(file_exists) then
           
           write(LVT_logunit,*) '[INFO] Reading ..',trim(fname_D)
           call read_LPRM_data(source, fname_D,smobs_D)
        endif
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(smobs_A(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
                 LPRM_AMSREsmobs(source)%smobs(c,r) = smobs_A(c+(r-1)*LVT_rc%lnc)
              endif
              if(smobs_D(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
                 LPRM_AMSREsmobs(source)%smobs(c,r) = smobs_D(c+(r-1)*LVT_rc%lnc)
              endif
              
           enddo
        enddo
     endif
  elseif(LPRM_AMSREsmobs(source)%version.eq."GES-DISC") then 
     LPRM_AMSREsmobs(source)%smobs = LVT_rc%udef
     
     if(LPRM_AMSREsmobs(source)%startmode.or.alarmCheck.or.&
          LVT_rc%resetFlag(source)) then 
        
        LVT_rc%resetFlag(source) = .false. 
        LPRM_AMSREsmobs(source)%startmode = .false. 
        
        call create_LPRM_AMSREsm_GESDISC_filename(&
             LPRM_AMSREsmobs(source)%odir, &
             LVT_rc%dyr(source), LVT_rc%dmo(source), &
             LVT_rc%dda(source), fname)
        
        inquire(file=trim(fname),exist=file_exists)
        
        if(file_exists) then
           
           write(LVT_logunit,*) '[INFO] Reading ..',trim(fname)
           call read_LPRM_GESDISC_data(source, fname,smobs)
        endif
        
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(smobs(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
                 LPRM_AMSREsmobs(source)%smobs(c,r) = smobs(c+(r-1)*LVT_rc%lnc)
              endif
           enddo
        enddo
     endif

  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, &
       LPRM_AMSREsmobs(source)%smobs,vlevel=1,units="m3/m3")

end subroutine readLPRM_AMSREsmObs


!BOP
! 
! !ROUTINE: read_LPRM_data
! \label(read_LPRM_data)
!
! !INTERFACE:
subroutine read_LPRM_data(source, fname, smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LVT_coreMod,  only : LVT_rc
  use LVT_logMod,   only : LVT_verify
  use map_utils,    only : latlon_to_ij
  use LPRM_AMSREsm_obsMod, only : LPRM_AMSREsmobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: source
  character (len=*)             :: fname
  real                          :: smobs_ip(LVT_rc%lnc*LVT_rc%lnr)


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the LPRM NETCDF file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  land surface temperature is below freezing, if rain is present, if 
!  RFI is present, if residual error is above 0.5 or if optical depth
!  is above 0.8. Finally the routine combines both the C-band and X-band
!  retrievals. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the LPRM AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LVT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  real                        :: rfi(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: tskin(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: rainf(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: optc(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: optx(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: smerrc(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: smerrx(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: smc(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: smx(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: sm_combined(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: sm_raw(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: sm_flags(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)
  real                        :: sm_cdf(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)


  real                        :: sm_data(LPRM_AMSREsmobs(source)%lprmnc* & 
       LPRM_AMSREsmobs(source)%lprmnr)
  logical*1                   :: sm_data_b(LPRM_AMSREsmobs(source)%lprmnc* & 
       LPRM_AMSREsmobs(source)%lprmnr)
  logical*1                   :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)

  integer                     :: c,r,i,j
  real                        :: rlat,rlon,ri,rj
  integer                     :: nid, rfiid
  integer                     :: smcombId,smcdfId
  integer                     :: tskinId, rainfId, optCId, optXid,smfId
  integer                     :: smerrCid, smerrXid, smcId, smXId
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LVT_verify(ios,'Error opening file '//trim(fname))
  
  ios = nf90_inq_varid(nid, 'RFI_Daily',rfiid)
  call LVT_verify(ios, 'Error nf90_inq_varid: RFI_Daily')
  
  ios = nf90_inq_varid(nid, 'Land_Surface_Temperature',tskinid)
  call LVT_verify(ios, 'Error nf90_inq_varid: Land_Surface_Temperature')
  
  ios = nf90_inq_varid(nid, 'Rainfall_Diagnostic',rainfid)
  call LVT_verify(ios, 'Error nf90_inq_varid: Rainfall_Diagnostic')
  
  ios = nf90_inq_varid(nid, 'Optical_Depth_from_C_band',optCid)
  call LVT_verify(ios, 'Error nf90_inq_varid: Optical_Depth_from_C_band')
  
  ios = nf90_inq_varid(nid, 'Optical_Depth_from_X_band',optXid)
  call LVT_verify(ios, 'Error nf90_inq_varid: Optical_Depth_from_X_band')
  
  ios = nf90_inq_varid(nid, 'Soil_Moisture_Error_from_C_band',smerrCid)
  call LVT_verify(ios, &
       'Error nf90_inq_varid: Soil_Moisture_Error_from_C_band')
  
  ios = nf90_inq_varid(nid, 'Soil_Moisture_Error_from_X_band',smerrXid)
  call LVT_verify(ios, &
       'Error nf90_inq_varid: Soil_Moisture_Error_from_X_band')
  
  ios = nf90_inq_varid(nid, 'SM_Flags',smfid)
  call LVT_verify(ios, 'Error nf90_inq_varid: SM_Flags')

  ios = nf90_inq_varid(nid, 'SM_Combined',smcombid)
  call LVT_verify(ios, 'Error nf90_inq_varid: SM_Combined')

  if(LPRM_AMSREsmobs(source)%rawdata.eq.0) then 
     ios = nf90_inq_varid(nid, 'SM_CDF',smcdfid)
     call LVT_verify(ios, 'Error nf90_inq_varid: SM_CDF')
  endif

  !values
  ios = nf90_get_var(nid,rfiid, rfi)
  call LVT_verify(ios, 'Error nf90_get_var: RFI_Daily')
  
  ios = nf90_get_var(nid,tskinid, tskin)
  call LVT_verify(ios, 'Error nf90_get_var: Land_Surface_Temperature')
  
  ios = nf90_get_var(nid, rainfid, rainf)
  call LVT_verify(ios, 'Error nf90_get_var: Rainfall_Diagnostic')
  
  ios = nf90_get_var(nid, optcid, optc)
  call LVT_verify(ios, 'Error nf90_get_var: Optical_Depth_from_C_band')
  
  ios = nf90_get_var(nid, optxid, optx)
  call LVT_verify(ios, 'Error nf90_get_var: Optical_Depth_from_X_band')
  
  ios = nf90_get_var(nid, smerrcid, smerrc)
  call LVT_verify(ios, 'Error nf90_get_var: Soil_Moisture_Error_from_C_band')
  
  ios = nf90_get_var(nid, smerrxid,smerrx)
  call LVT_verify(ios, 'Error nf90_get_var: Soil_Moisture_Error_from_X_band')
  
  ios = nf90_get_var(nid, smfid, sm_flags)
  call LVT_verify(ios, 'Error nf90_get_var: SM_Flags')

  ios = nf90_get_var(nid, smcombid, sm_combined)
  call LVT_verify(ios, 'Error nf90_get_var: SM_Combined')

  if(LPRM_AMSREsmobs(source)%rawdata.eq.0) then 
     ios = nf90_get_var(nid, smcdfid, sm_cdf)
     call LVT_verify(ios, 'Error nf90_get_var: SM_CDF')
  endif
  
  ios = nf90_close(ncid=nid)
  call LVT_verify(ios,'Error closing file '//trim(fname))

  smc = 1.0
  smx = 1.0
  sm_raw = -9999.0

  do r=1, LPRM_AMSREsmobs(source)%lprmnr
     do c=1, LPRM_AMSREsmobs(source)%lprmnc
!--------------------------------------------------------------------------
! Choose soil moisture retrievals only where the reported optical 
! depth values are between 0 and 0.8
!--------------------------------------------------------------------------
        if(optc(c,r)*0.01.gt.0.8.or.optc(c,r).lt.0) smc(c,r) = -1
        if(optx(c,r)*0.01.gt.0.8.or.optx(c,r).lt.0) smx(c,r) = -1
!--------------------------------------------------------------------------
! Reject retrievals when soil temperature is below freezing. 
!--------------------------------------------------------------------------
        if(tskin(c,r)*0.1.lt.273.15.or.tskin(c,r).lt.0) smc(c,r) = -1
        if(tskin(c,r)*0.1.lt.273.15.or.tskin(c,r).lt.0) smx(c,r) = -1
!--------------------------------------------------------------------------
! Reject retrievals when RFI is non zero
!--------------------------------------------------------------------------
        if(rfi(c,r).ne.0) smc(c,r) = -1
        if(rfi(c,r).ne.0.and.rfi(c,r).ne.9) smx(c,r) = -1
!--------------------------------------------------------------------------
! Reject retrievals when residual error is large (> 0.2)
!-------------------------------------------------------------------------- 
        if(smerrc(c,r)*0.01.gt.0.2.or.smerrc(c,r)*0.01.lt.0) smc(c,r) = -1
        if(smerrx(c,r)*0.01.gt.0.2.or.smerrx(c,r)*0.01.lt.0) smx(c,r) = -1
!--------------------------------------------------------------------------
! Reject retrievals when rain exists
!-------------------------------------------------------------------------- 
        if(rainf(c,r).ne.0) then 
           smc(c,r) = -1
           smx(c,r) = -1
        endif
!--------------------------------------------------------------------------
! Reject retrievals if flag is not set
! -99: Edge of the swath is masked (AMSR-E observations are disturbed along the edges)
!-10: Default (open water (oceans), no brightness temperature observations over land and land surface temperature below freezing)
!-9: No land surface temperature retrieval possible (Thomas Holmes developed some filtering of suspicious areas for retrieving LST)
!-5: Vegetation is too dense for regions where we use C-band
!-4: Vegetation is too dense for regions where we had to switch to X-band due to RFI in C-band
!-3: Both observations in C- and X-band frequencies are suspicious (happens sometimes in India)
!6: Areas where we used C-band (6.9 GHz) observations
!10: Areas where we used X-band (10.7 GHz) observations
!-------------------------------------------------------------------------- 
        if(sm_flags(c,r).lt.0) then 
           smc(c,r) = -1
           smx(c,r) = -1
        endif
!--------------------------------------------------------------------------
! Apply the QC flags 
!-------------------------------------------------------------------------- 
        if(smc(c,r).gt.0.and.smx(c,r).gt.0) then 
           if(LPRM_AMSREsmobs(source)%rawdata.eq.1) then 
              if(sm_combined(c,r).gt.0) then 
                 sm_raw(c,LPRM_AMSREsmobs(source)%lprmnr-r+1)  = &
                      sm_combined(c,r)/100.0
              else
                 sm_raw(c,LPRM_AMSREsmobs(source)%lprmnr-r+1) =LVT_rc%udef
              endif
           else
              if(sm_cdf(c,r).gt.0.and.sm_cdf(c,r).lt.50) then 
                 sm_combined(c,LPRM_AMSREsmobs(source)%lprmnr-r+1) =&
                      sm_cdf(c,r)/100.0
              else
                 sm_combined(c,LPRM_AMSREsmobs(source)%lprmnr-r+1) =LVT_rc%udef
              endif
           endif
        else
           if(LPRM_AMSREsmobs(source)%rawdata.eq.1) then 
              sm_raw(c,LPRM_AMSREsmobs(source)%lprmnr-r+1) =LVT_rc%udef
           else
              sm_combined(c,LPRM_AMSREsmobs(source)%lprmnr-r+1) = LVT_rc%udef
           endif
        endif
        if(LPRM_AMSREsmobs(source)%rawdata.eq.1) then 
           if(sm_raw(c,LPRM_AMSREsmobs(source)%lprmnr-r+1).gt.1) then 
              print*, c,r,&
                   sm_raw(c,LPRM_AMSREsmobs(source)%lprmnr-r+1)
           endif
        else
           if(sm_combined(c,LPRM_AMSREsmobs(source)%lprmnr-r+1).gt.1) then
              print*, c,r,&
                   sm_combined(c,LPRM_AMSREsmobs(source)%lprmnr-r+1)   
              stop
           endif
        endif
     enddo
  enddo
  if(LPRM_AMSREsmobs(source)%rawdata.eq.1) then 
     sm_combined = sm_raw
  endif

  do r=1, LPRM_AMSREsmobs(source)%lprmnr
     do c=1, LPRM_AMSREsmobs(source)%lprmnc
        sm_data(c+(r-1)*LPRM_AMSREsmobs(source)%lprmnc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LVT_rc%udef) then 
           sm_data_b(c+(r-1)*LPRM_AMSREsmobs(source)%lprmnc) = .true. 
        else
           sm_data_b(c+(r-1)*LPRM_AMSREsmobs(source)%lprmnc) = .false.
        endif
     enddo
  enddo
  
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
  call bilinear_interp(LVT_rc%gridDesc(:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       LPRM_AMSREsmobs(source)%lprmnc*LPRM_AMSREsmobs(source)%lprmnr, &
       LVT_rc%lnc*LVT_rc%lnr, &
       LPRM_AMSREsmobs(source)%rlat, LPRM_AMSREsmobs(source)%rlon, &
       LPRM_AMSREsmobs(source)%w11, LPRM_AMSREsmobs(source)%w12, &
       LPRM_AMSREsmobs(source)%w21, LPRM_AMSREsmobs(source)%w22, &
       LPRM_AMSREsmobs(source)%n11, LPRM_AMSREsmobs(source)%n12, &
       LPRM_AMSREsmobs(source)%n21, LPRM_AMSREsmobs(source)%n22, &
       LVT_rc%udef, ios)

#endif
  
end subroutine read_LPRM_data

!BOP
! 
! !ROUTINE: read_LPRM_GESDISC_data
! \label(read_LPRM_GESDISC_data)
!
! !INTERFACE:
subroutine read_LPRM_GESDISC_data(source, fname, smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LVT_coreMod,  only : LVT_rc
  use LVT_logMod,   only : LVT_verify
  use map_utils,    only : latlon_to_ij
  use LPRM_AMSREsm_obsMod, only : LPRM_AMSREsmobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: source
  character (len=*)             :: fname
  real                          :: smobs_ip(LVT_rc%lnc*LVT_rc%lnr)


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the LPRM NETCDF file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  land surface temperature is below freezing, if rain is present, if 
!  RFI is present, if residual error is above 0.5 or if optical depth
!  is above 0.8. Finally the routine combines both the C-band and X-band
!  retrievals. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the LPRM AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LVT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  real                        :: smerrc(LPRM_AMSREsmobs(source)%lprmnr,&
       LPRM_AMSREsmobs(source)%lprmnc)
  real                        :: smerrx(LPRM_AMSREsmobs(source)%lprmnr,&
       LPRM_AMSREsmobs(source)%lprmnc)
  real                        :: smc(LPRM_AMSREsmobs(source)%lprmnr,&
       LPRM_AMSREsmobs(source)%lprmnc)
  real                        :: smx(LPRM_AMSREsmobs(source)%lprmnr,&
       LPRM_AMSREsmobs(source)%lprmnc)
  real                        :: sm_combined(LPRM_AMSREsmobs(source)%lprmnc,&
       LPRM_AMSREsmobs(source)%lprmnr)

  real                        :: sm_data(LPRM_AMSREsmobs(source)%lprmnc* & 
       LPRM_AMSREsmobs(source)%lprmnr)
  logical*1                   :: sm_data_b(LPRM_AMSREsmobs(source)%lprmnc* & 
       LPRM_AMSREsmobs(source)%lprmnr)
  logical*1                   :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)

  integer                     :: c,r,i,j
  real                        :: rlat,rlon,ri,rj
  integer                     :: nid
  integer                     :: smcombId
  integer                     :: tskinId
  integer                     :: smerrCid, smerrXid, smcId, smXId
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  if(LPRM_AMSREsmobs(source)%channel.eq."C-band") then 
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios,'Error opening file '//trim(fname))

     ios = nf90_inq_varid(nid, 'sm_c_error',smerrCid)
     call LVT_verify(ios, &
          'Error nf90_inq_varid: sm_c_error')
     ios = nf90_get_var(nid, smerrcid, smerrc)
     call LVT_verify(ios, 'Error nf90_get_var: sm_c_error')
     
     ios = nf90_inq_varid(nid, 'soil_moisture_c',smcid)
     call LVT_verify(ios, 'Error nf90_inq_varid: soil_moisture_c')
     
     ios = nf90_get_var(nid, smcid, smc)
     call LVT_verify(ios, 'Error nf90_get_var: soil_moisture_c')

     ios = nf90_close(ncid=nid)
     call LVT_verify(ios,'Error closing file '//trim(fname))
     
     sm_combined = LVT_rc%udef

     do r=1, LPRM_AMSREsmobs(source)%lprmnr
        do c=1, LPRM_AMSREsmobs(source)%lprmnc
!--------------------------------------------------------------------------
! Reject retrievals when residual error is large (> 0.2)
!-------------------------------------------------------------------------- 
           if(smerrc(r,c)*0.01.gt.0.2.or.smerrc(r,c)*0.01.lt.0) smc(r,c) = -1
!--------------------------------------------------------------------------
! Apply the QC flags 
!-------------------------------------------------------------------------- 
           if(smc(r,c).gt.0) then 
              sm_combined(c,LPRM_AMSREsmobs(source)%lprmnr-r+1)  = &
                   smc(r,c)/100.0
           else
              sm_combined(c,LPRM_AMSREsmobs(source)%lprmnr-r+1) =LVT_rc%udef
           endif
        enddo
     enddo
  elseif(LPRM_AMSREsmobs(source)%channel.eq."X-band") then 
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios,'Error opening file '//trim(fname))

     ios = nf90_inq_varid(nid, 'sm_x_error',smerrXid)
     call LVT_verify(ios, &
          'Error nf90_inq_varid: sm_x_error')
     ios = nf90_get_var(nid, smerrxid,smerrx)
     call LVT_verify(ios, 'Error nf90_get_var: sm_x_error')

     ios = nf90_inq_varid(nid, 'soil_moisture_x',smxid)
     call LVT_verify(ios, 'Error nf90_inq_varid: soil_moisture_x')
     
     ios = nf90_get_var(nid, smxid, smx)
     call LVT_verify(ios, 'Error nf90_get_var: soil_moisture_x')

     ios = nf90_close(ncid=nid)
     call LVT_verify(ios,'Error closing file '//trim(fname))
     
     sm_combined = LVT_rc%udef

     do r=1, LPRM_AMSREsmobs(source)%lprmnr
        do c=1, LPRM_AMSREsmobs(source)%lprmnc
!--------------------------------------------------------------------------
! Reject retrievals when residual error is large (> 0.2)
!-------------------------------------------------------------------------- 
           if(smerrx(r,c)*0.01.gt.0.2.or.smerrx(r,c)*0.01.lt.0) smx(r,c) = -1
!--------------------------------------------------------------------------
! Apply the QC flags 
!-------------------------------------------------------------------------- 
           if(smx(r,c).gt.0) then 
              sm_combined(c,LPRM_AMSREsmobs(source)%lprmnr-r+1)  = &
                   smx(r,c)/100.0
           else
              sm_combined(c,LPRM_AMSREsmobs(source)%lprmnr-r+1) =LVT_rc%udef
           endif
        enddo
     enddo
  endif

  do r=1, LPRM_AMSREsmobs(source)%lprmnr
     do c=1, LPRM_AMSREsmobs(source)%lprmnc
        sm_data(c+(r-1)*LPRM_AMSREsmobs(source)%lprmnc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LVT_rc%udef) then 
           sm_data_b(c+(r-1)*LPRM_AMSREsmobs(source)%lprmnc) = .true. 
        else
           sm_data_b(c+(r-1)*LPRM_AMSREsmobs(source)%lprmnc) = .false.
        endif
     enddo
  enddo
  
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
  call bilinear_interp(LVT_rc%gridDesc(:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       LPRM_AMSREsmobs(source)%lprmnc*LPRM_AMSREsmobs(source)%lprmnr, &
       LVT_rc%lnc*LVT_rc%lnr, &
       LPRM_AMSREsmobs(source)%rlat, LPRM_AMSREsmobs(source)%rlon, &
       LPRM_AMSREsmobs(source)%w11, LPRM_AMSREsmobs(source)%w12, &
       LPRM_AMSREsmobs(source)%w21, LPRM_AMSREsmobs(source)%w22, &
       LPRM_AMSREsmobs(source)%n11, LPRM_AMSREsmobs(source)%n12, &
       LPRM_AMSREsmobs(source)%n21, LPRM_AMSREsmobs(source)%n22, &
       LVT_rc%udef, ios)

#endif
  
end subroutine read_LPRM_GESDISC_data

!BOP
! !ROUTINE: create_LPRM_AMSREsm_filename
! \label{create_LPRM_AMSREsm_filename}
! 
! !INTERFACE: 
subroutine create_LPRM_AMSREsm_filename(ndir, path,yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: path
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the LPRM AMSRE filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the LPRM AMSRE soil moisture directory
!  \item[path] name of the sensor path (A-ascending, D-descending)
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated LPRM filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/'//trim(fyr)//'-'//trim(fmo)//'/AMSR_L3_LPRMv05_' &
       //trim(path)//'_' &
       //trim(fyr)//trim(fmo)//trim(fda)//'T000000_' &
       //trim(fyr)//trim(fmo)//trim(fda)//'T235959D.v05.nc'
  
end subroutine create_LPRM_AMSREsm_filename


!BOP
! !ROUTINE: create_LPRM_AMSREsm_GESDISC_filename
! \label{create_LPRM_AMSREsm_GESDISC_filename}
! 
! !INTERFACE: 
subroutine create_LPRM_AMSREsm_GESDISC_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the LPRM AMSRE filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the LPRM AMSRE soil moisture directory
!  \item[path] name of the sensor path (A-ascending, D-descending)
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated LPRM filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/'//trim(fyr)//'/'//&
       '/LPRM-AMSR_E_L3_A_SOILM3_V002_' &
       //trim(fyr)//trim(fmo)//trim(fda)//'.nc'
  
end subroutine create_LPRM_AMSREsm_GESDISC_filename
