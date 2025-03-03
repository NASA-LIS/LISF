!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readLPRM_AMSRE_ANNdata
! \label{readLPRM_AMSRE_ANNdata}
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readLPRM_AMSRE_ANNdata(n,iomode)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_ANNMod
  use LPRM_AMSREsm_ANNdataMod
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: iomode
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the standard 
! LPRM soil moisture retrieval product. 
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character(len=LDT_CONST_PATH_LEN) :: fname_A, fname_D
  real              :: smobs_A(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real              :: smobs_D(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real              :: lat,lon
!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------

  timenow = float(LDT_rc%hr)*3600 + 60*LDT_rc%mn + LDT_rc%ss
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  LPRM_AMSREsmANNdata(n)%smobs = LDT_rc%udef
  smobs_A= LDT_rc%udef
  smobs_D= LDT_rc%udef
        
  if(LPRM_AMSREsmANNdata(n)%startmode.or.alarmCheck) then 
     
     LPRM_AMSREsmANNdata(n)%startmode = .false. 

     call create_LPRM_AMSREsm_ANN_filename(LPRM_AMSREsmANNdata(n)%odir, &
          'A',LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname_A)

     inquire(file=trim(fname_A),exist=file_exists)
     if(file_exists) then

        write(LDT_logunit,*) 'Reading ..',trim(fname_A)
        call read_LPRM_data_ANN(n, fname_A,smobs_A)
     endif

     call create_LPRM_AMSREsm_ANN_filename(LPRM_AMSREsmANNdata(n)%odir, &
          'D',LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname_D)

     inquire(file=trim(fname_D),exist=file_exists)
     if(file_exists) then

        write(LDT_logunit,*) 'Reading ..',trim(fname_D)
        call read_LPRM_data_ANN(n, fname_D,smobs_D)
     endif
     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(smobs_A(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
              LPRM_AMSREsmANNdata(n)%smobs(c,r) = smobs_A(c+(r-1)*LDT_rc%lnc(n))
           endif
           if(smobs_D(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
              LPRM_AMSREsmANNdata(n)%smobs(c,r) = smobs_D(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  endif

  call LDT_logSingleANNdata(n,&
       LPRM_AMSREsmANNdata(n)%smobs,pindex=1,iomode = iomode, &
       name = "SoilMoist",&
       units="m3/m3")

end subroutine readLPRM_AMSRE_ANNdata


!BOP
! 
! !ROUTINE: read_LPRM_data_ANN
! \label(read_LPRM_data_ANN)
!
! !INTERFACE:
subroutine read_LPRM_data_ANN(n, fname, smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_verify
  use map_utils,    only : latlon_to_ij
  use LPRM_AMSREsm_ANNdataMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname
  real                          :: smobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))


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
!  \item[smobs\_ip]    soil moisture data processed to the LDT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  real                        :: rfi(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: tskin(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: rainf(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: optc(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: optx(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: smerrc(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: smerrx(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: smc(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: smx(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: sm_combined(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: sm_flags(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)
  real                        :: sm_cdf(LPRM_AMSREsmANNdata(n)%lprmnc,&
       LPRM_AMSREsmANNdata(n)%lprmnr)


  real                        :: sm_data(LPRM_AMSREsmANNdata(n)%lprmnc* & 
       LPRM_AMSREsmANNdata(n)%lprmnr)
  logical*1                   :: sm_data_b(LPRM_AMSREsmANNdata(n)%lprmnc* & 
       LPRM_AMSREsmANNdata(n)%lprmnr)
  logical*1                   :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  integer                     :: c,r,i,j
  real                        :: rlat,rlon,ri,rj
  integer                     :: nid, rfiid
  integer                     :: smcombId,smcdfId
  integer                     :: tskinId, rainfId, optCId, optXid,smfId
  integer                     :: smerrCid, smerrXid, smcId, smXId
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LDT_verify(ios,'Error opening file '//trim(fname))
  
  ios = nf90_inq_varid(nid, 'RFI_Daily',rfiid)
  call LDT_verify(ios, 'Error nf90_inq_varid: RFI_Daily')
  
  ios = nf90_inq_varid(nid, 'Land_Surface_Temperature',tskinid)
  call LDT_verify(ios, 'Error nf90_inq_varid: Land_Surface_Temperature')
  
  ios = nf90_inq_varid(nid, 'Rainfall_Diagnostic',rainfid)
  call LDT_verify(ios, 'Error nf90_inq_varid: Rainfall_Diagnostic')
  
  ios = nf90_inq_varid(nid, 'Optical_Depth_from_C_band',optCid)
  call LDT_verify(ios, 'Error nf90_inq_varid: Optical_Depth_from_C_band')
  
  ios = nf90_inq_varid(nid, 'Optical_Depth_from_X_band',optXid)
  call LDT_verify(ios, 'Error nf90_inq_varid: Optical_Depth_from_X_band')
  
  ios = nf90_inq_varid(nid, 'Soil_Moisture_Error_from_C_band',smerrCid)
  call LDT_verify(ios, &
       'Error nf90_inq_varid: Soil_Moisture_Error_from_C_band')
  
  ios = nf90_inq_varid(nid, 'Soil_Moisture_Error_from_X_band',smerrXid)
  call LDT_verify(ios, &
       'Error nf90_inq_varid: Soil_Moisture_Error_from_X_band')
  
  ios = nf90_inq_varid(nid, 'SM_Flags',smfid)
  call LDT_verify(ios, 'Error nf90_inq_varid: SM_Flags')

  ios = nf90_inq_varid(nid, 'SM_Combined',smcombid)
  call LDT_verify(ios, 'Error nf90_inq_varid: SM_Combined')

  if(LPRM_AMSREsmANNdata(n)%rawdata.eq.0) then 
     ios = nf90_inq_varid(nid, 'SM_CDF',smcdfid)
     call LDT_verify(ios, 'Error nf90_inq_varid: SM_CDF')
  endif

  !values
  ios = nf90_get_var(nid,rfiid, rfi)
  call LDT_verify(ios, 'Error nf90_get_var: RFI_Daily')
  
  ios = nf90_get_var(nid,tskinid, tskin)
  call LDT_verify(ios, 'Error nf90_get_var: Land_Surface_Temperature')
  
  ios = nf90_get_var(nid, rainfid, rainf)
  call LDT_verify(ios, 'Error nf90_get_var: Rainfall_Diagnostic')
  
  ios = nf90_get_var(nid, optcid, optc)
  call LDT_verify(ios, 'Error nf90_get_var: Optical_Depth_from_C_band')
  
  ios = nf90_get_var(nid, optxid, optx)
  call LDT_verify(ios, 'Error nf90_get_var: Optical_Depth_from_X_band')
  
  ios = nf90_get_var(nid, smerrcid, smerrc)
  call LDT_verify(ios, 'Error nf90_get_var: Soil_Moisture_Error_from_C_band')
  
  ios = nf90_get_var(nid, smerrxid,smerrx)
  call LDT_verify(ios, 'Error nf90_get_var: Soil_Moisture_Error_from_X_band')
  
  ios = nf90_get_var(nid, smfid, sm_flags)
  call LDT_verify(ios, 'Error nf90_get_var: SM_Flags')

  ios = nf90_get_var(nid, smcombid, sm_combined)
  call LDT_verify(ios, 'Error nf90_get_var: SM_Combined')

  if(LPRM_AMSREsmANNdata(n)%rawdata.eq.0) then 
     ios = nf90_get_var(nid, smcdfid, sm_cdf)
     call LDT_verify(ios, 'Error nf90_get_var: SM_CDF')
  endif
  
  ios = nf90_close(ncid=nid)
  call LDT_verify(ios,'Error closing file '//trim(fname))

  smc = 1.0
  smx = 1.0

  do r=1, LPRM_AMSREsmANNdata(n)%lprmnr
     do c=1, LPRM_AMSREsmANNdata(n)%lprmnc
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
           if(LPRM_AMSREsmANNdata(n)%rawdata.eq.1) then 
              sm_combined(c,LPRM_AMSREsmANNdata(n)%lprmnr-r+1)  = &
                   sm_combined(c,r)/100.0
           else
              if(sm_cdf(c,r).gt.0.and.sm_cdf(c,r).lt.50) then 
                 sm_combined(c,LPRM_AMSREsmANNdata(n)%lprmnr-r+1) =&
                      sm_cdf(c,r)/100.0
              else
                 sm_combined(c,LPRM_AMSREsmANNdata(n)%lprmnr-r+1) =LDT_rc%udef
              endif
           endif
        else
           sm_combined(c,LPRM_AMSREsmANNdata(n)%lprmnr-r+1) = LDT_rc%udef
        endif
        if(sm_combined(c,LPRM_AMSREsmANNdata(n)%lprmnr-r+1).gt.1) then
           print*, c,r,&
                sm_combined(c,LPRM_AMSREsmANNdata(n)%lprmnr-r+1)   
           stop
        endif
     enddo
  enddo

  do r=1, LPRM_AMSREsmANNdata(n)%lprmnr
     do c=1, LPRM_AMSREsmANNdata(n)%lprmnc
        sm_data(c+(r-1)*LPRM_AMSREsmANNdata(n)%lprmnc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LDT_rc%udef) then 
           sm_data_b(c+(r-1)*LPRM_AMSREsmANNdata(n)%lprmnc) = .true. 
        else
           sm_data_b(c+(r-1)*LPRM_AMSREsmANNdata(n)%lprmnc) = .false.
        endif
     enddo
  enddo
  
!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  call bilinear_interp(LDT_rc%gridDesc(n,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       LPRM_AMSREsmANNdata(n)%lprmnc*LPRM_AMSREsmANNdata(n)%lprmnr, &
       LDT_rc%lnc(n)*LDT_rc%lnr(n), &
       LPRM_AMSREsmANNdata(n)%rlat, LPRM_AMSREsmANNdata(n)%rlon, &
       LPRM_AMSREsmANNdata(n)%w11, LPRM_AMSREsmANNdata(n)%w12, &
       LPRM_AMSREsmANNdata(n)%w21, LPRM_AMSREsmANNdata(n)%w22, &
       LPRM_AMSREsmANNdata(n)%n11, LPRM_AMSREsmANNdata(n)%n12, &
       LPRM_AMSREsmANNdata(n)%n21, LPRM_AMSREsmANNdata(n)%n22, &
       LDT_rc%udef, ios)

#endif
  
end subroutine read_LPRM_data_ANN

!BOP
! !ROUTINE: create_LPRM_AMSREsm_ANN_filename
! \label{create_LPRM_AMSREsm_ANN_filename}
! 
! !INTERFACE: 
subroutine create_LPRM_AMSREsm_ANN_filename(ndir, path,yr, mo,da, filename)
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
  
end subroutine create_LPRM_AMSREsm_ANN_filename
