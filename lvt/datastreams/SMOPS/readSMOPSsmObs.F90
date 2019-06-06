!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readSMOPSsmObs
! \label{readSMOPSsmObs}
!
! !INTERFACE: 
subroutine readSMOPSsmObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SMOPSsm_obsMod, only : SMOPSsmobs
  !use LVT_DAobsDataMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer,   intent(in) :: source
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!  11 Dec 2014: Sujay Kumar, Initial Specification
!  30 May 2018 : Mahdi Navari, Updated to 
!                support reading ASCAT Metop A and B & real time data 
!                support reading different version of the SMOPS data 
!                support binary QC flag as implemented in the LIS and LDT 
! ! 
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j, status
  integer           :: yr, mo, da, hr, mn, ss
  character*100     :: fname
  real              :: smobs(LVT_rc%lnc*LVT_rc%lnr)
  real              :: smobs_2d(LVT_rc%lnc,LVT_rc%lnr)
  real              :: relsmc(LVT_rc%lnc,LVT_rc%lnr)
  real              :: smobs_av(LVT_rc%lnc, LVT_rc%lnr)
  integer           :: count_av(LVT_rc%lnc, LVT_rc%lnr)
  type(ESMF_Time)   :: cTime, stTime, enTime

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
!MN : 86400 --> 3600 to read 6-hourly NRT data
  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 3600.0).eq.0) ! MN: 86400.0 --> 3600.0

  smobs= LVT_rc%udef
  smobs_2d = LVT_rc%udef

  if(SMOPSsmobs(source)%startmode.or.alarmCheck.or.&  
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false.      
     SMOPSsmobs(source)%startmode = .false. 

     if(LVT_rc%smoothObs.eq.1) then 
        smobs_av = 0.0
        count_av = 0

        call ESMF_TimeSet(cTime, yy = LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h  = LVT_rc%dhr(source), &
             m  = LVT_rc%dmn(source), & 
             s  = LVT_rc%dss(source), &
             calendar = LVT_calendar, & 
             rc = status)
        call LVT_verify(status)
        
        stTime = cTime - LVT_obsSmTwL
        enTime = cTime + LVT_obsSmTwL

        cTime = stTime
        do while (cTime.le.enTime) 
           call ESMF_TimeGet(cTime, yy = yr, mm=mo, dd=da,&
                h = hr, m = mn, s =ss, calendar=LVT_calendar, &
                rc=status)
           call LVT_verify(status)
           
           call create_SMOPSsm_filename(SMOPSsmobs(source)%odir, &
                SMOPSsmobs(source)%useRealtime , yr, mo, da, hr, fname)
           
           inquire(file=trim(fname),exist=file_exists)
           if(file_exists) then
              smobs = LVT_rc%udef
              write(LVT_logunit,*) '[INFO] Reading ',trim(fname)
              call read_SMOPS_data(source, fname,smobs)
           endif
           
           cTime = cTime + LVT_obsSmTwI
           
           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 if(smobs(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
                    smobs_av(c,r) = smobs_av(c,r) + smobs(c+(r-1)*LVT_rc%lnc)
                    count_av(c,r) = count_av(c,r) + 1
                 endif
              enddo
           enddo
        enddo
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(count_av(c,r).gt.0) then 
                 smobs_2d(c,r) = smobs_av(c,r)/count_av(c,r)
              else
                 smobs_2d(c,r) = LVT_rc%udef
              endif
           enddo
        enddo
        
     else
        call create_SMOPSsm_filename(SMOPSsmobs(source)%odir, &
             SMOPSsmobs(source)%useRealtime , &
             LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), &
             LVT_rc%dhr(source), fname)
        
        inquire(file=trim(fname),exist=file_exists)
        if(file_exists) then
           
           write(LVT_logunit,*) '[INFO] Reading ',trim(fname)
           call read_SMOPS_data(source, fname,smobs)
        endif
        
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(smobs(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
                 smobs_2d(c,r) = smobs(c+(r-1)*LVT_rc%lnc)
              endif
           enddo
        enddo

     endif
  else
     smobs_2d = LVT_rc%udef     
  endif

  relsmc = LVT_rc%udef

  if(LVT_rc%curr_pass.eq.1) then 
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(smobs_2d(c,r).ne.LVT_rc%udef) then 
              if(smobs_2d(c,r).gt.SMOPSsmobs(source)%maxv(c,r)) then
                 SMOPSsmobs(source)%maxv(c,r) = smobs_2d(c,r)
              endif
              
              if(smobs_2d(c,r).lt.SMOPSsmobs(source)%minv(c,r)) then
                 SMOPSsmobs(source)%minv(c,r) = smobs_2d(c,r)
              endif
           endif
        enddo
     enddo
  elseif(LVT_rc%curr_pass.eq.2) then 

     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc

           if(smobs_2d(c,r).ne.LVT_rc%udef.and.&
                (SMOPSsmobs(source)%maxv(c,r).ne.SMOPSsmobs(source)%minv(c,r))) then 
              relsmc(c,r) = (smobs_2d(c,r) - SMOPSsmobs(source)%minv(c,r))/&
                   (SMOPSsmobs(source)%maxv(c,r) - SMOPSsmobs(source)%minv(c,r))
           else
              relsmc(c,r) = LVT_rc%udef
           endif
              
        enddo
     enddo
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_RELSMC, source, &
       relsmc,vlevel=1,units="%")

  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, &
       smobs_2d,vlevel=1,units="m3/m3")

!  open(100,file='test.bin',form='unformatted')
!  write(100) smobs_2d
!  close(100)
!  stop

end subroutine readSMOPSsmObs


!BOP
! 
! !ROUTINE: read_SMOPS_data
! \label(read_SMOPS_data)
!
! !INTERFACE:
subroutine read_SMOPS_data(source, fname, smobs_ip)
! 
! !USES:   
#if (defined USE_GRIBAPI)
   use grib_api
#endif
  use LVT_coreMod,  only : LVT_rc
  use LVT_logMod
  use map_utils,    only : latlon_to_ij
  use SMOPSsm_obsMod, only : SMOPSsmobs
  use LVT_timeMgrMod,   only : LVT_date2time

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
!  This subroutine reads the SMOPS grib2 file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  the estimated error is above a predefined threshold (the recommeded
!  value is 5%). 
! MN: NOTE THE FOLLOWING INFORMATION IS VALID FOR SMOPS V1.3 AND V2 GRIB ID CHANGE IN VERSION 3
! The information from gribtab is: 
!{ 2, 0, 7, 1, 3, 210, "BLENDEDSM", "NOAA Blended Liquid Volumetric Soil Moisture (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 211, "AMSRESM", "NOAA AMSR-E Liquid Volumetric Soil Moisture (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 212, "SMOSSM", "SMOS Liquid Volumetric Soil Moisture (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 213, "ASCATSM", "ASCAT Liquid Volumetric Soil Moisture (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 214, "WINDSATSM", "NOAA WindSat Liquid Volumetric Soil Moisture (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 215, "RESSM1", "Reserved Liquid Volumetric Soil Moisture 1 (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 216, "RESSM2", "Reserved Liquid Volumetric Soil Moisture 2 (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 217, "BLENDEDSMHR", "Observation Hour for BLENDEDSM", "-"},
!{ 2, 0, 7, 1, 3, 218, "BLENDEDSMMN", "Observation Minute for BLENDEDSM", "-"},
!{ 2, 0, 7, 1, 3, 219, "AMSRESMHR", "Observation Hour for AMSRESM", "-"},
!{ 2, 0, 7, 1, 3, 220, "AMSRESMMN", "Observation Minute for AMSRESM", "-"},
!{ 2, 0, 7, 1, 3, 221, "SMOSSMHR", "Observation Hour for SMOSSM", "-"},
!{ 2, 0, 7, 1, 3, 222, "SMOSSMMN", "Observation Minute for SMOSSM", "-"},
!{ 2, 0, 7, 1, 3, 223, "ASCATSMHR", "Observation Hour for ASCATSM", "-"},
!{ 2, 0, 7, 1, 3, 224, "ASCATSMMN", "Observation Minute for ASCATSM", "-"},
!{ 2, 0, 7, 1, 3, 225, "WINDSATSMHR", "Observation Hour for WINDSATSM", "-"},
!{ 2, 0, 7, 1, 3, 226, "WINDSATSMMN", "Observation Minute for WINDSATSM", "-"},
!{ 2, 0, 7, 1, 3, 227, "RESSM1HR", "Observation Hour for RESSM1", "-"},
!{ 2, 0, 7, 1, 3, 228, "RESSM1MN", "Observation Minute for RESSM1", "-"},
!{ 2, 0, 7, 1, 3, 229, "RESSM2HR", "Observation Hour for RESSM2", "-"},
!{ 2, 0, 7, 1, 3, 230, "RESSM2MN", "Observation Minute for RESSM2", "-"},
!{ 2, 0, 7, 1, 3, 231, "BLENDEDSMQA", "NOAA Blended Liquid Volumetric Soil Moisture QA", "-"},
!{ 2, 0, 7, 1, 3, 232, "AMSRESMQA", "NOAA AMSR-E Liquid Volumetric Soil Moisture QA", "-"},
!{ 2, 0, 7, 1, 3, 233, "SMOSSMQA", "SMOS Liquid Volumetric Soil Moisture QA", "-"},
!{ 2, 0, 7, 1, 3, 234, "ASCATSMQA", "ASCAT Liquid Volumetric Soil Moisture QA", "-"},
!{ 2, 0, 7, 1, 3, 235, "WINDSATSMQA", "NOAA WindSat Liquid Volumetric Soil Moisture QA", "-"},
!{ 2, 0, 7, 1, 3, 236, "RESSM1QA", "Reserved Liquid Volumetric Soil Moisture 1 QA", "-"},
!{ 2, 0, 7, 1, 3, 237, "RESSM2QA", "Reserved Liquid Volumetric Soil Moisture 2 QA", "-"},

!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SMOPS AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LVT domain
! \end{description}
!
!
!EOP

! !DESCRIPTION:
!  This subroutine reads the SMOPS grib2 file and applies the data
!  quality flags to filter the data.
!
!  Quality flags are defined in:
!      NOAA NESDIS
!      CENTER FOR SATELLITE APPLICATIONS AND RESEARCH
!
!      SOIL MOISTURE OPERATIONAL PRODUCT SYSTEM (SMOPS)
!
!      ALGORITHM THEORETICAL BASIS DOCUMENT
!      Version 4.0
!
!  Found at http://www.ospo.noaa.gov/Products/land/smops/documents.html
!
!  The SMOPS QA flags are 16-bit (2-byte) integers, with the least
!  significant byte referred to as byte1 and the most significant byte
!  referred to as byte2.
!
!  bits: 15 | 14 | 13 | 12 | 11 | 10 | 9 | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 | 0
!                   byte2                    |           byte1
!
!  The arguments are:
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SMOPS AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LVT domain
! \end{description}
!
!EOP
   integer        :: param_ASCAT_A, param_ASCAT_A_qa
   integer        :: param_ASCAT_B, param_ASCAT_B_qa
   integer        :: param_SMOS, param_SMOS_qa
   integer        :: param_AMSR2, param_AMSR2_qa
   integer        :: param_SMAP, param_SMAP_qa
   integer*1,parameter :: err_threshold = 5 ! in percent
   integer*1,parameter :: AMSR2_accept = b'00000001'
   integer*2,parameter :: SMOS_accept1 = b'0000000000000000'
   integer*2,parameter :: SMOS_accept2 = b'0000000000000001'
   integer*2,parameter :: SMOS_accept3 = b'0000000000001000'
   integer*2,parameter :: SMOS_accept4 = b'0000000000001001'
   integer*2,parameter :: SMOS_accept5 = b'0000000010000000'
   integer*2,parameter :: SMOS_accept6 = b'0000000010000001'
   integer*2,parameter :: SMOS_accept7 = b'0000000010001000'
   integer*2,parameter :: SMOS_accept8 = b'0000000010001001'


!  INTEGER*2, PARAMETER :: FF = 255
!  real,    parameter  :: err_threshold = 5 ! in percent
!  real,    parameter  :: err_threshold = 100 ! in percent
! pre-3.0 versions
!  integer, parameter  :: param_ASCAT= 213, param_ASCAT_qa = 234
!  integer, parameter  :: param_windsat = 214, param_windsat_qa = 235
!  integer, parameter  :: param_SMOS = 212, param_SMOS_qa = 233

!  integer, parameter  :: param_ASCAT= 213, param_ASCAT_qa = 234
!  integer, parameter  :: param_windsat = 214, param_windsat_qa = 235
!  integer, parameter  :: param_SMOS = 212, param_SMOS_qa = 242
!SMAPSM
!  integer, parameter  :: param_SMAP = 218, param_SMAP_qa = 248
!NRTSMAP
!  integer, parameter  :: param_SMAP = 217, param_SMAP_qa = 247

!  real                :: sm_ascat(SMOPSsmobs(source)%rtsmopsnc*&
!       SMOPSsmobs(source)%rtsmopsnr)
!  real                :: sm_ascat_t(SMOPSsmobs(source)%rtsmopsnc*&
!  SMOPSsmobs(source)%rtsmopsnr)
!  real                :: sm_ascat_qa(SMOPSsmobs(source)%rtsmopsnc*&
!       SMOPSsmobs(source)%rtsmopsnr)
!  integer*2           :: sm_ascat_qa_t(SMOPSsmobs(source)%rtsmopsnc*&
!  SMOPSsmobs(source)%rtsmopsnr)
   real           :: sm_ASCAT_A(SMOPSsmobs(source)%smopsnc*&
      SMOPSsmobs(source)%smopsnr)
   real           :: sm_ASCAT_A_t(SMOPSsmobs(source)%smopsnc*&
      SMOPSsmobs(source)%smopsnr)
   real           :: sm_ASCAT_A_qa(SMOPSsmobs(source)%smopsnc*&
      SMOPSsmobs(source)%smopsnr)
   integer*2      :: sm_ASCAT_A_qa_t(SMOPSsmobs(source)%smopsnc*&
      SMOPSsmobs(source)%smopsnr)
   real           :: sm_ASCAT_B(SMOPSsmobs(source)%smopsnc*&
      SMOPSsmobs(source)%smopsnr)
   real           :: sm_ASCAT_B_t(SMOPSsmobs(source)%smopsnc*&
      SMOPSsmobs(source)%smopsnr)
   real           :: sm_ASCAT_B_qa(SMOPSsmobs(source)%smopsnc*&
      SMOPSsmobs(source)%smopsnr)
   integer*2      :: sm_ASCAT_B_qa_t(SMOPSsmobs(source)%smopsnc*&
      SMOPSsmobs(source)%smopsnr)

!  real                :: sm_windsat(SMOPSsmobs(source)%rtsmopsnc*&
!       SMOPSsmobs(source)%rtsmopsnr)
!  real                :: sm_windsat_t(SMOPSsmobs(source)%rtsmopsnc*&
!  SMOPSsmobs(source)%rtsmopsnr)
!  real                :: sm_windsat_qa(SMOPSsmobs(source)%rtsmopsnc*&
!       SMOPSsmobs(source)%rtsmopsnr)
!  integer*2           :: sm_windsat_qa_t(SMOPSsmobs(source)%rtsmopsnc*&
!  SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_smos(SMOPSsmobs(source)%smopsnc*&
       SMOPSsmobs(source)%smopsnr)
  real                :: sm_smos_t(SMOPSsmobs(source)%smopsnc*&
  SMOPSsmobs(source)%smopsnr)
  real                :: sm_smos_qa(SMOPSsmobs(source)%smopsnc*&
       SMOPSsmobs(source)%smopsnr)
  integer*2           :: sm_smos_qa_t(SMOPSsmobs(source)%smopsnc*&
  SMOPSsmobs(source)%smopsnr)
  real               :: sm_amsr2(SMOPSsmobs(source)%smopsnc*& 
	SMOPSsmobs(source)%smopsnr)
  real               :: sm_amsr2_t(SMOPSsmobs(source)%smopsnc*& 
	SMOPSsmobs(source)%smopsnr)
  real               :: sm_amsr2_qa(SMOPSsmobs(source)%smopsnc*&
	SMOPSsmobs(source)%smopsnr)
  integer*2          :: sm_amsr2_qa_t(SMOPSsmobs(source)%smopsnc*&
	SMOPSsmobs(source)%smopsnr)
  real                :: sm_smap(SMOPSsmobs(source)%smopsnc*&
       SMOPSsmobs(source)%smopsnr)
  real                :: sm_smap_t(SMOPSsmobs(source)%smopsnc*&
  SMOPSsmobs(source)%smopsnr)
  real                :: sm_smap_qa(SMOPSsmobs(source)%smopsnc*&
       SMOPSsmobs(source)%smopsnr)
  integer*2           :: sm_smap_qa_t(SMOPSsmobs(source)%smopsnc*&
  SMOPSsmobs(source)%smopsnr)
  real                :: sm_data(SMOPSsmobs(source)%smopsnc*&
       SMOPSsmobs(source)%smopsnr)
  logical*1           :: sm_data_b(SMOPSsmobs(source)%smopsnc*&
       SMOPSsmobs(source)%smopsnr)
  logical*1           :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)

  integer             :: c,r,i,j,kk
  integer             :: ftn,iret,igrib,nvars, num_message
  integer             :: param_num
  logical             :: var_found_ascat
  logical             :: var_found_windsat
  logical             :: var_found_smos
  logical             :: var_found_smap
  integer*2      :: qavalue
  real                :: err, ql, qaflags
   logical        :: var_found_AMSR2
   integer        :: updoy,yr1,mo1,da1,hr1,mn1,ss1
   real           :: upgmt
   real*8         :: timenow
   logical        :: smDataNotAvailable

   integer        :: ix, jx,c_s, c_e, r_s, r_e

#if (defined USE_GRIBAPI)
   smDataNotAvailable = .false.
   ! Set QA values to NESDIS SMOPS undefined value.
   sm_ASCAT_A_qa_t = 9999
   sm_ASCAT_B_qa_t = 9999
   sm_smos_qa_t    = 9999
   sm_amsr2_qa_t   = 9999
   sm_smap_qa_t    = 9999

   if ( SMOPSsmobs(source)%version == '1.3' ) then
      timenow = SMOPSsmobs(source)%version2_time - 1.0
   elseif ( SMOPSsmobs(source)%version == '2.0' ) then
      timenow = SMOPSsmobs(source)%version2_time
   elseif ( SMOPSsmobs(source)%version == '3.0' ) then
      timenow = SMOPSsmobs(source)%version3_time
   else
      yr1 = LVT_rc%yr
      mo1 = LVT_rc%mo
      da1 = LVT_rc%da
      hr1 = LVT_rc%hr
      mn1 = LVT_rc%mn
      ss1 = 0
      call LVT_date2time(timenow,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)
   endif

   if ( timenow < SMOPSsmobs(source)%version2_time ) then
      ! SMOPS version 1.3
      param_ASCAT_A = 213; param_ASCAT_A_qa = 234
      param_ASCAT_B = 214; param_ASCAT_B_qa = 235
      param_SMOS  = 212; param_SMOS_qa  = 233
      if ( SMOPSsmobs(source)%useSMOS.eq.1 ) then
         write(LVT_logunit,*) '[WARN] LVT does not process SMOS ' // &
            'from SMOPS version 1.3.'
         smDataNotAvailable = .true.
         smobs_ip = LVT_rc%udef
      endif
      if ( SMOPSsmobs(source)%useAMSR2.eq.1 ) then
         write(LVT_logunit,*) '[WARN] AMSR2 is not available ' // &
            'in SMOPS version 1.3'
         smDataNotAvailable = .true.
         smobs_ip = LVT_rc%udef
      endif
      if ( SMOPSsmobs(source)%useSMAP.eq.1 ) then
         write(LVT_logunit,*) '[WARN] SMAP is not available ' // &
            'in SMOPS version 1.3.'
         smDataNotAvailable = .true.
         smobs_ip = LVT_rc%udef
      endif
      write(LVT_logunit,*) '[INFO] Reading SMOPS dataset '//&
         'as SMOPS version 1.3'
   elseif ( timenow >= SMOPSsmobs(source)%version2_time .and. &
            timenow <  SMOPSsmobs(source)%version3_time ) then
      ! SMOPS version 2
      param_ASCAT_A = 213; param_ASCAT_A_qa = 234
      param_ASCAT_B = 214; param_ASCAT_B_qa = 235
      param_SMOS  = 212; param_SMOS_qa  = 233
      param_AMSR2 = 215; param_AMSR2_qa = 236
      if ( SMOPSsmobs(source)%useSMOS.eq.1 ) then
         write(LVT_logunit,*) '[WARN] LVT does not process SMOS ' // &
            'from SMOPS version 2.0.'
         smDataNotAvailable = .true.
         smobs_ip = LVT_rc%udef
      endif
      if ( SMOPSsmobs(source)%useAMSR2.eq.1 ) then
         write(LVT_logunit,*) '[WARN] LVT does not process AMSR2 ' // &
            'in SMOPS version 2.0.'
         smDataNotAvailable = .true.
         smobs_ip = LVT_rc%udef
      endif
      if(SMOPSsmobs(source)%useSMAP.eq.1 ) then
         write(LVT_logunit,*) '[WARN] SMAP is not availabe in SMOPS version 2'
         smDataNotAvailable = .true.
         smobs_ip = LVT_rc%udef
      endif
      write(LVT_logunit,*) '[INFO] Reading SMOPS dataset '//&
         'as SMOPS version 2.0'
   else ! ( timenow >= SMOPSsmobs(source)%version3_time ) then
      ! SMOPS version 3
      param_ASCAT_A = 213; param_ASCAT_A_qa = 243
      param_ASCAT_B = 214; param_ASCAT_B_qa = 244
      param_SMOS  = 212; param_SMOS_qa  = 242
      param_AMSR2 = 215; param_AMSR2_qa = 245
      param_SMAP  = 218; param_SMAP_qa  = 248
      write(LVT_logunit,*) '[INFO] Reading SMOPS dataset '//&
         'as SMOPS version 3.0'
   endif


   if ( smDataNotAvailable .eqv. .false. ) then


  call grib_open_file(ftn,trim(fname), 'r',iret)
  if(iret.ne.0) then 
     write(LVT_logunit,*) '[ERR] Could not open file: ',trim(fname)
     call LVT_endrun()
  endif
  call grib_multi_support_on

  do
     call grib_new_from_file(ftn,igrib,iret)
!     if ( iret /= 0 ) then
!        write(LVT_logunit,*) '[WARN] grib_new_from_file: failed in readSMOPSsmobs', trim(fname)
!        call grib_release(igrib,iret)
!        call LVT_verify(iret, 'error in grib_release in readSMOPSsmObs')
!        cycle
!     endif

     if ( iret == GRIB_END_OF_FILE ) then
        exit
     endif

!     call grib_count_in_file(ftn, num_message,iret)
!     print*, 'num_message', num_message
!     call grib_index_get_size(ftn, 'number',num_message,iret)
!     print*, 'num_message', num_message     
!     if (  num_message .lt.2 ) then
!        write(LVT_logunit,*) '[WARN] grib file is corrupted', trim(fname)
!        exit
!     endif

     call grib_get(igrib, 'parameterNumber',param_num, iret)
     if ( iret /= 0 ) then
        write(LVT_logunit,*) '[WARN] grib_get: parameterNumber failed in readSMOPSsmobs', trim(fname)
        cycle
     endif
     !call LVT_verify(iret, &
     !     'grib_get: parameterNumber failed in readSMOPSsmobs')
     var_found_ascat = .false. 
     if(SMOPSsmobs(source)%useASCAT.eq.1) then 
        if(param_num.eq.param_ASCAT_A) then 
           var_found_ascat = .true.
        endif
     endif
     
     if(var_found_ascat) then
        call grib_get(igrib, 'values',sm_ASCAT_A,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%smopsnr
           do c=1,SMOPSsmobs(source)%smopsnc
              sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                   sm_ASCAT_A(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%smopsnc)

           enddo
        enddo     
        
     endif

     var_found_ascat = .false. 
     if(SMOPSsmobs(source)%useASCAT.eq.1) then 
        if(param_num.eq.param_ASCAT_A_qa) then 
           var_found_ascat = .true.
        endif
     endif
 
!if (param_num .eq. 334 )then
!print * , param_num
!endif  
    
     if(var_found_ascat) then
        call grib_get(igrib, 'values',sm_ASCAT_A_qa,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%smopsnr
           do c=1,SMOPSsmobs(source)%smopsnc
              sm_ASCAT_A_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                   INT(sm_ASCAT_A_qa(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%smopsnc))
           enddo
        enddo       
     endif


     var_found_ascat = .false. 
     if(SMOPSsmobs(source)%useASCAT.eq.1) then 
        if(param_num.eq.param_ASCAT_B) then 
           var_found_ascat = .true.
        endif
     endif
     
     if(var_found_ascat) then
        call grib_get(igrib, 'values',sm_ASCAT_B,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%smopsnr
           do c=1,SMOPSsmobs(source)%smopsnc
              sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                   sm_ASCAT_B(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%smopsnc)

           enddo
        enddo     
        
     endif

     var_found_ascat = .false. 
     if(SMOPSsmobs(source)%useASCAT.eq.1) then 
        if(param_num.eq.param_ASCAT_B_qa) then 
           var_found_ascat = .true.
        endif
     endif
     
     if(var_found_ascat) then
        call grib_get(igrib, 'values',sm_ASCAT_B_qa,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%smopsnr
           do c=1,SMOPSsmobs(source)%smopsnc
              sm_ASCAT_B_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                   INT(sm_ASCAT_B_qa(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%smopsnc))
           enddo
        enddo       
     endif

!windsat
# if 0  
     var_found_windsat = .false. 
     if(SMOPSsmobs(source)%useWINDSAT.eq.1) then 
        if(param_num.eq.param_windsat) then 
           var_found_windsat = .true.
        endif
     endif
     
     if(var_found_windsat) then
        call grib_get(igrib, 'values',sm_windsat,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%smopsnr
           do c=1,SMOPSsmobs(source)%smopsnc
              sm_windsat_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                   sm_windsat(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%smopsnc)
           enddo
        enddo     
        
     endif

     var_found_windsat = .false. 
     if(SMOPSsmobs(source)%useWINDSAT.eq.1) then 
        if(param_num.eq.param_windsat_qa) then 
           var_found_windsat = .true.
        endif
     endif
     
     if(var_found_windsat) then
        call grib_get(igrib, 'values',sm_windsat_qa,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%smopsnr
           do c=1,SMOPSsmobs(source)%smopsnc
              sm_windsat_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                   INT(sm_windsat_qa(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%smopsnc))
           enddo
        enddo       
     endif
# endif 

!smos
     var_found_smos = .false. 
     if(SMOPSsmobs(source)%useSMOS.eq.1) then 
        if(param_num.eq.param_smos) then 
           var_found_smos = .true.
        endif
     endif
     
     if(var_found_smos) then
        call grib_get(igrib, 'values',sm_smos,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%smopsnr
           do c=1,SMOPSsmobs(source)%smopsnc
              sm_smos_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                   sm_smos(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%smopsnc)
           enddo
        enddo     
        
     endif

     var_found_smos = .false. 
     if(SMOPSsmobs(source)%useSMOS.eq.1) then 
        if(param_num.eq.param_smos_qa) then 
           var_found_smos = .true.
        endif
     endif
     
     if(var_found_smos) then
        call grib_get(igrib, 'values',sm_smos_qa,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%smopsnr
           do c=1,SMOPSsmobs(source)%smopsnc
              sm_smos_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                   INT(sm_smos_qa(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%smopsnc))
           enddo
        enddo       
     endif


        ! AMSR2
         var_found_amsr2 = .false.
         if(SMOPSsmobs(source)%useAMSR2.eq.1) then
            if(param_num.eq.param_amsr2) then
               var_found_amsr2 = .true.
            endif
         endif

         if(var_found_amsr2) then
            call grib_get(igrib, 'values',sm_amsr2,iret)
            call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
            !print *, 'sm_amsr2' , sm_amsr2
            do r=1,SMOPSsmobs(source)%smopsnr
               do c=1,SMOPSsmobs(source)%smopsnc
                  sm_amsr2_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                     sm_amsr2(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(source)%smopsnc)
                  if (sm_amsr2_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) .GT. 0.1 .and. &
                     sm_amsr2_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) .LT. 1) then
                     write(101,'(I5, 2x, I5, 2x, I8, 2x, G0, 2x)'), &
                        c, r, c+(r-1)*SMOPSsmobs(source)%smopsnc,      &
                        sm_amsr2_t(c+(r-1)*SMOPSsmobs(source)%smopsnc)
                  endif
               enddo
            enddo
         endif

         var_found_amsr2 = .false.
         if(SMOPSsmobs(source)%useAMSR2.eq.1) then
            if(param_num.eq.param_amsr2_qa) then
               var_found_amsr2 = .true.
            endif
         endif

         if(var_found_amsr2) then
            call grib_get(igrib, 'values',sm_amsr2_qa,iret)
            call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
            !print *, 'sm_amsr2_qa' , sm_amsr2_qa
            do r=1,SMOPSsmobs(source)%smopsnr
               do c=1,SMOPSsmobs(source)%smopsnc
                  sm_amsr2_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                     int(sm_amsr2_qa(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(source)%smopsnc))
                  !write(102,'(I5, 2x,I5, 2x, I8, 2x, f11.8,2x)'), &
                  !      c, r ,c+(r-1)*SMOPSsmobs(source)%smopsnc,      &
                  !      sm_amsr2_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc)
               enddo
            enddo
         endif


!smap
     var_found_smap = .false. 
     if(SMOPSsmobs(source)%useSMAP.eq.1) then 
        if(param_num.eq.param_smap) then 
           var_found_smap = .true.
        endif
     endif
     
     if(var_found_smap) then
        call grib_get(igrib, 'values',sm_smap,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%smopsnr
           do c=1,SMOPSsmobs(source)%smopsnc
              sm_smap_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                   sm_smap(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%smopsnc)
           enddo
        enddo     
        
     endif

     var_found_smap = .false. 
     if(SMOPSsmobs(source)%useSMAP.eq.1) then 
        if(param_num.eq.param_smap_qa) then 
           var_found_smap = .true.
        endif
     endif
     
     if(var_found_smap) then
        call grib_get(igrib, 'values',sm_smap_qa,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%smopsnr
           do c=1,SMOPSsmobs(source)%smopsnc
              sm_smap_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = &
                   INT(sm_smap_qa(c+((SMOPSsmobs(source)%smopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%smopsnc))
           enddo
        enddo       
     endif
     
     call grib_release(igrib,iret)
     call LVT_verify(iret, 'error in grib_release in readSMOPSsmObs')
  enddo

  call grib_close_file(ftn)

      ! Table 3.6.1 – SMOPS soil moisture product Quality Assessment (QA) bits.
      ! (d) ASCAT Soil Moisture Product QA
      !
      ! Byte |  Description
      ! -----------------------------------------------------------------------
      ! 0    |  Estimated Error in Soil Moisture. (Integer. Scale factor: 0.01)
      ! 1    |  Soil Moisture Quality (Integer, Scale factor: 0.01)
      !
      ! The retrievals are rejected when the estimated error is above
      ! a predefined threshold (the recommeded value is 5%).
      !
      ! Technically speaking, err_threshold should be 0.05 and
      ! err should be scaled by 0.01.  But below is consistent with the NESDIS
      ! documentation.
      !
      ! Note that I am assuming that Byte 0 above refers to the least
      ! significant byte and Byte 1 above refers to the most significant byte.
      ! Meaning estimated error is get_byte1, and quality is get_byte2.
  if(SMOPSsmobs(source)%useASCAT.eq.1) then 
     do r=1, SMOPSsmobs(source)%smopsnr
        do c=1, SMOPSsmobs(source)%smopsnc
           qavalue = sm_ASCAT_A_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc)
               if ( qavalue .ne. 9999) then
           !if(sm_ascat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc).ne.9999) then 
              !estimated error
              !err = ISHFT(sm_ascat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),-8)
              err = get_byte1(qavalue)
              !quality flag - not used currently
              !ql = IAND(sm_ascat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),FF)
              ql = get_byte2(qavalue)
              
              if(err.lt.err_threshold) then 
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .true. 
              else
                 sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false.
              endif
              if(sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(source)%smopsnc).lt.0.001) then
                 sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false.
              endif
              if(sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(source)%smopsnc).gt.1) then
                 sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
                !print*, 'ASCAT A', c, r ,sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(source)%smopsnc)
              endif

           else
              sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
              sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false. 
           endif
        enddo
     enddo

         do r=1, SMOPSsmobs(source)%smopsnr
            do c=1, SMOPSsmobs(source)%smopsnc
               qavalue = sm_ASCAT_B_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc)
               if ( qavalue .ne. 9999) then
                  !estimated error
                  err = get_byte1(qavalue)
                  !quality flag - not used currently
                  ql = get_byte2(qavalue)

                  if(err >= err_threshold) then
                     sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
                  endif
                  if(sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(source)%smopsnc).lt.0.001) then
                     sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
                  endif

                  if(sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(source)%smopsnc).gt.1) then
                    sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
                   !print*, 'ASCAT B', c, r ,sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(source)%smopsnc)
                 endif
               else
                  sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
               endif
            enddo
         enddo

         sm_data = sm_ASCAT_A_t
         where ( sm_ASCAT_B_t /= LVT_rc%udef )
            sm_data = sm_ASCAT_B_t
            sm_data_b = .true.
         endwhere

!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
     call neighbor_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
          SMOPSsmobs(source)%smopsnc*SMOPSsmobs(source)%smopsnr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          SMOPSsmobs(source)%rlat, SMOPSsmobs(source)%rlon, &
          SMOPSsmobs(source)%n11, LVT_rc%udef, iret)
  endif

#if 0 
  if(SMOPSsmobs(source)%useWINDSAT.eq.1) then 
     do r=1, SMOPSsmobs(source)%smopsnr
        do c=1, SMOPSsmobs(source)%smopsnc
           if(sm_windsat_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc).ne.9999) then 
              !estimated error
              err = ISHFT(sm_windsat_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc),-8)
              !quality flag - not used currently
              ql = IAND(sm_windsat_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc),FF)
              
              if(err.lt.err_threshold) then 
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .true. 
              else
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false.
                 sm_windsat_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
              endif
           else
              sm_windsat_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
              sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false. 
           endif
        enddo
     enddo
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
     call neighbor_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_windsat_t, smobs_b_ip, smobs_ip, &
          SMOPSsmobs(source)%smopsnc*SMOPSsmobs(source)%smopsnr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          SMOPSsmobs(source)%rlat, SMOPSsmobs(source)%rlon, &
          SMOPSsmobs(source)%n11, LVT_rc%udef, iret)
  endif
#endif 

      ! Table 3.6.1 – SMOPS soil moisture product Quality Assessment (QA) bits.
      ! (c) SMOS Soil Moisture Product QA
      !
      ! Byte 1:
      !
      ! Bit |  Description
      ! --------------------------------------------------------
      ! 0   |  Spare bit
      ! 1   |  1 = RFI for H pol above threshold, 0 = otherwise
      ! 2   |  1 = RFI for V pol above threshold, 0 = otherwise
      ! 3   |  Spare bit
      ! 4   |  1 = No products are generated, 0 = otherwise
      ! 5   |  1 = Retrieval values outside range, 0 = otherwise
      ! 6   |  1 = High retrieval DQX, 0 = otherwise
      ! 7   |  1 = Poor fit quality, 0 = otherwise
      !
      ! Byte 2:
      !
      ! Bit |  Description
      ! -------------------------------------------------------------
      ! 0   |  1 = Presence of other than nominal soil; 0 = otherwise
      ! 1   |  1 = Rocks; 0 = not rocks
      ! 2   |  1 = Moderate or strong topography; 0 = otherwise
      ! 3   |  1 = Open water; 0 = not open water
      ! 4   |  1 = Snow; 0 = not snow
      ! 5   |  1 = Forest; 0 = not forest
      ! 6   |  1 = Flood risk; 0 = no flood risk
      ! 7   |  1 = Urban area; 0 = not urban area
      !
      !
      ! From bytes 1 and 2, we will reject an observation if
      !
      !     bit 1 is 1 or bit 2 is 1 or bit 4 is 1 or bit 5 is 1 or
      !     bit 6 is 1 or bit 7 is 1 or bit 8 is 1 or bit 9 is 1 or
      !     bit 10 is 1 or bit 11 is 1 or bit 12 is 1 or bit 13 is 1 or
      !     bit 14 is 1 or bit 15 is 1
      !
      ! Thus we will accept an observation only if
      !
      !     bit 0 is (0|1) and bit 1 is 0 and bit 2 is 0 and
      !     bit 3 is (0|1) and bit 1 is 0 and bit 2 is 0 and
      !     bit 4 is 0 and bit 5 is 0 and bit 6 is 0 and bit 4 is 0 and
      !     bit 5 is 0 and bit 6 is 0 and bit 7 is 0 and bit 8 is 0 and
      !     bit 9 is 0 and bit 10 is 0 and bit 11 is 0 and bit 12 is 0 and
      !     bit 13 is 0 and bit 14 is 0 and bit 15 is 0
      !
      ! I.e., accept when bytes 1 and 2 are either b'0000000000000000' or
      ! b'0000000000000001' or b'0000000000001000' or b'0000000000001001';
      ! otherwise reject.


  if(SMOPSsmobs(source)%useSMOS.eq.1) then 
     do r=1, SMOPSsmobs(source)%smopsnr
        do c=1, SMOPSsmobs(source)%smopsnc
           qavalue = sm_smos_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc)
           if ( qavalue .ne. 9999 ) then
              if ( qavalue == SMOS_accept1 .or. &
                  qavalue == SMOS_accept2 .or. &
                  qavalue == SMOS_accept3 .or. &
                  qavalue == SMOS_accept4 .or. &
                  qavalue == SMOS_accept5 .or. &
                  qavalue == SMOS_accept6 .or. &
                  qavalue == SMOS_accept7 .or. &
                  qavalue == SMOS_accept8 ) then
                  sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .true.
              else
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false.
                 sm_smos_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
              endif
           else
              sm_smos_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
              sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false. 
           endif
        enddo
     enddo
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
     call neighbor_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_smos_t, smobs_b_ip, smobs_ip, &
          SMOPSsmobs(source)%smopsnc*SMOPSsmobs(source)%smopsnr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          SMOPSsmobs(source)%rlat, SMOPSsmobs(source)%rlon, &
          SMOPSsmobs(source)%n11, LVT_rc%udef, iret)
  endif

!  open(100,file='test_inp1.bin',form='unformatted')
!  write(100) sm_smap_t
!  close(100)


     ! Table 3.6.1 – SMOPS soil moisture product Quality Assessment (QA) bits.
      ! (e) AMSR2 Soil Moisture Layer QA
      !
      ! Byte 1:
      !
      ! Bit |  Description
      ! -------------------------------------------------------------------
      ! 0   |  0 = overall quality is not good; 1 = overall quality is good
      ! 1   |  1 = retrieval attempted but quality is not good; 0 = otherwise
      ! 2   |  1 = retrieval attempted but unsuccessful due to input
      !     |      brightness temperature data quality; 0 = otherwise
      ! 3   |  1 = retrieval attempted but unsuccessful due to the quality
      !     |      of other input data; 0 = otherwise
      ! 4   |  1 = retrieval not attempted; 0 = retrieval attempted
      ! 5   |  0 = not cold desert; 1 = cold desert
      ! 6   |  0 = not snow or rain; 1 = snow or rain
      ! 7   |  0 = not frozen ground; 1 = frozen ground
      !
      ! Byte 2:
      !
      ! Bit |  Description
      ! ----------------------------------------------
      ! 0   |  1: 0 ≤ GVF < 0.1; 0: otherwise
      ! 1   |  1: 0.1 ≤ GVF < 0.2; 0: otherwise
      ! 2   |  1: 0.2 ≤ GVF < 0.3; 0: otherwise
      ! 3   |  1: 0.3 ≤ GVF < 0.4; 0: otherwise
      ! 4   |  1: 0.4 ≤ GVF < 0.5; 0: otherwise
      ! 5   |  1: 0.5 ≤ GVF; 0: otherwise
      ! 6   |  1: overall input TB quality is good;
      !     |  0: overall input TB quality is not good
      ! 7   |  1 = real time NDVI; 0 = NDVI climate
      !
      ! Note that we are NOT considering byte 2.
      !
      ! From byte 1, we will reject an observation if
      !
      !    bit 0 is 0 or bit 1 is 1 or bit 2 is 1 or bit 3 is 1 or
      !    bit 4 is 1 or bit 5 is 1 or bit 6 is 1 or bit 7 is 1
      !
      ! Thus we will accept an observation only if
      !
      !    bit 0 is 1 and bit 1 is 0 and bit 2 is 0 and bit 3 is 0 and
      !    bit 4 is 0 and bit 5 is 0 and bit 6 is 0 and bit 7 is 0
      !
      ! I.e., accept when byte 1 is b'00000001'; otherwise reject.
      if(SMOPSsmobs(source)%useAMSR2.eq.1) then
         do r=1, SMOPSsmobs(source)%smopsnr
            do c=1, SMOPSsmobs(source)%smopsnc
               qavalue = sm_amsr2_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc)
               if ( qavalue .ne. 9999 ) then
                  qaflags = get_byte1(qavalue)
                  if (qaflags == AMSR2_accept ) then
                     sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .true.
                  else
                     sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false.
                     sm_amsr2_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
                  endif
               else
                  sm_amsr2_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
                  sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false.
               endif

               if (sm_amsr2_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) .GT. 0.1 .and. &
                  sm_amsr2_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) .LT. 1) then
                  write(103,'(I5, 2x, I5, 2x, I8, 2x, F10.4, 2x)'), &
                     c, r ,c+(r-1)*SMOPSsmobs(source)%smopsnc,   &
                     sm_amsr2_t(c+(r-1)*SMOPSsmobs(source)%smopsnc)
               endif
            enddo
         enddo
         !--------------------------------------------------------------------------
         ! Interpolate to the LVT running domain
         !--------------------------------------------------------------------------       
     call neighbor_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_amsr2_t, smobs_b_ip, smobs_ip, &
          SMOPSsmobs(source)%smopsnc*SMOPSsmobs(source)%smopsnr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          SMOPSsmobs(source)%rlat, SMOPSsmobs(source)%rlon, &
          SMOPSsmobs(source)%n11, LVT_rc%udef, iret)
    endif


  if(SMOPSsmobs(source)%useSMAP.eq.1) then 
     do r=1, SMOPSsmobs(source)%smopsnr
        do c=1, SMOPSsmobs(source)%smopsnc
               qavalue = sm_smap_qa_t(c+(r-1)*SMOPSsmobs(source)%smopsnc)
               if ( qavalue .ne. 9999 ) then
                  !estimated error
                  err = get_byte1(qavalue)
                  !quality flag - not used currently
                  ql = get_byte2(qavalue)

                  if(err.lt.err_threshold) then
                     sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .true.
                  else
                     sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false.
                     sm_smap_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
                  endif
               else
                  sm_smap_t(c+(r-1)*SMOPSsmobs(source)%smopsnc) = LVT_rc%udef
                  sm_data_b(c+(r-1)*SMOPSsmobs(source)%smopsnc) = .false.
               endif
        enddo
     enddo

!  open(100,file='test_inp.bin',form='unformatted')
!  write(100) sm_smap_t
!  close(100)
!  stop
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 

     call neighbor_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_smap_t, smobs_b_ip, smobs_ip, &
          SMOPSsmobs(source)%smopsnc*SMOPSsmobs(source)%smopsnr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          SMOPSsmobs(source)%rlat, SMOPSsmobs(source)%rlon, &
          SMOPSsmobs(source)%n11, LVT_rc%udef, iret)
  endif

!  print*, 'writing '
!  open(100,file='test_ip.bin',form='unformatted')
!  write(100) smobs_ip
!  close(100)
!  stop
  endif
#endif

contains 

integer*1 function get_byte1(i)
   implicit none
   integer*2, intent(in) :: i
   integer*2 :: j
   ! This function expects a 16-bit integer as input.  It returns
   ! the least significant byte (as an 8-bit integer), referred
   ! to as byte1 in the NESDIS documention.
   !
   ! For example,
   ! i = b'0000001000000001' <--> 00000010|00000001
   !                         <--> byte2|byte1
   !                         <--> 0x0201
   ! Here byte1 is b'00000001'; byte2 is b'00000010'
   j = iand(i, z'00ff')
   get_byte1 = j
end function get_byte1

integer*1 function get_byte2(i)
   implicit none
   integer*2, intent(in) :: i
   integer*2 :: j
   ! This function expects a 16-bit integer as input.  It returns
   ! the most significant byte (as an 8-bit integer), referred
   ! to as byte2 in the NESDIS documention.
   !
   ! For example,
   ! i = b'0000001000000001' <--> 00000010|00000001
   !                         <--> byte2|byte1
   !                         <--> 0x0201
   ! Here byte1 is b'00000001'; byte2 is b'00000010'
   j = ishft(i, -8)
   get_byte2 = j
end function get_byte2
end subroutine read_SMOPS_data

!BOP
! !ROUTINE: create_SMOPSsm_filename
! \label{create_SMOPSsm_filename}
! 
! !INTERFACE: 
subroutine create_SMOPSsm_filename(ndir, useRT, yr, mo, da, hr, filename)
! !USES:   
  use LVT_timeMgrMod,   only : LVT_date2time
  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da, hr
  character (len=*) :: ndir
   real*8         :: timenow , naming3_time
   integer        :: updoy
   real           :: upgmt
 integer           :: useRT
!  logical, save           :: first_time=.true.
! 
! !DESCRIPTION: 
!  This subroutine creates the SMOPS filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the SMOPS soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated RT SMOPS filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
     
  call LVT_date2time(naming3_time,updoy,upgmt,2017,10,5,0,0,0)
  call LVT_date2time(timenow,updoy,upgmt,yr,mo,da,hr,0,0)

  if(useRT.eq.1) then 
        if ( timenow >= naming3_time) then
           filename = trim(ndir)//'/'//'/NPR_SMOPS_CMAP_D' &
              //trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//'.gr2'     
        else
           filename = trim(ndir)//'/smops_d' &
              //trim(fyr)//trim(fmo)//trim(fda)//'_s'//trim(fhr)//'0000_cness.gr2'
        endif
  else

    filename = trim(ndir)//'/'//trim(fyr)//'/NPR_SMOPS_CMAP_D' &
       //trim(fyr)//trim(fmo)//trim(fda)//'.gr2'
  endif  
end subroutine create_SMOPSsm_filename




