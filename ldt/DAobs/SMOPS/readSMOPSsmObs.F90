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
! !ROUTINE: readSMOPSsmObs
! \label{readSMOPSsmObs}
!
! !REVISION HISTORY:
!  8 May 2013: Sujay Kumar, Initial Specification
!  22 Sep 2017: Mahdi Navari, filter added to screen the the ASCAT data
!               in the arid and semi-arid regions
!  9 Oct 2017: Mahdi Navari, added capability to read AMSR2 and SMAP
!
! !INTERFACE:
subroutine readSMOPSsmObs(n)
! !USES:
  use ESMF
  use LDT_coreMod,      only : LDT_rc
  use LDT_timeMgrMod,   only : LDT_get_julss
  use LDT_logMod,       only : LDT_logunit, LDT_getNextUnitNumber, &
                               LDT_releaseUnitNumber
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use SMOPSsm_obsMod,   only : SMOPSsmobs
  use map_utils

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
!
! !DESCRIPTION:
!
! This subroutine provides the data reader for the SMOPS
! soil moisture retrieval product. The data has many layers
! and the reader provides the options to select the
! desired layer(s).
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character(len=LDT_CONST_PATH_LEN)     :: fname
  real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations.
!-----------------------------------------------------------------------

  SMOPSsmobs(n)%smobs = LDT_rc%udef
  smobs = LDT_rc%udef

  call create_SMOPSsm_filename(SMOPSsmobs(n)%odir, &
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)

  inquire(file=trim(fname),exist=file_exists)
  if(file_exists) then

     write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
     call read_SMOPS_data(n, fname,smobs)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(smobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then
              SMOPSsmobs(n)%smobs(c,r) = smobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  else
     SMOPSsmobs(n)%smobs = LDT_rc%udef
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%soilmoist_obs,&
       SMOPSsmobs(n)%smobs,vlevel=1)

end subroutine readSMOPSsmObs


!BOP
!
! !ROUTINE: read_SMOPS_data
! \label(read_SMOPS_data)
!
! !INTERFACE:
subroutine read_SMOPS_data(n, fname, smobs_ip)
!
! !USES:
#if (defined USE_GRIBAPI)
   use grib_api
#endif
   use LDT_coreMod,      only : LDT_rc, LDT_domain
   use LDT_timeMgrMod,   only : LDT_date2time
   use LDT_logMod
   use map_utils,        only : latlon_to_ij
   use SMOPSsm_obsMod,   only : SMOPSsmobs

   ! MN : added to read the greenness, land cover and land mask
   use LDT_gfracMod,     only : LDT_gfrac_struc
   use LDT_LMLCMod
   use LDT_paramDataMod, only : LDT_LSMparam_struc

   implicit none
!
! !ARGUMENTS:
!
   integer                       :: n
   character (len=*)             :: fname
   real                          :: smobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   integer                       :: search_rad
!
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
!  \item[smobs\_ip]    soil moisture data processed to the LDT domain
! \end{description}
!
!EOP
   integer        :: param_ASCAT_A, param_ASCAT_A_qa
   integer        :: param_ASCAT_B, param_ASCAT_B_qa
   integer        :: param_SMOS, param_SMOS_qa
   integer        :: param_AMSR2, param_AMSR2_qa
   integer        :: param_SMAP, param_SMAP_qa
   integer*1, parameter :: err_threshold = 5 ! in percent
   ! EMK 21 Jul 2021...ISO Fortran restricts use of binary integer constants
   ! to DATA statements and as a principal argument to an intrinsic function.
   ! And gfortran 10 and 11 rejects their use in parameters.
   ! So we replace with arrays here.
   ! integer*1,parameter :: AMSR2_accept = b'00000001'
   ! integer*2,parameter :: SMOS_accept1 = b'0000000000000000'
   ! integer*2,parameter :: SMOS_accept2 = b'0000000000000001'
   ! integer*2,parameter :: SMOS_accept3 = b'0000000000001000'
   ! integer*2,parameter :: SMOS_accept4 = b'0000000000001001'
   ! integer*2,parameter :: SMOS_accept5 = b'0000000010000000'
   ! integer*2,parameter :: SMOS_accept6 = b'0000000010000001'
   ! integer*2,parameter :: SMOS_accept7 = b'0000000010001000'
   ! integer*2,parameter :: SMOS_accept8 = b'0000000010001001'
   integer*1 :: AMSR2_accept(1)
   integer*2 :: SMOS_accept(8)
   real           :: sm_ASCAT_A(SMOPSsmobs(n)%smopsnc*&
      SMOPSsmobs(n)%smopsnr)
   real           :: sm_ASCAT_A_t(SMOPSsmobs(n)%smopsnc*&
      SMOPSsmobs(n)%smopsnr)
   real           :: sm_ASCAT_A_qa(SMOPSsmobs(n)%smopsnc*&
      SMOPSsmobs(n)%smopsnr)
   integer*2      :: sm_ASCAT_A_qa_t(SMOPSsmobs(n)%smopsnc*&
      SMOPSsmobs(n)%smopsnr)
   real           :: sm_ASCAT_B(SMOPSsmobs(n)%smopsnc*&
      SMOPSsmobs(n)%smopsnr)
   real           :: sm_ASCAT_B_t(SMOPSsmobs(n)%smopsnc*&
      SMOPSsmobs(n)%smopsnr)
   real           :: sm_ASCAT_B_qa(SMOPSsmobs(n)%smopsnc*&
      SMOPSsmobs(n)%smopsnr)
   integer*2      :: sm_ASCAT_B_qa_t(SMOPSsmobs(n)%smopsnc*&
      SMOPSsmobs(n)%smopsnr)
   real           :: sm_smos(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   real           :: sm_smos_t(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   real           :: sm_smos_qa(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   integer*2      :: sm_smos_qa_t(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   real           :: sm_amsr2(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   real           :: sm_amsr2_t(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   real           :: sm_amsr2_qa(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   integer*2      :: sm_amsr2_qa_t(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   real           :: sm_smap(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   real           :: sm_smap_t(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   real           :: sm_smap_qa(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   integer*2      :: sm_smap_qa_t(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)

   real           :: sm_data(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   logical*1      :: sm_data_b(SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr)
   logical*1      :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))

   integer        :: c,r,i,j,kk
   integer        :: ftn,iret,igrib,nvars
   integer        :: param_num
   logical        :: var_found_ascat
   logical        :: var_found_smos
   logical        :: var_found_AMSR2
   logical        :: var_found_SMAP
   integer*2      :: qavalue
   integer*1      :: err, ql, qaflags
   integer        :: updoy,yr1,mo1,da1,hr1,mn1,ss1
   real           :: upgmt
   real*8         :: timenow
   logical        :: smDataNotAvailable

   integer        :: ix, jx,c_s, c_e, r_s, r_e

   ! EMK 21 Jul 2021...ISO Fortran restricts use of binary integer constants
   ! to DATA statements and as a principal argument to an intrinsic function.
   ! And gfortran 10 and 11 rejects their use in parameters.
   ! So we replace with arrays here.
   ! Old code for reference
   !integer*1,parameter :: AMSR2_accept = b'00000001'
   !integer*2,parameter :: SMOS_accept1 = b'0000000000000000'
   !integer*2,parameter :: SMOS_accept2 = b'0000000000000001'
   !integer*2,parameter :: SMOS_accept3 = b'0000000000001000'
   !integer*2,parameter :: SMOS_accept4 = b'0000000000001001'
   !integer*2,parameter :: SMOS_accept5 = b'0000000010000000'
   !integer*2,parameter :: SMOS_accept6 = b'0000000010000001'
   !integer*2,parameter :: SMOS_accept7 = b'0000000010001000'
   !integer*2,parameter :: SMOS_accept8 = b'0000000010001001'
   data AMSR2_accept / b'00000001' /
   data SMOS_accept / b'0000000000000000', &
        b'0000000000000001', &
        b'0000000000001000', &
        b'0000000000001001', &
        b'0000000010000000', &
        b'0000000010000001', &
        b'0000000010001000', &
        b'0000000010001001' /

#if (defined USE_GRIBAPI)
   smDataNotAvailable = .false.

   ! Set QA values to NESDIS SMOPS undefined value.
   sm_ASCAT_A_qa_t = 9999
   sm_ASCAT_B_qa_t = 9999
   sm_smos_qa_t    = 9999
   sm_amsr2_qa_t   = 9999
   sm_smap_qa_t    = 9999

   if ( SMOPSsmobs(n)%version == '1.3' ) then
      timenow = SMOPSsmobs(n)%version2_time - 1.0
   elseif ( SMOPSsmobs(n)%version == '2.0' ) then
      timenow = SMOPSsmobs(n)%version2_time
   elseif ( SMOPSsmobs(n)%version == '3.0' ) then
      timenow = SMOPSsmobs(n)%version3_time
   else
      yr1 = LDT_rc%yr
      mo1 = LDT_rc%mo
      da1 = LDT_rc%da
      hr1 = LDT_rc%hr
      mn1 = LDT_rc%mn
      ss1 = 0
      call LDT_date2time(timenow,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)
   endif

   if ( timenow < SMOPSsmobs(n)%version2_time ) then
      ! SMOPS version 1.3
      param_ASCAT_A = 213; param_ASCAT_A_qa = 234
      param_ASCAT_B = 214; param_ASCAT_B_qa = 235
      param_SMOS  = 212; param_SMOS_qa  = 233
      if ( SMOPSsmobs(n)%useSMOS.eq.1 ) then
         write(LDT_logunit,*) '[Warning] LDT does not process SMOS ' // &
            'from SMOPS version 1.3.'
         smDataNotAvailable = .true.
         smobs_ip = LDT_rc%udef
      endif
      if ( SMOPSsmobs(n)%useAMSR2.eq.1 ) then
         write(LDT_logunit,*) '[Warning] AMSR2 is not available ' // &
            'in SMOPS version 1.3'
         smDataNotAvailable = .true.
         smobs_ip = LDT_rc%udef
      endif
      if ( SMOPSsmobs(n)%useSMAP.eq.1 ) then
         write(LDT_logunit,*) '[Warning] SMAP is not available ' // &
            'in SMOPS version 1.3.'
         smDataNotAvailable = .true.
         smobs_ip = LDT_rc%udef
      endif
      write(LDT_logunit,*) '[MSG] Reading SMOPS dataset '//&
         'as SMOPS version 1.3'
   elseif ( timenow >= SMOPSsmobs(n)%version2_time .and. &
            timenow <  SMOPSsmobs(n)%version3_time ) then
      ! SMOPS version 2
      param_ASCAT_A = 213; param_ASCAT_A_qa = 234
      param_ASCAT_B = 214; param_ASCAT_B_qa = 235
      param_SMOS  = 212; param_SMOS_qa  = 233
      param_AMSR2 = 215; param_AMSR2_qa = 236
      if ( SMOPSsmobs(n)%useSMOS.eq.1 ) then
         write(LDT_logunit,*) '[Warning] LDT does not process SMOS ' // &
            'from SMOPS version 2.0.'
         smDataNotAvailable = .true.
         smobs_ip = LDT_rc%udef
      endif
      if ( SMOPSsmobs(n)%useAMSR2.eq.1 ) then
         write(LDT_logunit,*) '[Warning] LDT does not process AMSR2 ' // &
            'in SMOPS version 2.0.'
         smDataNotAvailable = .true.
         smobs_ip = LDT_rc%udef
      endif
      if(SMOPSsmobs(n)%useSMAP.eq.1 ) then
         write(LDT_logunit,*) '[Warning] SMAP is not availabe in SMOPS version 2'
         smDataNotAvailable = .true.
         smobs_ip = LDT_rc%udef
      endif
      write(LDT_logunit,*) '[MSG] Reading SMOPS dataset '//&
         'as SMOPS version 2.0'
   else ! ( timenow >= SMOPSsmobs(n)%version3_time ) then
      ! SMOPS version 3
      param_ASCAT_A = 213; param_ASCAT_A_qa = 243
      param_ASCAT_B = 214; param_ASCAT_B_qa = 244
      param_SMOS  = 212; param_SMOS_qa  = 242
      param_AMSR2 = 215; param_AMSR2_qa = 245
      param_SMAP  = 218; param_SMAP_qa  = 248
      write(LDT_logunit,*) '[MSG] Reading SMOPS dataset '//&
         'as SMOPS version 3.0'
   endif



   if ( smDataNotAvailable .eqv. .false. ) then

      call grib_open_file(ftn,trim(fname),'r',iret)
      if(iret.ne.0) then
         write(LDT_logunit,*) 'Could not open file: ',trim(fname)
         flush(LDT_logunit)
         call LDT_endrun()
      endif
      call grib_multi_support_on

      do
         call grib_new_from_file(ftn,igrib,iret)

         if ( iret == GRIB_END_OF_FILE ) then
            exit
         endif

         call grib_get(igrib, 'parameterNumber',param_num, iret)
         call LDT_verify(iret, &
            'grib_get: parameterNumber failed in readSMOPSsmobs')

         var_found_ascat = .false.
         if(SMOPSsmobs(n)%useASCAT.eq.1) then
            if(param_num.eq.param_ASCAT_A) then
               var_found_ascat = .true.
            endif
         endif

         if(var_found_ascat) then
            call grib_get(igrib, 'values',sm_ASCAT_A,iret)
            call LDT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
            do r=1,SMOPSsmobs(n)%smopsnr
               do c=1,SMOPSsmobs(n)%smopsnc
                  sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = &
                     sm_ASCAT_A(c+((SMOPSsmobs(n)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(n)%smopsnc)
               enddo
            enddo
         endif

         var_found_ascat = .false.
         if(SMOPSsmobs(n)%useASCAT.eq.1) then
            if(param_num.eq.param_ASCAT_A_qa) then
               var_found_ascat = .true.
            endif
         endif

         if(var_found_ascat) then
            call grib_get(igrib, 'values',sm_ASCAT_A_qa,iret)
            call LDT_warning(iret,'error in grib_get:values in readSMOPSsmObs')

            do r=1,SMOPSsmobs(n)%smopsnr
               do c=1,SMOPSsmobs(n)%smopsnc
                  sm_ASCAT_A_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = &
                     int(sm_ASCAT_A_qa(c+((SMOPSsmobs(n)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(n)%smopsnc))
               enddo
            enddo
         endif

         var_found_ascat = .false.
         if(SMOPSsmobs(n)%useASCAT.eq.1) then
            if(param_num.eq.param_ASCAT_B) then
               var_found_ascat = .true.
            endif
         endif

         if(var_found_ascat) then
            call grib_get(igrib, 'values',sm_ASCAT_B,iret)
            call LDT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
            do r=1,SMOPSsmobs(n)%smopsnr
               do c=1,SMOPSsmobs(n)%smopsnc
                  sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = &
                     sm_ASCAT_B(c+((SMOPSsmobs(n)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(n)%smopsnc)
               enddo
            enddo
         endif

         var_found_ascat = .false.
         if(SMOPSsmobs(n)%useASCAT.eq.1) then
            if(param_num.eq.param_ASCAT_B_qa) then
               var_found_ascat = .true.
            endif
         endif

         if(var_found_ascat) then
            call grib_get(igrib, 'values',sm_ASCAT_B_qa,iret)
            call LDT_warning(iret,'error in grib_get:values in readSMOPSsmObs')

            do r=1,SMOPSsmobs(n)%smopsnr
               do c=1,SMOPSsmobs(n)%smopsnc
                  sm_ASCAT_B_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = &
                     int(sm_ASCAT_B_qa(c+((SMOPSsmobs(n)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(n)%smopsnc))
               enddo
            enddo
         endif

         !smos
         var_found_smos = .false.
         if(SMOPSsmobs(n)%useSMOS.eq.1) then
            if(param_num.eq.param_smos) then
               var_found_smos = .true.
            endif
         endif

         if(var_found_smos) then
            call grib_get(igrib, 'values',sm_smos,iret)
            call LDT_warning(iret,'error in grib_get:values in readSMOPSsmObs')

            do r=1,SMOPSsmobs(n)%smopsnr
               do c=1,SMOPSsmobs(n)%smopsnc
                  sm_smos_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = &
                     sm_smos(c+((SMOPSsmobs(n)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(n)%smopsnc)
               enddo
            enddo
         endif

         var_found_smos = .false.
         if(SMOPSsmobs(n)%useSMOS.eq.1) then
            if(param_num.eq.param_smos_qa) then
               var_found_smos = .true.
            endif
         endif

         if(var_found_smos) then
            call grib_get(igrib, 'values',sm_smos_qa,iret)
            call LDT_warning(iret,'error in grib_get:values in readSMOPSsmObs')

            do r=1,SMOPSsmobs(n)%smopsnr
               do c=1,SMOPSsmobs(n)%smopsnc
                  sm_smos_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = &
                     int(sm_smos_qa(c+((SMOPSsmobs(n)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(n)%smopsnc))
               enddo
            enddo
         endif

         ! AMSR2
         var_found_amsr2 = .false.
         if(SMOPSsmobs(n)%useAMSR2.eq.1) then
            if(param_num.eq.param_amsr2) then
               var_found_amsr2 = .true.
            endif
         endif

         if(var_found_amsr2) then
            call grib_get(igrib, 'values',sm_amsr2,iret)
            call LDT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
            !print *, 'sm_amsr2' , sm_amsr2
            do r=1,SMOPSsmobs(n)%smopsnr
               do c=1,SMOPSsmobs(n)%smopsnc
                  sm_amsr2_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = &
                     sm_amsr2(c+((SMOPSsmobs(n)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(n)%smopsnc)
                  if (sm_amsr2_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) .GT. 0.1 .and. &
                     sm_amsr2_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) .LT. 1) then
                     write(101,'(I5, 2x, I5, 2x, I8, 2x, G0, 2x)')  &
                        c, r, c+(r-1)*SMOPSsmobs(n)%smopsnc,      &
                        sm_amsr2_t(c+(r-1)*SMOPSsmobs(n)%smopsnc)
                  endif
               enddo
            enddo
         endif

         var_found_amsr2 = .false.
         if(SMOPSsmobs(n)%useAMSR2.eq.1) then
            if(param_num.eq.param_amsr2_qa) then
               var_found_amsr2 = .true.
            endif
         endif

         if(var_found_amsr2) then
            call grib_get(igrib, 'values',sm_amsr2_qa,iret)
            call LDT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
            !print *, 'sm_amsr2_qa' , sm_amsr2_qa
            do r=1,SMOPSsmobs(n)%smopsnr
               do c=1,SMOPSsmobs(n)%smopsnc
                  sm_amsr2_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = &
                     int(sm_amsr2_qa(c+((SMOPSsmobs(n)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(n)%smopsnc))
                  !write(102,'(I5, 2x,I5, 2x, I8, 2x, f11.8,2x)')  &
                  !      c, r ,c+(r-1)*SMOPSsmobs(n)%smopsnc,      &
                  !      sm_amsr2_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc)
               enddo
            enddo
         endif
         !SMAP
         var_found_smap = .false.
         if(SMOPSsmobs(n)%useSMAP.eq.1) then
            if(param_num.eq.param_smap) then
               var_found_smap = .true.
            endif
         endif

         if(var_found_smap) then
            call grib_get(igrib, 'values',sm_smap,iret)
            call LDT_warning(iret,'error in grib_get:values in readSMOPSsmObs')

            do r=1,SMOPSsmobs(n)%smopsnr
               do c=1,SMOPSsmobs(n)%smopsnc
                  sm_smap_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = &
                     sm_smap(c+((SMOPSsmobs(n)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(n)%smopsnc)
               enddo
            enddo

         endif

         var_found_smap = .false.
         if(SMOPSsmobs(n)%useSMAP.eq.1) then
            if(param_num.eq.param_smap_qa) then
               var_found_smap = .true.
            endif
         endif

         if(var_found_smap) then
            call grib_get(igrib, 'values',sm_smap_qa,iret)
            call LDT_warning(iret,'error in grib_get:values in readSMOPSsmObs')

            do r=1,SMOPSsmobs(n)%smopsnr
               do c=1,SMOPSsmobs(n)%smopsnc
                  sm_smap_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = &
                     int(sm_smap_qa(c+((SMOPSsmobs(n)%smopsnr-r+1)-1)*&
                     SMOPSsmobs(n)%smopsnc))
               enddo
            enddo
         endif

         call grib_release(igrib,iret)
         call LDT_verify(iret, 'error in grib_release in readSMOPSsmObs')
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
      if(SMOPSsmobs(n)%useASCAT.eq.1) then
         do r=1, SMOPSsmobs(n)%smopsnr
            do c=1, SMOPSsmobs(n)%smopsnc
               qavalue = sm_ASCAT_A_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc)
               if ( qavalue .ne. 9999 ) then
                  !estimated error
                  err = get_byte1(qavalue)
                  !quality flag - not used currently
                  ql = get_byte2(qavalue)

                  if(err.lt.err_threshold) then
                     sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .true.
                  else
                     sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                     sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .false.
                  endif
                  if(sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(n)%smopsnc).lt.0.001) then
                     sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                     sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .false.
                  endif
               else
                  sm_ASCAT_A_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                  sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .false.
               endif
            enddo
         enddo

         do r=1, SMOPSsmobs(n)%smopsnr
            do c=1, SMOPSsmobs(n)%smopsnc
               qavalue = sm_ASCAT_B_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc)
               if ( qavalue .ne. 9999 ) then
                  !estimated error
                  err = get_byte1(qavalue)
                  !quality flag - not used currently
                  ql = get_byte2(qavalue)

                  if(err >= err_threshold) then
                     sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                  endif
                  if(sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(n)%smopsnc).lt.0.001) then
                     sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                  endif
               else
                  sm_ASCAT_B_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
               endif
            enddo
         enddo

         sm_data = sm_ASCAT_A_t
         where ( sm_ASCAT_B_t /= LDT_rc%udef )
            sm_data = sm_ASCAT_B_t
            sm_data_b = .true.
         endwhere

         !--------------------------------------------------------------------------
         ! Interpolate to the LDT running domain
         !--------------------------------------------------------------------------
         call neighbor_interp(LDT_rc%gridDesc(n,:),&
            sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
            SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr, &
            LDT_rc%lnc(n)*LDT_rc%lnr(n), &
            LDT_domain(n)%lat, LDT_domain(n)%lon,&
            SMOPSsmobs(n)%n11, LDT_rc%udef, iret)

         do r = 1, LDT_rc%lnr(n)
            do c = 1, LDT_rc%lnc(n)
               ! MN : add condition to check the GVF and land cover
               !      screen the ASCAT if GVF is less than 0.1 and
               !      land cover type is Barren or Sparsely Vegetated and
               !      the value of the pixel is larger than 0.5
               if ( allocated(LDT_gfrac_struc(n)%gfrac%value) ) then
                  if ( maxval(LDT_gfrac_struc(n)%gfrac%value(c,r,:)) &
                     .le. 0.1 .and. &
                     LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%bareclass) &
                     .ge. 0.5  ) then
                     smobs_ip(c+(r-1)*LDT_rc%lnc(n)) = LDT_rc%udef
                  endif
               else
                  !              write(LDT_logunit, *) '[INFO] Enabling greenness fraction ' // &
                  !                                    'will improve the SMOPS ASCAT filtering'
                  if ( LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%bareclass) &
                     .ge. 0.5  ) then
                     smobs_ip(c+(r-1)*LDT_rc%lnc(n)) = LDT_rc%udef
                  endif
               endif
            enddo
         enddo
      endif



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
      !
      if(SMOPSsmobs(n)%useSMOS.eq.1) then
         do r=1, SMOPSsmobs(n)%smopsnr
            do c=1, SMOPSsmobs(n)%smopsnc
               qavalue = sm_smos_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc)
               if ( qavalue .ne. 9999 ) then
                  if ( qavalue == SMOS_accept(1) .or. &
                     qavalue == SMOS_accept(2) .or. &
                     qavalue == SMOS_accept(3) .or. &
                     qavalue == SMOS_accept(4) .or. &
                     qavalue == SMOS_accept(5) .or. &
                     qavalue == SMOS_accept(6) .or. &
                     qavalue == SMOS_accept(7) .or. &
                     qavalue == SMOS_accept(8) ) then
                     sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .true.
                  else
                     sm_smos_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                     sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .false.
                  endif
               else
                  sm_smos_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                  sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .false.
               endif
            enddo
         enddo
         !--------------------------------------------------------------------------
         ! Interpolate to the LDT running domain
         !--------------------------------------------------------------------------
         call neighbor_interp(LDT_rc%gridDesc(n,:),&
            sm_data_b, sm_smos_t, smobs_b_ip, smobs_ip, &
            SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr, &
            LDT_rc%lnc(n)*LDT_rc%lnr(n), &
            LDT_domain(n)%lat, LDT_domain(n)%lon,&
            SMOPSsmobs(n)%n11, LDT_rc%udef, iret)
      endif





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
      !
      if(SMOPSsmobs(n)%useAMSR2.eq.1) then
         do r=1, SMOPSsmobs(n)%smopsnr
            do c=1, SMOPSsmobs(n)%smopsnc
               qavalue = sm_amsr2_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc)
               if ( qavalue .ne. 9999 ) then
                  qaflags = get_byte1(qavalue)
                  if (qaflags == AMSR2_accept(1) ) then
                     sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .true.
                  else
                     sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .false.
                     sm_amsr2_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                  endif
               else
                  sm_amsr2_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                  sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .false.
               endif

               if (sm_amsr2_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) .GT. 0.1 .and. &
                  sm_amsr2_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) .LT. 1) then
                  write(103,'(I5, 2x, I5, 2x, I8, 2x, F10.4, 2x)')  &
                     c, r ,c+(r-1)*SMOPSsmobs(n)%smopsnc,   &
                     sm_amsr2_t(c+(r-1)*SMOPSsmobs(n)%smopsnc)
               endif
            enddo
         enddo
         !--------------------------------------------------------------------------
         ! Interpolate to the LDT running domain
         !--------------------------------------------------------------------------
         call neighbor_interp(LDT_rc%gridDesc(n,:),&
            sm_data_b, sm_amsr2_t, smobs_b_ip, smobs_ip, &
            SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr, &
            LDT_rc%lnc(n)*LDT_rc%lnr(n), &
            LDT_domain(n)%lat, LDT_domain(n)%lon,&
            SMOPSsmobs(n)%n11, LDT_rc%udef, iret)



         ! MN print
         do r=1, LDT_rc%lnr(n)
            do c=1, LDT_rc%lnc(n)
               if (smobs_ip(c+(r-1)*LDT_rc%lnc(n)) .GT. 0.1 .and. &
                  smobs_ip(c+(r-1)*LDT_rc%lnc(n)) .LT. 1) then
                  write(104,'(I5, 2x, I5, 2x, I8, 2x, F10.4 ,2x)')  &
                     c, r , c+(r-1)*LDT_rc%lnc(n),         &
                     smobs_ip(c+(r-1)*LDT_rc%lnc(n))
               endif
            enddo
         enddo




      endif


      ! SMAP
      ! Not yet documented.  Applying ASCAT QA logic.
      if(SMOPSsmobs(n)%useSMAP.eq.1) then
         do r=1, SMOPSsmobs(n)%smopsnr
            do c=1, SMOPSsmobs(n)%smopsnc
               qavalue = sm_smap_qa_t(c+(r-1)*SMOPSsmobs(n)%smopsnc)
               if ( qavalue .ne. 9999 ) then
                  !estimated error
                  err = get_byte1(qavalue)
                  !quality flag - not used currently
                  ql = get_byte2(qavalue)

                  if(err.lt.err_threshold) then
                     sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .true.
                  else
                     sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .false.
                     sm_smap_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                  endif
               else
                  sm_smap_t(c+(r-1)*SMOPSsmobs(n)%smopsnc) = LDT_rc%udef
                  sm_data_b(c+(r-1)*SMOPSsmobs(n)%smopsnc) = .false.
               endif
            enddo
         enddo
         !--------------------------------------------------------------------------
         ! Interpolate to the LDT running domain
         !--------------------------------------------------------------------------
         call neighbor_interp(LDT_rc%gridDesc(n,:),&
            sm_data_b, sm_smap_t, smobs_b_ip, smobs_ip, &
            SMOPSsmobs(n)%smopsnc*SMOPSsmobs(n)%smopsnr, &
            LDT_rc%lnc(n)*LDT_rc%lnr(n), &
            LDT_domain(n)%lat, LDT_domain(n)%lon,&
            SMOPSsmobs(n)%n11, LDT_rc%udef, iret)
      endif

      !------------------------------------------------------------------------
      !  Remove pixel close to open water
      !------------------------------------------------------------------------
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            if (smobs_ip(c+(r-1)*LDT_rc%lnc(n)) .ne. LDT_rc%udef) then 
               search_rad = nint(SMOPSsmobs(n)%search_radius)
               
               c_s = max(1,c-search_rad)
               c_e = min(LDT_rc%lnc(n),c+search_rad)
               
               r_s = max(1,r-search_rad)
               r_e = min(LDT_rc%lnr(n),r+search_rad)
               
               do ix=c_s,c_e
                  do jx=r_s,r_e
                     if( LDT_LSMparam_struc(n)%landmask%value(ix,jx,1) .eq. 0 ) then
                        smobs_ip(c+(r-1)*LDT_rc%lnc(n)) = LDT_rc%udef
                     end if
                  enddo
               enddo
            endif
         enddo
      enddo
      
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
subroutine create_SMOPSsm_filename(ndir, yr, mo,da, filename)
! !USES:

  implicit none
! !ARGUMENTS:
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
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
  character (len=2) :: fmo,fda

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  filename = trim(ndir)//'/'//trim(fyr)//'/NPR_SMOPS_CMAP_D' &
       //trim(fyr)//trim(fmo)//trim(fda)//'.gr2'

end subroutine create_SMOPSsm_filename
