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
! !ROUTINE: readSMAPEOPL_SMObs
! \label{readSMAPEOPL_SMObs}
! 
! !REVISION HISTORY: 
!  6 Jun 2022: Yonghwan Kwon, Initial Specification
! 
! !INTERFACE: 
subroutine readSMAPEOPL_SMObs(n)
! !USES: 
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use LDT_DAobsDataMod
  use SMAPEOPLSMobsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for 
! SMAP_E_OPL soil moisture data
!
!EOP

  real*8            :: timenow
  real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  character*100     :: fname
  integer           :: mn_ind
  integer           :: mn, ss
  integer           :: doy
  character*8       :: yyyymmdd
  character*2       :: hh
  character*200     :: list_files
  character*100     :: smap_filename(10)
  integer           :: i
  integer           :: ftn, ierr
  real              :: gmt

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations.
!-----------------------------------------------------------------------
  SMAPEOPLsmobs(n)%smobs = LDT_rc%udef
  smobs = LDT_rc%udef

  if (LDT_rc%ts .gt. 3600) then
     write (LDT_logunit, *) '[ERR] Please set the LDT timestep to 1hr or less'
     write (LDT_logunit, *) '[ERR] This is required for SMAP_E_OPL data processing'
     call LDT_endrun()
  endif

  write (yyyymmdd, '(i4.4,2i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da
  write (hh, '(i2.2)') LDT_rc%hr

  list_files = 'ls '//trim(SMAPEOPLsmobs(n)%odir)// &
               '/ARFS_SM_*'//trim(yyyymmdd)//'T'//trim(hh) &
               //"*.dat > SMAP_filelist.dat"

  call system(trim(list_files))

  i = 1
  ftn = LDT_getNextUnitNumber()
  open (ftn, file="./SMAP_filelist.dat", &
        status='old', iostat=ierr)

  do while (ierr .eq. 0)
     read (ftn, '(a)', iostat=ierr) fname
     if (ierr .ne. 0) then
        exit
     endif
     mn_ind = index(fname, trim(yyyymmdd)//'T'//trim(hh))

     mn_ind = index(fname, trim(yyyymmdd)//'T'//trim(hh)) + 11
     read (fname(mn_ind:mn_ind + 1), '(i2.2)') mn
     ss = 0
     call LDT_tick(timenow, doy, gmt, LDT_rc%yr, LDT_rc%mo, LDT_rc%da, &
                   LDT_rc%hr, mn, ss, 0.0)

     smap_filename(i) = fname

     write (LDT_logunit, *) '[INFO] reading ', trim(smap_filename(i))

     call read_SMAPEOPLsm_data(n, smap_filename(i), &
                             SMAPEOPLsmobs(n)%smobs, timenow)

     i = i + 1
  enddo
  call LDT_releaseUnitNumber(ftn)

  call LDT_logSingleDAobs(n, LDT_DAobsData(n)%soilmoist_obs, &
                          SMAPEOPLsmobs(n)%smobs, vlevel=1)

end subroutine readSMAPEOPL_SMObs

!BOP
! 
! !ROUTINE: read_SMAPEOPLsm_data
! \label{read_SMAPEOPLsm_data}
!
! !INTERFACE:
subroutine read_SMAPEOPLsm_data(n, fname, smobs_inp, time)
! 
! !USES:

  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use SMAPEOPLSMobsMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n
  character (len=*)        :: fname
  real                     :: smobs_inp(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real*8                   :: time
!EOP
  integer,  parameter     :: nc=2560, nr=1920
  real*4                  :: sm_raw(SMAPEOPLsmobs(n)%nc,SMAPEOPLsmobs(n)%nr)
  real                    :: sm_in(SMAPEOPLsmobs(n)%nc*SMAPEOPLsmobs(n)%nr)
  real                    :: smobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1               :: sm_data_b(SMAPEOPLsmobs(n)%nc*SMAPEOPLsmobs(n)%nr)
  logical*1               :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer                 :: smid
  integer                 :: ios, nid
  integer                 :: c,r
  integer                 :: ftn1

  ftn1 = LDT_getNextUnitNumber()
  open(unit=ftn1,file=fname,form='unformatted',access='direct',recl=4*nc*nr,status='old')
  read(ftn1, rec=1) sm_raw
  close(1)
  call LDT_releaseUnitNumber(ftn1)

  do r=1,SMAPEOPLsmobs(n)%nr
     do c=1,SMAPEOPLsmobs(n)%nc
        if (sm_raw(c,r)>=0.and.sm_raw(c,r)<=1) then
           sm_in(c+(r-1)*SMAPEOPLsmobs(n)%nc) = sm_raw(c,r)
           sm_data_b(c+(r-1)*SMAPEOPLsmobs(n)%nc) = .true.
        else
           sm_in(c+(r-1)*SMAPEOPLsmobs(n)%nc) = LDT_rc%udef
           sm_data_b(c+(r-1)*SMAPEOPLsmobs(n)%nc) = .false.
        endif
     enddo
  enddo

  if(LDT_isLDTatAfinerResolution(n,0.0937500)) then
     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          sm_data_b, sm_in, smobs_b_ip, smobs_ip, &
          SMAPEOPLsmobs(n)%nc*SMAPEOPLsmobs(n)%nr, &
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          SMAPEOPLsmobs(n)%w11,SMAPEOPLsmobs(n)%w12,&
          SMAPEOPLsmobs(n)%w21,SMAPEOPLsmobs(n)%w22,&
          SMAPEOPLsmobs(n)%n11,SMAPEOPLsmobs(n)%n12,&
          SMAPEOPLsmobs(n)%n21,SMAPEOPLsmobs(n)%n22,LDT_rc%udef,ios)
  else
     call upscaleByAveraging(SMAPEOPLsmobs(n)%nc*SMAPEOPLsmobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_rc%udef, SMAPEOPLsmobs(n)%n11,&
          sm_data_b,sm_in, smobs_b_ip, smobs_ip)
  endif

!overwrite the input data
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(smobs_ip(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then
           smobs_inp(c,r) = &
                smobs_ip(c+(r-1)*LDT_rc%lnc(n))
        endif
     enddo
  enddo

end subroutine read_SMAPEOPLsm_data


