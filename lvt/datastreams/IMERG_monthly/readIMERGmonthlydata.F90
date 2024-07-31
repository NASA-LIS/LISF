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

!NOTE:  Currently only V07A IMERG Final Run monthly data are supported.
subroutine readIMERGmonthlydata(source)

  ! Imports
  use ESMF
  use IMERG_monthly_dataMod, only: imergmonthlydata
  use LVT_coreMod, only: LVT_rc, LVT_isAtAFinerResolution
  use LVT_histDataMod, only: LVT_logSingleDataStreamVar, &
       LVT_MOC_totalprecip
  use LVT_logMod, only: LVT_verify, LVT_logunit
  use LVT_timeMgrMod, only: LVT_calendar

  ! Defaults
  implicit none

  ! Arguments
  integer, intent(in) :: source

  ! Local variables
  character*255       :: filename
  logical             :: file_exists
  real                :: prcp_in(imergmonthlydata(source)%nc, &
       imergmonthlydata(source)%nr)
  real                :: prcp_in1(imergmonthlydata(source)%nc * &
       imergmonthlydata(source)%nr)
  logical*1           :: lb(imergmonthlydata(source)%nc * &
       imergmonthlydata(source)%nr)
  logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                :: prcp(LVT_rc%lnc*LVT_rc%lnr)
  real                :: prcp_final(LVT_rc%lnc, LVT_rc%lnr)
  integer             :: t, c, r
  integer             :: iret, ireaderr
  real                :: currTime
  logical             :: alarmCheck
  integer             :: yr1, mo1, da1, hr1, mn1, ss1
  integer             :: yr2, mo2, da2, hr2, mn2, ss2
  type(ESMF_Time)     :: time1
  type(ESMF_Time)     :: time2
  type(ESMF_TimeInterval) :: lis_ts
  integer :: status
  integer :: days_in_month

  ! External routines
  external :: create_IMERG_monthly_filename
  external :: read_imergmonthly_hdf
  external :: conserv_interp
  external :: upscaleByAveraging

  ! Initialize variables
  prcp = LVT_rc%udef
  prcp_final = LVT_rc%udef
  currTime = float(LVT_rc%dhr(source))*3600 + &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)

  ! Read IMERG at beginning of the month
  alarmCheck = .false.
  if (LVT_rc%dda(source) == 1 .and. &
       LVT_rc%dhr(source) == 0 .and. &
       LVT_rc%dmn(source) == 0 .and. &
       LVT_rc%dss(source) == 0) alarmCheck = .true.

  yr1 = LVT_rc%dyr(source)
  mo1 = LVT_rc%dmo(source)
  da1 = LVT_rc%dda(source)
  hr1 = LVT_rc%dhr(source)
  mn1 = LVT_rc%dmn(source)
  ss1 = 0
  call ESMF_TimeSet(time1, yy=yr1, mm=mo1, dd=da1, &
       h=hr1, m=mn1, s=ss1, calendar=LVT_calendar, rc=status)
  call LVT_verify(status)

  ! Use previous month file.
  call ESMF_TimeIntervalSet(lis_ts, mm=1, calendar=LVT_calendar, &
       rc=status)
  call LVT_verify(status)
  time2 = time1 - lis_ts
  call ESMF_TimeGet(time2, yy=yr2, mm=mo2, dd=da2, &
        h=hr2, m=mn2, s=ss2, calendar=LVT_calendar, &
        rc=status)
  call LVT_verify(status)

  if (alarmCheck) then
     call create_IMERG_monthly_filename(source, yr2, mo2, filename)
     inquire(file=trim(filename), exist=file_exists)

     if (file_exists) then
        write(LVT_logunit,*) '[INFO] Reading IMERG data ',trim(filename)
        call read_imergmonthly_hdf(filename, &
             imergmonthlydata(source)%nc, &
             imergmonthlydata(source)%nr, prcp_in, ireaderr)
        if (ireaderr .eq. 0) then

           ! Use budget-bilinear interpolation if IMERG data are at
           ! coarser resolution than the analysis grid; otherwise, use
           ! upscale averaging.
           prcp_in1 = LVT_rc%udef
           lb = .false.
           t = 1
           do r = 1, imergmonthlydata(source)%nr
              do c = 1,imergmonthlydata(source)%nc
                 prcp_in1(t) = prcp_in(c,r)
                 if (prcp_in1(t) .ge. 0) then
                    lb(t) = .true.
                 end if
                 t = t + 1
              end do ! c
           end do ! r

           if (LVT_isAtAFinerResolution( &
                imergmonthlydata(source)%datares)) then
              call conserv_interp(LVT_rc%gridDesc, lb, prcp_in1, &
                   lo, prcp, &
                   (imergmonthlydata(source)%nc * &
                   imergmonthlydata(source)%nr), &
                   (LVT_rc%lnc*LVT_rc%lnr), &
                   imergmonthlydata(source)%rlat, &
                   imergmonthlydata(source)%rlon, &
                   imergmonthlydata(source)%w112, &
                   imergmonthlydata(source)%w122, &
                   imergmonthlydata(source)%w212, &
                   imergmonthlydata(source)%w222, &
                   imergmonthlydata(source)%n112, &
                   imergmonthlydata(source)%n122, &
                   imergmonthlydata(source)%n212, &
                   imergmonthlydata(source)%n222, &
                   LVT_rc%udef, iret)
           else
              call upscaleByAveraging( &
                   (imergmonthlydata(source)%nc * &
                   imergmonthlydata(source)%nr), &
                   (LVT_rc%lnc*LVT_rc%lnr), LVT_rc%udef, &
                   imergmonthlydata(source)%n11, lb, &
                   prcp_in1, lo, prcp)
           endif
           write(LVT_logunit,*) '[INFO] Finished processing ', &
                trim(filename)
        else
           write(LVT_logunit,*) '[ERR] Read error with IMERG file ', &
                trim(filename)
           prcp = LVT_rc%udef
        endif
     else
        write(LVT_logunit,*) '[ERR] Missing IMERG file ', trim(filename)
        prcp = LVT_rc%udef
     end if ! file_exists

     do r = 1, LVT_rc%lnr
        do c = 1, LVT_rc%lnc
           prcp_final(c,r) = prcp(c+(r-1)*LVT_rc%lnc)
        end do ! c
     end do ! r
  end if ! alarmCheck

  ! Convert mm/hr to kg/m2s.
  do r = 1, LVT_rc%lnr
     do c = 1, LVT_rc%lnc
        if (prcp_final(c,r) .ge. 0) then
           prcp_final(c,r) = prcp_final(c,r) / 3600.
        else
           prcp_final(c,r) = LVT_rc%udef
        end if
     end do ! c
  end do ! r
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip, source, &
       prcp_final, vlevel=1, units='kg/m2s')

  ! Now convert from kg/m2s to kg/m2.  Need to calculate length of
  ! month in seconds.
  lis_ts = time1 - time2
  call ESMF_TimeIntervalGet(lis_ts, d=days_in_month, rc=status)
  call LVT_verify(status)
  do r = 1, LVT_rc%lnr
     do c = 1, LVT_rc%lnc
        if (prcp_final(c,r) >= 0) then
           prcp_final(c,r) = &
                prcp_final(c,r) * days_in_month * 86400 ! kg/m2
        else
           prcp_final(c,r) = LVT_rc%udef
        endif
     enddo ! c
  enddo ! r
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip, source, &
       prcp_final, vlevel=1, units='kg/m2')

end subroutine readIMERGmonthlydata

subroutine read_imergmonthly_hdf(filename, col, row, precipout, ireaderr)

  ! Imports
  use LVT_coreMod, only: LVT_rc
  use LVT_logMod, only: LVT_logunit

#if (defined USE_HDF5)
  use HDF5
#endif

  ! Defaults
  implicit none

  ! Arguments
  character(len=*), intent(in) :: filename
  integer, intent(in)          :: col, row
  integer, intent(out)         :: ireaderr
  real,    intent(out)         :: precipout(col,row)

  ! Local variables
  integer :: xsize, ysize
  character(len=40) :: dsetname = '/Grid/precipitation'
  real :: precipin(row,col)
  integer :: istatus
  integer :: i, j
#if (defined USE_HDF5)
  integer(HSIZE_T), dimension(2) :: dims
  integer(HID_T) :: fileid, dsetid
#endif

  precipin = LVT_rc%udef
  precipout = LVT_rc%udef

#if (defined USE_HDF5)
  xsize = col
  ysize = row
  dims(1) = xsize
  dims(2) = ysize

  ireaderr = 0

  ! Open Fortran interface
  call h5open_f(istatus)
  if (istatus .ne. 0) then
     write(LVT_logunit,*) 'Error opening HDF5 fortran interface'
     ireaderr = istatus
     return
  end if

  ! Open HDF5 file
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, fileid, istatus)
  if (istatus .ne. 0) then
     write(LVT_logunit,*) 'Error opening IMERG file', trim(filename)
     ireaderr = istatus
     call h5close_f(istatus) ! Close HDF5 interface
     return
  end if

  ! Open precip dataset
  call h5dopen_f(fileid, dsetname, dsetid, istatus)
  if (istatus .ne. 0) then
     write(LVT_logunit,*) 'Error opening IMERG dataset', trim(dsetname)
     ireaderr = istatus
     call h5fclose_f(fileid, istatus) ! Close HDF5 file
     call h5close_f(istatus) ! Close HDF5 interface
     return
  end if

  ! Read dataset
  call h5dread_f(dsetid, H5T_NATIVE_REAL, precipin, dims, istatus)
  if (istatus .ne. 0) then
     write(LVT_logunit,*) 'Error reading IMERG dataset', trim(dsetname)
     ireaderr = istatus
     call h5dclose_f(dsetid, istatus) ! Close dataset
     call h5fclose_f(fileid, istatus) ! Close HDF5 file
     call h5close_f(istatus) ! Close HDF5 interface
     return
  end if

  ! Put the real(1:,1:) on the precipout(0:,0:)
  ! precipin is (ysize,xsize) starting at (lon=-179.9,lat=-89.9)
  precipout(1:xsize,1:ysize) = transpose(precipin)

  ! Clean up.  Since the data have already been copied, we will assume
  ! any HDF5 errors from this point will have no impact on the data, and
  ! therefore we will ignore the status codes.
  call h5dclose_f(dsetid, istatus) ! Close HDF5 dataset
  call h5fclose_f(fileid, istatus) ! Close HDF5 file
  call h5close_f(istatus) ! Close HDF5 interface

#endif

end subroutine read_imergmonthly_hdf

subroutine create_IMERG_monthly_filename(source, yr, mo, filename)

  ! Imports
  use IMERG_monthly_dataMod, only: imergmonthlydata
  use LVT_logMod, only: LVT_logunit, LVT_endrun

  ! Defaults
  implicit none

  ! Arguments
  integer, intent(in) :: source, yr, mo
  character(len=*), intent(out) :: filename

  ! Local variables
  character*4   :: cyr
  character*2   :: cmo
  character*100 :: fstem, fext
  character*255 :: odir
  character*4   :: imVer

  write(cyr, '(I4.4)') yr
  write(cmo, '(I2.2)') mo
  odir = trim(imergmonthlydata(source)%odir)
  imVer = imergmonthlydata(source)%imergver

  ! FIXME:  Add support for Early and Late Runs
  if (imergmonthlydata(source)%imergprd == 'final') then
     fstem = '/3B-MO.MS.MRG.3IMERG.'
     fext = '.HDF5'
  else
     write(LVT_logunit,*) "[ERR] Invalid IMERG product option was chosen."
     write(LVT_logunit,*) "[ERR] Please choose either 'final'."
     call LVT_endrun()
  endif

  filename = trim(odir) // "/" // cyr // trim(fstem) // &
       cyr // cmo // "01-S000000-E235959." // cmo // "." // &
       trim(imVer) // fext

end subroutine create_IMERG_monthly_filename
