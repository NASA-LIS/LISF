!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: readLISlsmTEFFobs
! \label{readLISlsmTEFFobs}
!
! !INTERFACE:
subroutine readLISlsmTEFFobs(n)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
#if (defined USE_GRIBAPI)
  use grib_api
#endif
  use LDT_coreMod
  use LDT_DAobsDataMod
  use LDT_historyMod
  use LDT_logMod
  use LISlsmTEFF_obsMod,    only : lsmteffobs
  use LDT_timeMgrMod
!
! !DESCRIPTION:
!  This routine reads the soil temperature fields from a LIS model
!  simulation.
!
!EOP
  implicit none

  integer,   intent(in) :: n

  character*200    :: fname
  logical          :: file_exists
  real             :: teff_data(LDT_rc%lnc(n),LDT_rc%lnr(n))

  integer          :: t, index
  integer          :: ftn
  integer          :: iret
  real             :: topLev,botLev,param_num,trange
  integer          :: igrib,nvars
  real             :: teffvalue1d(lsmteffobs%nc*lsmteffobs%nr)
  real             :: Tsoil01value2d(lsmteffobs%nc, lsmteffobs%nr)
  real             :: Tsoil02value2d(lsmteffobs%nc, lsmteffobs%nr)
  real :: tsoil(lsmteffobs%nc, lsmteffobs%nr, 4)
  integer :: rc1
  integer          :: c,r
  character*20     :: vname
  integer          :: varid
  real             :: kk, cc_6am, cc_6pm    !parameters for calculating effective soil temperature
                                            !from soil layer temperature at 6 AM and 6 PM local time
  real       :: lon
  real       :: gmt
  real       :: lhour
  integer    :: zone
  integer, allocatable :: ntiles_pergrid(:)
  integer, allocatable :: str_tind(:)
  integer :: gid

  interface
     subroutine create_lsm_teff_output_filename(n, form, fname, odir, wstyle, wopt, &
                                         run_dd, map_proj, security_class,   &
                                         distribution_class, data_category,  &
                                         area_of_data, write_interval)
        integer,   intent(IN)        :: n
        character(len=*)             :: fname
        character(len=*)             :: form
        character(len=*)             :: odir
        character(len=*)             :: wstyle
        character(len=*)             :: wopt
        real, dimension(8), optional :: run_dd
        character(len=*),   optional :: map_proj
        character(len=*),   optional :: security_class
        character(len=*),   optional :: distribution_class
        character(len=*),   optional :: data_category
        character(len=*),   optional :: area_of_data
        character(len=*),   optional :: write_interval
     end subroutine create_lsm_teff_output_filename
  end interface


#if (defined USE_GRIBAPI)
  teff_data = LDT_rc%udef

  call create_lsm_teff_output_filename(lsmteffobs%nest,               &
                                  lsmteffobs%format,             &
                                  fname,                       &
                                  lsmteffobs%odir,               &
                                  lsmteffobs%wstyle,             &
                                  lsmteffobs%wopt,               &
                                  lsmteffobs%run_dd,             &
                                  lsmteffobs%map_proj,           &
                                  lsmteffobs%security_class,     &
                                  lsmteffobs%distribution_class, &
                                  lsmteffobs%data_category,      &
                                  lsmteffobs%area_of_data,       &
                                  lsmteffobs%write_interval)

  inquire(file=trim(fname),exist=file_exists)

  if(file_exists) then
     write(LDT_logunit,*) '[INFO] reading LSM output ',trim(fname)
     if(lsmteffobs%format.eq."binary") then
        write(LDT_logunit,*) '[ERR] DA preprocessing on the binary format is not '
        write(LDT_logunit,*) '[ERR] currently supported. Program stopping....'
        call LDT_endrun()

     elseif(lsmteffobs%format.eq."grib1") then
        write(LDT_logunit,*) '[ERR] DA preprocessing on the grib1 format is not '
        write(LDT_logunit,*) '[ERR] currently supported. Program stopping....'
        call LDT_endrun()

     elseif(lsmteffobs%format.eq."netcdf") then
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

        iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
        call LDT_verify(iret, 'Error opening file '//trim(fname))

!  The code looks for instantaneous variables. If it doesn't exist,
!  the time averaged data fields will be read in.

        iret = nf90_inq_varid(ftn, 'SoilTemp_inst', varid)
        vname = 'SoilTemp_inst'
        if(iret.ne.0) then
           vname = 'SoilTemp_tavg'
        endif
        iret = nf90_close(ftn)
        call LDT_verify(iret,'Error in nf90_close')

        if ((LDT_rc%lnc(n) .ne. lsmteffobs%nc) .or. &
             (LDT_rc%lnr(n) .ne. lsmteffobs%nr)) then
           write(LDT_logunit,*)'[ERR] Dimension mismatch for LIS data!'
           write(LDT_logunit,*)'[ERR] LDT_rc%lnc, LDT_rc%lnr = ', &
                LDT_rc%lnc(n), LDT_rc%lnr(n)
           write(LDT_Logunit,*)'[ERR] lsmteffobs%nc, lsmteffobs%nr = ', &
                lsmteffobs%nc, lsmteffobs%nr
           call LDT_endrun()
        end if
        tsoil = 0
        allocate(ntiles_pergrid(lsmteffobs%nc * lsmteffobs%nr))
        ntiles_pergrid = lsmteffobs%ntiles_pergrid ! Copy scalar to array
        allocate(str_tind(lsmteffobs%nc * lsmteffobs%nr))
        do gid = 1, lsmteffobs%nc * lsmteffobs%nr
           str_tind(gid) = ((gid - 1) * lsmteffobs%num_ens) + 1
        end do
        call read_LIStsoil_data_usaf(n, lsmteffobs%num_tiles, str_tind, &
             ntiles_pergrid, &
             lsmteffobs%num_ens, &
             fname, tsoil, rc1)
        if (rc1 .ne. 0) then
           write(LDT_logunit,*) '[ERR] Cannot read from ', trim(fname)
           call LDT_endrun()
        end if
        tsoil01value2d = tsoil(:,:,1)
        tsoil02value2d = tsoil(:,:,2)
        deallocate(ntiles_pergrid)
        deallocate(str_tind)

        kk = 1.007
        cc_6am = 0.246;    !Descending
        cc_6pm = 1.000;    !Ascending

        do r=1,lsmteffobs%nr
           do c=1, lsmteffobs%nc
              !calculate effective soil temperature
              lon = LDT_domain(n)%lon(c+(r-1)*lsmteffobs%nc)
              gmt = LDT_rc%hr
              call LDT_gmt2localtime(gmt, lon, lhour, zone)

              !if(lhour.eq.6) then ! Orig
              if(lhour > 4 .and. lhour < 8) then ! EMK for 3-hrly data

                 if (Tsoil01value2d(c,r).gt.273.15.and.Tsoil02value2d(c,r).gt.273.15) then
                    teffvalue1d(c+(r-1)*lsmteffobs%nc) = &
                          kk * (cc_6am * Tsoil01value2d(c,r) + (1 - cc_6am) * Tsoil02value2d(c,r))
                 else
                    teffvalue1d(c+(r-1)*lsmteffobs%nc) = LDT_rc%udef
                 endif
              !elseif(lhour.eq.18) then ! Orig
              elseif (lhour > 16 .and. lhour < 20) then ! EMK for 3-hrly data

                 if (Tsoil01value2d(c,r).gt.273.15.and.Tsoil02value2d(c,r).gt.273.15) then
                    teffvalue1d(c+(r-1)*lsmteffobs%nc) = &
                          kk * (cc_6pm * Tsoil01value2d(c,r) + (1 - cc_6pm) * Tsoil02value2d(c,r))
                 else
                    teffvalue1d(c+(r-1)*lsmteffobs%nc) = LDT_rc%udef
                 endif
              else
                 teffvalue1d(c+(r-1)*lsmteffobs%nc) = LDT_rc%udef
              endif

           enddo
        enddo

        call transformDataToLDTgrid_teff(n,teffvalue1d,teff_data)

#endif
     endif
  else
     write(LDT_logunit,*) '[WARN] LIS file '//trim(fname)
     write(LDT_logunit,*) '[WARN] not found ...'
     teff_data = LDT_rc%udef
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(1)%teff_obs,&
       teff_data,vlevel=1)

#endif
end subroutine readLISlsmTEFFobs


!BOP
!
! !ROUTINE: create_lsm_teff_output_filename
! \label{create_lsm_teff_output_filename}
!
! !INTERFACE:
subroutine create_lsm_teff_output_filename(n, form, fname, odir, wstyle, wopt, &
                                      run_dd, map_proj, security_class,   &
                                      distribution_class, data_category,  &
                                      area_of_data, write_interval)
! !USES:
   use LDT_coreMod,  only : LDT_rc
   use LDT_logMod

   implicit none
! !ARGUMENTS:
   integer,   intent(IN)        :: n
   character(len=*)             :: fname
   character(len=*)             :: form
   character(len=*)             :: odir
   character(len=*)             :: wstyle
   character(len=*)             :: wopt
   real, dimension(8), optional :: run_dd
   character(len=*),   optional :: map_proj
   character(len=*),   optional :: security_class
   character(len=*),   optional :: distribution_class
   character(len=*),   optional :: data_category
   character(len=*),   optional :: area_of_data
   character(len=*),   optional :: write_interval
!
! !DESCRIPTION:
!  Create the file name for the output data files. It creates both the GSWP
!  style of output filenames and the standard LIS style. The convention used
!  in LIS creates a filename in a hierarchical style (output directory,
!  model name, date, file extention)
!
!  2 level hierarchy
!  \begin{verbatim}
!   <output directory>/<model name>/LIS_HIST_<yyyymmddhhmnss>.<extension>
!  \end{verbatim}
!  3 level hierarchy
!  \begin{verbatim}
!   <output directory>/<model name>/<yyyymm>/LIS_HIST_<yyyymmddhhmnss>.<extension>
!  \end{verbatim}
!  4 level hierarchy
!  \begin{verbatim}
!   <output directory>/<model name>/<yyyy>/<yyyymm>/LIS_HIST_<yyyymmddhhmnss>.<extension>
!  \end{verbatim}
!  WMO convention
!  \begin{verbatim}
!   <output directory>/<AFWA Weather product style>
!  \end{verbatim}
!   A filename in the convention of weather products (such as): \newline
!   {\small
!   PS.AFWA\_SC.U\_DI.C\_DC.ANLYS\_GP.LIS\_GR.C0P25DEG\_AR.GLOBAL\_PA.03-HR-SUM\_DD.YYYYMMDD\_DT.HH00\_DF.GR1 \newline
!   }
!   where                             \newline
!    PS = Product source              \newline
!    SC = security classification     \newline
!    DI = distribution classification \newline
!    DC = data category               \newline
!    GP = generating process          \newline
!    GR = grid                        \newline
!    AR = area of data                \newline
!    PA = parameter                   \newline
!    DD = date                        \newline
!    DT = data time                   \newline
!    DF = data format                 \newline
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest
!   \item [fname]
!     the created file name.
!   \item [model\_name]
!    string describing the name of the model
!   \item [writeint]
!    output writing interval  of the model
!   \item [style]
!    style option as described above
!  \end{description}
!EOP
   character(len=8)        :: date
   character(len=10)       :: time
   character(len=5)        :: zone
   integer, dimension(8)   :: values
   character(len=20)       :: mname
   character(len=10)       :: cdate
   character(len=14)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=200)       :: dname
   character(len=200), save :: out_fname
   integer                  :: i, c

   mname = 'SURFACEMODEL'
   if(wstyle.eq."4 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn

      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mname)//'/'

      write(unit=cdate, fmt='(i4.4)') LDT_rc%yr
      dname = trim(dname)//trim(cdate)//'/'

      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da
      dname = trim(dname)//trim(cdate)

      out_fname = trim(dname)//'/LIS_HIST_'//trim(cdate1)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      out_fname = trim(out_fname)//trim(cdate)

      select case ( form )
      case ( "binary" )
         if(wopt.eq."1d tilespace") then
            out_fname = trim(out_fname)//'.ts4r'
         elseif(wopt.eq."2d gridspace") then
            out_fname = trim(out_fname)//'.gs4r'
         elseif(wopt.eq."1d gridspace") then
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_lsm_teff_output_filename -- '// &
              'Unrecognized output format')
         call LDT_endrun
      endselect
   elseif(wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn

      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mname)//'/'

      write(unit=cdate, fmt='(i4.4, i2.2)') LDT_rc%yr, LDT_rc%mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      out_fname = trim(out_fname)//trim(cdate)

      select case ( form )
      case ("binary")
         if(wopt.eq."1d tilespace") then
            out_fname = trim(out_fname)//'.ts4r'
         elseif(wopt.eq."2d gridspace") then
            out_fname = trim(out_fname)//'.gs4r'
         elseif(wopt.eq."1d gridspace") then
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_lsm_teff_output_filename -- '// &
              'Unrecognized form value')
         call LDT_endrun
      endselect
   elseif(wstyle.eq."2 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn

      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mname)//'/'

      out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      out_fname = trim(out_fname)//trim(cdate)

      select case ( form )
      case ("binary")
         if(wopt.eq."1d tilespace") then
            out_fname = trim(out_fname)//'.ts4r'
         elseif(wopt.eq."2d gridspace") then
            out_fname = trim(out_fname)//'.gs4r'
         elseif(wopt.eq."1d gridspace") then
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_lsm_teff_output_filename -- '// &
              'Unrecognized form value')
         call LDT_endrun
      endselect
   elseif(wstyle.eq."WMO convention") then
      if ( .not. present(run_dd)             .or. &
           .not. present(security_class)     .or. &
           .not. present(distribution_class) .or. &
           .not. present(data_category)      .or. &
           .not. present(area_of_data)       .or. &
           .not. present(write_interval) ) then
         call ldt_log_msg('ERR: create_lsm_teff_output_filename -- '// &
              'missing WMO convention identifiers')
         call LDT_endrun
      endif

      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, LDT_rc%da

      write(unit=cdate, fmt='(i2.2, i2.2)') LDT_rc%hr, LDT_rc%mn

      if(map_proj.eq."polar") then
         fproj = 'P'
         print *,"fres ",run_dd(6)
         if (run_dd(6) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(run_dd(6))
         else
            write(unit=fres, fmt='(i1)') nint(run_dd(6))
         endif
         fres2 = trim(fres)//'KM'
      elseif(map_proj.eq."lambert") then
         fproj = 'L'
         print *,"fres ",run_dd(6)
         write(unit=fres, fmt='(f3.0)') run_dd(6)
         if (run_dd(6) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(run_dd(6))
         else
            write(unit=fres, fmt='(i1)') nint(run_dd(6))
         endif
         fres2 = trim(fres)//'KM'
      elseif(map_proj.eq."mercator") then
         fproj = 'M'
         write(unit=fres, fmt='(i2.2)') run_dd(6)
         fres = trim(fres)//'KM'
      elseif(map_proj.eq."gaussian") then
         fproj = 'G'
         write(unit=fres, fmt='(i2.2)') run_dd(5)*100
         fres2 = '0P'//trim(fres)//'DEG'
      else
         fproj = 'C'
         write(unit=fres, fmt='(i10)') nint(run_dd(6)*100)
         read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
         c = 0
         do i=1,10
            if(fres1(i).ne.' '.and.c==0) c = i
         enddo
         if (run_dd(6) .lt. 0.1) then
            fres3 = '0P0'
         else
            fres3 = '0P'
         end if
         fres2 = fres3
         do i=c,10
            fres2 = trim(fres2)//trim(fres1(i))
         enddo
         fres2 = trim(fres2)//'DEG'
      endif

      out_fname = trim(odir)//'/'//&
           '/PS.AFWA_SC.'//trim(security_class)//&
           '_DI.'//trim(distribution_class)//&
           '_DC.'//trim(data_category)//&
           '_GP.LIS_GR.'//&
           trim(fproj)//trim(fres2)//&
           '_AR.'//trim(area_of_data)//&
           '_PA.'//trim(write_interval)//'-HR-SUM_DD.'//&
           trim(cdate1)//'_DT.'//trim(cdate)//'_DF'
      if (form == "netcdf") then
         out_fname = trim(out_fname) // ".nc"
      else if (form == "grib1") then
         out_fname = trim(out_fname) // ".GR1"
      else if (form == "grib2") then
         out_fname = trim(out_fname) // ".GR2"
      else
         write(LDT_logunit,*)'[ERR] Invalid LIS file format ', trim(form)
         call LDT_endrun()
      end if

   endif
   fname = out_fname
 end subroutine create_lsm_teff_output_filename

!BOP
!
! !ROUTINE: transformDataToLDTgrid_teff
! \label{transformDataToLDTgrid_teff}
!
! !INTERFACE:
 subroutine transformDataToLDTgrid_teff(n, teff_inp, teff_out)
! !USES:
   use LDT_coreMod
   use LISlsmTEFF_obsMod

   implicit none
! !ARGUMENTS:
   integer         :: n
   real            :: teff_inp(lsmteffobs%nc*lsmteffobs%nr)
   real            :: teff_out(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This routine interpolates or upscales the input data to
!  the LDT grid. If the input data is finer than the LDT
!  grid, the input data is upscaled. If the input data is
!  coarser, then it is interpolated to the LDT grid.
!
!EOP
   integer         :: ios
   integer         :: c,r
   logical*1       :: teff_data_b(lsmteffobs%nc*lsmteffobs%nr)
   real            :: teffobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   logical*1       :: teffobs_b_ip(lsmteffobs%nc*lsmteffobs%nr)

   do r=1,lsmteffobs%nr
      do c=1, lsmteffobs%nc
         if(teff_inp(c+(r-1)*lsmteffobs%nc).ne.LDT_rc%udef) then
            teff_data_b(c+(r-1)*lsmteffobs%nc) = .true.
         else
            teff_data_b(c+(r-1)*lsmteffobs%nc) = .false.
         endif
         if(teff_inp(c+(r-1)*lsmteffobs%nc).lt.0) then
            teff_inp(c+(r-1)*lsmteffobs%nc) = LDT_rc%udef
            teff_data_b(c+(r-1)*lsmteffobs%nc) = .false.
         endif
      enddo
   enddo

   if(LDT_isLDTatAfinerResolution(n,lsmteffobs%datares)) then

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!--------------------------------------------------------------------------
      call bilinear_interp(LDT_rc%gridDesc(n,:),&
           teff_data_b, teff_inp, teffobs_b_ip, teffobs_ip, &
           lsmteffobs%nc*lsmteffobs%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           lsmteffobs%w11, lsmteffobs%w12, &
           lsmteffobs%w21, lsmteffobs%w22, &
           lsmteffobs%n11, lsmteffobs%n12, &
           lsmteffobs%n21, lsmteffobs%n22, &
           LDT_rc%udef, ios)

      call neighbor_interp(LDT_rc%gridDesc(n,:),&
           teff_data_b, teff_inp, teffobs_b_ip, teffobs_ip, &
           lsmteffobs%nc*lsmteffobs%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           lsmteffobs%n11, LDT_rc%udef, ios)
   else
      call upscaleByAveraging(&
           lsmteffobs%nc*lsmteffobs%nr,&
           LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%udef, &
           lsmteffobs%n11,teff_data_b, teff_inp, teffobs_b_ip,teffobs_ip)

   endif

   do r=1,LDT_rc%lnr(n)
      do c=1,LDT_rc%lnc(n)
         if(teffobs_b_ip(c+(r-1)*LDT_rc%lnc(n))) then
            teff_out(c,r) = teffobs_ip(c+(r-1)*LDT_rc%lnc(n))
         else
            teff_out(c,r) = LDT_rc%udef
         endif
      enddo
   enddo

 end subroutine transformDataToLDTgrid_teff
