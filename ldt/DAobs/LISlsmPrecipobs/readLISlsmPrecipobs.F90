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
! !ROUTINE: readLISlsmPrecipobs
! \label{readLISlsmPrecipobs}
!
! !INTERFACE: 
! !REVISION HISTORY: 
! 2Dec2021: Mahdi Navari ; Initial Specification (based on readLISlsmSM)
subroutine readLISlsmPrecipobs(n)
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
  use LISlsmPrecip_obsMod,    only : lsmprecipobs
!
! !DESCRIPTION: 
!  This routine reads the total precipitation fields from a LIS model 
!  simulation.  
!
!EOP
  implicit none

  integer,   intent(in) :: n

  character*200    :: fname 
  logical          :: file_exists
  real             :: precip_data(LDT_rc%lnc(n),LDT_rc%lnr(n))
  
  integer          :: t, index
  integer          :: ftn
  integer          :: iret
  real             :: topLev,botLev,param_num,trange
  integer          :: igrib,nvars
  real             :: precipvalue1d(lsmprecipobs%nc*lsmprecipobs%nr)
  real             :: precipvalue2d(lsmprecipobs%nc, lsmprecipobs%nr)
  integer          :: c,r
  character*20     :: vname
  integer          :: varid

  interface
     subroutine create_lsm_output_fname(n, form, fname, odir, wstyle, wopt, &
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
     end subroutine create_lsm_output_fname
  end interface


#if (defined USE_GRIBAPI)
  precip_data = LDT_rc%udef

  call create_lsm_output_fname(lsmprecipobs%nest,               &
                                  lsmprecipobs%format,             &
                                  fname,                       &
                                  lsmprecipobs%odir,               &
                                  lsmprecipobs%wstyle,             &
                                  lsmprecipobs%wopt,               &
                                  lsmprecipobs%run_dd,             &
                                  lsmprecipobs%map_proj,           &
                                  lsmprecipobs%security_class,     &
                                  lsmprecipobs%distribution_class, &
                                  lsmprecipobs%data_category,      &
                                  lsmprecipobs%area_of_data,       &
                                  lsmprecipobs%write_interval)

  inquire(file=trim(fname),exist=file_exists)

  if(file_exists) then 
     write(LDT_logunit,*) '[INFO] reading LSM output ',trim(fname)
     if(lsmprecipobs%format.eq."binary") then 
        write(LDT_logunit,*) '[ERR] DA preprocessing on the binary format is not '
        write(LDT_logunit,*) '[ERR] currently supported. Program stopping....'
        call LDT_endrun()
  
     elseif(lsmprecipobs%format.eq."grib1") then 
        if(lsmprecipobs%wstyle.ne."WMO convention") then 
           write(LDT_logunit,*) '[ERR] LDT currently does not support this style of grib output'
           call LDT_endrun
        endif

        call grib_open_file(ftn,trim(fname),'r',iret)
        if(iret.ne.0) then 
           write(LDT_logunit,*) '[ERR] Could not open file: ',trim(fname)
           call LDT_endrun()
        endif
        call grib_multi_support_on

        do 
           call grib_new_from_file(ftn,igrib,iret)
           if(iret==GRIB_END_OF_FILE) then 
              exit
           endif
           
           call grib_get(igrib, 'indicatorOfParameter',param_num, iret)
           call LDT_verify(iret, &
                'grib_get: indicatorOfParameter failed in readLISlsmPrecipObs')

           call grib_get(igrib, 'topLevel',topLev, iret)
           call LDT_verify(iret, &
                'grib_get: topLevel failed in readLISlsmPrecipObs')

           call grib_get(igrib, 'bottomLevel',botLev, iret)
           call LDT_verify(iret, &
                'grib_get: bottomLevel failed in readLISlsmPrecipObs')

           call grib_get(igrib, 'timeRangeIndicator',trange, iret)
           call LDT_verify(iret, &
                'grib_get: timeRangeIndicator failed in readLISlsmPrecipObs')

!right now specifically geared for AFWA outputs.            
           if(param_num.eq.201.and.topLev.eq.0.and.botlev.eq.10.and.&
                trange.eq.1) then 

              call grib_get(igrib,'values',precipvalue1d,iret)
              call LDT_verify(iret,&
                   'grib_get: values failed in readLISlsmSMobs')
                           
              call transformPrecipDataToLDTgrid(n,precipvalue1d,precip_data)

           endif

           call grib_release(igrib,iret)
           call LDT_verify(iret, 'error in grib_release in readLISlsmSMObs')
        enddo

        call grib_close_file(ftn)

     elseif(lsmprecipobs%format.eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
        
        iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
        call LDT_verify(iret, 'Error opening file '//trim(fname))
        
!  The code looks for instantaneous variables. If it doesn't exist, 
!  the time averaged data fields will be read in. 

        iret = nf90_inq_varid(ftn, 'TotalPrecip_inst', varid)
        vname = 'TotalPrecip_inst'
        if(iret.ne.0) then 
           vname = 'TotalPrecip_tavg'
        endif

        if(lsmprecipobs%datares.eq.LDT_rc%gridDesc(n,10)) then 
           call LDT_readLISSingleNetcdfVar(n,ftn, vname,&
                1,lsmprecipobs%nc, lsmprecipobs%nr, precipvalue2d)
        else
           call LDT_readLISSingleNetcdfVar(n,ftn, vname,&
                1,lsmprecipobs%nc, lsmprecipobs%nr, precipvalue2d)
        endif
        
        iret = nf90_close(ftn)
        call LDT_verify(iret,'Error in nf90_close')
        
        do r=1,lsmprecipobs%nr
           do c=1, lsmprecipobs%nc
              precipvalue1d(c+(r-1)*lsmprecipobs%nc) = precipvalue2d(c,r)
           enddo
        enddo
        
        call transformPrecipDataToLDTgrid(n,precipvalue1d,precip_data)

#endif
     endif
  else
     write(LDT_logunit,*) '[WARN] LIS file '//trim(fname)
     write(LDT_logunit,*) '[WARN] not found ...'
     precip_data = LDT_rc%udef
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(1)%totalprecip_obs,&
       precip_data,vlevel=1)

#endif
end subroutine readLISlsmPrecipobs


!BOP
!
! !ROUTINE: create_lsm_output_fname
! \label{create_lsm_output_fname}
!
! !INTERFACE:
subroutine create_lsm_output_fname(n, form, fname, odir, wstyle, wopt, &
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
         call ldt_log_msg('ERR: create_lsm_output_fname -- '// &
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
         call ldt_log_msg('ERR: create_lsm_output_fname -- '// &
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
         call ldt_log_msg('ERR: create_lsm_output_fname -- '// &
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
         call ldt_log_msg('ERR: create_lsm_output_fname -- '// &
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
!         write(unit=fres, fmt='(f2.0)') run_dd(6)
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
           trim(cdate1)//'_DT.'//trim(cdate)//'_DF.GR1'
   endif
   fname = out_fname
 end subroutine create_lsm_output_fname

!BOP
! 
! !ROUTINE: transformPrecipDataToLDTgrid
! \label{trasnformDataToLDTgrid}
!
! !INTERFACE: 
 subroutine transformPrecipDataToLDTgrid(n, precip_inp, precip_out)
! !USES:    
   use LDT_coreMod
   use LISlsmPrecip_obsMod

   implicit none
! !ARGUMENTS: 
   integer         :: n 
   real            :: precip_inp(lsmprecipobs%nc*lsmprecipobs%nr)
   real            :: precip_out(LDT_rc%lnc(n),LDT_rc%lnr(n))
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
   logical*1       :: precip_data_b(lsmprecipobs%nc*lsmprecipobs%nr)
   real            :: precipobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   logical*1       :: precipobs_b_ip(lsmprecipobs%nc*lsmprecipobs%nr)

   do r=1,lsmprecipobs%nr
      do c=1, lsmprecipobs%nc
         if(precip_inp(c+(r-1)*lsmprecipobs%nc).ne.LDT_rc%udef) then 
            precip_data_b(c+(r-1)*lsmprecipobs%nc) = .true. 
         else
            precip_data_b(c+(r-1)*lsmprecipobs%nc) = .false.
         endif
         if(precip_inp(c+(r-1)*lsmprecipobs%nc).gt.1) then 
            precip_inp(c+(r-1)*lsmprecipobs%nc) = LDT_rc%udef
            precip_data_b(c+(r-1)*lsmprecipobs%nc) = .false.
         endif
      enddo
   enddo

   if(LDT_isLDTatAfinerResolution(n,lsmprecipobs%datares)) then 

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
      call bilinear_interp(LDT_rc%gridDesc(n,:),&
           precip_data_b, precip_inp, precipobs_b_ip, precipobs_ip, &
           lsmprecipobs%nc*lsmprecipobs%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           lsmprecipobs%w11, lsmprecipobs%w12, &
           lsmprecipobs%w21, lsmprecipobs%w22, &
           lsmprecipobs%n11, lsmprecipobs%n12, &
           lsmprecipobs%n21, lsmprecipobs%n22, &
           LDT_rc%udef, ios)

      call neighbor_interp(LDT_rc%gridDesc(n,:),&
           precip_data_b, precip_inp, precipobs_b_ip, precipobs_ip, &
           lsmprecipobs%nc*lsmprecipobs%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           lsmprecipobs%n11, LDT_rc%udef, ios)
   else
      call upscaleByAveraging(&
           lsmprecipobs%nc*lsmprecipobs%nr,&
           LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%udef, &
           lsmprecipobs%n11,precip_data_b, precip_inp, precipobs_b_ip,precipobs_ip)
      
   endif
   
   do r=1,LDT_rc%lnr(n)
      do c=1,LDT_rc%lnc(n)
         if(precipobs_b_ip(c+(r-1)*LDT_rc%lnc(n))) then 
            precip_out(c,r) = precipobs_ip(c+(r-1)*LDT_rc%lnc(n))
         else
            precip_out(c,r) = LDT_rc%udef
         endif
      enddo
   enddo
         
 end subroutine transformPrecipDataToLDTgrid
