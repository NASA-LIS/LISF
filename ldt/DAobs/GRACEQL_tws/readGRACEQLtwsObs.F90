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
! !ROUTINE: readGRACEQLtwsObs
! \label{readGRACEQLtwsObs}
! 
! !REVISION HISTORY: 
!  14 Feb 2018: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readGRACEQLtwsObs(n)
! !USES:   

  use LDT_coreMod
  use LDT_logMod
  use LDT_historyMod
  use LDT_DAobsDataMod
  use LDT_timeMgrMod
  use LDT_constantsMod
  use GRACEQLtws_obsMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine reads the GRACE QL data, computes the anomalies and 
! generates a new set of GRACE observations by incorporating the anomalies
! into the model generated terrestrial water storage observations. 
!
! TWS outputs from LIS is expected to be in units of mm.
!
!EOP
  character(len=LDT_CONST_PATH_LEN) :: fname,filename,gracefile
  integer               :: c,r,c1,r1,k,t,iret
  integer               :: ftn
  integer               :: yr,mo,da,hr
  integer               :: timeId, tId, tbId, lweId, scaleId,errId,err1Id
  real                  :: tws_data_ip(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real                  :: tws_data(GRACEQLtwsobs%nc,GRACEQLtwsobs%nr)
  real                  :: output_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real,    allocatable  :: output_basin_data(:)
  integer, allocatable  :: output_basin_data_count(:)
  real,    allocatable  :: var_merr(:)
  real,    allocatable  :: var_lerr(:)
  real,    allocatable  :: basin_merr(:)
  real,    allocatable  :: basin_lerr(:)
  integer, allocatable  :: nerr(:)
  real                  :: expdbm,expdbl
  real                  :: lat1,lon1,lat2,lon2,dist,betam, betal
  real                  :: output_err_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                  :: output_merr_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                  :: output_lerr_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1             :: output_bitmap(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1             :: output_err_bitmap(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                  :: input_data(GRACEQLtwsobs%gracenc*GRACEQLtwsobs%gracenr)
  logical*1             :: input_bitmap(GRACEQLtwsobs%gracenc*GRACEQLtwsobs%gracenr)
  integer               :: npts(GRACEQLtwsobs%gracenc,GRACEQLtwsobs%gracenr) 
  logical               :: file_exists
  real                  :: dt
  integer               :: cat_val, twsid
  integer               :: yyyy,mm,dd,hh
  character*10          :: ftime
  character*4           :: fyr
  logical               :: valid_data
  integer               :: md_nc
  integer               :: status
  integer               :: dayoffset

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
!
! At initial time, read the data into memory. At all other times, 
! simply index into the correct temporal location of the data. 
!
  if( GRACEQLtwsobs%startMode .or. GRACEQLtwsobs%yr.ne.LDT_rc%yr) then 
     GRACEQLtwsobs%startMode = .false. 
     GRACEQLtwsobs%yr = LDT_rc%yr

     ! Read 'raw' GRACE obs file:
     write(unit=fyr,fmt='(i4.4)') LDT_rc%yr
     gracefile = trim(GRACEQLtwsobs%gracefile)//'/RSWM_'//&
          trim(fyr)//'.nc'
     inquire(file=trim(gracefile),exist=file_exists) 
     if(file_exists) then 
        write(LDT_logunit,*) '[INFO] Reading GRACE data '//&
             trim(gracefile)
        call LDT_verify(nf90_open(path=trim(gracefile),&
             mode=nf90_nowrite, ncid = ftn), &
             'nf90_open failed in readGRACEQLtwsObs ')
        call LDT_verify(nf90_inq_dimid(ftn,'time',timeId),&
             'nf90_inq_dimid failed in readGRACEQLtwsObs')
        call LDT_verify(nf90_inquire_dimension(ftn,timeId,&
             len=GRACEQLtwsobs%tdims),&
             'nf90_inq_dimension failed in readGRACEQLtwsObs')

        if(allocated(GRACEQLtwsobs%tvals)) then 
           deallocate(GRACEQLtwsobs%tvals)
        endif
        if(allocated(GRACEQLtwsobs%lwe_thickness)) then 
           deallocate(GRACEQLtwsobs%lwe_thickness)
        endif
        if(allocated(GRACEQLtwsobs%twsavg)) then 
           deallocate(GRACEQLtwsobs%twsavg)
        endif

        allocate(GRACEQLtwsobs%tvals(GRACEQLtwsobs%tdims))
        allocate(GRACEQLtwsobs%lwe_thickness(GRACEQLtwsobs%gracenc, &
             GRACEQLtwsobs%gracenr,GRACEQLtwsobs%tdims))
        allocate(GRACEQLtwsobs%twsavg(GRACEQLtwsobs%gracenc, &
             GRACEQLtwsobs%gracenr))

        call LDT_verify(nf90_inq_varid(ftn,'TWS',lweId),&
             'nf90_inq_varid failed for TWS in readGRACEQLtwsObs')
        call LDT_verify(nf90_get_var(ftn,lweId,GRACEQLtwsobs%lwe_thickness),&
             'nf90_get_var failed for TWS in readGRACEQLtwsObs')
        
        call LDT_verify(nf90_inq_varid(ftn,'time',tId),&
             'nf90_inq_varid failed for time in readGRACEQLtwsObs')
        call LDT_verify(nf90_get_var(ftn,tId,GRACEQLtwsobs%tvals),&
             'nf90_get_var failed for time in readGRACEQLtwsObs')
        
        call LDT_verify(nf90_close(ftn))
                       
        npts = 0 
        GRACEQLtwsobs%twsavg = 0 
        
!        call getGRACEQLtimeoffset(LDT_rc%yr,LDT_rc%mo,LDT_rc%da,dayoffset)
!        print*, dayoffset
!        stop
        
        ! Loop over GRACE time dimensions:
        do k=1,GRACEQLtwsobs%tdims


           ! Within baseline GRACE data period for calculating TWS climatology:
           if(LDT_rc%yr.ge.GRACEQLtwsobs%b_syr.and.&
                LDT_rc%yr.le.GRACEQLtwsobs%b_eyr) then

              do r=1,GRACEQLtwsobs%gracenr
                 do c=1,GRACEQLtwsobs%gracenc
                    if(GRACEQLtwsobs%lwe_thickness(c,r,k).ne.32767.0) then 
                       GRACEQLtwsobs%twsavg(c,r) = GRACEQLtwsobs%twsavg(c,r) + & 
                            GRACEQLtwsobs%lwe_thickness(c,r,k)
                       npts(c,r) = npts(c,r) + 1
                    endif
                 enddo
              enddo              
           endif   ! Conditional for baseline averaging period
        enddo

        ! Calculate average and convert to mm.
        do r=1,GRACEQLtwsobs%gracenr
           do c=1,GRACEQLtwsobs%gracenc
              if( npts(c,r).gt.0 ) then   ! Convert to mm. 
                 GRACEQLtwsobs%twsavg(c,r) = GRACEQLtwsobs%twsavg(c,r)&
                                         /npts(c,r)
              else
                 GRACEQLtwsobs%twsavg(c,r) = LDT_rc%udef
              endif
           enddo
        enddo
        write(LDT_logunit,*) '[INFO] Finished reading GRACE data '//&
             trim(gracefile)


     else
        write(LDT_logunit,*) '[ERR] GRACE raw obs file '&
             //trim(gracefile)//& 
             'does not exist ...'
        call LDT_endrun()
     endif
  endif
  
  ! Read during the first pass for averaging:
  if(LDT_rc%pass_id.eq.1) then 
     
     call create_lsm_QL_twsoutput_filename(GRACEQLtwsobs%nest, &
          GRACEQLtwsobs%format,&
          fname,GRACEQLtwsobs%odir, GRACEQLtwsobs%wstyle, &
          GRACEQLtwsobs%wopt,'SURFACEMODEL')
     
     ! Average only between baseline years:
     if(LDT_rc%yr.ge.GRACEQLtwsobs%b_syr.and.&
          LDT_rc%yr.le.GRACEQLtwsobs%b_eyr) then

        inquire(file=trim(fname),exist=file_exists)
        if(file_exists) then 
           write(LDT_logunit,*) '[INFO] reading LSM output ',trim(fname)
        
           if(GRACEQLtwsobs%format.eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
              
              iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
              call LDT_verify(iret, 'Error opening file '//trim(fname))
              
              iret = nf90_inq_varid(ftn, "TWS_tavg", twsid)
              call LDT_verify(iret, 'Error in nf90_inq_varid: TWS_tavg')
              
              iret = nf90_get_var(ftn, twsid, tws_data,start = (/1,1/),&
                   count = (/GRACEQLtwsobs%nc,GRACEQLtwsobs%nr/))
              call LDT_verify(iret,'Error in nf90_get_var: TWS_tavg')

              iret = nf90_close(ftn)
              call LDT_verify(iret,'Error in nf90_close')
#endif
              call transformLISoutToGRACEQLgrid(n,tws_data,tws_data_ip)

              do r=1,LDT_rc%lnr(n)
                 do c=1,LDT_rc%lnc(n)
                    if(tws_data_ip(c,r).ne.LDT_rc%udef) then 
                       GRACEQLtwsobs%lisavg(c,r) = &
                            GRACEQLtwsobs%lisavg(c,r) + tws_data_ip(c,r)
                       GRACEQLtwsobs%nlisavg(c,r) = &
                            GRACEQLtwsobs%nlisavg(c,r) + 1
                    endif
                 enddo
              enddo

           endif
        else
           write(LDT_logunit,*) '[ERR] LIS file '//trim(fname)
           write(LDT_logunit,*) '[ERR] not found. Program stopping ...'
           call LDT_endrun()
        endif
     endif
     ! At the end of the first cycle. 

     if(LDT_rc%pass_id.eq.1.and.LDT_rc%endtime.eq.1) then 
        do r=1,LDT_rc%lnr(n)
           do c=1,LDT_rc%lnc(n)
              if(GRACEQLtwsobs%nlisavg(c,r).gt.0) then 
                 GRACEQLtwsobs%lisavg(c,r) = &
                      GRACEQLtwsobs%lisavg(c,r)/&
                      GRACEQLtwsobs%nlisavg(c,r) 
              else
                 GRACEQLtwsobs%lisavg(c,r) = LDT_rc%udef
              endif
           enddo
        enddo
     endif

  ! Second pass through files:
  elseif(LDT_rc%pass_id.eq.2) then 
!     call LDT_get_julhr(LDT_rc%yr,LDT_rc%mo,LDT_rc%da,&
!          LDT_rc%hr,LDT_rc%mn,LDT_rc%ss,currTime)
!     
!     dt = float((currTime-GRACEQLtwsobs%refTime))/24.0
!
!    if(GRACEQLtwsobs%datasource.eq."GRACE TWS Mascon 0.5 deg") then 
!       md_nc=360
!    elseif(GRACEQLtwsobs%datasource.eq."GRACE TWS Original 1 deg") then 
!       md_nc=180
!    endif
     md_nc = 360

     call getGRACEQLtimeoffset(LDT_rc%yr,LDT_rc%mo,LDT_rc%da,&
          GRACEQLtwsobs%tdims, GRACEQLtwsobs%tvals, k)
     
     input_data = LDT_rc%udef
           
     if(k.ne.-1) then 
        ! NO scaling performed:
        do r=1,GRACEQLtwsobs%gracenr
           do c=1,GRACEQLtwsobs%gracenc               
              if(GRACEQLtwsobs%twsavg(c,r).ne.LDT_rc%udef.or.&
                   GRACEQLtwsobs%lwe_thickness(c,r,k).ne.32767) then 
                 ! Convert GRACE 0 to 360dg domain to -180 to 180E:
                 if(c.le.md_nc) then 
                    input_data(c+(r-1)*GRACEQLtwsobs%gracenc+md_nc) = & 
                         (GRACEQLtwsobs%lwe_thickness(c,r,k)-&
                         GRACEQLtwsobs%twsavg(c,r))
                 else
                    input_data(c+(r-1)*GRACEQLtwsobs%gracenc-md_nc) = & 
                         (GRACEQLtwsobs%lwe_thickness(c,r,k)-&
                         GRACEQLtwsobs%twsavg(c,r))
                 endif
              endif
           enddo
        enddo
        
        input_bitmap = .false. 
        do t=1,GRACEQLtwsobs%gracenc*GRACEQLtwsobs%gracenr
           if(input_data(t).ne.LDT_rc%udef) then 
              input_bitmap(t) = .true.
           endif
        enddo
        
        call neighbor_interp(LDT_rc%gridDesc(n,:),&
             input_bitmap, input_data, &
             output_bitmap, output_data, & 
             GRACEQLtwsobs%gracenc*GRACEQLtwsobs%gracenr, & 
             LDT_rc%lnc(n)*LDT_rc%lnr(n),&
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             GRACEQLtwsobs%n111,&
             LDT_rc%udef, iret)
        
        
        output_err_data = LDT_rc%udef
        !
        ! Generate processed observations by adding the anomalies to the 
        ! LIS output data. 
        ! 
        if(GRACEQLtwsobs%process_basin_scale.eq.1) then 
           
           allocate(output_basin_data(GRACEQLtwsobs%basin_cat_max))
           allocate(output_basin_data_count(GRACEQLtwsobs%basin_cat_max))
           
           output_basin_data = 0.0
           output_basin_data_count = 0
           
           do r=1,LDT_rc%lnr(n)
              do c=1,LDT_rc%lnc(n)
                 if(GRACEQLtwsobs%basin_cat(c,r).ne.LDT_rc%udef) then 
                    cat_val = nint(GRACEQLtwsobs%basin_cat(c,r))
                    output_basin_data(cat_val) = &
                         output_basin_data(cat_val) + & 
                         output_data(c+(r-1)*LDT_rc%lnc(n))
                    output_basin_data_count(cat_val) = &
                         output_basin_data_count(cat_val) + 1 
                 endif
              enddo
           enddo
           
           do t=1,GRACEQLtwsobs%basin_cat_max
              output_basin_data(t) = &
                   output_basin_data(t)/&
                   float(output_basin_data_count(t))
           enddo
           
           do r=1,LDT_rc%lnr(n)
              do c=1,LDT_rc%lnc(n)
                 if(output_data(c+(r-1)*LDT_rc%lnc(n)).ne.LDT_rc%udef) then
                    if(GRACEQLtwsobs%basin_cat(c,r).ne.LDT_rc%udef) then
                       cat_val = nint(GRACEQLtwsobs%basin_cat(c,r))
                       output_data(c+(r-1)*LDT_rc%lnc(n)) = &
                            output_basin_data(cat_val)
                    endif
                 endif
              enddo
           enddo
           
           deallocate(output_basin_data)
           deallocate(output_basin_data_count)
        endif
        
        do r=1,LDT_rc%lnr(n)
           do c=1,LDT_rc%lnc(n)
              if(GRACEQLtwsobs%lisavg(c,r).gt.0.0 .and. & !ag
                   output_data(c+(r-1)*LDT_rc%lnc(n)).ne.LDT_rc%udef) then 
                 output_data(c+(r-1)*LDT_rc%lnc(n)) = &
                      GRACEQLtwsobs%lisavg(c,r) + & !ag
                      output_data(c+(r-1)*LDT_rc%lnc(n))*10.0 !to mm. 
              else
                 output_data(c+(r-1)*LDT_rc%lnc(n)) = &
                      LDT_rc%udef
              endif
           enddo
        enddo
        
        ! Write processed observation at the GRACE observation timestamp
        valid_data = .false. 
        do r=1,LDT_rc%lnr(n)
           do c=1,LDT_rc%lnc(n)
              if(output_data(c+(r-1)*LDT_rc%lnc(n)).ne.LDT_rc%udef) then
                 valid_data = .true.
              endif
           enddo
        enddo
        
        ! Write out processed GRACE output data:
        if(valid_data) then
           ftn = LDT_getNextUnitNumber()
           !            currTime = GRACEQLtwsobs%tvals(k)*24+GRACEQLtwsobs%reftime
           !            call LDT_julhr_date(currTime, yyyy,mm,dd,hh)
           
           write(unit=ftime,fmt='(i4.4,i2.2,i2.2)') &
                LDT_rc%yr, LDT_rc%mo, LDT_rc%da
           filename = trim(LDT_rc%odir)//'/GRACE_obs_'//trim(ftime)//'.bin'
           
           write(LDT_logunit,*) "[INFO] writing processed GRACE obs "//&
                trim(filename)
           open(ftn,file=trim(filename), form='unformatted')
           write(ftn) output_data
           write(ftn) output_err_data
           call LDT_releaseUnitNumber(ftn)
           
        endif
        
     endif
  endif
#endif

end subroutine readGRACEQLtwsObs

!BOP
!
! !ROUTINE: create_lsm_QL_twsoutput_filename
! \label{create_lsm_QL_twsoutput_filename}
!
! !INTERFACE:
subroutine create_lsm_QL_twsoutput_filename(n, form, fname, odir, wstyle, wopt,mname)
! !USES:
   use LDT_coreMod,  only : LDT_rc
   use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

   implicit none 
! !ARGUMENTS:
   integer,   intent(IN)        :: n
   character(len=*)             :: fname
   character(len=*)             :: form
   character(len=*)             :: odir
   character(len=*)             :: wstyle
   character(len=*)             :: wopt
   character(len=*)             :: mname !ag (21Dec2017)
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
   !character(len=20)       :: mname !ag (21Dec2017)
   character(len=10)       :: cdate
   character(len=14)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=LDT_CONST_PATH_LEN)       :: dname
   character(len=LDT_CONST_PATH_LEN), save :: out_fname
   integer                  :: i, c

   !ag (21Dec2017)
   !mname = 'SURFACEMODEL' 

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
         call ldt_log_msg('ERR: create_lsm_twsoutput_filename -- '// &
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
         call ldt_log_msg('ERR: create_lsm_twsoutput_filename -- '// &
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
         call ldt_log_msg('ERR: create_lsm_twsoutput_filename -- '// &
              'Unrecognized form value')
         call LDT_endrun 
      endselect

   elseif(wstyle.eq."WMO convention") then 
      write(LDT_logunit,*) '[WARN] WMO convention style not currently supported '
   endif

   fname = out_fname

 end subroutine create_lsm_QL_twsoutput_filename


 subroutine transformLISoutToGRACEQLgrid(n,tws_inp,tws_out)

  use LDT_coreMod
  use GRACEQLtws_obsMod

  integer,    intent(in) :: n
  real                   :: tws_inp(GRACEQLtwsobs%nc*GRACEQLtwsobs%nr)
  real                   :: tws_out(LDT_rc%lnc(n),LDT_rc%lnr(n))
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
  logical*1       :: tws_data_b(GRACEQLtwsobs%nc*GRACEQLtwsobs%nr)
  real            :: twsobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1       :: twsobs_b_ip(GRACEQLtwsobs%nc*GRACEQLtwsobs%nr)
  
   do r=1,GRACEQLtwsobs%nr
      do c=1, GRACEQLtwsobs%nc
         if(tws_inp(c+(r-1)*GRACEQLtwsobs%nc).ne.LDT_rc%udef) then 
            tws_data_b(c+(r-1)*GRACEQLtwsobs%nc) = .true. 
         else
            tws_data_b(c+(r-1)*GRACEQLtwsobs%nc) = .false.
         endif
      enddo
   enddo

   if(LDT_isLDTatAfinerResolution(n,GRACEQLtwsobs%datares)) then 

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
      call bilinear_interp(LDT_rc%gridDesc(n,:),&
           tws_data_b, tws_inp, twsobs_b_ip, twsobs_ip, &
           GRACEQLtwsobs%nc*GRACEQLtwsobs%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           GRACEQLtwsobs%w11, GRACEQLtwsobs%w12, &
           GRACEQLtwsobs%w21, GRACEQLtwsobs%w22, &
           GRACEQLtwsobs%n11, GRACEQLtwsobs%n12, &
           GRACEQLtwsobs%n21, GRACEQLtwsobs%n22, &
           LDT_rc%udef, ios)
   else
      call upscaleByAveraging(&
           GRACEQLtwsobs%nc*GRACEQLtwsobs%nr,&
           LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%udef, &
           GRACEQLtwsobs%n11,tws_data_b, tws_inp, twsobs_b_ip,twsobs_ip)
      
   endif
   
   do r=1,LDT_rc%lnr(n)
      do c=1,LDT_rc%lnc(n)
         if(twsobs_b_ip(c+(r-1)*LDT_rc%lnc(n))) then 
            tws_out(c,r) = twsobs_ip(c+(r-1)*LDT_rc%lnc(n))
         else
            tws_out(c,r) = LDT_rc%udef
         endif
      enddo
   enddo
  
 end subroutine transformLISoutToGRACEQLgrid


!BOP
! !ROUTINE: getGRACEQLtimeoffset
! \label{getGRACEQLtimeoffset}
! 
! !INTERFACE:
  subroutine getGRACEQLtimeoffset(yr,mo,da,tdims, tvals,dayoffset)

    implicit none
   
! !ARGUMENTS: 
    integer, intent(in)     :: yr
    integer, intent(in)     :: mo
    integer, intent(in)     :: da
    integer                 :: tdims
    real                    :: tvals(tdims)
    integer                 :: dayoffset
!
! !DESCRIPTION: 
!
! Returns the julian hour. In this convention, julian day began at
! midnight at the beginning of May 24, 1968. Interestingly, this
! convention was introduced by NASA for the space program. 
! 
! The arguments are: 
!  \begin{description}
!   \item [yr]
!     year
!   \item [mo]
!     month
!   \item [da]
!     day of the month
!   \item [hr]
!     hour of day
!   \item [mn]
!     minute
!   \item[ss]
!     second 
!   \item [dayoffset]
!     julian hour
!  \end{description}
!EOP
    integer      :: years
    integer      :: flagly
    integer      :: days
    integer      :: lyears
    integer      :: month(12) 
    integer      :: temphr
    integer      :: thours
    integer      :: hours
    integer      :: k

    integer                 :: jultime

    data month /31,28,31,30,31,30,31,31,30,31,30,31/

    jultime = 0 

    do k=1859,yr-1
       if((mod(k,4) .eq. 0 .and. mod(k, 100).ne.0) &!leap year
            .or.(mod(k,400) .eq.0)) then 
          jultime = jultime + 366
       else 
          jultime = jultime + 365
       endif
    enddo

    !add the offset from 1858/11/17
    jultime = jultime + 13 + 31

    if((mod(yr,4) .eq. 0 .and. mod(yr, 100).ne.0) &!leap year
         .or.(mod(yr,400) .eq.0)) then 
       month(2) = 29
    else
       month(2) = 28
    endif
    
    do k=1,mo-1
       jultime = jultime + month(mo)
    enddo

    jultime = jultime + da

    dayoffset = -1
    do k=1,tdims
       if(tvals(k).eq.jultime) then 
          dayoffset = k
          exit;
       endif
    enddo

  end subroutine getGRACEQLtimeoffset

