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
! !ROUTINE: readLISoutNatureRun
! \label{readLISoutNatureRun}
!
! !INTERFACE: 
subroutine readLISoutNatureRun(n)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
#if (defined USE_GRIBAPI)
  use grib_api
#endif
  use LDT_coreMod
  use LDT_obsSimMod
  use LDT_historyMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LISoutNatureRun_Mod
!
! !DESCRIPTION: 
!  This routine reads the soil moisture fields from a LIS model 
!  simulation.  
!
!EOP
  implicit none

  integer,   intent(in) :: n
  character(len=LDT_CONST_PATH_LEN)         :: fname
  integer               :: c,r,k
  logical               :: file_exists
  integer               :: ftn
  integer               :: iret
  integer               :: varid
  real                  :: value2d(LISoutNatureRunData%nc,&
       LISoutNatureRunData%nr,LDT_obsSim_struc%nVars)
  real                  :: value1d(LISoutNatureRunData%nc*&
       LISoutNatureRunData%nr,LDT_obsSim_struc%nVars)
  real                  :: var_data(LDT_rc%lnc(n),LDT_rc%lnr(n),&
       LDT_obsSim_struc%nVars)

! create LIS filename
  call create_lisout_naturerun_filename(n, &
       LISoutNatureRunData%mClass,&
       LISoutNatureRunData%format,&
       fname, &
       LISoutNatureRunData%odir,&
       LISoutNatureRunData%wstyle, &
       LISoutNatureRunData%wopt)

! read, subset and interpolate the data
  
  inquire(file=trim(fname),exist=file_exists)

  if(file_exists) then 
     write(LDT_logunit,*) '[INFO] reading LIS output ',trim(fname)
     if(LISoutNatureRunData%format.eq."binary") then 
        write(LDT_logunit,*) '[ERR] observation simulator in binary format is not '
        write(LDT_logunit,*) '[ERR] currently supported. Program stopping....'
        call LDT_endrun()
  
     elseif(LISoutNatureRunData%format.eq."grib1") then 
        write(LDT_logunit,*) '[ERR] observation simulator in grib1 format is not '
        write(LDT_logunit,*) '[ERR] currently supported. Program stopping....'
        call LDT_endrun()

     elseif(LISoutNatureRunData%format.eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
        
        iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
        call LDT_verify(iret, 'Error opening file '//trim(fname))
        
        do k=1,LDT_obsSim_struc%nVars
           iret = nf90_inq_varid(ftn, LDT_obsSim_struc%varNames(k), varid)

           if(LISoutNatureRunData%datares.eq.LDT_rc%gridDesc(n,10)) then 
              call LDT_readLISSingleNetcdfVar(n,ftn, &
                   LDT_obsSim_struc%varNames(k),&
                   1,LISoutNatureRunData%nc, LISoutNatureRunData%nr, &
                   value2d(:,:,k))
           else
              call LDT_readLISSingleNetcdfVar(n,ftn, &
                   LDT_obsSim_struc%varNames(k),&
                   1,LISoutNatureRunData%nc, LISoutNatureRunData%nr, &
                   value2d(:,:,k))
           endif
        enddo

        iret = nf90_close(ftn)
        call LDT_verify(iret,'Error in nf90_close')
 
        do k=1,LDT_obsSim_struc%nVars
           do r=1,LISoutNatureRunData%nr
              do c=1, LISoutNatureRunData%nc
                 value1d(c+(r-1)*LISoutNatureRunData%nc,k) = value2d(c,r,k)
              enddo
           enddo
           
           call convertNatureRunToLDTgrid(n,value1d(:,k),var_data(:,:,k))
        enddo

#endif
     endif
  else
     write(LDT_logunit,*) '[WARN] LIS file '//trim(fname)
     write(LDT_logunit,*) '[WARN] not found ...'
     var_data = LDT_rc%udef
  endif

! store the data
  do k=1,LDT_obsSim_struc%nVars
     call LDT_logNatureRunData(n, var_data(:,:,k),k)
  enddo

end subroutine readLISoutNatureRun

!BOP
!
! !ROUTINE: create_lisout_naturerun_filename
! \label{create_lisout_naturerun_filename}
!
! !INTERFACE:
subroutine create_lisout_naturerun_filename(n, &
     mclass,form, fname, odir, wstyle, wopt)
! !USES:
   use LDT_coreMod,  only : LDT_rc
   use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

   implicit none 
! !ARGUMENTS:
   integer,   intent(IN)        :: n
   character(len=*)             :: mClass
   character(len=*)             :: fname
   character(len=*)             :: form
   character(len=*)             :: odir
   character(len=*)             :: wstyle
   character(len=*)             :: wopt
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
   character(len=10)       :: cdate
   character(len=14)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=LDT_CONST_PATH_LEN) :: dname
   character(len=LDT_CONST_PATH_LEN), save :: out_fname
   integer                  :: i, c

   if(wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mClass)//'/'
      
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
         call ldt_log_msg('ERR: create_lisout_naturerun_filename -- '// &
              'Unrecognized output format')
         call LDT_endrun 
      endselect
   elseif(wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mClass)//'/'
      
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
         call ldt_log_msg('ERR: create_lisout_naturerun_filename -- '// &
              'Unrecognized form value')
         call LDT_endrun 
      endselect
   elseif(wstyle.eq."2 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mClass)//'/'
      
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
         call ldt_log_msg('ERR: create_lisout_naturerun_filename -- '// &
              'Unrecognized form value')
         call LDT_endrun 
      endselect
   endif
   fname = out_fname
 end subroutine create_lisout_naturerun_filename

!BOP
! 
! !ROUTINE: convertNatureRunToLDTgrid
! \label{convertNatureRunToLDTgrid}
!
! !INTERFACE: 
 subroutine convertNatureRunToLDTgrid(n, var_inp, var_out)
! !USES:    
   use LDT_coreMod
   use LISoutNatureRun_Mod

   implicit none
! !ARGUMENTS: 
   integer         :: n 
   real            :: var_inp(LISoutNatureRunData%nc*LISoutNatureRunData%nr)
   real            :: var_out(LDT_rc%lnc(n),LDT_rc%lnr(n))
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
   logical*1       :: var_data_b(LISoutNatureRunData%nc*LISoutNatureRunData%nr)
   real            :: varobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   logical*1       :: varobs_b_ip(LISoutNatureRunData%nc*LISoutNatureRunData%nr)

   do r=1,LISoutNatureRunData%nr
      do c=1, LISoutNatureRunData%nc
         if(var_inp(c+(r-1)*LISoutNatureRunData%nc).ne.LDT_rc%udef) then 
            var_data_b(c+(r-1)*LISoutNatureRunData%nc) = .true. 
         else
            var_data_b(c+(r-1)*LISoutNatureRunData%nc) = .false.
         endif
      enddo
   enddo

   if(LDT_isLDTatAfinerResolution(n,LISoutNatureRunData%datares)) then 

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
      call bilinear_interp(LDT_rc%gridDesc(n,:),&
           var_data_b, var_inp, varobs_b_ip, varobs_ip, &
           LISoutNatureRunData%nc*LISoutNatureRunData%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           LISoutNatureRunData%w11, LISoutNatureRunData%w12, &
           LISoutNatureRunData%w21, LISoutNatureRunData%w22, &
           LISoutNatureRunData%n11, LISoutNatureRunData%n12, &
           LISoutNatureRunData%n21, LISoutNatureRunData%n22, &
           LDT_rc%udef, ios)

      call neighbor_interp(LDT_rc%gridDesc(n,:),&
           var_data_b, var_inp, varobs_b_ip, varobs_ip, &
           LISoutNatureRunData%nc*LISoutNatureRunData%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           LISoutNatureRunData%n11, LDT_rc%udef, ios)
   else
      call upscaleByAveraging(&
           LISoutNatureRunData%nc*LISoutNatureRunData%nr,&
           LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%udef, &
           LISoutNatureRunData%n11,var_data_b, var_inp, varobs_b_ip,varobs_ip)
      
   endif
   
   do r=1,LDT_rc%lnr(n)
      do c=1,LDT_rc%lnc(n)
         if(varobs_b_ip(c+(r-1)*LDT_rc%lnc(n))) then 
            var_out(c,r) = varobs_ip(c+(r-1)*LDT_rc%lnc(n))
         else
            var_out(c,r) = LDT_rc%udef
         endif
      enddo
   enddo
         
 end subroutine convertNatureRunToLDTgrid
