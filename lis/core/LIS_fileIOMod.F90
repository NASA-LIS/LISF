!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module LIS_fileIOMod
!BOP
!
! !MODULE: LIS_fileIOMod
! 
! !DESCRIPTION: 
!   This module contains a number of routines useful for various file I/O
!   operations in LIS. The module 
!   provides routines to create output directories, filenames, that can 
!   be used in the model output routines. 
!   
! !REVISION HISTORY: 
!  08 Apr 2004    James Geiger Initial Specification
!  11 Oct 2018    Nargess Memarsadeghi, cleaned up and corrected LIS_create_output_directory
!  18 Oct 2019    David Mocko, corrected creation of sub-directories for "WMO convention"
!  26 Apr 2023    Eric Kemp, added new output naming style for 557 WW for
!                 streamflow output.
!  22 May 2023    Eric Kemp, added new output naming style for 557 WW for
!                 medium range forecasts.
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none 
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_create_output_filename   ! create an output filename
  public :: LIS_create_output_directory  ! create the output directory
  public :: LIS_create_restart_filename  ! create a restart filename
  public :: LIS_create_stats_filename    ! create a stats filename
  public :: LIS_create_dapert_filename   ! create a perturbations file
  public :: LIS_create_innov_filename    ! create an innovations filename
  public :: LIS_create_incr_filename     ! create an analysis increments filename
  public :: LIS_create_daspread_filename ! create an innovations filename
  public :: LIS_create_obs_filename      ! create an observations filename
  public :: LIS_create_gain_filename
  public :: LIS_putget                   !a generic read/write method. 
  public :: LIS_readData                 !generic method for reading surface parameters
  public :: LIS_readDomainConfigSpecs    !generic method to read configurable options
                                         !for a particular surface dataset 
  public :: LIS_checkDomainExtents       !checks if the data domain extents are 
                                         !contained in the LIS running domain
  public :: LIS_read_param
  public :: LIS_read_gparam
!EOP
  
!BOP
! 
!  !ROUTINE: LIS_readData
! \label{LIS_readData}
! 
! !INTERFACE: 
  interface LIS_readData
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure read2Ddata
! 
! !DESCRIPTION: 
!  Routine to read 2d data from a binary file, with direct access format. A special
!  routine is required to read the landcover data since it includes a distribution of 
!  vegetation types at each grid point. 
!  
!EOP
  end interface
!BOP
! 
!  !ROUTINE: LIS_create_dapert_filename
! \label{LIS_create_dapert_filename}
! 
! !INTERFACE: 
  interface LIS_create_dapert_filename
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure create_dapert_filename
     module procedure create_dapert_filename_withtime
! 
! !DESCRIPTION: 
!
!  
!EOP
  end interface
  
!BOP
! 
!  !ROUTINE: LIS_create_restart_filename
! \label{LIS_create_restart_filename}
! 
! !INTERFACE: 
  interface LIS_create_restart_filename
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure create_restart_filename
     module procedure create_restart_filename_withtime
! 
! !DESCRIPTION: 
! This routine puts together a restart filename, either based on the
! LIS current time or based on the time that is provided. 
!  
!EOP
  end interface

  
!BOP 
! 
! !ROUTINE: LIS_putget
! \label{LIS_putget}
! 
! !INTERFACE:
  interface LIS_putget
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure putget_int
     module procedure putget_real
! !DESCRIPTION: 
!
!     This is a generic method to read from or write to direct access files.
!     
!     method: \newline
!     - open file \newline
!     - if iofunc = r, read buffer from file, abort on error \newline
!     - if iofunc = w, write buffer to file, abort on error \newline
!     - if iofunc = anything else, abort   
!
!
! !REVISION HISTORY:
!
!     04 aug 1999 initial version ...........................mr gayno/dnxm  
!     02 nov 2005 incorporated into LIS                      sujay kumar     
!EOP
  end interface

!BOP
! 
!  !ROUTINE: LIS_read_param
! \label{LIS_read_param}
! 
! !INTERFACE: 
  interface LIS_read_param
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure readparam_real_2d
     module procedure readparam_real_2d_rc
     module procedure readparam_int_2d
! 
! !DESCRIPTION: 
!  Routine to read parameter data from a netcdf file.  
!  
!EOP
  end interface

!BOP
! 
!  !ROUTINE: LIS_read_gparam
! \label{LIS_read_gparam}
! 
! !INTERFACE: 
  interface LIS_read_gparam
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure readgparam_real_2d
     module procedure readgparam_real_2d_rc
! 
! !DESCRIPTION: 
!  Routine to read parameter data from a netcdf file.  
!  
!EOP
  end interface

  interface LIS_create_output_filename
     module procedure create_output_filename
     module procedure create_output_filename_expected

  end interface
contains

!BOP
!
! !ROUTINE: LIS_create_output_directory
! \label{LIS_create_output_directory}
!
! !INTERFACE:
subroutine LIS_create_output_directory(mname)
! !USES:
   use LIS_coreMod, only : LIS_rc
   use LIS_logMod,  only : LIS_log_msg, LIS_logunit
   implicit none 
! !ARGUMENTS:
   character(len=*)  :: mname

!
! !DESCRIPTION:  
!  Create the output directory for the output data files. The call creates
!  a hierarchy of directories in the following format, if the directory 
!  name is not specified. 
!
!  2 level hierarchy
!  \begin{verbatim}
!   <output directory>/<model name>
!  \end{verbatim}
!  3 level hierarchy
!  \begin{verbatim}
!  <output directory>/<model name>/<yrmo>
!  \end{verbatim}
!  4 level hierarchy
!  \begin{verbatim}
!  <output directory>/<model name>/<yr>/<yrmoda>
!  \end{verbatim}
!  
!  Once the directory name is created, the subroutine issues a 
!  system call to create the structure. 
! 
!  The arguments are: 
!  \begin{description}
!   \item [mname]
!     a string describing the name of the model
!  \end{description}
!
!EOP
   character(len=4) :: cdate
   character(len=8) :: cdate1
   character(len=LIS_CONST_PATH_LEN) :: out_dname
   integer            :: try,ios
#if (!defined AIX )
   integer            :: system 
#endif
   
   ! EMK...Calls to 'system' fail when using SGI MPT as the MPI implementation
   ! on Pleiades. We replace with a C wrapper function that calls the 'mkdir' 
   ! standard POSIX function. This requires defining the C wrapper function,
   ! and specifying new variables to pass to said C function.
   integer, external :: LIS_create_subdirs
   character(len=LIS_CONST_PATH_LEN+1) :: c_string

   if(LIS_rc%wstyle.eq."4 level hierarchy") then
      out_dname = trim(LIS_rc%odir)//'/'
      out_dname = trim(out_dname)//trim(mname)//'/'
      write(unit=cdate, fmt='(i4.4)') LIS_rc%yr
      out_dname = trim(out_dname)//trim(cdate)//'/'
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      out_dname = trim(out_dname)//trim(cdate1)
   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then 
      out_dname = trim(LIS_rc%odir)//'/'
      out_dname = trim(out_dname)//trim(mname)//'/'
      write(unit=cdate1, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
      out_dname = trim(out_dname)//trim(cdate1)
   elseif(LIS_rc%wstyle.eq."2 level hierarchy") then 
      out_dname = trim(LIS_rc%odir)//'/'
      out_dname = trim(out_dname)//trim(mname)//'/'
   elseif(LIS_rc%wstyle.eq."WMO convention") then
! If output style is "WMO convention", ensure that the below
! sub-directories are created before other parts of LIS will
! try to write datasets into those sub-directories. - Mocko
      out_dname = trim(LIS_rc%odir)
      if (trim(mname).eq."SURFACEMODEL") then
         continue
      elseif (trim(mname).eq."DAPERT") then
         out_dname = trim(LIS_rc%odir)//'/'
         out_dname = trim(out_dname)//trim(mname)//'/'
      elseif (trim(mname).eq."DAOBS") then
         out_dname = trim(LIS_rc%odir)//'/'
         out_dname = trim(out_dname)//trim(mname)//'/'
         write(unit=cdate1, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
         out_dname = trim(out_dname)//trim(cdate1)
      endif
   elseif(LIS_rc%wstyle.eq."557WW streamflow convention" .or. &
        LIS_rc%wstyle .eq. "557WW medium range forecast convention") then ! EMK
      out_dname = trim(LIS_rc%odir)
      if (trim(mname).eq."SURFACEMODEL") then
         continue
      elseif (trim(mname) .eq. "ROUTING") then
         continue
      elseif (trim(mname).eq."DAPERT") then
         out_dname = trim(LIS_rc%odir)//'/'
         out_dname = trim(out_dname)//trim(mname)//'/'
      elseif (trim(mname).eq."DAOBS") then
         out_dname = trim(LIS_rc%odir)//'/'
         out_dname = trim(out_dname)//trim(mname)//'/'
         write(unit=cdate1, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
         out_dname = trim(out_dname)//trim(cdate1)
      endif
   endif

#if ( defined AIX )
   call system('mkdir -p '//trim(out_dname))
#else
   ! EMK...Calls to 'system' fail when using SGI MPT as the MPI implementation
   ! on Pleiades. We replace with a C wrapper function that calls the 'mkdir' 
   ! standard POSIX function. 
   !         ios = system('mkdir -p '//trim(out_dname))
   c_string = trim(out_dname)
   ios = LIS_create_subdirs(len_trim(c_string),trim(c_string))
#endif

   if (ios .ne. 0) then
     write(LIS_logunit,*)'[ERR] creating directory ',trim(out_dname)
     flush(LIS_logunit)
   end if

end subroutine LIS_create_output_directory


!BOP
!
! !ROUTINE: create_output_filename
! \label{create_output_filename}
!
! !INTERFACE:
subroutine create_output_filename(n, fname, model_name, odir, writeint)
! !USES:
   use LIS_coreMod,  only : LIS_rc
   use LIS_logMod,   only : LIS_log_msg, LIS_endrun, LIS_logunit
!   use LIS_forecastMod

   implicit none 
! !ARGUMENTS:
   integer, intent(in) :: n
   character(len=*), intent(out)          :: fname
   character(len=*), intent(in), optional :: model_name ! needed for gswp run
   character(len=*), intent(in), optional :: odir ! needed for gswp run
   real, intent(in), optional             :: writeint ! output writing interval
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
   character(len=10)       :: ensString
   character(len=14)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=LIS_CONST_PATH_LEN)       :: dname
   character(len=LIS_CONST_PATH_LEN), save :: out_fname
   character(len=LIS_CONST_PATH_LEN)       :: odir_temp
   integer                  :: i, c

   character(len=8) :: initdate
   character(len=2) :: inithr
   character(len=3) :: fhr
   integer :: hr
   integer :: rc
   type(ESMF_Time) :: starttime, curtime
   type(ESMF_TimeInterval) :: deltatime

   if ( present(odir) ) then
      odir_temp = odir
   else
      odir_temp = LIS_rc%odir
   endif

   if(LIS_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(odir_temp)//'/'
      dname = trim(dname)//trim(model_name)//'/'
      
      write(unit=cdate, fmt='(i4.4)') LIS_rc%yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LIS_HIST_'//trim(cdate1)

!      if(LIS_rc%forecastMode.eq.1) then 
!         write(unit=ensString,fmt='(a2,i3.3)') '.e',LIS_forecast_struc(n)%iterId
!         out_fname = trim(out_fname)//trim(ensString)
!      endif
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      

      select case ( LIS_rc%wout )
      case ( "binary" )
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("distributed binary")
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif         
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call lis_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized output format')
         call LIS_endrun 
      endselect

   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(odir_temp)//'/'
      dname = trim(dname)//trim(model_name)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)
      
!      if(LIS_rc%forecastMode.eq.1) then 
!         write(unit=ensString,fmt='(a2,i3.3)') '.e',LIS_forecast_struc(n)%iterId
!         out_fname = trim(out_fname)//trim(ensString)
!      endif

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( LIS_rc%wout )
      case ("binary")
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'            
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ( "distributed binary" )
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'            
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif         
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call lis_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized LIS_rc%wout value')
         call LIS_endrun 
      endselect

   elseif(LIS_rc%wstyle.eq."2 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(odir_temp)//'/'
      dname = trim(dname)//trim(model_name)//'/'
      
      out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)
      
!      if(LIS_rc%forecastMode.eq.1) then 
!         write(unit=ensString,fmt='(a2,i3.3)') '.e',LIS_forecast_struc(n)%iterId
!         out_fname = trim(out_fname)//trim(ensString)
!      endif

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( LIS_rc%wout )
      case ("binary")
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'            
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ( "distributed binary" )
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'            
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif         
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call lis_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized LIS_rc%wout value')
         call LIS_endrun 
      endselect

   elseif(LIS_rc%wstyle.eq."WMO convention") then 
      write(unit=fint,fmt='(i2.2)') nint(writeint)/3600
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      
      write(unit=cdate, fmt='(i2.2, i2.2)') LIS_rc%hr, LIS_rc%mn
      
      if(LIS_rc%lis_map_proj.eq."polar") then 
         fproj = 'P'
         print *,"fres ",LIS_rc%gridDesc(n, 9)
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."lambert") then 
         fproj = 'L'
         print *,"fres ",LIS_rc%gridDesc(n, 9)
         write(unit=fres, fmt='(f2.0)') LIS_rc%gridDesc(n, 9)
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."mercator") then 
         fproj = 'M'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)
         fres = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."gaussian") then 
         fproj = 'G'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)*100        
         fres = '0P'//trim(fres)//'DEG'
      else
         fproj = 'C'
         write(unit=fres, fmt='(i10)') nint(LIS_rc%gridDesc(n,10)*100)
         read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
         c = 0 
         do i=1,10
            if(fres1(i).ne.' '.and.c==0) c = i
         enddo
         if (LIS_rc%gridDesc(n,10) .lt. 0.1) then
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
      
      dname = trim(odir_temp)//&
           '/PS.AFWA_SC.'//trim(LIS_rc%security_class)//&
           '_DI.'//trim(LIS_rc%distribution_class)//&
           '_DC.'//trim(LIS_rc%data_category)//'_GP.LIS_GR.'//&
           trim(fproj)//trim(fres2)//'_AR.'//trim(LIS_rc%area_of_data)//&
           '_PA.'//trim(fint)//'-HR-SUM_DD.'//&
           trim(cdate1)//'_DT.'//trim(cdate)//'_DF'
      select case (LIS_rc%wout)
      case ("binary")
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(dname)//'.DAT'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(dname)//'.DAT'
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.DAT'
         endif
      case ("grib1")
         out_fname = trim(dname)//'.GR1'
      case ("netcdf")
         out_fname = trim(dname)//'.nc'
      case ("grib2")
         out_fname = trim(dname)//'.GR2'
      case default            
      end select
   elseif(LIS_rc%wstyle.eq."557WW streamflow convention") then ! EMK
      write(unit=fint,fmt='(i2.2)') nint(writeint)/3600
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, LIS_rc%da

      write(unit=cdate, fmt='(i2.2, i2.2)') LIS_rc%hr, LIS_rc%mn

      if(LIS_rc%lis_map_proj.eq."polar") then
         fproj = 'P'
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."lambert") then
         fproj = 'L'
         write(unit=fres, fmt='(f2.0)') LIS_rc%gridDesc(n, 9)
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."mercator") then
         fproj = 'M'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)
         fres = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."gaussian") then
         fproj = 'G'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)*100
         fres = '0P'//trim(fres)//'DEG'
      else
         fproj = 'C'
         write(unit=fres, fmt='(i10)') nint(LIS_rc%gridDesc(n,10)*100)
         read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
         c = 0
         do i=1,10
            if(fres1(i).ne.' '.and.c==0) c = i
         enddo
         if (LIS_rc%gridDesc(n,10) .lt. 0.1) then
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

      !dname = trim(odir_temp)//&
      !     '/PS.AFWA_SC.'//trim(LIS_rc%security_class)//&
      !     '_DI.'//trim(LIS_rc%distribution_class)//&
      !     '_DC.'//trim(LIS_rc%data_category)//'_GP.LIS_GR.'//&
      !     trim(fproj)//trim(fres2)//'_AR.'//trim(LIS_rc%area_of_data)//&
      !     '_PA.'//trim(fint)//'-HR-SUM_DD.'//&
      !     trim(cdate1)//'_DT.'//trim(cdate)//'_DF'
      dname = trim(odir_temp) &
           //'/PS.557WW' &
           //'_SC.'//trim(LIS_rc%security_class) &
           //'_DI.'//trim(LIS_rc%distribution_class)
      if (LIS_rc%lsm .eq. "Noah.3.9") then
         dname = trim(dname) &
           //'_GP.'//'LIS-NOAH'
      else if (LIS_rc%lsm .eq. "Noah-MP.4.0.1") then
         dname = trim(dname) &
           //'_GP.'//'LIS-NOAHMP'
      else if (LIS_rc%lsm .eq. "JULES.5.0") then
         dname = trim(dname) &
           //'_GP.'//'LIS-JULES'
      else
         write(LIS_logunit,*) &
              '[ERR] Invalid Land surface model for 557WW streamflow convention ', &
              trim(LIS_rc%lsm)
         call LIS_endrun()
      end if
      if (LIS_rc%routingmodel .eq. "HYMAP2 router") then
         dname = trim(dname) &
           //'-HYMAP'
      else if (LIS_rc%routingmodel .eq. "RAPID router") then
         dname = trim(dname) &
           //'-RAPID'
      else if (LIS_rc%routingmodel .eq. "none") then
         continue
      else
         write(LIS_logunit,*) &
              '[ERR] Invalid Routing model for 557WW streamflow convention ', &
              trim(LIS_rc%routingmodel)
         call LIS_endrun()
      end if

      dname = trim(dname) &
           //'_GR.'//trim(fproj)//trim(fres2) &
           //'_AR.'//trim(LIS_rc%area_of_data) &
           //'_PA.'//trim(model_name) &
           //'_DD.'//trim(cdate1) &
           //'_DT.'//trim(cdate) &
           //'_DF'

      select case (LIS_rc%wout)
      case ("binary")
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(dname)//'.DAT'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(dname)//'.DAT'
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.DAT'
         endif
      case ("grib1")
         out_fname = trim(dname)//'.GR1'
      case ("netcdf")
         out_fname = trim(dname)//'.NC' ! EMK
      case ("grib2")
         out_fname = trim(dname)//'.GR2'
      case default            
      end select

   elseif(LIS_rc%wstyle.eq."557WW medium range forecast convention") then ! EMK
      write(unit=fint,fmt='(i2.2)') nint(writeint)/3600
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, LIS_rc%da

      write(unit=cdate, fmt='(i2.2, i2.2)') LIS_rc%hr, LIS_rc%mn

      if(LIS_rc%lis_map_proj.eq."polar") then
         fproj = 'P'
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."lambert") then
         fproj = 'L'
         write(unit=fres, fmt='(f2.0)') LIS_rc%gridDesc(n, 9)
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."mercator") then
         fproj = 'M'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)
         fres = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."gaussian") then
         fproj = 'G'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)*100
         fres = '0P'//trim(fres)//'DEG'
      else
         fproj = 'C'
         write(unit=fres, fmt='(i10)') nint(LIS_rc%gridDesc(n,10)*100)
         read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
         c = 0
         do i=1,10
            if(fres1(i).ne.' '.and.c==0) c = i
         enddo
         if (LIS_rc%gridDesc(n,10) .lt. 0.1) then
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

      dname = trim(odir_temp) &
           //'/PS.557WW' &
           //'_SC.'//trim(LIS_rc%security_class) &
           //'_DI.'//trim(LIS_rc%distribution_class)
      if (LIS_rc%lsm .eq. "Noah.3.9") then
         dname = trim(dname) &
           //'_GP.'//'LIS-MR-NOAH'
      else if (LIS_rc%lsm .eq. "Noah-MP.4.0.1") then
         dname = trim(dname) &
           //'_GP.'//'LIS-MR-NOAHMP'
      else if (LIS_rc%lsm .eq. "JULES.5.0") then
         dname = trim(dname) &
           //'_GP.'//'LIS-MR-JULES'
      else
         write(LIS_logunit,*) &
              '[ERR] Invalid Land surface model for ', &
              '557WW medium range forecast convention ', &
              trim(LIS_rc%lsm)
         call LIS_endrun()
      end if
      if (LIS_rc%routingmodel .eq. "HYMAP2 router") then
         dname = trim(dname) &
           //'-HYMAP'
      else if (LIS_rc%routingmodel .eq. "RAPID router") then
         dname = trim(dname) &
           //'-RAPID'
      else if (LIS_rc%routingmodel .eq. "none") then
         continue
      else
         write(LIS_logunit,*) &
              '[ERR] Invalid Routing model for ', &
              '557WW medium range forecast convention ', &
              trim(LIS_rc%routingmodel)
         call LIS_endrun()
      end if

      write(unit=initdate, fmt='(i4.4, i2.2, i2.2)') &
           LIS_rc%syr, LIS_rc%smo, LIS_rc%sda
      write(unit=inithr, fmt='(i2.2)') LIS_rc%shr
      call ESMF_TimeSet(starttime, &
           yy=LIS_rc%syr, mm=LIS_rc%smo, dd=LIS_rc%sda, &
           h=LIS_rc%shr, m=LIS_rc%smn, s=LIS_rc%sss, rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         write(LIS_logunit,*)'[ERR] Cannot set starttime object!'
         call LIS_endrun()
      end if
      call ESMF_TimeSet(curtime, &
           yy=LIS_rc%yr, mm=LIS_rc%mo, dd=LIS_rc%da, &
           h=LIS_rc%hr, m=LIS_rc%mn, s=LIS_rc%ss, rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         write(LIS_logunit,*)'[ERR] Cannot set curtime object!'
         call LIS_endrun()
      end if
      deltatime = curtime - starttime
      call ESMF_TimeIntervalGet(deltatime, h=hr, rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         write(LIS_logunit,*)'[ERR] Cannot get hr from deltatime!'
         call LIS_endrun()
      end if
      write(unit=fhr, fmt='(i3.3)') hr
      dname = trim(dname) &
           //'_GR.'//trim(fproj)//trim(fres2) &
           //'_AR.'//trim(LIS_rc%area_of_data) &
           //'_PA.'//trim(model_name) &
           //'_DD.'//trim(initdate) &
           //'_CY.'//trim(inithr) &
           //'_FH.'//trim(fhr)

      select case (LIS_rc%wout)
      case ("binary")
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(dname)//'_DF.DAT'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(dname)//'_DF.DAT'
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'_DF.DAT'
         endif
      case ("grib1")
         out_fname = trim(dname)//'_DF.GR1'
      case ("netcdf")
         out_fname = trim(dname)//'_DF.NC'
      case ("grib2")
         out_fname = trim(dname)//'_DF.GR2'
      case default            
      end select

   endif
   fname = out_fname
 end subroutine create_output_filename



!BOP
!
! !ROUTINE: create_output_filename_expected
! \label{create_output_filename_expected}
!
! !INTERFACE:
subroutine create_output_filename_expected(n, fname, wout, flag, model_name, odir,&
     writeint)
! !USES:
   use LIS_coreMod
   use LIS_surfaceModelDataMod
   use LIS_logMod
   use LIS_timeMgrMod

   implicit none 
! !ARGUMENTS:
   integer, intent(in) :: n
   character(len=*), intent(out)          :: fname
   character(len=*), intent(in)           :: wout
   logical         , intent(in)           :: flag 
   character(len=*), intent(in), optional :: model_name ! needed for gswp run
   character(len=*), intent(in), optional :: odir ! needed for gswp run
   real, intent(in), optional             :: writeint ! output writing interval
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
   character(len=8) :: initdate
   character(len=2) :: inithr
   character(len=3) :: fhr
   integer :: rc
   type(ESMF_Time) :: starttime, curtime
   type(ESMF_TimeInterval) :: deltatime
   character(len=LIS_CONST_PATH_LEN)       :: dname
   character(len=LIS_CONST_PATH_LEN), save :: out_fname
   character(len=LIS_CONST_PATH_LEN)       :: odir_temp
   integer                  :: i, c
   integer                  :: yr, mo, da, hr, mn, ss
   type(ESMF_Time)          :: currTime, outTime
   type(ESMF_TimeInterval)  :: outTS
   integer                  :: status

   if ( present(odir) ) then
      odir_temp = odir
   else
      odir_temp = LIS_rc%odir
   endif   

   call ESMF_TimeIntervalSet(outTS, s=nint(LIS_sfmodel_struc(n)%outInterval), &
        rc=status)
! Find the next output time (only supported upto 1 day intervals)

   if(LIS_sfmodel_struc(n)%outInterval.eq.86400) then 
      
      call ESMF_TimeSet(currTime, yy = LIS_rc%yr, &
           mm = LIS_rc%mo, &
           dd = LIS_rc%da, &
           h  = 0, &
           m  = 0, & 
           s  = 0, &
           calendar = LIS_calendar, & 
           rc = status)
      call LIS_verify(status,'error in ESMF_TimeSet:currTime in LIS_fileIOMod')
      outTime = currTime + outTS
      call ESMF_TimeGet(outTime, yy = yr, &
           mm = mo, &
           dd = da, &
           h  = hr, &
           m  = mn, & 
           s  = ss, &
           calendar = LIS_calendar, & 
           rc = status)
      call LIS_verify(status,'error in ESMF_TimeSet:currTime in LIS_fileIOMod')
   elseif(LIS_sfmodel_struc(n)%outInterval.eq.10800) then 
      
      call ESMF_TimeSet(currTime, yy = LIS_rc%yr, &
           mm = LIS_rc%mo, &
           dd = LIS_rc%da, &
           h  = 3*(int(real(LIS_rc%hr)/3.0)),&
           m  = 0, & 
           s  = 0, &
           calendar = LIS_calendar, & 
           rc = status)
      call LIS_verify(status,'error in ESMF_TimeSet:currTime in LIS_fileIOMod')
      outTime = currTime + outTS
      call ESMF_TimeGet(outTime, yy = yr, &
           mm = mo, &
           dd = da, &
           h  = hr, &
           m  = mn, & 
           s  = ss, &
           calendar = LIS_calendar, & 
           rc = status)
      call LIS_verify(status,'error in ESMF_TimeSet:currTime in LIS_fileIOMod')
   elseif(LIS_sfmodel_struc(n)%outInterval.eq.3600) then 
      
      call ESMF_TimeSet(currTime, yy = LIS_rc%yr, &
           mm = LIS_rc%mo, &
           dd = LIS_rc%da, &
           h  = LIS_rc%hr, &
           m  = 0, & 
           s  = 0, &
           calendar = LIS_calendar, & 
           rc = status)
      call LIS_verify(status,'error in ESMF_TimeSet:currTime in LIS_fileIOMod')
      outTime = currTime + outTS
      call ESMF_TimeGet(outTime, yy = yr, &
           mm = mo, &
           dd = da, &
           h  = hr, &
           m  = mn, & 
           s  = ss, &
           calendar = LIS_calendar, & 
           rc = status)
      call LIS_verify(status,'error in ESMF_TimeSet:currTime in LIS_fileIOMod')
   elseif(LIS_sfmodel_struc(n)%outInterval.eq.1800) then 
      
      call ESMF_TimeSet(currTime, yy = LIS_rc%yr, &
           mm = LIS_rc%mo, &
           dd = LIS_rc%da, &
           h  = LIS_rc%hr, &
           m  = LIS_rc%mn, & 
           s  = 0, &
           calendar = LIS_calendar, & 
           rc = status)
      call LIS_verify(status,'error in ESMF_TimeSet:currTime in LIS_fileIOMod')
      outTime = currTime + outTS
      call ESMF_TimeGet(outTime, yy = yr, &
           mm = mo, &
           dd = da, &
           h  = hr, &
           m  = mn, & 
           s  = ss, &
           calendar = LIS_calendar, & 
           rc = status)
      call LIS_verify(status,'error in ESMF_TimeSet:currTime in LIS_fileIOMod')
   else
      write(LIS_logunit,*) '[ERR] This output interval is not supported for '
      write(LIS_logunit,*) '[ERR] for computing expected output locations'
      call LIS_endrun
   endif

   if(LIS_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           yr, mo, &
           da, hr,mn
      
      dname = trim(odir_temp)//'/'
      dname = trim(dname)//trim(model_name)//'/'
      
      write(unit=cdate, fmt='(i4.4)') yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') yr, mo, da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LIS_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( wout )
      case ( "binary" )
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'            
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ( "distributed binary" )
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'            
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif         
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call lis_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized output format')
         call LIS_endrun 
      endselect
   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           yr, mo, &
           da, hr,mn
      
      dname = trim(odir_temp)//'/'
      dname = trim(dname)//trim(model_name)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') yr, mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( wout )
      case ("binary")
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'            
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ( "distributed binary" )
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'            
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif         
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call lis_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized wout value')
         call LIS_endrun 
      endselect
   elseif(LIS_rc%wstyle.eq."2 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           yr, mo, &
           da, hr,mn
      
      dname = trim(odir_temp)//'/'
      dname = trim(dname)//trim(model_name)//'/'
      
      out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( wout )
      case ("binary")
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'            
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ( "distributed binary" )
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'            
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif                  
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call lis_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized wout value')
         call LIS_endrun 
      endselect
   elseif(LIS_rc%wstyle.eq."WMO convention") then 
      write(unit=fint,fmt='(i2.2)') nint(writeint)/3600
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
           yr, mo, da
      
      write(unit=cdate, fmt='(i2.2, i2.2)') hr, mn
      
      if(LIS_rc%lis_map_proj.eq."polar") then 
         fproj = 'P'
         print *,"fres ",LIS_rc%gridDesc(n, 9)
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."lambert") then 
         fproj = 'L'
         print *,"fres ",LIS_rc%gridDesc(n, 9)
         write(unit=fres, fmt='(f2.0)') LIS_rc%gridDesc(n, 9)
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."mercator") then 
         fproj = 'M'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)
         fres = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."gaussian") then 
         fproj = 'G'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)*100        
         fres = '0P'//trim(fres)//'DEG'
      else
         fproj = 'C'
         write(unit=fres, fmt='(i10)') nint(LIS_rc%gridDesc(n,10)*100)
         read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
         c = 0 
         do i=1,10
            if(fres1(i).ne.' '.and.c==0) c = i
         enddo
         if (LIS_rc%gridDesc(n,10) .lt. 0.1) then
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
      
      dname = trim(odir_temp)//&
           '/PS.AFWA_SC.'//trim(LIS_rc%security_class)//&
           '_DI.'//trim(LIS_rc%distribution_class)//&
           '_DC.'//trim(LIS_rc%data_category)//'_GP.LIS_GR.'//&
           trim(fproj)//trim(fres2)//'_AR.'//trim(LIS_rc%area_of_data)//&
           '_PA.'//trim(fint)//'-HR-SUM_DD.'//&
           trim(cdate1)//'_DT.'//trim(cdate)//'_DF'
      select case (wout)
      case ("binary")
         if(LIS_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(dname)//'.DAT'
         elseif(LIS_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(dname)//'.DAT'
         elseif(LIS_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.DAT'
         endif
      case ("grib1")
         out_fname = trim(dname)//'.GR1'
      case ("netcdf")
         out_fname = trim(dname)//'.nc'
      case ("grib2")
         out_fname = trim(dname)//'.GR2'
      case default            
      end select

   elseif(LIS_rc%wstyle.eq."557WW streamflow convention") then ! EMK
      write(unit=fint,fmt='(i2.2)') nint(writeint)/3600
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
           yr, mo, da

      write(unit=cdate, fmt='(i2.2, i2.2)') hr, mn

      if(LIS_rc%lis_map_proj.eq."polar") then
         fproj = 'P'
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."lambert") then
         fproj = 'L'
         write(unit=fres, fmt='(f2.0)') LIS_rc%gridDesc(n, 9)
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."mercator") then
         fproj = 'M'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)
         fres = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."gaussian") then
         fproj = 'G'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)*100
         fres = '0P'//trim(fres)//'DEG'
      else
         fproj = 'C'
         write(unit=fres, fmt='(i10)') nint(LIS_rc%gridDesc(n,10)*100)
         read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
         c = 0
         do i=1,10
            if(fres1(i).ne.' '.and.c==0) c = i
         enddo
         if (LIS_rc%gridDesc(n,10) .lt. 0.1) then
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

      !dname = trim(odir_temp)//&
      !     '/PS.AFWA_SC.'//trim(LIS_rc%security_class)//&
      !     '_DI.'//trim(LIS_rc%distribution_class)//&
      !     '_DC.'//trim(LIS_rc%data_category)//'_GP.LIS_GR.'//&
      !     trim(fproj)//trim(fres2)//'_AR.'//trim(LIS_rc%area_of_data)//&
      !     '_PA.'//trim(fint)//'-HR-SUM_DD.'//&
      !     trim(cdate1)//'_DT.'//trim(cdate)//'_DF'
      dname = trim(odir_temp) &
           //'/PS.557WW' &
           //'_SC.'//trim(LIS_rc%security_class) &
           //'_DI.'//trim(LIS_rc%distribution_class)
      if (LIS_rc%lsm .eq. "Noah.3.9") then
         dname = trim(dname) &
           //'_GP.'//'LIS-NOAH'
      else if (LIS_rc%lsm .eq. "Noah-MP.4.0.1") then
         dname = trim(dname) &
           //'_GP.'//'LIS-NOAHMP'
      else if (LIS_rc%lsm .eq. "JULES.5.0") then
         dname = trim(dname) &
           //'_GP.'//'LIS-JULES'
      else
         write(LIS_logunit,*) &
              '[ERR] Invalid Land surface model for 557WW streamflow convention ', &
              trim(LIS_rc%lsm)
         call LIS_endrun()
      end if
      if (LIS_rc%routingmodel .eq. "HYMAP2 router") then
         dname = trim(dname) &
           //'-HYMAP'
      else if (LIS_rc%routingmodel .eq. "RAPID router") then
         dname = trim(dname) &
           //'-RAPID'
      else if (LIS_rc%routingmodel .eq. "none") then
         continue
      else
         write(LIS_logunit,*) &
              '[ERR] Invalid Routing model for 557WW streamflow convention ', &
              trim(LIS_rc%routingmodel)
         call LIS_endrun()
      end if

      dname = trim(dname) &
           //'_GR.'//trim(fproj)//trim(fres2) &
           //'_AR.'//trim(LIS_rc%area_of_data) &
           //'_PA.'//trim(model_name) &
           //'_DD.'//trim(cdate1) &
           //'_DT.'//trim(cdate) &
           //'_DF'

      select case (wout)
      case ("binary")
         if(LIS_rc%wopt.eq."1d tilespace") then
            out_fname = trim(dname)//'.DAT'
         elseif(LIS_rc%wopt.eq."2d gridspace") then
            out_fname = trim(dname)//'.DAT'
         elseif(LIS_rc%wopt.eq."1d gridspace") then
            out_fname = trim(out_fname)//'.DAT'
         endif
      case ("grib1")
         out_fname = trim(dname)//'.GR1'
      case ("netcdf")
         out_fname = trim(dname)//'.NC'
      case ("grib2")
         out_fname = trim(dname)//'.GR2'
      case default
      end select

   elseif(LIS_rc%wstyle.eq. &
        "557WW medium range forecast convention") then ! EMK
      write(unit=fint,fmt='(i2.2)') nint(writeint)/3600
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
           yr, mo, da

      write(unit=cdate, fmt='(i2.2, i2.2)') hr, mn

      if(LIS_rc%lis_map_proj.eq."polar") then
         fproj = 'P'
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."lambert") then
         fproj = 'L'
         write(unit=fres, fmt='(f2.0)') LIS_rc%gridDesc(n, 9)
         if (LIS_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LIS_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LIS_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."mercator") then
         fproj = 'M'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)
         fres = trim(fres)//'KM'
      elseif(LIS_rc%lis_map_proj.eq."gaussian") then
         fproj = 'G'
         write(unit=fres, fmt='(i2.2)') LIS_rc%gridDesc(n, 9)*100
         fres = '0P'//trim(fres)//'DEG'
      else
         fproj = 'C'
         write(unit=fres, fmt='(i10)') nint(LIS_rc%gridDesc(n,10)*100)
         read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
         c = 0
         do i=1,10
            if(fres1(i).ne.' '.and.c==0) c = i
         enddo
         if (LIS_rc%gridDesc(n,10) .lt. 0.1) then
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

      dname = trim(odir_temp) &
           //'/PS.557WW' &
           //'_SC.'//trim(LIS_rc%security_class) &
           //'_DI.'//trim(LIS_rc%distribution_class)
      if (LIS_rc%lsm .eq. "Noah.3.9") then
         dname = trim(dname) &
           //'_GP.'//'LIS-MR-NOAH'
      else if (LIS_rc%lsm .eq. "Noah-MP.4.0.1") then
         dname = trim(dname) &
           //'_GP.'//'LIS-MR-NOAHMP'
      else if (LIS_rc%lsm .eq. "JULES.5.0") then
         dname = trim(dname) &
           //'_GP.'//'LIS-MR-JULES'
      else
         write(LIS_logunit,*) &
              '[ERR] Invalid Land surface model for ', &
              '557WW medium range forecast convention ', &
              trim(LIS_rc%lsm)
         call LIS_endrun()
      end if
      if (LIS_rc%routingmodel .eq. "HYMAP2 router") then
         dname = trim(dname) &
           //'-HYMAP'
      else if (LIS_rc%routingmodel .eq. "RAPID router") then
         dname = trim(dname) &
           //'-RAPID'
      else if (LIS_rc%routingmodel .eq. "none") then
         continue
      else
         write(LIS_logunit,*) &
              '[ERR] Invalid Routing model for ', &
              '557WW medium range forecast convention ', &
              trim(LIS_rc%routingmodel)
         call LIS_endrun()
      end if

      write(unit=initdate, fmt='(i4.4, i2.2, i2.2)') &
           LIS_rc%syr, LIS_rc%smo, LIS_rc%sda
      write(unit=inithr, fmt='(i2.2)') LIS_rc%shr
      call ESMF_TimeSet(starttime, &
           yy=LIS_rc%syr, mm=LIS_rc%smo, dd=LIS_rc%sda, &
           h=LIS_rc%shr, m=LIS_rc%smn, s=LIS_rc%sss, rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         write(LIS_logunit,*)'[ERR] Cannot set starttime object!'
         call LIS_endrun()
      end if
      call ESMF_TimeSet(curtime, &
           yy=LIS_rc%yr, mm=LIS_rc%mo, dd=LIS_rc%da, &
           h=LIS_rc%hr, m=LIS_rc%mn, s=LIS_rc%ss, rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         write(LIS_logunit,*)'[ERR] Cannot set curtime object!'
         call LIS_endrun()
      end if
      deltatime = curtime - starttime
      call ESMF_TimeIntervalGet(deltatime, h=hr, rc=rc)
      if (rc .ne. ESMF_SUCCESS) then
         write(LIS_logunit,*)'[ERR] Cannot get hr from deltatime!'
         call LIS_endrun()
      end if
      write(unit=fhr, fmt='(i3.3)') hr
      dname = trim(dname) &
           //'_GR.'//trim(fproj)//trim(fres2) &
           //'_AR.'//trim(LIS_rc%area_of_data) &
           //'_PA.'//trim(model_name) &
           //'_DD.'//trim(initdate) &
           //'_CY.'//trim(inithr) &
           //'_FH.'//trim(fhr)

      select case (wout)
      case ("binary")
         if(LIS_rc%wopt.eq."1d tilespace") then
            out_fname = trim(dname)//'_DF.DAT'
         elseif(LIS_rc%wopt.eq."2d gridspace") then
            out_fname = trim(dname)//'_DF.DAT'
         elseif(LIS_rc%wopt.eq."1d gridspace") then
            out_fname = trim(out_fname)//'_DF.DAT'
         endif
      case ("grib1")
         out_fname = trim(dname)//'_DF.GR1'
      case ("netcdf")
         out_fname = trim(dname)//'_DF.NC'
      case ("grib2")
         out_fname = trim(dname)//'_DF.GR2'
      case default
      end select

   endif
   fname = out_fname
 end subroutine create_output_filename_expected


!BOP
!
! !ROUTINE: create_dapert_filename
! \label{create_dapert_filename}
!
! !INTERFACE:
subroutine create_dapert_filename(n, fname)
! !USES:
   use LIS_coreMod,  only : LIS_rc
   use LIS_logMod,   only : LIS_log_msg, LIS_endrun

   implicit none 
! !ARGUMENTS:
   integer, intent(in) :: n
   character(len=*), intent(out)          :: fname

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
   character(len=LIS_CONST_PATH_LEN)       :: dname
   character(len=LIS_CONST_PATH_LEN), save :: out_fname
   integer                  :: i, c

   if(LIS_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//'DAPERT/'
      
      write(unit=cdate, fmt='(i4.4)') LIS_rc%yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LIS_DAPERT_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.bin'
   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//'DAPERT/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_DAPERT_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.bin'
   elseif(LIS_rc%wstyle.eq."2 level hierarchy".or.&
        LIS_rc%wstyle.eq."WMO convention" .or. &
        LIS_rc%wstyle.eq."557WW streamflow convention" .or. &
        LIS_rc%wstyle.eq."557WW medium range forecast convention") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//'DAPERT/'
      
      out_fname = trim(dname)//'LIS_DAPERT_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.bin'
   endif
   fname = out_fname
 end subroutine create_dapert_filename

!BOP
!
! !ROUTINE: create_dapert_filename_withtime
! \label{create_dapert_filename_withtime}
!
! !INTERFACE:
subroutine create_dapert_filename_withtime(n, fname, yr, mo, da, hr, mn, ss)
! !USES:
   use LIS_coreMod,  only : LIS_rc
   use LIS_logMod,   only : LIS_log_msg, LIS_endrun

   implicit none 
! !ARGUMENTS:
   integer, intent(in) :: n
   character(len=*), intent(out)          :: fname
   integer                                :: yr, mo, da, hr, mn, ss

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
   character(len=LIS_CONST_PATH_LEN)       :: dname
   character(len=LIS_CONST_PATH_LEN), save :: out_fname
   integer                  :: i, c

   if(LIS_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           yr, mo, &
           da, hr,mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//'DAPERT/'
      
      write(unit=cdate, fmt='(i4.4)') yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') yr, mo, da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LIS_DAPERT_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.bin'
   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           yr, mo, &
           da, hr,mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//'DAPERT/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') yr, mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_DAPERT_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.bin'
   elseif(LIS_rc%wstyle.eq."2 level hierarchy".or.&
        LIS_rc%wstyle.eq."WMO convention" .or. &
        LIS_rc%wstyle.eq."557WW streamflow convention" .or. &
        LIS_rc%wstyle.eq."557WW medium range forecast convention") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           yr, mo, &
           da, hr,mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//'DAPERT/'
      
      out_fname = trim(dname)//'LIS_DAPERT_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.bin'
   endif
   fname = out_fname
 end subroutine create_dapert_filename_withtime

!BOP
!
! !ROUTINE: create_restart_filename
! \label{create_restart_filename}
!
! !INTERFACE:
 subroutine create_restart_filename(n, fname,dir_name,model_name,wformat)
! !USES:
   use LIS_coreMod,    only : LIS_rc

   implicit none 
! !ARGUMENTS:
   integer, intent(in)           :: n
   character(len=*), intent(out) :: fname
   character(len=*), intent(in)  :: dir_name 
   character(len=*), intent(in)  :: model_name 
   character(len=*), optional    :: wformat

! !DESCRIPTION:  
!  Create the file name for the restart data files.  The convention used
!  in LIS creates a restart filename in the following format. 
!
!  \begin{verbatim}
!  <output directory>/EXP<expno>/<dir name>/<yr>/<yrmoda>/<yrmodahrmn>.<extn>
!  \end{verbatim}
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest
!   \item [fname]
!     the created file name. 
!   \item [model\_name]
!    string describing the name of the model 
!   \item [wformat]
!    optional argument specifying the restart file format. 
!  \end{description}
! 
!EOP
   character(len=10)        :: cdate
   character(len=12)        :: cdate1
   character(len=LIS_CONST_PATH_LEN)       :: dname
   character(len=LIS_CONST_PATH_LEN)       :: out_fname
   character*50             :: wformat_temp
   integer                  :: yr, mo, da, hr, mn, ss

   if ( present(wformat) ) then
      wformat_temp = wformat
   else
      wformat_temp = LIS_rc%wout
   endif

   yr = LIS_rc%yr
   mo = LIS_rc%mo
   da = LIS_rc%da
   hr = LIS_rc%hr
   mn = LIS_rc%mn
   ss = LIS_rc%ss

   if(LIS_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr, mo, da, hr,mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(dir_name)//'/'
      
      write(unit=cdate, fmt='(i4.4)') yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
           yr, mo, da
      dname = trim(dname)//trim(cdate)//'/'

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      select case (wformat_temp)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
!         if(LIS_rc%wopt.eq."1d tilespace") then 
!            out_fname = trim(out_fname)//'.ts4r'
!         elseif(LIS_rc%wopt.eq."2d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         elseif(LIS_rc%wopt.eq."1d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select
      fname = out_fname

   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr,mo,da,hr,mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(dir_name)//'/'
      
      write(unit=cdate, fmt='(i4.4,i2.2)') yr, mo
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      select case (wformat_temp)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
!         if(LIS_rc%wopt.eq."1d tilespace") then 
!            out_fname = trim(out_fname)//'.ts4r'
!         elseif(LIS_rc%wopt.eq."2d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         elseif(LIS_rc%wopt.eq."1d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select
      fname = out_fname
      
   elseif(LIS_rc%wstyle.eq."2 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr,mo,da,hr,mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(dir_name)//'/'
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      select case (wformat_temp)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
!         if(LIS_rc%wopt.eq."1d tilespace") then 
!            out_fname = trim(out_fname)//'.ts4r'
!         elseif(LIS_rc%wopt.eq."2d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         elseif(LIS_rc%wopt.eq."1d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select
      fname = out_fname
      
   elseif(LIS_rc%wstyle.eq."WMO convention" .or. &
        LIS_rc%wstyle.eq."557WW streamflow convention" .or. &
        LIS_rc%wstyle.eq."557WW medium range forecast convention") then 

      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr,mo,da,hr,mn
      
      dname = LIS_rc%odir
         
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'/LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      
      select case (wformat_temp)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
!         if(LIS_rc%wopt.eq."1d tilespace") then 
!            out_fname = trim(dname)//'.DAT'
!         elseif(LIS_rc%wopt.eq."2d gridspace") then 
!            out_fname = trim(dname)//'.DAT'
!         elseif(LIS_rc%wopt.eq."1d gridspace") then 
!            out_fname = trim(out_fname)//'.DAT'
!         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select

      fname = out_fname
   endif
 end subroutine create_restart_filename

!BOP
!
! !ROUTINE: create_restart_filename_withtime
! \label{create_restart_filename_withtime}
!
! !INTERFACE:
 subroutine create_restart_filename_withtime(n, fname,dir_name,model_name,&
      yr,mo,da,hr,mn,ss,wformat)
! !USES:
   use LIS_coreMod,    only : LIS_rc

   implicit none 
! !ARGUMENTS:
   integer, intent(in)           :: n
   character(len=*), intent(out) :: fname
   character(len=*), intent(in)  :: dir_name 
   character(len=*), intent(in)  :: model_name 
   integer         , intent(in)  :: yr, mo, da, hr, mn, ss
   character(len=*), optional    :: wformat

! !DESCRIPTION:  
!  Create the file name for the restart data files.  The convention used
!  in LIS creates a restart filename in the following format. 
!
!  \begin{verbatim}
!  <output directory>/EXP<expno>/<dir name>/<yr>/<yrmoda>/<yrmodahrmn>.<extn>
!  \end{verbatim}
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest
!   \item [fname]
!     the created file name. 
!   \item [model\_name]
!    string describing the name of the model 
!   \item [wformat]
!    optional argument specifying the restart file format. 
!  \end{description}
! 
!EOP
   character(len=10)        :: cdate
   character(len=12)        :: cdate1
   character(len=LIS_CONST_PATH_LEN)       :: dname
   character(len=LIS_CONST_PATH_LEN)       :: out_fname
   character*50             :: wformat_temp


   if ( present(wformat) ) then
      wformat_temp = wformat
   else
      wformat_temp = LIS_rc%wout
   endif

   if(LIS_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr, mo, da, hr,mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(dir_name)//'/'
      
      write(unit=cdate, fmt='(i4.4)') yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
           yr, mo, da
      dname = trim(dname)//trim(cdate)//'/'

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      select case (wformat_temp)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
!         if(LIS_rc%wopt.eq."1d tilespace") then 
!            out_fname = trim(out_fname)//'.ts4r'
!         elseif(LIS_rc%wopt.eq."2d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         elseif(LIS_rc%wopt.eq."1d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select
      fname = out_fname

   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr,mo,da,hr,mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(dir_name)//'/'
      
      write(unit=cdate, fmt='(i4.4,i2.2)') yr, mo
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      select case (wformat_temp)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
!         if(LIS_rc%wopt.eq."1d tilespace") then 
!            out_fname = trim(out_fname)//'.ts4r'
!         elseif(LIS_rc%wopt.eq."2d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         elseif(LIS_rc%wopt.eq."1d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select
      fname = out_fname
      
   elseif(LIS_rc%wstyle.eq."2 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr,mo,da,hr,mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(dir_name)//'/'
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      select case (wformat_temp)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
!         if(LIS_rc%wopt.eq."1d tilespace") then 
!            out_fname = trim(out_fname)//'.ts4r'
!         elseif(LIS_rc%wopt.eq."2d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         elseif(LIS_rc%wopt.eq."1d gridspace") then 
!            out_fname = trim(out_fname)//'.gs4r'
!         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select
      fname = out_fname
      
   elseif(LIS_rc%wstyle.eq."WMO convention" .or. &
        LIS_rc%wstyle.eq."557WW streamflow convention" .or. &
        LIS_rc%wstyle.eq."557WW medium range forecast convention") then 

      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr,mo,da,hr,mn
      
      dname = LIS_rc%odir
         
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'/LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      
      select case (wformat_temp)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
!         if(LIS_rc%wopt.eq."1d tilespace") then 
!            out_fname = trim(dname)//'.DAT'
!         elseif(LIS_rc%wopt.eq."2d gridspace") then 
!            out_fname = trim(dname)//'.DAT'
!         elseif(LIS_rc%wopt.eq."1d gridspace") then 
!            out_fname = trim(out_fname)//'.DAT'
!         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select

      fname = out_fname
   endif
 end subroutine Create_restart_filename_withtime

!BOP
!
! !ROUTINE: LIS_create_stats_filename
! \label{LIS_create_stats_filename}
! 
! !INTERFACE:
subroutine LIS_create_stats_filename(n, fname, mname)
! !USES:
   use LIS_coreMod,  only : LIS_rc
  
   implicit none 
! !ARGUMENTS:
   integer, intent(in)           :: n
   character(len=*), intent(out) :: fname
   character(len=*), intent(in)  :: mname
!
! !DESCRIPTION:  
!  Create the file name for the stats files.  The convention used
!  in LIS creates a restart filename in the following format. 
!
!  \begin{verbatim}
!  <output directory>/EXP<expno>/<model name>stats.dat
!  \end{verbatim}
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest
!   \item [fname]
!     the created file name. 
!   \item [mname]
!    string describing the name of the model
!   \item [style]
!    output directory style
!  \end{description}
! 
!
!EOP
   character(len=LIS_CONST_PATH_LEN) :: out_fname
   character*100      :: cdate

   write(unit=cdate, fmt='(a2,i2.2)') '.d',n
   
   out_fname = trim(LIS_rc%odir)//'/' &
        //trim(mname)//trim(cdate)//'.stats'
   
   fname = out_fname

 end subroutine LIS_create_stats_filename


!BOP
!
! !ROUTINE: LIS_create_innov_filename
! \label{LIS_create_innov_filename}
! 
! !INTERFACE:
subroutine LIS_create_innov_filename(n, k, fname, mname)
   use LIS_coreMod,  only : LIS_rc
   use LIS_logMod,   only : LIS_log_msg, LIS_endrun

   implicit none 
! !ARGUMENTS:
   integer, intent(in)            :: n
   integer, intent(in)            :: k
   character(len=*), intent(out)  :: fname
   character(len=*), intent(in)   :: mname
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
   character(len=10)       :: cda
   character(len=14)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=LIS_CONST_PATH_LEN)       :: dname
   character(len=LIS_CONST_PATH_LEN), save :: out_fname
   integer                  :: i, c

   if(LIS_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4)') LIS_rc%yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_innov'
      
      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_innov'
      
      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   elseif(LIS_rc%wstyle.eq."2 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      out_fname = trim(dname)//'LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_innov'

      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   elseif(LIS_rc%wstyle.eq."WMO convention" .or. &
        LIS_rc%wstyle.eq."557WW streamflow convention" .or. &
        LIS_rc%wstyle.eq."557WW medium range forecast convention") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//'/'
      
      out_fname = trim(dname)//'LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_innov'
      
      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   endif
   fname = out_fname

        
 end subroutine LIS_create_innov_filename

!BOP
!
! !ROUTINE: LIS_create_incr_filename
! \label{LIS_create_incr_filename}
! 
! !INTERFACE:
subroutine LIS_create_incr_filename(n, k, fname, mname)
   use LIS_coreMod,  only : LIS_rc
   use LIS_logMod,   only : LIS_log_msg, LIS_endrun

   implicit none 
! !ARGUMENTS:
   integer, intent(in)            :: n
   integer, intent(in)            :: k
   character(len=*), intent(out)  :: fname
   character(len=*), intent(in)   :: mname
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
   character(len=10)       :: cda
   character(len=14)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=LIS_CONST_PATH_LEN)       :: dname
   character(len=LIS_CONST_PATH_LEN), save :: out_fname
   integer                  :: i, c

   if(LIS_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4)') LIS_rc%yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_incr'
      
      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_incr'
      
      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   elseif(LIS_rc%wstyle.eq."2 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      out_fname = trim(dname)//'LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_incr'

      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   elseif(LIS_rc%wstyle.eq."WMO convention" .or. &
        LIS_rc%wstyle.eq."557WW streamflow convention" .or. &
        LIS_rc%wstyle.eq."557WW medium range forecast convention") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//'/'
      
      out_fname = trim(dname)//'LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_incr'
      
      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   endif
   fname = out_fname

        
 end subroutine LIS_create_incr_filename

!BOP
!
! !ROUTINE: LIS_create_daspread_filename
! \label{LIS_create_daspread_filename}
! 
! !INTERFACE:
subroutine LIS_create_daspread_filename(n, k, fname, mname)
   use LIS_coreMod,  only : LIS_rc
   use LIS_logMod,   only : LIS_log_msg, LIS_endrun

   implicit none 
! !ARGUMENTS:
   integer, intent(in) :: n
   integer, intent(in) :: k
   character(len=*), intent(out)          :: fname
   character(len=*), intent(in)  :: mname
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
   character(len=10)       :: cda
   character(len=14)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=LIS_CONST_PATH_LEN)       :: dname
   character(len=LIS_CONST_PATH_LEN), save :: out_fname
   integer                  :: i, c

   if(LIS_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4)') LIS_rc%yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_spread'
      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_spread'
      
      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   elseif(LIS_rc%wstyle.eq."2 level hierarchy".or.&
        LIS_rc%wstyle.eq."WMO convention" .or. &
        LIS_rc%wstyle.eq."557WW streamflow convention" .or. &
        LIS_rc%wstyle.eq."557WW medium range forecast convention") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      out_fname = trim(dname)//'LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_spread'
      
      write(unit=cda, fmt='(a2,i2.2)') '.a',k      
      out_fname = trim(out_fname)//trim(cda)

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   endif
   fname = out_fname

        
 end subroutine LIS_create_daspread_filename


!BOP
! !ROUTINE: LIS_create_obs_filename
! \label{LIS_create_obs_filename}
! 
! !INTERFACE: 
subroutine LIS_create_obs_filename(n, fname, mname)
! !USES:
   use LIS_coreMod,   only : LIS_rc
  
   implicit none 
! !ARGUMENTS:
   integer, intent(in)           :: n
   character(len=*), intent(out) :: fname
   character(len=*), intent(in)  :: mname
!
! !DESCRIPTION:  
!  Create the file name for the observation files.  The convention used
!  in LIS creates a restart filename in the following format. 
!
!  \begin{verbatim}
!  <output directory>/EXP<expno>/<model name><timestamp>.obs
!  \end{verbatim}
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest
!   \item [fname]
!     the created file name. 
!   \item [mname]
!    string describing the name of the model 
!  \end{description}
! 
!
!EOP
   character(len=LIS_CONST_PATH_LEN) :: out_fname
   character*100      :: cdate
   character(len=12)       :: cdate1

   write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') LIS_rc%yr, LIS_rc%mo, &
        LIS_rc%da, LIS_rc%hr,LIS_rc%mn

   write(unit=cdate, fmt='(a2,i2.2)') '.d',n

   out_fname = trim(LIS_rc%odir)//'/'//&
               trim(mname)//'/'//trim(mname)//trim(cdate)//'.'//&
               trim(cdate1)//'.obs'

   fname = out_fname
        
 end subroutine LIS_create_obs_filename

  

!BOP
! !ROUTINE: putget_int
!
! !INTERFACE:
! !Private name: call using putget()
    subroutine putget_int ( buffer, iofunc, file_name, &
         routine_name, imax, jmax )
! !USES: 
      use LIS_logMod,    only : LIS_abort, LIS_endrun, & 
           LIS_getNextUnitNumber, LIS_releaseUnitNumber

      implicit none
! !ARGUMENTS: 
      integer,           intent(in)     :: imax
      integer,           intent(in)     :: jmax
      integer,           intent(inout)  :: buffer (imax , jmax)
      character*1,       intent(in)     :: iofunc
      character(len=*),  intent(in)     :: file_name
      character(len=*),  intent(in)     :: routine_name
! !DESCRIPTION: 
! 
!  put or get any size integer array routine
!
!  The arguments are: 
!  \begin{description}
!   \item [buffer]
!      integer buffer array
!   \item [iofunc]
!      I/O function (r=read, w=write)
!   \item [file\_name]
!      input file name
!   \item [routine\_name]
!      name of the calling routine
!   \item [imax]
!      grid dimension of the buffer array - East/West direction
!   \item [jmax]
!      grid dimension of the buffer array - North/South direction
!  \end{description}
!
!  The methods invoked are: 
!  \begin{description}
!   \item[LIS\_abort](\ref{LIS_abort}) \newline
!     abort the program in case of a fatal error
!  \end{description}
!EOP

      character*9                   :: cstat
      character*100                 :: message     (20)
      integer                       :: rec_length
      integer                       :: istat
      integer                       :: istat1
      integer                       :: ftn

!     ------------------------------------------------------------------
!     executable code starts here ... open file, abort on error
!     ------------------------------------------------------------------

      rec_length = imax * jmax * 4

!YDT 9/27/07, create dir for write in case it does not exit
      if( iofunc .eq. 'w' ) &
        call system('mkdir -p `dirname '//trim(file_name)//'`')
      

      ftn = LIS_getNextUnitNumber()
      open( ftn, file=trim(file_name), form='unformatted', &
           access='direct', recl=rec_length, iostat=istat )
      if( istat .eq. 0 ) then 

!     ------------------------------------------------------------------
!     read from file, abort on error
!     ------------------------------------------------------------------

         if( iofunc .eq. 'r' )then
            read( ftn, rec=1, iostat=istat ) buffer
            if( istat .ne. 0 ) then 
               call LIS_releaseUnitNumber(ftn)
               write(cstat,'(i9)',iostat=istat1) istat
               message(1) = 'program:  LIS'
               message(2) = '  routine:  '//trim(routine_name)
               message(3) = '  error reading file '//trim(file_name)
               if( istat1 .eq. 0 )then
                  message(4) = '  status = '//trim(cstat)
               endif
               call LIS_abort( message )
               call LIS_endrun
            endif
!     ------------------------------------------------------------------
!     write to file, abort on error
!     ------------------------------------------------------------------

         elseif( iofunc .eq. 'w' )then
            write( ftn, rec=1, iostat=istat ) buffer
            if( istat .ne. 0 ) then 
               call LIS_releaseUnitNumber(ftn)
               write(cstat,'(i9)',iostat=istat1) istat
               message(1) = 'program:  LIS'
               message(2) = '  routine:  '//trim(routine_name)
               message(3) = '  error writing file '//trim(file_name)
               if( istat1 .eq. 0 )then
                  message(4) = '  status = '//trim(cstat)
               endif
               call LIS_abort( message )
               call LIS_endrun
            endif

!     ------------------------------------------------------------------
!     else abort due to invalid iofunc value
!     ------------------------------------------------------------------

         else
            call LIS_releaseUnitNumber(ftn)
            message(1) = 'program:  LIS'
            message(2) = '  routine:  '//trim(routine_name)
            message(3) = '  invalid iofunc value ='//trim(iofunc)
            call LIS_abort( message )
            call LIS_endrun
         endif

         call LIS_releaseUnitNumber(ftn)
      
         return

!     ------------------------------------------------------------------
!     error handling
!     ------------------------------------------------------------------
      else
         call LIS_releaseUnitNumber(ftn)
         write(cstat,'(i9)',iostat=istat1) istat
         message(1) = 'program:  LIS'
         message(2) = '  routine:  '//trim(routine_name)
         message(3) = '  error opening file '//trim(file_name)
         if( istat1 .eq. 0 )then
            message(4) = '  status = '//trim(cstat)
         endif
         call LIS_abort( message )
         stop
      endif
    end subroutine putget_int

!BOP
! 
! !ROUTINE: putget_real
! \label{putget_real}
!     
! !INTERFACE:
! !Private name: call using putget()
    subroutine putget_real ( buffer, iofunc, file_name, &
         routine_name, imax, jmax )
! !USES: 
      use LIS_logMod,    only : LIS_abort, & 
           LIS_getNextUnitNumber, LIS_releaseUnitNumber
      implicit none
      
      integer,       intent(in)     :: imax
      integer,       intent(in)     :: jmax
      real,         intent(inout)   :: buffer      (imax , jmax)
      character(len=*), intent(in)  :: file_name
      character*1,   intent(in)     :: iofunc
      character(len=*),  intent(in) :: routine_name

! !DESCRIPTION: 
! 
!  put or get any size real array routine
!
!  The arguments are: 
!  \begin{description}
!   \item [buffer]
!      integer buffer array
!   \item [iofunc]
!      I/O function (r=read, w=write)
!   \item [file\_name]
!      input file name
!   \item [routine\_name]
!      name of the calling routine
!   \item [imax]
!      grid dimension of the buffer array - East/West direction
!   \item [jmax]
!      grid dimension of the buffer array - North/South direction
!  \end{description}
!
!
!  The methods invoked are: 
!  \begin{description}
!   \item[LIS\_abort](\ref{LIS_abort}) \newline
!     abort the program in case of a fatal error
!  \end{description}
!EOP      

      character*100                 :: message     (20)
      integer                       :: rec_length
      character*9                   :: cstat
      integer                       :: istat
      integer                       :: istat1
      integer                       :: ftn

!     ------------------------------------------------------------------
!     executable code starts here ... open file, abort on error
!     ------------------------------------------------------------------
!YDT 9/27/07, create dir for write in case it does not exit
      if( iofunc .eq. 'w' ) &
        call system('mkdir -p `dirname '//trim(file_name)//'`')

      rec_length = imax * jmax * 4

      ftn = LIS_getNextUnitNumber()
      open( ftn, file=trim(file_name), form='unformatted', &
           access='direct', recl=rec_length, iostat=istat )
      if( istat .eq. 0 ) then 
         
!     ------------------------------------------------------------------
!     read from file, abort on error
!     ------------------------------------------------------------------
  
         if( iofunc .eq. 'r' )then
            read( ftn, rec=1, iostat=istat ) buffer
            if( istat .ne. 0 ) then 
               call LIS_releaseUnitNumber(ftn)
               write(cstat,'(i9)',iostat=istat1) istat
               message(1) = 'program:  LIS'
               message(2) = '  routine:  '//trim(routine_name)
               message(3) = '  error reading file '//trim(file_name)
               if( istat1 .eq. 0 )then
                  message(4) = '  status = '//trim(cstat)
               endif
               call LIS_abort( message )
            endif
     
!     ------------------------------------------------------------------
!     write to file, abort on error
!     ------------------------------------------------------------------
     
         elseif( iofunc .eq. 'w' )then
            write( ftn, rec=1, iostat=istat ) buffer
            if( istat .ne. 0 ) then 
               call LIS_releaseUnitNumber(ftn)
               write(cstat,'(i9)',iostat=istat1) istat
               message(1) = 'program:  LIS'
               message(2) = '  routine:  '//trim(routine_name)
               message(3) = '  error writing file '//trim(file_name)
               if( istat1 .eq. 0 )then
                  message(4) = '  status = '//trim(cstat)
               endif
               call LIS_abort( message )
            endif
     
!     ------------------------------------------------------------------
!     else abort due to invalid iofunc value
!     ------------------------------------------------------------------
            
         else
            call LIS_releaseUnitNumber(ftn)
            message(1) = 'program:  LIS'
            message(2) = '  routine:  '//trim(routine_name)
            message(3) = '  invalid iofunc value ='//trim(iofunc)
            call LIS_abort( message )
         endif
      else
         
         call LIS_releaseUnitNumber(ftn)
         write(cstat,'(i9)',iostat=istat1) istat
         message(1) = 'program:  LIS'
         message(2) = '  routine:  '//trim(routine_name)
         message(3) = '  error opening file '//trim(file_name)
         if( istat1 .eq. 0 )then
            message(4) = '  status = '//trim(cstat)
         endif
         call LIS_abort( message )
      endif
      
      call LIS_releaseUnitNumber(ftn)
      
      return
      
    end subroutine putget_real


!BOP
!
! !ROUTINE: LIS_create_gain_filename
! \label{LIS_create_gain_filename}
! 
! !INTERFACE:
subroutine LIS_create_gain_filename(n, fname, mname)
   use LIS_coreMod,  only : LIS_rc
   use LIS_logMod,   only : LIS_log_msg, LIS_endrun

   implicit none 
! !ARGUMENTS:
   integer, intent(in) :: n
   character(len=*), intent(out)          :: fname
   character(len=*), intent(in)  :: mname
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
   character(len=LIS_CONST_PATH_LEN)       :: dname
   character(len=LIS_CONST_PATH_LEN), save :: out_fname
   integer                  :: i, c

   if(LIS_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4)') LIS_rc%yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_gain'
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   elseif(LIS_rc%wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') LIS_rc%yr, LIS_rc%mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_gain'
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   elseif(LIS_rc%wstyle.eq."2 level hierarchy".or.&
        LIS_rc%wstyle.eq."WMO convention" .or. &
        LIS_rc%wstyle.eq."557WW streamflow convention" .or. &
        LIS_rc%wstyle.eq."557WW medium range forecast convention") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn
      
      dname = trim(LIS_rc%odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      out_fname = trim(dname)//'LIS_DA_'//trim(mname)//'_'//trim(cdate1)//&
           '_gain'
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      out_fname = trim(out_fname)//'.nc'
   endif
   fname = out_fname

        
 end subroutine LIS_create_gain_filename

!BOP
! !ROUTINE: read2Ddata
!  \label{read2Ddata}
! 
! !INTERFACE: 
 subroutine read2Ddata(n, ftn, gridDesc, array)
! !USES: 
   use LIS_coreMod
   use LIS_logMod

   implicit none
! !ARGUMENTS: 
   integer,  INTENT(IN)           :: n
   integer,  INTENT(IN)           :: ftn
   real,     INTENT(INOUT)        :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
   real,     INTENT(IN)           :: gridDesc(6)
!
! !DESCRIPTION: 
!   This routine retrieves a 2-d data from a binary direct access file. 
!EOP  
   integer  :: c,r,line,istat
   real     :: data_in(LIS_rc%gnc(n),LIS_rc%gnr(n))
   
   do r=1,LIS_rc%gnr(n)
      do c=1,LIS_rc%gnc(n)
         line = c+(r-1)*LIS_rc%gnc(n)
         read(ftn,rec=line,iostat=istat) data_in(c,r)
      enddo
   enddo

   array(:,:) = data_in(&
        LIS_ews_halo_ind(n,LIS_localPet+1):&         
        LIS_ewe_halo_ind(n,LIS_localPet+1), &
        LIS_nss_halo_ind(n,LIS_localPet+1): &
        LIS_nse_halo_ind(n,LIS_localPet+1))

#if 0 
   if(LIS_rc%param_proj.eq."gaussian") then ! gaussian 
      line1 = gaussian_find_row(LIS_rc%gridDesc(n,4))  -   &
           gaussian_find_row(LIS_rc%gridDesc(n,44)) + 1
      line2 = gaussian_find_col(LIS_rc%gridDesc(n,5))  -   &
           gaussian_find_col(LIS_rc%gridDesc(n,45)) + 1
      
      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            glnc = line2+c-1
            glnr = line1+r-1
            line = (glnr-1)*nint(LIS_rc%gridDesc(n,42))+glnc
            read(ftn,rec=line,iostat=istat) array(c,r)
            if( istat .ne. 0 ) then
               message(1) = 'program:  LIS'
               message(2) = '  routine:  read2DData'
               message(3) = '  iostat != 0'
               call LIS_abort( message )
               call LIS_endrun
            endif
         enddo
      enddo
   elseif(LIS_rc%param_proj.eq."latlon") then 

      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            call ij_to_latlon(LIS_domain(n)%lisproj,float(c),float(r),&
                 rlat(c,r),rlon(c,r))
         enddo
      enddo
      
      nc_dom = nint((gridDesc(4)-gridDesc(2))/(gridDesc(5)))+1
      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            line1 = nint((rlat(c,r)-gridDesc(1))/gridDesc(6))+1
            line2 = nint((rlon(c,r)-gridDesc(2))/gridDesc(5))+1
            line = (line1-1)*nc_dom + line2
            read(ftn,rec=line) array(c,r)
         enddo
      enddo
   elseif(LIS_rc%param_proj.eq."polar") then !ps

#if 0
      call map_set(PROJ_PS,LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,5), &
                   LIS_rc%gridDesc(n,8)*1000.0,                       &
                   LIS_rc%gridDesc(n,11),LIS_rc%gridDesc(n,10),0.0,   &
                   LIS_rc%lnc(n),LIS_rc%lnr(n),proj)

      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            call ij_to_latlon(proj,float(c),float(r),rlat(c,r),rlon(c,r))
         enddo
      enddo

      call map_set(PROJ_PS,gridDesc(1),gridDesc(2),  &
                   gridDesc(6)*1000.0,               &
                   gridDesc(4),gridDesc(3),0.0,      &
                   int(gridDesc(7)),int(gridDesc(8)),&
                   proj)
     
      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            call latlon_to_ij(proj,rlat(c,r),rlon(c,r),ctmp,rtmp)

            line1 = nint(rtmp)
            line2 = nint(ctmp)
            line = (line1-1)*gridDesc(7)+line2
            read(ftn,rec=line) array(c,r) 
         enddo
      enddo
#endif

   elseif(LIS_rc%param_proj.eq."UTM") then !utm
!rlat/rlon used here to store northing and easting
      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            rlat(c,r) = LIS_rc%gridDesc(n,4)+(r-1)*LIS_rc%gridDesc(n,9)
            rlon(c,r) = LIS_rc%gridDesc(n,5)+(c-1)*LIS_rc%gridDesc(n,9)
         enddo
      enddo
      
      nc_dom = gridDesc(4)

      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            line1 = nint((rlat(c,r)-gridDesc(2))/gridDesc(6))+1
            line2 = nint((rlon(c,r)-gridDesc(3))/gridDesc(6))+1
            line = (line1-1)*nc_dom +line2
           read(ftn,rec=line) array(c,r) 
         enddo
      enddo
   else
      print*, 'This parameter projection is not supported...'
      print*, 'Program stopping ....'
      stop
   endif
#endif
 end subroutine read2Ddata


!BOP
! !ROUTINE: LIS_readDomainConfigSpecs
! \label{LIS_readDomainConfigSpecs}
! 
! !INTERFACE: 
 subroutine LIS_readDomainConfigSpecs(segment_name, domain_info)
! !USES: 
   use LIS_coreMod,  only : LIS_rc, LIS_config
   use LIS_logMod,   only : LIS_verify
! !ARGUMENTS: 
   character(len=*),  intent(in)    :: segment_name
   real,              intent(inout) :: domain_info(:,:)
!   real,              intent(inout) :: domain_info(LIS_rc%nnest,6)
! 
! !DESCRIPTION: 
!   This subroutine reads the relevant section of the lis.config file 
!   that defines the domain specifications for a particular surface 
!   parameter dataset, based on the map projection used to define 
!   the surface parameter dataset. 
!EOP
   integer                         :: i, rc

   if(LIS_rc%param_proj.eq."latlon") then 
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" lower left lat:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,1),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' lower left lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" lower left lon:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,2),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' lower left lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" upper right lat:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,3),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' upper right lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" upper right lon:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,4),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' upper right lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" resolution (dx):",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,5),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' resolution (dx):')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" resolution (dy):",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,6),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' resolution (dy):')
      enddo
      
   elseif(LIS_rc%param_proj.eq."gaussian") then !gaussian
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" first grid point lat:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,1),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' first grid point lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" first grid point lon:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,2),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' first grid point lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" last grid point lat:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,3),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' last grid point lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" last grid point lon:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,4),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' last grid point lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" resolution dlon:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,5),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' resolution dlon:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" number of lat circles:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,6),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' number of lat circles:')
      enddo
   elseif(LIS_rc%param_proj.eq."polar") then !polar
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" lower left lat:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,1),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' lower left lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" lower left lon:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,2),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' lower left lon:')
      enddo
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" true lat:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,3),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' true lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" standard lon:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,4),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' standard lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" orientation:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,5),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' orientation:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" resolution:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,6),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' resolution:')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" x-dimension size:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,7),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' x-dimension size:')
      enddo

      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" y-dimension size:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,8),rc=rc)
         call LIS_verify(rc,'please specify '//trim(segment_name)//' y-dimension size:')
      enddo
   elseif(LIS_rc%param_proj.eq."UTM") then !UTM
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" UTM zone:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,1),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' UTM zone: not defined')
      enddo
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" northing of SW corner:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,2),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' northing of SW corner: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" easting of SW corner:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,3),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' easting of SW corner: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" x-dimension size:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,4),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' x-dimension size: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" y-dimension size:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,5),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' y-dimension size: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" resolution:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,6),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' resolution: not defined')
      enddo

   elseif(LIS_rc%param_proj.eq."hrap") then !HRAP
      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" lower left lat:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,1),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' lower left lat: not defined')
      enddo

      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" lower left lon:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,2),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' lower left lon: not defined')
      enddo

      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" true lat:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,3),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' true lat: not defined')
      enddo

      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" standard lon:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,4),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' standard lon: not defined')
      enddo

      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" orientation:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,5),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' orientation: not defined')
      enddo

      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" resolution:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,6),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' resolution: not defined')
      enddo

      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" x-dimension size:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,7),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' x-dimension size: not defined')
      enddo

      call ESMF_ConfigFindLabel(LIS_config,trim(segment_name)//" y-dimension size:",rc=rc)
      do i=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,domain_info(i,8),rc=rc)
         call LIS_verify(rc,trim(segment_name)//' y-dimension size: not defined')
      enddo

   else
      print*, 'Reading configuration settings for this projection ',trim(LIS_rc%param_proj)
      print*, 'is not supported'
      print*, 'Program stopping ....'
      stop
   endif
 end subroutine LIS_readDomainConfigSpecs


!BOP
! !ROUTINE: LIS_checkDomainExtents
! \label{LIS_checkDomainExtents}
! 
! !INTERFACE: 
 subroutine LIS_checkDomainExtents(n, data_gridDesc)
! !USES: 
   use LIS_coreMod,  only : LIS_rc, LIS_domain
   use LIS_logMod,   only : LIS_logunit, LIS_endrun
   use map_utils,    only : ij_to_latlon
! !ARGUMENTS: 
   integer              :: n 
   real                 :: data_gridDesc(6)
! 
! !DESCRIPTION: 
!   This subroutine checks to see if the domain extents of the data is
!   within the LIS running domain. 
!EOP
   real                            :: max_lat, max_lon
   real                            :: min_lat, min_lon
   real                            :: rlat, rlon
   integer                         :: c,r
   max_lat = -10000
   max_lon = -10000
   min_lat = 10000 
   min_lon = 10000
   
   do r=1,LIS_rc%lnr(n)
      do c=1,LIS_rc%lnc(n)
         call ij_to_latlon(LIS_domain(n)%lisproj,float(c),float(r),&
              rlat,rlon)
         
         if(rlat.gt.max_lat) max_lat = rlat
         if(rlon.gt.max_lon) max_lon = rlon
         if(rlat.lt.min_lat) min_lat = rlat
         if(rlon.lt.min_lon) min_lon = rlon
      enddo
   enddo
   
   if((min_lat.lt.data_gridDesc(1)).or.&
        (min_lon.lt.data_gridDesc(2)).or.&
        (max_lat.gt.data_gridDesc(3)).or.&
        (max_lon.gt.data_gridDesc(4))) then 
      
      write(LIS_logunit,*) '[ERR]  The STATSGO data is specified only for the CONUS..'
      write(LIS_logunit,*) '[ERR] The running domain is outside the STATSGO boundaries..'
      write(LIS_logunit,*) '[ERR] Stopping program....'
      call LIS_endrun
   endif
   
 end subroutine LIS_checkDomainExtents


!BOP
! 
! !ROUTINE: readparam_real_2d
! \label{readparam_real_2d}
!
! !INTERFACE:
  subroutine readparam_real_2d(n,pname,array)
! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,        only : LIS_rc, LIS_localPet,&
       LIS_domain, LIS_ews_ind, LIS_ewe_ind,&
       LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
       LIS_nss_halo_ind, LIS_nse_halo_ind
  use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_endrun, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  character(len=*)       :: pname
  real,    intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
! 
! !DESCRIPTION:
!  Reads a parameter field from the input LIS parameter file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[pname]
!    name of the parameter field
!   \item[array]
!    retrieved parameter value
!   \end{description}
!
!EOP      

  integer :: ios1
  integer :: ios,nid,paramid,ncId, nrId
  integer :: nc,nr,c,r
  logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=trim(LIS_rc%paramfile(n)), exist=file_exists)
  if(file_exists) then 

     ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in readparam_real_2d')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in readparam_real_2d')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in readparam_real_2d')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in readparam_real_2d')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in readparam_real_2d')

     ios = nf90_inq_varid(nid,trim(pname),paramid)
     call LIS_verify(ios,trim(pname)//' field not found in the LIS param file')

     ios = nf90_get_var(nid,paramid,array,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))            
     call LIS_verify(ios,'Error in nf90_get_var in readparam_real_2d')
     
     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in readparam_real_2d')

  else
     write(LIS_logunit,*) '[ERR] '//trim(pname)//' map: ',&
          trim(LIS_rc%paramfile(n)), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif

#endif

end subroutine readparam_real_2d


!BOP
! 
! !ROUTINE: readparam_real_2d_rc
! \label{readparam_real_2d_rc}
!
! !INTERFACE:
  subroutine readparam_real_2d_rc(n,pname,array,rc)
! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,        only : LIS_rc, LIS_localPet,&
       LIS_domain, LIS_ews_ind, LIS_ewe_ind,&
       LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
       LIS_nss_halo_ind, LIS_nse_halo_ind
  use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_endrun, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  character(len=*)       :: pname
  real,    intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                :: rc
! 
! !DESCRIPTION:
!  Reads a parameter field from the input LIS parameter file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[pname]
!    name of the parameter field
!   \item[array]
!    retrieved parameter value
!   \end{description}
!
!EOP      

  integer :: ios1
  integer :: ios,nid,paramid,ncId, nrId
  integer :: nc,nr,c,r
  logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  rc = 1
  inquire(file=trim(LIS_rc%paramfile(n)), exist=file_exists)
  if(file_exists) then 
     ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in readparam_real_2d')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in readparam_real_2d')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in readparam_real_2d')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in readparam_real_2d')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in readparam_real_2d')

     ios = nf90_inq_varid(nid,trim(pname),paramid)
     if(ios.ne.0) then 
        rc = 1
     else
        ios = nf90_get_var(nid,paramid,array,&
             start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
             LIS_nss_halo_ind(n,LIS_localPet+1)/),&
             count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
        call LIS_verify(ios,'Error in nf90_get_var in readparam_real_2d')
        
        ios = nf90_close(nid)
        call LIS_verify(ios,'Error in nf90_close in readparam_real_2d')
        
        rc = 0 
     endif
  else
     write(LIS_logunit,*) '[ERR] '//trim(pname)//' map: ',&
          trim(LIS_rc%paramfile(n)), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
     rc = 1
  endif

#endif

end subroutine readparam_real_2d_rc

!BOP
! 
! !ROUTINE: readgparam_real_2d
! \label{readgparam_real_2d}
!
! !INTERFACE:
  subroutine readgparam_real_2d(n,pname,array)
! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,        only : LIS_rc
  use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_endrun, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  character(len=*)       :: pname
  real,    intent(inout) :: array(LIS_rc%gnc(n),LIS_rc%gnr(n))
! 
! !DESCRIPTION:
!  Reads a parameter field from the input LIS parameter file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[pname]
!    name of the parameter field
!   \item[array]
!    retrieved parameter value
!   \end{description}
!
!EOP      

  integer :: ios1
  integer :: ios,nid,paramid,ncId, nrId
  integer :: nc,nr,c,r
  real    :: param(LIS_rc%gnc(n),LIS_rc%gnr(n))
  logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=trim(LIS_rc%paramfile(n)), exist=file_exists)
  if(file_exists) then 

     ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in readparam_real_2d')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in readparam_real_2d')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in readparam_real_2d')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in readparam_real_2d')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in readparam_real_2d')

     ios = nf90_inq_varid(nid,trim(pname),paramid)
     call LIS_verify(ios,trim(pname)//' field not found in the LIS param file')

     ios = nf90_get_var(nid,paramid,param)
     call LIS_verify(ios,'Error in nf90_get_var in readparam_real_2d')
     
     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in readparam_real_2d')

     array(:,:) = param(:,:)

  else
     write(LIS_logunit,*) '[ERR] '//trim(pname)//' map: ',&
          trim(LIS_rc%paramfile(n)), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif

#endif

end subroutine readgparam_real_2d


!BOP
! 
! !ROUTINE: readgparam_real_2d_rc
! \label{readgparam_real_2d_rc}
!
! !INTERFACE:
  subroutine readgparam_real_2d_rc(n,pname,array,rc)
! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,        only : LIS_rc
  use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_endrun, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  character(len=*)       :: pname
  real,    intent(inout) :: array(LIS_rc%gnc(n),LIS_rc%gnr(n))
  integer                :: rc
! 
! !DESCRIPTION:
!  Reads a parameter field from the input LIS parameter file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[pname]
!    name of the parameter field
!   \item[array]
!    retrieved parameter value
!   \item[rc]
!    file read return code
!   \end{description}
!
!EOP      

  integer :: ios1
  integer :: ios,nid,paramid,ncId, nrId
  integer :: nc,nr,c,r
  real    :: param(LIS_rc%gnc(n),LIS_rc%gnr(n))
  logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  rc = 1
  inquire(file=trim(LIS_rc%paramfile(n)), exist=file_exists)
  if(file_exists) then 

     ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in readparam_real_2d_rc')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in readparam_real_2d_rc')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in readparam_real_2d_rc')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in readparam_real_2d_rc')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in readparam_real_2d_rc')

     ios = nf90_inq_varid(nid,trim(pname),paramid)
     if(ios.ne.0) then 
        rc = 1
     else
        ios = nf90_get_var(nid,paramid,param)
        call LIS_verify(ios,'Error in nf90_get_var in readparam_real_2d_rc')
        
        ios = nf90_close(nid)
        call LIS_verify(ios,'Error in nf90_close in readparam_real_2d_rc')
        
        array(:,:) = param(:,:)
        rc = 0
     endif
  else
     write(LIS_logunit,*) '[ERR] '//trim(pname)//' map: ',&
          trim(LIS_rc%paramfile(n)), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
     rc = 1
  endif

#endif

end subroutine readgparam_real_2d_rc



!BOP
! 
! !ROUTINE: readparam_int_2d
! \label{readparam_int_2d}
!
! !INTERFACE:
  subroutine readparam_int_2d(n,pname,array)
! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,        only : LIS_rc, LIS_localPet,&
       LIS_domain, LIS_ews_ind, LIS_ewe_ind,&
       LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
       LIS_nss_halo_ind, LIS_nse_halo_ind
  use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_endrun, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  character(len=*)       :: pname
  integer,    intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
! 
! !DESCRIPTION:
!  Reads a parameter field from the input LIS parameter file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[pname]
!    name of the parameter field
!   \item[array]
!    retrieved parameter value
!   \end{description}
!
!EOP      

  integer :: ios1
  integer :: ios,nid,paramid,ncId, nrId
  integer :: nc,nr,c,r
  logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=trim(LIS_rc%paramfile(n)), exist=file_exists)
  if(file_exists) then 

     write(LIS_logunit,*)'[INFO] Reading '//trim(pname)//' map '
     ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in readparam_int_2d')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in readparam_int_2d')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in readparam_int_2d')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in readparam_int_2d')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in readparam_int_2d')

     ios = nf90_inq_varid(nid,trim(pname),paramid)
     call LIS_verify(ios,trim(pname)//' field not found in the LIS param file')

     ios = nf90_get_var(nid,paramid,array,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))          
     call LIS_verify(ios,'Error in nf90_get_var in readparam_int_2d')
     
     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in readparam_int_2d')

  else
     write(LIS_logunit,*) '[ERR] '//trim(pname)//' map: ',&
          trim(LIS_rc%paramfile(n)), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif

#endif

 end subroutine readparam_int_2d


end module LIS_fileIOMod
