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
module LDT_fileIOMod
!BOP
!
! !MODULE: LDT_fileIOMod
! 
! !DESCRIPTION: 
!   This module contains a number of routines useful for various file I/O
!   operations in LDT. The module provides routines to create output directories, 
!   filenames, that are be used in the model output routines. 
!   
! !REVISION HISTORY: 
!  02 Oct 2008  Sujay Kumar;  Initial Specification
!  18 Feb 2015  KR Arsenault;  Added additional parameter reading capability
! 
! !USES: 
  use LDT_coreMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none 
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_create_output_directory
  public :: LDT_create_output_filename   ! create an output filename  
  public :: LDT_create_stats_filename    ! create a stats filename
  public :: LDT_create_daobs_filename
  public :: LDT_create_restart_filename  ! create a restart filename
  public :: LDT_transform_paramgrid      ! transform parameter file grid to LIS target grid
  public :: readLISdata                  ! generic method for reading surface parameters
  public :: LDT_readDomainConfigSpecs    ! generic method to read configurable options
                                         !  for a particular surface dataset 
  public :: LDT_checkDomainExtents       ! checks if the data domain extents are 
                                         !  contained in the LDT running domain
  public :: LDT_read_param               ! read LDT-processed parameters
!EOP

!BOP
! 
!  !ROUTINE: readLISdata
! \label{readLISdata}
! 
! !INTERFACE: 
  interface readLISdata
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure readparam_real_2d
     module procedure readparam_int_2d
! 
! !DESCRIPTION: 
!  Generic routine to read parameter data.  
!  
!EOP
  end interface
  
!BOP
! 
!  !ROUTINE: LDT_read_param
! \label{LDT_read_param}
! 
! !INTERFACE: 
  interface LDT_read_param
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure readldtparam_real_2d
!     module procedure readldtparam_int_2d
! 
! !DESCRIPTION: 
!  Routine to read parameter data from a netcdf file.  
!  
!EOP
  end interface

!BOP
! 
!  !ROUTINE: LDT_create_restart_filename
! \label{LDT_create_restart_filename}
! 
! !INTERFACE: 
  interface LDT_create_restart_filename
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure create_restart_filename_withtime
! 
! !DESCRIPTION: 
! This routine puts together a restart filename, either based on the
! LIS current time or based on the time that is provided. 
!  
!EOP
  end interface

contains

!BOP
!
! !ROUTINE: LDT_create_output_directory
! \label{LDT_create_output_directory}
!
! !INTERFACE:
subroutine LDT_create_output_directory(mname,dir_name)
! !USES:
   use LDT_coreMod, only : LDT_rc
   use LDT_logMod,  only : LDT_log_msg
   implicit none 
! !ARGUMENTS:
   character(len=*)  :: mname
   character(len=*), optional   :: dir_name
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
!   \item [dir\_name]
!     name of the directory to override the above format
!  \end{description}
!
!EOP
   character(len=4) :: cdate
   character(len=8) :: cdate1
   character(len=200) :: out_dname
   integer            :: try,ios
#if (!defined AIX )
   integer            :: system 
#endif
   
   ! EMK...Calls to 'system' fail when using SGI MPT as the MPI implementation
   ! on Pleiades. We replace with a C wrapper function that calls the 'mkdir' 
   ! standard POSIX function. This requires defining the C wrapper function,
   ! and specifying new variables to pass to said C function.
   integer, external :: LDT_create_subdirs 
   character(len=201) :: c_string  

   if(LDT_rc%wstyle.eq."4 level hierarchy") then 

      out_dname = trim(LDT_rc%odir)//'/'
      
      out_dname = trim(out_dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4)') LDT_rc%yr
      out_dname = trim(out_dname)//trim(cdate)//'/'
    
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da
      out_dname = trim(out_dname)//trim(cdate1)

      if ( present(dir_name) ) then
         dir_name = out_dname
      else
#if ( defined AIX )
         call system('mkdir -p '//trim(out_dname),ios)
#else
! EMK...Calls to 'system' fail when using SGI MPT as the MPI implementation
! on Pleiades. We replace with a C wrapper function that calls the 'mkdir' 
! standard POSIX function. 
!            ios = system('mkdir -p '//trim(out_dname))
            c_string = trim(out_dname)
            ios = LDT_create_subdirs(len_trim(c_string),trim(c_string))
#endif
      endif
   elseif(LDT_rc%wstyle.eq."3 level hierarchy") then 
      out_dname = trim(LDT_rc%odir)//'/'
      
      out_dname = trim(out_dname)//trim(mname)//'/'

      write(unit=cdate1, fmt='(i4.4, i2.2)') LDT_rc%yr, LDT_rc%mo
      out_dname = trim(out_dname)//trim(cdate1)

      if ( present(dir_name) ) then
         dir_name = out_dname
      else
#if ( defined AIX )
         call system('mkdir -p '//trim(out_dname),ios)
#else
! EMK...Calls to 'system' fail when using SGI MPT as the MPI implementation
! on Pleiades. We replace with a C wrapper function that calls the 'mkdir' 
! standard POSIX function. 
!         ios = system('mkdir -p '//trim(out_dname))
         c_string = trim(out_dname)
         ios = LDT_create_subdirs(len_trim(c_string),trim(c_string))
#endif
      endif      
   elseif(LDT_rc%wstyle.eq."2 level hierarchy") then 
      out_dname = trim(LDT_rc%odir)//'/'
      
      out_dname = trim(out_dname)//trim(mname)//'/'

      if ( present(dir_name) ) then
         dir_name = out_dname
      else
#if ( defined AIX )
         call system('mkdir -p '//trim(out_dname),ios)
#else
! EMK...Calls to 'system' fail when using SGI MPT as the MPI implementation
! on Pleiades. We replace with a C wrapper function that calls the 'mkdir' 
! standard POSIX function. 
!         ios = system('mkdir -p '//trim(out_dname))
         c_string = trim(out_dname)
         ios = LDT_create_subdirs(len_trim(c_string),trim(c_string))
      endif      
#endif
   elseif(LDT_rc%wstyle.eq."WMO convention") then 
      out_dname = trim(LDT_rc%odir)

      if ( present(dir_name) ) then
         dir_name = out_dname
      else
#if ( defined AIX )
         call system('mkdir -p '//trim(out_dname))
#else
! EMK...Calls to 'system' fail when using SGI MPT as the MPI implementation
! on Pleiades. We replace with a C wrapper function that calls the 'mkdir' 
! standard POSIX function. 
!         ios = system('mkdir -p '//trim(out_dname))
         c_string = trim(out_dname)
         ios = LDT_create_subdirs(len_trim(c_string),trim(c_string))
#endif
      endif      
   endif

 end subroutine LDT_create_output_directory

!BOP
!
! !ROUTINE: LDT_create_output_filename
! \label{LDT_create_output_filename}
!
! !INTERFACE:
subroutine LDT_create_output_filename(n, fname, model_name, odir, writeint)
! !USES:
   use LDT_coreMod,  only : LDT_rc
   use LDT_logMod,   only : LDT_log_msg, LDT_logunit, LDT_endrun

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
!   <output directory>/<model name>/LDT_HIST_<yyyymmddhhmnss>.<extension>
!  \end{verbatim}
!  3 level hierarchy
!  \begin{verbatim}
!   <output directory>/<model name>/<yyyymm>/LDT_HIST_<yyyymmddhhmnss>.<extension>
!  \end{verbatim}
!  4 level hierarchy
!  \begin{verbatim}
!   <output directory>/<model name>/<yyyy>/<yyyymm>/LDT_HIST_<yyyymmddhhmnss>.<extension>
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
   character(len=LDT_CONST_PATH_LEN)       :: dname
   character(len=LDT_CONST_PATH_LEN), save :: out_fname
   character(len=LDT_CONST_PATH_LEN)       :: odir_temp
   integer                  :: i, c

   if ( present(odir) ) then
      odir_temp = odir
   else
      odir_temp = LDT_rc%odir
   endif

   if(LDT_rc%wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir_temp)//'/'
      dname = trim(dname)//trim(model_name)//'/'
      
      write(unit=cdate, fmt='(i4.4)') LDT_rc%yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LDT_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( LDT_rc%wout )
      case ( "binary" )
         if(LDT_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LDT_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LDT_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized output format')
         call LDT_endrun 
      endselect
   elseif(LDT_rc%wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir_temp)//'/'
      dname = trim(dname)//trim(model_name)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') LDT_rc%yr, LDT_rc%mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LDT_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( LDT_rc%wout )
      case ("binary")
         if(LDT_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LDT_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LDT_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized LDT_rc%wout value')
         call LDT_endrun 
      endselect
   elseif(LDT_rc%wstyle.eq."2 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir_temp)//'/'
      dname = trim(dname)//trim(model_name)//'/'
      
      out_fname = trim(dname)//'LDT_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( LDT_rc%wout )
      case ("binary")
         if(LDT_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(LDT_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(LDT_rc%wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized LDT_rc%wout value')
         call LDT_endrun 
      endselect
   elseif(LDT_rc%wstyle.eq."WMO convention") then 
      write(unit=fint,fmt='(i2.2)') nint(writeint)/3600
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, LDT_rc%da
      
      write(unit=cdate, fmt='(i2.2, i2.2)') LDT_rc%hr, LDT_rc%mn
      
      if(LDT_rc%lis_map_proj(n).eq."polar") then 
         fproj = 'P'
         write(LDT_logunit,*) "[INFO] fres ",LDT_rc%gridDesc(n, 9)
         if (LDT_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LDT_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LDT_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LDT_rc%lis_map_proj(n).eq."lambert") then 
         fproj = 'L'
         write(LDT_logunit,*)"[INFO] fres ",LDT_rc%gridDesc(n, 9)
!         write(unit=fres, fmt='(f2.0)') LDT_rc%gridDesc(n, 9)
         write(unit=fres, fmt='(f3.0)') LDT_rc%gridDesc(n, 9)
         if (LDT_rc%gridDesc(n, 9) .ge. 10.) then
            write(unit=fres, fmt='(i2)') nint(LDT_rc%gridDesc(n, 9))
         else
            write(unit=fres, fmt='(i1)') nint(LDT_rc%gridDesc(n, 9))
         endif
         fres2 = trim(fres)//'KM'
      elseif(LDT_rc%lis_map_proj(n).eq."mercator") then 
         fproj = 'M'
         write(unit=fres, fmt='(i2.2)') LDT_rc%gridDesc(n, 9)
         fres = trim(fres)//'KM'
      elseif(LDT_rc%lis_map_proj(n).eq."gaussian") then 
         fproj = 'G'
         write(unit=fres, fmt='(i2.2)') LDT_rc%gridDesc(n, 9)*100        
         fres = '0P'//trim(fres)//'DEG'
      else
         fproj = 'C'
         write(unit=fres, fmt='(i10)') nint(LDT_rc%gridDesc(n,10)*100)
         read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
         c = 0 
         do i=1,10
            if(fres1(i).ne.' '.and.c==0) c = i
         enddo
         if (LDT_rc%gridDesc(n,10) .lt. 0.1) then
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
           '/PS.AFWA_SC.'//trim(LDT_rc%security_class)//&
           '_DI.'//trim(LDT_rc%distribution_class)//&
           '_DC.'//trim(LDT_rc%data_category)//'_GP.LDT_GR.'//&
           trim(fproj)//trim(fres2)//'_AR.'//trim(LDT_rc%area_of_data)//&
           '_PA.'//trim(fint)//'-HR-SUM_DD.'//&
           trim(cdate1)//'_DT.'//trim(cdate)//'_DF'
      select case (LDT_rc%wout)
      case ("binary")
         if(LDT_rc%wopt.eq."1d tilespace") then 
            out_fname = trim(dname)//'.DAT'
         elseif(LDT_rc%wopt.eq."2d gridspace") then 
            out_fname = trim(dname)//'.DAT'
         elseif(LDT_rc%wopt.eq."1d gridspace") then 
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
   endif
   fname = out_fname
 end subroutine LDT_create_output_filename


!BOP
!
! !ROUTINE: LDT_create_stats_filename
! \label{LDT_create_stats_filename}
! 
! !INTERFACE:
subroutine LDT_create_stats_filename(n, fname, mname, style)

! !USES:
   use LDT_coreMod,  only : LDT_rc
  
   implicit none 
! !ARGUMENTS:
   integer, intent(in)           :: n
   character(len=*), intent(out) :: fname
   character(len=*), intent(in)  :: mname
   integer, intent(in), optional :: style ! output directory style
!
! !DESCRIPTION:  
!  Create the file name for the stats files.  The convention used
!  in LIS creates a restart filename in the following format. 
!
!  <output directory>/EXP<expno>/<model name>stats.dat
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
   character(len=200) :: out_fname
   character*100      :: cdate
   integer            :: style_temp 

   if(.not. PRESENT(style)) then 
      style_temp = 1
   else
      style_temp = style
   endif

   if(style_temp.eq.1) then 

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(LDT_rc%odir)//'/EXP'//trim(LDT_rc%expcode)&
           //'/'//trim(mname)//trim(cdate)//'.stats'
      
      fname = trim(out_fname)
   elseif(style_temp.eq.4) then 
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(LDT_rc%odir)//'/'//trim(mname)//trim(cdate)//'.stats'
      
      fname = trim(out_fname)
   endif
 end subroutine LDT_create_stats_filename

!BOP
!
! !ROUTINE: LDT_create_daobs_filename
! \label{LDT_create_daobs_filename}
! 
! !INTERFACE:
subroutine LDT_create_daobs_filename(n, fname)
! !USES:
   use LDT_coreMod,  only : LDT_rc
  
   implicit none 
! !ARGUMENTS:
   integer, intent(in)           :: n
   character(len=*), intent(out) :: fname
!
! !DESCRIPTION:  
!  Create the file name for the daobs files.  The convention used
!  in LIS creates a restart filename in the following format. 
!
!  <output directory>/EXP<expno>/<model name>daobs.dat
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
   character(len=200) :: out_fname
   character*100      :: cdate
   integer            :: style_temp 

   
   write(unit=cdate, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
        LDT_rc%yr, LDT_rc%mo, &
        LDT_rc%da, LDT_rc%hr,LDT_rc%mn
   
   
   out_fname = trim(LDT_rc%odir)//'/EXP'//trim(LDT_rc%expcode)&
        //'/DAOBS/'//trim(cdate)//'.1gs4r'
   
   fname = trim(out_fname)
   
 end subroutine LDT_create_daobs_filename

!BOP
! !ROUTINE: LDT_readDomainConfigSpecs
! \label{LDT_readDomainConfigSpecs}
! 
! !INTERFACE: 
 subroutine LDT_readDomainConfigSpecs(segment_name, param_proj, domain_info)

! !USES: 
   use ESMF
   use LDT_coreMod,  only : LDT_rc, LDT_config
   use LDT_logMod,   only : LDT_verify, LDT_logunit, LDT_endrun

! !ARGUMENTS: 
   character(len=*),  intent(in)    :: segment_name
   character(50),     intent(in)    :: param_proj
   real,              intent(inout) :: domain_info(LDT_rc%nnest,20)
! 
! !DESCRIPTION: 
!   This subroutine reads the relevant section of the lis.config file 
!   that defines the domain specifications for a particular surface 
!   parameter dataset, based on the map projection used to define 
!   the surface parameter dataset. 
!EOP
   integer  :: i, rc

   domain_info = 0.

! ---
  select case ( param_proj )

    case ( "latlon" )  

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" lower left lat:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,4),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' lower left lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" lower left lon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,5),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' lower left lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" upper right lat:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,7),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' upper right lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" upper right lon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,8),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' upper right lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" resolution (dx):",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,9),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' resolution (dx):')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" resolution (dy):",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,10),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' resolution (dy):')
      enddo

      do i=1,LDT_rc%nnest
      ! glpnc:
        domain_info(i,2) = nint((domain_info(i,8)-domain_info(i,5))/domain_info(i,9) ) + 1 
      ! glpnr:
        domain_info(i,3) = nint((domain_info(i,7)-domain_info(i,4))/domain_info(i,10)) + 1
        domain_info(i,6)  = 128
        domain_info(i,11) = 64
        domain_info(i,20) = 64
      enddo

! -- Added Lambert parameter domain (KRA)

    case ( "lambert" )

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" lower left lat:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,4),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' lower left lat:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" lower left lon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,5),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' lower left lon:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" true lat1:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,10),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' true lat1:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" true lat2:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,7),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' true lat2:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" standard lon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,11),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' standard lon:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" resolution:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,8),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' resolution:')
         domain_info(i,9) = domain_info(i,8)   ! set dx = dy for now
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" x-dimension size:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,2),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' x-dimension size:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" y-dimension size:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,3),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' y-dimension size:')
      enddo

      do i=1,LDT_rc%nnest
         domain_info(i,6) = 8.0
      enddo

! ---
    case ( "gaussian" )  

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" first grid point lat:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,4),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' first grid point lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" first grid point lon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,5),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' first grid point lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" last grid point lat:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,7),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' last grid point lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" last grid point lon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,8),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' last grid point lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" resolution dlon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,9),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' resolution dlon:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" number of lat circles:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,3),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' number of lat circles:')
         ! "NUMBER OF LAT CIRCLES" represents the number of latitudes from 
         ! a pole to the equator.
         ! LDT wants the total number of latitudes, hence multiply by 2.
         domain_info(i,3) = 2*domain_info(i,3)
      enddo

! ---
    case ( "polar" )  
      write(LDT_logunit,*) '[ERR] Polar stereographic case not supported currently'
      call LDT_endrun
#if 0       
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" lower left lat:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,1),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' lower left lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" lower left lon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,2),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' lower left lon:')
      enddo
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" true lat:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,3),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' true lat:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" standard lon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,4),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' standard lon:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" orientation:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,5),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' orientation:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" resolution:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,6),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' resolution:')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" x-dimension size:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,7),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' x-dimension size:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" y-dimension size:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,8),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' y-dimension size:')
      enddo
#endif
! ---
    case ( "hrap" )
      write(LDT_logunit,*) " Reading in hrap domain for parameter: ",trim(segment_name)
      return
 
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" lower left lat:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,1),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' lower left lat:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" lower left lon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,2),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' lower left lon:')
      enddo
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" true lat:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,3),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' true lat:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" standard lon:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,4),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' standard lon:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" orientation:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,5),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' orientation:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" resolution:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,6),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' resolution:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" x-dimension size:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,7),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' x-dimension size:')
      enddo

      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" y-dimension size:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,8),rc=rc)
         call LDT_verify(rc,'please specify '//trim(segment_name)//' y-dimension size:')
      enddo

! ---

    case ( "UTM" )
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" UTM zone:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,1),rc=rc)
         call LDT_verify(rc,trim(segment_name)//' UTM zone: not defined')
      enddo
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" northing of SW corner:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,2),rc=rc)
         call LDT_verify(rc,trim(segment_name)//' northing of SW corner: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" easting of SW corner:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,3),rc=rc)
         call LDT_verify(rc,trim(segment_name)//' easting of SW corner: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" x-dimension size:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,4),rc=rc)
         call LDT_verify(rc,trim(segment_name)//' x-dimension size: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" y-dimension size:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,5),rc=rc)
         call LDT_verify(rc,trim(segment_name)//' y-dimension size: not defined')
      enddo
      
      call ESMF_ConfigFindLabel(LDT_config,trim(segment_name)//" resolution:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,domain_info(i,6),rc=rc)
         call LDT_verify(rc,trim(segment_name)//' resolution: not defined')
      enddo

!- Non-supported projections (at this time):
    case default 
      write(LDT_logunit,*) "[ERR] This type of run domain or parameter projection: ",&
          trim(param_proj)
      write(LDT_logunit,*) " is not supported.  Please change to one of these:"
      write(LDT_logunit,*) " -- latlon, lambert, polar, gaussian, mercator, hrap."
      write(LDT_logunit,*) " Program stopping ..."
      call LDT_endrun
  end select

 end subroutine LDT_readDomainConfigSpecs


!BOP
! !ROUTINE: LDT_transform_paramgrid
!  \label{LDT_transform_paramgrid}
! 
! !INTERFACE: 
!
 subroutine LDT_transform_paramgrid(n, gridtransform_opt, param_gridDesc, &
      numinpts, numtiles, array_in, li,               &
      numoutpts, array_out, lo )
!
! !USES: 
   use LDT_coreMod,  only : LDT_rc, LDT_domain, LDT_localPet
   use LDT_logMod,   only : LDT_logunit, LDT_endrun

   implicit none
!
! !ARGUMENTS: 
   integer,      intent(IN) :: n
   character(50),intent(IN) :: gridtransform_opt
   real,         intent(IN) :: param_gridDesc(20)
   integer,      intent(IN) :: numinpts
   integer,      intent(IN) :: numtiles
   logical*1,    intent(IN) :: li(numinpts)
   real,         intent(IN) :: array_in(numinpts)
   integer,      intent(IN) :: numoutpts
   real,      intent(INOUT) :: array_out(numoutpts,numtiles)
   logical*1, intent(INOUT) :: lo(numoutpts,numtiles)
!
! !DESCRIPTION: 
!   This routine transforms spatially the original parameter grid
!    to the target LIS output grid.
!EOP  
!
   integer :: iret

!- Grid transform arrays:
   integer, allocatable     :: n11(:)     ! Map array for aggregating methods

   integer, allocatable     :: n113(:)    ! Map array for nearest neighbor interp

   integer, allocatable     :: n111(:)    ! Map array for bilinear interp
   integer, allocatable     :: n121(:)
   integer, allocatable     :: n211(:)
   integer, allocatable     :: n221(:)
   real, allocatable        :: w111(:),w121(:)
   real, allocatable        :: w211(:),w221(:)

   integer, allocatable     :: n112(:,:)  ! Map array for budget/conservative interp
   integer, allocatable     :: n122(:,:)
   integer, allocatable     :: n212(:,:)
   integer, allocatable     :: n222(:,:)
   real, allocatable        :: w112(:,:),w122(:,:)
   real, allocatable        :: w212(:,:),w222(:,:)

! __________________________________________________________________

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
    if( gridtransform_opt == "average" .or. &
        gridtransform_opt == "mode"    .or. &
        gridtransform_opt == "tile" )  then

       allocate( n11(numinpts) )

    !- Create mapping between parameter domain and LIS grid domain:
       call upscaleByAveraging_input( param_gridDesc, &
                     LDT_rc%gridDesc(n,:), numinpts, numoutpts, n11 )
    end if

 !- Select grid spatial transform option:
    select case ( gridtransform_opt )

   !- Aggregate by calculating average of each output gridcell:
      case ( "average" )
        write(LDT_logunit,*) "#Regridding: Applying average to input parameter"
        call upscaleByAveraging( numinpts, numoutpts, LDT_rc%udef, n11,  &
                                 li, array_in, lo(:,1), array_out(:,1) )

   !- Select dominant type using the mode statistic:
      case ( "mode" )
        write(LDT_logunit,*) "#Regridding: Applying mode to input parameter"
        call upscaleByMode( numinpts, numoutpts, LDT_rc%udef, n11, &
                            li, array_in, lo(:,1), array_out(:,1) )

   !- Aggregate via tiling by counting up pixels within gridcell:
      case ( "tile" )
        write(LDT_logunit,*) "#Regridding: Applying tiling to input parameter"
        call upscaleByCnt( numinpts, numoutpts, numtiles, LDT_rc%udef, n11, &
                           li, array_in, lo, array_out )


   !- Select neighboring point:
      case ( "neighbor" )
         allocate( n113(numoutpts) )

         write(LDT_logunit,*) "#Regridding: Applying nearest neighbor to input parameter"
         call neighbor_interp_input( n, param_gridDesc, n113 )

         call neighbor_interp( LDT_rc%gridDesc(n,:), li, array_in,   &
                       lo(:,1), array_out(:,1), numinpts, numoutpts, &
                       LDT_domain(n)%lat, LDT_domain(n)%lon, &
                       n113, LDT_rc%udef, iret )

      case ( "bilinear" )
         allocate( n111(numoutpts), n121(numoutpts), n211(numoutpts), n221(numoutpts) )
         allocate( w111(numoutpts), w121(numoutpts), w211(numoutpts), w221(numoutpts) )

         write(LDT_logunit,*) "#Regridding: Applying bilinear interp to input parameter"
         call bilinear_interp_input( n, param_gridDesc, &
              n111, n121, n211, n221, w111, w121, w211, w221)

         call bilinear_interp(LDT_rc%gridDesc(n,:), li, array_in, lo(:,1), array_out(:,1), &
              numinpts, numoutpts, LDT_domain(n)%lat, LDT_domain(n)%lon,&
              w111, w121, w211, w221, n111, n121, n211, n221,&
              LDT_rc%udef, iret)

      case ( "budget-bilinear" )
         allocate( n112(numoutpts,25), n122(numoutpts,25), n212(numoutpts,25), n222(numoutpts,25) )
         allocate( w112(numoutpts,25), w122(numoutpts,25), w212(numoutpts,25), w222(numoutpts,25) )

         write(LDT_logunit,*) "#Regridding: Applying budget-bilinear interp to input parameter"
         call conserv_interp_input( n, param_gridDesc, &
              n112, n122, n212, n222, w112, w122, w212, w222 )

         call conserv_interp( LDT_rc%gridDesc(n,:), li, array_in, lo(:,1), array_out(:,1), &
              numinpts, numoutpts, LDT_domain(n)%lat, LDT_domain(n)%lon,&
              w112, w122, w212, w222, n112, n122, n212, n222, &
              LDT_rc%udef,iret)

   !- When no transform is performed:
      case ( "none" )
         write(LDT_logunit,*) "[INFO] No aggregation applied for parameter file ... "
         array_out(:,1) = array_in(:)

      case default
         write(LDT_logunit,*) "[ERR] This spatial transformation option ("//trim(gridtransform_opt)//") "
         write(LDT_logunit,*) " is not currently supported."
         write(LDT_logunit,*) " Program stopping ..."
         call LDT_endrun
    end select

 !- Deallocate arrays:
    if ( gridtransform_opt == "average" .or.  &
         gridtransform_opt == "mode"    .or.  &
         gridtransform_opt == "tile" ) then
       deallocate( n11 )

    elseif ( gridtransform_opt == "neighbor" ) then
       deallocate( n113 )

    elseif ( gridtransform_opt == "bilinear" ) then
       deallocate( n111, n121, n211, n221, w111, w121, w211, w221 )

    elseif ( gridtransform_opt == "budget-bilinear" ) then
       deallocate( n112, n122, n212, n222, w112, w122, w212, w222 )
    endif


 end subroutine LDT_transform_paramgrid


!BOP
! !ROUTINE: readparam_real_2d
!  \label{readparam_real_2d}
! 
! !INTERFACE: 
 subroutine readparam_real_2d( n, ftn, param_proj, gridtransform_opt, &
                      readparam_gridDesc, numtiles, array )

! !USES: 
   use LDT_coreMod, only : LDT_rc, LDT_domain, LDT_localPet
   use LDT_logMod,  only : LDT_logunit, LDT_endrun
   use LDT_gridmappingMod

   implicit none
! !ARGUMENTS: 
   integer,  intent(IN)     :: n 
   integer,  intent(IN)     :: ftn
   character(50),intent(IN) :: param_proj
   character(50),intent(IN) :: gridtransform_opt
   real,     intent(IN)     :: readparam_gridDesc(20)
   integer,  intent(IN)     :: numtiles
   real,     intent(INOUT)  :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),numtiles)
!
! !DESCRIPTION: 
!   This routine retrieves a 2D real data from a binary direct access file. 
!EOP  
   logical :: file_exists
   integer :: c, r, t, i, line, iret
   integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
   integer :: mi                       ! Total number of input param grid array points
   integer :: mo                       ! Total number of output LIS grid array points
   real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)
   real,    allocatable  :: read_inputparm(:,:)   ! Read input parameter
   real,    allocatable  :: gi(:)      ! input parameter 1d grid
   logical*1,allocatable :: li(:)      ! input logical mask (to match gi)
   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n),numtiles)  ! output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n),numtiles)  ! output logical mask (to match go1)
! __________________________________________________________________________________

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

!- Map Parameter Grid Info to LIS Target Grid/Projection Info --
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, param_proj, readparam_gridDesc(:), &
               glpnc, glpnr, subpnc, subpnr, subparam_gridDesc,    &
               lat_line, lon_line )

! -------------------------------------------------------------------

  array = LDT_rc%udef

! ---  Parameter Projection  ---
  select case ( param_proj )

!- Lat/lon:
   case ( "latlon" ) 

  !- Initialize parameter read-in array:
     allocate( read_inputparm(subpnc, subpnr) )
     read_inputparm = LDT_rc%udef

     mi = subpnc*subpnr
     allocate( li(mi), gi(mi) )
     li = .false.
     gi = LDT_rc%udef
     mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
     lo1 = .false.

! -------------------------------------------------------------------
!   READ IN 2-D FIELDS 
! -------------------------------------------------------------------
     line = 0
     do r = 1, subpnr
        do c = 1, subpnc
           line = (lat_line(c,r)-1)*glpnc + lon_line(c,r)
           read(ftn,rec=line) read_inputparm(c,r)
        enddo
     enddo

  !- Assign 2-D array to 1-D for grid transformation routines:
     i = 0
     do r = 1, subpnr
        do c = 1, subpnc;  i = i + 1
           gi(i) = read_inputparm(c,r)
           if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
        enddo
     enddo

!- Select grid spatial transform option:
   if( gridtransform_opt == "none" ) then
      write(LDT_logunit,*) " No aggregation applied for parameter file ... "
      array(:,:,1) = read_inputparm(:,:)
   
   else
  !- Transform parameter from original grid to LIS output grid:
     call LDT_transform_paramgrid( n, gridtransform_opt, &
              subparam_gridDesc, mi, numtiles, gi, li, mo, go1, lo1 )

   endif
 
!- Convert 1D to 2D grid arrays:
   if( gridtransform_opt .ne. "none" ) then
      i = 0
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n);  i = i + 1
            do t = 1, numtiles
               array(c,r,t) = go1(i,t)
            end do
         enddo
      enddo
   end if

!- Deallocate arrays:
   deallocate( li, gi )
   deallocate( lat_line, lon_line )
   deallocate( read_inputparm )

! ----------------------------------------------------------------------

#if 0
!- Gaussian:
   case ( "gaussian" ) 

!KRA -- Need to modify bottom two lines:
      line1 = gaussian_find_row(LDT_rc%gridDesc(n,4))  -   &
!KRA           gaussian_find_row(LDT_rc%gridDesc(n,44)) + 1
              gaussian_find_row(LDT_rc%gridDesc(n,1)) + 1
      line2 = gaussian_find_col(LDT_rc%gridDesc(n,5))  -   &
!KRA           gaussian_find_col(LDT_rc%gridDesc(n,45)) + 1
              gaussian_find_col(LDT_rc%gridDesc(n,2)) + 1

      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            pnc = line2+c-1
            pnr = line1+r-1
!KRA            line = (pnr-1)*nint(LDT_rc%gridDesc(n,42))+pnc
            line = (pnr-1)*nint(360/LDT_rc%gridDesc(n,5))+pnc

         !- Read in 2D parameter data:
            read(ftn,rec=line) array(c,r)

         enddo
      enddo

   case ( "UTM" ) 
!rlat/rlon used here to store northing and easting
      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            rlat(c,r) = LDT_rc%gridDesc(n,4)+(r-1)*LDT_rc%gridDesc(n,9)
            rlon(c,r) = LDT_rc%gridDesc(n,5)+(c-1)*LDT_rc%gridDesc(n,9)
         enddo
      enddo
      nc_dom = param_gridDesc(4)

      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            line1 = nint((rlat(c,r)-param_gridDesc(2))/param_gridDesc(6))+1
            line2 = nint((rlon(c,r)-param_gridDesc(3))/param_gridDesc(6))+1
            line = (line1-1)*nc_dom + line2

         !- Read in 2D parameter data:
            read(ftn,rec=line) array(c,r) 
         enddo
      enddo
#endif
    case default
      write(LDT_logunit,*) "[ERR] This parameter projection is not supported (in readparam_real_2d)"
      write(LDT_logunit,*) " Program stopping ..."
      call LDT_endrun
   end select 

 end subroutine readparam_real_2d


!BOP
! !ROUTINE: readparam_int_2d
!  \label{readparam_int_2d}
! 
! !INTERFACE: 
 subroutine readparam_int_2d(n, ftn, param_proj, gridtransform_opt, &
                           readparam_gridDesc, numtiles, array)

! !USES: 
   use LDT_coreMod,  only : LDT_rc, LDT_domain, LDT_localPet
   use LDT_logMod,   only : LDT_logunit, LDT_endrun
   use LDT_gridmappingMod

   implicit none
! !ARGUMENTS: 
   integer,  intent(IN)     :: n 
   integer,  intent(IN)     :: ftn
   character(50),intent(IN) :: param_proj
   character(50),intent(IN) :: gridtransform_opt
   real,     intent(IN)     :: readparam_gridDesc(20)
   integer,  intent(IN)     :: numtiles
   integer,  intent(INOUT)  :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),numtiles)
!
! !DESCRIPTION: 
!   This routine retrieves a 2D integer data from a binary direct access file. 
!EOP  
   logical :: file_exists
   integer :: c, r, t, i, line, iret
   integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
   integer :: mi                       ! Total number of input param grid array points
   integer :: mo                       ! Total number of output LIS grid array points
   real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)
   real,    allocatable  :: read_inputparm(:,:)   ! Read input parameter
   integer, allocatable  :: read_inputparm_int(:,:)   ! Read input parameter
   real,    allocatable  :: gi(:)      ! input parameter 1d grid
   logical*1,allocatable :: li(:)      ! input logical mask (to match gi)
   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n),numtiles)  ! output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n),numtiles)  ! output logical mask (to match go)
! __________________________________________________________________________________

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

!- Map Parameter Grid Info to LIS Target Grid/Projection Info --
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, param_proj, readparam_gridDesc(:), &
               glpnc, glpnr, subpnc, subpnr, subparam_gridDesc,    &
               lat_line, lon_line )

! -------------------------------------------------------------------
  array = LDT_rc%udef

! ---  Parameter Projection  ---
  select case ( param_proj )

!- Lat/lon:
   case ( "latlon" ) 

  !- Initialize parameter read-in array:
     allocate( read_inputparm_int(subpnc, subpnr) )
     allocate( read_inputparm(subpnc, subpnr) )
     read_inputparm = LDT_rc%udef

     mi = subpnc*subpnr
     allocate( li(mi), gi(mi) )
     li = .false.
     gi = LDT_rc%udef
     mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
     lo1 = .false.

! -------------------------------------------------------------------
!     READ IN 2-D FIELDS (NON-TILED OPTIONS)
! -------------------------------------------------------------------
     line = 0
     do r = 1, subpnr
        do c = 1, subpnc
           line = (lat_line(c,r)-1)*glpnc + lon_line(c,r)
           read(ftn,rec=line) read_inputparm_int(c,r)
        enddo
     enddo
     read_inputparm = read_inputparm_int

  !- Assign 2-D array to 1-D for aggregation routines:
     i = 0
     do r = 1, subpnr
        do c = 1, subpnc;  i = i + 1
           gi(i) = read_inputparm(c,r)
           if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
        enddo
     enddo

!- Select grid spatial transform option:
   if( gridtransform_opt == "none" ) then
      write(LDT_logunit,*) " No aggregation applied for parameter file ... "
      array(:,:,1) = read_inputparm(:,:)

   else
  !- Transform parameter from original grid to LIS output grid:
     call LDT_transform_paramgrid( n, gridtransform_opt, &
              subparam_gridDesc, mi, numtiles, gi, li, mo, go1, lo1 )

   endif

!- Convert 1D to 2D grid arrays:
   if ( gridtransform_opt .ne. "none" ) then
      i = 0
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n);  i = i + 1
            do t = 1, numtiles
               array(c,r,t) = nint(go1(i,t))
            enddo
         enddo
      enddo
   end if

!- Deallocate arrays:
   deallocate( li, gi )
   deallocate( lat_line, lon_line )
   deallocate( read_inputparm )
   deallocate( read_inputparm_int )

! ----------------------------------------------------------------------

#if 0
!- Gaussian:
   case ( "gaussian" ) 

!KRA -- Need to modify bottom two lines:
      line1 = gaussian_find_row(LDT_rc%gridDesc(n,4))  -   &
!KRA           gaussian_find_row(LDT_rc%gridDesc(n,44)) + 1
              gaussian_find_row(LDT_rc%gridDesc(n,1)) + 1
      line2 = gaussian_find_col(LDT_rc%gridDesc(n,5))  -   &
!KRA           gaussian_find_col(LDT_rc%gridDesc(n,45)) + 1
              gaussian_find_col(LDT_rc%gridDesc(n,2)) + 1

      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            pnc = line2+c-1
            pnr = line1+r-1
!KRA            line = (pnr-1)*nint(LDT_rc%gridDesc(n,42))+pnc
            line = (pnr-1)*nint(360/LDT_rc%gridDesc(n,5))+pnc

         !- Read in 2D parameter data:
            read(ftn,rec=line) array(c,r)

         enddo
      enddo

   case ( "UTM" ) 
!rlat/rlon used here to store northing and easting
      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            rlat(c,r) = LDT_rc%gridDesc(n,4)+(r-1)*LDT_rc%gridDesc(n,9)
            rlon(c,r) = LDT_rc%gridDesc(n,5)+(c-1)*LDT_rc%gridDesc(n,9)
         enddo
      enddo
      nc_dom = param_gridDesc(4)

      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            line1 = nint((rlat(c,r)-param_gridDesc(2))/param_gridDesc(6))+1
            line2 = nint((rlon(c,r)-param_gridDesc(3))/param_gridDesc(6))+1
            line = (line1-1)*nc_dom + line2

         !- Read in 2D parameter data:
            read(ftn,rec=line) array(c,r) 
         enddo
      enddo
#endif
    case default
      write(LDT_logunit,*) "[ERR] This parameter projection is not supported (in readparam_int_2d)"
      write(LDT_logunit,*) " Program stopping ..."
      call LDT_endrun
   end select 

 end subroutine readparam_int_2d


!BOP
! 
! !ROUTINE: readldtparam_real_2d
! \label{readldtparam_real_2d}
!
! !INTERFACE:
  subroutine readldtparam_real_2d(n,pname,array)
! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_coreMod, only : LDT_rc, LDT_localPet,&
          LDT_domain, LDT_ews_ind, LDT_ewe_ind,&
          LDT_nss_ind, LDT_nse_ind, LDT_ews_halo_ind,LDT_ewe_halo_ind, &
          LDT_nss_halo_ind, LDT_nse_halo_ind
  use LDT_logMod,  only : LDT_logunit, LDT_endrun, LDT_verify
  use LDT_paramDataMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in)    :: n
  character(len=*)       :: pname
  real,    intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))
! 
! !DESCRIPTION:
!  Reads a parameter field from the input LIS/LDT parameter file
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
  integer :: ios1
  integer :: ios,nid,paramid,ncId, nrId
  integer :: nc,nr,c,r
  real    :: param(LDT_rc%gnc(n),LDT_rc%gnr(n))
  logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
!  inquire(file=trim(LDT_rc%paramfile(n)), exist=file_exists)
  inquire(file=trim(LDT_LSMparam_struc(n)%param_filename), exist=file_exists)
  if(file_exists) then

     write(LDT_logunit,*)'Reading '//trim(pname)//' map ...'
!     ios = nf90_open(path=trim(LDT_rc%paramfile(n)),&
     ios = nf90_open(path=trim(LDT_LSMparam_struc(n)%param_filename),&
          mode=NF90_NOWRITE,ncid=nid)
     call LDT_verify(ios,'Error in nf90_open in readldtparam_real_2d')

     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LDT_verify(ios,'Error in nf90_inq_dimid in readldtparam_real_2d')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LDT_verify(ios,'Error in nf90_inq_dimid in readldtparam_real_2d')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LDT_verify(ios,'Error in nf90_inquire_dimension in readldtparam_real_2d')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LDT_verify(ios,'Error in nf90_inquire_dimension in readldtparam_real_2d')

     ios = nf90_inq_varid(nid,trim(pname),paramid)
     call LDT_verify(ios,trim(pname)//' field not found in the LDT param file')

     ios = nf90_get_var(nid,paramid,param)
     call LDT_verify(ios,'Error in nf90_get_var in readldtparam_real_2d')

     ios = nf90_close(nid)
     call LDT_verify(ios,'Error in nf90_close in readldtparam_real_2d')

     array(:,:) = &
          param(LDT_ews_halo_ind(n,LDT_localPet+1):&
          LDT_ewe_halo_ind(n,LDT_localPet+1), &
          LDT_nss_halo_ind(n,LDT_localPet+1): &
          LDT_nse_halo_ind(n,LDT_localPet+1))

     write(LDT_logunit,*)'[INFO] Successfully read '//trim(pname)//' map '
  else
     write(LDT_logunit,*)'[ERR] '//trim(pname)//' map: ',&
          trim(LDT_LSMparam_struc(n)%param_filename), ' does not exist'
     write(LDT_logunit,*) ' Program stopping ...'
     call LDT_endrun
  endif

#endif

end subroutine readldtparam_real_2d


!BOP
! !ROUTINE: LDT_checkDomainExtents
! \label{LDT_checkDomainExtents}
! 
! !INTERFACE: 
 subroutine LDT_checkDomainExtents(n, data_gridDesc)

! !USES: 
   use LDT_coreMod,  only : LDT_rc, LDT_domain
   use LDT_logMod,   only : LDT_logunit, LDT_endrun
   use map_utils,    only : ij_to_latlon

! !ARGUMENTS: 
   integer, intent(in) :: n 
   real,    intent(in) :: data_gridDesc(20)
! 
! !DESCRIPTION: 
!   This subroutine checks to see if the domain extents of the data 
!    are within the LDT running domain. 
!EOP
   real     :: max_lat, max_lon
   real     :: min_lat, min_lon
   real     :: rlat, rlon
   integer  :: c,r

   max_lat = -10000
   max_lon = -10000
   min_lat = 10000 
   min_lon = 10000
   
   do r=1,LDT_rc%lnr(n)
      do c=1,LDT_rc%lnc(n)
         call ij_to_latlon(LDT_domain(n)%ldtproj,float(c),float(r),&
              rlat,rlon)
         if(rlat.gt.max_lat) max_lat = rlat
         if(rlon.gt.max_lon) max_lon = rlon
         if(rlat.lt.min_lat) min_lat = rlat
         if(rlon.lt.min_lon) min_lon = rlon
      enddo
   enddo

   if((min_lat.lt.data_gridDesc(4)).or.&
      (min_lon.lt.data_gridDesc(5)).or.&
      (max_lat.gt.data_gridDesc(7)).or.&
      (max_lon.gt.data_gridDesc(8))) then 
     write(LDT_logunit,*) "[ERR] The parameter data grid is specified only for the CONUS ..."
     write(LDT_logunit,*) " The LIS running domain is outside the input parameter grid boundaries ..."
     call LDT_endrun
   endif
   
 end subroutine LDT_checkDomainExtents

!BOP
!
! !ROUTINE: create_restart_filename_withtime
! \label{create_restart_filename_withtime}
!
! !INTERFACE:
 subroutine create_restart_filename_withtime(n, rootdir,dir_name,model_name,&
      yr,mo,da,hr,mn,ss, wstyle, wformat,fname)
! !USES:

   implicit none 
! !ARGUMENTS:
   integer, intent(in)           :: n
   character(len=*), intent(in)  :: rootdir
   character(len=*), intent(in)  :: dir_name 
   character(len=*), intent(in)  :: model_name 
   integer         , intent(in)  :: yr, mo, da, hr, mn, ss
   character(len=*), intent(in)  :: wstyle
   character(len=*), intent(in)  :: wformat
   character(len=*), intent(out) :: fname


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
   character(len=100)       :: dname
   character(len=200)       :: out_fname


   if(wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr, mo, da, hr,mn
      
      dname = trim(rootdir)//'/'
      dname = trim(dname)//trim(dir_name)//'/'
      
      write(unit=cdate, fmt='(i4.4)') yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
           yr, mo, da
      dname = trim(dname)//trim(cdate)//'/'

      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      select case (wformat)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select
      fname = out_fname

   elseif(wstyle.eq."3 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr,mo,da,hr,mn
      
      dname = trim(rootdir)//'/'
      dname = trim(dname)//trim(dir_name)//'/'
      
      write(unit=cdate, fmt='(i4.4,i2.2)') yr, mo
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      select case (wformat)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select
      fname = out_fname
      
   elseif(wstyle.eq."2 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr,mo,da,hr,mn
      
      dname = trim(rootdir)//'/'
      dname = trim(dname)//trim(dir_name)//'/'
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      select case (wformat)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
      case ("grib1")
         out_fname = trim(out_fname)//'.GR1'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.GR2'
      case default            
      end select
      fname = out_fname
      
   elseif(wstyle.eq."WMO convention") then 

      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2,i2.2)') &
           yr,mo,da,hr,mn
      
      dname = rootdir
         
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n
      
      out_fname = trim(dname)//'/LIS_RST'//&
           '_'//trim(model_name)//'_'//cdate1//trim(cdate)
      
      select case (wformat)
      case ("binary")
         out_fname = trim(out_fname)//'.bin'
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


end module LDT_fileIOMod
