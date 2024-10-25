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
!BOP
! 
! !MODULE: LVT_fileIOMod
! \label(LVT_fileIOMod)
!
! !INTERFACE:
module LVT_fileIOMod
! 
! !USES:   
  implicit none 
  PRIVATE
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This module contains a number of routines useful for various file I/O
!   operations in LIS. The module 
!   provides routines to create output directories, filenames, that are 
!   be used in the model output routines. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008  Sujay Kumar;  Initial Specification
! 
!EOP
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LVT_create_output_filename   ! create an output filename  
  public :: LVT_create_output_directory  ! create the output directory
  public :: LVT_create_daobs_filename
  public :: LVT_convertParamDataToLocalDomain

!EOP


!BOP
! 
! !ROUTINE: LVT_create_output_filename
! \label{LVT_create_output_filename}
!
! !INTERFACE:
  interface LVT_create_output_filename 
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure create_output_filename 
     module procedure create_output_filename_with_timestamp
!
! !DESCRIPTION: 
!
!EOP
  end interface

!BOP
! 
! !ROUTINE: LVT_convertParamDataToLocalDomain
! \label{LVT_convertParamDataToLocalDomain}
!
! !INTERFACE:
  interface LVT_convertParamDataToLocalDomain
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure convertParam_int
     module procedure convertParam_real
!
! !DESCRIPTION: 
!  Routine to subset the data read from the parameter attributes
!  file supplied by LDT to the LVT domain being run. It is assumed
!  that the parameter data domain is at the same spatial resolution 
!  and map projection as that of the LVT domain. 
!
!EOP
  end interface

contains

!BOP
! 
! !ROUTINE: LVT_create_output_directory
! \label{LVT_create_output_directory}
!
! !INTERFACE:
subroutine LVT_create_output_directory(mname,dir_name,style)
! 
! !USES:
   use LVT_coreMod, only : LVT_rc
   use LVT_logMod,  only : LVT_log_msg, LVT_endrun
   implicit none 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:  
!  Create the output directory for the output data files. The call creates
!  a hierarchy of directories in the following format, if the directory 
!  name is not specified. 
!
!  style option 1: 
!  \begin{verbatim}
!  <output directory>/EXP<expno>/<model name>/<yr>/<yrmoda>
!  \end{verbatim}
!  style option 2: 
!  \begin{verbatim}
!   <output directory>/EXP<expno>/<model name>
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
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS:
   character(len=*)  :: mname
   character(len=*), optional   :: dir_name
   character(len=*), intent(IN), optional :: style
!EOP
   character(len=4) :: cdate
   character(len=8) :: cdate1
   character(len=200) :: out_dname
   character(len=50)  :: style_temp
   
   if(PRESENT(style)) then 
      style_temp = style
   else
      style_temp = "4 level hierarchy"
   endif
   
   if(style_temp.eq."4 level hierarchy") then

      out_dname = trim(LVT_rc%odir)//'/'
      
      out_dname = trim(out_dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4)') LVT_rc%yr
      out_dname = trim(out_dname)//trim(cdate)//'/'
      
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') LVT_rc%yr, LVT_rc%mo, LVT_rc%da
      out_dname = trim(out_dname)//trim(cdate1)
      
      if ( present(dir_name) ) then
         dir_name = trim(out_dname)
      else
         call system("mkdir -p "//trim(out_dname))
      endif
   elseif(style_temp.eq."3 level hierarchy") then 
      out_dname = trim(LVT_rc%odir)//'/'
      
      out_dname = trim(out_dname)//trim(mname)//'/'

      write(unit=cdate1, fmt='(i4.4, i2.2)') LVT_rc%yr, LVT_rc%mo
      out_dname = trim(out_dname)//trim(cdate1)

      if ( present(dir_name) ) then
         dir_name = trim(out_dname)
      else
         call system("mkdir -p "//trim(out_dname))
      endif      
   elseif(style_temp.eq."2 level hierarchy") then 
      out_dname = trim(LVT_rc%odir)//'/'
      
      out_dname = trim(out_dname)//trim(mname)//'/'

      if ( present(dir_name) ) then
         dir_name = trim(out_dname)
      else
         call system("mkdir -p "//trim(out_dname))
      endif      
   elseif((style_temp.eq."WMO convention").or.  &
          (style_temp.eq."WMO convention (AFW OPS)")) then 
      out_dname = trim(LVT_rc%odir)

      if ( present(dir_name) ) then
         dir_name = trim(out_dname)
      else
         call system("mkdir -p "//trim(out_dname))
      endif
   else
      call lvt_log_msg('ERR: LVT_create_output_directory --')
      call lvt_log_msg('  Unrecognized LIS output naming style:')
      call lvt_log_msg('      '//style_temp)
      call LVT_endrun 
   endif

 end subroutine LVT_create_output_directory

!BOP
! 
! !ROUTINE: create_output_filename
! \label{create_output_filename}
!
! !INTERFACE:
subroutine create_output_filename(n, source, fname, model_name, writeint, &
     wout, style,odir)
! 
! !USES:
   use LVT_coreMod,  only : LVT_rc, LVT_LIS_rc
   use LVT_logMod,   only : LVT_log_msg, LVT_endrun

   implicit none 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:  
!  Create the file name for LIS output data files. It creates both the GSWP
!  style of output filenames and the standard LIS style. The convention used
!  in LIS creates a filename in the following default format (style==1) 
!
!   <output directory>/EXP<expno>/<model name>/<yr>/<yrmoda>/<yrmodahrmn>.<extension>
!  Style option ==2 corresponds to the following style: 
!  
!   <output directory>/EXP<expno>/<model name>.<yrmodahrmn>.<extension>
!  Style option ==3 :
!   <output directory>/EXP<expno>/<model name>.<extension>
! Style option ==4 :
!   <output directory>/<AFWA Weather product style>
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
!   \item [wout]
!    output writing format
!   \item [style]
!    style option as described above
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
!
! !ARGUMENTS:
   integer, intent(in) :: n
   integer, intent(in) :: source
   character(len=*), intent(out)          :: fname
   character(len=*), intent(in), optional :: model_name ! needed for gswp run
   integer, intent(in), optional          :: writeint ! output writing interval
   character(len=*), intent(in), optional :: wout ! output format
   character(len=*), intent(in), optional :: style ! output directory style
   character(len=*), intent(in), optional :: odir
! 
!EOP
   character(len=8)        :: date
   character(len=10)       :: time
   character(len=5)        :: zone
   integer, dimension(8)   :: values
 
   character(len=10)       :: cdate
   character(len=12)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=200)       :: dname
   character(len=200), save :: out_fname
   character(len=50)        :: style_temp
   character(len=100)        :: odir_temp
   integer                  :: i, c

   if(.not.PRESENT(odir)) then 
      odir_temp = LVT_rc%odir
   else
      odir_temp = odir
   endif

   if(.not. PRESENT(style)) then 
      style_temp = "4 level hierarchy"
   else
      style_temp = style
   endif
   ! added by Shugong
   if(LVT_rc%lis_version == 6) then
     write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
          LVT_rc%dyr(source), LVT_rc%dmo(source), &
          LVT_rc%dda(source), LVT_rc%dhr(source),LVT_rc%dmn(source)
     
     dname = trim(LVT_rc%odir)//'/EXP'//trim(adjustl(LVT_rc%expcode))//'/'
     dname = trim(dname)//trim(LVT_rc%lsm)//'/'
     
     write(unit=cdate, fmt='(i4.4)') LVT_rc%dyr(source)
     dname = trim(dname)//trim(cdate)//'/'
     
     write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') LVT_rc%dyr(source), &
          LVT_rc%dmo(source), LVT_rc%dda(source)
     dname = trim(dname)//trim(cdate)
     
     out_fname = trim(dname)//'/'//cdate1
     
     write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
     out_fname = trim(out_fname)//trim(cdate)
     out_fname = trim(out_fname)//'.gs4r'
     fname = trim(out_fname) 
  else ! LIS version is 7 or even higher number
     if(style_temp.eq."4 level hierarchy") then 
        write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
             LVT_rc%dyr(source), LVT_rc%dmo(source), &
             LVT_rc%dda(source), LVT_rc%dhr(source),LVT_rc%dmn(source)

        dname = trim(odir_temp)//'/'
        dname = trim(dname)//trim(model_name)//'/'

        write(unit=cdate, fmt='(i4.4)') LVT_rc%dyr(source)
        dname = trim(dname)//trim(cdate)//'/'

        write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
             LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source)
        dname = trim(dname)//trim(cdate)

        out_fname = trim(dname)//'/LIS_HIST_'//cdate1

        write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
        out_fname = trim(out_fname)//trim(cdate)

        if(present(wout)) then 
           select case ( wout )
           case ( "binary")
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif
           case ( "distributed binary")
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif
           case ( "grib1" )
              out_fname = trim(out_fname)//'.grb'
           case ( "netcdf")
              out_fname = trim(out_fname)//'.nc'
           case ( "ascii" )
              out_fname = trim(out_fname)//'.txt'
           case default
              call lvt_log_msg('ERR: create_output_filename -- '// &
                   'Unrecognized output format')
              call LVT_endrun 
           endselect
        endif
     elseif(style_temp.eq."3 level hierarchy") then 
        write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
             LVT_rc%dyr(source), LVT_rc%dmo(source), &
             LVT_rc%dda(source), LVT_rc%dhr(source),LVT_rc%dmn(source)

        dname = trim(odir_temp)//'/'
        dname = trim(dname)//trim(model_name)//'/'

        write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') LVT_rc%dyr(source), &
             LVT_rc%dmo(source)
        dname = trim(dname)//trim(cdate)

        out_fname = trim(dname)//'/LIS_HIST_'//cdate1

        write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
        out_fname = trim(out_fname)//trim(cdate)

        if(present(wout)) then 
           select case ( wout )
           case ( "binary" )
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif
           case ( "distributed binary")
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif              
           case ( "grib1" )
              out_fname = trim(out_fname)//'.grb'
           case ( "netcdf" )
              out_fname = trim(out_fname)//'.nc'
           case ( "grib2" )
              out_fname = trim(out_fname)//'.gr2'
           case default
              call lvt_log_msg('ERR: create_output_filename -- '// &
                   'Unrecognized output format')
              call LVT_endrun 
           endselect
        endif
     elseif(style_temp.eq."2 level hierarchy") then 
        write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
             LVT_rc%dyr(source), LVT_rc%dmo(source), &
             LVT_rc%dda(source), LVT_rc%dhr(source),LVT_rc%dmn(source)

        dname = trim(odir_temp)//'/'
        dname = trim(dname)//trim(model_name)//'/'

        out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)

        write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
        out_fname = trim(out_fname)//trim(cdate)

        if(present(wout)) then 
           select case ( wout )
           case ( "binary")
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif
           case ( "distributed binary")
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif                            
           case ( "grib1" )
              out_fname = trim(out_fname)//'.grb'
           case ( "netcdf" )
              out_fname = trim(out_fname)//'.nc'
           case ( "grib2" )
              out_fname = trim(out_fname)//'.gr2'
           case default
              call lvt_log_msg('ERR: create_output_filename -- '// &
                   'Unrecognized output format')
              call LVT_endrun 
           endselect
        endif
     elseif(style_temp.eq."WMO convention") then 
        write(unit=fint,fmt='(i2.2)') writeint/3600
        write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
             LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source)

        write(unit=cdate, fmt='(i2.2, i2.2)') LVT_rc%dhr(source), LVT_rc%dmn(source)

        if(LVT_rc%domain.eq."polar") then 
           fproj = 'P'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."lambert") then 
           fproj = 'L'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."polar") then 
           fproj = 'M'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."gaussian") then 
           fproj = 'G'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)*100        
           fres = '0P'//trim(fres)//'DEG'
        else
           fproj = 'C'
           write(unit=fres, fmt='(i10)') nint(LVT_LIS_rc(source)%gridDesc(10)*100)
           read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
           c = 0 
           do i=1,10
              if(fres1(i).ne.' '.and.c==0) c = i
           enddo
           ! EMK...Make code consistent with LIS
!           fres3 = '0P'
           if (LVT_LIS_rc(source)%gridDesc(10) .lt. 0.1) then
              fres3 = '0P0'
           else
              fres3 = '0P'
           end if
           fres2 = trim(fres3)
           do i=c,10
              fres2 = trim(fres2)//trim(fres1(i))
           enddo
           fres2 = trim(fres2)//'DEG'
        endif

        ! EMK TEST...Remove date directory
!        dname = trim(odir_temp)//'/'//trim(cdate1)//&
        dname = trim(odir_temp)//'/'//&

             '/PS.AFWA_SC.' &
             //trim(LVT_rc%security_class)//'_DI.' &
             //trim(LVT_rc%distribution_class)//'_DC.' &
!             //trim(LVT_rc%data_category)//'_GP.LIS_GR.' &
             //'ANLYS'//'_GP.LIS_GR.' &
             //trim(fproj)//trim(fres2)//'_AR.'//trim(LVT_rc%area_of_data)//&
             '_PA.'//trim(fint)//'-HR-SUM_DD.'// &
             trim(cdate1)//'_DT.'//trim(cdate)//'_DF'
        select case (LVT_LIS_rc(source)%format)
        case ( "binary" )
           out_fname = trim(dname)//'.DAT'
        case ( "grib1" )
           out_fname = trim(dname)//'.GR1'
        case ( "netcdf" )
           out_fname = trim(dname)//'.nc'
        case ( "grib2" )
           out_fname = trim(dname)//'.GR2'
        case default            
        end select
     elseif(style_temp.eq."WMO convention (AFW OPS)") then 
        write(unit=fint,fmt='(i2.2)') writeint/3600
        write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
             LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source)

        write(unit=cdate, fmt='(i2.2, i2.2)') LVT_rc%dhr(source), LVT_rc%dmn(source)

        if(LVT_rc%domain.eq."polar") then 
           fproj = 'P'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."lambert") then 
           fproj = 'L'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."polar") then 
           fproj = 'M'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."gaussian") then 
           fproj = 'G'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)*100        
           fres = '0P'//trim(fres)//'DEG'
        else
           fproj = 'C'
           write(unit=fres, fmt='(i10)') nint(LVT_LIS_rc(source)%gridDesc(10)*100)
           read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
           c = 0 
           do i=1,10
              if(fres1(i).ne.' '.and.c==0) c = i
           enddo
           fres3 = '0P'
           fres2 = trim(fres3)
           do i=c,10
              fres2 = trim(fres2)//trim(fres1(i))
           enddo
           fres2 = trim(fres2)//'DEG'
        endif

        dname = trim(odir_temp)//'/'//trim(cdate1)//&
!        dname = trim(odir_temp)//&
             '/PS.AFWA_SC.'//trim(LVT_rc%security_class)//&
             '_DI.'//trim(LVT_rc%distribution_class)//& 
             '_GP.LIS_GR.'//&
             trim(fproj)//trim(fres2)//'_AR.'//trim(LVT_rc%area_of_data)//&
             '_PA.LIS_DD.'// &
             trim(cdate1)//'_DT.'//trim(cdate)//'_DF'
        select case (LVT_LIS_rc(source)%format)
        case ( "binary" )
           out_fname = trim(dname)//'.DAT'
        case ( "grib1" )
           out_fname = trim(dname)//'.GR1'
        case ( "netcdf" )
           out_fname = trim(dname)//'.nc'
        case ( "grib2" )
           out_fname = trim(dname)//'.GR2'
        case default            
        end select
     else
        call lvt_log_msg('ERR: LVT_create_output_filename --')
        call lvt_log_msg('  Unrecognized LIS output naming style:')
        call lvt_log_msg('      '//style_temp)
        call LVT_endrun 
     endif
     fname = trim(out_fname)
   endif 
 end subroutine Create_output_filename

!BOP
! 
! !ROUTINE: create_output_filename_with_timestamp
! \label{create_output_filename_with_timestamp}
!
! !INTERFACE:
subroutine create_output_filename_with_timestamp(&
     n, source, fname, yr,mo,da,hr,mn,ss,&
     model_name, writeint, &
     wout, style,odir)
! 
! !USES:
   use LVT_coreMod,  only : LVT_rc, LVT_LIS_rc
   use LVT_logMod,   only : LVT_log_msg, LVT_endrun

   implicit none 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:  
!  Create the file name for LIS output data files. It creates both the GSWP
!  style of output filenames and the standard LIS style. The convention used
!  in LIS creates a filename in the following default format (style==1) 
!
!   <output directory>/EXP<expno>/<model name>/<yr>/<yrmoda>/<yrmodahrmn>.<extension>
!  Style option ==2 corresponds to the following style: 
!  
!   <output directory>/EXP<expno>/<model name>.<yrmodahrmn>.<extension>
!  Style option ==3 :
!   <output directory>/EXP<expno>/<model name>.<extension>
! Style option ==4 :
!   <output directory>/<AFWA Weather product style>
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
!   \item [wout]
!    output writing format
!   \item [style]
!    style option as described above
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
!
! !ARGUMENTS:
   integer, intent(in) :: n
   integer, intent(in) :: source
   character(len=*), intent(out)          :: fname
   integer,  intent(in) :: yr, mo,da,hr,mn,ss
   character(len=*), intent(in), optional :: model_name ! needed for gswp run
   integer, intent(in), optional          :: writeint ! output writing interval
   character(len=*), intent(in), optional :: wout ! output format
   character(len=*), intent(in), optional :: style ! output directory style
   character(len=*), intent(in), optional :: odir

! 
!EOP
   character(len=8)        :: date
   character(len=10)       :: time
   character(len=5)        :: zone
   integer, dimension(8)   :: values
 
   character(len=10)       :: cdate
   character(len=12)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=200)       :: dname
   character(len=200), save :: out_fname
   character(len=50)        :: style_temp
   character(len=100)        :: odir_temp
   integer                  :: i, c

   if(.not.PRESENT(odir)) then 
      odir_temp = LVT_rc%odir
   else
      odir_temp = odir
   endif

   if(.not. PRESENT(style)) then 
      style_temp = "4 level hierarchy"
   else
      style_temp = style
   endif
   ! added by Shugong
   if(LVT_rc%lis_version == 6) then
     write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
          yr, mo, da, hr, mn
     
     dname = trim(LVT_rc%odir)//'/EXP'//trim(adjustl(LVT_rc%expcode))//'/'
     dname = trim(dname)//trim(LVT_rc%lsm)//'/'
     
     write(unit=cdate, fmt='(i4.4)') yr
     dname = trim(dname)//trim(cdate)//'/'
     
     write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') yr,mo,da
     dname = trim(dname)//trim(cdate)
     
     out_fname = trim(dname)//'/'//cdate1
     
     write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
     out_fname = trim(out_fname)//trim(cdate)
     out_fname = trim(out_fname)//'.gs4r'
     fname = trim(out_fname) 
  else ! LIS version is 7 or even higher number
     if(style_temp.eq."4 level hierarchy") then 
        write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
             yr,mo,da,hr,mn

        dname = trim(odir_temp)//'/'
        dname = trim(dname)//trim(model_name)//'/'

        write(unit=cdate, fmt='(i4.4)') yr
        dname = trim(dname)//trim(cdate)//'/'

        write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
             yr,mo,da
        dname = trim(dname)//trim(cdate)

        out_fname = trim(dname)//'/LIS_HIST_'//cdate1

        write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
        out_fname = trim(out_fname)//trim(cdate)

        if(present(wout)) then 
           select case ( wout )
           case ( "binary")
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif
           case ( "distributed binary")
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif                            
           case ( "grib1" )
              out_fname = trim(out_fname)//'.grb'
           case ( "netcdf")
              out_fname = trim(out_fname)//'.nc'
           case ( "ascii" )
              out_fname = trim(out_fname)//'.txt'
           case default
              call lvt_log_msg('ERR: create_output_filename -- '// &
                   'Unrecognized output format')
              call LVT_endrun 
           endselect
        endif
     elseif(style_temp.eq."3 level hierarchy") then 
        write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
             yr,mo,da,hr,mn

        dname = trim(odir_temp)//'/'
        dname = trim(dname)//trim(model_name)//'/'

        write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') yr,mo
        dname = trim(dname)//trim(cdate)

        out_fname = trim(dname)//'/LIS_HIST_'//cdate1

        write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
        out_fname = trim(out_fname)//trim(cdate)

        if(present(wout)) then 
           select case ( wout )
           case ( "binary" )
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif
           case ( "distributed binary")
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif                            
           case ( "grib1" )
              out_fname = trim(out_fname)//'.grb'
           case ( "netcdf" )
              out_fname = trim(out_fname)//'.nc'
           case ( "grib2" )
              out_fname = trim(out_fname)//'.gr2'
           case default
              call lvt_log_msg('ERR: create_output_filename -- '// &
                   'Unrecognized output format')
              call LVT_endrun 
           endselect
        endif
     elseif(style_temp.eq."2 level hierarchy") then 
        write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
             yr,mo,da,hr,mn

        dname = trim(odir_temp)//'/'
        dname = trim(dname)//trim(model_name)//'/'

        out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)

        write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
        out_fname = trim(out_fname)//trim(cdate)

        if(present(wout)) then 
           select case ( wout )
           case ( "binary")
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif
           case ( "distributed binary")
              if(LVT_LIS_rc(source)%wopt.eq."1d tilespace") then 
                 out_fname = trim(out_fname)//'.ts4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              elseif(LVT_LIS_rc(source)%wopt.eq."2d ensemble gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'                 
              elseif(LVT_LIS_rc(source)%wopt.eq."1d gridspace") then 
                 out_fname = trim(out_fname)//'.gs4r'
              endif                            
           case ( "grib1" )
              out_fname = trim(out_fname)//'.grb'
           case ( "netcdf" )
              out_fname = trim(out_fname)//'.nc'
           case ( "grib2" )
              out_fname = trim(out_fname)//'.gr2'
           case default
              call lvt_log_msg('ERR: create_output_filename -- '// &
                   'Unrecognized output format')
              call LVT_endrun 
           endselect
        endif
     elseif(style_temp.eq."WMO convention") then 
        write(unit=fint,fmt='(i2.2)') writeint/3600
        write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)')        &
             yr,mo,da

        write(unit=cdate, fmt='(i2.2, i2.2)') hr, mn

        if(LVT_rc%domain.eq."polar") then 
           fproj = 'P'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."lambert") then 
           fproj = 'L'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."polar") then 
           fproj = 'M'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."gaussian") then 
           fproj = 'G'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)*100        
           fres = '0P'//trim(fres)//'DEG'
        else
           fproj = 'C'
           write(unit=fres, fmt='(i10)') nint(LVT_LIS_rc(source)%gridDesc(10)*100)
           read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
           c = 0 
           do i=1,10
              if(fres1(i).ne.' '.and.c==0) c = i
           enddo
           fres3 = '0P'
           fres2 = trim(fres3)
           do i=c,10
              fres2 = trim(fres2)//trim(fres1(i))
           enddo
           fres2 = trim(fres2)//'DEG'
        endif

        dname = trim(odir_temp)//'/'//trim(cdate1)//&
!        dname = trim(odir_temp)//&
             '/PS.AFWA_SC.'//trim(LVT_rc%security_class)//&
             '_DI.'//trim(LVT_rc%distribution_class)//& 
             '_DC.'//&
             trim(LVT_rc%data_category)//'_GP.LIS_GR.'//&
             trim(fproj)//trim(fres2)//'_AR.'//trim(LVT_rc%area_of_data)//&
             '_PA.'//trim(fint)//'-HR-SUM_DD.'// &
             trim(cdate1)//'_DT.'//trim(cdate)//'_DF'
        select case (LVT_LIS_rc(source)%format)
        case ( "binary" )
           out_fname = trim(dname)//'.DAT'
        case ( "grib1" )
           out_fname = trim(dname)//'.GR1'
        case ( "netcdf" )
           out_fname = trim(dname)//'.nc'
        case ( "grib2" )
           out_fname = trim(dname)//'.GR2'
        case default            
        end select
     elseif(style_temp.eq."WMO convention (AFW OPS)") then 
        write(unit=fint,fmt='(i2.2)') writeint/3600
        write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
             yr, mo, da

        write(unit=cdate, fmt='(i2.2, i2.2)') hr, mn

        if(LVT_rc%domain.eq."polar") then 
           fproj = 'P'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."lambert") then 
           fproj = 'L'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."polar") then 
           fproj = 'M'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)
           fres = trim(fres)//'KM'
        elseif(LVT_rc%domain.eq."gaussian") then 
           fproj = 'G'
           write(unit=fres, fmt='(i2.2)') LVT_LIS_rc(source)%gridDesc(9)*100        
           fres = '0P'//trim(fres)//'DEG'
        else
           fproj = 'C'
           write(unit=fres, fmt='(i10)') nint(LVT_LIS_rc(source)%gridDesc(10)*100)
           read(unit=fres,fmt='(10a1)') (fres1(i),i=1,10)
           c = 0 
           do i=1,10
              if(fres1(i).ne.' '.and.c==0) c = i
           enddo
           fres3 = '0P'
           fres2 = trim(fres3)
           do i=c,10
              fres2 = trim(fres2)//trim(fres1(i))
           enddo
           fres2 = trim(fres2)//'DEG'
        endif

        dname = trim(odir_temp)//'/'//trim(cdate1)//&
!        dname = trim(odir_temp)//&
             '/PS.AFWA_SC.'//trim(LVT_rc%security_class)//&
             '_DI.'//trim(LVT_rc%distribution_class)//& 
             '_GP.LIS_GR.'//&
             trim(fproj)//trim(fres2)//'_AR.'//trim(LVT_rc%area_of_data)//&
             '_PA.LIS_DD.'// &
             trim(cdate1)//'_DT.'//trim(cdate)//'_DF'
        select case (LVT_LIS_rc(source)%format)
        case ( "binary" )
           out_fname = trim(dname)//'.DAT'
        case ( "grib1" )
           out_fname = trim(dname)//'.GR1'
        case ( "netcdf" )
           out_fname = trim(dname)//'.nc'
        case ( "grib2" )
           out_fname = trim(dname)//'.GR2'
        case default            
        end select
     else
        call lvt_log_msg('ERR: LVT_create_output_filename --')
        call lvt_log_msg('  Unrecognized LIS output naming style:')
        call lvt_log_msg('      '//style_temp)
        call LVT_endrun 
     endif
     fname = trim(out_fname)
   endif 
 end subroutine Create_output_filename_with_timestamp

!BOP
! 
! !ROUTINE: LVT_create_daobs_filename
! \label{LVT_create_daobs_filename}
!
! !INTERFACE:
subroutine LVT_create_daobs_filename(n, fname)
! 
! !USES:
   use LVT_coreMod,  only : LVT_rc
  
   implicit none 
!
! !INPUT PARAMETERS: 
   integer, intent(in)           :: n
! 
! !OUTPUT PARAMETERS:
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
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

   character(len=200) :: out_fname
   character*100      :: cdate, cdate1

   
   write(unit=cdate, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
        LVT_rc%yr, LVT_rc%mo, &
        LVT_rc%da, LVT_rc%hr,LVT_rc%mn
   write(unit=cdate1, fmt='(i4.4, i2.2)') &
        LVT_rc%yr, LVT_rc%mo
   
   out_fname = trim(LVT_rc%odir)//'/' & 
        //'/DAOBS/'//trim(cdate1)//'/'//trim(cdate)//'.1gs4r'
   
   fname = trim(out_fname)
   
 end subroutine LVT_create_daobs_filename

!BOP
! 
! !ROUTINE: convertParam_real
! \label{convertParam_real}
!
! !INTERFACE: 
subroutine convertParam_real(k,pdata, ldata)
! !USES: 
  use LVT_coreMod
  use map_utils

! !ARGUMENTS: 
  integer,     intent(in)    :: k
  real,        intent(in)    :: pdata(LVT_rc%pnc,LVT_rc%pnr)
  real,        intent(inout) :: ldata(LVT_LIS_rc(k)%lnc,LVT_LIS_rc(k)%lnr)
!
! !DESCRIPTION:
!  This routine subsets the floating point data 
!  read from the parameter attributes
!  file supplied by LDT to the LVT domain being run. It is assumed
!  that the parameter data domain is at the same spatial resolution 
!  and map projection as that of the LVT domain. 
!  
!EOP  
  real                :: rlat, rlon,ctmp, rtmp
  integer             :: c,r

! Using the map_utils routines could lead to roundoff errors. We use
! simpler calculations for latlon projection. 

  do r=1,LVT_LIS_rc(k)%lnr
     do c=1,LVT_LIS_rc(k)%lnc

        call ij_to_latlon(LVT_LIS_domain(k)%proj, float(c), float(r), &
             rlat, rlon)
        
        call latlon_to_ij(LVT_domain%lvtparamproj, rlat, rlon, &
             ctmp, rtmp)
        ldata(c,r) = pdata(nint(ctmp), nint(rtmp))
        
     enddo
  enddo
   
end subroutine convertParam_real

!BOP
! 
! !ROUTINE: convertParam_int
! \label{convertParam_int}
!
! !INTERFACE: 
subroutine convertParam_int(k, pdata, ldata)
! !USES: 
  use LVT_coreMod
  use map_utils

! !ARGUMENTS: 
  integer,     intent(in)    :: k
  integer,     intent(in)    :: pdata(LVT_rc%pnc,LVT_rc%pnr)
  integer,     intent(inout) :: ldata(LVT_LIS_rc(k)%lnc,LVT_LIS_rc(k)%lnr)
!
! !DESCRIPTION:
!  This routine subsets the floating point data 
!  read from the parameter attributes
!  file supplied by LDT to the LVT domain being run. It is assumed
!  that the parameter data domain is at the same spatial resolution 
!  and map projection as that of the LVT domain. 
!  
!EOP  
  real                :: rlat, rlon,ctmp, rtmp
  integer             :: c,r

  do r=1,LVT_LIS_rc(k)%lnr
     do c=1,LVT_LIS_rc(k)%lnc
        call ij_to_latlon(LVT_LIS_domain(k)%proj, float(c), float(r), &
             rlat, rlon)
        
        call latlon_to_ij(LVT_domain%lvtparamproj, rlat, rlon, &
             ctmp, rtmp)
        ldata(c,r) = pdata(nint(ctmp), nint(rtmp))
        
     enddo
  enddo
   
end subroutine convertParam_int


end module LVT_fileIOMod
