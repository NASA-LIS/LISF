!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module LIS_forecastMod
!BOP
!
! !MODULE: LIS_forecastMod
!
! !DESCRIPTION:
! 
!  \subsubsection{Overview}
!  This module provides routines for reading and modifying forecasting data.
!
! !REVISION HISTORY:
!
!  19 Feb 2016: Sujay Kumar; Initial implementation
!
  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_numerRecipesMod
  use LIS_ran2_gasdev

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_forecast_init  ! initializes data structures and memory
  public :: LIS_sample_forecastDate
  public :: LIS_forecast_writerestart
  public :: LIS_forecast_readrestart
  public :: LIS_get_iteration_index
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_forecast_struc
!EOP
  type, public :: forecast_type_dec
     character*50            :: alg   
     integer                 :: nIterations
     integer                 :: st_iterId
     integer                 :: iterId
     character*50            :: startMode
     character(len=LIS_CONST_PATH_LEN) :: rstfilename
     integer, allocatable    :: seed(:,:)

  end type forecast_type_dec
  
  type(forecast_type_dec), allocatable :: LIS_forecast_struc(:)
!EOP

contains

!BOP
! 
! !ROUTINE: LIS_forecast_init
! \label{LIS_forecast_init}
! 
! !INTERFACE:
  subroutine LIS_forecast_init()
! !USES:
    use LIS_coreMod
! 
! !DESCRIPTION:
!
! Allocates memory for data structures for supporting forecasting
! instances.
!
!EOP
    implicit none
    integer                  :: n, i,k
    integer                  :: rc
    
    allocate(LIS_forecast_struc(LIS_rc%nnest))
    LIS_forecast_struc(:)%st_iterId  =  1

    call ESMF_ConfigFindLabel(LIS_config,"Forecast forcing source mode:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,LIS_forecast_struc(n)%alg,rc=rc)
       call LIS_verify(rc,"Forecast forcing source mode: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Forecast forcing start mode:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,LIS_forecast_struc(n)%startMode,rc=rc)
       call LIS_verify(rc,"Forecast forcing start mode: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Forecast forcing restart filename:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,LIS_forecast_struc(n)%rstfilename,rc=rc)
       call LIS_verify(rc,"Forecast forcing restart filename: not defined")
    enddo

    do n=1,LIS_rc%nnest
       call forecastalginit(trim(LIS_forecast_struc(n)%alg)//char(0))
    enddo

    LIS_rc%forecastMode = 1
    
  end subroutine LIS_forecast_init

!BOP
!
! !ROUTINE: LIS_sample_forecastDate
! \label{LIS_sample_forecastDate}
!
! !INTERFACE:
  subroutine LIS_sample_forecastDate(n,kk, k, yr,mo,da)
! !ARGUMENTS
    integer                 :: n 
    integer                 :: kk
    integer                 :: k
    integer                 :: yr, mo, da
! 
! !DESCRIPTION: 
!  This subroutine returns a "sampled" date from the 
!  forcing archive. Based on the input date, the routine
!  first randomly samples the year from the archive period.  
!  It then samples a day that is within +/- user defined 
!  time window length in that year. This random sampling
!  is done once in every 10 days. At other times, the code
!  simply returns consecutive dates using the randomly
!  sampled date as the temporary reference. 
!
!EOP

    
    call forecastalgsampledate(trim(LIS_forecast_struc(n)%alg)//char(0), &
         n, kk,k, yr, mo, da)
    
  end subroutine LIS_sample_forecastDate

  subroutine LIS_forecast_writerestart

    integer         :: n 
    integer         :: ios
    integer         :: ftn    
    character(len=LIS_CONST_PATH_LEN) :: out_dname
    character(len=LIS_CONST_PATH_LEN) :: rst_fname
    character*4     :: fiter
    integer, external :: LIS_create_subdirs 
    character(len=201) :: c_string  

    n = 1
    
    if(LIS_masterproc) then 
       
       out_dname = trim(LIS_rc%odir)//'/FCST/'
#if ( defined AIX )
       call system('mkdir -p '//trim(out_dname),ios)
#else
       c_string = trim(out_dname)
       ios = LIS_create_subdirs(len_trim(c_string),trim(c_string))
       if (ios .ne. 0) then
          write(LIS_logunit,*)'[ERR] problem creating directory ', &
               trim(c_string)
          flush(LIS_logunit)
       end if
#endif
       rst_fname = trim(LIS_rc%odir)//'/FCST/LIS_RST_FCST_'
       write(unit=fiter,fmt='(i4.4)')  LIS_forecast_struc(n)%iterId
       rst_fname = trim(rst_fname)//trim(fiter)//'.bin'
       
       ftn = LIS_getNextUnitNumber()
       
       open(ftn,file=rst_fname,status='unknown',form='unformatted')
       write(ftn) LIS_forecast_struc(n)%iterId
       write(ftn) LIS_forecast_struc(n)%seed       
       call LIS_releaseUnitNumber(ftn)
       
    endif
    
  end subroutine LIS_forecast_writerestart

  subroutine LIS_forecast_readrestart

    integer         :: n 
    integer         :: ftn    
    integer         :: iterid
    character(len=LIS_CONST_PATH_LEN) :: out_dname
    character(len=LIS_CONST_PATH_LEN) :: rst_fname    

    n = 1
    rst_fname = LIS_forecast_struc(n)%rstfilename
    
    ftn = LIS_getNextUnitNumber()
    
    open(ftn,file=rst_fname,status='unknown',form='unformatted')
    read(ftn) iterId
    LIS_forecast_struc(n)%st_iterId = iterid + 1
    read(ftn) LIS_forecast_struc(n)%seed
    call LIS_releaseUnitNumber(ftn)

  end subroutine LIS_forecast_readrestart
  
!BOP
!
! !ROUTINE: LIS_get_iteration_index
! \label{LIS_get_iteration_index}
!
! !INTERFACE:
integer function LIS_get_iteration_index(n, k, index1, mfactor) result(kk)
! !USES:
! None

   implicit none

! !ARGUMENTS:
   integer, intent(in)  :: n, k, index1, mfactor
!
! !DESCRIPTION:
!  Returns the iteration index used to offset into the metdata arrays
!  when running in forecase mode; else returns 1.
!
! Note that the arithmetic assumes that sub-grid tiling is disabled
! when running in forecast mode.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!     nest index
!  \item[k]
!     tile index
!  \item[index1]
!     grid-cell index corresponding to k
!  \item[mfactor]
!     nensem/nIter
!  \item[kk]
!     iteration index
!  \end{description}
!
!EOP

   if ( LIS_rc%forecastMode == 1 ) then
      kk = k-(index1-1)*LIS_rc%nensem(n)/mfactor
   else
      kk = 1
   endif

end function LIS_get_iteration_index

end module LIS_forecastMod
