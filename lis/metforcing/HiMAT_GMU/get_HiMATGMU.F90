!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_HiMATGMU
! \label{get_HiMATGMU}
!
! !REVISION HISTORY:
! 28 July 2017: Sujay Kumar;  Data and Code implementation
!
! !INTERFACE:
subroutine get_HiMATGMU(n, findex)

! !USES:
  use LIS_coreMod, only     : LIS_rc, LIS_domain
  use LIS_timeMgrMod, only  : LIS_tick, LIS_get_nstep
  use LIS_logMod,      only : LIS_logunit, LIS_endrun
  use HiMATGMU_forcingMod, only : HiMATGMU_struc
  use LIS_constantsMod,    only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n 
  integer, intent(in) :: findex

! !DESCRIPTION:
!  Opens, reads, and interpolates hourly HiMAT GMU forcing. 
!  At the beginning of a simulation, the code reads the most
!  recent past data (nearest the hour interval), and the nearest
!  future data. These two datasets are used to temporally 
!  interpolate the data to the current model timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[lis\_time]
!    Current LIS Time 
!  \item[HiMATGMU\_file\_time]
!    End boundary time of HiMAT GMU file 
!  \item[file\_name]
!    HiMAT GMU filename - passed back to getstg 
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the HiMAT GMU data times
!  \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    Computes the neighbor, weights for bilinear interpolation
!  \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    Computes the neighbor, weights for conservative interpolation
!  \item[HiMATGMUfile](\ref{HiMATGMUfile}) \newline
!    Puts together appropriate file name for 1-hour intervals
!  \item[read\_HiMATGMU](\ref{read_HiMATGMU}) \newline
!    Interpolates HiMAT GMU data to LIS grid
!  \end{description}
!EOP
   
!== Local Variables =======================
    integer :: ferror_HiMATGMU           ! Error flags for precip data sources
!    integer :: endtime_HiMATGMU         ! 1=get a new file 
    integer :: order

    real*8  :: HiMATGMU_file_time1  
    real*8  :: HiMATGMU_file_time2  
    character(len=LIS_CONST_PATH_LEN) :: file_name      

    integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
    integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
    real    :: gmt1, gmt2,ts1,ts2                    
    real    :: gridDesci(LIS_rc%nnest,50)

!=== End Variable Definition =======================

!    endtime_HiMATGMU = 0

!-- Determine LIS's current time and the time of the data file:
    yr1 = LIS_rc%yr
    mo1 = LIS_rc%mo
    da1 = LIS_rc%da
    hr1 = LIS_rc%hr    
    mn1 = 0
    ss1 = 0
    ts1 = 0
    call LIS_tick( HiMATGMU_file_time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

!-- data product time; end accumulation time data
    yr2 = LIS_rc%yr     
    mo2 = LIS_rc%mo
    da2 = LIS_rc%da
    hr2 = LIS_rc%hr+1  
    mn2 = 0
    ss2 = 0
    ts2 = 0
    call LIS_tick( HiMATGMU_file_time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

!-- Ensure that data is found during first time step
    if ( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1) then 
!         endtime_HiMATGMU = 1
         LIS_rc%rstflag(n) = 0
    endif

    ferror_HiMATGMU = 0
    order = 2

    if( LIS_rc%ts < HiMATGMU_struc(n)%ts ) then
      if ( LIS_rc%time > HiMATGMU_struc(n)%HiMATGMUtime ) then

        call HiMATGMUfile( file_name, HiMATGMU_struc(n)%HiMATGMUdir, yr2, mo2, da2, hr2 )
        write(LIS_logunit,*) '[INFO] Getting new HiMAT GMU precip data: ', trim(file_name)
        call read_HiMATGMU ( n, file_name, findex, order, ferror_HiMATGMU )
        HiMATGMU_struc(n)%HiMATGMUtime = HiMATGMU_file_time2
      endif

    elseif( LIS_rc%ts == HiMATGMU_struc(n)%ts ) then

       call HiMATGMUfile( file_name, HiMATGMU_struc(n)%HiMATGMUdir, yr1, mo1, da1, hr1 )
       write(LIS_logunit,*) '[INFO] Getting new HiMAT GMU precip data: ', trim(file_name)
       call read_HiMATGMU ( n, file_name, findex, order, ferror_HiMATGMU )
       HiMATGMU_struc(n)%HiMATGMUtime = HiMATGMU_file_time1

    else
      write(LIS_logunit,*) "[ERR] ERR MSG: HiMAT GMU READER CANNOT HANDLE LIS "
      write(LIS_logunit,*) "[ERR]   TIMESTEP > 3600 secs -- AT THIS TIME.  STOPPING ..."
      call LIS_endrun
    endif

end subroutine get_HiMATGMU

