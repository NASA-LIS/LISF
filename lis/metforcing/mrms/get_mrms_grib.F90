!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_mrms_grib
! \label{get_mrms_grib}
!
! !REVISION HISTORY:
! 25 May 2006: Kristi Arsenault;  Data and Code implementation
!            : Jon Case: Written for MRMS binary data (pre-operational)
! 5 September 2017: Jessica Erlingis; update for operational GRIB2 format
! 22 February 2019: Jessica Erlingis; updated to include masking option
!
! !INTERFACE:
subroutine get_mrms_grib(n, findex)

! !USES:
  use LIS_coreMod, only     : LIS_rc
  use LIS_logMod, only      : LIS_logunit
  use LIS_timeMgrMod, only  : LIS_tick, LIS_get_nstep
  use mrms_grib_forcingMod, only : mrms_grib_struc

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n 
  integer, intent(in) :: findex

! !DESCRIPTION:
!  Opens, reads, and interpolates hourly MRMS GRIB2 forcing. 
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
!  \item[mrms\_file\_time]
!    End boundary time of MRMS file 
!  \item[file\_name]
!    MRMS filename - passed back to getmrms_grib 
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick})\\
!    determines the MRMS data times
!  \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \\
!    Computes the neighbor, weights for bilinear interpolation
!  \item[conserv\_interp\_input](\ref{conserv_interp_input}) \\
!    Computes the neighbor, weights for conservative interpolation
!  \item[mrmsgribfile](\ref{mrmsgribfile})\\
!    Puts together appropriate file name for 1-hour intervals
!  \item[read\_mrms_grib](\ref{read_mrms_grib})\\
!    Interpolates MRMS data to LIS grid
!  \end{description}
!EOP
   
!== Local Variables =======================
  integer :: ferror_mrms_grib              ! Error flags for precip data sources
  integer :: endtime_mrms_grib             ! 1=get a new file 
  integer :: order

  real*8  :: lis_time, mrms_grib_file_time ! Current LIS Time and end boundary time for MRMS file
  character(150) :: file_name               ! Filename variables for precip data sources

  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  real    :: gmt1, gmt2,ts1,ts2                    
                                             !  for precip data sources
  real    :: gridDesci(LIS_rc%nnest,50)
! J.Case -- local variables
  character*50 :: dom
  real         :: res

!=== End Variable Definition =======================

  endtime_mrms_grib = 0

!-- Determine LIS's current time and the time of the MRMS file:
! - LIS's current date and time
  yr1 = LIS_rc%yr  
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0      ! Time increment in seconds

  call LIS_tick( lis_time, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )   

!-- MRMS product time; end accumulation time data
  yr2 = LIS_rc%yr     
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = LIS_rc%hr+1  ! Advance MRMS by one hour increment
  mn2 = 0
  ss2 = 0
  ts2 = 0

  call LIS_tick( mrms_grib_file_time, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

!-- Determine if MRMS Grid LIS_domain has changed (before opening file)   

    mrms_grib_struc(n)%ncol = 7000
    mrms_grib_struc(n)%nrow = 3500

!   -- Reinitialize the weights and neighbors
    gridDesci(n,:) = 0.0
    gridDesci(n,1) = 0                       ! Projection type (Lat/Lon)
    gridDesci(n,2) = mrms_grib_struc(n)%ncol ! X-dir amount of points
    gridDesci(n,3) = mrms_grib_struc(n)%nrow ! y-dir amount of points
    gridDesci(n,4) =   54.995                ! Starting latitude point
    gridDesci(n,5) = -129.995                ! Starting longitude point
    gridDesci(n,6) = 128                     ! (not used)
    gridDesci(n,7) =   20.005                ! Ending latitude point
    gridDesci(n,8) =  -60.005                ! Ending longitude point 
    gridDesci(n,9) =  0.01                   ! spatial resolution in W-E dirn (deg)
    gridDesci(n,10) = 0.01                   ! spatial resolution in S-N dirn (deg)
    gridDesci(n,20) = 64                     ! N-S ordering (number divisible by 32; same as in NLDAS2)

    dom=LIS_rc%lis_map_proj
    res=LIS_rc%gridDesc(n,9)

! J.Case (8/22/2012) -- Use interpolation only if LIS resolution is < MRMS resolution.
!                       Otherwise, use upscaling by default.
!
! J. Erlingis (9/26/2018) -- Add neighbor option if resolutions match. Fix arguments for bilinear
!                            and budget bilinear

    if ( (((dom.eq."latlon").or.(dom.eq."gaussian")) .and. (res.eq.0.01)) .or. & !Case where resolutions match
      (((dom.eq."mercator").or.(dom.eq."lambert").or.                       &
      (dom.eq."polar").or.(dom.eq."UTM")) .and. (res.eq.1.0)) ) then

! === NEAREST NEIGHBOR INTERPOLATION ===

      call neighbor_interp_input(n, gridDesci(n,:), &
        mrms_grib_struc(n)%n113)

    elseif ( (((dom.eq."latlon").or.(dom.eq."gaussian")) .and. (res.lt.0.01)) .or. &
      (((dom.eq."mercator").or.(dom.eq."lambert").or.                       &
      (dom.eq."polar").or.(dom.eq."UTM")) .and. (res.lt.1.0)) ) then

      if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then
! === BILINEAR INTERPOLATION ==== 
        call bilinear_interp_input(n, gridDesci(n,:),&
         mrms_grib_struc(n)%n111,mrms_grib_struc(n)%n121,&
          mrms_grib_struc(n)%n211,mrms_grib_struc(n)%n221,&
          mrms_grib_struc(n)%w111,mrms_grib_struc(n)%w121,&
          mrms_grib_struc(n)%w211,mrms_grib_struc(n)%w221)


      elseif ( trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") then 
! === BUDGET BILINEAR INTERPOLATION ==== 
        call conserv_interp_input(n, gridDesci(n,:),&
          mrms_grib_struc(n)%n112,mrms_grib_struc(n)%n122,&
          mrms_grib_struc(n)%n212,mrms_grib_struc(n)%n222,&
          mrms_grib_struc(n)%w112,mrms_grib_struc(n)%w122,&
          mrms_grib_struc(n)%w212,mrms_grib_struc(n)%w222)

      endif

      else !! use upscaling
! === UPSCALING ==== 
        call upscaleByAveraging_input( gridDesci(n,:), LIS_rc%gridDesc(n,:), &
          mrms_grib_struc(n)%ncol*mrms_grib_struc(n)%nrow, &
          LIS_rc%lnc(n)*LIS_rc%lnr(n), &
          mrms_grib_struc(n)%n11)
              
      endif ! dom projection check for interpolation


!-- Ensure that data are found during first time step
    if ( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1) then 
      endtime_mrms_grib = 1
      LIS_rc%rstflag(n) = 0
    endif

!-- Check for and get MRMS Precipitation data
    if ( lis_time > mrms_grib_struc(n)%mrms_grib_time )  endtime_mrms_grib = 1

!-- Get new second file (time2 data)
    if ( endtime_mrms_grib == 1 ) then  

!   -- Determine and return filename of MRMS file 
      call mrms_gribfile( file_name, mrms_grib_struc(n)%mrms_grib_dir, yr2, mo2, da2, hr2 )

      write(LIS_logunit,*) 'Getting new MRMS precip data:: ', trim(file_name)
      ferror_mrms_grib = 0
!       print *, "FERROR (get_mrms) :: ", ferror_mrms
      order = 2
      call read_mrms_grib( n, file_name, findex, order, yr2, mo2, da2, ferror_mrms_grib )

!   -- Assign latest MRMS file time to stored MRMS time variable
      mrms_grib_struc(n)%mrms_grib_time = mrms_grib_file_time

    endif  

  return

end subroutine get_mrms_grib
