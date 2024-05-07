!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module micromet_module
!BOP
!
! !MODULE: micromet_module
!
! !DESCRIPTION:
!  This module contains the key MicroMet routines for each
!  of the different meteorological forcing fields (e.g.,
!  temperature, etc.).
!
!  The routines in this module include: 
!
!  \begin{description}
!   \item[micromet_relhum]
!    micromet_relhum  ...
!   \end{description}
!
! !REVISION HISTORY:
!  01 Mar 2006: G. Liston and K. Elder; MicroMet routine authors
!  16 Feb 2021: K. Arsenault; Added MicroMet routines to LIS
!  03 Feb 2022: K. Arsenault; Added MicroMet's solar routines
! 
!EOP
  implicit none

contains

!BOP
!
! !ROUTINE: micromet_relhum
! \label{micromet_relhum}
!
! !DESCRIPTION:  
!  This routine calculates relative humidity (rh; %) for given
!  specific humidity (sphm; kg/kg) and pressure (pres; Pa)
!
! !INTERFACE:
  subroutine micromet_relhum(nx, ny, tair, sphm, pres, rh)

    implicit none

    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real, intent(in)  :: tair(nx,ny)
    real, intent(in)  :: sphm(nx,ny)
    real, intent(in)  :: pres(nx,ny)
    real, intent(out) :: rh(nx,ny)

    integer :: i, j
    real :: Tf       ! Temperature at freezing point
    real :: e        ! vapor pressure 
    real :: es       ! saturation vapor pressure
    real :: A, B, C  ! Coeffs for sat vapor pressure over water

! ______________________________________________________________

   Tf = 273.16   

! Coeffs for saturation vapor pressure over water (Buck 1981).
!   Note: temperatures for Buck`s equations are in deg C, and
!   vapor pressures are in mb.  Do the adjustments so that the
!   calculations are done with temperatures in K, and vapor
!   pressures in Pa.
! Over water.
    A = 6.1121 * 100.0
    B = 17.502
    C = 240.97
! Over ice.
!    A = 6.1115 * 100.0
!    B = 22.452
!    C = 272.55

  do j = 1, ny
   do i = 1, nx

! Compute the saturation water vapor pressure for given air temperature.
!    CALL SATVAPOR(es,tair,Tf)
     es = A * exp((B * (Tair(i,j) - Tf))/(C + (Tair(i,j) - Tf)))

! Compute the water vapor pressure from the specific humidity (kg/kg)
!   and pressure (Pa).
     e = sphm(i,j) * pres(i,j) / (0.622 + 0.378 * sphm(i,j))

! Relative humidity (%).
     rh(i,j) = 100.0 * (e / es)
     rh(i,j) = min(100.0,rh(i,j))
     rh(i,j) = max(0.0,rh(i,j))
   
   enddo
  enddo

  end subroutine micromet_relhum

!BOP
!
! !ROUTINE: micromet_snowfall
! \label{micromet_snowfall}
!
! !DESCRIPTION:  
!  This routine discriminates the amount of snowfall
!   based on total precipitation, 2m air temperature
!   and the snow fraction method based on the option
!   selected in snowmodel.par
!
! !INTERFACE:
  subroutine micromet_snowfall(nx, ny, snowfall_frac, prec, tair, sprec)

    implicit none

    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real, intent(in)  :: snowfall_frac  ! Snowfall method option
    real, intent(in)  :: prec(nx,ny)    ! Total precip
    real, intent(in)  :: tair(nx,ny)    ! 2m Air temp
    real, intent(out) :: sprec(nx,ny)   ! Snowfall

    real :: Tair_C(nx,ny)   ! Temperature in Celsius
    real :: Tf              ! Temperature at freezing point

    integer :: i,j
    real :: snowfall_frac_1, snowfall_frac_2, snowfall_frac_3
    real :: Tair_C_center
    real :: slope, b
! ______________________________________________________________

   Tf = 273.15

   do j = 1, ny
    do i = 1, nx
 
     Tair_C(i,j) = Tair(i,j) - Tf 

! Option 1: Auer (1974); a rain-snow threshold at +2.0 C.
     if (snowfall_frac.eq.1.0) then

       if (Tair_C(i,j).lt.2.0) then
         snowfall_frac_1 = 1.0
        else
         snowfall_frac_1 = 0.0
       endif
       sprec(i,j) = snowfall_frac_1 * prec(i,j)

! Option 2: Dai, A. (2008): Temperature and pressure dependenxe of the
!   rain-snow phase transition over land and ocean, Geophys.
!   Res. Lett., 35, L12802, doi:10.1029/2008GL033295.
! In this implementation, G. Liston clipped Dai's values to 
!   0 and 1.
     elseif (snowfall_frac.eq.2.0) then
       snowfall_frac_2 = - 0.482292 * &
              (tanh(0.7205 * (Tair_C(i,j) - 1.1662)) - 1.0223)
       if (Tair_C(i,j).lt.-4.0) then
         snowfall_frac_2 = 1.0
       elseif (Tair_C(i,j).gt.6.0) then
         snowfall_frac_2 = 0.0
       endif
       sprec(i,j) = snowfall_frac_2 * prec(i,j)

! Option 3: Glen's linear approximation to Dai (2008). This plots right
!   over the top of Dai, between frac = 0.1 and 0.9.
     elseif (snowfall_frac.eq.3.0) then

! Define where you want the center temperature to be when frac = 0.5.
       Tair_C_center = 1.1
! Define the slope of the line.
       slope = -0.30
! Calculate the intercept (the 0.5 is the frac for Tair_C_center).
       b = 0.5 - slope * Tair_C_center

! Solve the equation in the form y = m*x + b
       snowfall_frac_3 = slope * Tair_C(i,j) + b
       snowfall_frac_3 = max(0.0,snowfall_frac_3)
       snowfall_frac_3 = min(1.0,snowfall_frac_3)
       sprec(i,j) = snowfall_frac_3 * prec(i,j)

     endif

    enddo
   enddo

 end subroutine micromet_snowfall


!BOP
!
! !ROUTINE: micromet_wind
! \label{micromet_wind}
!
! !DESCRIPTION:  
!  This routine calculates the wind distribution, for both
!   U- and V-wind gridded fields, and the gridded wind speed
!   and wind direction fields.
!  Options selected in snowmodel.par
!
! !INTERFACE:
  subroutine micromet_wind( nx, ny, deltax, deltay, &
                         uwind_grid,vwind_grid,slopewt,curvewt,& 
                         curvature,slope_az,terrain_slope,windspd_grid,& 
                         winddir_grid,windspd_flag,winddir_flag,windspd_min,&
                         vegtype,forest_LAI,calc_subcanopy_met,vegsnowd_xy,&
                         topo,wind_lapse_rate,curve_len_scale)

    use LIS_mpiMod
    use snowmodel_lsmMod, only : snowmodel_struc

    implicit none

    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real,    intent(in) :: deltax      ! grid increment in x
    real,    intent(in) :: deltay      ! grid increment in y

    ! Output wind fields
    real,    intent(inout) :: uwind_grid(nx,ny)    ! output U-wind values
    real,    intent(inout) :: vwind_grid(nx,ny)    ! output V-wind values
    real,    intent(inout) :: winddir_grid(nx,ny)  ! output wind direction
    real,    intent(inout) :: windspd_grid(nx,ny)  ! output wind speed
    real,    intent(inout) :: windspd_flag
    real,    intent(inout) :: winddir_flag

    ! Input parameters
    real,    intent(in) :: curvature(nx,ny)
    real,    intent(in) :: slope_az(nx,ny)
    real,    intent(in) :: terrain_slope(nx,ny)
    real,    intent(in) :: vegtype(nx,ny)
    real,    intent(in) :: vegsnowd_xy(nx,ny)
    real,    intent(in) :: topo(nx,ny)

    ! User-specified inputs:
    real,    intent(in) :: slopewt     ! wind model slope weight
    real,    intent(in) :: curvewt     ! wind model curvature weight
    real,    intent(in) :: windspd_min ! Min wind speed threshold (set by user)
    real,    intent(in) :: curve_len_scale   ! length scale for curvature calculation
    real,    intent(in) :: wind_lapse_rate
    real,    intent(in) :: calc_subcanopy_met
    real,    intent(in) :: forest_LAI(5)       
   
    integer :: i, j
    integer :: ierr
    integer :: loops_windwt_smoother
    real    :: topo_ref_grid(nx,ny)
    real    :: delta_topo,alfa1,alfa2
    real    :: u_sum, v_sum
    real    :: pi,deg2rad,rad2deg

! Initialize any local fields:
    topo_ref_grid = 0.0

! Define the required constants.
    pi = 2.0 * acos(0.0)
    deg2rad = pi / 180.0
    rad2deg = 180.0 / pi

! If desired, impose a wind speed increase with elevation.  Here
!   the wind_lapse_rate = the wind speed increase per 1-km elevation
!   gain.  The adjustment is patterned after the precipitation-
!   elevation adjustment.  
!   Note: topo_ref_grid here comes from the precipitation adjustment.
    if (wind_lapse_rate.ne.0.0) then
      alfa1 = (wind_lapse_rate - 1.0) / (1.0 + wind_lapse_rate)
      ! Convert to m-1.
      alfa1 = alfa1 / 1000.0
      do j=1,ny
        do i=1,nx
          delta_topo = topo(i,j) - topo_ref_grid(i,j)
          ! Impose some limits to the adjustment.
          delta_topo = min(delta_topo,1800.0)
          alfa2 = alfa1 * delta_topo
          uwind_grid(i,j) = uwind_grid(i,j) * (1.0 + alfa2)/(1.0 - alfa2)
          vwind_grid(i,j) = vwind_grid(i,j) * (1.0 + alfa2)/(1.0 - alfa2)
        enddo
      enddo
    endif

! Convert these u and v components to speed and directions.
    do j=1,ny
      do i=1,nx
! Some compilers do not allow both u and v to be 0.0 in
!   the atan2 computation.
        if (abs(uwind_grid(i,j)).lt.1e-10) then
           uwind_grid(i,j) = 1e-10
        endif
        winddir_grid(i,j) = rad2deg * atan2(uwind_grid(i,j),vwind_grid(i,j))
        if (winddir_grid(i,j).ge.180.0) then
          winddir_grid(i,j) = winddir_grid(i,j) - 180.0
        else
          winddir_grid(i,j) = winddir_grid(i,j) + 180.0
        endif

        windspd_grid(i,j) = sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
      enddo
    enddo

! Modify the wind speed and direction according to simple
!   wind-topography relationships.
    call topo_mod_winds(nx,ny,winddir_grid,slopewt,curvewt,&
     &  windspd_grid,uwind_grid,vwind_grid,curvature,slope_az,&
     &  terrain_slope,vegtype,forest_LAI,calc_subcanopy_met,&
     &  vegsnowd_xy,curve_len_scale,deltax,deltay)

! Avoid problems of zero (low) winds (for example, turbulence
!   theory, log wind profile, etc., says that we must have some
!   wind.  Thus, some equations blow up when the wind speed gets
!   very small).
    do j=1,ny
      do i=1,nx
        if (windspd_grid(i,j).lt.windspd_min) then
          windspd_grid(i,j) = windspd_min
          uwind_grid(i,j) = (- windspd_grid(i,j)) * &
     &        sin(deg2rad*winddir_grid(i,j))
          vwind_grid(i,j) = (- windspd_grid(i,j)) * &
     &        cos(deg2rad*winddir_grid(i,j))
        endif
      enddo
    enddo

! Find the maximum wind speed in the domain, and the
!   domain-averaged wind direction.
    windspd_flag = 0.0
    u_sum = 0.0
    v_sum = 0.0
    do j=1,ny
      do i=1,nx
        windspd_flag = max(windspd_flag,windspd_grid(i,j))
        u_sum = u_sum + uwind_grid(i,j)
        v_sum = v_sum + vwind_grid(i,j)
      enddo
    enddo

#if (defined SPMD)
     call MPI_Barrier(LIS_MPI_COMM, ierr)
     call MPI_ALLREDUCE(u_sum, snowmodel_struc(1)%usum_glb, 1,&
          MPI_REAL, MPI_SUM,&
          LIS_mpi_comm, ierr)
     call MPI_Barrier(LIS_MPI_COMM, ierr)
     call MPI_ALLREDUCE(v_sum, snowmodel_struc(1)%vsum_glb, 1,&
          MPI_REAL, MPI_SUM,&
          LIS_mpi_comm, ierr)
     call MPI_Barrier(LIS_MPI_COMM, ierr)
     call MPI_ALLREDUCE(windspd_flag, snowmodel_struc(1)%windspdflg_glb, 1,&
          MPI_REAL, MPI_MAX,&
          LIS_mpi_comm, ierr)

     u_sum = snowmodel_struc(1)%usum_glb
     v_sum = snowmodel_struc(1)%vsum_glb
     windspd_flag = snowmodel_struc(1)%windspdflg_glb
#endif

! NOTE: winddir_flag never gets used in any subsequent 
!      SnowModel submodels, like snowtran3d
     u_sum = u_sum / real(nx*ny)
     v_sum = v_sum / real(nx*ny)

! Some compilers do not allow both u and v to be 0.0 in
!   the atan2 computation.
     if (abs(u_sum).lt.1e-10) then
       u_sum = 1e-10
     endif
     winddir_flag = rad2deg * atan2(u_sum,v_sum)
     if (winddir_flag.ge.180.0) then
       winddir_flag = winddir_flag - 180.0
     else
       winddir_flag = winddir_flag + 180.0
     endif

 end subroutine micromet_wind


!BOP
!
! !ROUTINE: micromet_solar
! \label{micromet_solar}
!
! !DESCRIPTION:  
!  Calculate solar radiation and estiamte a cloud cover
!  for SWdown inputs to SnowModel.
!
!  Options selected in snowmodel.par
!
! !INTERFACE:
! subroutine micromet_solar( nx,ny,xhour,J_day,topo,rh_grid,Tair_grid,&
!       xlat_grid,Qsi_grid,slope_az,terrain_slope,dt,vegtype,&
!       forest_LAI,T_lapse_rate,Td_lapse_rate,&
!       calc_subcanopy_met,gap_frac,cloud_frac_factor,UTC_flag,&
!       xlon_grid,cloud_frac_grid)

!    use LIS_mpiMod
!    use snowmodel_lsmMod, only : snowmodel_struc

!    implicit none

!    integer, intent(in) :: nx
!    integer, intent(in) :: ny
!    real,    intent(in) :: xhour                 ! model decimal hour
!    integer, intent(in) :: J_day                 ! model day of year 

    ! Input parameters and fields
!    real,    intent(in) :: rh_grid(nx,ny)    
!    real,    intent(in) :: Tair_grid(nx,ny)    
!    real,    intent(in) :: xlat_grid(nx,ny)    
!    real,    intent(in) :: xlon_grid(nx,ny)
!    real,    intent(in) :: topo(nx,ny)    
!    real,    intent(in) :: slope_az(nx,ny)
!    real,    intent(in) :: terrain_slope(nx,ny)
!    real,    intent(in) :: vegtype(nx,ny)    
!    real,    intent(in) :: cloud_frac_grid(nx,ny)    
!    real,    intent(in) :: dt
!    real,    intent(in) :: forest_LAI
!    real,    intent(in) :: T_lapse_rate
!    real,    intent(in) :: Td_lapse_rate
!    real,    intent(in) :: calc_subcanopy_met
!    real,    intent(in) :: gap_frac
!    real,    intent(in) :: cloud_frac_factor
!    real,    intent(in) :: UTC_flag

    ! Output SWdown fields
!    real,    intent(inout) :: Qsi_grid(nx,ny)    ! output SWdown values

! end subroutine micromet_solar

end module micromet_module

