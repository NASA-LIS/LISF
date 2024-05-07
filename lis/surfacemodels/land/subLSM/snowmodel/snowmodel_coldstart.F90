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
! !ROUTINE: SnowModel_coldstart
! \label{SnowModel_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!  05 Aug 2022: Kristi Arsenault; Added SnowModel coldstart routine
!
! !INTERFACE:
subroutine SnowModel_coldstart(mtype)

! USES:
   use LIS_coreMod
   use LIS_logMod
   use LIS_timeMgrMod, only: LIS_date2time

   use SnowModel_lsmMod
   use snowmodel_inc
   use snowmodel_vars
!
! !DESCRIPTION:
!
!  This routine initializes the SnowModel state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
    implicit none
    integer :: mtype
    integer :: t, l, n

    do n=1, LIS_rc%nnest
       if( trim(LIS_rc%startcode) .eq. "coldstart" ) then

          write(LIS_logunit,*) "[INFO] SnowModel_coldstart -- initializing SnowModel states"

          ! Initialize snow fields 
          snowmodel_struc(n)%sm(:)%sden = 0.
          snowmodel_struc(n)%sm(:)%sublim = 0.
          snowmodel_struc(n)%sm(:)%swemelt = 0.
          snowmodel_struc(n)%sm(:)%runoff = 0.

          ! Initialize other state fields used 
          snowmodel_struc(n)%sm(:)%swe_depth = 0.
          snowmodel_struc(n)%sm(:)%snow_depth = 0.
          snowmodel_struc(n)%sm(:)%snow_d = 0.
          snowmodel_struc(n)%sm(:)%canopy_int = 0.
          snowmodel_struc(n)%sm(:)%soft_snow_d = 0.
          snowmodel_struc(n)%sm(:)%ro_snow_grid = 0.
          snowmodel_struc(n)%sm(:)%ro_soft_snow_old = 0.
          snowmodel_struc(n)%sm(:)%snow_d_init = 0.
          snowmodel_struc(n)%sm(:)%swe_depth_old = 0.
          snowmodel_struc(n)%sm(:)%canopy_int_old = 0.
          snowmodel_struc(n)%sm(:)%topo = 0.
          snowmodel_struc(n)%sm(:)%sum_sprec = 0.

          ! SnowModel "initialize" routine code (from preprocess_code.f90)
          do i=1,nx
            do j=1,ny
              ! Fill the summing arrays.
              sum_runoff(i,j) = 0.0
              sum_prec(i,j) = 0.0
              sum_sprec(i,j) = 0.0
              sum_qsubl(i,j) = 0.0
              sum_trans(i,j) = 0.0
              sum_unload(i,j) = 0.0
              sum_Qcs(i,j) = 0.0
              sum_glacmelt(i,j) = 0.0
              sum_swemelt(i,j) = 0.0
              sum_d_canopy_int(i,j) = 0.0
              sum_sfcsublim(i,j) = 0.0

              ! Define the initial snow-depth distributions.
              snow_d_init(i,j) = snow_d_init_const
              snow_d(i,j) = snow_d_init(i,j)
              snow_depth(i,j) = snow_d_init(i,j)
              canopy_int(i,j) = 0.0
              soft_snow_d(i,j) = snow_d(i,j)
              ro_snow_grid(i,j) = ro_snow
              swe_depth(i,j) = snow_d(i,j) * ro_snow_grid(i,j) / ro_water
              ro_soft_snow_old(i,j) = 50.0
              swe_depth_old(i,j) = swe_depth(i,j)
              canopy_int_old(i,j) = canopy_int(i,j)

              ! Initialize the multi-layer snowpack arrays.
              KK(i,j) = 0
              tslsnowfall(i,j) = tsls_threshold

            enddo
          enddo

          do i=1,nx
            do j=1,ny
              do k=1,nz_max
                snod_layer(i,j,k) = 0.0
                swed_layer(i,j,k) = 0.0
                ro_layer(i,j,k) = ro_snow
                T_old(i,j,k) = 273.15
                gamma(i,j,k) = 0.138 - 1.01 * (ro_layer(i,j,k)/1000.0) + &
                               3.233 * (ro_layer(i,j,k)/1000.0)**2
                diam_layer(i,j,k) = 0.5 / 1000.0
              enddo
            enddo
          enddo

          if (topoflag.eq.1.0) then
            do i=1,nx
              do j=1,ny
                topo(i,j) = topo_land(i,j) + snow_d(i,j)
              enddo
            enddo
          elseif (topoflag.eq.0.0) then
            do i=1,nx
              do j=1,ny
                topo(i,j) = topo_land(i,j)
              enddo
            enddo
          endif

       endif
    enddo

end subroutine SnowModel_coldstart
