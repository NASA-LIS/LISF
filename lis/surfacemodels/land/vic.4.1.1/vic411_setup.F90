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
!BOP
!
! !ROUTINE: vic411_setup
! \label{vic411_setup}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 01 Jun 2012; Shugong Wang, Retrieve vegetation classes of grid and
!                            send to C routine to support COMPUTE_TREELINE
! 14 Aug 2013; Shugong Wang, modifed for LIS-7
! !INTERFACE:
subroutine vic411_setup()
! !USES:
   use LIS_coreMod,   only : LIS_rc, LIS_domain
   use vic411_lsmMod, only : vic411_struc
   use LIS_LMLCMod,   only : LIS_LMLC         ! added by Shugong Wang
! !ARGUMENTS: 
!
! !DESCRIPTION:
!  This routine initialized VIC's internal data structures.
! 
!EOP

  implicit none

  integer :: n, ngrids, npatch, i, index
  real, allocatable, dimension(:) :: lat_array, lon_array
  integer, allocatable, dimension(:) :: t2gindex, vegclasses
  integer :: vtscheme, VNT

  real*8 :: vic411_find_min_tfactor
  real*8 :: vic411_get_max_snow_temp
  ! added by Shugong Wang 
  integer :: nvegs, row, col , k ,j                                    ! number of vegetation classes
  real(kind=4) , allocatable, dimension(:,:) :: veg_fracs   ! array of fraction for each vegetatoin class
  real(kind=4), allocatable, dimension(:) :: veg_fracs_all  
  do n=1,LIS_rc%nnest
     ngrids = LIS_rc%ngrid(n)
     npatch = LIS_rc%npatch(n,LIS_rc%lsm_index)
     nvegs  = LIS_rc%nvegtypes            ! added by Shugong

     allocate(lat_array(npatch))
     allocate(lon_array(npatch))
     allocate(t2gindex(npatch))
     allocate(vegclasses(npatch))
     allocate(veg_fracs(npatch, nvegs)) ! added by Shugong
     allocate(veg_fracs_all(npatch * nvegs))
     do i = 1, npatch
        index = LIS_domain(n)%tile(i)%index
        lat_array(i) = LIS_domain(n)%grid(index)%lat
        lon_array(i) = LIS_domain(n)%grid(index)%lon
        t2gindex(i) = LIS_domain(n)%tile(i)%index
        vegclasses(i) = LIS_domain(n)%tile(i)%vegt
        row = LIS_domain(n)%tile(i)%row                      ! added by Shugong 
        col = lIS_domain(n)%tile(i)%col                      ! added by Shugong
        j = 1
        do k = 1, LIS_rc%nsurfacetypes
           if ( k == LIS_rc%waterclass   .or. &
                k == LIS_rc%wetlandclass .or. &
                k == LIS_rc%snowclass    .or. &
                k == LIS_rc%glacierclass ) then
                ! skip water points
                cycle
            else
               veg_fracs(i, j) = LIS_LMLC(n)%landcover(col, row, k)
               j = j + 1
           endif
        enddo
     enddo

     j = 1
     do i=1, npatch
        do k=1, nvegs
          veg_fracs_all(j) = veg_fracs(i, k)
          j = j + 1  
        enddo
     enddo  

     vtscheme = vic411_struc(n)%veg_tiling_scheme
     VNT = vic411_struc(n)%NT
     call setup_vic411(ngrids, npatch, lat_array, lon_array,       &
                       t2gindex, vegclasses, vtscheme, VNT,        &
                       veg_fracs_all,                              &
                       trim(vic411_struc(n)%global_param)//char(0))

     do i = 1, npatch
        vic411_struc(n)%vic(i)%min_Tfactor = vic411_find_min_tfactor(n,i)
     enddo
     vic411_struc(n)%MAX_SNOW_TEMP = vic411_get_max_snow_temp(n)

     deallocate(lat_array)
     deallocate(lon_array)
     deallocate(t2gindex)
     deallocate(vegclasses)
  enddo

end subroutine vic411_setup
