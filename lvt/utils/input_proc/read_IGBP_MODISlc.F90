!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "../../make/LVT_misc.h"
!BOP
!
! !ROUTINE: read_IGBP_MODISlc
!  \label{read_IGBP_MODISlc}
!
! !REVISION HISTORY:
! 03 Sep 2004: Sujay Kumar; Initial Specification
! 11 Mar 2011: David Mocko, Added modified IGBP MODIS landcover
!
! !INTERFACE:
subroutine read_IGBP_MODISlc(n, num_types, fgrd, maskarray, sfctype)

! !USES:
  use preprocMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_types
  real, intent(inout) :: fgrd(LVT_rc%lnc(n),LVT_rc%lnr(n),LVT_rc%nt)
  real, intent(inout) :: maskarray(LVT_rc%lnc(n),LVT_rc%lnr(n))
  real, intent(inout) :: sfctype(LVT_rc%lnc(n),LVT_rc%lnr(n),LVT_rc%nt)
!
! !DESCRIPTION:
!  This subroutine reads the modified IGBP MODIS landcover data
!  and returns the distribution of vegetation in each grid cell,
!  in a lat/lon projection.  The data has 20 vegetation types.
!  The dataset was modified by NCEP from the original Boston
!  University 1km global vegetation class map based on MODIS
!  data by Mark Friedl and his colleagues at Boston University.
!  It is based on the IGBP classification but excludes the
!  permanent wetland (class 11) and cropland/natural vegetation
!  mosaic (class 14).  Now 3 new classes of tundra have been added
!  by the Land Team/EMC/NCEP led by Ken Mitchell.  Lake is also added.
!  More detail can be found here:
!  ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/ldas/noahlsm/README
!
!  (A) New IGBP\_MODIS\_BU+tundra Landcover Class Legend:
!   1 Evergreen Needleleaf Forest
!   2 Evergreen Broadleaf Forest
!   3 Deciduous Needleleaf Forest
!   4 Deciduous Broadleaf Forest
!   5 Mixed Forests
!   6 Closed Shrublands
!   7 Open Shrublands
!   8 Woody Savannas
!   9 Savannas
!  10 Grasslands
!  11 Permanent Wetland
!  12 Croplands
!  13 Urban and Built-Up
!  14 Cropland/Natural Vegetation Mosaic
!  15 Snow and Ice
!  16 Barren or Sparsely Vegetated
!  17 Ocean
!  18 Wooded Tundra
!  19 Mixed Tundra
!  20 Bare Ground Tundra
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[fgrd]
!    fraction of grid covered by each vegetation type
!   \end{description}
!EOP      

  real, allocatable :: veg(:,:,:)
  real, allocatable :: tsum(:,:)

  real :: isum
  integer :: ftn
  integer :: c,r
  integer :: t, ierr, ios1
  logical :: file_exists
! ____________________________________________________

  fgrd(:,:,:) = 0.0

  LVT_rc%bareclass = 16
  LVT_rc%urbanclass = 13
  LVT_rc%snowclass = 15
  LVT_rc%waterclass = 17
  LVT_rc%wetlandclass = 11
  LVT_rc%glacierclass = 15

  ftn = 100
  allocate(veg(LVT_rc%lnc(n),LVT_rc%lnr(n),LVT_rc%nt), stat=ierr)
  veg = 0.0
  call LVT_verify(ierr,'Error allocating veg')

  inquire(file=trim(LVT_rc%vfile(n)),exist=file_exists) 
  if(.not.file_exists) then 
     write(*,*) 'landcover map: ',trim(LVT_rc%vfile(n)), ' does not exist'
     write(*,*) 'program stopping ...'
     stop
  endif
  write(*,*) 'Reading landcover data ',trim(LVT_rc%vfile(n))
  open(ftn,file=trim(LVT_rc%vfile(n)),status='old',form='unformatted',&
       access ='direct',recl=4,iostat=ios1)

  if (LVT_rc%vfile_form(n).eq.0) then ! lc data is not tiled
     allocate(tsum(LVT_rc%lnc(n),LVT_rc%lnr(n)), stat=ierr)
     tsum = 0.0
     call LVT_verify(ierr,'Error allocating tsum')

     call LVT_readData(n,ftn,LVT_rc%lc_proj,LVT_rc%lc_gridDesc(n,:),tsum)

     do r=1,LVT_rc%lnr(n)
        do c=1,LVT_rc%lnc(n) 
#if ( defined INC_WATER_PTS )
!kludged to include water points
! The right way to do this would be read an appropriate mask and veg 
! file with water points
           if ( tsum(c,r) .le. 0 ) then 
              tsum(c,r) = LVT_rc%waterclass
           endif
!kludge
           if (nint(tsum(c,r)) .ne. LVT_rc%udef ) then 
              veg(c,r,NINT(tsum(c,r))) = 1.0
           endif
#else
           if (nint(tsum(c,r)) .ne. LVT_rc%waterclass ) then 
              veg(c,r,NINT(tsum(c,r))) = 1.0
           endif
#endif
        enddo
     enddo
     deallocate(tsum)

  else !data is tiled
     call LVT_readData(n, ftn, LVT_rc%lc_proj, LVT_rc%lc_gridDesc(n,:), veg)

  endif

  do r=1,LVT_rc%lnr(n)
     do c=1,LVT_rc%lnc(n) 
        isum=0.0
        do t=1,LVT_rc%nt 
#if ( defined INC_WATER_PTS )
           isum=isum+veg(c,r,t)
#else
           if(t.ne.LVT_rc%waterclass) then 
              isum=isum+veg(c,r,t)  !recompute ISUM without water points
           endif
#endif
        enddo
        do t=1,LVT_rc%nt 
           fgrd(c,r,t)=0.0
#if ( defined INC_WATER_PTS )
!kluge
           if(isum.gt.0) then 
              fgrd(c,r,t)=veg(c,r,t)/isum
           else
              fgrd(c,r,17) = 1.0
           endif
#else
           if(isum.gt.0) fgrd(c,r,t)=veg(c,r,t)/isum
#endif
        enddo
#if ( ! defined INC_WATER_PTS )
        if(LVT_rc%waterclass.gt.0) &
             fgrd(c,r,LVT_rc%waterclass) = 0.0
#endif
     end do
  enddo

  deallocate(veg)

  close(ftn)
  
end subroutine read_IGBP_MODISlc
