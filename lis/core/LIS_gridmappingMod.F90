!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_gridmappingMod
!BOP
!
! !MODULE: LIS_gridmappingMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to map parameter fields to 
!   the target LIS output grid.
!
! !REVISION HISTORY: 
!  10 Jul 2012;  Kristi Arsenault;  Initial Specification
!  10 Feb 2014;  Kristi Arsenault;  Updated routines to include add. options
!  10 Jan 2015;  Kristi Arsenault;  Added to LIS-7
! 
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_RunDomainPts        ! Input Parameter Grid Array and points
!EOP

contains

!BOP
! 
! !ROUTINE: LIS_RunDomainPts
! \label{LIS_RunDomainPts}
! 
! !INTERFACE: 
   subroutine LIS_RunDomainPts( n, param_proj, param_grid, &
                  glpnc, glpnr, subparam_nc, subparam_nr,  &
                  subparam_gridDesc, lat_line, lon_line )

! !USES: 
   use LIS_coreMod,  only : LIS_rc, LIS_domain
   use LIS_logMod,   only : LIS_logunit, LIS_endrun
   use map_utils


! !ARGUMENTS: 
   integer, intent(in)      :: n
   character(50),intent(in) :: param_proj
   real,    intent(in)      :: param_grid(20)
   integer, intent(out)     :: glpnc, glpnr
   real,    intent(out)     :: subparam_gridDesc(20)
   integer, intent(out)     :: subparam_nc, subparam_nr
   integer, allocatable, intent(out) :: lat_line(:,:)
   integer, allocatable, intent(out) :: lon_line(:,:)

! !DESCRIPTION: 
!  This subroutine creates the input parameter grid information
!  required for reading in just a subsetted or the full global
!  parameter domain.
!
! REVISION HISTORY:
!  13FEB2014 -- K.R. Arsenault: Initial Specification
! 
!EOP      
   type(proj_info)  :: subset_paramproj
   integer          :: c, r, i
   real             :: lisdom_max_lat, lisdom_max_lon
   real             :: lisdom_min_lat, lisdom_min_lon
   real*8           :: lisdom_xres_ll, lisdom_yres_ll
   real*8           :: lisdom_xres_ur, lisdom_yres_ur
   real*8           :: plon_search, plat_search 
   real             :: subparm_lllat_ext, subparm_lllon_ext
   real             :: subparm_urlat_ext, subparm_urlon_ext
   real,allocatable :: rlat(:,:)
   real,allocatable :: rlon(:,:)
   real             :: rmin
   real,allocatable :: lats(:)
   real             :: diff_lon
   integer          :: nlats
! _____________________________________________________

   glpnr = 0
   glpnc = 0
   subparam_nc = 0 
   subparam_nr = 0
   subparam_gridDesc = 0.

   subparm_lllat_ext = -9999.
   subparm_urlat_ext = -9999.
   subparm_lllon_ext = -9999.
   subparm_urlon_ext = -9999.

   if(allocated(lat_line)) then
      deallocate(lat_line)
   end if
   if(allocated(lon_line)) then
      deallocate(lon_line)
   end if

! -------------------------------------------------------------
! Calculate extents for subsetted parameter file areas 
!  (for reading in files)
! -------------------------------------------------------------

!- Calculate extents for run domain:
   allocate(rlat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
   allocate(rlon(LIS_rc%lnc(n),LIS_rc%lnr(n)))
   rlat = 0.; rlon = 0.
   lisdom_max_lat = -10000
   lisdom_max_lon = -10000
   lisdom_min_lat = 10000
   lisdom_min_lon = 10000
   do r = 1, LIS_rc%lnr(n)
      do c = 1, LIS_rc%lnc(n)
         call ij_to_latlon( LIS_domain(n)%lisproj,float(c),float(r),&
                            rlat(c,r),rlon(c,r) )
         if(rlat(c,r).gt.lisdom_max_lat) lisdom_max_lat = rlat(c,r)
         if(rlon(c,r).gt.lisdom_max_lon) lisdom_max_lon = rlon(c,r)
         if(rlat(c,r).lt.lisdom_min_lat) lisdom_min_lat = rlat(c,r)
         if(rlon(c,r).lt.lisdom_min_lon) lisdom_min_lon = rlon(c,r)
      enddo
   enddo

!- Calculate LIS run domain corner point resolutions: 
   lisdom_xres_ll = abs(rlon(2,1)-rlon(1,1))
   lisdom_yres_ll = abs(rlat(1,2)-rlat(1,1))
   lisdom_xres_ur = abs(rlon(LIS_rc%lnc(n),  LIS_rc%lnr(n))  &
                  - rlon(LIS_rc%lnc(n)-1,LIS_rc%lnr(n)) )
   lisdom_yres_ur = abs(rlat(LIS_rc%lnc(n),  LIS_rc%lnr(n))  &
                  - rlat(LIS_rc%lnc(n),  LIS_rc%lnr(n)-1) )

! -------------------------------------------------------------
!  Select Parameter File Projection Type
! -------------------------------------------------------------

  select case ( param_proj )

! ------------------------------------------------

!= LAT-LON Geographic Coordinate System (GCS):
   case ( "latlon" )

  !- Calculate total points for global (or "full") parameter file domains:
     glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
     glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9) ) + 1
     !glpnc = nint(diff_lon(param_grid(8),param_grid(5))/param_grid(9) ) + 1

  !- LIS RUN DOMAIN GRID INFORMATION: 
  !   Used to determine parameter domain to be read in
     select case ( LIS_rc%lis_map_proj )

       case( "latlon" )

      ! Parameter grid resolution SAME as LIS run grid resolution:
        if( param_grid(9) == (LIS_rc%gridDesc(n,9)) ) then
          subparm_lllat_ext = lisdom_min_lat 
          subparm_lllon_ext = lisdom_min_lon 
          subparm_urlat_ext = lisdom_max_lat 
          subparm_urlon_ext = lisdom_max_lon 

      ! Parameter grid resolution < as LIS run grid resolution:
        elseif( param_grid(9) < (LIS_rc%gridDesc(n,9)) ) then
          subparm_lllat_ext = lisdom_min_lat - (lisdom_yres_ll/2.) + &
                             (param_grid(10)/2.0)
          subparm_lllon_ext = lisdom_min_lon - (lisdom_xres_ll/2.) + &
                             (param_grid(9)/2.0)
          subparm_urlat_ext = lisdom_max_lat + (lisdom_yres_ur/2.) - &
                             (param_grid(10)/2.0)
          subparm_urlon_ext = lisdom_max_lon + (lisdom_xres_ur/2.) - &
                             (param_grid(9)/2.0)

      ! Parameter grid resolution > as LIS run grid resolution:
        elseif( param_grid(9) > (LIS_rc%gridDesc(n,9)) ) then
          subparm_lllat_ext = param_grid(4)
          subparm_urlat_ext = param_grid(7)
          subparm_lllon_ext = param_grid(5)
          subparm_urlon_ext = param_grid(8)
        endif
 
   !- Lambert conformal LIS run domain:
      case( "lambert" )

      ! Parameter grid resolution <= as LIS run grid resolution:
        if( param_grid(9) <= (LIS_rc%gridDesc(n,9)/100.)) then

       !- Locate parameter longitude extents:
          do i = 1, nint(param_grid(2))
             plon_search = param_grid(5) + (i-1)*param_grid(9)
             if( plon_search > (lisdom_min_lon-lisdom_xres_ll) .and. &
                 plon_search <= lisdom_min_lon ) then
               subparm_lllon_ext = plon_search
             endif
             if( plon_search < (lisdom_max_lon+lisdom_xres_ur) .and. &
                 plon_search >=(lisdom_max_lon+lisdom_xres_ur-param_grid(9)*2.0) ) then
               subparm_urlon_ext = plon_search
             endif
          end do
       !- Locate parameter latitude extents:
          do i = 1, nint(param_grid(3))
             plat_search = param_grid(4) + (i-1)*param_grid(10)
             if( plat_search > (lisdom_min_lat-lisdom_yres_ll) .and. &
                 plat_search <= lisdom_min_lat ) then
               subparm_lllat_ext = plat_search
             endif
             if( plat_search <= (lisdom_max_lat+lisdom_yres_ur)  .and. &
                 plat_search >  lisdom_max_lat ) then
               subparm_urlat_ext = plat_search
             endif
          end do

      ! Parameter grid resolution > as LIS run grid resolution:
        elseif( param_grid(9) > (LIS_rc%gridDesc(n,9)/100.)) then

       !- Locate parameter longitude extents:
          do i = 1, nint(param_grid(2))
             plon_search = param_grid(5) + (i-1)*param_grid(9)
             if( lisdom_min_lon > ( plon_search-(param_grid(9)/2) ) .and. &
                 lisdom_min_lon < ( plon_search+(param_grid(9)/2))) then
               subparm_lllon_ext = plon_search
             endif
             if( lisdom_max_lon > ( plon_search-(param_grid(9)/2) ) .and. &
                 lisdom_max_lon < ( plon_search+(param_grid(9)/2))) then
               subparm_urlon_ext = plon_search
             endif
          end do

       !- Locate parameter latitude extents:
          do i = 1, nint(param_grid(3))
             plat_search = param_grid(4) + (i-1)*param_grid(10)
             if( lisdom_min_lat > ( plat_search-(param_grid(10)/2)) .and. &
                 lisdom_min_lat < ( plat_search+(param_grid(10)/2)) )  then
               subparm_lllat_ext = plat_search
             endif
             if( lisdom_max_lat > ( plat_search-(param_grid(10)/2)) .and. &
                 lisdom_max_lat < ( plat_search+(param_grid(10)/2)) )  then
               subparm_urlat_ext = plat_search
             endif
          end do
        endif  ! End of Lambert run domain

   !- Gaussian LIS run domain:
!      case( "gaussian" )

   !- Mercator LIS run domain:
      case( "mercator" )
! Parameter grid resolution <= as LIS run grid resolution:
        if( param_grid(9) <= (LIS_rc%gridDesc(n,9)/100.)) then

       !- Locate parameter longitude extents:
          do i = 1, nint(param_grid(2))
             plon_search = param_grid(5) + (i-1)*param_grid(9)
             if( plon_search > (lisdom_min_lon-lisdom_xres_ll) .and. &
                 plon_search <= lisdom_min_lon ) then
               subparm_lllon_ext = plon_search
             endif
             if( plon_search < (lisdom_max_lon+lisdom_xres_ur) .and. &
                 plon_search >=(lisdom_max_lon+lisdom_xres_ur-param_grid(9)*2.0) ) then
               subparm_urlon_ext = plon_search
             endif
          end do
       !- Locate parameter latitude extents:
          do i = 1, nint(param_grid(3))
             plat_search = param_grid(4) + (i-1)*param_grid(10)
             if( plat_search > (lisdom_min_lat-lisdom_yres_ll) .and. &
                 plat_search <= lisdom_min_lat ) then
               subparm_lllat_ext = plat_search
             endif
             if( plat_search <= (lisdom_max_lat+lisdom_yres_ur)  .and. &
                 plat_search >  lisdom_max_lat ) then
               subparm_urlat_ext = plat_search
             endif
          end do

      ! Parameter grid resolution > as LIS run grid resolution:
        elseif( param_grid(9) > (LIS_rc%gridDesc(n,9)/100.)) then

       !- Locate parameter longitude extents:
          do i = 1, nint(param_grid(2))
             plon_search = param_grid(5) + (i-1)*param_grid(9)
             if( lisdom_min_lon > ( plon_search-(param_grid(9)/2) ) .and. &
                 lisdom_min_lon < ( plon_search+(param_grid(9)/2))) then
               subparm_lllon_ext = plon_search
             endif
             if( lisdom_max_lon > ( plon_search-(param_grid(9)/2) ) .and. &
                 lisdom_max_lon < ( plon_search+(param_grid(9)/2))) then
               subparm_urlon_ext = plon_search
             endif
          end do

       !- Locate parameter latitude extents:
          do i = 1, nint(param_grid(3))
             plat_search = param_grid(4) + (i-1)*param_grid(10)
             if( lisdom_min_lat > ( plat_search-(param_grid(10)/2)) .and. &
                 lisdom_min_lat < ( plat_search+(param_grid(10)/2)) )  then
               subparm_lllat_ext = plat_search
             endif
             if( lisdom_max_lat > ( plat_search-(param_grid(10)/2)) .and. &
                 lisdom_max_lat < ( plat_search+(param_grid(10)/2)) )  then
               subparm_urlat_ext = plat_search
             endif
          end do
        endif  ! End of mercator run domain

   !- Polar Stereographic LIS run domain:
!      case( "polar" )

     end select

  !- Estimate final number of subsetted parameter points:
     subparam_nc = nint((subparm_urlon_ext - subparm_lllon_ext)/param_grid(9) ) + 1
     !subparam_nc = nint(diff_lon(subparm_urlon_ext,subparm_lllon_ext)/param_grid(9) ) + 1
     subparam_nr = nint((subparm_urlat_ext - subparm_lllat_ext)/param_grid(10)) + 1

  !- Set subsetted parameter grid array inputs:
     subparam_gridDesc(1)  = 0.    ! Latlon
     subparam_gridDesc(2)  = float(subparam_nc)
     subparam_gridDesc(3)  = float(subparam_nr)
     subparam_gridDesc(4)  = subparm_lllat_ext
     subparam_gridDesc(5)  = subparm_lllon_ext
     subparam_gridDesc(6)  = 128.
     subparam_gridDesc(7)  = subparm_urlat_ext
     subparam_gridDesc(8)  = subparm_urlon_ext
     subparam_gridDesc(9)  = param_grid(9)
     subparam_gridDesc(10) = param_grid(10)
     subparam_gridDesc(11) = 64.
     subparam_gridDesc(20) = 64.

  !- Double check subsetted and global parameter number of rows and columns:
     if( subparam_nc > glpnc .or. subparam_nr > glpnr ) then
       write(LIS_logunit,*) "[ERR]: The number of *subsetted* row or column points"
       write(LIS_logunit,*) "          EXCEEDS the total *global* row or column points "
       write(LIS_logunit,*) "          for the input parameter file, respectively."
       write(LIS_logunit,*) " "
       if( subparam_nr > glpnr ) then
          write(LIS_logunit,*) "  - The subsetted, global row points ::",subparam_nr,">",glpnr
       elseif( subparam_nc > glpnc ) then
          write(LIS_logunit,*) "  - The subsetted, global col points ::",subparam_nc,">",glpnc
       endif
       write(LIS_logunit,*) " "
       write(LIS_logunit,*) "  Please check your specified upper-right corner longitude and "
       write(LIS_logunit,*) "   latitude extents, especially if you have entered a gridcell "
       write(LIS_logunit,*) "   resolution that does not divide *nicely* into a 1-deg gridcell."
       write(LIS_logunit,*) "  Program stopping ..."
       call LIS_endrun
     endif

  !- Set up map_set parameter array ...
     call map_set( PROJ_LATLON, subparam_gridDesc(4), subparam_gridDesc(5), &
                   0.0, subparam_gridDesc(9), subparam_gridDesc(10), 0.0,   &
                   subparam_nc, subparam_nr, subset_paramproj )

     allocate( lat_line(subparam_nc,subparam_nr) )
     allocate( lon_line(subparam_nc,subparam_nr) )
     lat_line = 0.; lon_line = 0. 
     deallocate(rlat,rlon)
     allocate( rlat(subparam_nc,subparam_nr) )
     allocate( rlon(subparam_nc,subparam_nr) )
     rlat = 0.; rlon = 0.

  !- Extract LIS-target domain lat and long 2-d arrays:
     do r = 1, subparam_nr
        do c = 1, subparam_nc
           call ij_to_latlon( subset_paramproj, float(c), float(r),&
                              rlat(c,r), rlon(c,r) )
           lat_line(c,r) = nint((rlat(c,r)-param_grid(4))/param_grid(10))+1
           lon_line(c,r) = nint((rlon(c,r)-param_grid(5))/param_grid(9))+1
           !lon_line(c,r) = nint(diff_lon(rlon(c,r),param_grid(5))/param_grid(9))+1
        enddo
     enddo
     deallocate(rlat,rlon)

! ------------------------------------------------
!  Future Parameter File Projection Types:
! ------------------------------------------------
!= Lambert conformal projection:
!   case ( "lambert" )

! ------------------------------------------------
!= Mercator projection:
!   case ( "mercator" )

! ------------------------------------------------
!= Guassian grid system:

   case ( "gaussian" )

      !glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1
      glpnc = nint(diff_lon(param_grid(8),param_grid(5))/param_grid(9)) + 1
      glpnr = param_grid(3)

    ! Assume that parameter grid resolution SAME as LIS run grid resolution
      subparm_lllat_ext = lisdom_min_lat
      subparm_lllon_ext = lisdom_min_lon
      subparm_urlat_ext = lisdom_max_lat
      subparm_urlon_ext = lisdom_max_lon

  !- Estimate final number of subsetted parameter points:
      subparam_nc = LIS_rc%lnc(n)
      subparam_nr = LIS_rc%lnr(n) 

  !- Set subsetted parameter grid array inputs:
     subparam_gridDesc(1)  = 4  
     subparam_gridDesc(2)  = float(subparam_nc)
     subparam_gridDesc(3)  = float(subparam_nr)
     subparam_gridDesc(4)  = subparm_lllat_ext
     subparam_gridDesc(5)  = subparm_lllon_ext
     subparam_gridDesc(6)  = 128.
     subparam_gridDesc(7)  = subparm_urlat_ext
     subparam_gridDesc(8)  = subparm_urlon_ext
     subparam_gridDesc(9)  = param_grid(9)
     subparam_gridDesc(10) = param_grid(10)
     subparam_gridDesc(11) = 64.
     subparam_gridDesc(20) = 64.
  
  !- Set up map_set parameter array ...
!    subset_paramproj = LIS_domain(n)%lisproj
     call map_set( PROJ_GAUSS, subparam_gridDesc(4), subparam_gridDesc(5), &
                   subparam_gridDesc(9), subparam_gridDesc(3), subparam_gridDesc(4), &
                   subparam_gridDesc(7), subparam_nc, subparam_nr, subset_paramproj )

     allocate( lat_line(subparam_nc,subparam_nr) )
     allocate( lon_line(subparam_nc,subparam_nr) )
     lat_line = 0.; lon_line = 0.

  !- Extract LIS-target domain lat and long 2-d arrays:
     rmin = 1
     if (abs(rlat(1,1)-param_grid(4)).gt.0.001) then
       nlats = int(param_grid(3))
       allocate(lats(nlats))
       call gaussian_comp_lats(nlats, lats)
       do r = 1, nlats
         if (abs(lats(r)-subset_paramproj%truelat1).le.0.001) then
           rmin = r
           exit
         endif
       enddo
       deallocate(lats)
     endif
     do r = 1, subparam_nr
        do c = 1, subparam_nc
           lat_line(c,r) = rmin + r - 1
           !lon_line(c,r) = nint((rlon(c,r)-param_grid(5))/param_grid(9))+1
           lon_line(c,r) = nint( diff_lon(rlon(c,r),param_grid(5)) / &
                                 param_grid(9) ) + 1
        enddo
     enddo
     deallocate(rlat,rlon)

!      glnr = gaussian_find_row(output_gridDesc(4)) -  &
!!KRA           gaussian_find_row(output_gridDesc(44)) + 1
!              gaussian_find_row(output_gridDesc(1)) + 1

!      glnc = gaussian_find_col(output_gridDesc(5)) -  &
!!KRA           gaussian_find_col(output_gridDesc(45)) + 1
!              gaussian_find_col(output_gridDesc(2)) + 1

!      line1 = gaussian_find_row(LIS_rc%gridDesc(n,4))  -   &
!           gaussian_find_row(LIS_rc%gridDesc(n,44)) + 1
!      line2 = gaussian_find_col(LIS_rc%gridDesc(n,5))  -   &
!           gaussian_find_col(LIS_rc%gridDesc(n,45)) + 1
!
!      do r=1,LIS_rc%lnr(n)
!         do c=1,LIS_rc%lnc(n)
!            glnc = line2+c-1
!            glnr = line1+r-1
!            line = (glnr-1)*nint(LIS_rc%gridDesc(n,42))+glnc
!            read(ftn,rec=line,iostat=istat) array(c,r)
!            if( istat .ne. 0 ) then
!               message(1) = 'program:  LIS'
!               message(2) = '  routine:  read2DData'
!               message(3) = '  iostat != 0'
!               call LIS_abort( message )
!               call LIS_endrun
!            endif
!         enddo
!      enddo

! ------------------------------------------------
!= Polar stereographic and HRAP projections:
!   case ( "polar", "hrap" )
!      call map_set(PROJ_PS,LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,5), &
!                   LIS_rc%gridDesc(n,8)*1000.0,                       &
!                   LIS_rc%gridDesc(n,11),LIS_rc%gridDesc(n,10),0.0,   &
!                   LIS_rc%lnc(n),LIS_rc%lnr(n),proj)

!      do r=1,LIS_rc%lnr(n)
!         do c=1,LIS_rc%lnc(n)
!            call ij_to_latlon(proj,float(c),float(r),rlat(c,r),rlon(c,r))
!         enddo
!      enddo

!      call map_set(PROJ_PS,LIS_rc%lc_gridDesc(n,1),LIS_rc%lc_gridDesc(n,2),  &
!                   LIS_rc%lc_gridDesc(n,6)*1000.0,                           &
!                   LIS_rc%lc_gridDesc(n,4),LIS_rc%lc_gridDesc(n,3),0.0,      &
!                   int(LIS_rc%lc_gridDesc(n,7)),int(LIS_rc%lc_gridDesc(n,8)),&
!                   proj)

!      do r=1,LIS_rc%lnr(n)
!         do c=1,LIS_rc%lnc(n)
!            call latlon_to_ij(proj,rlat(c,r),rlon(c,r),ctmp,rtmp)
!            line1 = nint(rtmp)
!            line2 = nint(ctmp)
!         enddo
!      enddo

! ------------------------------------------------
!= Universal Tranverse Mercator coordinate system:
!   case ( "UTM" )

!  rlat/rlon used here to store northing and easting
!      do r=1,LIS_rc%lnr(n)
!         do c=1,LIS_rc%lnc(n)
!            rlat(c,r) = LIS_rc%gridDesc(n,4)+(r-1)*LIS_rc%gridDesc(n,9)
!            rlon(c,r) = LIS_rc%gridDesc(n,5)+(c-1)*LIS_rc%gridDesc(n,9)
!         enddo
!      enddo

!      glpnc = param_grid(4)
!      glpnr = 1
!      do r=1,LIS_rc%lnr(n)
!         do c=1,LIS_rc%lnc(n)
!            line1 = nint((rlat(c,r)-gridDesc(2))/gridDesc(6))+1
!            line2 = nint((rlon(c,r)-gridDesc(3))/gridDesc(6))+1
!         enddo
!      enddo

! ------------------------------------------------

   case default
      write(*,*) "This parameter projection(",trim(param_proj),") is not supported ..." 
      write(*,*) "Program stopping ... "
      call LIS_endrun

  end select


  end subroutine LIS_RunDomainPts

! __________________________________________________________________________


end module LIS_gridmappingMod

