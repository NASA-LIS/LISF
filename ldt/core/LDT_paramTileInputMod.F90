!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_paramTileInputMod
!BOP
!
! !MODULE: LDT_paramTileInputMod
! 
! !DESCRIPTION: 
!   The routines in this module support parameter class and/or bin
!   fraction estimations for a gridcell.
!
!
! !REVISION HISTORY: 
!  30 Jul 2012   Kristi Arsenault  Initial Specification
!  23 Mar 2017   Shugong Wang  Add param_index_fgrdcalc_pft for JULES 
! 
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: param_index_fgrdcalc  ! Index field grid fraction estimates for tiling
  public :: param_1dbin_areacalc  ! 1 Continuous field's grid fraction estimates for tiling
  public :: param_2dbin_areacalc  ! 2+ Continuous fields' grid fraction estimates for tiling
  public :: param_index_fgrdcalc_pft  !JULES pft grid fraction estimates for tiling

!EOP

contains

!BOP
! 
! !ROUTINE: param_index_fgrdcalc
! \label{param_index_fgrdcalc}
! 
! !INTERFACE: 
  subroutine param_index_fgrdcalc( n, param_proj, spatial_transform, &
                              water_class, num_bins, varcnt, fgrd )

! !USES: 
   use LDT_coreMod, only : LDT_rc, LDT_domain
   use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun

! !ARGUMENTS: 
   integer,      intent(in) :: n
   character(50),intent(in) :: param_proj
   character(50),intent(in) :: spatial_transform
   integer,      intent(in) :: water_class
   integer,      intent(in) :: num_bins
   real,         intent(in) :: varcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
   real,         intent(out):: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
! 
! !DESCRIPTION: 
!  This subroutine calculates the grid fraction (fgrd) for a 3-D
!   tiled parameter field.
!
! REVISION HISTORY:
!  28JUL2012 -- K.R. Arsenault: Initial Specification
! 
!EOP      
   integer :: c, r, t
   real    :: isum
! __________________________________________________________________

!   select case ( param_proj )
!    case ( "latlon" ) 

   fgrd = 0.0

!- Estimate fraction of grid (fgrid) represented by class/type::
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         isum = 0.0

         if( LDT_rc%inc_water_pts ) then
      !- COUNT total number of pixels for each type:  INCLUDING WATER PIXELS
            isum = sum( varcnt(c,r,1:num_bins), &
                   mask=varcnt(c,r,1:num_bins).ne.LDT_rc%udef )
         else
      !- COUNT total number of pixels for each type:  EXCLUDING WATER PIXELS
            do t=1,num_bins
               if( t /= water_class ) then
                 if( varcnt(c,r,t) /= LDT_rc%udef ) then
                     isum = isum + varcnt(c,r,t)
                 end if
               end if
            end do
         endif

      !- Calculate grid fraction for each type:  fgrd(c,r,t)
         do t = 1, num_bins

          ! Variable count ...
            if( isum > 0. ) then

            ! Include water points
              if( LDT_rc%inc_water_pts ) then
                 if(varcnt(c,r,t) /= LDT_rc%udef) fgrd(c,r,t) = varcnt(c,r,t)/isum

            ! Exclude water points
              else
                if( t == water_class ) then
                   fgrd(c,r,t) = 0.
                else
                   if(varcnt(c,r,t)+1 > LDT_rc%udef) then
                      fgrd(c,r,t) = varcnt(c,r,t)/isum
                   endif
                endif
              endif

          ! No variable count at all ... 
            else
               fgrd(c,r,t) = 0.0
            end if

         !- If no subgrid info:
            if( isum <= 0 .and. ( spatial_transform == "none" .or. &
                spatial_transform == "mode") ) then
               fgrd(c,r,water_class) = 1.0
            endif
         enddo

      end do ! end col loop
   end do    ! end row loop

!   case default 
!      print*, 'This parameter projection is not supported...'
!      print*, 'Program stopping ....'
!      stop

!  end select

  end subroutine param_index_fgrdcalc

!BOP
! 
! !ROUTINE: param_index_fgrdcalc_pft
! \label{param_index_fgrdcalc_pft}
! 
! !INTERFACE: 
  subroutine param_index_fgrdcalc_pft( n, param_proj, spatial_transform, &
                                       water_class, num_bins, varcnt, fpft)

! !USES: 
   use LDT_coreMod, only : LDT_rc, LDT_domain
   use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun

! !ARGUMENTS: 
   integer,      intent(in) :: n
   character(50),intent(in) :: param_proj
   character(50),intent(in) :: spatial_transform
   integer,      intent(in) :: water_class
   integer,      intent(in) :: num_bins
   real,         intent(in) :: varcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
!   real,         intent(out):: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
   real,         intent(out):: fpft(LDT_rc%lnc(n),LDT_rc%lnr(n),9)
! 
! !DESCRIPTION: 
!  This subroutine calculates the grid fraction (fgrd) for a 3-D
!   tiled parameter field.
!
! REVISION HISTORY:
!  28JUL2012 -- K.R. Arsenault: Initial Specification
!  23Mar2017 -- Shugng Wang: Revision for JULES PFT
!   
!EOP      
   integer :: c, r, t, igbp_type, k
   real    :: isum, isum_pft
   real    :: tmp_cpft(9), igbp_cnt
! __________________________________________________________________
   real :: table(9,20)=reshape((/        &
     ! Broadleaf Needleleaf  C3 Grass  C4 Grass  Shrub Urban Water Bare soil  Ice    IGBP class
     0.0     , 70.0     ,  20.0   ,  0.0    ,  0.0 , 0.0   , 0.0   , 10.0    ,  0.0   ,        & ! Evergreen needleleaf
     85.0    , 0.0      ,  0.0    ,  10.0   ,  0.0 , 0.0   , 0.0   , 5.0     ,  0.0   ,        & ! Evergreen broadleaf
     0.0     , 65.0     ,  25.0   ,  0.0    ,  0.0 , 0.0   , 0.0   , 10.0    ,  0.0   ,        & ! Deciduous needleleaf
     60.0    , 0.0      ,  5.0    ,  10.0   ,  5.0 , 0.0   , 0.0   , 20.0    ,  0.0   ,        & ! Deciduous broadleaf
     35.0    , 35.0     ,  20.0   ,  0.0    ,  0.0 , 0.0   , 0.0   , 10.0    ,  0.0   ,        & ! Mixed forest
     0.0     , 0.0      ,  25.0   ,  0.0    ,  60.0, 0.0   , 0.0   , 15.0    ,  0.0   ,        & ! Close shrub
     0.0     , 0.0      ,  5.0    ,  10.0   ,  35.0, 0.0   , 0.0   , 50.0    ,  0.0   ,        & ! Open shrub
     50.0    , 0.0      ,  15.0   ,  0.0    ,  25.0, 0.0   , 0.0   , 10.0    ,  0.0   ,        & ! Woody savanna
     20.0    , 0.0      ,  0.0    ,  75.0   ,  0.0 , 0.0   , 0.0   , 5.0     ,  0.0   ,        & ! Savanna
     0.0     , 0.0      ,  66.0   ,  15.7   ,  4.9 , 0.0   , 0.0   , 13.5    ,  0.0   ,        & ! Grassland
     0.0     , 0.0      ,  80.0   ,  0.0    ,  0.0 , 0.0   , 20.0  , 0.0     ,  0.0   ,        & ! Permanent wetland
     0.0     , 0.0      ,  75.0   ,  5.0    ,  0.0 , 0.0   , 0.0   , 20.0    ,  0.0   ,        & ! Cropland
     0.0     , 0.0      ,  0.0    ,  0.0    ,  0.0 , 100.0 , 0.0   , 0.0     ,  0.0   ,        & ! Urban
     5.0     , 5.0      ,  55.0   ,  15.0   ,  10.0, 0.0   , 0.0   , 10.0    ,  0.0   ,        & ! Cropland/natural mosaic
     0.0     , 0.0      ,  0.0    ,  0.0    ,  0.0 , 0.0   , 0.0   , 0.0     ,  100.0 ,        & ! Snow and ice
     0.0     , 0.0      ,  0.0    ,  0.0    ,  0.0 , 0.0   , 0.0   , 100.0   ,  0.0   ,        & ! Barren
     0.0     , 0.0      ,  0.0    ,  0.0    ,  0.0 , 0.0   , 100.0 , 0.0     ,  0.0   ,        & ! Inland water
     0.0     , 0.0      ,  0.0    ,  0.0    ,  0.0 , 0.0   , 0.0   , 100.0   ,  0.0   ,        & ! Wooded Tundra
     0.0     , 0.0      ,  0.0    ,  0.0    ,  0.0 , 0.0   , 0.0   , 100.0   ,  0.0   ,        & ! Mixed Tundra
     0.0     , 0.0      ,  0.0    ,  0.0    ,  0.0 , 0.0   , 0.0   , 100.0   ,  0.0/), (/9, 20/))! Bared ground Tundra


!- Estimate fraction of grid (fgrid) represented by class/type::
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
        isum = 0.0
        isum_pft = 0   
        tmp_cpft(:) = 0
        do t=1, num_bins
          igbp_type = t
          igbp_cnt  =varcnt(c,r,t)
          do k=1,9
            tmp_cpft(k) = tmp_cpft(k) + igbp_cnt*table(k, igbp_type)
          enddo
        enddo
        isum_pft = sum(tmp_cpft)

        do k=1, 9
          if(isum_pft >0) then
            fpft(c,r,k) = tmp_cpft(k)/isum_pft
          else
            fpft(c,r,k) = LDT_rc%udef
          endif
        enddo
        
        ! if landice only covers a portion of the grid box,
        ! add the portion to baresoil because landice is 
        ! exclusive with soil. landice should fully occupy a gridbox.
        if(fpft(c,r,9)>0 .and. fpft(c,r,9)<1.0) then
          if(fpft(c,r,9)<0.5) then !land ice covers less than a half of the grid
            fpft(c,r,8) = fpft(c,r,8) + fpft(c,r,9)
            fpft(c,r,9) = 0.0
          else ! set the grid into land ice
            fpft(c,r,1:8) = 0 
            fpft(c,r,9) = 1.0
          endif
        endif
      end do ! end col loop
   end do    ! end row loop

  end subroutine param_index_fgrdcalc_pft
  

!BOP
!
! !ROUTINE: param_1dbin_areacalc
! \label{param_1dbin_areacalc}
!
! !INTERFACE:
  subroutine param_1dbin_areacalc( n, num_bins, inpts, outpts, &
                         n11, undef, invar1d, fgrd, areaave )
! !USES:
   use LDT_coreMod, only : LDT_rc, LDT_domain
   use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun

! !ARGUMENTS:
   integer,      intent(in) :: n
   integer,      intent(in) :: num_bins
   integer,      intent(in) :: inpts
   integer,      intent(in) :: outpts
   integer,      intent(in) :: n11(inpts)
   real,         intent(in) :: undef
   real,         intent(in) :: invar1d(inpts)
   real,         intent(out):: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
   real,         intent(out):: areaave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
!
! !DESCRIPTION:
!  This subroutine calculates the grid fraction (fgrd), making a 3-D
!   tiled parameter field, along the bin averages for the field,
!
! !REFERENCE:
! - Currently based on the University of Washington VIC Elevation 
!    band calculation approach:
!  http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/PrepElevBand.shtml
!
! REVISION HISTORY:
!  28JUL2012 -- K.R. Arsenault: Initial Specification
!
!EOP
   integer :: c, r, i, j, z
   integer :: bin
   real    :: zmin, zmax, zmean
   real    :: max_spread
   real    :: bin_width
   real    :: bin_llim, bin_ulim
   real    :: pts_perbin(outpts,num_bins)
   real    :: area_binsum(outpts,num_bins)
   integer :: pixels_pergrid(outpts)
   integer :: z_cnt(outpts)
   real, allocatable :: z_ij_gridcell(:,:)

! __________________________________________________________________

   pts_perbin = 0.
   area_binsum = 0.0
   fgrd = 0.0
   areaave = 0.0
   pixels_pergrid = 0

!- Estimate number of pixels per gridcell (coarse domain):
   do i = 1, inpts
      if( n11(i) .ne. 0 ) then
         pixels_pergrid(n11(i)) = pixels_pergrid(n11(i)) + 1
      endif
   enddo

   allocate( z_ij_gridcell( outpts, maxval(pixels_pergrid(:))) )
   z_ij_gridcell = undef
   z_cnt = 0 
   do i = 1, inpts
      if( n11(i) .ne. 0 ) then
        z_cnt(n11(i)) = z_cnt(n11(i)) + 1
        z_ij_gridcell(n11(i),z_cnt(n11(i))) = invar1d(i)
      end if
   end do  ! input grid loop

   do j = 1, outpts  ! output grid loop

   !- Calculate gridcell stats for variable:
   !- Screen when having all undefined values in gridcell:
      if( minval(z_ij_gridcell(j,:),mask=z_ij_gridcell(j,:).ne.undef) > 10e20 .or. &
          maxval(z_ij_gridcell(j,:),mask=z_ij_gridcell(j,:).ne.undef) > 10e20 .or. &
          sum(z_ij_gridcell(j,:),mask=z_ij_gridcell(j,:).ne.undef) > 10e20 ) then
         zmin  = 0.
         zmax  = 0.
         zmean = 0.
      else
         zmin = minval(z_ij_gridcell(j,:), mask=z_ij_gridcell(j,:).ne.undef)
         zmax = maxval(z_ij_gridcell(j,:), mask=z_ij_gridcell(j,:).ne.undef)
         zmean= sum(z_ij_gridcell(j,:), mask=z_ij_gridcell(j,:).ne.undef) / &
                count( mask=z_ij_gridcell(j,:).ne.undef )
      endif
      max_spread = max( (zmax-zmean),(zmean-zmin) )
      bin_width = (2*max_spread) / float(num_bins)

   !- Designate fine-scale values to a bin:
      do bin = 1, num_bins
         bin_llim = (zmean-max_spread) + float(bin-1)*bin_width
         bin_ulim = (zmean-max_spread) + (float(bin)*bin_width)
         do z = 1, pixels_pergrid(j)
            if( z_ij_gridcell(j,z) >= bin_llim  .and. &
                z_ij_gridcell(j,z) <  bin_ulim  .and. & 
                z_ij_gridcell(j,z) .ne. undef ) then

               pts_perbin(j,bin) = pts_perbin(j,bin) + 1.0
               area_binsum(j,bin) = area_binsum(j,bin) + z_ij_gridcell(j,z)

            elseif( bin == num_bins  .and. &
                    z_ij_gridcell(j,z) == bin_ulim .and. & 
                    z_ij_gridcell(j,z) .ne. undef ) then
               pts_perbin(j,bin) = pts_perbin(j,bin) + 1.0
               area_binsum(j,bin) = area_binsum(j,bin) + z_ij_gridcell(j,z)
         
            endif
         end do  
      end do ! end bin loop
   end do    ! end output grid loop
   deallocate( z_ij_gridcell )

!- Convert 2D to 3D grid arrays:
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         i = i + 1
         do bin = 1, num_bins
         !- If no counts occur for bin:
            if( pts_perbin(i,bin) == 0. ) then
               fgrd(c,r,bin)    = 0.
               areaave(c,r,bin) = 0.
         !- If counts occur for bin:
            else
               fgrd(c,r,bin)    = pts_perbin(i,bin)/sum(pts_perbin(i,:))
               areaave(c,r,bin) = area_binsum(i,bin)/pts_perbin(i,bin)
            endif
         end do
      enddo
   enddo

  end subroutine param_1dbin_areacalc

!BOP
!
! !ROUTINE: param_2dbin_areacalc
! \label{param_2dbin_areacalc}
!
! !INTERFACE:
  subroutine param_2dbin_areacalc( n, num_tiles, inpts, outpts, n11, &
                         undef, invar1, invar2, outvar1, outvar2, fgrd )
! !USES:
   use LDT_coreMod, only : LDT_rc, LDT_domain
   use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun
   use LDT_numericalMethodsMod, only : LDT_quicksort_1arr

! !ARGUMENTS:
   integer,      intent(in) :: n
   integer,      intent(in) :: num_tiles
   integer,      intent(in) :: inpts
   integer,      intent(in) :: outpts
   integer,      intent(in) :: n11(inpts)
   real,         intent(in) :: undef
   real,         intent(in) :: invar1(inpts)
   real,         intent(in) :: invar2(inpts)
   real,         intent(out):: outvar1(LDT_rc%lnc(n),LDT_rc%lnr(n),num_tiles)
   real,         intent(out):: outvar2(LDT_rc%lnc(n),LDT_rc%lnr(n),num_tiles)
   real,         intent(out):: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_tiles)
!
! !DESCRIPTION:
!  This subroutine calculates dominant bin values for combined 
!   continuous arrays and estimates the combine grid fraction (fgrd).
!
! !REFERENCE:
! -- Sujay Kumar, James Geiger,and KR Arsenault 
!
! !NOTES:
! -- Current set up for this routine only supports datasets ranging 0-1.
! -- Some bin-calculation features hardwired (for testing) for now, 
!     but will be updated with more options for future releases.
!
! REVISION HISTORY:
!  06AUG2012 -- K.R. Arsenault: Initial Specification
!
!EOP
   integer :: c, r, i, j, z, i2
   integer :: bin1, bin2
   real    :: isum
   real    :: zmin1, zmax1
   real    :: zmin2, zmax2
   real    :: max_spread1, max_spread2
   real    :: bin_width1, bin_width2
   real    :: bin_llim1, bin_ulim1
   real    :: bin_llim2, bin_ulim2
   real    :: bin_mid1, bin_mid2
   integer :: pixels_pergrid(outpts)
   integer :: z_cnt(outpts)
   real,allocatable :: z_ij_gridcell1(:,:)
   real,allocatable :: z_ij_gridcell2(:,:)

   integer, parameter :: num_bins = 10
   real    :: pts_perbin(outpts,num_bins,num_bins)
   real    :: combo_binarray(num_bins*num_bins)

   real    :: fgrd_1d(outpts,num_tiles)
   real    :: var1_1d(outpts,num_tiles)
   real    :: var2_1d(outpts,num_tiles)
! __________________________________________________________________

   pts_perbin = 0.
   combo_binarray = 0.
   fgrd_1d = 0.;  fgrd = 0.
   var1_1d = 0.;  var2_1d = 0.
   outvar1 = 0.;  outvar2 = 0.
   pixels_pergrid = 0

!- Estimate number of pixels per gridcell (coarse domain):
   do i = 1, inpts
      if( n11(i) .ne. 0 ) then
         pixels_pergrid(n11(i)) = pixels_pergrid(n11(i)) + 1
      endif
   enddo

   allocate( z_ij_gridcell1(outpts, maxval(pixels_pergrid(:))) )
   allocate( z_ij_gridcell2(outpts, maxval(pixels_pergrid(:))) )

   z_ij_gridcell1 = undef
   z_ij_gridcell2 = undef
   z_cnt = 0
   do i = 1, inpts
      if( n11(i) .ne. 0 ) then
        z_cnt(n11(i)) = z_cnt(n11(i)) + 1
        z_ij_gridcell1(n11(i),z_cnt(n11(i))) = invar1(i)
        z_ij_gridcell2(n11(i),z_cnt(n11(i))) = invar2(i)
      endif
   end do  ! input grid loop

   do j = 1, outpts  ! output grid loop

   !- Set min, max and bin-width values:
      zmin1 = 0; zmax1 = 1.0
      bin_width1 = (zmax1 - zmin1) / float(num_bins)

   !- First, loop over array of bins for each field:
      do bin1 = 1, num_bins
         bin_llim1 = zmin1 + float(bin1-1)*bin_width1
         bin_ulim1 = zmin1 + (float(bin1)*bin_width1)

         do bin2 = 1, num_bins
            bin_llim2 = zmin1 + float(bin2-1)*bin_width1
            bin_ulim2 = zmin1 + (float(bin2)*bin_width1)

         !- Second, loop over each field's pixels for same gridcell:
            do z = 1, pixels_pergrid(j)
               if( z_ij_gridcell1(j,z) >= bin_llim1 .and. &
                   z_ij_gridcell1(j,z) <  bin_ulim1 .and. &
                   z_ij_gridcell1(j,z) .ne. undef   .and. &
                   z_ij_gridcell2(j,z) >= bin_llim2 .and. &
                   z_ij_gridcell2(j,z) <  bin_ulim2 .and. &
                   z_ij_gridcell2(j,z) .ne. undef ) then
                 pts_perbin(j,bin1,bin2) = pts_perbin(j,bin1,bin2) + 1.0
               endif
            end do ! End pixel loop count loop

         end do  ! End second bin loop
      end do     ! End first bin loop

   !- Assign 2-D bin matrix to 1D array:
      i = 0;  combo_binarray = 0.
      do bin1 = 1, num_bins
         do bin2 = 1, num_bins; i = i + 1
            combo_binarray(i) = pts_perbin(j,bin1,bin2)
         end do  ! End second bin loop
      end do     ! End first bin loop

   !- Sort to locate dominant combinations:
      call LDT_quicksort_1arr( num_bins*num_bins, combo_binarray )

   !- Estimate fractions or areas of grid for combined fields:
      isum = sum(combo_binarray((num_bins*num_bins-num_bins+1):(num_bins*num_bins)))
      do bin1 = 1, num_bins
         bin_mid1 = (zmin1 + float(bin1)*bin_width1)-(bin_width1*0.5)
         do bin2 = 1, num_bins
            bin_mid2 = (zmin1 + float(bin2)*bin_width1)-(bin_width1*0.5)
            i2 = 0
!            do i = (num_bins*num_bins), (num_bins*num_bins-num_bins+1), -1
            do i = (num_bins*num_bins), (num_bins*num_bins-num_tiles+1), -1
               i2 = i2 + 1
            !- Dominant bin match:
               if( combo_binarray(i) > 0. ) then
                  if(pts_perbin(j,bin1,bin2) == combo_binarray(i)) then
                    fgrd_1d(j,i2) = combo_binarray(i) / isum
                    var1_1d(j,i2) = bin_mid1
                    var2_1d(j,i2) = bin_mid2
                    cycle 
                  endif
            !- No counts for bin:
               else
                  fgrd_1d(j,i2) = 0.
               end if
            end do
         end do  ! end bin loops
      end do
   end do   ! end output grid loop
   deallocate( z_ij_gridcell1, z_ij_gridcell2 )

!- Convert 2D to 3D grid arrays:
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         i = i + 1
!         do bin1 = 1, num_bins
         do bin1 = 1, num_tiles
         !- If no counts occur for bin:
            if( fgrd_1d(i,bin1) == 0. ) then
               fgrd(c,r,bin1)    = 0.
               outvar1(c,r,bin1) = LDT_rc%udef
               outvar2(c,r,bin1) = LDT_rc%udef
         !- If counts occur for bin:
            else
               fgrd(c,r,bin1)    = fgrd_1d(i,bin1)
               outvar1(c,r,bin1) = var1_1d(i,bin1)
               outvar2(c,r,bin1) = var2_1d(i,bin1)
            endif
         end do
      enddo
   enddo

  end subroutine param_2dbin_areacalc

end module LDT_paramTileInputMod

