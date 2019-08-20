!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_paramOptCheckMod
!BOP
!
! !MODULE: LDT_paramOptCheckMod
!
! !DESCRIPTION:
!   This module contains a number of routines that check the LSM parameter
!    file option inputs from the ldt.config file.
!
! !REVISION HISTORY:
!  17 Jul 2012: KR Arsenault;  Initial Specification
!
! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_paramDataMod

  implicit none

contains

!BOP
! !ROUTINE: LDT_gridOptChecks
! \label{LDT_gridOptChecks}
!
! !INTERFACE:
 subroutine LDT_gridOptChecks( n, param_name, spatial_transform, &
                               param_proj, param_gridres )
! !USES:
   use LDT_coreMod,  only : LDT_rc, LDT_config
   use LDT_logMod,   only : LDT_verify, LDT_logunit, LDT_endrun

   implicit none

! !ARGUMENTS:
   integer,          intent(in) :: n
   character(len=*), intent(in) :: param_name
   character(50),    intent(in) :: spatial_transform
   character(50),    intent(in) :: param_proj
   real,             intent(in) :: param_gridres
!
! !DESCRIPTION:
!  Check resolutions of parameter and LIS grids and spatial transform.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[param_name]
!     parameter type
!   \item[spatial_transform]
!     Spatial transform of the parameter
!   \item[param_proj]
!     Parameter grid projection
!   \item[param_gridres]
!     Parameter grid resolution     
!   \end{description}
!EOP
   real :: lis_gridres
! ______________________________________________________________

!- Check for LIS domain projection option support:
   select case ( LDT_rc%lis_map_proj )

     case ( "latlon" )
     ! supported ...
       lis_gridres = LDT_rc%gridDesc(n,9)
     case ( "lambert" )
     ! supported ...
       lis_gridres = LDT_rc%gridDesc(n,8)
     case ( "gaussian" )
     ! supported
       lis_gridres = LDT_rc%gridDesc(n,9)
     case ( "hrap" )
     ! supported ...
       lis_gridres = LDT_rc%gridDesc(n,9)
     case ( "polar" )
     ! supported ...
       lis_gridres = LDT_rc%gridDesc(n,9)
     case ( "ease V2" )
     ! supported ...
       lis_gridres = LDT_rc%gridDesc(n,10)
     case ( "mercator" )
        write(LDT_logunit,*) "Grid_Check for LIS map projection: The ",trim(LDT_rc%lis_map_proj)
        write(LDT_logunit,*) "          projection for the targeted LIS run domain has not been"
        write(LDT_logunit,*) "          fully implemented or tested ..."
        write(LDT_logunit,*) " Program stopping ..."
        call LDT_endrun
     case ( "UTM" )
        write(LDT_logunit,*) "Grid_Check for LIS map projection: The ",trim(LDT_rc%lis_map_proj)
        write(LDT_logunit,*) "          projection for the targeted LIS run domain has not been"
        write(LDT_logunit,*) "          fully implemented or tested ..."
        write(LDT_logunit,*) " Program stopping ..."
        call LDT_endrun
   end select

! -------------------------------------------------------------

!- Check for parameter projection option support:
   select case ( param_proj )
     case ( "latlon" )
     ! supported ...
     case ( "hrap" )
     ! supported ...
     case ( "gaussian" )
     ! supported
     case ( "mercator" )
        write(LDT_logunit,*) "Grid_Check for param_proj: The ",trim(param_name)," parameter "
        write(LDT_logunit,*) "          projection is currently not supported. The following parameter"
        write(LDT_logunit,*) "          data grid projections include: "
        write(LDT_logunit,*) "       -- latlon, gaussian, hrap "
        write(LDT_logunit,*) " Program stopping ..."
        call LDT_endrun
     case ( "lambert" )
        write(LDT_logunit,*) "Grid_Check for param_proj: The ",trim(param_name)," parameter "
        write(LDT_logunit,*) "          projection is currently not supported. The following parameter"
        write(LDT_logunit,*) "          data grid projections include: "
        write(LDT_logunit,*) "       -- latlon, gaussian, hrap "
        write(LDT_logunit,*) " Program stopping ..."
        call LDT_endrun
     case ( "polar" )
        write(LDT_logunit,*) "Grid_Check for param_proj: The ",trim(param_name)," parameter "
        write(LDT_logunit,*) "          projection is currently not supported. The following parameter"
        write(LDT_logunit,*) "          data grid projections include: "
        write(LDT_logunit,*) "       -- latlon, gaussian, hrap "
        write(LDT_logunit,*) " Program stopping ..."
        call LDT_endrun
     case ( "UTM" )
        write(LDT_logunit,*) "Grid_Check for param_proj: The ",trim(param_name)," parameter "
        write(LDT_logunit,*) "          projection is currently not supported. The following parameter"
        write(LDT_logunit,*) "          data grid projections include: "
        write(LDT_logunit,*) "       -- latlon, gaussian, hrap "
        write(LDT_logunit,*) " Program stopping ..."
        call LDT_endrun
   end select

! -------------------------------------------------------------

!- Spatial transform checks:
! - Current checks below only work for latlon projections.  Future updates to come ...

   select case( spatial_transform )

     case ( "none" )
        if( param_gridres /= (lis_gridres/LDT_rc%lis_map_resfactor) ) then
           write(LDT_logunit,*) "Grid_Check for 'none':  The ",trim(param_name)," grid resolution"
           write(LDT_logunit,*) "           differs from the target LIS grid resolution. "
           write(LDT_logunit,*) "[WARN] For this spatial transform option, the parameter and LIS grid"
           write(LDT_logunit,*) "           resolutions must be the same."
        end if
     case ( "mode" )
        if( param_gridres >= lis_gridres/LDT_rc%lis_map_resfactor ) then
           write(LDT_logunit,*) "Grid_Check for 'mode':  The ",trim(param_name)," grid resolution"
           write(LDT_logunit,*) "           equals or is greater than the target LIS grid resolution."
           write(LDT_logunit,*) "[WARN] For this spatial transform option, the parameter resolution "
           write(LDT_logunit,*) "           must be smaller than the LIS grid resolution."
        end if
     case ( "average" )
        if( param_gridres >= lis_gridres/LDT_rc%lis_map_resfactor ) then
           write(LDT_logunit,*) "Grid_Check for 'average':  The ",trim(param_name)," grid resolution"
           write(LDT_logunit,*) "           equals or is greater than the target LIS grid resolution."
           write(LDT_logunit,*) "[WARN] For this spatial transform option, the parameter resolution "
           write(LDT_logunit,*) "           must be smaller than the LIS grid resolution."
        end if
    case ( "neighbor" )
        if( param_gridres < lis_gridres/LDT_rc%lis_map_resfactor ) then
           write(LDT_logunit,*) "Grid_Check for 'neighbor':  The ",trim(param_name)," grid resolution"
           write(LDT_logunit,*) "           is less than the target LIS grid resolution."
           write(LDT_logunit,*) "[WARN] For this spatial transform option, the parameter resolution "
           write(LDT_logunit,*) "           should be greater than or equal to the LIS grid resolution."
        end if
    case ( "bilinear" )
        if( param_gridres < lis_gridres/LDT_rc%lis_map_resfactor ) then
           write(LDT_logunit,*) "Grid_Check for 'bilinear':  The ",trim(param_name)," grid resolution"
           write(LDT_logunit,*) "           is less than the target LIS grid resolution."
           write(LDT_logunit,*) "[WARN] For this spatial transform option, the parameter resolution "
           write(LDT_logunit,*) "           should be greater than or equal to the LIS grid resolution."
        end if
     case ( "budget-bilinear" )
        if( param_gridres < lis_gridres/LDT_rc%lis_map_resfactor ) then
           write(LDT_logunit,*) "Grid_Check for 'budget-bilinear':  The ",trim(param_name)," grid resolution"
           write(LDT_logunit,*) "           is less than the target LIS grid resolution."
           write(LDT_logunit,*) "[WARN] For this spatial transform option, the parameter resolution "
           write(LDT_logunit,*) "           should be greater than or equal to the LIS grid resolution."
        end if
     case ( "tile" )
        if( trim(param_name) == "Landcover" ) then
           if( param_gridres > lis_gridres/LDT_rc%lis_map_resfactor ) then
             write(LDT_logunit,*) "Grid_Check for 'tile':  The ",trim(param_name)," grid resolution"
             write(LDT_logunit,*) "           is greater than the target LIS grid resolution. "
             write(LDT_logunit,*) "[WARN] For this spatial transform option, the parameter resolution must"
             write(LDT_logunit,*) "           be smaller or equal than the LIS grid resolution."
           endif
           if( param_gridres == lis_gridres/LDT_rc%lis_map_resfactor .and. &
               param_gridres == 0.01 ) then
             write(LDT_logunit,*) "Grid_Check for 'tile':  The ",trim(param_name)," grid resolution"
             write(LDT_logunit,*) "           is 1KM and equal to the target LIS grid resolution. "
             write(LDT_logunit,*) "[WARN] Please select 'none' for this situation."
           endif

        else
           if( param_gridres >= lis_gridres/LDT_rc%lis_map_resfactor ) then
             write(LDT_logunit,*) "Grid_Check for 'tile':  The ",trim(param_name)," grid resolution"
             write(LDT_logunit,*) "           equals or is greater than the target LIS grid resolution."
             write(LDT_logunit,*) "[WARN] For this spatial transform option, the parameter resolution "
             write(LDT_logunit,*) "           must be smaller than the LIS grid resolution."
           endif
        end if
     case default
        write(LDT_logunit,*) "Grid_Check: Spatial transform upscaling or downscaling option not "
        write(LDT_logunit,*) "             supported or recognized.  Please try one of the "
        write(LDT_logunit,*) "             following options (based on parameter chosen): "
        write(LDT_logunit,*) "          -- none "
        write(LDT_logunit,*) "          -- mode (upscaling) "
        write(LDT_logunit,*) "          -- average (upscaling) "
        write(LDT_logunit,*) "          -- tile (upscaling) "
        write(LDT_logunit,*) "          -- neighbor (upscaling/downscaling)"
        write(LDT_logunit,*) "          -- bilinear (downscaling)"
        write(LDT_logunit,*) "          -- budget-bilinear (downscaling) "
        write(LDT_logunit,*) "     Note:  Not all options work yet for all parameter types. "
        write(LDT_logunit,*) " Stopping ..."
        call LDT_endrun

   end select

 end subroutine LDT_gridOptChecks


!BOP
! !ROUTINE: LDT_LMLCOptChecks
! \label{LDT_LMLCOptChecks}
!
! !INTERFACE:
 subroutine LDT_LMLCOptChecks( param_name, lc_type, nt, spatial_transform )

! !USES:
   use LDT_coreMod,  only : LDT_rc, LDT_config
   use LDT_logMod,   only : LDT_verify, LDT_logunit, LDT_endrun

   implicit none

! !ARGUMENTS:
   character(len=*), intent(in) :: param_name
   character(50),    intent(in) :: lc_type
   integer,          intent(in) :: nt
   character(50),    intent(in) :: spatial_transform

   integer :: n
!
! !DESCRIPTION:
!   This subroutine checks options of landcover and landmask fields.
!EOP

   select case( param_name )

     case ( "Landcover" )

   !- Check on landcover index value ("vertical levels"):
     select case( lc_type )

      case( "UMD" )      
         if( nt .lt. 14 .or. nt .ge. 15 ) then
            write(LDT_logunit,*) "Param_Check: UMD has a minimum of 14 types (includes water)."
            write(LDT_logunit,*) "             Value is currently at: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
         endif
       case( "USGS" )
         if( nt .lt. 24 .or. nt .ge. 25 ) then 
            write(LDT_logunit,*) "Param_Check: USGS has a minimum of 24 types (includes water)."
            write(LDT_logunit,*) "             Value is currently at: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
         endif
       case("ECOCLIMAP2") 
         if( nt .lt. 12 .or. nt .ge. 13 ) then 
            write(LDT_logunit,*) "Param_Check: ECOCLIMAP2 has a minimum of 12 types (includes water)."
            write(LDT_logunit,*) "             Value is currently at: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
         endif
        case( "IGBPNCEP" )
          if( nt .lt. 20 .or. nt .ge. 21 ) then            
            write(LDT_logunit,*) "Param_Check: IGBPNCEP has a minimum of 20 types (includes water)."
            write(LDT_logunit,*) "             Value is currently at: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif
        case( "Bondville" )
          if( nt .lt. 20 .or. nt .ge. 21 ) then            
            write(LDT_logunit,*) "Param_Check: Bondville has a minimum of 20 types (includes water)."
            write(LDT_logunit,*) "             Value is currently at: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif
        case( "JULES_PFT" )
          if( nt .lt. 10) then            
            write(LDT_logunit,*) "Param_Check: UKMO_IGBP_PFT has a minimum of 10 types (includes water)."
            write(LDT_logunit,*) "             Value is currently at: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif
        case( "MOSAIC" )
          if( nt .lt. 7 .or. nt .ge. 8 ) then
            write(LDT_logunit,*) "Param_Check: MOSAIC has a minimum of 7 types (includes water)."
            write(LDT_logunit,*) "             Value is currently at: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif
        case( "ISA" )
          if( nt .lt. 13 .or. nt .ge. 14 ) then
            write(LDT_logunit,*) "Param_Check: ISA has a minimum of 13 types (includes water)."
            write(LDT_logunit,*) "             Value is currently at: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif
        case( "CLM45" )
          if( nt .lt. 36 .or. nt .ge. 37 ) then
            write(LDT_logunit,*) "Param_Check: CLM45 has a minimum of 36 types (includes water)."
            write(LDT_logunit,*) "             Value is currently at: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif
#if 0
        case( "CONSTANT" )
          if( nt .lt. 2 ) then
            write(LDT_logunit,*) "Param_Check: CONSTANT has a minimum of 2 types (includes water)."
            write(LDT_logunit,*) "             Value is currently at: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif
#endif
        case default 
          write(LDT_logunit,*) "Param_Check: Land classification selected, ",trim(lc_type),", currently not "
          write(LDT_logunit,*) "             supported or not entered correctly. Please try: "
          write(LDT_logunit,*) "          -- UMD "
          write(LDT_logunit,*) "          -- USGS "
          write(LDT_logunit,*) "          -- IGBPNCEP "
          write(LDT_logunit,*) "          -- JULES_PFT "
          write(LDT_logunit,*) "          -- MOSAIC "
          write(LDT_logunit,*) "          -- ISA "
          write(LDT_logunit,*) "          -- CLM45 "
          write(LDT_logunit,*) "          -- CONSTANT "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

   !- Spatial transform checks:
      select case( spatial_transform ) 

        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "mode" )
          write(LDT_logunit,*) "Param_Check: Dominant (mode) vegetation type being selected ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Neighboring vegetation type being selected ..."
        case ( "tile" )
          write(LDT_logunit,*) "Param_Check: Vegetation tiles are being created ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- tile "
          write(LDT_logunit,*) "              -- mode "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun

      end select

 !- Regional mask:
     case ( "Regional mask" )

   !- Spatial transform checks:
      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "mode" )
          write(LDT_logunit,*) "Param_Check: Dominant (mode) region index value being selected ..."
        case ( "average" )
          write(LDT_logunit,*) "Param_Check: Average region index value being selected ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Neighboring index value being selected ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- mode "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) "              -- average "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

   end select   ! end variable type

 end subroutine LDT_LMLCOptChecks


!BOP
! !ROUTINE: LDT_soilsOptChecks
! \label{LDT_soilsOptChecks}
!
! !INTERFACE:
 subroutine LDT_soilsOptChecks( n, param_name, soil_class, spatial_transform )

! !USES:
   use LDT_coreMod,  only : LDT_rc, LDT_config
   use LDT_logMod,   only : LDT_verify, LDT_logunit, LDT_endrun

   implicit none
! !ARGUMENTS:
   integer,          intent(in) :: n    ! nest
   character(len=*), intent(in) :: param_name
   character(50),    intent(in) :: soil_class
   character(50),    intent(in) :: spatial_transform
!
! !DESCRIPTION:
!   This subroutine checks options of the soil parameter 
!   fields.
!
   integer :: nt
   character(50) :: soilclass_check
!EOP

   nt = LDT_LSMparam_struc(n)%texture%num_bins
   if( index(soil_class,"STATSGO") > 0 ) then
      soilclass_check = "STATSGO"
   else
      soilclass_check = soil_class
   endif

   select case( param_name )

     case ( "Soil texture" )

   !- Check number of index values ("vertical levels"):
!      select case( soil_class )
      select case( soilclass_check )
        case( "STATSGO" )
          if( nt .lt. 16 ) then
            write(LDT_logunit,*) "Param_Check: STATSGO has a minimum of 16 types (includes water) for tiling."
            write(LDT_logunit,*) "             The current value is: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif

        case( "FAO" )

        case( "CONSTANT" )
          if( nt .lt. 16 ) then
            write(LDT_logunit,*) "Param_Check: USING CONSTANT SOIL TEXTURE VALUE WITH STATSGO AS INDEX RANGE ... "
            write(LDT_logunit,*) "      -- STATSGO has a minimum of 16 types (includes water) for tiling."
            write(LDT_logunit,*) "             The current value is: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif

        case( "ZOBLER_GFS")
          if (nt.ne.9) then 
            write(LDT_logunit,*) "Param_Check: ZOBLER_GFS has 9 types (includes water) for tiling."
            write(LDT_logunit,*) "             The current value is: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif

        case( "ISRIC")
          if (nt.ne.13) then 
            write(LDT_logunit,*) "Param_Check: ISRIC has 13 types (includes water) for tiling."
            write(LDT_logunit,*) "             The current value is: ",nt
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
          endif
        case( "Special" )

        case default
          write(LDT_logunit,*) "Param_Check: Soil classification, ",trim(soilclass_check),", selected "
          write(LDT_logunit,*) "             currently not supported or not entered correctly. Please try: "
          write(LDT_logunit,*) "          -- STATSGO (can include merged STATSGO+FAO)"
          write(LDT_logunit,*) "          -- FAO (which only includes FAO)"
          write(LDT_logunit,*) "          -- CONSTANT (which uses the STATSGO index range)"
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

   !- Spatial transform checks:
      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "mode" )
          write(LDT_logunit,*) "Param_Check: Dominant (mode) soil texture class being selected ..."
        case ( "tile" )
          write(LDT_logunit,*) "Param_Check: Tiled soil texture class being selected ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Nearest neighbor of soil texture class being selected ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- mode "
          write(LDT_logunit,*) "              -- tile "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

     case ( "Soil color" )

!      if( LDT_LSMparam_struc(n)%color%num_bins .lt. 9 ) then
!         write(LDT_logunit,*) "Param_Check: Soil color has a minimum of 9 types."
!         write(LDT_logunit,*) "             To process correctly, please change the value,",&
!                                             LDT_LSMparam_struc(n)%color%vlevels,","
!         write(LDT_logunit,*) " Stopping ..."
!         call LDT_endrun
!      endif

      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "average" )
          write(LDT_logunit,*) "Param_Check: Average soil color being estimated (not recommended) ..."
        case ( "mode" )
          write(LDT_logunit,*) "Param_Check: Dominant (mode) soil color class being selected ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Neighboring soil color class being selected ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- mode "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

     case ( "Slope type" )

!      if( LDT_LSMparam_struc(n)%slopetype%vlevels .lt. 9 ) then
!         write(LDT_logunit,*) "Param_Check: Slope type has 9 classes (including water)."
!         write(LDT_logunit,*) "             To process correctly, please change the value,",&
!                                            LDT_LSMparam_struc(n)%slopetype%vlevels
!         write(LDT_logunit,*) " Stopping ..."
!         call LDT_endrun
!      endif

       select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "mode" )
          write(LDT_logunit,*) "Param_Check: Dominant (mode) slope type being selected ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Estimate neareast neighbor slope type values ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options for slope type: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- mode "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

     case ( "Soils" )

      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "average" )
          write(LDT_logunit,*) "Param_Check: Average soil fraction being estimated ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Neighboring soil fraction being estimated ..."
        case ( "tile" )
          write(LDT_logunit,*) "Param_Check: Creating tiled soil fractions information ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- average "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) "              -- tile "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select
    
   end select

 end subroutine LDT_soilsOptChecks

!BOP
! !ROUTINE: LDT_albedoOptChecks
! \label{LDT_albedoOptChecks}
!
! !INTERFACE:
 subroutine LDT_albedoOptChecks( param_name, param_proj, spatial_transform )

! !USES:
   use LDT_coreMod, only : LDT_rc, LDT_config
   use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun

   implicit none
! !ARGUMENTS:
   character(len=*), intent(in) :: param_name
   character(50),    intent(in) :: param_proj
   character(50),    intent(in) :: spatial_transform
!
! !DESCRIPTION:
!   This subroutine checks options of the albedo fields.
!EOP
   select case( param_name )

     case ( "Albedo" )

      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "average" )
          write(LDT_logunit,*) "Param_Check: Average albedo values to be calculated ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Estimate neareast neighbor snow-free albedo values ..."
        case ( "bilinear" )
          write(LDT_logunit,*) "Param_Check: Estimate bilinearly interpolated snow-free albedo values ..."
        case ( "budget-bilinear" )
          write(LDT_logunit,*) "Param_Check: Estimate conservative interpolated snow-free values ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- average "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) "              -- bilinear "
          write(LDT_logunit,*) "              -- budget-bilinear "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

     case ( "Max snow albedo" )

      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "average" )
          write(LDT_logunit,*) "Param_Check: Average max snow albedo values to be calculated ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Estimate neareast neighbor max snow alb values ..."
        case ( "bilinear" )
          write(LDT_logunit,*) "Param_Check: Estimate bilinearly interpolated max snow alb values ..."
        case ( "budget-bilinear" )
          write(LDT_logunit,*) "Param_Check: Estimate conservative interpolated max snow alb values ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- average "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) "              -- bilinear "
          write(LDT_logunit,*) "              -- budget-bilinear "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
       end select

   end select
 end subroutine LDT_albedoOptChecks

!BOP
! !ROUTINE: LDT_gfracOptChecks
! \label{LDT_gfracOptChecks}
!
! !INTERFACE:
 subroutine LDT_gfracOptChecks( param_name, param_proj, spatial_transform )

! !USES:
   use LDT_coreMod, only : LDT_rc, LDT_config
   use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun

   implicit none
! !ARGUMENTS:
   character(len=*), intent(in) :: param_name
   character(50),    intent(in) :: param_proj
   character(50),    intent(in) :: spatial_transform
!
! !DESCRIPTION:
!   This subroutine checks options of the greenness fraction fields.
!EOP
   select case( param_name )

     case ( "Greenness" )

      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "average" )
          write(LDT_logunit,*) "Param_Check: Average greenness fraction values to be calculated ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Estimate neareast neighbor greenness fraction values ..."
        case ( "bilinear" )
          write(LDT_logunit,*) "Param_Check: Estimate bilinearly interpolated greenness fraction values ..."
        case ( "budget-bilinear" )
          write(LDT_logunit,*) "Param_Check: Estimate conservative interpolated greenness fraction values ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- average "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) "              -- bilinear "
          write(LDT_logunit,*) "              -- budget-bilinear "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

   end select  ! end param case
 end subroutine LDT_gfracOptChecks

!BOP
! !ROUTINE: LDT_topoOptChecks
! \label{LDT_topoOptChecks}
!
! !INTERFACE:
 subroutine LDT_topoOptChecks( param_name, param_proj, spatial_transform )

! !USES:
   use LDT_coreMod, only : LDT_rc, LDT_config
   use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun

   implicit none
! !ARGUMENTS:
   character(len=*), intent(in) :: param_name
   character(50),    intent(in) :: param_proj
   character(50),    intent(in) :: spatial_transform
!
! !DESCRIPTION:
!   This subroutine checks options of topographic fields.
!EOP
   select case( param_name )

     case ( "Topography" )

      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "average" )
          write(LDT_logunit,*) "Param_Check: Average topographic values to be calculated ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Neighboring topographic values to be calculated ..."
        case ( "bilinear" )
          write(LDT_logunit,*) "Param_Check: Using bilinear interpolation for topographic values ..."
        case ( "budget-bilinear" )
          write(LDT_logunit,*) "Param_Check: Using budget-bilinear interpolation for topographic values ..."
        case ( "tile" )
          write(LDT_logunit,*) "Param_Check: Tiling topographic values into fractions/averages ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "            -- none "
          write(LDT_logunit,*) "            -- average "
          write(LDT_logunit,*) "            -- neighbor "
          write(LDT_logunit,*) "            -- bilinear "
          write(LDT_logunit,*) "            -- budget-bilinear "
          write(LDT_logunit,*) "            -- tile "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

   end select
 end subroutine LDT_topoOptChecks

!BOP
! !ROUTINE: LDT_laisaiOptChecks
! \label{LDT_laisaiOptChecks}
!
! !INTERFACE:
 subroutine LDT_laisaiOptChecks( param_name, param_proj, spatial_transform )

! !USES:
   use LDT_coreMod, only : LDT_rc, LDT_config
   use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun

   implicit none
! !ARGUMENTS:
   character(len=*), intent(in) :: param_name
   character(50),    intent(in) :: param_proj
   character(50),    intent(in) :: spatial_transform
!
! !DESCRIPTION:
!   This subroutine checks options of LAI/SAI fields.
!EOP
   select case( param_name )

     case ( "LAI/SAI" )

      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "average" )
          write(LDT_logunit,*) "Param_Check: Average LAI/SAI values to be calculated ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Neighboring LAI/SAI values to be calculated ..."
        case ( "bilinear" )
          write(LDT_logunit,*) "Param_Check: Using bilinear interpolation for LAI/SAI values ..."
        case ( "budget-bilinear" )
          write(LDT_logunit,*) "Param_Check: Using budget-bilinear interpolation for LAI/SAI values ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- average "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) "              -- bilinear "
          write(LDT_logunit,*) "              -- budget-bilinear "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

    end select
 end subroutine LDT_laisaiOptChecks

!BOP
! !ROUTINE: LDT_noahparmsOptChecks
! \label{LDT_tbotOptChecks}
!
! !INTERFACE:
 subroutine LDT_noahparmsOptChecks( n, param_name, param_proj, &
      spatial_transform )

! !USES:
   use LDT_coreMod, only : LDT_rc, LDT_config
   use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun

   implicit none
! !ARGUMENTS:
   integer,          intent(in) :: n
   character(len=*), intent(in) :: param_name
   character(50),    intent(in) :: param_proj
   character(50),    intent(in) :: spatial_transform
!
! !DESCRIPTION:
!   This subroutine checks options of the bottom temperature field.
!EOP
   select case( param_name )

     case ( "Bottom temperature" )

      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "average" )
          write(LDT_logunit,*) "Param_Check: Average bottom temperature values to be calculated ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Estimate neareast neighbor bottom temperature values ..."
        case ( "bilinear" )
          write(LDT_logunit,*) "Param_Check: Estimate bilinearly interpolated bottom temperature values ..."
        case ( "budget-bilinear" )
          write(LDT_logunit,*) "Param_Check: Estimate conservative interpolated bottom temperature values ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- average "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) "              -- bilinear "
          write(LDT_logunit,*) "              -- budget-bilinear (conservative)"
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

   end select  ! end param case
 end subroutine LDT_noahparmsOptChecks


!BOP
! !ROUTINE: LDT_climateOptChecks
! \label{LDT_climateOptChecks}
!
! !INTERFACE:
 subroutine LDT_climateOptChecks( param_name, spatial_transform )

! !USES:
   use LDT_coreMod, only : LDT_rc, LDT_config
   use LDT_logMod,  only : LDT_verify, LDT_logunit, LDT_endrun

   implicit none
! !ARGUMENTS:
   character(len=*), intent(in) :: param_name
   character(50),    intent(in) :: spatial_transform
!
! !DESCRIPTION:
!   This subroutine checks options of the forcing climatology fields.
!EOP
   select case( param_name )

     case ( "Climate PPT" )

      select case( spatial_transform )
        case ( "none" )
          write(LDT_logunit,*) "Param_Check: NO spatial transform (read from native parameter grid)"
        case ( "average" )
          write(LDT_logunit,*) "Param_Check: Average Climate PPT values to be calculated ..."
        case ( "neighbor" )
          write(LDT_logunit,*) "Param_Check: Nearest neighbor Climate PPT values to be selected ..."
        case ( "bilinear" )
          write(LDT_logunit,*) "Param_Check: Bilinear interp used to calculate climate PPT values ..."
        case ( "budget-bilinear" )
          write(LDT_logunit,*) "Param_Check: Budget-bilinear interp used to calculate climate PPT values ..."
        case default
          write(LDT_logunit,*) "Param_Check: Spatial transform option selected currently not "
          write(LDT_logunit,*) "              supported or not entered correctly. Please try one "
          write(LDT_logunit,*) "              of the following options: "
          write(LDT_logunit,*) "              -- none "
          write(LDT_logunit,*) "              -- average "
          write(LDT_logunit,*) "              -- neighbor "
          write(LDT_logunit,*) "              -- bilinear "
          write(LDT_logunit,*) "              -- budget-bilinear "
          write(LDT_logunit,*) " Stopping ..."
          call LDT_endrun
      end select

   end select  ! end param case

 end subroutine LDT_climateOptChecks

end module LDT_paramOptCheckMod
