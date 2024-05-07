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
! !ROUTINE: read_NALCMS_SM_IGBPNCEP_lc
! \label{read_NALCMS_SM_IGBPNCEP_lc}
!
! !REVISION HISTORY:
!  16Jul2020: Kristi Arsenault; Added SnowModel topo-veg reader
!  03Jan2022 Kristi Arsenault; Added SM NALCMS landcover as separate field
!  13Jul2022 Kristi Arsenault; Added SM NALCMS landcover with IGBP-NCEP classes
!
! !INTERFACE:
subroutine read_NALCMS_SM_IGBPNCEP_lc( n, num_types, fgrd, maskarray )

! !USES:
  use LDT_coreMod
  use LDT_paramDataMod
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(inout) :: num_types
  real,    intent(inout):: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
  real,    intent(inout):: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine reads in SnowModel's NALCMS landcover-derived
!   SM vegtype classes and converts to the IGBP/NCEP-modified classes
!   for use by Noah and NoahMP LSMs.
!
!  Final mapped IGBP-NCEP classes include,
!  (A) New IGBP_MODIS_BU+tundra Landcover Class Legend:
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
! _________________________
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_types]
!   number of bins (or bands)
!  \item[fgrd]
!   output grid fractions for landcover bins (or bands)
!  \item[maskarray]
!   landcover average for each bin (or band)
!EOP

  integer :: ftn
  integer :: c, r, t
  logical :: file_exists
  real    :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: vegtype2(LDT_rc%lnc(n),LDT_rc%lnr(n))  ! Mapped IGBP/NCEP classes
  real    :: vegcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
  real, allocatable :: read_topomap(:,:)   ! Read input parameters
  real, allocatable :: read_vegmap(:,:)    ! Read input parameters
! ____________________________________________________________________________________

  ! Assign additional land cover types, including generic water points: 
  LDT_rc%wetlandclass = 11
  LDT_rc%urbanclass   = 13
  LDT_rc%snowclass    = 15
  LDT_rc%glacierclass = 15
  LDT_rc%bareclass    = 16
  LDT_rc%waterclass   = 17

  vegcnt   = 0.
  vegtype  = float(LDT_rc%waterclass)
  vegtype2  = float(LDT_rc%waterclass)
  maskarray= 0.0
  fgrd     = 0.0

  inquire(file=trim(LDT_rc%vfile(n)), exist=file_exists)
  if(.not.file_exists) then
     write(LDT_logunit,*) "[ERR] The SM-based NALCMS landcover map: ",&
           trim(LDT_rc%vfile(n))," not found."
     call LDT_endrun
  endif

  select case ( LDT_rc%lc_gridtransform(n) )
    case( "none" ) 
      write(LDT_logunit,*) "[INFO] Reading NALCMS (SnowModel-based) landcover file: ",&
            trim(LDT_rc%vfile(n))
  case default
     write(LDT_logunit,*) "[ERR] This NALCMS (SnowModel-based) landcover reader currently"
     write(LDT_logunit,*) "  only supports domains the same as the set LIS run domain." 
     call LDT_endrun
  end select

  ftn = LDT_getNextUnitNumber()
  open( ftn, file=LDT_rc%vfile(n), &
        form="unformatted", access='direct',status='old', &
        recl=4*LDT_rc%gnc(n)*LDT_rc%gnr(n) )

  ! Read in topographic map
  allocate(read_topomap(LDT_rc%gnc(n),LDT_rc%gnr(n)))
  read(ftn,rec=1) ((read_topomap(c,r),c=1,LDT_rc%gnc(n)),r=1,LDT_rc%gnr(n))
  deallocate(read_topomap)

  ! Read in landcover map
  allocate(read_vegmap(LDT_rc%gnc(n),LDT_rc%gnr(n)))
  read(ftn,rec=2) ((read_vegmap(c,r),c=1,LDT_rc%gnc(n)),r=1,LDT_rc%gnr(n))
  vegtype(1:LDT_rc%lnc(n),1:LDT_rc%lnr(n)) = &
        read_vegmap((LDT_ews_halo_ind(n,LDT_localPet+1)):(LDT_ewe_halo_ind(n,LDT_localPet+1)),&
                     (LDT_nss_halo_ind(n,LDT_localPet+1)):(LDT_nse_halo_ind(n,LDT_localPet+1)))
  deallocate(read_vegmap)

  write(LDT_logunit,*) "[INFO] SnowModel-NALCMS: Mapping SnowModel classes to IGBP-NCEP"

  ! Transform 2D landcover field to 3D
  do r = 1, LDT_rc%lnr(n)
     do c = 1, LDT_rc%lnc(n)

        ! Map SnowModel classes to the IGBP/NCEP-modified classes:
        !  SnowModel Classes   ==>> IGBP/NCEP-modified classes 
        if ( vegtype(c,r) == 1.0 )  vegtype2(c,r) = 1.0
        if ( vegtype(c,r) == 2.0 )  vegtype2(c,r) = 4.0
        if ( vegtype(c,r) == 3.0 )  vegtype2(c,r) = 5.0
!        if ( vegtype(c,r) == 4.0 )  vegtype2(c,r) = 1.0
!        if ( vegtype(c,r) == 5.0 )  vegtype2(c,r) = 1.0
        if ( vegtype(c,r) == 6.0 )  vegtype2(c,r) = 6.0
!        if ( vegtype(c,r) == 7.0 )  vegtype2(c,r) = 1.0
        if ( vegtype(c,r) == 8.0 )  vegtype2(c,r) = 7.0
!        if ( vegtype(c,r) == 9.0 )  vegtype2(c,r) = 1.0
        if ( vegtype(c,r) == 10.0 )  vegtype2(c,r) = 18.0
!        if ( vegtype(c,r) == 11.0 )  vegtype2(c,r) = 1.0
        if ( vegtype(c,r) == 12.0 )  vegtype2(c,r) = 10.0
!        if ( vegtype(c,r) == 13.0 )  vegtype2(c,r) = 1.0
        if ( vegtype(c,r) == 14.0 )  vegtype2(c,r) = 19.0
!        if ( vegtype(c,r) == 15.0 )  vegtype2(c,r) = 1.0
!        if ( vegtype(c,r) == 16.0 )  vegtype2(c,r) = 1.0
        if ( vegtype(c,r) == 17.0 )  vegtype2(c,r) = 11.0
        if ( vegtype(c,r) == 18.0 )  vegtype2(c,r) = 16.0
        if ( vegtype(c,r) == 19.0 )  vegtype2(c,r) = 17.0
        if ( vegtype(c,r) == 20.0 )  vegtype2(c,r) = 15.0
        if ( vegtype(c,r) == 21.0 )  vegtype2(c,r) = 13.0
        if ( vegtype(c,r) == 22.0 )  vegtype2(c,r) = 12.0
!        if ( vegtype(c,r) == 23.0 )  vegtype2(c,r) = 1.0
        if ( vegtype(c,r) == 24.0 )  vegtype2(c,r) = 17.0

        if( vegtype(c,r) == 4.0 .or. &
            vegtype(c,r) == 5.0 .or. &
            vegtype(c,r) == 7.0 .or. &
            vegtype(c,r) == 9.0 .or. &
            vegtype(c,r) == 11.0 .or. &
            vegtype(c,r) == 13.0 .or. &
            vegtype(c,r) == 15.0 .or. &
            vegtype(c,r) == 16.0 .or. &
            vegtype(c,r) == 23.0 ) then

            print *, "[WARN] We have SM class vegtypes in this domain that are "
            print *, " not mapped to an IGBP/NCEP class !!! "
            print *, c, r, vegtype(c,r)
        endif

        if( vegtype2(c,r) .le. 0 ) then
           vegtype2(c,r) = float(LDT_rc%waterclass)
        endif
        if( (nint(vegtype2(c,r)) .ne. LDT_rc%waterclass ) .and. &
            (nint(vegtype2(c,r)) .ne. LDT_rc%udef)) then
           vegcnt(c,r,NINT(vegtype2(c,r))) = 1.0
        endif
     enddo
  end do

  ! Estimate fraction of grid (fgrid) represented by vegetation type::
  call param_index_fgrdcalc( n, LDT_rc%lc_proj, LDT_rc%lc_gridtransform(n), &
       LDT_rc%waterclass, LDT_rc%nt, vegcnt, fgrd )
 
! -------------------------------------------------------------------
!    CREATE OR READ-IN (OR IMPOSE) LAND MASK FILE AND CREATE
!    SURFACE MAP
! -------------------------------------------------------------------

  ! "READ-IN" land mask file, if user-specified:
  if( LDT_rc%mask_type(n) == "readin" ) then
     write(LDT_logunit,*) "[WARN] The SnowModel-based NALCMS landcover map reader "
     write(LDT_logunit,*) "  has not been tested with reading in and imposing an"
     write(LDT_logunit,*) "  outside landmask at this time. Please use with caution."
     call read_maskfile( n, vegtype2, fgrd, maskarray )

  ! "CREATE" land mask and surface type fields (user-specified):
  elseif( LDT_rc%mask_type(n) == "create" ) then
     call create_maskfile( n, LDT_rc%nt, LDT_rc%lc_gridtransform(n), &
                 vegtype2, vegcnt, maskarray )
  end if

  ! Temporary option added for NoahMP401+SM coupling: (KRA)
  if( LDT_rc%allmaskland == 1 ) then
     LDT_LSMparam_struc(n)%landmask%value = 1
     LDT_LSMparam_struc(n)%dommask%value = 1
  endif
  ! Temporary option (KRA)

! --

  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "[INFO] Done reading SnowModel-based NALCMS landcover file"

end subroutine read_NALCMS_SM_IGBPNCEP_lc
