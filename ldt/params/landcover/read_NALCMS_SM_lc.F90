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
! !ROUTINE: read_NALCMS_SM_lc
! \label{read_NALCMS_SM_lc}
!
! !REVISION HISTORY:
!  16Jul2020: Kristi Arsenault; Added SnowModel topo-veg reader
!  03Jan2022 Kristi Arsenault; Added SM NALCMS landcover as separate field
!
! !INTERFACE:
subroutine read_NALCMS_SM_lc( n, num_types, fgrd, maskarray )

! !USES:
  use LDT_coreMod
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
!  This subroutine retrieves SnowModel's NALCMS landcover data files.
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
  real    :: vegcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
  real, allocatable :: read_topomap(:,:)   ! Read input parameters
  real, allocatable :: read_vegmap(:,:)   ! Read input parameters

! ____________________________________________________________________________________

  ! Assign additional land cover types, including generic water points: 
  LDT_rc%wetlandclass = 9
  LDT_rc%urbanclass   = 21
  LDT_rc%snowclass    = 19
  LDT_rc%glacierclass = 20
  LDT_rc%bareclass    = 18
  LDT_rc%waterclass   = 24

  vegcnt   = 0.
  vegtype  = float(LDT_rc%waterclass)
  maskarray= 0.0
  fgrd     = 0.0

  inquire(file=trim(LDT_rc%vfile(n)), exist=file_exists)
  if(.not.file_exists) then
     write(LDT_logunit,*) "[ERR] The landcover map: ",&
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

  ! Transform 2D landcover field to 3D
  do r = 1, LDT_rc%lnr(n)
     do c = 1, LDT_rc%lnc(n)
        if( vegtype(c,r) .le. 0 ) then
           vegtype(c,r) = float(LDT_rc%waterclass)
        endif
        if( (nint(vegtype(c,r)) .ne. LDT_rc%waterclass ) .and. &
            (nint(vegtype(c,r)) .ne. LDT_rc%udef)) then
           vegcnt(c,r,NINT(vegtype(c,r))) = 1.0
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
     call read_maskfile( n, vegtype, fgrd, maskarray )

  ! "CREATE" land mask and surface type fields (user-specified):
  elseif( LDT_rc%mask_type(n) == "create" ) then
     call create_maskfile( n, LDT_rc%nt, LDT_rc%lc_gridtransform(n), &
                 vegtype, vegcnt, maskarray )
  end if

! --

  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "[INFO] Done reading SnowModel-based NALCMS landcover file"

end subroutine read_NALCMS_SM_lc
