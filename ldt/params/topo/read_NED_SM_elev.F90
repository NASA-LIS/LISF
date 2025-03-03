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
! !ROUTINE: read_NED_SM_elev
! \label{read_NED_SM_elev}
!
! !REVISION HISTORY:
!  16Jul2020: Kristi Arsenault; Added SnowModel topo-veg reader
!  20Dec2021: Kristi Arsenault; Added SM NED elevation as separate field
!
! !INTERFACE:
subroutine read_NED_SM_elev( n, num_bins, fgrd, elevave )

! !USES:
  use LDT_coreMod
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_1dbin_areacalc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_bins
  real,    intent(out):: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out):: elevave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine retrieves SnowModel's NED elevation data files.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of bins (or bands)
!  \item[fgrd]
!   output grid fractions for elevation bins (or bands)
!  \item[elevave]
!   elevation average for each bin (or band)
!EOP

  integer :: ftn
  integer :: c, r, t
  logical :: file_exists
  real    :: elev(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
  real, allocatable :: read_topomap(:,:)   ! Read input parameters

! ____________________________________________________________________________________

  elev = 0.
  fgrd = 0.
  elevave = 0.

  inquire(file=trim(LDT_rc%elevfile(n)), exist=file_exists)
  if(.not.file_exists) then
     write(LDT_logunit,*) "[ERR] Elevation map ",&
            trim(LDT_rc%elevfile(n))," not found."
     call LDT_endrun
  endif

  select case ( LDT_rc%topo_gridtransform(n) )
    case( "none" ) 
      write(LDT_logunit,*) "[INFO] Reading NED (SnowModel-based) elevation file: ",&
            trim(LDT_rc%elevfile(n))
  case default
     write(LDT_logunit,*) "[ERR] This NED (SnowModel-based) elevation reader currently"
     write(LDT_logunit,*) "  only supports domains the same as the set LIS run domain." 
     call LDT_endrun
  end select

  ftn = LDT_getNextUnitNumber()
  open( ftn, file=LDT_rc%elevfile(n), &
        form="unformatted", access='direct',status='old', &
        recl=4*LDT_rc%gnc(n)*LDT_rc%gnr(n) )

  allocate(read_topomap(LDT_rc%gnc(n),LDT_rc%gnr(n)))

  ! Read in topographic map
  read(ftn,rec=1) ((read_topomap(c,r),c=1,LDT_rc%gnc(n)),r=1,LDT_rc%gnr(n))

  elev(1:LDT_rc%lnc(n),1:LDT_rc%lnr(n),1) = &
        read_topomap((LDT_ews_halo_ind(n,LDT_localPet+1)):(LDT_ewe_halo_ind(n,LDT_localPet+1)),&
                     (LDT_nss_halo_ind(n,LDT_localPet+1)):(LDT_nse_halo_ind(n,LDT_localPet+1)))
  deallocate(read_topomap)

  ! Transform 2D elevation field to 3D
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        do t = 1, num_bins
           ! Single elevation layer, write to first bin:
           if( t == 1 ) then
             fgrd(c,r,t)    = 1.0
             elevave(c,r,t) = elev(c,r,1)
           endif
        end do
     enddo
  enddo
 
  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "[INFO] Done reading SnowModel-based NED elevation file"

end subroutine read_NED_SM_elev
