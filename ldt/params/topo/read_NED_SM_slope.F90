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
! !ROUTINE: read_NED_SM_slope
! \label{read_NED_SM_slope}
!
! !REVISION HISTORY:
!  16Jul2020: Kristi Arsenault; Added SnowModel topo-veg reader
!  20Dec2021: Kristi Arsenault; Added SM NED slope as separate field
!
! !INTERFACE:
subroutine read_NED_SM_slope( n, num_bins, fgrd, slopeave )

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
  real,    intent(out):: slopeave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine retrieves SnowModel's NED slope data files.
!  The code in this routine is from Glen Liston's MicroMet code
!   described in Liston and Elder (2006) MicroMet paper.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of bins (or bands)
!  \item[fgrd]
!   output grid fractions for slope bins (or bands)
!  \item[slopeave]
!   slope average for each bin (or band)
!EOP

  integer :: ftn
  integer :: c, r, t
  logical :: file_exists
  real, allocatable :: read_topomap(:,:)   ! Read input parameters
!  real    :: elev(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
  real    :: elev(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: slope(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! For slope calculation:
  integer :: i,j
  real    :: pi, rad2deg
  real    :: deltax,deltay,deltaxy
  real    :: dzdx(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: dzdy(LDT_rc%lnc(n),LDT_rc%lnr(n))
! ____________________________________________________________________________________

  elev = 0.
  slope = 0.
  fgrd = 0.
  slopeave = 0.

  inquire(file=trim(LDT_rc%elevfile(n)), exist=file_exists)
  if(.not.file_exists) then
     write(LDT_logunit,*) "[ERR] Elevation/slope map ",&
           trim(LDT_rc%elevfile(n))," not found."
     call LDT_endrun
  endif

  select case ( LDT_rc%topo_gridtransform(n) )
    case( "none" ) 
      write(LDT_logunit,*) "[INFO] Reading NED (SnowModel-based) slope file: ",&
            trim(LDT_rc%elevfile(n))
  case default
     write(LDT_logunit,*) "[ERR] This NED (SnowModel-based) slope reader currently"
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

  elev(1:LDT_rc%lnc(n),1:LDT_rc%lnr(n)) = &
        read_topomap((LDT_ews_halo_ind(n,LDT_localPet+1)):(LDT_ewe_halo_ind(n,LDT_localPet+1)),&
                     (LDT_nss_halo_ind(n,LDT_localPet+1)):(LDT_nse_halo_ind(n,LDT_localPet+1)))

  deallocate(read_topomap)


  ! From MicroMet calculations:
  pi = 2.0 * acos(0.0)
  rad2deg = 180.0 / pi

  ! Compute the average grid increment.
  deltax = LDT_rc%topo_gridDesc(n,8) * 1000.  ! Convert to meters
  deltay = LDT_rc%topo_gridDesc(n,9) * 1000.  ! Convert to meters
  deltaxy = 0.5 * (deltax + deltay)

  ! Find dzdx.
  do j=1,LDT_rc%lnr(n) ! ny
     dzdx(1,j) = (elev(2,j) - elev(1,j)) / deltax
     do i=2,LDT_rc%lnc(n)-1 ! nx
        dzdx(i,j) = (elev(i+1,j) - elev(i-1,j)) / (2.0 * deltax)
     enddo
     dzdx(LDT_rc%lnc(n),j) = (elev(LDT_rc%lnc(n),j) - elev(LDT_rc%lnc(n)-1,j)) / deltax
  enddo

  ! Find dzdy.
  do i=1,LDT_rc%lnc(n) !nx
     dzdy(i,1) = (elev(i,2) - elev(i,1)) / deltay
     do j=2,LDT_rc%lnr(n)-1 ! ny
        dzdy(i,j) = (elev(i,j+1) - elev(i,j-1)) / (2.0 * deltay)
     enddo
     dzdy(i,LDT_rc%lnr(n)) = (elev(i,LDT_rc%lnr(n)) - elev(i,LDT_rc%lnr(n)-1)) / deltay
  enddo

  ! Calculate the terrain slope and azimuth.
  do i=1,LDT_rc%lnc(n)
     do j=1,LDT_rc%lnr(n)

        ! Some compilers will not allow dzdx and dzdy to both be 0.0 in
        !   the atan2 computation.
        !if (abs(dzdx(i,j)).lt.1e-10) dzdx(i,j) = 1e-10
        if (abs(dzdy(i,j)).lt.1e-10) then 
           dzdy(i,j) = 1e-10
        endif

        ! Compute the slope of the terrain.
!        terrain_slope(i,j) = rad2deg * &
        slope(i,j,1) = rad2deg * &
            atan(sqrt(dzdx(i,j)*dzdx(i,j) + dzdy(i,j)*dzdy(i,j)))
     enddo
  enddo

  ! Transform 2D slope field to 3D
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        do t = 1, num_bins
           ! Single slope layer, write to first bin:
           if( t == 1 ) then
             fgrd(c,r,t)    = 1.0
             slopeave(c,r,t) = slope(c,r,1)
           endif
        end do
     enddo
  enddo
 
  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "[INFO] Done reading SnowModel-based NED slope file"

end subroutine read_NED_SM_slope
