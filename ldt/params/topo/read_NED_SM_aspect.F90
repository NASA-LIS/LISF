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
! !ROUTINE: read_NED_SM_aspect
! \label{read_NED_SM_aspect}
!
! !REVISION HISTORY:
!  16Jul2020: Kristi Arsenault; Added SnowModel topo-veg reader
!  28Dec2021: Kristi Arsenault; Added SM NED aspect as separate field
!
! !INTERFACE:
subroutine read_NED_SM_aspect( n, num_bins, fgrd, aspectave )

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
  real,    intent(out):: aspectave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine retrieves SnowModel's NED elevation data files.
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
!   output grid fractions for aspect bins (or bands)
!  \item[aspectave]
!   aspect average for each bin (or band)
!EOP

  integer :: ftn
  integer :: c, r, t
  logical :: file_exists
  real, allocatable :: read_topomap(:,:)   ! Read input parameters
  real    :: elev(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: aspect(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! For aspect calculation:
  integer :: i,j
  real    :: pi, rad2deg
  real    :: deltax,deltay,deltaxy
  real    :: dzdx(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: dzdy(LDT_rc%lnc(n),LDT_rc%lnr(n))
! ____________________________________________________________________________________

  elev = 0.
  aspect = 0.
  fgrd = 0.
  aspectave = 0.

  inquire(file=trim(LDT_rc%elevfile(n)), exist=file_exists)
  if(.not.file_exists) then
     write(LDT_logunit,*) "[ERR] Elevation/aspect map ",&
           trim(LDT_rc%elevfile(n))," not found."
     call LDT_endrun
  endif

  select case ( LDT_rc%topo_gridtransform(n) )
    case( "none" ) 
      write(LDT_logunit,*) "[INFO] Reading NED (SnowModel-based) aspect file: ",&
            trim(LDT_rc%elevfile(n))
  case default
     write(LDT_logunit,*) "[ERR] This NED (SnowModel-based) aspect reader currently"
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

  ! Calculate the terrain aspect and azimuth.
  do i=1,LDT_rc%lnc(n)
     do j=1,LDT_rc%lnr(n)

        ! Some compilers will not allow dzdx and dzdy to both be 0.0 in
        !   the atan2 computation.
        !if (abs(dzdx(i,j)).lt.1e-10) dzdx(i,j) = 1e-10
        if (abs(dzdy(i,j)).lt.1e-10) then 
           dzdy(i,j) = 1e-10
        endif

        ! Compute the slope azimuth, making sure that north has zero
        !  azimuth.  Also note that for the Ryan wind rotation, the
        !  azimuth values must range from 0 to 360.
!        slope_az(i,j) = rad2deg * &
        aspect(i,j,1) = rad2deg * &
              (3.0 / 2.0 * pi - atan2(dzdy(i,j),dzdx(i,j)))
!        if (slope_az(i,j).ge.360.0) then
!           slope_az(i,j) = slope_az(i,j) - 360.0
        if (aspect(i,j,1).ge.360.0) then
           aspect(i,j,1) = aspect(i,j,1) - 360.0
        endif

     enddo
  enddo


  ! Transform 2D aspect field to 3D
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        do t = 1, num_bins
           ! Single aspect layer, write to first bin:
           if( t == 1 ) then
             fgrd(c,r,t)    = 1.0
             aspectave(c,r,t) = aspect(c,r,1)
           endif
        end do
     enddo
  enddo
 
  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "[INFO] Done reading SnowModel-based NED aspect file"

end subroutine read_NED_SM_aspect
