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
! !ROUTINE: read_NED_SM_curvature
! \label{read_NED_SM_curvature}
!
! !REVISION HISTORY:
!  16Jul2020: Kristi Arsenault; Added SnowModel topo-veg reader
!  28Dec2021: Kristi Arsenault; Added SM NED curvature as separate field
!
! !INTERFACE:
subroutine read_NED_SM_curvature( n, curv )

! !USES:
  use LDT_coreMod
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_1dbin_areacalc
#if ( defined SPMD )
  use mpi
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real,    intent(out):: curv(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!
!  This subroutine retrieves SnowModel's NED elevation data files.
!  The code in this routine is from Glen Liston's MicroMet code
!   described in Liston and Elder (2006) MicroMet paper.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[curv]
!   curvature average 
!EOP

  integer :: ftn
  integer :: c, r, t
  integer :: ierr
  logical :: file_exists
  real, allocatable :: read_topomap(:,:)   ! Read input parameters
  real    :: elev(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: curvature(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! For curvature calculation:
  integer :: i, j
  real    :: pi, rad2deg
  real    :: deltax,deltay,deltaxy
  real    :: curve_max, curve_max_all
  integer :: inc

! The curvature is used as part of the wind model.  Define a length
!   scale that the curvature calculation will be performed on.  This
!   has units of meters, and should be approximately one-half the
!   wavelength of the topographic features within the domain.
!!! (recommended default value: curve_len_scale = 500.0)
  real, parameter ::  curve_len_scale = 300.0

! ____________________________________________________________________________________

  elev = 0.
  curvature = 0.
  curv = 0.

  inquire(file=trim(LDT_rc%elevfile(n)), exist=file_exists)
  if(.not.file_exists) then
     write(LDT_logunit,*) "[ERR] Elevation/curvature map ",&
            trim(LDT_rc%elevfile(n))," not found."
     call LDT_endrun
  endif

  select case ( LDT_rc%topo_gridtransform(n) )
    case( "none" ) 
      write(LDT_logunit,*) "[INFO] Reading NED (SnowModel-based) curvature file: ",&
            trim(LDT_rc%elevfile(n))
  case default
     write(LDT_logunit,*) "[ERR] This NED (SnowModel-based) curvature reader currently"
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

  ! Convert the length scale to an appropriate grid increment.
  inc = max(1,nint(curve_len_scale/deltaxy))

  ! Compute the curvature.
  do i=1,LDT_rc%lnc(n)     ! nx
     do j=1,LDT_rc%lnr(n)  ! ny

        curvature(i,j,1) = (4.0 * elev(i,j) - &
          elev(max(1,i-inc),max(1,j-inc)) - &
          elev(min(LDT_rc%lnc(n),i+inc),min(LDT_rc%lnr(n),j+inc)) - &
          elev(min(LDT_rc%lnc(n),i+inc),max(1,j-inc)) - &
          elev(max(1,i-inc),min(LDT_rc%lnr(n),j+inc))) / &
          (sqrt(2.0) * 16.0 * real(inc) * deltaxy) + &
          (4.0 * elev(i,j) - &
          elev(min(LDT_rc%lnc(n),i+inc),j) - elev(max(1,i-inc),j) - &
          elev(i,min(LDT_rc%lnr(n),j+inc)) - elev(i,max(1,j-inc))) / &
          (16.0 * real(inc) * deltaxy)
     enddo
  enddo

  ! Scale the curvature such that the max abs(curvature) has a value
  !   of abs(0.5).  Include a 1 mm curvature in curve_max to prevent
  !   divisions by zero in flat terrain where the curvature is zero.
  curve_max = 0.0 + 0.001
  do j=1,LDT_rc%lnr(n)
     do i=1,LDT_rc%lnc(n)
        curve_max = max(curve_max,abs(curvature(i,j,1)))
     enddo
  enddo

#if (defined SPMD)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(curve_max, curve_max_all, 1,&
           MPI_REAL, MPI_MAX,&
           mpi_comm_world, ierr)
  curve_max = curve_max_all
#endif

  do j=1,LDT_rc%lnr(n)
     do i=1,LDT_rc%lnc(n)
        curvature(i,j,1) = curvature(i,j,1) / (2.0 * curve_max)
     enddo
  enddo

  curv = curvature

#if 0
  ! Transform 2D curvature field to 3D
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        do t = 1, num_bins
           ! Single curvature layer, write to first bin:
           if( t == 1 ) then
             fgrd(c,r,t)    = 1.0
             curvave(c,r,t) = curvature(c,r,1)
           endif
        end do
     enddo
  enddo
#endif
 
  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "[INFO] Done reading SnowModel-based NED curvature file"

end subroutine read_NED_SM_curvature
