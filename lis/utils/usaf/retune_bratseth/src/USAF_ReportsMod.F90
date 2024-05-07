!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!
!
! MODULE: USAF_ReportsMod
!
! This contains code for storing select observation information and
! corresponding (interpolated) background and analysis values.
!
! REVISION HISTORY:
! 26 Oct 2020:  Eric Kemp.  Initial specification.
!
module USAF_ReportsMod

  ! Defaults
  implicit none
  private

  ! Public structure
  type Reports
     private
     integer :: nobs
     character(len=10), allocatable :: networks(:)
     character(len=10), allocatable :: platforms(:)
     character(len=10), allocatable :: yyyymmddhh(:)
     real, allocatable :: latitudes(:)
     real, allocatable :: longitudes(:)
     real, allocatable :: O(:)
     real, allocatable :: B(:)
     real, allocatable :: A(:)
  end type Reports
  public :: Reports

  ! Interface for constructor
  interface createReports
     module procedure newReports
  end interface createReports
  public :: createReports

  ! Public methods
  public :: newReports
  public :: destroyReports
  public :: appendToReports
  public :: getReport
  public :: getNobs
  public :: getListOfNetworksAndPlatforms
  public :: bcast_reports
contains

  !---------------------------------------------------------------------------
  ! Constructor
  function newReports(maxReports) result(this)

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: maxReports

    ! Result
    type(Reports) :: this

    ! Initialize structure
    this%nobs = 0
    allocate(this%networks(maxReports))
    allocate(this%platforms(maxReports))
    allocate(this%yyyymmddhh(maxReports))
    allocate(this%latitudes(maxReports))
    allocate(this%longitudes(maxReports))
    allocate(this%O(maxReports))
    allocate(this%B(maxReports))
    allocate(this%A(maxReports))

    this%networks(:) = "NULL"
    this%platforms(:) = "NULL"
    this%yyyymmddhh(:) = "NULL"
    this%latitudes(:) = 0
    this%longitudes(:) = 0
    this%O(:) = 0
    this%B(:) = 0
    this%A(:) = 0

  end function newReports

  !---------------------------------------------------------------------------
  ! Destructor
  subroutine destroyReports(this)

    ! Defaults
    implicit none

    ! Arguments
    type(Reports), intent(inout) :: this

    this%nobs = 0
    deallocate(this%networks)
    deallocate(this%platforms)
    deallocate(this%yyyymmddhh)
    deallocate(this%latitudes)
    deallocate(this%longitudes)
    deallocate(this%O)
    deallocate(this%B)
    deallocate(this%A)

  end subroutine destroyReports

  !---------------------------------------------------------------------------
  subroutine appendToReports(this, network, platform, yyyymmddhh, &
       latitude, longitude, O, B, A)

    ! Defaults
    implicit none

    ! Arguments
    type(Reports), intent(inout) :: this
    character(len=10), intent(in) :: network
    character(len=10), intent(in) :: platform
    character(len=10), intent(in) :: yyyymmddhh
    real, intent(in) :: latitude
    real, intent(in) :: longitude
    real, intent(in) :: O
    real, intent(in) :: B
    real, intent(in) :: A

    ! Local variables
    integer :: nobs

    ! Sanity check
    nobs = this%nobs
    if (nobs .eq. size(this%platforms, 1)) then
       write(6,*)'ERROR, not enough memory to append report!'
       write(6,*)'Current allocation size is ', size(this%platforms,1)
       write(6,*)'Increase allocation and try again!'
       stop 1
    end if

    ! Assign the values
    nobs = nobs + 1
    this%nobs = nobs
    this%networks(nobs) = network
    this%platforms(nobs) = platform
    this%yyyymmddhh(nobs) = yyyymmddhh
    this%latitudes(nobs) = latitude
    this%longitudes(nobs) = longitude
    this%O(nobs) = O
    this%B(nobs) = B
    this%A(nobs) = A

  end subroutine appendToReports

  !---------------------------------------------------------------------------
  ! General getter
  subroutine getReport(this, index, network, platform, yyyymmddhh, &
       latitude, longitude, O, B, A)

    ! Defaults
    implicit none

    ! Arguments
    type(Reports), intent(in) :: this
    integer, intent(in) :: index
    character(len=10), optional, intent(out) :: network
    character(len=10), optional, intent(out) :: platform
    character(len=10), optional, intent(out) :: yyyymmddhh
    real, optional, intent(out) :: latitude
    real, optional, intent(out) :: longitude
    real, optional, intent(out) :: O
    real, optional, intent(out) :: B
    real, optional, intent(out) :: A

    ! Local variable
    integer :: i

    ! Make sure at least one data member is requested
    if ( .not. present(network) .and. &
         .not. present(platform) .and. &
         .not. present(yyyymmddhh) .and. &
         .not. present(latitude) .and. &
         .not. present(longitude) .and. &
         .not. present(O) .and. &
         .not. present(B) .and. &
         .not. present(A)) then
       write(6,*)'ERROR, called getReport without requesting data member!'
       write(6,*)'Need to return network, platform, yyyymmddhh, ', &
            'latitude, longitude, O, B, and/or A'
       stop 1
    end if

    ! Sanity check station index
    i = index
    if (i .gt. this%nobs) then
       write(6,*)'ERROR, invalid station number!'
       write(6,*)'Asked for i = ', i, ' but max number is ', this%nobs
       stop 1
    end if

    if (present(network)) then
       network = this%networks(i)
    end if
    if (present(platform)) then
       platform = this%platforms(i)
    end if
    if (present(yyyymmddhh)) then
       yyyymmddhh = this%yyyymmddhh(i)
    end if
    if (present(latitude)) then
       latitude = this%latitudes(i)
    end if
    if (present(longitude)) then
       longitude = this%longitudes(i)
    end if
    if (present(O)) then
       O = this%O(i)
    end if
    if (present(B)) then
       B = this%B(i)
    end if
    if (present(A)) then
       A = this%A(i)
    end if

  end subroutine getReport

  !---------------------------------------------------------------------------
  integer function getNobs(this)
    implicit none
    type(Reports), intent(in) :: this
    getNobs = this%nobs
  end function getNobs

  !---------------------------------------------------------------------------
  subroutine getListOfNetworksAndPlatforms(this, networks, platforms)
    implicit none
    type(Reports), intent(in) :: this
    character(len=10), allocatable, intent(inout) :: networks(:)
    character(len=10), allocatable, intent(inout) :: platforms(:)
    integer :: i
    if (this%nobs .eq. 0) return
    allocate(networks(this%nobs))
    allocate(platforms(this%nobs))
    networks(:) = "NULL"
    platforms(:) = "NULL"
    do i = 1, this%nobs
       networks(i) = this%networks(i)
       platforms(i) = this%platforms(i)
    end do ! i
  end subroutine getListOfNetworksAndPlatforms

  subroutine bcast_reports(this, myid, ierr)
    use mpi
    implicit none
    type(Reports), intent(inout) :: this
    integer, intent(in) :: myid
    integer, intent(inout) :: ierr
    integer :: nobs
    ierr = 0
    nobs = 0
    if (myid .eq. 0) then
       nobs = getNobs(this)
    end if
    call mpi_bcast(nobs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (nobs .eq. 0) return
    this%nobs = nobs
    if (ierr .ne. MPI_SUCCESS) return
    call mpi_bcast(this%networks, nobs, MPI_CHARACTER, 0, MPI_COMM_WORLD, &
         ierr)
    if (ierr .ne. MPI_SUCCESS) return
    call mpi_bcast(this%platforms, nobs, MPI_CHARACTER, 0, MPI_COMM_WORLD, &
         ierr)
    if (ierr .ne. MPI_SUCCESS) return
    call mpi_bcast(this%yyyymmddhh, nobs, MPI_CHARACTER, 0, &
         MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) return
    call mpi_bcast(this%latitudes, nobs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) return
    call mpi_bcast(this%longitudes, nobs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) return
    call mpi_bcast(this%O, nobs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) return
    call mpi_bcast(this%B, nobs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) return
    call mpi_bcast(this%A, nobs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) return
  end subroutine bcast_reports
end module USAF_ReportsMod
