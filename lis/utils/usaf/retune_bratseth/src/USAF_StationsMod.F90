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
! MODULE: USAF_stationsMod
!
! This contains code for storing station latitudes, longitudes, and platforms.
!
! REVISION HISTORY:
! 26 Oct 2020:  Eric Kemp.  Initial specification.
!
module USAF_stationsMod

  ! Defaults
  implicit none
  private

  ! Public structure
  type Stations
     private
     integer :: nstations
     character(len=10), allocatable :: platforms(:)
     real, allocatable :: latitudes(:)
     real, allocatable :: longitudes(:)
  end type Stations
  public :: Stations

  ! Interface for constructor
  interface createStations
     module procedure newStations
  end interface createStations
  public :: createStations

  ! Public methods
  public :: newStations
  public :: destroyStations
  public :: appendToStations
  public :: getLatLonOfPlatform
  public :: getNstations
  public :: getListOfPlatforms
  public :: great_circle_distance

contains

  !---------------------------------------------------------------------------
  ! Constructor
  function newStations(max_stations) result(this)

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: max_stations

    ! Result
    type(Stations) :: this

    ! Initialize structure
    this%nstations = 0
    allocate(this%platforms(max_stations))
    allocate(this%latitudes(max_stations))
    allocate(this%longitudes(max_stations))

    this%platforms(:) = "NULL"
    this%latitudes(:) = 0
    this%longitudes(:) = 0

  end function newStations

  !---------------------------------------------------------------------------
  ! Destructor
  subroutine destroyStations(this)

    ! Defaults
    implicit none

    ! Arguments
    type(Stations), intent(inout) :: this

    this%nstations = 0
    deallocate(this%platforms)
    deallocate(this%latitudes)
    deallocate(this%longitudes)
  end subroutine destroyStations

  !---------------------------------------------------------------------------
  ! Add new station to the structure.
  subroutine appendToStations(this, platform, latitude, longitude)

    ! Defaults
    implicit none

    ! Arguments
    type(Stations), intent(inout) :: this
    character(len=10), intent(in) :: platform
    real, intent(in) :: latitude
    real, intent(in) :: longitude

    ! Local variables
    integer :: nstations
    integer :: j

    ! See if platform is already in structure
    do j = 1, this%nstations
       if (this%platforms(j) .eq. platform) return
    end do ! j

    nstations = this%nstations
    if (nstations .eq. size(this%platforms,1)) then
       write(6,*)'ERROR, not enough memory for appending station!'
       stop 1
    end if

    nstations = nstations + 1
    this%nstations = nstations
    this%platforms(nstations) = platform
    this%latitudes(nstations) = latitude
    this%longitudes(nstations) = longitude

  end subroutine appendToStations

  !---------------------------------------------------------------------------
  ! Get lat and lon of specified platform
  subroutine getLatLonOfPlatform(this, platform, latitude, longitude)

    ! Defaults
    implicit none

    ! Arguments
    type(Stations), intent(in) :: this
    character(len=10), intent(in) :: platform
    real, intent(out) :: latitude
    real, intent(out) :: longitude

    ! Local variables
    logical :: found
    integer :: i

    ! Find the station
    found = .false.
    do i = 1, this%nstations
       if (this%platforms(i) .eq. platform) then
          found = .true.
          exit
       end if
    end do ! i
    if (.not. found) then
       write(6,*)"ERROR, cannot find platform ", platform
       stop 1
    end if

    latitude = this%latitudes(i)
    longitude = this%longitudes(i)

  end subroutine getLatLonOfPlatform

  !---------------------------------------------------------------------------
  function getNstations(this) result(Nstations)
    implicit none
    type(Stations),intent(in) :: this
    integer :: Nstations
    Nstations = this%Nstations
  end function getNstations

  !---------------------------------------------------------------------------
  subroutine getListOfPlatforms(this, platforms)
    implicit none
    type(Stations), intent(in) :: this
    character(len=10), allocatable, intent(out) :: platforms(:)
    integer :: i
    if (this%nstations .eq. 0) return
    allocate(platforms(this%nstations))
    platforms(:) = "NULL"
    do i = 1, this%nstations
       platforms(i) = this%platforms(i)
    end do ! i
  end subroutine getListOfPlatforms

  !---------------------------------------------------------------------------
  real function great_circle_distance(lat1, lon1, lat2, lon2)

    ! Defaults
    implicit none

    ! Arguments
    real, intent(in) :: lat1, lon1, lat2, lon2

    ! Local variables
    double precision :: radlat1, radlon1, radlat2, radlon2
    double precision :: pi
    double precision :: deg2rad
    double precision :: lon_abs_diff
    double precision :: central_angle
    double precision :: term1, term2, term3

    pi = 4d0*atan(1d0)
    deg2rad = pi / 180d0
    radlat1 = dble(lat1)*deg2rad
    radlon1 = dble(lon1)*deg2rad
    radlat2 = dble(lat2)*deg2rad
    radlon2 = dble(lon2)*deg2rad
    lon_abs_diff = abs(radlon1 - radlon2)
    term1 = cos(radlat2)*sin(lon_abs_diff)
    term1 = term1*term1
    term2 = (cos(radlat1)*sin(radlat2)) - &
         (sin(radlat1)*cos(radlat2)*cos(lon_abs_diff))
    term2 = term2*term2
    term3 = (sin(radlat1)*sin(radlat2)) + &
         (cos(radlat1)*cos(radlat2)*cos(lon_abs_diff))
    central_angle = atan2( sqrt(term1 + term2) , term3 )
    great_circle_distance = real(6381000d0 * central_angle)
  end function great_circle_distance

end module USAF_stationsMod
