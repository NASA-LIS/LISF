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
! MODULE: USAF_GridHashMod
!
! This contains code for storing and retrieving observation data in a
! "gridded hash table," a 2d array of linked lists storing the observations
! in each grid box.
!
! REVISION HISTORY:
! 26 Oct 2020:  Eric Kemp.  Initial specification.
!
module USAF_GridHashMod

  ! Defaults
  implicit none
  private

  ! Private structure -- a node in a linked list
  type ObindexNode
     private
     integer :: obindex
     type(ObindexNode), pointer :: next
  end type ObindexNode

  ! Public structure -- Defines a hash table, where the "hash" is the
  ! 2d grid coordinate.
  type GridHash
     private
     integer :: imax
     integer :: jmax
     real :: dlon
     real :: dlat
     type(ObindexNode), pointer :: lists(:, :)
  end type GridHash
  public :: GridHash

  ! Interface for constructor
  interface createGridHash
     module procedure newGridHash
  end interface createGridHash
  public :: createGridHash

  ! Public methods
  public :: newGridHash
  public :: destroyGridHash
  public :: insertIntoGridHash
  public :: getObindexVectorFromGridHash
  public :: createIJForGridHash

contains

  ! ** Constructor **
  function newGridHash(imax, jmax, dlon, dlat) result(this)

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: imax
    integer, intent(in) :: jmax
    real, intent(in) :: dlon
    real, intent(in) :: dlat

    ! Result
    type(GridHash) :: this

    ! Local variables
    integer :: i, j

    this%imax = imax
    this%jmax = jmax
    allocate(this%lists(imax, jmax))

    this%dlon = dlon
    this%dlat = dlat

    do j = 1, jmax
       do i = 1, imax
          this%lists(i, j)%obindex = 0
          nullify(this%lists(i, j)%next)
       end do ! i
    end do ! j

  end function newGridHash

  ! ** Destructor **
  subroutine destroyGridHash(this)

    ! Defaults
    implicit none

    ! Arguments
    type(GridHash), intent(inout) :: this

    ! Local variables
    integer :: imax, jmax
    type(ObindexNode), pointer :: node, next, first
    integer :: i, j

    imax = this%imax
    jmax = this%jmax
    nullify(node, next, first)

    do j = 1, jmax
       do i = 1, imax
          ! Very first node in list must be preserved until we deallocate
          ! the array. But all nodes beyond that can be deallocated one by
          ! one.
          first => this%lists(i, j)
          first%obindex = 0
          if (.not. associated(first%next)) cycle
          node => first%next
          nullify(first%next)
          do
             if (associated(node%next)) then
                next => node%next
                deallocate(node)
                node => next
             else ! Last node in list
                deallocate(node)
                nullify(node)
                nullify(next)
                exit
             end if
          end do
       end do ! i
    end do ! j

    deallocate(this%lists)

  end subroutine destroyGridHash

  ! ** Insert a new data member into the hash **
  subroutine insertIntoGridHash(this, i, j, obindex)

    ! Defaults
    implicit none

    ! Arguments
    type(GridHash), intent(inout) :: this
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: obindex

    ! Local variables
    type(ObindexNode), pointer :: node, new_node

    nullify(node, new_node)
    node => this%lists(i, j)
    do
       ! List member has no data.  Just assign.
       if (node%obindex .eq. 0) then
          node%obindex = obindex
          exit
       else
          ! At current end of list.  Append.
          if (.not. associated(node%next)) then
             allocate(new_node)
             new_node%obindex = obindex
             nullify(new_node%next)
             node%next => new_node
             exit
          else
             ! Not at end of list.
             node => node%next
          end if
       end if
    end do

  end subroutine insertIntoGridHash

  ! ** Get all obindex values in current GridHash i, j
  subroutine getObindexVectorFromGridHash(this, i, j, nobs, obindexVector)

    ! Defaults
    implicit none

    ! Arguments
    type(GridHash), intent(in) :: this
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(out) :: nobs
    integer, allocatable, intent(out) :: obindexVector(:)

    ! Local variables
    type(ObindexNode), pointer :: node
    integer :: k

    nullify(node)

    ! See how many obs are in this grid box.
    nobs = 0
    if (this%lists(i,j)%obindex .eq. 0) return
    node => this%lists(i, j)
    nobs = nobs + 1
    do
       if (.not. associated(node%next)) then
          exit
       else
          node => node%next
          nobs = nobs + 1
       end if
    end do

    ! Now collect the index values
    allocate(obindexVector(nobs))
    node => this%lists(i, j)
    k = 1
    obindexVector(k) = node%obindex
    do
       if (.not. associated(node%next)) then
          exit
       else
          node => node%next
          k = k + 1
          obindexVector(k) = node%obindex
       end if
    end do

  end subroutine getObindexVectorFromGridHash

  ! ** Create i,j based on lat/lon
  subroutine createIJForGridHash(this, lat, lon, i, j)

    ! Defaults
    implicit none

    ! Arguments
    type(GridHash), intent(in) :: this
    real, intent(in) :: lat
    real, intent(in) :: lon
    integer, intent(out) :: i
    integer, intent(out) :: j

    ! Local variables
    real :: tmp_lat, tmp_lon

    tmp_lat = lat + 90

    tmp_lon = lon
    if (tmp_lon .lt. 0) then
       tmp_lon = tmp_lon + 360
    end if

    i = int(tmp_lon / this%dlon) + 1
    j = int(tmp_lat / this%dlat) + 1

    if (j .lt. 1 .or. i .lt. 1) then
       !print*, 'lon, tmp_lon, j: ', lon, tmp_lon, j
       !print*, 'lat, tmp_lat, i: ', lat, tmp_lat, i
    end if

  end subroutine createIJForGridHash

end module USAF_GridHashMod
