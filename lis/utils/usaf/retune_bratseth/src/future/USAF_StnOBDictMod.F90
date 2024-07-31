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
! MODULE: USAF_StnOBDictMod
!
! This contains cone for storing and retrieving observations in a "dictionary"
! i.e., a hash table indexed by station names.
!
! REVISION HISTORY:
! 21 Dec 2020: Eric Kemp. Initial specification.
!
!------------------------------------------------------------------------------

module USAF_StnOBDictMod

  ! Defaults
  implicit none
  private

  ! Private structure -- a node with a O,B pair in a linked list
  type OBNode_t
     private
     real :: O  ! Observed value
     real :: B  ! Background value
     type(OBNode_t), pointer :: nextOBNode
  end type OBNode_t

  ! Private structure -- a node with a key and a O,B linked list, itself within
  ! a linked list
  type KeyNode_t
     private
     character(len=10) :: key
     type(OBNode_t), pointer :: firstOBNode
     type(KeyNode_t), pointer :: nextKeyNode
  end type KeyNode_t

  ! Public class storing lists of O,B pairs of different stations, stored by
  ! station name.
  type StnOBDict_t
     private
     logical :: empty
     type(KeyNode_t), pointer :: firstKeyNode
   contains
     procedure :: new => USAF_stnOBDict_new
     procedure :: destroy => USAF_stnOBDict_destroy
     procedure :: insert => USAF_stnOBDict_insert
     procedure :: getNumKeys => USAF_stnOBDict_getNumKeys
     procedure :: getMeanOMBs => USAF_stnOBDict_getMeanOMBs
  end type StnOBDict_t
  public :: StnOBDict_t

  ! Public methods
  public :: USAF_stnOBDict_new
  public :: USAF_stnOBDict_destroy
  public :: USAF_stnOBDict_insert
  public :: USAF_stnOBDict_getNumKeys
  public :: USAF_stnOBDict_getMeanOMBs

contains

  ! Constructor
  subroutine USAF_stnOBDict_new(this)
    implicit none
    class(StnOBDict_t), intent(inout) :: this
    this%empty = .true.
    nullify(this%firstKeyNode)
  end subroutine USAF_stnOBDict_new

  ! Destructor
  subroutine USAF_stnOBDict_destroy(this)

    ! Defaults
    implicit none

    ! Arguments
    class(stnOBDict_t), intent(inout) :: this

    ! Locals
    type(KeyNode_t), pointer :: firstKeyNode, keyNode, nextKeyNode

    if (this%empty) then
       ! This should only occur if the stnOBDict_t was created but never
       ! populated.
       return
    end if

    nullify(firstKeyNode, keyNode, nextKeyNode)

    firstKeyNode => this%firstKeyNode
    call delete_obnode_list(firstKeyNode)

    if (associated(firstKeyNode%nextKeyNode)) then
       keyNode => firstKeyNode%nextKeyNode
       nullify(firstKeyNode%nextKeyNode)
       do
          call delete_obnode_list(keyNode)
          if (associated(keyNode%nextKeyNode)) then
             nextKeyNode => keyNode%nextKeyNode
             deallocate(keyNode)
             keyNode => nextKeyNode
          else ! Last keyNode in list
             deallocate(keyNode)
             nullify(keyNode)
             nullify(nextKeyNode)
             exit
          end if
       end do
    end if
    deallocate(firstKeyNode)
    nullify(firstKeyNode)

    this%empty = .true.
  end subroutine USAF_stnOBDict_destroy

  ! Private subroutine. Delete OBnode linked list for a particular key.
  subroutine delete_obnode_list(keyNode)

    ! Defaults
    implicit none

    ! Arguments
    type(KeyNode_t), pointer, intent(inout) :: keyNode

    ! Locals
    type(OBNode_t), pointer :: firstOBNode, OBnode, nextOBNode

    nullify(firstOBNode, OBnode, nextOBNode)

    ! Loop through all OBNode_t linked list members
    firstOBNode => keyNode%firstOBnode
    if (associated(firstOBnode%nextOBNode)) then
       OBnode => firstOBnode%nextOBNode
       nullify(firstOBNode%nextOBNode)
       do
          if (associated(OBnode%nextOBNode)) then
             nextOBnode => OBnode%nextOBNode
             deallocate(OBnode)
             OBnode => nextOBNode
          else ! Last OBNode in list
             deallocate(OBnode)
             nullify(OBnode)
             nullify(nextOBnode)
             exit
          end if
       end do
    end if
    deallocate(firstOBnode)

  end subroutine delete_obnode_list

  ! Insert new O-B pair into the tnOBDict_t object
  subroutine USAF_stnOBDict_insert(this, key, O, B)

    use esmf

    ! Defaults
    implicit none

    ! Arguments
    class(StnOBDict_t), intent(inout) :: this
    character(len=10), intent(in) :: key
    real, intent(in) :: O
    real, intent(in) :: B

    ! Locals
    type(KeyNode_t), pointer :: keyNode, newKeyNode

    nullify(keyNode, newKeyNode)

    if (this%empty) then
       this%empty = .false.
       call create_keynode_pointer(this%firstKeyNode, key, O, B)
    else
       ! Search the linked list of keyNodes to see if we already have this
       ! key
       keyNode => this%firstKeyNode
       do
          if (trim(keyNode%key) .eq. trim(key)) then
             !call ESMF_LogWrite("Adding new OB pair for "//trim(key), &
             !     ESMF_LOGMSG_INFO)
             !call ESMF_LogFlush()
             call append_to_ob_list(keyNode, O, B)
             exit
          else if (associated(keyNode%nextKeyNode)) then
             keyNode => keyNode%nextKeyNode
          else
             call create_keynode_pointer(newKeyNode, key, O, B)
             keyNode%nextKeyNode => newKeyNode
             exit
          end if
       end do
    end if
  end subroutine USAF_stnOBDict_insert

  ! Private subroutine.  Creates and populates a single KeyNode pointer
  subroutine create_keynode_pointer(keyNode, key, O, B)
    implicit none
    type(KeyNode_t), pointer, intent(inout) :: keyNode
    character(len=10), intent(in) :: key
    real, intent(in) :: O
    real, intent(in) :: B
    allocate(keyNode)
    keyNode%key = key
    allocate(keyNode%firstOBNode)
    keyNode%firstOBNode%O = O
    keyNode%firstOBNode%B = B
    nullify(keyNode%firstOBNode%nextOBNode)
    nullify(keyNode%nextKeyNode)
  end subroutine create_keynode_pointer

  ! Private subroutine.  Appends O,B pair into linked list associated with
  ! provided keyNode_t object.
  subroutine append_to_ob_list(keyNode, O, B)

    ! Defaults
    implicit none

    ! Arguments
    type(KeyNode_t), pointer, intent(inout) :: keyNode
    real, intent(in) :: O
    real, intent(in) :: B

    ! Locals
    type(OBNode_t), pointer :: OBNode

    nullify(OBNode)

    if (.not. associated(keyNode%firstOBNode)) then
       allocate(keyNode%firstOBNode)
       keyNode%firstOBNode%O = O
       keyNode%firstOBNode%B = B
       nullify(keyNode%firstOBNode%nextOBNode)
    else
       OBNode => keyNode%firstOBNode
       do
          if (associated(OBNode%nextOBNode)) then
             OBNode => OBNode%nextOBNode
          else
             allocate(OBNode%nextOBNode)
             OBNode%nextOBNode%O = O
             OBNode%nextOBNode%B = B
             nullify(OBNode%nextOBNode%nextOBNode)
             exit
          end if
       end do
    end if
  end subroutine append_to_ob_list

  ! Public method.  Return total number of keys (gage names)
  function USAF_StnOBDict_getNumKeys(this) result(numKeys)
    implicit none
    class(StnOBDict_t), intent(in) :: this
    integer :: numKeys
    class(KeyNode_t), pointer :: keyNode
    numKeys = 0
    if (this%empty) return
    keyNode => this%firstKeyNode
    numKeys = 1
    do
       if (associated(keyNode%nextKeyNode)) then
          keyNode => keyNode%nextKeyNode
          numKeys = numKeys + 1
       else
          exit
       end if
    end do
  end function USAF_StnOBDict_getNumKeys

  ! Public method.  Calculate and return mean OMBs for each station, along
  ! with list of stations.
  subroutine USAF_StnOBDict_getMeanOMBs(this, nstns, stationNames, meanOMBs)

    ! Defaults
    implicit none

    ! Arguments
    class(StnOBDict_t), intent(in) :: this
    integer, intent(out) :: nstns
    character(len=10), allocatable, intent(out) :: stationNames(:)
    real, allocatable, intent(out) :: meanOMBs(:)

    ! Locals
    real, allocatable :: O_values(:)
    real, allocatable :: B_values(:)
    integer, allocatable :: counts(:)
    type(KeyNode_t), pointer :: keyNode
    type(obNode_t), pointer :: obNode
    integer :: k

    nullify(keyNode, obNode)

    nstns = this%getNumKeys()
    allocate(stationNames(nstns))
    stationNames = "NULL"
    allocate(meanOMBs(nstns))
    meanOMBs = 0
    allocate(O_values(nstns))
    O_values = 0
    allocate(B_values(nstns))
    B_values = 0
    allocate(counts(nstns))
    counts = 0

    keyNode => this%firstKeyNode
    k = 1
    do
       stationNames(k) = keyNode%key
       obNode => keyNode%firstOBNode
       do
          O_values(k) = O_values(k) + obNode%O
          B_values(k) = B_values(k) + obNode%B
          counts(k) = counts(k) + 1
          if (.not. associated(obNode%nextOBNode)) exit
          obNode => obNode%nextOBNode
       end do
       if (.not. associated(keyNode%nextKeyNode)) exit
       keyNode => keyNode%nextKeyNode
       k = k + 1
    end do

    do k = 1, nstns
       if (counts(k) > 0) then
          meanOMBs(k) = (O_values(k) - B_values(k)) / real(counts(k))
       end if
    end do

    deallocate(O_values)
    deallocate(B_values)
    deallocate(counts)

  end subroutine USAF_StnOBDict_getMeanOMBs
end module USAF_StnOBDictMod
