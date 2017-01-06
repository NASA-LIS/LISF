! $Id$
!
! Earth System Modeling Framework
! Copyright 2002-2016, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!
!==============================================================================
#define FILENAME "beta_NUOPC_Base.F90"
!==============================================================================

#define RECONCILE_MEMORY_DEBUG_off

module beta_NUOPC_Base

  !-----------------------------------------------------------------------------
  ! Generic Coupler Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC

  implicit none
  
  private
  
  public beta_NUOPC_GetStateMemberLists
  public beta_NUOPC_AddNamespace

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
!BOP
! !IROUTINE: beta_NUOPC_GetStateMemberLists - Build lists of information of State members
! !INTERFACE:
  recursive subroutine beta_NUOPC_GetStateMemberLists(state, StandardNameList, &
    ConnectedList, DomainList, NamespaceList, itemNameList, fieldList, rc)
! !ARGUMENTS:
    type(ESMF_State),       intent(in)            :: state
    character(ESMF_MAXSTR), pointer, optional     :: StandardNameList(:)
    character(ESMF_MAXSTR), pointer, optional     :: ConnectedList(:)
    character(ESMF_MAXSTR), pointer, optional     :: DomainList(:)
    character(ESMF_MAXSTR), pointer, optional     :: NamespaceList(:)
    character(ESMF_MAXSTR), pointer, optional     :: itemNameList(:)
    type(ESMF_Field),       pointer, optional     :: fieldList(:)
    integer,                intent(out), optional :: rc
! !DESCRIPTION:
!   Construct lists containing the StandardNames, field names, and connected 
!   status of the fields in {\tt state}. Return this information in the
!   list arguments. Recursively parse through nested States.
!
!   All pointer arguments present must enter this method unassociated. On 
!   return, the deallocation of an associated pointer becomes the responsibility
!   of the caller.
!
!   The arguments are:
!   \begin{description}
!   \item[state]
!     The {\tt ESMF\_State} object to be queried.
!   \item[{[StandardNameList]}]
!     If present, return a list of the "StandardName" attribute of each member.
!   \item[{[ConnectedList]}]
!     If present, return a list of the "Connected" attribute of each member.
!   \item[{[DomainList]}]
!     If present, return a list of the "Domain" attribute of each member.
!   \item[{[NamespaceList]}]
!     If present, return a list of the namespace of each member.
!   \item[{[itemNameList]}]
!     If present, return a list of each member name.
!   \item[{[fieldList]}]
!     If present, return a list of the member fields.
!   \item[{[rc]}]
!     Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    integer           :: item, itemCount, fieldCount, stat, i
    type(ESMF_Field)  :: field
    character(ESMF_MAXSTR), allocatable     :: ll_itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable  :: stateitemtypeList(:)
    type(ESMF_State)                        :: nestedState
    character(ESMF_MAXSTR), pointer         :: l_StandardNameList(:)
    character(ESMF_MAXSTR), pointer         :: l_itemNameList(:)
    character(ESMF_MAXSTR), pointer         :: l_ConnectedList(:)
    character(ESMF_MAXSTR), pointer         :: l_DomainList(:)
    character(ESMF_MAXSTR), pointer         :: l_NamespaceList(:)
    type(ESMF_Field),       pointer         :: l_fieldList(:)
    character(ESMF_MAXSTR)                  :: domain
    character(ESMF_MAXSTR)                  :: namespace
    
    if (present(rc)) rc = ESMF_SUCCESS
    
    call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=FILENAME)) &
      return  ! bail out
          
    if (itemCount > 0) then
      allocate(ll_itemNameList(itemCount))
      allocate(stateitemtypeList(itemCount))
      call ESMF_StateGet(state, itemNameList=ll_itemNameList, &
        itemtypeList=stateitemtypeList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=FILENAME)) &
        return  ! bail out
        
      fieldCount = 0  ! reset
      do item=1, itemCount
        if (stateitemtypeList(item) == ESMF_STATEITEM_FIELD) then
          fieldCount = fieldCount + 1
        else if (stateitemtypeList(item) == ESMF_STATEITEM_STATE) then
          ! recursively parse the nested state
          nullify(l_StandardNameList)
          call ESMF_StateGet(state, itemName=ll_itemNameList(item), &
            nestedState=nestedState, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
          call beta_NUOPC_GetStateMemberLists(nestedState, l_StandardNameList, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
          if (associated(l_StandardNameList)) then
            fieldCount = fieldCount + size(l_StandardNameList)
            deallocate(l_StandardNameList)
          endif
        endif
      enddo
      
      if (present(StandardNameList)) then
        if (associated(StandardNameList)) then
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="StandardNameList must enter unassociated", &
            line=__LINE__, &
            file=FILENAME, &
            rcToReturn=rc)
          return  ! bail out
        else
          allocate(StandardNameList(fieldCount), stat=stat)
          if (ESMF_LogFoundAllocError(stat, msg="allocating StandardNameList", &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
        endif
      endif
      
      if (present(itemNameList)) then
        if (associated(itemNameList)) then
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="itemNameList must enter unassociated", &
            line=__LINE__, &
            file=FILENAME, &
            rcToReturn=rc)
          return  ! bail out
        else
          allocate(itemNameList(fieldCount), stat=stat)
          if (ESMF_LogFoundAllocError(stat, msg="allocating itemNameList", &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
        endif
      endif

      if (present(ConnectedList)) then
        if (associated(ConnectedList)) then
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="ConnectedList must enter unassociated", &
            line=__LINE__, &
            file=FILENAME, &
            rcToReturn=rc)
          return  ! bail out
        else
          allocate(ConnectedList(fieldCount), stat=stat)
          if (ESMF_LogFoundAllocError(stat, msg="allocating ConnectedList", &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
        endif
      endif

      if (present(DomainList)) then
        if (associated(DomainList)) then
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="DomainList must enter unassociated", &
            line=__LINE__, &
            file=FILENAME, &
            rcToReturn=rc)
          return  ! bail out
        else
          allocate(DomainList(fieldCount), stat=stat)
          if (ESMF_LogFoundAllocError(stat, msg="allocating DomainList", &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
        endif
      endif

      if (present(NamespaceList)) then
        if (associated(NamespaceList)) then
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="NamespaceList must enter unassociated", &
            line=__LINE__, &
            file=FILENAME, &
            rcToReturn=rc)
          return  ! bail out
        else
          allocate(NamespaceList(fieldCount), stat=stat)
          if (ESMF_LogFoundAllocError(stat, msg="allocating NamespaceList", &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
        endif
      endif

      if (present(fieldList)) then
        if (associated(fieldList)) then
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="fieldList must enter unassociated", &
            line=__LINE__, &
            file=FILENAME, &
            rcToReturn=rc)
          return  ! bail out
        else
          allocate(fieldList(fieldCount), stat=stat)
          if (ESMF_LogFoundAllocError(stat, msg="allocating fieldList", &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
        endif
      endif

      fieldCount = 1  ! reset

      do item=1, itemCount
        call ESMF_AttributeGet(state, name="Domain", value=domain, &
          defaultvalue="0", convention="NUOPC", purpose="Instance", &
          attnestflag=ESMF_ATTNEST_ON, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=FILENAME)) &
          return  ! bail out
        call NUOPC_GetAttribute(state, name="Namespace", &
          value=namespace, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=FILENAME)) &
          return  ! bail out
        if (stateitemtypeList(item) == ESMF_STATEITEM_FIELD) then
          call ESMF_StateGet(state, itemName=ll_itemNameList(item), &
            field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
          if (present(StandardNameList)) then
            call NUOPC_GetAttribute(field, name="StandardName", &
              value=StandardNameList(fieldCount), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=FILENAME)) &
              return  ! bail out
          endif
          if (present(itemNameList)) then
            itemNameList(fieldCount)=ll_itemNameList(item)
          endif
          if (present(ConnectedList)) then
            call NUOPC_GetAttribute(field, name="Connected", &
              value=ConnectedList(fieldCount), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=FILENAME)) &
              return  ! bail out
          endif
          if (present(DomainList)) then
            DomainList(fieldCount)=trim(domain)
          endif
          if (present(NamespaceList)) then
            NamespaceList(fieldCount)=trim(namespace)
          endif
          if (present(fieldList)) then
            fieldList(fieldCount)=field
          endif
          fieldCount = fieldCount + 1
        else if (stateitemtypeList(item) == ESMF_STATEITEM_STATE) then
          ! recursively parse the nested state
          nullify(l_StandardNameList)
          nullify(l_itemNameList)
          nullify(l_ConnectedList)
          nullify(l_DomainList)
          nullify(l_NamespaceList)
          nullify(l_fieldList)
          call ESMF_StateGet(state, itemName=ll_itemNameList(item), &
            nestedState=nestedState, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
          call beta_NUOPC_GetStateMemberLists(nestedState, &
            StandardNameList=l_StandardNameList, &
            itemNameList=l_itemNameList, &
            ConnectedList=l_ConnectedList, &
            DomainList=l_DomainList, &
            NamespaceList=l_NamespaceList, &
            fieldList=l_fieldList, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
          if (associated(l_StandardNameList)) then
            do i=1, size(l_StandardNameList)
              if (present(StandardNameList)) then
                StandardNameList(fieldCount) = l_StandardNameList(i)
              endif
              if (present(itemNameList)) then
                itemNameList(fieldCount) = l_itemNameList(i)
              endif
              if (present(ConnectedList)) then
                ConnectedList(fieldCount) = l_ConnectedList(i)
              endif
              if (present(DomainList)) then
                DomainList(fieldCount) = l_DomainList(i)
              endif
              if (present(NamespaceList)) then
                if (l_NamespaceList(i) /= "[UNLABELED_NESTED_STATE]") then
                  NamespaceList(fieldCount) = trim(namespace)//":"// &
                    trim(l_NamespaceList(i))
                else
                  NamespaceList(fieldCount) = trim(namespace)
                endif
              endif
              if (present(fieldList)) then
                fieldList(fieldCount) = l_fieldList(i)
              endif
              fieldCount = fieldCount + 1
            enddo
            deallocate(l_StandardNameList)
            deallocate(l_itemNameList)
            deallocate(l_ConnectedList)
            deallocate(l_DomainList)
            deallocate(l_NamespaceList)
            deallocate(l_fieldList)
          endif
        endif
      enddo
        
      deallocate(ll_itemNameList)
      deallocate(stateitemtypeList)
    endif
    
  end subroutine
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: beta_NUOPC_AddNamespace - Add a namespace and or domain to a State
! !INTERFACE:
  subroutine beta_NUOPC_AddNamespace(state, namespace, domain, &
    nestedStateName, nestedState, rc)
! !ARGUMENTS:
    type(ESMF_State), intent(inout)         :: state
    character(len=*), intent(in),  optional :: namespace
    character(len=*), intent(in),  optional :: domain
    character(len=*), intent(in),  optional :: nestedStateName
    type(ESMF_State), intent(out), optional :: nestedState
    integer,          intent(out), optional :: rc
! !DESCRIPTION:
!   Add a namespace to {\tt state}. Namespaces are implemented via nested
!   states. This creates a nested state inside of {\tt state}. The nested state
!   is returned as {\tt nestedState}. If provided, {\tt nestedStateName} will
!   be used to name the newly created nested state. The default name of the
!   nested state is equal to {\tt namespace}.
!
!   The arguments are:
!   \begin{description}
!   \item[state]
!     The {\tt ESMF\_State} object to which the namespace is added.
!   \item[namespace]
!     The namespace string.
!   \item[domain]
!     The domain identifier string.
!   \item[{[nestedStateName]}]
!     Name of the nested state. Defaults to {\tt namespace}.
!   \item[{[nestedState]}]
!     Optional return of the newly created nested state.
!   \item[{[rc]}]
!     Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    type(ESMF_State)        :: nestedS
    character(len=80)       :: nestedSName
    character(ESMF_MAXSTR)  :: s_namespace

    if (present(rc)) rc = ESMF_SUCCESS

    if (.not.present(namespace) .and. &
        .not.present(domain)) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg="Missing namespace or domain", &
        line=__LINE__, &
        file=FILENAME, &
        rcToReturn=rc)
      return  ! bail out
    endif

    if (present(nestedStateName)) then
      nestedSName = trim(nestedStateName)
    else
      nestedSName = trim(namespace)
    endif

    nestedS = ESMF_StateCreate(name=nestedSName, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    call beta_NUOPC_InitAttributesState(nestedS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    if (present(namespace)) then
      call NUOPC_SetAttribute(nestedS, name="Namespace", &
        value=trim(namespace), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
    else
      call NUOPC_SetAttribute(nestedS, name="Namespace", &
        value="[UNLABELED_NESTED_STATE]", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
    endif

    if (present(domain)) then
      call NUOPC_SetAttribute(nestedS, name="Domain", &
        value=trim(domain), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
    else
      call NUOPC_SetAttribute(nestedS, name="Domain", &
        value="0", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
    endif

    call ESMF_StateAdd(state, (/nestedS/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    if (present(nestedState)) &
      nestedState = nestedS

  end subroutine
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
!BOPI
! !IROUTINE: beta_NUOPC_InitAttribute - Initialize the NUOPC State Attributes
! !INTERFACE:
  ! call using generic interface: NUOPC_InitAttributes
  subroutine beta_NUOPC_InitAttributesState(state, rc)
! !ARGUMENTS:
    type(ESMF_state)                      :: state
    integer,      intent(out), optional   :: rc
! !DESCRIPTION:
!   Add the standard NUOPC State AttPack hierarchy to the State.
!
!   The highest level in the AttPack hierarchy will have convention="NUOPC" and
!   purpose="Instance".
!EOPI
  !-----------------------------------------------------------------------------
    ! local variables
    character(ESMF_MAXSTR)            :: attrList(3)

    if (present(rc)) rc = ESMF_SUCCESS

    ! Set up a customized list of Attributes to be added to the Fields
    attrList(1) = "Namespace"           ! namespace of this State
    attrList(2) = "FieldTransferPolicy" ! indicates to connectors to transfer/mirror fields:
                                        !    one of transferNone, transferAll
    attrList(3) = "Domain"              ! domain identifier of this State

    ! add Attribute packages
    call ESMF_AttributeAdd(state, convention="NUOPC", purpose="Instance", &
      attrList=attrList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    ! set Attributes to defaults
    call ESMF_AttributeSet(state, attrList(2), "transferNone", &
        convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

  end subroutine
  !-----------------------------------------------------------------------------

end module
