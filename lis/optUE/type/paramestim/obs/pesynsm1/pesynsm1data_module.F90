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
! !MODULE: pesynsm1data_module
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle synthetic soil moisture test data generated 
!   from a previous LIS simulation
!   
! !REVISION HISTORY: 
!  02Feb08    Sujay Kumar;   Initial Specification
! 
module pesynsm1data_module
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

  PUBLIC :: pesynsm1data_setup
contains
!BOP
! 
! !ROUTINE: pesynsm1data_setup
! \label{pesynsm1data_setup}
! 
! !INTERFACE: 
  subroutine pesynsm1data_setup(Obj_Space)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecGrid
    use LIS_logMod, only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obj_Space(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for soil moisture
!    assimilation
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                   ::  n 
    integer                   ::  status
    integer                   ::  nobstypes
    integer                   ::  i 
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field), allocatable ::  obsField(:)
    character(len=LIS_CONST_PATH_LEN) ::  synsmobsdir
    character*100,  allocatable   ::  vname(:)
    character(len=LIS_CONST_PATH_LEN) ::  objspaceAttribFile(LIS_rc%nnest)
    integer                   ::  ftn

    
    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"Syn SM data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,synsmobsdir,&
            rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obj_Space(n),"Data Directory",&
            synsmobsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obj_Space(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Syn SM observations attributes file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,objspaceAttribFile(n),rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Syn SM number of observation types:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,nobstypes,&
            rc=status)
       call LIS_verify(status)
   
       allocate(obsField(nobstypes))
       allocate(vname(nobstypes))

       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(objspaceAttribFile(n)),status='old')
       read(ftn,*)
       
       do i=1,nobstypes
          read(ftn,fmt='(a100)') vname(i)
          obsField(i) = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecGrid(n), &
               name = trim(vname(i)), rc=status)
          call LIS_verify(status)

          call ESMF_StateAdd(Obj_Space(n),(/obsField(i)/),rc=status)
          call LIS_verify(status)
       enddo
       call LIS_releaseUnitNumber(ftn)
       deallocate(vname)
    enddo

    write(LIS_logunit,*) 'Created the States to hold the observations data'
    
  end subroutine pesynsm1data_setup
  
end module pesynsm1data_module
