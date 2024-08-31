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
! !MODULE: Empty_obsMod
! 
! !DESCRIPTION: 
!   
! !REVISION HISTORY: 
!  09 Jul 09    Sujay Kumar;   Initial Specification
! 
module Empty_obsMod
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: Empty_obs_setup
  public :: read_Empty_em_obsdata
  public :: write_EmptyObsdata
  public :: reset_Empty_obsdata
  public :: Empty_setupobspred
  public :: Empty_getpeobspred

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------

contains
!BOP
! 
! !ROUTINE: Empty_obs_setup
! \label{Empty_obs_setup}
! 
! !INTERFACE: 
  subroutine Empty_obs_setup(Obs_State)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecGrid
    use LIS_timeMgrMod, only : LIS_calendar
    use LIS_logMod, only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obs_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!
!   The arguments are: 
!   \begin{description}
!    \item[Obj\_Space]   observation/Objective space object 
!   \end{description}
!EOP
    integer                   ::  n 
    character(len=LIS_CONST_PATH_LEN) ::  emptyobsdir
    integer             :: status 

!    allocate(Empty_obs_struc(LIS_rc%nnest))
    emptyobsdir='Empty'
    !Read lis.config entries
    do n=1,LIS_rc%nnest
       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            emptyobsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),'Data Update Status',&
            .false., rc=status)
       call LIS_verify(status)
    enddo

  end subroutine Empty_obs_setup

  subroutine read_Empty_em_obsdata(Obj_Space)
     ! !ARGUMENTS: 
    type(ESMF_State)    :: Obj_Space

    !Empty--used to support MCSIM for which there are no obs to deal with

  end subroutine read_Empty_em_obsdata

  subroutine write_EmptyObsdata(Obj_Space)
    ! !ARGUMENTS: 
    type(ESMF_State)    :: Obj_Space


    !Empty--used to support MCSIM for which there are no obs to deal with

  end subroutine write_EmptyObsdata

  subroutine reset_Empty_obsdata()

    !Empty--used to support MCSIM for which there are no obs to deal with

  end subroutine reset_Empty_obsdata

subroutine Empty_setupobspred(OBSPred)
! !USES:
  use LIS_coreMod,      only : LIS_vecTile
  use LIS_logMod,       only : LIS_verify

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: OBSPred
!
! !DESCRIPTION:
!  
! 
!EOP
  integer                :: n
  type(ESMF_ArraySpec)   :: realarrspec
  type(ESMF_Field)       :: emField
  integer                :: status

  n = 1

!!$  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
!!$       rc=status)
!!$  call LIS_verify(status)
!!$
!!$!   call ESMF_ArraySpecSet(realarrspec, rank=2,typekind=ESMF_TYPEKIND_R4,&
!!$!        rc=status)
!!$!   call LIS_verify(status)
!!$
!  emField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecTile(n), &
!      name="Empty", rc=status)
!  call LIS_verify(status)
  
!  call ESMF_StateAdd(OBSPred,emField,rc=status)
!  call LIS_verify(status)

end subroutine Empty_setupobspred

subroutine Empty_getpeobspred(Obj_Func)
! !USES:
!!  use LIS_coreMod, only : LIS_rc, LIS_domain
!!  use LIS_soilsMod,  only : LIS_soils
!!  use LIS_logMod,       only : LIS_verify
!!!  use CRTM2_EMMod  !, only           : Empty_ctl
!!use Empty_Mod, only : cmem3_struc
!!
!!  implicit none

!ARGUMENTS: 
  type(ESMF_State)       :: Obj_Func

!!!
! !DESCRIPTION:
!!!  
!!! 
!EOP
!!$  integer                :: n
!!$  type(ESMF_Field)       :: emField
!!$!  real, allocatable          :: em_data(:,:)
!!$  real, allocatable          :: em_data(:)  ! just focusing on one channel now
!!$  integer                :: t,col,row
!!$  integer                :: i
!!$  integer                :: status
!!$  integer                :: k
!!$  integer, parameter     :: cnrs_numchannels=7  !19V,19H,22V,37V,37H,85V,85H
!!$
!!$  n = 1
!!$  call ESMF_StateGet(Obj_Func,"Emissivity estimate",emField,rc=status)
!!$  call LIS_verify(status)
!!$
!!$  call ESMF_FieldGet(emField,localDE=0,farray=em_data,rc=status)
!!$  call LIS_verify(status)

!   allocate(em_data(LIS_rc%nch(n),cnrs_numchannels))
!   allocate(em_data(LIS_rc%nch(n)))

!   do t=1,LIS_rc%nch(n)
!   enddo

end subroutine Empty_getpeobspred

end module Empty_obsMod
