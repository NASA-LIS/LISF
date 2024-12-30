!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#define FILENAME "lis_nuopc_flags"
#define MODNAME "lis_nuopc_flags"
#include "LIS_NUOPC_Macros.h"

!> @file LIS_NUOPC_Flags.F90 LIS NUOPC Cap configuration flags
module LIS_NUOPC_Flags
!BOP
!
! !MODULE: LIS_NUOPC_Flags
!
! !DESCRIPTION:
!   This module contains configuration setting flags
!
! !REVISION HISTORY:
!  2022May02    Dan Rosen  Added
!
! !USES:
  use ESMF, only: ESMF_UtilStringUpperCase, ESMF_SUCCESS
  IMPLICIT NONE

  PRIVATE

!-----------------------------------------------------------------------------
! !FLAG TYPES AND VALUES
!-----------------------------------------------------------------------------

!> @cond IGNORE_TYPE_FLAGS
  type missingval_flag
    sequence
    private
      integer :: opt
  end type missingval_flag

  type field_init_flag
    sequence
    private
      integer :: opt
  end type field_init_flag
!> @endcond

  type(missingval_flag), parameter ::        &
    MISSINGVAL_ERROR  = missingval_flag(-1), &
    MISSINGVAL_IGNORE = missingval_flag(0),  &
    MISSINGVAL_FAIL   = missingval_flag(1),  &
    MISSINGVAL_SKPCPY = missingval_flag(2)

  type(field_init_flag), parameter ::       &
    FLD_INIT_ERROR   = field_init_flag(-1), &
    FLD_INIT_ZERO    = field_init_flag(0),  &
    FLD_INIT_MODEL   = field_init_flag(1),  &
    FLD_INIT_FILLV   = field_init_flag(2),  &
    FLD_INIT_IMPORT  = field_init_flag(3)

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------

  public missingval_flag
  public field_init_flag
  public MISSINGVAL_ERROR
  public MISSINGVAL_IGNORE
  public MISSINGVAL_FAIL
  public MISSINGVAL_SKPCPY
  public FLD_INIT_ERROR
  public FLD_INIT_ZERO
  public FLD_INIT_MODEL
  public FLD_INIT_FILLV
  public FLD_INIT_IMPORT

  public operator(==), assignment(=)

!-----------------------------------------------------------------------------
! !INTERFACE DEFINITIONS:
!-----------------------------------------------------------------------------

!> @cond IGNORE_INTERFACES
  interface operator (==)
    module procedure missingval_eq
    module procedure field_init_eq
  end interface

  interface assignment (=)
    module procedure missingval_toString
    module procedure missingval_frString
    module procedure field_init_toString
    module procedure field_init_frString
  end interface
!> @endcond

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "missingval_eq"
  function missingval_eq(val1, val2)
    logical missingval_eq
    type(missingval_flag), intent(in) :: val1, val2
    missingval_eq = (val1%opt == val2%opt)
  end function missingval_eq

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "missingval_toString"
  subroutine missingval_toString(string, val)
    character(len=*), intent(out) :: string
    type(missingval_flag), intent(in) :: val
    if (val == MISSINGVAL_IGNORE) then
      write(string,'(a)') 'MISSINGVAL_IGNORE'
    elseif (val == MISSINGVAL_FAIL) then
      write(string,'(a)') 'MISSINGVAL_FAIL'
    elseif (val == MISSINGVAL_SKPCPY) then
      write(string,'(a)') 'MISSINGVAL_SKPCPY'
    else
      write(string,'(a)') 'MISSINGVAL_ERROR'
    endif
  end subroutine missingval_toString

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "missingval_frString"
  subroutine missingval_frString(val, string)
    type(missingval_flag), intent(out) :: val
    character(len=*), intent(in) :: string
    character(len=32) :: ustring
    integer :: rc
    ustring = ESMF_UtilStringUpperCase(string, rc=rc)
    if (rc .ne. ESMF_SUCCESS) then
      val = MISSINGVAL_ERROR
    elseif (ustring .eq. 'MISSINGVAL_IGNORE') then
      val = MISSINGVAL_IGNORE
    elseif (ustring .eq. 'MISSINGVAL_FAIL') then
      val = MISSINGVAL_FAIL
    elseif (ustring .eq. 'MISSINGVAL_SKPCPY') then
      val = MISSINGVAL_SKPCPY
    else
      val = MISSINGVAL_ERROR
    endif
  end subroutine missingval_frString

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "field_init_eq"
  function field_init_eq(val1, val2)
    logical field_init_eq
    type(field_init_flag), intent(in) :: val1, val2
    field_init_eq = (val1%opt == val2%opt)
  end function field_init_eq

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "field_init_toString"
  subroutine field_init_toString(string, val)
    character(len=*), intent(out) :: string
    type(field_init_flag), intent(in) :: val
    if (val == FLD_INIT_ZERO) then
      write(string,'(a)') 'FLD_INIT_ZERO'
    elseif (val == FLD_INIT_MODEL) then
      write(string,'(a)') 'FLD_INIT_MODEL'
    elseif (val == FLD_INIT_FILLV) then
      write(string,'(a)') 'FLD_INIT_FILLV'
    elseif (val == FLD_INIT_IMPORT) then
      write(string,'(a)') 'FLD_INIT_IMPORT'
    else
      write(string,'(a)') 'FLD_INIT_ERROR'
    endif
  end subroutine field_init_toString

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "field_init_frString"
  subroutine field_init_frString(val, string)
    type(field_init_flag), intent(out) :: val
    character(len=*), intent(in) :: string
    character(len=16) :: ustring
    integer :: rc
    ustring = ESMF_UtilStringUpperCase(string, rc=rc)
    if (rc .ne. ESMF_SUCCESS) then
      val = FLD_INIT_ERROR
    elseif (ustring .eq. 'FLD_INIT_ZERO') then
      val = FLD_INIT_ZERO
    elseif (ustring .eq. 'FLD_INIT_MODEL') then
      val = FLD_INIT_MODEL
    elseif (ustring .eq. 'FLD_INIT_FILLV') then
      val = FLD_INIT_FILLV
    elseif (ustring .eq. 'FLD_INIT_MISSING') then
      val = FLD_INIT_FILLV
    elseif (ustring .eq. 'FLD_INIT_IMPORT') then
      val = FLD_INIT_IMPORT
    else
      val = FLD_INIT_ERROR
    endif
  end subroutine field_init_frString

  !-----------------------------------------------------------------------------

end module
