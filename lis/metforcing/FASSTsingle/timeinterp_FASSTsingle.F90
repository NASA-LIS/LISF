!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: timeinterp_FASSTsingle
! \label{timeinterp_FASSTsingle}
! 
!
! !REVISION HISTORY:
! 13 Apr 2007: Bailing Li, Initial Specification
! 15 Oct 2010: David Mocko, Updated for FASST single-point test case
!
! !INTERFACE:
subroutine timeinterp_FASSTsingle(n,findex)
  ! !USES:
  use ESMF
  use LIS_logMod, only           : LIS_logunit,LIS_verify,LIS_endrun
  use LIS_coreMod, only          : LIS_rc,LIS_domain
  use LIS_metforcingMod, only   : LIS_FORC_Base_State
  use FASSTsingle_forcingMod, only : FASSTsingle_struc
  use LIS_FORC_AttributesMod

  implicit none
  ! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex

  ! !DESCRIPTION:
  ! Temporally interpolates the FASSTsingle forcing data to
  ! the model timestep.  All variables except precipitation
  ! is linearly interpolated. 
  ! 
  !EOP
  integer :: t,f
  integer :: index1
  integer :: tid
  !      write(LIS_logunit,*) 'starting timeinterp_FASSTsingle'

  do t = 1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     do f = 1,35
        FASSTsingle_struc(n)%metm(f) =                             &
             FASSTsingle_struc(n)%metm_back(f)
     enddo
  enddo

end subroutine timeinterp_FASSTsingle

