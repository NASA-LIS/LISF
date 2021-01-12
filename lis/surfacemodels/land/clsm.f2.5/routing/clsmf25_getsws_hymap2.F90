!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: clsmf25_getsws_hymap2
!  \label{clsmf25_getsws_hymap2}
!
! !REVISION HISTORY:
! 12 Sep 2019: Augusto Getirana; implementation of two-way coupling
!
! !INTERFACE:
subroutine clsmf25_getsws_hymap2(n)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_routingMod, only : LIS_runoff_state
  use LIS_logMod
  use LIS_historyMod
  use clsmf25_lsmMod, only : clsmf25_struc

  implicit none
! !ARGUMENTS: 
  integer,  intent(in)   :: n 
!
! !DESCRIPTION:
!   This routine defines the surface water storage variables in the LSM
!   to be updated based on feedback from HYMAP2
!  
!EOP
  integer                :: t
  integer                :: c,r
  integer                :: status
  integer                :: enable2waycpl
  
  call ESMF_AttributeGet(LIS_runoff_state(n),"2 way coupling",&
       enable2waycpl, rc=status)
  call LIS_verify(status)
  
  if(enable2waycpl==1) then
     write(LIS_logunit,*) '[ERR] Two-way coupling between CLSM and HYMAP2'
     write(LIS_logunit,*) '[ERR] is not currently supported'
     call LIS_endrun()
  endif

end subroutine clsmf25_getsws_hymap2
