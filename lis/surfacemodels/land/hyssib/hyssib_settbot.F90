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
! !ROUTINE: hyssib_settbot
! \label{hyssib_settbot}
!
! !REVISION HISTORY:
! 03 Jun 2005: David Mocko, Conversion from NOAH to HY-SSiB
! 05 Sep 2007: Chuck Alonge, Updates for LIS 5.0 compliance
!
! !INTERFACE:
subroutine hyssib_settbot
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_fileIOMod, only : LIS_read_param
  use hyssib_lsmMod      

!
! !DESCRIPTION:
!  This subroutine retrieves HY-SSiB bottom soil temperature
!  Currently a static data is being used. If a climatology is available, 
!  it can be used by setting the appropriate frequency of the climatology
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
!   retrieves the bottom temperature data
!  \end{description}
!EOP
  implicit none

  integer :: i,n
  real, allocatable :: placetbot1(:,:)
 
  do n=1,LIS_rc%nnest
!-----------------------------------------------------------------------
! Read in bottom temperature fields and adjust for elevation differnce
! with either Eta (NLDAS) or GDAS (GLDAS) correction datasets. 
!-----------------------------------------------------------------------
     allocate(placetbot1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     call LIS_read_param(n, "TBOT",placetbot1)
     
     do i = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
        if(placetbot1(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
             LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)&
             .ne.-9999.00) then
           hyssib_struc(n)%hyssib(i)%tempbot = &
                placetbot1(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
        endif
     enddo
     deallocate(placetbot1)
  enddo
end subroutine hyssib_settbot

