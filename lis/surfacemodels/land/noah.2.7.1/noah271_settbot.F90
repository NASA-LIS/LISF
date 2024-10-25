!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: noah271_settbot
! \label{noah271_settbot}
!
! !REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added Noah2.5 LSM, Initial Code
!  13 Oct 2003: Sujay Kumar; Domain independent modifications
!
! !INTERFACE:
subroutine noah271_settbot(mtype)
! !USES:
  use LIS_coreMod
  use LIS_tbotAdjustMod, only : LIS_initTmnUpdateTile
  use LIS_fileIOMod

  use noah271_lsmMod      
!
! !DESCRIPTION:
!  This subroutine retrieves bottom boundary temperature for Noah2.7.1 LSM.
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
  integer :: mtype
  real, allocatable :: placetbot1(:,:)

  do n=1,LIS_rc%nnest
!-----------------------------------------------------------------------
! Read in bottom temperature fields and adjust for elevation difference
! with either Eta (NLDAS) or GDAS (GLDAS) correction datasets. 
!-----------------------------------------------------------------------
     allocate(placetbot1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     call LIS_read_param(n,"TBOT",placetbot1)

     do i = 1, LIS_rc%npatch(n,mtype)
        if(placetbot1(LIS_surface(n,mtype)%tile(i)%col,LIS_surface(n,mtype)%tile(i)%row)&
             .ne.-9999.00) then
           noah271_struc(n)%noah(i)%tempbot = &
                placetbot1(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
        endif

        !Initialize variables for dynamic deep soil temperature
        call LIS_initTmnUpdateTile(n,i,noah271_struc(n)%noah(i)%tempbot)

     enddo
     deallocate(placetbot1)
  enddo

end subroutine noah271_settbot
