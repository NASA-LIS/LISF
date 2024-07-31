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
! !ROUTINE: noah36_settbot
! \label{noah36_settbot}
!
! !REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added Noah2.5 LSM, Initial Code
!  13 Oct 2003: Sujay Kumar; Domain independent modifications
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  23 May 2014: David Mocko, changes for Noah3.3 in LIS7.0
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
!
! !INTERFACE:
subroutine noah36_settbot(mtype)
! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_surface
  use LIS_logMod,     only : LIS_logunit
  use LIS_tbotAdjustMod, only : LIS_initTmnUpdateTile
  use LIS_fileIOMod,  only : LIS_read_param
  use noah36_lsmMod      
!
! !DESCRIPTION:
!  This subroutine retrieves bottom boundary temperature for Noah-3.6 LSM.
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
  integer :: mtype

  integer :: i,n
  real, allocatable :: placetbot1(:,:)

  do n=1,LIS_rc%nnest
!-----------------------------------------------------------------------
! Read in bottom temperature fields and adjust for elevation difference
! with either Eta (NLDAS) or GDAS (GLDAS) correction datasets. 
!-----------------------------------------------------------------------
     allocate(placetbot1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     if (noah36_struc(n)%fixedtbot.ne.0.0) then
        write(LIS_logunit,*) '[INFO] Fixing Noah-3.6 deep soil temperature: ',&
                              noah36_struc(n)%fixedtbot
     else
        call LIS_read_param(n,"TBOT",placetbot1)
     endif

     do i = 1, LIS_rc%npatch(n,mtype)
        if (noah36_struc(n)%fixedtbot.ne.0.0) then
           noah36_struc(n)%noah(i)%tempbot = noah36_struc(n)%fixedtbot
        else
           if (placetbot1(LIS_surface(n,mtype)%tile(i)%col,                   &
                          LIS_surface(n,mtype)%tile(i)%row).ne.-9999.00) then
              noah36_struc(n)%noah(i)%tempbot =                        &
                            placetbot1(LIS_surface(n,mtype)%tile(i)%col,      &
                                       LIS_surface(n,mtype)%tile(i)%row)
           else
              noah36_struc(n)%noah(i)%tempbot = 285.0
              write(LIS_logunit,*) '[INFO] Noah-3.6 bottom temperature -- ',  &
                                   'not defined for point ',i
              write(LIS_logunit,*) '[INFO] Noah-3.6 bottom temperature -- ',  &
                                   'set to default value of 285.0'
           endif
        endif

        !Initialize variables for dynamic deep soil temperature
        call LIS_initTmnUpdateTile(n,i,noah36_struc(n)%noah(i)%tempbot)

     enddo
     deallocate(placetbot1)
  enddo

end subroutine noah36_settbot
