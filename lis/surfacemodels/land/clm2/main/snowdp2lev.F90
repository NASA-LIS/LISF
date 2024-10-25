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
! Purpose: 
! create snow layers and interfaces given snow depth
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: snowdp2lev.F90,v 1.5 2004/05/07 22:18:36 jim Exp $
!-----------------------------------------------------------------------
subroutine snowdp2lev(n)

  use LIS_coreMod, only : LIS_rc
  use LIS_precisionMod
  use clm2_varpar, only : nlevsoi, nlevsno, nlevlak
  use clm2_lsmMod, only : clm2_struc
  implicit none

! ------------------- local variables -----------------------------
  integer, intent(in) :: n 
  integer i,k    !indices
! -----------------------------------------------------------------

  do k = 1,LIS_rc%ntiles(n)
     clm2_struc(n)%clm(k)%dz(-nlevsno+1:0) = 1.e36 
     clm2_struc(n)%clm(k)%z (-nlevsno+1:0) = 1.e36 
     clm2_struc(n)%clm(k)%zi(-nlevsno:-1)  = 1.e36 
     if (.not. clm2_struc(n)%clm(k)%lakpoi) then  !not lake
        if (clm2_struc(n)%clm(k)%snowdp < 0.01) then
           clm2_struc(n)%clm(k)%snl = 0
           clm2_struc(n)%clm(k)%dz(-nlevsno+1:0) = 0.
           clm2_struc(n)%clm(k)%z (-nlevsno+1:0) = 0.
           clm2_struc(n)%clm(k)%zi(-nlevsno+0:0) = 0.
        else
           if ((clm2_struc(n)%clm(k)%snowdp >= 0.01) .AND. (clm2_struc(n)%clm(k)%snowdp <= 0.03)) then
              clm2_struc(n)%clm(k)%snl = -1
              clm2_struc(n)%clm(k)%dz(0)  = clm2_struc(n)%clm(k)%snowdp
           else if ((clm2_struc(n)%clm(k)%snowdp > 0.03) .AND. (clm2_struc(n)%clm(k)%snowdp <= 0.04)) then
              clm2_struc(n)%clm(k)%snl = -2
              clm2_struc(n)%clm(k)%dz(-1) = clm2_struc(n)%clm(k)%snowdp/2.
              clm2_struc(n)%clm(k)%dz( 0) = clm2_struc(n)%clm(k)%dz(-1)
           else if ((clm2_struc(n)%clm(k)%snowdp > 0.04) .AND. (clm2_struc(n)%clm(k)%snowdp <= 0.07)) then
              clm2_struc(n)%clm(k)%snl = -2
              clm2_struc(n)%clm(k)%dz(-1) = 0.02
              clm2_struc(n)%clm(k)%dz( 0) = clm2_struc(n)%clm(k)%snowdp - clm2_struc(n)%clm(k)%dz(-1)
           else if ((clm2_struc(n)%clm(k)%snowdp > 0.07) .AND. (clm2_struc(n)%clm(k)%snowdp <= 0.12)) then
              clm2_struc(n)%clm(k)%snl = -3
              clm2_struc(n)%clm(k)%dz(-2) = 0.02
              clm2_struc(n)%clm(k)%dz(-1) = (clm2_struc(n)%clm(k)%snowdp - 0.02)/2.
              clm2_struc(n)%clm(k)%dz( 0) = clm2_struc(n)%clm(k)%dz(-1)
           else if ((clm2_struc(n)%clm(k)%snowdp > 0.12) .AND. (clm2_struc(n)%clm(k)%snowdp <= 0.18)) then
              clm2_struc(n)%clm(k)%snl = -3
              clm2_struc(n)%clm(k)%dz(-2) = 0.02
              clm2_struc(n)%clm(k)%dz(-1) = 0.05
              clm2_struc(n)%clm(k)%dz( 0) = clm2_struc(n)%clm(k)%snowdp - clm2_struc(n)%clm(k)%dz(-2) - clm2_struc(n)%clm(k)%dz(-1)
           else if ((clm2_struc(n)%clm(k)%snowdp > 0.18) .AND. (clm2_struc(n)%clm(k)%snowdp <= 0.29)) then
              clm2_struc(n)%clm(k)%snl = -4
              clm2_struc(n)%clm(k)%dz(-3) = 0.02
              clm2_struc(n)%clm(k)%dz(-2) = 0.05
              clm2_struc(n)%clm(k)%dz(-1) = (clm2_struc(n)%clm(k)%snowdp - &
                                            clm2_struc(n)%clm(k)%dz(-3) - &
                                            clm2_struc(n)%clm(k)%dz(-2))/2.
              clm2_struc(n)%clm(k)%dz( 0) = clm2_struc(n)%clm(k)%dz(-1)
           else if ((clm2_struc(n)%clm(k)%snowdp > 0.29) .and. &
                    (clm2_struc(n)%clm(k)%snowdp <= 0.41)) then
              clm2_struc(n)%clm(k)%snl = -4
              clm2_struc(n)%clm(k)%dz(-3) = 0.02
              clm2_struc(n)%clm(k)%dz(-2) = 0.05
              clm2_struc(n)%clm(k)%dz(-1) = 0.11
              clm2_struc(n)%clm(k)%dz( 0) = clm2_struc(n)%clm(k)%snowdp - &
                                           clm2_struc(n)%clm(k)%dz(-3) - &
                                           clm2_struc(n)%clm(k)%dz(-2) - &
                                           clm2_struc(n)%clm(k)%dz(-1)
           else if ((clm2_struc(n)%clm(k)%snowdp > 0.41) .and. &
                    (clm2_struc(n)%clm(k)%snowdp <= 0.64)) then
              clm2_struc(n)%clm(k)%snl = -5
              clm2_struc(n)%clm(k)%dz(-4) = 0.02
              clm2_struc(n)%clm(k)%dz(-3) = 0.05
              clm2_struc(n)%clm(k)%dz(-2) = 0.11
              clm2_struc(n)%clm(k)%dz(-1) = (clm2_struc(n)%clm(k)%snowdp - &
                                            clm2_struc(n)%clm(k)%dz(-4) - &
                                            clm2_struc(n)%clm(k)%dz(-3) - &
                                            clm2_struc(n)%clm(k)%dz(-2))/2.
              clm2_struc(n)%clm(k)%dz( 0) = clm2_struc(n)%clm(k)%dz(-1)
           else if (clm2_struc(n)%clm(k)%snowdp > 0.64) then 
              clm2_struc(n)%clm(k)%snl = -5
              clm2_struc(n)%clm(k)%dz(-4) = 0.02
              clm2_struc(n)%clm(k)%dz(-3) = 0.05
              clm2_struc(n)%clm(k)%dz(-2) = 0.11
              clm2_struc(n)%clm(k)%dz(-1) = 0.23
              clm2_struc(n)%clm(k)%dz( 0) = clm2_struc(n)%clm(k)%snowdp - &
                                           clm2_struc(n)%clm(k)%dz(-4) - &
                                           clm2_struc(n)%clm(k)%dz(-3) - &
                                           clm2_struc(n)%clm(k)%dz(-2) - &
                                           clm2_struc(n)%clm(k)%dz(-1)
           endif
           do i = 0, clm2_struc(n)%clm(k)%snl+1, -1
              clm2_struc(n)%clm(k)%z(i)    = clm2_struc(n)%clm(k)%zi(i) - 0.5*clm2_struc(n)%clm(k)%dz(i)
              clm2_struc(n)%clm(k)%zi(i-1) = clm2_struc(n)%clm(k)%zi(i) - clm2_struc(n)%clm(k)%dz(i)
           enddo
        endif 
     else   !lake points
        clm2_struc(n)%clm(k)%snl = 0
        clm2_struc(n)%clm(k)%dz(-nlevsno+1:0) = 0.
        clm2_struc(n)%clm(k)%z (-nlevsno+1:0) = 0.
        clm2_struc(n)%clm(k)%zi(-nlevsno+0:0) = 0.
     endif
  end do

end subroutine snowdp2lev
