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

subroutine CombineSnowLayers (clm) 

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Combine snow layers that are less than a minimum thickness or mass
!
! Method:
! If the snow element thickness or mass is less than a prescribed minimum,
! then it is combined with a neighboring element.  The subroutine 
! clm_combo.f90 then executes the combination of mass and energy.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: CombineSnowLayers.F90,v 1.6 2004/11/24 22:56:21 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, only : istsoil
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm      !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer i         ! do loop index 
  integer j         ! node index
  integer k         ! do loop index
  integer l         ! node index
  integer msn_old   ! number of snow layer 1 (top) to msn0 (bottom)
  integer mssi      ! node index
  integer neibor    ! adjacent node selected for combination
  real(r8) dzmin(5) ! minimum of snow layer 1 (top) to msn0 (bottom)
  real(r8) zwice    ! total ice mass in snow
  real(r8) zwliq    ! total liquid water in snow

  data dzmin /0.010, 0.015, 0.025, 0.055, 0.115/

!----End Variable List--------------------------------------------------

!
! Check the mass of ice lens of snow, when the total is less than a small value,
! combine it with the underlying neighbor.
!

  msn_old = clm%snl
  do j = msn_old+1, 0
     if(clm%h2osoi_ice(j) <= .1)then

        if (clm%itypwat == istsoil) then
         clm%h2osoi_liq(j+1) = clm%h2osoi_liq(j+1) + clm%h2osoi_liq(j)
         clm%h2osoi_ice(j+1) = clm%h2osoi_ice(j+1) + clm%h2osoi_ice(j)
        else if (clm%itypwat /= istsoil .and. j /= 0) then
         clm%h2osoi_liq(j+1) = clm%h2osoi_liq(j+1) + clm%h2osoi_liq(j)
         clm%h2osoi_ice(j+1) = clm%h2osoi_ice(j+1) + clm%h2osoi_ice(j)
        endif

!
! shift all elements above this down one.
!

        if(j > clm%snl+1 .AND. clm%snl < -1)then
           do i =  j, clm%snl+2, -1
              clm%t_soisno(i) = clm%t_soisno(i-1)
              clm%h2osoi_liq(i) = clm%h2osoi_liq(i-1)
              clm%h2osoi_ice(i) = clm%h2osoi_ice(i-1)
              clm%dz(i) = clm%dz(i-1)
           enddo
        endif
        clm%snl = clm%snl + 1
     endif
  enddo

  if(clm%snl == 0)then
     clm%h2osno = 0.
     clm%snowdp = 0.
     return
  else
     clm%h2osno = 0.
     clm%snowdp = 0.
     zwice = 0.
     zwliq = 0.
     do j = clm%snl + 1, 0
        clm%h2osno = clm%h2osno + clm%h2osoi_ice(j) + clm%h2osoi_liq(j)
        clm%snowdp = clm%snowdp + clm%dz(j)
        zwice = zwice + clm%h2osoi_ice(j)
        zwliq = zwliq + clm%h2osoi_liq(j)
     enddo
  endif

!
! Check the snow depth
!

  if(clm%snowdp < 0.01)then       !!! all snow gone 
     clm%snl = 0
     clm%h2osno = zwice
     if(clm%h2osno <= 0.) clm%snowdp = 0.

!
! The liquid water assumed ponding on soil surface.
!

     if (clm%itypwat == istsoil) then
      clm%h2osoi_liq(1) = clm%h2osoi_liq(1) + zwliq
     endif

     return
  else                        !!! snow layers combined

!
! Two or more layers 
!

     if(clm%snl < -1)then
        msn_old = clm%snl
        mssi = 1
        do i = msn_old+1, 0

!
! If top node is removed, combine with bottom neighbor.
!

           if(clm%dz(i) < dzmin(mssi))then
              if(i == clm%snl+1)then
                 neibor = i + 1

!
! If the bottom neighbor is not snow, combine with the top neighbor.
!

              else if(i == 0)then
                 neibor = i - 1

!
! If none of the above special cases apply, combine with the thinnest neighbor
!

              else
                 neibor = i + 1
                 if((clm%dz(i-1)+clm%dz(i)) < (clm%dz(i+1)+clm%dz(i))) neibor = i-1
              endif

!
! Node l and j are combined and stored as node j.
!

              if(neibor > i)then
                 j = neibor
                 l = i
              else
                 j = i
                 l = neibor
              endif

              call Combo (clm%dz(j),         clm%h2osoi_liq(j), clm%h2osoi_ice(j), &
                          clm%t_soisno(j),   clm%dz(l),         clm%h2osoi_liq(l), &
                          clm%h2osoi_ice(l), clm%t_soisno(l) )

!
! Now shift all elements above this down one.
!

              if(j-1 > clm%snl+1) then
                 do k = j-1, clm%snl+2, -1
                    clm%t_soisno(k) = clm%t_soisno(k-1)
                    clm%h2osoi_ice(k) = clm%h2osoi_ice(k-1)
                    clm%h2osoi_liq(k) = clm%h2osoi_liq(k-1)
                    clm%dz(k) = clm%dz(k-1)
                 enddo
              endif

              clm%snl = clm%snl + 1

              if(clm%snl >= -1) EXIT

!
! The layer thickness is greater than the prescribed minimum value                      
!

           else
              mssi = mssi + 1 
           endif
        enddo

     endif

!
! Reset the node depth and the depth of layer interface
!

     do k = 0, clm%snl+1, -1
        clm%z(k) = clm%zi(k) - 0.5*clm%dz(k)
        clm%zi(k-1) = clm%zi(k) - clm%dz(k)
     enddo

  endif                       !!! snow layers combined 

end subroutine CombineSnowLayers 
