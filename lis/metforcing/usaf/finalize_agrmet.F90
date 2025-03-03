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
! !ROUTINE: finalize_agrmet
! \label{finalize_agrmet}
!
! !REVISION HISTORY:
! 25Oct2005; Sujay Kumar, Initial Code
! 31 MAR 2010 Add release of additional resolutions' gridded data structures
!             .....see in code comments for where...........Michael Shaw/WXE
! 13 APR 2012 Test for offline mode and deallocate only valid allocatables
!             ....................................Chris Franks/16WS/WXE/SEMS
! 9 May 2013  Add deallocations for cmorph
!             ......................................Ryan Ruhge/16WS/WXE/SEMS
! !INTERFACE:    
subroutine finalize_agrmet
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_pluginIndices, only : LIS_retroId
  use AGRMET_forcingMod, only : agrmet_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for AGRMET forcing. 
!
!EOP 
  implicit none
  integer :: n 

#if 0
  do n=1,LIS_rc%nnest
     deallocate(agrmet_struc(n)%smask1)
     deallocate(agrmet_struc(n)%smask2)
     if ( LIS_rc%runmode .ne. LIS_retroId ) then

       if(LIS_rc%gridDesc(n,4).ge.0.or.LIS_rc%gridDesc(n,1).eq.0.and.LIS_rc%gridDesc(n,4).ge.0.and.LIS_rc%gridDesc(n,7).ge.0) then
          deallocate(agrmet_struc(n)%rlat1_nh)
          deallocate(agrmet_struc(n)%rlon1_nh)
          deallocate(agrmet_struc(n)%n111_nh)
          deallocate(agrmet_struc(n)%n121_nh)
          deallocate(agrmet_struc(n)%n211_nh)
          deallocate(agrmet_struc(n)%n221_nh)
          deallocate(agrmet_struc(n)%w111_nh)
          deallocate(agrmet_struc(n)%w121_nh)
          deallocate(agrmet_struc(n)%w211_nh)
          deallocate(agrmet_struc(n)%w221_nh)
          deallocate(agrmet_struc(n)%rlat2_nh)
          deallocate(agrmet_struc(n)%rlon2_nh)
          deallocate(agrmet_struc(n)%n112_nh)
! Michael Shaw - If there's a different geoprecip
! from native, structures were created, so release them
          if(agrmet_struc(n)%imaxgp /= agrmet_struc(n)%imaxnative)then
             deallocate(agrmet_struc(n)%rlat1_nh2)
             deallocate(agrmet_struc(n)%rlon1_nh2)
             deallocate(agrmet_struc(n)%n111_nh2)
             deallocate(agrmet_struc(n)%n121_nh2)
             deallocate(agrmet_struc(n)%n211_nh2)
             deallocate(agrmet_struc(n)%n221_nh2)
             deallocate(agrmet_struc(n)%w111_nh2)
             deallocate(agrmet_struc(n)%w121_nh2)
             deallocate(agrmet_struc(n)%w211_nh2)
             deallocate(agrmet_struc(n)%w221_nh2)
             deallocate(agrmet_struc(n)%rlat2_nh2)
             deallocate(agrmet_struc(n)%rlon2_nh2)
             deallocate(agrmet_struc(n)%n112_nh2)
          endif
! Michael Shaw - If there's a different ssmi(/s) 
! from native, structures were created, so release them
          if(agrmet_struc(n)%imaxsmi /= agrmet_struc(n)%imaxnative)then
             deallocate(agrmet_struc(n)%rlat1_nh3)
             deallocate(agrmet_struc(n)%rlon1_nh3)
             deallocate(agrmet_struc(n)%n111_nh3)
             deallocate(agrmet_struc(n)%n121_nh3)
             deallocate(agrmet_struc(n)%n211_nh3)
             deallocate(agrmet_struc(n)%n221_nh3)
             deallocate(agrmet_struc(n)%w111_nh3)
             deallocate(agrmet_struc(n)%w121_nh3)
             deallocate(agrmet_struc(n)%w211_nh3)
             deallocate(agrmet_struc(n)%w221_nh3)
             deallocate(agrmet_struc(n)%rlat2_nh3)
             deallocate(agrmet_struc(n)%rlon2_nh3)
             deallocate(agrmet_struc(n)%n112_nh3)
          endif
       elseif(LIS_rc%gridDesc(n,4).le.0.or.LIS_rc%gridDesc(n,1).eq.0.and.LIS_rc%gridDesc(n,4).ge.0.and.LIS_rc%gridDesc(n,7).le.0) then
          deallocate(agrmet_struc(n)%rlat1_sh)
          deallocate(agrmet_struc(n)%rlon1_sh)
          deallocate(agrmet_struc(n)%n111_sh)
          deallocate(agrmet_struc(n)%n121_sh)
          deallocate(agrmet_struc(n)%n211_sh)
          deallocate(agrmet_struc(n)%n221_sh)
          deallocate(agrmet_struc(n)%w111_sh)
          deallocate(agrmet_struc(n)%w121_sh)
          deallocate(agrmet_struc(n)%w211_sh)
          deallocate(agrmet_struc(n)%w221_sh)
          deallocate(agrmet_struc(n)%rlat2_sh)
          deallocate(agrmet_struc(n)%rlon2_sh)
          deallocate(agrmet_struc(n)%n112_sh)
! Michael Shaw - If there's a different geoprecip
! from native, structures were created, so release them
          if(agrmet_struc(n)%imaxgp /= agrmet_struc(n)%imaxnative)then
             deallocate(agrmet_struc(n)%rlat1_sh2)
             deallocate(agrmet_struc(n)%rlon1_sh2)
             deallocate(agrmet_struc(n)%n111_sh2)
             deallocate(agrmet_struc(n)%n121_sh2)
             deallocate(agrmet_struc(n)%n211_sh2)
             deallocate(agrmet_struc(n)%n221_sh2)
             deallocate(agrmet_struc(n)%w111_sh2)
             deallocate(agrmet_struc(n)%w121_sh2)
             deallocate(agrmet_struc(n)%w211_sh2)
             deallocate(agrmet_struc(n)%w221_sh2)
             deallocate(agrmet_struc(n)%rlat2_sh2)
             deallocate(agrmet_struc(n)%rlon2_sh2)
             deallocate(agrmet_struc(n)%n112_sh2)
          endif
! Michael Shaw - If there's a different ssmi(/s) 
! from native, structures were created, so release them
          if(agrmet_struc(n)%imaxsmi /= agrmet_struc(n)%imaxnative)then
             deallocate(agrmet_struc(n)%rlat1_sh3)
             deallocate(agrmet_struc(n)%rlon1_sh3)
             deallocate(agrmet_struc(n)%n111_sh3)
             deallocate(agrmet_struc(n)%n121_sh3)
             deallocate(agrmet_struc(n)%n211_sh3)
             deallocate(agrmet_struc(n)%n221_sh3)
             deallocate(agrmet_struc(n)%w111_sh3)
             deallocate(agrmet_struc(n)%w121_sh3)
             deallocate(agrmet_struc(n)%w211_sh3)
             deallocate(agrmet_struc(n)%w221_sh3)
             deallocate(agrmet_struc(n)%rlat2_sh3)
             deallocate(agrmet_struc(n)%rlon2_sh3)
             deallocate(agrmet_struc(n)%n112_sh3)
          endif
       else
          deallocate(agrmet_struc(n)%rlat1_nh)
          deallocate(agrmet_struc(n)%rlon1_nh)
          deallocate(agrmet_struc(n)%n111_nh)
          deallocate(agrmet_struc(n)%n121_nh)
          deallocate(agrmet_struc(n)%n211_nh)
          deallocate(agrmet_struc(n)%n221_nh)
          deallocate(agrmet_struc(n)%w111_nh)
          deallocate(agrmet_struc(n)%w121_nh)
          deallocate(agrmet_struc(n)%w211_nh)
          deallocate(agrmet_struc(n)%w221_nh)
          deallocate(agrmet_struc(n)%rlat2_nh)
          deallocate(agrmet_struc(n)%rlon2_nh)
          deallocate(agrmet_struc(n)%n112_nh)
          deallocate(agrmet_struc(n)%rlat1_sh)
          deallocate(agrmet_struc(n)%rlon1_sh)
          deallocate(agrmet_struc(n)%n111_sh)
          deallocate(agrmet_struc(n)%n121_sh)
          deallocate(agrmet_struc(n)%n211_sh)
          deallocate(agrmet_struc(n)%n221_sh)
          deallocate(agrmet_struc(n)%w111_sh)
          deallocate(agrmet_struc(n)%w121_sh)
          deallocate(agrmet_struc(n)%w211_sh)
          deallocate(agrmet_struc(n)%w221_sh)
          deallocate(agrmet_struc(n)%rlat2_sh)
          deallocate(agrmet_struc(n)%rlon2_sh)
          deallocate(agrmet_struc(n)%n112_sh)
! Michael Shaw - If there's a different geoprecip
! from native, structures were created, so release them
          if(agrmet_struc(n)%imaxgp /= agrmet_struc(n)%imaxnative)then
             deallocate(agrmet_struc(n)%rlat1_nh2)
             deallocate(agrmet_struc(n)%rlon1_nh2)
             deallocate(agrmet_struc(n)%n111_nh2)
             deallocate(agrmet_struc(n)%n121_nh2)
             deallocate(agrmet_struc(n)%n211_nh2)
             deallocate(agrmet_struc(n)%n221_nh2)
             deallocate(agrmet_struc(n)%w111_nh2)
             deallocate(agrmet_struc(n)%w121_nh2)
             deallocate(agrmet_struc(n)%w211_nh2)
             deallocate(agrmet_struc(n)%w221_nh2)
             deallocate(agrmet_struc(n)%rlat2_sh2)
             deallocate(agrmet_struc(n)%rlon2_sh2)
             deallocate(agrmet_struc(n)%n112_sh2)
             deallocate(agrmet_struc(n)%rlat1_sh2)
             deallocate(agrmet_struc(n)%rlon1_sh2)
             deallocate(agrmet_struc(n)%n111_sh2)
             deallocate(agrmet_struc(n)%n121_sh2)
             deallocate(agrmet_struc(n)%n211_sh2)
             deallocate(agrmet_struc(n)%n221_sh2)
             deallocate(agrmet_struc(n)%w111_sh2)
             deallocate(agrmet_struc(n)%w121_sh2)
             deallocate(agrmet_struc(n)%w211_sh2)
             deallocate(agrmet_struc(n)%w221_sh2)
          endif
! Michael Shaw - If there's a different ssmi(/s) 
! from native, structures were created, so release them
          if(agrmet_struc(n)%imaxsmi /= agrmet_struc(n)%imaxnative)then
             deallocate(agrmet_struc(n)%rlat1_nh3)
             deallocate(agrmet_struc(n)%rlon1_nh3)
             deallocate(agrmet_struc(n)%n111_nh3)
             deallocate(agrmet_struc(n)%n121_nh3)
             deallocate(agrmet_struc(n)%n211_nh3)
             deallocate(agrmet_struc(n)%n221_nh3)
             deallocate(agrmet_struc(n)%w111_nh3)
             deallocate(agrmet_struc(n)%w121_nh3)
             deallocate(agrmet_struc(n)%w211_nh3)
             deallocate(agrmet_struc(n)%w221_nh3)
             deallocate(agrmet_struc(n)%rlat2_nh3)
             deallocate(agrmet_struc(n)%rlon2_nh3)
             deallocate(agrmet_struc(n)%n112_nh3)
             deallocate(agrmet_struc(n)%rlat1_sh3)
             deallocate(agrmet_struc(n)%rlon1_sh3)
             deallocate(agrmet_struc(n)%n111_sh3)
             deallocate(agrmet_struc(n)%n121_sh3)
             deallocate(agrmet_struc(n)%n211_sh3)
             deallocate(agrmet_struc(n)%n221_sh3)
             deallocate(agrmet_struc(n)%w111_sh3)
             deallocate(agrmet_struc(n)%w121_sh3)
             deallocate(agrmet_struc(n)%w211_sh3)
             deallocate(agrmet_struc(n)%w221_sh3)
             deallocate(agrmet_struc(n)%rlat2_sh3)
             deallocate(agrmet_struc(n)%rlon2_sh3)
             deallocate(agrmet_struc(n)%n112_sh3)
          endif
       endif
       deallocate(agrmet_struc(n)%clippd)
       deallocate(agrmet_struc(n)%cliprc)
       deallocate(agrmet_struc(n)%clirtn)
       
       deallocate(agrmet_struc(n)%clippd1)
       deallocate(agrmet_struc(n)%cliprc1)
       deallocate(agrmet_struc(n)%clirtn1)

       deallocate(agrmet_struc(n)%clippd2)
       deallocate(agrmet_struc(n)%cliprc2)
       deallocate(agrmet_struc(n)%clirtn2)
  
       deallocate(agrmet_struc(n)%mrgp)
     
       deallocate(agrmet_struc(n)%rlatcmor)
       deallocate(agrmet_struc(n)%rloncmor)
       deallocate(agrmet_struc(n)%n112cmor)
       deallocate(agrmet_struc(n)%n122cmor)
       deallocate(agrmet_struc(n)%n212cmor)
       deallocate(agrmet_struc(n)%n222cmor)
       deallocate(agrmet_struc(n)%w112cmor)
       deallocate(agrmet_struc(n)%w122cmor)
       deallocate(agrmet_struc(n)%w212cmor)
       deallocate(agrmet_struc(n)%w222cmor)
     
     else 
       deallocate(agrmet_struc(n)%rlat_glb)
       deallocate(agrmet_struc(n)%rlon_glb)
       deallocate(agrmet_struc(n)%n11_glb)
       deallocate(agrmet_struc(n)%n12_glb)
       deallocate(agrmet_struc(n)%n21_glb)
       deallocate(agrmet_struc(n)%n22_glb)
       deallocate(agrmet_struc(n)%w11_glb)
       deallocate(agrmet_struc(n)%w12_glb)
       deallocate(agrmet_struc(n)%w21_glb)
       deallocate(agrmet_struc(n)%w22_glb)
     end if
  enddo
  deallocate(agrmet_struc)
#endif
end subroutine finalize_agrmet
