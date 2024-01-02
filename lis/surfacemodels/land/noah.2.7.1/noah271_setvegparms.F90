!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
!BOP
!
! !ROUTINE: noah271_setvegparms
! \label{noah271_setvegparms}
!
! !REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added Noah2.5 LSM, Initial Code
!  13 Oct 2003: Sujay Kumar; Domain independent modifications
!  25 Jul 2005: Sujay Kumar; Removed the dependency to SIB classes
!  04 Oct 2005: Matthew Garcia; additions for Noah Unified
!  04 Oct 2005: Matthew Garcia; additions for Distributed Noah Router
!
! !INTERFACE:
subroutine noah271_setvegparms(mtype)
! !USES:
  use LIS_coreMod
  use noah271_lsmMod      

! !DESCRIPTION:
!  This subroutine retrieves Noah2.7.1 vegetation parameters.
!  The current implementation uses a table-based lookup based
!  on vegetation classes to initialize the following parameters
!  \begin{verbatim}
!   nroot  - number of root zones
!   rsmin  - minimum canopy resistance
!   rgl    - solar radiation term in canopy resistance function
!   hs     - vapor pressure deficit term
!   snup   - threshold snow depth
!   z0     - roughness length
!   lai    - leaf area index
!  \end{verbatim}
!EOP

  implicit none
  integer :: mtype
  integer :: n,i,j,m
  real, allocatable:: value(:,:)

  real, allocatable :: temp(:,:)
  integer       :: c,r
  integer       :: nvegtypes

  do m=1,LIS_rc%nnest
! Set the section of the vegetation parameter table to read for Noah 2.7.1
     nvegtypes = 0
     if (LIS_rc%lcscheme.eq."UMD") nvegtypes = 13
     if (LIS_rc%lcscheme.eq."USGS") nvegtypes = 24
     if (LIS_rc%lcscheme.eq."MODIS") nvegtypes = 20
     if (LIS_rc%lcscheme.eq."IGBPNCEP") nvegtypes = 20

     !if (LIS_rc%uselcmap(m).eq."ECOCLIMAP2") nvegtypes = 13

     allocate(value(nvegtypes,noah271_struc(m)%nvegp))
!-----------------------------------------------------------------------
! Set Noah2.7.1 vegetation type at tile from the LIS domain
!-----------------------------------------------------------------------
     do n=1,LIS_rc%npatch(m,mtype)
        noah271_struc(m)%noah(n)%vegt = LIS_surface(m,mtype)%tile(n)%vegt
     enddo
     
!-----------------------------------------------------------------------
! Get Vegetation Parameters for Noah2.7.1 Model in Tile Space
! Read in the Noah2.7.1 Static Vegetation Parameter Files
!-----------------------------------------------------------------------
     open(unit=11,file=noah271_struc(m)%vfile,status='old')
     do j=1,noah271_struc(m)%nvegp
        read(11,*)
        read(11,*)(value(i,j),i=1,nvegtypes)
     enddo
     close(11)

     do i=1,LIS_rc%npatch(m,mtype)
        noah271_struc(m)%noah(i)%nroot = value(LIS_surface(m,mtype)%tile(i)%vegt,1)
        noah271_struc(m)%noah(i)%rsmin = value(LIS_surface(m,mtype)%tile(i)%vegt,2)
        noah271_struc(m)%noah(i)%rgl = value(LIS_surface(m,mtype)%tile(i)%vegt,3)
        noah271_struc(m)%noah(i)%hs = value(LIS_surface(m,mtype)%tile(i)%vegt,4)
        noah271_struc(m)%noah(i)%snup = value(LIS_surface(m,mtype)%tile(i)%vegt,5)
        noah271_struc(m)%noah(i)%z0 = value(LIS_surface(m,mtype)%tile(i)%vegt,6)
        noah271_struc(m)%noah(i)%lai = value(LIS_surface(m,mtype)%tile(i)%vegt,7)
     enddo
     deallocate(value)
     
!begin hack
#if 0 
     n=1
     allocate(temp(LIS_rc%lnc(n),LIS_rc%lnr(n)))

     temp = -9999.0
     do i=1,LIS_rc%ntiles(n)
        c = LIS_domain(n)%tile(i)%col
        r = LIS_domain(n)%tile(i)%row
        temp(c,r) = noah271_struc(n)%noah(i)%rsmin
     enddo
     open(100,file='rsmin_def.1gd4r',form='unformatted',access='direct',&
          recl=LIS_rc%lnc(n)*LIS_rc%lnr(n)*4)
     write(100,rec=1) temp
     close(100)

     temp = -9999.0
     do i=1,LIS_rc%ntiles(n)
        c = LIS_domain(n)%tile(i)%col
        r = LIS_domain(n)%tile(i)%row
        temp(c,r) = noah271_struc(n)%noah(i)%z0
     enddo
     open(100,file='z0_def.1gd4r',form='unformatted',access='direct',&
          recl=LIS_rc%lnc(n)*LIS_rc%lnr(n)*4)
     write(100,rec=1) temp
     close(100)
#endif
#if 0 
     allocate(z0(LIS_rc%lnc(m),LIS_rc%lnr(m)))
     open(111,file='z0_calib.bin',form='unformatted')
     read(111) z0
     close(111)
     
     do i=1,LIS_rc%npatch(m,mtype)
        noah271_struc(m)%noah(i)%z0 = &
             z0(LIS_surface(m,mtype)%tile(i)%col, LIS_surface(m,mtype)%tile(i)%row)
     enddo
     deallocate(z0)
#endif
!end hack

  enddo
end subroutine noah271_setvegparms
