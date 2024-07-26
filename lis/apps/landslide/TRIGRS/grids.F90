!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

! SY: File indentation performed for integration into LIS
! By Rex L. Baum, 1 April 2004
module grids
        integer,allocatable:: pf2(:),indx(:),nxt(:),nv(:),nvu(:)
        integer,allocatable:: dsctr(:),dsc(:),zo(:), itemp(:) ! Added itemp() 2/13/2013 
        real,allocatable:: rikzero(:)
        real,allocatable::  rik(:),rik1(:),ri(:),rizero(:),pf1(:)
        real,allocatable::  pf1_Slope(:) ! SY: For assisting rain assignment from 2-d LIS grid to 1-d TRIGRS array
        real,allocatable:: ri_Mat(:,:) ! SY
        real,allocatable:: temp(:),ro(:),wf(:),ir(:),tfg(:)
        real,allocatable:: zmax(:),slo(:),depth(:)
        real,allocatable:: zfmin(:),fsmin(:),pmin(:) 
        real,allocatable:: elev(:),wtab(:) ! added 4/21/2010, wtab() added 4/19/11
        character (len=4):: grxt ! Added grxt 4/26/2010
        real,allocatable::  pf1_Mat(:,:) ! SY
end module grids
