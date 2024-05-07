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
! !ROUTINE: canhtset
! \label{canhtset}
!
! !REVISION HISTORY:
!  15 Nov 2002: Jon Gottschalck; Initial code
!
! !INTERFACE:
subroutine canhtset (n)
! !USES:
  use clm2_lsmMod          ! CLM tile variables
  use clm2_varpar     , only : maxpatch
  use LIS_coreMod, only  : LIS_rc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
! !DESCRIPTION:
!  This subroutine reads in canopy height information into CLM. 
!  Expects a file with the top and bottom canopy heights for each 
!  veg type. 
!
!EOP

  integer :: i,t
  real    :: pthtop(maxpatch),pthbot(maxpatch)  ! Canopy top and bottom height for the 13 UMD vegetation types


!---------------------------------------------------------------------
! Open canopy heights file and read into temporary variables
!---------------------------------------------------------------------
  open(57,file=clm2_struc(n)%clm_chtfile,status='old')
  read(57,*) (pthtop(i),i=1,maxpatch)
  read(57,*) (pthbot(i),i=1,maxpatch)
  close(57)
!---------------------------------------------------------------------
! Assign canopy top and height values for each 
! vegetation type to CLM tiles
!---------------------------------------------------------------------
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(clm2_struc(n)%clm(t)%itypveg.eq.0) then !noveg
        clm2_struc(n)%clm(t)%htop = 0.0
        clm2_struc(n)%clm(t)%hbot = 0.0
     else
        clm2_struc(n)%clm(t)%htop = pthtop(clm2_struc(n)%clm(t)%itypveg)
        clm2_struc(n)%clm(t)%hbot = pthbot(clm2_struc(n)%clm(t)%itypveg)
     endif
  enddo

end subroutine canhtset


