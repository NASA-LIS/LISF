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
! !ROUTINE: snowmodel_finalize
! \label{snowmodel_finalize}
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel 
!
! !INTERFACE:
subroutine snowmodel_finalize()
! !USES:
  use LIS_logMod,   only : LIS_logunit
  use snowmodel_lsmMod
!
! !DESCRIPTION:
!  
!  This routine cleans up the allocated memory structures in SnowModel 
!  
!EOP
  implicit none

  integer :: t

  ! Print a banner when the model run is finished.
  write(LIS_logunit,*)&
       & 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
  write(LIS_logunit,*)&
       & '                 The SnowModel Run Has Finished                '
  write(LIS_logunit,*)&
       & 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'


end subroutine snowmodel_finalize


