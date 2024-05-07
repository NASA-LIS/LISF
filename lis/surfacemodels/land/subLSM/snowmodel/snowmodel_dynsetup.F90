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
! !ROUTINE: snowmodel_dynsetup
! \label{snowmodel_dynsetup}
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel 
!
! !INTERFACE:
subroutine snowmodel_dynsetup(n)
! !USES: 
  use LIS_coreMod, only : LIS_rc,LIS_domain, LIS_surface
  use LIS_logMod,  only : LIS_logunit, LIS_endrun
!  use LIS_snowMod, only : LIS_snow_struc
!  use LIS_vegDataMod, only : LIS_gfrac, LIS_lai, LIS_roughness
  use snowmodel_lsmMod, only : snowmodel_struc
  use LIS_timeMgrMod, only : LIS_date2time,LIS_tick
! 
! !DESCRIPTION: 
!  This routine sets up the time-dependent variables in SnowModel.
! 
!EOP   

  implicit none
  integer, intent(in) :: n

  write(LIS_logunit,*) '[INFO] Call to the SnowModel dynamic setup routine ...'

end subroutine snowmodel_dynsetup
