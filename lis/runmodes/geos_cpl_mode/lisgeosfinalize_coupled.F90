!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine lisgeosfinalize_coupled()
! 
!BOP
! !ROUTINE: lisgeosfinalize_coupled
! 
! !REVISION HISTORY: 
! 9 Nov 2026  Ronnie Abolafia-Rosenzweig  Initial Specification
! !USES: 
  use LISGEOSGridCompMod,     only : LISGEOS_import
  use LIS_coreMod,           only : LIS_rc
  use LIS_timeMgrMod,        only : LIS_timemgr_set
  use LIS_surfaceModelMod,   only : LIS_surfaceModel_output,         &
                                    LIS_surfaceModel_writerestart
  use LIS_logMod,            only : LIS_verify, LIS_logunit
  use LIS_irrigationMod,     only : LIS_irrigation_output ! EMK
!EOP
 
  implicit none 
! !ARGUMENTS:   
  !integer, intent(in)   :: n 
!
! !DESCRIPTION: 
!  This routine defines the set of steps required from LIS to finalize
!  a coupled LIS-GEOS run. The routine writes LIS output and restart files
!  valid at the end of the run.
! 
!EOP  

  integer :: n 

  n = 1
  call LIS_timemgr_set(LIS_rc,LISGEOS_import(n)%yr,&
                       LISGEOS_import(n)%mo,LISGEOS_import(n)%da,&
                       LISGEOS_import(n)%hr,LISGEOS_import(n)%mn,&
                       LISGEOS_import(n)%ss,0,0.0)

  LIS_rc%endtime = 1

  do n = 1, LIS_rc%nnest
     call LIS_surfaceModel_output(n)
     call LIS_surfaceModel_writerestart(n)
     call LIS_irrigation_output(n) ! EMK
  enddo

  write(unit=LIS_logunit,fmt=*) 'LIS cycle completed'

 flush(LIS_logunit)

end subroutine lisgeosfinalize_coupled
