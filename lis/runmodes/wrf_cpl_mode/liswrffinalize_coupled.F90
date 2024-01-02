!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine liswrffinalize_coupled()
! 
!BOP
! !ROUTINE: liswrffinalize_coupled
! 
! !REVISION HISTORY: 
! 14 Aug 2014  James Geiger  Initial Specification
! !USES: 
  use LISWRFGridCompMod,     only : LISWRF_import
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
!  a coupled LIS-WRF run. The routine writes LIS output and restart files
!  valid at the end of the run.
! 
!EOP  

  integer :: n 

  n = 1
  call LIS_timemgr_set(LIS_rc,LISWRF_import(n)%yr,&
                       LISWRF_import(n)%mo,LISWRF_import(n)%da,&
                       LISWRF_import(n)%hr,LISWRF_import(n)%mn,&
                       LISWRF_import(n)%ss,0,0.0)

  LIS_rc%endtime = 1

  do n = 1, LIS_rc%nnest
     call LIS_surfaceModel_output(n)
     call LIS_surfaceModel_writerestart(n)
     call LIS_irrigation_output(n) ! EMK
  enddo

  write(unit=LIS_logunit,fmt=*) 'LIS cycle completed'

 flush(LIS_logunit)

end subroutine liswrffinalize_coupled
