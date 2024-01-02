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
!BOP
!
! !ROUTINE: read_MODIS_RT_sai
! \label{read_MODIS_RT_sai}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!
! !INTERFACE:
subroutine read_MODIS_RT_sai(n, time1, array)
! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_localPet
  use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_endrun
  use LIS_vegDataMod,    only : LIS_sai
  use LIS_fileIOMod,     only : LIS_readData
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n
  type(ESMF_Time)     :: time1, newtime
  type(ESMF_TimeInterval) :: deltaT
  real, intent(inout) :: array(LIS_rc%ntiles(n))

! !DESCRIPTION:
!  This subroutine retrieves the LAI data for the 
!  specified date and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[mo]
!   time index (month or quarter)
!  \item[array]
!   output field with the retrieved LAI 
!  \end{description}
!
!EOP     
  integer*1, allocatable      :: tmparr(:,:)
  character(len=LIS_CONST_PATH_LEN) :: filename2
  integer :: yr, mo, dy 
  integer :: t, it, rc
  integer :: ftn
  character * 4 :: cyr
  character * 2 :: cmo, cdy 
  real        :: sai_gridDesc(6), sai
  logical     :: file_exists

! MODIS RT LAI file example: 
!  input/LS_PARAMETERS/MODIS/MCD15A2.005/2004/lai.03.13.1gd1c

  allocate(tmparr(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  call ESMF_TimeGet(time1,yy=yr, mm=mo, dd=dy)
  Do it=1, 8  ! search forward for 8 days 
  write(cyr, '(I4.4)') yr
  write(cmo, '(I2.2)') mo
  write(cdy, '(I2.2)') dy
  filename2 = trim(LIS_sai(n)%saifile)//'/'//cyr//'/'//'sai.'//cmo//'.'//cdy//'.1gd1c' 

  inquire(file=trim(filename2), exist=file_exists)
    if(.not.file_exists) then 
       call ESMF_TimeIntervalSet(deltaT,d=it,rc=rc)
       newtime = time1 + deltaT
       call ESMF_TimeGet(newtime, yy=yr, mm=mo, dd=dy)
    else 
       exit
    end if
  End Do

  if(.not.file_exists) then   ! did not find after 8 steps 
     write(LIS_logunit,*) 'SAI map ',trim(filename2),' not found'
     write(LIS_logunit,*) 'Program stopping ...'
     call LIS_endrun
  endif

  ftn = LIS_getNextUnitNumber()
  open(ftn, file=filename2, access='direct',status='old', &
       form="unformatted", recl=1)

  sai_gridDesc(1) = LIS_rc%gridDesc(n,34)
  sai_gridDesc(2) = LIS_rc%gridDesc(n,35)
  sai_gridDesc(3) = LIS_rc%gridDesc(n,37)
  sai_gridDesc(4) = LIS_rc%gridDesc(n,38)
  sai_gridDesc(5) = LIS_rc%gridDesc(n,39)
  sai_gridDesc(6) = LIS_rc%gridDesc(n,40)

  call read_lai_1gd1c(n,ftn,sai_gridDesc,tmparr)

  call LIS_releaseUnitNumber(ftn)

  do t=1,LIS_rc%ntiles(n)
     sai = tmparr(LIS_domain(n)%tile(t)%col,LIS_domain(n)%tile(t)%row) * 0.1
     ! reset undef to 0 
     if (sai .LT.0 .or. sai .GT. 10) sai = 0 
     array(t) = sai 
  enddo
  deallocate(tmparr)
  write(LIS_logunit,*)'Read SAI File ',trim(filename2)


end subroutine read_MODIS_RT_sai

