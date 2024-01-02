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
! !ROUTINE: create_gdasT1534filename
!  \label{create_gdasT1534filename}
!
! !REVISION HISTORY:
!  20 June 2014: Sujay Kumar; initial implementation
! !INTERFACE:
subroutine create_gdasT1534filename(name, gdasdir, yr, mo, da, hr)
! !USES:

  use LIS_timeMgrMod,     only : LIS_tick

  implicit none
! !ARGUMENTS: 
  character(len=*)             :: name
  character(len=*), intent(in) :: gdasdir
  integer                      :: yr
  integer                      :: mo
  integer                      :: da
  integer                      :: hr
!  
! !DESCRIPTION:
!
!EOP
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss
  real    :: ts1
  integer :: ghh, gff
  integer :: remainder
  integer :: doy
  real    :: gmt
  real*8  :: dumbtime
  character(len=2) :: initcode, fcstcode
  character(len=13) :: fdir
  character(len=22) :: fsubs
!=== End Variable Definition ===============

!-----------------------------------------------------------------
!  Make variables for the time used to create the file
!  We don't want these variables being passed out
!-----------------------------------------------------------------
  dumbtime = 0.
  doy = 1
  gmt = 0.
  uyr = yr
  umo = mo
  uda = da
  uhr = hr
  umn = 0
  uss = 0
  ts1 = -1.0*60*60 !roll back date.

  ghh = hr 
  if ( ghh .EQ. 0 ) then
     call LIS_tick(dumbtime, doy, gmt, uyr, umo, uda, uhr, umn, uss, ts1)
     ghh = 24
  endif

  gff = modulo(ghh,6) 
  if ( gff .EQ. 0 ) gff =  6
  ghh = ghh - gff
 
  write(initcode,'(i2.2)') ghh
  write(fcstcode,'(i2.2)') gff

  write(UNIT=fdir, fmt='(a5, i4, i2.2, i2.2)') 'gdas.', uyr, umo, uda

  fsubs = 'gdas1.t' // initcode // 'z.sfluxgrbf' // fcstcode

  name = trim(gdasdir) // '/' // fdir // '/' // fsubs

end subroutine create_gdasT1534filename
