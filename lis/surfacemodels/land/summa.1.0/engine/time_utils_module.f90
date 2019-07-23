! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module time_utils_module
USE nrtype
implicit none
private
public::extractTime
public::compjulday
public::compcalday
public::elapsedSec
contains


 ! ******************************************************************************************
 ! public subroutine extractTime: extract year/month/day/hour/minute/second from units string
 ! ******************************************************************************************
 subroutine extractTime(refdate,iyyy,im,id,ih,imin,dsec,err,message)
 implicit none
 ! dummy variables
 character(*),intent(in)    :: refdate             ! units string (time since...)
 integer(i4b),intent(out)   :: iyyy,im,id,ih,imin  ! time (year/month/day/hour/minute)
 real(dp),intent(out)       :: dsec                ! seconds
 integer(i4b),intent(out)   :: err                 ! error code
 character(*),intent(out)   :: message             ! error message
 ! local variables
 integer(i4b)               :: n                   ! length of the string
 integer(i4b)               :: istart,iend         ! position in string
 ! iniitalize error control
 err=0; message="extractTime/"

 ! get the length of the string
 n      = len_trim(refdate)
 ! move to a position in string past the time units (seconds since , days since , hours since )
 istart = index(refdate,'since')  ! get the index at the beginning of the word "since"
 if (istart>0) then ! if the word "since" exists
  iend   = index(refdate(istart:n)," ")
  istart = istart+iend
 else
  istart=1
 end if

 ! get the year
 call extract(refdate(istart:n),"-",iend,iyyy,err,message); if (err/=0) return
 if(iyyy < 1900)then; err=20; message=trim(message)//'year < 1900'; return; end if
 if(iyyy > 2100)then; err=20; message=trim(message)//'year > 2100'; return; end if
 ! get the month
 istart=istart+iend
 call extract(refdate(istart:n),"-",iend,im,err,message);   if (err/=0) return
 if(im <  1)then; err=20; message=trim(message)//'month < 1'; return; end if
 if(im > 12)then; err=20; message=trim(message)//'month > 12'; return; end if
 ! get the day
 istart=istart+iend
 call extract(refdate(istart:n)," ",iend,id,err,message);   if (err/=0) return
 if(id <  1)then; err=20; message=trim(message)//'day < 1'; return; end if
 if(id > 31)then; err=20; message=trim(message)//'day > 31'; return; end if
 ! check if we are at the end of the string
 if (istart+(iend-2)==n) then
  ih=0; imin=0; dsec=0._dp; return
 end if

 ! get the hour (":" at end of hour)
 istart = istart+iend
 if(istart > len_trim(refdate))then; err=20; message=trim(message)//'string does not include hours'; return; end if
 call extract(refdate(istart:n),":",iend,ih,err,message);   if (err/=0) return
 if(ih <  0)then; err=20; message=trim(message)//'hour < 0'; return; end if
 if(ih > 24)then; err=20; message=trim(message)//'hour > 24'; return; end if
 ! get the minute (":" at end of minute)
 istart = istart+iend
 if(istart > len_trim(refdate))then; err=20; message=trim(message)//'string does not include minutes'; return; end if
 call extract(refdate(istart:n),":",iend,imin,err,message); if (err/=0) return
 if(imin <  0)then; err=20; message=trim(message)//'minute < 0'; return; end if
 if(imin > 60)then; err=20; message=trim(message)//'minute > 60'; return; end if

 ! get the second
 istart = istart+iend
 if(istart > len_trim(refdate)) return
 iend   = index(refdate(istart:n)," ")
 read(refdate(istart:n),*) dsec
 !write(*,'(a,i4,1x,4(i2,1x))') 'refdate: iyyy, im, id, ih, imin = ', iyyy, im, id, ih, imin

 contains


  ! ******************************************************************************************
  ! internal subroutine extract: extract substring
  ! ******************************************************************************************
  subroutine extract(substring,cdelim,iend,itemp,err,message)
  implicit none
  ! input
  character(*),intent(in)     :: substring  ! sub-string to process
  character(len=1),intent(in) :: cdelim     ! string delimiter
  ! output
  integer(i4b),intent(out)    :: iend       ! index at the end of desired string
  integer(i4b),intent(out)    :: itemp      ! output date
  integer(i4b),intent(out)    :: err        ! error code
  character(*),intent(out)    :: message    ! error message
  ! initialize error code and message
  err=0; message="extract/"
  ! identify end-point of string
  iend = index(substring,cdelim)
  ! if sub-string does not exist, assume end is at end of string
  if (iend==0) iend=len_trim(substring)+1
  ! convert string to integer
  read(substring(1:iend-1),*,iostat=err) itemp
  ! read error
  if (err/=0) then
   err=20; message=trim(message)//"unexpected characters [string='"//trim(substring)//"']"; return
  end if
  end subroutine extract

 end subroutine extractTime


 ! ***************************************************************************************
 ! public subroutine compjulday: convert date to julian day (units of days)
 ! ***************************************************************************************
 subroutine compjulday(iyyy,mm,id,ih,imin,dsec,&  ! input
                       juldayss,err,message)      ! output
 USE multiconst,only:secprday,secprhour,secprmin  ! seconds in an (day, hour, minute)
 implicit none
 ! input variables
 integer(i4b),intent(in)   :: iyyy,mm,id   ! year, month, day
 integer(i4b),intent(in)   :: ih,imin      ! hour, minute
 real(dp),intent(in)       :: dsec         ! seconds
 ! output
 real(dp),intent(out)      :: juldayss
  integer(i4b),intent(out) :: err          ! error code
  character(*),intent(out) :: message      ! error message
 ! local variables
 integer(i4b)              :: julday       ! julian day
 integer(i4b),parameter    :: igreg=15+31*(10+12*1582)  !IGREG = 588829
 integer(i4b)              :: ja,jm,jy
 real(dp)                  :: jfrac        ! fraction of julian day

 ! initialize errors
 err=0; message="juldayss"

 ! compute julian day
 jy=iyyy
 if (jy.eq.0) then; err=10; message=trim(message)//"noYearZero/"; return; end if
 if (jy.lt.0) jy=jy+1
 if (mm.gt.2) then
  jm=mm+1
 else
  jy=jy-1
  jm=mm+13
 end if
 julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
 if (id+31*(mm+12*iyyy).ge.IGREG) then
  ja=int(0.01*jy)
  julday=julday+2-ja+int(0.25*ja)
 end if

 ! compute fraction of the day
 jfrac = (real(ih,kind(dp))*secprhour + real(imin,kind(dp))*secprmin + dsec) / secprday

 ! and return the julian day, expressed in fraction of a day
 juldayss = real(julday,kind(dp)) + jfrac

 end subroutine compjulday

 ! ***************************************************************************************
 ! public subroutine compgregcal: convert julian day (units of days) to calendar date
 ! source: https://en.wikipedia.org/wiki/Julian_day#Julian_or_Gregorian_calendar_from_Julian_day_number
 ! ***************************************************************************************

 subroutine compcalday(julday,                              & !input
                       iyyy,mm,id,ih,imin,dsec,err,message)   !output
 USE multiconst,only:secprmin  ! seconds in an (day, hour, minute)
 implicit none

 ! input variables	
 real(dp), intent(in)          :: julday       ! julian day

 ! output varibles
 integer(i4b), intent(out)     :: iyyy         ! year
 integer(i4b), intent(out)     :: mm           ! month
 integer(i4b), intent(out)     :: id           ! day
 integer(i4b), intent(out)     :: ih           ! hour
 integer(i4b), intent(out)     :: imin         ! minute
 real(dp),     intent(out)     :: dsec         ! seconds
 integer(i4b), intent(out)     :: err          ! error code
 character(*), intent(out)     :: message      ! error message

 ! local parameters
 integer(i4b),parameter       :: y = 4716
 integer(i4b),parameter       :: j = 1401
 integer(i4b),parameter       :: m = 2
 integer(i4b),parameter       :: n = 12
 integer(i4b),parameter       :: r = 4
 integer(i4b),parameter       :: p = 1461
 integer(i4b),parameter       :: v = 3
 integer(i4b),parameter       :: u = 5
 integer(i4b),parameter       :: s = 153
 integer(i4b),parameter       :: w = 2
 integer(i4b),parameter       :: b = 274277
 integer(i4b),parameter       :: c = -38
 real(dp),parameter           :: hr_per_day = 24.0_dp
 real(dp),parameter           :: min_per_hour = 60.0_dp

 ! local variables
 integer(i4b)          :: f,e,g,h                            ! various step variables from wikipedia
 integer(i4b)          :: step_1a,step_1b,step_1c,step_1d    ! temporary variables for calendar calculations
 real(dp)              :: frac_day  ! fractional day 
 real(dp)              :: remainder ! remainder of modulus operation

 ! initialize errors
 err=0; message="compcalday"
 if(julday<=0)then;err=10;message=trim(message)//"no negative julian days/"; return; end if

 ! step 1
 step_1a = 4*int(julday)+b
 step_1b = step_1a/146097
 step_1c = step_1b*3
 step_1d = step_1c/4

 f = int(julday)+j+step_1d+c

 ! step 2
 e = r * f + v

 ! step 3
 g = mod(e,p)/r

 ! step 4
 h = u * g + w

 ! find day
 id = (mod(h,s))/u + 1

 ! find month
 mm = mod(h/s+m,n)+1

 ! find year
 iyyy = (e/p)-y + (n+m-mm)/n

 ! now find hour,min,second

 frac_day = julday - floor(julday)
 ih = floor((frac_day+1e-9)*hr_per_day)

 remainder = (frac_day+1e-9)*hr_per_day - ih
 imin = floor(remainder*min_per_hour)

 remainder = remainder*min_per_hour - imin
 dsec = nint(remainder*secprmin)

 end subroutine compcalday
 
 ! ***************************************************************************************
 ! public function elapsedSec: calculate difference of two time marks obtained by date_and_time()
 ! *************************************************************************************** 
 function elapsedSec(startTime, endTime)
 USE multiconst,only            :  secprday,secprhour,secprmin        ! seconds in an (day, hour, minute)
 integer(i4b),intent(in)        :: startTime(8),endTime(8)            ! state time and end time
 real(dp)                       :: elapsedSec                         ! elapsed time in seconds
 ! local variables
 integer(i4b)                   :: elapsedDay                         ! elapsed full days
 integer(i4b)                   :: yy                                 ! index of year
 ! number of days of each month
 integer(i4b)                   :: days1(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
 integer(i4b)                   :: days2(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
 
 ! calculate the elapsed time smaller than a day
 elapsedSec = (endTime(8)-startTime(8))*.001_dp + (endTime(7)-startTime(7)) + (endTime(6)-startTime(6))*secprmin + (endTime(5)-startTime(5))*secprhour

 ! check if the run is within the same day otherwise calculate how many days
 if (endTime(1) > startTime(1) .or. endTime(2) > startTime(2) .or. endTime(3) > startTime(3)) then

  elapsedDay = 0
  ! diffenece in year
  do yy = startTime(1), endTime(1) - 1
   elapsedDay = elapsedDay + 365
   if ((mod(yy,4)==0 .and. .not. mod(yy,100)==0) .or. (mod(yy,400)==0)) elapsedDay = elapsedDay + 1
  end do
  if ((mod(startTime(1),4)==0 .and. .not. mod(startTime(1),100)==0) .or. (mod(startTime(1),400)==0)) days1(2) = 29
  if ((mod(endTime(1),4)==0 .and. .not. mod(endTime(1),100)==0) .or. (mod(endTime(1),400)==0)) days2(2) = 29
  ! difference in month 
  if (startTime(2) > 1) elapsedDay = elapsedDay - sum(days1(1:(startTime(2)-1)))
  elapsedDay = elapsedDay - startTime(3) 
  ! difference in day
  if (endTime(2) > 1) elapsedDay = elapsedDay + sum(days2(1:(endTime(2)-1)))
  elapsedDay = elapsedDay + endTime(3)
  ! convert to seconds
  elapsedSec = elapsedSec + elapsedDay * secprday
 end if 
 end function  

end module time_utils_module
