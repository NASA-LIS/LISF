!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_rhoneAGG
!   \label{get_rhoneAGG}
!
! !REVISION HISTORY:
!  5 Nov 2003: Dave Mocko, Initial Specification
!
! !INTERFACE:
subroutine get_rhoneAGG(n,findex)
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use LIS_timeMgrMod,     only : LIS_get_nstep, LIS_tick
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use rhoneAGG_forcingMod,   only : rhoneAGG_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!
! !DESCRIPTION:
!  Opens, reads, and interpolates RHONEAGG forcing.
!
!    TIME1 = most recent past data
!    TIME2 = nearest future data
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!      Determines RHONEAGG data times
!  \item[rhoneAGGfile](\ref{rhoneAGGfile}) \newline
!      Puts together appropriate file name for 3 hour intervals
!  \item[readrhoneAGG](\ref{readrhoneAGG}) \newline
!      Interpolates RHONEAGG data to LIS grid
!  \end{description}
!EOP

  integer, parameter :: ndays = 10 ! # days to look back for forcing data
  integer :: try, ferror
  integer :: c,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real*8 :: time1,time2,dumbtime1,dumbtime2
  real*8 :: timenow
  character*80 :: name
  real :: gmt1,gmt2,ts1,ts2
  integer :: movetime       ! 1=move time 2 data into time 1
  integer :: nforce         ! RHONEAGG forcing file time, # forcing variables
  integer :: nstep

  nstep = LIS_get_nstep(LIS_rc, n)
!-----------------------------------------------------------------------
! Determine the correct number of forcing variables
!-----------------------------------------------------------------------
  nforce = rhoneAGG_struc(n)%nmif
  
  rhoneAGG_struc(n)%findtime1=0
  rhoneAGG_struc(n)%findtime2=0
  movetime=0
!-----------------------------------------------------------------------
! Determine Required RHONEAGG Data Times
! (The previous hour & the future hour)
!-----------------------------------------------------------------------
  yr1=LIS_rc%yr              !Time now
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr1=LIS_rc%yr              !Previous Hour
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=3*((LIS_rc%hr)/3)
  mn1=0
  ss1=0
  ts1=0
  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr2=LIS_rc%yr              !Next Hour
  mo2=LIS_rc%mo
  da2=LIS_rc%da
  hr2=3*((LIS_rc%hr)/3)
  mn2=0
  ss2=0
  ts2=3*60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
  
  if (timenow.gt.rhoneAGG_struc(n)%rhoneAGGtime2) then
     movetime=1
     rhoneAGG_struc(n)%findtime2=1
  endif
  
  if ((nstep.eq.0).or.(nstep.eq.1).or.(LIS_rc%rstflag(n).eq.1)) then
     rhoneAGG_struc(n)%findtime1=1
     rhoneAGG_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
  LIS_rc%shortflag=2         !Time averaged SW
  LIS_rc%longflag=2          !Time averaged LW
  
  !-----------------------------------------------------------------------
  ! Establish rhoneAGGtime1
  !-----------------------------------------------------------------------
  if (rhoneAGG_struc(n)%findtime1.eq.1) then
     order = 1
     ferror = 0
     try = 0
     ts1 = -24*60*60
     do
        if (ferror.ne.0) then
           exit
        endif
        try = try+1
        call rhoneAGGfile(name,rhoneAGG_struc(n)%rhoneAGGdir,yr1,mo1,da1,hr1,     &
             rhoneAGG_struc(n)%ncold)
        call readrhoneAGG(n, order,findex, name,LIS_rc%tscount(n),ferror)
        if (ferror.eq.1) then
!-----------------------------------------------------------------------
! successfully retrieved forcing data
!-----------------------------------------------------------------------
           rhoneAGG_struc(n)%rhoneAGGtime1=time1
        else
!-----------------------------------------------------------------------
! ferror still=0, so roll back one day
!-----------------------------------------------------------------------
           call LIS_tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        endif
        if (try.gt.ndays) then
           write(LIS_logunit,*)'ERROR: RHONEAGG data gap exceeds 10 days on file 1'
           call LIS_endrun
        endif
     enddo
  endif
  if (movetime.eq.1) then
     rhoneAGG_struc(n)%rhoneAGGtime1=rhoneAGG_struc(n)%rhoneAGGtime2
     rhoneAGG_struc(n)%findtime2 =1
     do f=1,nforce
        do c=1,LIS_rc%ngrid(n)
           rhoneAGG_struc(n)%metdata1(f,c)=rhoneAGG_struc(n)%metdata2(f,c)
        enddo
     enddo
  endif
  
!-----------------------------------------------------------------------
! Establish rhoneAGGtime2
!-----------------------------------------------------------------------
  if (rhoneAGG_struc(n)%findtime2.eq.1) then
     order = 2
     ferror = 0
     try = 0
     ts2 = -24*60*60
     do
        if (ferror.ne.0) exit
        try = try+1
        call rhoneAGGfile(name,rhoneAGG_struc(n)%rhoneAGGdir,yr2,mo2,da2,hr2,     &
             rhoneAGG_struc(n)%ncold)
        call readrhoneAGG(n, order,findex, name,LIS_rc%tscount(n),ferror)
        if (ferror.eq.1) then
!-----------------------------------------------------------------------
! successfully retrieved forcing data
!-----------------------------------------------------------------------
           rhoneAGG_struc(n)%rhoneAGGtime2=time2
        else
!-----------------------------------------------------------------------
! ferror still=0, so roll back one day
!-----------------------------------------------------------------------
           call LIS_tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        endif
        if (try.gt.ndays) then
           print *,'ERROR: RHONEAGG data gap exceeds 10 days on file 2'
           call LIS_endrun
        endif
     enddo
  endif
  
84 format('now',i4,4i3,2x,'pvt ',a22,' nxt ',a22)
  
  return
end subroutine get_rhoneAGG

!BOP
! !ROUTINE: rhoneAGGfile
! \label{rhoneAGGfile}
!
! !DESCRIPTION:
!  This subroutine puts together RHONEAGG file name
!
! !INTERFACE:
subroutine rhoneAGGfile(name,rhoneAGGdir,yr,mo,da,hr,ncold)

  implicit none
      
! !INPUT PARAMETERS:
  character*40, intent(in) :: rhoneAGGdir
  integer, intent(in)      :: yr,mo,da,hr,ncold
! !OUTPUT PARAMETERS:
  character*80, intent(out) :: name
!EOP
  integer uyr,umo,uda,uhr,i,c,ii,jj
  character(len=2) :: initcode
  character*1 fbase(80),fsubs(80)
  character*1 ftime(10),fdir(8)
  character(LEN=100) :: temp
  
  ii = ncold
  jj = 6

!-----------------------------------------------------------------------
! Make variables for the time used to create the file
! We don't want these variables being passed out
!-----------------------------------------------------------------------
  uyr=yr
  umo=mo
  uda=da
  uhr = 3*(hr/3)            !hour needs to be a multiple of 3 hours
!-----------------------------------------------------------------------
!  Determine initcode for the hour of the forecast file
!  If the time is 12 or later the file is time stamped
!  with the next day.  So check for that first
!-----------------------------------------------------------------------
  if (uhr.lt.3) then
     initcode = '00'
  elseif (uhr.lt.6) then
     initcode = '03'
  elseif (uhr.lt.9) then
     initcode = '06'
  elseif (uhr.lt.12) then
     initcode = '09'
  elseif (uhr.lt.15) then
     initcode = '12'
  elseif (uhr.lt.18) then
     initcode = '15'
  elseif (uhr.lt.21) then
     initcode = '18'
  elseif (uhr.lt.24) then
     initcode = '21'
  endif

  write(UNIT=temp,FMT='(A40)') rhoneAGGdir
  read(UNIT=temp,FMT='(80A1)') (fbase(i),i=1,80)

  write(UNIT=temp,FMT='(a1,i4,i2,a1)') '/',uyr,umo,'/'
  read(UNIT=temp,FMT='(8A1)') fdir
  do i=1,8
     if (fdir(i).eq.(' ')) fdir(i)='0'
  enddo
  
  write(UNIT=temp,FMT='(i4,i2,i2,a2)') uyr,umo,uda,initcode
  read(UNIT=temp,FMT='(10A1)') ftime
  do i=1,10
     if (ftime(i).eq.(' ')) ftime(i)='0'
  enddo
  
  if (ncold.eq.360) then
     write(UNIT=temp,FMT='(A6)') '.RHONE'
     read(UNIT=temp,FMT='(80A1)') (fsubs(i),i=1,6)
  else
     write(UNIT=temp,FMT='(A6)') '.RHONE'
     read(UNIT=temp,FMT='(80A1)') (fsubs(i),i=1,6)
  endif
  c=0
  do i=1,80
     if ((fbase(i).eq.(' ')).and.(c.eq.0)) c=i-1
  enddo
  
  if (ncold.eq.360) then
     write(UNIT=temp,FMT='(80a1)')(fbase(i),i=1,c),(fdir(i),i=1,8),&
          (ftime(i),i=1,10),(fsubs(i),i=1,6)
  else
     write(UNIT=temp,FMT='(80a1)')(fbase(i),i=1,c),(fdir(i),i=1,8),&
          (ftime(i),i=1,10),(fsubs(i),i=1,6)
  endif
  
  read(UNIT=temp,FMT='(a80)') name
  return
end subroutine rhoneAGGfile

