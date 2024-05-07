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
! !ROUTINE: AGRMET_read_pcpclimodata
! \label{AGRMET_read_pcpclimodata}
! 
! !REVISION HISTORY:
!
!    04 dec 97  initial version................capt andrus/dnxm(agromet)
!    31 mar 99  changed retrieval of data from files to new unix based
!               system.  increased size of arrays and loops to full 
!               hemisphere .............................. mr moore/dnxm
!     8 oct 99  ported to ibm sp-2, updated prolog, incorporated
!               FORTRAN 90 features.................capt hidalgo/agrmet
!    21 feb 01  reformatted diagnostic prints.  added intent attribute
!               to arguments...............................mr gayno/dnxm
!    10 jun 02  updated comments to reflect switch from rtneph to
!               cdfs2......................................mr gayno/dnxm
!    1  nov 05  Sujay Kumar, Adopted in LIS
!   26  Dec 07  Marv Freimund, Simplify filename creation
!
! !INTERFACE:
subroutine AGRMET_read_pcpclimodata(n)
! !USES: 
  use LIS_coreMod, only : LIS_rc
  use LIS_LMLCMod, only : LIS_LMLC
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod, only : LIS_logunit, LIS_endrun
  use AGRMET_forcingMod, only :agrmet_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in)         :: n 
!
! !DESCRIPTION: 
! 
!    This routine retrieves monthly climo data from file and interpolate to the
!    current day.   
!
!    method: \newline
!
!    - retrieve climo for current month. \newline
!    - if date is 15th, current months climo is valid, no more work   
!      needs to be done.  \newline
!    - else   \newline
!    - determine if need to interpolate to previous month,  \newline
!      or to the next month. \newline
!    - read appropriate month from file and interpolate \newline
!      values to the current day. \newline
! 
!  The arguments and variables are: 
!  \begin{description}
!   \item[n] 
!     index of the nest
!   \item[days]
!     number of days in each month (for leap year)
!   \item[newday]
!    number of days to interpolate to
!   \item[totday]
!    total days in the current, or previous 
!   \item[i,j]
!    loop indices
!   \item[quad9r]
!    constant = 9999.0
!  \end{description}
!
!    remarks: \newline
!    
!    NOTE: although we use cdfs2 cloud analysis data, we use rtneph
!          climo information in the cdfs2 precip estimate.  this is
!          because we do not have a climo database based on cdfs2
!          as of Jun, 2002.  ideally, once we develop a large
!          cdfs2 database, we should create a new climo.
! 
!   The routines invoked are: 
!   \begin{description}
!    \item[LIS\_isAlarmRinging](\ref{LIS_isAlarmRinging}) \newline
!     check to see if it time to read a new climo data
!    \item[getclimortnfilename](\ref{getclimortnfilename}) \newline
!     filename for the climatological rtneph percent cloud cover 
!    \item[getclimoprcfilename](\ref{getclimoprcfilename}) \newline
!     filename for the precip used for estimate
!    \item[getclimoppdfilename](\ref{getclimoppdfilename}) \newline
!     filename for the precip-per-precip day amount
!    \item[AGRMET\_getcli](\ref{AGRMET_getcli}) \newline
!     retrieves the climo data for any of the above datasets
!    \end{description}
!EOP
  logical                     :: midmonth
  logical                     :: pcpclimoAlarmCheck
  character*100               :: filename
  integer                     :: days(12)
  integer                     :: i,j
  integer                     :: newday
  integer                     :: totday
  real                        :: ratio
  integer                     :: mo1, mo2
  real, parameter             :: quad9r = 9999.0   
  data days    / 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /  

  midmonth = .true. 
  pcpclimoAlarmCheck = LIS_isAlarmRinging(LIS_rc,&
       agrmet_struc(n)%pcpclimoAlarmTime,&
       "monthly",midmonth)
  
  if((agrmet_struc(n)%cdfs2swch.eq.1).or.&
       (agrmet_struc(n)%clswch.eq.1)) then         
     mo1 = LIS_rc%mo
     ! JVG: I have changed the valid day to 16 to be consistent with 
     ! LIS' time manager.
     if(LIS_rc%da.ne.16) then 
        if(LIS_rc%da.lt.16) then 
           if(mo1.eq.1) then 
              mo2 = 12
           else
              mo2 = mo1-1
           endif
           newday = 16-LIS_rc%da
           totday = days(mo2)
        else
           if(mo1.eq.12) then 
              mo2 = 1
           else
              mo2 = mo1+1
           endif
           newday = LIS_rc%da-16
           totday = days(mo2)
        endif
     else
        mo2 = mo1
        ! Chris' problem was due to a bad ratio value which results from
        ! newday and totday being undefined at da=16 (was da=15).  This will 
        ! apply full weight to agrmet_struc(n)%cliprc1(i,j).
        newday = 0
        totday = days(mo2)
     endif
     if(pcpclimoAlarmCheck) then                 
        call getclimortnfilename(filename,agrmet_struc(n)%climodir,mo1)
        call AGRMET_getcli(n, filename, 1, agrmet_struc(n)%clirtn1)
        call getclimortnfilename(filename,agrmet_struc(n)%climodir,mo2)
        call AGRMET_getcli(n, filename, 1, agrmet_struc(n)%clirtn2)
     endif
     do j=1,LIS_rc%lnr(n)
        do i=1,LIS_rc%lnc(n)
           if(LIS_LMLC(n)%landmask(i,j).ne.0) then 
              ratio = float(newday)/float(totday)                    
              agrmet_struc(n)%clirtn(i,j) = (ratio*&
                   (agrmet_struc(n)%clirtn2(i,j) - &
                   agrmet_struc(n)%clirtn1(i,j))+ &
                   agrmet_struc(n)%clirtn1(i,j))
           else
              agrmet_struc(n)%clirtn(i,j) = quad9r
           endif
        enddo
     enddo
        
     if(pcpclimoAlarmCheck) then                 
        call getclimoprcfilename(filename,agrmet_struc(n)%climodir,mo1)
        call AGRMET_getcli(n, filename, 1, agrmet_struc(n)%cliprc1)
        call getclimoprcfilename(filename,agrmet_struc(n)%climodir,mo2)
        call AGRMET_getcli(n, filename, 1, agrmet_struc(n)%cliprc2)
     endif
     do j=1,LIS_rc%lnr(n)
        do i=1,LIS_rc%lnc(n)
           if(LIS_LMLC(n)%landmask(i,j).ne.0) then 
              ratio = float(newday)/float(totday)                    
              agrmet_struc(n)%cliprc(i,j) = (ratio* &
                   (agrmet_struc(n)%cliprc2(i,j) - &
                   agrmet_struc(n)%cliprc1(i,j))+ &
                   agrmet_struc(n)%cliprc1(i,j))
              if(agrmet_struc(n)%cliprc(i,j).lt.0) then 
                 write(LIS_logunit,*) 'ERROR - NEGATIVE PRECIP CLIMO',     &
                                      agrmet_struc(n)%cliprc(i,j) , ratio, &
                                      agrmet_struc(n)%cliprc2(i,j),        &
                                      agrmet_struc(n)%cliprc1(i,j),        &
                                      agrmet_struc(n)%cliprc1(i,j)
                 write(LIS_logunit,*) 'Stopping program '
                 call LIS_endrun
              endif              
           else
              agrmet_struc(n)%cliprc(i,j) = quad9r
           endif
        enddo
     enddo
     if(agrmet_struc(n)%cdfs2swch.eq.1) then 
        if(pcpclimoAlarmCheck) then                 
           call getclimoppdfilename(filename,agrmet_struc(n)%climodir,mo1)
           call AGRMET_getcli(n, filename, 1, agrmet_struc(n)%clippd1)
           call getclimoppdfilename(filename,agrmet_struc(n)%climodir,mo2)
           call AGRMET_getcli(n, filename, 1, agrmet_struc(n)%clippd2)
        endif
        do j=1,LIS_rc%lnr(n)
           do i=1,LIS_rc%lnc(n)
              if(LIS_LMLC(n)%landmask(i,j).ne.0) then 
                 ratio = float(newday)/float(totday)                    
                 agrmet_struc(n)%clippd(i,j) = (ratio*&
                      (agrmet_struc(n)%clippd2(i,j) - &
                      agrmet_struc(n)%clippd1(i,j))+ &
                      agrmet_struc(n)%clippd1(i,j)) 
              else
                 agrmet_struc(n)%clippd(i,j) = quad9r
              endif
           enddo
        enddo
     endif
  end if
end subroutine AGRMET_read_pcpclimodata

!BOP
! 
! !ROUTINE: getclimortnfilename
!  \label{getclimortnfilename}
! 
! !INTERFACE: 
subroutine getclimortnfilename(fname,dir,mo)

  implicit none
! !ARGUMENTS:   
  character(*)        :: fname
  character(*)        :: dir
  integer, intent(in) :: mo

! 
! !DESCRIPTION: 
!  This routines generates the name of the climatalogical rtneph
!  percent cloud cover file, by 
!  appending the root directory and the hemisphere information. 
!  The name of the file is expected to be: 
!  <dir>/clim\_rtn\_<mon>\_<hh>, where mon is the 3-character
!  name of the month and hh is 'nh' or 'sh', for
!  northern and southern hemisphere, respectively. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[fname]
!    created filename
!   \item[dir]
!    path to the directory containing the climo file
!   \item[mo]
!    integer value of month (1-12)
!  \end{description}
!EOP

  character(3), parameter :: months(12) = (/'jan','feb','mar','apr','may','jun',&
                                            'jul','aug','sep','oct','nov','dec'/)

  fname = trim(dir) // '/cli_rtn_' // months(mo) // '.1gd4r'

end subroutine getclimortnfilename

!BOP
! 
! !ROUTINE: getclimoprcfilename
! \label{getclimoprcfilename}
! 
! !INTERFACE: 
subroutine getclimoprcfilename(fname,dir,mo)

  implicit none
! !ARGUMENTS:  
  character(*)        :: fname
  character(*)        :: dir
  integer, intent(in) :: mo
! 
! !DESCRIPTION: 
!  This routines generates the name of the 3-hour climo precip 
!  file, by appending the root directory and the hemisphere information. 
!  The name of the file is expected to be: 
!  <dir>/clim\_prc\_<mon>\_<hh>, where mon is the 3-character
!  name of the month and hh is 'nh' or 'sh', for
!  northern and southern hemisphere, respectively. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[fname]
!    created filename
!   \item[dir]
!    path to the directory containing the climo file
!   \item[mo]
!    integer value of month (1-12)
!  \end{description}
!EOP

  character(3), parameter :: months(12) = (/'jan','feb','mar','apr','may','jun',&
                                            'jul','aug','sep','oct','nov','dec'/)

  fname = trim(dir) // '/cli_prec_' // months(mo) // '.1gd4r'

end subroutine getclimoprcfilename

!BOP
! 
! !ROUTINE: getclimoppdfilename
! \label{getclimoppdfilename}
! 
! !INTERFACE: 
subroutine getclimoppdfilename(fname,dir,mo)

  implicit none
! !ARGUMENTS:   
  character(*)        :: fname
  character(*)        :: dir
  integer, intent(in) :: mo
! 
! !DESCRIPTION: 
!  This routines generates the name of the precip-per-precip day amount
!  file, by appending the root directory and the hemisphere information. 
!  The name of the file is expected to be: 
!  <dir>/clim\_ppd\_<mon>\_<hh>, where mon is the 3-character
!  name of the month and hh is 'nh' or 'sh', for
!  northern and southern hemisphere, respectively. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[fname]
!    created filename
!   \item[dir]
!    path to the directory containing the climo file
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[mo]
!    integer value of month (1-12)
!  \end{description}
!EOP
  character(3), parameter :: months(12) = (/'jan','feb','mar','apr','may','jun',&
                                            'jul','aug','sep','oct','nov','dec'/)

  fname = trim(dir) // '/cli_ppd_' // months(mo) //'.1gd4r'

end subroutine getclimoppdfilename
