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
! !ROUTINE: get_gldas
!  \label{get_gldas}
!
! !REVISION HISTORY:
!  19 Sept 2008: Sujay Kumar: Initial Implementation
!
! !INTERFACE:
subroutine get_gldas(n,findex)
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use LIS_timeMgrMod,     only : LIS_tick, LIS_get_nstep
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use gldas_forcingMod,    only : gldas_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, GLDAS forcing. 
!  The temporal frequency of the GLDAS data is 3 hours. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 3 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
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
!    determines the GLDAS data times
!  \item[gldasfile](\ref{gldasfile}) \newline
!    Puts together a time stamped filename for GLDAS data. 
!  \item[read\_gldas](\ref{read_gldas}) \newline
!    reads and spatially interpolates GLDAS data
!  \end{description}
!EOP

  integer :: t,f
  integer :: ferror
  integer :: try
  integer, parameter :: ndays = 10  ! # of days to look back for forcing data
  integer :: order     ! 1 indicates lesser interpolation boundary time
                       ! 2 indicates greater interpolation boundary time
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  integer :: movetime  ! 1=move time2 into time1
  real*8  :: timenow, time1, time2
  real*8  :: dumbtime1, dumbtime2
  real    :: gmt1, gmt2,ts1,ts2
  character(len=LIS_CONST_PATH_LEN) :: name
  integer :: nstep

  gldas_struc(n)%findtime1=0
  gldas_struc(n)%findtime2=0
  movetime=0
  nstep = LIS_get_nstep(LIS_rc,n)

!-----------------------------------------------------------------
! Determine required GLDAS data times 
! (previous assimilation, current & future assimilation hours)
! The adjustment of the hour and the direction will be done
! in the subroutines that generate the names
!-----------------------------------------------------------------
  yr1=LIS_rc%yr    !Time now
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0        
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr1=LIS_rc%yr    !Previous Hour
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=3*((LIS_rc%hr)/3)
  mn1=0
  ss1=0
  ts1=0
  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr2=LIS_rc%yr    !Next Hour
  mo2=LIS_rc%mo
  da2=LIS_rc%da
  hr2=3*((LIS_rc%hr)/3)
  mn2=0
  ss2=0
  ts2=3*60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
  if(timenow.ge.gldas_struc(n)%gldastime2) then
     movetime=1
     gldas_struc(n)%findtime2=1
  endif
  
  if ( nstep.eq.0 .or. nstep.eq.1 .or.LIS_rc%rstflag(n).eq.1 ) then 
     gldas_struc(n)%findtime1=1
     gldas_struc(n)%findtime2=1
     movetime=0
        LIS_rc%rstflag(n) = 0
     endif
     LIS_rc%shortflag=2            !Time averaged SW
     LIS_rc%longflag=2             !Time averaged LW 
     
!-------------------------------------------------------------------
! Establish gldastime1
!-------------------------------------------------------------------
     if (gldas_struc(n)%findtime1==1) then  
        order=1   
        ferror = 0
        try = 0
        ts1 = -24*60*60
        do
           if ( ferror /= 0 ) then
              exit
           end if
           try = try+1
           call gldasfile(name,gldas_struc(n)%gldasdir,&
                yr1,mo1,da1,hr1,gldas_struc(n)%ncold)
           write(LIS_logunit,*) 'Reading ',trim(name)
           call read_gldas(order,n,findex, name,ferror, try)
           if ( ferror == 1 ) then 
!-------------------------------------------------------------------
! successfully retrieved forcing data
!-------------------------------------------------------------------
              gldas_struc(n)%gldastime1=time1
           else  
!-------------------------------------------------------------------
! ferror still=0, so roll back one day
!-------------------------------------------------------------------
              call LIS_tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
           end if
           if ( try > ndays ) then 
              print *, 'ERROR: GLDAS data gap exceeds 10 days on file 1'
              call LIS_endrun
           end if
        end do
     endif
     if(movetime.eq.1) then
        gldas_struc(n)%gldastime1=gldas_struc(n)%gldastime2
        gldas_struc(n)%findtime2=1 
        do f=1,LIS_rc%met_nf(findex)
           do t=1,LIS_rc%ngrid(n)
              gldas_struc(n)%metdata1(f,t)=gldas_struc(n)%metdata2(f,t)
           enddo
        enddo
     endif
     if(gldas_struc(n)%findtime2.eq.1) then 
        order=2   
        ferror = 0
        try = 0
        ts2 = -24*60*60
        do
           if ( ferror /= 0 ) exit
           try = try+1
           call gldasfile(name,gldas_struc(n)%gldasdir,&
                yr2,mo2,da2,hr2,gldas_struc(n)%ncold)
           write(LIS_logunit,*) 'Reading ',trim(name)
           call read_gldas(order,n,findex,name,ferror,try)
           if ( ferror == 1 ) then 
!-------------------------------------------------------------------
! successfully retrieved forcing data
!-------------------------------------------------------------------
              gldas_struc(n)%gldastime2=time2
           else  
!-------------------------------------------------------------------
! ferror still=0, so roll back one day
!-------------------------------------------------------------------
              call LIS_tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
           end if
           if ( try > ndays ) then 
              print *, 'ERROR: GLDAS data gap exceeds 10 days on file 2'
              call LIS_endrun
           end if
        end do
     endif
     
84   format('now',i4,4i3,2x,'pvt ',a22,' nxt ',a22)
!  if ((LIS_rc%gridchange==1).and.(gldas_struc(n)%ncold==288)) then
!     LIS_rc%gridchange=0
!  endif
  return 

end subroutine get_gldas


!BOP
! !ROUTINE: gldasfile
! \label{gldasfile}
!
! !INTERFACE:
subroutine gldasfile(name, gldasdir, yr, mo, da, hr, ncold)
  
  implicit none
! !ARGUMENTS: 
  character(len=*), intent(in)  :: gldasdir
  integer, intent(in)       :: yr,mo,da,hr,ncold
  character(len=*), intent(out) :: name
! !DESCRIPTION:
!   This subroutine puts together GLDAS file name for 
!   3 hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[name]
!   name of the timestamped GLDAS file
!  \item[gldasdir]
!    Name of the GLDAS directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!  \item[ncold]
!   Number of columns (used as to determine the GLDAS resolution)
!  \end{description}
!
!EOP
  integer uyr,umo,uda,uhr,i,c,ii,jj
  character(len=2)  :: initcode
  character(len=10) :: ftime
  character(len=4)  :: fdir
  character(len=17), parameter :: fsuffix = '.force.mosaic.grb'
  
  ii = ncold
  jj = 181

!-------------------------------------------------------------------  
! Make variables for the time used to create the file
! We don't want these variables being passed out
!-------------------------------------------------------------------
  uyr=yr
  umo=mo
  uda=da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
!-------------------------------------------------------------------
!  Determine initcode for the hour of the forecast file
!  If the time is 12 or later the file is time stamped
!  with the next day.  So check for that first
!-------------------------------------------------------------------

  if(uhr<3)then
     initcode = '00'   
  elseif(uhr<6)then
     initcode = '03'
  elseif(uhr<9)then
     initcode = '06'
  elseif(uhr<12)then
     initcode = '09'
  elseif(uhr<15)then
     initcode = '12'
  elseif(uhr<18)then
     initcode = '15'
  elseif(uhr<21)then
     initcode = '18'
  elseif(uhr<24)then
     initcode = '21'
  endif
  
  write(UNIT=fdir,FMT='(i4)') uyr
  
  write(UNIT=ftime,FMT='(i4, i2.2, i2.2, a2)') uyr, umo, uda, initcode

  name = trim(gldasdir) // '/' // fdir // '/' // ftime // fsuffix
  
  return

end subroutine gldasfile



 
