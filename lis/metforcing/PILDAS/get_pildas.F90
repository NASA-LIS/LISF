!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_pildas
! \label{get_pildas}
!  
!
! !REVISION HISTORY:
!  25 Apr 2013: Sujay Kumar, initial version
! 14 Jul 2016: Mahdi Navari - Modified for PILDAS
!
! !INTERFACE:
subroutine get_pildas(n,findex)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_metforcingMod
  use pildas_forcingMod

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n 
  integer, intent(in)          :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly PILDAS forcing.
!  At the beginning of a simulation, the code reads the most
!  recent past data (nearest 1 hour interval), and the nearest
!  future data.  These two datasets are used to temporally
!  interpolate the data to the current model timestep.
!  The strategy for missing data is to go backwards up
!  to 10 days to get forcing at the same time of day.
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
!    call to advance or retract time
!  \item[pildasfile](\ref{pildasfile}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_pildas](\ref{read_pildas}) \newline
!    call to read the PILDAS data and perform spatial interpolation 
!  \end{description}
!EOP
  integer :: c,f,ferror,try
  integer :: order
  real*8  :: time1,time2,timenow
  real*8  :: dtime1, dtime2
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  character*80 :: name
  real :: gmt1,gmt2,ts1,ts2
  integer:: movetime     ! 1=move time 2 data into time 1  

  try=-999

  pildas_struc(n)%findtime1=0
  pildas_struc(n)%findtime2=0
  movetime=0

  yr1 = LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=LIS_rc%ss
  ts1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  if(LIS_rc%ts.gt.3600) then 
     write(LIS_logunit,*) 'WARNING: The model timestep is > forcing data timestep'
     write(LIS_logunit,*) 'LIS does not support this mode currently'
     write(LIS_logunit,*) 'Program stopping ...'
     call LIS_endrun()
  endif

  if(mod(nint(LIS_rc%ts),3600).eq.0) then 
     if(timenow.ge.pildas_struc(n)%pildastime2) then 
        yr1=LIS_rc%yr
        mo1=LIS_rc%mo
        da1=LIS_rc%da
        hr1=LIS_rc%hr
        mn1=30
        ss1=LIS_rc%ss
        ts1=-3600.0
        
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        
        yr2=LIS_rc%yr    !next hour
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=LIS_rc%hr
        mn2=30
        ss2=LIS_rc%ss
        ts2=0.0
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        movetime = 1
        pildas_struc(n)%findtime2 = 1
     endif
  else
     if(timenow.ge.pildas_struc(n)%pildastime2) then 
        yr1=LIS_rc%yr
        mo1=LIS_rc%mo
        da1=LIS_rc%da
        if(LIS_rc%mn.ge.30) then 
           hr1=LIS_rc%hr+1
        else
           hr1=LIS_rc%hr
        endif
        mn1=30
        ss1=0
        ts1=-60*60
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        yr2=LIS_rc%yr    !next hour
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        if(LIS_rc%mn.ge.30) then 
           hr2=LIS_rc%hr+1
        else
           hr2=LIS_rc%hr
        endif
        mn2=30
        ss2=0
        ts2=0
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        movetime = 1
        pildas_struc(n)%findtime2 = 1
     endif
  endif
    
  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then    !beginning of the run	
     pildas_struc(n)%findtime1=1
     pildas_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
  if(movetime.eq.1) then
     pildas_struc(n)%pildastime1 = pildas_struc(n)%pildastime2
     pildas_struc(n)%gmt1 = pildas_struc(n)%gmt2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           pildas_struc(n)%metdata1(f,c)=pildas_struc(n)%metdata2(f,c)
        enddo
     enddo
  endif    !end of movetime=1
  
  if(pildas_struc(n)%findtime1.eq.1) then

     ferror=0
     try=0  
     ts1=-60*60*24

     do 
        if ( ferror /= 0 ) exit
        try=try+1
        
        call pildasfile(name,pildas_struc(n)%pildasdir,&
             pildas_struc(n)%version,&
             yr1,mo1,da1,hr1,mn1)
        
        write(unit=LIS_logunit,fmt=*)'getting file1.. ',name
        order = 1
        call read_pildas(n,findex,order,name,ferror)
        if(ferror.ge.1) then
           pildas_struc(n)%pildastime1=time1
           pildas_struc(n)%gmt1=gmt1
        endif
        call LIS_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(*,*)'error: Pildas data gap exceeds 10 days on file 1'
           stop
        endif
     enddo
  endif   !end of LIS_rc%findtime=1	   	


  if(pildas_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        
        call pildasfile(name,pildas_struc(n)%pildasdir,&
             pildas_struc(n)%version,&
             yr2,mo2,da2,hr2,mn2)
        
        write(unit=LIS_logunit,fmt=*)'getting file2.. ',name
        order = 2
        call read_pildas(n,findex,order,name,ferror)
        if(ferror.ge.1) then
           pildas_struc(n)%pildastime2=time2
           pildas_struc(n)%gmt2=gmt2
        endif
        call LIS_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(*,*)'error: Pildas data gap exceeds 10 days on file 2'
           stop
        endif
     enddo
  endif
end subroutine get_pildas

!BOP
! !ROUTINE: pildasfile
! \label{pildasfile}
!
! !INTERFACE:
subroutine pildasfile(name,pildasdir,version,yr,mo,da,hr,mn)
  
  implicit none
! !ARGUMENTS: 
  character(len=*), intent(out) :: name
  character(len=*), intent(in)  :: pildasdir
  integer, intent(in)           :: version
  integer, intent(in)           :: yr,mo,da,hr,mn
! !DESCRIPTION:
!  This subroutine puts together PILDAS file name for 
!  1 hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[name]
!   name of the timestamped GEOS5 forecast file
!  \item[pildasdir]
!    Name of the GEOS5 forecast data directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!  \item[mn]
!   minute
!  \end{description}
!
!EOP
  character*3 :: months(12)
  character*4 :: ftime1
  character*8 :: ftime2
  character*4 :: ftime3
  character*3 :: fversion
  character*2 :: fmo

  data months /'jan','feb','mar','apr','may','jun','jul','aug','sep','oct',&
       'nov','dec'/

  write(unit=ftime1, fmt='(i4.4)') yr
  write(unit=ftime2, fmt='(i4.4,i2.2,i2.2)') yr,mo,da
  write(unit=ftime3, fmt='(i2.2,i2.2)') hr,mn
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fversion,fmt='(i3.3)') version

  name = trim(pildasdir)//'/Y'//trim(ftime1)//'/M'//trim(fmo)//&
       '/PILDAS-1_forcng_F'//trim(fversion)//'V001_'//&
       trim(ftime2)//'_'//trim(ftime3)//'00.nc4'

end subroutine pildasfile
