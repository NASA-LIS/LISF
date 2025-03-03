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
! !ROUTINE: get_geos5fcst
! \label{get_geos5fcst}
!
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
! 
! !INTERFACE:
subroutine get_geos5fcst(n,findex)
! !USES:
  use LIS_coreMod,          only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,       only : LIS_tick
  use LIS_metforcingMod,    only : LIS_forc
  use LIS_logMod,           only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,     only : LIS_CONST_PATH_LEN
  use geos5fcst_forcingMod, only : geos5fcst_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates hourly, GEOS5 forecast forcing. 
!  At the beginning of a simulation, the code reads the most recent
!  past data (nearest hourly interval), and the nearest future data.
!  These two datasets are used to temporally interpolate the data to
!  the current model timestep.   The strategy for missing data is to
!  go backwards up to 10 days to get forcing at the same time of day.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the GEOS5 data times
!  \item[geos5fcstfile](\ref{geos5fcstfile}) \newline
!    Puts together appropriate timestamped GEOS5 data
!  \item[read\_geos5fcst](\ref{read_geos5fcst}) \newline
!    Reads and Interpolates GEOS5 forecast data to the LIS grid
!  \end{description}
!
!EOP

  integer       :: c,f,m,ferror,try
  integer       :: order
  real*8        :: time1,time2,timenow
  real*8        :: dtime1, dtime2
  integer       :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer       :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  character(len=LIS_CONST_PATH_LEN) :: name
  real          :: gmt1,gmt2,ts1,ts2
  integer       :: movetime     ! 1=move time 2 data into time 1

  try=-999

  geos5fcst_struc(n)%findtime1=0
  geos5fcst_struc(n)%findtime2=0
  movetime=0

  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=LIS_rc%ss
  ts1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
 
  if(LIS_rc%ts.gt.3600) then 
     write(LIS_logunit,*) '[WARN] The model timestep is > forcing data timestep'
     write(LIS_logunit,*) '[WARN] LIS does not support this mode currently'
     write(LIS_logunit,*) '[WARN] Program stopping ...'
     call LIS_endrun()
  endif

  if(timenow.gt.geos5fcst_struc(n)%fcsttime2) then 
     yr1=LIS_rc%yr
     mo1=LIS_rc%mo
     da1=LIS_rc%da
     hr1=LIS_rc%hr
     mn1=30 
     ss1=0
     ts1=-60*60
     
     call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
     
     yr2=LIS_rc%yr    !next hour
     mo2=LIS_rc%mo
     da2=LIS_rc%da
     hr2=LIS_rc%hr
     mn2=30
     ss2=0
     ts2=0
     call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
     movetime = 1
     geos5fcst_struc(n)%findtime2 = 1
  elseif(timenow.eq.geos5fcst_struc(n)%fcsttime2) then 
     yr1 = LIS_rc%yr
     mo1=LIS_rc%mo
     da1=LIS_rc%da
     hr1=LIS_rc%hr
     mn1=30
     ss1=0
     ts1=0
     
     call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
     
     yr2=LIS_rc%yr    !next hour
     mo2=LIS_rc%mo
     da2=LIS_rc%da
     hr2=LIS_rc%hr
     mn2=30
     ss2=0
     ts2=60*60
     call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
     
     movetime = 1
     geos5fcst_struc(n)%findtime2 = 1
  endif
    
  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then    !beginning of the run	
     geos5fcst_struc(n)%findtime1=1
     geos5fcst_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
  if(movetime.eq.1) then
     geos5fcst_struc(n)%fcsttime1=&
          geos5fcst_struc(n)%fcsttime2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           geos5fcst_struc(n)%metdata1(f,:,c) = geos5fcst_struc(n)%metdata2(f,:,c)
        enddo
     enddo
  endif    !end of movetime=1

  
  if(geos5fcst_struc(n)%findtime1.eq.1) then

     ferror=0
     try=0  
     ts1=-60*60*24

     do 
        if ( ferror /= 0 ) exit
        try=try+1
        
        do m=1,geos5fcst_struc(n)%max_ens_members
           call geos5fcstfile(name,geos5fcst_struc(n)%geos5fcstdir,&
                yr1,mo1,da1,hr1,mn1,m)
        
           write(LIS_logunit,*)'[INFO] getting file1.. ',trim(name)
           order = 1
           call read_geos5fcst(n,m,findex,order,name,ferror)
           if(ferror.ge.1) &
                geos5fcst_struc(n)%fcsttime1=time1
        enddo
        call LIS_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(*,*)'error: GEOS5fcst data gap exceeds 10 days on file 1'
           stop
        endif
     enddo
  endif   !end of LIS_rc%findtime=1	   	


  if(geos5fcst_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        
        do m=1,geos5fcst_struc(n)%max_ens_members
           call geos5fcstfile(name,geos5fcst_struc(n)%geos5fcstdir,&
                yr2,mo2,da2,hr2,mn2,m)
        
           write(LIS_logunit,*)'[INFO] getting file2.. ',trim(name)
           order = 2
           call read_geos5fcst(n,m,findex,order,name,ferror)
           if(ferror.ge.1) then
              geos5fcst_struc(n)%fcsttime2=time2
           endif
        enddo
        call LIS_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(*,*)'error: GEOS5fcst data gap exceeds 10 days on file 2'
           stop
        endif
     enddo
  endif   

end subroutine get_geos5fcst

!BOP
! !ROUTINE: geos5fcstfile
! \label{geos5fcstfile}
!
! !INTERFACE:
subroutine geos5fcstfile(name,geos5fcstdir,yr,mo,da,hr,mn,ens_id)
  
  implicit none
! !ARGUMENTS: 
  character(len=*), intent(out) :: name
  character(len=*), intent(in)  :: geos5fcstdir
  integer, intent(in)           :: yr,mo,da,hr,mn
  integer                       :: ens_id
! !DESCRIPTION:
!  This subroutine puts together GEOS5FCST file name for 
!  1 hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[name]
!   name of the timestamped GEOS5 forecast file
!  \item[geos5fcstdir]
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
!  \item[ens\_id]
!   id of the ensemble member
!  \end{description}
!
!EOP
  character*3 :: month
  character*6 :: ftime1
  character*8 :: ftime2
  character*4 :: ftime3
  character*2 :: ftime4
  character*2 :: fens

  write(unit=ftime1, fmt='(i4.4,i2.2)') yr,mo
  write(unit=ftime2, fmt='(i4.4,i2.2,i2.2)') yr,mo,da
  write(unit=ftime3, fmt='(i2.2,i2.2)') hr,mn  
  if(mo.eq.1) then 
     month = 'jan'
     ftime4 = '01'
  elseif(mo.eq.2) then
     month = 'jan'
     ftime4 = '31'
  elseif(mo.eq.3) then
     month = 'mar'
     ftime4 = '02'    
  elseif(mo.eq.4) then
     month = 'apr'
     ftime4 = '01'    
  elseif(mo.eq.5) then
     month = 'may'
     ftime4 = '01'    
  elseif(mo.eq.6) then
     month = 'may'
     ftime4 = '31'    
  elseif(mo.eq.7) then
     month = 'jun'
     ftime4 = '30'    
  elseif(mo.eq.8) then
     month = 'jul'
     ftime4 = '30'    
  elseif(mo.eq.9) then
     month = 'aug'
     ftime4 = '29'    
  elseif(mo.eq.10) then
     month = 'oct'
     ftime4 = '03'    
  elseif(mo.eq.11) then
     month = 'nov'
     ftime4 = '02'    
  elseif(mo.eq.12) then
     month = 'dec'
     ftime4 = '02'    
  endif

  write(unit=fens, fmt='(i2)') ens_id

  name = trim(geos5fcstdir)//'/'//trim(ftime1)//'/ens'//trim(adjustl(fens))//'/'&
       //trim(month)//trim(ftime4)//'.geosgcm_surfh.'//trim(ftime2)//'_'//&
       trim(ftime3)//'z.nc4'

end subroutine geos5fcstfile

