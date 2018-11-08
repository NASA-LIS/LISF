!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: get_nldas1
! \label{get_nldas1}
!
!
! !REVISION HISTORY:
!  27 Apr 2000: Initial Specification
!  02 Feb 2004: Sujay Kumar ; Initial Version in LDT
!  06 Mar 2007: Kristi Arsenault; Implemented EDAS-Elevation Correction Var
!  15 Feb 2012: K. Arsenault; Accommodate GES DISC, NCEP filename conventions
!  14 Mar 2014: David Mocko: Updated filename options for config file
! 
! !INTERFACE:
subroutine get_nldas1(n,findex)
! !USES:
  use LDT_coreMod,        only : LDT_rc, LDT_domain
  use LDT_timeMgrMod,     only : LDT_tick
  use LDT_metforcingMod,  only : LDT_forc
  use LDT_logMod,         only : LDT_logunit, LDT_endrun
  use nldas1_forcingMod,  only : nldas1_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates hourly, NLDAS-1 forcing. 
!  At the beginning of a simulation, the code reads the most recent
!  past data (nearest hourly interval), and the nearest future data.
!  These two datasets are used to temporally interpolate the data to
!  the current model timestep.  The strategy for missing data is to
!  go backwards up to 10 days to get forcing at the same time of day.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    determines the NLDAS-1 data times
!  \item[get\_nldas1\_filename](\ref{get_nldas1_filename}) \newline
!    Puts together appropriate timestamped filename
!  \item[read\_nldas1](\ref{read_nldas1}) \newline
!    Interpolates NLDAS-1 data to LDT grid
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

!=== End Variable Definition =============================================
  try=-999
  
!====Assumption will be not to find or move any data
  nldas1_struc(n)%findtime1=0
  nldas1_struc(n)%findtime2=0
  movetime=0

!=== Determine Required NLDAS-1 Data Times (The previous hour and the future hour)
  yr1=LDT_rc%yr
  mo1=LDT_rc%mo
  da1=LDT_rc%da
  hr1=LDT_rc%hr
  mn1=LDT_rc%mn
  ss1=0
  ts1=0
  call LDT_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
 
  if(LDT_rc%ts.gt.3600) then 
     write(LDT_logunit,*) '[WARN] The model timestep is > forcing data timestep'
     write(LDT_logunit,*) ' LDT does not support this mode currently'
     write(LDT_logunit,*) 'Program stopping ...'
     call LDT_endrun()
  endif

  if(mod(nint(LDT_rc%ts),3600).eq.0) then 
     if(timenow.ge.nldas1_struc(n)%nldas1time2) then 
        yr1 = LDT_rc%yr
        mo1=LDT_rc%mo
        da1=LDT_rc%da
        hr1=LDT_rc%hr
        mn1=0
        ss1=0
        ts1=-60*60
        
        call LDT_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        
        yr2=LDT_rc%yr    !next hour
        mo2=LDT_rc%mo
        da2=LDT_rc%da
        hr2=LDT_rc%hr
        mn2=0
        ss2=0
        ts2=0
        call LDT_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        movetime = 1
        nldas1_struc(n)%findtime2 = 1
     endif
  else
     if(timenow.gt.nldas1_struc(n)%nldas1time2) then 
        yr1 = LDT_rc%yr
        mo1=LDT_rc%mo
        da1=LDT_rc%da
        hr1=LDT_rc%hr
        mn1=0
        ss1=0
        ts1=0
        
        call LDT_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

        yr2=LDT_rc%yr    !next hour
        mo2=LDT_rc%mo
        da2=LDT_rc%da
        hr2=LDT_rc%hr
        mn2=0
        ss2=0
        ts2=60*60
        call LDT_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

        movetime = 1
        nldas1_struc(n)%findtime2 = 1
     endif
  endif

  if(LDT_rc%tscount(n).eq.1 .or.LDT_rc%rstflag(n).eq.1  ) then    !beginning of the run	
     nldas1_struc(n)%findtime1=1
     nldas1_struc(n)%findtime2=1
     movetime=0
     LDT_rc%rstflag(n) = 0
  endif
  
  if(movetime.eq.1) then
     nldas1_struc(n)%nldas1time1=nldas1_struc(n)%nldas1time2
     do f=1,LDT_rc%met_nf(findex)
        do c=1,LDT_rc%ngrid(n)
           LDT_forc(n,findex)%metdata1(f,c)=LDT_forc(n,findex)%metdata2(f,c)
        enddo
     enddo
  endif    !end of movetime=1
  
  if(nldas1_struc(n)%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts1=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        call get_nldas1_filename(n,name,nldas1_struc(n)%nldas1dir,&
             yr1,mo1,da1,doy1,hr1)
        write(unit=LDT_logunit,fmt=*)'Getting NLDAS-1 file1.. ',trim(name)
        order = 1
        call read_nldas1(n,findex, order, name,ferror)
        if(ferror.eq.1) then 
           nldas1_struc(n)%nldas1time1=time1
        endif
        call LDT_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(*,*)'error: NLDAS-1 data gap exceeds 10 days on file 1'
           stop
        endif
     enddo
!=== end of data search
  endif   !end of LDT_rc%findtime=1	   	


  if(nldas1_struc(n)%findtime2.eq.1) then  
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        call get_nldas1_filename(n,name,nldas1_struc(n)%nldas1dir,&
             yr2,mo2,da2,doy2,hr2)
        write(unit=LDT_logunit,fmt=*)'Getting NLDAS-1 file2.. ',trim(name)
        order = 2
        call read_nldas1(n,findex, order, name,ferror)
        if(ferror.eq.1) then 
           nldas1_struc(n)%nldas1time2=time2
        endif
        call LDT_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(*,*)'error: NLDAS-1 data gap exceeds 10 days on file 2'
           stop
        endif
     enddo
     !=== end of data search
  endif   ! end of findtime2=1

end subroutine get_nldas1
