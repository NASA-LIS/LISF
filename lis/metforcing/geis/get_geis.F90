!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: get_geis
! \label{get_geis}
!
!
! !REVISION HISTORY:
! 15 Nov 2023: Sujay Kumar; Initial Version in LIS
! 
! !INTERFACE:
subroutine get_geis(n,findex)
! !USES:
  use LIS_coreMod, only        : LIS_rc, LIS_domain
  use LIS_timeMgrMod,     only : LIS_tick
  use LIS_metforcingMod,  only : LIS_forc
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use geis_forcingMod,  only : geis_struc
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates hourly, GLBOAL EIS forcing. 
!  At the beginning of a simulation, the code reads the most recent
!  past data (nearest hourly interval), and the nearest future data.
!  These two datasets are used to temporally interpolate the data to
!  the current model timestep.  The strategy for missing data is to
!  go backwards up to 10 days to get forcing at the same time of day.
!
!EOP

  integer :: c,f,ferror,try
  integer :: order
  real*8  :: time1,time2,timenow
  real*8  :: dtime1, dtime2
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  character(len=LIS_CONST_PATH_LEN) :: name
  real    :: gmt1,gmt2,ts1,ts2
  integer :: movetime     ! 1=move time 2 data into time 1  
  integer :: kk           ! Forecast member index

!=== End Variable Definition =============================================
  try=-999

  kk = 1
  
!====Assumption will be not to find or move any data
  geis_struc(n)%findtime1=0
  geis_struc(n)%findtime2=0
  movetime=0

!=== Determine Required NLDAS-2 Data Times (The previous hour and the future hour)
  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
 
  if(LIS_rc%ts.gt.3600) then 
     write(LIS_logunit,*) '[ERR] The model timestep is > forcing data timestep'
     write(LIS_logunit,*) '[ERR] LIS does not support this mode currently'
     write(LIS_logunit,*) '[ERR] Program stopping ...'
     call LIS_endrun()
  endif

  if(mod(nint(LIS_rc%ts),3600).eq.0) then 
     if(timenow.ge.geis_struc(n)%geistime2) then 
        yr1 = LIS_rc%yr
        mo1=LIS_rc%mo
        da1=LIS_rc%da
        hr1=LIS_rc%hr
        mn1=0
        ss1=0
        ts1=-60*60
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        
        yr2=LIS_rc%yr    !next hour
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=LIS_rc%hr
        mn2=0
        ss2=0
        ts2=0
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        movetime = 1
        geis_struc(n)%findtime2 = 1
     endif
  else
     if(timenow.ge.geis_struc(n)%geistime2) then 
        yr1 = LIS_rc%yr
        mo1=LIS_rc%mo
        da1=LIS_rc%da
        hr1=LIS_rc%hr
        mn1=0
        ss1=0
        ts1=0
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

        yr2=LIS_rc%yr    !next hour
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=LIS_rc%hr
        mn2=0
        ss2=0
        ts2=60*60
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

        movetime = 1
        geis_struc(n)%findtime2 = 1
     endif
  endif

  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then    !beginning of the run	
     geis_struc(n)%findtime1=1
     geis_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
  
  if(movetime.eq.1) then
     geis_struc(n)%geistime1=geis_struc(n)%geistime2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           geis_struc(n)%metdata1(:,f,c)=geis_struc(n)%metdata2(:,f,c)
        enddo
     enddo
  endif    !end of movetime=1
  
  if(geis_struc(n)%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts1=-60*60*24

     do 
        if ( ferror /= 0 ) exit
        try=try+1

        call get_geis_filename(n,kk,findex, &
             name,geis_struc(n)%geisdir,&
             yr1,mo1,da1,hr1)
        
        write(unit=LIS_logunit,fmt=*)'[INFO] getting file1.. ',trim(name)
        order = 1
        call read_geis(n,kk,findex,order,mo1,name,ferror)


        if(ferror.ge.1) geis_struc(n)%geistime1=time1
        call LIS_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(*,*)'error: GLOBAL EIS data gap exceeds 10 days on file 1'
           stop
        endif
     enddo
!=== end of data search
  endif   !end of LIS_rc%findtime=1	   	


  if(geis_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.

     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        
        call get_geis_filename(n,kk,findex,&
             name,geis_struc(n)%geisdir,&
             yr2,mo2,da2,hr2)

        write(unit=LIS_logunit,fmt=*)'[INFO] getting file2a.. ',trim(name)
        order = 2

        call read_geis(n,kk,findex,order,mo2,name,ferror)
        
        if(ferror.ge.1) then
           geis_struc(n)%geistime2=time2
        endif
        call LIS_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(*,*)'error: GLOBAL EIS data gap exceeds 10 days on file 2'
           stop
        endif
     enddo
     !=== end of data search
  endif   ! end of findtime2=1

end subroutine get_geis


!BOP
! !ROUTINE: get_geis_filename
! \label{get_geis_filename}
!
! !REVISION HISTORY:
!  15 Nov 2023: Sujay Kumar, initial specification
!
! !INTERFACE:
 subroutine get_geis_filename(n,kk,findex,filename,odir,yr,mo,da,hr)

! !USES:
   use LIS_coreMod
   use LIS_logMod

   implicit none
! !ARGUMENTS:
   integer, intent(in)        :: n
   integer, intent(in)        :: kk
   integer, intent(in)        :: findex
   character(len=*), intent(out) :: filename
   character(len=*), intent(in)   :: odir
   integer, intent(in)        :: yr,mo,da,hr

! !DESCRIPTION:
!   This subroutine puts together file name for
!   the Global EIS forcing data
!
!  The arguments are:
!  \begin{description}
!  \item[odir]
!    Name of the NLDAS-2 directory
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[filename]
!   name of the timestamped Global EIS data
!  \end{description}
!
!EOP

   character(len=6) :: fdir
   character(len=10) :: ftime

   !=== end variable definition =============================================

   !=== put together filename

   
   write(UNIT=fdir,fmt='(i4,i2.2)') yr, mo 
   write(UNIT=ftime,fmt='(i4,i2.2,i2.2,i2.2)') yr, mo, da, hr

   filename = trim(odir) // '/' // fdir // '/MERRA2_CERES_IMERG_' // &
        trim(ftime)//'00.d01.nc' 

end subroutine get_geis_filename
