!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_geos
! \label{get_geos}
!  
!
! !REVISION HISTORY:
!  27 Apr 2000: Brian Cosgrove; Added correction for use of old shortwave
!               data with opposite sign convention from recent shortwave data.
!               Added capability to use time averaged shortwave & longwave data
!               Altered times which are passed into ZTERP--used to be GMT1 
!               and GMT2, now they are LDAS%ETATIME1 and LDAS%ETATIME2
!   5 Nov 2001: Urszula Jambor; Reset tiny negative SW values to zero. 
!  10 Feb 2003: Sujay Kumar: Initial Specification in LIS.    
!   2 Dec 2008: Yudong Tian, added support for GEOS5 from 00Z2MAR2007
!
! !INTERFACE:
subroutine get_geos(n, findex)
! !USES:
  use LIS_coreMod,        only  : LIS_rc, LIS_domain
  use LIS_timeMgrMod,     only  : LIS_get_nstep, LIS_tick
  use LIS_logMod,         only  : LIS_logunit, LIS_endrun
  use geos_forcingMod,    only  : geos_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, 1 degree 
!  GEOS forcing. At the beginning of a simulation, the code 
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
!    call to advance or retract time
!  \item[geosfile](\ref{geosfile}) \newline
!    Puts together appropriate file name for 3 hour intervals
!  \item[read\_geos](\ref{read_geos}) \newline
!    call to read the GEOS data and perform spatial interpolation 
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer :: try, ferror
  integer :: c,t,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real*8  :: time1,time2,dumbtime1,dumbtime2
  real*8  :: timenow
  character*80 :: name
  character*40 :: elevfile, fpart1, fpart2
  real :: gmt1,gmt2,ts1,ts2
  integer :: movetime      ! 1=move time 2 data into time 1
  integer :: nstep
  real :: gridDesci(50)

     nstep = LIS_get_nstep(LIS_rc,n)

     geos_struc(n)%findtime1=0
     geos_struc(n)%findtime2=0
     movetime=0
!-------------------------------------------------------------------
! Determine Required GEOS Data Times 
! (The previous hour & the future hour)
!-------------------------------------------------------------------
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
     if(timenow.ge.geos_struc(n)%geostime2) then
        movetime=1
        geos_struc(n)%findtime2=1
     endif
     
     if ( nstep.eq.0 .or. nstep.eq.1 .or.LIS_rc%rstflag(n).eq.1 ) then 
        geos_struc(n)%findtime1=1
        geos_struc(n)%findtime2=1
        movetime=0
        LIS_rc%rstflag(n) = 0
     endif
     LIS_rc%shortflag=2            !Time averaged SW
     LIS_rc%longflag=2             !Time averaged LW 
     
     if(time1>geos_struc(n)%griduptime1.and.&
          time1 < geos_struc(n)%griduptime2.and.&
          geos_struc(n)%gridchange1) then 
        write(LIS_logunit,*) 'Time change..., Switching to GEOS4'
        geos_struc(n)%ncold = 288
        geos_struc(n)%mi = geos_struc(n)%ncold*geos_struc(n)%nrold
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
        gridDesci = 0
        gridDesci(1) = 0
        gridDesci(2) = geos_struc(n)%ncold
        gridDesci(3) = geos_struc(n)%nrold
        gridDesci(4) = -90.000
        gridDesci(7) = 90.000
        gridDesci(5) = -180.000
        gridDesci(6) = 128
        gridDesci(8) = 181.000
        gridDesci(9) = 1.254355
        gridDesci(10) = 1.000
        gridDesci(20) = 0
        if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
           call bilinear_interp_input(n,gridDesci,&
                geos_struc(n)%n111,geos_struc(n)%n121,&
                geos_struc(n)%n211,geos_struc(n)%n221,&
                geos_struc(n)%w111,geos_struc(n)%w121,&
                geos_struc(n)%w211,geos_struc(n)%w221)
        elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
           call bilinear_interp_input(n,gridDesci,&
                geos_struc(n)%n111,geos_struc(n)%n121,&
                geos_struc(n)%n211,geos_struc(n)%n221,&
                geos_struc(n)%w111,geos_struc(n)%w121,&
                geos_struc(n)%w211,geos_struc(n)%w221)
           call conserv_interp_input(n,gridDesci,&
                geos_struc(n)%n112,geos_struc(n)%n122,&
                geos_struc(n)%n212,geos_struc(n)%n222,&
                geos_struc(n)%w112,geos_struc(n)%w122,&
                geos_struc(n)%w212,geos_struc(n)%w222)
        endif
        
        geos_struc(n)%gridchange1 = .false.
        
        if ( trim(LIS_rc%met_ecor(findex)).ne."none") then 
           write(LIS_logunit,*) 'Elevation correction for GEOS forcing is not '
           write(LIS_logunit,*) 'Implemented... Program stopping ...'
           call LIS_endrun()
           LIS_rc%gridchange(n) = 1
           elevfile = LIS_rc%elevfile(n)
           c = index(elevfile,"geos3")
           fpart1 = elevfile(1:c+3)
           fpart2 = elevfile(c+5:40)
           LIS_rc%elevfile(n) = trim(fpart1) // "4" // trim(fpart2)
           write(LIS_logunit,*) 'Use newer elevation difference file: ', LIS_rc%elevfile(n)
           write(LIS_logunit,*) 'Transitioned from GEOS3 to GEOS4 grid dimensions.'
           call get_geos4_diff(n)
        endif
     elseif(time1 >=geos_struc(n)%griduptime2.and.&
          geos_struc(n)%gridchange2) then 

        write(LIS_logunit,*) 'Time change..., Switching to GEOS5'
        geos_struc(n)%ncold = 540
        geos_struc(n)%nrold = 361
        geos_struc(n)%mi = geos_struc(n)%ncold*geos_struc(n)%nrold
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
        gridDesci = 0
        gridDesci(1) = 0
        gridDesci(2) = geos_struc(n)%ncold
        gridDesci(3) = geos_struc(n)%nrold
        gridDesci(4) = -90.000
        gridDesci(5) = -180.000
        gridDesci(6) = 128
        gridDesci(7) = 90.000
        gridDesci(8) = 179.333333 
        gridDesci(9) = 0.66666666667 
        gridDesci(10) = 0.5
        gridDesci(20) = 0
        if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
           call bilinear_interp_input(n,gridDesci,&
                geos_struc(n)%n111,geos_struc(n)%n121,&
                geos_struc(n)%n211,geos_struc(n)%n221,&
                geos_struc(n)%w111,geos_struc(n)%w121,&
                geos_struc(n)%w211,geos_struc(n)%w221)
        elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
           call bilinear_interp_input(n,gridDesci,&
                geos_struc(n)%n111,geos_struc(n)%n121,&
                geos_struc(n)%n211,geos_struc(n)%n221,&
                geos_struc(n)%w111,geos_struc(n)%w121,&
                geos_struc(n)%w211,geos_struc(n)%w221)
           call conserv_interp_input(n,gridDesci,&
                geos_struc(n)%n112,geos_struc(n)%n122,&
                geos_struc(n)%n212,geos_struc(n)%n222,&
                geos_struc(n)%w112,geos_struc(n)%w122,&
                geos_struc(n)%w212,geos_struc(n)%w222)
        endif

        geos_struc(n)%gridchange2 = .false.

     endif
!-------------------------------------------------------------------
! Establish geostime1
!-------------------------------------------------------------------
     if (geos_struc(n)%findtime1==1) then  
        order=1   
        ferror = 0
        try = 0
        ts1 = -24*60*60
        do
           if ( ferror /= 0 ) then
              exit
           end if
           try = try+1
           call geosfile(name,geos_struc(n)%geosdir,&
                yr1,mo1,da1,hr1,geos_struc(n)%ncold)
           call read_geos(order,n, findex,&
                name,LIS_rc%tscount(n),ferror)
           if ( ferror == 1 ) then 
!-------------------------------------------------------------------
! successfully retrieved forcing data
!-------------------------------------------------------------------
              geos_struc(n)%geostime1=time1
           else  
!-------------------------------------------------------------------
! ferror still=0, so roll back one day
!-------------------------------------------------------------------
              call LIS_tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
           end if
           if ( try > ndays ) then 
              print *, 'ERROR: GEOS data gap exceeds 10 days on file 1'
              call LIS_endrun
           end if
        end do
     endif
     if(movetime.eq.1) then
        geos_struc(n)%geostime1=geos_struc(n)%geostime2
        geos_struc(n)%findtime2=1 
        do f=1,LIS_rc%met_nf(findex)
           do t=1,LIS_rc%ngrid(n)
              geos_struc(n)%metdata1(f,t)=geos_struc(n)%metdata2(f,t)
           enddo
        enddo
     endif
     if(geos_struc(n)%findtime2.eq.1) then 
        order=2   
        ferror = 0
        try = 0
        ts2 = -24*60*60
        do
           if ( ferror /= 0 ) exit
           try = try+1
           call geosfile(name,geos_struc(n)%geosdir,yr2,mo2,da2,hr2,geos_struc(n)%ncold)
           call read_geos(order,n,findex, &
                name,LIS_rc%tscount(n),ferror)
           if ( ferror == 1 ) then 
!-------------------------------------------------------------------
! successfully retrieved forcing data
!-------------------------------------------------------------------
              geos_struc(n)%geostime2=time2
           else  
!-------------------------------------------------------------------
! ferror still=0, so roll back one day
!-------------------------------------------------------------------
              call LIS_tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
           end if
           if ( try > ndays ) then 
              print *, 'ERROR: GEOS data gap exceeds 10 days on file 2'
              call LIS_endrun
           end if
        end do
     endif
     
84   format('now',i4,4i3,2x,'pvt ',a22,' nxt ',a22)
!  if ((LIS_rc%gridchange==1).and.(geos_struc(n)%ncold==288)) then
!     LIS_rc%gridchange=0
!  endif
  return 

end subroutine get_geos


!BOP
! !ROUTINE: geosfile
! \label{geosfile}
!
! !INTERFACE:
subroutine geosfile(name,geosdir,yr,mo,da,hr,ncold)
  
  implicit none
! !ARGUMENTS: 
  character(len=*), intent(in)  :: geosdir
  integer, intent(in)           :: yr,mo,da,hr,ncold
  character(len=*), intent(out) :: name
! !DESCRIPTION:
!   This subroutine puts together GEOS file name for 
!   3 hour file intervals.
! 
!  The arguments are:
!  \begin{description}
!  \item[geosdir]
!    Name of the GEOS directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!  \item[ncold]
!   Number of columns (used as to determine the GEOS resolution)
!  \item[name]
!   name of the timestamped GEOS file
!  \end{description}
!
!EOP
  integer uyr,umo,uda,uhr,i,c,ii,jj
  character(len=2) :: initcode
  character*1 fbase(80),fsubs(80)
  character*1 ftime(10),fdir(8)
  
  character(LEN=100) :: temp
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
  
  write(UNIT=temp,FMT='(A40)') geosdir 
  read(UNIT=temp,FMT='(80A1)') (fbase(i),i=1,80)
  
  write(UNIT=temp,FMT='(a1,i4,i2,a1)') '/',uyr,umo,'/'
  read(UNIT=temp,FMT='(8A1)') fdir
  do i=1,8
     if(fdir(i).eq.(' ')) fdir(i)='0'
  enddo
  
  write(UNIT=temp,FMT='(i4,i2,i2,a2)') uyr,umo,uda,initcode
  read(UNIT=temp,FMT='(10A1)') ftime
  do i=1,10
     if(ftime(i).eq.(' ')) ftime(i)='0'
  enddo
  
  if(ncold==360) then 
     write(UNIT=temp,FMT='(A8)') '.GEOS323'
     read(UNIT=temp,FMT='(80A1)') (fsubs(i),i=1,8)
  elseif(ncold==288) then 
     write(UNIT=temp,FMT='(A6)') '.GEOS4'
     read(UNIT=temp,FMT='(80A1)') (fsubs(i),i=1,6)
  else if(ncold==540) then 
     write(UNIT=temp,FMT='(A6)') '.GEOS5'
     read(UNIT=temp,FMT='(80A1)') (fsubs(i),i=1,6)
  endif
  c=0
  do i=1,80
     if(fbase(i).eq.(' ').and.c.eq.0) c=i-1 
  enddo
  
  if (ncold==360) then       
     write(UNIT=temp,FMT='(80a1)') (fbase(i),i=1,c),(fdir(i),i=1,8), &
          (ftime(i),i=1,10),(fsubs(i),i=1,8)
  else !covers both GEOS4 and GEOS5
     write(UNIT=temp,FMT='(80a1)') (fbase(i),i=1,c),(fdir(i),i=1,8), &
          (ftime(i),i=1,10),(fsubs(i),i=1,6)
  endif
  
  read(UNIT=temp, FMT='(a80)') name
  return

end subroutine geosfile

subroutine get_geos4_diff(n)
   use LIS_coreMod, only : LIS_rc, LIS_domain
   use LIS_logMod,  only : LIS_logunit, LIS_verify

   implicit none
   integer, intent(in) :: n
   real, allocatable, dimension(:,:) :: elev
   integer :: ierr, i, j

   print*, 'get_geos4_diff needs to be fixed'
   stop
#if 0 
   allocate(elev(LIS_rc%lnc(n),LIS_rc%lnr(n)), stat=ierr)
   call LIS_verify(ierr,'Error allocating elev diff.')

   elev = 0.0

   call readelev(n, LIS_rc%elevsrc, elev)
        
   write(LIS_logunit,*) 'MSG: get_geos4_diff -- done reading elevation '// &
                    'difference file'

!   do i=1,LIS_rc%ntiles
!      tile(i)%elev = elev(tile(i)%col, tile(i)%row-lis_tnroffset)
!   enddo
   do j = 1, LIS_rc%lnr(n)
      do i = 1, LIS_rc%lnc(n)
         if ( elev(i,j) == -9999.0 ) then
            elev(i,j) = 0.0
         endif
         LIS_domain(n)%grid(LIS_domain(n)%gindex(i,j))%elev=elev(i,j)
      enddo
   enddo

   deallocate(elev,stat=ierr)
   call LIS_verify(ierr,'Error deallocating elev diff.')
#endif
end subroutine get_geos4_diff
