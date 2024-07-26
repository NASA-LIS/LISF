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
! !ROUTINE: read_Loobos
! \label{read_Loobos}
! 
! !REVISION HISTORY:
! 06 Oct 2010: David Mocko, Updated for Bondville test case
! 20 Feb 2018: Shugong Wang, Updated for Loobos test case 
! 
! !INTERFACE:
subroutine read_Loobos(n,ftn,findex,order,itime)
! !USES:
  use LIS_logMod, only            : LIS_logunit,LIS_endrun
  use LIS_coreMod, only           : LIS_rc,LIS_domain
  use LIS_metforcingMod,     only : LIS_forc
  use LIS_timeMgrMod, only        : LIS_date2time,LIS_tick
  use Loobos_forcingMod, only : Loobos_struc
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: ftn
  integer, intent(in) :: itime
  integer, intent(in) :: findex
  integer, intent(in) :: order
!
! !DESCRIPTION:
!  For the given time, reads parameters from the correct Loobos
!  station data (ASCII), transforms into LIS forcing parameters, and
!  interpolates to the LIS domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[ftn]
!    unit number for the Loobos station data
!  \item[itime]
!    index for reading current time or reading a future time
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[normalize\_stnwts](\ref{normalize_stnwts}) \newline
!    renormalizes the station weights accounting for
!    missing data
!  \item[interp\_stndata](\ref{interp_stndata}) \newline
!    spatially interpolates the station data onto the LIS grid.
!  \end{description}
!EOP

  integer :: i,c,r,f,count1, k 
  real    :: rain(Loobos_struc(n)%nstns)
  real    :: snow(Loobos_struc(n)%nstns)
  !real    :: pcp(Loobos_struc(n)%nstns)
  real    :: psurf(Loobos_struc(n)%nstns)
  real    :: tair(Loobos_struc(n)%nstns)
  real    :: qair(Loobos_struc(n)%nstns)
  real    :: swdown(Loobos_struc(n)%nstns)
  real    :: lwdown(Loobos_struc(n)%nstns)
  real    :: u(Loobos_struc(n)%nstns)
  real    :: v(Loobos_struc(n)%nstns)
  real    :: varfield(9,LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real    :: varfield1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real*8  :: listime,loobos_time
  real    :: lisgmt,loobos_gmt
  integer :: lisdoy,loobos_doy
  integer :: loobos_yr,loobos_mon,loobos_day,loobos_hr,loobos_min,loobos_sec
  real    :: loobos_tick
  logical :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: Loobos_filename
  character(len=500) :: line

  !      write(LIS_logunit,*) 'starting read_Loobos'
  do i = 1,Loobos_struc(n)%nstns
     ! Generate the Loobos filename and see if it exists
     Loobos_filename = trim(Loobos_struc(n)%Loobosfile)
     write(LIS_logunit,*) 'Reading Loobos file: ',             &
          trim(Loobos_filename)
     inquire(file=Loobos_filename,exist=file_exists)
     if (file_exists) then
        !            write(LIS_logunit,*) 'File is open!!'
        open(ftn,file=Loobos_filename,form='formatted',        &
             status='old')
        ! segment to skip the header
        do k=1,5
           read(ftn,'(a)') line
        enddo
        
        ! Actual data reading section
        do
           read(ftn,'(a)') line
           read(line,40) loobos_yr,loobos_mon,loobos_day,loobos_hr,loobos_min,         &
                swdown(i),lwdown(i), rain(i), snow(i), tair(i),u(i),psurf(i),qair(i) 
           v(i)   = 0
           !pcp(i) = rain(i) + snow(i) 
           loobos_sec = 0
           !loobos_tick = 21600
           loobos_tick = 0 
           call LIS_date2time(loobos_time,loobos_doy,loobos_gmt,loobos_yr,loobos_mon,loobos_day,loobos_hr,loobos_min,loobos_sec)
           call LIS_date2time(listime    ,    lisdoy,    lisgmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da,LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
           ! Convert local solar time of Loobos forcing data into GMT time
           ! as in LIS; add 1 hours to Loobos to get GMT.  When Loobos
           ! +1 hours is greater than or equal to the LIS time (depending if
           ! we are interpolating between half-hourly data), then stop reading
           ! the Loobos data and use that line as forcing at the LIS time.
           call LIS_tick(loobos_time,loobos_doy,loobos_gmt,loobos_yr,loobos_mon,       &
                loobos_day,loobos_hr,loobos_min,loobos_sec,loobos_tick)
           if ((loobos_time.ge.listime).and.(itime.eq.1)) exit
           if ((loobos_time.gt.listime).and.(itime.eq.2)) exit
        enddo
        close(ftn)
     else
        write(LIS_logunit,*) 'Filename not found; stopping program'
        call LIS_endrun
     endif
  enddo

!40 format(i4,4i3,8f17.10)
40 format(i4,4i3, f8.1, f7.1, e14.2, e14.2, f10.3, f10.3, f13.1, e12.2)

  call normalize_stnwts(psurf,Loobos_struc(n)%nstns,           &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%undef,  &
       Loobos_struc(n)%stnwt)
  call interp_stndata(Loobos_struc(n)%stnwt,                   &
       Loobos_struc(n)%undef,psurf,varfield(1,:),          &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%nstns)

  call normalize_stnwts(tair,Loobos_struc(n)%nstns,            &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%undef,  &
       Loobos_struc(n)%stnwt)
  call interp_stndata(Loobos_struc(n)%stnwt,                   &
       Loobos_struc(n)%undef,tair,varfield(2,:),           &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%nstns)

  call normalize_stnwts(qair,Loobos_struc(n)%nstns,            &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%undef,  &
       Loobos_struc(n)%stnwt)
  call interp_stndata(Loobos_struc(n)%stnwt,                   &
       Loobos_struc(n)%undef,qair,varfield(3,:),           &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%nstns)

  call normalize_stnwts(u,Loobos_struc(n)%nstns,               &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%undef,  &
       Loobos_struc(n)%stnwt)
  call interp_stndata(Loobos_struc(n)%stnwt,                   &
       Loobos_struc(n)%undef,u,varfield(4,:),              &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%nstns)

  call normalize_stnwts(v,Loobos_struc(n)%nstns,               &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%undef,  &
       Loobos_struc(n)%stnwt)
  call interp_stndata(Loobos_struc(n)%stnwt,                   &
       Loobos_struc(n)%undef,v,varfield(5,:),              &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%nstns)

  call normalize_stnwts(swdown,Loobos_struc(n)%nstns,          &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%undef,  &
       Loobos_struc(n)%stnwt)
  call interp_stndata(Loobos_struc(n)%stnwt,                   &
       Loobos_struc(n)%undef,swdown,varfield(6,:),         &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%nstns)

  call normalize_stnwts(lwdown,Loobos_struc(n)%nstns,          &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%undef,  &
       Loobos_struc(n)%stnwt)
  call interp_stndata(Loobos_struc(n)%stnwt,                   &
       Loobos_struc(n)%undef,lwdown,varfield(7,:),         &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%nstns)

  call normalize_stnwts(rain,Loobos_struc(n)%nstns,             &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%undef,  &
       Loobos_struc(n)%stnwt)
  call interp_stndata(Loobos_struc(n)%stnwt,                   &
       Loobos_struc(n)%undef,rain,varfield(8,:),            &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%nstns)
  
  call normalize_stnwts(snow,Loobos_struc(n)%nstns,             &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%undef,  &
       Loobos_struc(n)%stnwt)
  call interp_stndata(Loobos_struc(n)%stnwt,                   &
       Loobos_struc(n)%undef,snow,varfield(9,:),            &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%nstns)

  do f = 1,8
     count1 = 0
     do r = 1,LIS_rc%lnr(n)
        do c = 1,LIS_rc%lnc(n)
           varfield1(c,r) = varfield(f,c+count1)
        enddo
        count1 = count1 + LIS_rc%lnc(n)
     enddo

     do r = 1,LIS_rc%lnr(n)
        do c = 1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r).ne.-1) then
              if(order.eq.1) then 
                 Loobos_struc(n)%metdata1(f,LIS_domain(n)%gindex(c,r)) = varfield1(c,r)
              elseif(order.eq.2) then 
                 Loobos_struc(n)%metdata2(f,LIS_domain(n)%gindex(c,r)) = varfield1(c,r)
              endif
           endif
        enddo
     enddo
  enddo

end subroutine read_Loobos

