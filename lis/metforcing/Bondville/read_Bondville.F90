!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_Bondville
! \label{read_Bondville}
! 
! !REVISION HISTORY:
! 06 Oct 2010: David Mocko, Updated for Bondville test case
! 26 Oct 2018: David Mocko, Updated for Noah-MP-4.0.1 HRLDAS test case
!
! !INTERFACE:
subroutine read_Bondville(n,ftn,findex,order,itime)
! !USES:
  use LIS_logMod, only            : LIS_logunit,LIS_endrun
  use LIS_coreMod, only           : LIS_rc,LIS_domain
  use LIS_metforcingMod,     only : LIS_forc
  use LIS_timeMgrMod, only        : LIS_date2time,LIS_tick
  use Bondville_forcingMod, only : Bondville_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: ftn
  integer, intent(in) :: itime
  integer, intent(in) :: findex
  integer, intent(in) :: order
!
! !DESCRIPTION:
!  For the given time, reads parameters from the correct Bondville
!  station data (ASCII), transforms into LIS forcing parameters, and
!  interpolates to the LIS domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[ftn]
!    unit number for the Bondville station data
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

  integer :: i,c,r,f,count1
  real    :: pcp(Bondville_struc(n)%nstns)
  real    :: psurf(Bondville_struc(n)%nstns)
  real    :: tair(Bondville_struc(n)%nstns)
  real    :: qair(Bondville_struc(n)%nstns)
  real    :: swdown(Bondville_struc(n)%nstns)
  real    :: lwdown(Bondville_struc(n)%nstns)
  real    :: u(Bondville_struc(n)%nstns)
  real    :: v(Bondville_struc(n)%nstns)
  real    :: varfield(8,LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real    :: varfield1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real*8  :: listime,bontime
  real    :: lisgmt,bongmt
  integer :: lisdoy,bondoy
  integer :: bonyr,bonmon,bonday,bonhr,bonmin,bonsec
  real    :: bontick
  logical :: file_exists
  character*80       :: Bondville_filename
  character(len=500) :: line

  ! write(LIS_logunit,*) 'starting read_Bondville'
  do i = 1,Bondville_struc(n)%nstns
     ! Generate the Bondville filename and see if it exists
     Bondville_filename = trim(Bondville_struc(n)%Bondvillefile)
     write(LIS_logunit,*) 'Reading Bondville file: ',             &
          trim(Bondville_filename)
     inquire(file=Bondville_filename,exist=file_exists)
     if (file_exists) then
        ! write(LIS_logunit,*) 'File is open!!'
        open(ftn,file=Bondville_filename,form='formatted',        &
             status='old')

! If MP=0, use the older "bondville.dat" forcing file.
!    This option is to support older testcases.
! If MP.ne.0, use the new Noah-MP-4.0.1 HRLDAS "bondville.dat" file.
!
! Segment to skip the header
        if (Bondville_struc(n)%MP.eq.0) then
           do
              read(ftn,'(a)') line
              if (line(1:9).eq.'<Forcing>') exit
           enddo
        else
           do
              read(ftn,'(a)') line
              if (line(1:4).eq.'yyyy') exit
           enddo
           read(ftn,'(a)') line
        endif
! Actual data reading section
        do
           read(ftn,'(a)') line
           if (Bondville_struc(n)%MP.eq.0) then
              read(line,40) bonyr,bonmon,bonday,bonhr,bonmin,         &
                   u(i),v(i),tair(i),qair(i),psurf(i),                &
                   swdown(i),lwdown(i),pcp(i)
           else
              read(line,50) bonyr,bonmon,bonday,bonhr,bonmin,         &
                            u(i),tair(i),qair(i),psurf(i),            &
                            swdown(i),lwdown(i),pcp(i)
! Convert data into expected units
              v(i) = 0.0
              tair(i) = tair(i) + 273.15
              pcp(i) = pcp(i) / 1800.0 * 25.4
           endif
           ! write(LIS_logunit,*) bonyr,bonmon,bonday,bonhr,bonmin,      &
           !                      u(i),v(i),tair(i),qair(i),psurf(i),    &
           !                      swdown(i),lwdown(i),pcp(i)
           bonsec = 0
           bontick = 21600
           call LIS_date2time(bontime,bondoy,bongmt,bonyr,bonmon,  &
                bonday,bonhr,bonmin,bonsec)
           call LIS_date2time(listime,lisdoy,lisgmt,LIS_rc%yr,     &
                LIS_rc%mo,LIS_rc%da,LIS_rc%hr,       &
                LIS_rc%mn,LIS_rc%ss)
           ! Convert local solar time of Bondville forcing data into GMT time
           ! as in LIS; add 6 hours to Bondville to get GMT.  When Bondville
           ! +6 hours is greater than or equal to the LIS time (depending if
           ! we are interpolating between half-hourly data), then stop reading
           ! the Bondville data and use that line as forcing at the LIS time.
           call LIS_tick(bontime,bondoy,bongmt,bonyr,bonmon,       &
                         bonday,bonhr,bonmin,bonsec,bontick)
           ! write(LIS_logunit,*) 'read_time: ',bontime,bonyr,bonmon,&
           !                                    bonday,bonhr,bonmin
           ! write(LIS_logunit,*) 'LIS_time: ',listime,LIS_rc%yr,    &
           !                LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn
           if ((bontime.ge.listime).and.(itime.eq.1)) exit
           if ((bontime.gt.listime).and.(itime.eq.2)) exit
        enddo
        close(ftn)
     else
        write(LIS_logunit,*) 'Filename not found; stopping program'
        call LIS_endrun
     endif
  enddo

40 format(i4,4i3,8f17.10)
50 format(i4,4i3,f6.2,2f6.1,2f6.0,f5.0,f6.2)

  call normalize_stnwts(psurf,Bondville_struc(n)%nstns,           &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%undef,  &
       Bondville_struc(n)%stnwt)
  call interp_stndata(Bondville_struc(n)%stnwt,                   &
       Bondville_struc(n)%undef,psurf,varfield(1,:),          &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%nstns)

  call normalize_stnwts(tair,Bondville_struc(n)%nstns,            &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%undef,  &
       Bondville_struc(n)%stnwt)
  call interp_stndata(Bondville_struc(n)%stnwt,                   &
       Bondville_struc(n)%undef,tair,varfield(2,:),           &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%nstns)

  call normalize_stnwts(qair,Bondville_struc(n)%nstns,            &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%undef,  &
       Bondville_struc(n)%stnwt)
  call interp_stndata(Bondville_struc(n)%stnwt,                   &
       Bondville_struc(n)%undef,qair,varfield(3,:),           &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%nstns)

  call normalize_stnwts(u,Bondville_struc(n)%nstns,               &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%undef,  &
       Bondville_struc(n)%stnwt)
  call interp_stndata(Bondville_struc(n)%stnwt,                   &
       Bondville_struc(n)%undef,u,varfield(4,:),              &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%nstns)

  call normalize_stnwts(v,Bondville_struc(n)%nstns,               &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%undef,  &
       Bondville_struc(n)%stnwt)
  call interp_stndata(Bondville_struc(n)%stnwt,                   &
       Bondville_struc(n)%undef,v,varfield(5,:),              &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%nstns)

  call normalize_stnwts(swdown,Bondville_struc(n)%nstns,          &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%undef,  &
       Bondville_struc(n)%stnwt)
  call interp_stndata(Bondville_struc(n)%stnwt,                   &
       Bondville_struc(n)%undef,swdown,varfield(6,:),         &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%nstns)

  call normalize_stnwts(lwdown,Bondville_struc(n)%nstns,          &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%undef,  &
       Bondville_struc(n)%stnwt)
  call interp_stndata(Bondville_struc(n)%stnwt,                   &
       Bondville_struc(n)%undef,lwdown,varfield(7,:),         &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%nstns)

  call normalize_stnwts(pcp,Bondville_struc(n)%nstns,             &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%undef,  &
       Bondville_struc(n)%stnwt)
  call interp_stndata(Bondville_struc(n)%stnwt,                   &
       Bondville_struc(n)%undef,pcp,varfield(8,:),            &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%nstns)

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
                 Bondville_struc(n)%metdata1(f,LIS_domain(n)%gindex(c,r)) = varfield1(c,r)
              elseif(order.eq.2) then 
                 Bondville_struc(n)%metdata2(f,LIS_domain(n)%gindex(c,r)) = varfield1(c,r)
              endif
           endif
        enddo
     enddo
  enddo

end subroutine read_Bondville

