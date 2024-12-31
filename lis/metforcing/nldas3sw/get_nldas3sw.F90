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
! !ROUTINE: get_nldas3sw
! \label{get_nldas3sw}
!
! !REVISION HISTORY:
! 27 Dec 2024: David Mocko, Initial Specification
!                           (derived from get_nldas20.F90)
!
! !INTERFACE:
subroutine get_nldas3sw(n,findex)
! !USES:
  use LIS_coreMod, only         : LIS_rc
  use LIS_timeMgrMod,     only  : LIS_tick
  use LIS_metforcingMod,  only  : LIS_forc
  use LIS_logMod,         only  : LIS_logunit,LIS_endrun
  use LIS_constantsMod,   only  : LIS_CONST_PATH_LEN
  use nldas3sw_forcingMod, only : nldas3sw_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates hourly 4-km CERES SWdown forcing.
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
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the forcing data times
!  \item[ceres\_nldas3swfile](\ref{ceres_nldas3swfile}) \newline
!    Puts together appropriate timestamped CERES filename
!  \item[read\_nldas3sw](\ref{read_nldas3sw}) \newline
!    Reads and Interpolates CERES SWdown data to LIS grid
!  \end{description}
!
!EOP
  integer :: c,f,ferrora,ferror,try
  integer :: order
  real*8  :: time1,time2,timenow
  real*8  :: dtime1, dtime2
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  character(len=LIS_CONST_PATH_LEN) :: name_a
  real    :: gmt1,gmt2,ts1,ts2
  integer :: movetime       ! 1=move time 2 data into time 1

  external :: ceres_nldas3swfile
  external :: read_nldas3sw

!=== End Variable Definition ===========================================
  try = -999

!=== Assumption will be not to find or move any data
  nldas3sw_struc(n)%findtime1 = 0
  nldas3sw_struc(n)%findtime2 = 0
  movetime = 0

!=== Determine Required CERES Data Times (The previous hour and the future hour)
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  if (LIS_rc%ts.gt.3600) then
     write(LIS_logunit,*)                                              &
          "[ERR] The model timestep is > forcing data timestep"
     write(LIS_logunit,*)                                              &
          "[ERR] LIS does not support this mode currently"
     write(LIS_logunit,*) "[ERR] Program stopping ..."
     call LIS_endrun()
  endif

  if (mod(nint(LIS_rc%ts),3600).eq.0) then
     if (timenow.ge.nldas3sw_struc(n)%nldas3swtime2) then
        yr1 = LIS_rc%yr
        mo1 = LIS_rc%mo
        da1 = LIS_rc%da
        hr1 = LIS_rc%hr
        mn1 = 0
        ss1 = 0
        ts1 = -60*60     !previous hour
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

        yr2 = LIS_rc%yr
        mo2 = LIS_rc%mo
        da2 = LIS_rc%da
        hr2 = LIS_rc%hr
        mn2 = 0
        ss2 = 0
        ts2 = 0          !current hour
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        movetime = 1
        nldas3sw_struc(n)%findtime2 = 1
     endif
  else
     if (timenow.ge.nldas3sw_struc(n)%nldas3swtime2) then
        yr1 = LIS_rc%yr
        mo1 = LIS_rc%mo
        da1 = LIS_rc%da
        hr1 = LIS_rc%hr
        mn1 = 0
        ss1 = 0
        ts1 = 0          !current hour
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

        yr2 = LIS_rc%yr
        mo2 = LIS_rc%mo
        da2 = LIS_rc%da
        hr2 = LIS_rc%hr
        mn2 = 0
        ss2 = 0
        ts2 = 60*60      !next hour
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

        movetime = 1
        nldas3sw_struc(n)%findtime2 = 1
     endif
  endif

! Beginning of the run
  if ((LIS_rc%tscount(n).eq.1).or.(LIS_rc%rstflag(n).eq.1)) then
     nldas3sw_struc(n)%findtime1 = 1
     nldas3sw_struc(n)%findtime2 = 1
     movetime = 0
     LIS_rc%rstflag(n) = 0
  endif

  if (movetime.eq.1) then
     nldas3sw_struc(n)%nldas3swtime1 = nldas3sw_struc(n)%nldas3swtime2
     do f = 1,LIS_rc%met_nf(findex)
        do c = 1,LIS_rc%ngrid(n)
           nldas3sw_struc(n)%metdata1(f,c) =                           &
                             nldas3sw_struc(n)%metdata2(f,c)
        enddo
     enddo
  endif                     !end of movetime=1

! The following looks back 10 days, at the same hour to fill data gaps.
  if (nldas3sw_struc(n)%findtime1.eq.1) then
     ferrora = 0
     ferror = 0
     try = 0
     ts1 = -60*60*24
     do
        if (ferror.ne.0) exit
        try = try + 1
!- Obtaining CERES SWdown file:
        call ceres_nldas3swfile(n,1,findex,name_a,                     &
                nldas3sw_struc(n)%nldas3swfordir,yr1,mo1,da1,doy1,hr1)
        write(unit=LIS_logunit,fmt=*)                                  &
                "[INFO] getting CERES file: ",trim(name_a)
        order = 1
        call read_nldas3sw(n,1,findex,order,hr1,name_a,ferrora)

        ferror = ferrora
        if (ferror.ge.0) nldas3sw_struc(n)%nldas3swtime1 = time1
        call LIS_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if (try.gt.11) then
           write(LIS_logunit,*)                                        &
                "[ERR] CERES data gap exceeds 10 days on file 1"
           write(LIS_logunit,*) "[ERR] Program stopping ..."
           call LIS_endrun()
        endif
     enddo
!=== end of data search
  endif                     !end of LIS_rc%findtime=1

  if (nldas3sw_struc(n)%findtime2.eq.1) then
! The following looks back 10 days, at the same hour to fill data gaps.
     ferrora = 0
     ferror = 0
     try = 0
     ts2 = -60*60*24
     do
        if (ferror.ne.0) exit
        try = try+1

 !- Obtaining CERES SWdown file:
        call ceres_nldas3swfile(n,1,findex,name_a,                     &
                nldas3sw_struc(n)%nldas3swfordir,yr2,mo2,da2,doy2,hr2)
        write(unit=LIS_logunit,fmt=*)                                  &
                "[INFO] getting CERES file: ",trim(name_a)
        order = 2
        call read_nldas3sw(n,1,findex,order,hr2,name_a,ferrora)

        ferror = ferrora
        if (ferror.ge.0) then
           nldas3sw_struc(n)%nldas3swtime2 = time2
        endif
        call LIS_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if (try.gt.11) then
           write(LIS_logunit,*)                                        &
                "[ERR] CERES data gap exceeds 10 days on file 2"
           write(LIS_logunit,*) "[ERR] Program stopping ..."
           call LIS_endrun()
        endif
     enddo
!=== end of data search
  endif                     ! end of findtime2=1

end subroutine get_nldas3sw

