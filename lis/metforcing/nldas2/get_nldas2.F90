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
! !ROUTINE: get_nldas2
! \label{get_nldas2}
!
!
! !REVISION HISTORY:
! 02 Feb 2004: Sujay Kumar; Initial Version in LIS
! 22 Aug 2007: Chuck Alonge; Updated for NLDAS2 Forcing
! 22 Jan 2012: K. Arsenault; Accommodate GESDISC, NCEP filename conventions
! 14 Mar 2014: David Mocko: Trim writing of filenames to log file
! 
! !INTERFACE:
subroutine get_nldas2(n,findex)
! !USES:
  use LIS_coreMod, only        : LIS_rc, LIS_domain
  use LIS_timeMgrMod,     only : LIS_tick
  use LIS_metforcingMod,  only : LIS_forc
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use nldas2_forcingMod,  only : nldas2_struc
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates hourly, NLDAS-2 forcing (NARR based). 
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
!    determines the NLDAS-2 data times
!  \item[gesdisc\_nldas2filea](\ref{gesdisc_nldas2filea}) \newline
!    Puts together appropriate timestamped GES DISC filename - a data
!  \item[gesdisc\_nldas2fileb](\ref{gesdisc_nldas2fileb}) \newline
!    Puts together appropriate timestamped GES DISC filename - b data
!  \item[ncep\_nldas2filea](\ref{ncep_nldas2filea}) \newline
!    Puts together appropriate timestamped NCEP/EMC filename - a data
!  \item[ncep\_nldas2fileb](\ref{ncep_nldas2fileb}) \newline
!    Puts together appropriate timestamped NCEP/EMC filename - b data
!  \item[read\_nldas2a](\ref{read_nldas2a}) \newline
!    Reads and Interpolates NLDAS-2 A data to LIS grid
!  \item[read\_nldas2b](\ref{read_nldas2b}) \newline
!    Reads and Interpolates NLDAS-2 B data to LIS grid
!  \end{description}
!EOP

  integer :: c,f,ferrora,ferrorb,ferror,try
  integer :: order
  integer :: readbfile
  real*8  :: time1,time2,timenow
  real*8  :: dtime1, dtime2
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  character(len=LIS_CONST_PATH_LEN) :: name_a,name_b
  real    :: gmt1,gmt2,ts1,ts2
  integer :: movetime     ! 1=move time 2 data into time 1  
  integer :: kk           ! Forecast member index

!=== End Variable Definition =============================================
  try=-999
  
!====Check to see if b-file needs to be opened
  readbfile = 0
  if(nldas2_struc(n)%model_level_data .gt. 0 .or. & 
     nldas2_struc(n)%model_dswrf_data .gt. 0 .or. &
     nldas2_struc(n)%model_level_press .gt. 0 .or. &
     nldas2_struc(n)%model_pcp_data .gt. 0 ) then     
      readbfile = 1
  endif

!====Assumption will be not to find or move any data
  nldas2_struc(n)%findtime1=0
  nldas2_struc(n)%findtime2=0
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
     if(timenow.ge.nldas2_struc(n)%nldas2time2) then 
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
        nldas2_struc(n)%findtime2 = 1
     endif
  else
     if(timenow.ge.nldas2_struc(n)%nldas2time2) then 
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
        nldas2_struc(n)%findtime2 = 1
     endif
  endif

  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then    !beginning of the run	
     nldas2_struc(n)%findtime1=1
     nldas2_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
  
  if(movetime.eq.1) then
     nldas2_struc(n)%nldas2time1=nldas2_struc(n)%nldas2time2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           nldas2_struc(n)%metdata1(:,f,c)=nldas2_struc(n)%metdata2(:,f,c)
        enddo
     enddo
  endif    !end of movetime=1
  
  if(nldas2_struc(n)%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferrora=0
     ferrorb=0
     ferror=0
     try=0  
     ts1=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        !- Obtaining NLDAS-2 File-A:
        do kk= nldas2_struc(n)%st_iterid, nldas2_struc(n)%en_iterid
           if( trim(nldas2_struc(n)%nldas2_filesrc)== "GES-DISC") then    ! GES-DISC
              call gesdisc_nldas2filea(n,kk,findex,&
                    name_a,nldas2_struc(n)%nldas2dir,&
                    yr1,mo1,da1,doy1,hr1)
           elseif( trim(nldas2_struc(n)%nldas2_filesrc) == "NCEP") then   ! NCEP
              call ncep_nldas2filea(n,kk,findex, &
                    name_a,nldas2_struc(n)%nldas2dir,&
                    yr1,mo1,da1,hr1)
           endif
           write(unit=LIS_logunit,fmt=*)'[INFO] getting file1a.. ',trim(name_a)
           order = 1
           call read_nldas2a(n,kk,findex,order,mo1,name_a,ferrora)
        enddo

        if( readbfile .gt. 0) then
          do kk= nldas2_struc(n)%st_iterid, nldas2_struc(n)%en_iterid
             if( trim(nldas2_struc(n)%nldas2_filesrc) == "GES-DISC" ) then  ! GES-DISC
                call gesdisc_nldas2fileb(n,kk,findex, &
                      name_b,nldas2_struc(n)%nldas2dir,&
                      yr1,mo1,da1,doy1,hr1)
             elseif( trim(nldas2_struc(n)%nldas2_filesrc) == "NCEP" ) then  ! NCEP
                call ncep_nldas2fileb(n,kk,findex,&
                      name_b,nldas2_struc(n)%nldas2dir,&
                      yr1,mo1,da1,hr1)
             endif
             write(unit=LIS_logunit,fmt=*)'[INFO] getting file1b.. ',trim(name_b)
             call read_nldas2b(n,kk,findex,order,name_b,ferrorb)
          enddo
        else
           ferrorb = 1
        endif
        ferror = ferrora + ferrorb
        if(ferror.ge.1) nldas2_struc(n)%nldas2time1=time1
        call LIS_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(*,*)'error: NLDAS-2 data gap exceeds 10 days on file 1'
           stop
        endif
     enddo
!=== end of data search
  endif   !end of LIS_rc%findtime=1	   	


  if(nldas2_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferrora=0
     ferrorb=0
     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1

     !- Obtaining NLDAS-2 File-A:
        do kk= nldas2_struc(n)%st_iterid, nldas2_struc(n)%en_iterid
          if( trim(nldas2_struc(n)%nldas2_filesrc) == "GES-DISC" ) then   ! GES-DISC
            call gesdisc_nldas2filea(n,kk,findex, &
                  name_a,nldas2_struc(n)%nldas2dir,&
                  yr2,mo2,da2,doy2,hr2)
          elseif( trim(nldas2_struc(n)%nldas2_filesrc) == "NCEP" ) then   ! NCEP
            call ncep_nldas2filea(n,kk,findex,&
                  name_a,nldas2_struc(n)%nldas2dir,&
                  yr2,mo2,da2,hr2)
          endif
          write(unit=LIS_logunit,fmt=*)'[INFO] getting file2a.. ',trim(name_a)
          order = 2
          call read_nldas2a(n,kk,findex,order,mo2,name_a,ferrora)
        end do

        if( readbfile .gt. 0) then
          do kk= nldas2_struc(n)%st_iterid, nldas2_struc(n)%en_iterid
             if( trim(nldas2_struc(n)%nldas2_filesrc) == "GES-DISC" ) then   ! GES-DISC
                call gesdisc_nldas2fileb(n,kk,findex,&
                      name_b,nldas2_struc(n)%nldas2dir,&
                      yr2,mo2,da2,doy2,hr2)
             elseif( trim(nldas2_struc(n)%nldas2_filesrc) == "NCEP" ) then   ! NCEP
                call ncep_nldas2fileb(n,kk,findex,&
                      name_b,nldas2_struc(n)%nldas2dir,&
                      yr2,mo2,da2,hr2)
             end if
             write(unit=LIS_logunit,fmt=*)'[INFO] getting file2b.. ',trim(name_b)
             call read_nldas2b(n,kk,findex,order,name_b,ferrorb)
          enddo
        else
           ferrorb = 1
        endif
        ferror = ferrora + ferrorb
        if(ferror.ge.1) then
           nldas2_struc(n)%nldas2time2=time2
        endif
        call LIS_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(*,*)'error: NLDAS-2 data gap exceeds 10 days on file 2'
           stop
        endif
     enddo
     !=== end of data search
  endif   ! end of findtime2=1

end subroutine get_nldas2
