!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center 
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE:  get_chirps
!  \label{get_chirps}
!
! !REVISION HISTORY:
!  26 Jan 2007: Hiroko Kato; Initial Specification adopted from LDT
!  15 Jul 2015: K. R. Arsenault;  Adapted for CHIRPS
!  29 Oct 2025: J. Nattala; Updated to support both CHIRPS v2.0 and v3.0;
!                           added actual arg. version (return_chirps_filename)
! 
! !INTERFACE:
subroutine get_chirps(n,findex)
! !USES:
  use LDT_coreMod,         only : LDT_rc
  use LDT_metforcingMod,   only : LDT_forc
  use LDT_timeMgrMod,      only : LDT_get_nstep, LDT_tick
  use LDT_logMod,          only : LDT_logunit, LDT_endrun, LDT_verify
  use chirps_forcingMod,   only : chirps_struc
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates daily, 0.05 or 0.25 deg 
!  CHIRPS precipitation forcing. At the beginning of a simulation, 
!  the code reads the most recent past data (nearest 3 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day
!.
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    call to advance or retract time
!  \item[read\_chirps](\ref{read_chirps}) \newline
!    call to read the CHIRPS data and perform spatial interpolation 
!  \end{description}
!EOP  
   
  integer :: openfile2_flag
  integer :: ferror_chirps     ! Error flags for precip data sources
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real*8  :: time1,time2
  real    :: gmt1,gmt2,ts1,ts2

  character(len=LDT_CONST_PATH_LEN) :: chirps_filename
  logical       :: file_exists

! ___________________________________________________________________

  openfile2_flag = 0

  !=== Determine Required Data Times (The previous hour & the future hour)
  yr1=LDT_rc%yr    ! Current forcing time
  mo1=LDT_rc%mo
  da1=LDT_rc%da
  hr1=0
  mn1=0
  ss1=0
  ts1=0
  call LDT_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  yr2=LDT_rc%yr    ! Next forcing time
  mo2=LDT_rc%mo
  da2=LDT_rc%da
  hr2=0
  mn2=0
  ss2=0
  ts2=24*3600
  call LDT_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

! Beginning of the run:
  if( LDT_get_nstep(LDT_rc,n) == 1 .or. chirps_struc(n)%reset_flag ) then
     write(LDT_logunit,*) "[INFO] Valid CHIRPS data value reached ... "
     openfile2_flag = 1
     chirps_struc(n)%chirpstime2 = time2
     chirps_struc(n)%reset_flag = .false.
  endif

  !=== Check if time interval boundary was crossed
  if( LDT_rc%nts(n) >= 86400. ) then   ! LDT timestep - daily or greater
     if( LDT_rc%time >= chirps_struc(n)%chirpstime2 ) then  ! ORIG
       write(LDT_logunit,*) "[INFO] Valid CHIRPS data value reached ... "
       openfile2_flag = 1
       chirps_struc(n)%chirpstime2 = time2
     endif

  elseif( LDT_rc%nts(n) < 86400. ) then  ! LDT timesteps < daily
     if( LDT_rc%time > chirps_struc(n)%chirpstime2 ) then       ! UPDATED
       write(LDT_logunit,*) "[INFO] Valid CHIRPS data value reached ... "
       openfile2_flag = 1
       chirps_struc(n)%chirpstime2 = time2
     endif
  endif

  ! Open Next CHIRPS file:
  if( openfile2_flag == 1 ) then
     call return_chirps_filename( chirps_filename, chirps_struc(n)%directory, &
       yr1, mo1, chirps_struc(n)%xres, chirps_struc(n)%version )

     ! Determine if the netcdf file exists:
     inquire( file=trim(chirps_filename), exist=file_exists)
     if( file_exists ) then
        write(LDT_logunit,*) "[INFO] Opening the CHIRPS file: ",trim(chirps_filename)
        call read_chirps( n, findex, chirps_filename, LDT_rc%yr, &
                           LDT_rc%mo, LDT_rc%da, ferror_chirps )
     else
        write(LDT_logunit,*) "[WARN] CHIRPS file missing: ",trim(chirps_filename)
        write(LDT_logunit,*) "[WARN] Make sure to include first portion of filename."
        write(LDT_logunit,*) " for example, 'chirps-v3.0' or 'chirp'. "
        write(LDT_logunit,*) "[WARN] For now, relying only on baseforcing precipitation ..."
        if( LDT_rc%nmetforc == 1 ) then
          write(LDT_logunit,*)"[ERR] No other underlying forcings and no available CHIRPS precipitation files."
          write(LDT_logunit,*)" Program halting ..."
          call LDT_endrun
        endif
     endif
  endif

  ! IF first timestep, Assign Metforcing Data 2 to Metforcing Data 1:
  if( LDT_get_nstep(LDT_rc,n) == 1 ) then
     LDT_forc(n,findex)%metdata1(:,:) = LDT_forc(n,findex)%metdata2(:,:)
  endif

end subroutine get_chirps


!BOP
! !ROUTINE: return_chirps_filename
! \label{return_chirps_filename}
!
! !INTERFACE:
subroutine return_chirps_filename( filename, dirpath, yr, mo, res, version )

  use LDT_logMod, only : LDT_logunit
  implicit none

! !ARGUMENTS: 
  character(len=1024), intent(out)  :: filename
  character(len=*), intent(in)     :: dirpath
  integer, intent(in)              :: yr, mo, version
  real, intent(in) :: res

! !DESCRIPTION:
!   This subroutine puts together the CHIRPS file name
!
!  The arguments are:
!  \begin{description}
!  \item[dirpath]
!    Name of the CHIRPS directory pathname
!  \item[yr]
!    year 
!  \item[mo]
!    month
!  \item[res]
!    CHIRPS dataset resolution
!  \item[version]
!    CHIRPS dataset version (2 or 3)
!  \item[filename]
!    name of the timestamped CHIRPS file
!  \end{description}
!
!EOP

   character*4  :: cyr
   character*2  :: cmo
   character*3  :: cres

!=== End Variable Definition ===============

   write(unit=cyr, fmt='(i4.4)')  yr
   write(unit=cmo, fmt='(i2.2)')  mo

   if( res == 0.05 ) then
      cres = "p05"
   elseif( res == 0.25 ) then
      cres = "p25"
   endif

  !=== Assemble CHIRPS filename:
  !    e.g., chirps-v2.0.1981.days_p05.nc  or  chirps-v3.0.2025.days_p05.nc
  filename = trim(dirpath)//"."//trim(cyr)//".days_"//trim(cres)//".nc"

end subroutine return_chirps_filename
