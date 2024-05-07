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
! !ROUTINE: AGRMET_smiest
!  \label{AGRMET_smiest}
!
! !REVISION HISTORY:
!
!    04 dec 97  initial version................capt andrus/dnxm(agromet)
!    09 may 00  modified code to read in values from ported version
!               of smiedr..............................capt hidalgo/dnxm
! 29Jul2005 Sujay Kumar, Initial Code in LIS
! 20Dec2007 Marv Freimund, Simplify filename creation and correct bug
! 11Mar2010 Changed program names in messages to LIS.  Left argument
!           to LIS_alert with old program name to keep alert numbers
!           unique............................Chris Franks/16WS/WXE/SEMS
! 31 MAR 2010 Add handling of multiple resolutions ...see in code comments
!             for where........................Michael Shaw/WXE
!
! !INTERFACE:    
subroutine AGRMET_smiest( n,j3hr, quad9r, ra, razero,alert_number,imax,jmax)
                                                ! Michael Shaw - added imax and 
                                                ! jmax; could probably avoid this
                                                ! with current version of code,
                                                ! though was used to ensure proper
                                                ! dimensioning on multi processors
                                                ! in previous versions 
! !USES:
  use AGRMET_forcingMod, only : agrmet_struc
  use LIS_coreMod,   only : LIS_rc, LIS_masterproc
  use LIS_fileIOMod, only : LIS_putget
  use LIS_timeMgrMod, only : LIS_julhr_date
  use LIS_logMod, only : lis_alert, LIS_logunit

  implicit none
! !ARGUMENTS: 
  integer,       intent(in)         :: n,imax,jmax
  integer                           :: ii,jj
  integer,       intent(in)         :: j3hr
  real,          intent(out)        :: ra(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                              :: ra_temp(1,imax,jmax)
  integer,       intent(in)         :: razero  
  real,          intent(in)         :: quad9r   
  integer                           :: alert_number 

!
! !DESCRIPTION:
!
!    to make an ssmi based estimate array.  
!
!    \textbf{Method}
!    
!    - retrieve rain rates for the current hemisphere from file. \newline
!    - the files retrieved will be 3hrly amounts (mm) \newline
!    - if the razero flag equals 0, reset the ssmi zeros to quad9r \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[n]  index of nest
!   \item[alert\_number]  alert message number
!   \item[exists]     a logical that indicates whether or not a file exists
!   \item[hemi]       hemisphere (1=nh, 2=sh)
!   \item[ifil]       input ssmi data file name
!   \item[message]    array of alert messages
!   \item[quad9r]     parameter for 9999.0 value  
!   \item[ra]         ssmi based estimated precip amount (mm/3hr)
!   \item[razero]     flag to use zero ssmi amounts (1=use zeros)
!   \item[use\_zeros]  logical variable that is true if ssmi zeros are
!                      to be used
!  \end{description}
!
!
!  The routines invoked are: 
!  \begin{description}
!  \item[julhr\_date] (\ref{LIS_julhr_date}) \newline
!   converts julian hour to a date format
!  \item[agrmet\_ssmiprec\_filename](\ref{agrmet_ssmiprec_filename}) \newline
!   generates the SSM/I filename
!  \item[LIS\_putget](\ref{LIS_putget}) \newline
!   read the SSM/I data
!  \item[LIS\_alert](\ref{LIS_alert}) \newline
!   print out an alert message in case of missing data
!  \end{description}
!EOP
  real            :: ra_tmp(2,imax,jmax)
  character*100   :: ifil
  character*100   :: message(20)
  character*30    :: routine_name
  logical         :: exists
  logical         :: use_zeros
  integer         :: yr,mo,da,hr
  integer         :: hemi,first,last
  integer         :: ip,o,p
  real            :: udef

  data routine_name     / 'AGRMET_smiest' /
!     ------------------------------------------------------------------
!     executable code begins here...set logical variable for use of 
!     ssmi zero values.
!     ------------------------------------------------------------------

  ra = quad9r
  udef = 9999
  ra_tmp = udef
  use_zeros = .false. ! EMK BUG FIX
  if (razero .eq. 1) use_zeros = .true.

! Michael Shaw - Similar to geoest, adding dependency on global or hemispheric land mask
! of first and lass "hemispheres" of mask, e.g.

  if (agrmet_struc(n)%global_or_hemi.eq.0) then
     first=agrmet_struc(n)%shemi
     last=agrmet_struc(n)%nhemi
  else
     first=1
     last=1
  endif
!     ------------------------------------------------------------------
!     determine the file containing the data for julhr.  check if it
!     exists.
!     ------------------------------------------------------------------
  do hemi=first,last
     call LIS_julhr_date( j3hr, yr,mo,da,hr)
     call agrmet_ssmiprec_filename(ifil,agrmet_struc(n)%agrmetdir,&
          agrmet_struc(n)%ssmidir,agrmet_struc(n)%use_timestamp,&
          hemi,yr,mo,da,hr,n,imax,jmax)

     inquire(file = trim(ifil), exist = exists)
     if (.not. exists) then
        write(LIS_logunit,*) ' '
        write(LIS_logunit,*)'[WARN] precip/smiedr:  error opening file'
        write(LIS_logunit,*)'[WARN]  SSMI data file ', trim(ifil),' does not exist.'
        write(LIS_logunit,*)'[WARN]  SSMI estimates will not be used in ',&
             'precip analysis.'
        write(LIS_logunit,*) ' '
        message   =' '
        message(1)='program:  LIS'
        message(2)='  routine:  AGRMET_smiest'
        message(3)='  SSMI data file '//trim(ifil)//' does not exist.'
        message(4)='  SSMI estimates will not be used in '//&
             'precip analysis.'
        alert_number = alert_number + 1
        if(LIS_masterproc) & 
             call lis_alert('precip              ', alert_number, message )

        ra_tmp(hemi,:,:) = udef
     else
        write(LIS_logunit,*) '[INFO] READING ',trim(ifil)
        call LIS_putget( ra_tmp(hemi,:,:), 'r', ifil, routine_name, &
             imax, jmax)
     endif
     
  enddo

  if ( .not. use_zeros ) then
     write(LIS_logunit,*)'[INFO] SSMI ZEROS NOT USED'
     where ( ra_tmp .eq. 0.0 ) ra_tmp = udef 
  endif

  ip = 1  
  call interp_agrmetvar(n,ip,ra_tmp,udef,ra,imax,jmax)

!     ------------------------------------------------------------------
!     loop through all points and adjust the zero values as necessary
!     ------------------------------------------------------------------
  
  return
end subroutine AGRMET_smiest

!BOP
! 
! !ROUTINE: agrmet_ssmiprec_filename
! \label{agrmet_ssmiprec_filename}
! 
! !INTERFACE: 
subroutine agrmet_ssmiprec_filename(fname,rootdir,dir,use_timestamp,&
     hemi,yr,mo,da,hr,n,imax,jmax) ! Distinguish what grid - legacy 
                                   ! from previous versions to ensure
                                   ! using correct files on multiprocessors

  use AGRMET_forcingMod, only : agrmet_struc
  implicit none
! !ARGUMENTS: 
  character(*)        :: fname
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: hemi,n,imax,jmax
  integer, intent(in) :: yr,mo,da,hr

! 
! !DESCRIPTION: 
!  This routines generates the name of the SSM/I file to be read, 
!  with the appropriate hemisphere and time stamps. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[fname]
!    created filename
!   \item[rootdir]
!    path to the root directory containing the data
!   \item[dir]
!    name of the subdirectory containing the data
!   \item[use\_timestamp]
!    flag to indicate whether the directories 
!    should be timestamped or not
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[yr]
!    4 digit year
!   \item[mo]
!    integer value of month (1-12)
!   \item[da]
!    day of the month
!   \item[hr]
!    hour of the day
!  \end{description}
!EOP

  character( 2), parameter :: fhemi(2) = (/'nh','sh'/)
  character(10)            :: ftime1, ftime2

  write (UNIT=ftime2, FMT='(i4, i2.2, i2.2, i2.2)') yr, mo, da, hr

  if(imax == 512)then
  if (use_timestamp .eq. 1) then 
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/ssmira_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  else
     fname = trim(rootdir) //  '/'   // trim(dir) // '/ssmira_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  endif
  elseif(imax == 1024)then
  if (use_timestamp .eq. 1) then 
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/ssmira_' // &
             fhemi(hemi) // '.16.03hr.' // ftime2 // '.dat'
  else
     fname = trim(rootdir) //  '/'   // trim(dir) // '/ssmira_' // &
             fhemi(hemi) // '.16.03hr.' // ftime2 // '.dat'
  endif
  elseif(imax == 1440)then
  if (use_timestamp .eq. 1) then
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // '/smiedr/ssmira_' // &
             '0p25deg' // '.03hr.' // ftime2 // '.dat'
  else
     fname = trim(rootdir) //  '/'   // trim(dir) // '/ssmira_' // &
             '0p25deg' // '.03hr.' // ftime2 // '.dat'
  endif
  endif

end subroutine agrmet_ssmiprec_filename
