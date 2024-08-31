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
! !ROUTINE: AGRMET_fgfill
! \label{AGRMET_fgfill}
! 
! !REVISION HISTORY:
!      7 oct 97  initial version............ssgt mccormick,mr moore/sysm
!     28 feb 99  changed array structure to supergrid. corrected 
!                temporal interpolation. removed equivalencing and 
!                simplified array structures.  eliminated z6 and z6m
!                ..........................................mr moore/dnxm
!     22 Jul 99  ported to IBM SP-2.  added intent attributes to all 
!                arguments.  changed to grid specific arrays to
!                be dynamically allocatable................mr gayno/dnxm
!     01 aug 02  added read of data_levels namelist........mr gayno/dnxm
!     17 Jun 17  Added GRIB ground height and 2-m T, RH; removed
!                land mask.................................Eric Kemp/GSFC
! !INTERFACE: 
! EMK...Removed land mask
subroutine AGRMET_fgfill( n, alt, sfctmp, sfcrlh, sfcspd, sfcprs, &
     jul, lastmp, lasrlh, lasspd, lasprs, &
     step, order, minwnd, imax,jmax, &
     kprs1, agr_tmp_c, agr_hgt_c, agr_rh_c, &
     agr_tmp_sfc_c, agr_hgt_sfc_c, agr_rh_sfc_c,&
     agr_wspd_c,agr_pres_c,&
     agr_tmp_p, agr_hgt_p, agr_rh_p, &
     agr_tmp_sfc_p, agr_hgt_sfc_p, agr_rh_sfc_p, &
     agr_wspd_p,agr_pres_p,pathfb)

! !USES:
  use LIS_coreMod,  only : LIS_rc , LIS_localPet
  use LIS_topoMod,  only : LIS_topo

  implicit none
! !ARGUMENTS:
  integer,   intent(in)           :: n 
  integer,   intent(in)           :: imax
  integer,   intent(in)           :: jmax
  integer,   intent(in)           :: order
  real,      intent(in)           :: alt (imax,jmax)
  real,      intent(inout)        :: sfcprs ( 6,imax,jmax )
  real,      intent(inout)        :: sfcrlh ( 6,imax,jmax )
  real,      intent(inout)        :: sfcspd ( 6,imax,jmax )
  real,      intent(inout)        :: sfctmp ( 6,imax,jmax )
  integer,   intent(in)           :: jul  
  real,      intent(inout)        :: lasprs ( imax,jmax )
  real,      intent(inout)        :: lasrlh ( imax,jmax )
  real,      intent(inout)        :: lasspd ( imax,jmax )
  real,      intent(inout)        :: lastmp ( imax,jmax )
  integer,   intent(inout)        :: step
  real,      intent(in)           :: minwnd
  integer,     intent(in)         :: kprs1
  real,        intent(in)    :: agr_hgt_c( imax,jmax,kprs1 )
  real,        intent(in)    :: agr_rh_c ( imax,jmax,kprs1 )
  real,        intent(in)    :: agr_tmp_c( imax,jmax,kprs1 )
  real,        intent(in)    :: agr_hgt_sfc_c( imax,jmax)
  real,        intent(in)    :: agr_rh_sfc_c ( imax,jmax)
  real,        intent(in)    :: agr_tmp_sfc_c( imax,jmax)
  real,        intent(in)    :: agr_wspd_c( imax,jmax)
  real,        intent(in)    :: agr_pres_c( imax,jmax)
  real,        intent(in)    :: agr_hgt_p ( imax,jmax,kprs1 )
  real,        intent(in)    :: agr_rh_p  ( imax,jmax,kprs1 )
  real,        intent(in)    :: agr_tmp_p ( imax,jmax,kprs1 )
  real,        intent(in)    :: agr_hgt_sfc_p ( imax,jmax)
  real,        intent(in)    :: agr_rh_sfc_p  ( imax,jmax)
  real,        intent(in)    :: agr_tmp_sfc_p ( imax,jmax)
  real,        intent(in)    :: agr_wspd_p( imax,jmax)
  real,        intent(in)    :: agr_pres_p( imax,jmax)
  character(len=*), intent(in) :: pathfb


! !DESCRIPTION:
!     to determine a first guess field for the current time which will
!     later be blended with observations to make a final analysis.
!   
!     \textbf{Method} \newline
!     - the first time through the 6-hour cycle (t+1) \newline
!       - read data\_levels namelist to get the isobaric levels
!         for the first guess data. \newline
!       - retrieve the first guess data for the end
!         of the 6 hour cycle (t+6). (e.g. 12z data for a 12z run.) \newline
!       - retrieve the first guess data for the start of 
!         the 6 hour cycle (t+0). (e.g. 06z data for a 12z run.)
!       - time interpolate to the current julian hour using the
!         retrieved first guess data from t+0 and t+6. \newline
!     - for all other hours in the 6-hour cycle (t+2 to t+6): \newline
!       - time interpolate to the current hour using the barnes
!         analysis from the previous hour and first guess data from 
!         the end of the 6 hour cycle. \newline
!    
!  The arguments and variables are: 
!  \begin{description}
!  \item[land]
!   array of land mass points for agrmet grid
!  \item[jul]
!    current julian hour
!  \item[sfcprs]
!    first guess pressure fields on the AGRMET grid interpolated to 
!    the current time
!  \item[sfcrlh]
!    first guess rh fields on the AGRMET grid interpolated to 
!    the current time
!  \item[sfcspd]
!    first guess wind speed fields on the AGRMET grid interpolated to 
!    the current time   
!  \item[sfctmp]
!    first guess temperature fields on the AGRMET grid interpolated to 
!    the current time
!  \item[lasprs]
!    first guess pressure fields for ending 6 hour time
!  \item[lasrlh]
!    first guess rh fields for ending 6 hour time
!  \item[lasspd]
!    first guess wind speed fields for ending 6 hour time
!  \item[lastmp]
!    first guess temperature fields for ending 6 hour time
!  \item[step]
!    time loop counter used to calculate time weights
!  \item[order]
!   flag to indicate which data is to be read
!  \item[minwnd]
!    minimum allowable windspeed on the agrmet grid   
!  \item[imax]
!    east/west dimension of agrmet grid
!  \item[jmax]
!    north/south dimension of agrmet grid
!  \item[pathfb]
!    directory path of the surface calculations
!  \item[kprs1]
!    number of isobaric levels
!  \item[agr\_tmp\_c]
!    temperature on the AGRMET grid for the current time
!  \item[agr\_hgt\_c]
!    geopotential heights on the AGRMET grid for the current time
!  \item[agr\_rh\_c] 
!    relative humidity on the AGRMET grid for the current time
!  \item[agr\_wspd\_c]
!    wind speeds on the AGRMET grid for the current time
!  \item[agr\_tmp\_p]
!    temperature on the AGRMET grid from 6 hours ago
!  \item[agr\_hgt\_p]
!    geopotential heights on the AGRMET grid from 6 hours ago
!  \item[agr\_rh\_p] 
!    relative humidity on the AGRMET grid from 6 hours ago
!  \item[agr\_wspd\_p]
!    wind speeds on the AGRMET grid from 6 hours ago
!  \item[action]
!    file i/o action - reading or opening
!  \item[ctime]
!    character equivalent of variable time
!  \item[i]
!   loop counter
!  \item[julm6]
!   current julian time of run minus 6 hours 
!  \item[nmlfil]
!   name, including path, of data\_levels namelist
!  \item[prslvls]
!   isobaric levels in fldbld data
!  \item[time]
!   zulu hour of data\_levels file
!  \item[wt1,wt2]
!   weight factors used for time smoothing prior hour's
!   first guess data
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[AGRMET\_sfcval](\ref{AGRMET_sfcval}) \newline
!    vertically interpolate the first guess fields
!  \end{description}
!EOP
  
  integer                         :: i
  integer                         :: j
  integer                         :: julm6
  integer                         :: kprs
  integer                         :: prslvls (30)
  integer                         :: time
  real                            :: wt1
  real                            :: wt2
  character*2                     :: ctime
  integer                         :: step_p
  integer                         :: step_c
  data prslvls / 1000,975,950,925,900,850,800,750,700,650,600,550,500,&
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
  kprs = 13
!     ------------------------------------------------------------------
!     executable code starts here ...
!     if this is the first time thru this routine (t+1), call sfcval
!     to read in first guess data for current time and 6 hours ago.
!     but first, read in data_levels namelist to get the first guess
!     isobaric levels for this time period.
!     ------------------------------------------------------------------

  if( step .eq. 1 .or. order.eq.2)then
     
!     kprs    = 0
!     prslvls = 0
     
     time = mod(jul,24)

     write(ctime,'(i2.2)') time
     
 !    nmlfil = trim(pathfb) // "/data_levels." // ctime // 'z'

 !    write(LIS_logunit,*)' '
 !    write(LIS_logunit,*)'- OPENING ', trim(nmlfil)

 !    action = 'opening'
 !    open(22, file=trim(nmlfil), iostat=istat)
 !    action = 'reading'
 !    read(22, nml=fg_info, iostat = istat)
 !    close(22)
     
!     ------------------------------------------------------------------
!       call sfcval to calculate read in fldbld data and vertically 
!       interpolate it to the agrmet grid for the final hour of
!       the 6-hour cycle (t+6) (end period).
!     ------------------------------------------------------------------
     ! EMK...Removed land mask
!     call AGRMET_sfcval( alt, jul, &!land, &
!          minwnd, lasprs, lasrlh, &
!          lasspd, lastmp, imax, jmax, &
!          kprs, prslvls, kprs1, &
!          agr_tmp_c, agr_hgt_c, agr_rh_c, &
!          agr_tmp_sfc_c, agr_hgt_sfc_c, agr_rh_sfc_c, &
!          agr_wspd_c, agr_pres_c)
     call USAF_sfcval( alt, &
          minwnd, lasprs, lasrlh, &
          lasspd, lastmp, imax, jmax, &
          kprs, prslvls, kprs1, &
          agr_tmp_c, agr_hgt_c, agr_rh_c, &
          agr_tmp_sfc_c, agr_hgt_sfc_c, agr_rh_sfc_c, &
          agr_wspd_c, agr_pres_c)


!      print*,'EMK: LIS_localPet, maxval(agr_rh_c) = ', &
!           LIS_localPet, maxval(agr_rh_c)
!      print*,'EMK: LIS_localPet, minval(agr_rh_c) = ', &
!           LIS_localPet, minval(agr_rh_c)

!     ------------------------------------------------------------------
!       call sfcval again to load the data from the beginning of the
!       6-hour period (t+0) (start period).  but first, read in 
!       data_levels namelist to get the first guess isobaric levels
!       for this time period.
!     ------------------------------------------------------------------

     julm6 = jul - 6 
     
 !    kprs    = 0
 !    prslvls = 0
     
     time = mod(julm6,24)
     write(ctime,'(i2.2)') time
     
 !    nmlfil = trim(pathfb) // "/data_levels." // ctime // 'z'
 !    
 !    write(LIS_logunit,*)' '
 !    write(LIS_logunit,*)'- OPENING ',trim(nmlfil)
 !    
 !    action = 'opening'
 !    open(22, file=trim(nmlfil), iostat=istat)
 !    action = 'reading'
 !    read(22, nml=fg_info, iostat=istat)
 !    close(22)

     ! EMK...Removed land mask
!     call AGRMET_sfcval( alt, julm6, &!land, 
!          minwnd, sfcprs(1,:,:), &
!          sfcrlh(1,:,:), sfcspd(1,:,:), sfctmp(1,:,:), imax, jmax, &
!          kprs, prslvls, kprs1, &
!          agr_tmp_p, agr_hgt_p, agr_rh_p, &
!          agr_tmp_sfc_p, agr_hgt_sfc_p, agr_rh_sfc_p, &
!          agr_wspd_p, agr_pres_p)
     call USAF_sfcval( alt, &
          minwnd, sfcprs(1,:,:), &
          sfcrlh(1,:,:), sfcspd(1,:,:), sfctmp(1,:,:), imax, jmax, &
          kprs, prslvls, kprs1, &
          agr_tmp_p, agr_hgt_p, agr_rh_p, &
          agr_tmp_sfc_p, agr_hgt_sfc_p, agr_rh_sfc_p, &
          agr_wspd_p, agr_pres_p)


!      print*,'EMK: LIS_localPet, maxval(agr_rh_p) = ', &
!           LIS_localPet, maxval(agr_rh_p)
!      print*,'EMK: LIS_localPet, minval(agr_rh_p) = ', &
!           LIS_localPet, minval(agr_rh_p)
     
  endif

!     ------------------------------------------------------------------
!     time-interpolate the first guess fields to the current julian
!     hour.  if this is the hour 1 (t+1) of the 6 hour cycle, 
!     then interpolate between the first guess data at start (t+0) and
!     end periods (t+6). for all subsequent times, the 
!     interpolation is between the first guess at t+6
!     and the previous barnes analysis.  if this is the last
!     time through (t+6), the previous barnes analysis is used equally
!     with the end period first guess.  this method of time 
!     interpolation, which incorporates observations,
!     was found to give better results than using only 
!     the pure first guess data.
!     ------------------------------------------------------------------

  if(step.eq.1) then
     step_p = 1
     step_c = 1
  else
     step_c = step
     step_p = step-1
  endif

  if( step .eq. 6 ) step = 5

  wt1 = 6.0 - float(step)
  wt2 = 1.0 / (wt1 + 1.0)

  do j = 1, jmax
     do i = 1, imax
!        if( land(i,j) .gt. 0 )then
           sfctmp(step_c,i,j) = wt2 * ( wt1 * sfctmp(step_p,i,j) +&
                lastmp(i,j) )
           sfcrlh(step_c,i,j) = wt2 * ( wt1 * sfcrlh(step_p,i,j) +&
                lasrlh(i,j) )
           sfcspd(step_c,i,j) = wt2 * ( wt1 * sfcspd(step_p,i,j) +&
                lasspd(i,j) )
           sfcprs(step_c,i,j) = wt2 * ( wt1 * sfcprs(step_p,i,j) +&
                lasprs(i,j) )
!        else
!           sfctmp(step_c,i,j) = -1.0
!           sfcrlh(step_c,i,j) = -1.0
!           sfcspd(step_c,i,j) = -1.0
!           sfcprs(step_c,i,j) = -1.0
!           
!        endif
     enddo
  enddo
  step = step + 1  

!  print*,'EMK: LIS_localPet, step_c, maxval(sfcrlh) = ', &
!       LIS_localPet, step_c,maxval(sfcrlh(step_c,:,:))
!  print*,'EMK: LIS_localPet, step_c, minval(sfcrlh) = ', &
!       LIS_localPet, step_c, minval(sfcrlh(step_c,:,:))

  return  
end subroutine AGRMET_fgfill
