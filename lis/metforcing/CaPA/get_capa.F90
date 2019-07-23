!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_capa
! \label{get_capa}
!
! !REVISION HISTORY:
! 23 Nov 2009: Sujay Kumar, Initial Specification
! 
! !INTERFACE:
subroutine get_capa(n, findex)
! !USES:
  use LIS_coreMod,     only : LIS_rc
  use LIS_timeMgrMod,  only : LIS_tick, LIS_get_nstep
  use LIS_logMod,      only : LIS_logunit
  use capa_forcingMod, only : capa_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 6-hrly, CAPA forcing. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 6 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the CAPA data times
!  \item[capafile](\ref{capafile}) \newline
!    Puts together appropriate file name for 6 hour intervals
!  \item[read\_capa](\ref{read_capa}) \newline
!    Interpolates CAPA data to LIS grid
!  \end{description}
!EOP
   

  integer :: ferror_capa ! Error flags for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy5, yr5, mo5, da5, hr5, mn5, ss5
  integer :: endtime_capa   
  real    :: ts1, ts5
  real*8  :: ctime,ftime_capa
  real    :: gmt1,gmt5       
  character*80 :: name 

!------------------------------------------------------------------------
! Determine required observed precip data times 
! (current, accumulation end time)
! Model current time
!------------------------------------------------------------------------
  yr1 = LIS_rc%yr  !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick( ctime, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )   
!------------------------------------------------------------------------ 
! CAPA product end time
!------------------------------------------------------------------------
  yr5 = LIS_rc%yr  !end accumulation time data
  mo5 = LIS_rc%mo
  da5 = LIS_rc%da
  hr5 = 6*(LIS_rc%hr/6)
  mn5 = 0
  ss5 = 0
  ts5 = 6*60*60
  call LIS_tick( ftime_capa, doy5, gmt5, yr5, mo5, da5, hr5, mn5, ss5, ts5 )

!------------------------------------------------------------------------
! Ensure that data is found during first time step
!------------------------------------------------------------------------
  if ( LIS_get_nstep(LIS_rc,n).eq. 1 & 
       .or.LIS_rc%rstflag(n) .eq. 1) then 
     endtime_capa = 1
     LIS_rc%rstflag(n) = 0
  endif
!------------------------------------------------------------------------
! Check for and get CAPA CPC Precipitation data
!------------------------------------------------------------------------

   name=''
   if ( ctime > capa_struc(n)%capatime ) endtime_capa = 1
   if ( endtime_capa == 1 ) then  !get new time2 data
        ferror_capa = 0
        call capafile( name, capa_struc(n)%capadir, yr5, mo5, da5, hr5 )
        write(LIS_logunit,*) 'Getting new CAPA CPC precip data',name
        call read_capa( n, findex, name, ferror_capa, hr5 )
        capa_struc(n)%capatime = ftime_capa
   endif  !need new time2

 return
end subroutine get_capa

!BOP
! !ROUTINE: capafile
! \label{capafile}
!
! !INTERFACE:
subroutine capafile( name, capadir, yr, mo, da, hr)

  implicit none
! !ARGUMENTS: 
  character(len=*)   :: name
  character(len=*)   :: capadir
  integer            :: yr, mo, da, hr
! !DESCRIPTION:
!   This subroutine puts together CAPA file name for 
!   6 hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[capadir]
!    Name of the CAPA directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[name]
!   name of the timestamped CAPA file
!  \end{description}
!
!EOP
  character(len=6)      :: cdate1
  character(len=10)      :: cdate2
 
  write(unit=cdate1, fmt='(i4.4,I2.2)') yr, mo
  write(unit=cdate2, fmt='(i4.4,3i2.2)') yr, mo, da, hr

  name = trim(capadir)//'/'//cdate1//'/pr_an_ps_'//cdate2//'.grb'

end subroutine capafile
