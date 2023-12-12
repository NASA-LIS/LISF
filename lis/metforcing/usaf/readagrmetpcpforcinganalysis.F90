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
!
! !ROUTINE: readagrmetpcpforcinganalysis
! \label{readagrmetpcpforcinganalysis}
!
! !REVISION HISTORY:
! 29Jul2005; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readagrmetpcpforcinganalysis(n,findex, order)
! !USES:
  use LIS_coreMod, only         : LIS_rc,LIS_domain, LIS_masterproc
  use LIS_timeMgrMod, only      : LIS_julhr_date, LIS_get_julhr
  use LIS_logMod, only          : LIS_logunit, LIS_endrun
  use LIS_fileIOMod, only       : LIS_putget
  use AGRMET_forcingMod, only : agrmet_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: order
!
! !DESCRIPTION:
!  This routine calls the previously generated 
!  AGRMET precipitation analysis for the 
!  current time.
!
!  The arguments and variables are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!  \item[n]
!    index of the nest
!  \item[hemi]
!    index of hemisphere loops
!  \item[c,r,i,t]
!    looping and indexing variables
!  \item[use\_twelve]
!   flag to use the 12-hourly precip amounts or the 6 hourly 
!   amounts 
!  \item[julbeg]
!    starting julian hour
!  \item[julend]
!    ending julian hour
!  \item[varfield]
!    interpolated variable
!  \item[yr1,mo1,da1,hr1,mn1,ss1]
!    time/date specific variables 
!  \item[alert\_number]
!    alert message number
!  \item[gi]
!    input merged precip field
!  \item[quad9r]
!    undefined value
!  \item[ip]
!    interpolation option
!  \item[c,r,k]
!   looping and indexing variables
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[get\_agrpcp\_readtime](\ref{get_agrpcp_readtime}) \newline
!    computes the reading time for the precip processing.
!  \item[julhr\_date] (\ref{LIS_julhr_date}) \newline
!    converts julian hour to a date format
!  \item[interp\_agrmetvar](\ref{interp_agrmetvar}) \newline
!    spatial interpolation of an AGRMET variable to LIS grid
!  \item[AGRMET\_fillgaps](\ref{AGRMET_fillgaps}) \newline
!    fills the gaps in the interpolated field due to mismatches in 
!    LIS and AGRMET masks
!  \item[LIS\_putget](\ref{LIS_putget}) \newline
!    reads the data from the previously generated analysis
!  \end{description}
!EOP
  real                  :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                  :: gi(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
  integer               :: julhr
  character*10          :: date10_03
  integer               :: yr1,mo1,da1,hr1
  integer               :: hemi
  character*4           :: fyr
  character*2           :: fmo,fda
  character*3           :: chemi(2)
  character*100         :: ifil
  character*30          :: routine_name
  real                  :: udef
  integer               :: c,r
  logical               :: exists
  integer               :: ip

  data chemi  / '_nh','_sh' /
  data routine_name     / 'readagrmetpcpforcinganalysis' /

!read every 3 hours, the nearest next analysis. 
     call get_agrpcp_readtime(LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn,julhr)
     
     do hemi=1,2
        call LIS_julhr_date(julhr,yr1,mo1,da1,hr1)
        call AGRMET_julhr_date10( julhr, date10_03 )
        
        write(unit=fyr,fmt='(I4.4)') yr1
        write(unit=fmo,fmt='(I2.2)') mo1
        write(unit=fda,fmt='(I2.2)') da1
        
        if(agrmet_struc(n)%use_timestamp.eq.1) then 
           ifil =  trim(agrmet_struc(n)%agrmetdir)//'/'&
                //trim(fyr)//trim(fmo)// &
                trim(fda)//'/'//trim(agrmet_struc(n)%mrgpcpdir)//&
                '/premrgp'//chemi(hemi) &
                //'.03hr.'//date10_03
        else
           ifil =  trim(agrmet_struc(n)%agrmetdir)//'/'&
                //trim(agrmet_struc(n)%mrgpcpdir)//&
                '/premrgp'//chemi(hemi) &
                //'.03hr.'//date10_03
        endif

        write(LIS_logunit,*)'[INFO] READING precip ',ifil
        inquire(file=ifil,exist=exists)
        if(exists) then      
           call LIS_putget(gi(hemi,:,:), 'r', ifil, &
                routine_name, agrmet_struc(n)%imax, agrmet_struc(n)%jmax )
        else
           write(LIS_logunit,*) '[ERR] premrg file does not exist'
           write(LIS_logunit,*) ifil
           call LIS_endrun()
        endif
     enddo
     
     ip = 1
     udef = 9999.0
     varfield = -9999.0
     call interp_agrmetvar(n,ip,gi,udef,varfield,agrmet_struc(n)%imax, agrmet_struc(n)%jmax) !,agrmet_struc(n)%shemi,agrmet_struc(n)%nhemi)
     call AGRMET_fillgaps(n,ip,varfield)
     
     do c =1, LIS_rc%lnc(n)
        do r = 1,LIS_rc%lnr(n)
           if (LIS_domain(n)%gindex(c,r).ne. -1) then 
              agrmet_struc(n)%metdata2(8,LIS_domain(n)%gindex(c,r)) = varfield(c,r)
           else
              agrmet_struc(n)%metdata2(8,LIS_domain(n)%gindex(c,r)) = -9999.0
           endif
        end do
     end do
end subroutine readagrmetpcpforcinganalysis


!BOP
! 
! !ROUTINE: get_agrpcp_readtime
! \label{get_agrpcp_readtime}
! 
! !INTERFACE:
subroutine get_agrpcp_readtime(yr,mo,da,hr,mn,julhr)        
! !USES: 
  use LIS_timeMgrMod, only : LIS_tick, LIS_get_julhr

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: yr
  integer, intent(in) :: mo
  integer, intent(in) :: da
  integer, intent(in) :: hr
  integer, intent(in) :: mn
  integer, intent(inout) :: julhr
! 
! !DESCRIPTION: 
!  This routine gets the julian hour to read the AGRMET
!  precip from, based on the current input time. 
! 
!  The arguments are:
!  \begin{description}
!  \item[yr]
!    the current year
!  \item[mo]
!    the current month
!  \item[da]
!    the current day
!  \item[hr]
!    the current hour
!  \item[julbeg]
!    output AGRMET precip reading time
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick] (\ref{LIS_tick}) \newline
!    computes previous 6 hour time. 
!  \item[get\_julhr] (\ref{LIS_get_julhr}) \newline
!    converts the date to a julian hour
!  \end{description}
!EOP
  
  integer          :: yr1,mo1,da1,hr1,mn1,ss1
  real             :: curr_hr
  real*8           :: time1
  integer          :: doy1
  real             :: gmt1,ts1

  yr1 = yr
  mo1 = mo
  da1 = da
  mn1 = 0 
  ss1 = 0 

  curr_hr = float(hr)+float(mn)/60.0

  if(curr_hr.eq.0.or.curr_hr.eq.3.or.curr_hr.eq.6&
       .or.curr_hr.eq.9.or.curr_hr.eq.12.&
       .or.curr_hr.eq.15.or.curr_hr.eq.18.or.curr_hr.eq.21) then 
     yr1 = yr
     mo1 = mo
     da1 = da
     hr1 = hr
  else
    yr1 = yr  !next assimilation/forecast hour
    mo1 = mo
    da1 = da
    hr1 = 3*((hr/3))
    mn1 = 0
    ss1 = 0
    ts1 = 3*60*60 
    call LIS_tick( time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
 endif
 call LIS_get_julhr(yr1,mo1,da1,hr1,mn1,ss1,julhr)

end subroutine get_agrpcp_readtime
