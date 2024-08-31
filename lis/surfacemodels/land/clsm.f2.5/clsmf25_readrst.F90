!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: clsmf25_readrst
! \label{clsmf25_readrst}
!
! !REVISION HISTORY:
!  16 Dec 2005: Sujay Kumar, Initial Specification
! 
! !INTERFACE:
subroutine clsmf25_readrst
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_timeMgrMod
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use clsmf25_constants
  use clsmf25_lsmMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for Catchment.  This
!  includes all relevant water/energy storages, tile information.
!
!  The routines invoked are: 
! \begin{description}
! \item[drv\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!  reads a variable from the restart file
! \end{description}

!EOP
  implicit none      
  real :: ts,mc
  integer :: l,t, k,n
  integer :: ftn,status
  integer :: nc,nr,ntiles
  real, allocatable :: temp(:)
  logical           :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: filen
  character*20      :: wformat
  integer           :: yr,mo,da,hr,mn,ss,doy
  real*8            :: time
  real              :: gmt
  logical           :: read_restart

!=== End Variable Definition =============================================
  wformat = "netcdf"
!-------------------------------------------------------------------------
! Read Active Archive File
!-------------------------------------------------------------------------
  do n=1,LIS_rc%nnest

     read_restart = .false. 

     if(LIS_rc%startcode.eq."restart") then 
        read_restart = .true. 
     endif

     if(LIS_rc%runmode.eq."ensemble smoother") then 
        if(LIS_rc%iterationId(n).gt.1) then 
           read_restart = .true. 

!           if(clsmf25_struc(n)%rstInterval.ne.2592000) then 
!              write(LIS_logunit,*) 'restart interval must be set to a month' 
!              write(LIS_logunit,*) 'when running in the ensemble smoother mode'
!              call LIS_endrun()
!           endif

           if(clsmf25_struc(n)%rstInterval.eq.2592000) then 
           !create the restart filename based on the timewindow start time
              call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                   dd=da,calendar=LIS_calendar,rc=status)
              hr = 0 
              mn = 0 
              ss = 0 
              call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,(-1)*LIS_rc%ts)
           else
              call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                   dd=da,calendar=LIS_calendar,rc=status)
              hr = 0 
              mn = 0 
              ss = 0 
           endif

           call LIS_create_restart_filename(n,filen,'SURFACEMODEL','CLSMF25', &
                yr,mo,da,hr,mn,ss, wformat=wformat)
           clsmf25_struc(n)%rfile = filen
        endif
     endif

     if(read_restart) then 
        allocate(temp(LIS_rc%npatch(n,LIS_rc%lsm_index)))

        inquire(file=clsmf25_struc(n)%rfile,exist=file_exists) 

        if(.not.file_exists) then 
           write(LIS_logunit,*) '[ERR] CLSM F2.5 restart file ',               &
                trim(clsmf25_struc(n)%rfile),' does not exist. '
           write(LIS_logunit,*) 'Program stopping ...'
           call LIS_endrun()
           
        endif
        write(LIS_logunit,*)                        &
             '[INFO] CLSM F2.5 restart file used: ',clsmf25_struc(n)%rfile
        
        if(wformat.eq."binary") then 
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=clsmf25_struc(n)%rfile,form='unformatted')        
           read(ftn) nc,nr,ntiles  !time, veg class, no. tiles

           if(nc.ne.LIS_rc%gnc(n).or.nr.ne.LIS_rc%gnr(n))then
              write(LIS_logunit,*) "[ERR] ",trim(clsmf25_struc(n)%rfile),&
                   "grid space mismatch - Catchment run halted."
              call LIS_endrun
           endif
!------------------------------------------------------------------------
! Transfer Restart tile space to LIS tile space
!------------------------------------------------------------------------
           if(ntiles.ne.LIS_rc%glbnpatch(n,LIS_rc%lsm_index))then           
              write(LIS_logunit,*)'[ERR] Restart tile space mismatch, halting ...',&
                   ntiles,LIS_rc%glbnpatch(n,LIS_rc%lsm_index)
              call LIS_endrun
           endif
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status=nf90_open(path=trim(clsmf25_struc(n)%rfile),&
                mode=NF90_NOWRITE,ncid=ftn)
           call LIS_verify(status,&
                'Error opening file '//trim(clsmf25_struc(n)%rfile))
#endif
        endif

        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="TC1",wformat=wformat)
        clsmf25_struc(n)%cat_progn%tc1 = temp
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="TC2",wformat=wformat)
        clsmf25_struc(n)%cat_progn%tc2 = temp
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="TC4",wformat=wformat)
        clsmf25_struc(n)%cat_progn%tc4 = temp
        
        ! skip tsnow -- this is needed to read Sarith's restart file.
        ! This variable is not necessary for a restart.  It is a forcing
        ! variable -- snow fall rate.
        ! LIS' restart file does not contain it.
        !call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp)

        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="QA1",wformat=wformat)
        clsmf25_struc(n)%cat_progn%qa1 = temp
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="QA2",wformat=wformat)
        clsmf25_struc(n)%cat_progn%qa2 = temp
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="QA4",wformat=wformat)
        clsmf25_struc(n)%cat_progn%qa4 = temp
        
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="CAPAC",wformat=wformat)
        clsmf25_struc(n)%cat_progn%capac = temp

        do l=1,N_gt
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="GHT",dim=l,vlevels=N_gt,wformat=wformat)
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              clsmf25_struc(n)%cat_progn(t)%ght(l) = temp(t)
           enddo
        enddo

        do l=1,N_snow
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="WESN",dim=l,vlevels=N_snow,wformat=wformat)
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              clsmf25_struc(n)%cat_progn(t)%wesn(l) = temp(t)
           enddo
        enddo
        
        do l=1,N_snow
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="HTSN",dim=l,vlevels=N_snow,wformat=wformat)
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              clsmf25_struc(n)%cat_progn(t)%htsn(l) = temp(t)
           enddo
        enddo
        
        do l=1,N_snow
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="SNDZ",dim=l,vlevels=N_snow,wformat=wformat)
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              clsmf25_struc(n)%cat_progn(t)%sndz(l) = temp(t)
           enddo
        enddo

        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="CATDEF",wformat=wformat)
        clsmf25_struc(n)%cat_progn%catdef = temp
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="RZEXC",wformat=wformat)
        clsmf25_struc(n)%cat_progn%rzexc = temp
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,temp,&
             varname="SRFEXC",wformat=wformat)
        clsmf25_struc(n)%cat_progn%srfexc = temp

        if(wformat.eq."binary") then 
           call LIS_releaseUnitNumber(ftn)        
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in clsmf25_readrst')
#endif     
        endif
        deallocate(temp)
     elseif(LIS_rc%startcode.eq."coldstart") then 
        
        clsmf25_struc(n)%cat_progn%tc1 = clsmf25_struc(n)%initST
        clsmf25_struc(n)%cat_progn%tc2 = clsmf25_struc(n)%initST
        clsmf25_struc(n)%cat_progn%tc4 = clsmf25_struc(n)%initST
        
        clsmf25_struc(n)%cat_progn%qa1 = 0.002
        clsmf25_struc(n)%cat_progn%qa2 = 0.002
        clsmf25_struc(n)%cat_progn%qa4 = 0.002
        
        clsmf25_struc(n)%cat_progn%capac = 0
        
        ts = clsmf25_struc(n)%initST
        mc = clsmf25_struc(n)%initSM
        do k=1,N_gt
           clsmf25_struc(n)%cat_progn%ght(k) =            &
                (ts-273.16)*                                  &
                ( 2.4e6*0.55*clsmf25_struc(n)%cat_param%dzgt(k)   &
                +   0.5*0.45*clsmf25_struc(n)%cat_param%dzgt(k)*4.185e6)
        enddo
        
        do k=1,N_snow
           
           clsmf25_struc(n)%cat_progn%wesn(k) = 0.
           clsmf25_struc(n)%cat_progn%htsn(k) = 0.
           clsmf25_struc(n)%cat_progn%sndz(k) = 0.
           
        end do

        clsmf25_struc(n)%cat_progn%catdef = &
             max( (clsmf25_struc(n)%cat_param%poros-mc), 0.05) * &
             clsmf25_struc(n)%cat_param%dzpr
         
        clsmf25_struc(n)%cat_progn%rzexc  = 0.
        clsmf25_struc(n)%cat_progn%srfexc = 0.
        
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           if (clsmf25_struc(n)%cat_progn(t)%tc1.le.0.0) then
              write (LIS_logunit,*) 'Negative TC1 at tile ',t
              call LIS_endrun()
           endif
           if (clsmf25_struc(n)%cat_progn(t)%tc2.le.0.0) then
              write (LIS_logunit,*) 'Negative negative TC2 at tile ',t
              call LIS_endrun()
           endif
           if (clsmf25_struc(n)%cat_progn(t)%tc4.le.0.0) then
              write (LIS_logunit,*) 'Negative negative TC4 at tile ',t
              call LIS_endrun()
           endif
        enddo
     endif
!     write(LIS_logunit,*) 'restart ',clsmf25_struc(n)%cat_progn(327)%catdef
  enddo
end subroutine clsmf25_readrst
