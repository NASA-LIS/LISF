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
! !ROUTINE: noah36_readrst
! \label{noah36_readrst}
!
! !REVISION HISTORY:
!  05 Sep 2001: Brian Cosgrove; Modified code to use Dag Lohmann's NOAA
!               initial conditions if necessary.  This is controlled with
!               local variable NOAAIC.  Normally set to 0 in this subroutine
!               but set to 1 if want to use Dag's NOAA IC's.  Changed output
!               directory structure, and commented out if-then check so that
!               directory is always made.
!  28 Apr 2002: Kristi Arsenault; Added Noah2.5 LSM into LDAS
!  28 May 2002: Kristi Arsenault; For STARTCODE=4, corrected SNEQV values  
!                and put SMC, SH2O, STC limit for GDAS and GEOS forcing.
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  10 Jun 2012: Sujay Kumar, added support for netcdf formats
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
!
! !INTERFACE:
subroutine noah36_readrst()
! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber,&
       LIS_verify
  use module_sfcdif_wrf_36, only: MYJSFCINIT
  use noah36_lsmMod
  use LIS_tbotAdjustMod, only: LIS_readTmnUpdateRestart
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
!
! !DESCRIPTION:
!  This program reads restart files for Noah-3.6.  This
!  includes all relevant water/energy storages and tile information. 
!  The following is the list of variables specified in the Noah-3.6 
!  restart file: 
!  \begin{verbatim}
!   nc,nr,ntiles    - grid and tile space dimensions 
!   t1           - Noah-3.6 skin temperature
!   cmc          - Noah-3.6 canopy moisture storage
!   snowh        - Noah-3.6 snow depth
!   sneqv        - Noah-3.6 snow water equivalent
!   stc          - Noah-3.6 soil temperature (for each layer)
!   smc          - Noah-3.6 soil moisture (for each layer) 
!   sh2o         - Noah-3.6 liquid only soil moisture (for each layer)
!   ch           - Noah-3.6 heat/moisture exchange coefficient
!   cm           - Noah-3.6 momentum exchange coefficient
!  \end{verbatim}
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!  reads a variable from the restart file
! \item[noah36\_coldstart](\ref{noah36_coldstart}) \newline
!   initializes the Noah-3.6 state variables
! \end{description}
!EOP
  implicit none      

  integer           :: t,l
  integer           :: nc,nr,npatch
  integer           :: n
  integer           :: ftn
  integer           :: status
  real, allocatable :: tmptilen(:)
  logical           :: file_exists
  character*20      :: wformat

  wformat="netcdf"

  do n=1,LIS_rc%nnest

     if(LIS_rc%startcode.eq."coldstart") then  ! coldstart
        call noah36_coldstart(LIS_rc%lsm_index)
        ! Set initial q1 to a negative value - D. Mocko
        ! Will be set to q2 for the first timestep in noah36_main
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           noah36_struc(n)%noah(t)%q1 = -0.5
        enddo
     elseif(LIS_rc%startcode.eq."restart") then  ! start from restart
        allocate(tmptilen(LIS_rc%npatch(n,LIS_rc%lsm_index)))
!-------------------------------------------------------------------------
! Read Active Archive File
!-------------------------------------------------------------------------
        inquire(file=noah36_struc(n)%rfile,exist=file_exists) 

        if(.not.file_exists) then 
           write(LIS_logunit,*) '[ERR] Noah-3.6 restart file ',              &
                noah36_struc(n)%rfile,' does not exist '
           write(LIS_logunit,*) '[ERR] Program stopping ...'
           call LIS_endrun()
           
        endif
        write(LIS_logunit,*)                        &
             '[INFO] Noah-3.6 restart file used: ',noah36_struc(n)%rfile
        
        if(wformat.eq."binary") then 
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=noah36_struc(n)%rfile,form='unformatted')        
           read(ftn) nc,nr,npatch  !time, veg class, no. tiles
           
!------------------------------------------------------------------------
!   Check for Grid Space Conflict 
!------------------------------------------------------------------------
           if(nc.ne.LIS_rc%gnc(n) .or. nr.ne.LIS_rc%gnr(n))then
              write(LIS_logunit,*) '[ERR]',noah36_struc(n)%rfile,              &
                   'grid space mismatch - Noah-3.6 halted'
              call LIS_endrun
           endif
!------------------------------------------------------------------------
! Transfer Restart tile space to LIS tile space
!------------------------------------------------------------------------
           if(npatch.ne.LIS_rc%glbnpatch_red(n,LIS_rc%lsm_index))then           
              write(LIS_logunit,*) '[ERR] restart tile space mismatch, halting..'
              call LIS_endrun
           endif
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status=nf90_open(path=noah36_struc(n)%rfile,&
                mode=NF90_NOWRITE,ncid=ftn)
           call LIS_verify(status,&
                'Error opening file '//noah36_struc(n)%rfile)
#endif
        endif

        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,noah36_struc(n)%noah%t1, &
             varname="T1",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,noah36_struc(n)%noah%cmc,&
             varname="CMC", wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,noah36_struc(n)%noah%snowh,&
             varname="SNOWH",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,noah36_struc(n)%noah%sneqv,&
             varname="SNEQV",wformat=wformat)
        
        do l=1,noah36_struc(n)%nslay
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen, varname="STC",&
                dim=l,vlevels = noah36_struc(n)%nslay, &
                wformat=wformat)
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              noah36_struc(n)%noah(t)%stc(l) = tmptilen(t)
           enddo
        enddo
        
        do l=1,noah36_struc(n)%nslay
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen, varname="SMC",&
                dim=l,vlevels = noah36_struc(n)%nslay, &
                wformat=wformat)
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              noah36_struc(n)%noah(t)%smc(l) = tmptilen(t)
           enddo
        enddo
        
        do l=1,noah36_struc(n)%nslay
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen, varname="SH2O",&
                dim=l,vlevels = noah36_struc(n)%nslay, &
                wformat=wformat)
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              if(tmptilen(t).lt.0.0001)  tmptilen(t) = 0.0001
              noah36_struc(n)%noah(t)%sh2o(l) = tmptilen(t)
           enddo
        enddo
        
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index, noah36_struc(n)%noah%ch,&
             varname="CH",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index, noah36_struc(n)%noah%cm,&
             varname="CM",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index, noah36_struc(n)%noah%z0,&
             varname="Z0",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index, noah36_struc(n)%noah%emiss,&
             varname="EMISS",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index, noah36_struc(n)%noah%albedo,&
             varname="ALBEDO",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index, noah36_struc(n)%noah%snotime1,&
             varname="SNOTIME1",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index, noah36_struc(n)%noah%q1,&
             varname="Q1",wformat=wformat)

        !Read in variables for dynamic deep soil temperature
        call LIS_readTmnUpdateRestart(n,ftn,wformat)
        
        if(wformat.eq."binary") then 
           call LIS_releaseUnitNumber(ftn)        
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in noah36_readrst')
#endif     
        endif
        deallocate(tmptilen)
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           noah36_struc(n)%noah(t)%z0_old = noah36_struc(n)%noah(t)%z0
        enddo
        
     endif
  enddo
  
  do n=1,LIS_rc%nnest
     if (noah36_struc(n)%sfcdifoption.eq.1) then
        call MYJSFCINIT(noah36_struc(n)%ztmax2, &
             noah36_struc(n)%dzeta2, &
             noah36_struc(n)%psih2,  &
             noah36_struc(n)%psim2)
     endif
  enddo
end subroutine noah36_readrst
