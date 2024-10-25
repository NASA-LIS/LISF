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
! !ROUTINE: cable_readrst
! \label{cable_readrst}
!
! !REVISION HISTORY:
!  27 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!  13 Sep 2011: Claire Carouge (ccc), CABLE LSM improvements
!  10 Jun 2012: Sujay Kumar, added support for netcdf formats
!  23 May 2013: David Mocko, latest CABLE v1.4b version for LIS6.2
!
! !INTERFACE:
subroutine cable_readrst
! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber,&
       LIS_verify
  use cable_dimensions,   only : ms,msn,ncp,ncs
  use cable_lsmMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
!
! !DESCRIPTION:
!  This subroutine initializes the state variables in the CABLE LSM.
!  It will either read the restart files if one is specified, or it
!  will initialize the variables by calling the coldstart routine.
!  All relevant water/energy storages and tile info is initialized.
!
!
!  The following is the list of variables specified in the
!  CABLE restart file:
!   \begin{verbatim}
!   nc,nr,ntiles    - grid and tile space dimensions
!   cansto       - CABLE canopy water storage (kg m-2)
!   rtsoil       - CABLE turbulent resistance for soil (units?)
!   ssdnn        - CABLE overall snow density (kg/m3)
!   snowd        - CABLE snow liquid water equivalent depth (kg m-2)
!   osnowd       - CABLE snowd from previous timestep (kg m-2)
!   snage        - CABLE snow age (seconds?)
!   isflag       - CABLE snow layer scheme flag (0 = no or little snow, 1 = snow)
!   wbice        - CABLE soil ice - "ms" soil layers (kg m-2)
!   tggsn        - CABLE snow temperature - "msn" snow layers (K)
!   ssdn         - CABLE snow density - "msn" snow layers (kg/m3)
!   smass        - CABLE snow mass - "msn" snow layers (kg m-2)
!   albsoilsn    - CABLE soil+snow albedo - "msn" snow layers (0-1)
!   wb           - CABLE soil moisture - "ms" soil layers (kg m-2)
!   tgg          - CABLE soil temperature - "ms" soil layers (K)
!   cplant       - CABLE plant carbon - "ncp" vegetation carbon stores (g C/m2)
!   csoil        - CABLE soil carbon - "ncs" soil carbon stores (g C/m2)
!   gammzz       - CABLE heat capacity for each soil layer.
!   ktau         - CABLE timestep counter
!   \end{verbatim}
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!    Reads a variable from the restart file
!   \item[cable\_coldstart](\ref{cable_coldstart}) \newline
!    Initializes the CABLE state variables
!  \end{description}
!EOP
  implicit none

  integer :: n,t,l
  integer :: nc,nr,npatch
  real, allocatable :: tmptilen(:)
  logical      :: file_exists
  integer      :: ftn
  integer      :: status
  character*20 :: wformat

  wformat="netcdf"
  if(LIS_rc%startcode.eq."restart") then 
     do n=1,LIS_rc%nnest
        allocate(tmptilen(LIS_rc%npatch(n,LIS_rc%lsm_index)))
!-------------------------------------------------------------------------
! Read Active Archive File
!-------------------------------------------------------------------------
        inquire(file=cable_struc(n)%rfile,exist=file_exists) 

        if(.not.file_exists) then 
           write(LIS_logunit,*) 'CABLE restart file ',               &
                trim(cable_struc(n)%rfile),' does not exist '
           write(LIS_logunit,*) 'Program stopping ...'
           call LIS_endrun()
           
        endif
        write(LIS_logunit,*)                        &
             'CABLE restart file used: ',cable_struc(n)%rfile
        
        if(wformat.eq."binary") then 
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=cable_struc(n)%rfile,form='unformatted')        
           read(ftn) nc,nr,npatch  !time, veg class, no. tiles
           
!------------------------------------------------------------------------
!   Check for Grid Space Conflict 
!------------------------------------------------------------------------
           if(nc.ne.LIS_rc%gnc(n) .or. nr.ne.LIS_rc%gnr(n))then
              write(LIS_logunit,*) cable_struc(n)%rfile,                 &
                   'grid space mismatch - CABLE halted'
              call LIS_endrun
           endif
!------------------------------------------------------------------------
! Transfer Restart tile space to LIS tile space
!------------------------------------------------------------------------
           if(npatch.ne.LIS_rc%glbnpatch_red(n,LIS_rc%lsm_index))then           
              write(LIS_logunit,*) 'restart tile space mismatch, halting..'
              call LIS_endrun
           endif
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status=nf90_open(path=trim(cable_struc(n)%rfile),&
                mode=NF90_NOWRITE,ncid=ftn)
           call LIS_verify(status,&
                'Error opening file '//trim(cable_struc(n)%rfile))
#endif
        endif
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             cable_struc(n)%cable%cansto,varname="CANSTO",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             cable_struc(n)%cable%rtsoil,varname="RTSOIL",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             cable_struc(n)%cable%ssdnn,varname="SSDNN",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             cable_struc(n)%cable%snowd,varname="SNOWD",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             cable_struc(n)%cable%osnowd,varname="OSNOWD",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             cable_struc(n)%cable%snage,varname="SNAGE",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             cable_struc(n)%cable%isflag,varname="ISFLAG",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             cable_struc(n)%cable%wbtot,varname="WBTOT",wformat=wformat)
        do l = 1,ms
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
                varname="WBICE",dim=l,vlevels=ms,wformat=wformat)
           do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              cable_struc(n)%cable(t)%wbice(l) = tmptilen(t)
           enddo
        enddo
        do l = 1,msn
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
                varname="TGGSN",dim=l,vlevels=msn,wformat=wformat)
           do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              cable_struc(n)%cable(t)%tggsn(l) = tmptilen(t)
           enddo
        enddo
        do l = 1,msn
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
                varname="SSDN",dim=l,vlevels=msn,wformat=wformat)
           do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              cable_struc(n)%cable(t)%ssdn(l) = tmptilen(t)
           enddo
        enddo
        do l = 1,msn
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
                varname="SMASS",dim=l,vlevels=msn,wformat=wformat)
           do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              cable_struc(n)%cable(t)%smass(l) = tmptilen(t)
           enddo
        enddo
        do l = 1,msn
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
                varname="ALBSOILSN",dim=l,vlevels=msn,wformat=wformat)
           do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              cable_struc(n)%cable(t)%albsoilsn(l) = tmptilen(t)
           enddo
        enddo
        do l = 1,ms
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
                varname="WB",dim=l,vlevels=ms,wformat=wformat)
           do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              cable_struc(n)%cable(t)%wb(l) = tmptilen(t)
           enddo
        enddo
        do l = 1,ms
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
                varname="TGG",dim=l,vlevels=ms,wformat=wformat)
           do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              cable_struc(n)%cable(t)%tgg(l) = tmptilen(t)
           enddo
        enddo
        do l = 1,ncp
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
                varname="CPLANT",dim=l,vlevels=ncp,wformat=wformat)
           do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              cable_struc(n)%cable(t)%cplant(l) = tmptilen(t)
           enddo
        enddo
        do l = 1,ncs
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
                varname="CSOIL",dim=l,vlevels=ncs,wformat=wformat)
           do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              cable_struc(n)%cable(t)%csoil(l) = tmptilen(t)
           enddo
        enddo
        do l = 1,ms
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
                varname="GAMMZZ",dim=l,vlevels=ms,wformat=wformat)
           do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              cable_struc(n)%cable(t)%gammzz(l) = tmptilen(t)
           enddo
        enddo
        tmptilen = 0.0
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
             varname="KTAU",wformat=wformat)
        cable_struc(n)%ktau = tmptilen(1)

        if(wformat.eq."binary") then 
           call LIS_releaseUnitNumber(ftn)        
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in noah33_readrst')
#endif     
        endif
        deallocate(tmptilen)
     end do

     ! If no restart file is use, call coldstart to provide initial values
  elseif (LIS_rc%startcode.eq."coldstart") then
     call cable_coldstart
  endif

end subroutine cable_readrst
