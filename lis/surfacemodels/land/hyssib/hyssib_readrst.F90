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
! !ROUTINE: hyssib_readrst
! \label{hyssib_readrst}
!
! !REVISION HISTORY:
! 21 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
! 25 Aug 2007: Chuck Alonge, Updates for LIS 5.0 Compliance
!
! !INTERFACE:
subroutine hyssib_readrst
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod
  use hyssib_lsmMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
!
! !DESCRIPTION:
!  This program reads restart files for HY-SSiB.  This includes all
!  relevant water/energy storages, tile information, and time information.
!  The following is the list of variables specified in the Hyssib 
!  restart file: 
!  \begin{verbatim}
!   nc,nr,ntiles    - grid and tile space dimensions 
!   tc           - hyssib canopy temperature
!   tg           - hyssib ground temperature
!   tsn          - hyssib snow on ground temperature
!   td           - hyssib deep soil temperature
!   www(3)       - hyssib soil wetness (for each layer)
!   capac(2)     - hyssib liq equiv. water&snow on canopy/ground 
!   snow(2)      - hyssib snow storage (swe) on snow/ground
!   sgfg         - hyssib density of bulk snow layer
!   sdens        - hyssib momentum exchange coefficient
!  \end{verbatim}
!
!  The routines invoked are: 
! \begin{description}
! \item[drv\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!  reads a variable from the restart file
! \end{description}
!EOP
  implicit none

  INTEGER           :: L,N   ! Loop counters
  INTEGER           :: NC,NR,npatch
  integer           :: ftn
  logical           :: file_exists
  character*20      :: wformat
  integer           :: status
  real, allocatable :: tmptile(:)
!=== End Variable Definition ===========================================
  wformat = "netcdf"

  do n=1,LIS_rc%nnest
!-----------------------------------------------------------------------
! Read Active Archive File
!-----------------------------------------------------------------------
     if (LIS_rc%startcode.eq."restart") then
        allocate(tmptile(LIS_rc%npatch(n,LIS_rc%lsm_index)))
        inquire(file=hyssib_struc(n)%rfile,exist=file_exists) 

        if(.not.file_exists) then 
           write(LIS_logunit,*) 'HySSIB restart file ',               &
                trim(hyssib_struc(n)%rfile),' does not exist '
           write(LIS_logunit,*) 'Program stopping ...'
           call LIS_endrun()
           
        endif
        write(LIS_logunit,*)                        &
             'HySSIB restart file used: ',hyssib_struc(n)%rfile
        
        if(wformat.eq."binary") then 
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=hyssib_struc(n)%rfile,form='unformatted')        
           read(ftn) nc,nr,npatch  !time, veg class, no. tiles
           
!------------------------------------------------------------------------
!   Check for Grid Space Conflict 
!------------------------------------------------------------------------
           if(nc.ne.LIS_rc%gnc(n) .or. nr.ne.LIS_rc%gnr(n))then
              write(LIS_logunit,*) hyssib_struc(n)%rfile,                 &
                   'grid space mismatch - Noah3.3 halted'
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
           status=nf90_open(path=trim(hyssib_struc(n)%rfile),&
                mode=NF90_NOWRITE,ncid=ftn)
           call LIS_verify(status,&
                'Error opening file '//trim(hyssib_struc(n)%rfile))
#endif
        endif
        
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,hyssib_struc(n)%hyssib%tc,&
             varname="TC",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,hyssib_struc(n)%hyssib%tg,&
             varname="TG",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,hyssib_struc(n)%hyssib%tsn,&
             varname="TSN",wformat=wformat) 
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,hyssib_struc(n)%hyssib%td,&
             varname="TD",wformat=wformat)

        do l=1,3
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptile,&
             varname="WWW",dim=l,vlevels =3,wformat=wformat)
            hyssib_struc(n)%hyssib%www(l) = tmptile
        enddo

        do l=1,2
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptile,&
             varname="CAPAC",dim=l,vlevels =2,wformat=wformat)
           hyssib_struc(n)%hyssib%capac(l) = tmptile
        enddo

        do l=1,2
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,tmptile,&
             varname="SNOW",dim=l,vlevels =2,wformat=wformat)
           hyssib_struc(n)%hyssib%snow(l) = tmptile
        enddo

        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index, hyssib_struc(n)%hyssib%sgfg,&
             varname="SGFG",wformat=wformat)
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index, hyssib_struc(n)%hyssib%sdens,&
             varname="SDENS",wformat=wformat)

        if(wformat.eq."binary") then 
           call LIS_releaseUnitNumber(ftn)        
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in hyssib_readrst')
#endif     
        endif
        deallocate(tmptile)
     endif
  enddo
end subroutine hyssib_readrst
