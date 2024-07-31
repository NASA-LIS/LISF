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
! !ROUTINE: mos_readrestart
!  \label{mos_readrestart}
!
! !REVISION HISTORY:
!  19 Jan 2001: Brian Cosgrove; Initial Code
!  12 Dec. 2002: Brian Cosgrove; Fixed usage of Wiltpoint variable.  Before,
!               Wiltpoint1 and Wiltpoint2 were used in calculation of
!               root zone soil moisture availability...now, only Wiltpoint2
!               is used since wiltpoint1 is not the correct wilting point
!               needed for the calculation.
!  23 Jan 2003: Urszula Jambor; Switch index order of GEOS forcing
!               array.  Snow is 12, soil wetness is 13.
!  25 Sep 2007: Sujay Kumar, Upgraded for LIS 5.0
!  10 Jun 2012: Sujay Kumar, added support for netcdf formats
! 
! !INTERFACE:
subroutine mos_readrestart
! !USES:
  use LIS_coreMod, only : LIS_rc
  use mos_lsmMod      
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for Mosaic.  This
!  includes all relevant water/energy storages and tile information. 
!  The following is the list of variables specified in the Mosaic 
!  restart file: 
!  \begin{verbatim}
!   nc,nr,ntiles    - grid and tile space dimensions 
!   ct           - mosaic canopy temperature
!   qa           - mosaic canopy humidity
!   ics          - mosaic interception canopy storage
!   snow         - mosaic snow depth
!   sot          - mosaic deep soil temperature
!   soWet        - mosaic soil wetness
!  \end{verbatim}
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!  reads a variable from the restart file
! \end{description}
!EOP
  implicit none      
  integer :: l,t 
  integer :: nc,nr,ntiles 
  integer :: n
  integer :: ftn
  logical :: file_exists
  integer :: status
  character*20 :: wformat
  real, allocatable :: tmptilen(:)
  wformat = "netcdf"
!-------------------------------------------------------------------------
! Read Active Archive File
!-------------------------------------------------------------------------
  if(LIS_rc%startcode.eq."restart")then
     do n=1,LIS_rc%nnest
        
        inquire(file=mos_struc(n)%mos_rfile,exist=file_exists)

        if(.not.file_exists) then 
            write(LIS_logunit,*) 'Mosaic restart file ',               &
                mos_struc(n)%mos_rfile,' does not exist '
           write(LIS_logunit,*) 'Program stopping ...'
           call LIS_endrun()
           
        endif
        write(LIS_logunit,*)                        &
             'Mosaic restart file used: ',mos_struc(n)%mos_rfile
        
        if(wformat.eq."binary") then 
           ftn = LIS_getNextUnitNumber()
           
           open(ftn,file=mos_struc(n)%mos_rfile,form='unformatted')
           write(LIS_logunit,*)'mosaic restart file used: ',mos_struc(n)%mos_rfile
           read(ftn) nc,nr,ntiles  !veg class, no. columns, no. rows, no. tiles
!------------------------------------------------------------------------
!   Check for Grid Space Conflict 
!------------------------------------------------------------------------
           if(nc.ne.LIS_rc%gnc(n) .or. nr.ne.LIS_rc%gnr(n))then
              write(LIS_logunit,*)mos_struc(n)%mos_rfile,'grid space mismatch - mosaic halted'
              call LIS_endrun
           endif
           if(ntiles.ne.LIS_rc%glbnpatch(n,LIS_rc%lsm_index))then           
              write(LIS_logunit,*)'restart tile space mismatch, halting..'
              call LIS_endrun
           endif
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status=nf90_open(path=mos_struc(n)%mos_rfile,&
                mode=NF90_NOWRITE,ncid=ftn)
           call LIS_verify(status,&
                'Error opening file '//mos_struc(n)%mos_rfile)
#endif
        endif
!------------------------------------------------------------------------
! Transfer Restart tile space to LIS tile space
!------------------------------------------------------------------------
        
        call LIS_readvar_restart(ftn,n,mos_struc(n)%mos%ct,&
             varname="CT",wformat=wformat)
        call LIS_readvar_restart(ftn,n,mos_struc(n)%mos%qa,&
             varname="QA",wformat=wformat)
        call LIS_readvar_restart(ftn,n,mos_struc(n)%mos%ics,&
             varname="ICS",wformat=wformat)             
        call LIS_readvar_restart(ftn,n,mos_struc(n)%mos%snow,&
             varname="SNOW",wformat=wformat)                          
        call LIS_readvar_restart(ftn,n,mos_struc(n)%mos%SoT,&
             varname="SoT",wformat=wformat)  
        allocate(tmptilen(LIS_rc%npatch(n,LIS_rc%lsm_index)))
        do l=1,3
           call LIS_readvar_restart(ftn,n,tmptilen,varname="SoWET",&
                dim=l,vlevels=3,wformat=wformat)
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              mos_struc(n)%mos(t)%SoWET(l) = tmptilen(t)
           enddo
        enddo
        deallocate(tmptilen)
        if(wformat.eq."binary") then 
           close(ftn)        
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in mos_readrst')
#endif     
        endif
     enddo
  endif
end subroutine mos_readrestart






