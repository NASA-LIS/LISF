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
! !ROUTINE: mos_writerst
! \label{mos_writerst}
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
subroutine mos_writerst(n)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only  : LIS_isAlarmRinging
  use mos_lsmMod, only    : mos_struc
  use LIS_logMod
  use lis_fileIOMod, only : LIS_create_output_directory, &
       LIS_create_restart_filename
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none      
! !ARGUMENTS: 
  integer, intent(in) :: n 
!
! !DESCRIPTION:
!  This program writes restart files for Mosaic.  This
!  includes all relevant water/energy storage and tile information
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[mos\_dump\_restart](\ref{mos_dump_restart}) \newline
!   writes the Mosaic variables into the restart file
! \end{description}
!EOP

  character(len=LIS_CONST_PATH_LEN) :: filen
  logical       :: alarmCheck
  integer       :: status
  integer       :: ftn 
  character*20  :: wformat
  character*3         :: fnest

  write(fnest,'(i3.3)') n    

  alarmCheck = LIS_isAlarmRinging(LIS_rc, &
       "Mosaic restart alarm "//trim(fnest))
  wformat = "netcdf"
  if(alarmCheck .or.(LIS_rc%endtime ==1)) then 

     if(LIS_masterproc) then
        call LIS_create_output_directory('SURFACEMODEL')
        call LIS_create_restart_filename(n, filen,'SURFACEMODEL',&
             'MOS',wformat=wformat)
        if(wformat.eq."binary") then 
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=filen,status='unknown',form='unformatted')
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3)
           status = nf90_create(path=filen,cmode=nf90_clobber,&
                ncid = ftn)
           call LIS_verify(status,'Error in nf90_open in mos_writerst')
#endif
#if (defined USE_NETCDF4)
           status = nf90_create(path=trim(filen),cmode=nf90_hdf5,&
               ncid = ftn)
           call LIS_verify(status,'Error in nf90_open in mos_writerst')
#endif
           
        endif
     endif
     call mos_dump_restart(ftn,n,wformat)

     if(LIS_masterproc) then
        if(wformat.eq."binary") then 
           call LIS_releaseUnitNumber(ftn)
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in mos_writerst')
#endif           
        endif
        write(LIS_logunit,*)'mosaic archive restart written: ',trim(filen)
     endif
  endif
end subroutine mos_writerst

!BOP
! 
! !ROUTINE: mos_dump_restart
! \label{mos_dump_restart}
!
! !REVISION HISTORY:
!  25 Sep 2007: Sujay Kumar, Initial Specification
! !INTERFACE:
subroutine mos_dump_restart(ftn,n,wformat)
! !USES: 
  use LIS_coreMod, only  : LIS_rc, LIS_masterproc
  use mos_lsmMod   
  use LIS_historyMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: ftn 
  character(len=*), intent(in) :: wformat
!
! !DESCRIPTION:
!  This routine gathers the necessary restart variables and performs
!  the actual write statements to create the restart files.
!
!  The arguments are: 
!  \begin{description}
!   \item[n] 
!    index of the nest
!   \item[ftn]
!    unit number for the restart file
!   \item[wformat]
!    restart file format (binary/netcdf)
!  \end{description}
!
!  The following is the list of variables written in the Mosaic
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
! \item[LIS\_writeGlobalHeader\_restart](\ref{LIS_writeGlobalHeader_restart}) \newline
!  writes the global header information 
! \item[LIS\_writeHeader\_restart](\ref{LIS_writeHeader_restart}) \newline
!  writes the header information for a variable
! \item[LIS\_closeHeader\_restart](\ref{LIS_closeHeader_restart}) \newline
!  close the header
! \item[LIS\_writevar\_restart](\ref{LIS_writevar_restart}) \newline
!  writes a variable to the restart file
! \end{description}
!
!EOP
  integer             :: l,t
  integer             :: ctId, qaId, icsId, snowId,SoTId,SoWETid
  integer             :: dimID(11)
  real, allocatable   :: tmptilen(:)

  allocate(tmptilen(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  call LIS_writeGlobalHeader_restart(ftn,n,LIS_rc%lsm_index,&
       "Mosaic", dim1=3, dim2 = 1, dimID=dimID,&
       output_format="netcdf")
  
  call LIS_writeHeader_restart(ftn,n,dimID,ctId,"CT",&
       "CT","-",vlevels=1,valid_min=0.0,valid_max=0.0)
  call LIS_writeHeader_restart(ftn,n,dimID,qaId,"QA",&
       "QA","-",vlevels=1,valid_min=0.0,valid_max=0.0)
  call LIS_writeHeader_restart(ftn,n,dimID,icsId,"ICS",&
       "ICS","-",vlevels=1,valid_min=0.0,valid_max=0.0)
  call LIS_writeHeader_restart(ftn,n,dimID,snowId,"SNOW",&
       "SNOW","-",vlevels=1,valid_min=0.0,valid_max=0.0)
  call LIS_writeHeader_restart(ftn,n,dimID,soTId,"SoT",&
       "SoT","K",vlevels=1,valid_min=0.0,valid_max=0.0)
  call LIS_writeHeader_restart(ftn,n,dimID,SoWETId,"SoWET",&
       "SoWET","-",vlevels=3,valid_min=0.0,valid_max=0.0,&
       var_flag="dim1")

  call LIS_closeHeader_restart(ftn,n,LIS_rc%lsm_index,&
       dimID,mos_struc(n)%rstInterval)

  call LIS_writevar_restart(ftn,n,mos_struc(n)%mos%ct,&
       varid=ctId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,mos_struc(n)%mos%qa,&
       varid=qaId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,mos_struc(n)%mos%ics,&
       varid=icsId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,mos_struc(n)%mos%snow,&
       varid=snowId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,mos_struc(n)%mos%SoT,&
       varid=soTId,dim=1,wformat=wformat)

  do l=1,3
     tmptilen = 0.0
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = mos_struc(n)%mos(t)%SoWET(l)
     enddo
     call LIS_writevar_restart(ftn,n,tmptilen,&
          varid=sowetId,dim=l,wformat=wformat)
  enddo  
  deallocate(tmptilen)

end subroutine mos_dump_restart
