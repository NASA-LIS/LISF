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
! !ROUTINE: hyssib_writerst
! \label{hyssib_writerst}
!
! !REVISION HISTORY:
! 21 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
! 25 Aug 2007: Chuck Alonge, Updates for LIS 5.0 Compliance
!
! RESTART FILE FORMAT(fortran sequential binary):
!  YR,MO,DA,HR,MN,SS,VCLASS,NCH !Restart time,Veg class,no.tiles, no.soil lay
!  TILE(NCH)%COL        !Grid Col of Tile
!  TILE(NCH)%ROW        !Grid Row of Tile
!  TILE(NCH)%FGRD       !Fraction of Grid covered by Tile
!  TILE(NCH)%VEGT       !Vegetation Type of Tile
!  HYSSIB(NCH)%STATES   !Model States in Tile Space
!
! !INTERFACE:
subroutine hyssib_writerst(n)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod
  use lis_fileIOMod, only : LIS_create_output_directory, &
                              LIS_create_restart_filename
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use hyssib_lsmMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n

!
! !DESCRIPTION:
!  This program writes restart files for Hyssib.  This
!  includes all relevant water/energy storage and tile information
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[hyssib\_dump\_restart](\ref{hyssib_dump_restart}) \newline
!   writes the hyssib variables into the restart file
! \end{description}
!EOP

  character(len=LIS_CONST_PATH_LEN) :: filen
  logical       :: alarmCheck
  integer       :: status
  integer       :: ftn
  character*20  :: wformat
  character*3   :: fnest

  write(fnest,'(i3.3)') n

  wformat ="netcdf"

  alarmCheck = LIS_isAlarmRinging(LIS_rc, &
       "Hyssib restart alarm "//trim(fnest))
  if(alarmCheck .or. (LIS_rc%endtime ==1)) then 

     if ( LIS_masterproc ) then
        call LIS_create_output_directory('SURFACEMODEL')
        call LIS_create_restart_filename(n,filen,'SURFACEMODEL','HYSSIB',&
             wformat=wformat)
        if(wformat.eq."binary") then 
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=filen,status='unknown',form='unformatted')
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3)
           status = nf90_create(path=filen,cmode=nf90_clobber,&
                ncid = ftn)
           call LIS_verify(status,'Error in nf90_open in hyssib_writerst')
#endif
#if (defined USE_NETCDF4)
           status = nf90_create(path=trim(filen),cmode=nf90_hdf5,&
               ncid = ftn)
           call LIS_verify(status,'Error in nf90_open in hyssib_writerst')
#endif           
        endif
     endif
     call hyssib_dump_restart(n,ftn,wformat)
     
     if ( LIS_masterproc ) then
        if(wformat.eq."binary") then 
           call LIS_releaseUnitNumber(ftn)
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in hyssib_writerst')
#endif
        endif
        write(LIS_logunit,*) 'Noah3.3 archive restart written: ',trim(filen)
     endif
  end if
end subroutine hyssib_writerst

!BOP
!
! !ROUTINE: hyssib_dump_restart
! \label{hyssib_dump_restart}
!
! !REVISION HISTORY:
! 24 Aug 2004: James Geiger, Initial Specification
! 27 Sep 2005: Chuck Alonge, Updates for LIS 5.0
!
! !INTERFACE:
subroutine hyssib_dump_restart(n, ftn,wformat)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_historyMod
  use hyssib_lsmMod
      
  implicit none      
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: ftn
  character(len=*)    :: wformat
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
!  The following is the list of variables written in the Hyssib 
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
  integer :: l
  integer :: dimID(11)
  integer :: tcId, tgId, tsnId, tdId, wwwId, capacId, snowId, sgfgId, sdensId

  call LIS_writeGlobalHeader_restart(ftn,n,LIS_rc%lsm_index,&
       "HySSIB", dim1=3, dim2= 3, dim3=2,dimID=dimID)
  
  call LIS_writeHeader_restart(ftn,n,dimID,tcid, "TC", &
       "TC","K",vlevels=1,valid_min = 220.0, valid_max=330.0)
  call LIS_writeHeader_restart(ftn,n,dimID,tgid, "TG", &
       "TG","-",vlevels=1,valid_min = 0.0, valid_max=0.0)
  call LIS_writeHeader_restart(ftn,n,dimID,tsnid, "TSN", &
       "TSN","-",vlevels=1,valid_min = 0.0, valid_max=0.0)
  call LIS_writeHeader_restart(ftn,n,dimID,tsnid, "TD", &
       "TD","-",vlevels=1,valid_min = 0.0, valid_max=0.0)
  call LIS_writeHeader_restart(ftn,n,dimID,tsnid, "WWW", &
       "WWW","-",vlevels=3,valid_min = 0.0, valid_max=0.0,&
       var_flag = "dim1")
  call LIS_writeHeader_restart(ftn,n,dimID,tsnid, "CAPAC", &
       "CAPAC","-",vlevels=2,valid_min = 0.0, valid_max=0.0,&
       var_flag = "dim3")
  call LIS_writeHeader_restart(ftn,n,dimID,tsnid, "SNOW", &
       "SNOW","-",vlevels=2,valid_min = 0.0, valid_max=0.0,&
       var_flag = "dim3")
  call LIS_writeHeader_restart(ftn,n,dimID,tsnid, "SGFG", &
       "SGFG","-",vlevels=1,valid_min = 0.0, valid_max=0.0)
  call LIS_writeHeader_restart(ftn,n,dimID,tsnid, "SDENS", &
       "SDENS","-",vlevels=1,valid_min = 0.0, valid_max=0.0)

  call LIS_closeHeader_restart(ftn,n,LIS_rc%lsm_index,dimID,&
       hyssib_struc(n)%rstInterval)

  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,hyssib_struc(n)%hyssib%tc,&
       varid=tcId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,hyssib_struc(n)%hyssib%tg,&
       varid=tgId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,hyssib_struc(n)%hyssib%tsn,&
       varid=tsnId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,hyssib_struc(n)%hyssib%td,&
       varid=tdId,dim=1,wformat=wformat)
  do l=1,3
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
          hyssib_struc(n)%hyssib%www(l),&
          varid=wwwId,dim=l,wformat=wformat)
  enddo
  do l=1,2
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
          hyssib_struc(n)%hyssib%capac(l),&
          varid=capacId,dim=l,wformat=wformat)
  enddo
  do l=1,2
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
          hyssib_struc(n)%hyssib%snow(l),&
          varid=snowId,dim=l,wformat=wformat)
  enddo
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
       hyssib_struc(n)%hyssib%sgfg,&
       varid=sgfgId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
       hyssib_struc(n)%hyssib%sdens,&
       varid=sdensId,dim=1,wformat=wformat)
  
end subroutine hyssib_dump_restart
