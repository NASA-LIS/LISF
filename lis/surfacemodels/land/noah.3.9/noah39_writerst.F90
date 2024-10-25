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
! !ROUTINE: noah39_writerst
! \label{noah39_writerst}
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
!  14 Jun 2003: Sujay Kumar; Separated the write restart from the original 
!                code
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  10 Jun 2012: Sujay Kumar, added support for netcdf formats
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
!
! !INTERFACE:
subroutine noah39_writerst(n)
! !USES:
  use LIS_coreMod,   only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod,    only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_verify
  use LIS_fileIOMod, only : LIS_create_output_directory, &
                              LIS_create_restart_filename
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use noah39_lsmMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none       
! !ARGUMENTS: 
  integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for Noah-3.9.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[noah39\_dump\_restart](\ref{noah39_dump_restart}) \newline
!   writes the Noah-3.9 variables into the restart file
! \end{description}
!EOP
!  real :: curr_time
  character(len=LIS_CONST_PATH_LEN) :: filen
  logical       :: alarmCheck
  integer       :: ftn
  integer       :: status
  character*20  :: wformat
  character*3        :: fnest

  write(fnest,'(i3.3)') n
  alarmCheck = LIS_isAlarmRinging(LIS_rc,"Noah39 restart alarm "//trim(fnest))
       
  wformat = "netcdf"

  if(alarmCheck .or.(LIS_rc%endtime ==1)) then 
     if ( LIS_masterproc ) then
        call LIS_create_output_directory('SURFACEMODEL')
        call LIS_create_restart_filename(n,filen,'SURFACEMODEL','NOAH39',&
             wformat=wformat)
        if(wformat.eq."binary") then 
           ftn = LIS_getNextUnitNumber()       
           open(ftn,file=filen,status='unknown',form='unformatted')
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF4)
           status = nf90_create(path=filen,cmode=nf90_hdf5,&
               ncid = ftn)
           call LIS_verify(status,'Error in nf90_open in noah39_writerst')
#endif
#if (defined USE_NETCDF3)
           status = nf90_create(path=filen,cmode=nf90_clobber,&
                ncid = ftn)
           call LIS_verify(status,'Error in nf90_open in noah39_writerst')
#endif
        endif
     endif
     
     call noah39_dump_restart(n,ftn,wformat)
     
     if ( LIS_masterproc ) then
        if(wformat.eq."binary") then 
           call LIS_releaseUnitNumber(ftn)
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in noah39_writerst')
#endif
        endif
        write(LIS_logunit,*) '[INFO] Noah-3.9 archive restart written: ',trim(filen)
     endif
  end if
end subroutine noah39_writerst

!BOP
! 
! !ROUTINE: noah39_dump_restart
! \label{noah39_dump_restart}
!
! !REVISION HISTORY:
!  24 Aug 2004: James Geiger, Initial Specification
! !INTERFACE:
subroutine noah39_dump_restart(n, ftn, wformat)

! !USES:
   use LIS_coreMod,    only : LIS_rc, LIS_masterproc
   use LIS_historyMod
   use LIS_logMod,     only : LIS_logunit
   use noah39_lsmMod
   use LIS_tbotAdjustMod, only : LIS_writeTmnUpdateRestart

   implicit none

   integer, intent(in) :: ftn
   integer, intent(in) :: n
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
!  The following is the list of variables written in the Noah-3.9
!  restart file: 
!  \begin{verbatim}
!   nc,nr,ntiles    - grid and tile space dimensions 
!   t1           - Noah-3.9 skin temperature
!   cmc          - Noah-3.9 canopy moisture storage
!   snowh        - Noah-3.9 snow depth
!   sneqv        - Noah-3.9 snow water equivalent
!   stc          - Noah-3.9 soil temperature (for each layer)
!   smc          - Noah-3.9 soil moisture (for each layer) 
!   sh2o         - Noah-3.9 liquid only soil moisture (for each layer)
!   ch           - Noah-3.9 heat/moisture exchange coefficient
!   cm           - Noah-3.9 momentum exchange coefficient
!   z0           - Noah-3.9 roughness length
!   emiss        - Noah-3.9 emissivity
!   albedo       - Noah-3.9 albedo
!   snotime1     - Noah-3.9 snow timing variable
!   q1           - Noah-3.9 mixing ratio
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
   integer :: l,t
   integer :: dimID(11)
   integer :: t1Id, cmcid, snowhid, sneqvid, smcid, sh2oId, stcId
   integer :: cmId, chId,snotime1Id,emissId,z0Id,albId,q1Id
   integer :: xsda_qfx_Id 
   integer :: hfx_phy_Id
   integer :: qfx_phy_Id
   integer :: xqnorm_Id
   real, allocatable :: tmptilen(:)

   allocate(tmptilen(LIS_rc%npatch(n,LIS_rc%lsm_index)))

   call LIS_writeGlobalHeader_restart(ftn,n,LIS_rc%lsm_index,&
        "Noah version 3.9", dimID, dim1=noah39_struc(n)%nslay, &
        dim2= noah39_struc(n)%nslay,&
        output_format="netcdf")

   if ( wformat == "netcdf" ) then
      call LIS_writeHeader_restart(ftn,n,dimID,t1id, "T1", &
           "skin temperature",&
           "K",vlevels=1,valid_min = 220.0, valid_max=330.0)

      call LIS_writeHeader_restart(ftn,n,dimID,cmcid, "CMC", &
           "canopy moisture content",&
           "-",vlevels=1,valid_min = 0.0, valid_max=0.0)
      call LIS_writeHeader_restart(ftn,n,dimID,snowhid, "SNOWH",&
           "snow depth",&
           "m",vlevels=1,valid_min = 0.0, valid_max=0.0)
      call LIS_writeHeader_restart(ftn,n,dimID,sneqvid, "SNEQV",&
           "snow water equivalent",&
           "m",vlevels=1,valid_min = 0.0, valid_max=0.0)
      call LIS_writeHeader_restart(ftn,n,dimID,stcid, "STC",&
           "soil temperature",&
           "K",vlevels=noah39_struc(n)%nslay,valid_min = 220.0, valid_max=330.0,&
           var_flag = "dim2")
      call LIS_writeHeader_restart(ftn,n,dimID,smcid, "SMC",&
           "total soil moisture",&
           "m^3 m-3",vlevels=noah39_struc(n)%nslay,valid_min =0.0, valid_max=0.5,&
           var_flag = "dim1")
      call LIS_writeHeader_restart(ftn,n,dimID,sh2oid, "SH2O",&
           "liquid soil moisture",&
           "m^3 m-3",vlevels=noah39_struc(n)%nslay,valid_min =0.0, valid_max=0.5,&
           var_flag = "dim1")
      call LIS_writeHeader_restart(ftn,n,dimID,cmid, "CM",&
           "exchange coefficient for momentum",&
           "-",vlevels=1,valid_min =0.0, valid_max=0.0)
      call LIS_writeHeader_restart(ftn,n,dimID,chid, "CH",&
           "exchange coefficient for heat",&
           "-",vlevels=1,valid_min =0.0, valid_max=0.0)
      call LIS_writeHeader_restart(ftn,n,dimID,z0id, "Z0",&
           "roughness length",&
           "m",vlevels=1,valid_min =0.0, valid_max=0.0)
      call LIS_writeHeader_restart(ftn,n,dimID,emissid, "EMISS",&
           "surface emissivity",&
           "-",vlevels=1,valid_min =0.0, valid_max=0.0)
      call LIS_writeHeader_restart(ftn,n,dimID,albid, "ALBEDO",&
           "surface albedo",&
           "-",vlevels=1,valid_min =0.0, valid_max=0.0)
      call LIS_writeHeader_restart(ftn,n,dimID,snotime1id, "SNOTIME1",&
           "Snotime1",&
           "-",vlevels=1,valid_min =0.0, valid_max=0.0)
      call LIS_writeHeader_restart(ftn,n,dimID,q1id, "Q1",&
           "Q1",&
           "-",vlevels=1,valid_min =0.0, valid_max=0.0)
      if(noah39_struc(n)%opt_fasdas .eq. 1) then 
        call LIS_writeHeader_restart(ftn,n,dimID,xsda_qfx_Id, "XSDA_QFX",&
             "XSDA_QFX",&
             "-",vlevels=1,valid_min =0.0, valid_max=0.0)
      endif 

      if(noah39_struc(n)%opt_fasdas .eq. 1) then
        call LIS_writeHeader_restart(ftn,n,dimID,hfx_phy_Id, "HFX_PHY",&
             "HFX_PHY",&
             "-",vlevels=1,valid_min =0.0, valid_max=0.0)
      endif

      if(noah39_struc(n)%opt_fasdas .eq. 1) then
        call LIS_writeHeader_restart(ftn,n,dimID,qfx_phy_Id, "QFX_PHY",&
             "QFX_PHY",&
             "-",vlevels=1,valid_min =0.0, valid_max=0.0)
      endif

      if(noah39_struc(n)%opt_fasdas .eq. 1) then
        call LIS_writeHeader_restart(ftn,n,dimID,xqnorm_Id, "XQNORM",&
             "XQNORM",&
             "-",vlevels=1,valid_min =0.0, valid_max=0.0)
      endif

      call LIS_closeHeader_restart(ftn,n,LIS_rc%lsm_index,dimID,noah39_struc(n)%rstInterval)
   endif

   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%t1,&
        varid=t1Id,dim=1,wformat=wformat)

   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%cmc,&
        varid=cmcId,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%snowh,&
        varid=snowhId,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%sneqv,&
        varid=sneqvId,dim=1,wformat=wformat)
   do l=1,noah39_struc(n)%nslay
      tmptilen = 0.0
      do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         tmptilen(t) = noah39_struc(n)%noah(t)%stc(l)         
      enddo
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen, &
           varid=stcId,dim=l,wformat=wformat)
   enddo
   do l=1,noah39_struc(n)%nslay
      tmptilen = 0.0
      do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         tmptilen(t) = noah39_struc(n)%noah(t)%smc(l)         
      enddo
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
           varid=smcId,dim=l,wformat=wformat)
   enddo

   do l=1,noah39_struc(n)%nslay
      tmptilen = 0.0
      do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         tmptilen(t) = noah39_struc(n)%noah(t)%sh2o(l)         
      enddo
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
           varid=sh2oId,dim=l,wformat=wformat)
   enddo

   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%ch,&
        varid=chId,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%cm,&
        varid=cmId,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%z0,&
        varid=z0Id,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%emiss,&
        varid=emissId,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%albedo,&
        varid=albId,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%snotime1,&
        varid=snotime1Id,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%q1,&
        varid=q1Id,dim=1,wformat=wformat)

   if(noah39_struc(n)%opt_fasdas .eq. 1) then 
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%xsda_qfx,&
          varid=xsda_qfx_Id,dim=1,wformat=wformat)
   endif

   if(noah39_struc(n)%opt_fasdas .eq. 1) then
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%hfx_phy,&
          varid=hfx_phy_Id,dim=1,wformat=wformat)
   endif

   if(noah39_struc(n)%opt_fasdas .eq. 1) then
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%qfx_phy,&
          varid=qfx_phy_Id,dim=1,wformat=wformat)
   endif

   if(noah39_struc(n)%opt_fasdas .eq. 1) then
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,noah39_struc(n)%noah%xqnorm,&
          varid=xqnorm_Id,dim=1,wformat=wformat)
   endif

   !Write extra data for dynamic deep soil temperature option
   call LIS_writeTmnUpdateRestart(n,ftn,dimID,wformat)

   deallocate(tmptilen)
end subroutine noah39_dump_restart
