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
! !ROUTINE: cable_writerst
! \label{cable_writerst}
!
! !REVISION HISTORY:
!  27 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!  13 Sep 2011: Claire Carouge (ccc), CABLE LSM improvements
!  10 Jun 2012: Sujay Kumar, added support for netcdf formats
!  23 May 2013: David Mocko, latest CABLE v1.4b version for LIS6.2
!
! !INTERFACE:
subroutine cable_writerst(n)
! !USES:
  use LIS_coreMod,   only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod,    only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_verify
  use LIS_fileIOMod, only : LIS_create_output_directory, &
                              LIS_create_restart_filename
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use cable_lsmMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
!
! !DESCRIPTION:
!  This subroutine writes the state variables to restart files for CABLE.
!  All relevant water/energy storages and tile info is written.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!    Creates a timestamped directory for the restart files
!   \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!    Creates a timestamped restart filename
!   \item[cable\_dump\_restart](\ref{cable_dump_restart}) \newline
!    Writes the CABLE variables into the restart file
!  \end{description}
!EOP
  character(len=LIS_CONST_PATH_LEN) :: filen
  logical       :: alarmCheck
  integer       :: status
  integer       :: ftn 
  character*20  :: wformat
  character*3   :: fnest
  
  write(fnest,'(i3.3)') n
  
  alarmCheck = LIS_isAlarmRinging(LIS_rc,"CABLE restart alarm "//trim(fnest))
       
  wformat = "netcdf"

  if(alarmCheck .or. (LIS_rc%endtime ==1)) then 

     if ( LIS_masterproc ) then
        call LIS_create_output_directory('SURFACEMODEL')
        call LIS_create_restart_filename(n,filen,'SURFACEMODEL','CABLE',&
             wformat=wformat)
        if(wformat.eq."binary") then 
           ftn = LIS_getNextUnitNumber()       
           open(ftn,file=filen,status='unknown',form='unformatted')
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3)
           status = nf90_create(path=filen,cmode=nf90_clobber,&
                ncid = ftn)
           call LIS_verify(status,'Error in nf90_open in cable_writerst')
#endif
#if (defined USE_NETCDF4)
           status = nf90_create(path=trim(filen),cmode=nf90_hdf5,&
               ncid = ftn)
           call LIS_verify(status,'Error in nf90_open in cable_writerst')
#endif
        endif
     endif

     call cable_dump_restart(n,ftn,wformat)

     
     if ( LIS_masterproc ) then
        if(wformat.eq."binary") then 
           call LIS_releaseUnitNumber(ftn)
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in cable_writerst')
#endif
        endif
        write(LIS_logunit,*) 'CABLE archive restart written: ',trim(filen)
     endif
  end if

end subroutine cable_writerst

!BOP
!
! !ROUTINE: cable_dump_restart
! \label{cable_dump_restart}
!
! !REVISION HISTORY:
!  27 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!
! !INTERFACE:
subroutine cable_dump_restart(n,ftn,wformat)
        
! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_historyMod
  use LIS_logMod,     only : LIS_logunit
  use cable_dimensions,   only : ms,msn,ncp,ncs
  use cable_lsmMod
  
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: ftn
  integer, intent(in) :: n
  character(len=*), intent(in) :: wformat
!
! !DESCRIPTION:
!  This subroutine gathers the necessary restart variables and performs
!  the actual write statements to create the CABLE restart files.
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
!  The following is the list of variables written in the
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
!   ktau         - CABLE time step counter
!   \end{verbatim}
!
!  The routines invoked are:
!  \begin{description}
! \item[LIS\_writeGlobalHeader\_restart](\ref{LIS_writeGlobalHeader_restart}) \newline
!  writes the global header information 
! \item[LIS\_writeHeader\_restart](\ref{LIS_writeHeader_restart}) \newline
!  writes the header information for a variable
! \item[LIS\_closeHeader\_restart](\ref{LIS_closeHeader_restart}) \newline
!  close the header
! \item[LIS\_writevar\_restart](\ref{LIS_writevar_restart}) \newline
!  writes a variable to the restart file
!  \end{description}
!EOP
  integer :: l,t
  integer :: dimID(11)
  integer :: canstoId,rtsoilId, ssdnnId, snowdId,osnowdId
  integer :: snageId,isflagId, wbtotId,wbiceId,tggsnId
  integer :: ssdnId, smassId,albsoilsnId,wgId,tggId
  integer :: cplantId,csoilId,gammzzId,wbId
  integer :: ktauId
  real, allocatable :: tmptilen(:)

  allocate(tmptilen(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  call LIS_writeGlobalHeader_restart(ftn,n,LIS_rc%lsm_index,&
        "CABLE", dim1=ms, dim2= ms, dim3=msn,dim4 =ncp, &
        dim5=ncs, dimID=dimID)
   
   call LIS_writeHeader_restart(ftn,n,dimID,canstoid, "CANSTO", &
        "CANSTO","-",vlevels=1,valid_min =0.0, valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,rtsoilid, "RTSOIL", &
        "RTSOIL","-",vlevels=1,valid_min =0.0, valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,ssdnnid, "SSDNN", &
        "SSDNN","-",vlevels=1,valid_min =0.0, valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,snowdid, "SNOWD", &
        "SNOWD","-",vlevels=1,valid_min =0.0, valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,osnowdid, "OSNOWD", &
        "OSNOWD","-",vlevels=1,valid_min =0.0, valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,snageid, "SNAGE", &
        "SNAGE","-",vlevels=1,valid_min =0.0, valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,isflagid, "ISFLAG", &
        "ISFLAG","-",vlevels=1,valid_min =0.0, valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,wbtotid, "WBTOT", &
        "WBTOT","-",vlevels=1,valid_min =0.0, valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,wbiceid, "WBICE", &
        "WBICE","-",vlevels=ms,valid_min =0.0, valid_max=0.0,&
        var_flag ="dim1")
   call LIS_writeHeader_restart(ftn,n,dimID,tggsnid, "TGGSN", &
        "TGGSN","-",vlevels=ms,valid_min =0.0, valid_max=0.0,&
        var_flag ="dim3")
   call LIS_writeHeader_restart(ftn,n,dimID,ssdnid, "SSDN", &
        "SSDN","-",vlevels=ms,valid_min =0.0, valid_max=0.0,&
        var_flag ="dim3")
   call LIS_writeHeader_restart(ftn,n,dimID,smassid, "SMASS", &
        "SMASS","-",vlevels=ms,valid_min =0.0, valid_max=0.0,&
        var_flag ="dim3")
   call LIS_writeHeader_restart(ftn,n,dimID,albsoilsnid, "ALBSOILSN", &
        "ALBSOILSN","-",vlevels=ms,valid_min =0.0, valid_max=0.0,&
        var_flag ="dim3")
   call LIS_writeHeader_restart(ftn,n,dimID,wbid, "WB", &
        "WB","-",vlevels=ms,valid_min =0.0, valid_max=0.0,&
        var_flag ="dim1")
   call LIS_writeHeader_restart(ftn,n,dimID,tggid, "TGG", &
        "TGG","-",vlevels=ms,valid_min =0.0, valid_max=0.0,&
        var_flag ="dim1")
   call LIS_writeHeader_restart(ftn,n,dimID,cplantid, "CPLANT", &
        "CPLANT","-",vlevels=ncp,valid_min =0.0, valid_max=0.0,&
        var_flag ="dim4")
   call LIS_writeHeader_restart(ftn,n,dimID,csoilid, "CSOIL", &
        "CSOIL","-",vlevels=ncs,valid_min =0.0, valid_max=0.0,&
        var_flag ="dim5")
   call LIS_writeHeader_restart(ftn,n,dimID,gammzzid, "GAMMZZ", &
        "GAMMZZ","-",vlevels=ms,valid_min =0.0, valid_max=0.0,&
        var_flag ="dim1")
   call LIS_writeHeader_restart(ftn,n,dimID,ktauid, "KTAU", &
        "KTAU","-",vlevels=1,valid_min =0.0, valid_max=0.0)

   call LIS_closeHeader_restart(ftn,n,LIS_rc%lsm_index,dimID,&
        cable_struc(n)%rstInterval)

  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,cable_struc(n)%cable%cansto,&
        varid=canstoId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,cable_struc(n)%cable%rtsoil,&
        varid=rtsoilId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,cable_struc(n)%cable%ssdnn,&
        varid=ssdnnId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,cable_struc(n)%cable%snowd,&
       varid=snowdId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,cable_struc(n)%cable%osnowd,&
        varid=osnowdId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,cable_struc(n)%cable%snage,&
        varid=snageId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,cable_struc(n)%cable%isflag,&
       varid=isflagId,dim=1,wformat=wformat)
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,cable_struc(n)%cable%wbtot,&
       varid=wbtotId,dim=1,wformat=wformat)
  tmptilen = 0.0
  do l = 1,ms
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = cable_struc(n)%cable(t)%wbice(l)
     enddo
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
          varid=wbiceId,dim=l,wformat=wformat)
  enddo
  tmptilen = 0.0
  do l = 1,msn
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = cable_struc(n)%cable(t)%tggsn(l)
     enddo
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
          varid=tggsnId,dim=l,wformat=wformat)
  enddo
  tmptilen = 0.0
  do l = 1,msn
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = cable_struc(n)%cable(t)%ssdn(l)
     enddo
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
          varid=ssdnId,dim=l,wformat=wformat)
  enddo
  tmptilen = 0.0
  do l = 1,msn
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = cable_struc(n)%cable(t)%smass(l)
     enddo
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
          varid=smassId,dim=l,wformat=wformat)
  enddo
  tmptilen = 0.0
  do l = 1,msn
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = cable_struc(n)%cable(t)%albsoilsn(l)
     enddo
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
          varid=albsoilsnId,dim=l,wformat=wformat)
  enddo
  tmptilen = 0.0
  do l = 1,ms
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = cable_struc(n)%cable(t)%wb(l)
     enddo
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
          varid=wbId,dim=l,wformat=wformat)
  enddo
  tmptilen = 0.0
  do l = 1,ms
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = cable_struc(n)%cable(t)%tgg(l)
     enddo
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
          varid=tggId,dim=l,wformat=wformat)
  enddo
  tmptilen = 0.0
  do l = 1,ncp
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = cable_struc(n)%cable(t)%cplant(l)
     enddo
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
          varid=cplantId,dim=l,wformat=wformat)
  enddo
  tmptilen = 0.0
  do l = 1,ncs
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = cable_struc(n)%cable(t)%csoil(l)
     enddo
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
          varid=csoilId,dim=l,wformat=wformat)
  enddo
  tmptilen = 0.0
  do l = 1,ms
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tmptilen(t) = cable_struc(n)%cable(t)%gammzz(l)
     enddo
     call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
          varid=gammzzId,dim=l,wformat=wformat)
  enddo
  tmptilen = 0.0
  tmptilen = cable_struc(n)%ktau
  call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen,&
        varid=ktauId,dim=1,wformat=wformat)

  deallocate(tmptilen)
end subroutine cable_dump_restart
