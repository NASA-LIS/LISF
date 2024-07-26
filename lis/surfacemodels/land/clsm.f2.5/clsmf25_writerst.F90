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
! !ROUTINE: clsmf25_writerst
! \label{clsmf25_writerst}
!
! !REVISION HISTORY:
! 18 Dec 2005: Sujay Kumar, Initial Specification
!
!
! !INTERFACE:
subroutine clsmf25_writerst(n)
! !USES:

  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_fileIOMod, only : LIS_create_output_directory, &
                            LIS_create_restart_filename
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use clsmf25_lsmMod      
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

! !ARGUMENTS: 
  implicit none       
  integer, intent(in) :: n 
!
! !DESCRIPTION:
!  This program writes restart files for catchment.  This
!  includes all relevant water/energy storage and tile information
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[cat\_dump\_restart](\ref{clsmf25_dump_restart}) \newline
!   writes the catchment variables into the restart file
! \end{description}
!EOP
  character(len=LIS_CONST_PATH_LEN) :: filen
  logical       :: alarmCheck
  integer       :: ftn 
  integer       :: status
  character*20  :: wformat
  logical       :: check_flag
  character*3 :: fnest

  write(fnest,'(i3.3)') n

  alarmCheck = LIS_isAlarmRinging(LIS_rc, &
       "CLSM F2.5 restart alarm "//trim(fnest))
  
  wformat = "netcdf"
  if(LIS_rc%runmode.eq."ensemble smoother") then 
     check_flag = alarmCheck
  else
     check_flag = alarmCheck.or.(LIS_rc%endtime==1)
  endif

  if(check_flag) then 
     if ( LIS_masterproc ) then
        call LIS_create_output_directory('SURFACEMODEL')
        call LIS_create_restart_filename(n,filen,'SURFACEMODEL','CLSMF25', &
             wformat=wformat)
        if(wformat.eq."binary") then 
           ftn = LIS_getNextUnitNumber()           
           open(ftn,file=filen,status='unknown',form='unformatted')
        elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF3)
           status = nf90_create(path=filen,cmode=nf90_clobber,&
                ncid = ftn)
           call LIS_verify(status,'Error in nf90_open in clsmf25_writerst')
#endif
#if (defined USE_NETCDF4)
           status = nf90_create(path=trim(filen),cmode=nf90_hdf5,&
                ncid = ftn)
           call LIS_verify(status,'Error in nf90_open in clsmf25_writerst')
#endif
        endif
     endif
     
     call clsmf25_dump_restart(n,ftn,wformat)
     
     if ( LIS_masterproc ) then
        if(wformat.eq."binary") then 
           call LIS_releaseUnitNumber(ftn)
        elseif(wformat.eq."netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in clsmf25_writerestart')
#endif
        endif
        write(LIS_logunit,*)                                           &
           'Catchment Fortuna-2.5 archive restart written: ',trim(filen)
     endif
  end if
  return

end subroutine clsmf25_writerst

!BOP
! 
! !ROUTINE: clsmf25_dump_restart
!  \label{clsmf25_dump_restart}
!
! !REVISION HISTORY:
!  24 Aug 2004: James Geiger, Initial Specification
! !INTERFACE:
subroutine clsmf25_dump_restart(n, ftn, wformat)

! !USES:
   use LIS_coreMod, only : LIS_rc, LIS_masterproc
   use LIS_historyMod
   use clsmf25_constants
   use clsmf25_lsmMod

   implicit none
! 
! !ARGUMENTS: 
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
!EOP
   integer :: t,l
   real, allocatable :: temp(:)
   integer :: dimID(11)
   integer :: tc1Id, tc2Id, tc4Id, qa1Id, qa2Id, qa4Id
   integer :: capacId, ghtId, wesnId, htsnId, sndzId
   integer :: catdefId, rzexcId, srfexcId
   
   allocate(temp(LIS_rc%npatch(n,LIS_rc%lsm_index)))
   call LIS_writeGlobalHeader_restart(ftn,n,LIS_rc%lsm_index,&
        "Catchment",dimID, dim1=3,dim2=N_gt,&
        dim3=N_snow,output_format="netcdf")

   call LIS_writeHeader_restart(ftn,n,dimID,tc1Id,"TC1",&
        "TC1","K",vlevels=1,valid_min=0.0,valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,tc2Id,"TC2",&
        "TC2","K",vlevels=1,valid_min=0.0,valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,tc4Id,"TC4",&
        "TC4","K",vlevels=1,valid_min=0.0,valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,qa1Id,"QA1",&
        "QA1","-",vlevels=1,valid_min=0.0,valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,qa2Id,"QA2",&
        "QA2","-",vlevels=1,valid_min=0.0,valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,qa4Id,"QA4",&
        "QA4","-",vlevels=1,valid_min=0.0,valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,capacId,"CAPAC",&
        "CAPAC","-",vlevels=1,valid_min=0.0,valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,ghtId,"GHT",&
        "GHT","-",vlevels=N_gt,valid_min=0.0,valid_max=0.0,&
        var_flag="dim2")
   call LIS_writeHeader_restart(ftn,n,dimID,wesnId,"WESN",&
        "WESN","-",vlevels=N_snow,valid_min=0.0,valid_max=0.0,&
        var_flag = "dim3")
   call LIS_writeHeader_restart(ftn,n,dimID,htsnId,"HTSN",&
        "HTSN","-",vlevels=N_snow,valid_min=0.0,valid_max=0.0,&
        var_flag = "dim3")
   call LIS_writeHeader_restart(ftn,n,dimID,sndzId,"SNDZ",&
        "SNDZ","-",vlevels=N_snow,valid_min=0.0,valid_max=0.0,&
        var_flag = "dim3")
   call LIS_writeHeader_restart(ftn,n,dimID,catdefId,"CATDEF",&
        "CATDEF","-",vlevels=1,valid_min=0.0,valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,rzexcId,"RZEXC",&
        "RZEXC","-",vlevels=1,valid_min=0.0,valid_max=0.0)
   call LIS_writeHeader_restart(ftn,n,dimID,srfexcId,"SRFEXC",&
        "SRFEXC","-",vlevels=1,valid_min=0.0,valid_max=0.0)

   call LIS_closeHeader_restart(ftn,n,LIS_rc%lsm_index,dimID,&
        clsmf25_struc(n)%rstInterval)
   
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
        clsmf25_struc(n)%cat_progn%tc1,varid=tc1Id,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
        clsmf25_struc(n)%cat_progn%tc2,varid=tc2Id,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
        clsmf25_struc(n)%cat_progn%tc4,varid=tc4Id,dim=1,wformat=wformat)

   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
        clsmf25_struc(n)%cat_progn%qa1,varid=qa1Id,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
        clsmf25_struc(n)%cat_progn%qa2,varid=qa2Id,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
        clsmf25_struc(n)%cat_progn%qa4,varid=qa4Id,dim=1,wformat=wformat)

   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
        clsmf25_struc(n)%cat_progn%capac,varid=capacId,dim=1,wformat=wformat)

   do l=1,N_gt
      do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         temp(t) = clsmf25_struc(n)%cat_progn(t)%ght(l)
      enddo
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,temp,&
           varid=ghtId,dim=l,wformat=wformat)
   enddo

   do l=1,N_snow
      do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         temp(t) = clsmf25_struc(n)%cat_progn(t)%wesn(l)
      enddo
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,temp,&
           varid=wesnId,dim=l,wformat=wformat)
   enddo

   do l=1,N_snow
      do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         temp(t) = clsmf25_struc(n)%cat_progn(t)%htsn(l)
      enddo
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,temp,&
           varid=htsnId,dim=l,wformat=wformat)
   enddo

   do l=1,N_snow
      do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         temp(t) = clsmf25_struc(n)%cat_progn(t)%sndz(l)
      enddo
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,temp,&
           varid=sndzId,dim=l,wformat=wformat)
   enddo
   
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
        clsmf25_struc(n)%cat_progn%catdef,varid=catdefId,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
        clsmf25_struc(n)%cat_progn%rzexc,varid=rzexcId,dim=1,wformat=wformat)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
        clsmf25_struc(n)%cat_progn%srfexc,varid=srfexcId,dim=1,wformat=wformat)

   deallocate(temp)
end subroutine clsmf25_dump_restart

