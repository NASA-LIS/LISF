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
! !ROUTINE: read_WUS_UCLAsnow
! \label{read_WUS_UCLAsnow}
!
! !REVISION HISTORY:
!  08 Jun 2022: Sujay Kumar; Initial version
!
! !INTERFACE: 
subroutine read_WUS_UCLAsnow(n, k, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use map_utils
  use LIS_pluginIndices
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use WUS_UCLAsnowMod, only : WUS_UCLAsnow_struc
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the WUS_UCLAsnow observations and prepares them for DA
!  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  integer                :: ftn,status
  character(len=LIS_CONST_PATH_LEN) :: sndobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname
  logical                :: alarmCheck
  integer                :: t,c,r,i,j,p,jj
  integer                :: sndId,ierr
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: sndfield, pertField
  logical*1              :: sndobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  logical*1, allocatable :: snd_data_b(:)
  real                   :: sndobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real, allocatable      :: var(:,:)
  real, allocatable      :: snd1d(:)
  real                   :: snd_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer                :: fnd
  logical                :: file_exists
    
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sndobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  sndobs = LIS_rc%udef
  data_upd = .false. 
!-------------------------------------------------------------------------
!   Read both ascending and descending passes at 0Z and then store
!   the overpass time as 1.30AM for the descending pass and 1.30PM 
!   for the ascending pass. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "WUS_UCLAsnow read alarm")
  
  if(alarmCheck.or.WUS_UCLAsnow_struc(n)%startMode) then
     WUS_UCLAsnow_struc(n)%startMode = .false.
     sndobs = LIS_rc%udef
     
     if(WUS_UCLAsnow_struc(n)%obsMap) then             
        call create_WUS_UCLAsnow_filename(sndobsdir, &
             LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname)
        
        inquire(file=fname,exist=file_exists)
        
        if(file_exists) then
           
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
           allocate(snd_data_b(WUS_UCLAsnow_struc(n)%nc*WUS_UCLAsnow_struc(n)%nr))
           allocate(var(WUS_UCLAsnow_struc(n)%nc,WUS_UCLAsnow_struc(n)%nr))
           allocate(snd1d(WUS_UCLAsnow_struc(n)%nc*WUS_UCLAsnow_struc(n)%nr))
           
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
           ierr = nf90_open(path=fname,mode=NF90_NOWRITE,ncid=ftn)
           call LIS_verify(ierr,'error opening WUS UCLA file')
           
           ierr = nf90_inq_varid(ftn,'SD_Post',sndId)
           call LIS_verify(ierr, 'nf90_inq_varid failed for SD_Post in read_WUS_UCLAsnow')
           
           ierr = nf90_get_var(ftn,sndId, var,&
                count=(/WUS_UCLAsnow_struc(n)%nc,WUS_UCLAsnow_struc(n)%nr,1/),&
                start=(/WUS_UCLAsnow_struc(n)%c_off,WUS_UCLAsnow_struc(n)%r_off,1/))
           call LIS_verify(ierr, 'nf90_get_var failed for SD_Post')
           
           ierr = nf90_close(ftn)
           call LIS_verify(ierr)
#endif
           
           snd1d = LIS_rc%udef
           snd_data_b = .false.
           
           do r=1, WUS_UCLAsnow_struc(n)%nr
              do c=1, WUS_UCLAsnow_struc(n)%nc
                 if(var(c,r).ge.0) then 
                    snd1d(c+(r-1)*WUS_UCLAsnow_struc(n)%nc) = &
                         var(c,r)
                    snd_data_b(c+(r-1)*WUS_UCLAsnow_struc(n)%nc)=.true.
                 endif
              enddo
           enddo
           deallocate(var)
!--------------------------------------------------------------------------
! Interpolate to the observation grid
!-------------------------------------------------------------------------- 

           call upscaleByAveraging(&
                WUS_UCLAsnow_struc(n)%nc*WUS_UCLAsnow_struc(n)%nr,&
                LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
                LIS_rc%udef,&
                WUS_UCLAsnow_struc(n)%n11,&
                snd_data_b,snd1d,&
                sndobs_b_ip,sndobs)
           
           
           deallocate(snd1d)
           deallocate(snd_data_b)
           
        endif
        
     endif

     call ESMF_StateGet(OBS_State,"Observation01",sndfield,&
             rc=status)
     call LIS_verify(status, 'Error: StateGet Observation01')
     
     call ESMF_FieldGet(sndfield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status, 'Error: FieldGet')
     
     fnd = 0 
     
     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              
              if(sndobs(c+(r-1)*LIS_rc%obs_lnc(k)).ge.0) then 
                 fnd = 1
                 exit
              endif
           endif
        enddo
     enddo
     
     obsl = LIS_rc%udef 
     if(fnd.ne.0) then 
        do r=1, LIS_rc%obs_lnr(k)
           do c=1, LIS_rc%obs_lnc(k)
              if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                 obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                      sndobs(c+(r-1)*LIS_rc%obs_lnc(k))
              endif
           enddo
        enddo
     endif
     
!-------------------------------------------------------------------------
!  Apply LSM based quality control and screening of observations
!-------------------------------------------------------------------------     
     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
          //trim(LIS_wusUCLAobsId)//char(0),n, k, OBS_state)
     
     snd_current = LIS_rc%udef
     call LIS_checkForValidObs(n,k,obsl,fnd,sndobs)

     if(fnd.eq.0) then 
        data_upd_flag_local = .false. 
     else
        data_upd_flag_local = .true. 
     endif
     
#if (defined SPMD)
     call MPI_ALLGATHER(data_upd_flag_local,1, &
          MPI_LOGICAL, data_upd_flag(:),&
          1, MPI_LOGICAL, LIS_mpi_comm, status)
    data_upd = any(data_upd_flag)
#else
    data_upd = data_upd_flag_local
#endif
     
     if(data_upd) then
        do t=1,LIS_rc%obs_ngrid(k)
           gid(t) = t
           if(obsl(t).ne.-9999.0) then 
              assimflag(t) = 1
              if(obsl(t).gt.10.0) then
                  !print *,'WARNING: OBSL > 10.0m; OBSL=',obsl(t)
                  assimflag(t) = 0 !Do not assimilate values > 10m
              endif
           else
              assimflag(t) = 0
           endif
        enddo
        
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .true. , rc=status)
        call LIS_verify(status)
        
        if(LIS_rc%obs_ngrid(k).gt.0) then 
           call ESMF_AttributeSet(sndField,"Grid Number",&
                gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(sndField,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
        endif

     else
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .false., rc=status)
        call LIS_verify(status)     
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif
end subroutine read_WUS_UCLAsnow


!BOP
! !ROUTINE: create_WUS_UCLAsnow_filename
! \label{create_WUS_UCLAsnow_filename}
! 
! !INTERFACE: 
subroutine create_WUS_UCLAsnow_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the WUS_UCLAsnow data filename
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the WUS_UCLAsnow directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated WUS_UCLAsnow filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/'//fyr//&
       '/WUS_UCLA_SR_v01_agg_16_SWE_SCA_SD_POST_' &
       //fyr//fmo//fda//'.nc4'
  
end subroutine create_WUS_UCLAsnow_filename


