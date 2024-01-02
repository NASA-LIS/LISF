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
! !ROUTINE: read_SNODAS
! \label{read_SNODAS}
!
! !REVISION HISTORY:
!  09 Apr 2021: Sujay Kumar; Initial version
!
! !INTERFACE: 
subroutine read_SNODAS(n, k, OBS_State, OBS_Pert_State)
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
  use SNODAS_Mod, only : SNODAS_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the SNODAS observations and prepares them for DA
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
  integer*2, allocatable :: var(:,:)
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
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SNODAS read alarm")
  
  if(alarmCheck.or.SNODAS_struc(n)%startMode) then 

     SNODAS_struc(n)%startMode = .false.
     sndobs = LIS_rc%udef
        
     call create_SNODAS_filename(sndobsdir, &
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname)
     
     inquire(file=trim(fname),exist=file_exists)

     if(file_exists) then

        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
        allocate(snd_data_b(SNODAS_struc(n)%nc*SNODAS_struc(n)%nr))
        allocate(var(SNODAS_struc(n)%nc,SNODAS_struc(n)%nr))
        allocate(snd1d(SNODAS_struc(n)%nc*SNODAS_struc(n)%nr))
        
        ftn = LIS_getNextUnitNumber()
        open(ftn,file=trim(fname),form='unformatted',&
             access='direct',convert='big_endian',&
             recl=SNODAS_struc(n)%nc*SNODAS_struc(n)%nr*2)
        read(ftn,rec=1) var
        call LIS_releaseUnitNumber(ftn)


        snd1d = LIS_rc%udef
        snd_data_b = .false.
        
        do r=1, SNODAS_struc(n)%nr
           do c=1, SNODAS_struc(n)%nc
              if(var(c,SNODAS_struc(n)%nr-r+1).ge.0) then 
                 snd1d(c+(r-1)*SNODAS_struc(n)%nc) = &
                      float(var(c,SNODAS_struc(n)%nr-r+1))/1000.0
                 snd_data_b(c+(r-1)*SNODAS_struc(n)%nc)=.true.
              endif
           enddo
        enddo
        deallocate(var)
!--------------------------------------------------------------------------
! Interpolate to the observation grid
!-------------------------------------------------------------------------- 

        call upscaleByAveraging(&
             SNODAS_struc(n)%nc*SNODAS_struc(n)%nr,&
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
             LIS_rc%udef,&
             SNODAS_struc(n)%n11,&
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
       //trim(LIS_SNODASobsId)//char(0),n, k, OBS_state)

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
#endif
  data_upd = .false.
  do p=1,LIS_npes
     data_upd = data_upd.or.data_upd_flag(p)
  enddo

  if(data_upd) then 
     do t=1,LIS_rc%obs_ngrid(k)
        gid(t) = t
        if(obsl(t).ne.-9999.0) then 
           assimflag(t) = 1
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
end subroutine read_SNODAS


!BOP
! !ROUTINE: create_SNODAS_filename
! \label{create_SNODAS_filename}
! 
! !INTERFACE: 
subroutine create_SNODAS_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the SNODAS data filename
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the SNODAS directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated SNODAS filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//&
       '/us_ssmv11036tS__T0001TTNATS' &
       //trim(fyr)//trim(fmo)//trim(fda)//'05HP001.dat'
end subroutine create_SNODAS_filename


