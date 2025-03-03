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
! !ROUTINE: read_SMOSNESDISsm
! \label{read_SMOSNESDISsm}
!
! !REVISION HISTORY:
!  17 Jun 2010: Sujay Kumar; Updated for use with LPRM AMSRE Version 5. 
!  20 Sep 2012: Sujay Kumar; Updated to the NETCDF version of the data. 
!
! !INTERFACE: 
subroutine read_SMOSNESDISsm(n, k, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use map_utils
  use LIS_pluginIndices
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use SMOSNESDISsm_Mod, only : SMOSNESDISsm_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the AMSRE soil moisture observations 
!  from NETCDF files and applies the spatial masking for dense
!  vegetation, rain and RFI. The data is then rescaled
!  to the land surface model's climatology using rescaling
!  algorithms. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  real, parameter        ::  minssdev = 0.05
  real, parameter        ::  maxssdev = 0.11
  real,  parameter       :: MAX_SM_VALUE=0.45, MIN_SM_VALUE=0.0001
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: smfield, pertField

  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  real                   :: obs_unsc(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: smobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: smtime(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: sm_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: dt
  real                   :: lon
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
  integer                :: fnd
  real, allocatable      :: ssdev(:)
  integer                :: lis_julss
  real                   :: smvalue
  real                   :: model_delta(LIS_rc%obs_ngrid(k))
  real                   :: obs_delta(LIS_rc%obs_ngrid(k))
  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 
  obs_unsc = LIS_rc%udef
!-------------------------------------------------------------------------
!   Read both ascending and descending passes at 0Z and then store
!   the overpass time as 1.30AM for the descending pass and 1.30PM 
!   for the ascending pass. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SMOSNESDIS read alarm")
  
  if(alarmCheck.or.SMOSNESDISsm_struc(n)%startMode) then 
     SMOSNESDISsm_struc(n)%startMode = .false.

     SMOSNESDISsm_struc(n)%smobs = LIS_rc%udef
     SMOSNESDISsm_struc(n)%smtime = -1

     smobs = LIS_rc%udef
     smtime = LIS_rc%udef

     
     call create_SMOSNESDISsm_filename(smobsdir, &
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, LIS_rc%hr,  fname)

     inquire(file=fname,exist=file_exists)

     if(file_exists) then 
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
        call read_SMOSNESDIS_data(n,k,fname,smobs,smtime)
     else
        write(LIS_logunit,*) '[WARN] Missing SMOSNESDIS ',trim(fname)
     endif

     SMOSNESDISsm_struc(n)%smobs  = LIS_rc%udef
     SMOSNESDISsm_struc(n)%smtime = -1

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           grid_index = LIS_obs_domain(n,k)%gindex(c,r)
           if(grid_index.ne.-1) then 
              if(smobs(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then 
                 SMOSNESDISsm_struc(n)%smobs(c,r) = &
                      smobs(c+(r-1)*LIS_rc%obs_lnc(k))                 
                 SMOSNESDISsm_struc(n)%smtime(c,r) = &
                      smtime(c+(r-1)*LIS_rc%obs_lnc(k))
              endif
           endif
        enddo
     enddo

  endif
  
  
  call ESMF_StateGet(OBS_State,"Observation01",smfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
  fnd = 0 
  sm_current = LIS_rc%udef
 
! dt is not defined as absolute value of the time difference to avoid
! double counting of the data in assimilation. 

  call LIS_get_timeoffset_sec(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
       LIS_rc%hr, LIS_rc%mn, LIS_rc%ss, lis_julss)

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           if(SMOSNESDISsm_struc(n)%smtime(c,r).ge.0) then 
              dt = (lis_julss-SMOSNESDISsm_struc(n)%smtime(c,r))
              if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
                 sm_current(c,r) = & 
                      SMOSNESDISsm_struc(n)%smobs(c,r)
                 fnd = 1
              endif           
           endif
        endif
     enddo
  enddo

  obsl = LIS_rc%udef 
  if(fnd.ne.0) then 
     do r=1, LIS_rc%obs_lnr(k)
        do c=1, LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=sm_current(c,r)
           endif
        enddo
     enddo
     
  endif

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     

  if(fnd.ne.0) then
     ! Store the unscaled obs (ie, before the rescaling)
     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = &
                   sm_current(c,r)
           end if
        end do
     end do
        
     call LIS_rescale_with_CDF_matching(    &
          n,k,                              & 
          SMOSNESDISsm_struc(n)%nbins,         & 
          SMOSNESDISsm_struc(n)%ntimes,         & 
          MAX_SM_VALUE,                        & 
          MIN_SM_VALUE,                        & 
          SMOSNESDISsm_struc(n)%model_xrange,  &
          SMOSNESDISsm_struc(n)%obs_xrange,    &
          SMOSNESDISsm_struc(n)%model_cdf,     &
          SMOSNESDISsm_struc(n)%obs_cdf,       &
          sm_current)

  endif

  !-------------------------------------------------------------------------
  !  Apply LSM based QC and screening of observations
  !-------------------------------------------------------------------------  
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_SMOSNESDISsmobsId)//char(0),n,k,OBS_state)

  call LIS_checkForValidObs(n,k,obsl,fnd,sm_current)

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
        call ESMF_AttributeSet(smField,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)
        
        call ESMF_AttributeSet(smField,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(smfield, "Unscaled Obs",&
             obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error in setting Unscaled Obs attribute')

     endif

     if(SMOSNESDISsm_struc(n)%useSsdevScal.eq.1) then
        call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')
        
        allocate(ssdev(LIS_rc%obs_ngrid(k)))
        ssdev = SMOSNESDISsm_struc(n)%ssdev_inp 

        if(SMOSNESDISsm_struc(n)%ntimes.eq.1) then 
           jj = 1
        else
           jj = LIS_rc%mo
        endif
        do t=1,LIS_rc%obs_ngrid(k)
           if(SMOSNESDISsm_struc(n)%obs_sigma(t,jj).gt.0) then 
              ssdev(t) = ssdev(t)*SMOSNESDISsm_struc(n)%model_sigma(t,jj)/&
                   SMOSNESDISsm_struc(n)%obs_sigma(t,jj)
              if(ssdev(t).lt.minssdev) then 
                 ssdev(t) = minssdev
              endif
           endif
        enddo
        
        if(LIS_rc%obs_ngrid(k).gt.0) then 
           call ESMF_AttributeSet(pertField,"Standard Deviation",&
                ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
        endif
        deallocate(ssdev)
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif

end subroutine read_SMOSNESDISsm

!BOP
! 
! !ROUTINE: read_SMOSNESDIS_data
! \label{read_SMOSNESDIS_data}
!
! !INTERFACE:
subroutine read_SMOSNESDIS_data(n, k, fname, smobs_ip, smtime_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use SMOSNESDISsm_Mod, only : SMOSNESDISsm_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                          :: smtime_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the SMOS NESDIS binary file and applies the data
!  quality flags to filter the data. !\normalsize

!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTSMOSNESDIS AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP
  integer*2           :: sm_raw(SMOSNESDISsm_struc(n)%nc,&
       SMOSNESDISsm_struc(n)%nr)
  integer*2           :: sm_raw_time(SMOSNESDISsm_struc(n)%nc,&
       SMOSNESDISsm_struc(n)%nr)
  real                :: sm_data(SMOSNESDISsm_struc(n)%nc*&
       SMOSNESDISsm_struc(n)%nr)
  integer*1           :: sm_data_hr(SMOSNESDISsm_struc(n)%nc,&
       SMOSNESDISsm_struc(n)%nr)
  integer*1           :: sm_data_mn(SMOSNESDISsm_struc(n)%nc,&
       SMOSNESDISsm_struc(n)%nr)
  integer*1           :: sm_data_hr1(SMOSNESDISsm_struc(n)%nc*&
       SMOSNESDISsm_struc(n)%nr)
  integer*1           :: sm_data_mn1(SMOSNESDISsm_struc(n)%nc*&
       SMOSNESDISsm_struc(n)%nr)
  real                :: sm_time(SMOSNESDISsm_struc(n)%nc*&
       SMOSNESDISsm_struc(n)%nr)
  logical*1           :: sm_data_b(SMOSNESDISsm_struc(n)%nc*&
       SMOSNESDISsm_struc(n)%nr)
  logical*1           :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer             :: hr_val, mn_val
  integer             :: julss
  integer             :: c,r,i,j,kk,ios
  integer             :: ftn,iret,igrib,nvars
  integer             :: param_num
  logical             :: var_found
  real                :: err, ql

  
  ftn = LIS_getNextUnitNumber()
  open(ftn,file=trim(fname),form='unformatted',access='direct',&
       recl=SMOSNESDISsm_struc(n)%nc*SMOSNESDISsm_struc(n)%nr*2,&
       convert='little_endian')
  read(ftn,rec=1) sm_raw
  read(ftn,rec=2) sm_raw_time
  call LIS_releaseUnitNumber(ftn)
  
  do r=1,SMOSNESDISsm_struc(n)%nr
     do c=1,SMOSNESDISsm_struc(n)%nc
        if((sm_raw(c,r)*0.0001.lt.0.01).or.sm_raw_time(c,r).lt.0) then 
           
           sm_data_hr(c,r) = -1
           sm_data_mn(c,r) = -1
        else
           sm_data_hr(c,r) = sm_raw_time(c,r)/100
           sm_data_mn(c,r) = sm_raw_time(c,r) - &
                sm_data_hr(c,r)*100
        endif
     enddo
  enddo

  do r=1,SMOSNESDISsm_struc(n)%nr
     do c=1,SMOSNESDISsm_struc(n)%nc
        if((sm_raw(c,r)*0.0001.lt.0.01).or.sm_raw_time(c,r).lt.0) then 
           sm_data(c+((SMOSNESDISsm_struc(n)%nr-r+1)-1)*&
                SMOSNESDISsm_struc(n)%nc) = LIS_rc%udef
           sm_data_hr1(c+((SMOSNESDISsm_struc(n)%nr-r+1)-1)*&
                SMOSNESDISsm_struc(n)%nc) = sm_data_hr(c,r)
           sm_data_mn1(c+((SMOSNESDISsm_struc(n)%nr-r+1)-1)*&
                SMOSNESDISsm_struc(n)%nc) = sm_data_mn(c,r)
           sm_data_b(c+((SMOSNESDISsm_struc(n)%nr-r+1)-1)*&
                SMOSNESDISsm_struc(n)%nc) = .false.
        else
           sm_data(c+((SMOSNESDISsm_struc(n)%nr-r+1)-1)*&
                SMOSNESDISsm_struc(n)%nc) = &
                sm_raw(c,r)*0.0001
           sm_data_hr1(c+((SMOSNESDISsm_struc(n)%nr-r+1)-1)*&
                SMOSNESDISsm_struc(n)%nc) =  sm_data_hr(c,r)
           sm_data_mn1(c+((SMOSNESDISsm_struc(n)%nr-r+1)-1)* &
                SMOSNESDISsm_struc(n)%nc) =  sm_data_mn(c,r)
           sm_data_b(c+((SMOSNESDISsm_struc(n)%nr-r+1)-1)*&
                SMOSNESDISsm_struc(n)%nc) = .true. 

        endif
     enddo
  enddo

  sm_time = LIS_rc%udef

  do r=1, SMOSNESDISsm_struc(n)%nr
     do c=1, SMOSNESDISsm_struc(n)%nc
        hr_val = sm_data_hr1(c+&
             (r-1)*SMOSNESDISsm_struc(n)%nc)
        mn_val =  sm_data_mn1(c+&
             (r-1)*SMOSNESDISsm_struc(n)%nc)
        call LIS_get_timeoffset_sec(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
             hr_val, mn_val, 0, julss)
        sm_time(c+(r-1)*SMOSNESDISsm_struc(n)%nc) = julss
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       SMOSNESDISsm_struc(n)%nc*SMOSNESDISsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       SMOSNESDISsm_struc(n)%rlat, SMOSNESDISsm_struc(n)%rlon, &
       SMOSNESDISsm_struc(n)%n11,  LIS_rc%udef, ios)

  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_time, smobs_b_ip, smtime_ip, &
       SMOSNESDISsm_struc(n)%nc*SMOSNESDISsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       SMOSNESDISsm_struc(n)%rlat, SMOSNESDISsm_struc(n)%rlon, &
       SMOSNESDISsm_struc(n)%n11,  LIS_rc%udef, ios)

end subroutine read_SMOSNESDIS_data


!BOP
! !ROUTINE: create_SMOSNESDISsm_filename
! \label{create_SMOSNESDISsm_filename}
! 
! !INTERFACE: 
subroutine create_SMOSNESDISsm_filename(ndir, yr, mo,da, hr,  filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da,hr
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the SMOSNESDIS filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the SMOSNESDIS soil moisture data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated SMOS NESDIS filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  logical           :: oper_flag
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  oper_flag = .false.

  if(yr.gt.2011) then
     if(yr.eq.2011) then 
        if(mo.eq.12) then 
           oper_flag = .true. 
        elseif(mo.eq.11) then 
           if(da.eq.30) then 
              oper_flag = .true. 
           endif
        endif
     else
        oper_flag = .true. 
     endif
  endif

  
  if(oper_flag) then 
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/SMOS_OPER_MIR_SMUDP2_SoilMoisture_'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.dat'
  else
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/SMOS_REPR_MIR_SMUDP2_SoilMoisture_'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.dat'
  endif

end subroutine create_SMOSNESDISsm_filename




