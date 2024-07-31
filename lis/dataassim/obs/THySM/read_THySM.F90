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
! !ROUTINE: read_THySM
! \label{read_THySM}
!
! !REVISION HISTORY:
!  28 Mar 2021    Sujay Kumar; initial specification
!
! !INTERFACE: 
subroutine read_THySM(n, k, OBS_State, OBS_Pert_State)
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
  use THySM_Mod, only : THySM_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
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
!  reads the THySM observations 
!  and rescales the data to the model soil moisture data
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
  real,  parameter       :: MAX_SM_VALUE=10.0, MIN_SM_VALUE=0.0001
  integer                :: status
  integer                :: ftn
  integer                :: latid,lonid,smid
  real                   :: lat(THySM_struc(n)%nr)
  real                   :: lon(THySM_struc(n)%nc)
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p,jj,c1,r1
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: smfield, pertField
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  real                   :: obs_unsc(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: sm_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: dt
  real                   :: lon_value
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
  integer                :: fnd
  real, allocatable      :: ssdev(:)
  integer                :: lis_julss
  real                   :: smvalue
  real                   :: sm_file(THySM_struc(n)%nc,THySM_struc(n)%nr)
  real                   :: sm_inp(THySM_struc(n)%nc*THySM_struc(n)%nr)
  logical*1              :: sm_b_inp(THySM_struc(n)%nc*THySM_struc(n)%nr)
  real                   :: sm_out(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  logical*1              :: sm_b_out(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
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
  sm_out = LIS_rc%udef

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "THySM read alarm")

  if(alarmCheck.or.THySM_struc(n)%startMode) then 
     THySM_struc(n)%startMode = .false.

     call create_THySM_filename(smobsdir, &
          'PM',&
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname)
     
     inquire(file=fname,exist=file_exists)
     
     if(file_exists) then 
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        status = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=ftn)
        call LIS_verify(status,'Error opening file '//trim(fname))
        
        status = nf90_inq_varid(ftn, 'Band1',smid)
        call LIS_verify(status, 'Error nf90_inq_varid: sm')
        
        !values
        status = nf90_get_var(ftn, smid, sm_file, &
             start=(/THySM_struc(n)%c1,THySM_struc(n)%r1/),&
             count=(/THySM_struc(n)%nc,THySM_struc(n)%nr/))
        call LIS_verify(status, 'Error nf90_get_var: sm')
                
        status = nf90_close(ncid=ftn)
        call LIS_verify(status,'Error closing file '//trim(fname))
#endif

        sm_inp = LIS_rc%udef
        sm_b_inp  = .false. 

        do r=1,THySM_struc(n)%nr
           do c=1,THySM_struc(n)%nc 
              if(.not.isNaN(sm_file(c,r)).and.sm_file(c,r).gt.0) then 
                 sm_inp(c+(r-1)*THySM_struc(n)%nc) = & 
                      sm_file(c,r)
                 sm_b_inp(c+(r-1)*THySM_struc(n)%nc) = & 
                      .true. 
              endif
           enddo
        enddo

        call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
             sm_b_inp, sm_inp, sm_b_out, sm_out, &
             THySM_struc(n)%nc*THySM_struc(n)%nr, &
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             THySM_struc(n)%rlat, THySM_struc(n)%rlon, &
             THySM_struc(n)%n11,LIS_rc%udef, status)
        
     else
        write(LIS_logunit,*) '[WARN] Missing THySM file: ',trim(fname)
     endif
     
     THySM_struc(n)%smobs  = LIS_rc%udef
     THySM_struc(n)%smtime = -1
     
     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           grid_index = LIS_obs_domain(n,k)%gindex(c,r)
           if(grid_index.ne.-1) then 
              if(sm_out(c+(r-1)*LIS_rc%obs_lnc(k)).gt.0) then 
                 THySM_struc(n)%smobs(c,r) = &
                         sm_out(c+(r-1)*LIS_rc%obs_lnc(k))                 
                 lon_value = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                 lhour = 6.0
                 call LIS_localtime2gmt (gmt,lon_value,lhour,zone)
                 THySM_struc(n)%smtime(c,r) = gmt
              endif
           endif
        enddo
     enddo
     call create_THySM_filename(smobsdir, &
          'AM',&
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname)
     
     inquire(file=fname,exist=file_exists)

     if(file_exists) then 
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        status = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=ftn)
        call LIS_verify(status,'Error opening file '//trim(fname))
        
        status = nf90_inq_varid(ftn, 'Band1',smid)
        call LIS_verify(status, 'Error nf90_inq_varid: sm')
        
        !values
        status = nf90_get_var(ftn, smid, sm_file, &
             start=(/THySM_struc(n)%c1,THySM_struc(n)%r1/),&
             count=(/THySM_struc(n)%nc,THySM_struc(n)%nr/))
        call LIS_verify(status, 'Error nf90_get_var: sm')
                
        status = nf90_close(ncid=ftn)
        call LIS_verify(status,'Error closing file '//trim(fname))
#endif
        sm_inp = -9999.0
        sm_b_inp = .false.
        
        do r=1,THySM_struc(n)%nr
           do c=1,THySM_struc(n)%nc 
              if(.not.isNaN(sm_file(c,r)).and.sm_file(c,r).gt.0) then 
                 sm_inp(c+(r-1)*THySM_struc(n)%nc) = & 
                      sm_file(c,r)
                 sm_b_inp(c+(r-1)*THySM_struc(n)%nc) = & 
                      .true. 
              endif
           enddo
        enddo

        call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
             sm_b_inp, sm_inp, sm_b_out, sm_out, &
             THySM_struc(n)%nc*THySM_struc(n)%nr, &
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             THySM_struc(n)%rlat, THySM_struc(n)%rlon, &
             THySM_struc(n)%n11,LIS_rc%udef, status)
        
     else
        write(LIS_logunit,*) '[WARN] Missing THySM file: ',trim(fname)
     endif
     
     
     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           grid_index = LIS_obs_domain(n,k)%gindex(c,r)
           if(grid_index.ne.-1) then 
              if(sm_out(c+(r-1)*LIS_rc%obs_lnc(k)).gt.0.and.&
                   THySM_struc(n)%smobs(c,r).eq.-9999.0) then 
                 THySM_struc(n)%smobs(c,r) = &
                         sm_out(c+(r-1)*LIS_rc%obs_lnc(k))                 
                 lon_value = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                 lhour = 18.0
                 call LIS_localtime2gmt (gmt,lon_value,lhour,zone)
                 THySM_struc(n)%smtime(c,r) = gmt
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

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1.and.&
             THySM_struc(n)%smtime(c,r).gt.0) then 
           grid_index = c+(r-1)*LIS_rc%obs_lnc(k)

           dt = (LIS_rc%gmt - THySM_struc(n)%smtime(c,r))*3600.0
           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
              sm_current(c,r) = & 
                   THySM_struc(n)%smobs(c,r)
              if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                 obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = &
                      sm_current(c,r)
              endif
              if(sm_current(c,r).ne.LIS_rc%udef) then 
                 fnd = 1
              endif
           endif
        endif
     enddo
  enddo

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     

  if(LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then
     call LIS_rescale_with_CDF_matching(     &
          n,k,                               & 
          THySM_struc(n)%nbins,         & 
          THySM_struc(n)%ntimes,        & 
          MAX_SM_VALUE,                      & 
          MIN_SM_VALUE,                      & 
          THySM_struc(n)%model_xrange,  &
          THySM_struc(n)%obs_xrange,    &
          THySM_struc(n)%model_cdf,     &
          THySM_struc(n)%obs_cdf,       &
          sm_current)
  elseif(LIS_rc%dascaloption(k).eq."Linear scaling".and.fnd.ne.0) then
     call LIS_rescale_with_linear_scaling(    &
          n,                                   & 
          k,                                   & 
          THySM_struc(n)%nbins,         & 
          THySM_struc(n)%ntimes,        & 
          THySM_struc(n)%obs_xrange,    &
          THySM_struc(n)%obs_cdf,       &
          sm_current)
     
  endif
  
  obsl = LIS_rc%udef 
  do r=1, LIS_rc%obs_lnr(k)
     do c=1, LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           obsl(LIS_obs_domain(n,k)%gindex(c,r))=sm_current(c,r)
        endif
     enddo
  enddo
  write(100) sm_current
  !-------------------------------------------------------------------------
  !  Apply LSM based QC and screening of observations
  !-------------------------------------------------------------------------  
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_THySMId)//char(0),n,k,OBS_state)

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

     if(LIS_rc%dascaloption(k).eq."CDF matching") then 
        if(THySM_struc(n)%useSsdevScal.eq.1) then
           call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                rc=status)
           call LIS_verify(status, 'Error: StateGet Observation01')
           
           allocate(ssdev(LIS_rc%obs_ngrid(k)))
           ssdev = THySM_struc(n)%ssdev_inp 
           
           if(THySM_struc(n)%ntimes.eq.1) then 
              jj = 1
           else
              jj = LIS_rc%mo
           endif
           do t=1,LIS_rc%obs_ngrid(k)
              if(THySM_struc(n)%obs_sigma(t,jj).gt.0) then 
                 ssdev(t) = ssdev(t)*THySM_struc(n)%model_sigma(t,jj)/&
                      THySM_struc(n)%obs_sigma(t,jj)
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
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif
end subroutine read_THySM


!BOP
! !ROUTINE: create_THySM_filename
! \label{create_THySM_filename}
! 
! !INTERFACE: 
subroutine create_THySM_filename(ndir, overpass, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: overpass
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the THySM filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the THySM data directory
!  \item[overpass]  AM or PM overpass
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated THySM filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  filename=trim(ndir)//'/'//trim(fyr)//&
       '/SMAP-HYB-1KM-DAILY_'//trim(fyr)//'.'//&
       trim(fmo)//'.'//trim(fda)//'_'//trim(overpass)//'.nc'
  
end subroutine create_THySM_filename





