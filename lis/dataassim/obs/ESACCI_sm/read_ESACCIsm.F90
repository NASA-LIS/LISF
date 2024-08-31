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
! !ROUTINE: read_ESACCIsm
! \label{read_ESACCIsm}
!
! !REVISION HISTORY:
!  01 Oct 2012: Sujay Kumar, Initial Specification
!  13 Jul 2016: Sujay Kumar, Updated the code to support DA in observation space
!
! !INTERFACE: 
subroutine read_ESACCIsm(n,k,  OBS_State, OBS_Pert_State)
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
  use ESACCI_sm_Mod, only : ESACCI_sm_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
! 
! This subroutine provides the data reader for the ESACCI
! soil moisture retrieval product. The routine also applies
! online bias correction and LSM based quality control. The 
! processed data is packaged into an ESMF State object for 
! later use with DA. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[k] index of the data assimilation instance
!  \item[OBS\_State] observations state
!  \item[OBS\_Pert\_State] observations perturbation state
!  \end{description}
!
!EOP
  real, parameter        ::  minssdev = 0.01
  real, parameter        ::  maxssdev = 0.11
  real,  parameter       :: MAX_SM_VALUE=0.45, MIN_SM_VALUE=0.0001
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN)          :: smobsdir, fname
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: smfield, pertfield

  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  real                   :: obs_unsc(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: smobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: sm_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: dt
  integer                :: fnd
  real                   :: lon
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
  real                   :: smvalue
  real, allocatable      :: ssdev(:)
  real                   :: model_delta(LIS_rc%obs_ngrid(k))
  real                   :: obs_delta(LIS_rc%obs_ngrid(k))
  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status, 'ESMF_AttributeGet failed in read_ESACCIsm')
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status, 'ESMF_AttributeGet failed in read_ESACCIsm')

  data_upd = .false. 
  obs_unsc = LIS_rc%udef
!-------------------------------------------------------------------------
!   Read the data at 0Z and store it. The reference time is defined
!   also as 0Z. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "ESACCI read alarm")
  
  if(alarmCheck.or.ESACCI_sm_struc(n)%startMode) then 
     ESACCI_sm_struc(n)%startMode = .false.

     ESACCI_sm_struc(n)%smobs = LIS_rc%udef
     smobs = LIS_rc%udef
     ESACCI_sm_struc(n)%smtime = -1
     
     call create_ESACCIsm_filename(smobsdir, &
          ESACCI_sm_struc(n)%version, ESACCI_sm_struc(n)%sensor,&
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname)

     inquire(file=fname,exist=file_exists)

     if(file_exists) then 
        
        write(LIS_logunit,*) 'Reading ',trim(fname)
        call read_ESACCI_data(n,k,fname,ESACCI_sm_struc(n)%version,smobs) ! NT: include version in reading

     endif

     ESACCI_sm_struc(n)%smobs  = LIS_rc%udef
     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           grid_index = LIS_obs_domain(n,k)%gindex(c,r)
           if(grid_index.ne.-1) then 
              if(smobs(c+(r-1)*LIS_rc%obs_lnc(k)).gt.0) then             
                 ESACCI_sm_struc(n)%smobs(c,r) = &
                      smobs(c+(r-1)*LIS_rc%obs_lnc(k))                 
                 lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                 
                 lhour = 12.0
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 ESACCI_sm_struc(n)%smtime(c,r) = gmt
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
! Assimilate at 12 z localtime

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           grid_index = c+(r-1)*LIS_rc%obs_lnc(k)

           dt = (LIS_rc%gmt - ESACCI_sm_struc(n)%smtime(c,r))*3600.0
           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
              sm_current(c,r) = & 
                   ESACCI_sm_struc(n)%smobs(c,r)
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
  if(LIS_rc%dascaloption(k).ne."none".and.fnd.ne.0) then  
     call LIS_rescale_with_CDF_matching(    &
          n,k,                              & 
          ESACCI_sm_struc(n)%nbins,         & 
          ESACCI_sm_struc(n)%ntimes,        & 
          MAX_SM_VALUE,                     & 
          MIN_SM_VALUE,                     & 
          ESACCI_sm_struc(n)%model_xrange,  &
          ESACCI_sm_struc(n)%obs_xrange,    &
          ESACCI_sm_struc(n)%model_cdf,     &
          ESACCI_sm_struc(n)%obs_cdf,       &
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
!-------------------------------------------------------------------------
!  Apply LSM-based QC and screening of observations
!-------------------------------------------------------------------------     

  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_ESACCIsmobsId)//char(0),n, k,OBS_state)

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
     call LIS_verify(status,&
          'ESMF_AttributeSet: Data Update Status failed in read_ESACCIsm')

     if(LIS_rc%obs_ngrid(k).gt.0) then 
        call ESMF_AttributeSet(smField,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status,&
             'ESMF_AttributeSet: Grid Number failed in read_ESACCIsm')          
  
        call ESMF_AttributeSet(smField,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status,&
             'ESMF_AttributeSet: Assimilation Flag failed in read_ESACCIsm')

        call ESMF_AttributeSet(smfield, "Unscaled Obs",&
             obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error in setting Unscaled Obs attribute')      
     endif

     if(ESACCI_sm_struc(n)%useSsdevScal.eq.1.and.&
          ESACCI_sm_struc(n)%ntimes.gt.1) then 

        call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')

        allocate(ssdev(LIS_rc%obs_ngrid(k)))
        ssdev = ESACCI_sm_struc(n)%ssdev_inp
        if(ESACCI_sm_struc(n)%ntimes.eq.1) then 
           jj = 1
        else
           jj = LIS_rc%mo
        endif

        do t=1,LIS_rc%obs_ngrid(k)
           if(ESACCI_sm_struc(n)%obs_sigma(t,jj).gt.0) then 
              ssdev(t) = ssdev(t)*ESACCI_sm_struc(n)%model_sigma(t,jj)/&
                   ESACCI_sm_struc(n)%obs_sigma(t,jj)
              if(ssdev(t).gt.maxssdev) ssdev(t) = maxssdev
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
     call LIS_verify(status,&
          'ESMF_AttributeSet: Data Update Status failed in read_ESACCIsm')     
  endif

end subroutine read_ESACCIsm


!BOP
! 
! !ROUTINE: read_ESACCI_data
! \label{read_ESACCI_data}
!
! !INTERFACE:
subroutine read_ESACCI_data(n, k, fname, version, smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,   only : LIS_verify
  use ESACCI_sm_Mod, only : ESACCI_sm_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real	                        :: version

! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the ESACCI NETCDF file and applies the data
!  quality flags to filter the data.The data quality flags provided
!  with the data is applied by excluding data masked for dense vegetation,
!  temperature below zero and lack of convergence of the algorithm. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the ESACCI AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  integer      :: sm(ESACCI_sm_struc(n)%ecvnc,ESACCI_sm_struc(n)%ecvnr)
  real         :: sm1(ESACCI_sm_struc(n)%ecvnc,ESACCI_sm_struc(n)%ecvnr)
  integer      :: flag(ESACCI_sm_struc(n)%ecvnc,ESACCI_sm_struc(n)%ecvnr)

  real         :: sm_combined(ESACCI_sm_struc(n)%ecvnc,ESACCI_sm_struc(n)%ecvnr)
  real         :: sm_data(ESACCI_sm_struc(n)%ecvnc*ESACCI_sm_struc(n)%ecvnr)
  logical*1    :: sm_data_b(ESACCI_sm_struc(n)%ecvnc*ESACCI_sm_struc(n)%ecvnr)
  logical*1    :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  integer      :: c,r,i,j
  real         :: rlat,rlon,ri,rj
  integer      :: nid
  integer      :: smid, flagid
  integer      :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LIS_verify(ios,'Error opening file '//trim(fname))
  
  ios = nf90_inq_varid(nid, 'sm',smid)
  call LIS_verify(ios, 'Error nf90_inq_varid: sm')
  
  ios = nf90_inq_varid(nid, 'flag',flagid)
  call LIS_verify(ios, 'Error nf90_inq_varid: flag')
  
  !values
  if(version .lt. 3) then
      ios = nf90_get_var(nid, smid, sm)
  else
      ios = nf90_get_var(nid, smid, sm1)
  endif
  call LIS_verify(ios, 'Error nf90_get_var: sm')
  
  ios = nf90_get_var(nid, flagid,flag)
  call LIS_verify(ios, 'Error nf90_get_var: flag')
  
  ios = nf90_close(ncid=nid)
  call LIS_verify(ios,'Error closing file '//trim(fname))

  do r=1, ESACCI_sm_struc(n)%ecvnr
     do c=1, ESACCI_sm_struc(n)%ecvnc
!------------------------------------------------------------------------
! All data flagged for snow coverage or temperature below zero (flag=1), 
! dense vegetation (flag=2) and no convergence in the ESACCI algorithm 
! (flag =3) and undefined values are masked out. 
!------------------------------------------------------------------------
        if(version .lt. 3) then
            if(flag(c,r).ne.0.or.sm(c,r).le.0) then 
               sm_combined(c,ESACCI_sm_struc(n)%ecvnr-r+1) = LIS_rc%udef
            else
              sm_combined(c,ESACCI_sm_struc(n)%ecvnr-r+1) = sm(c,r)*0.0001
            endif
        else
            if(flag(c,r).ne.0.or.sm1(c,r).le.0) then ! NT: checking sm and sm1 separately
               sm_combined(c,ESACCI_sm_struc(n)%ecvnr-r+1) = LIS_rc%udef
            else
              sm_combined(c,ESACCI_sm_struc(n)%ecvnr-r+1) = sm1(c,r)
           endif
        endif
     enddo
  enddo

  do r=1, ESACCI_sm_struc(n)%ecvnr
     do c=1, ESACCI_sm_struc(n)%ecvnc
        sm_data(c+(r-1)*ESACCI_sm_struc(n)%ecvnc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LIS_rc%udef) then 
           sm_data_b(c+(r-1)*ESACCI_sm_struc(n)%ecvnc) = .true. 
        else
           sm_data_b(c+(r-1)*ESACCI_sm_struc(n)%ecvnc) = .false.
        endif
        if(sm_combined(c,r).gt.0.5) then 
           sm_combined(c,r) = LIS_rc%udef
           sm_data_b(c+(r-1)*ESACCI_sm_struc(n)%ecvnc) = .false.
        endif
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the DA observation space
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       ESACCI_sm_struc(n)%ecvnc*ESACCI_sm_struc(n)%ecvnr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       ESACCI_sm_struc(n)%rlat, ESACCI_sm_struc(n)%rlon, &
       ESACCI_sm_struc(n)%n11,LIS_rc%udef, ios)

#endif
  
end subroutine read_ESACCI_data

!BOP
! !ROUTINE: create_ESACCIsm_filename
! \label{create_ESACCIsm_filename}
! 
! !INTERFACE: 
subroutine create_ESACCIsm_filename(ndir, version, sensor, yr, mo,da, filename)
! !USES:   
  use LIS_logMod

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  real              :: version
  character(len=*)  :: sensor
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the ESACCI filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the ESACCI soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated ESACCI filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  character (len=3) :: cversion3
  character (len=4) :: cversion4
  character (len=8) :: sensortxt !20220622 Pang 
 
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  if (version .lt. 10.0) then
     write(unit=cversion3, fmt='(f3.1)') version
  else
     ! For future version, version number > 10, e.g., 10.2
     write(unit=cversion4, fmt='(f4.1)') version
  endif
 
  if(sensor == 'passive') then
      sensortxt = 'PASSIVE'
  elseif(sensor == 'active') then
      sensortxt = 'ACTIVE'
  elseif(sensor == 'combined') then
      sensortxt = 'COMBINED'
  else
      write(LIS_logunit,*) "[ERR] Invalid ESA CCI soil moisture sensor type was chosen."
      write(LIS_logunit,*) "[ERR] Please choose either 'passive', 'active', or 'combined'."
      call LIS_endrun
  endif
 
  if(version.eq.1) then 
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/ESACCI-L3S_SOILMOISTURE-SSMV-MERGED-' &
          //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv00.1.nc'
     
  elseif(version.eq.2) then 
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-' & 
          //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv02.0.nc'
  elseif(version.eq.2.2) then 
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-' & 
          //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv02.2.nc'
  else
     ! NT: for versions after 2.2
     if (version .lt. 10.0) then
         filename = trim(ndir)//'/'//trim(fyr)//&
              '/ESACCI-SOILMOISTURE-L3S-SSMV-'//trim(sensortxt)//'-' &  !20220622 Pang
              //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv0'//cversion3//'.nc'
     else
         filename = trim(ndir)//'/'//trim(fyr)//&
              '/ESACCI-SOILMOISTURE-L3S-SSMV-'//trim(sensortxt)//'-' & !20220622 Pang
              //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv'//cversion4//'.nc'
     endif
  endif
end subroutine create_ESACCIsm_filename



