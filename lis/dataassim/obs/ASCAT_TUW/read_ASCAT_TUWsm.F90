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
! !ROUTINE: read_ASCAT_TUWsm
! \label{read_ASCAT_TUWsm}
!
! !REVISION HISTORY:
!  17 Jun 2010: Sujay Kumar; Updated for use with LPRM AMSRE Version 5. 
!  20 Sep 2012: Sujay Kumar; Updated to the NETCDF version of the data. 
!
! !INTERFACE: 
subroutine read_ASCAT_TUWsm(n, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use map_utils
  use LIS_pluginIndices
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use ASCAT_TUWsm_Mod, only : ASCAT_TUWsm_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
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
  real,  parameter       :: MAX_SM_VALUE=0.45, MIN_SM_VALUE=0.0001
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: smobsdir, fname_A, fname_D
  logical                :: alarmCheck, file_exists,file_status
  integer                :: t,c,r,i,j,p
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: smfield

  integer                :: gid(LIS_rc%ngrid(n))
  integer                :: assimflag(LIS_rc%ngrid(n))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: smobs(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                   :: sm_current(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                   :: dt
  real                   :: lon
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
  integer                :: fnd
  integer                :: lis_julss
  real                   :: smvalue
  real                   :: model_delta(LIS_rc%ngrid(n))
  real                   :: obs_delta(LIS_rc%ngrid(n))
  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  file_status = .false.
  data_upd = .false. 
!-------------------------------------------------------------------------
!   Read both ascending and descending passes at 0Z and then store
!   the overpass time as 1.30AM for the descending pass and 1.30PM 
!   for the ascending pass. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "ASCAT TUW read alarm")
  
  if(alarmCheck.or.ASCAT_TUWsm_struc(n)%startMode) then 
     ASCAT_TUWsm_struc(n)%startMode = .false.

     smobs = LIS_rc%udef
     
     call create_ASCAT_TUWsm_filenames(smobsdir, &
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, LIS_rc%hr, fname_A, fname_D)

     call read_ASCATTUW_data(n,fname_A, fname_D,smobs, &
          file_status)
       
  else
     smobs = LIS_rc%udef
  endif

  if(file_status)  then 
     call ESMF_StateGet(OBS_State,"Observation01",smfield,&
          rc=status)
     call LIS_verify(status, 'Error: StateGet Observation01')
     
     call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status, 'Error: FieldGet')
     
     
     obsl = LIS_rc%udef 
     do r=1, LIS_rc%lnr(n)
        do c=1, LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              obsl(LIS_domain(n)%gindex(c,r))=smobs(c+(r-1)*LIS_rc%lnc(n))
           endif
        enddo
     enddo
     
  
!lsm based qc
     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
          //trim(LIS_ASCAT_TUWsmobsId)//char(0),n, OBS_state)
     
     call ESMF_StateGet(OBS_State,"Observation01",smField,&
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
     fnd = 0 
     sm_current = LIS_rc%udef
     
     do r =1,LIS_rc%lnr(n)
        do c =1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r) .ne. -1)then
              sm_current(c,r) = obsl(LIS_domain(n)%gindex(c,r))
              if(sm_current(c,r).ne.LIS_rc%udef) then 
                 fnd = 1
              endif
           end if
        end do
     end do
     
!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     

     if(fnd.ne.0) then        
        call LIS_rescale_with_CDF_matching(    &
             n,                                   & 
             ASCAT_TUWsm_struc(n)%nbins,         & 
             MAX_SM_VALUE,                        & 
             MIN_SM_VALUE,                        & 
             ASCAT_TUWsm_struc(n)%model_xrange,  &
             ASCAT_TUWsm_struc(n)%obs_xrange,    &
             ASCAT_TUWsm_struc(n)%model_cdf,     &
             ASCAT_TUWsm_struc(n)%obs_cdf,       &
             sm_current)
        
     endif
     
     fnd = 0 
     data_upd_flag_local = .false.   
     do r =1,LIS_rc%lnr(n)
        do c =1,LIS_rc%lnc(n)
           if(sm_current(c,r).ne.LIS_rc%udef) then
              fnd = 1
           endif
        enddo
     enddo
     
     call ESMF_StateGet(OBS_State,"Observation01",smField,&
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
     obsl = LIS_rc%udef  
     
     if(fnd.eq.0) then
        obsl = LIS_rc%udef
     else
        do r =1,LIS_rc%lnr(n)
           do c =1,LIS_rc%lnc(n)
              if (LIS_domain(n)%gindex(c,r) .ne. -1)then
                 obsl(LIS_domain(n)%gindex(c,r))=sm_current(c,r)
              end if
           end do
        end do
     endif
     
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
        do t=1,LIS_rc%ngrid(n)
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
        
        if(LIS_rc%ngrid(n).gt.0) then 
           call ESMF_AttributeSet(smField,"Grid Number",&
                gid,itemCount=LIS_rc%ngrid(n),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(smField,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
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
  
end subroutine read_ASCAT_TUWsm

!BOP
! 
! !ROUTINE: read_ASCATTUW_data
! \label{read_ASCATTUW_data}
!
! !INTERFACE:
subroutine read_ASCATTUW_data(n, fname_A, fname_D, sm_data, data_status)
! 
! !USES:   
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use map_utils
  use ASCAT_TUWsm_Mod, only : ASCAT_TUWsm_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname_A
  character (len=*)             :: fname_D
  real                          :: sm_data(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical                       :: data_status

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the ASCATTUW grib2 file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  the estimated error is above a predefined threshold (the recommeded
!  value is 5%). 
! 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the ASCATTUW AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP
  real, parameter  :: err_threshold = 0.14
  logical          :: file_exists
  integer          :: ftn
  integer          :: n_data
  real, allocatable    :: lat(:)
  real, allocatable    :: lon(:)
  real, allocatable    :: sds(:)
  real, allocatable    :: err(:)
  real             :: col,row
  integer          :: c,r,i,t,iret
  logical*1        :: sm_data_b(ASCAT_TUWsm_struc(n)%nc*ASCAT_TUWsm_struc(n)%nr)
  logical*1        :: smobs_b_ip(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real             :: sm_ascat(ASCAT_TUWsm_struc(n)%nc*ASCAT_TUWsm_struc(n)%nr)

  sm_ascat = LIS_rc%udef
  inquire(file=fname_A, exist=file_exists)
  data_status = .false. 
  if(file_exists) then 
     data_status = .true. 
     ftn = LIS_getNextUnitNumber()
     write(LIS_logunit,*) 'Reading ',trim(fname_A)
     open(ftn,file=fname_A,form='unformatted')
     read(ftn) n_data

     if(n_data.gt.0) then 

        allocate(lat(n_data))
        allocate(lon(n_data))
        allocate(sds(n_data))
        allocate(err(n_data))

        read(ftn) lon
        read(ftn) lat
        read(ftn) sds
        read(ftn) err

        do i=1,n_data
           call latlon_to_ij(ASCAT_TUWsm_struc(n)%ascattuwproj,&
                lat(i),lon(i),col,row)
           c = nint(col)
           r = nint(row)
           if(err(i).le.err_threshold) then 
              sm_ascat(c+(r-1)*ASCAT_TUWsm_struc(n)%nc) = sds(i)
           endif
        enddo

        deallocate(lat)
        deallocate(lon)
        deallocate(sds)
        deallocate(err)
     endif
     call LIS_releaseUnitNumber(ftn)
  endif

  inquire(file=fname_D, exist=file_exists)

  if(file_exists) then 
     data_status = .true. 
     ftn = LIS_getNextUnitNumber()
     write(LIS_logunit,*) 'Reading ',fname_D
     open(ftn,file=fname_D,form='unformatted')
     read(ftn) n_data

     if(n_data.gt.0) then 

        allocate(lat(n_data))
        allocate(lon(n_data))
        allocate(sds(n_data))
        allocate(err(n_data))

        read(ftn) lon
        read(ftn) lat
        read(ftn) sds
        read(ftn) err

        do i=1,n_data
           call latlon_to_ij(ASCAT_TUWsm_struc(n)%ascattuwproj,&
                lat(i),lon(i),col,row)
           c = nint(col)
           r = nint(row)
           
           if(err(i).le.err_threshold) then 
              sm_ascat(c+(r-1)*ASCAT_TUWsm_struc(n)%nc) = sds(i)
           endif
        enddo

        deallocate(lat)
        deallocate(lon)
        deallocate(sds)
        deallocate(err)
     endif
     call LIS_releaseUnitNumber(ftn)
  endif

  sm_data_b = .false. 
  do t=1,ASCAT_TUWsm_struc(n)%nc*ASCAT_TUWsm_struc(n)%nr
     if(sm_ascat(t).ne.LIS_rc%udef) then 
        sm_data_b(t) = .true. 
     endif
  enddo
  
  call neighbor_interp(LIS_rc%gridDesc(n,:),&
       sm_data_b, sm_ascat, smobs_b_ip, sm_data,&
       ASCAT_TUWsm_struc(n)%nc*ASCAT_TUWsm_struc(n)%nr,&
       LIS_rc%lnc(n)*LIS_rc%lnr(n), & 
       LIS_domain(n)%lat, LIS_domain(n)%lon,&
       ASCAT_TUWsm_struc(n)%n11, LIS_rc%udef, iret)

!  open(100,file='test_ip.bin',form='unformatted')
!  write(100) sm_data
!  close(100)

end subroutine read_ASCATTUW_data


!BOP
! !ROUTINE: create_ASCAT_TUWsm_filenames
! \label{create_ASCAT_TUWsm_filenames}
! 
! !INTERFACE: 
subroutine create_ASCAT_TUWsm_filenames(ndir, yr, mo,da,hr,fname_A,fname_D)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: fname_A
  character(len=*)  :: fname_D
  integer           :: yr, mo, da,hr
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the RT SMOPS filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the RT SMOPS soil moisture data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated RT SMOPS  filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
 
  fname_A = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'/'//&
       'SDS_'//trim(fyr)//trim(fmo)//trim(fda)//&
       '.'//trim(fhr)//'_A.bin'
  fname_D = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'/'//&
       'SDS_'//trim(fyr)//trim(fmo)//trim(fda)//&
       '.'//trim(fhr)//'_D.bin'
  
end subroutine create_ASCAT_TUWsm_filenames




