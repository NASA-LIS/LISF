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
! !ROUTINE: read_GCOMW_AMSR2L3sm
! \label{read_GCOMW_AMSR2L3sm}
!
! !REVISION HISTORY:
!  12 Jan 15: Sujay Kumar; Initial version
!
! !INTERFACE: 
subroutine read_GCOMW_AMSR2L3sm(n, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use map_utils
  use LIS_pluginIndices
  use LIS_constantsMod,    only : LIS_CONST_PATH_LEN
  use GCOMW_AMSR2L3sm_Mod, only : GCOMW_AMSR2L3sm_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the AMSR2 soil moisture observations 
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
  real, parameter        ::  minssdev = 0.01
  real, parameter        ::  maxssdev = 0.11
  real                   :: MAX_SM_VALUE, MIN_SM_VALUE
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: smobsdir, fname_A, fname_D
  logical                :: alarmCheck
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: smfield, pertField

  integer                :: gid(LIS_rc%ngrid(n))
  integer                :: assimflag(LIS_rc%ngrid(n))
  real                   :: obs_unsc(LIS_rc%ngrid(n))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: smobs_A(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                   :: smobs_D(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                   :: sm_current(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                   :: dt
  real                   :: lon
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
  integer                :: fnd
  real                   :: smvalue
  real, allocatable      :: ssdev(:)
  real                   :: model_delta(LIS_rc%ngrid(n))
  real                   :: obs_delta(LIS_rc%ngrid(n))
  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  MAX_SM_VALUE=0.45
  MIN_SM_VALUE=0.0001

  data_upd = .false. 
  obs_unsc = LIS_rc%udef
!-------------------------------------------------------------------------
!   Read both ascending and descending passes at 0Z and then store
!   the overpass time as 1.30AM for the descending pass and 1.30PM 
!   for the ascending pass. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "AMSR2(GCOMW) read alarm")
  
  if(alarmCheck.or.GCOMW_AMSR2L3sm_struc(n)%startMode) then 
     GCOMW_AMSR2L3sm_struc(n)%startMode = .false.

     GCOMW_AMSR2L3sm_struc(n)%smobs = LIS_rc%udef
     smobs_A = LIS_rc%udef
     smobs_D = LIS_rc%udef             
     GCOMW_AMSR2L3sm_struc(n)%smtime = -1

     call create_GCOMW_AMSR2L3sm_filename(smobsdir, 'A',&
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname_A)

     call create_GCOMW_AMSR2L3sm_filename(smobsdir, 'D',&
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname_D)

     call read_AMSR2_data(n,fname_A, smobs_A)
     call read_AMSR2_data(n,fname_D, smobs_D)

     GCOMW_AMSR2L3sm_struc(n)%smobs  = LIS_rc%udef
     GCOMW_AMSR2L3sm_struc(n)%smtime = -1


     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           grid_index = LIS_domain(n)%gindex(c,r)
           if(grid_index.ne.-1) then 

              if(smobs_A(c+(r-1)*LIS_rc%lnc(n)).ne.-9999.0) then               
                 GCOMW_AMSR2L3sm_struc(n)%smobs(c,r) = &
                      smobs_A(c+(r-1)*LIS_rc%lnc(n))                 
                 lon = LIS_domain(n)%grid(grid_index)%lon
                 lhour = 13.5 
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 GCOMW_AMSR2L3sm_struc(n)%smtime(c,r) = gmt
              endif

              if(smobs_D(c+(r-1)*LIS_rc%lnc(n)).ne.-9999.0) then               
                 GCOMW_AMSR2L3sm_struc(n)%smobs(c,r) = &
                      smobs_D(c+(r-1)*LIS_rc%lnc(n))

                 lon = LIS_domain(n)%grid(grid_index)%lon
                 lhour = 1.5
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 GCOMW_AMSR2L3sm_struc(n)%smtime(c,r) = gmt
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

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(GCOMW_AMSR2L3sm_struc(n)%smtime(c,r).ge.0) then 
              dt = (LIS_rc%gmt - GCOMW_AMSR2L3sm_struc(n)%smtime(c,r))*3600.0
              if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
                 sm_current(c,r) = & 
                      GCOMW_AMSR2L3sm_struc(n)%smobs(c,r)
                 fnd = 1
              endif
           endif
        endif
     enddo
  enddo

  obsl = LIS_rc%udef 
  if(fnd.ne.0) then 
     do r=1, LIS_rc%lnr(n)
        do c=1, LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              obsl(LIS_domain(n)%gindex(c,r))=sm_current(c,r)
           endif
        enddo
     enddo
  endif

  !lsm based qc
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),n, OBS_state)

  call ESMF_StateGet(OBS_State,"Observation01",smField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)

  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           sm_current(c,r) = obsl(LIS_domain(n)%gindex(c,r))
           obs_unsc(LIS_domain(n)%gindex(c,r)) = &
                obsl(LIS_domain(n)%gindex(c,r))
        end if
     end do
  end do

  !-------------------------------------------------------------------------
  !  Transform data to the LSM climatology using a CDF-scaling approach
  !-------------------------------------------------------------------------     

  if(GCOMW_AMSR2L3sm_struc(n)%scal.ne.0.and.fnd.ne.0) then        
     call LIS_rescale_with_CDF_matching(    &
          n,                                   & 
          GCOMW_AMSR2L3sm_struc(n)%nbins,         & 
          GCOMW_AMSR2L3sm_struc(n)%ntimes,        & 
          MAX_SM_VALUE,                        & 
          MIN_SM_VALUE,                        & 
          GCOMW_AMSR2L3sm_struc(n)%model_xrange,  &
          GCOMW_AMSR2L3sm_struc(n)%obs_xrange,    &
          GCOMW_AMSR2L3sm_struc(n)%model_cdf,     &
          GCOMW_AMSR2L3sm_struc(n)%obs_cdf,       &
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

        call ESMF_AttributeSet(smfield, "Unscaled Obs",&
             obs_unsc, itemCount=LIS_rc%ngrid(n), rc=status)
        call LIS_verify(status, 'Error in setting Unscaled Obs attribute')
     endif

     if(GCOMW_AMSR2L3sm_struc(n)%useSsdevScal.eq.1.and.&
          GCOMW_AMSR2L3sm_struc(n)%ntimes.gt.1) then 

        call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')

        allocate(ssdev(LIS_rc%ngrid(n)))
        ssdev = GCOMW_AMSR2L3sm_struc(n)%ssdev_inp
        if(GCOMW_AMSR2L3sm_struc(n)%ntimes.eq.1) then 
           jj = 1
        else
           jj = LIS_rc%mo
        endif

        do t=1,LIS_rc%ngrid(n)
           if(GCOMW_AMSR2L3sm_struc(n)%obs_sigma(t,jj).gt.0) then 
              ssdev(t) = ssdev(t)*GCOMW_AMSR2L3sm_struc(n)%model_sigma(t,jj)/&
                   GCOMW_AMSR2L3sm_struc(n)%obs_sigma(t,jj)
              if(ssdev(t).gt.maxssdev) ssdev(t) = maxssdev
              if(ssdev(t).lt.minssdev) then 
                 ssdev(t) = minssdev
              endif
           endif
        enddo

        if(LIS_rc%ngrid(n).gt.0) then 
           call ESMF_AttributeSet(pertField,"Standard Deviation",&
                ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
           call LIS_verify(status)
        endif
        deallocate(ssdev)
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif
end subroutine read_GCOMW_AMSR2L3sm

!BOP
! 
! !ROUTINE: read_AMSR2_data
! \label{read_AMSR2_data}
!
! !INTERFACE:
subroutine read_AMSR2_data(n, fname,smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod
  use map_utils,    only : latlon_to_ij
  use GCOMW_AMSR2L3sm_Mod, only : GCOMW_AMSR2L3sm_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname
  real                          :: smobs_ip(LIS_rc%lnc(n)*LIS_rc%lnr(n))


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the LPRM NETCDF file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  land surface temperature is below freezing, if rain is present, if 
!  RFI is present, if residual error is above 0.5 or if optical depth
!  is above 0.8. Finally the routine combines both the C-band and X-band
!  retrievals. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the LPRM AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  real                        :: sm(1,GCOMW_AMSR2L3sm_struc(n)%amsr2nc,&
       GCOMW_AMSR2L3sm_struc(n)%amsr2nr)
  real                        :: time(GCOMW_AMSR2L3sm_struc(n)%amsr2nc,&
       GCOMW_AMSR2L3sm_struc(n)%amsr2nr)
  integer                     :: smtime(GCOMW_AMSR2L3sm_struc(n)%amsr2nc,&
       GCOMW_AMSR2L3sm_struc(n)%amsr2nr)
  real                        :: sm_combined(GCOMW_AMSR2L3sm_struc(n)%amsr2nc,&
       GCOMW_AMSR2L3sm_struc(n)%amsr2nr)

  real                        :: sm_data(GCOMW_AMSR2L3sm_struc(n)%amsr2nc*&
       GCOMW_AMSR2L3sm_struc(n)%amsr2nr)
  logical*1                   :: sm_data_b(GCOMW_AMSR2L3sm_struc(n)%amsr2nc*&
       GCOMW_AMSR2L3sm_struc(n)%amsr2nr)
  logical*1                   :: smobs_b_ip(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  logical                     :: file_exists
  integer                     :: c,r,c1,i,j
  real                        :: rlat,rlon,ri,rj
  integer                     :: nid
  integer                     :: smId
  integer                     :: timeid
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  smtime = -1.0
  sm_combined = LIS_rc%udef

  inquire(file=fname, exist=file_exists) 
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading ',trim(fname)
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(fname))
     
     ios = nf90_inq_varid(nid, 'Geophysical Data',smid)
     call LIS_verify(ios, 'Error nf90_inq_varid: Geophysical Data')
     
     ios = nf90_inq_varid(nid, 'Time Information',timeId)
     call LIS_verify(ios, 'Error nf90_inq_varid: Time Information')
     
     !values
     ios = nf90_get_var(nid, smid, sm)
     call LIS_verify(ios, 'Error nf90_get_var: SM_Flags')
     
     ios = nf90_get_var(nid, timeid, time)
     call LIS_verify(ios, 'Error nf90_get_var: time')
     
     ios = nf90_close(ncid=nid)
     call LIS_verify(ios,'Error closing file '//trim(fname))
     
     do r=1, GCOMW_AMSR2L3sm_struc(n)%amsr2nr
        do c=1, GCOMW_AMSR2L3sm_struc(n)%amsr2nc
           c1 = c + GCOMW_AMSR2L3sm_struc(n)%amsr2nc/2 - 1
           if(c1.gt.GCOMW_AMSR2L3sm_struc(n)%amsr2nc) then 
              c1 = c - GCOMW_AMSR2L3sm_struc(n)%amsr2nc/2 -1
           endif
           if(time(c,r).lt.0.or.sm(1,c,r).le.0.001) then 
              sm_combined(c1,GCOMW_AMSR2L3sm_struc(n)%amsr2nr-r+1) = LIS_rc%udef
           else
              sm_combined(c1,GCOMW_AMSR2L3sm_struc(n)%amsr2nr-r+1) = &
                   (sm(1,c,r)*0.1)/100.0
              smtime(c1,GCOMW_AMSR2L3sm_struc(n)%amsr2nr-r+1) = time(c,r)
           endif
        enddo
     enddo
  endif

  do r=1, GCOMW_AMSR2L3sm_struc(n)%amsr2nr
     do c=1, GCOMW_AMSR2L3sm_struc(n)%amsr2nc
        sm_data(c+(r-1)*GCOMW_AMSR2L3sm_struc(n)%amsr2nc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LIS_rc%udef) then 
           sm_data_b(c+(r-1)*GCOMW_AMSR2L3sm_struc(n)%amsr2nc) = .true. 
        else
           sm_data_b(c+(r-1)*GCOMW_AMSR2L3sm_struc(n)%amsr2nc) = .false.
        endif
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  if(LIS_isatAfinerResolution(n,GCOMW_AMSR2L3sm_struc(n)%datares)) then   
     call bilinear_interp(LIS_rc%gridDesc(n,:),&
          sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
          GCOMW_AMSR2L3sm_struc(n)%amsr2nc*GCOMW_AMSR2L3sm_struc(n)%amsr2nr, &
          LIS_rc%lnc(n)*LIS_rc%lnr(n), &
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          GCOMW_AMSR2L3sm_struc(n)%w11, GCOMW_AMSR2L3sm_struc(n)%w12, &
          GCOMW_AMSR2L3sm_struc(n)%w21, GCOMW_AMSR2L3sm_struc(n)%w22, &
          GCOMW_AMSR2L3sm_struc(n)%n11, GCOMW_AMSR2L3sm_struc(n)%n12, &
          GCOMW_AMSR2L3sm_struc(n)%n21, GCOMW_AMSR2L3sm_struc(n)%n22, &
          LIS_rc%udef, ios)
  else
     call upscaleByAveraging(&
          GCOMW_AMSR2L3sm_struc(n)%amsr2nc*GCOMW_AMSR2L3sm_struc(n)%amsr2nr, &
          LIS_rc%lnc(n)*LIS_rc%lnr(n), LIS_rc%udef, &
          GCOMW_AMSR2L3sm_struc(n)%n11,sm_data_b, &
          sm_data, smobs_b_ip,smobs_ip)
  endif
#endif
  
end subroutine read_AMSR2_data


!BOP
! !ROUTINE: create_GCOMW_AMSR2L3sm_filename
! \label{create_GCOMW_AMSR2L3sm_filename}
! 
! !INTERFACE: 
subroutine create_GCOMW_AMSR2L3sm_filename(ndir, path, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: path
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the AMSR2 cmd based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the AMSR2 soil moisture directory
!  \item[path] name of the sensor path (A-ascending, D-descending)
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated AMSR2 filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  if(path.eq.'A') then 
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
          '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
          '_01D_EQMA_L3SGSMCHA1110100.h5'         
  else     
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
          '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
          '_01D_EQMD_L3SGSMCHA1110100.h5'      
  endif
  
end subroutine create_GCOMW_AMSR2L3sm_filename




