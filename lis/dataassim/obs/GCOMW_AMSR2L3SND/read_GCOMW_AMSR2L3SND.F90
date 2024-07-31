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
! !ROUTINE: read_GCOMW_AMSR2L3SND
! \label{read_GCOMW_AMSR2L3SND}
!
! !REVISION HISTORY:
!  12 Jan 15: Sujay Kumar; Initial version
!
! !INTERFACE: 
subroutine read_GCOMW_AMSR2L3SND(n, k, OBS_State, OBS_Pert_State)
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
  use GCOMW_AMSR2L3SND_Mod, only : GCOMW_AMSR2L3SND_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
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
  integer,         parameter    :: IMSnc = 1500,IMSnr = 375
  integer                :: ftn,status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: sndobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname_A, fname_D,imsfile, MODISfile,fname
  logical                :: alarmCheck
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: sndfield, pertField
  logical*1              :: lo(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: sndobs_A(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: sndobs_D(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: snd_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: dt
  real                   :: lon
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
  integer                :: fnd
  real                   :: sndvalue
  real, allocatable      :: ssdev(:)
  logical                :: file_exists
  real                   :: IMSdata(IMSnc*IMSnr)
  logical*1              :: IMS_li(IMSnc*IMSnr)
  real                   :: IMSdata_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: MODISdata_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                     :: ios

    
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sndobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 
!-------------------------------------------------------------------------
!   Read both ascending and descending passes at 0Z and then store
!   the overpass time as 1.30AM for the descending pass and 1.30PM 
!   for the ascending pass. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "AMSR2(GCOMW) read alarm")
  
  if(alarmCheck.or.GCOMW_AMSR2L3SND_struc(n)%startMode) then 
     
     GCOMW_AMSR2L3SND_struc(n)%IMSdata_obs = -9999.0
     GCOMW_AMSR2L3SND_struc(n)%MODISdata_obs = -9999.0

     GCOMW_AMSR2L3SND_struc(n)%startMode = .false.

     GCOMW_AMSR2L3SND_struc(n)%sndobs = LIS_rc%udef

     if(GCOMW_AMSR2L3SND_struc(n)%bc_version.eq.1) then 

        sndobs_A = LIS_rc%udef
        call create_GCOMW_AMSR2L3SND_BC_filename(sndobsdir, &
             LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname) 
        call read_AMSR2snd_bc_data(n,k, fname, sndobs_A)       

        GCOMW_AMSR2L3SND_struc(n)%sndobs  = LIS_rc%udef
        GCOMW_AMSR2L3SND_struc(n)%sndtime = -1
        
        
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              grid_index = LIS_obs_domain(n,k)%gindex(c,r)
              if(grid_index.ne.-1) then 
                 
                 if(sndobs_A(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then        
                    GCOMW_AMSR2L3SND_struc(n)%sndobs(c,r) = &
                         sndobs_A(c+(r-1)*LIS_rc%obs_lnc(k))                 
                    lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))

                    lhour = 14.0
                    call LIS_localtime2gmt (gmt,lon,lhour,zone)
                    GCOMW_AMSR2L3SND_struc(n)%sndtime(c,r) = gmt
                 endif
                 
              endif
           enddo
        enddo
        if(GCOMW_AMSR2L3SND_struc(n)%usr_input_mask.eq.1) then 
           
           do r=1,LIS_rc%obs_lnr(k)
              do c=1,LIS_rc%obs_lnc(k)
                 grid_index = LIS_obs_domain(n,k)%gindex(c,r)
                 if(grid_index.ne.-1) then 
                    if(GCOMW_AMSR2L3SND_struc(n)%input_mask(c,r).gt.0) then
                       GCOMW_AMSR2L3SND_struc(n)%sndobs(c,r) = LIS_rc%udef
                    endif
                 endif
              enddo
           enddo
        endif
        
     else
        sndobs_A = LIS_rc%udef
        sndobs_D = LIS_rc%udef             
        GCOMW_AMSR2L3SND_struc(n)%sndtime = -1
        
        call create_GCOMW_AMSR2L3SND_filename(sndobsdir, 'A',&
             LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname_A)
        
        call create_GCOMW_AMSR2L3SND_filename(sndobsdir, 'D',&
             LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname_D)
        
        call read_AMSR2snd_data(n,k, fname_A, sndobs_A)
        call read_AMSR2snd_data(n,k, fname_D, sndobs_D)
        
        GCOMW_AMSR2L3SND_struc(n)%sndobs  = LIS_rc%udef
        GCOMW_AMSR2L3SND_struc(n)%sndtime = -1
        
        
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              grid_index = LIS_obs_domain(n,k)%gindex(c,r)
              if(grid_index.ne.-1) then 
                 
                 if(sndobs_A(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then        
                    GCOMW_AMSR2L3SND_struc(n)%sndobs(c,r) = &
                         sndobs_A(c+(r-1)*LIS_rc%obs_lnc(k))/100.0                 
                    lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                    lhour = 13.5 
                    call LIS_localtime2gmt (gmt,lon,lhour,zone)
                    GCOMW_AMSR2L3SND_struc(n)%sndtime(c,r) = gmt
                 endif
                 
                 if(sndobs_D(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then          
                    GCOMW_AMSR2L3SND_struc(n)%sndobs(c,r) = &
                         sndobs_D(c+(r-1)*LIS_rc%obs_lnc(k))/100.0
                    
                    lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                    lhour = 1.5
                    call LIS_localtime2gmt (gmt,lon,lhour,zone)
                    GCOMW_AMSR2L3SND_struc(n)%sndtime(c,r) = gmt
                 endif
              endif
           enddo
        enddo
        
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              grid_index = LIS_obs_domain(n,k)%gindex(c,r)
              if(grid_index.ne.-1) then 
                 if(GCOMW_AMSR2L3SND_struc(n)%usr_input_mask.eq.1) then 
                    if(GCOMW_AMSR2L3SND_struc(n)%input_mask(c,r).gt.0) then
                       GCOMW_AMSR2L3SND_struc(n)%sndobs(c,r) = LIS_rc%udef
                    endif
                 endif
              endif
           enddo
        enddo
     endif

!-------------------------------------------------------------------------
! Process IMS data if it exists
!------------------------------------------------------------------------- 
     if(GCOMW_AMSR2L3SND_struc(n)%useIMS.eq.1) then 
        call create_IMS_filename_AMSR2(imsfile,&
             GCOMW_AMSR2L3SND_struc(n)%IMSdir,&
             LIS_rc%yr,LIS_rc%doy)  
        inquire(file=imsfile,exist=file_exists)
        
        if(file_exists) then 
           
           ftn = LIS_getNextUnitNumber()
           
           write(LIS_logunit,*) '[INFO] Reading ',trim(imsfile)
           open(ftn,file=trim(imsfile),form='unformatted')
           read(ftn) IMSdata
           close(ftn)
           
           call LIS_releaseUnitNumber(ftn)
           
           ims_li  = .false.
           do c=1,IMSnc*IMSnr
              if(IMSdata(c).ge.0) then 
                 ims_li(c) = .true. 
              endif
           enddo
           
!-------------------------------------------------------------------------
! Interpolate the IMS data to the observations grid. It is assumed that
! IMS data is at a coarser resolution than the AMSR2 observation space. 
!------------------------------------------------------------------------- 

           call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
                ims_li, IMSdata, lo, IMSdata_ip,&
                GCOMW_AMSR2L3SND_struc(n)%ims_mi,&
                LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
                GCOMW_AMSR2L3SND_struc(n)%ims_rlat, &
                GCOMW_AMSR2L3SND_struc(n)%ims_rlon, &
                GCOMW_AMSR2L3SND_struc(n)%ims_w11, &
                GCOMW_AMSR2L3SND_struc(n)%ims_w12, &
                GCOMW_AMSR2L3SND_struc(n)%ims_w21, &
                GCOMW_AMSR2L3SND_struc(n)%ims_w22, &
                GCOMW_AMSR2L3SND_struc(n)%ims_n11, &
                GCOMW_AMSR2L3SND_struc(n)%ims_n12, &
                GCOMW_AMSR2L3SND_struc(n)%ims_n21, &
                GCOMW_AMSR2L3SND_struc(n)%ims_n22, &
                LIS_rc%udef, ios)
           
           do r=1,LIS_rc%obs_lnr(k)
              do c=1,LIS_rc%obs_lnc(k)
                 if(IMSdata_ip(c+LIS_rc%obs_lnc(k)*(r-1)).lt.0.0) then 
                    GCOMW_AMSR2L3SND_struc(n)%sndobs(&
                         c,r) = LIS_rc%udef
                 elseif(IMSdata_ip(c+LIS_rc%obs_lnc(k)*(r-1)).eq.0.0) then 
                    GCOMW_AMSR2L3SND_struc(n)%sndobs(&
                         c,r) = 0.0
                 endif
                 if(IMSdata_ip(c+LIS_rc%obs_lnc(k)*(r-1)).gt.0.0.and.&
                      GCOMW_AMSR2L3SND_struc(n)%sndobs(&
                      c,r).eq.0) then 
                    GCOMW_AMSR2L3SND_struc(n)%sndobs(&
                         c,r) = 0.0
                 endif
                 if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                    GCOMW_AMSR2L3SND_struc(n)%IMSdata_obs(&
                         LIS_obs_domain(n,k)%gindex(c,r)) = & 
                         IMSdata_ip(c+LIS_rc%obs_lnc(k)*(r-1))
                 endif
                 
              enddo
           enddo
           
        endif
     endif

!-------------------------------------------------------------------------
! Process MODIS data if it exists
!------------------------------------------------------------------------- 
     if(GCOMW_AMSR2L3SND_struc(n)%useMODIS.eq.1) then 
        call get_MOD10C1_filename_AMSR2(MODISfile,&
             GCOMW_AMSR2L3SND_struc(n)%MODISdir)
        
        inquire(file=MODISfile,exist=file_exists) 
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(MODISfile)
           call getMOD10data_AMSR2(n,k,MODISfile,MODISdata_ip)
           
           do r=1,LIS_rc%obs_lnr(k)
              do c=1,LIS_rc%obs_lnc(k)
                 if(MODISdata_ip(c+LIS_rc%obs_lnc(k)*(r-1)).lt.0.0) then 
                    GCOMW_AMSR2L3SND_struc(n)%sndobs(c,r) = &
                         LIS_rc%udef
                 elseif(MODISdata_ip(c+LIS_rc%obs_lnc(k)*(r-1)).eq.0.0) then 
                    GCOMW_AMSR2L3SND_struc(n)%sndobs(c,r) = &
                         0.0
                 endif
                 
                 if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                    GCOMW_AMSR2L3SND_struc(n)%MODISdata_obs(&
                         LIS_obs_domain(n,k)%gindex(c,r)) = & 
                         MODISdata_ip(c+LIS_rc%obs_lnc(k)*(r-1))
                 endif
                 
              enddo
           enddo
        endif
        
     endif

  endif
  

  call ESMF_StateGet(OBS_State,"Observation01",sndfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')

  call ESMF_FieldGet(sndfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')

  fnd = 0 
  snd_current = LIS_rc%udef

  ! dt is not defined as absolute value of the time difference to avoid
  ! double counting of the data in assimilation. 

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           if(GCOMW_AMSR2L3SND_struc(n)%sndtime(c,r).ge.0) then 
              dt = (LIS_rc%gmt - GCOMW_AMSR2L3SND_struc(n)%sndtime(c,r))*3600.0
              if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
                 snd_current(c,r) = & 
                      GCOMW_AMSR2L3SND_struc(n)%sndobs(c,r)
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
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=snd_current(c,r)
           endif
        enddo
     enddo
  endif

!-------------------------------------------------------------------------
!  Apply LSM based quality control and screening of observations
!-------------------------------------------------------------------------     
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_GCOMW_AMSR2L3SNDobsId)//char(0),n, k, OBS_state)

  snd_current = LIS_rc%udef
  call LIS_checkForValidObs(n,k,obsl,fnd,snd_current)

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
        
        call ESMF_AttributeSet(sndfield, "IMS data",&
             GCOMW_AMSR2L3SND_struc(n)%IMSdata_obs, &
             itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error in setting IMS attribute')

        call ESMF_AttributeSet(sndfield, "MODIS data",&
             GCOMW_AMSR2L3SND_struc(n)%MODISdata_obs, &
             itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error in setting IMS attribute')
     endif

  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif
end subroutine read_GCOMW_AMSR2L3SND

!BOP
! 
! !ROUTINE: read_AMSR2snd_data
! \label{read_AMSR2snd_data}
!
! !INTERFACE:
subroutine read_AMSR2snd_data(n, k, fname,sndobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod
  use map_utils,    only : latlon_to_ij
  use GCOMW_AMSR2L3SND_Mod, only : GCOMW_AMSR2L3SND_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: sndobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


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
!  \item[sndobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  real                        :: snd(1,GCOMW_AMSR2L3SND_struc(n)%amsr2nc,&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr)
  real                        :: time(GCOMW_AMSR2L3SND_struc(n)%amsr2nc,&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr)
  integer                     :: sndtime(GCOMW_AMSR2L3SND_struc(n)%amsr2nc,&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr)
  real                        :: snd_combined(GCOMW_AMSR2L3SND_struc(n)%amsr2nc,&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr)

  real                        :: snd_data(GCOMW_AMSR2L3SND_struc(n)%amsr2nc*&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr)
  logical*1                   :: snd_data_b(GCOMW_AMSR2L3SND_struc(n)%amsr2nc*&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr)
  logical*1                   :: sndobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  logical                     :: file_exists
  integer                     :: c,r,c1,i,j
  real                        :: rlat,rlon,ri,rj
  integer                     :: nid
  integer                     :: sndId
  integer                     :: timeid
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  sndtime = -1.0
  snd_combined = LIS_rc%udef

  inquire(file=fname, exist=file_exists) 
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading ',trim(fname)
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(fname))
     
     ios = nf90_inq_varid(nid, 'Geophysical Data',sndid)
     call LIS_verify(ios, 'Error nf90_inq_varid: Geophysical Data')
     
     ios = nf90_inq_varid(nid, 'Time Information',timeId)
     call LIS_verify(ios, 'Error nf90_inq_varid: Time Information')
     
     !values
     ios = nf90_get_var(nid, sndid, snd)
     call LIS_verify(ios, 'Error nf90_get_var: SND_Flags')
     
     ios = nf90_get_var(nid, timeid, time)
     call LIS_verify(ios, 'Error nf90_get_var: time')
     
     ios = nf90_close(ncid=nid)
     call LIS_verify(ios,'Error closing file '//trim(fname))
     
     do r=1, GCOMW_AMSR2L3SND_struc(n)%amsr2nr
        do c=1, GCOMW_AMSR2L3SND_struc(n)%amsr2nc
           c1 = c + GCOMW_AMSR2L3SND_struc(n)%amsr2nc/2 - 1
           if(c1.gt.GCOMW_AMSR2L3SND_struc(n)%amsr2nc) then 
              c1 = c - GCOMW_AMSR2L3SND_struc(n)%amsr2nc/2 -1
           endif
           if(time(c,r).lt.0.or.snd(1,c,r).le.0.0001) then 
              snd_combined(c1,GCOMW_AMSR2L3SND_struc(n)%amsr2nr-r+1) = &
                   LIS_rc%udef
           else
              snd_combined(c1,GCOMW_AMSR2L3SND_struc(n)%amsr2nr-r+1) = &
                   (snd(1,c,r)*0.1)*10.0 !cm to mm
              sndtime(c1,GCOMW_AMSR2L3SND_struc(n)%amsr2nr-r+1) = time(c,r)
           endif
        enddo
     enddo
  endif

  do r=1, GCOMW_AMSR2L3SND_struc(n)%amsr2nr
     do c=1, GCOMW_AMSR2L3SND_struc(n)%amsr2nc
        snd_data(c+(r-1)*GCOMW_AMSR2L3SND_struc(n)%amsr2nc) = snd_combined(c,r)
!--------------------------------------------------------------------------
! flag all values less than 10 mm
!--------------------------------------------------------------------------
!        if(snd_combined(c,r).ne.LIS_rc%udef) then 
!        if(snd_combined(c,r).ge.10) then 
        if(snd_combined(c,r).gt.50) then 
           snd_data_b(c+(r-1)*GCOMW_AMSR2L3SND_struc(n)%amsr2nc) = .true. 
        else
           snd_data_b(c+(r-1)*GCOMW_AMSR2L3SND_struc(n)%amsr2nc) = .false.
        endif
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the observation grid
!-------------------------------------------------------------------------- 
  call bilinear_interp(LIS_rc%obs_gridDesc(k,:),snd_data_b,snd_data,&
       sndobs_b_ip,sndobs_ip,&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nc*GCOMW_AMSR2L3SND_struc(n)%amsr2nr,&
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
       GCOMW_AMSR2L3SND_struc(n)%rlat,GCOMW_AMSR2L3SND_struc(n)%rlon, &
       GCOMW_AMSR2L3SND_struc(n)%w11,&
       GCOMW_AMSR2L3SND_struc(n)%w12,&
       GCOMW_AMSR2L3SND_struc(n)%w21,&
       GCOMW_AMSR2L3SND_struc(n)%w22,&
       GCOMW_AMSR2L3SND_struc(n)%n11,&
       GCOMW_AMSR2L3SND_struc(n)%n12,&
       GCOMW_AMSR2L3SND_struc(n)%n21,&
       GCOMW_AMSR2L3SND_struc(n)%n22,&
       LIS_rc%udef,ios)

#endif
  
end subroutine read_AMSR2snd_data


!BOP
! 
! !ROUTINE: read_AMSR2snd_bc_data
! \label{read_AMSR2snd_bc_data}
!
! !INTERFACE:
subroutine read_AMSR2snd_bc_data(n, k, fname,sndobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod
  use map_utils,    only : latlon_to_ij
  use GCOMW_AMSR2L3SND_Mod, only : GCOMW_AMSR2L3SND_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: sndobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


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
!  \item[sndobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  real                        :: snd(GCOMW_AMSR2L3SND_struc(n)%amsr2nc,&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr)
  real                        :: snd_combined(GCOMW_AMSR2L3SND_struc(n)%amsr2nc,&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr)

  real                        :: snd_data(GCOMW_AMSR2L3SND_struc(n)%amsr2nc*&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr)
  logical*1                   :: snd_data_b(GCOMW_AMSR2L3SND_struc(n)%amsr2nc*&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nr)
  logical*1                   :: sndobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  logical                     :: file_exists
  integer                     :: c,r,c1,i,j
  real                        :: rlat,rlon,ri,rj
  integer                     :: nid
  integer                     :: sndId
  integer                     :: timeid
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  snd_combined = LIS_rc%udef
  snd = LIS_rc%udef

  inquire(file=fname, exist=file_exists) 
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading ',trim(fname)
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(fname))
     
     ios = nf90_inq_varid(nid, 'AMSR2_snowdepth',sndid)
     call LIS_verify(ios, 'Error nf90_inq_varid: Geophysical Data')
     
     !values
     ios = nf90_get_var(nid, sndid, snd)
     call LIS_verify(ios, 'Error nf90_get_var: SND_Flags')
     
     ios = nf90_close(ncid=nid)
     call LIS_verify(ios,'Error closing file '//trim(fname))
     
  endif

  do r=1, GCOMW_AMSR2L3SND_struc(n)%amsr2nr
     do c=1, GCOMW_AMSR2L3SND_struc(n)%amsr2nc
        snd_data(c+(r-1)*GCOMW_AMSR2L3SND_struc(n)%amsr2nc) = snd(c,r)
!--------------------------------------------------------------------------
! flag all values less than 10 mm
!--------------------------------------------------------------------------
!        if(snd_combined(c,r).ne.LIS_rc%udef) then 
!        if(snd_combined(c,r).ge.10) then 
        if(snd_data(c+(r-1)*GCOMW_AMSR2L3SND_struc(n)%amsr2nc).gt.50) then 
           snd_data_b(c+(r-1)*GCOMW_AMSR2L3SND_struc(n)%amsr2nc) = .true. 
        else
           snd_data_b(c+(r-1)*GCOMW_AMSR2L3SND_struc(n)%amsr2nc) = .false.
        endif
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the observation grid
!-------------------------------------------------------------------------- 
  call bilinear_interp(LIS_rc%obs_gridDesc(k,:),snd_data_b,snd_data,&
       sndobs_b_ip,sndobs_ip,&
       GCOMW_AMSR2L3SND_struc(n)%amsr2nc*GCOMW_AMSR2L3SND_struc(n)%amsr2nr,&
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
       GCOMW_AMSR2L3SND_struc(n)%rlat,GCOMW_AMSR2L3SND_struc(n)%rlon, &
       GCOMW_AMSR2L3SND_struc(n)%w11,&
       GCOMW_AMSR2L3SND_struc(n)%w12,&
       GCOMW_AMSR2L3SND_struc(n)%w21,&
       GCOMW_AMSR2L3SND_struc(n)%w22,&
       GCOMW_AMSR2L3SND_struc(n)%n11,&
       GCOMW_AMSR2L3SND_struc(n)%n12,&
       GCOMW_AMSR2L3SND_struc(n)%n21,&
       GCOMW_AMSR2L3SND_struc(n)%n22,&
       LIS_rc%udef,ios)

#endif
  
end subroutine read_AMSR2snd_bc_data


!BOP
! !ROUTINE: create_GCOMW_AMSR2L3SND_filename
! \label{create_GCOMW_AMSR2L3SND_filename}
! 
! !INTERFACE: 
subroutine create_GCOMW_AMSR2L3SND_filename(ndir, path, yr, mo,da, filename)
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
 
  if(yr.ge.2015) then 
     if(yr.eq.2015.and.mo.le.2) then 
        if(path.eq.'A') then 
           filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
                '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
                '_01D_EQMA_L3SGSNDHA1100100.h5'
        else     
           filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
                '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
                '_01D_EQMD_L3SGSNDHA1100100.h5'
        endif
     else
        if(path.eq.'A') then 
           filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
                '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
                '_01D_EQMA_L3SGSNDHG2210210.h5'
        else     
           filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
                '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
                '_01D_EQMD_L3SGSNDHG2210210.h5'
        endif
        
     endif
  else

     if(path.eq.'A') then 
        filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
             '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
             '_01D_EQMA_L3SGSNDHA1100100.h5'
     else     
        filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
             '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
             '_01D_EQMD_L3SGSNDHA1100100.h5'
     endif
  endif

end subroutine create_GCOMW_AMSR2L3SND_filename

!BOP
! !ROUTINE: create_GCOMW_AMSR2L3SND_BC_filename
! \label{create_GCOMW_AMSR2L3SND_BC_filename}
! 
! !INTERFACE: 
subroutine create_GCOMW_AMSR2L3SND_BC_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
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
 
  filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
       '/AMSR2_BC_'//trim(fyr)//trim(fmo)//trim(fda)//'.nc4'

end subroutine create_GCOMW_AMSR2L3SND_BC_filename





!BOP
!
! !ROUTINE: create_IMS_filename_AMSR2
! \label{create_IMS_filename_AMSR2}
! 
! !INTERFACE: 
subroutine create_IMS_filename_AMSR2(name, ndir, yr, doy)
  
  implicit none
! !ARGUMENTS: 
  character(len=*) :: name, ndir
  integer           :: yr, doy
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped IMS filename
!  
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the IMS filename
!  \item[ndir] name of the IMS root directory
!  \item[yr]  current year
!  \item[doy] current day of the year
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  name = trim(ndir)//'/'//trim(fyr)//'/ims'//trim(fyr)//trim(fdoy)//'.bin'
    
end subroutine create_IMS_filename_AMSR2

!BOP
!
! !ROUTINE: getMOD10data_AMSR2
! 
! !INTERFACE: 
subroutine getMOD10data_AMSR2(n,k,name,tmp_obsl)
! !USES: 
  use LIS_coreMod,only : LIS_rc,LIS_domain
  use LIS_logMod, only : LIS_logunit, LIS_verify
  use GCOMW_AMSR2L3SND_Mod, only : GCOMW_AMSR2L3SND_struc
  
  implicit none
! 
! !DESCRIPTION: 
!  This routine retrievs the MODIS sca observations. The data is read
!  from a file in HDF-EOS format followed by the application of the 
!  confidence index and cloud cover flags. Finally the data is upscaled
!  to the AMSR2 observation grid
!EOP

#if (defined USE_HDFEOS2)
#include "hdf.f90"
#endif
  
  integer              :: n 
  integer              :: k
  character(len=*) :: name
  real                 :: tmp_obsl(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

#if (defined USE_HDFEOS2)
  integer, parameter   :: modis_nc=7200, modis_nr=3600
  integer              :: local_nc, local_nr

  !declare the hdf-eos library functions 
  integer              :: gdopen,gdattach,gddefboxreg,gdreginfo
  integer              :: gdextreg,gddetach,gdclose
  character*50         :: grid_name,ps_name,ci_name,pc_name,qa_name
  integer              :: ntype,rank,dims(2),size
  integer              :: file_id,grid_id,region_id,ret
  integer*1,allocatable    :: ps(:),ci(:),pc(:),qa(:)
  real*8               :: upleftpt(2),lowrightpt(2)
  real*8               :: cornerlon(2),cornerlat(2)

  integer              :: sfstart, sfselect, sfrdata, sfend
  integer              :: s_id1, s_id2, s_id3, sd_index, istat
  integer              :: stride(2), edges(2), start(2)
  integer*1, allocatable   :: data_map(:,:),cloud_map(:,:)
  integer*1, allocatable   :: cloud_pers_map(:,:)
  integer*1, allocatable   :: gf_snow_map(:,:)
  real,      allocatable   :: sca1(:)
  integer              :: c1,r1,c2,r2,c3,r3,c,r,gid
  
  !Grid and field names
  grid_name ="MOD_CMG_Snow_5km"
  ps_name   ="Day_CMG_Snow_Cover"
  ci_name   ="Day_CMG_Clear_Index"  ! MLW collection 6
  ! ci_name   ="Day_CMG_Confidence_Index" ! MLW collection 5
  pc_name   ="Day_CMG_Cloud_Obscured"
  qa_name   ="Snow_Spatial_QA"
  
  !open the hdf file
  
  file_id = gdopen(trim(name),DFACC_READ)
  if (file_id.eq.-1)then
     write(LIS_logunit,*)"[ERR] Failed to open hdf file",name
     return
  end if
  !  write(LIS_logunit,*) 'opened file',file_id
  
  !get the grid id
  grid_id = gdattach(file_id,grid_name)
  if (grid_id.eq.-1)then
     write(LIS_logunit,*)"[ERR] Failed to attach grid: ",grid_name,name
     ret = gdclose(file_id)
     return
  end if
  !  write(LIS_logunit,*) 'gdattach',grid_id
  
  !get the LIS domain
  cornerlat(1)=GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(4)
  cornerlon(1)=GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(5)
  cornerlat(2)=GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(7)
  cornerlon(2)=GCOMW_AMSR2L3SND_struc(n)%mod_gridDesci(8)
  
  region_id = gddefboxreg(grid_id,cornerlon,cornerlat)
  if (region_id <0)then
     write(LIS_logunit,*)"[ERR] Failed to obtain region id: gddefboxreg"
     ret = gdclose(file_id)
     return
  end if
  !  write(LIS_logunit,*) 'gddefboxreg',region_id
  
  ! find the dimensions of the subregion:dims(2)
  ret = gdreginfo(grid_id,region_id,ps_name,ntype,rank,&
       dims,size,upleftpt,lowrightpt)
  if (ret <0)then
     write(LIS_logunit,*)"[ERR] Failed to get region info: gdreginfo"
     ret = gdclose(file_id)
     return
  end if
  
  !    write(LIS_logunit,*) 'gdreginfo',ret
  
  ! get the percent snow (ps) cover field
  allocate(ps(dims(1)*dims(2)))
  ret = gdextreg(grid_id,region_id,ps_name,ps)
  if (ret <0)then
     write(LIS_logunit,*)"[ERR] Failed to get the ps field"
     ret=gddetach(grid_id)
     ret=gdclose(file_id)
     return
  end if
  
  ! get the confidence index (ci) field
  allocate(ci(dims(1)*dims(2)))
  ret = gdextreg(grid_id,region_id,ci_name,ci)
  if (ret <0)then
     write(LIS_logunit,*)"[ERR] Failed to get the CI field"
     ret=gddetach(grid_id)
     ret=gdclose(file_id)
     return
  end if
  
  ! get the percent cloud (pc) cover field
  allocate(pc(dims(1)*dims(2)))
  ret = gdextreg(grid_id,region_id,pc_name,pc)
  if (ret <0)then
     write(LIS_logunit,*)"[ERR] Failed to get the cloud field"
     ret=gddetach(grid_id)
     ret=gdclose(file_id)
     return
  end if
  
  ! get the qa flag of the grid
  allocate(qa(dims(1)*dims(2)))
  ret = gdextreg(grid_id,region_id,qa_name,qa)
  if (ret <0)then
     write(LIS_logunit,*)"[ERR] Failed to get the QA field"
     ret=gddetach(grid_id)
     ret=gdclose(file_id)
     return
  end if
  
  call calculate_mod10sca_AMSR2(n, k, dims(1), dims(2), ps,ci,pc, tmp_obsl)
  
  deallocate(ps)
  deallocate(ci)
  deallocate(pc)
  deallocate(qa)
  
  ret=gddetach(grid_id)
  if (ret <0)then
     write(LIS_logunit,*)"[ERR] Failed to detach grid_id: ",grid_id
  end if
  ret=gdclose(file_id)
  ! write(LIS_logunit,*) 'gdclose',file_id
  if (ret <0)then
     write(LIS_logunit,*)"[ERR] Failed to close file: ",file_id
  end if

#endif
end subroutine getMOD10data_AMSR2



!BOP
! 
! !ROUTINE: calculate_mod10sca_AMSR2
! \label{calculate_mod10sca_AMSR2}
!
! !INTERFACE: 
subroutine calculate_mod10sca_AMSR2(n, k, dim1, dim2, ps,ci,pc,go)
! !USES: 
  use GCOMW_AMSR2L3SND_Mod, only : GCOMW_AMSR2L3SND_struc
  use LIS_coreMod, only : LIS_rc, LIS_domain

  implicit none
! !ARGURMENTS: 
  integer,  intent(IN) :: n 
  integer,  intent(IN) :: k
  integer,  intent(IN) :: dim1
  integer,  intent(IN) :: dim2
  integer*1            :: ps(dim1*dim2)
  integer*1            :: ci(dim1*dim2)
  integer*1            :: pc(dim1*dim2)
  real                 :: go(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
! 
! !DESCRIPTION: 
!  This routine processes the raw sca observations by eliminating
!  the cloud covered points and points of low CI. This routine also
!  upscales the MODIS sca data to the AMSR2 grid.
! 
!EOP

  logical*1      :: lo(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  logical*1, allocatable :: li(:)

  integer        :: i
  real           :: sca(dim1*dim2)
  real, allocatable  :: sca1(:)

  integer        :: iret
  integer        :: c,r

  do i = 1,dim1*dim2
    if ((ps(i) >= 0.0001 .and. ps(i)<=100).and.&
         (pc(i).le.10.and.pc(i).ge.0001))then
       sca(i)=real(ps(i))*100.0/real(ci(i))
    else
      sca(i)=-9999.0
    end if
  end do
 
  
  allocate(sca1(dim1*dim2))
  allocate(li(dim1*dim2))
  
  li = .false. 
  do r=1,dim2
     do c=1,dim1
        sca1(c+(r-1)*dim1) = sca(c+(dim2-r+1-1)*dim1)
        if(sca1(c+(r-1)*dim1).ne.-9999.0) then 
           li(c+(r-1)*dim1) = .true.              
           
        endif
     enddo
  enddo
  
  call bilinear_interp(LIS_rc%obs_gridDesc(k,:),li,sca1,&
       lo,go,&
       GCOMW_AMSR2L3SND_struc(n)%mod_mi,&
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
       GCOMW_AMSR2L3SND_struc(n)%mod_rlat,&
       GCOMW_AMSR2L3SND_struc(n)%mod_rlon, &
       GCOMW_AMSR2L3SND_struc(n)%mod_w11,&
       GCOMW_AMSR2L3SND_struc(n)%mod_w12,&
       GCOMW_AMSR2L3SND_struc(n)%mod_w21,&
       GCOMW_AMSR2L3SND_struc(n)%mod_w22,&
       GCOMW_AMSR2L3SND_struc(n)%mod_n11,&
       GCOMW_AMSR2L3SND_struc(n)%mod_n12,&
       GCOMW_AMSR2L3SND_struc(n)%mod_n21,&
       GCOMW_AMSR2L3SND_struc(n)%mod_n22,&
       LIS_rc%udef,iret)

!  call upscaleByAveraging(GCOMW_AMSR2L3SND_struc(n)%mod_mi,&
!       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),LIS_rc%udef, &
!       GCOMW_AMSR2L3SND_struc(n)%mod_n11,li, sca1, lo,go)
  deallocate(sca1)
  deallocate(li)
  
end subroutine calculate_mod10sca_AMSR2

!BOP
! !ROUTINE: get_MOD10C1_filename_AMSR2
! \label{get_MOD10C1_filename_AMSR2}
!
! !REVISION HISTORY:
!  10 Sep 08    Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine get_MOD10C1_filename_AMSR2(name, ndir)
! !USES:   
  use LIS_coreMod,only : LIS_rc
  use LIS_logMod, only : LIS_logunit

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: name
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
! 
! This routine creates the name of the MODIS sca observation 
! based on the current model date. 
! 
! The arguments are:
! \begin{description}
!  \item[name] Name of the MODIS sca observation file
!  \item[ndir] Name of the MODIS sca root directory
! \end{description}
!
!EOP
 
  integer           :: yr_lp, doy

  character (len=4) :: fyr
  character (len=3) :: fdoy

  integer :: DaysOfPrev(12)
  integer :: DaysOfPrev1(12)=(/ 0,31,59,90,120,151,181,212,243,273,304,334 /)
  integer :: DaysOfPrev2(12)=(/ 0,31,60,91,121,152,182,213,244,274,305,335 /)
  
  write(unit=fyr, fmt='(i4.4)') LIS_rc%yr
  
  do yr_lp = 1904, 2096, 4
    if ( LIS_rc%yr == yr_lp ) then;  DaysOfPrev=DaysOfPrev2
    else;   DaysOfPrev=DaysOfPrev1
    end if
  end do
    
  doy=LIS_rc%da+DaysOfPrev(LIS_rc%mo)
  write(unit=fdoy, fmt='(i3.3)') doy

  name = trim(ndir)//'/'//trim(fyr)//'/'//'MOD10C1.A'&
            //trim(fyr)//trim(fdoy)//'.006.hdf' ! MLW read MOD10C1 collection 6
end subroutine get_MOD10C1_filename_AMSR2


