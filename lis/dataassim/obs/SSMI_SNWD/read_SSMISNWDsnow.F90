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
! !ROUTINE: read_SSMISNWDsnow
! \label{read_SSMISNWDsnow}
!
! !REVISION HISTORY:
!  16 Oct 2012: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_SSMISNWDsnow(n, k, OBS_State,OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use LIS_pluginIndices, only : LIS_SSMISNWDsnowobsId
  use SSMISNWDsnow_Mod, only : SSMISNWDsnow_struc
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  This routine reads and processes SSMI snow depth observations. 
!  The data is read at 0z every day and is kept in memory. At 
!  each timestep, a subset of data is chosen for use in DA if 
!  the local time of the grid point is 6AM (personal 
!  communication with George Riggs). If snow cover data from 
!  IMS or MODIS exists, they are used as additional cross-checks
!  of the data. The non-zero snow depth data from SSMI is only chosen if 
!  IMS and MODIS datasets also indicate the presence of snow. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[k]    index of the assimilation instance
!  \item[OBS\_State] observations state
!  \item[OBS\_Pert\_State] observations perturbation state
!  \end{description}
!
!EOP
  integer,         parameter    :: IMSnc = 1500,IMSnr = 375
  type(ESMF_Field)              :: snowField,pertfield
  logical                       :: alarmCheck
  logical                       :: data_upd, file_exists
  logical                       :: dataflag(LIS_npes)
  logical                       :: dataflag_local
  integer                       :: c,r, p, t
  character(len=LIS_CONST_PATH_LEN) :: obsdir, ssmi_filename, imsfile, MODISfile
  real, allocatable             :: snwd_field(:,:)
  real                          :: tsnow(SSMISNWDsnow_struc(n)%nc*SSMISNWDsnow_struc(n)%nr)
  logical*1                     :: li(SSMISNWDsnow_struc(n)%nc*SSMISNWDsnow_struc(n)%nr)
  integer                       :: ftn
  real                          :: lon, lhour
  real                          :: gmt
  real                          :: dt
  integer                       :: zone
  integer                       :: grid_index
  real                          :: ssdev(LIS_rc%obs_ngrid(k))
  logical*1                     :: lo(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real,             pointer     :: obsl(:)
  integer                       :: gid(LIS_rc%obs_ngrid(k))
  integer                       :: assimflag(LIS_rc%obs_ngrid(k))
  integer                       :: status, iret, ierr
  real                          :: IMSdata(IMSnc*IMSnr)
  logical*1                     :: IMS_li(IMSnc*IMSnr)
  real                          :: IMSdata_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                          :: MODISdata_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                       :: fnd
  real                          :: snwd_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))


  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status,'Error in AttributeGet: Data Directory')

!-------------------------------------------------------------------------
!   Read the data at 0z daily. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SSMI snow depth read alarm")

  if(alarmCheck.or.SSMISNWDsnow_struc(n)%startMode) then 

     SSMISNWDsnow_struc(n)%IMSdata_obs = -9999.0
     SSMISNWDsnow_struc(n)%MODISdata_obs = -9999.0

     SSMISNWDsnow_struc(n)%startMode = .false.
     
     SSMISNWDsnow_struc(n)%snwd = LIS_rc%udef

     call SSMIsnow_filename3(ssmi_filename,obsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)       

     inquire(file=ssmi_filename,exist=file_exists)
     if(file_exists) then 

        write(LIS_logunit,*)  '[INFO] Reading SSMI SNWD data ',trim(ssmi_filename)
        
        allocate(snwd_field(SSMISNWDsnow_struc(n)%nc, &
             SSMISNWDsnow_struc(n)%nr))
        
        ftn = LIS_getNextUnitNumber()
        open(ftn,file=ssmi_filename,form='unformatted')
        read(ftn) snwd_field
        call LIS_releaseUnitNumber(ftn)

        do r=1,SSMISNWDsnow_struc(n)%nr
           do c=1,SSMISNWDsnow_struc(n)%nc 
              tsnow(c+(r-1)*SSMISNWDsnow_struc(n)%nc) = snwd_field(c,r)
           enddo
        enddo

        li  = .false.
        do c=1,SSMISNWDsnow_struc(n)%mi
!-------------------------------------------------------------------------
! assume that 10mm is the threshold of detecting snow for passive microwave
! sensors
!-------------------------------------------------------------------------
           if(tsnow(c).ge.10) then 
              li(c) = .true. 
           endif
        enddo
        
!-------------------------------------------------------------------------
! use neighbor search approach to interpolate the SSMI data to the 
! observation grid used in assimilation. 
!-------------------------------------------------------------------------
        call neighbor_interp(LIS_rc%obs_gridDesc(k,:),li,tsnow,&
             lo,SSMISNWDsnow_struc(n)%snwd,&
             SSMISNWDsnow_struc(n)%mi,LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
             SSMISNWDsnow_struc(n)%rlat,SSMISNWDsnow_struc(n)%rlon, &
             SSMISNWDsnow_struc(n)%n11,LIS_rc%udef,iret)

        deallocate(snwd_field)

        SSMISNWDsnow_struc(n)%snwdtime = -1

!-------------------------------------------------------------------------
! Store the GMT corresponding to 6AM localtime at each grid point
!-------------------------------------------------------------------------
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                 grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
                 lon = LIS_obs_domain(n,k)%lon(grid_index)
                 
                 lhour = 6.0
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 SSMISNWDsnow_struc(n)%snwdtime(c,r) = gmt

              endif
           enddo
        enddo
        
!-------------------------------------------------------------------------
! Process IMS data if it exists
!------------------------------------------------------------------------- 
        if(SSMISNWDsnow_struc(n)%useIMS.eq.1) then 
           call create_IMS_filename_SSMI(imsfile,SSMISNWDsnow_struc(n)%IMSdir,&
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
! Upscale the IMS data to the observations grid. It is assumed that
! IMS data is at a finer resolution than the SSMI observation space. 
!------------------------------------------------------------------------- 
              call upscaleByAveraging(SSMISNWDsnow_struc(n)%ims_mi,&
                   LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),LIS_rc%udef, &
                   SSMISNWDsnow_struc(n)%ims_n11,ims_li, IMSdata, lo,IMSdata_ip)
              
              do r=1,LIS_rc%obs_lnr(k)
                 do c=1,LIS_rc%obs_lnc(k)
                    if(IMSdata_ip(c+LIS_rc%obs_lnc(k)*(r-1)).lt.0.0) then 
                       SSMISNWDsnow_struc(n)%snwd(&
                            c+LIS_rc%obs_lnc(k)*(r-1)) = LIS_rc%udef
                    elseif(IMSdata_ip(c+LIS_rc%obs_lnc(k)*(r-1)).eq.0.0) then 
                       SSMISNWDsnow_struc(n)%snwd(&
                            c+LIS_rc%obs_lnc(k)*(r-1)) = 0.0
                    endif
                    if(IMSdata_ip(c+LIS_rc%obs_lnc(k)*(r-1)).gt.0.0.and.&
                         SSMISNWDsnow_struc(n)%snwd(&
                         c+LIS_rc%obs_lnc(k)*(r-1)).eq.0) then 
                       SSMISNWDsnow_struc(n)%snwd(&
                            c+LIS_rc%obs_lnc(k)*(r-1)) = 0.0
                    endif
                    if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                       SSMISNWDsnow_struc(n)%IMSdata_obs(&
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
        if(SSMISNWDsnow_struc(n)%useMODIS.eq.1) then 
           
           call get_MOD10C1_filename_SSMI(MODISfile,&
                SSMISNWDsnow_struc(n)%MODISdir)
           
           inquire(file=MODISfile,exist=file_exists) 
           if(file_exists) then 
              write(LIS_logunit,*) '[INFO] Reading ',trim(MODISfile)
              call getMOD10data_SSMI(n,k,MODISfile,MODISdata_ip)

              do r=1,LIS_rc%obs_lnr(k)
                 do c=1,LIS_rc%obs_lnc(k)
                    if(MODISdata_ip(c+LIS_rc%obs_lnc(k)*(r-1)).lt.0.0) then 
                       SSMISNWDsnow_struc(n)%snwd(c+LIS_rc%obs_lnc(k)*(r-1)) = &
                            LIS_rc%udef
                    elseif(MODISdata_ip(c+LIS_rc%obs_lnc(k)*(r-1)).eq.0.0) then 
                       SSMISNWDsnow_struc(n)%snwd(c+LIS_rc%obs_lnc(k)*(r-1)) = &
                            0.0
                    endif
                    
                    if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                       SSMISNWDsnow_struc(n)%MODISdata_obs(&
                            LIS_obs_domain(n,k)%gindex(c,r)) = & 
                            MODISdata_ip(c+LIS_rc%obs_lnc(k)*(r-1))
                    endif
                 enddo
              enddo
           endif
           
        endif

               
     endif
  endif

  call ESMF_StateGet(OBS_State,"Observation01",snowfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(snowfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
  obsl = LIS_rc%udef 
 
!-------------------------------------------------------------------------
!  Update the OBS_State by subsetting to the local grid time  
!-------------------------------------------------------------------------     

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
           
           dt = (LIS_rc%gmt - SSMISNWDsnow_struc(n)%snwdtime(c,r))*3600.0
           lon = LIS_obs_domain(n,k)%lon(grid_index)

           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = & 
                   SSMISNWDsnow_struc(n)%snwd(grid_index)
           endif
        endif
     enddo
  enddo
 
  dataflag_local = .false. 

!-------------------------------------------------------------------------
!  Apply LSM based quality control and screening of observations
!-------------------------------------------------------------------------     
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_SSMISNWDsnowobsId)//char(0), & 
       n, k, OBS_state)

  snwd_current = LIS_rc%udef
  call LIS_checkForValidObs(n,k,obsl,fnd,snwd_current)

  if(fnd.eq.0) then 
     dataflag_local = .false. 
  else
     dataflag_local = .true. 
  endif

#if (defined SPMD)
  call MPI_ALLGATHER(dataflag_local,1, MPI_LOGICAL, dataflag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, ierr)
#endif
  data_upd = .false.
  
  do p=1,LIS_npes
     data_upd = data_upd.or.dataflag(p)
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
          .true., rc=status)
     call LIS_verify(status, 'Error: AttributeSet in Data Update Status')
     

     call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
          rc=status)
     call LIS_verify(status, 'ESMF_StateGet for Observation01 for OBS_Pert_State failed in read_SSMISNWDsnow')
     
     if(LIS_rc%obs_ngrid(k).gt.0) then 


        ssdev = SSMISNWDsnow_struc(n)%ssdev 
        do t=1,LIS_rc%obs_ngrid(k)
           if(obsl(t).ne.-9999.0) then 
              ssdev(t) =  SSMISNWDsnow_struc(n)%ssdev !+ 0.05*obsl(t)
           endif
        enddo

        call ESMF_AttributeSet(pertField,"Standard Deviation",&
             ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(snowfield,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status,'Error: AttributeSet in Grid Number')
        
        call ESMF_AttributeSet(snowfield,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')
        
        call ESMF_AttributeSet(snowfield, "IMS data",&
             SSMISNWDsnow_struc(n)%IMSdata_obs, &
             itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error in setting IMS attribute')

        call ESMF_AttributeSet(snowfield, "MODIS data",&
             SSMISNWDsnow_struc(n)%MODISdata_obs, &
             itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error in setting IMS attribute')
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status, "Error: AttributeSet Data Update Status")
     return
  end if
  
end subroutine read_SSMISNWDsnow

!BOP
!
! !ROUTINE: SSMIsnow_filename3
! \label{SSMIsnow_filename3}
! 
! !INTERFACE: 
subroutine SSMIsnow_filename3(name, ndir, yr, mo,da)
  
  implicit none
! !ARGUMENTS: 
  character(len=*)  :: name
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped SSMI filename
!  
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the NESDIS AMSRE soil moisture filename
!  \item[ndir] name of the NESDIS AMSRE soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da 

  name = trim(ndir)//'/'//trim(fyr)//'/SSMI_F08_LL_'//trim(fyr)//trim(fmo)//trim(fda)//'_SD.bin'
    
end subroutine SSMIsnow_filename3



!BOP
!
! !ROUTINE: create_IMS_filename_SSMI
! \label{create_IMS_filename_SSMI}
! 
! !INTERFACE: 
subroutine create_IMS_filename_SSMI(name, ndir, yr, doy)
  
  implicit none
! !ARGUMENTS: 
  character(len=*)  :: name
  integer           :: yr, doy
  character (len=*) :: ndir
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
    
end subroutine create_IMS_filename_SSMI

!BOP
!
! !ROUTINE: getMOD10data_SSMI
! 
! !INTERFACE: 
subroutine getMOD10data_SSMI(n,k,name,tmp_obsl)
! !USES: 
  use LIS_coreMod,only : LIS_rc,LIS_domain
  use LIS_logMod, only : LIS_logunit, LIS_verify
  use SSMISNWDsnow_Mod, only : SSMISNWDsnow_struc
  
  implicit none
! 
! !DESCRIPTION: 
!  This routine retrievs the MODIS sca observations. The data is read
!  from a file in HDF-EOS format followed by the application of the 
!  confidence index and cloud cover flags. Finally the data is upscaled
!  to the SSMI observation grid
!EOP

#if (defined USE_HDFEOS2)
#include "hdf.f90"
#endif
  
  integer              :: n 
  integer              :: k
  character(len=*)     :: name
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
  ci_name   ="Day_CMG_Confidence_Index"
  pc_name   ="Day_CMG_Cloud_Obscured"
  qa_name   ="Snow_Spatial_QA"
  
  !open the hdf file
  
  file_id = gdopen(trim(name),DFACC_READ)
  if (file_id.eq.-1)then
     write(LIS_logunit,*)"[ERR] Failed to open hdf file",trim(name)
     return
  end if
  !  write(LIS_logunit,*) 'opened file',file_id
  
  !get the grid id
  grid_id = gdattach(file_id,grid_name)
  if (grid_id.eq.-1)then
     write(LIS_logunit,*)"[ERR] Failed to attach grid: ",grid_name,trim(name)
     ret = gdclose(file_id)
     return
  end if
  !  write(LIS_logunit,*) 'gdattach',grid_id
  
  !get the LIS domain
  cornerlat(1)=SSMISNWDsnow_struc(n)%mod_gridDesci(4)
  cornerlon(1)=SSMISNWDsnow_struc(n)%mod_gridDesci(5)
  cornerlat(2)=SSMISNWDsnow_struc(n)%mod_gridDesci(7)
  cornerlon(2)=SSMISNWDsnow_struc(n)%mod_gridDesci(8)
  
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
  
  call calculate_mod10sca_SSMI(n, k, dims(1), dims(2), ps,ci,pc, tmp_obsl)
  
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
end subroutine getMOD10data_SSMI



!BOP
! 
! !ROUTINE: calculate_mod10sca_SSMI
! \label{calculate_mod10sca_SSMI}
!
! !INTERFACE: 
subroutine calculate_mod10sca_SSMI(n, k, dim1, dim2, ps,ci,pc,go)
! !USES: 
  use SSMISNWDsnow_Mod, only : SSMISNWDsnow_struc
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
!  upscales the MODIS sca data to the SSMI grid.
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
  call upscaleByAveraging(SSMISNWDsnow_struc(n)%mod_mi,&
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),LIS_rc%udef, &
       SSMISNWDsnow_struc(n)%mod_n11,li, sca1, lo,go)
  deallocate(sca1)
  deallocate(li)
  
end subroutine calculate_mod10sca_SSMI

!BOP
! !ROUTINE: get_MOD10C1_filename_SSMI
! \label{get_MOD10C1_filename_SSMI}
!
! !REVISION HISTORY:
!  10 Sep 08    Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine get_MOD10C1_filename_SSMI(name, ndir)
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
            //trim(fyr)//trim(fdoy)//'.005.hdf'
end subroutine get_MOD10C1_filename_SSMI


