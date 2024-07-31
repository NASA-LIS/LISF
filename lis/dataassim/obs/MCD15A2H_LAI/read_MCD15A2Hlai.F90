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
! !ROUTINE: read_MCD15A2Hlai
! \label{read_MCD15A2Hlai}
!
! !REVISION HISTORY:
!  16 Jun 2020    Wanshu Nie; initial specification
!
! !INTERFACE: 
subroutine read_MCD15A2Hlai(n, k, OBS_State, OBS_Pert_State)
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
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use MCD15A2Hlai_Mod, only : MCD15A2Hlai_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the MCD15A2H LAI observations from NETCDF files.

! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[k] number of observation state
!  \item[OBS\_State] observations state
!  \item[OBS\_Pert\_State] observation perturbations state
!  \end{description}
!
!EOP
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: laiobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname1,fname2, climofile1, climofile2
  integer                :: cyr, cmo, cda, chr,cmn,css,cdoy
  real                   :: wt1, wt2,ts
  integer                :: count
  real                   :: cgmt
  real*8                 :: time
  logical                :: alarmCheck, file_exists,dataCheck
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: laifield, pertField
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: laiobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                :: fnd
  real                   :: timenow

  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       laiobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "MCD15A2H LAI read alarm")

  if(MCD15A2Hlai_struc(n)%tsmooth.eq.1) then 

     if(alarmCheck.or.MCD15A2Hlai_struc(n)%startMode) then 
        MCD15A2Hlai_struc(n)%startMode = .false.
        
        cyr = LIS_rc%yr
        cmo = LIS_rc%mo
        cda = LIS_rc%da
        cdoy = LIS_rc%doy
        chr = 0 
        cmn = 0 
        css = 0 
        ts = -86400.0
        
        file_exists = .false.
        count = 0
        do while(.not.file_exists.and.count.lt.8) 
           call create_MCD15A2Hlai_filename(laiobsdir, &
                MCD15A2Hlai_struc(n)%version, cyr, cdoy, fname1,climofile1)

           inquire(file=fname1,exist=file_exists)          
           if(file_exists) then 
              call LIS_tick(MCD15A2Hlai_struc(n)%time1,cdoy,cgmt,cyr,cmo,cda, &
                   chr,cmn,css,0.0)
              exit; 
           else
              !go back a day till 8 days
              call LIS_tick(MCD15A2Hlai_struc(n)%time1,cdoy,cgmt,cyr,cmo,cda, &
                   chr,cmn,css,ts)
              count = count + 1
           endif
        enddo
        if(count.ne.8) then 
           !Find the next data location
           if(cdoy.lt.361) then
              call LIS_tick(MCD15A2Hlai_struc(n)%time2,cdoy,cgmt,cyr,cmo,cda, &
                   chr,cmn,css,(-8.0)*ts)
           else
             if(mod(cyr,4).ne.0) then
                call LIS_tick(MCD15A2Hlai_struc(n)%time2,cdoy,cgmt,cyr,cmo,cda, &
                     chr,cmn,css,(-5.0)*ts)
             else
                call LIS_tick(MCD15A2Hlai_struc(n)%time2,cdoy,cgmt,cyr,cmo,cda, &
                     chr,cmn,css,(-6.0)*ts)
             endif
           endif
           call create_MCD15A2Hlai_filename(laiobsdir, &
                MCD15A2Hlai_struc(n)%version,cyr, cdoy, fname2,climofile2)
           
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname1)
           call read_MCD15A2H_LAI_data(n,k,fname1,climofile1,&
                MCD15A2Hlai_struc(n)%laiobs1)
           
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname2)
           call read_MCD15A2H_LAI_data(n,k, fname2,climofile2,&
                MCD15A2Hlai_struc(n)%laiobs2)
           MCD15A2Hlai_struc(n)%fnd = 1
        else
           MCD15A2Hlai_struc(n)%fnd = 0 
           MCD15A2Hlai_struc(n)%laiobs1 = LIS_rc%udef
           MCD15A2Hlai_struc(n)%laiobs2 = LIS_rc%udef
        endif
     endif
  else
     if(alarmCheck.or.MCD15A2Hlai_struc(n)%startMode) then 
        MCD15A2Hlai_struc(n)%startMode = .false.
        
        call create_MCD15A2Hlai_filename(laiobsdir, &
             MCD15A2Hlai_struc(n)%version, LIS_rc%yr, LIS_rc%doy, fname1, climofile1)
        
        inquire(file=fname1,exist=file_exists)          
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname1)
           call read_MCD15A2H_LAI_data(n,k, fname1,climofile1,laiobs)
           fnd = 1
        else
           fnd = 0 
           write(LIS_logunit,*) '[WARN] Missing LAI file: ',trim(fname1)
        endif
     else
        fnd = 0 
        laiobs = LIS_rc%udef
     endif
  endif
     
  if(MCD15A2Hlai_struc(n)%tsmooth.eq.1) then 
     !interpolate between two data points every day
     timenow = float(LIS_rc%hr)*3600 + 60.0*LIS_rc%mn + LIS_rc%ss
     alarmCheck = (mod(timenow, 86400.0).eq.0)
     
     if(alarmCheck) then 
        call LIS_tick(time,cdoy,cgmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
             LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,0.0)

        wt2 = (time - MCD15A2Hlai_struc(n)%time1)/&
             (MCD15A2Hlai_struc(n)%time2-MCD15A2Hlai_struc(n)%time1)
        wt1 = 1.0 - wt2
        
        if(MCD15A2Hlai_struc(n)%fnd.eq.1) then 
           do t=1,LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)
              if(MCD15A2Hlai_struc(n)%laiobs1(t).ne.-9999.0.and.&
                   MCD15A2Hlai_struc(n)%laiobs2(t).ne.-9999.0) then 
                 laiobs(t) = (MCD15A2Hlai_struc(n)%laiobs1(t)*wt1 + & 
                      MCD15A2Hlai_struc(n)%laiobs2(t)*wt2)
              else
                 laiobs(t) = LIS_rc%udef
              endif
           enddo
           fnd = 1
        else
           laiobs = LIS_rc%udef
           fnd = 0 
        endif
     endif

  endif

  dataCheck = .false.
  if(alarmCheck) then 
     if(MCD15A2Hlai_struc(n)%tsmooth.eq.1) then 
        if(MCD15A2Hlai_struc(n)%fnd.ne.0) then 
           dataCheck = .true. 
           fnd = 1
        endif
     else
        if(fnd.eq.1) then 
           dataCheck = .true. 
        endif
     endif
  else
     fnd = 0 
     dataCheck = .false.
  endif

  if(dataCheck) then 
        
     call ESMF_StateGet(OBS_State,"Observation01",laifield,&
          rc=status)
     call LIS_verify(status, 'Error: StateGet Observation01')
     
     call ESMF_FieldGet(laifield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status, 'Error: FieldGet')
     
     obsl = LIS_rc%udef 
     do r=1, LIS_rc%obs_lnr(k)
        do c=1, LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                   laiobs(c+(r-1)*LIS_rc%obs_lnc(k))
           endif
        enddo
     enddo
     
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
           call ESMF_AttributeSet(laifield,"Grid Number",&
                gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(laifield,"Assimilation Flag",&
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
end subroutine read_MCD15A2Hlai

!BOP
! 
! !ROUTINE: read_MCD15A2H_LAI_data
! \label{read_MCD15A2H_LAI_data}
!
! !INTERFACE:
subroutine read_MCD15A2H_LAI_data(n, k, fname, climofile, laiobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use MCD15A2Hlai_Mod, only : MCD15A2Hlai_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  character (len=*)             :: climofile
  real                          :: laiobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real*8                        :: cornerlat(2), cornerlon(2)
  character*3                   :: fdoy

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the MCD15A2H LAI file and applies the data
!  quality flags to filter the data. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[k]            number of observation state
!  \item[k]            number of observation state
!  \item[fname]        name of the MCD15A2H LAI file
!  \item[climofile]    Generated MCD152AH LAI climatology file
!  \item[laiobs\_ip]   MCD15A2H LAI data processed to the LIS domain
!  \end{description}
!
!
!EOP

!--------------Wanshu -----------------------
  integer,  parameter     :: nc=86400, nr=43200
  integer                 :: lat_off, lon_off
  integer                 :: lai(MCD15A2Hlai_struc(n)%nc,MCD15A2Hlai_struc(n)%nr)
  integer                 :: flag(MCD15A2Hlai_struc(n)%nc,MCD15A2Hlai_struc(n)%nr)
  real                    :: lai_flagged(MCD15A2Hlai_struc(n)%nc,MCD15A2Hlai_struc(n)%nr)
  real                    :: lai_in(MCD15A2Hlai_struc(n)%nc*MCD15A2Hlai_struc(n)%nr)
  logical*1               :: lai_data_b(MCD15A2Hlai_struc(n)%nc*MCD15A2Hlai_struc(n)%nr)
  logical*1               :: laiobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                    :: laiobs_climo_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                 :: c,r,t
  integer                 :: nid
  integer                 :: laiid, flagid
  integer                 :: ios
  
  integer, dimension(nf90_max_var_dims) :: dimIDs
  integer                                :: numLons, numLats
  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LIS_verify(ios,'Error opening file '//trim(fname))
  
  ios = nf90_inq_varid(nid, 'Lai_500m',laiid)
  call LIS_verify(ios, 'Error nf90_inq_varid: Lai_500m')
  
  ios = nf90_inq_varid(nid, 'FparLai_QC',flagid)
  call LIS_verify(ios, 'Error nf90_inq_varid: flag')

  !values
  
  cornerlat(1)=MCD15A2Hlai_struc(n)%gridDesci(4)
  cornerlon(1)=MCD15A2Hlai_struc(n)%gridDesci(5)
  cornerlat(2)=MCD15A2Hlai_struc(n)%gridDesci(7)
  cornerlon(2)=MCD15A2Hlai_struc(n)%gridDesci(8)
  
  lai_data_b = .false.
  
  lat_off = nint((cornerlat(1)+89.9979167)/0.00416667)+1
  lon_off = nint((cornerlon(1)+179.9979167)/0.00416667)+1


  ios = nf90_get_var(nid, laiid, lai, &
       start=(/lon_off,lat_off/), &
       count=(/MCD15A2Hlai_struc(n)%nc,MCD15A2Hlai_struc(n)%nr/)) 
  
  call LIS_verify(ios, 'Error nf90_get_var: Lai_500m')
  
  ios = nf90_get_var(nid, flagid, flag, &
       start=(/lon_off,lat_off/), &
       count=(/MCD15A2Hlai_struc(n)%nc,MCD15A2Hlai_struc(n)%nr/))
  
  call LIS_verify(ios, 'Error nf90_get_var: flag')
  
  ios = nf90_close(ncid=nid)
  call LIS_verify(ios,'Error closing file '//trim(fname))
  
  do r=1, MCD15A2Hlai_struc(n)%nr
     do c=1, MCD15A2Hlai_struc(n)%nc

        if(MCD15A2Hlai_struc(n)%qcflag.eq.1) then !apply QC flag

          if(lai(c,r).gt.0.and.lai(c,r).le.100) then
             if (MOD(flag(c,r),2) ==0.and.flag(c,r).le.62) then
                lai_flagged(c,r) =&
                   lai(c,r)*0.1
             else
               lai_flagged(c,r) = LIS_rc%udef
             endif
          else
            lai_flagged(c,r) = LIS_rc%udef
          endif

        else  ! no QC flag applied                

           if(lai(c,r).gt.0.and.lai(c,r).le.100) then
              lai_flagged(c,r) =&
                   lai(c,r)*0.1
           else
              lai_flagged(c,r) = LIS_rc%udef
           endif
        endif
     end do
  end do


  do r=1, MCD15A2Hlai_struc(n)%nr
     do c=1, MCD15A2Hlai_struc(n)%nc
        lai_in(c+(r-1)*MCD15A2Hlai_struc(n)%nc) = lai_flagged(c,r)
        if(lai_flagged(c,r).ne.LIS_rc%udef) then
           lai_data_b(c+(r-1)*MCD15A2Hlai_struc(n)%nc) = .true.
        else
           lai_data_b(c+(r-1)*MCD15A2Hlai_struc(n)%nc) = .false.
        endif
     enddo
  enddo

  if(LIS_rc%obs_gridDesc(k,10).le.0.00416667) then 
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
          lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
          MCD15A2Hlai_struc(n)%nc*MCD15A2Hlai_struc(n)%nr, &
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          MCD15A2Hlai_struc(n)%rlat,MCD15A2Hlai_struc(n)%rlon,&
          MCD15A2Hlai_struc(n)%w11,MCD15A2Hlai_struc(n)%w12,&
          MCD15A2Hlai_struc(n)%w21,MCD15A2Hlai_struc(n)%w22,&
          MCD15A2Hlai_struc(n)%n11,MCD15A2Hlai_struc(n)%n12,&
          MCD15A2Hlai_struc(n)%n21,MCD15A2Hlai_struc(n)%n22,LIS_rc%udef,ios)
  else
     call upscaleByAveraging(MCD15A2Hlai_struc(n)%nc*MCD15A2Hlai_struc(n)%nr,&
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          LIS_rc%udef, MCD15A2Hlai_struc(n)%n11,&
          lai_data_b,lai_in, laiobs_b_ip, laiobs_ip)

  endif
  
  if(MCD15A2Hlai_struc(n)%climofill.eq.1) then 
  
     write(LIS_logunit,*) '[INFO] Opening climo file ',trim(climofile)
     ios = nf90_open(path=trim(climofile),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(climofile))
     
     ios = nf90_inq_varid(nid, 'Lai_500m',laiid)
     call LIS_verify(ios, 'Error nf90_inq_varid: Lai_500m')
     
     cornerlat(1)=MCD15A2Hlai_struc(n)%gridDesci(4)
     cornerlon(1)=MCD15A2Hlai_struc(n)%gridDesci(5)
     cornerlat(2)=MCD15A2Hlai_struc(n)%gridDesci(7)
     cornerlon(2)=MCD15A2Hlai_struc(n)%gridDesci(8)
     
     lai_data_b = .false.
     
     lat_off = nint((cornerlat(1)+89.9979167)/0.00416667)+1
     lon_off = nint((cornerlon(1)+179.9979167)/0.00416667)+1

     ios = nf90_get_var(nid, laiid, lai, &
          start=(/lon_off,lat_off/), &
          count=(/MCD15A2Hlai_struc(n)%nc,MCD15A2Hlai_struc(n)%nr/)) 
     
     call LIS_verify(ios, 'Error nf90_get_var: Lai_500m')
     
     ios = nf90_close(ncid=nid)
     call LIS_verify(ios,'Error closing file '//trim(fname))
     
     do r=1, MCD15A2Hlai_struc(n)%nr
        do c=1, MCD15A2Hlai_struc(n)%nc
           
           if(lai(c,r).gt.0.and.lai(c,r).le.100) then
              lai_flagged(c,r) = lai(c,r)*0.1
           else
              lai_flagged(c,r) = LIS_rc%udef
           endif
           
        end do
     end do

     do r=1, MCD15A2Hlai_struc(n)%nr
        do c=1, MCD15A2Hlai_struc(n)%nc
           lai_in(c+(r-1)*MCD15A2Hlai_struc(n)%nc) = lai_flagged(c,r)
           if(lai_flagged(c,r).ne.LIS_rc%udef) then
              lai_data_b(c+(r-1)*MCD15A2Hlai_struc(n)%nc) = .true.
           else
              lai_data_b(c+(r-1)*MCD15A2Hlai_struc(n)%nc) = .false.
           endif
        enddo
     enddo

     if(LIS_rc%obs_gridDesc(k,10).le.0.00416667) then 
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
        call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
             lai_data_b, lai_in, laiobs_b_ip, laiobs_climo_ip, &
             MCD15A2Hlai_struc(n)%nc*MCD15A2Hlai_struc(n)%nr, &
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             MCD15A2Hlai_struc(n)%rlat,MCD15A2Hlai_struc(n)%rlon,&
             MCD15A2Hlai_struc(n)%w11,MCD15A2Hlai_struc(n)%w12,&
             MCD15A2Hlai_struc(n)%w21,MCD15A2Hlai_struc(n)%w22,&
             MCD15A2Hlai_struc(n)%n11,MCD15A2Hlai_struc(n)%n12,&
             MCD15A2Hlai_struc(n)%n21,MCD15A2Hlai_struc(n)%n22,LIS_rc%udef,ios)
     else
        call upscaleByAveraging(&
             MCD15A2Hlai_struc(n)%nc*MCD15A2Hlai_struc(n)%nr,&
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             LIS_rc%udef, MCD15A2Hlai_struc(n)%n11,&
             lai_data_b,lai_in, laiobs_b_ip, laiobs_climo_ip)
        
     endif
     
     do t=1,LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)
        if(laiobs_ip(t).eq.-9999.0.and.laiobs_climo_ip(t).ne.-9999.0) then 
           laiobs_ip(t) = laiobs_climo_ip(t)
        endif
     enddo
  endif
#endif

end subroutine read_MCD15A2H_LAI_data


!BOP
! !ROUTINE: create_MCD15A2Hlai_filename
! \label{create_MCD15A2Hlai_filename}
! 
! !INTERFACE: 
subroutine create_MCD15A2Hlai_filename(ndir, version, yr, doy, filename, climofile)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: climofile
  character(len=*)  :: version
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the MCD15A2H LAI filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the MCD15A2H LAI data directory
!  \item[version] version of the MCD15A2H LAI data
!  \item[yr]  current year
!  \item[doy]  current day of the year
!  \item[filename] Generated MCD15A2H LAI filename
!  \item[climofile] Generated MCD152AH LAI climatology file
!  \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  if(version.eq."006") then
     filename = trim(ndir)//'/'//trim(fyr)//'/MCD15A2H.006_LAI_'//&
          trim(fyr)//trim(fdoy)//'.nc4'
  endif

  climofile = trim(ndir)//'/MCD15A2H.006_LAI_YYYY'//&
       trim(fdoy)//'.nc4'

end subroutine create_MCD15A2Hlai_filename





