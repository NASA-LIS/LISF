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
! !ROUTINE: read_GLASSlai
! \label{read_GLASSlai}
!
! !REVISION HISTORY:
!  21 Dec 2017    Sujay Kumar; initial specification
!
! !INTERFACE: 
subroutine read_GLASSlai(n, k, OBS_State, OBS_Pert_State)
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
  use GLASSlai_Mod, only : GLASSlai_struc

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
!  !NOTES: The subroutine only supports V4 of the data. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: laiobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname1,fname2
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

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "GLASS LAI read alarm")

  if(GLASSlai_struc(n)%tsmooth.eq.1) then 

     if(alarmCheck.or.GLASSlai_struc(n)%startMode) then 
        GLASSlai_struc(n)%startMode = .false.
        
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
           call create_GLASSlai_filename(laiobsdir, &
                GLASSlai_struc(n)%source, cyr, cdoy, fname1)
           
           inquire(file=fname1,exist=file_exists)          
           if(file_exists) then 
              call LIS_tick(GLASSlai_struc(n)%time1,cdoy,cgmt,cyr,cmo,cda, &
                   chr,cmn,css,0.0)
              exit; 
           else
              !go back a day till 8 days
              call LIS_tick(GLASSlai_struc(n)%time1,cdoy,cgmt,cyr,cmo,cda, &
                   chr,cmn,css,ts)
              count = count + 1
           endif
        enddo
        if(count.ne.8) then 
           !Find the next data location
           call LIS_tick(GLASSlai_struc(n)%time2,cdoy,cgmt,cyr,cmo,cda, &
                chr,cmn,css,(-8.0)*ts)
           call create_GLASSlai_filename(laiobsdir, &
                GLASSlai_struc(n)%source,cyr, cdoy, fname2)
           
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname1)
           call read_GLASS_LAI_data(n,k, fname1,GLASSlai_struc(n)%laiobs1)
           
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname2)
           call read_GLASS_LAI_data(n,k, fname2,GLASSlai_struc(n)%laiobs2)
           GLASSlai_struc(n)%fnd = 1
        else
           GLASSlai_struc(n)%fnd = 0 
           GLASSlai_struc(n)%laiobs1 = LIS_rc%udef
           GLASSlai_struc(n)%laiobs2 = LIS_rc%udef
        endif
     endif
  else
     if(alarmCheck.or.GLASSlai_struc(n)%startMode) then 
        GLASSlai_struc(n)%startMode = .false.
        
        call create_GLASSlai_filename(laiobsdir, &
             GLASSlai_struc(n)%source, LIS_rc%yr, LIS_rc%doy, fname1)
        
        inquire(file=fname1,exist=file_exists)          
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname1)
           call read_GLASS_LAI_data(n,k, fname1,laiobs)
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
     
  if(GLASSlai_struc(n)%tsmooth.eq.1) then 
     !interpolate between two data points every day
     timenow = float(LIS_rc%hr)*3600 + 60.0*LIS_rc%mn + LIS_rc%ss
     alarmCheck = (mod(timenow, 86400.0).eq.0)
     
     if(alarmCheck) then 
        call LIS_tick(time,cdoy,cgmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
             LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,0.0)
        wt2 = (time - GLASSlai_struc(n)%time1)/&
             (GLASSlai_struc(n)%time2-GLASSlai_struc(n)%time1)
        wt1 = 1.0 - wt2
        
        if(GLASSlai_struc(n)%fnd.eq.1) then 
           do t=1,LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)
              if(GLASSlai_struc(n)%laiobs1(t).ne.-9999.0.and.&
                   GLASSlai_struc(n)%laiobs2(t).ne.-9999.0) then 
                 laiobs(t) = (GLASSlai_struc(n)%laiobs1(t)*wt1 + & 
                      GLASSlai_struc(n)%laiobs2(t)*wt2)
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
     if(GLASSlai_struc(n)%tsmooth.eq.1) then 
        if(GLASSlai_struc(n)%fnd.ne.0) then 
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
end subroutine read_GLASSlai

!BOP
! 
! !ROUTINE: read_GLASS_LAI_data
! \label{read_GLASS_LAI_data}
!
! !INTERFACE:
subroutine read_GLASS_LAI_data(n, k, fname, laiobs_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use GLASSlai_Mod, only : GLASSlai_struc

  implicit none
#if (defined USE_HDFEOS2)
#include "hdf.f90"
#endif
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: laiobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real*8                        :: cornerlat(2), cornerlon(2)


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the GLASS LAI file and applies the data
!  quality flags to filter the data. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTGLASS AMSR-E file
!  \item[laiobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDFEOS2)
  integer,  parameter     :: nc=7200, nr=3600
  integer                 :: gdopen,gdattach,gdrdfld
  integer                 :: gddetach,gdclose
  integer                 :: file_id,grid_id,region_id,iret,c,r,c1,r1
  character*50            :: grid_name,lai_name
  integer                 :: start(2), edge(2), stride(2)
  integer*2, allocatable  :: lai_raw_avhrr(:)
  integer*1, allocatable  :: lai_raw_modis(:)
  integer                 :: lat_off, lon_off
  real                    :: lai_in(GLASSlai_struc(n)%nc*GLASSlai_struc(n)%nr)
  logical*1               :: lai_data_b(GLASSlai_struc(n)%nc*GLASSlai_struc(n)%nr)
  logical*1               :: laiobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  real                    :: test(nc,nr)
  if(GLASSlai_struc(n)%source.eq."AVHRR") then 
     grid_name ="GLASS01B02"
  elseif(GLASSlai_struc(n)%source.eq."MODIS") then 
     grid_name ="GLASS01B01"
  endif

  file_id = gdopen(trim(fname),DFACC_READ)
  if (file_id.eq.-1)then
     write(LIS_logunit,*) "[ERR] Failed to open hdf file",trim(fname)
  end if
  
  lai_name = "LAI"

  grid_id = gdattach(file_id,grid_name)

  start(1)=0  !hdfeos lib uses 0-based count
  start(2)=0
  edge(1)=nc
  edge(2)=nr
  stride(1)=1
  stride(2)=1
  
  cornerlat(1)=GLASSlai_struc(n)%gridDesci(4)
  cornerlon(1)=GLASSlai_struc(n)%gridDesci(5)
  cornerlat(2)=GLASSlai_struc(n)%gridDesci(7)
  cornerlon(2)=GLASSlai_struc(n)%gridDesci(8)

  if(GLASSlai_struc(n)%source.eq."AVHRR") then 
     allocate(lai_raw_avhrr(nc*nr))
     iret = gdrdfld(grid_id,lai_name,start,stride,edge,lai_raw_avhrr)
     
     lai_data_b = .false. 
     
     lat_off = nint((cornerlat(1)+89.975)/0.05)+1
     lon_off = nint((cornerlon(1)+179.975)/0.05)+1
     
     do r=1,GLASSlai_struc(n)%nr
        do c=1,GLASSlai_struc(n)%nc
           c1 = c + lon_off
           r1 = nr - (r + lat_off) + 1
           
           if(lai_raw_avhrr(c1+(r1-1)*nc).gt.0.and.&
                lai_raw_avhrr(c1+(r1-1)*nc).ne.2550) then 
              lai_in(c+(r-1)*GLASSlai_struc(n)%nc) =&
                   lai_raw_avhrr(c1+(r1-1)*nc)*0.01
              lai_data_b(c+(r-1)*GLASSlai_struc(n)%nc) =  .true. 
           else
              lai_in(c+(r-1)*GLASSlai_struc(n)%nc) = -9999.0
              lai_data_b(c+(r-1)*GLASSlai_struc(n)%nc) = .false. 
           endif
        enddo
     enddo
     deallocate(lai_raw_avhrr)
  elseif(GLASSlai_struc(n)%source.eq."MODIS") then 
     allocate(lai_raw_modis(nc*nr))
     iret = gdrdfld(grid_id,lai_name,start,stride,edge,lai_raw_modis)
     lai_data_b = .false. 
     
     lat_off = nint((cornerlat(1)+89.975)/0.05)+1
     lon_off = nint((cornerlon(1)+179.975)/0.05)+1
     
     do r=1,GLASSlai_struc(n)%nr
        do c=1,GLASSlai_struc(n)%nc
           c1 = c + lon_off
           r1 = nr - (r + lat_off) + 1
           
           if(lai_raw_modis(c1+(r1-1)*nc).gt.0.and.&
                lai_raw_modis(c1+(r1-1)*nc).ne.2550) then 
              lai_in(c+(r-1)*GLASSlai_struc(n)%nc) =&
                   lai_raw_modis(c1+(r1-1)*nc)*0.1
              lai_data_b(c+(r-1)*GLASSlai_struc(n)%nc) =  .true. 
           else
              lai_in(c+(r-1)*GLASSlai_struc(n)%nc) = -9999.0
              lai_data_b(c+(r-1)*GLASSlai_struc(n)%nc) = .false. 
           endif
        enddo
     enddo
     deallocate(lai_raw_modis)

  endif


  iret=gddetach(grid_id)
  iret=gdclose(file_id)

  if(LIS_rc%obs_gridDesc(k,10).le.0.05) then 
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
          lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
          GLASSlai_struc(n)%nc*GLASSlai_struc(n)%nr, &
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          GLASSlai_struc(n)%rlat,GLASSlai_struc(n)%rlon,&
          GLASSlai_struc(n)%w11,GLASSlai_struc(n)%w12,&
          GLASSlai_struc(n)%w21,GLASSlai_struc(n)%w22,&
          GLASSlai_struc(n)%n11,GLASSlai_struc(n)%n12,&
          GLASSlai_struc(n)%n21,GLASSlai_struc(n)%n22,LIS_rc%udef,iret)
  else
     call upscaleByAveraging(GLASSlai_struc(n)%nc*GLASSlai_struc(n)%nr,&
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          LIS_rc%udef, GLASSlai_struc(n)%n11,&
          lai_data_b,lai_in, laiobs_b_ip, laiobs_ip)
  endif
  
#endif

end subroutine read_GLASS_LAI_data




!BOP
! !ROUTINE: create_GLASSlai_filename
! \label{create_GLASSlai_filename}
! 
! !INTERFACE: 
subroutine create_GLASSlai_filename(ndir, source, yr, doy, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: source
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the GLASS LAI filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the GLASS LAI data directory
!  \item[yr]  current year
!  \item[mo]  current doy
!  \item[filename] Generated GLASS LAI filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  if(source.eq."AVHRR") then 
     filename = trim(ndir)//'/'//trim(fyr)//'/GLASS01B02.V04.A'//&
          trim(fyr)//trim(fdoy)//'.hdf'
  elseif(source.eq."MODIS") then 
     filename = trim(ndir)//'/'//trim(fyr)//'/GLASS01B01.V04.A'//&
          trim(fyr)//trim(fdoy)//'.hdf'

  endif
end subroutine create_GLASSlai_filename





