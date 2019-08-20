!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_GLASSalbedo
! \label{read_GLASSalbedo}
!
! !REVISION HISTORY:
!  21 Dec 2017    Sujay Kumar; initial specification
!
! !INTERFACE: 
subroutine read_GLASSalbedo(n, k, OBS_State, OBS_Pert_State)
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
  use GLASSalbedo_Mod, only : GLASSalbedo_struc

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
  character*100          :: albedoobsdir
  character*100          :: fname1,fname2
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: albedofield, pertField
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: wt1, wt2,ts
  integer                :: count
  integer                :: cyr, cmo, cda, chr,cmn,css,cdoy
  real                   :: cgmt
  real*8                 :: time
  real                   :: albedoobs_bs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: albedoobs_ws(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                :: fnd
  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       albedoobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "GLASS Albedo read alarm")

  if(alarmCheck.or.GLASSalbedo_struc(n)%startMode) then 
     GLASSalbedo_struc(n)%startMode = .false.

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
        call create_GLASSalbedo_filename(albedoobsdir, &
             GLASSalbedo_struc(n)%source, &
             cyr, cdoy, fname1)
     
        inquire(file=fname1,exist=file_exists)          
        if(file_exists) then 
           exit;
        else
           !go back a day till 8 days
           call LIS_tick(GLASSalbedo_struc(n)%time1,cdoy,cgmt,cyr,cmo,cda, &
                chr,cmn,css,ts)
           count = count + 1
        endif
     enddo
     if(count.ne.8) then 
        !Find the next data location
        
        call LIS_tick(GLASSalbedo_struc(n)%time2,cdoy,cgmt,cyr,cmo,cda, &
             chr,cmn,css,(-8.0)*ts)
        call create_GLASSalbedo_filename(albedoobsdir, &
             GLASSalbedo_struc(n)%source, &
             cyr, cdoy, fname2)

        write(LIS_logunit,*) '[INFO] Reading ',trim(fname1)
        call read_GLASS_ALBEDO_data(n,k, fname1,&
             GLASSalbedo_struc(n)%source, &
             GLASSalbedo_struc(n)%obs_bs1,&
             GLASSalbedo_struc(n)%obs_ws1)
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname2)
        call read_GLASS_ALBEDO_data(n,k, fname2,&
             GLASSalbedo_struc(n)%source, &
             GLASSalbedo_struc(n)%obs_bs2,&
             GLASSalbedo_struc(n)%obs_ws2)
        GLASSalbedo_struc(n)%fnd = 1
     else
        GLASSalbedo_struc(n)%fnd = 0 
        GLASSalbedo_struc(n)%obs_bs1 = LIS_rc%udef
        GLASSalbedo_struc(n)%obs_ws1 = LIS_rc%udef
        GLASSalbedo_struc(n)%obs_bs2 = LIS_rc%udef
        GLASSalbedo_struc(n)%obs_ws2 = LIS_rc%udef
     endif
  endif
  !interpolate in time between the two datapoints. 
  call LIS_tick(time,cdoy,cgmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
       LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,0.0)
  wt1 = (time - GLASSalbedo_struc(n)%time1)/&
       (GLASSalbedo_struc(n)%time2-GLASSalbedo_struc(n)%time1)
  wt2 = 1.0 - wt1
  
  if(GLASSalbedo_struc(n)%fnd.eq.1) then 
     albedoobs_bs = (GLASSalbedo_struc(n)%obs_bs1*wt1 + & 
          GLASSalbedo_struc(n)%obs_bs2*wt2)
     albedoobs_ws = (GLASSalbedo_struc(n)%obs_ws1*wt1 + & 
          GLASSalbedo_struc(n)%obs_ws2*wt2)
     fnd = 1
  else
     albedoobs_bs = LIS_rc%udef
     albedoobs_ws = LIS_rc%udef
     fnd = 0 
  endif
  
  if(GLASSalbedo_struc(n)%fnd.ne.0) then 
     call ESMF_StateGet(OBS_State,"Observation01",albedofield,&
          rc=status)
     call LIS_verify(status, 'Error: StateGet Observation01')
     
     call ESMF_FieldGet(albedofield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status, 'Error: FieldGet')
     
     obsl = LIS_rc%udef 
     do r=1, LIS_rc%obs_lnr(k)
        do c=1, LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                   albedoobs_bs(c+(r-1)*LIS_rc%obs_lnc(k))
           endif
        enddo
     enddo

     call ESMF_StateGet(OBS_State,"Observation02",albedofield,&
          rc=status)
     call LIS_verify(status, 'Error: StateGet Observation02')
     
     call ESMF_FieldGet(albedofield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status, 'Error: FieldGet')
     
     obsl = LIS_rc%udef 
     do r=1, LIS_rc%obs_lnr(k)
        do c=1, LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                   albedoobs_ws(c+(r-1)*LIS_rc%obs_lnc(k))
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
           call ESMF_AttributeSet(albedofield,"Grid Number",&
                gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(albedofield,"Assimilation Flag",&
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
end subroutine read_GLASSalbedo

!BOP
! 
! !ROUTINE: read_GLASS_ALBEDO_data
! \label{read_GLASS_ALBEDO_data}
!
! !INTERFACE:
subroutine read_GLASS_ALBEDO_data(n, k, fname, source, albedoobs_bs_ip, &
     albedoobs_ws_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use GLASSalbedo_Mod, only : GLASSalbedo_struc

  implicit none
#if (defined USE_HDF4) 
#include "hdf.f90"
#endif
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  character (len=*)             :: source
  real                          :: albedoobs_bs_ip(LIS_rc%obs_lnc(k)*&
       LIS_rc%obs_lnr(k))
  real                          :: albedoobs_ws_ip(LIS_rc%obs_lnc(k)*&
       LIS_rc%obs_lnr(k))
  real*8                        :: cornerlat(2), cornerlon(2)


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the GLASS ALBEDO file and applies the data
!  quality flags to filter the data. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTGLASS AMSR-E file
!  \item[albedoobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDF4)
  integer,  parameter     :: nc=7200, nr=3600
  integer                 :: gdopen,gdattach,gdrdfld
  integer                 :: gddetach,gdclose
  integer                 :: file_id,grid_id,region_id,iret,c,r,c1,r1
  character*50            :: grid_name,albedo_bs_name, albedo_ws_name
  integer                 :: start(2), edge(2), stride(2)
  integer*2, allocatable  :: albedo_raw_avhrr(:)
  integer*2, allocatable  :: albedo_raw_modis(:)
  integer*2, allocatable  :: albedo_qc(:)
  integer                 :: lat_off, lon_off
  real                    :: albedo_bs_in(GLASSalbedo_struc(n)%nc*&
       GLASSalbedo_struc(n)%nr)
  real                    :: albedo_ws_in(GLASSalbedo_struc(n)%nc*&
       GLASSalbedo_struc(n)%nr)
  logical*1               :: albedo_data_b(GLASSalbedo_struc(n)%nc*&
       GLASSalbedo_struc(n)%nr)
  logical*1               :: albedoobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  if(source.eq."AVHRR") then 
     grid_name ="GLASS02B05"
  elseif(source.eq."MODIS") then 
     grid_name ="GLASS02B06"
  endif

  file_id = gdopen(trim(fname),DFACC_READ)
  if (file_id.eq.-1)then
     write(LIS_logunit,*) "[ERR] Failed to open hdf file",fname
  end if
  
  albedo_bs_name = "ABD_BSA_VIS"
  albedo_ws_name = "ABD_WSA_VIS"

  grid_id = gdattach(file_id,grid_name)

  start(1)=0  !hdfeos lib uses 0-based count
  start(2)=0
  edge(1)=nc
  edge(2)=nr
  stride(1)=1
  stride(2)=1
  
  cornerlat(1)=GLASSalbedo_struc(n)%gridDesci(4)
  cornerlon(1)=GLASSalbedo_struc(n)%gridDesci(5)
  cornerlat(2)=GLASSalbedo_struc(n)%gridDesci(7)
  cornerlon(2)=GLASSalbedo_struc(n)%gridDesci(8)

  if(GLASSalbedo_struc(n)%qcflag.eq.1) then 

     allocate(albedo_qc(nc*nr))
     
     iret = gdrdfld(grid_id,"QC_VIS",&
          start, stride, edge, albedo_qc)           
! The 0-1  bit in the QC flag denotes the overall quality. 
! 0 indicates a good value and 1 indicates acceptable value
  endif

  if(GLASSalbedo_struc(n)%source.eq."AVHRR") then 
!black sky albedo /direct
     allocate(albedo_raw_avhrr(nc*nr))
     iret = gdrdfld(grid_id,albedo_bs_name,start,stride,edge,albedo_raw_avhrr)
     
     albedo_data_b = .false. 
     
     lat_off = nint((cornerlat(1)+89.975)/0.05)+1
     lon_off = nint((cornerlon(1)+179.975)/0.05)+1
     
     do r=1,GLASSalbedo_struc(n)%nr
        do c=1,GLASSalbedo_struc(n)%nc
           c1 = c + lon_off
           r1 = nr - (r + lat_off) + 1
           
           if(albedo_raw_avhrr(c1+(r1-1)*nc).ge.0) then 
              albedo_bs_in(c+(r-1)*GLASSalbedo_struc(n)%nc) =&
                   albedo_raw_avhrr(c1+(r1-1)*nc)*0.0001
              albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) =  .true. 
              if(GLASSalbedo_struc(n)%qcflag.eq.1) then 
                 if(ibits(albedo_qc(c1+(r1-1)*nc),0,1).gt.1) then 
                    albedo_bs_in(c+(r-1)*GLASSalbedo_struc(n)%nc) = -9999.0
                    albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) = .false. 
                 endif                 
              endif
           else
              albedo_bs_in(c+(r-1)*GLASSalbedo_struc(n)%nc) = -9999.0
              albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) = .false. 
           endif
        enddo
     enddo
!white sky albedo /diffuse

     iret = gdrdfld(grid_id,albedo_ws_name,start,stride,edge,albedo_raw_avhrr)
     
     albedo_data_b = .false. 
     
     lat_off = nint((cornerlat(1)+89.975)/0.05)+1
     lon_off = nint((cornerlon(1)+179.975)/0.05)+1
     
     do r=1,GLASSalbedo_struc(n)%nr
        do c=1,GLASSalbedo_struc(n)%nc
           c1 = c + lon_off
           r1 = nr - (r + lat_off) + 1
           
           if(albedo_raw_avhrr(c1+(r1-1)*nc).ge.0) then 
              albedo_ws_in(c+(r-1)*GLASSalbedo_struc(n)%nc) =&
                   albedo_raw_avhrr(c1+(r1-1)*nc)*0.0001
              albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) =  .true.
              if(GLASSalbedo_struc(n)%qcflag.eq.1) then 
                 if(ibits(albedo_qc(c1+(r1-1)*nc),0,1).gt.1) then 
                    albedo_ws_in(c+(r-1)*GLASSalbedo_struc(n)%nc) = -9999.0
                    albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) = .false. 
                 endif                 
              endif
 
           else
              albedo_ws_in(c+(r-1)*GLASSalbedo_struc(n)%nc) = -9999.0
              albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) = .false. 
           endif
        enddo
     enddo

  elseif(GLASSalbedo_struc(n)%source.eq."MODIS") then 
!black sky albedo
     allocate(albedo_raw_modis(nc*nr))
     iret = gdrdfld(grid_id,albedo_bs_name,start,stride,edge,albedo_raw_modis)
     
     albedo_data_b = .false. 
     
     lat_off = nint((cornerlat(1)+89.975)/0.05)+1
     lon_off = nint((cornerlon(1)+179.975)/0.05)+1
     
     do r=1,GLASSalbedo_struc(n)%nr
        do c=1,GLASSalbedo_struc(n)%nc
           c1 = c + lon_off
           r1 = nr - (r + lat_off) + 1
           
           if(albedo_raw_modis(c1+(r1-1)*nc).ge.0) then 
              albedo_bs_in(c+(r-1)*GLASSalbedo_struc(n)%nc) =&
                   albedo_raw_modis(c1+(r1-1)*nc)*0.0001
              albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) =  .true. 
              
              if(GLASSalbedo_struc(n)%qcflag.eq.1) then 
                 if(ibits(albedo_qc(c1+(r1-1)*nc),0,1).gt.1) then 
                    print*, ibits(albedo_qc(c1+(r1-1)*nc),0,1),albedo_bs_in(c+(r-1)*GLASSalbedo_struc(n)%nc)
                    albedo_bs_in(c+(r-1)*GLASSalbedo_struc(n)%nc) = -9999.0
                    albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) = .false. 
                 endif                 
              endif
           else
              albedo_bs_in(c+(r-1)*GLASSalbedo_struc(n)%nc) = -9999.0
              albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) = .false. 
              
           endif
        enddo
     enddo

!white sky albedo
     iret = gdrdfld(grid_id,albedo_ws_name,start,stride,edge,albedo_raw_modis)
     
     albedo_data_b = .false. 
     
     lat_off = nint((cornerlat(1)+89.975)/0.05)+1
     lon_off = nint((cornerlon(1)+179.975)/0.05)+1
     
     do r=1,GLASSalbedo_struc(n)%nr
        do c=1,GLASSalbedo_struc(n)%nc
           c1 = c + lon_off
           r1 = nr - (r + lat_off) + 1
           
           if(albedo_raw_modis(c1+(r1-1)*nc).ge.0) then 
              albedo_ws_in(c+(r-1)*GLASSalbedo_struc(n)%nc) =&
                   albedo_raw_modis(c1+(r1-1)*nc)*0.0001
              albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) =  .true. 
              if(GLASSalbedo_struc(n)%qcflag.eq.1) then 
                 if(ibits(albedo_qc(c1+(r1-1)*nc),0,1).gt.1) then 
                    albedo_ws_in(c+(r-1)*GLASSalbedo_struc(n)%nc) = -9999.0
                    albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) = .false. 
                 endif                 
              endif
           else
              albedo_ws_in(c+(r-1)*GLASSalbedo_struc(n)%nc) = -9999.0
              albedo_data_b(c+(r-1)*GLASSalbedo_struc(n)%nc) = .false. 
           endif
        enddo
     enddo

     deallocate(albedo_raw_modis)

  endif

  if(GLASSalbedo_struc(n)%qcflag.eq.1) then 
     deallocate(albedo_qc)
  endif
!  open(100,file='test_inp.bin',form='unformatted')
!  write(100) albedo_ws_in
!  write(100) albedo_bs_in
!  close(100)
!  stop

  iret=gddetach(grid_id)
  iret=gdclose(file_id)

  if(LIS_rc%obs_gridDesc(k,10).le.0.05) then 
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
          albedo_data_b, albedo_bs_in, albedoobs_b_ip, albedoobs_bs_ip, &
          GLASSalbedo_struc(n)%nc*GLASSalbedo_struc(n)%nr, &
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          GLASSalbedo_struc(n)%rlat,GLASSalbedo_struc(n)%rlon,&
          GLASSalbedo_struc(n)%w11,GLASSalbedo_struc(n)%w12,&
          GLASSalbedo_struc(n)%w21,GLASSalbedo_struc(n)%w22,&
          GLASSalbedo_struc(n)%n11,GLASSalbedo_struc(n)%n12,&
          GLASSalbedo_struc(n)%n21,GLASSalbedo_struc(n)%n22,LIS_rc%udef,iret)

     call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
          albedo_data_b, albedo_ws_in, albedoobs_b_ip, albedoobs_ws_ip, &
          GLASSalbedo_struc(n)%nc*GLASSalbedo_struc(n)%nr, &
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          GLASSalbedo_struc(n)%rlat,GLASSalbedo_struc(n)%rlon,&
          GLASSalbedo_struc(n)%w11,GLASSalbedo_struc(n)%w12,&
          GLASSalbedo_struc(n)%w21,GLASSalbedo_struc(n)%w22,&
          GLASSalbedo_struc(n)%n11,GLASSalbedo_struc(n)%n12,&
          GLASSalbedo_struc(n)%n21,GLASSalbedo_struc(n)%n22,LIS_rc%udef,iret)
  else
     call upscaleByAveraging(GLASSalbedo_struc(n)%nc*GLASSalbedo_struc(n)%nr,&
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          LIS_rc%udef, GLASSalbedo_struc(n)%n11,&
          albedo_data_b,albedo_bs_in, albedoobs_b_ip, albedoobs_bs_ip)
     call upscaleByAveraging(GLASSalbedo_struc(n)%nc*GLASSalbedo_struc(n)%nr,&
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          LIS_rc%udef, GLASSalbedo_struc(n)%n11,&
          albedo_data_b,albedo_ws_in, albedoobs_b_ip, albedoobs_ws_ip)
  endif
  
!  open(100,file='test_out.bin',form='unformatted')
!  write(100) albedoobs_bs_ip
!  write(100) albedoobs_ws_ip
!  close(100)
!  stop

#endif

end subroutine read_GLASS_ALBEDO_data




!BOP
! !ROUTINE: create_GLASSalbedo_filename
! \label{create_GLASSalbedo_filename}
! 
! !INTERFACE: 
subroutine create_GLASSalbedo_filename(ndir, source, yr, doy, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: source
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the GLASS ALBEDO filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the GLASS ALBEDO data directory
!  \item[yr]  current year
!  \item[mo]  current doy
!  \item[filename] Generated GLASS ALBEDO filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  if(source.eq."AVHRR") then 
     filename = trim(ndir)//'/'//trim(fyr)//'/GLASS02B05.V04.A'//&
          trim(fyr)//trim(fdoy)//'.hdf'
  elseif(source.eq."MODIS") then 
     filename = trim(ndir)//'/'//trim(fyr)//'/GLASS02B06.V04.A'//&
          trim(fyr)//trim(fdoy)//'.hdf'
  endif

end subroutine create_GLASSalbedo_filename





