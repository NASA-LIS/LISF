!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readJULESObs
! \label{readJULESObs}
!
! !INTERFACE: 
subroutine readJULESObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_logMod
  use LVT_histDataMod
  use JULES_obsMod 
  use map_utils

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!   This subroutine reads and extracts variables from JULES
!   output data for analysis within LVT. The variable names
!   are assumed to be in ALMA convention. 
!
! !NOTES: 
!   The code employs a neighbor lookup strategy to map the 
!   JULES data to the LVT analysis grid. No reprojection
!   is performed. 
!     
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  8 July 2015: Sujay Kumar, Initial Specification
! 
!EOP

  character*500           :: filename
  logical                 :: file_exists
  integer                 :: k,nid, ios
  integer                 :: timeid, tId, xId, yId, soilId
  integer                 :: latid, lonid
  integer                 :: rainfId, snowfId, qleId, qhId, qtauId, evapId, smcId, stcId, pstarId, tstarId,gppId

  real,  allocatable      :: rainf(:,:)
  real,  allocatable      :: snowf(:,:)
  real,  allocatable      :: qle(:,:)
  real,  allocatable      :: qh(:,:)
  real,  allocatable      :: qtau(:,:)
  real,  allocatable      :: evap(:,:)
  real,  allocatable      :: smc(:,:,:)
  real,  allocatable      :: stc(:,:,:)
  real,  allocatable      :: tstar(:,:)
  real,  allocatable      :: pstar(:,:)
  real,  allocatable      :: gpp(:,:)
  integer                 :: c,r,t,kk, tindex
  integer                 :: num_secs
  integer                 :: yr, mo, da, hr, mn, ss
  type(ESMF_Time)         :: currTime
  type(ESMF_TimeInterval) :: ts
  integer                 :: status
  integer                 :: stn_row, stn_col
  real                    :: col,row
  integer                 :: Rconst   !specific gas constant for dry air
  
#if (defined USE_NETCDF3 || defined USE_NETCDF4)  

!------------------------------------------------------------
!  check if the JULES file exists. If not, throw an error
!------------------------------------------------------------   
  filename = JULESobs(source)%odir
  inquire(file=trim(filename),exist=file_exists) 
  
  if(.not.file_exists) then 
     write(LVT_logunit,*) '[INFO] Finished reading JULES data '
     write(LVT_logunit,*) '[ERR] JULES file not found'
     write(LVT_logunit,*) '[ERR] Program stopping ..'
     call LVT_endrun()
  endif

  if(JULESobs(source)%startmode) then 
     JULESobs(source)%startMode = .false. 
     write(LVT_logunit,*) '[INFO] Reading JULES file ',trim(filename)
     
     ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios, 'Error opening file '//trim(filename))
     
     ! dimensions
     ios = nf90_inq_dimid(nid,'x',xId)
     call LVT_verify(ios, 'Error nf90_inq_dimid: x')
     
     ios = nf90_inquire_dimension(nid,xId, len=JULESobs(source)%nx)
     call LVT_verify(ios, 'Error nf90_inquire_dimension: x')
     
     ios = nf90_inq_dimid(nid,'y',yId)
     call LVT_verify(ios, 'Error nf90_inq_dimid: y')
     
     ios = nf90_inquire_dimension(nid,yId, len=JULESobs(source)%ny)
     call LVT_verify(ios, 'Error nf90_inquire_dimension: y')

     ios = nf90_inq_dimid(nid,'soil', soilId)
     if(ios.ne.0) then 
        write(LVT_logunit,*) '[WARN] soil dimension is missing in the JULES file'
     else
        JULESobs(source)%nsoil = 1
     endif
     ios = nf90_inquire_dimension(nid,soilId, len=JULESobs(source)%nsoil)
     !        call LVT_verify(ios, 'Error nf90_inquire_dimension: soil')

     ios = nf90_inq_dimid(nid,'time',tId)
     call LVT_verify(ios, 'Error nf90_inq_dimid: t')

     ios = nf90_inquire_dimension(nid,tId, len=JULESobs(source)%ntimes)
     call LVT_verify(ios, 'Error nf90_inquire_dimension: time')

     allocate(JULESobs(source)%time_val(JULESobs(source)%ntimes))
     allocate(JULESobs(source)%lat(JULESobs(source)%nx,JULESobs(source)%ny))
     allocate(JULESobs(source)%lon(JULESobs(source)%nx,JULESobs(source)%ny))

     if(LVT_MOC_RAINF(source).ge.1) then 
        allocate(JULESobs(source)%rainf_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%ntimes))
     endif
     if(LVT_MOC_SNOWF(source).ge.1) then 
        allocate(JULESobs(source)%snowf_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%ntimes))
     endif
     if(LVT_MOC_QLE(source).ge.1) then 
        allocate(JULESobs(source)%qle_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%ntimes))
     endif
     if(LVT_MOC_QH(source).ge.1) then 
        allocate(JULESobs(source)%qh_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%ntimes))
     endif
     if(LVT_MOC_QTAU(source).ge.1) then 
        allocate(JULESobs(source)%qtau_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%ntimes))
     endif
     if(LVT_MOC_EVAP(source).ge.1) then 
        allocate(JULESobs(source)%evap_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%ntimes))
     endif
     if(LVT_MOC_SOILMOIST(source).ge.1) then 
        allocate(JULESobs(source)%smc_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%nsoil,JULESobs(source)%ntimes))
     endif
     if(LVT_MOC_SOILTEMP(source).ge.1) then 
        allocate(JULESobs(source)%stc_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%nsoil,JULESobs(source)%ntimes))
     endif
     if(LVT_MOC_PSURFFORC(source).ge.1) then 
        allocate(JULESobs(source)%pstar_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%ntimes))
     endif
     if(LVT_MOC_AVGSURFT(source).ge.1) then 
        allocate(JULESobs(source)%tstar_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%ntimes))
     endif
     if(LVT_MOC_GPP(source).ge.1) then 
        allocate(JULESobs(source)%gpp_jules(JULESobs(source)%nx,JULESobs(source)%ny,JULESobs(source)%ntimes))
     endif


     !values
     ios = nf90_inq_varid(nid,'latitude',latid)
     call LVT_verify(ios, 'Error nf90_inq_varid: latitude')

     ios = nf90_get_var(nid,latid, JULESobs(source)%lat)
     call LVT_verify(ios, 'Error nf90_get_var: latitude')

     ios = nf90_inq_varid(nid,'longitude',lonid)
     call LVT_verify(ios, 'Error nf90_inq_varid: longitude')

     ios = nf90_get_var(nid,lonid, JULESobs(source)%lon)
     call LVT_verify(ios, 'Error nf90_get_var: longitude')

     ios = nf90_inq_varid(nid,'time',timeid)
     call LVT_verify(ios, 'Error nf90_inq_varid: time')

     ios = nf90_get_var(nid,timeid, JULESobs(source)%time_val)
     call LVT_verify(ios, 'Error nf90_get_var: time')


     if(LVT_MOC_RAINF(source).ge.1) then 
        !rainf
        ios = nf90_inq_varid(nid,'Rainf',rainfid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,rainfid, JULESobs(source)%rainf_jules, &
                start=(/1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Rainf')
        else
            JULESobs(source)%rainf_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_SNOWF(source).ge.1) then 
        !snowf
        ios = nf90_inq_varid(nid,'Snowf',snowfid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,snowfid, JULESobs(source)%snowf_jules, &
                start=(/1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Snowf')
        else
            JULESobs(source)%snowf_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_QLE(source).ge.1) then 
        !qle
        ios = nf90_inq_varid(nid,'Qle',qleid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,qleid, JULESobs(source)%qle_jules, &
                start=(/1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Qle')
        else
            JULESobs(source)%qle_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_QH(source).ge.1) then 
        !qh
        ios = nf90_inq_varid(nid,'Qh',qhid)
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,qhid, JULESobs(source)%qh_jules, &
                start=(/1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Qh')
        else
            JULESobs(source)%qh_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_QTAU(source).ge.1) then 
        !qtau
        ios = nf90_inq_varid(nid,'Qtau',qtauid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,qtauid, JULESobs(source)%qtau_jules, &
                start=(/1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Qtau')
        else
            JULESobs(source)%qtau_jules = LVT_rc%udef
        endif
     endif

     if(LVT_MOC_PSURFFORC(source).ge.1) then 
        !pstar
        ios = nf90_inq_varid(nid,'pstar',pstarid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,pstarid, JULESobs(source)%pstar_jules, &
                start=(/1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: pstar')
        else
            JULESobs(source)%pstar_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_AVGSURFT(source).ge.1) then 
        !tstar
        ios = nf90_inq_varid(nid,'tstar',tstarid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,tstarid, JULESobs(source)%tstar_jules, &
                start=(/1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: tstar')
        else
            JULESobs(source)%tstar_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_GPP(source).ge.1) then 
        !gpp
        ios = nf90_inq_varid(nid,'GPP',gppid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,gppid, JULESobs(source)%gpp_jules, &
                start=(/1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: GPP')
        else
            JULESobs(source)%gpp_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_EVAP(source).ge.1) then 
        !evap
        ios = nf90_inq_varid(nid,'Evap',evapid)
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,evapid, JULESobs(source)%evap_jules, &
                start=(/1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Evap')
        else
            JULESobs(source)%evap_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_SOILMOIST(source).ge.1) then 
        !soil moisture
        ios = nf90_inq_varid(nid,'SoilMoist',smcid)
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,smcid, JULESobs(source)%smc_jules, &
                start=(/1,1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%nsoil,JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: SoilMoist')
        else
            JULESobs(source)%smc_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_SOILTEMP(source).ge.1) then 
        !soil temperature
        ios = nf90_inq_varid(nid,'SoilTemp',stcid)
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,stcid, JULESobs(source)%stc_jules, &
                start=(/1,1,1,1/), &
                count=(/JULESobs(source)%nx,JULESobs(source)%ny,&
                JULESobs(source)%nsoil,JULESobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: SoilTemp')
        else
            JULESobs(source)%stc_jules = LVT_rc%udef
        endif
     endif

     ios = nf90_close(nid)
     call LVT_verify(ios, 'Error in nf90_close')
  end if

  allocate(rainf(LVT_rc%lnc, LVT_rc%lnr))
  allocate(snowf(LVT_rc%lnc, LVT_rc%lnr))
  allocate(qle(LVT_rc%lnc, LVT_rc%lnr))
  allocate(qh(LVT_rc%lnc, LVT_rc%lnr))
  allocate(qtau(LVT_rc%lnc, LVT_rc%lnr))
  allocate(evap(LVT_rc%lnc, LVT_rc%lnr))
  allocate(smc(LVT_rc%lnc, LVT_rc%lnr,JULESobs(source)%nsoil))
  allocate(stc(LVT_rc%lnc, LVT_rc%lnr,JULESobs(source)%nsoil))
  allocate(tstar(LVT_rc%lnc, LVT_rc%lnr))
  allocate(pstar(LVT_rc%lnc, LVT_rc%lnr))
  allocate(gpp(LVT_rc%lnc, LVT_rc%lnr))

  rainf = LVT_rc%udef
  snowf = LVT_rc%udef
  qle = LVT_rc%udef
  qh  = LVT_rc%udef
  qtau  = LVT_rc%udef
  evap  = LVT_rc%udef
  smc = LVT_rc%udef
  stc = LVT_rc%udef
  pstar = LVT_rc%udef
  tstar = LVT_rc%udef
  gpp = LVT_rc%udef

!------------------------------------------------------------
!  Find the time offset for current time. 
!------------------------------------------------------------   
  call ESMF_TimeSet(currTime, yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
       h = LVT_rc%dhr(source), m = LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'error in ESMF_TimeSet in readJULESobs')
  
  ts = currTime - JULESobs(source)%refTime
  
  call ESMF_TimeIntervalGet(ts, s=num_secs,rc=status)
  
  tindex = -1
  do t=1, JULESobs(source)%ntimes
     if(JULESobs(source)%time_val(t).eq.num_secs) then 
        tindex = t
        exit
     endif
  enddo
  
!------------------------------------------------------------
!   map the JULES data to the LVT grid           
!------------------------------------------------------------     
#if 0 
  if(tindex.gt.0) then 
     do c=1,JULESobs(source)%nx
        do r=1,JULESobs(source)%ny
           call latlon_to_ij(LVT_domain%lvtproj, &
                JULESobs(source)%lat(c,r), &
                JULESobs(source)%lon(c,r),&
                col,row)
           stn_col = nint(col)
           stn_row = nint(row)
           if(LVT_MOC_QLE(source).ge.1) then 
           !qle
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 qle(stn_col,stn_row) = &
                      JULESobs(source)%qle_jules(c,r,tindex)
              endif
           endif
           if(LVT_MOC_QH(source).ge.1) then 
           !qh
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 qh(stn_col,stn_row) = &
                      JULESobs(source)%qh_jules(c,r,tindex)
              endif
           endif
           if(LVT_MOC_QTAU(source).ge.1) then 
           !qtau
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 qtau(stn_col,stn_row) = &
                      JULESobs(source)%qtau_jules(c,r,tindex)
                 !print*,qtau(stn_col,stn_row)
              endif
           endif

           if(LVT_MOC_PSURFFORC(source).ge.1) then 
           !pstar
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 pstar(stn_col,stn_row) = &
                      JULESobs(source)%pstar_jules(c,r,tindex)
              endif
           endif

           if(LVT_MOC_AVGSURFT(source).ge.1) then 
           !tstar
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 tstar(stn_col,stn_row) = &
                      JULESobs(source)%tstar_jules(c,r,tindex)
              endif
           endif

           if(LVT_MOC_GPP(source).ge.1) then 
           !gpp
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 gpp(stn_col,stn_row) = &
                      JULESobs(source)%gpp_jules(c,r,tindex)


           if(LVT_MOC_EVAP(source).ge.1) then 
           !evap
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 evap(stn_col,stn_row) = &
                      JULESobs(source)%evap_jules(c,r,tindex)
              endif
           endif
           if(LVT_MOC_SOILMOIST(source).ge.1) then 
           !soil moisture
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 smc(stn_col,stn_row,:) = &
                      JULESobs(source)%smc_jules(c,r,:,tindex)
              endif
           endif
           if(LVT_MOC_SOILTEMP(source).ge.1) then 
           !soil temperature
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 stc(stn_col,stn_row,:) = &
                      JULESobs(source)%stc_jules(c,r,:,tindex)
              endif
           endif
        enddo
     enddo
  endif
#endif
  if(tindex.gt.0) then 
     do c=1,JULESobs(source)%nx
        do r=1,JULESobs(source)%ny
           if(LVT_MOC_QLE(source).ge.1) then 
           !qle
              qle(:,:) = & 
                   JULESobs(source)%qle_jules(c,r,tindex)
           endif
           if(LVT_MOC_QH(source).ge.1) then 
              qh(:,:) = &
                   JULESobs(source)%qh_jules(c,r,tindex)
           endif
           if(LVT_MOC_QTAU(source).ge.1) then 
           !qtau
              qtau(:,:) = & 
                   JULESobs(source)%qtau_jules(c,r,tindex)
           endif
           if(LVT_MOC_PSURFFORC(source).ge.1) then 
           !pstar
              pstar(:,:) = & 
                   JULESobs(source)%pstar_jules(c,r,tindex)
           endif
           if(LVT_MOC_AVGSURFT(source).ge.1) then 
           !tstar
              tstar(:,:) = & 
                   JULESobs(source)%tstar_jules(c,r,tindex)
           endif
           if(LVT_MOC_GPP(source).ge.1) then 
           !gpp
              gpp(:,:) = & 
                   JULESobs(source)%gpp_jules(c,r,tindex)
           endif

           if(LVT_MOC_EVAP(source).ge.1) then 
              evap(:,:) = &
                   JULESobs(source)%evap_jules(c,r,tindex)
           endif
           if(LVT_MOC_SOILMOIST(source).ge.1) then 
              do k=1,JULESobs(source)%nsoil
                 smc(:,:,k) = &
                      JULESobs(source)%smc_jules(c,r,k,tindex)
              enddo
           endif
           if(LVT_MOC_SOILTEMP(source).ge.1) then 
           !soil temperature
              do k=1,JULESobs(source)%nsoil
                 stc(:,:,k) = &
                      JULESobs(source)%stc_jules(c,r,k,tindex)
              enddo
           endif
        enddo
     enddo
  endif
#endif
     
!------------------------------------------------------------
! log the processed variables to LVT 
!------------------------------------------------------------

  call LVT_logSingleDataStreamVar(LVT_MOC_RAINF,source,&
       rainf,vlevel=1,units="kg/m2s")

  call LVT_logSingleDataStreamVar(LVT_MOC_SNOWF,source,&
       snowf,vlevel=1,units="kg/m2s")

  call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,&
       qle,vlevel=1,units="W/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_QH,source,&
       qh,vlevel=1,units="W/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_QTAU,source,&
       qtau,vlevel=1,units="kg/m/s2")

  call LVT_logSingleDataStreamVar(LVT_MOC_PSURFFORC,source,&
       pstar,vlevel=1,units="Pa")

  call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT,source,&
       tstar,vlevel=1,units="K")
 
  !Convert tau to friction velocity
  do r=1, LVT_rc%lnr
      do c=1, LVT_rc%lnc
          if(qtau(c,r).ne.LVT_rc%udef) then 
              ! tau/density of air
              
              Rconst = 287.05
              qtau(c,r) = SQRT(qtau(c,r)/((pstar(c,r)/(Rconst*tstar(c,r)))))  
          else
              qtau(c,r) = LVT_rc%udef
          endif
      enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_QTAU,source,&
       qtau,vlevel=1,units="m/s")

  !GPP
  call LVT_logSingleDataStreamVar(LVT_MOC_GPP,source,&
       gpp,vlevel=1,units="kg/m2s")
  !Convert gpp to umolC m-2 s-1
  !Assuming the molar mass of carbon 12.0107g/mol
  do r=1, LVT_rc%lnr
    do c=1, LVT_rc%lnc
        
        if(gpp(c,r).ne.LVT_rc%udef) then
            gpp(c,r) = (gpp(c,r)/0.0120107)*1000000.0
        else
            gpp(c,r) = LVT_rc%udef
        endif
    enddo
  enddo
  
  call LVT_logSingleDataStreamVar(LVT_MOC_GPP,source,&
      gpp,vlevel=1,units="umol/m2s")

  call LVT_logSingleDataStreamVar(LVT_MOC_EVAP,source,&
       evap,vlevel=1,units="kg/m2s")

  do k=1,JULESobs(source)%nsoil
     call LVT_logSingleDataStreamVar(LVT_MOC_SoilMoist,source,&
          smc(:,:,k),vlevel=k,units="kg/m2")
  enddo

  do k=1,JULESobs(source)%nsoil
     do r=1, LVT_rc%lnr
        do c=1, LVT_rc%lnc
           if(smc(c,r,k).ne.LVT_rc%udef) then 
              ! density of water*layer depth
              if(k.eq.1) then
                 smc(c,r,k) = smc(c,r,k)/(1000.0*0.1) 
              elseif(k.eq.2) then 
                 smc(c,r,k) = smc(c,r,k)/(1000.0*0.25)
              elseif(k.eq.3) then 
                 smc(c,r,k) = smc(c,r,k)/(1000.0*0.65)
              elseif(k.eq.4) then 
                 smc(c,r,k) = smc(c,r,k)/(1000.0*2.0)
              endif
           else
              smc(c,r,k) = LVT_rc%udef
           endif
        enddo
     enddo
  enddo

  do k=1,JULESobs(source)%nsoil
     call LVT_logSingleDataStreamVar(LVT_MOC_SoilMoist,source,&
          smc(:,:,k),vlevel=k,units="m3/m3")
  enddo

  do k=1,JULESobs(source)%nsoil
     call LVT_logSingleDataStreamVar(LVT_MOC_SoilTemp,source,&
          stc(:,:,k),vlevel=k,units="K")
  enddo
  
  deallocate(rainf)
  deallocate(snowf)
  deallocate(qle)
  deallocate(qh)
  deallocate(qtau)
  deallocate(evap)
  deallocate(smc)
  deallocate(stc)
  deallocate(pstar)
  deallocate(tstar)
  deallocate(gpp)

end subroutine readJULESObs

