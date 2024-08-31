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
! !ROUTINE: readJULES2DObs
! \label{readJULES2DObs}
!
! !INTERFACE: 
subroutine readJULES2DObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_logMod
  use LVT_histDataMod
  use JULES2D_obsMod
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
!  08 Jul 2015: Sujay Kumar, Initial Specification
!  26 Apr 2018: Abheera Hazra, JULES native resolution mapping 
!                                and interpolation to LVT grid
! 
!EOP

  character*500           :: filename
  logical                 :: file_exists
  integer                 :: k,nid, ios,l,iret
  integer                 :: timeid, tId, xId, yId, soilId
  integer                 :: latid, lonid
  integer                 :: rainfId, snowfId, qleId, qhId, qtauId, evapId, smcId, stcId, pstarId, tstarId,gppId,ndviId
  integer		  :: albedoId, albedodirvisId, albedodifvisId, albedodirnirId, albedodifnirId 
  real,  allocatable      :: rainf(:,:), rainfJ(:,:), rainfJ1(:)
  real,  allocatable      :: snowf(:,:), snowfJ(:,:), snowfJ1(:)
  real,  allocatable      :: qle(:,:), qleJ(:,:), qleJ1(:)
  real,  allocatable      :: qh(:,:), qhJ(:,:), qhJ1(:)
  real,  allocatable      :: qtau(:,:), qtauJ(:,:), qtauJ1(:)
  real,  allocatable      :: evap(:,:), evapJ(:,:), evapJ1(:)
  real,  allocatable      :: smc(:,:,:),smc1(:,:), smcJ(:,:,:), smcJ1(:,:),smcJ2(:) 
  real,  allocatable      :: stc(:,:,:),stc1(:,:), stcJ(:,:,:), stcJ1(:,:),stcJ2(:)
  real,  allocatable      :: tstar(:,:),tstarJ(:,:), tstarJ1(:)
  real,  allocatable      :: pstar(:,:), pstarJ(:,:), pstarJ1(:)
  real,  allocatable      :: ndvi(:,:), ndviJ(:,:), ndviJ1(:)
  real,  allocatable      :: albedo(:,:), albedoJ(:,:), albedoJ1(:)
  real,  allocatable      :: albedodirvis(:,:), albedodirvisJ(:,:), albedodirvisJ1(:)
  real,  allocatable      :: albedodifvis(:,:), albedodifvisJ(:,:), albedodifvisJ1(:)
  real,  allocatable      :: albedodirnir(:,:), albedodirnirJ(:,:), albedodirnirJ1(:)
  real,  allocatable      :: albedodifnir(:,:), albedodifnirJ(:,:), albedodifnirJ1(:)
  real,  allocatable      :: gpp(:,:), gppJ(:,:), gppJ1(:)
  real	 		  :: dx,dy,latmin,latmax,lonmin,lonmax
  integer                 :: c,r,t,kk, tindex, i  
  integer                 :: num_secs
  integer                 :: yr, mo, da, hr, mn, ss
  type(ESMF_Time)         :: currTime
  type(ESMF_TimeInterval) :: ts
  integer                 :: status
  integer                 :: stn_row, stn_col
  real                    :: col,row
  integer                 :: Rconst   !specific gas constant for dry air
  logical*1, allocatable  :: lo(:)
  integer   		  :: ila,ilo,cnt1,cnt2,nla,nlo  
  real, allocatable 	  :: latJ(:,:),lonJ(:,:),mskla(:),msklo(:)
  integer, allocatable	  :: laJ(:,:),loJ(:,:) 

#if (defined USE_NETCDF3 || defined USE_NETCDF4)  

!------------------------------------------------------------
!  check if the JULES file exists. If not, throw an error
!------------------------------------------------------------   
										
						
  filename = JULES2Dobs(source)%odir
  inquire(file=trim(filename),exist=file_exists) 
  
 
  if(.not.file_exists) then 
     write(LVT_logunit,*) '[INFO] Finished reading JULES data '
     write(LVT_logunit,*) '[ERR] JULES file not found'
     write(LVT_logunit,*) '[ERR] Program stopping ..'
     call LVT_endrun()
  endif


  if(JULES2Dobs(source)%startmode) then 
!     JULES2Dobs(source)%startMode = .false.
     write(LVT_logunit,*) '[INFO] Reading JULES file ',trim(filename)
     
     ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios, 'Error opening file '//trim(filename))
     
     ! dimensions
     ios = nf90_inq_dimid(nid,'x',xId)
     call LVT_verify(ios, 'Error nf90_inq_dimid: x')
     
     ios = nf90_inquire_dimension(nid,xId, len=JULES2Dobs(source)%nx)
     call LVT_verify(ios, 'Error nf90_inquire_dimension: x')
     
     ios = nf90_inq_dimid(nid,'y',yId)
     call LVT_verify(ios, 'Error nf90_inq_dimid: y')
     
     ios = nf90_inquire_dimension(nid,yId, len=JULES2Dobs(source)%ny)
     call LVT_verify(ios, 'Error nf90_inquire_dimension: y')

     ios = nf90_inq_dimid(nid,'soil', soilId)
     if(ios.ne.0) then 
        write(LVT_logunit,*) '[WARN] soil dimension is missing in the JULES file'
     else
        JULES2Dobs(source)%nsoil = 1
     endif
     ios = nf90_inquire_dimension(nid,soilId, len=JULES2Dobs(source)%nsoil)
     !        call LVT_verify(ios, 'Error nf90_inquire_dimension: soil')

     ios = nf90_inq_dimid(nid,'time',tId)
     call LVT_verify(ios, 'Error nf90_inq_dimid: t')

     ios = nf90_inquire_dimension(nid,tId, len=JULES2Dobs(source)%ntimes)
     call LVT_verify(ios, 'Error nf90_inquire_dimension: time')

     allocate(JULES2Dobs(source)%time_val(JULES2Dobs(source)%ntimes))
     allocate(JULES2Dobs(source)%lat(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny))
     allocate(JULES2Dobs(source)%lon(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny))

     if(LVT_MOC_RAINF(source).ge.1) then 
        allocate(JULES2Dobs(source)%rainf_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_SNOWF(source).ge.1) then 
        allocate(JULES2Dobs(source)%snowf_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_QLE(source).ge.1) then 
        allocate(JULES2Dobs(source)%qle_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_QH(source).ge.1) then 
        allocate(JULES2Dobs(source)%qh_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_QTAU(source).ge.1) then 
        allocate(JULES2Dobs(source)%qtau_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_EVAP(source).ge.1) then 
        allocate(JULES2Dobs(source)%evap_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_SOILMOIST(source).ge.1) then 
        allocate(JULES2Dobs(source)%smc_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%nsoil,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_SOILTEMP(source).ge.1) then 
        allocate(JULES2Dobs(source)%stc_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%nsoil,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_PSURFFORC(source).ge.1) then 
        allocate(JULES2Dobs(source)%pstar_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_AVGSURFT(source).ge.1) then 
        allocate(JULES2Dobs(source)%tstar_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
!     if(LVT_MOC_LAI(source).ge.1) then 
!        allocate(JULES2Dobs(source)%lai_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
!     endif
     if(LVT_MOC_NDVI(source).ge.1) then 
        allocate(JULES2Dobs(source)%ndvi_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
!     if(LVT_MOC_SNOWCOVER(source).ge.1) then 
!        allocate(JULES2Dobs(source)%snowm_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
!     endif
     if(LVT_MOC_ALBEDO(source).ge.1) then 
        allocate(JULES2Dobs(source)%albedo_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_VISDIRALBEDO(source).ge.1) then 
        allocate(JULES2Dobs(source)%albedodirvis_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_VISDIFALBEDO(source).ge.1) then 
        allocate(JULES2Dobs(source)%albedodifvis_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_NIRDIRALBEDO(source).ge.1) then 
        allocate(JULES2Dobs(source)%albedodirnir_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif
     if(LVT_MOC_NIRDIFALBEDO(source).ge.1) then 
        allocate(JULES2Dobs(source)%albedodifnir_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif

!     if(LVT_MOC_EMISSFORC(source).ge.1) then 
!        allocate(JULES2Dobs(source)%emis_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
!     endif
!     if(LVT_MOC_RFLO(source).ge.1) then 
!        allocate(JULES2Dobs(source)%ndvi_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
!     endif
     if(LVT_MOC_GPP(source).ge.1) then 
        allocate(JULES2Dobs(source)%gpp_jules(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,JULES2Dobs(source)%ntimes))
     endif


     !values
     ios = nf90_inq_varid(nid,'latitude',latid)
     call LVT_verify(ios, 'Error nf90_inq_varid: latitude')

     ios = nf90_get_var(nid,latid, JULES2Dobs(source)%lat)
     call LVT_verify(ios, 'Error nf90_get_var: latitude')

     ios = nf90_inq_varid(nid,'longitude',lonid)
     call LVT_verify(ios, 'Error nf90_inq_varid: longitude')

     ios = nf90_get_var(nid,lonid, JULES2Dobs(source)%lon)
     call LVT_verify(ios, 'Error nf90_get_var: longitude')

     ios = nf90_inq_varid(nid,'time',timeid)
     call LVT_verify(ios, 'Error nf90_inq_varid: time')

     ios = nf90_get_var(nid,timeid, JULES2Dobs(source)%time_val)
     call LVT_verify(ios, 'Error nf90_get_var: time')



     if(LVT_MOC_RAINF(source).ge.1) then 
        !rainf
        ios = nf90_inq_varid(nid,'Rainf',rainfid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,rainfid, JULES2Dobs(source)%rainf_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Rainf')
        else
            JULES2Dobs(source)%rainf_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_SNOWF(source).ge.1) then 
        !snowf
        ios = nf90_inq_varid(nid,'Snowf',snowfid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,snowfid, JULES2Dobs(source)%snowf_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Snowf')
        else
            JULES2Dobs(source)%snowf_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_QLE(source).ge.1) then 
        !qle
        ios = nf90_inq_varid(nid,'Qle',qleid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,qleid, JULES2Dobs(source)%qle_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Qle')
        else
            JULES2Dobs(source)%qle_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_QH(source).ge.1) then 
        !qh
        ios = nf90_inq_varid(nid,'Qh',qhid)
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,qhid, JULES2Dobs(source)%qh_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Qh')
        else
            JULES2Dobs(source)%qh_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_QTAU(source).ge.1) then 
        !qtau
        ios = nf90_inq_varid(nid,'Qtau',qtauid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,qtauid, JULES2Dobs(source)%qtau_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Qtau')
        else
            JULES2Dobs(source)%qtau_jules = LVT_rc%udef
        endif
     endif

     if(LVT_MOC_PSURFFORC(source).ge.1) then 
        !pstar
        ios = nf90_inq_varid(nid,'pstar',pstarid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,pstarid, JULES2Dobs(source)%pstar_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: pstar')
        else
            JULES2Dobs(source)%pstar_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_AVGSURFT(source).ge.1) then 
        !tstar
        ios = nf90_inq_varid(nid,'AvgSurfT',tstarid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,tstarid, JULES2Dobs(source)%tstar_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: tstar')
        else
            JULES2Dobs(source)%tstar_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_GPP(source).ge.1) then 
        !gpp
        ios = nf90_inq_varid(nid,'GPP',gppid)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,gppid, JULES2Dobs(source)%gpp_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: GPP')
        else
            JULES2Dobs(source)%gpp_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_EVAP(source).ge.1) then 
        !evap
        ios = nf90_inq_varid(nid,'Evap',evapid)
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,evapid, JULES2Dobs(source)%evap_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: Evap')
        else
            JULES2Dobs(source)%evap_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_SOILMOIST(source).ge.1) then 
        !soil moisture
        ios = nf90_inq_varid(nid,'SoilMoist',smcid)
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,smcid, JULES2Dobs(source)%smc_jules, &
                start=(/1,1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%nsoil,JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: SoilMoist')
        else
            JULES2Dobs(source)%smc_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_SOILTEMP(source).ge.1) then 
        !soil temperature
        ios = nf90_inq_varid(nid,'SoilTemp',stcid)
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,stcid, JULES2Dobs(source)%stc_jules, &
                start=(/1,1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%nsoil,JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: SoilTemp')
        else
            JULES2Dobs(source)%stc_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_NDVI(source).ge.1) then 
        !NDVI
        ios = nf90_inq_varid(nid,'ndvi',ndviId)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,ndviId, JULES2Dobs(source)%ndvi_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: ndvi')
        else
            JULES2Dobs(source)%ndvi_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_ALBEDO(source).ge.1) then 
        !ALBEDO
        ios = nf90_inq_varid(nid,'Albedo_Total',albedoId)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,albedoId, JULES2Dobs(source)%albedo_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: albedo')
        else
            JULES2Dobs(source)%albedo_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_VISDIFALBEDO(source).ge.1) then 
        !ALBEDO_DIF_VIS
        ios = nf90_inq_varid(nid,'AlbedoDiffuseVisible',albedodifvisId)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,albedodifvisId, JULES2Dobs(source)%albedodifvis_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: albedodifvis')
        else
            JULES2Dobs(source)%albedodifvis_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_VISDIRALBEDO(source).ge.1) then 
        !ALBEDO_DIR_VIS
        ios = nf90_inq_varid(nid,'AlbedoDirectVisible',albedodirvisId)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,albedodirvisId, JULES2Dobs(source)%albedodirvis_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: albedodirvis')
        else
            JULES2Dobs(source)%albedodirvis_jules = LVT_rc%udef
        endif
     endif

     if(LVT_MOC_NIRDIFALBEDO(source).ge.1) then 
        !ALBEDO_DIF_NIR
        ios = nf90_inq_varid(nid,'AlbedoDiffuseNIR',albedodifnirId)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,albedodifnirId, JULES2Dobs(source)%albedodifnir_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: albedodifnir')
        else
            JULES2Dobs(source)%albedodifnir_jules = LVT_rc%udef
        endif
     endif
     if(LVT_MOC_NIRDIRALBEDO(source).ge.1) then 
        !ALBEDO_DIR_NIR
        ios = nf90_inq_varid(nid,'AlbedoDirectNIR',albedodirnirId)
        
        if(ios.eq.0) then 
           ios = nf90_get_var(nid,albedodirnirId, JULES2Dobs(source)%albedodirnir_jules, &
                start=(/1,1,1/), &
                count=(/JULES2Dobs(source)%nx,JULES2Dobs(source)%ny,&
                JULES2Dobs(source)%ntimes/))
           call LVT_verify(ios, 'Error nf90_get_var: albedodirnir')
        else
            JULES2Dobs(source)%albedodirnir_jules = LVT_rc%udef
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
  allocate(ndvi(LVT_rc%lnc, LVT_rc%lnr))
  allocate(albedo(LVT_rc%lnc, LVT_rc%lnr))
  allocate(albedodirvis(LVT_rc%lnc, LVT_rc%lnr))
  allocate(albedodifvis(LVT_rc%lnc, LVT_rc%lnr))
  allocate(albedodirnir(LVT_rc%lnc, LVT_rc%lnr))
  allocate(albedodifnir(LVT_rc%lnc, LVT_rc%lnr))
  allocate(smc(LVT_rc%lnc, LVT_rc%lnr,JULES2Dobs(source)%nsoil))
  allocate(stc(LVT_rc%lnc, LVT_rc%lnr,JULES2Dobs(source)%nsoil))
  allocate(tstar(LVT_rc%lnc, LVT_rc%lnr))
  allocate(pstar(LVT_rc%lnc, LVT_rc%lnr))
  allocate(gpp(LVT_rc%lnc, LVT_rc%lnr))
  allocate(smc1(LVT_rc%lnc,LVT_rc%lnr)) 
  allocate(stc1(LVT_rc%lnc,LVT_rc%lnr)) 
  allocate(mskla(JULES2Dobs(source)%nx))			
  allocate(msklo(JULES2Dobs(source)%nx))			
  allocate(latJ(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny))       
  allocate(lonJ(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny))	
  allocate(laJ(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny))       
  allocate(loJ(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny))	

  mskla = LVT_rc%udef						
  msklo = LVT_rc%udef						
  rainf = LVT_rc%udef
  snowf = LVT_rc%udef
  qle = LVT_rc%udef
  qh  = LVT_rc%udef
  qtau  = LVT_rc%udef
  evap  = LVT_rc%udef
  ndvi  = LVT_rc%udef
  albedo  = LVT_rc%udef
  albedodirvis  = LVT_rc%udef
  albedodifvis  = LVT_rc%udef
  albedodirnir  = LVT_rc%udef
  albedodifnir  = LVT_rc%udef
  smc = LVT_rc%udef
  smc1 = LVT_rc%udef
  stc = LVT_rc%udef
  stc1 = LVT_rc%udef
  pstar = LVT_rc%udef
  tstar = LVT_rc%udef
  gpp = LVT_rc%udef
  latJ(:,:)=JULES2Dobs(source)%lat
  lonJ(:,:)=JULES2Dobs(source)%lon

!------------------------------------------------------------
!  Find the time offset for current time. 
!------------------------------------------------------------   

  call ESMF_TimeSet(currTime, yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
       h = LVT_rc%dhr(source), m = LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'error in ESMF_TimeSet in readJULESobs')
  
  ts = currTime - JULES2Dobs(source)%refTime
  
  call ESMF_TimeIntervalGet(ts, s=num_secs,rc=status)
  
  tindex = -1
  do t=1, JULES2Dobs(source)%ntimes
     if(JULES2Dobs(source)%time_val(t).eq.num_secs) then 
        tindex = t
        exit
     endif
  enddo

!------------------------------------------------------------
!   maps the JULES2D data at its native resolution  (Abheera)  
!------------------------------------------------------------     


  if(JULES2Dobs(source)%startMode .and.JULES2Dobs(source)%nx.gt.1 .and. JULES2Dobs(source)%ny.ge.1 ) then 

     JULES2Dobs(source)%startMode = .false.
     allocate(JULES2Dobs(source)%latJ1(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny))
     allocate(JULES2Dobs(source)%lonJ1(JULES2Dobs(source)%nx,JULES2Dobs(source)%ny))	
     allocate(JULES2Dobs(source)%latc(JULES2Dobs(source)%nx))			
     allocate(JULES2Dobs(source)%lonc(JULES2Dobs(source)%nx))			


     ila=1
     ilo=1


!---------------------------------------------------------------------
!!! extracts the longitudes and finds out the number of longitudes
!---------------------------------------------------------------------

     do c=1,JULES2Dobs(source)%nx
     	if(c.eq.1) then
	   JULES2Dobs(source)%lonc(1) = MINVAL(lonJ)
	   msklo(1) = MINVAL(lonJ)
	   where(lonJ .eq. JULES2Dobs(source)%lonc(1))
	      JULES2Dobs(source)%lonJ1=c
           end where	
	else if(c.gt.1)then
	   ilo=ilo+1
	   if((ANY(lonJ.gt.MAXVAL(msklo(1:ilo-1))).eqv..true.))then
	      JULES2Dobs(source)%lonc(ilo) = MINVAL(lonJ,mask=lonJ.gt.MAXVAL(msklo(1:ilo-1)))
	      msklo(ilo) = JULES2Dobs(source)%lonc(ilo)
	      where(lonJ .eq. JULES2Dobs(source)%lonc(ilo))
	     	 JULES2Dobs(source)%lonJ1=c
              end where	
	      JULES2Dobs(source)%clo=ilo
	   endif	     	     
	endif	   	   
     enddo

     write(LVT_logunit,*) 'Number of JULES Longitudes ',JULES2Dobs(source)%clo

!----------------------------------------------------------------
!!! extracts the latitudes and finds out the number of latitudes
!----------------------------------------------------------------

     do r=1,JULES2Dobs(source)%nx
	if(r.eq.1) then
	   JULES2Dobs(source)%latc(1) = MINVAL(latJ)
	   mskla(1) = MINVAL(latJ)
	   where(latJ .eq. JULES2Dobs(source)%latc(1))
	      JULES2Dobs(source)%latJ1=r
           end where	 
	else if(r.gt.1)then
	   ila=ila+1
	   if(ANY(latJ.gt.MAXVAL(mskla(1:ila-1))).eqv..true.)then
	      JULES2Dobs(source)%latc(ila) = MINVAL(latJ,mask=latJ.gt.MAXVAL(mskla(1:ila-1)))
	      mskla(ila) = JULES2Dobs(source)%latc(ila)
	      where(latJ .eq. JULES2Dobs(source)%latc(ila))
	     	  JULES2Dobs(source)%latJ1=r
              end where	
	      JULES2Dobs(source)%cla=ila
	   endif
	endif
     enddo

     write(LVT_logunit,*) 'Number of JULES Latitudes ',JULES2Dobs(source)%cla

     allocate(JULES2Dobs(source)%latf(JULES2Dobs(source)%cla))
     allocate(JULES2Dobs(source)%lonf(JULES2Dobs(source)%clo))
  
     JULES2Dobs(source)%latf(1:JULES2Dobs(source)%cla)=JULES2Dobs(source)%latc(1:JULES2Dobs(source)%cla)
     JULES2Dobs(source)%lonf(1:JULES2Dobs(source)%clo)=JULES2Dobs(source)%lonc(1:JULES2Dobs(source)%clo)

!------------------------------------------------------------
!!!   Run function bilinear_interp_input         
!------------------------------------------------------------    
     allocate(JULES2DObs(source)%rlat(LVT_rc%lnc*LVT_rc%lnr))
     allocate(JULES2DObs(source)%rlon(LVT_rc%lnc*LVT_rc%lnr))
     allocate(JULES2DObs(source)%n11(LVT_rc%lnc*LVT_rc%lnr))
     allocate(JULES2DObs(source)%n12(LVT_rc%lnc*LVT_rc%lnr))
     allocate(JULES2DObs(source)%n21(LVT_rc%lnc*LVT_rc%lnr))
     allocate(JULES2DObs(source)%n22(LVT_rc%lnc*LVT_rc%lnr))
     allocate(JULES2DObs(source)%w11(LVT_rc%lnc*LVT_rc%lnr))
     allocate(JULES2DObs(source)%w12(LVT_rc%lnc*LVT_rc%lnr))
     allocate(JULES2DObs(source)%w21(LVT_rc%lnc*LVT_rc%lnr))
     allocate(JULES2DObs(source)%w22(LVT_rc%lnc*LVT_rc%lnr))
    
     JULES2DObs(source)%gridDesci    = 0
     JULES2DObs(source)%gridDesci(1) = 0
     JULES2DObs(source)%gridDesci(2) = JULES2Dobs(source)%clo
     JULES2DObs(source)%gridDesci(3) = JULES2Dobs(source)%cla
     JULES2DObs(source)%gridDesci(4) = MINVAL(JULES2Dobs(source)%latf)
     JULES2DObs(source)%gridDesci(5) = MINVAL(JULES2Dobs(source)%lonf)
     JULES2DObs(source)%gridDesci(7) = MAXVAL(JULES2Dobs(source)%latf)
     JULES2DObs(source)%gridDesci(8) = MAXVAL(JULES2Dobs(source)%lonf)
     JULES2DObs(source)%gridDesci(6) = 128
     JULES2DObs(source)%gridDesci(9) = 0.5
     JULES2DObs(source)%gridDesci(10) = 0.5
     JULES2DObs(source)%gridDesci(20) = 64

     call bilinear_interp_input(JULES2DObs(source)%gridDesci,LVT_rc%gridDesc,  &
           LVT_rc%lnc*LVT_rc%lnr,                            &
           JULES2DObs(source)%rlat, JULES2DObs(source)%rlon,     &
           JULES2DObs(source)%n11, JULES2DObs(source)%n12,       &
           JULES2DObs(source)%n21, JULES2DObs(source)%n22,       &
           JULES2DObs(source)%w11, JULES2DObs(source)%w12,       &
           JULES2DObs(source)%w21, JULES2DObs(source)%w22)


  end if

!------------------------------------------------------------
!!!   Allocate variables to JULES2D native grids          
!------------------------------------------------------------     

  if(tindex.gt.0.and.JULES2Dobs(source)%nx.gt.1.and.JULES2Dobs(source)%ny.ge.1) then

     laJ(:,:)=JULES2Dobs(source)%latJ1(:,:)
     loJ(:,:)=JULES2Dobs(source)%lonJ1(:,:)
     nla=JULES2Dobs(source)%cla
     nlo=JULES2Dobs(source)%clo

     allocate(rainfJ(nlo,nla))
     allocate(snowfJ(nlo,nla))
     allocate(qleJ(nlo,nla))
     allocate(qhJ(nlo,nla))
     allocate(qtauJ(nlo,nla))
     allocate(evapJ(nlo,nla))
     allocate(ndviJ(nlo,nla))
     allocate(albedoJ(nlo,nla))
     allocate(albedodirvisJ(nlo,nla))
     allocate(albedodifvisJ(nlo,nla))
     allocate(albedodirnirJ(nlo,nla))
     allocate(albedodifnirJ(nlo,nla))
     allocate(smcJ(nlo,nla,JULES2Dobs(source)%nsoil))
     allocate(stcJ(nlo,nla,JULES2Dobs(source)%nsoil))
     allocate(tstarJ(nlo,nla))   
     allocate(pstarJ(nlo,nla))
     allocate(gppJ(nlo,nla))

     allocate(rainfJ1(nlo*nla))
     allocate(snowfJ1(nlo*nla))
     allocate(qleJ1(nlo*nla))
     allocate(qhJ1(nlo*nla))
     allocate(qtauJ1(nlo*nla))
     allocate(evapJ1(nlo*nla))
     allocate(ndviJ1(nlo*nla))
     allocate(albedoJ1(nlo*nla))
     allocate(albedodirvisJ1(nlo*nla))
     allocate(albedodifvisJ1(nlo*nla))
     allocate(albedodirnirJ1(nlo*nla))
     allocate(albedodifnirJ1(nlo*nla))
     allocate(smcJ1(nlo*nla,JULES2Dobs(source)%nsoil))
     allocate(smcJ2(nlo*nla))
     allocate(stcJ1(nlo*nla,JULES2Dobs(source)%nsoil))
     allocate(stcJ2(nlo*nla))
     allocate(tstarJ1(nlo*nla))   
     allocate(pstarJ1(nlo*nla))
     allocate(gppJ1(nlo*nla))

     rainfJ = LVT_rc%udef
     snowfJ = LVT_rc%udef
     qleJ = LVT_rc%udef
     qhJ = LVT_rc%udef
     qtauJ = LVT_rc%udef
     evapJ = LVT_rc%udef
     ndviJ = LVT_rc%udef
     albedoJ = LVT_rc%udef
     albedodirvisJ = LVT_rc%udef
     albedodifvisJ = LVT_rc%udef
     albedodirnirJ = LVT_rc%udef
     albedodifnirJ = LVT_rc%udef
     smcJ = LVT_rc%udef
     stcJ = LVT_rc%udef
     tstarJ = LVT_rc%udef
     pstarJ = LVT_rc%udef
     gppJ = LVT_rc%udef

     rainfJ1 = LVT_rc%udef
     snowfJ1 = LVT_rc%udef
     qleJ1 = LVT_rc%udef
     qhJ1 = LVT_rc%udef
     qtauJ1 = LVT_rc%udef
     evapJ1 = LVT_rc%udef
     ndviJ1 = LVT_rc%udef
     albedoJ1 = LVT_rc%udef
     albedodirvisJ1 = LVT_rc%udef
     albedodifvisJ1 = LVT_rc%udef
     albedodirnirJ1 = LVT_rc%udef
     albedodifnirJ1 = LVT_rc%udef
     smcJ1 = LVT_rc%udef
     stcJ1 = LVT_rc%udef
     tstarJ1 = LVT_rc%udef
     pstarJ1 = LVT_rc%udef
     gppJ1 = LVT_rc%udef

  endif

!------------------------------------------------------------
!   map the JULES data to relevant JULES grids          
!------------------------------------------------------------     

 
#if 0 
  if(tindex.gt.0) then 
     do c=1,JULES2Dobs(source)%nx
        do r=1,JULES2Dobs(source)%ny
           call latlon_to_ij(LVT_domain%lvtproj, &
                JULES2Dobs(source)%lat(c,r), &
                JULES2Dobs(source)%lon(c,r),&
                col,row)
           stn_col = nint(col)
           stn_row = nint(row)
           if(LVT_MOC_QLE(source).ge.1) then 
           !qle
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 qle(stn_col,stn_row) = &
                      JULES2Dobs(source)%qle_jules(c,r,tindex)
              endif
           endif
               if(LVT_MOC_QH(source).ge.1) then 
           !qh
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 qh(stn_col,stn_row) = &
                      JULES2Dobs(source)%qh_jules(c,r,tindex)
              endif
           endif
           if(LVT_MOC_QTAU(source).ge.1) then 
           !qtau
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 qtau(stn_col,stn_row) = &
                      JULES2Dobs(source)%qtau_jules(c,r,tindex)
                 !print*,qtau(stn_col,stn_row)
              endif
           endif

           if(LVT_MOC_PSURFFORC(source).ge.1) then 
           !pstar
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 pstar(stn_col,stn_row) = &
                      JULES2Dobs(source)%pstar_jules(c,r,tindex)
              endif
           endif

           if(LVT_MOC_AVGSURFT(source).ge.1) then 
           !tstar
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 tstar(stn_col,stn_row) = &
                      JULES2Dobs(source)%tstar_jules(c,r,tindex)
              endif
           endif

           if(LVT_MOC_GPP(source).ge.1) then 
           !gpp
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 gpp(stn_col,stn_row) = &
                      JULES2Dobs(source)%gpp_jules(c,r,tindex)


           if(LVT_MOC_EVAP(source).ge.1) then 
           !evap
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 evap(stn_col,stn_row) = &
                      JULES2Dobs(source)%evap_jules(c,r,tindex)
              endif
           endif
           if(LVT_MOC_SOILMOIST(source).ge.1) then 
           !soil moisture
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 smc(stn_col,stn_row,:) = &
                      JULES2Dobs(source)%smc_jules(c,r,:,tindex)
              endif
           endif
           if(LVT_MOC_SOILTEMP(source).ge.1) then 
           !soil temperature
              if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
                   stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
                 stc(stn_col,stn_row,:) = &
                      JULES2Dobs(source)%stc_jules(c,r,:,tindex)
              endif
           endif
        enddo
     enddo
  endif
#endif 

  if(tindex.gt.0.and.JULES2Dobs(source)%nx.gt.1.and.JULES2Dobs(source)%ny.ge.1) then 
     do c=1,JULES2Dobs(source)%nx
        do r=1,JULES2Dobs(source)%ny
           if(LVT_MOC_QLE(source).ge.1) then 
           !qle
              qleJ(loJ(c,r),laJ(c,r)) = & 
                   JULES2Dobs(source)%qle_jules(c,r,tindex)
           endif
           if(LVT_MOC_QH(source).ge.1) then 
              qhJ(loJ(c,r),laJ(c,r)) = &
                   JULES2Dobs(source)%qh_jules(c,r,tindex)
           endif
           if(LVT_MOC_QTAU(source).ge.1) then 
           !qtau
              qtauJ(loJ(c,r),laJ(c,r)) = & 
                   JULES2Dobs(source)%qtau_jules(c,r,tindex)
           endif
           if(LVT_MOC_PSURFFORC(source).ge.1) then 
           !pstar
              pstarJ(loJ(c,r),laJ(c,r)) = & 
                   JULES2Dobs(source)%pstar_jules(c,r,tindex)
           endif
           if(LVT_MOC_AVGSURFT(source).ge.1) then 
           !tstar
              tstarJ(loJ(c,r),laJ(c,r)) = & 
                   JULES2Dobs(source)%tstar_jules(c,r,tindex)
           endif
           if(LVT_MOC_GPP(source).ge.1) then 
           !gpp
              gppJ(loJ(c,r),laJ(c,r)) = & 
                   JULES2Dobs(source)%gpp_jules(c,r,tindex)
           endif
	   !evap
           if(LVT_MOC_EVAP(source).ge.1) then 
              evapJ(loJ(c,r),laJ(c,r)) = &
                   JULES2Dobs(source)%evap_jules(c,r,tindex)
           endif
	   !ndvi
           if(LVT_MOC_NDVI(source).ge.1) then 
              ndviJ(loJ(c,r),laJ(c,r)) = &
                   JULES2Dobs(source)%ndvi_jules(c,r,tindex)
           endif
	   !albedo
           if(LVT_MOC_ALBEDO(source).ge.1) then 
              albedoJ(loJ(c,r),laJ(c,r)) = &
                   JULES2Dobs(source)%albedo_jules(c,r,tindex)
           endif
	   !albedo_dir_vis
           if(LVT_MOC_VISDIRALBEDO(source).ge.1) then 
              albedodirvisJ(loJ(c,r),laJ(c,r)) = &
                   JULES2Dobs(source)%albedodirvis_jules(c,r,tindex)
           endif
	   !albedo_dif_vis
           if(LVT_MOC_VISDIFALBEDO(source).ge.1) then 
              albedodifvisJ(loJ(c,r),laJ(c,r)) = &
                   JULES2Dobs(source)%albedodifvis_jules(c,r,tindex)
           endif
	   !albedo_dir_nir
           if(LVT_MOC_NIRDIRALBEDO(source).ge.1) then 
              albedodirnirJ(loJ(c,r),laJ(c,r)) = &
                   JULES2Dobs(source)%albedodirnir_jules(c,r,tindex)
           endif
	   !albedo_dif_nir
           if(LVT_MOC_NIRDIFALBEDO(source).ge.1) then 
              albedodifnirJ(loJ(c,r),laJ(c,r)) = &
                   JULES2Dobs(source)%albedodifnir_jules(c,r,tindex)
           endif
	   !soil moist
           if(LVT_MOC_SoilMoist(source).ge.1) then
              do k=1,JULES2Dobs(source)%nsoil
                   smcJ(loJ(c,1),laJ(c,1),k) = &
		      JULES2Dobs(source)%smc_jules(c,r,k,tindex)
	      enddo
           endif
           if(LVT_MOC_SOILTEMP(source).ge.1) then 
           !soil temperature
              do k=1,JULES2Dobs(source)%nsoil
                   stcJ(loJ(c,r),laJ(c,r),k) = &
                      JULES2Dobs(source)%stc_jules(c,r,k,tindex)
              enddo
           endif
        enddo
     enddo

!-------------------------------------------
!JULES 2D to nc*nr
!-------------------------------------------
     do r=1,nla
        do c=1,nlo
           if(LVT_MOC_QLE(source).ge.1) then 
           !qle
              qleJ1(c+(r-1)*nlo)=qleJ(c,r)
           endif
           if(LVT_MOC_QH(source).ge.1) then 
           !qh
              qhJ1(c+(r-1)*nlo)=qhJ(c,r)
           endif
           if(LVT_MOC_QTAU(source).ge.1) then 
           !qtau
              qtauJ1(c+(r-1)*nlo)=qtauJ(c,r)
           endif
           if(LVT_MOC_PSURFFORC(source).ge.1) then 
           !pstar
              pstarJ1(c+(r-1)*nlo)=pstarJ(c,r)
           endif
           if(LVT_MOC_AVGSURFT(source).ge.1) then 
           !tstar
              tstarJ1(c+(r-1)*nlo)=tstarJ(c,r)
           endif
           if(LVT_MOC_GPP(source).ge.1) then 
           !GPP
              gppJ1(c+(r-1)*nlo)=gppJ(c,r)
           endif
           if(LVT_MOC_EVAP(source).ge.1) then 
           !evap
              evapJ1(c+(r-1)*nlo)=evapJ(c,r)
           endif
           if(LVT_MOC_NDVI(source).ge.1) then 
           !ndvi
              ndviJ1(c+(r-1)*nlo)=ndviJ(c,r)
           endif
           if(LVT_MOC_ALBEDO(source).ge.1) then 
           !albedo
              albedoJ1(c+(r-1)*nlo)=albedoJ(c,r)
           endif
           if(LVT_MOC_VISDIRALBEDO(source).ge.1) then 
           !albedo_dir_vis
              albedodirvisJ1(c+(r-1)*nlo)=albedodirvisJ(c,r)
           endif
           if(LVT_MOC_VISDIFALBEDO(source).ge.1) then 
           !albedo_dif_vis
              albedodifvisJ1(c+(r-1)*nlo)=albedodifvisJ(c,r)
           endif
           if(LVT_MOC_NIRDIRALBEDO(source).ge.1) then 
           !albedo_dir_nir
              albedodirnirJ1(c+(r-1)*nlo)=albedodirnirJ(c,r)
           endif
           if(LVT_MOC_NIRDIFALBEDO(source).ge.1) then 
           !albedo_dif_nir
              albedodifnirJ1(c+(r-1)*nlo)=albedodifnirJ(c,r)
           endif
           if(LVT_MOC_SoilMoist(source).ge.1) then 
           !soil_moist
              smcJ1(c+(r-1)*nlo,1:JULES2Dobs(source)%nsoil)=smcJ(c,r,1:JULES2Dobs(source)%nsoil) 
           endif
           if(LVT_MOC_SOILTEMP(source).ge.1) then 
           !soil_temp
              stcJ1(c+(r-1)*nlo,1:JULES2Dobs(source)%nsoil)=stcJ(c,r,1:JULES2Dobs(source)%nsoil)
           endif
        enddo
     enddo

!------------------------------------------------------------
! Interpolate and log processed data to LVT 
!------------------------------------------------------------


!    call interp_jules2var(source,nlo,nla,rainfJ1,rainf)
     call LVT_logSingleDataStreamVar(LVT_MOC_RAINF,source,&
       rainf,vlevel=1,units="kg/m2s")

!    call interp_jules2var(source,nlo,nla,snowfJ1,snowf)
     call LVT_logSingleDataStreamVar(LVT_MOC_SNOWF,source,&
       snowf,vlevel=1,units="kg/m2s")

     if(LVT_MOC_QLE(source).ge.1) then
        call interp_jules2var(source,nlo,nla,qleJ1,qle)
     endif
     call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,&
        qle,vlevel=1,units="W/m2")


     if(LVT_MOC_QH(source).ge.1) then
        call interp_jules2var(source,nlo,nla,qhJ1,qh)
     endif
     call LVT_logSingleDataStreamVar(LVT_MOC_QH,source,&
        qh,vlevel=1,units="W/m2")


     if(LVT_MOC_EVAP(source).ge.1) then
        call interp_jules2var(source,nlo,nla,evapJ1,evap)
     endif
     call LVT_logSingleDataStreamVar(LVT_MOC_EVAP,source,&
        evap,vlevel=1,units="kg/m2s")


     if(LVT_MOC_NDVI(source).ge.1) then
        call interp_jules2var(source,nlo,nla,ndviJ1,ndvi)
     endif
     call LVT_logSingleDataStreamVar(LVT_MOC_NDVI,source,&
        ndvi,vlevel=1,units="-")


     if(LVT_MOC_ALBEDO(source).ge.1) then
        call interp_jules2var(source,nlo,nla,albedoJ1,albedo)
     endif
     call LVT_logSingleDataStreamVar(LVT_MOC_ALBEDO,source,&
        albedo,vlevel=1,units="-")


     if(LVT_MOC_VISDIRALBEDO(source).ge.1) then
        call interp_jules2var(source,nlo,nla,albedodirvisJ1,albedodirvis)
     endif
     call LVT_logSingleDataStreamVar(LVT_MOC_VISDIRALBEDO,source,&
        albedodirvis,vlevel=1,units="-")

     if(LVT_MOC_VISDIFALBEDO(source).ge.1) then
        call interp_jules2var(source,nlo,nla,albedodifvisJ1,albedodifvis)
     endif
     call LVT_logSingleDataStreamVar(LVT_MOC_VISDIFALBEDO,source,&
        albedodifvis,vlevel=1,units="-")


    if(LVT_MOC_NIRDIRALBEDO(source).ge.1) then
        call interp_jules2var(source,nlo,nla,albedodirnirJ1,albedodirnir)
     endif
     call LVT_logSingleDataStreamVar(LVT_MOC_NIRDIRALBEDO,source,&
        albedodirnir,vlevel=1,units="-")

     if(LVT_MOC_NIRDIFALBEDO(source).ge.1) then
        call interp_jules2var(source,nlo,nla,albedodifnirJ1,albedodifnir)
     endif
     call LVT_logSingleDataStreamVar(LVT_MOC_NIRDIFALBEDO,source,&
        albedodifnir,vlevel=1,units="-")


     if(LVT_MOC_AVGSURFT(source).ge.1) then
        call interp_jules2var(source,nlo,nla,tstarJ1,tstar)
     endif
     call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT,source,&
        tstar,vlevel=1,units="K")


     do k=1,JULES2Dobs(source)%nsoil
        if(LVT_MOC_SoilMoist(source).ge.1) then 
           smcJ2=smcJ1(:,k)
           call interp_jules2var(source,nlo,nla,smcJ2,smc1)
           smc(:,:,k)=smc1(:,:)
        endif
        call LVT_logSingleDataStreamVar(LVT_MOC_SoilMoist,source,&
           smc(:,:,k),vlevel=k,units="kg/m2")
     enddo

     do k=1,JULES2Dobs(source)%nsoil
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

     do k=1,JULES2Dobs(source)%nsoil
        call LVT_logSingleDataStreamVar(LVT_MOC_SoilMoist,source,&
           smc(:,:,k),vlevel=k,units="m3/m3")
     enddo

     do k=1,JULES2Dobs(source)%nsoil
        if(LVT_MOC_SoilTemp(source).ge.1) then 
           stcJ2=stcJ1(:,k)
           call interp_jules2var(source,nlo,nla,stcJ2,stc1)
           stc(:,:,k)=stc1(:,:)
        endif
        call LVT_logSingleDataStreamVar(LVT_MOC_SoilTemp,source,&
           stc(:,:,k),vlevel=k,units="K")
     enddo
  endif
#endif


     
!------------------------------------------------------------
! deallocating
!------------------------------------------------------------

  deallocate(rainf)
  deallocate(snowf)
  deallocate(qle)
  deallocate(qh)
  deallocate(qtau)
  deallocate(evap)
  deallocate(ndvi)
  deallocate(albedo)
  deallocate(albedodirvis)
  deallocate(albedodifvis)
  deallocate(albedodirnir)
  deallocate(albedodifnir)
  deallocate(smc)
  deallocate(stc)
  deallocate(pstar)
  deallocate(tstar)
  deallocate(gpp)
  deallocate(smc1)
  deallocate(stc1)
  deallocate(mskla)
  deallocate(msklo)
  deallocate(latJ)
  deallocate(lonJ)
  deallocate(laJ)
  deallocate(loJ)

  if(tindex .gt.0.and.JULES2Dobs(source)%nx.gt.1.and.JULES2Dobs(source)%ny.ge.1) then 
     deallocate(rainfJ)
     deallocate(rainfJ1)
     deallocate(snowfJ)
     deallocate(snowfJ1)
     deallocate(qleJ)
     deallocate(qleJ1)
     deallocate(qhJ)
     deallocate(qhJ1)
     deallocate(qtauJ)
     deallocate(qtauJ1)
     deallocate(evapJ)
     deallocate(evapJ1)
     deallocate(ndviJ)
     deallocate(ndviJ1)
     deallocate(albedoJ)
     deallocate(albedoJ1)
     deallocate(albedodirvisJ)
     deallocate(albedodirvisJ1)
     deallocate(albedodifnirJ)
     deallocate(albedodifnirJ1)
     deallocate(smcJ)
     deallocate(smcJ1)
     deallocate(smcJ2)
     deallocate(stcJ)
     deallocate(stcJ1)
     deallocate(stcJ2)
     deallocate(pstarJ)
     deallocate(pstarJ1)
     deallocate(tstarJ)
     deallocate(tstarJ1)
     deallocate(gppJ)  
     deallocate(gppJ1)
  endif

end subroutine readJULES2DObs

!BOP
!
! !ROUTINE: interp_jules2var
!  \label{interp_nldas2var}
!
! !INTERFACE:
  subroutine interp_jules2var(source, nc,nr,var_input,var_output)
!
! !USES:
    use LVT_coreMod,    only : LVT_rc
    use JULES2D_obsMod      
    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This subroutine spatially interpolates the JULES variable to the
!   model (LIS) output grid and resolution, using a bilinear interpolation
!   approach.
!
!   The arguments are:
!   \begin{description}
!    \item[nc]      number of columns in the input (JULES) grid
!    \item[nr]      number of rows in the input (JULES) grid
!    \item[var_input] input variable to be interpolated
!    \item[lb]        input bitmap (true//false)
!    \item[var_output] resulting interpolated field
!   \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
!BOP
!
! !ARGUMENTS:
    integer            :: source 
    integer            :: nc
    integer            :: nr
    real               :: var_input(nc*nr)
    logical*1          :: lb(nc*nr)
    real               :: var_output(LVT_rc%lnc, LVT_rc%lnr)
    !EOP
    integer            :: iret
    integer            :: c,r
    logical*1          :: lo(LVT_rc%lnc*LVT_rc%lnr)
    real               :: go(LVT_rc%lnc*LVT_rc%lnr)

     var_output = LVT_rc%udef
     lb = .false.
     do r = 1,nr
        do c = 1,nc
           if (var_input(c+(r-1)*nc).ne.LVT_rc%udef) then
              lb(c+(r-1)*nc) = .true.
           endif
        enddo
     enddo

!
     call bilinear_interp(LVT_rc%gridDesc,lb,var_input,               &
        lo,go,nc*nr,LVT_rc%lnc*LVT_rc%lnr,          &
        JULES2Dobs(source)%rlat,JULES2Dobs(source)%rlon,            &
        JULES2Dobs(source)%w11,JULES2Dobs(source)%w12,              &
        JULES2Dobs(source)%w21,JULES2Dobs(source)%w22,              &
        JULES2Dobs(source)%n11,JULES2Dobs(source)%n12,              &
        JULES2Dobs(source)%n21,JULES2Dobs(source)%n22,              &
        LVT_rc%udef,iret)

     do r = 1,LVT_rc%lnr
        do c = 1,LVT_rc%lnc
           var_output(c,r) = go(c+(r-1)*LVT_rc%lnc)
        enddo
     enddo

  end subroutine interp_jules2var




