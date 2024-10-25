!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_NUWRFprocMod
!BOP
!
! !MODULE: LDT_NUWRFprocMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manage NUWRF related
!   preprocessing
!
! !REVISION HISTORY: 
!  13 Feb 2014:  Sujay Kumar;  Initial Specification
! 
  use ESMF
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_NUWRFprocInit   
  public :: LDT_procDataForREAL

!EOP
  type, public :: nuwrfdec
     character(len=LDT_CONST_PATH_LEN)        :: LIShistfile
     character(len=LDT_CONST_PATH_LEN)        :: realfile
     type(LDT_paramEntry) :: avgsurft
  end type nuwrfdec

  type(nuwrfdec), allocatable :: LDT_NUWRF_struc(:)
contains

!BOP
! !ROUTINE: LDT_NUWRFprocInit
! \label{LDT_NUWRFprocInit}
!
! !INTERFACE: 
  subroutine LDT_NUWRFprocInit()

! !USES:
    use LDT_coreMod, only : LDT_rc, LDT_config
! 
! !DESCRIPTION: 
!  This subroutine reads the configuration attributes related to 
!  the NUWRF preprocessing mode. 
! 
!EOP
    integer   :: n 
    integer   :: rc
    integer   :: k
! ____________________________________________

    allocate(LDT_NUWRF_struc(LDT_rc%nnest))
    call ESMF_ConfigFindLabel(LDT_config,&
         "LIS history file for land state initialization:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            LDT_NUWRF_struc(n)%LIShistfile,rc=rc)
       call LDT_verify(rc,&
            "LIS history file for land state initialization: not defined")
    enddo

    call ESMF_ConfigFindLabel(LDT_config,&
         "Processed NUWRF file for input to real:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            LDT_NUWRF_struc(n)%realfile,rc=rc)
       call LDT_verify(rc,&
            "Processed NUWRF file for input to real: not defined")
    enddo

  end subroutine LDT_NUWRFprocInit

!BOP
! !ROUTINE: LDT_procDataForREAL
! \label{LDT_procDataForREAL}
!
! !INTERFACE: 
  subroutine LDT_procDataForREAL(n)

! !USES:
    use LDT_coreMod, only : LDT_rc, LDT_config
    use LDT_logMod,  only : LDT_logunit, LDT_verify
! !ARGUMENTS: 
    integer,   intent(in) :: n
!
! !DESCRIPTION: 
!  This subroutine generates the blended file for use with the
!  real data initialization program (real.exe) in NUWRF. The 
!  routine reads the LDT generated parameter file (that was
!  used in a LIS spinup and a LIS history file to generate the 
!  blended file. 
!
!  The following variables are read from the parameter file: 
!   Landcover
!   Landcover fraction
!   Soiltype
!   Terrain elevation
!   Greenness monthly climatology
!   Greenness min/max fields
!
!  The following variables are read from the given LIS history file
!   Snow cover (instantaneous)
!   Snow depth (instantaneous)
!   SWE (instantaneous)
!   Canopy interception (instantaneous)
!   Albedo (instantaneous)
!   Bottom temperature (instantaneous)
!   Average surface temperature (time averaged)
!   Green vegetation fraction (instantaneous)
!   Total soil moisture (instantaneous)
!   Liquid fraction of soil moisture (instantaneous)
!   Soil temperature (instantaneous)
!   Relative soil moisture (instantaneous)
! 
!EOP

    integer               :: ftn_real, ftn_lis, ftn_ldt
    integer               :: dimID(3)
    integer               :: ncId, nrId, vegId, soilId,elevId,maskId
    integer               :: smpId,nsmp,stpId,nstp,slpId,nslp
    integer               :: monthId, nmonth,nrsmc,relsmcId
    integer               :: shdminId, shdmaxId,lu1id
    integer               :: lc1Id,slt1Id,elev1Id,mask1Id
    integer               :: shdmin1Id, shdmax1Id
    integer               :: scaId,sca1Id,sweid,swe1id
    integer               :: snodid,snod1id,avgsurftid,avgsurft1id
    integer               :: canopintid,canopint1id,albid,alb1id
    integer               :: gvfid,gvf1id,rsmcid,rsmc1id
    integer               :: gvfinstid, gvfinst1id,tbotid,tbot1Id
    integer               :: smcId,smc1Id,sh2oId,sh2o1Id,stcid,stc1id
    integer               :: mxalbId, mxalb1Id,mxsnalb1Id
    integer               :: nc,nr,nsfctypes,nslttypes
    integer               :: ios,iret
    integer               :: lcid,textureId,calbId,calb1Id
    real                  :: maxv, maxt
    integer               :: c,r,t
    integer               :: varid
    character*20          :: units
    character*100         :: standard_name
    real                  :: scale_factor
    real                  :: offset
    real                  :: vmin
    real                  :: vmax
    character*50          :: lcscheme
    integer               :: nsoillayers
    real, allocatable         :: lyrthk(:)
    real, allocatable         :: lc(:,:,:),texture(:,:,:),mxsnalb(:,:)
    real, allocatable         :: mask(:,:),ltype(:,:), stype(:,:),elev(:,:)
    real, allocatable         :: shdmin(:,:),shdmax(:,:),relsmc(:,:,:)
    real, allocatable         :: sca(:,:),snod(:,:),swe(:,:),avgsurft(:,:)
    real, allocatable         :: canopint(:,:),alb(:,:), gvf_inst(:,:),tbot(:,:)
    real, allocatable         :: smc(:,:,:), stc(:,:,:),sh2o(:,:,:),gvf(:,:,:)
    real, allocatable         :: albclimo(:,:,:)

! read from LSM parameter file, process dominant fields, write 
! to real input file
#if(defined USE_NETCDF3 || defined USE_NETCDF4)        

!write the merged input file for real.exe
#if(defined USE_NETCDF3)
     iret=nf90_create(path=trim(LDT_NUWRF_struc(n)%realfile),&
          cmode=nf90_clobber, ncid=ftn_real)
#else
    iret=nf90_create(path=trim(LDT_NUWRF_struc(n)%realfile),&
         cmode=nf90_netcdf4, ncid=ftn_real)
#endif

    call LDT_verify(iret,'creating netcdf file failed in LDT_paramProcWrite')
    
    call LDT_verify(nf90_def_dim(ftn_real,'east_west',LDT_rc%gnc(n),dimID(1)))
    call LDT_verify(nf90_def_dim(ftn_real,'north_south',LDT_rc%gnr(n),dimID(2)))
    call LDT_verify(nf90_def_var(ftn_real,'Landmask',&
         nf90_float, dimids=dimID(1:2),varid=mask1Id),&
         'nf90_def_var failed for Landmask')
    call LDT_verify(nf90_def_var(ftn_real,'Landcover',&
         nf90_float, dimids=dimID(1:2),varid=lc1Id),&
         'nf90_def_var failed for Landcover')
    call LDT_verify(nf90_def_var(ftn_real,'Soiltype',&
         nf90_float, dimids=dimID(1:2),varid=slt1Id),&
         'nf90_def_var failed for Soiltype')
    call LDT_verify(nf90_def_var(ftn_real,'Elevation',&
         nf90_float, dimids=dimID(1:2),varid=elev1Id),&
         'nf90_def_var failed for Elevation')
    call LDT_verify(nf90_def_var(ftn_real,'Shdmax',&
         nf90_float, dimids=dimID(1:2),varid=shdmax1Id),&
         'nf90_def_var failed for Shdmax')
    call LDT_verify(nf90_def_var(ftn_real,'Shdmin',&
         nf90_float, dimids=dimID(1:2),varid=shdmin1Id),&
         'nf90_def_var failed for Shdmin')
    call LDT_verify(nf90_def_var(ftn_real,'Mxsnowalbedo',&
         nf90_float, dimids=dimID(1:2),varid=mxsnalb1Id),&
         'nf90_def_var failed for Mxsnowalbedo')
    call LDT_verify(nf90_def_var(ftn_real,'Snowcover',&
         nf90_float, dimids=dimID(1:2),varid=sca1Id),&
         'nf90_def_var failed for Snowcover')
    call LDT_verify(nf90_def_var(ftn_real,'Snowdepth',&
         nf90_float, dimids=dimID(1:2),varid=snod1Id),&
         'nf90_def_var failed for Snowdepth')
    call LDT_verify(nf90_def_var(ftn_real,'SWE',&
         nf90_float, dimids=dimID(1:2),varid=swe1Id),&
         'nf90_def_var failed for SWE')
    call LDT_verify(nf90_def_var(ftn_real,'CanopInt',&
         nf90_float, dimids=dimID(1:2),varid=canopint1Id),&
         'nf90_def_var failed for CanopInt')
    call LDT_verify(nf90_def_var(ftn_real,'AvgSurfT',&
         nf90_float, dimids=dimID(1:2),varid=avgsurft1Id),&
         'nf90_def_var failed for AvgSurfT')
    call LDT_verify(nf90_def_var(ftn_real,'Albedo_inst',&
         nf90_float, dimids=dimID(1:2),varid=alb1Id),&
         'nf90_def_var failed for Albedo_inst')
    call LDT_verify(nf90_def_var(ftn_real,'Tempbot',&
         nf90_float, dimids=dimID(1:2),varid=tbot1Id),&
         'nf90_def_var failed for Tempbot')
    call LDT_verify(nf90_def_var(ftn_real,'Greenness_inst',&
         nf90_float, dimids=dimID(1:2),varid=gvfinst1Id),&
         'nf90_def_var failed for Greenness_inst')

    write(LDT_logunit,*) 'Reading '//trim(LDT_LSMparam_struc(n)%param_filename)
    call LDT_verify(nf90_open(path=LDT_LSMparam_struc(n)%param_filename, &
         mode=NF90_NOWRITE, ncid = ftn_ldt),&
         'Error opening file '//trim(LDT_LSMparam_struc(n)%param_filename))

    ios = nf90_inq_dimid(ftn_ldt,"east_west",ncId)
    call LDT_verify(ios,'Error in nf90_inq_dimid in LDT_procDataForREAL')
    
    ios = nf90_inq_dimid(ftn_ldt,"north_south",nrId)
    call LDT_verify(ios,'Error in nf90_inq_dimid in LDT_procDataForREAL')

    ios = nf90_inq_dimid(ftn_ldt,"sfctypes",vegId)
    call LDT_verify(ios,'Error in nf90_inq_dimid in LDT_procDataForREAL')

    ios = nf90_inq_dimid(ftn_ldt,"soiltypes",soilId)
    call LDT_verify(ios,'Error in nf90_inq_dimid in LDT_procDataForREAL')

    ios = nf90_inq_dimid(ftn_ldt,"month",monthId)
    call LDT_verify(ios,'Error in nf90_inq_dimid in LDT_procDataForREAL')

    ios = nf90_inquire_dimension(ftn_ldt,ncId, len=nc)
    call LDT_verify(ios,'Error in nf90_inquire_dimension in LDT_procDataForREAL')
    
    ios = nf90_inquire_dimension(ftn_ldt,nrId, len=nr)
    call LDT_verify(ios,'Error in nf90_inquire_dimension in LDT_procDataForREAL')
    
     ios = nf90_inquire_dimension(ftn_ldt,vegId, len=nsfctypes)
     call LDT_verify(ios,'Error in nf90_inquire_dimension in LDT_procDataForREAL')

     ios = nf90_inquire_dimension(ftn_ldt,monthId, len=nmonth)
     call LDT_verify(ios,'Error in nf90_inquire_dimension in LDT_procDataForREAL')

     ios = nf90_inquire_dimension(ftn_ldt,soilId, len=nslttypes)
     call LDT_verify(ios,'Error in nf90_inquire_dimension in LDT_procDataForREAL')

    if(nc.ne.LDT_rc%gnc(n).or.nr.ne.LDT_rc%gnr(n)) then 
       write(LDT_logunit,*) 'The domain dimensions in the Processed LDT parameter file'
       write(LDT_logunit,*) 'do not match the domain dimensions in the LDT'
       write(LDT_logunit,*) 'config file. Program stopping...'
       call LDT_endrun()
    endif
    
    allocate(mask(nc,nr))

    ios = nf90_inq_varid(ftn_ldt,'LANDMASK',maskid)
    call LDT_verify(ios,'LANDMASK field not found in the LIS param file')
    
    ios = nf90_get_var(ftn_ldt,maskid,mask)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_ldt,maskid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,mask1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(lc(nc,nr,nsfctypes))
    
    call LDT_verify(nf90_def_dim(ftn_real,'Nvegtypes',nsfctypes,dimID(3)))
    call LDT_verify(nf90_def_var(ftn_real,'Landusef',&
         nf90_float, dimids=dimID(1:3),varid=lu1Id),&
         'nf90_def_var failed for Landusef')

    ios = nf90_inq_varid(ftn_ldt,'LANDCOVER',lcid)
    call LDT_verify(ios,'LANDCOVER field not found in the LIS param file')
    
    ios = nf90_get_var(ftn_ldt,lcid,lc)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')
    
    call getVarAttribs(ftn_ldt,lcid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,lc1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,lu1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(ltype(nc,nr))

    ltype = LDT_rc%udef
    do r=1,nr
       do c=1,nc
          maxt = -1
          maxv = -10.0

          do t=1,nsfctypes
             if(lc(c,r,t).gt.maxv) then 
                maxv = lc(c,r,t)
                maxt = t
             endif
          enddo
          ltype(c,r) = maxt
       enddo
    enddo

!    allocate(texture(nc,nr,nsfctypes))
    allocate(texture(nc,nr,nslttypes))
    
    ios = nf90_inq_varid(ftn_ldt,'TEXTURE',textureid)
    call LDT_verify(ios,'TEXTURE field not found in the LIS param file')
    
    ios = nf90_get_var(ftn_ldt,textureid,texture)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    allocate(stype(nc,nr))
    stype = LDT_rc%udef
    do r=1,nr
       do c=1,nc
          maxt = -1
          maxv = -10.0
!          do t=1,nsfctypes
          do t=1,nslttypes
             if(texture(c,r,t).gt.maxv) then 
                maxv = texture(c,r,t)
                maxt = t
             endif
          enddo
          stype(c,r) = maxt    
       enddo
    enddo

    call getVarAttribs(ftn_ldt,textureid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,slt1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(elev(nc,nr))

    ios = nf90_inq_varid(ftn_ldt,'ELEVATION',elevid)
    call LDT_verify(ios,'ELEVATION field not found in the LIS param file')
    
    ios = nf90_get_var(ftn_ldt,elevid,elev)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_ldt,elevid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,elev1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    
    allocate(shdmax(nc,nr))
    ios = nf90_inq_varid(ftn_ldt,'SHDMAX',shdmaxid)
    call LDT_verify(ios,'SHDMAX field not found in the LIS param file')
    
    ios = nf90_get_var(ftn_ldt,shdmaxid,shdmax)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_ldt,shdmaxid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,shdmax1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(shdmin(nc,nr))
    ios = nf90_inq_varid(ftn_ldt,'SHDMIN',shdminid)
    call LDT_verify(ios,'SHDMIN field not found in the LIS param file')
    
    ios = nf90_get_var(ftn_ldt,shdminid,shdmin)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_ldt,shdminid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,shdmin1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(gvf(nc,nr,nmonth))
    call LDT_verify(nf90_def_dim(ftn_real,'month',nmonth,dimID(3)))
    call LDT_verify(nf90_def_var(ftn_real,'Greenness',&
         nf90_float, dimids=dimID(1:3),varid=gvf1Id),&
         'nf90_def_var failed for Greenness')

    ios = nf90_inq_varid(ftn_ldt,'GREENNESS',gvfid)
    call LDT_verify(ios,'GREENNESS field not found in the LIS param file')
    
    ios = nf90_get_var(ftn_ldt,gvfid,gvf)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_ldt,gvfid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,gvf1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(albclimo(nc,nr,nmonth))
!    call LDT_verify(nf90_def_dim(ftn_real,'month',nmonth,dimID(3)))
    call LDT_verify(nf90_def_var(ftn_real,'Albedo',&
         nf90_float, dimids=dimID(1:3),varid=calb1Id),&
         'nf90_def_var failed for Albedo')

    ios = nf90_inq_varid(ftn_ldt,'ALBEDO',calbid)
    call LDT_verify(ios,'ALBEDO field not found in the LIS param file')
    
    ios = nf90_get_var(ftn_ldt,calbid,albclimo)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_ldt,calbid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,calb1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(mxsnalb(nc,nr))
    ios = nf90_inq_varid(ftn_ldt,'MXSNALBEDO',mxalbid)
    call LDT_verify(ios,'MXSNALBEDO field not found in the LIS param file')
    
    ios = nf90_get_var(ftn_ldt,mxalbid,mxsnalb)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_ldt,mxalbid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,mxsnalb1Id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

! read from the LIS history file, 
    
    write(LDT_logunit,*) 'Reading LIS history file '//trim(LDT_NUWRF_struc(n)%LIShistfile)
    call LDT_verify(nf90_open(path=LDT_NUWRF_struc(n)%LIShistfile, &
         mode=NF90_NOWRITE, ncid = ftn_lis),&
         'Error opening file '//trim(LDT_NUWRF_struc(n)%LIShistfile))
    
    allocate(sca(nc,nr))
    ios = nf90_inq_varid(ftn_lis,'Snowcover_inst',scaid)
    call LDT_verify(ios,'Snowcover_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,scaid,sca)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,scaid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,sca1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(snod(nc,nr))
    ios = nf90_inq_varid(ftn_lis,'SnowDepth_inst',snodid)
    call LDT_verify(ios,'SnowDepth_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,snodid,snod)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,snodid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,snod1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)


    allocate(swe(nc,nr))
    ios = nf90_inq_varid(ftn_lis,'SWE_inst',sweid)
    call LDT_verify(ios,'SWE_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,sweid,swe)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,sweid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,swe1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    
    allocate(canopint(nc,nr))
    ios = nf90_inq_varid(ftn_lis,'CanopInt_inst',canopintid)
    call LDT_verify(ios,'CanopInt_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,canopintid,canopint)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,canopintid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,canopint1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(avgsurft(nc,nr))
    ios = nf90_inq_varid(ftn_lis,'AvgSurfT_tavg',avgsurftid)
    call LDT_verify(ios,'AvgSurfT_tavg field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,avgsurftid,avgsurft)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,avgsurftid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,avgsurft1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(alb(nc,nr))
    ios = nf90_inq_varid(ftn_lis,'Albedo_inst',albid)
    call LDT_verify(ios,'Albedo_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,albid,alb)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,albid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,alb1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(tbot(nc,nr))
    ios = nf90_inq_varid(ftn_lis,'Tempbot_inst',tbotid)
    call LDT_verify(ios,'Tempbot_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,tbotid,tbot)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,tbotid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,tbot1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    allocate(gvf_inst(nc,nr))
    ios = nf90_inq_varid(ftn_lis,'Greenness_inst',gvfinstid)
    call LDT_verify(ios,'Greenness_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,gvfinstid,gvf_inst)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,gvfinstid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,gvfinst1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    ios = nf90_inq_dimid(ftn_lis,"SoilMoist_profiles",smpId)
    call LDT_verify(ios,'Error in nf90_inq_dimid in LDT_procDataForREAL')

    ios = nf90_inquire_dimension(ftn_lis,smpId, len=nsmp)
    call LDT_verify(ios,'Error in nf90_inquire_dimension in LDT_procDataForREAL')

    call LDT_verify(nf90_def_dim(ftn_real,'SoilMoist_profiles',nsmp,dimID(3)))
    call LDT_verify(nf90_def_var(ftn_real,'SoilMoist',&
         nf90_float, dimids=dimID(1:3),varid=smc1Id),&
         'nf90_def_var failed for SoilMoist')

    allocate(smc(nc,nr,nsmp))
    ios = nf90_inq_varid(ftn_lis,'SoilMoist_inst',smcid)
    call LDT_verify(ios,'SoilMoist_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,smcid,smc)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,smcid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,smc1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)


    ios = nf90_inq_dimid(ftn_lis,"SoilTemp_profiles",stpId)
    call LDT_verify(ios,'Error in nf90_inq_dimid in LDT_procDataForREAL')

    ios = nf90_inquire_dimension(ftn_lis,stpId, len=nstp)
    call LDT_verify(ios,'Error in nf90_inquire_dimension in LDT_procDataForREAL')

    call LDT_verify(nf90_def_dim(ftn_real,'SoilTemp_profiles',nsmp,dimID(3)))
    call LDT_verify(nf90_def_var(ftn_real,'SoilTemp',&
         nf90_float, dimids=dimID(1:3),varid=stc1Id),&
         'nf90_def_var failed for SoilMoist')

    allocate(stc(nc,nr,nstp))
    ios = nf90_inq_varid(ftn_lis,'SoilTemp_inst',stcid)
    call LDT_verify(ios,'SoilTemp_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,stcid,stc)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,stcid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,stc1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)


    ios = nf90_inq_dimid(ftn_lis,"SmLiqFrac_profiles",slpId)
    call LDT_verify(ios,'Error in nf90_inq_dimid in LDT_procDataForREAL')

    ios = nf90_inquire_dimension(ftn_lis,slpId, len=nslp)
    call LDT_verify(ios,'Error in nf90_inquire_dimension in LDT_procDataForREAL')

    call LDT_verify(nf90_def_dim(ftn_real,'SoilWet_profiles',nslp,dimID(3)))
    call LDT_verify(nf90_def_var(ftn_real,'SoilWet',&
         nf90_float, dimids=dimID(1:3),varid=sh2o1Id),&
         'nf90_def_var failed for SoilWet')

    allocate(sh2o(nc,nr,nslp))
    ios = nf90_inq_varid(ftn_lis,'SmLiqFrac_inst',sh2oid)
    call LDT_verify(ios,'SmLiqFrac_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,sh2oid,sh2o)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,sh2oid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,sh2o1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)

    ios = nf90_inq_dimid(ftn_lis,"RelSMC_profiles",relsmcId)
    call LDT_verify(ios,'Error in nf90_inq_dimid in LDT_procDataForREAL')

    ios = nf90_inquire_dimension(ftn_lis,relsmcId, len=nrsmc)
    call LDT_verify(ios,'Error in nf90_inquire_dimension in LDT_procDataForREAL')
    call LDT_verify(nf90_def_dim(ftn_real,'RelSMC_profiles',nrsmc,dimID(3)))
    call LDT_verify(nf90_def_var(ftn_real,'RelSMC',&
         nf90_float, dimids=dimID(1:3),varid=rsmc1Id),&
         'nf90_def_var failed for RelSMC')

    allocate(relsmc(nc,nr,nrsmc))
    ios = nf90_inq_varid(ftn_lis,'RelSMC_inst',rsmcid)
    call LDT_verify(ios,'RelSMC_inst field not found in LIS history file')
    
    ios = nf90_get_var(ftn_lis,rsmcid,relsmc)
    call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')

    call getVarAttribs(ftn_lis,rsmcid,units,standard_name,scale_factor,&
         offset,vmin,vmax)
    call setVarAttribs(ftn_real,rsmc1id,units,standard_name,scale_factor,&
         offset,vmin,vmax)


    ios = nf90_get_att(ftn_ldt,NF90_GLOBAL,"LANDCOVER_SCHEME",lcscheme)
    call LDT_verify(ios,'Error in nf90_get_att:LANDCOVER_SCHEME in LDT_procDataForREAL')

    ios = nf90_put_att(ftn_real,NF90_GLOBAL,"landcover_scheme",lcscheme)
    call LDT_verify(ios,'Error in nf90_put_att in LDT_procDataForREAL')

    ios = nf90_get_att(ftn_lis,NF90_GLOBAL,"NUM_SOIL_LAYERS",nsoillayers)
    call LDT_verify(ios,'Error in nf90_get_att:NUM_SOIL_LAYERS in LDT_procDataForREAL')

    ios = nf90_put_att(ftn_real,NF90_GLOBAL,"num_soil_layers",nsoillayers)
    call LDT_verify(ios,'Error in nf90_put_att in LDT_procDataForREAL')

    allocate(lyrthk(nsoillayers))

    ios = nf90_get_att(ftn_lis,NF90_GLOBAL,"SOIL_LAYER_THICKNESSES",lyrthk)
    call LDT_verify(ios,'Error in nf90_get_att:SOIL_LAYER_THICKNESSES in LDT_procDataForREAL')

    ios = nf90_put_att(ftn_real,NF90_GLOBAL,"soil_layer_thicknesses",lyrthk)
    call LDT_verify(ios,'Error in nf90_put_att in LDT_procDataForREAL')

    call LDT_verify(nf90_put_att(ftn_real,NF90_GLOBAL,"missing_value", -9999.0))          

    deallocate(lyrthk)

!write to realfile    
    call LDT_verify(nf90_enddef(ftn_real))

    call LDT_verify(nf90_put_var(ftn_real,mask1Id,mask,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for landmask ')

    call LDT_verify(nf90_put_var(ftn_real,lc1Id,ltype,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for landcover type ')

    call LDT_verify(nf90_put_var(ftn_real,lu1Id,lc,(/1,1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n),nsfctypes/)),&
         'nf90_put_var failed for landcover fraction ')

    call LDT_verify(nf90_put_var(ftn_real,slt1Id,stype,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for soil type ')

    call LDT_verify(nf90_put_var(ftn_real,elev1Id,elev,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for elevation ')

    call LDT_verify(nf90_put_var(ftn_real,shdmax1Id,shdmax,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for shdmax ')

    call LDT_verify(nf90_put_var(ftn_real,shdmin1Id,shdmin,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for shdmin ')

    call LDT_verify(nf90_put_var(ftn_real,sca1Id,sca,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for sca ')

    call LDT_verify(nf90_put_var(ftn_real,snod1Id,snod,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for snod ')

    call LDT_verify(nf90_put_var(ftn_real,swe1Id,swe,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for swe ')

    call LDT_verify(nf90_put_var(ftn_real,canopint1Id,canopint,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for canopint ')

    call LDT_verify(nf90_put_var(ftn_real,avgsurft1Id,avgsurft,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for avgsurft ')

    call LDT_verify(nf90_put_var(ftn_real,alb1Id,alb,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for albedo_inst ')

    call LDT_verify(nf90_put_var(ftn_real,calb1Id,albclimo,(/1,1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n),nmonth/)),&
         'nf90_put_var failed for alb ')

    call LDT_verify(nf90_put_var(ftn_real,tbot1Id,tbot,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for Tempbot ')

    call LDT_verify(nf90_put_var(ftn_real,mxsnalb1Id,mxsnalb,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for maxsnoalbedo')

    call LDT_verify(nf90_put_var(ftn_real,gvfinst1Id,gvf_inst,(/1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_var failed for greeness_inst')

    call LDT_verify(nf90_put_var(ftn_real,gvf1Id,gvf,(/1,1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n),nmonth/)),&
         'nf90_put_var failed for greeness')

    call LDT_verify(nf90_put_var(ftn_real,smc1Id,smc,(/1,1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n),nsmp/)),&
         'nf90_put_var failed for smc ')

    call LDT_verify(nf90_put_var(ftn_real,sh2o1Id,sh2o,(/1,1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n),nslp/)),&
         'nf90_put_var failed for sh2o ')

    call LDT_verify(nf90_put_var(ftn_real,stc1Id,stc,(/1,1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n),nstp/)),&
         'nf90_put_var failed for stc ')

    call LDT_verify(nf90_put_var(ftn_real,rsmc1Id,relsmc,(/1,1,1/),&
         (/LDT_rc%gnc(n),LDT_rc%gnr(n),nrsmc/)),&
         'nf90_put_var failed for relsmc ')

    call LDT_verify(nf90_close(ftn_real))
    call LDT_verify(nf90_close(ftn_lis))
    call LDT_verify(nf90_close(ftn_ldt))

    write(LDT_logunit,*) 'Completed the processsing for real data initialization'
    deallocate(lc)
    deallocate(texture)
    deallocate(mask)
    deallocate(ltype)
    deallocate(stype)
    deallocate(elev)
    deallocate(shdmin)
    deallocate(shdmax)
    deallocate(sca)
    deallocate(snod)
    deallocate(swe)
    deallocate(avgsurft)
    deallocate(canopint)
    deallocate(alb)
    deallocate(albclimo)
    deallocate(mxsnalb)
    deallocate(gvf_inst)
    deallocate(tbot)
    deallocate(smc)
    deallocate(stc)
    deallocate(sh2o)
    deallocate(gvf)
    deallocate(relsmc)
    
#endif    
  end subroutine LDT_procDataForREAL

!BOP
! 
! !ROUTINE: getVarAttribs
! \label{getVarAttribs}
!
! !INTERFACE: 
  subroutine  getVarAttribs(ftn,varid,units,standard_name,scale_factor,&
       offset,vmin,vmax)
    
    implicit none
! !ARGUMENTS: 
    integer               :: ftn
    integer               :: varid
    character(len=*)      :: units
    character(len=*)      :: standard_name
    real                  :: scale_factor
    real                  :: offset
    real                  :: vmin
    real                  :: vmax
!
! !DESCRIPTION: 
! 
!  This subroutine extracts the attributes (units, standard name, scale
!  factor, offset, min and max values) associated with a certain variable. 
!EOP
#if(defined USE_NETCDF3 || defined USE_NETCDF4)        
    call LDT_verify(nf90_get_att(ftn,varid,"units",units),&
         'nf90_get_att failed in getVarAttribs:LDT_NUWRFprocMod')
    call LDT_verify(nf90_get_att(ftn,varid,"standard_name",standard_name),&
         'nf90_get_att failed in getVarAttribs:LDT_NUWRFprocMod')
    call LDT_verify(nf90_get_att(ftn,varid,"scale_factor",scale_factor),&
         'nf90_get_att failed in getVarAttribs:LDT_NUWRFprocMod')
    call LDT_verify(nf90_get_att(ftn,varid,"add_offset",offset),&
         'nf90_get_att failed in getVarAttribs:LDT_NUWRFprocMod')
    call LDT_verify(nf90_get_att(ftn,varid,"vmin",vmin),&
         'nf90_get_att failed in getVarAttribs:LDT_NUWRFprocMod')
    call LDT_verify(nf90_get_att(ftn,varid,"vmax",vmax),&
         'nf90_get_att failed in getVarAttribs:LDT_NUWRFprocMod')
#endif
    
  end subroutine getVarAttribs

!BOP
! 
! !ROUTINE: setVarAttribs
! \label{setVarAttribs}
!
! !INTERFACE: 
  subroutine  setVarAttribs(ftn,varid,units,standard_name,scale_factor,&
       offset,vmin,vmax)
    
    implicit none
! !ARGUMENTS: 
    integer               :: ftn
    integer               :: varid
    character(len=*)      :: units
    character(len=*)      :: standard_name
    real                  :: scale_factor
    real                  :: offset
    real                  :: vmin
    real                  :: vmax
!
! !DESCRIPTION: 
! 
!  This subroutine sets the attributes (units, standard name, scale
!  factor, offset, min and max values) associated with a certain variable. 
!EOP
#if(defined USE_NETCDF3 || defined USE_NETCDF4)        
    call LDT_verify(nf90_put_att(ftn,varid,"units",units),&
         'nf90_put_att units failed in setVarAttribs:LDT_NUWRFprocMod')
    call LDT_verify(nf90_put_att(ftn,varid,"standard_name",standard_name),&
         'nf90_put_att standard_name failed in setVarAttribs:LDT_NUWRFprocMod')
    call LDT_verify(nf90_put_att(ftn,varid,"scale_factor",scale_factor),&
         'nf90_put_att scale_factor failed in setVarAttribs:LDT_NUWRFprocMod')
    call LDT_verify(nf90_put_att(ftn,varid,"add_offset",offset),&
         'nf90_put_att add_offset failed in setVarAttribs:LDT_NUWRFprocMod')
    call LDT_verify(nf90_put_att(ftn,varid,"vmin",vmin),&
         'nf90_put_att vmin failed in setVarAttribs:LDT_NUWRFprocMod')
    call LDT_verify(nf90_put_att(ftn,varid,"vmax",vmax),&
         'nf90_put_att vmax failed in setVarAttribs:LDT_NUWRFprocMod')
#endif
    
  end subroutine setVarAttribs

end module LDT_NUWRFprocMod

