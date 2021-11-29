!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module alltypes_irrigationMod
!BOP
!
! !MODULE: alltypes_irrigationMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  11 Nov 2012: Sujay Kumar; Initial implementation
!  25 Feb 2020: Jessica Erlingis; Update irrigation window
!  20 Nov 2020: Hiroko Beaudoing; Consolidated common routines for sprinkler,
!               drip, and flood schemes and incorporated irrigation
!               type map and calendar information as well as a routine for 
!               mapping croptypes to specific irrigation type (developped by 
!               Matt). Updated IrrigScale determination.
!  29 Oct 2021: Sarith Mahanama; Added mapping croptypes to irrigation types.
!
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_irrigationMod
  
  implicit none

  PRIVATE
  
  PUBLIC  :: alltypes_irrigation_init
  PUBLIC  :: alltypes_irrigation_updates

contains
  
  subroutine alltypes_irrigation_init(irrigState)

    type(ESMF_State) :: irrigState(LIS_rc%nnest)

    integer              :: n 
    integer              :: rc, status
    integer              :: i
    type(ESMF_ArraySpec) :: arrspec1
    type(ESMF_Field)     :: irrigRateField, irrigFracField
    type(ESMF_Field)     :: irrigRootDepthField, irrigScaleField
    type(ESMF_Field)     :: irrigTypeField
    real,  allocatable   :: irrigFrac(:)
    real,  allocatable   :: irrigRootdepth(:)
    real,  allocatable   :: irrigScale(:)
    real,  allocatable   :: irrigType(:)
    real, pointer        :: frac(:),scale(:),rootdepth(:),irrigRate(:)
    real, pointer        :: itype(:)
    character*100        :: maxrootdepthfile
    character*50         :: irrigtypetocrop
    integer              :: nlctypes ! non-crop land cover types
! __________________________________________________________


    do n=1,LIS_rc%nnest
       allocate(irrigFrac(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(irrigRootDepth(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(irrigScale(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(irrigType(LIS_rc%npatch(n,LIS_rc%lsm_index)))   ! HKB
 
       write(LIS_logunit,*) "[INFO] Setting up irrigation method ... "

     ! Read irrigation fraction (or "intensity") input:
       call read_irrigFrac(n, irrigFrac)

     ! Read crop type maximum root depth file (combined landcover/crop classifications):
!HKB the file is the same for Sprinkler, Drip, and Flood
       call ESMF_ConfigGetAttribute(LIS_config,maxrootdepthfile,&
            label="Irrigation max root depth file:",rc=rc)
       call LIS_verify(rc,&
            'Irrigation max root depth file: option not specified in the config file')

       call read_irrigRootdepth( n, maxrootdepthfile, irrigRootDepth, nlctypes )

!HKB ! Read irrigation type frac and assign a type for the tile by opted method:
       call ESMF_ConfigGetAttribute(LIS_config,irrigtypetocrop,&
            label="Irrigation type to crop mapping method:",rc=rc)
       call LIS_verify(rc,&
            'Irrigation type to crop mapping method: option not specified in the config file')
       call get_irrigType(n, irrigtypetocrop, nlctypes, irrigType)

!HKB ! Read plant/harvest dates if opted:
       if ( LIS_irrig_struc(n)%cropcalendar .ne. "none" ) then
         call read_cropcalendar(n, LIS_irrig_struc(n)%cropseasons, nlctypes)
       endif 

     ! Calculate irrigation scale:
       call compute_irrigScale(n, irrigFrac, irrigScale)

       call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status, &
            "ESMF_ArraySpecSet failed in alltypes_irrigation_init")

       irrigRateField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation rate", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in alltypes_irrigation_init")

       call ESMF_FieldGet(irrigRateField,localDE=0,&
            farrayPtr=irrigRate,rc=status)
       call LIS_verify(status,'ESMF_FieldGet failed for irrigrate ')
       irrigRate = 0.0

       call ESMF_StateAdd(irrigState(n),(/irrigRateField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for irrigRate failed in alltypes_irrigation_init")

       irrigFracField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation frac", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in alltypes_irrigation_init")

       call ESMF_FieldGet(irrigFracField,localDE=0,&
            farrayPtr=frac,rc=status)
       call LIS_verify(status,'ESMF_FieldGet failed for IrrigFrac')
       
       frac = irrigFrac

       call ESMF_StateAdd(irrigState(n),(/irrigFracField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for irrigFrac failed in alltypes_irrigation_init")
       deallocate(irrigFrac)

       irrigRootdepthField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation max root depth", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in alltypes_irrigation_init")

       call ESMF_FieldGet(irrigRootdepthField,localDE=0,&
            farrayPtr=rootdepth,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet failed for root depth')
       rootdepth=irrigRootDepth
       
       call ESMF_StateAdd(irrigState(n),(/irrigRootdepthField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for max root depth failed in alltypes_irrigation_init")
       deallocate(irrigRootDepth)

       irrigScaleField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation scale",rc=status)
       call LIS_verify(status, &
            'ESMF_FieldCreate failed in alltypes_irrigation_init')
       
       call ESMF_FieldGet(irrigScaleField,localDE=0,&
            farrayPtr = scale,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet failed for irrigation scale')

       scale = irrigScale
       call ESMF_StateAdd(irrigState(n),(/irrigScaleField/),rc=status)
       call LIS_verify(status,&
            'ESMF_StateAdd for irrigation scale failed in alltypes_irrigation_init')
       deallocate(irrigScale)

!HKB added below
       irrigTypeField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation type",rc=status)
       call LIS_verify(status, &
            'ESMF_FieldCreate failed in alltypes_irrigation_init')
       
       call ESMF_FieldGet(irrigTypeField,localDE=0,&
            farrayPtr = itype ,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet failed for irrigation type')

       itype = irrigType
       call ESMF_StateAdd(irrigState(n),(/irrigTypeField/),rc=status)
       call LIS_verify(status,&
            'ESMF_StateAdd for irrigation type failed in alltypes_irrigation_init')
       deallocate(irrigType)

    enddo

  end subroutine alltypes_irrigation_init
  

  subroutine alltypes_irrigation_updates(n, irrigState)

    use LIS_FORC_AttributesMod 
    use LIS_histDataMod
    use LIS_metforcingMod, only : LIS_FORC_State    

    implicit none

    integer, intent(in) :: n 
    type(ESMF_State)    :: irrigState
    
! moved to lis.config
!    real, parameter     :: otimes = 6.0 ! local trigger check start time [hour]
!    real, parameter     :: otimee = 10.0 !local trigger check end time [hour]
    integer             :: t,gid
    real                :: ltime
    integer             :: chhr, lhr
    integer             :: status
    integer             :: ftn

    type(ESMF_Field)    :: irrigRateField,prcpField
    type(ESMF_Field)    :: irrigTypeField
    real,    pointer    :: prcp(:)
    real,    pointer    :: irrigRate(:)
    real,    pointer    :: irrigType(:)
    real,  allocatable  :: irrigAmt(:)    ! local variable

    integer             :: tmpval

    real                :: timestep, shift_otimess, shift_otimese
    real                :: shift_otimeds, shift_otimede
    real                :: shift_otimefs, shift_otimefe

    allocate(irrigAmt(LIS_rc%npatch(n,LIS_rc%lsm_index)))
    call ESMF_StateGet(irrigState,&
         "Irrigation rate",&
         irrigRateField,rc=status)
    call LIS_verify(status,&
         'ESMF_StateGet failed for Irrigation rate')
    
    call ESMF_FieldGet(irrigRateField,localDE=0,&
         farrayPtr=irrigrate,rc=status)
    call LIS_verify(status,'ESMF_FieldGet failed for irrigrate ')

    call ESMF_StateGet(irrigState,&
         "Irrigation type",&
         irrigTypeField,rc=status)
    call LIS_verify(status,&
         'ESMF_StateGet failed for Irrigation type')
    
    call ESMF_FieldGet(irrigTypeField,localDE=0,&
         farrayPtr=irrigtype,rc=status)
    call LIS_verify(status,'ESMF_FieldGet failed for irrigtype ')

    call ESMF_StateGet(LIS_FORC_State(n),&
         trim(LIS_FORC_Rainf%varname(1)),prcpField,&
         rc=status)
    call LIS_verify(status,&
         'ESMF_StateGet failed for rainf in alltypes_irrigation')

    
    call ESMF_FieldGet(prcpField,localDE=0, farrayPtr=prcp,rc=status)
    call LIS_verify(status,&
         'ESMF_FieldGet failed for rainf in alltypes_irrigation')

    irrigAmt = 0.0

    timestep = LIS_rc%ts

    ! Adjust bounds by timestep to account for the fact that LIS_rc%hr, etc.
    ! will represents the END of the integration timestep window
    ! Sprinkler
    shift_otimess = LIS_irrig_struc(n)%sprinkler_start + (timestep/3600.)
    shift_otimese = (LIS_irrig_struc(n)%sprinkler_start + &
                     LIS_irrig_struc(n)%sprinkler_duration) + (timestep/3600.)
    ! Drip
    shift_otimeds = LIS_irrig_struc(n)%drip_start + (timestep/3600.)
    shift_otimede = (LIS_irrig_struc(n)%drip_start + &
                     LIS_irrig_struc(n)%drip_duration) + (timestep/3600.)
    ! Flood
    shift_otimefs = LIS_irrig_struc(n)%flood_start + (timestep/3600.)
    shift_otimefe = (LIS_irrig_struc(n)%flood_start + &
                     LIS_irrig_struc(n)%flood_duration) + (timestep/3600.)


    ! HKB check 
    ftn = LIS_getNextUnitNumber()
    open(ftn,file="tileirrtype2.txt",status='unknown',form='formatted')

    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

       gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
       chhr = nint(24.0*(LIS_domain(n)%grid(gid)%lon/360.0))
       if((LIS_domain(n)%grid(gid)%lon.lt.0.0).and.&
            (abs(mod(LIS_domain(n)%grid(gid)%lon,15.0)).ge.0.0001)) &
          chhr = chhr -1
       lhr = LIS_rc%hr+chhr
       if(lhr.ge.24) lhr = lhr-24
       if(lhr.lt.0) lhr = lhr+24
       
       ltime = real(lhr)+real(LIS_rc%mn)/60.0+real(LIS_rc%ss)/3600.0
!HKB apply water during specified hours for the irrigation type
       if ( irrigType(t) .eq. 1 ) then  ! sprinkler, add to rain
         if((ltime.ge.shift_otimess).and.(ltime.lt.shift_otimese)) then 
          prcp(t)=prcp(t)+irrigRate(t)
          irrigAmt(t) = irrigRate(t)
         end if
       elseif ( irrigType(t) .eq. 2 ) then  ! drip
         if((ltime.ge.shift_otimeds).and.(ltime.lt.shift_otimede)) then 
          irrigAmt(t) = irrigRate(t)
         end if
       elseif ( irrigType(t) .eq. 3 ) then  ! flood
         if((ltime.ge.shift_otimefs).and.(ltime.lt.shift_otimefe)) then 
          irrigAmt(t) = irrigRate(t)
         end if
       endif
    
     ! Write out irrigated amount of water to separate directory:
       if( LIS_MOC_IRRIGATEDWATER == 1 ) then
         call LIS_diagnoseIrrigationOutputVar(n,t,LIS_MOC_IRRIGATEDWATER,&
               value=irrigAmt(t),unit="kg m-2 s-1",direction="-",vlevel=1)
       else
         write(LIS_logunit,*) "[ERR] 'Irrigated water:' was not selected. Program stopping ..."
         call LIS_endrun
       endif
     ! HKB check 
       write(ftn,fmt='(i8,2f10.4,f3.0)') &
       t,LIS_domain(n)%grid(gid)%lon,LIS_domain(n)%grid(gid)%lat,irrigType(t)

    enddo
    call LIS_releaseUnitNumber(ftn)  !HKB
    deallocate(irrigAmt)
    
  end subroutine alltypes_irrigation_updates


  subroutine read_irrigFrac(n,frac)

    use LIS_fileIOMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer,  intent(in) :: n 
    real                 :: frac(LIS_rc%npatch(n,LIS_rc%lsm_index))

    integer              :: t,col,row
    integer              :: nid,ios,status,fracId
    logical              :: file_exists    
    real,  allocatable   :: l_frac(:,:)
    real,  allocatable   :: glb_frac(:,:)
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 

       allocate(l_frac(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       
       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in the lis input netcdf file')

       write(LIS_logunit,*) "[INFO] Reading in the irrigation fraction field ... "
       
       allocate(glb_frac(LIS_rc%gnc(n),LIS_rc%gnr(n)))

       ios = nf90_inq_varid(nid,'IRRIGFRAC',fracId)
       call LIS_verify(ios,'nf90_inq_varid failed for IRRIGFRAC')
       
       ios = nf90_get_var(nid,fracId, glb_frac)
       call LIS_verify(ios,'nf90_get_var failed for in alltypes_irrigationMod')

       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in alltypes_irrigationMod')
       
       l_frac(:,:) = glb_frac(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1):&
            LIS_nse_halo_ind(n,LIS_localPet+1))
       deallocate(glb_frac)

       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
          row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
          frac(t) = l_frac(col,row)
       enddo

       deallocate(l_frac)
    else
       write(LIS_logunit,*) "[ERR] Irrigation frac map: ",&
             LIS_rc%paramfile(n),"[ERR] does not exist."
       write(LIS_logunit,*) "[ERR] Program stopping ..."
       call LIS_endrun
    endif
#endif
  end subroutine read_irrigFrac


  subroutine read_irrigRootdepth(n, rdfile, rootdepth, nlctypes)

!  25 Nov 2020: Hiroko Beaudoing; updated CROPTYPE array content 
!               and structure changes.
!
    use LIS_fileIOMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer,    intent(in) :: n 
    character(len=*)       :: rdfile
    real                   :: rootdepth(LIS_rc%npatch(n,LIS_rc%lsm_index))
    integer,intent(out)    :: nlctypes
    integer                :: total_vegtypes
    real, allocatable      :: rootd(:)
    integer                :: ftn
    integer                :: t,j,col,row
    integer                :: nid, ios, status, croptypeId, cropdimid
    logical                :: file_exists    
    real,   allocatable    :: l_croptype(:,:,:)
    real,   allocatable    :: glb_croptype(:,:,:)
    integer                :: ncroptypes
    integer                :: vegt
! __________________________________________________________________________
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

 !- Read in LDT input crop classification information (done here for now):
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
     ! Read in LDT-generated netcdf file information:
       write(LIS_logunit,*)"[INFO] Reading crop classification information ..."
       ios = nf90_open(path=LIS_rc%paramfile(n),&
                       mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_irrigRootdepth')
 
       ios = nf90_get_att(nid, NF90_GLOBAL, 'CROPCLASS_SCHEME', LIS_rc%cropscheme)
       call LIS_verify(ios,'Error in nf90_get_att in read_irrigRootdepth')
 
       ios = nf90_get_att(nid, NF90_GLOBAL, 'CROPCLASS_NUMBER', LIS_rc%numbercrops)
       call LIS_verify(ios,'Error in nf90_get_att in read_irrigRootdepth')
       write(LIS_logunit,*)"[INFO] Read in crop classfication: ",trim(LIS_rc%cropscheme),&
                          ", with the number of crop types:",LIS_rc%numbercrops
       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in read_irrigRootdepth')
    endif

  ! Estimate total crop types, added to landcover scheme class total:
    select case ( LIS_rc%lcscheme )
     case( "UMD", "UMD+CROPMAP" )
       total_vegtypes = 13 + LIS_rc%numbercrops
     case( "UMD+MIRCA" )
       total_vegtypes = 14 + LIS_rc%numbercrops
     !case( "IGBP", "IGBPNCEP")
     case( "IGBP", "IGBPNCEP", "IGBP+MIRCA", "IGBPNCEP+MIRCA" )
       total_vegtypes = 20 + LIS_rc%numbercrops
     case( "USGS" )
       total_vegtypes = 24 + LIS_rc%numbercrops
     case default
       write(LIS_logunit,*) "[ERR] The landcover scheme, ",trim(LIS_rc%lcscheme),","
       write(LIS_logunit,*) "[ERR] is not supported for irrigation. "
       write(LIS_logunit,*) " Stopping program ... "
       call LIS_endrun()
    end select

    nlctypes = total_vegtypes - LIS_rc%numbercrops ! number of land cover types

 !- Read the max root depth table file:
    inquire(file=rdfile,exist=file_exists)
    if(file_exists) then 
       write(LIS_logunit,*) "[INFO] Reading in the max root depth file: ",trim(rdfile)
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=rdfile,status='old')
       allocate( rootd(total_vegtypes) )
       read(ftn,*) (rootd(j),j=1,total_vegtypes)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) "[ERR] Max root depth file, ",trim(rdfile),", not found."
       write(LIS_logunit,*) "[ERR] Stopping program ..."
       call LIS_endrun()
    endif

 !- Read in crop type map file (specified in LIS parameter input file)
 !- CROPTYPE is now 3D array (croptypes,lat,lon) if old parameter file
 !  with 2D CROPTYPE or no croptypes dimension, need to rerun LDT -HKB
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       ios = nf90_open(path=LIS_rc%paramfile(n),&
                       mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_irrigRootdepth')

       write(LIS_logunit,*) "[INFO] Reading in the crop type field ... "
       
       ! LIS_rc%numbercrops = ncroptypes dimension normally, unless  
       ! single or constant croptype is assinged, so get ncroptypes info
       ios = nf90_inq_dimid(nid, "croptypes", cropdimid)
       call LIS_verify(ios,'nf90_inq_dimid failed for CROPTYPE, NEED new lis_input')
       ios = nf90_inquire_dimension(nid, cropdimid, len = ncroptypes)
       call LIS_verify(ios,'nf90_inquire_dimension failed for CROPTYPE')

       allocate(l_croptype(LIS_rc%lnc(n),LIS_rc%lnr(n),ncroptypes))
       allocate(glb_croptype(LIS_rc%gnc(n),LIS_rc%gnr(n),ncroptypes))
       
       ios = nf90_inq_varid(nid,'CROPTYPE',croptypeId)
       call LIS_verify(ios,'nf90_inq_varid failed for CROPTYPE')
       
       ios = nf90_get_var(nid, croptypeId, glb_croptype)
       call LIS_verify(ios,'nf90_get_var failed for CROPTYPE')

       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in read_irrigRootdepth')
       
       l_croptype(:,:,:) = glb_croptype(&
          LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1),:)

       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
          row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
          vegt = LIS_domain(n)%tile(t)%vegt  !veg index 1~total_vegtype(eg 1-46)
          ! root depth is set to zero for forest, water, bare soil etc, per
          ! land cover class, so assign all tiles 
!          if ( vegt .gt. nlctypes ) then     !crops
!            j = vegt - nlctypes  !crop index nlctypes~total_vegtype(eg 21-46)
!            if(l_croptype(col,row,j).gt.0) then 
!             !rootdepth(t) = rootd(nint(l_croptype(col,row))) 
             rootdepth(t) = rootd(vegt)
!            else
!             rootdepth(t) = 0 
!            endif
!          else
!           rootdepth(t) = 0 
!          endif
       enddo
       deallocate( rootd )
       deallocate(glb_croptype)
       deallocate(l_croptype)

    else
       write(LIS_logunit,*) "[ERR] The irrigation croptype map: ",&
             LIS_rc%paramfile(n)," does not exist."
       write(LIS_logunit,*) "[ERR] Program stopping ..."
       call LIS_endrun
    endif
#endif

  end subroutine read_irrigRootdepth

  subroutine read_cropcalendar(n, cropseasons, nlctypes)
! Reading PLANT and HARVEST dates if applicable. 

    use LIS_fileIOMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer,  intent(in) :: n 
    integer,  intent(in) :: cropseasons
    integer,  intent(in) :: nlctypes      !land cover index for veg 
    integer              :: t,col,row,j,i
    integer              :: nid,ios,status,pdayId,hdayId
    logical              :: file_exists    
    real,  allocatable   :: l_frac_p(:,:,:)
    real,  allocatable   :: l_frac_h(:,:,:)
    real,  allocatable   :: glb_frac_p(:,:,:,:)
    real,  allocatable   :: glb_frac_h(:,:,:,:)
    integer              :: vegt
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 

       allocate(l_frac_p(LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%numbercrops))
       allocate(l_frac_h(LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%numbercrops))
       
       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in the lis input netcdf file')

       write(LIS_logunit,*) "[INFO] Reading in the plant/harvest days.. "
       
       allocate(glb_frac_p(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%numbercrops,cropseasons))
       allocate(glb_frac_h(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%numbercrops,cropseasons))

       ios = nf90_inq_varid(nid,'PLANTDAY',pdayId)
       call LIS_verify(ios,'nf90_inq_varid failed for PLANTDAY')
       
       ios = nf90_get_var(nid,pdayId, glb_frac_p)
       call LIS_verify(ios,'nf90_get_var failed for in alltypes_irrigationMod')

       ios = nf90_inq_varid(nid,'HARVESTDAY',hdayId)
       call LIS_verify(ios,'nf90_inq_varid failed for HARVESTDAY')
       
       ios = nf90_get_var(nid,hdayId, glb_frac_h)
       call LIS_verify(ios,'nf90_get_var failed for in alltypes_irrigationMod')
       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in alltypes_irrigationMod')
       
       CROPSEASON : do j=1,cropseasons
           l_frac_p(:,:,:) = glb_frac_p(&
               LIS_ews_halo_ind(n,LIS_localPet+1):&         
               LIS_ewe_halo_ind(n,LIS_localPet+1),&
               LIS_nss_halo_ind(n,LIS_localPet+1):&
               LIS_nse_halo_ind(n,LIS_localPet+1),:,j)
           l_frac_h(:,:,:) = glb_frac_h(&
               LIS_ews_halo_ind(n,LIS_localPet+1):&         
               LIS_ewe_halo_ind(n,LIS_localPet+1),&
               LIS_nss_halo_ind(n,LIS_localPet+1):&
               LIS_nse_halo_ind(n,LIS_localPet+1),:,j)

           CROPTILE : do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
              row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
              vegt = LIS_domain(n)%tile(t)%vegt  !veg index (eg 1-46)
              if ( vegt .gt. nlctypes ) then     !crops
                i = vegt - nlctypes  !crop index (eg 21-46 -> 1-26)
                LIS_irrig_struc(n)%plantDay(t,j) = l_frac_p(col,row,i)
                LIS_irrig_struc(n)%harvestDay(t,j) = l_frac_h(col,row,i)
              else
                LIS_irrig_struc(n)%plantDay(t,j) = LIS_rc%udef
                LIS_irrig_struc(n)%harvestDay(t,j) = LIS_rc%udef
              endif
           enddo CROPTILE
        enddo CROPSEASON

       deallocate(l_frac_p)
       deallocate(glb_frac_p)
       deallocate(l_frac_h)
       deallocate(glb_frac_h)
    else
       write(LIS_logunit,*) "[ERR] Crop Plant/Harvest days map: ",&
             LIS_rc%paramfile(n),"[ERR] does not exist."
       write(LIS_logunit,*) "[ERR] Program stopping ..."
       call LIS_endrun
    endif
#endif
  end subroutine read_cropcalendar

  subroutine get_irrigType(n,irrigtypetocrop,nlctypes,itype)

!HKB new routine for reading in grid level irrigType field and mapping to 
!crop tiles using A) Matt's algorithm, B) dominant, or C) single type
!Added reading in county/country input fields 
! 

    use LIS_fileIOMod
    use LIS_constantsMod,  ONLY : radius => LIS_CONST_REARTH, pi => LIS_CONST_PI
    use getCropIrrigTypes, ONLY : MattA => matt_algorithm
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer,  intent(in) :: n 
    character*50,intent(in) :: irrigtypetocrop
    integer, intent(in)  :: nlctypes ! non-crop land cover types
    real,intent(inout)   :: itype(LIS_rc%npatch(n,LIS_rc%lsm_index))

    integer              :: t,col,row,j
    integer              :: nid,ios,status,itypeId
    integer              :: nirrigtypes, irrigdimid
    integer              :: nsfctypes,   sfcdimid, lcoverid
    integer              :: countyId, countryId
    logical              :: file_exists
    real,  allocatable   :: l_croptype (:,:,:),s_croptype (:,:,:) 
    real,  allocatable   :: g_croptype (:,:,:),t_croptype (:,:,:)
    real,  allocatable   :: l_itype(:,:,:),    s_itype    (:,:,:)
    real,  allocatable   :: glb_itype(:,:,:)
    real,  allocatable   :: l_country(:,:)
    real,  allocatable   :: glb_country(:,:)
    real,  allocatable   :: l_county(:,:)
    real,  allocatable   :: glb_county(:,:)
    real,  allocatable   :: temp(:)
    real,  allocatable   :: PREFTYPE(:,:,:)
    real,  allocatable   :: cell_area (:,:)
    integer              :: iindex
    real                 :: tempval
    integer              :: vegt
    integer              :: ss,ccc
    type(MattA)          :: MA_global, MA_USA
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

 !- IRRIGTYPE is now 3D array (irrigtypes,lat,lon) for fractions of types,
 !  if old parameter file with 2D IRRIGTYPE index map, it will crash with the
 !  error message and need to rerun LDT -HKB
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       
       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in the lis input netcdf file')

       write(LIS_logunit,*) "[INFO] Reading in the irrigation type field ... "
       write(LIS_logunit,*) "[INFO] type determination method is ",irrigtypetocrop
       
       ios = nf90_inq_dimid(nid, "irrigtypes", irrigdimid)
       call LIS_verify(ios,'nf90_inq_dimid failed for IRRIGTYPE, NEED new lis_input')
       ios = nf90_inquire_dimension(nid, irrigdimid, len = nirrigtypes)
       call LIS_verify(ios,'nf90_inquire_dimension failed for IRRIGTYPE')

       allocate(l_itype(LIS_rc%lnc(n),LIS_rc%lnr(n),nirrigtypes))
       allocate(s_itype(LIS_rc%lnc(n),LIS_rc%lnr(n),nirrigtypes))
       allocate(glb_itype(LIS_rc%gnc(n),LIS_rc%gnr(n),nirrigtypes))
       allocate(temp(nirrigtypes))

       ios = nf90_inq_varid(nid,'IRRIGTYPE',itypeId)
       call LIS_verify(ios,'nf90_inq_varid failed for IRRIGTYPE')
       
       ios = nf90_get_var(nid, itypeId, glb_itype)
       call LIS_verify(ios,'nf90_get_var failed for in alltypes_irrigationMod')

       !-- real in LANDCOVER fractions to populate CROPTYPES (lnc,lnc,numbercrops)
       ios = nf90_inq_dimid(nid, "sfctypes", sfcdimid)
       call LIS_verify(ios,'nf90_inq_dimid failed for LANDCOVER, NEED new lis_input')
       ios = nf90_inquire_dimension(nid, sfcdimid, len = nsfctypes)
       call LIS_verify(ios,'nf90_inquire_dimension failed for LANDCOVER')

       allocate(l_croptype (LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%numbercrops))
       allocate(s_croptype (LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%numbercrops))
       allocate(g_croptype (LIS_rc%gnc(n),LIS_rc%gnr(n),nsfctypes))
       allocate(t_croptype (LIS_rc%lnc(n),LIS_rc%lnr(n),nsfctypes))
       ios = nf90_inq_varid(nid,'LANDCOVER',lcoverId)
       call LIS_verify(ios,'nf90_inq_varid failed for LANDCOVER')

       ios = nf90_get_var(nid, lcoverId, g_croptype)
       call LIS_verify(ios,'nf90_get_var failed for in alltypes_irrigationMod')
       
       !-- read in country and county fields       
       allocate(glb_country(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       allocate(glb_county(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       allocate(l_country(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(l_county(LIS_rc%lnc(n),LIS_rc%lnr(n)))

       ios = nf90_inq_varid(nid,'COUNTRY',countryId)
       call LIS_verify(ios,'nf90_inq_varid failed for COUNTRY')
       
       ios = nf90_get_var(nid, countryId, glb_country)
       call LIS_verify(ios,'nf90_get_var failed for in COUNTRY')

       ios = nf90_inq_varid(nid,'COUNTY',countyId)
       call LIS_verify(ios,'nf90_inq_varid failed for COUNTY')
       
       ios = nf90_get_var(nid, countyId, glb_county)
       call LIS_verify(ios,'nf90_get_var failed for in COUNTY')

       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in alltypes_irrigationMod')

       ! Grid to tile mapping CROPTYPES fractions
       t_croptype (:,:,:) = g_croptype (&
          LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1):&
          LIS_nse_halo_ind(n,LIS_localPet+1),:)
       l_croptype = t_croptype (:,:, nsfctypes - LIS_rc%numbercrops + 1: nsfctypes)
       s_croptype = l_croptype
       
       ! Grid to tile mapping index 1=Sprinkler, 2=Drip, 3=Floodg
       ! Note irrigType values are non-missing regardless of croptype or
       ! irrigFrac, based on country/state distribution
       l_itype(:,:,:) = glb_itype(&
          LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1):&
          LIS_nse_halo_ind(n,LIS_localPet+1),:)
       s_itype = l_itype
       
       if (irrigtypetocrop .eq. "distribute") then
 
        allocate(PREFTYPE(LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nsurfacetypes))
        PREFTYPE = -1
        l_county(:,:) = glb_county(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1):&
            LIS_nse_halo_ind(n,LIS_localPet+1))
        l_country(:,:) = glb_country(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1):&
            LIS_nse_halo_ind(n,LIS_localPet+1))
         ss = maxval(l_county)/1000
         ccc = maxval(l_county)-ss*1000
         
         allocate (cell_area (LIS_rc%lnc(n),LIS_rc%lnr(n)))
         call get_area (n,cell_area)

         call MA_global%init_thres (ncrops=LIS_rc%numbercrops,nitypes=nirrigtypes, global = .true.)
         ! process by country irrigtypes and populate PREFTYPE
         call MA_global%git(LIS_rc%gnc(n),LIS_rc%gnr(n), cell_area,l_country,l_croptype, l_itype, PREFTYPE)
         ! if ITYPE_MIN_FRAC/CTYPE_AREA_TOL parameters differ between country and US county
         !    implementation
         call MA_USA%init_thres (ncrops=LIS_rc%numbercrops,nitypes=nirrigtypes)

         ! again with an optional argument. 
         ! rerun with usa option and update PREFTYPE using irrigation fraction vy county
         !     in the US and overwrite PREFTYPE over the US.
         ! initaialize l_croptype, l_itype and PREFTYPE in COUNTY grid cells again
         do t = 1, LIS_rc%nsurfacetypes
            
            where (l_county > 0.)
               PREFTYPE (:,:,t) = -1.               
            endwhere
            if(t <= LIS_rc%numbercrops) then
               where (l_county > 0.)
                  l_croptype (:,:,t) = s_croptype (:,:,t)
               endwhere
            endif
            if (t <= nirrigtypes) then
               where (l_county > 0.)
                  l_itype (:,:,t) = s_itype (:,:,t)
               endwhere
            endif
         end do
         
         call MA_USA%git(LIS_rc%gnc(n),LIS_rc%gnr(n),cell_area, l_county, l_croptype, l_itype, PREFTYPE, usa = .true.)
         !do t = 21, 46
         !   do j=1, LIS_rc%lnr(n)
         !      write(800, '(1440i3)') NINT(preftype(:,j,t))
         !   end do
         !end do
         deallocate (cell_area)
         !stop
       endif

       TILE_LOOP: do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          
          col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
          row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
          vegt= LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt

          select case (irrigtypetocrop)
           case("single")
             ! Simple case: single irrigation type (0 or 1)
             ! it will be overwritten with latter type if more than one type
              do j=1,nirrigtypes
               if ( l_itype(col,row,j).gt.0 ) then
                 itype(t) = j * 1.0
               else
                 itype(t) = 0         ! no irrigation
               endif
              enddo   ! j
           case("dominant")  ! find the max fraction
              temp = l_itype(col,row,:)
              if ( maxval(temp).gt.0 ) then
                 iindex = 0
                 tempval = 0.0
                 do j=1,nirrigtypes
                    if ( temp(j) > tempval ) then
                      tempval = temp(j)
                      iindex = j
                    end if
                 end do
                 itype(t) = iindex * 1.0
              else
                 itype(t) = 0         ! no irrigation
              endif
           case("distribute")

              itype(t) = PREFTYPE(col,row,vegt)
              if ( itype(t).lt.0 .and. itype(t).ne.LIS_rc%udef .and. vegt.gt.nlctypes) then
                 print*,'invalid entry',PREFTYPE(col,row,vegt),col,row,vegt
              endif

           case default
              write(LIS_logunit,*) "[ERR] Irrigation type is not supported for ", &
                trim(irrigtypetocrop)
              call LIS_endrun()
          end select

       enddo TILE_LOOP

       deallocate(l_itype, s_itype)
       deallocate(glb_itype)
       deallocate(temp)
       deallocate(l_country)
       deallocate(glb_country)
       deallocate(l_county)
       deallocate(glb_county)
       deallocate(l_croptype, t_croptype, g_croptype, s_croptype)
       if (allocated(PREFTYPE)) deallocate(PREFTYPE)
    else
       write(LIS_logunit,*) "[ERR] Irrigation type map: ",&
             LIS_rc%paramfile(n),"[ERR] does not exist."
       write(LIS_logunit,*) "[ERR] Program stopping ..."
       call LIS_endrun
    endif
#endif
  end subroutine get_irrigType

  subroutine compute_irrigScale(n, irrigFrac, irrigScale)
! Compute scale to be applied to irrigation amount when the grid total
! crop fraction is less than the irrigation intensity.
! Irrigation is expanded to non-crop, non-forest,
! non-baresoil/urban tiles if intensity exceeds grid total crop fraction.
! In the latter case, scaled irrigation is applied to grassland first,
! then further applied over the rest of tiles equally if the intensity
! exceeds grassland fraction as well.

 !- Inputs/Outputs:
    integer      :: n 
    real         :: irrigFrac(LIS_rc%npatch(n,LIS_rc%lsm_index))
    real         :: irrigScale(LIS_rc%npatch(n,LIS_rc%lsm_index))
  
 !- Local:
    integer      :: t,gid,vegt
    real         :: irrpix,excess
    integer      :: crop1,crop2,grass,shrub1,shrub2
    real, allocatable :: crppix(:)
    real, allocatable :: grasspix(:)
    real, allocatable :: restpix(:)

    allocate(crppix(LIS_rc%ngrid(n)))
    allocate(grasspix(LIS_rc%ngrid(n)))
    allocate(restpix(LIS_rc%ngrid(n)))

    crppix    = 0.0
    grasspix  = 0.0
    restpix   = 0.0

!------------------------------------------------------------------------
! WARNING: The following code is valid only for the no-tiling or the 
! vegetation only tiling. The fgrd values are not valid when multiple
! modes of tiling are turned on. 
! --HKB this WARNIG is no longer true??
!------------------------------------------------------------------------
   select case ( LIS_rc%lcscheme ) 
     case( "UMD" ) 
       crop1  = 11
       crop2  = 11
       grass  = 10 
       shrub1 = 6
       shrub2 = 9
     case( "IGBP", "IGBPNCEP", "MODIS" )
       crop1  = 12
       crop2  = 14
       grass  = 10 
       shrub1 = 6
       shrub2 = 9
     case( "USGS" ) 
       crop1  = 2
       crop2  = 6
       grass  = 7 
       shrub1 = 8
       shrub2 = 10
     case( "UMDCROPMAP" )  
       crop1  = 14   ! Barley
       crop2  = 32   ! Wheat
       grass  = 10
       shrub1 = 6
       shrub2 = 9 
     case( "UMD+MIRCA" )  ! TEMPORARILY ADDED HERE (KRA)
       crop1  = 15   ! Wheat
       crop2  = 40   ! Others annual
       grass  = 10
       shrub1 = 6
       shrub2 = 9
     case( "IGBP+MIRCA", "IGBPNCEP+MIRCA" )
       crop1  = 21   ! Wheat
       crop2  = 46   ! Others annual
       grass  = 10
       shrub1 = 6
       shrub2 = 9
     case default
       write(LIS_logunit,*) "[ERR] The landcover scheme, ",trim(LIS_rc%lcscheme),","
       write(LIS_logunit,*) "[ERR] is not supported for concurrent irrigation."
       write(LIS_logunit,*) " Stopping program ... "
       call LIS_endrun()
   end select
   
   do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
   !- If not water or urban class:
      if(LIS_domain(n)%tile(t)%vegt.ne.LIS_rc%waterclass.and.&
          LIS_domain(n)%tile(t)%vegt.ne.LIS_rc%urbanclass) then 

      !- Crop tiles:
         if(LIS_domain(n)%tile(t)%vegt.ge.crop1.and.&
             LIS_domain(n)%tile(t)%vegt.le.crop2) then 
            crppix(gid) = crppix(gid) + LIS_domain(n)%tile(t)%fgrd*&
                 LIS_domain(n)%tile(t)%pens
      !- Grass tiles:
         elseif(LIS_domain(n)%tile(t)%vegt.eq.grass) then 
            grasspix(gid) = LIS_domain(n)%tile(t)%fgrd*&
                 LIS_domain(n)%tile(t)%pens
      !- Shrub tiles:
         elseif(LIS_domain(n)%tile(t)%vegt.ge.shrub1 .and.&
                LIS_domain(n)%tile(t)%vegt.le.shrub2) then 
            restpix(gid) = restpix(gid) + LIS_domain(n)%tile(t)%fgrd* & 
                 LIS_domain(n)%tile(t)%pens
         endif
      endif
   enddo
   
   irrigScale = 1.0
   do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
      vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
    ! Irrigation pixels:
      irrpix = irrigFrac(t)*0.01
      
      if(vegt.ne.LIS_rc%waterclass) then 
         if(irrpix < crppix(gid)) then 
            if(vegt.ge.crop1.and.vegt.le.crop2) then 
               irrigScale(t) = irrpix/crppix(gid)
            else
               irrigScale(t) = 0.0
            endif
         else
            excess = irrpix -crppix(gid)
            if(excess.gt.grasspix(gid)) then 
               if(vegt.ge.shrub1.and.vegt.le.shrub2) then 
                  irrigScale(t) = (excess - grasspix(gid)) /restpix(gid)
               elseif(vegt.eq.grass) then !grass
                  irrigScale(t) = 1.0
               elseif(vegt.ge.crop1.and.vegt.le.crop2) then !crop
                  irrigScale(t) = 1.0
               else
                  irrigScale(t) = 0.0
               endif
            elseif(excess.lt.grasspix(gid)) then 
               if(vegt.eq.grass) then 
                  irrigScale(t) = excess/grasspix(gid)
               elseif(vegt.ge.crop1.and.vegt.le.crop2) then 
                  irrigScale(t) = 1.0
               else
                  irrigScale(t) = 0.0
               endif
            else
               if(vegt.eq.grass) then 
                  irrigScale(t) =1.0
               elseif(vegt.ge.crop1.and.vegt.le.crop2) then 
                  irrigScale(t) =1.0
               else
                  irrigScale(t) =0.0
               endif
            endif
         endif
      endif
   enddo
#if 0        
       irrigScale =1.0
       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
          vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
          irrpix = irrigFrac(t)*0.01
          if(vegt.ne.LIS_rc%waterclass) then 
             if(irrpix < crppix(gid)) then 
                if(vegt.ge.14) then 
                   irrigScale(t) = irrpix/crppix(gid)
                else
                   irrigScale(t) = 0.0
                endif
             else
                excess = irrpix -crppix(gid)
                if(excess.gt.grasspix(gid)) then 
                   if(vegt.gt.5.and.vegt.lt.10) then 
                      irrigScale(t) = (excess - grasspix(gid)) /restpix(gid)
                   elseif(vegt.eq.10) then !grass
                      irrigScale(t) = 1.0
                   elseif(vegt.ge.14) then !crop
                      irrigScale(t) = 1.0
                   else
                      irrigScale(t) = 0.0
                   endif
                elseif(excess.lt.grasspix(gid)) then 
                   if(vegt.eq.10) then 
                      irrigScale(t) = excess/grasspix(gid)
                   elseif(vegt.ge.14) then 
                      irrigScale(t) = 1.0
                   else
                      irrigScale(t) = 0.0
                   endif
                else
                   if(vegt.eq.10) then 
                      irrigScale(t) =1.0
                   elseif(vegt.ge.14) then 
                      irrigScale(t) =1.0
                   else
                      irrigScale(t) =0.0
                   endif
                endif
             endif
          endif
       enddo
#endif

       deallocate(crppix)
       deallocate(grasspix)
       deallocate(restpix)

     end subroutine compute_irrigScale
     
! ----------------------------------------------------------------

  SUBROUTINE get_area (nest, area)

    use map_utils,        only : ij_to_latlon
    use LIS_constantsMod, ONLY : radius => LIS_CONST_REARTH, pi => LIS_CONST_PI
    
    implicit none

    integer, intent (in)                    :: nest
    real, dimension (:,:),   intent (inout) :: area
    integer                                 :: i,j
    real                                    :: lat_ll, lat_ur , lat_ul, lat_lr, c, r
    real                                    :: lon_ll, lon_ur , lon_ul, lon_lr

    area    = 0.

    do j = 1, LIS_rc%lnr(nest)
       do i = 1, LIS_rc%lnc(nest)
           
          r = float (j)
          c = float (i)

          select case (LIS_domain(nest)%lisproj%code)
          case (0)
             ! lat/lon
             call ij_to_latlon(LIS_domain(nest)%lisproj,c, r, lat_ll, lon_ll)
             area (i,j) = area_latlon (nest, lat_ll)
             
          case (3)
             ! Lambert conical follows WPS
             area (i,j) = area_wps (nest)
             
          case DEFAULT
             ! Area of a polygon areaint.m from Matlab
             call ij_to_latlon(LIS_domain(nest)%lisproj,c-0.5, r-0.5, lat_ll, lon_ll) ! SW corner (A)
             call ij_to_latlon(LIS_domain(nest)%lisproj,c-0.5, r+0.5, lat_ul, lon_ul) ! NW corner (B)         
             call ij_to_latlon(LIS_domain(nest)%lisproj,c+0.5, r+0.5, lat_ur, lon_ur) ! NE corner (C)        
             call ij_to_latlon(LIS_domain(nest)%lisproj,c+0.5, r-0.5, lat_lr, lon_lr) ! SE corner (D)
             area (i,j) = areaint((/lat_ll, lat_ul, lat_ur, lat_lr/), (/lon_ll, lon_ul, lon_ur, lon_lr/))
             
          END select
                    
       end do
    end do
    
  contains
    
    ! ----------------------------------------------------------------

    real function area_latlon (nest, lat)

      implicit none

      real, intent (in)    :: lat
      integer, intent (in) :: nest
      
      area_latlon = radius * radius * &
                 (sin(d2r(lat + 0.5*LIS_rc%gridDesc(nest,9))) - &
                  sin(d2r(lat - 0.5*LIS_rc%gridDesc(nest,9))))* &
                  d2r(LIS_rc%gridDesc(nest,10))/1000./1000.    ! [km2]
      
    end function area_latlon

    ! ----------------------------------------------------------------

    real function area_wps (nest)

      implicit none
      integer , intent (in) :: nest
      integer           :: rc
      real              :: DX, DY,  MSFTX, MSFTY

      MSFTY = 1.
      MSFTX = 1.

      DX = LIS_rc%gridDesc(nest,8)
      DY = LIS_rc%gridDesc(nest,9)
      area_wps = DX*DY/MSFTX/MSFTY
      
    end function area_wps
    
    ! ----------------------------------------------------------------
    
    real FUNCTION areaint (lat, lon)
            
      ! simplified from Matlab's areaint.m to compute area of a single polygon
      ! AREAINT Surface area of polygon on sphere 
      !   A = AREAINT(LAT,LON) calculates the spherical surface area of the
      !   polygon specified by the input vectors LAT, LON.  LAT and LON are in
      !   degrees.  The calculation uses a line integral approach.  The output,
      !   A, is the surface area fraction covered by the polygon on a unit
      !   sphere.
      
      implicit none
      real, intent(in), dimension(:)    :: lat, lon
      real, allocatable , dimension (:) :: latc, lonc, colat, az, integrands 
      real                              :: lat0,lon0,dlat,dlon,a,deltas,daz,colat2
      integer                           :: n, i
      
      n = size (lat) + 1
      allocate (latc (1:n))
      allocate (lonc (1:n))
      allocate (colat(1:n))
      allocate (az   (1:n))
      
      latc(1:n-1) = lat
      lonc(1:n-1) = lon
      latc(n)     = lat(1)
      lonc(n)     = lon(1)
      lat0 = 0.
      lon0 = 0.
      
      ! greatcircle distance, and greatcircle azimuth wrt 0.,0 (Matlab's distance.m)
      ! ----------------------------------------------------------------------------
      
      do i = 1,n
       
         latc(i) = d2r(latc(i))
         lonc(i) = d2r(lonc(i))
         dlat     = latc(i) - lat0
         dlon     = lonc(i) - lon0
         
         ! haversine
         a        = (sin(dlat/2.))**2 + cos(lat0)*cos(latc(i))*(sin(dlon/2.))**2
         if(a < 0.) a =0.
         if(a > 1.) a =1.
         
         colat(i) = 2.*atan2(sqrt(a),sqrt(1.-a))         
         az(i)    = atan2(cos(latc(i)) * sin(lonc(i)-lon0),  &
              cos(lat0) * sin(latc(i)) - sin(lat0) * cos(latc(i))* cos(lonc(i)-lon0))
         ! wrap az to the range 0-2pi
         az(i)    = az(i) - 2.*pi*floor(az(i)/2./pi)
         
      end do
      
      n = n -1
      allocate (integrands (1:n))
      
      do i = 1, n
         
         ! Calculate step sizes
         daz = az(i+1) - az(i)
         ! wrap to -pi <= daz <=pi
         daz = daz - 2.*pi*floor((daz+pi)/2./pi) 
         
         ! Determine average surface distance for each step
         deltas = (colat (i+1) - colat (i))/2.
         colat2 = colat(i) + deltas
         
         ! Integral over azimuth is 1-cos(colatitudes)
         integrands (i) = (1. - cos(colat2)) * daz
      end do
      
      areaint = abs (sum (integrands))/4./pi
      areaint = MIN (areaint, 1. - areaint)
      deallocate (integrands, latc, lonc, colat, az)
      
    end FUNCTION areaint

    ! ----------------------------------------------------------------
    
    function d2r (degree) result(rad)
      
      ! degrees to radians
      real,intent(in) :: degree
      real :: rad
      
      rad = degree*PI/180.
      
    end function d2r
  end SUBROUTINE get_area
   
end module alltypes_irrigationMod
