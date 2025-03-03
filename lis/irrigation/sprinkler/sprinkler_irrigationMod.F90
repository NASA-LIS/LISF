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
module sprinkler_irrigationMod
!BOP
!
! !MODULE: sprinkler_irrigationMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  11 Nov 2012: Sujay Kumar; Initial implementation
!  25 Feb 2020: Jessica Erlingis; Update irrigation window
!  14 Apr 2021: Wanshu Nie; Add support for GW/SW irrigation partitioning
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
  
  PUBLIC  :: sprinkler_irrigation_init
  PUBLIC  :: sprinkler_irrigation_updates

contains
  
  subroutine sprinkler_irrigation_init(irrigState)

    type(ESMF_State) :: irrigState(LIS_rc%nnest)

    integer              :: n 
    integer              :: rc, status
    type(ESMF_ArraySpec) :: arrspec1
    type(ESMF_Field)     :: irrigRateField, irrigFracField
    type(ESMF_Field)     :: irrigRootDepthField, irrigScaleField
    real,  allocatable   :: irrigFrac(:)
    real,  allocatable   :: irrigRootdepth(:)
    real,  allocatable   :: irrigScale(:)
    real, pointer        :: frac(:),scale(:),rootdepth(:),irrigrate(:)
    character(len=LIS_CONST_PATH_LEN) :: maxrootdepthfile

    type(ESMF_Field)    :: irriggwratioField
    real,  allocatable  :: irriggwratio(:)
    real,  pointer      :: gwratio(:)
! __________________________________________________________


    do n=1,LIS_rc%nnest
       allocate(irrigFrac(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(irrigRootDepth(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(irrigScale(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(irriggwratio(LIS_rc%npatch(n,LIS_rc%lsm_index))) 

       write(LIS_logunit,*) "[INFO] Running the 'Sprinkler' irrigation method ... "

     ! Read irrigation fraction (or "intensity") input:
       call read_irrigFrac(n, irrigFrac)
       call read_irriggwratio(n,irriggwratio)

     ! Read crop type maximum root depth file (combined landcover/crop classifications):
       call ESMF_ConfigGetAttribute(LIS_config,maxrootdepthfile,&
            label="Sprinkler irrigation max root depth file:",&
            rc=rc)
       call LIS_verify(rc,&
            'Sprinkler irrigation max root depth file: option not specified in the config file')

       call read_irrigRootdepth( n, maxrootdepthfile, irrigRootDepth )

     ! Calculate irrigation scale:
       call compute_irrigScale(n, irrigFrac, irrigScale)

       call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status, &
            "ESMF_ArraySpecSet failed in sprinkler_irrigation_init")

       irrigRateField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation rate", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in sprinkler_irrigation_init")

       call ESMF_FieldGet(irrigRateField,localDE=0,&
            farrayPtr=irrigrate,rc=status)
       call LIS_verify(status,'ESMF_FieldGet failed for irrigrate ')
       irrigrate = 0.0

       call ESMF_StateAdd(irrigState(n),(/irrigRateField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for irrigRate failed in sprinkler_irrigation_init")

       irrigFracField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation frac", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in sprinkler_irrigation_init")

       call ESMF_FieldGet(irrigFracField,localDE=0,&
            farrayPtr=frac,rc=status)
       call LIS_verify(status,'ESMF_FieldGet failed for IrrigFrac')
       
       frac = irrigFrac

       call ESMF_StateAdd(irrigState(n),(/irrigFracField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for irrigFrac failed in sprinkler_irrigation_init")
       deallocate(irrigFrac)

       irriggwratioField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Groundwater irrigation ratio", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in sprinkler_irrigation_init")

       call ESMF_FieldGet(irriggwratioField,localDE=0,&
            farrayPtr=gwratio,rc=status)
       call LIS_verify(status,'ESMF_FieldGet failed for irriggwratio')

       gwratio = irriggwratio

       call ESMF_StateAdd(irrigState(n),(/irriggwratioField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for irriggwratio failed in sprinkler_irrigation_init")
       deallocate(irriggwratio)

       irrigRootdepthField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation max root depth", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in sprinkler_irrigation_init")

       call ESMF_FieldGet(irrigRootdepthField,localDE=0,&
            farrayPtr=rootdepth,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet failed for root depth')
       rootdepth=irrigRootDepth
       
       call ESMF_StateAdd(irrigState(n),(/irrigRootdepthField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for max root depth failed in sprinkler_irrigation_init")
       deallocate(irrigRootDepth)

       irrigScaleField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation scale",rc=status)
       call LIS_verify(status, &
            'ESMF_FieldCreate failed in sprinkler_irrigation_init')
       
       call ESMF_FieldGet(irrigScaleField,localDE=0,&
            farrayPtr = scale,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet failed for irrigation scale')

       scale = irrigScale
       call ESMF_StateAdd(irrigState(n),(/irrigScaleField/),rc=status)
       call LIS_verify(status,&
            'ESMF_StateAdd for irrigation scale failed in sprinkler_irrigation_init')
       deallocate(irrigScale)
    enddo

  end subroutine sprinkler_irrigation_init
  

  subroutine sprinkler_irrigation_updates(n, irrigState)

    use LIS_FORC_AttributesMod 
    use LIS_histDataMod
    use LIS_metforcingMod, only : LIS_FORC_State    

    implicit none

    integer, intent(in) :: n 
    type(ESMF_State)    :: irrigState
    
    real, parameter     :: otimes = 6.0 ! local trigger check start time [hour]
    real, parameter     :: otimee = 10.0 !local trigger check end time [hour]
    integer             :: t,gid
    real                :: ltime
    integer             :: chhr, lhr
    integer             :: status

    type(ESMF_Field)    :: irrigRateField,prcpField
    real,    pointer    :: prcp(:)
    real,  allocatable  :: irrigAmt(:)
    real,    pointer    :: irrigRate(:)

    integer             :: tmpval

    real                 :: timestep, shift_otimes, shift_otimee

    allocate(irrigAmt(LIS_rc%npatch(n,LIS_rc%lsm_index)))
    call ESMF_StateGet(irrigState,&
         "Irrigation rate",&
         irrigRateField,rc=status)
    call LIS_verify(status,&
         'ESMF_StateGet failed for Irrigation rate')
    
    call ESMF_FieldGet(irrigRateField,localDE=0,&
         farrayPtr=irrigrate,rc=status)
    call LIS_verify(status,'ESMF_FieldGet failed for irrigrate ')

    call ESMF_StateGet(LIS_FORC_State(n),&
         trim(LIS_FORC_Rainf%varname(1)),prcpField,&
         rc=status)
    call LIS_verify(status,&
         'ESMF_StateGet failed for rainf in sprinkler_irrigation')

    
    call ESMF_FieldGet(prcpField,localDE=0, farrayPtr=prcp,rc=status)
    call LIS_verify(status,&
         'ESMF_FieldGet failed for rainf in sprinkler_irrigation')

    irrigAmt = 0.0

    timestep = LIS_rc%ts

    ! Adjust bounds by timestep to account for the fact that LIS_rc%hr, etc. will
    ! represents the END of the integration timestep window

    shift_otimes = otimes + (timestep/3600.)
    shift_otimee = otimee + (timestep/3600.)

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
       if((ltime.ge.shift_otimes).and.(ltime.lt.shift_otimee)) then           
          tmpval = LIS_domain(n)%tile(t)%index
          prcp(t)=prcp(t)+irrigRate(t)
          irrigAmt(t) = irrigRate(t)
       endif
    
     ! Write out irrigated amount of water to separate directory:
       if( LIS_MOC_IRRIGATEDWATER == 1 ) then
         call LIS_diagnoseIrrigationOutputVar(n,t,LIS_MOC_IRRIGATEDWATER,&
               value=irrigAmt(t),unit="kg m-2 s-1",direction="-",vlevel=1)
       else
         write(LIS_logunit,*) "[ERR] 'Irrigated water:' was not selected. Program stopping ..."
         call LIS_endrun
       endif

    enddo
    deallocate(irrigAmt)
    
  end subroutine sprinkler_irrigation_updates


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
    real,  allocatable   :: glb_frac1(:,:)
    
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
       call LIS_verify(ios,'nf90_get_var failed for in sprinkler_irrigationMod')

       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in sprinkler_irrigationMod')
       
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

  subroutine read_irriggwratio(n,gwratio)

    use LIS_fileIOMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer,  intent(in) :: n
    real                 :: gwratio(LIS_rc%npatch(n,LIS_rc%lsm_index))

    integer              :: t,col,row
    integer              :: nid,ios,status,gwratioId
    logical              :: file_exists
    real,  allocatable   :: l_gwratio(:,:)
    real,  allocatable   :: glb_gwratio(:,:)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    if ( LIS_rc%irrigation_GWabstraction .ne. 0 .and. &
         LIS_rc%irrigation_SourcePartition .ne. 0 ) then
       inquire(file=LIS_rc%paramfile(n), exist=file_exists)
       if(file_exists) then

          allocate(l_gwratio(LIS_rc%lnc(n),LIS_rc%lnr(n)))

          ios = nf90_open(path=LIS_rc%paramfile(n),&
             mode=NF90_NOWRITE,ncid=nid)
          call LIS_verify(ios,'Error in nf90_open in the lis input netcdf file')

          write(LIS_logunit,*) "[INFO] Reading in the groundwater irrigation ratio field ... "

          allocate(glb_gwratio(LIS_rc%gnc(n),LIS_rc%gnr(n)))

          ios = nf90_inq_varid(nid,'irriggwratio',gwratioId)
          call LIS_verify(ios,'nf90_inq_varid failed for irriggwratio')

          ios = nf90_get_var(nid,gwratioId, glb_gwratio)
          call LIS_verify(ios,'nf90_get_var failed for in sprinkler_irrigationMod')

          ios = nf90_close(nid)
          call LIS_verify(ios,'nf90_close failed in sprinkler_irrigationMod')
          l_gwratio(:,:) = glb_gwratio(&
             LIS_ews_halo_ind(n,LIS_localPet+1):&
             LIS_ewe_halo_ind(n,LIS_localPet+1),&
             LIS_nss_halo_ind(n,LIS_localPet+1):&
             LIS_nse_halo_ind(n,LIS_localPet+1))
          deallocate(glb_gwratio)

          do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
             col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
             row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
             gwratio(t) = l_gwratio(col,row)
          enddo

          deallocate(l_gwratio)
       else
          write(LIS_logunit,*) "[ERR] Groundwater irrigation ratio map: ",&
             LIS_rc%paramfile(n),"[ERR] does not exist."
          write(LIS_logunit,*) "Program stopping ..."
          call LIS_endrun
       endif
    else
       gwratio = LIS_rc%udef
    endif
#endif
  end subroutine read_irriggwratio

  subroutine read_irrigRootdepth(n, rdfile, rootdepth)

    use LIS_fileIOMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer,    intent(in) :: n 
    character(len=*)       :: rdfile
    real                   :: rootdepth(LIS_rc%npatch(n,LIS_rc%lsm_index))
    integer                :: total_vegtypes
!    integer, parameter     :: nt = 32 ! used to be hardcoded 
    real, allocatable      :: rootd(:)
    integer                :: ftn
    integer                :: t,j,col,row
    integer                :: nid, ios, status, croptypeId
    logical                :: file_exists    
    real,   allocatable    :: l_croptype(:,:)
    real,   allocatable    :: glb_croptype(:,:)
    real,   allocatable    :: glb_croptype1(:,:)
! __________________________________________________________________________
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

 !- Read in LDT input crop classification information (done here for now):
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
     ! Read in LDT-generated netcdf file information:
       write(LIS_logunit,*)"[INFO] Reading crop classification information ..."
       ios = nf90_open(path=LIS_rc%paramfile(n),&
                       mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in LIS_irrigation_init')
 
       ios = nf90_get_att(nid, NF90_GLOBAL, 'CROPCLASS_SCHEME', LIS_rc%cropscheme)
       call LIS_verify(ios,'Error in nf90_get_att in LIS_irrigation_init')
 
       ios = nf90_get_att(nid, NF90_GLOBAL, 'CROPCLASS_NUMBER', LIS_rc%numbercrops)
       call LIS_verify(ios,'Error in nf90_get_att in LIS_irrigation_init')
       write(LIS_logunit,*)"[INFO] Read in crop classfication: ",trim(LIS_rc%cropscheme),&
                          ", with the number of crop types:",LIS_rc%numbercrops
       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in sprinkler_irrigationMod')
    endif

  ! Estimate total crop types, added to landcover scheme class total:
    select case ( LIS_rc%lcscheme )
     case( "UMD" )
       total_vegtypes = 13 + LIS_rc%numbercrops
     case( "IGBP", "IGBPNCEP" )
       total_vegtypes = 20 + LIS_rc%numbercrops
     case( "USGS" )
       total_vegtypes = 24 + LIS_rc%numbercrops
     case default
       write(LIS_logunit,*) "[ERR] The landcover scheme, ",trim(LIS_rc%lcscheme),","
       write(LIS_logunit,*) "[ERR] is not supported for irrigation. Stopping program ... "
       call LIS_endrun()
    end select
  ! Assign default 32 UMD+CROPMAP for now, due to indexing for max root depth input files:
    total_vegtypes = 13 + LIS_rc%numbercrops

    allocate(l_croptype(LIS_rc%lnc(n),LIS_rc%lnr(n)))

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
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       ios = nf90_open(path=LIS_rc%paramfile(n),&
                       mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in the lis input netcdf file')

       write(LIS_logunit,*) "[INFO] Reading in the crop type field ... "
       
       allocate(glb_croptype(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       
       ios = nf90_inq_varid(nid,'CROPTYPE',croptypeId)
       call LIS_verify(ios,'nf90_inq_varid failed for CROPTYPE')
       
       ios = nf90_get_var(nid, croptypeId, glb_croptype)
       call LIS_verify(ios,'nf90_get_var failed for CROPTYPE')

       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in sprinkler_irrigationMod')
       
       l_croptype(:,:) = glb_croptype(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1): &
            LIS_nse_halo_ind(n,LIS_localPet+1))
       deallocate(glb_croptype)

       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
          row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
          
          if(l_croptype(col,row).gt.0) then 
             rootdepth(t) = rootd(nint(l_croptype(col,row)))
          else
             rootdepth(t) = 0 
          endif
       enddo
       deallocate( rootd )

    else
       write(LIS_logunit,*) "[ERR] The irrigation croptype map: ",&
             LIS_rc%paramfile(n)," does not exist."
       write(LIS_logunit,*) "[ERR] Program stopping ..."
       call LIS_endrun
    endif
    deallocate(l_croptype)
#endif

  end subroutine read_irrigRootdepth


  subroutine compute_irrigScale(n, irrigFrac, irrigScale)

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
     case default
       write(LIS_logunit,*) "[ERR] The landcover scheme, ",trim(LIS_rc%lcscheme),","
       write(LIS_logunit,*) "[ERR] is not supported for irrigation. Stopping program ... "
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
! logic for the more detailed (32 category map)
!             if(LIS_domain(n)%tile(t)%vegt.ge.14) then !crop tiles
!                crppix(gid) = crppix(gid) + LIS_domain(n)%tile(t)%fgrd
!             elseif(LIS_domain(n)%tile(t)%vegt.eq.10) then !grassland
!                grasspix(gid) = LIS_domain(n)%tile(t)%fgrd
!             elseif(LIS_domain(n)%tile(t)%vegt.gt.5.and.&
!                  LIS_domain(n)%tile(t)%vegt.lt.10) then !shrubs
!                restpix(gid) = restpix(gid) + LIS_domain(n)%tile(t)%fgrd
!             endif
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

end module sprinkler_irrigationMod
