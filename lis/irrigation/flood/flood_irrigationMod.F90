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
module flood_irrigationMod
!BOP
!
! !MODULE: flood_irrigationMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  11 Nov 2012: Sujay Kumar; Initial implementation
!  18 Jun 2014: Ben Zaitchik; Modified for flood
!
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
  
  PUBLIC  :: flood_irrigation_init
  PUBLIC  :: flood_irrigation_updates

contains
  
  subroutine flood_irrigation_init(irrigState)


    type(ESMF_State) :: irrigState(LIS_rc%nnest)

    integer              :: n 
    integer              :: rc, status
    type(ESMF_ArraySpec) :: arrspec1
    type(ESMF_Field)     :: irrigRateField, irrigFracField
    type(ESMF_Field)     :: irrigRootDepthField, irrigScaleField
    real,  pointer       :: irrigrate(:)
    real,  pointer       :: irrigFrac(:), frac(:)
    real,  pointer       :: irrigRootdepth(:), rootdepth(:)
    real,  pointer       :: irrigScale(:),scale(:)
    character(len=LIS_CONST_PATH_LEN) :: maxrootdepthfile

    do n=1,LIS_rc%nnest
       allocate(irrigFrac(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(irrigRootDepth(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(irrigScale(LIS_rc%npatch(n,LIS_rc%lsm_index)))

       write(LIS_logunit,*) " Running the 'Flood' irrigation method ... "

       call ESMF_ConfigGetAttribute(LIS_config,maxrootdepthfile,&
            label="Flood irrigation max root depth file:",&
            rc=rc)
       call LIS_verify(rc,&
            'Flood irrigation max root depth file: option not specified in the config file')

       call read_irrigfrac(n, irrigFrac)
       call read_irrigrootdepth(n, maxrootdepthfile, irrigRootDepth)

       call compute_irrigScale(n,irrigFrac, irrigScale)

       call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status, &
            "ESMF_ArraySpecSet failed in flood_irrigation_init")

       irrigRateField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation rate", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in flood_irrigation_init")
       
       call ESMF_StateAdd(irrigState(n),(/irrigRateField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for irrigRate failed in flood_irrigation_init")
       
       irrigFracField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation frac", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in flood_irrigation_init")

       call ESMF_FieldGet(irrigFracField,localDE=0,&
            farrayPtr=frac,rc=status)
       call LIS_verify(status,'ESMF_FieldGet failed for IrrigFrac')
       
       frac = irrigFrac

       call ESMF_StateAdd(irrigState(n),(/irrigFracField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for irrigFrac failed in flood_irrigation_init")
       deallocate(irrigFrac)

       irrigRootdepthField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation max root depth", rc=status)
       call LIS_verify(status, &
            "ESMF_FieldCreate failed in flood_irrigation_init")

       call ESMF_FieldGet(irrigRootdepthField,localDE=0,&
            farrayPtr=rootdepth,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet failed for root depth')
       rootdepth=irrigRootDepth
       
       call ESMF_StateAdd(irrigState(n),(/irrigRootdepthField/),rc=status)
       call LIS_verify(status,&
            "ESMF_StateAdd for max root depth failed in flood_irrigation_init")
       deallocate(irrigRootDepth)

       irrigScaleField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1,&
            name="Irrigation scale",rc=status)
       call LIS_verify(status, &
            'ESMF_FieldCreate failed in flood_irrigation_init')
       
       call ESMF_FieldGet(irrigScaleField,localDE=0,&
            farrayPtr = scale,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet failed for irrigation scale')

       scale = irrigScale
       call ESMF_StateAdd(irrigState(n),(/irrigScaleField/),rc=status)
       call LIS_verify(status,&
            'ESMF_StateAdd for irrigation scale failed in flood_irrigation_init')
       deallocate(irrigScale)
    enddo


  end subroutine flood_irrigation_init
  

  subroutine flood_irrigation_updates(n, irrigState)

    use LIS_FORC_AttributesMod 
    use LIS_histDataMod
    use LIS_metforcingMod, only : LIS_FORC_State    

    implicit none

    integer, intent(in) :: n 
    type(ESMF_State)    :: irrigState
    
    real, parameter     :: otimes = 6.0 ! local trigger check start time [hour]
    real, parameter     :: irrhrf = 0.5 !duration of flood irrigation [hour]
    integer             :: t,gid
    real                :: ltime,otimee
    integer             :: chhr, lhr
    integer             :: status

    type(ESMF_Field)    :: irrigRateField,prcpField
    real,    pointer    :: prcp(:)
    real                :: irrigAmt(LIS_rc%npatch(n,LIS_rc%lsm_index))
    real,    pointer    :: irrigRate(:)

    call ESMF_StateGet(irrigState,&
         "Irrigation rate",&
         irrigRateField,rc=status)
    call LIS_verify(status,&
         'ESMF_StateGet failed for Irrigation rate')
    
    call ESMF_FieldGet(irrigRateField,localDE=0,&
         farrayPtr=irrigrate,rc=status)
    call LIS_verify(status,'ESMF_FieldGet failed for irrigrate ')

!    call ESMF_StateGet(LIS_FORC_State(n),&
!         trim(LIS_FORC_Rainf%varname(1)),prcpField,&
!         rc=status)
!    call LIS_verify(status,&
!         'ESMF_StateGet failed for rainf in flood_irrigation')

    
!    call ESMF_FieldGet(prcpField,localDE=0, farrayPtr=prcp,rc=status)
!    call LIS_verify(status,&
!         'ESMF_FieldGet failed for rainf in flood_irrigation')

    irrigAmt = 0.0
    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

       gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
       chhr = nint(24.0*(LIS_domain(n)%grid(gid)%lon/360.0))
       if((LIS_domain(n)%grid(gid)%lon.lt.0.0).and.&
            (abs(mod(LIS_domain(n)%grid(gid)%lon,15.0)).ge.0.0001)) &
          chhr = chhr -1
       lhr = LIS_rc%hr +chhr
       if(lhr.ge.24) lhr = lhr-24
       if(lhr.lt.0) lhr = lhr+24
       
       ltime = real(lhr)+real(LIS_rc%mn)/60.0+real(LIS_rc%ss)/3600.0
       otimee = otimes + irrhrf
       if((ltime.ge.otimes).and.(ltime.lt.otimee)) then           
          irrigAmt(t) = irrigRate(t)
       endif
       call LIS_diagnoseIrrigationOutputVar(n,t,LIS_MOC_IRRIGATEDWATER,&
            value=irrigAmt(t),unit="kg m-2 s-1",direction="-",vlevel=1)
    enddo
    
  end subroutine flood_irrigation_updates


  subroutine read_irrigFrac(n,frac)
    use LIS_fileIOMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer,      intent(in) :: n 
    real                     :: frac(LIS_rc%npatch(n,LIS_rc%lsm_index))

    integer                  :: t,col,row
    integer                  :: nid,ios,status,fracId
    logical                  :: file_exists    
    real                     :: l_frac(LIS_rc%lnc(n),LIS_rc%lnr(n))
    real,         pointer    :: glb_frac(:,:)
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 

       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in the lis input netcdf file')

       write(LIS_logunit,*) " Reading in the irrigation fraction field ... "
       
       allocate(glb_frac(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       
       ios = nf90_inq_varid(nid,'IRRIGFRAC',fracId)
       call LIS_verify(ios,'nf90_inq_varid failed for IRRIGFRAC')
       
       ios = nf90_get_var(nid,fracId, glb_frac)
       call LIS_verify(ios,'nf90_get_var failed for in flood_irrigationMod')

       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in flood_irrigationMod')
       
       l_frac(:,:) = glb_frac(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1): &
            LIS_nse_halo_ind(n,LIS_localPet+1))
       deallocate(glb_frac)

       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
          row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
          
          frac(t) = l_frac(col,row)
       enddo

    else
       write(LIS_logunit,*) 'irrigation frac map: ',&
            LIS_rc%paramfile(n), ' does not exist'
       write(LIS_logunit,*) 'program stopping ...'
       call LIS_endrun
    endif
#endif
  end subroutine read_irrigFrac


  subroutine read_irrigRootdepth(n,rdfile, rootdepth)
    use LIS_fileIOMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer,      intent(in) :: n 
    character(len=*)         :: rdfile
    real                     :: rootdepth(LIS_rc%npatch(n,LIS_rc%lsm_index))

    integer,   parameter     :: nt = 32 !hardcoded for now
    integer                  :: ftn
    real                     :: rootd(nt)
    integer                  :: t,j,col,row
    integer                  :: nid,ios,status,croptypeId
    logical                  :: file_exists    
    real                     :: l_croptype(LIS_rc%lnc(n),LIS_rc%lnr(n))
    real,         pointer    :: glb_croptype(:,:)
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=rdfile,exist=file_exists)
    if(file_exists) then 
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=rdfile,status='old')
       read(ftn,*) (rootd(j),j=1,nt)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) 'Max root depth file ',trim(rdfile), ' not found'
       call LIS_endrun()
    endif

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 

       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in the lis input netcdf file')

       write(LIS_logunit,*) " Reading in the crop type field ... "
       
       allocate(glb_croptype(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       
       ios = nf90_inq_varid(nid,'CROPTYPE',croptypeId)
       call LIS_verify(ios,'nf90_inq_varid failed for CROPTYPE')
       
       ios = nf90_get_var(nid,croptypeId, glb_croptype)
       call LIS_verify(ios,'nf90_get_var failed for CROPTYPE')

       ios = nf90_close(nid)
       call LIS_verify(ios,'nf90_close failed in flood_irrigationMod')
       
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

    else
       write(LIS_logunit,*) 'irrigation croptype map: ',&
            LIS_rc%paramfile(n), ' does not exist'
       write(LIS_logunit,*) 'program stopping ...'
       call LIS_endrun
    endif
#endif
  end subroutine read_irrigRootdepth

  subroutine compute_irrigScale(n,irrigFrac, irrigScale)

    integer                :: n 
    real                   :: irrigFrac(LIS_rc%npatch(n,LIS_rc%lsm_index))
    real                   :: irrigScale(LIS_rc%npatch(n,LIS_rc%lsm_index))
    
    integer                :: t,gid,vegt
    real                   :: crppix(LIS_rc%ngrid(n))
    real                   :: grasspix(LIS_rc%ngrid(n))
    real                   :: restpix(LIS_rc%ngrid(n))
    real                   :: irrpix,excess
    integer                :: crop1,crop2,grass,shrub1,shrub2
    crppix    = 0.0
    grasspix  = 0.0
    restpix   = 0.0

!------------------------------------------------------------------------
! WARNING: The following code is valid only for the no-tiling or the 
! vegetation only tiling. The fgrd values are not valid when multiple
! modes of tiling are turned on. 
!------------------------------------------------------------------------
    if(LIS_rc%lcscheme.eq."UMD") then !UMD
      crop1 = 11
      crop2 = 11
      grass = 10 
      shrub1 = 6
      shrub2 = 9
   elseif(LIS_rc%lcscheme.eq."MODIS".or.LIS_rc%lcscheme.eq."IGBPNCEP") then 
      crop1 = 12
      crop2 = 14
      grass = 10 
      shrub1 = 6
      shrub2 = 9
   elseif(LIS_rc%lcscheme.eq."USGS") then !UMD
      crop1 = 2
      crop2 = 6
      grass = 7 
      shrub1 = 8
      shrub2 = 10
   else
      write(LIS_logunit,*) 'The landcover scheme ',trim(LIS_rc%lcscheme)
      write(LIS_logunit,*) 'is not supported for irrigation '
      call LIS_endrun()
   endif
   
   do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
      if(LIS_domain(n)%tile(t)%vegt.ne.LIS_rc%waterclass.and.&
           LIS_domain(n)%tile(t)%vegt.ne.LIS_rc%urbanclass) then 
         if(LIS_domain(n)%tile(t)%vegt.ge.crop1.and.&
              LIS_domain(n)%tile(t)%vegt.le.crop2) then !crop tiles
            crppix(gid) = crppix(gid) + LIS_domain(n)%tile(t)%fgrd* & 
                 LIS_domain(n)%tile(t)%pens
         elseif(LIS_domain(n)%tile(t)%vegt.eq.grass) then !grassland
            grasspix(gid) = LIS_domain(n)%tile(t)%fgrd*& 
                 LIS_domain(n)%tile(t)%pens
         elseif(LIS_domain(n)%tile(t)%vegt.ge.shrub1.and.&
              LIS_domain(n)%tile(t)%vegt.le.shrub2) then !shrubs
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
   
   irrigScale =1.0
   do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
      vegt = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
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

  end subroutine compute_irrigScale

end module flood_irrigationMod
