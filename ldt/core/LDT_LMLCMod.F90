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
module LDT_LMLCMod
!BOP
!
! !MODULE: LDT_LMLCMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read various sources of
!  landmask and landcover data. 
! 
!  \subsubsection{Overview}
!  This module provides routines for reading and modifying landmask 
!  and landcover data. 
!
! !REVISION HISTORY:
!
!  18 Jul 2008: Sujay Kumar; Initial implementation
!  18 Jul 2013: KR Arsenault; Expanded options
!  30 Nov 2018: David Mocko; Added Bondville landcover classification
!
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_SurfaceTypeMod, only: LDT_assign_lcsfctype
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_LMLC_init            ! initializes data structures and memory
  public :: LDT_LMLC_writeHeader
  public :: LDT_LMLC_writeData
  
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------

!BOP 
! 
! !ROUTINE: LDT_LMLC_writeHeader
! \label{LDT_LMLC_writeHeader}
! 
! !INTERFACE:
  interface LDT_LMLC_writeHeader
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure LMLC_writeHeader_LIS
     module procedure LMLC_writeHeader_LISHydro
! 
! !DESCRIPTION:
! This interface provides routines for writing NETCDF header both 
! in the standard preprocessing mode for LIS 
! as well as in the LISHydro(WRFHydro) 
! preprocessing mode.
!EOP 
  end interface

!BOP 
! 
! !ROUTINE: LDT_LMLC_writeHeader
! \label{LDT_LMLC_writeHeader}
! 
! !INTERFACE:
  interface LDT_LMLC_init
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure LMLC_init_LIS
     module procedure LMLC_init_LISHydro
! 
! !DESCRIPTION:
! This interface provides routines for writing NETCDF header both 
! in the standard preprocessing mode for LIS as well as LISHydro(WRFHydro) 
! preprocessing mode.
!EOP 
  end interface

!BOP 
! 
! !ROUTINE: LDT_LMLC_writeHeader
! \label{LDT_LMLC_writeHeader}
! 
! !INTERFACE:
  interface LDT_LMLC_writeData
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure LMLC_writeData_LIS
     module procedure LMLC_writeData_LISHydro
! 
! !DESCRIPTION:
! This interface provides routines for writing NETCDF header both 
! in standard preprocessing mode for LIS as well as LISHydro(WRFHydro) 
! preprocessing mode.
!EOP 
  end interface


!EOP
contains

!BOP
! 
! !ROUTINE: LMLC_init_LIS
! \label{LMLC_init_LIS}
! 
! !INTERFACE:
  subroutine LMLC_init_LIS()

! !USES:
    use ESMF
    use LDT_coreMod,  only : LDT_rc, LDT_config, LDT_domain
    use LDT_logMod,   only : LDT_verify, LDT_logunit
    use LDT_fileIOMod,only : LDT_readDomainConfigSpecs
    use LDT_paramOptCheckMod, only: LDT_LMLCOptChecks, LDT_gridOptChecks
! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! landmask and landcover datasets in the standard 
! preprocessing mode for LIS
!
!EOP
   implicit none
   integer :: rc
   integer :: n, c, r
   real, allocatable  :: landcover_fgrd(:,:,:)
   type(LDT_fillopts) :: landcover
! _________________________________________________


  write(LDT_logunit,*)" [INFO] - - - - - - - - - Landcover/Landmask Parameters - - - - - - - - - - -"

  ! Initialize gridcell water fraction (default value 0.5):
    LDT_rc%gridcell_water_frac = 0.5    

    allocate(LDT_rc%lc_gridDesc(LDT_rc%nnest,20))
    allocate(LDT_rc%mask_gridDesc(LDT_rc%nnest,20))
    allocate(LDT_rc%reg_gridDesc(LDT_rc%nnest,20))
    
    call ESMF_ConfigFindLabel(LDT_config,"Water fraction cutoff value:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%gridcell_water_frac(n),rc=rc)
       call LDT_verify(rc,'"Water fraction cutoff value: not specified')
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Create or readin landmask:",rc=rc)
    do n = 1, LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%mask_type(n),rc=rc)
       call LDT_verify(rc,"Create or read-in landmask option: not specified")
    enddo
    
    do n=1,LDT_rc%nnest

    !- Read-in land mask config entries:
       if( trim(LDT_rc%mask_type(n)) == "readin" ) then

         call ESMF_ConfigFindLabel(LDT_config,"Landmask file:",rc=rc)
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%mfile(n),rc=rc)
         call LDT_verify(rc,'Landmask file: not specified')

         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%mask_proj,&
              label="Landmask map projection:",rc=rc)
         call LDT_verify(rc,'Landmask map projection: option not specified in the config file')

         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%mask_gridtransform(n),&
              label="Landmask spatial transform:",rc=rc)
         call LDT_verify(rc,'Landmask spatial transform: option not specified in the config file')

         call LDT_readDomainConfigSpecs("Landmask", LDT_rc%mask_proj, LDT_rc%mask_gridDesc)

       end if
    enddo

    ! Option to assign all landmask and domainmask points to land:
    ! 0 - "off"; 1 - "on", where all points are set to 1
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%allmaskland,&
         label="Set all land and domain masks to 1:",&
         default=0,rc=rc)
    call LDT_verify(rc,'Set all land and domain masks to 1: not defined')

!-- Landcover dataset file and option inputs:

    call ESMF_ConfigFindLabel(LDT_config,"Landcover file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%vfile(n),rc=rc)
       call LDT_verify(rc,'Landcover file: not specified')
    enddo

    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%lc_proj,&
         label="Landcover map projection:",rc=rc)
    call LDT_verify(rc,'Landcover map projection: option not specified in the config file')

    call ESMF_ConfigFindLabel(LDT_config,"Landcover spatial transform:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%lc_gridtransform(n),rc=rc)
       call LDT_verify(rc,'Landcover spatial transform: not specified')
    enddo

  ! Don't need to read in "Native" grid extents/resolution, just for LIS inputs:
    do n=1,LDT_rc%nnest
       if( index(LDT_LSMparam_struc(n)%landcover%source,"Native").eq.0 .and. &
            index(LDT_LSMparam_struc(n)%landcover%source,"CONSTANT").eq.0 .and.&
            index(LDT_LSMparam_struc(n)%landcover%source,"MCD12Q1").eq.0) then
          
         call LDT_readDomainConfigSpecs("Landcover", LDT_rc%lc_proj, LDT_rc%lc_gridDesc)
         if( LDT_rc%lc_proj == "latlon" ) then
           call LDT_gridOptChecks( n, "Landcover", &
                LDT_rc%lc_gridtransform(n), &
                LDT_rc%lc_proj,&
                LDT_rc%lc_gridDesc(n,9) )
         endif
       endif
    enddo

    landcover%filltype = "none"
    call ESMF_ConfigGetAttribute(LDT_config, landcover%filltype, &
         label="Landcover fill option:",rc=rc)
    call LDT_verify(rc,"Landcover fill option: option not specified in the config file")

    ! Check if landcover fill is turned on when Source = 'CONSTANT':
    do n=1,LDT_rc%nnest
       if( LDT_LSMparam_struc(n)%landcover%source == "CONSTANT" .and. &
           landcover%filltype == "none" ) then
          write(LDT_logunit,*)"[ERR] 'CONSTANT' landcover source type selected"
          write(LDT_logunit,*)"[ERR]  and fill type set to 'none'. If wanting a"
          write(LDT_logunit,*)"[ERR]  constant land class for your LIS run, then" 
          write(LDT_logunit,*)"[ERR]  a landcover type needs to be set using"
          write(LDT_logunit,*)"[ERR]  the 'fill' options. You can specify 'neighbor'"
          write(LDT_logunit,*)"[ERR]  and for the 'fill value' specify the land class"
          write(LDT_logunit,*)"[ERR]  you need.  LDT run ending ... "
          call LDT_endrun
       endif
    enddo

    if( landcover%filltype == "neighbor" .or. landcover%filltype == "average" ) then
      write(LDT_logunit,*) "[INFO] Parameter-Mask Agreement Selected for Landcover (",&
                           trim(landcover%filltype),")"
      call ESMF_ConfigGetAttribute(LDT_config, landcover%fillradius, &
           label="Landcover fill radius:",rc=rc)
      call LDT_verify(rc,"Landcover fill radius: option not specified in the config file")

      call ESMF_ConfigGetAttribute(LDT_config, landcover%fillvalue, &
           label="Landcover fill value:",rc=rc)
      call LDT_verify(rc,"Landcover fill value: option not specified in the config file")

    elseif( landcover%filltype == "none" ) then
      write(LDT_logunit,*) " -- 'NONE' Parameter-Mask Agreement Option Selected for Landcover"
    else
         write(LDT_logunit,*) "[ERR] Fill option for Landcover is not valid: ",trim(landcover%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none or neighbor "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
    end if

!-- Regional mask dataset file and option inputs:
    if( LDT_LSMparam_struc(1)%regmask%selectOpt > 0 ) then

       LDT_rc%cliplandmask = .false.
       call ESMF_ConfigFindLabel(LDT_config,"Clip landmask with regional mask:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%cliplandmask(n),rc=rc)
          call LDT_verify(rc,'Clip landmask with regional mask: not specified')
       enddo

       call ESMF_ConfigFindLabel(LDT_config,"Regional mask file:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%regfile(n),rc=rc)
          call LDT_verify(rc,'Regional mask file: not specified')
       enddo
       
       call ESMF_ConfigFindLabel(LDT_config,"Regional mask spatial transform:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%reg_gridtransform(n),rc=rc)
          call LDT_verify(rc,'Regional mask spatial transform: option not specified in the config file')
       enddo

       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%reg_proj,&
            label="Regional mask map projection:",rc=rc)
       call LDT_verify(rc,'Regional mask map projection: option not specified in the config file')

       call LDT_readDomainConfigSpecs("Regional mask", LDT_rc%reg_proj, LDT_rc%reg_gridDesc)
    end if


!-- Allocate landcover and landmask file arrays:
    do n = 1, LDT_rc%nnest
       allocate(LDT_LSMparam_struc(n)%dommask%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            (LDT_LSMparam_struc(n)%landmask%num_bins)))

       allocate(LDT_LSMparam_struc(n)%landmask%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            (LDT_LSMparam_struc(n)%landmask%num_bins)))
       
       allocate(LDT_LSMparam_struc(n)%landmask2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            LDT_LSMparam_struc(n)%landmask%num_bins))
       
       allocate(LDT_LSMparam_struc(n)%watermask%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            LDT_LSMparam_struc(n)%landmask%num_bins))
       
     ! Fill in landmask placeholder parameter entries:
       call populate_param_attribs( "MASK2", &
            "Mask used for parameter consistency", "-",     &
            LDT_LSMparam_struc(n)%landmask, &
            LDT_LSMparam_struc(n)%landmask2 )
       
     ! Fill in water mask parameter entries:
       call populate_param_attribs( "WATERMASK", &
            "Water Mask used for water points", "-",     &
            LDT_LSMparam_struc(n)%landmask, &
            LDT_LSMparam_struc(n)%watermask )

       allocate(landcover_fgrd( LDT_rc%lnc(n),LDT_rc%lnr(n),&
            LDT_LSMparam_struc(n)%landcover%num_bins) )
       
       allocate(LDT_LSMparam_struc(n)%landcover%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            LDT_LSMparam_struc(n)%sfctype%vlevels))
       
       if( LDT_LSMparam_struc(n)%regmask%selectOpt == 1 ) then
          allocate(LDT_LSMparam_struc(n)%regmask%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%regmask%num_bins))
       endif

    !- Assign landmask source, in relation to original source or landcover:
       LDT_rc%mask_source(n) = LDT_LSMparam_struc(n)%landmask%source
       if( trim(LDT_rc%mask_type(n)) == "create" .or. &
            trim(LDT_LSMparam_struc(n)%landcover%source) == "ISA" ) then
          LDT_rc%mask_source(n) = LDT_LSMparam_struc(n)%landcover%source
       endif

    
!-- Perform checks on available landcover/landmask map input options: 

       if( LDT_LSMparam_struc(n)%regmask%selectOpt > 0 ) then
          call LDT_gridOptChecks( n, "Regional mask", &
               LDT_rc%reg_gridtransform(n), &
               LDT_rc%reg_proj, LDT_rc%reg_gridDesc(n,9) ) 
       endif

       call LDT_LMLCOptChecks( "Landcover", &
            LDT_rc%lc_type(n), LDT_rc%nt, LDT_rc%lc_gridtransform(n) )
       
       
!-- Read land cover and mask files:

       write(LDT_logunit,*) "[INFO] Reading landcover values"

    !- Landcover:
       call readlandcover(&
            trim(LDT_LSMparam_struc(n)%landcover%source)//char(0),&
            n, LDT_LSMparam_struc(n)%landcover%num_bins, &
            landcover_fgrd,                              &
            LDT_LSMparam_struc(n)%landmask%value  )
       
       LDT_LSMparam_struc(n)%landcover%vlevels = &
            LDT_LSMparam_struc(n)%sfctype%num_bins
       LDT_LSMparam_struc(n)%landcover%value = 0.
       LDT_LSMparam_struc(n)%landcover%value(:,:,1:LDT_LSMparam_struc(n)%landcover%num_bins) = &
            landcover_fgrd(:,:,1:LDT_LSMparam_struc(n)%landcover%num_bins)

    !- Copy readin/created landmask to place holder:
       LDT_LSMparam_struc(n)%landmask2%value(:,:,:) = &
            LDT_LSMparam_struc(n)%landmask%value(:,:,:)
       
    !- Reverse landmask to be a watermask for lake/open-water map fills:
       LDT_LSMparam_struc(n)%watermask%value(:,:,:) = 0.
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             if( LDT_LSMparam_struc(n)%landmask%value(c,r,1) == 1 ) then
                LDT_LSMparam_struc(n)%watermask%value(c,r,1) = 0.
             elseif( LDT_LSMparam_struc(n)%landmask%value(c,r,1) == 0 ) then
                LDT_LSMparam_struc(n)%watermask%value(c,r,1) = 1.
             endif
          enddo
       enddo

     ! Fill where parameter values are missing compared to land/water mask:
       if( landcover%filltype == "neighbor" ) then
          write(LDT_logunit,*) "[INFO] Checking/filling mask values for: ", &
                             trim(LDT_LSMparam_struc(n)%landcover%short_name)
          write(fill_logunit,*) "Checking/filling mask values for: ", &
                             trim(LDT_LSMparam_struc(n)%landcover%short_name)
          landcover%watervalue = LDT_rc%waterclass

          ! EMK...Make sure filling ignores water points.
          call LDT_discreteParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                LDT_rc%lc_gridtransform(n),                        &
                LDT_LSMparam_struc(n)%landcover%num_bins,          &
                LDT_LSMparam_struc(n)%landcover%value, landcover%watervalue, &
                LDT_LSMparam_struc(n)%landmask2%value,             &
                landcover%filltype, landcover%fillvalue, &
                landcover%fillradius , force_exclude_water=.true.)
       endif

    !- Surface type field update:
       call LDT_assign_lcsfctype( n, &
            LDT_LSMparam_struc(n)%sfctype%num_bins,   & 
            LDT_LSMparam_struc(n)%landcover%num_bins, &
            LDT_LSMparam_struc(n)%sfctype%value,      &
            LDT_LSMparam_struc(n)%landmask%value,     &
            LDT_LSMparam_struc(n)%landcover%value )
       write(LDT_logunit,*) "[INFO] Finished assigning LSM surface types"

    ! __________________________________________________

    !- Regional mask:
       if( LDT_LSMparam_struc(n)%regmask%selectOpt > 0 ) then
         call readregmask(&
              trim(LDT_LSMparam_struc(n)%regmask%source)//char(0),&
              n, LDT_LSMparam_struc(n)%regmask%value )

         LDT_LSMparam_struc(n)%regmask%vlevels = &
             LDT_LSMparam_struc(n)%regmask%num_bins

      !- Use the Regional land mask to "clip" the main landmask:
         if( LDT_rc%cliplandmask(n) .eqv. .true. ) then
           write(LDT_logunit,*) "[INFO] -- Using Regional Mask to Clip Landmask"
           do r = 1, LDT_rc%lnr(n)
             do c = 1, LDT_rc%lnc(n)
                if( LDT_LSMparam_struc(n)%regmask%value(c,r,1) == 0 ) then
                  LDT_LSMparam_struc(n)%landmask%value(c,r,1) = 0.
                endif
             enddo
           enddo
        endif

     endif
     LDT_LSMparam_struc(n)%dommask%value(:,:,:) = &
          LDT_LSMparam_struc(n)%landmask%value(:,:,:)

     deallocate( landcover_fgrd )
   enddo

   write(LDT_logunit,*) "[INFO] Finished reading landmask and landcover data"

 end subroutine LMLC_init_LIS

!BOP
! 
! !ROUTINE: LMLC_init_LISHydro
! \label{LMLC_init_LISHydro}
! 
! !INTERFACE:
  subroutine LMLC_init_LISHydro(flag)

! !USES:
    use ESMF
    use LDT_coreMod,  only : LDT_rc, LDT_config, LDT_domain
    use LDT_logMod,   only : LDT_verify, LDT_logunit
    use LDT_fileIOMod,only : LDT_readDomainConfigSpecs
    use LDT_paramOptCheckMod, only: LDT_LMLCOptChecks, LDT_gridOptChecks
! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! landmask and landcover datasets in the preprocessing
! mode for LIShydro
!
!EOP
   implicit none
   integer :: rc
   integer :: n, c, r, t
   real, allocatable  :: landcover_fgrd(:,:,:)
   type(LDT_fillopts) :: landcover
   integer            :: flag

   integer             :: nc,nr,nlctypes, ntxtypes
   real, allocatable   :: luindex(:,:),sctdom(:,:)
   real                :: maxv, domv
   real, allocatable   :: landcover1(:,:,:), texture1(:,:,:)
! _________________________________________________


  write(LDT_logunit,*)" [INFO] - - - - - - - - - Landcover/Landmask Parameters - - - - - - - - - - -"

  ! Initialize gridcell water fraction (default value 0.5):
    LDT_rc%gridcell_water_frac = 0.5    

    allocate(LDT_rc%lc_gridDesc(LDT_rc%nnest,20))
    allocate(LDT_rc%mask_gridDesc(LDT_rc%nnest,20))
    allocate(LDT_rc%reg_gridDesc(LDT_rc%nnest,20))
    
    call ESMF_ConfigFindLabel(LDT_config,"Water fraction cutoff value:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%gridcell_water_frac(n),rc=rc)
       call LDT_verify(rc,'"Water fraction cutoff value: not specified')
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Create or readin landmask:",rc=rc)
    do n = 1, LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%mask_type(n),rc=rc)
       call LDT_verify(rc,"Create or read-in landmask option: not specified")
    enddo
    
    do n=1,LDT_rc%nnest

    !- Read-in land mask config entries:
       if( trim(LDT_rc%mask_type(n)) == "readin" ) then

         call ESMF_ConfigFindLabel(LDT_config,"Landmask file:",rc=rc)
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%mfile(n),rc=rc)
         call LDT_verify(rc,'Landmask file: not specified')

         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%mask_proj,&
              label="Landmask map projection:",rc=rc)
         call LDT_verify(rc,'Landmask map projection: option not specified in the config file')

         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%mask_gridtransform(n),&
              label="Landmask spatial transform:",rc=rc)
         call LDT_verify(rc,'Landmask spatial transform: option not specified in the config file')

         call LDT_readDomainConfigSpecs("Landmask", LDT_rc%mask_proj, LDT_rc%mask_gridDesc)

       end if
    enddo

!-- Landcover dataset file and option inputs:

    call ESMF_ConfigFindLabel(LDT_config,"Landcover file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%vfile(n),rc=rc)
       call LDT_verify(rc,'Landcover file: not specified')
    enddo

    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%lc_proj,&
         label="Landcover map projection:",rc=rc)
    call LDT_verify(rc,'Landcover map projection: option not specified in the config file')

    call ESMF_ConfigFindLabel(LDT_config,"Landcover spatial transform:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%lc_gridtransform(n),rc=rc)
       call LDT_verify(rc,'Landcover spatial transform: not specified')
    enddo

  ! Don't need to read in "Native" grid extents/resolution, just for LIS inputs:
    do n=1,LDT_rc%nnest
       if( index(LDT_LSMparam_struc(n)%landcover%source,"Native").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%landcover%source,"CONSTANT").eq.0 ) then
         call LDT_readDomainConfigSpecs("Landcover", LDT_rc%lc_proj, LDT_rc%lc_gridDesc)
         if( LDT_rc%lc_proj == "latlon" ) then
           call LDT_gridOptChecks( n, "Landcover", &
                LDT_rc%lc_gridtransform(n), &
                LDT_rc%lc_proj,&
                LDT_rc%lc_gridDesc(n,9) )
         endif
       endif
    enddo

    landcover%filltype = "none"
    call ESMF_ConfigGetAttribute(LDT_config, landcover%filltype, &
         label="Landcover fill option:",rc=rc)
    call LDT_verify(rc,"Landcover fill option: option not specified in the config file")

    ! Check if landcover fill is turned on when Source = 'CONSTANT':
    do n=1,LDT_rc%nnest
       if( LDT_LSMparam_struc(n)%landcover%source == "CONSTANT" .and. &
           landcover%filltype == "none" ) then
          write(LDT_logunit,*)"[ERR] 'CONSTANT' landcover source type selected"
          write(LDT_logunit,*)"[ERR]  and fill type set to 'none'. If wanting a"
          write(LDT_logunit,*)"[ERR]  constant land class for your LIS run, then" 
          write(LDT_logunit,*)"[ERR]  a landcover type needs to be set using"
          write(LDT_logunit,*)"[ERR]  the 'fill' options. You can specify 'neighbor'"
          write(LDT_logunit,*)"[ERR]  and for the 'fill value' specify the land class"
          write(LDT_logunit,*)"[ERR]  you need.  LDT run ending ... "
          call LDT_endrun
       endif
    enddo

    if( landcover%filltype == "neighbor" .or. landcover%filltype == "average" ) then
      write(LDT_logunit,*) "[INFO] Parameter-Mask Agreement Selected for Landcover (",&
                           trim(landcover%filltype),")"
      call ESMF_ConfigGetAttribute(LDT_config, landcover%fillradius, &
           label="Landcover fill radius:",rc=rc)
      call LDT_verify(rc,"Landcover fill radius: option not specified in the config file")

      call ESMF_ConfigGetAttribute(LDT_config, landcover%fillvalue, &
           label="Landcover fill value:",rc=rc)
      call LDT_verify(rc,"Landcover fill value: option not specified in the config file")

    elseif( landcover%filltype == "none" ) then
      write(LDT_logunit,*) " -- 'NONE' Parameter-Mask Agreement Option Selected for Landcover"
    else
         write(LDT_logunit,*) "[ERR] Fill option for Landcover is not valid: ",trim(landcover%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none or neighbor "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
    end if

!-- Regional mask dataset file and option inputs:
    if( LDT_LSMparam_struc(1)%regmask%selectOpt > 0 ) then

       LDT_rc%cliplandmask = .false.
       call ESMF_ConfigFindLabel(LDT_config,"Clip landmask with regional mask:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%cliplandmask(n),rc=rc)
          call LDT_verify(rc,'Clip landmask with regional mask: not specified')
       enddo

       call ESMF_ConfigFindLabel(LDT_config,"Regional mask file:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%regfile(n),rc=rc)
          call LDT_verify(rc,'Regional mask file: not specified')
       enddo
       
       call ESMF_ConfigFindLabel(LDT_config,"Regional mask spatial transform:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%reg_gridtransform(n),rc=rc)
          call LDT_verify(rc,'Regional mask spatial tranform: option not specified in the config file')
       enddo

       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%reg_proj,&
            label="Regional mask map projection:",rc=rc)
       call LDT_verify(rc,'Regional mask map projection: option not specified in the config file')

       call LDT_readDomainConfigSpecs("Regional mask", LDT_rc%reg_proj, LDT_rc%reg_gridDesc)
    end if


!-- Allocate landcover and landmask file arrays:
    do n = 1, LDT_rc%nnest
       allocate(LDT_LSMparam_struc(n)%dommask%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            (LDT_LSMparam_struc(n)%landmask%num_bins)))

       allocate(LDT_LSMparam_struc(n)%landmask%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            (LDT_LSMparam_struc(n)%landmask%num_bins)))
       
       allocate(LDT_LSMparam_struc(n)%landmask2%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            LDT_LSMparam_struc(n)%landmask%num_bins))
       
       allocate(LDT_LSMparam_struc(n)%watermask%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            LDT_LSMparam_struc(n)%landmask%num_bins))
       allocate(LDT_LSMparam_struc(n)%luindex%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            1))       

     ! Fill in landmask placeholder parameter entries:
       call populate_param_attribs( "MASK2", &
            "Mask used for parameter consistency", "-",     &
            LDT_LSMparam_struc(n)%landmask, &
            LDT_LSMparam_struc(n)%landmask2 )
       
     ! Fill in water mask parameter entries:
       call populate_param_attribs( "WATERMASK", &
            "Water Mask used for water points", "-",     &
            LDT_LSMparam_struc(n)%landmask, &
            LDT_LSMparam_struc(n)%watermask )

       allocate(landcover_fgrd( LDT_rc%lnc(n),LDT_rc%lnr(n),&
            LDT_LSMparam_struc(n)%landcover%num_bins) )
       
       allocate(LDT_LSMparam_struc(n)%landcover%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            LDT_LSMparam_struc(n)%sfctype%vlevels))
       
       if( LDT_LSMparam_struc(n)%regmask%selectOpt == 1 ) then
          allocate(LDT_LSMparam_struc(n)%regmask%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%regmask%num_bins))
       endif

    !- Assign landmask source, in relation to original source or landcover:
       LDT_rc%mask_source(n) = LDT_LSMparam_struc(n)%landmask%source
       if( trim(LDT_rc%mask_type(n)) == "create" .or. &
            trim(LDT_LSMparam_struc(n)%landcover%source) == "ISA" ) then
          LDT_rc%mask_source(n) = LDT_LSMparam_struc(n)%landcover%source
       endif

    
!-- Perform checks on available landcover/landmask map input options: 

       if( LDT_LSMparam_struc(n)%regmask%selectOpt > 0 ) then
          call LDT_gridOptChecks( n, "Regional mask", &
               LDT_rc%reg_gridtransform(n), &
               LDT_rc%reg_proj, LDT_rc%reg_gridDesc(n,9) ) 
       endif

       call LDT_LMLCOptChecks( "Landcover", &
            LDT_rc%lc_type(n), LDT_rc%nt, LDT_rc%lc_gridtransform(n) )
       
       
!-- Read land cover and mask files:

       write(LDT_logunit,*) "[INFO] Reading landcover values"

    !- Landcover:
       call readlandcover(&
            trim(LDT_LSMparam_struc(n)%landcover%source)//char(0),&
            n, LDT_LSMparam_struc(n)%landcover%num_bins, &
            landcover_fgrd,                              &
            LDT_LSMparam_struc(n)%landmask%value  )
       
       LDT_LSMparam_struc(n)%landcover%vlevels = &
            LDT_LSMparam_struc(n)%sfctype%num_bins
       LDT_LSMparam_struc(n)%landcover%value = 0.
       LDT_LSMparam_struc(n)%landcover%value(:,:,1:LDT_LSMparam_struc(n)%landcover%num_bins) = &
            landcover_fgrd(:,:,1:LDT_LSMparam_struc(n)%landcover%num_bins)
       
    !- Copy readin/created landmask to place holder:
       LDT_LSMparam_struc(n)%landmask2%value(:,:,:) = &
            LDT_LSMparam_struc(n)%landmask%value(:,:,:)
       
    !- Reverse landmask to be a watermask for lake/open-water map fills:
       LDT_LSMparam_struc(n)%watermask%value(:,:,:) = 0.
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             if( LDT_LSMparam_struc(n)%landmask%value(c,r,1) == 1 ) then
                LDT_LSMparam_struc(n)%watermask%value(c,r,1) = 0.
             elseif( LDT_LSMparam_struc(n)%landmask%value(c,r,1) == 0 ) then
                LDT_LSMparam_struc(n)%watermask%value(c,r,1) = 1.
             endif
          enddo
       enddo

     ! Fill where parameter values are missing compared to land/water mask:
       if( landcover%filltype == "neighbor" ) then
          write(LDT_logunit,*) "[INFO] Checking/filling mask values for: ", &
                             trim(LDT_LSMparam_struc(n)%landcover%short_name)
          write(fill_logunit,*) "Checking/filling mask values for: ", &
                             trim(LDT_LSMparam_struc(n)%landcover%short_name)
          landcover%watervalue = LDT_rc%waterclass

          ! EMK...Make sure filling ignores water points.
          call LDT_discreteParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                LDT_rc%lc_gridtransform(n),                        &
                LDT_LSMparam_struc(n)%landcover%num_bins,          &
                LDT_LSMparam_struc(n)%landcover%value, landcover%watervalue, &
                LDT_LSMparam_struc(n)%landmask2%value,             &
                landcover%filltype, landcover%fillvalue, &
                landcover%fillradius , force_exclude_water=.true.)
       endif

    !- Surface type field update:
       call LDT_assign_lcsfctype( n, &
            LDT_LSMparam_struc(n)%sfctype%num_bins,   & 
            LDT_LSMparam_struc(n)%landcover%num_bins, &
            LDT_LSMparam_struc(n)%sfctype%value,      &
            LDT_LSMparam_struc(n)%landmask%value,     &
            LDT_LSMparam_struc(n)%landcover%value )
       write(LDT_logunit,*) "[INFO] Finished assigning LSM surface types"

    ! __________________________________________________

    ! _____________________LU_INDEX_____________________________
       nlctypes= LDT_LSMparam_struc(n)%landcover%num_bins  !size(landcover,3)
       !print*, "nlctypes:", nlctypes
       allocate(luindex(LDT_rc%lnc(n),LDT_rc%lnr(n)))
       allocate(landcover1(LDT_rc%lnc(n),LDT_rc%lnr(n),nlctypes))
       landcover1 = LDT_LSMparam_struc(n)%landcover%value

       do r = 1,LDT_rc%lnr(n)    ! nr
           do c = 1, LDT_rc%lnc(n) ! nc
              maxv = 0.0
              domv = -1
                do  t = 1, nlctypes
                  if (landcover1(c,r,t).gt.maxv) then
                    maxv = landcover1(c,r,t)
                    domv = t
                  endif
                enddo
                luindex(c,r) = domv

                if(domv.eq.-1.and.sum(landcover1(c,r,:)).eq.0) then 
                   luindex(c,r) = LDT_rc%waterclass
                   landcover1(c,r,LDT_rc%waterclass) = 1.0
                endif
             enddo
          enddo
          
          LDT_LSMparam_struc(n)%landcover%value = landcover1
          deallocate(landcover1)

        LDT_LSMparam_struc(n)%luindex%value(:,:,1) =luindex
        LDT_LSMparam_struc(n)%luindex%short_name    = "LU_INDEX"

    !-----------------------------------------------------


    !- Regional mask:
       if( LDT_LSMparam_struc(n)%regmask%selectOpt > 0 ) then
         call readregmask(&
              trim(LDT_LSMparam_struc(n)%regmask%source)//char(0),&
              n, LDT_LSMparam_struc(n)%regmask%value )

         LDT_LSMparam_struc(n)%regmask%vlevels = &
             LDT_LSMparam_struc(n)%regmask%num_bins

      !- Use the Regional land mask to "clip" the main landmask:
         if( LDT_rc%cliplandmask(n) .eqv. .true. ) then
           write(LDT_logunit,*) "[INFO] -- Using Regional Mask to Clip Landmask"
           do r = 1, LDT_rc%lnr(n)
             do c = 1, LDT_rc%lnc(n)
                if( LDT_LSMparam_struc(n)%regmask%value(c,r,1) == 0 ) then
                  LDT_LSMparam_struc(n)%landmask%value(c,r,1) = 0.
                endif
             enddo
           enddo
        endif

     endif
     LDT_LSMparam_struc(n)%dommask%value(:,:,:) = &
          LDT_LSMparam_struc(n)%landmask%value(:,:,:)

     deallocate( landcover_fgrd )
   enddo

   write(LDT_logunit,*) "[INFO] Finished reading landmask and landcover data"

  end subroutine LMLC_init_LISHydro

!BOP
! 
! !ROUTINE: LMLC_writeHeader_LIS
! \label{LMLC_writeHeader_LIS}
!
! !INTERFACE:
  subroutine LMLC_writeHeader_LIS(n,ftn,dimID)
! !USES: 
    use LDT_coreMod, only : LDT_rc
! 
! !DESCRIPTION: 
!
! This subroutine writes the NetCDF headers for the landcover
! and landmask fields in the standard preprocessing mode for LIS
!
!EOP

    integer      :: n 
    integer      :: ftn
    integer      :: dimID(3)
    integer      :: tdimID(3)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)
   
    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LDT_LSMparam_struc(n)%dommask)

    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LDT_LSMparam_struc(n)%landmask)

    if( LDT_LSMparam_struc(n)%regmask%selectOpt > 0 ) then
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%regmask)
    endif

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"LANDCOVER_SCHEME", &
         LDT_rc%lc_type(n)))

  ! Attributes serving Noah-MP only (at this time):
    if ((LDT_rc%lsm.eq."Noah-MP.3.6").or.(LDT_rc%lsm.eq."Noah-MP.4.0.1")) then
      select case( LDT_rc%lc_type(n) ) 
       case( "IGBPNCEP" ) 
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              20))
       case( "USGS" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              27))
       case( "USGS-RUC" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              28))
       case( "MODI-RUC" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              21))
       case( "NLCD40" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              40))
       case( "UMD" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              13))
       case( "NALCMS_SM" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              24))
       case( "NALCMS_SM_IGBPNCEP" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              20))
       case( "Bondville" ) 
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              20))
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"LANDCOVER_SCHEME", &
              "IGBPNCEP"))
      end select
    endif

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"BARESOILCLASS", &
         LDT_rc%bareclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"URBANCLASS", &
         LDT_rc%urbanclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SNOWCLASS", &
         LDT_rc%snowclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"WATERCLASS", &
         LDT_rc%waterclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"WETLANDCLASS", &
         LDT_rc%wetlandclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"GLACIERCLASS", &
         LDT_rc%glacierclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"CROPCLASS", &
         LDT_rc%cropclass1))

! - Enter number of vegetation land use types only:
!   (no wetland, snow/ice, water classes included at this time)
    select case( LDT_rc%lc_type(n) )
      case ( "UMD" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             13))
!             14))
      case ( "IGBP" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             14))
      case ( "JULES_PFT")
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             9))
      case ( "IGBPNCEP" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             17))
!             20))
      case ( "ECOCLIMAP2" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             12))
      case ( "USGS" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             22))
!             24))
      case ( "MOSAIC" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             6))
      case ( "ISA" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             13))
      case ( "CLM45" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             36))
      case ( "NALCMS_SM" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             24))
      case ( "NALCMS_SM_IGBPNCEP" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             17))
!             20))
      case ( "Bondville" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             17))
      case ( "CONSTANT" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             LDT_LSMparam_struc(n)%landcover%num_bins))
      case default
         write(LDT_logunit,*) ' Land cover scheme not currently recognized.' 
         write(LDT_logunit,*) ' Please select one of the following: ' 
         write(LDT_logunit,*) ' -- UMD ' 
         write(LDT_logunit,*) ' -- IGBP ' 
         write(LDT_logunit,*) ' -- IGBPNCEP ' 
         write(LDT_logunit,*) ' -- USGS ' 
         write(LDT_logunit,*) ' -- MOSAIC ' 
         write(LDT_logunit,*) ' -- ISA ' 
         write(LDT_logunit,*) ' -- CLM45 ' 
         write(LDT_logunit,*) ' -- NALCMS_SM '
         write(LDT_logunit,*) ' -- NALCMS_SM_IGBPNCEP '
         write(LDT_logunit,*) ' -- CONSTANT ' 
         write(LDT_logunit,*) ' ... program stopping. ' 
         call LDT_endrun
    end select

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"LANDMASK_SOURCE", &
         LDT_rc%mask_source(1)))

#endif
  end subroutine LMLC_writeHeader_LIS


!BOP
! !ROUTINE: LMLC_writeHeader_LISHydro
! \label{LMLC_writeHeader_LISHydro}
!
! !INTERFACE: 
 subroutine LMLC_writeHeader_LISHydro(n,ftn,dimID, flag)
! !USES:    
    use LDT_coreMod, only : LDT_rc
! 
! !DESCRIPTION: 
!
! This subroutine writes the NetCDF headers for the landcover
! and landmask fields in the LISHydro preprocessing mode
!
!EOP
    integer      :: n 
    integer      :: ftn
    integer      :: dimID(4)    !changed dimID(3) to dimID(4)
    integer      :: tdimID(4)
    integer      :: varid
    integer      :: luindexId
    integer      :: flag

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)
    tdimID(3) = dimID(4)   !assigned 'Time' dim here
   
       ! LU_INDEX field attributes: !!
    call LDT_verify(nf90_def_var(ftn, &
            "LU_INDEX", &
            nf90_float, (/dimID(1:2),dimID(4)/),LDT_LSMparam_struc(n)%luindexId), &
            'nf90_def_var failed for LDT_LSMparam_struc(n)%luindex')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%luindexId, &
            "standard_name","Dominant category"),&
            'nf90_put_att failed for luindex:standard_name')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%luindexId, &
            "units","-"),&
            'nf90_put_att failed for luindex:units')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%luindexId, &
            "scale_factor",1.0),&
            'nf90_put_att failed for luindex:scale_factor')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%luindexId, &
            "add_offset",0.0),&
            'nf90_put_att failed for luindex:add_offset')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%luindexId, &
            "vmin",0.0),&
            'nf90_put_att failed for luindex:vmin')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%luindexId, &
            "vmax",0.0),&
            'nf90_put_att failed for luindex:vmax')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%luindexId, &
            "num_bins",20),&
            'nf90_put_att failed for luindex:vmax')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%luindexId, &
            "stagger","M"),&
            'nf90_put_att failed for luindex:stagger')


    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LDT_LSMparam_struc(n)%dommask,flag)

    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LDT_LSMparam_struc(n)%landmask,flag)

    if( LDT_LSMparam_struc(n)%regmask%selectOpt > 0 ) then
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%regmask,flag)
    endif


    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"LANDCOVER_SCHEME", &
         LDT_rc%lc_type(n)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MMINLU", &
         "MODIFIED_IGBP_MODIS_NOAH"))

  ! Attributes serving Noah-MP only (at this time):
    if( LDT_rc%lsm == "Noah-MP.3.6" ) then
      select case( LDT_rc%lc_type(n) ) 
       case( "IGBPNCEP" ) 
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUM_LAND_CAT", &
              20))
       case( "USGS" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              27))
       case( "USGS-RUC" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              28))
       case( "MODI-RUC" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              21))
       case( "NLCD40" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              40))
       case( "UMD" )
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_LANDCATS", &
              13))
      end select
    endif

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"BARESOILCLASS", &
         LDT_rc%bareclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ISURBAN", &
         LDT_rc%urbanclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ISICE", &
         LDT_rc%snowclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ISWATER", &
         LDT_rc%waterclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"WETLANDCLASS", &
         LDT_rc%wetlandclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"GLACIERCLASS", &
         LDT_rc%glacierclass))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"grid_id", &
         LDT_rc%nnest))

! - Enter number of vegetation land use types only:
!   (no wetland, snow/ice, water classes included at this time)
    select case( LDT_rc%lc_type(n) )
      case ( "UMD" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             13))
!             14))
      case ( "IGBP" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             14))
      case ( "JULES_PFT")
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             9))
      case ( "IGBPNCEP" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             17))
!             20))
      case ( "ECOCLIMAP2" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             12))
      case ( "USGS" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             22))
!             24))
      case ( "MOSAIC" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             6))
      case ( "ISA" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             13))
      case ( "NALCMS_SM" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             24))
      case ( "NALCMS_SM_IGBPNCEP" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             17))
!             20))
      case ( "CONSTANT" )
        call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMVEGTYPES", &
             LDT_LSMparam_struc(n)%landcover%num_bins))
      case default
         write(LDT_logunit,*) ' Land cover scheme not currently recognized.' 
         write(LDT_logunit,*) ' Please select one of the following: ' 
         write(LDT_logunit,*) ' -- UMD ' 
         write(LDT_logunit,*) ' -- IGBP ' 
         write(LDT_logunit,*) ' -- IGBPNCEP ' 
         write(LDT_logunit,*) ' -- USGS ' 
         write(LDT_logunit,*) ' -- MOSAIC ' 
         write(LDT_logunit,*) ' -- ISA ' 
         write(LDT_logunit,*) ' -- NALCMS_SM '
         write(LDT_logunit,*) ' -- NALCMS_SM_IGBPNCEP '
         write(LDT_logunit,*) ' -- CONSTANT ' 
         write(LDT_logunit,*) ' ... program stopping. ' 
         call LDT_endrun
    end select

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"LANDMASK_SOURCE", &
         LDT_rc%mask_source(1)))

#endif
  end subroutine LMLC_writeHeader_LISHydro

!BOP
! 
! !ROUTINE: LMLC_writeData_LIS
! \label{LMLC_writeData_LIS}
! 
! !INTERFACE:
  subroutine LMLC_writeData_LIS(n,ftn)
! !USES:
    use LDT_coreMod
!
! !DESCRIPTION: 
! 
! This subroutine writes the landcover
! and landmask fields in NetCDF file in the 
! standard preprocessing mode for LIS
!
! 
!EOP
    integer  :: n 
    integer  :: ftn
    integer  :: ierr

!    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%landcover)

    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%dommask)

    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%landmask)
    
    if( LDT_LSMparam_struc(n)%regmask%selectOpt.eq.1 ) then
      call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%regmask)
    endif

  end subroutine LMLC_writeData_LIS

!BOP
! 
! !ROUTINE: LMLC_writeData_LISHydro
! \label{LMLC_writeData_LISHydro}
! 
! !INTERFACE:
  subroutine LMLC_writeData_LISHydro(n,ftn,flag)
! !USES: 
    use LDT_coreMod
!
! !DESCRIPTION: 
! 
! This subroutine writes the landcover
! and landmask fields in NetCDF file in the 
! standard preprocessing mode for LISHydro
!
! 
!EOP
    integer  :: n 
    integer  :: ftn
    integer  :: ierr
    integer  :: nc,nr
    integer  :: flag

    nr=size(LDT_LSMparam_struc(n)%luindex%value,1)
    nc=size(LDT_LSMparam_struc(n)%luindex%value,2)
    !print *, "nr and nc", nr, nc
    call LDT_verify(nf90_put_var(ftn,&
            LDT_LSMparam_struc(n)%luindexId, LDT_LSMparam_struc(n)%luindex%value(:,:,1),&
            (/1,1/),(/nr, nc/)),&
           'nf90_put_att failed for luindex')

!    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%landcover)

    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%dommask)

    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%landmask)
    
    if( LDT_LSMparam_struc(n)%regmask%selectOpt.eq.1 ) then
      call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%regmask)
    endif

  end subroutine LMLC_writeData_LISHydro
end module LDT_LMLCMod
