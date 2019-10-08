!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
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

!EOP
contains

!BOP
! 
! !ROUTINE: LDT_LMLC_init
! \label{LDT_LMLC_init}
! 
! !INTERFACE:
  subroutine LDT_LMLC_init()

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
! landmask and landcover datasets
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

  end subroutine LDT_LMLC_init


  subroutine LDT_LMLC_writeHeader(n,ftn,dimID)
    
    use LDT_coreMod, only : LDT_rc

    integer      :: n 
    integer      :: ftn
    integer      :: dimID(3)
    integer      :: tdimID(3)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)
   
!    call LDT_verify(nf90_def_dim(ftn,'sfctypes',&
!         LDT_LSMparam_struc(n)%sfctype%vlevels,tdimID(3)))

!    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
!         LDT_LSMparam_struc(n)%landcover)

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
         write(LDT_logunit,*) ' -- CONSTANT ' 
         write(LDT_logunit,*) ' ... program stopping. ' 
         call LDT_endrun
    end select

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"LANDMASK_SOURCE", &
         LDT_rc%mask_source(1)))

#endif
  end subroutine LDT_LMLC_writeHeader

  subroutine LDT_LMLC_writeData(n,ftn)

    use LDT_coreMod

    integer  :: n 
    integer  :: ftn
    integer  :: ierr

!    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%landcover)

    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%dommask)

    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%landmask)
    
    if( LDT_LSMparam_struc(n)%regmask%selectOpt.eq.1 ) then
      call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%regmask)
    endif

  end subroutine LDT_LMLC_writeData

end module LDT_LMLCMod
