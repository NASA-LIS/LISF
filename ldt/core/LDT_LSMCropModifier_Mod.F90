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
module LDT_LSMCropModifier_Mod
!BOP
!
! !MODULE: LDT_LSMCropModifier_Mod
!
! !DESCRIPTION:
!  This module is designed to provide options to modify land surface model
!   (LSM) vegetation and soil parameters, where needed, if crop model or
!   parameter data information is incorporated or "assimilated".
!
!  \subsubsection{Overview}
!   Options are provided here to modify land surface model
!   (LSM) vegetation and soil parameters, where needed, if crop model or
!   parameter data information is incorporated.
!
! !REVISION HISTORY:
!  14 Jan 2014: K. Arsenault: Added to modify LSM parameters with crop data 
! ___________________________________________________________

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_LSMCrop_readParamSpecs
  public :: LDT_LSMCropMod_init     ! allocates memory for required structures
  public :: LDT_LSMCropMod_writeHeader
  public :: LDT_LSMCropMod_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_LSMCrop_struc

  type, public :: LSMCrop_type_dec
     integer           :: selectOpt
     character*50      :: crop_gridtransform
     character*50      :: crop_proj
     character(len=LDT_CONST_PATH_LEN)     :: croptfile
   ! LSM/Crop-specific entries:
!     logical            :: assimcropinfo
     character*20      :: assign_cropvalue
     character*20      :: config_croptype
     character*20      :: crop_classification
     character(len=LDT_CONST_PATH_LEN)     :: croplib_dir

     type(LDT_paramEntry) :: croptype    ! Crop type land cover
  end type LSMCrop_type_dec

  type(LSMCrop_type_dec), allocatable :: LDT_LSMCrop_struc(:)

!BOP 
! 
! !ROUTINE: LDT_LSMCropMod_writeHeader 
! \label{LDT_LSMCropMod_writeHeader}
! 
! !INTERFACE:
  interface LDT_LSMCropMod_writeHeader
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure LDT_LSMCropMod_writeHeader_LIS
     module procedure LDT_LSMCropMod_writeHeader_LISHydro
! 
! !DESCRIPTION:
! This interface provides routines for writing NETCDF header both 
! in LIS preprocessing requirements as well as LISHydro(WRFHydro) 
! preprocessing requiremetns. A dummy argument call "flagX" was added 
! to overload the LISHydro procedue.
!EOP 
  end interface

contains

  subroutine LDT_LSMCrop_readParamSpecs
    
    character*100     :: source
    integer           :: rc
    integer           :: n
    
    allocate(LDT_LSMCrop_struc(LDT_rc%nnest))
    LDT_LSMCrop_struc(:)%selectOpt = 0

    call ESMF_ConfigFindLabel(LDT_config,"Crop type data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       if( LDT_rc%assimcropinfo(n) ) then
         if(rc.ne.0) then 
            LDT_LSMCrop_struc(n)%selectOpt = 0
            call LDT_warning(rc,"WARNING: Crop type data source: not defined")
            return
         else
            if(source.eq."none") then 
               LDT_LSMCrop_struc(n)%selectOpt = 0
            else
               LDT_LSMCrop_struc(n)%selectOpt = 1
            endif
         endif
       endif
       call LDT_set_param_attribs(rc,LDT_LSMCrop_struc(n)%croptype,&
            "CROPTYPE",source)
    enddo

  end subroutine LDT_LSMCrop_readParamSpecs

!BOP
! 
! !ROUTINE: LDT_LSMCropMod_init
! \label{LDT_LSMCropMod_init}
! 
! !INTERFACE:
  subroutine LDT_LSMCropMod_init

    ! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify

    !
    ! !DESCRIPTION:
    !
    ! Allocates memory for data structures for reading 
    ! the LSM parameters that are impacted by crop parameter presence.
    ! 
    !  The routines invoked are: 
    !  \begin{description}
    !   \item[LSMCropModsetup](\ref{LSMParmssetup}) \newline
    !    calls the registry to invoke the LSMCropMod setup methods. 
    !  \end{description}
    !
    !EOP
    implicit none

    integer  :: n
    integer  :: c,r,t,k,i,j
    integer  :: rc
    integer  :: crop_index
    real     :: landcover_fgrd
    real     :: crop_fgrd
    logical  :: crop_select
    logical  :: croptype_select
    character(20) :: croptemp
    type(LDT_fillopts) :: croptype

   ! _____________________________________________________

    if(LDT_LSMCrop_struc(1)%selectOpt.eq.1) then 
       crop_select = .false.
       croptype_select = .false.

       do n=1,LDT_rc%nnest
          if( LDT_rc%assimcropinfo(n) .eqv. .true. ) crop_select = .true.
       enddo

       !- Turn on Incorporating Crop Information:
       if( crop_select ) then
          write(LDT_logunit,*)" - - - - - - Updating LSM Parameters with Crop Information - - - - - - -"

          !     do n = 1, LDT_rc%nnest
          !        allocate(sosdata(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          !        allocate(lgpdata(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          !     enddo

      !-- Crop Parameters:
          allocate ( LDT_rc%numcrop(LDT_rc%nnest) )

       !- Crop type:
          do n=1,LDT_rc%nnest
             if(LDT_LSMCrop_struc(n)%croptype%selectOpt.eq.1) then
                croptype_select = .true.
             !- Allocate croptype values:
                allocate(LDT_LSMCrop_struc(n)%croptype%value(&
                     LDT_rc%lnc(n),LDT_rc%lnr(n),&
                     LDT_LSMCrop_struc(n)%croptype%num_bins))
                LDT_LSMCrop_struc(n)%croptype%value = -9999.   ! perhaps change to 0. later
             endif
             call setCropParmsFullnames(n,"croptype",LDT_LSMCrop_struc(n)%croptype%source)
          enddo   ! End nest loop

          if( croptype_select ) then
             call ESMF_ConfigFindLabel(LDT_config,"Crop type file:", rc=rc)
             do n=1,LDT_rc%nnest
                call ESMF_ConfigGetAttribute(LDT_config,LDT_LSMCrop_struc(n)%croptfile,rc=rc)
             enddo
             call ESMF_ConfigFindLabel(LDT_config,"Crop map spatial transform:",&
                  rc=rc)
             do n=1,LDT_rc%nnest
                call ESMF_ConfigGetAttribute(LDT_config,LDT_LSMCrop_struc(n)%crop_gridtransform,&
                     rc=rc)
                call LDT_verify(rc,'Crop map spatial transform: option not specified in the config file')
             enddo
          endif

          !- Crop classification:
          LDT_LSMCrop_struc(:)%crop_classification = "none"
          call ESMF_ConfigFindLabel(LDT_config,"Crop classification:",rc=rc)
          do n=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,LDT_LSMCrop_struc(n)%crop_classification,rc=rc)
             call LDT_verify(rc,'Crop classification: not specified')
          enddo

          !- Assign number of crop types based on classification selected:
          do n=1,LDT_rc%nnest
             select case( LDT_LSMCrop_struc(n)%crop_classification )
             case( "CROPMAP" )
                LDT_rc%numcrop(n) = 19
             case( "FAOSTAT01" )
                LDT_rc%numcrop(n) = 18
             case( "FAOSTAT05" )
                LDT_rc%numcrop(n) = 175
             case default
                write(LDT_logunit,*) "[ERR] THE CROP CLASSIFICATION TYPE, ",&
                     trim(LDT_LSMCrop_struc(n)%crop_classification),", IS NOT RECOGNIZED."
                write(LDT_logunit,*) "  Please enter one of the following options: "
                write(LDT_logunit,*) "   -- CROPMAP   "
                write(LDT_logunit,*) "   -- FAOSTAT01 "
                write(LDT_logunit,*) "   -- FAOSTAT05 "
                write(LDT_logunit,*) "  Stopping ..." 
                call LDT_endrun
             end select
             write(LDT_logunit,*) " -- Number of crop types for, ", &
                  trim(LDT_LSMCrop_struc(n)%crop_classification),", is :: ",LDT_rc%numcrop(n)
          enddo   ! End Nest Loop

          LDT_LSMCrop_struc(:)%croplib_dir = "none"
          call ESMF_ConfigFindLabel(LDT_config,"Crop library directory:",rc=rc)
          do n=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,LDT_LSMCrop_struc(n)%croplib_dir,rc=rc)
             call LDT_verify(rc,'Crop library directory: not specified')
          enddo
          LDT_LSMCrop_struc(:)%assign_cropvalue = "none"
          call ESMF_ConfigFindLabel(LDT_config,"Assign crop value type:",rc=rc)
          do n=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,LDT_LSMCrop_struc(n)%assign_cropvalue,rc=rc)
             call LDT_verify(rc,'Assign crop value type: not specified')
          enddo
          LDT_LSMCrop_struc(:)%config_croptype = "none"
          call ESMF_ConfigFindLabel(LDT_config,"Default crop type:",rc=rc)
          do n=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,&
                  LDT_LSMCrop_struc(n)%config_croptype,rc=rc)
             call LDT_verify(rc,'Default crop type: not specified')
          enddo

       !- Append crop information to landcover map:
          do n = 1, LDT_rc%nnest

          !- Determine which landcover classification "crop" index:
             select case ( LDT_rc%lc_type(n) )
             case ( "UMD" )      
                LDT_rc%cropclass1 = 11   ! One dominant crop type
                LDT_rc%cropclass2 = 11   ! Just assigned again here
                LDT_rc%grassclass = 10
             case ( "IGBP", "IGBPNCEP" ) ! 12 (primary), 14
                LDT_rc%cropclass1 = 12   ! Croplands
                LDT_rc%cropclass2 = 14   ! Cropland/natural vegetation mosaic
                LDT_rc%grassclass = 10   ! Grassland
             case ( "USGS" )     ! 2-6, primary (2)
                LDT_rc%cropclass1 = 2    ! Dryland cropland and pasture
                LDT_rc%cropclass2 = 3    ! Irrigated cropland and pasture
                LDT_rc%cropclass3 = 5    ! Cropland/grassland mosaic
                LDT_rc%cropclass4 = 6    ! Cropland/woodland mosaic
                LDT_rc%grassclass = 7    ! Grassland
             case ( "ISA" )      
                LDT_rc%cropclass1 = 11
                LDT_rc%cropclass2 = 11
             case ( "MOSAIC" )   ! 9, if following SiB2 classification
                LDT_rc%cropclass1 = 9
                LDT_rc%cropclass2 = 9
                !    case default ! non-supported options
             end select

          !- Read in Cropmap file or assign "Constant" crop information:
             if( croptype_select ) then

                if( LDT_LSMCrop_struc(n)%crop_gridtransform == "tile" ) then
                   LDT_LSMCrop_struc(n)%croptype%vlevels = &
                        LDT_LSMCrop_struc(n)%croptype%num_bins
                elseif( LDT_LSMCrop_struc(n)%crop_gridtransform == "mode" ) then
                   LDT_LSMCrop_struc(n)%croptype%vlevels = 1
                endif

             !- Read crop type file and "blend" with output landcover map:
                if( LDT_LSMCrop_struc(n)%croptype%source .ne. "CONSTANT" ) then

                   call readcroptype(&
                        trim(LDT_LSMCrop_struc(n)%croptype%source)//char(0),&
                        n, LDT_LSMCrop_struc(n)%croptype%num_bins, &
                        LDT_LSMCrop_struc(n)%croptype%value )

                 ! Ensure consistency of dominant crop map with landcover cropland:
                   croptemp = "maize"   
                   call assigncroptype( n, LDT_LSMCrop_struc(n)%crop_classification,   &
                        LDT_rc%numcrop(n), croptemp, crop_index ) ! default - most common crop
                !  like in Ozdogan etal (2010)
                   do r = 1, LDT_rc%lnr(n)
                      do c = 1, LDT_rc%lnc(n)
                       ! Landcover map crop(s') gridcell fractions:
                         if( LDT_rc%cropclass1 .ne. LDT_rc%cropclass2 ) then  ! Account for two croptypes
                           landcover_fgrd = &   
                               LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%cropclass1) &
                             + LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%cropclass2)
                         else   ! Account only for one croptype
                           landcover_fgrd = &   
                               LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%cropclass1) 
                         endif                       

                       ! Crop map gridcell fraction:
                         crop_fgrd = LDT_LSMCrop_struc(n)%croptype%value(c,r,1)

                       ! If landcover indicates crop, but cropmap does not: Assign to maize type
                         if( landcover_fgrd > 0. .and. crop_fgrd <= 0. ) then
                            LDT_LSMCrop_struc(n)%croptype%value(c,r,1) = float(crop_index)
                      
                       ! If landcover has no crop, but cropmap does: Assign as undefined
                         elseif( landcover_fgrd == 0. .and. crop_fgrd > 0. ) then
                            LDT_LSMCrop_struc(n)%croptype%value(c,r,1) = LDT_rc%udef
                         endif

                      enddo
                   enddo

             !- Assign crop type file:
                else
                   if( LDT_LSMCrop_struc(n)%assign_cropvalue == "single") then
                      call assigncroptype( n, LDT_LSMCrop_struc(n)%crop_classification,   &
                           LDT_rc%numcrop(n), LDT_LSMCrop_struc(n)%config_croptype, & 
                           crop_index )

                   !- Assign the crop layer for landcover-based crop class:
                      do r = 1, LDT_rc%lnr(n) 
                         do c = 1, LDT_rc%lnc(n)
!                            do t = 1, LDT_rc%nt
                               landcover_fgrd = &
                                    LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%cropclass1)
                               if( landcover_fgrd > 0. ) then
                                  LDT_LSMCrop_struc(n)%croptype%value(c,r,1) = crop_index
                               endif
!                            enddo
                         enddo
                      enddo

                   else
                      print *, " You've selected CONSTANT but need to select "
                      print*, "'single' crop value assignment !! "
                   endif
                endif  ! End Crop Type Map or Assignment Condition

             endif     ! End crop type logical selection
          end do   ! End nest loop

       end if  ! Ending the LSM-Crop modifier check
    end if

  end subroutine LDT_LSMCropMod_init


  subroutine LDT_LSMCropMod_writeHeader_LIS(n,ftn,dimID)

    use LDT_coreMod, only : LDT_rc

    integer      :: n
    integer      :: ftn
    integer      :: dimID(3)
    integer      :: tdimID(3)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   if( LDT_rc%assimcropinfo(n) ) then

     tdimID(1) = dimID(1)
     tdimID(2) = dimID(2)

     if( LDT_LSMCrop_struc(n)%croptype%selectOpt > 0 ) then
       call LDT_verify(nf90_def_dim(ftn,'croptypes',&
            LDT_LSMCrop_struc(n)%croptype%vlevels,tdimID(3)))
 
       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_LSMCrop_struc(n)%croptype)
     end if
 
     call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"CROPCLASS_SCHEME", &
          LDT_LSMCrop_struc(n)%crop_classification))
 
     call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"CROPCLASS_NUMBER", &
          LDT_rc%numcrop(n)))
   endif
#endif
   
 end subroutine LDT_LSMCropMod_writeHeader_LIS

  subroutine LDT_LSMCropMod_writeHeader_LISHydro(n,ftn,dimID,flag)

    use LDT_coreMod, only : LDT_rc

    integer      :: n
    integer      :: ftn
    integer      :: dimID(4)
    integer      :: tdimID(4)
    integer :: flag, flagn

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   if( LDT_rc%assimcropinfo(n) ) then

     tdimID(1) = dimID(1)
     tdimID(2) = dimID(2)
     tdimID(4) = dimID(4)

     if( LDT_LSMCrop_struc(n)%croptype%selectOpt > 0 ) then
       call LDT_verify(nf90_def_dim(ftn,'croptypes',&
            LDT_LSMCrop_struc(n)%croptype%vlevels,tdimID(3)))
 
       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_LSMCrop_struc(n)%croptype,flagn)
     end if
 
     call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"CROPCLASS_SCHEME", &
          LDT_LSMCrop_struc(n)%crop_classification))
 
     call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"CROPCLASS_NUMBER", &
          LDT_rc%numcrop(n)))
   endif
#endif
   
  end subroutine LDT_LSMCropMod_writeHeader_LISHydro


  subroutine LDT_LSMCropMod_writeData(n,ftn)

    use LDT_coreMod, only : LDT_rc

    integer  :: n
    integer  :: ftn

   if( LDT_rc%assimcropinfo(n) ) then

     if(LDT_LSMCrop_struc(n)%croptype%selectOpt.eq.1) then
        call LDT_writeNETCDFdata(n,ftn,LDT_LSMCrop_struc(n)%croptype)
     endif

   endif

  end subroutine LDT_LSMCropMod_writeData


end module LDT_LSMCropModifier_Mod

