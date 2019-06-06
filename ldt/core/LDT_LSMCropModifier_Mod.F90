!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
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
!  14 May 2019: K. Arsenault: Expanded crop support for tiling
! ____________________________________________________________________________

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
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
     character*140     :: croptfile
   ! LSM/Crop-specific entries:
     character*20      :: assign_cropvalue
     character*20      :: config_croptype
     character*100     :: croplib_dir

     type(LDT_paramEntry) :: croptype    ! Crop type land cover
  end type LSMCrop_type_dec

  type(LSMCrop_type_dec), allocatable :: LDT_LSMCrop_struc(:)

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
            call LDT_warning(rc,"WARNING: Crop type data source: not defined in config file")
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
    integer  :: c,r,t,k,i,i2,j
    integer  :: rc
    integer  :: crop_index
    real     :: landcover_fgrd
    real     :: crop_fgrd
    real     :: ratio
    real     :: vegtype_cropfrac1, vegtype_cropfrac2
    real     :: croptypes_sumfrac
    real     :: sumtemp, fractemp   ! FOR TESTING PURPOSES ONLY
    logical  :: crop_select
    logical  :: croptype_select
    character(20) :: croptemp
    integer  :: domcroptype
    real, allocatable    :: croparray(:)

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

       !- Crop type:
          do n=1,LDT_rc%nnest
             if(LDT_LSMCrop_struc(n)%croptype%selectOpt.eq.1) then
                croptype_select = .true.
             !- Allocate croptype values:
                LDT_LSMCrop_struc(n)%croptype%num_bins = LDT_rc%numcrop(n) 
                allocate(LDT_LSMCrop_struc(n)%croptype%value(&
                     LDT_rc%lnc(n),LDT_rc%lnr(n),&
                     LDT_LSMCrop_struc(n)%croptype%num_bins))
                LDT_LSMCrop_struc(n)%croptype%value = -9999.   ! perhaps change to 0. later
             endif
             ! Set crop type parameter derived type:
             call setCropParmsFullnames(n,"croptype",LDT_LSMCrop_struc(n)%croptype%source)
          enddo   ! End nest loop

          ! Read in croptype map config entries:
          if( croptype_select ) then
             call ESMF_ConfigFindLabel(LDT_config,"Crop type file:", rc=rc)
             do n=1,LDT_rc%nnest
                call ESMF_ConfigGetAttribute(LDT_config,LDT_LSMCrop_struc(n)%croptfile,rc=rc)
                call LDT_verify(rc,'Crop type file: option not specified in the config file')
             enddo
             call ESMF_ConfigFindLabel(LDT_config,"Crop map spatial transform:",&
                  rc=rc)
             do n=1,LDT_rc%nnest
                call ESMF_ConfigGetAttribute(LDT_config,LDT_LSMCrop_struc(n)%crop_gridtransform,&
                     rc=rc)
                call LDT_verify(rc,'Crop map spatial transform: option not specified in the config file')
             enddo
          endif

          ! Read in the other crop library and type config file entries:
          LDT_LSMCrop_struc(:)%croplib_dir = "none"
          call ESMF_ConfigFindLabel(LDT_config,"Crop library directory:",rc=rc)
          do n=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,LDT_LSMCrop_struc(n)%croplib_dir,rc=rc)
             call LDT_verify(rc,'Crop library directory: not specified in the config file')
          enddo
          LDT_LSMCrop_struc(:)%assign_cropvalue = "none"
          call ESMF_ConfigFindLabel(LDT_config,"Assign crop value type:",rc=rc)
          do n=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,LDT_LSMCrop_struc(n)%assign_cropvalue,rc=rc)
             call LDT_verify(rc,'Assign crop value type: not specified in the config file')
          enddo
          LDT_LSMCrop_struc(:)%config_croptype = "none"
          call ESMF_ConfigFindLabel(LDT_config,"Default crop type:",rc=rc)
          do n=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,&
                  LDT_LSMCrop_struc(n)%config_croptype,rc=rc)
             call LDT_verify(rc,'Default crop type: not specified in the config file')
          enddo

       !- Append crop information to landcover map:
          do n = 1, LDT_rc%nnest

             write(LDT_logunit,*)"[INFO] CURRENT LSM-Crop Tile Module accounts for"
             write(LDT_logunit,*)"[INFO]  single LANDCOVER crop type for now. Will"
             write(LDT_logunit,*)"[INFO]  expand to include other CROP/IRRIG classes ..."

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
                   print *, "[WARN] CURRENTLY CROPTILE IS NOT SUPPORTED IN LIS ..."
                   LDT_LSMCrop_struc(n)%croptype%vlevels = &
                        LDT_LSMCrop_struc(n)%croptype%num_bins
                elseif( LDT_LSMCrop_struc(n)%crop_gridtransform == "mode" ) then
                   print *, "[WARN] NEED TO TRANSFORM THE MODE TO WRITING OUT"
                   print *, "[WARN] DOMINANT CROP TYPE IN A CROP TILE ..."
                   LDT_LSMCrop_struc(n)%croptype%vlevels = 1
                elseif( LDT_LSMCrop_struc(n)%crop_gridtransform == "average" ) then
                   print *, "[WARN] CURRENTLY LIS IS SETUP FOR DOMINANT CROPTYPE"
                   LDT_LSMCrop_struc(n)%croptype%vlevels = 1
                endif

             !- Read crop type file and "blend" with output landcover map:
                if( LDT_LSMCrop_struc(n)%croptype%source .ne. "CONSTANT" ) then

                  call readcroptype(&
                        trim(LDT_LSMCrop_struc(n)%croptype%source)//char(0),&
                        n, LDT_LSMCrop_struc(n)%croptype%num_bins, &
                        LDT_LSMCrop_struc(n)%croptype%value )

                   ! CHECKS:
!                   print *, "LSMCrop module: ", LDT_LSMCrop_struc(n)%croptype%value(218,128,:)
!                   print *, "LSMCrop sum:    ", &
!                    sum(LDT_LSMCrop_struc(n)%croptype%value(218,128,:), &
!                     mask=LDT_LSMCrop_struc(n)%croptype%value(218,128,:).ne.LDT_rc%udef) 
!                   print *, "VEGETATION FRACS: ",LDT_LSMparam_struc(n)%landcover%value(218,128,:)
!                   print *, "VEG-CROP TYPE FRAC: ",LDT_LSMparam_struc(n)%landcover%value(218,128,LDT_rc%cropclass1)
!                   print *, "SURFACE TYPES: ",LDT_LSMparam_struc(n)%sfctype%value(218,128,:)
!                   print *, "SURFACE-CROP TYPE: ",LDT_LSMparam_struc(n)%sfctype%value(218,128,LDT_rc%cropclass1)

                  ! Perform user-options on croptype array and merge with vegtype:
                  do r = 1, LDT_rc%lnr(n)
                    do c = 1, LDT_rc%lnc(n)

                       ! Compare vegetation crop class fraction with sum of croptype fractions:
                       vegtype_cropfrac1 = LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%cropclass1)
                       croptypes_sumfrac = sum(LDT_LSMCrop_struc(n)%croptype%value(c,r,:), &
                                           mask=LDT_LSMCrop_struc(n)%croptype%value(c,r,:).ne.LDT_rc%udef)

                       ! Re-normalize the crop fraction totals (user option):
                       if( croptypes_sumfrac > 1.0 ) then
                         ratio = 1.0 / croptypes_sumfrac
                         do i = 1, LDT_rc%numcrop(n)
                            if( LDT_LSMCrop_struc(n)%croptype%value(c,r,i).ne.LDT_rc%udef ) then
                               LDT_LSMCrop_struc(n)%croptype%value(c,r,i) = &
                                   LDT_LSMCrop_struc(n)%croptype%value(c,r,i) * ratio
                            endif
                            croptypes_sumfrac = sum(LDT_LSMCrop_struc(n)%croptype%value(c,r,:), &
                                      mask=LDT_LSMCrop_struc(n)%croptype%value(c,r,:).ne.LDT_rc%udef)
                         enddo
                       endif

                       ! If croptypes not present, but vegetation class crop present:
                       if( croptypes_sumfrac == 0. .and. vegtype_cropfrac1 > 0. ) then

                          LDT_LSMCrop_struc(n)%croptype%value(c,r,1) = &    ! Replace with crop selected by user
                              LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%cropclass1)

                          LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%cropclass1) = 0.
                          LDT_LSMparam_struc(n)%sfctype%value(c,r,LDT_rc%cropclass1) = 0.

                       ! If croptypes present, but vegetation class crop not:
                       elseif( croptypes_sumfrac > 0. .and. vegtype_cropfrac1 == 0. ) then
                          LDT_LSMCrop_struc(n)%croptype%value(c,r,:) = 0.
     
                       ! Both croptypes and veg-crop class present:
                       elseif( croptypes_sumfrac > 0. .and. vegtype_cropfrac1 > 0. ) then

                         ratio = vegtype_cropfrac1 / croptypes_sumfrac 
!                         write(500,*) c, r, croptypes_sumfrac, vegtype_cropfrac1, ratio
                         do i = 1, LDT_rc%numcrop(n)                  
                            if( LDT_LSMCrop_struc(n)%croptype%value(c,r,i).ne.LDT_rc%udef ) then
                               fractemp = LDT_LSMCrop_struc(n)%croptype%value(c,r,i)
                               LDT_LSMCrop_struc(n)%croptype%value(c,r,i) = &
                                   LDT_LSMCrop_struc(n)%croptype%value(c,r,i) * ratio
!                               write(501,*) c, r, i, LDT_LSMCrop_struc(n)%croptype%value(c,r,i), fractemp
                            endif
                         enddo
                         LDT_LSMparam_struc(n)%landcover%value(c,r,LDT_rc%cropclass1) = 0.
                         LDT_LSMparam_struc(n)%sfctype%value(c,r,LDT_rc%cropclass1) = 0.

!                         sumtemp=sum(LDT_LSMCrop_struc(n)%croptype%value(c,r,:), &
!                              mask=LDT_LSMCrop_struc(n)%croptype%value(c,r,:).ne.LDT_rc%udef)
!                         if( sumtemp > 1.05 .or. sumtemp < 0.95 ) then
!                           write(502,*) c,r, sumtemp
!                         endif
                       endif

                       ! Double-check for CROPMAP classification with UMD vegtype:
                       if( LDT_rc%crop_classification(n) .eq. "CROPMAP" .and. &  
                            LDT_rc%lc_type(n) .ne. "UMD" ) then
                          write(LDT_logunit,*)"[ERR] CROPMAP Classification only works "
                          write(LDT_logunit,*)"  with the UMD landcover classification."
                          call LDT_endrun

                       ! Populate crop tiles within expanded landcover/surfacetype arrays:         
                       else
                         i2 = 0
                         do i = 1, LDT_rc%numcrop(n)
                           i2 = LDT_rc%nt + i    ! Add onto vegtypes

                           if( LDT_LSMCrop_struc(n)%croptype%value(c,r,i).ne.LDT_rc%udef ) then
                             LDT_LSMparam_struc(n)%landcover%value(c,r,i2) = &
                                 LDT_LSMCrop_struc(n)%croptype%value(c,r,i) 
                             ! Ensure updates to surfacetype array is updated
                             if( LDT_LSMCrop_struc(n)%croptype%value(c,r,i) > 0.0 ) then
                               LDT_LSMparam_struc(n)%sfctype%value(c,r,i2) = 1.0
                             endif
                           endif
!                           ! TEMPORARY CHCECK FOR FINAL LANDCOVER == 1.0
!                           croptypes_sumfrac = sum(LDT_LSMparam_struc(n)%landcover%value(c,r,:), &
!                                     mask=LDT_LSMparam_struc(n)%landcover%value(c,r,:).ne.LDT_rc%udef)
!                           if( croptypes_sumfrac > 1.05 .or. croptypes_sumfrac < 0.95 ) then
!                              write(503,*) c,r,croptypes_sumfrac
!                           endif
                         enddo
                       endif

                    end do
                  end do

!                  print *, "LSMCrop module(2): ", LDT_LSMCrop_struc(n)%croptype%value(218,128,:)
!                  print *, "LSMCrop sum(2):    ", &
!                   sum(LDT_LSMCrop_struc(n)%croptype%value(218,128,:), &
!                    mask=LDT_LSMCrop_struc(n)%croptype%value(218,128,:).ne.LDT_rc%udef)

                  ! Convert CROPTYPE to single dominant type for now, LIS - required!
                  !  This is for testing purposes; will be replaced with the vegtiles in LIS ...
                  allocate(croparray(LDT_LSMCrop_struc(n)%croptype%num_bins))
                  croparray = LDT_rc%udef
                  do r = 1, LDT_rc%lnr(n)
                    do c = 1, LDT_rc%lnc(n)

                      if( maxval(LDT_LSMCrop_struc(n)%croptype%value(c,r,:)) == LDT_rc%udef) then
                         LDT_LSMCrop_struc(n)%croptype%value(c,r,1) = LDT_rc%udef
                      else
                       croptypes_sumfrac = sum(LDT_LSMCrop_struc(n)%croptype%value(c,r,:), &
                                            mask=LDT_LSMCrop_struc(n)%croptype%value(c,r,:).ne.LDT_rc%udef)

                       if( croptypes_sumfrac == 0 ) then 
                          LDT_LSMCrop_struc(n)%croptype%value(c,r,1) = LDT_rc%udef
                       else
                         croparray(:) = LDT_LSMCrop_struc(n)%croptype%value(c,r,:)
                        ! For final CROPTYPE output field -- designate as dominant type 
                        !  for now -- as required by LIS at this time... 
                         if( LDT_rc%crop_classification(n) == "CROPMAP" ) then
                           domcroptype = &
                             maxloc(croparray(:),1,mask=croparray.ne.LDT_rc%udef) + 13.0
                         else
                           domcroptype = &
                             maxloc(croparray(:),1,mask=croparray.ne.LDT_rc%udef) + LDT_rc%nt
                           LDT_LSMCrop_struc(n)%croptype%value(c,r,1) = domcroptype
                         endif
                       endif

                      endif
                    end do  
                  end do


!!  OLDER VERSION OF THE CODE ... SUPPORTED SINGLE CROP TYPE ONLY - DOMINANT TYPE ...
#if 0
                 ! Ensure consistency of dominant crop map with landcover cropland:
                   croptemp = "maize"   
                   call assigncroptype( n, LDT_rc%crop_classification(n),   &
                        LDT_rc%numcrop(n), croptemp, crop_index ) ! default - most common crop

                   print *, "Crop_index:: ",crop_index
                   print *, " Crop class 1 and 2 :: ",LDT_rc%cropclass1, LDT_rc%cropclass2

                !  Like in Ozdogan et al (2010):
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
#endif
!! OLDER CODE, BASED ON CROPMAP - OZDOGAN ET AL(2010) RULES

             !- Assign single crop type value:
                else
                   if( LDT_LSMCrop_struc(n)%assign_cropvalue == "single" ) then
                      call assigncroptype( n, LDT_rc%crop_classification(n),   &
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
                      print *, "[WARN] You've selected CONSTANT, so you need to select "
                      print*,  "'single' crop value assignment for this case. "
                   endif
                endif  ! End Crop Type Map or Assignment Condition

             endif     ! End crop type logical selection
          end do       ! End nest loop

       end if   ! Ending the LSM-Crop modifier check
    end if

  end subroutine LDT_LSMCropMod_init


  subroutine LDT_LSMCropMod_writeHeader(n,ftn,dimID)

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
          LDT_rc%crop_classification(n)))
 
     call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"CROPCLASS_NUMBER", &
          LDT_rc%numcrop(n)))
   endif
#endif
   
  end subroutine LDT_LSMCropMod_writeHeader


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

