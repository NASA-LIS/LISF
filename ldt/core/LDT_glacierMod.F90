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
module LDT_glacierMod
!BOP
!
! !MODULE: LDT_glacierMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read various sources of
!  soil parameter data. 
! 
!  \subsubsection{Overview}
!  This module provides routines for reading and manipulating various 
!  parameters related to glaciers

! !REVISION HISTORY:
!
!  30 Mar 2018: Sujay Kumar; Initial implementation
!
  use ESMF
#if ( defined SPMD )
  use mpi
#endif
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_coreMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_SurfaceTypeMod
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_glacier_readParamSpecs
  public :: LDT_glacier_init  ! initializes data structures and memory
  public :: LDT_glacier_writeHeader
  public :: LDT_glacier_writeData
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_glacier_struc

  type, public :: glacier_type_dec
     character*50         :: mask_proj
     character*50         :: mask_gridtransform
     real                 :: mask_gridDesc(20)
  end type glacier_type_dec

  type(glacier_type_dec), allocatable :: LDT_glacier_struc(:)

!EOP

 contains

   subroutine LDT_glacier_readParamSpecs

     character*100    :: source
     integer          :: rc
     integer          :: n,k
     
     logical          :: check
     rc = 0

     allocate(LDT_glacier_struc(LDT_rc%nnest))
     
     check = .false. 

     do k = 1, LDT_rc%nsf_model_types
        if( LDT_rc%sf_model_type_name_select(k) == "Glacier" ) then
           check  = .true. 
        endif
     enddo
     if(check) then 
        call ESMF_ConfigFindLabel(LDT_config,"Glacier fraction cutoff value:",rc=rc)
        do n=1,LDT_rc%nnest
           call ESMF_ConfigGetAttribute(LDT_config,&
                LDT_rc%gridcell_glacier_frac(n),rc=rc)
           call LDT_verify(rc,'"Glacier fraction cutoff value: not specified')
        enddo
     endif

   end subroutine LDT_glacier_readParamSpecs

!BOP
! 
! !ROUTINE: LDT_glacier_init
! \label{LDT_glacier_init}
! 
! !INTERFACE:
  subroutine LDT_glacier_init()

! !USES:
    use LDT_coreMod,   only : LDT_rc, LDT_config
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_paramOptCheckMod, only: LDT_gridOptChecks!,&
!         LDT_glacierOptChecks

! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! glacier datasets and reads the glacier data based on the 
! choice of options specified in the ldt configuration. 
! 
!
!EOP
    implicit none
    integer  :: n, i, c, r,k
    integer  :: rc
    character*50        :: mask_proj
    real, allocatable   :: mask_gridDesc(:,:)
    type(LDT_fillopts)  :: glaciermask
    logical             :: check_data, check_modis_data
    integer             :: nl_start, nl_end
    real,    allocatable :: glacier_fgrd(:,:,:)
    integer, allocatable :: useExternalMask(:)
    character(len=LDT_CONST_PATH_LEN)        :: GLmaskname
! _____________________________________________________________________________
    

    write(LDT_logunit,*)" - - - - - - - - - - Glacier Parameters - - - - - - - - - - - - - -"

    allocate(mask_gridDesc(LDT_rc%nnest,20))

    do n=1, LDT_rc%nnest

       if(LDT_LSMparam_struc(n)%glaciermask%selectOpt.eq.1) then
          LDT_LSMparam_struc(n)%glaciermask%vlevels = &
               LDT_LSMparam_struc(n)%glaciermask%num_bins
          allocate(LDT_LSMparam_struc(n)%glaciermask%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%glaciermask%num_bins))
       endif
    enddo
    
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%glaciermask%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo

    check_modis_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%glaciermask%source.eq."MODIS_Native") then
          check_modis_data = .true. 
       endif
    enddo
    if(check_data) then 
       if(.not.check_modis_data) then 
          call ESMF_ConfigFindLabel(LDT_config,"Glacier mask map:", rc=rc)
          do i=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%glaciermask(i),rc=rc)
          enddo
          call ESMF_ConfigGetAttribute(LDT_config,mask_proj,&
               label="Glacier mask map projection:",rc=rc)
          call LDT_verify(rc,'Glacier mask map projection: option not specified in the config file')
          
          call ESMF_ConfigFindLabel(LDT_config,"Glacier mask spatial transform:",&
               rc=rc)
          do i=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,&
                  LDT_glacier_struc(i)%mask_gridtransform,&
                  label="Glacier mask spatial transform:",rc=rc)
             call LDT_verify(rc,'Glacier mask spatial transform: option not specified in the config file')
          enddo
          do n=1, LDT_rc%nnest
             
             if(LDT_LSMparam_struc(n)%glaciermask%selectOpt.eq.1) then
                call LDT_readDomainConfigSpecs("Glacier mask", &
                     mask_proj, mask_gridDesc)
                call LDT_gridOptChecks( n, "Glacier mask", &
                     LDT_glacier_struc(n)%mask_gridtransform, &
                  mask_proj, mask_gridDesc(n,9) )
!       call LDT_glacierOptChecks(n, "Glacier mask", &
!            mask_proj, LDT_glacier_struc(n)%mask_gridtransform )
             endif
             
             LDT_glacier_struc(n)%mask_proj = mask_proj
             LDT_glacier_struc(n)%mask_gridDesc(:) = mask_gridDesc(n,:)
          enddo
       endif
       do n=1,LDT_rc%nnest
          !if the source is MODIS_Native, then simply remap the existing
          !categories
          if(LDT_LSMparam_struc(n)%glaciermask%source.ne."MODIS_Native") then 
             
             call readglaciermask(&
                  trim(LDT_LSMparam_struc(n)%glaciermask%source)//char(0),&
                  n,LDT_LSMparam_struc(n)%glaciermask%value)
          else
             LDT_LSMparam_struc(n)%glaciermask%value = 0.0
          endif

       ! Fill where parameter values are missing compared to land/water mask:
!      glaciermask%filltype = "neighbor"
!      if(glaciermask%filltype == "neighbor" ) then
!         write(LDT_logunit,*) "Checking/filling mask values for: ", &
!              trim(LDT_LSMparam_struc(n)%glaciermask%short_name)
!         write(fill_logunit,*) "Checking/filling mask values for: ", &
!              trim(LDT_LSMparam_struc(n)%glaciermask%short_name)
!         glaciermask%watervalue = LDT_rc%udef
       !            glaciermask%filltype = "neighbor"
       !            glaciermask%fillvalue = 4.
       !            glaciermask%fillradius = 2.
!         call LDT_discreteParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
!              LDT_glacier_struc(n)%mask_gridtransform,         &
!              LDT_LSMparam_struc(n)%glaciermask%num_bins,                &
!              LDT_LSMparam_struc(n)%glaciermask%value, glaciermask%watervalue, &
!              LDT_LSMparam_struc(n)%landmask2%value,               &
!              glaciermask%filltype, glaciermask%fillvalue, &
!              glaciermask%fillradius )
!      endif
       enddo

!Now apply the additional mask if specified. 
       allocate(useExternalMask(LDT_rc%nnest))
       useExternalMask = 0
       call ESMF_ConfigFindLabel(LDT_config,"Apply additional external glacier mask:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,useExternalMask(i),rc=rc)
       enddo

       do n=1, LDT_rc%nnest
          if(useExternalMask(n).eq.1) then 
             
             call ESMF_ConfigFindLabel(LDT_config,"External glacier mask file:",rc=rc)
             call ESMF_ConfigGetAttribute(LDT_config,GLmaskname,rc=rc)
             call LDT_verify(rc,'External glacier mask file: not specified')

             call applyExtraGLmask(n,GLmaskname)

          endif
       enddo
       deallocate(useExternalMask)
       
       do k = 1, LDT_rc%nsf_model_types
          if( LDT_rc%sf_model_type_name_select(k) == "Glacier" ) then
             
         !- Assign Open water surface types:
             do n = 1,LDT_rc%nnest
                allocate(glacier_fgrd( LDT_rc%lnc(n),LDT_rc%lnr(n),&
                     LDT_LSMparam_struc(n)%sfctype%num_bins))
                glacier_fgrd = 0
! May need revision - will this work if multiple patches are turned on?
! land + lake + glacier
                nl_start = LDT_LSMparam_struc(n)%landcover%num_bins+1
                nl_end   = LDT_LSMparam_struc(n)%landcover%num_bins+&
                     LDT_LSMparam_struc(n)%glaciermask%num_bins
                
                glacier_fgrd(:,:,nl_start) = &
                     LDT_LSMparam_struc(n)%glaciermask%value(:,:,LDT_LSMparam_struc(n)%glaciermask%num_bins)
            !- Update surface type inforamtion:
                call LDT_assign_glaciersfctype( n, &
                     LDT_LSMparam_struc(n)%landcover%num_bins, &
                     LDT_LSMparam_struc(n)%sfctype%num_bins,   &
                     LDT_LSMparam_struc(n)%sfctype%value,      &
                     LDT_LSMparam_struc(n)%dommask%value,     &
                     LDT_LSMparam_struc(n)%landcover%value,    &
                     glacier_fgrd )
                
                write(LDT_logunit,*) "[INFO] Finished assigning glacier surface types"
             enddo
             
          endif
       enddo

    endif

 end subroutine LDT_glacier_init

 subroutine applyExtraGLmask(n,GLmaskname)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

   integer           :: n 
   character(len=*)  :: GLmaskname

   integer           :: ftn
   integer           :: varid
   integer           :: c,r
   real              :: mask_value(LDT_rc%lnc(n),LDT_rc%lnr(n))
   
! Currently the code uses the SWE fields from the LSM output 
! as way to determine Glacier locations

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
   call LDT_verify(nf90_open(path=GLmaskname,&
        mode=nf90_nowrite,ncid=ftn),&
        'nf90_open file failed in applyExtraGLmask')
   
   call LDT_verify(nf90_inq_varid(ftn,"SWE_tavg",varid),&
        'nf90_inq_varid failed in applyExtraGLmask')

   call LDT_verify(nf90_get_var(ftn,varid,mask_value),&
        'nf90_get_var failed in applyExtraGLmask')
   call LDT_verify(nf90_close(ftn))

   do r=1,LDT_rc%lnr(n)
      do c=1,LDT_rc%lnc(n)
         if(mask_value(c,r).gt.1990) then 
            LDT_LSMparam_struc(n)%glaciermask%value(c,r,1) = 1.0
         endif
      enddo
   enddo
            
#endif

 end subroutine applyExtraGLmask


 subroutine LDT_glacier_writeHeader(n,ftn,dimID)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer    :: n 
    integer    :: ftn
    integer    :: dimID(3)
    integer    :: tdimID(3)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)

    if(LDT_LSMparam_struc(n)%glaciermask%selectOpt.eq.1) then
       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_LSMparam_struc(n)%glaciermask)
    endif
#endif
  end subroutine LDT_glacier_writeHeader

  subroutine LDT_glacier_writeData(n,ftn)

    use LDT_coreMod, only : LDT_rc

    integer      :: n 
    integer      :: ftn
    integer      :: ierr

    if(LDT_LSMparam_struc(n)%glaciermask%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%glaciermask)
    endif

  end subroutine LDT_glacier_writeData
  
end module LDT_glacierMod
