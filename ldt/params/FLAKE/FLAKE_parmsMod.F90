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
module FLAKE_parmsMod
!BOP
!
! !MODULE: FLAKE_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read lake fraction
!  and type data.
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the
!  lake depth data and generates lake surface fraction.
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  15 Apr 2013: Kristi Arsenault;  Modified for lake datasets.
!  21 Jul 2021: Eric Kemp; Moved read_FLake_lakedepth into module to better
!               handle optional parameter.
  !
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_SurfaceTypeMod, only: LDT_assign_lakesfctype
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: FLAKEparms_init    !allocates memory for required structures
  public :: FLAKEparms_writeHeader
  public :: FLAKEparms_writeData
  public :: read_FLake_lakedepth ! EMK

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: FLAKE_struc

  type, public :: flake_type_dec

     character*50         :: lakeparms_proj
     character*50         :: lakeparms_gridtransform
     character(len=LDT_CONST_PATH_LEN)        :: inlandwaterfile
     character(len=LDT_CONST_PATH_LEN)        :: lakedepthfile
     character(len=LDT_CONST_PATH_LEN)        :: lakedepthQCfile
     character*50         :: inlandwater_gridtransform

     ! -  Lake model-specific:
     type(LDT_paramEntry) :: inlandwatertype  ! Inland water body types (-)
     type(LDT_paramEntry) :: inlandwaterfrac  ! Inland water type fraction (-)
     type(LDT_paramEntry) :: lakedepth        ! Lake depth (m)
     type(LDT_paramEntry) :: lakefrac         ! Lake fraction (-)
     type(LDT_paramEntry) :: lakedepthQC      ! Lake depth quality flag (-)
     type(LDT_paramEntry) :: lakewindfetch    ! Lake wind fetch (m)
     type(LDT_paramEntry) :: lakeseddepth     ! Lake sediment depth (m)
     type(LDT_paramEntry) :: lakesedtemp      ! Lake sediment temperature (K)

  end type flake_type_dec

  type(flake_type_dec), allocatable :: FLAKE_struc(:)

contains

!BOP
!
! !ROUTINE: FLAKEparms_init
! \label{FLAKEparms_init}
!
! !INTERFACE:
  subroutine FLAKEparms_init

! !USES:
!   use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
!   use LDT_paramOptCheckMod, only:  FLAKEOptChecks
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading
! the lake fraction datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[lakesetup](\ref{lakesetup}) \newline
!    calls the registry to invoke the lake setup methods.
!  \end{description}
!
!EOP
   implicit none
   integer   :: n, k, c ,r, i
   integer   :: rc
   integer   :: nl_start, nl_end
   real      :: temp
   integer   :: numlaketypes, totaltypes
   real, allocatable :: lake_fgrd(:,:,:)
   logical   :: lakesfctype_present
   logical   :: lake_select, lakedepth_select, inlandwatertype_select
   logical   :: lakedepthqc_select, lakewindfetch_select
   logical   :: lakesedtemp_select, lakeseddepth_select

   type(LDT_fillopts) :: lakedepth
   type(LDT_fillopts) :: inlandwater
! _____________________________________________________________

   allocate(FLAKE_struc(LDT_rc%nnest))

   do n=1,LDT_rc%nnest
      ! - Lake model parameters:
      call set_param_attribs(FLAKE_struc(n)%inlandwatertype,"INLANDWATERTYPE")
      call set_param_attribs(FLAKE_struc(n)%lakedepth,"LAKEDEPTH")
      call set_param_attribs(FLAKE_struc(n)%lakedepthqc,"LAKEDEPTHQC")
      call set_param_attribs(FLAKE_struc(n)%lakewindfetch,"LAKEWINDFETCH")
      call set_param_attribs(FLAKE_struc(n)%lakeseddepth,"LAKESEDIMDEPTH")
      call set_param_attribs(FLAKE_struc(n)%lakesedtemp,"LAKESEDIMTEMP")
   enddo

   lakesfctype_present = .false.
   lake_select = .false.
   lakedepth_select= .false.
   lakedepthqc_select= .false.
   inlandwatertype_select = .false.
   lakewindfetch_select = .false.
   lakeseddepth_select = .false.
   lakesedtemp_select = .false.

   do n=1,LDT_rc%nnest
      if( FLAKE_struc(n)%lakeseddepth%selectOpt.gt.0 ) then
         lakeseddepth_select = .true.
      endif
   enddo
   do n=1,LDT_rc%nnest
      if(FLAKE_struc(n)%lakewindfetch%selectOpt.gt.0 ) then
         lakewindfetch_select = .true.
      endif
   enddo
   do n=1,LDT_rc%nnest
      if(FLAKE_struc(n)%inlandwatertype%selectOpt.gt.0 ) then
         inlandwatertype_select = .true.
      endif
   enddo
   do n=1,LDT_rc%nnest
      if( FLAKE_struc(n)%lakedepth%selectOpt.gt.0 ) then
         lakedepth_select = .true.
      endif
   enddo
   do n=1,LDT_rc%nnest
      if( FLAKE_struc(n)%lakedepthqc%selectOpt.gt.0 ) then
         lakedepthqc_select = .true.
      endif
   enddo
   do n=1,LDT_rc%nnest
      if(FLAKE_struc(n)%lakesedtemp%selectOpt.gt.0 ) then
         lakesedtemp_select = .true.
      endif
   enddo
   do n=1,LDT_rc%nnest
      if( FLAKE_struc(n)%inlandwatertype%selectOpt.gt.0 .or. &
          FLAKE_struc(n)%lakedepth%selectOpt.gt.0 ) then
         lake_select = .true.
      endif
   enddo

   if( lake_select ) then
      write(LDT_logunit,*) &
           " - - - - - - - - - Lake Parameters - - - - - - - - -"
   endif

   do i=1,LDT_rc%nsf_model_types
      if( LDT_rc%sf_model_type_name_select(i) == "Lake" ) then
         lakesfctype_present = .true.
      endif
   enddo
   if( .not. lakesfctype_present .and. lakedepth_select ) then
      write(LDT_logunit,*) "[ERR] 'Lake' surface type is NOT selected but lake"
      write(LDT_logunit,*) &
           "    depth field is.  To use the lake depth information,"
      write(LDT_logunit,*) &
           "    increase the number of surface model types to at least"
      write(LDT_logunit,*) &
           "    '2' and include 'Lake' in the 'Land surface model:' array."
      write(LDT_logunit,*) " Stopping ..."
      call LDT_endrun
   endif
   if( lakesfctype_present .and. .not. lakedepth_select ) then
      write(LDT_logunit,*) "[ERR] The 'Lake' surface type is selected but lake"
      write(LDT_logunit,*) &
           "    depth field is NOT.  To use the lake depth information,"
      write(LDT_logunit,*) &
           "    select 'Lake depth:' and point to an actual lake"
      write(LDT_logunit,*) "    depth file (e.g., FLake)."
      write(LDT_logunit,*) " Stopping ..."
      call LDT_endrun
   endif

   FLAKE_struc(:)%lakedepthfile   = "no_file"
   FLAKE_struc(:)%lakedepthQCfile = "no_file"

! ------- READ IN LDT.CONFIG FILE ENTRIES -------

!- Inland water classifications map:
    if( inlandwatertype_select ) then

       call ESMF_ConfigFindLabel(LDT_config,"Inland waterbody data source:", &
            rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               FLAKE_struc(n)%inlandwatertype%source,rc=rc)
          call LDT_verify(rc,'Inland waterbody data source: not specified')

          if( FLAKE_struc(n)%inlandwatertype%source == "GLWD" ) then
            FLAKE_struc(n)%inlandwatertype%num_bins = 12
          else
             write(LDT_logunit) &
                  " I don't recognize that Inlandwater body data source:",&
                     trim(FLAKE_struc(n)%inlandwatertype%source)
            call LDT_endrun()
          endif
       enddo


       call ESMF_ConfigFindLabel(LDT_config,"Inland waterbody type map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               FLAKE_struc(n)%inlandwaterfile,rc=rc)
          call LDT_verify(rc,'Inland waterbody type map: not specified')
       enddo

       call ESMF_ConfigFindLabel(LDT_config, &
            "Inland waterbody spatial transform:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               FLAKE_struc(n)%inlandwater_gridtransform,&
               rc=rc)
          call LDT_verify(rc,'Inland waterbody spatial transform: &
               option not specified in the config file')
       enddo
    endif

 !- Lake depth (in meters):
    if( lakedepth_select ) then
       call ESMF_ConfigFindLabel(LDT_config,"Lake depth map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               FLAKE_struc(n)%lakedepthfile,rc=rc)
          call LDT_verify(rc,'Lake depth map: not specified')
       enddo
    endif
    if( lakedepthqc_select ) then
       call ESMF_ConfigFindLabel(LDT_config,"Lake depth QC map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               FLAKE_struc(n)%lakedepthQCfile,rc=rc)
          call LDT_verify(rc,'Lake depth QC map: not specified')
       enddo
    end if

  ! Static FLake parameters:
    do n=1,LDT_rc%nnest
       if( FLAKE_struc(n)%lakewindfetch%selectOpt.gt.0 ) then
          FLAKE_struc(n)%lakewindfetch%vlevels = &
               FLAKE_struc(n)%lakewindfetch%num_times
          allocate(FLAKE_struc(n)%lakewindfetch%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               FLAKE_struc(n)%lakewindfetch%num_bins))
       endif
       if( FLAKE_struc(n)%lakeseddepth%selectOpt.gt.0 ) then
          FLAKE_struc(n)%lakeseddepth%vlevels = &
               FLAKE_struc(n)%lakeseddepth%num_times
          allocate(FLAKE_struc(n)%lakeseddepth%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               FLAKE_struc(n)%lakeseddepth%num_bins))
       endif
       if( FLAKE_struc(n)%lakesedtemp%selectOpt.gt.0 ) then
          FLAKE_struc(n)%lakesedtemp%vlevels = &
               FLAKE_struc(n)%lakesedtemp%num_times
          allocate(FLAKE_struc(n)%lakesedtemp%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               FLAKE_struc(n)%lakesedtemp%num_bins))
       endif
    enddo

  ! Typical wind fetch (m):
    if( lakewindfetch_select ) then
       call ESMF_ConfigFindLabel(LDT_config,&
            "Lake wind fetch value:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config, temp,rc=rc)
          FLAKE_struc(n)%lakewindfetch%value = temp
          call LDT_verify(rc,'Lake wind fetch value: not defined')
       enddo
    endif

  ! Thermally active layer depth of bottom sediments (m)
    if( lakeseddepth_select ) then
       call ESMF_ConfigFindLabel(LDT_config,&
            "Lake bottom sediments depth value:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config, temp, rc=rc)
          FLAKE_struc(n)%lakeseddepth%value = temp
          call LDT_verify(rc,'Lake bottom sediments depth value: not defined')
       enddo
    endif

    ! Outer edge temperature (K) of the thermally active layer of the bottom
    ! sediments
    if(lakesedtemp_select ) then
       call ESMF_ConfigFindLabel(LDT_config,&
            "Lake bottom sediments temperature value:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config, temp, rc=rc)
          FLAKE_struc(n)%lakesedtemp%value = temp
          call LDT_verify(rc, &
               'Lake bottom sediments temperature value: not defined')
       enddo
    endif

 !- Read lake parameter projection and grid info:
    if( lakedepth_select .or. &
        inlandwatertype_select ) then

       call ESMF_ConfigFindLabel(LDT_config,"Lake params spatial transform:",&
            rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config, &
               FLAKE_struc(n)%lakeparms_gridtransform,&
               rc=rc)
          call LDT_verify(rc, &
          'Lake params spatial transform: option not specified in config file')
       enddo

    endif


! ------- Setting Number of Bins/Tiles for Lake Types -------

   do n = 1, LDT_rc%nnest

   !== Inland water body types:
      if( FLAKE_struc(n)%inlandwatertype%selectOpt.gt.0 ) then

         allocate(FLAKE_struc(n)%inlandwatertype%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              FLAKE_struc(n)%inlandwatertype%num_bins))

         allocate(FLAKE_struc(n)%inlandwaterfrac%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              FLAKE_struc(n)%inlandwatertype%num_bins))

        ! Fill in derived lake fraction parameter entries:
        ! ( input_parmattribs -> output_parmattribs )
         call populate_param_attribs( "INLANDWATERFRAC", &
              "Inland water type fraction", "-",     &
              FLAKE_struc(n)%inlandwatertype, &
              FLAKE_struc(n)%inlandwaterfrac )
      endif

      !== Lake depth parameters:
      if( FLAKE_struc(n)%lakedepthQC%selectOpt == 1 .and. &
           FLAKE_struc(n)%lakedepth%selectOpt == 0 ) then
         write(LDT_logunit,*) "Lake depth option must be turned 'ON' for "
         write(LDT_logunit,*) " Lake depth QC option to be read and used. "
         write(LDT_logunit,*) "Program stopping ..."
         call LDT_endrun
      end if

      if( FLAKE_struc(n)%lakedepth%num_bins .ne. &
           FLAKE_struc(n)%lakedepthQC%num_bins ) then
         write(LDT_logunit,*) &
              " Number of bins for lake depth QC field MUST MATCH "
         write(LDT_logunit,*)"  the number of bins of the lake depth field. "
         write(LDT_logunit,*)"Program stopping ..."
         call LDT_endrun
      end if

      if( FLAKE_struc(n)%lakedepth%selectOpt.gt.0 ) then

         numlaketypes = FLAKE_struc(n)%lakedepth%num_bins

         FLAKE_struc(n)%lakedepth%vlevels = numlaketypes
!              FLAKE_struc(n)%lakedepth%num_bins

         allocate(FLAKE_struc(n)%lakedepth%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              FLAKE_struc(n)%lakedepth%num_bins))

         allocate(FLAKE_struc(n)%lakefrac%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              FLAKE_struc(n)%lakedepth%num_bins) )

        ! Fill in derived lake fraction parameter entries:
        ! ( input_parmattribs -> output_parmattribs )
         call populate_param_attribs( "LAKEFRAC", "Lake fraction", "-", &
              FLAKE_struc(n)%lakedepth,  &
              FLAKE_struc(n)%lakefrac )

         allocate(lake_fgrd( LDT_rc%lnc(n),LDT_rc%lnr(n),&
              LDT_LSMparam_struc(n)%sfctype%num_bins))

      endif

    ! Lake depth QC parameters:
      if( FLAKE_struc(n)%lakedepthQC%selectOpt.gt.0 ) then
         allocate(FLAKE_struc(n)%lakedepthQC%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              FLAKE_struc(n)%lakedepthQC%num_bins))
      endif


! ------- Read in Lake Parameters -------

    !- Read inland water body type map:
       if( FLAKE_struc(n)%inlandwatertype%selectOpt == 1 ) then
          write(LDT_logunit,*) "Reading inland water body type map: "&
                               //trim(FLAKE_struc(n)%inlandwaterfile)
          call read_GLWD_inlandwatertype( &
               n, FLAKE_struc(n)%inlandwatertype%num_bins,  &
               FLAKE_struc(n)%inlandwatertype%value, &
               FLAKE_struc(n)%inlandwaterfrac%value  )

       endif  ! End lake map-based fraction calculation

    !- Read Lake depth (and QC) map:
       if( FLAKE_struc(n)%lakedepth%selectOpt == 1 ) then
          write(LDT_logunit,*) "Reading lake depth: "//&
               trim(FLAKE_struc(n)%lakedepthfile)

       !- Reading Lake Depth with QC File:
          if( FLAKE_struc(n)%lakedepthQC%selectOpt == 1 ) then

            write(LDT_logunit,*) "Reading lake depth QC data: "//&
                  trim(FLAKE_struc(n)%lakedepthQCfile)
            call read_FLake_lakedepth( &
               n, FLAKE_struc(n)%lakedepth%num_bins,  &
               FLAKE_struc(n)%lakedepth%value, &
               FLAKE_struc(n)%lakefrac%value,  &
               FLAKE_struc(n)%lakedepthQC%value )

       !- Reading Lake Depth with NO QC File:
          else
             call read_FLake_lakedepth( &
                    n, FLAKE_struc(n)%lakedepth%num_bins, &
                    FLAKE_struc(n)%lakedepth%value,  &
                    FLAKE_struc(n)%lakefrac%value  )
          endif


        ! ------- Add to Surface Type Bins/Tiles -------

          nl_start = LDT_LSMparam_struc(n)%landcover%num_bins+1
          nl_end   = LDT_LSMparam_struc(n)%landcover%num_bins+&
                     FLAKE_struc(n)%lakefrac%num_bins

          FLAKE_struc(n)%lakefrac%vlevels = &
                 FLAKE_struc(n)%lakefrac%num_bins

          lake_fgrd(:,:,nl_start) = &
            FLAKE_struc(n)%lakefrac%value(:,:,FLAKE_struc(n)%lakefrac%num_bins)

       !- Update surface type inforamtion:
          call LDT_assign_lakesfctype( n, &
                  LDT_LSMparam_struc(n)%landcover%num_bins, &
                  LDT_LSMparam_struc(n)%sfctype%num_bins,   &
                  LDT_LSMparam_struc(n)%sfctype%value,      &
                  LDT_LSMparam_struc(n)%dommask%value,     &
                  LDT_LSMparam_struc(n)%landcover%value,    &
                  lake_fgrd )

          write(LDT_logunit,*) "Finished assigning lake surface types"


       !- Set general mask values to 1 wherever valid lake points:
          totaltypes = LDT_LSMparam_struc(n)%sfctype%num_bins
          do r = 1, LDT_rc%lnr(n)
             do c = 1, LDT_rc%lnc(n)
                if( FLAKE_struc(n)%lakedepth%value(c,r,1) > 0. )then
                   LDT_LSMparam_struc(n)%dommask%value(c,r,1) = 1.
                endif

                do k = 1, LDT_rc%nsf_model_types
                   if( LDT_rc%sf_model_type_name_select(k) == "Openwater" &
                        .and. &
                        LDT_LSMparam_struc(n)%sfctype%value(c,r,totaltypes) &
                        == 5 ) then
                      LDT_LSMparam_struc(n)%dommask%value(c,r,1) = 1.
                   endif
                enddo

             enddo
          enddo

          deallocate(lake_fgrd)

       endif  ! End lakedepth-based fraction calculation

    enddo

  end subroutine FLAKEparms_init


  subroutine FLAKEparms_writeHeader(n,ftn,dimID,monthId)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer   :: n
    integer   :: ftn
    integer   :: dimID(3)
    integer   :: monthId
    integer   :: tdimID(3)

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    if( FLAKE_struc(n)%inlandwatertype%selectOpt.gt.0 ) then

       FLAKE_struc(n)%inlandwatertype%vlevels = &
           FLAKE_struc(n)%inlandwatertype%num_bins
       FLAKE_struc(n)%inlandwaterfrac%vlevels = &
           FLAKE_struc(n)%inlandwatertype%num_bins

       call LDT_verify(nf90_def_dim(ftn,'inlandwatertypes',&
                FLAKE_struc(n)%inlandwatertype%num_bins,tdimID(3)))

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
                FLAKE_struc(n)%inlandwatertype)

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
                FLAKE_struc(n)%inlandwaterfrac)

    endif

    if( FLAKE_struc(n)%lakedepth%selectOpt.gt.0 ) then

       FLAKE_struc(n)%lakedepth%vlevels = &
           FLAKE_struc(n)%lakedepth%num_bins
!       FLAKE_struc(n)%lakefrac%vlevels = &
!           LDT_LSMparam_struc(n)%sfctype%num_bins

!       call LDT_verify(nf90_def_dim(ftn,'laketiles',&
!            LDT_LSMparam_struc(n)%sfctype%num_bins,tdimID(3)))

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            FLAKE_struc(n)%lakefrac)

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
               FLAKE_struc(n)%lakedepth)

       if( FLAKE_struc(n)%lakedepthQC%selectOpt.gt.0 ) &
          call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
               FLAKE_struc(n)%lakedepthQC)
    endif

    if( FLAKE_struc(n)%lakewindfetch%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            FLAKE_struc(n)%lakewindfetch)

    if( FLAKE_struc(n)%lakeseddepth%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            FLAKE_struc(n)%lakeseddepth)

    if( FLAKE_struc(n)%lakesedtemp%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            FLAKE_struc(n)%lakesedtemp)

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMLAKES", &
       FLAKE_struc(n)%lakedepth%num_bins ))

#endif
  end subroutine FLAKEparms_writeHeader

  subroutine FLAKEparms_writeData(n,ftn)

    integer     :: n
    integer     :: ftn

    if( FLAKE_struc(n)%inlandwatertype%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdata(n,ftn,FLAKE_struc(n)%inlandwatertype)

    if( FLAKE_struc(n)%inlandwatertype%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdata(n,ftn,FLAKE_struc(n)%inlandwaterfrac)

    if( FLAKE_struc(n)%lakedepth%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdata(n,ftn,FLAKE_struc(n)%lakefrac)

    if( FLAKE_struc(n)%lakedepth%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdata(n,ftn,FLAKE_struc(n)%lakedepth)

    if( FLAKE_struc(n)%lakedepthQC%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdata(n,ftn,FLAKE_struc(n)%lakedepthQC)

    if( FLAKE_struc(n)%lakewindfetch%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdata(n,ftn,FLAKE_struc(n)%lakewindfetch)

    if( FLAKE_struc(n)%lakeseddepth%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdata(n,ftn,FLAKE_struc(n)%lakeseddepth)

    if( FLAKE_struc(n)%lakesedtemp%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdata(n,ftn,FLAKE_struc(n)%lakesedtemp)

  end subroutine FLAKEparms_writeData

!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name)

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name

! ____________________________________________________


   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = 1
   paramEntry%selectOpt = 1
   paramEntry%source = "FLAKE"
   paramEntry%units ="none"
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(short_name)

  end subroutine set_param_attribs

!BOP
!
! !ROUTINE: read_FLake_lakedepth
!  \label{read_FLake_lakedepth}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  16 Apr 2013: K. Arsenault; Modified for lake maps
!  21 Jul 2021: Eric Kemp; Moved into FLAKE_parmsMod to better handle
!               optional argument.
!
! !INTERFACE:
  subroutine read_FLake_lakedepth(n, num_bins, lakedepth, &
       lakefrac, lakedepthQC )

! !USES:
    use LDT_coreMod,       only : LDT_rc, LDT_domain
    use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
         LDT_releaseUnitNumber, LDT_endrun
    use LDT_gridmappingMod
    use LDT_fileIOMod,     only : LDT_transform_paramgrid
    use LDT_paramTileInputMod, only: param_1dbin_areacalc, &
         param_index_fgrdcalc

!EOP
    implicit none

! !ARGUMENTS:
    integer, intent(in) :: n
    integer, intent(in) :: num_bins
    real, intent(inout) :: lakedepth(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
    real, intent(inout) :: lakefrac(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
    real, intent(inout), optional :: &
         lakedepthQC(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
!
! !DESCRIPTION:
!  This subroutine retrieves the lake depth for each gridcell
!   and returns the values in a latlon projection.
!
!  The Global database provides the external parameters fields for the
!  parameterization of lakes in atmospheric modeling.  It combines depth
!  information for the individual lakes from different sources with a map.
!  For mapping, the raster map of ECOCLIMAP2 dataset for ecosystems was used.
!  For some large lakes the bathymetry is included.  Additionally, the software
!  to project the lake-related information accurately onto an atmospheric
!  model grid is provided.
!
!  The global gridded datasets (in GlobalLake.tar.gz) contain the following
!  information on the geographical grid with the resolution of 30 arc sec.
!  (approx. 1 km):
!
!   1) the distributed mean lake depth values OR bathymetry data, and
!   2) the distributed flagm (QC) map to estimate reliability of the lake depth
!     information in every pixel of the grid:
!
!    = 0 - no inland water,
!    = 1 - the lake was not recognized by the automating mapping software,
!    = 2 - the lake depth value was missing in the dataset for individual lakes,
!    = 3 - the real depth value was used,
!    = 4 - a river.
!
!  * Note:  For flags 1 and 2, a default lake depth value of 10 m is assigned
! ** Note:  For flag 4, 3 m is set for lake depth.
!
!  Ref:  Kourzeneva, E., E. Martin, Y. Batrak and P.L.E. Moigne, 2012:
!        Climate data for parameterisation of lakes in Numerical Weather
!        Prediction models, Tellus A 2012, 64: 17226.
!        DOI: 10.3402/tellusa.v64i0.17226
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of bins (or bands)
!  \item[lakedepth]
!   output field FLake lake depth
!  \item[lakefrac]
!   output grid fractions for lake tiles
!  \item[lakedepthQC]
!   optional output for lake depth QC output
!  \end{description}
!
!EOP
    integer :: ftn1, ftn2
    logical :: file_exists
    logical :: qc_file_present   ! Flag to indicate that QC file is present
    integer :: c, r, t, i, ilon, ilat, ilon2, ilat2 ! loop indexes
    integer :: NorthPix, SouthPix, WestPix, EastPix ! The coordinates in pixels
    integer :: ErCode

    integer, parameter :: NLonB = 43200, & ! No. of longitude pixels of bitmap
         NLatB = 21600    ! Number of latitude pixels of bitmap
    integer, parameter :: PixSize = 30   ! Pixel size of the bitmap in arcsecs
    integer(1), dimension(NLonB) :: LonPix1 ! Status data pixels along latitude
    integer(2), dimension(NLonB) :: LonPix2 ! Depth data pixels along latitude
    real    :: IN_yres, IN_xres

    integer :: glpnc, glpnr          ! Global total columns and rows
    integer :: subpnc, subpnr        ! Parameter subsetted columns and rows
    real    :: param_gridDesc(20)    ! Input parameter grid desc array
    real    :: subparam_gridDesc(20) ! Subsetted Input parameter grid desc
    integer, allocatable  :: lat_line(:,:), lon_line(:,:)
    integer :: mi           ! Total number of input param grid array points
    integer :: mo           ! Total number of output LIS grid array points
    integer, allocatable  :: n11(:) ! array that maps the location of each
                                    ! input grid point in the output grid.
    real,    allocatable  :: gi1(:), gi2(:) ! input parameter 1d grid
    logical*1, allocatable :: li1(:), li2(:) ! input logical mask (to match gi)
    real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output lis 1d grid
    real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output lis 1d grid
    real      :: go3(LDT_rc%lnc(n)*LDT_rc%lnr(n),1)
    logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output logical mask
                                                   !  (matching go1)
    logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output logical mask
                                                   !  (matching go2)
    !  real      :: lo3(LDT_rc%lnc(n)*LDT_rc%lnr(n),1)

    real, dimension(:,:), allocatable       :: ReadDepth
    ! Depth of lake(s) on the bitmap in the circumscribed rectangular
    integer(1), dimension(:,:), allocatable :: ReadStatus
    ! Status of lake(s) on the bitmap in the circumscribed rectangular

    real  :: lakedepth2d(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real  :: lakedepthQC2d(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real  :: lakedepthcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
    real  :: total_cnt(LDT_rc%lnc(n)*LDT_rc%lnr(n))

    character(50) :: projection

! __________________________________________________________

    projection = "latlon"
    IN_yres = 1.0/120.0
    IN_xres = 1.0/120.0

!- Set parameter grid array inputs:
    param_gridDesc(1)  = 0.          ! Latlon
    param_gridDesc(2)  = NLonB       ! input_cols
    param_gridDesc(3)  = NlatB       ! input_rows
    param_gridDesc(4)  = -90.0  + (IN_yres/2) ! LL lat (-89.9960000S)
    param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon (-179.9960000W)
    param_gridDesc(6)  = 128
    param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat (89.99570000N)
    param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon (179.9960000W)
    param_gridDesc(9)  = IN_yres     ! dy: 0.0083333
    param_gridDesc(10) = IN_xres     ! dx: 0.0083333
    param_gridDesc(20) = 64

    inquire(file=trim(FLAKE_struc(n)%lakedepthfile), exist=file_exists)
    if (.not. file_exists) then
       write(LDT_logunit,*) "Lake depth map (", &
            trim(FLAKE_struc(n)%lakedepthfile),") not found."
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
    endif
    inquire(file=trim(FLAKE_struc(n)%lakedepthQCfile), exist=qc_file_present)
    if (.not. qc_file_present) then
       write(LDT_logunit,*) "Lake depth QC map (", &
            trim(FLAKE_struc(n)%lakedepthQCfile),") not present."
       write(LDT_logunit,*) " No QC applied to lake depth map ..."
    endif

    write(LDT_logunit,*) "[INFO] Reading FLake lake depth files"

    lakedepth = 10.  ! Default lake depth value for FLake model
    lakefrac  = 0.
    if ( qc_file_present )  lakedepthQC = 0

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

   !- Map Parameter Grid Info to LIS Target Grid/Projection Info --
    subparam_gridDesc = 0.
    call LDT_RunDomainPts( n, projection, param_gridDesc(:), &
         glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, &
         lon_line )

!  WestPix: CALL Coor2Num(West, 1, WestPix, ErCode)
    WestPix = INT((180.+MIN(subparam_gridDesc(5),179.999))*3600/PixSize)+1
!  EastPix: CALL Coor2Num(East, 1, EastPix, ErCode)
    EastPix = INT((180.+MIN(subparam_gridDesc(8),179.999))*3600/PixSize)+1
!  SouthPix: CALL Coor2Num(South, 2, SouthPix, ErCode)
    SouthPix=  &
         NLatB-(INT((90.+MIN(subparam_gridDesc(4),89.999))*3600/PixSize)+1)+1
!  NorthPix: CALL Coor2Num(North, 2, NorthPix, ErCode)
    NorthPix = &
         NLatB-(INT((90.+MIN(subparam_gridDesc(7),89.999))*3600/PixSize)+1)+1

    allocate (ReadDepth(subpnc,subpnr), stat=ErCode)
    IF (ErCode.NE.0) then
       write(LDT_logunit,*) "[ERR] Can't allocate the array <<ReadDepth>>"
       call LDT_endrun()
    end if
    allocate( ReadStatus(subpnc,subpnr), stat=ErCode )
    IF (ErCode.NE.0) then
       write(LDT_logunit,*) "[ERR] Can't allocate the array <<ReadStatus>>"
       call LDT_endrun()
    end if
    ReadDepth = LDT_rc%udef
    if ( qc_file_present )  ReadStatus = -9

! -------------------------------------------------------------------
!    READ IN LAKE DEPTH AND QC LAKE DEPTH FILES/INFO
! -------------------------------------------------------------------

    ftn1 = LDT_getNextUnitNumber()
    if ( qc_file_present )  ftn2 = LDT_getNextUnitNumber()

    open(ftn1, file=FLAKE_struc(n)%lakedepthfile, &
         form='unformatted', access='direct',&
         convert="little_endian", recl=NLonB*2)

    if ( qc_file_present ) then
       open(ftn2, file=FLAKE_struc(n)%lakedepthQCfile, &
            form='unformatted', access='direct',&
            convert="little_endian", recl=NLonB)
    endif

    ilat2 = 0
    do ilat = SouthPix, NorthPix, -1
       ilat2 = ilat2 + 1
       ! Read actual lake depth file:
       read(ftn1,REC=ilat) LonPix2
       ! Read lake depth QC-file:
       if( qc_file_present ) then
          read(ftn2,REC=ilat) LonPix1
       endif

       ilon2 = 0
       do ilon = WestPix, EastPix
          ilon2 = ilon2 + 1
          ! Read-in Lake Depth 2-D Array:
          ReadDepth(ilon2,ilat2)=LonPix2(ilon)/10.
          ! Read-in QC Lake Depth 2-D Array:
          if ( qc_file_present ) then
             ReadStatus(ilon2,ilat2)=LonPix1(ilon)
          endif
       end do
    end do

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

    mi = subpnc*subpnr
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    allocate( gi1(mi), li1(mi), gi2(mi), li2(mi), n11(mi) )
    gi1 = LDT_rc%udef;  gi2 = 0.   ! LDT_rc%udef
    li1 = .false.; li2 = .false.
    lo1 = .false.; lo2 = .false.;! lo3 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
    i = 0
    do r = 1, subpnr
       do c = 1, subpnc;  i = i + 1
          gi1(i) = ReadDepth(c,r)
          if ( gi1(i) .ne. LDT_rc%udef ) li1(i) = .true.

          if ( qc_file_present ) then
             gi2(i) = float(ReadStatus(c,r))
             if( gi2(i) .ne. 0. ) li2(i) = .true.
             !if( gi2(i).ne. LDT_rc%udef ) li2(i) = .true.
          endif
       enddo
    enddo
    deallocate( ReadDepth )
    if ( qc_file_present) deallocate( ReadStatus )

    !- Create mapping between parameter domain and LIS grid domain:
    call upscaleByAveraging_input( subparam_gridDesc, &
         LDT_rc%gridDesc(n,:), mi, mo, n11 )


!- Transform parameter grid to LIS run domain:
    select case ( FLAKE_struc(n)%lakeparms_gridtransform )

 !- Transforming 2-D lake depth field:
    !case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )
    case( "none", "neighbor", "average" )

   !- Transform parameter from original grid to LIS output grid:
       lo1 = .false.
       call LDT_transform_paramgrid(n, FLAKE_struc(n)%lakeparms_gridtransform,&
            subparam_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )

   !- Convert 1D to 2D grid arrays:
       lakedepth2d = LDT_rc%udef
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             lakedepth2d(c,r) = go1(i)
          enddo
       enddo

    !- Estimate lake fraction from lake depth map:
       go3 = 0.
       total_cnt = 0.
       do i = 1, mi
          if ( li1(i) ) then
             !- Count depths:
             if ( n11(i).ne.0 .and. gi1(i) > 0 ) then
                total_cnt(n11(i)) = total_cnt(n11(i)) + 1.
                go3(n11(i),1) = go3(n11(i),1) + 1.0
             endif
          endif
       enddo
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             if ( total_cnt(i) .ne. 0  ) then
                lakefrac(c,r,1) = go3(i,1)/total_cnt(i)
             else
                lakefrac(c,r,1) = 0.
             endif
          enddo
       enddo

 !- 3D Tile Case:
    case( "tile" )

       write(LDT_logunit,*)  &
            " [ERR] FLake Lake Depth 'tile' option is disabled... "
       write(LDT_logunit,*)  &
            " [ERR] Stopping ... (please select 'average' for the time being)"
       write(LDT_logunit,*)
       call LDT_endrun

   !- Create mapping between parameter domain and LIS grid domain:
       call param_1dbin_areacalc( n, num_bins, mi, mo, n11, &
            0., gi1, lakefrac, lakedepth )

 !- All Other Cases:
    case default
       write(LDT_logunit,*) " This lake depth spatial transform, ",&
            trim(FLAKE_struc(n)%lakeparms_gridtransform), &
            ", is not available at this time ..."
       write(LDT_logunit,*) " Program stopping ..."
       call LDT_endrun

    end select
    deallocate( gi1, li1 )

!- Apply QC file to 2-D lake depth map:
    if ( qc_file_present ) then

  !- Select primary QC value (per gridcell):
       call upscaleByMode( mi, mo, LDT_rc%udef, n11, li2, gi2, &
            lo2, go2 )

  !- Convert 1D to 2D grid arrays:
       lakedepthQC2d = LDT_rc%udef
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             if (go2(i) == LDT_rc%udef) go2(i) = 0.
             lakedepthQC2d(c,r) = go2(i)    ! Real lake depth used

        !- Using QC values can modify the final output lake depths:
             ! Lake not orig. recognized; filled
             if ( lakedepthQC2d(c,r) == 1 ) lakedepth2d(c,r) = 10.0
             ! Lake depth value missing; filled
             if ( lakedepthQC2d(c,r) == 2 ) lakedepth2d(c,r) = 10.0
             ! River points; filled
             if ( lakedepthQC2d(c,r) == 4 ) lakedepth2d(c,r) =  3.0
          enddo
       enddo
    endif
    deallocate( gi2, li2 )

!- Bring 2-D Array to 3-D lake depth tile space:
    if ( FLAKE_struc(n)%lakeparms_gridtransform == "none"     .or. &
         FLAKE_struc(n)%lakeparms_gridtransform == "neighbor" .or. &
         FLAKE_struc(n)%lakeparms_gridtransform == "average" ) then
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             !- Single lake layer, write to first bin:
             !            lakefrac(c,r,1)  = 1.0
             lakedepth(c,r,1) = lakedepth2d(c,r)
             if ( qc_file_present ) then
                lakedepthQC(c,r,1)= lakedepthQC2d(c,r)
             endif
          enddo
       enddo
    end if

    call LDT_releaseUnitNumber(ftn1)
    if ( qc_file_present ) then
       call LDT_releaseUnitNumber(ftn2)
    endif

    write(LDT_logunit,*) "[INFO] Done reading FLake lake depth files."

  end subroutine read_FLake_lakedepth

end module FLAKE_parmsMod
