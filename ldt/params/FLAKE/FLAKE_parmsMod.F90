!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
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
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
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

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: FLAKE_struc

  type, public :: flake_type_dec

     character*50         :: lakeparms_proj
     character*50         :: lakeparms_gridtransform
     character*140        :: inlandwaterfile
     character*140        :: lakedepthfile
     character*140        :: lakedepthQCfile
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
      write(LDT_logunit,*)" - - - - - - - - - Lake Parameters - - - - - - - - -"
   endif

   do i=1,LDT_rc%nsf_model_types
      if( LDT_rc%sf_model_type_name_select(i) == "Lake" ) then
         lakesfctype_present = .true.
      endif
   enddo
   if( .not. lakesfctype_present .and. lakedepth_select ) then
      write(*,*) "[ERR] 'Lake' surface type is NOT selected but lake"
      write(*,*) "    depth field is.  To use the lake depth information,"
      write(*,*) "    increase the number of surface model types to at least"
      write(*,*) "    '2' and include 'Lake' in the 'Land surface model:' array."
      write(*,*) " Stopping ..."
      call LDT_endrun
   endif
   if( lakesfctype_present .and. .not. lakedepth_select ) then
      write(*,*) "[ERR] The 'Lake' surface type is selected but lake"
      write(*,*) "    depth field is NOT.  To use the lake depth information,"
      write(*,*) "    select 'Lake depth:' and point to an actual lake"
      write(*,*) "    depth file (e.g., FLake)."
      write(*,*) " Stopping ..."
      call LDT_endrun
   endif

   FLAKE_struc(:)%lakedepthfile   = "no_file"
   FLAKE_struc(:)%lakedepthQCfile = "no_file"

! ------- READ IN LDT.CONFIG FILE ENTRIES -------

!- Inland water classifications map: 
    if( inlandwatertype_select ) then

       call ESMF_ConfigFindLabel(LDT_config,"Inland waterbody data source:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               FLAKE_struc(n)%inlandwatertype%source,rc=rc)
          call LDT_verify(rc,'Inland waterbody data source: not specified')

          if( FLAKE_struc(n)%inlandwatertype%source == "GLWD" ) then
            FLAKE_struc(n)%inlandwatertype%num_bins = 12
          else
            print *, " I don't recognize that Inlandwater body data source:",&
                     trim(FLAKE_struc(n)%inlandwatertype%source)
            stop
          endif
       enddo


       call ESMF_ConfigFindLabel(LDT_config,"Inland waterbody type map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               FLAKE_struc(n)%inlandwaterfile,rc=rc)
          call LDT_verify(rc,'Inland waterbody type map: not specified')
       enddo

       call ESMF_ConfigFindLabel(LDT_config,"Inland waterbody spatial transform:",rc=rc)
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

  ! Outer edge temperature (K) of the thermally active layer of the bottom sediments
    if(lakesedtemp_select ) then
       call ESMF_ConfigFindLabel(LDT_config,&
            "Lake bottom sediments temperature value:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config, temp, rc=rc)
          FLAKE_struc(n)%lakesedtemp%value = temp
          call LDT_verify(rc,'Lake bottom sediments temperature value: not defined')
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
          call LDT_verify(rc,'Lake params spatial transform: option not specified in the config file')
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
         write(LDT_logunit,*)" Number of bins for lake depth QC field MUST MATCH "
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
                   if( LDT_rc%sf_model_type_name_select(k) == "Openwater" .and. &
                       LDT_LSMparam_struc(n)%sfctype%value(c,r,totaltypes) == 5 ) then
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

end module FLAKE_parmsMod
