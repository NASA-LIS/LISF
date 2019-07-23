!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_topoMod
!BOP
!
! !MODULE: LDT_topoMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read various sources of
!  topography data
! 
!  \subsubsection{Overview}
!  This module provides routines for reading and modifying topography data.
!
! !REVISION HISTORY:
!
!  18 Jul 2008: Sujay Kumar; Initial implementation
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_logMod
  use LDT_paramDataMod
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_topo_init  ! initializes data structures and memory
  public :: LDT_topo_writeHeader
  public :: LDT_topo_writeData
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------

!EOP

contains

!BOP
! 
! !ROUTINE: LDT_topo_init
! \label{LDT_topo_init}
! 
! !INTERFACE:
  subroutine LDT_topo_init()
! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_paramOptCheckMod, only: LDT_topoOptChecks, LDT_gridOptChecks

! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! topo datasets
!
!EOP
    implicit none
    integer :: n, i
    integer :: rc
    logical :: check_data, topo_select
    type(LDT_fillopts)  :: elev
    type(LDT_fillopts)  :: slope
    type(LDT_fillopts)  :: aspect
! ______________________________________________________________
    
    topo_select = .false. 

    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%elevation%selectOpt.eq.1 .or. &
            LDT_LSMparam_struc(n)%slope%selectOpt.eq.1   .or. &
            LDT_LSMparam_struc(n)%aspect%selectOpt.eq.1 ) then 
          topo_select = .true. 
       endif
    enddo
    if(topo_select) &
       write(LDT_logunit,*)" - - - - - - - - - Topographic Parameters - - - - - - - - - - - -"

    allocate(LDT_rc%topo_gridDesc(LDT_rc%nnest,20))

    call ESMF_ConfigFindLabel(LDT_config,"Elevation number of bands:",rc=rc)
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%elevation%selectOpt.eq.1) then
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_LSMparam_struc(n)%elevation%num_bins,rc=rc)
          call LDT_verify(rc,"Elevation number of bands: not defined")
       endif
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Slope number of bands:",rc=rc)
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%slope%selectOpt.eq.1) then
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_LSMparam_struc(n)%slope%num_bins,rc=rc)
          call LDT_verify(rc,"Slope number of bands: not defined")
       endif
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Aspect number of bands:",rc=rc)
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%aspect%selectOpt.eq.1) then
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_LSMparam_struc(n)%aspect%num_bins,rc=rc)
          call LDT_verify(rc,"Aspect number of bands: not defined")
       endif
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Curvature number of bands:",rc=rc)
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%curvature%selectOpt.eq.1) then
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_LSMparam_struc(n)%curvature%num_bins,rc=rc)
          call LDT_verify(rc,"Curvature number of bands: not defined")
       endif
    enddo

       
    do n=1, LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%elevation%selectOpt.eq.1) then
          allocate(LDT_LSMparam_struc(n)%elevation%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%elevation%num_bins))
          
          allocate(LDT_LSMparam_struc(n)%elevfgrd%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%elevation%num_bins))
       endif
       if(LDT_LSMparam_struc(n)%slope%selectOpt.eq.1) then
          allocate(LDT_LSMparam_struc(n)%slope%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%slope%num_bins))
          
          allocate(LDT_LSMparam_struc(n)%slopefgrd%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%slope%num_bins))
       endif
       if(LDT_LSMparam_struc(n)%aspect%selectOpt.eq.1) then
          allocate(LDT_LSMparam_struc(n)%aspect%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%aspect%num_bins))
          
          allocate(LDT_LSMparam_struc(n)%aspectfgrd%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%aspect%num_bins))

       endif
       if(LDT_LSMparam_struc(n)%curvature%selectOpt.eq.1) then
          allocate(LDT_LSMparam_struc(n)%curvature%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%curvature%num_bins))
       endif
    enddo

! - Elevation:
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%elevation%selectOpt.eq.1) then
          check_data = .true.
       endif
    enddo

    if(check_data) then 
      call ESMF_ConfigFindLabel(LDT_config,"Elevation map:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elevfile(i),rc=rc)
      enddo

      elev%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, elev%filltype, &
           label="Elevation fill option:",rc=rc)
      call LDT_verify(rc,"Elevation fill option: option not specified in the config file")

      if( elev%filltype == "neighbor" .or. elev%filltype == "average" ) then
         call ESMF_ConfigGetAttribute(LDT_config, elev%fillradius, &
              label="Elevation fill radius:",rc=rc)
         call LDT_verify(rc,"Elevation fill radius: option not specified in the config file")

         call ESMF_ConfigGetAttribute(LDT_config, elev%fillvalue, &
              label="Elevation fill value:",rc=rc)
         call LDT_verify(rc,"Elevation fill value: option not specified in the config file")
      elseif( elev%filltype == "none" ) then
         write(LDT_logunit,*) " -- 'NONE' Parameter-Mask Agreement Option Selected for Elevation"
      else
         write(LDT_logunit,*) "[ERR] Fill option for Elevation is not valid: ",trim(elev%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
      end if

    ! Set units and full names:
      do n=1,LDT_rc%nnest
         LDT_LSMparam_struc(n)%elevation%units="m"     
         call setTopoParmsFullnames( n, "elevation", &
                 LDT_LSMparam_struc(n)%elevation%source )
      enddo

    endif
    
! - Slope:
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%slope%selectOpt.eq.1) then
          check_data = .true.
       endif
    enddo
    if( check_data) then
      call ESMF_ConfigFindLabel(LDT_config,"Slope map:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%slfile(i),rc=rc)
      enddo

      slope%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, slope%filltype, &
           label="Slope fill option:",rc=rc)
      call LDT_verify(rc,"Slope fill option: option not specified in the config file")

      if( slope%filltype == "neighbor" .or. slope%filltype == "average" ) then
         call ESMF_ConfigGetAttribute(LDT_config, slope%fillradius, &
              label="Slope fill radius:",rc=rc)
         call LDT_verify(rc,"Slope fill radius: option not specified in the config file")

         call ESMF_ConfigGetAttribute(LDT_config, slope%fillvalue, &
              label="Slope fill value:",rc=rc)
         call LDT_verify(rc,"Slope fill value: option not specified in the config file")
      elseif( slope%filltype == "none" ) then
         write(LDT_logunit,*) " -- 'NONE' Parameter-Mask Agreement Option Selected for Slope"
      else
         write(LDT_logunit,*) "[ERR] Fill option for Slope is not valid: ",trim(slope%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
      end if

    ! Set units and full names:
      do n=1,LDT_rc%nnest
         LDT_LSMparam_struc(n)%slope%units="-"
         call setTopoParmsFullnames( n, "slope", &
                 LDT_LSMparam_struc(n)%slope%source )
      enddo

    endif
    
! - Aspect:
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%aspect%selectOpt.eq.1) then
          check_data = .true.
       endif
    enddo

    if( check_data ) then
      call ESMF_ConfigFindLabel(LDT_config,"Aspect map:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%aspfile(i),rc=rc)
      enddo

      aspect%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, aspect%filltype, &
           label="Aspect fill option:",rc=rc)
      call LDT_verify(rc,"Aspect fill option: option not specified in the config file")

      if( aspect%filltype == "neighbor" .or. aspect%filltype == "average" ) then
         call ESMF_ConfigGetAttribute(LDT_config, aspect%fillradius, &
              label="Aspect fill radius:",rc=rc)
         call LDT_verify(rc,"Aspect fill radius: option not specified in the config file")

         call ESMF_ConfigGetAttribute(LDT_config, aspect%fillvalue, &
              label="Aspect fill value:",rc=rc)
         call LDT_verify(rc,"Aspect fill value: option not specified in the config file")
      elseif( aspect%filltype == "none" ) then
         write(LDT_logunit,*) " -- 'NONE' Parameter-Mask Agreement Option Selected for Aspect"
      else
         write(LDT_logunit,*) "[ERR] Fill option for Aspect is not valid: ",trim(aspect%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
      end if
    ! Set units and full names:
      do n=1,LDT_rc%nnest
         LDT_LSMparam_struc(n)%aspect%units="-"
         call setTopoParmsFullnames( n, "aspect", &
                 LDT_LSMparam_struc(n)%aspect%source )
      enddo

    endif
    
! - Curvature:
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%curvature%selectOpt.eq.1) then
          check_data = .true.
       endif
    enddo
    if(check_data) then
      call ESMF_ConfigFindLabel(LDT_config,"Curvature map:",rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%curvfile(i),rc=rc)
      enddo
    endif
    
! - Read in geographic input data:
    if(topo_select) then 

      call ESMF_ConfigFindLabel(LDT_config,"Topography spatial transform:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%topo_gridtransform(n),&
              rc=rc)
         call LDT_verify(rc,'Topography spatial transform: option not specified in the config file')
    
       ! Check tiled settings:
         if( LDT_LSMparam_struc(n)%elevation%num_bins > 1 .or. &
             LDT_LSMparam_struc(n)%slope%num_bins     > 1 .or. &
             LDT_LSMparam_struc(n)%aspect%num_bins    > 1 ) then
            if( LDT_rc%topo_gridtransform(n) .ne. 'tile' ) then
               write(LDT_logunit,*) "[ERR]  Number of elevation, slope or aspect"
               write(LDT_logunit,*) "    bins (or tiles) > 1, but 'tile' spatial transform "
               write(LDT_logunit,*) "    not selected.  Program stopping ..."
               call LDT_endrun
            endif
         endif

       ! Don't need to read in "Native" grid extents/resolution, just for LIS inputs
         if( index(LDT_LSMparam_struc(n)%elevation%source,"Native").eq.0 .and. &
             index(LDT_LSMparam_struc(n)%elevation%source,"CONSTANT").eq.0 ) then
            call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%topo_proj,&
                 label="Topography map projection:",rc=rc)
            call LDT_verify(rc,'Topography map projection: option not specified in the config file')

            call LDT_readDomainConfigSpecs("Topography", LDT_rc%topo_proj, LDT_rc%topo_gridDesc)
            if( LDT_rc%topo_proj == "latlon" ) then
              call LDT_gridOptChecks( n, "Topography", LDT_rc%topo_gridtransform(n), &
                                      LDT_rc%topo_proj, LDT_rc%topo_gridDesc(n,9) )
            endif
         endif
      enddo
    endif

!-- Loop over each nest:
    do n = 1, LDT_rc%nnest

    !- Topography parameter input option checks:
       if( (LDT_LSMparam_struc(n)%elevation%selectOpt== 1 )  .or. &
           (LDT_LSMparam_struc(n)%slope%selectOpt    == 1 )  .or. & 
           (LDT_LSMparam_struc(n)%aspect%selectOpt   == 1 ) )  then

         call LDT_topoOptChecks( "Topography", LDT_rc%topo_proj, LDT_rc%topo_gridtransform(n) )
       endif

    !- Elevation Parameter (Tiling available):
       if( LDT_LSMparam_struc(n)%elevation%selectOpt.eq.1 ) then

        ! Fill in derived elevation fraction parameter entries:
        ! ( input_parmattribs -> output_parmattribs ) 
          call populate_param_attribs( "ELEVFGRD", &
                                  "Elevation Area Fraction", "-",  &
                                  LDT_LSMparam_struc(n)%elevation, &
                                  LDT_LSMparam_struc(n)%elevfgrd )

          call readelev(trim(LDT_LSMparam_struc(n)%elevation%source)//char(0),&
               n, LDT_LSMparam_struc(n)%elevation%num_bins, &
               LDT_LSMparam_struc(n)%elevfgrd%value,        &
               LDT_LSMparam_struc(n)%elevation%value)

       ! Fill where parameter values are missing compared to land/water mask:
         if( elev%filltype == "neighbor" .or. elev%filltype == "average" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%elevation%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%elevation%short_name)
            if( LDT_LSMparam_struc(n)%elevation%source == "SRTM" ) then
               elev%watervalue = 0.
            else  ! e.g., GTOPO30
               elev%watervalue = LDT_rc%udef
            endif
         !- Tiled output data fields:
            if( LDT_rc%topo_gridtransform(n) == "tile" ) then
               call LDT_contTileParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                  LDT_rc%topo_gridtransform(n),                         &
                  LDT_LSMparam_struc(n)%elevation%num_bins,             &
                  LDT_LSMparam_struc(n)%elevation%value,                &
                  LDT_LSMparam_struc(n)%elevfgrd%value, elev%watervalue,&
                  LDT_LSMparam_struc(n)%landmask2%value, elev%filltype, &
                  elev%fillvalue, elev%fillradius )
         !- Non-tiled output data fields:
            else
               call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                  LDT_rc%topo_gridtransform(n),                               &
                  LDT_LSMparam_struc(n)%elevation%num_bins,                &
                  LDT_LSMparam_struc(n)%elevation%value, elev%watervalue,  &
                  LDT_LSMparam_struc(n)%landmask2%value, elev%filltype,    &
                  elev%fillvalue, elev%fillradius )
            endif
         endif 
       endif

    !- Slope Parameter (Tiling available):
       if( LDT_LSMparam_struc(n)%slope%selectOpt.eq.1 ) then

        ! Fill in derived slope fraction parameter entries:
        ! ( input_parmattribs -> output_parmattribs ) 
          call populate_param_attribs( "SLOPEFGRD", &
                                  "Slope Area Fraction", "-",  &
                                  LDT_LSMparam_struc(n)%slope, &
                                  LDT_LSMparam_struc(n)%slopefgrd )

          call readslope( trim(LDT_LSMparam_struc(n)%slope%source)//char(0), &
               n, LDT_LSMparam_struc(n)%slope%num_bins, &
               LDT_LSMparam_struc(n)%slopefgrd%value,   &
               LDT_LSMparam_struc(n)%slope%value )

       ! Fill where parameter values are missing compared to land/water mask:
         if( slope%filltype == "neighbor" .or. slope%filltype == "average" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%slope%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%slope%short_name)

            if( LDT_LSMparam_struc(n)%slope%source == "SRTM" ) then
               slope%watervalue = 0.
            else  ! e.g., GTOPO30
               slope%watervalue = LDT_rc%udef
            endif
         !- Tiled output data fields:
            if( LDT_rc%topo_gridtransform(n) == "tile" ) then
               call LDT_contTileParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                  LDT_rc%topo_gridtransform(n),                              &
                  LDT_LSMparam_struc(n)%slope%num_bins,                   &
                  LDT_LSMparam_struc(n)%slope%value,                      &
                  LDT_LSMparam_struc(n)%slopefgrd%value, slope%watervalue,&
                  LDT_LSMparam_struc(n)%landmask2%value, slope%filltype,  &
                  slope%fillvalue, slope%fillradius )
         !- Non-tiled output data fields:
            else
               call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                  LDT_rc%topo_gridtransform(n),                              &
                  LDT_LSMparam_struc(n)%slope%num_bins,                   &
                  LDT_LSMparam_struc(n)%slope%value, slope%watervalue,    &
                  LDT_LSMparam_struc(n)%landmask2%value, slope%filltype,  &
                  slope%fillvalue, slope%fillradius )
            endif
         endif
       endif

    !- Aspect Parameter (Tiling available):
       if( LDT_LSMparam_struc(n)%aspect%selectOpt.eq.1 ) then

        ! Fill in derived elevation fraction parameter entries:
        ! ( input_parmattribs -> output_parmattribs ) 
          call populate_param_attribs( "ASPECTFGRD", &
                                  "Aspect Area Fraction", "-",  &
                                  LDT_LSMparam_struc(n)%aspect, &
                                  LDT_LSMparam_struc(n)%aspectfgrd )

          call readaspect( trim(LDT_LSMparam_struc(n)%aspect%source)//char(0), &
               n, LDT_LSMparam_struc(n)%aspect%num_bins,  &
               LDT_LSMparam_struc(n)%aspectfgrd%value,    &
               LDT_LSMparam_struc(n)%aspect%value)

        ! Fill where parameter values are missing compared to land/water mask:
          if( aspect%filltype == "neighbor" .or. aspect%filltype == "average" ) then
             write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                trim(LDT_LSMparam_struc(n)%aspect%short_name)
             write(fill_logunit,*) "Checking/filling mask values for: ", &
                                trim(LDT_LSMparam_struc(n)%aspect%short_name)

             if( LDT_LSMparam_struc(n)%aspect%source == "SRTM" ) then
                aspect%watervalue = 3.1416
             else  ! e.g., GTOPO30
                aspect%watervalue = LDT_rc%udef
             endif
             if( LDT_rc%topo_gridtransform(n) == "tile" ) then
               call LDT_contTileParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                   LDT_rc%topo_gridtransform(n),                                &
                   LDT_LSMparam_struc(n)%aspect%num_bins,                    &
                   LDT_LSMparam_struc(n)%aspect%value,                       &
                   LDT_LSMparam_struc(n)%aspectfgrd%value, aspect%watervalue,&
                   LDT_LSMparam_struc(n)%landmask2%value, aspect%filltype,   &
                   aspect%fillvalue, aspect%fillradius )
          !- Non-tiled output data fields:
             else
               call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n),&
                   LDT_rc%topo_gridtransform(n),                               &
                   LDT_LSMparam_struc(n)%aspect%num_bins,                   &
                   LDT_LSMparam_struc(n)%aspect%value, aspect%watervalue,   &
                   LDT_LSMparam_struc(n)%landmask2%value, aspect%filltype,  &
                   aspect%fillvalue, aspect%fillradius )
             endif
          endif
       endif

       if(LDT_LSMparam_struc(n)%curvature%selectOpt.eq.1) then
          call readcurv(&
               trim(LDT_LSMparam_struc(n)%curvature%source)//char(0),n, &
               LDT_LSMparam_struc(n)%curvature%value)
       endif       
    enddo
  end subroutine LDT_topo_init


  subroutine LDT_topo_writeHeader(n,ftn,dimID)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer   :: n 
    integer   :: ftn
    integer   :: dimID(3)
    integer   :: tdimID(3)

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

 !- Write elevation headers:
    if( LDT_LSMparam_struc(n)%elevation%selectOpt.eq.1 ) then
      LDT_LSMparam_struc(n)%elevation%vlevels = LDT_LSMparam_struc(n)%elevation%num_bins
      LDT_LSMparam_struc(n)%elevfgrd%vlevels = LDT_LSMparam_struc(n)%elevfgrd%num_bins

      call LDT_verify(nf90_def_dim(ftn,'elevbins',&
           LDT_LSMparam_struc(n)%elevation%num_bins,tdimID(3)))

      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%elevfgrd)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%elevation)
    endif
 !- Write slope headers:
    if( LDT_LSMparam_struc(n)%slope%selectOpt.eq.1 ) then
      LDT_LSMparam_struc(n)%slope%vlevels = LDT_LSMparam_struc(n)%slope%num_bins
      LDT_LSMparam_struc(n)%slopefgrd%vlevels = LDT_LSMparam_struc(n)%slopefgrd%num_bins
      call LDT_verify(nf90_def_dim(ftn,'slopebins',&
           LDT_LSMparam_struc(n)%slope%num_bins,tdimID(3)))

      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%slopefgrd)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%slope)
    endif

 !- Write aspect headers:
    if( LDT_LSMparam_struc(n)%aspect%selectOpt.eq.1 ) then
      LDT_LSMparam_struc(n)%aspect%vlevels = LDT_LSMparam_struc(n)%aspect%num_bins
      LDT_LSMparam_struc(n)%aspectfgrd%vlevels = LDT_LSMparam_struc(n)%aspectfgrd%num_bins
      call LDT_verify(nf90_def_dim(ftn,'aspectbins',&
           LDT_LSMparam_struc(n)%aspect%num_bins,tdimID(3)))

      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%aspectfgrd)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%aspect)
    endif

 !- Write curvature headers:
    if(LDT_LSMparam_struc(n)%curvature%selectOpt.eq.1) then
      call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
           LDT_LSMparam_struc(n)%curvature)
    endif
#endif

  end subroutine LDT_topo_writeHeader

  subroutine LDT_topo_writeData(n,ftn)

    use LDT_coreMod, only : LDT_rc

    integer  :: n 
    integer  :: ftn
    
    if(LDT_LSMparam_struc(n)%elevation%selectOpt.eq.1) then
      call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%elevfgrd)
      call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%elevation)
    endif
    if(LDT_LSMparam_struc(n)%slope%selectOpt.eq.1) then
      call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%slope)
      call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%slopefgrd)
    endif
    if(LDT_LSMparam_struc(n)%aspect%selectOpt.eq.1) then
      call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%aspect)
      call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%aspectfgrd)
    endif
    if(LDT_LSMparam_struc(n)%curvature%selectOpt.eq.1) then
      call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%curvature)
    endif

  end subroutine LDT_topo_writeData

end module LDT_topoMod
