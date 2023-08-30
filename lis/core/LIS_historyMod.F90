!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
#include "LIS_NetCDF_inc.h"
module LIS_historyMod
!BOP
!
! !MODULE: LIS_historyMod
! 
! !DESCRIPTION: 
!  The code in this file provides interfaces to manage LIS model output 
!  in different file formats. 
!  
! \subsubsection{Overview}
! The module provides generic interfaces to manage output from different
! land surface models. Currently LIS supports three data formats. 
! \begin{itemize}
! \item{binary: in binary, big-endian.}
! \item{grib1: WMO GRIB-1 format}
! \item{grib2: WMO GRIB-1 format}
! \item{netcdf: NCAR unidata netcdf-4 format}
! \end{itemize}
! This module also provides generic interfaces for writing and reading
! model restart files. LIS uses the following convention on the 
! filename extensions: 
! \begin{itemize}
! \item{.1gd4r: binary output, 1 correspondes to one variable, g represents
! gridded data, d represnts direct access, 4r represents 4 byte real}
! \item{.grb: files in GRIB-1 format}
! \item{.grb2: files in GRIB-2 format}
! \item{.nc: files in netcdf-4 format}
! \item{*rst: restart files. Currently all restart files are written in
! netcdf-4 format}
! \end{itemize}
! This module also manages any data aggregation from computing nodes
! when parallel processing is employed. 
!
! !NOTES: LIS history writer employs grib-api library for writing grib
!  output and NETCDF-4 for writing NETCDF output. If the support for these
!  data formats are expected, these libraries must 
!  be downloaded externally and built before LIS compilation. 
! 
!  GRIB-API: http://www.ecmwf.int/products/data/software/grib\_api.html \newline
!  NETCDF-4: http://www.unidata.ucar.edu/software/netcdf/   \newline
!  
! !REVISION HISTORY:
!  10 Feb 2004: Sujay Kumar; Initial Specification
!   4 Jul 2008: Sujay Kumar; Redesigned with generic routines that can handle
!                   model output in different formats for all LSMs
!  25 Jan 2012: Sujay Kumar; The grib writers were modified to use grib-api
!  28 Feb 2012: Sujay Kumar; Updated NETCDF interfaces to use NETCDF-4 library
!  28 Jan 2014: David Mocko; Updates for AFWA GRIB-1 files using GRIBAPI library
!  21 Feb 2014: David Mocko; More updates for matching AFWA GRIB-1 files
!                   to the LIS-6 GRIB output
!  28 Mar 2014: David Mocko; More refinements/fixes to AFWA GRIB-1 CONFIGS
!   1 Apr 2015: Hiroko Beaudoing; Added GRIB-2 routines
!  15 May 2015: Hiroko Beaudoing; Added nsoillayers2, lyrthk2 for when soil
!                   moisture and temperature having different number of layers
!                   used in GRIB1 & GRIB2 format
!  18 Oct 2018: David Mocko: Check lis.config entry for option to turn off
!                   writing ASCII stats files with netCDF output format
!
! !USES: 
  use LIS_coreMod
  use LIS_histDataMod
  use LIS_timeMgrMod
  use LIS_logMod

#if ( defined USE_GRIBAPI)
  use grib_api
#endif

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  use LIS_mpiMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_writevar_bin     ! write a variable into a binary file
  public :: LIS_writevar_grib1   ! write a variable into a GRIB-1 file
  public :: LIS_writevar_grib2   ! write a variable into a GRIB-2 file
  public :: LIS_writevar_netcdf  ! write a variable into a netcdf file
  public :: LIS_writevar_restart ! writes a variable into a restart file
  public :: LIS_readvar_restart  ! reads a variable from a restart file
  public :: LIS_writevar_reduced_tilespace! writes a variable in reduced tilespace (ngrid)
  public :: LIS_readvar_reduced_tilespace ! reads a reduced tilespace variable
  public :: LIS_readvar_gridded  ! read a variable from a gridded file
  public :: LIS_writevar_gridded ! write a variable in grid space to a file
  public :: LIS_tile2grid        ! convert a tilespace variable to gridspace
  public :: LIS_grid2tile        ! convert gridspace to tilespace
  public :: LIS_patch2tile       ! convert a patchspace variable to tilespace
  public :: LIS_grid2patch       ! convert a gridspace variable to patchspace
  public :: LIS_writevar_spread  ! writes ensemble spared to a gridded file
  public :: LIS_writevar_incr    ! writes analysis increments to a gridded file
  public :: LIS_writeModelOutput ! writes model output based on the selected
  public :: LIS_writeRoutingModelOutput                                           ! format and list of variables
  public :: LIS_gather_gridded_output      ! gather the 1d tiled output variable into a 2d gridded array
  public :: LIS_gather_tiled_vector_output ! gather the 1d tiled output variables in tile space
  public :: LIS_gather_gridded_vector_output ! gather the 1d tiled output variable into a 1d gridded array
  public :: LIS_gather_tiled_vector_withhalo_output
  public :: LIS_writeGlobalHeader_restart
  public :: LIS_writeHeader_restart
  public :: LIS_closeHeader_restart
  public :: LIS_convertVarToLocalSpace
  public :: LIS_gather_1dgrid_to_2dgrid
  public :: LIS_scatter_global_to_local_grid
  public :: LIS_gather_2d_local_to_global
!EOP

!BOP 
! 
! !ROUTINE: LIS_writevar_bin
! \label{LIS_writevar_bin}
! 
! !INTERFACE:
  interface LIS_writevar_bin
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure writevar_bin_real
     module procedure writevar_bin_real_direct
     module procedure writevar_bin_withstats_real     
! 
! !DESCRIPTION:
! This interface provides routines for writing variables (real)
! in a binary file. The aggregation from decomposed tile spaces in 
! each processor to their respective grid spaces and the aggregations 
! of these grid spaces from each processor are also performed by 
! this routine. The interface also provides options to write some 
! diagnostic statistics about the variable being written. 
!
!EOP 
  end interface

!BOP 
! 
! !ROUTINE: LIS_writevar_restart
! \label{LIS_writevar_restart}
!
! !INTERFACE:
  interface LIS_writevar_restart
! !PRIVATE MEMBER FUNCTIONS:

     module procedure writevar_restart_tile_int
     module procedure writevar_restart_tile_real
     module procedure writevar_restart_patch_int
     module procedure writevar_restart_patch_real
!     module procedure writevar_restart_char
! 
! !DESCRIPTION: 
! This interface provides routines for writing variables (integer or real)
! in a restart file. By definition, the restart files are written in the
! model tile space. The aggregations from tile spaces of each processors 
! are also performed by this routine.  
!
!EOP
  end interface

!BOP
! !ROUTINE: LIS_grid2tile
! \label{LIS_grid2tile}
!
! !INTERFACE:
  interface LIS_grid2tile
     module procedure grid2tile
     module procedure grid2tile_ens
  end interface
  
!BOP
! 
! !ROUTINE: LIS_writevar_reduced_tilespace
! \label{LIS_writevar_reduced_tilespace}
! 
! !INTERFACE: 
  interface LIS_writevar_reduced_tilespace
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure writevar_tile_noensemble_real
! 
! !DESCRIPTION: 
! This interface provides routines for writing variables in a
! reduced tilespace (tilespace without ensembles). 
! 
!EOP
  end interface

!BOP
! 
! !ROUTINE: LIS_readvar_reduced_tilespace
! \label{LIS_readvar_reduced_tilespace}
! 
! !INTERFACE: 
  interface LIS_readvar_reduced_tilespace
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure readvar_tile_noensemble_real
! 
! !DESCRIPTION: 
! This interface provides routines for writing variables in a
! reduced tilespace (tilespace without ensembles). 
! 
!EOP
  end interface
!BOP 
! 
! !ROUTINE: LIS_readvar_restart
! \label{LIS_readvar_restart}
!
! !INTERFACE:
  interface LIS_readvar_restart
! !PRIVATE MEMBER FUNCTIONS:
     module procedure readvar_restart_tile_int
     module procedure readvar_restart_tile_real
     module procedure readvar_restart_patch_real
     module procedure readvar_restart_patch_int
!     module procedure readvar_restart_char

! !DESCRIPTION:
! This interface provides routines for reading variables (integer or real)
! from a restart file. The restart files are read by the master processor. 
! The domain decomposition of the ``global'' tile space on the master processor
! to individual processors are also performed by this routine. 
!EOP
  end interface

!BOP 
! 
! !ROUTINE: LIS_readvar_gridded
! \label{LIS_readvar_gridded}
!
! !INTERFACE:
  interface LIS_readvar_gridded
! !PRIVATE MEMBER FUNCTIONS:

     module procedure readvar_1dgridded_real
     module procedure readvar_2dgridded_real
     module procedure readvar_1dgridded_fromvector_real

! !DESCRIPTION:
! This interface provides routines for reading variables (real)
! from a gridded (output) file into a 1d or 2d local gridded space. 
! The gridded files are read by the master processor. 
! The domain decomposition of the ``global'' grid space on the master processor
! to individual processors are also performed by this routine. 
!EOP
  end interface

!BOP 
! 
! !ROUTINE: LIS_writevar_gridded
! \label{LIS_writevar_gridded}
!
! !INTERFACE:
  interface LIS_writevar_gridded
! !PRIVATE MEMBER FUNCTIONS:

     module procedure writevar_gridded_real
     module procedure writevar_gridded_real_withstats

! !DESCRIPTION:
! This interface provides routines for writing variables (real)
! in grid space to a file. 
!EOP
  end interface


!BOP 
! 
! !ROUTINE: LIS_writevar_grib1
! \label{LIS_writevar_grib1}
!
! !INTERFACE:
  interface LIS_writevar_grib1
! !PRIVATE MEMBER FUNCTIONS:
     module procedure writevar_grib1_withstats_real
! 
! !DESCRIPTION:
! This interface provides routines for writing variables (real)
! in a GRIB-1 file. The aggregation from decomposed tile spaces in 
! each processor to their respective grid spaces and the aggregations 
! of these grid spaces from each processor are also performed by 
! this routine. The interface also provides options to write some 
! diagnostic statistics about the variable being written. 
!
!EOP
  end interface
!BOP 
! 
! !ROUTINE: LIS_writevar_grib2
! \label{LIS_writevar_grib2}
!
! !INTERFACE:
  interface LIS_writevar_grib2
! !PRIVATE MEMBER FUNCTIONS:
     module procedure writevar_grib2_withstats_real
! 
! !DESCRIPTION:
! This interface provides routines for writing variables (real)
! in a GRIB-2 file. The aggregation from decomposed tile spaces in 
! each processor to their respective grid spaces and the aggregations 
! of these grid spaces from each processor are also performed by 
! this routine. The interface also provides options to write some 
! diagnostic statistics about the variable being written. 
!
!EOP
  end interface
!BOP
! 
! !ROUTINE: LIS_grid2patch
! \label{LIS_grid2patch}
! 
! !INTERFACE: 
  interface LIS_grid2patch
! !PRIVATE MEMBER FUNCTIONS: 
!     module procedure grid2patch
     module procedure grid2patch_local
     module procedure grid2patch_global
     module procedure grid2patch_global_ens
!
! !DESCRIPTON: 
!
!EOP
  end interface


!BOP 
! 
! !ROUTINE: LIS_writevar_netcdf
! \label{LIS_writevar_netcdf}
! 
! !INTERFACE:
  interface LIS_writevar_netcdf
! !PRIVATE MEMBER FUNCTIONS:

     module procedure writevar_netcdf_withstats_real
! 
! !DESCRIPTION:
! This interface provides routines for writing variables (2d or 3d)
! in a netcdf file. The aggregation from decomposed tile spaces in 
! each processor to their respective grid spaces and the aggregations 
! of these grid spaces from each processor are also performed by 
! this routine.
!
!EOP
  end interface

!BOP 
! 
! !ROUTINE: LIS_gather_gridded_output
! \label{LIS_gather_gridded_output}
! 
! !INTERFACE:
  interface LIS_gather_gridded_output
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure gather_gridded_output_tile
     module procedure gather_gridded_output_patch
! 
! !DESCRIPTION:
! This interface provides routines for writing variables (real)
! in a binary file. The aggregation from decomposed tile spaces in 
! each processor to their respective grid spaces and the aggregations 
! of these grid spaces from each processor are also performed by 
! this routine. The interface also provides options to write some 
! diagnostic statistics about the variable being written. 
!
!EOP 
  end interface


!BOP 
! 
! !ROUTINE: LIS_tile2grid
! \label{LIS_tile2grid}
! 
! !INTERFACE:
  interface LIS_tile2grid
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure tile2grid_local
     module procedure tile2grid_local_ens
     module procedure tile2grid_global_ens
     module procedure tile2grid_global_noens
! 
! !DESCRIPTION:
! This interface provides routines for converting a tile
! space variable to the grid space 
!
!EOP 
  end interface

contains

!BOP
! !ROUTINE: LIS_writeModelOutput
! \label{LIS_writeModelOutput}
! 
! !INTERFACE: 
  subroutine LIS_writeModelOutput(n, lsmoutfile, lsmstatsfile, &
       sopen, nsoillayers, outInterval, lyrthk, &
       nsoillayers2, group, model_name, lyrthk2)
! !USES:
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
! !ARGUMENTS: 
    integer,   intent(in)   :: n 
    character(len=*),   intent(in)   :: lsmoutfile
    character(len=*),   intent(in)   :: lsmstatsfile
    logical,   intent(in)   :: sopen
    real,      intent(in)   :: outInterval
    integer,   intent(in)   :: nsoillayers
    real,      intent(in)   :: lyrthk(nsoillayers)
    integer,   intent(in)   :: nsoillayers2
    integer,   intent(in),optional :: group
    character(len=*), intent(in),optional :: model_name
    real,      intent(in),optional   :: lyrthk2(nsoillayers2)
! 
! !DESCRIPTION: 
!   This subroutine invokes the routine to write LSM output in the selected
!    data format (binary/grib1/netcdf) and using the selected list of variables. 
!    Further, the variables are also written as instantaneous, time averaged,
!    or accumulated based on the user specifications. 
!   
!   The arguments are: 
!   \begin{description}
!    \item[n]  index of the nest \newline
!    \item[lsmoutfile]  name of the LSM history file \newline
!    \item[lsmstatsfile] name of the LSM history stats file  \newline
!    \item[sopen] flag determining whether to open the LIS history stats file \newline
!    \item[outInterval]   history output frequency \newline
!    \item[nsoillayers]  Number of soil layers \newline
!    \item[lyrthk]   Thickness of soil layers \newline
!  \end{description}
! 
!   The routines invoked are: 
!   \begin{description}
!   \item[writeBinaryOutput](\ref{writeBinaryOutput}) \newline
!     writes the history files in binary format
!   \item[writeGribOutput](\ref{writeGribOutput})\\
!     writes the history files in GRIB-1 or GRIB-2 format
!   \item[writeNetcdfOutput](\ref{writeNetcdfOutput}) \newline
!     writes the history files in NETCDF format
!   \item[LIS\_resetOutputVars](\ref{LIS_resetOutputVars}) \newline
!     resets the time averaged varibles for the next output. 
!   \end{description}
!EOP

    integer :: ftn, ftn_stats, ftnp(LIS_npes)
    integer :: iret
    integer :: group_temp
    character*100 :: mname_temp,temp1
    character(len=LIS_CONST_PATH_LEN) :: lsmoutfilep
    character*1    :: fproc(4)

    if(.NOT.PRESENT(group)) then 
       group_temp = 1
    else
       group_temp = group
    endif
    
    if(.NOT.PRESENT(model_name)) then 
       mname_temp = "model_not_specified"
    else
       mname_temp = model_name
    endif

    if(LIS_masterproc) then 
       if ( LIS_rc%sout ) then
          if ( sopen ) then
             if(LIS_rc%startcode.eq."restart") then 
                open(65+n+10*group_temp,file=lsmstatsfile,&
                   form='formatted', position='append')
             else
                open(65+n+10*group_temp,file=lsmstatsfile,&
                   form='formatted')
             endif
          endif

          write(65+n+10*group_temp,*)              
          write(65+n+10*group_temp,996)&
             '       Statistical Summary of LSM output for:  ', & 
             LIS_rc%mo,'/',LIS_rc%da,'/',LIS_rc%yr,LIS_rc%hr,':',&
             LIS_rc%mn,':',LIS_rc%ss
          996    format(a47,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
          write(65+n+10*group_temp,*)
          write(65+n+10*group_temp,997)
          997    format(t27,'Mean',t41,'Stdev',t56,'Min',t70,'Max')
       endif
    endif
    ftn = 12
    ftn_stats = 65 + n +10*group_temp

    call LIS_rescaleCount(n,group_temp)

    write(LIS_logunit,*)'[INFO] Writing surface model output to:  ', &
         trim(lsmoutfile) ! EMK

    if(LIS_rc%wout.eq."binary") then 
       if(LIS_masterproc) then 
          open(ftn,file=lsmoutfile,form='unformatted')
       endif
       call writeBinaryOutput(n,group_temp,ftn,ftn_stats)
       if(LIS_masterproc) then 
          close(ftn)
       endif
    elseif(LIS_rc%wout.eq."distributed binary") then
       ftnp(LIS_localPet+1) = LIS_getNextUnitNumber()

       write(temp1,'(i4.4)') LIS_localPet
       read(temp1,fmt='(4a1)') fproc
       
       lsmoutfilep = trim(lsmoutfile)//'.'//fproc(1)//fproc(2)//fproc(3)//fproc(4)
       open(ftnp(LIS_localPet+1),file=trim(lsmoutfilep),&
            form='unformatted')

       write(ftnp(LIS_localPet+1)) LIS_ews_ind(n,LIS_localPet+1)
       write(ftnp(LIS_localPet+1)) LIS_ewe_ind(n,LIS_localPet+1)
       write(ftnp(LIS_localPet+1)) LIS_nss_ind(n,LIS_localPet+1)
       write(ftnp(LIS_localPet+1)) LIS_nse_ind(n,LIS_localPet+1)
       
       call writeBinaryDistributedOutput(n,group_temp,ftnp(LIS_localPet+1))

       call LIS_releaseUnitNumber(ftnp(LIS_localPet+1))
       
    elseif(LIS_rc%wout.eq."grib1") then 
#if(defined USE_GRIBAPI)
       if(LIS_masterproc) then 
          call grib_open_file(ftn,lsmoutfile,'w',iret)
          call LIS_verify(iret, 'failed to open grib file '//lsmoutfile)

       endif
       call writeGribOutput(n,group_temp,ftn,ftn_stats,outInterval,&
                             nsoillayers,lyrthk,nsoillayers2)
       if(LIS_masterproc) then 
          call grib_close_file(ftn,iret)
          call LIS_verify(iret, 'failed to close grib file'//lsmoutfile)
       endif
#endif          
    elseif(LIS_rc%wout.eq."grib2") then 
#if(defined USE_GRIBAPI)
       if(LIS_masterproc) then 
          call grib_open_file(ftn,lsmoutfile,'w',iret)
          call LIS_verify(iret, 'failed to open grib2 file '//lsmoutfile)

       endif
       if(.NOT.PRESENT(lyrthk2)) then 
       call writeGribOutput(n,group_temp,ftn,ftn_stats,outInterval,&
                             nsoillayers,lyrthk,nsoillayers2)
       else
       call writeGribOutput(n,group_temp,ftn,ftn_stats,outInterval,&
                             nsoillayers,lyrthk,nsoillayers2,lyrthk2)
       endif
       if(LIS_masterproc) then 
          call grib_close_file(ftn,iret)
          call LIS_verify(iret, 'failed to close grib2file'//lsmoutfile)
       endif
#endif          
    elseif(LIS_rc%wout.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
       if(LIS_masterproc) then 
#if (defined USE_NETCDF4)
          iret = nf90_create(path=lsmoutfile,cmode=nf90_hdf5,&
               ncid = ftn)
          call LIS_verify(iret,'creating netcdf file failed in LIS_historyMod')
#endif
#if (defined USE_NETCDF3)
          iret = nf90_create(path=lsmoutfile,cmode=nf90_clobber,&
               ncid = ftn)
          call LIS_verify(iret,'creating netcdf file failed in LIS_historyMod')
#endif
       endif
       call writeNetcdfOutput(n,group_temp,ftn,ftn_stats, outInterval, &
            nsoillayers, lyrthk, mname_temp)
       if(LIS_masterproc) then 
          iret = nf90_close(ftn)
       endif
#endif
    endif
    ! After writing reset the variables
    call LIS_resetOutputVars(n,group_temp)
  end subroutine LIS_writeModelOutput


!BOP
! !ROUTINE: LIS_writeRoutingModelOutput
! \label{LIS_writeRoutingModelOutput}
! 
! !INTERFACE: 
  subroutine LIS_writeRoutingModelOutput(n, routingoutfile, routingstatsfile, &
       sopen, nsoillayers, outInterval, lyrthk, &
       nsoillayers2, group, model_name, lyrthk2)
! !USES:
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

! !ARGUMENTS: 
    integer,   intent(in)   :: n 
    character(len=*),   intent(in)   :: routingoutfile
    character(len=*),   intent(in)   :: routingstatsfile
    logical,   intent(in)   :: sopen
    real,      intent(in)   :: outInterval
    integer,   intent(in)   :: nsoillayers
    real,      intent(in)   :: lyrthk(nsoillayers)
    integer,   intent(in)   :: nsoillayers2
    integer,   intent(in),optional :: group
    character(len=*), intent(in),optional :: model_name
    real,      intent(in),optional   :: lyrthk2(nsoillayers2)
! 
! !DESCRIPTION: 
!   This subroutine invokes the routine to write routing output in the selected
!    data format (binary/grib1/netcdf) and using the selected list of variables. 
!    Further, the variables are also written as instantaneous, time averaged,
!    or accumulated based on the user specifications. 
!   
!   The arguments are: 
!   \begin{description}
!    \item[n]  index of the nest \newline
!    \item[routingoutfile]  name of the ROUTING history file \newline
!    \item[routingstatsfile] name of the ROUTING history stats file  \newline
!    \item[sopen] flag determining whether to open the LIS history stats file \newline
!    \item[outInterval]   history output frequency \newline
!    \item[nsoillayers]  Number of soil layers \newline
!    \item[lyrthk]   Thickness of soil layers \newline
!  \end{description}
! 
!   The routines invoked are: 
!   \begin{description}
!   \item[writeBinaryOutput](\ref{writeBinaryOutput}) \newline
!     writes the history files in binary format
!   \item[writeGribOutput](\ref{writeGribOutput})\\
!     writes the history files in GRIB-1 or GRIB-2 format
!   \item[writeNetcdfOutput](\ref{writeNetcdfOutput}) \newline
!     writes the history files in NETCDF format
!   \item[LIS\_resetOutputVars](\ref{LIS_resetOutputVars}) \newline
!     resets the time averaged varibles for the next output. 
!   \end{description}
!EOP

    integer :: ftn, ftn_stats, ftnp(LIS_npes)
    integer :: iret
    integer :: group_temp
    character*100 :: mname_temp, temp1
    character(len=LIS_CONST_PATH_LEN) :: routingoutfilep
    character*1   :: fproc(4)

    if(.NOT.PRESENT(group)) then 
       group_temp = 1
    else
       group_temp = group
    endif
    
    if(.NOT.PRESENT(model_name)) then 
       mname_temp = "model_not_specified"
    else
       mname_temp = model_name
    endif

    if(LIS_masterproc) then 
       if ( LIS_rc%sout ) then
          if ( sopen ) then
             if(LIS_rc%startcode.eq."restart") then 
                open(65+n+10*group_temp,file=routingstatsfile,&
                   form='formatted', position='append')
             else
                open(65+n+10*group_temp,file=routingstatsfile,&
                   form='formatted')
             endif
          endif

          write(65+n+10*group_temp,*)              
          write(65+n+10*group_temp,996)&
             '       Statistical Summary of ROUTING output for:  ', & 
             LIS_rc%mo,'/',LIS_rc%da,'/',LIS_rc%yr,LIS_rc%hr,':',&
             LIS_rc%mn,':',LIS_rc%ss
          996    format(a51,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
          write(65+n+10*group_temp,*)
          write(65+n+10*group_temp,997)
          997    format(t27,'Mean',t41,'Stdev',t56,'Min',t70,'Max')
       endif
    endif
    ftn = 12
    ftn_stats = 65 + n +10*group_temp

    call LIS_rescaleCount(n,group_temp)

    write(LIS_logunit,*)'[INFO] Writing routing model output to:  ', &
         trim(routingoutfile) ! EMK

    if(LIS_rc%wout.eq."binary") then 
       write(LIS_logunit,*)'[ERR] binary routing model outputs are not supported'
       call LIS_endrun()
    elseif(LIS_rc%wout.eq."distributed binary") then

       ftnp(LIS_localPet+1) = LIS_getNextUnitNumber()

       write(temp1,'(i4.4)') LIS_localPet
       read(temp1,fmt='(4a1)') fproc
       
       routingoutfilep = trim(routingoutfile)//'.'//fproc(1)//fproc(2)//fproc(3)//fproc(4)
       open(ftnp(LIS_localPet+1),file=trim(routingoutfilep),&
            form='unformatted')

       write(ftnp(LIS_localPet+1)) LIS_ews_ind(n,LIS_localPet+1)
       write(ftnp(LIS_localPet+1)) LIS_ewe_ind(n,LIS_localPet+1)
       write(ftnp(LIS_localPet+1)) LIS_nss_ind(n,LIS_localPet+1)
       write(ftnp(LIS_localPet+1)) LIS_nse_ind(n,LIS_localPet+1)
       
       call writeRoutingBinaryDistributedOutput(n,group_temp,ftnp(LIS_localPet+1))

       call LIS_releaseUnitNumber(ftnp(LIS_localPet+1))       
       
    elseif(LIS_rc%wout.eq."grib1") then 
#if(defined USE_GRIBAPI)
       write(LIS_logunit,*)'[ERR] grib1 routing model outputs are not supported'
       call LIS_endrun()
#endif          
    elseif(LIS_rc%wout.eq."grib2") then 
#if(defined USE_GRIBAPI)
       write(LIS_logunit,*)'[ERR] grib2 routing model outputs are not supported'
       call LIS_endrun()
#endif          
    elseif(LIS_rc%wout.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
       if(LIS_masterproc) then 
#if (defined USE_NETCDF4)
          iret = nf90_create(path=routingoutfile,cmode=nf90_hdf5,&
               ncid = ftn)
          call LIS_verify(iret,'creating netcdf file failed in LIS_historyMod')
#endif
#if (defined USE_NETCDF3)
          iret = nf90_create(path=routingoutfile,cmode=nf90_clobber,&
               ncid = ftn)
          call LIS_verify(iret,'creating netcdf file failed in LIS_historyMod')
#endif
       endif
       call writeRoutingNetcdfOutput(n,group_temp,ftn,ftn_stats, outInterval, &
            nsoillayers, lyrthk, mname_temp)
       if(LIS_masterproc) then 
          iret = nf90_close(ftn)
       endif
#endif
    endif
    ! After writing reset the variables
    call LIS_resetOutputVars(n,group_temp)
  end subroutine LIS_writeRoutingModelOutput

!BOP
! 
! !ROUTINE: writeBinaryOutput
! \label{writeBinaryOutput}
! 
! !INTERFACE: 
  subroutine writeBinaryOutput(n, group, ftn, ftn_stats)
!  !USES: 

! !ARGUMENTS: 
    implicit none
    integer,   intent(in)   :: n 
    integer,   intent(in)   :: group
    integer,   intent(in)   :: ftn
    integer,   intent(in)   :: ftn_stats
! 
! !DESCRIPTION: 
!  This routine writes a binary output file based on the list of selected 
!  output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest \newline
!    \item[group] output group (LSM, routing, RTM) \newline
!    \item[ftn] file unit for the output file \newline
!    \item[ftn\_stats] file unit for the output statistics file \newline
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[writeSingleBinaryVar](\ref{writeSingleBinaryVar}) \newline
!     writes a single variable into a flat binary file. 
!   \end{description}
!EOP
    type(LIS_metadataEntry), pointer :: dataEntry

    if(group.eq.1) then !LSM output
       dataEntry => LIS_histData(n)%head_lsm_list
    elseif(group.eq.2) then !ROUTING
       dataEntry => LIS_histData(n)%head_routing_list
    elseif(group.eq.3) then !RTM
       dataEntry => LIS_histData(n)%head_rtm_list
    elseif(group.eq.4) then !Irrigation
       dataEntry => LIS_histData(n)%head_irrig_list
    endif


    do while ( associated(dataEntry) )
       call writeSingleBinaryVar(ftn,ftn_stats,n,dataEntry)
       dataEntry => dataEntry%next
    enddo
    
  end subroutine writeBinaryOutput


!BOP
! 
! !ROUTINE: writeBinaryDistributedOutput
! \label{writeBinaryDistributedOutput}
! 
! !INTERFACE: 
  subroutine writeBinaryDistributedOutput(n, group, ftn)
!  !USES: 

! !ARGUMENTS: 
    implicit none
    integer,   intent(in)   :: n 
    integer,   intent(in)   :: group
    integer,   intent(in)   :: ftn

! 
! !DESCRIPTION: 
!  This routine writes a binary distributed
!  output file based on the list of selected output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest \newline
!    \item[group] output group (LSM, routing, RTM) \newline
!    \item[ftn] file unit for the output file \newline
!    \item[ftn\_stats] file unit for the output statistics file \newline
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[writeSingleBinaryVar](\ref{writeSingleBinaryVar}) \newline
!     writes a single variable into a flat binary file. 
!   \end{description}
!EOP
    type(LIS_metadataEntry), pointer :: dataEntry

    if(group.eq.1) then !LSM output
       dataEntry => LIS_histData(n)%head_lsm_list
    elseif(group.eq.2) then !ROUTING
       dataEntry => LIS_histData(n)%head_routing_list
    elseif(group.eq.3) then !RTM
       dataEntry => LIS_histData(n)%head_rtm_list
    elseif(group.eq.4) then !Irrigation
       dataEntry => LIS_histData(n)%head_irrig_list
    endif


    do while ( associated(dataEntry) )
       call writeSingleBinaryDistributedVar(ftn,n,dataEntry)
       dataEntry => dataEntry%next
    enddo
    
  end subroutine writeBinaryDistributedOutput

!BOP
! 
! !ROUTINE: writeRoutingBinaryDistributedOutput
! \label{writeRoutingBinaryDistributedOutput}
! 
! !INTERFACE: 
  subroutine writeRoutingBinaryDistributedOutput(n, group, ftn)
!  !USES: 

! !ARGUMENTS: 
    implicit none
    integer,   intent(in)   :: n 
    integer,   intent(in)   :: group
    integer,   intent(in)   :: ftn

! 
! !DESCRIPTION: 
!  This routine writes a binary distributed output file
!  based on the list of selected routing output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest \newline
!    \item[group] output group (LSM, routing, RTM) \newline
!    \item[ftn] file unit for the output file \newline
!    \item[ftn\_stats] file unit for the output statistics file \newline
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[writeSingleBinaryVar](\ref{writeSingleBinaryVar}) \newline
!     writes a single variable into a flat binary file. 
!   \end{description}
!EOP
    type(LIS_metadataEntry), pointer :: dataEntry

    dataEntry => LIS_histData(n)%head_routing_list

    do while ( associated(dataEntry) )
       call writeRoutingSingleBinaryDistributedVar(ftn,n,dataEntry)
       dataEntry => dataEntry%next
    enddo
    
  end subroutine writeRoutingBinaryDistributedOutput
  
!BOP
!
! !ROUTINE: writeSingleBinaryVar
! \label{writeSingleBinaryVar}
! 
! !INTERFACE:
  subroutine writeSingleBinaryVar(ftn, ftn_stats, n, dataEntry)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer :: ftn
    integer :: ftn_stats
    integer :: n 
    type(LIS_metadataEntry), pointer :: dataEntry
! 
! !DESCRIPTION: 
!  This routine writes a single variable to a binary file
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[dataEntry] object containing the specifications for 
!                     the output variable
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_writevar\_bin](\ref{LIS_writevar_bin})
!     writes a variable into a flat binary file. 
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \end{description}
!EOP
    integer :: k,i

    real, allocatable    :: var(:), meanv(:), stdv(:)
    integer :: t, m, c

    if(dataEntry%selectOpt.eq.1) then 

       do t=1,LIS_rc%ntiles(n)
          m = LIS_domain(n)%tile(t)%sftype
          do k=1,dataEntry%vlevels
             if(dataEntry%count(t,k).gt.0) then 
                if(dataEntry%timeAvgOpt.eq.3) then  !do nothing
                   continue       
                elseif(dataEntry%timeAvgOpt.eq.2.or.dataEntry%timeAvgOpt.eq.1) then
                   dataEntry%modelOutput(1,t,k) = dataEntry%modelOutput(1,t,k)/&
                        dataEntry%count(t,k)
                else !do nothing
                   continue
                endif
             else
                dataEntry%modelOutput(1,t,k) = LIS_rc%udef
             endif
          enddo
       enddo
       
       do k=1,dataEntry%vlevels
          ! accumulated values
          ! time-averaged values and instantaneous values
          if(dataEntry%timeAvgOpt.eq.2) then 
             call LIS_writevar_bin(ftn,ftn_stats,n,&
                  dataEntry%modelOutput(1,:,k),&
                  trim(dataEntry%short_name)//'('//&
                  trim(dataEntry%units)//')',dataEntry%form)
             call LIS_writevar_bin(ftn,ftn_stats,n,&
                  dataEntry%modelOutput(2,:,k),&
                  trim(dataEntry%short_name)//'('//&
                  trim(dataEntry%units)//')',dataEntry%form)
             ! time-averaged or instantaneous values
          else
             call LIS_writevar_bin(ftn,ftn_stats,n,&
                  dataEntry%modelOutput(1,:,k),&
                  trim(dataEntry%short_name)//'('//&
                  trim(dataEntry%units)//')',dataEntry%form)
          endif
          
          if ( dataEntry%minMaxOpt.ne.0 ) then
             call LIS_writevar_bin(ftn,ftn_stats,n,&
                  dataEntry%minimum(:,k),&
                  trim(dataEntry%short_name)//'Min'//'('//&
                  trim(dataEntry%units)//')',dataEntry%form)
             call LIS_writevar_bin(ftn,ftn_stats,n,&
                  dataEntry%maximum(:,k),&
                  trim(dataEntry%short_name)//'Max'//'('//&
                  trim(dataEntry%units)//')',dataEntry%form)
          endif
          if( dataEntry%stdOpt.ne.0) then 
             
             allocate(var(LIS_rc%ntiles(n)))
             allocate(meanv(LIS_rc%ngrid(n)))
             allocate(stdv(LIS_rc%ngrid(n)))
             
             var(:) = dataEntry%modelOutput(1,:,k)
             
             meanv = 0                 
             do i=1,LIS_rc%ntiles(n), LIS_rc%nensem(n)
                c=LIS_domain(n)%tile(i)%index
                do m=1, LIS_rc%nensem(n)
                   t = i+m-1
                   meanv(c) = meanv(c) + var(t)*LIS_domain(n)%tile(t)%fgrd*&
                        LIS_domain(n)%tile(t)%pens
                enddo
             enddo
             
             stdv = 0                 
             do i=1,LIS_rc%ntiles(n), LIS_rc%nensem(n)
                c=LIS_domain(n)%tile(i)%index
                do m=1,LIS_rc%nensem(n)
                   t = i+m-1
                   stdv(c) = stdv(c) + LIS_domain(n)%tile(t)%fgrd*&
                        LIS_domain(n)%tile(t)%pens*&
                        (var(t)-meanv(c))**2
                enddo
             enddo
             stdv = sqrt(stdv)
             
             call LIS_writevar_gridded(ftn,ftn_stats,n,stdv,&
                  trim(dataEntry%short_name)//'Std'//'('//&
                  trim(dataEntry%units)//')',dataEntry%form)
             
             deallocate(var)
             deallocate(meanv)
             deallocate(stdv)
          endif
       enddo
    endif
  end subroutine writeSingleBinaryVar

!BOP
!
! !ROUTINE: writeSingleBinaryDistributedVar
! \label{writeSingleBinaryDistributedVar}
! 
! !INTERFACE:
  subroutine writeSingleBinaryDistributedVar(ftn, n, dataEntry)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer :: ftn
    integer :: n 
    type(LIS_metadataEntry), pointer :: dataEntry
! 
! !DESCRIPTION: 
!  This routine writes a single variable to a distributed binary file
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[dataEntry] object containing the specifications for 
!                     the output variable
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_writevar\_bin](\ref{LIS_writevar_bin})
!     writes a variable into a flat binary file. 
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \end{description}
!EOP
    integer :: k,i

    real, allocatable    :: var(:), meanv(:), stdv(:)
    character*100        :: mvar,units
    integer :: t, m, c

    if(dataEntry%selectOpt.eq.1) then 

       do t=1,LIS_rc%ntiles(n)
          m = LIS_domain(n)%tile(t)%sftype
          do k=1,dataEntry%vlevels
             if(dataEntry%count(t,k).gt.0) then 
                if(dataEntry%timeAvgOpt.eq.3) then  !do nothing
                   continue       
                elseif(dataEntry%timeAvgOpt.eq.2.or.dataEntry%timeAvgOpt.eq.1) then
                   dataEntry%modelOutput(1,t,k) = dataEntry%modelOutput(1,t,k)/&
                        dataEntry%count(t,k)
                else !do nothing
                   continue
                endif
             else
                dataEntry%modelOutput(1,t,k) = LIS_rc%udef
             endif
          enddo
       enddo

       
       if(dataEntry%timeAvgOpt.eq.0) then
          mvar = trim(dataEntry%short_name)//'_inst'
       elseif(dataEntry%timeAvgOpt.eq.1) then
          mvar = trim(dataEntry%short_name)//'_tavg'
       elseif(dataEntry%timeAvgOpt.eq.3) then
          mvar = trim(dataEntry%short_name)//'_acc'
       endif
       units = dataEntry%units
       
       write(ftn) dataEntry%vlevels
       write(ftn) mvar
       write(ftn) units
       
       do k=1,dataEntry%vlevels
          ! accumulated values
          ! time-averaged values and instantaneous values
          if(dataEntry%timeAvgOpt.eq.2) then 
             call writevar_dist_bin(ftn,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%form)
             call writevar_dist_bin(ftn,n,&
                  dataEntry%modelOutput(2,:,k),&
                  dataEntry%form)
             ! time-averaged or instantaneous values
          elseif(dataEntry%timeAvgOpt.eq.3) then 
             call writevar_dist_bin(ftn,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%form)
          elseif(dataEntry%timeAvgOpt.eq.0) then 
             call writevar_dist_bin(ftn,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%form)             
          elseif(dataEntry%timeAvgOpt.eq.1) then 
             call writevar_dist_bin(ftn,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%form)
          endif
          
       enddo
    endif
  end subroutine writeSingleBinaryDistributedVar

!BOP
!
! !ROUTINE: writeRoutingSingleBinaryDistributedVar
! \label{writeRoutingSingleBinaryDistributedVar}
! 
! !INTERFACE:
  subroutine writeRoutingSingleBinaryDistributedVar(ftn, n, dataEntry)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer :: ftn
    integer :: n 
    type(LIS_metadataEntry), pointer :: dataEntry
! 
! !DESCRIPTION: 
!  This routine writes a single routing variable to a distributed binary file
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[dataEntry] object containing the specifications for 
!                     the output variable
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_writevar\_bin](\ref{LIS_writevar_bin})
!     writes a variable into a flat binary file. 
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \end{description}
!EOP
    integer :: k,i

    real, allocatable    :: var(:), meanv(:), stdv(:)
    character*100        :: mvar,units
    integer :: t, m, c

    if(dataEntry%selectOpt.eq.1) then 

       do t=1,LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)
          do k=1,dataEntry%vlevels
             if(dataEntry%count(t,k).gt.0) then 
                if(dataEntry%timeAvgOpt.eq.3) then  !do nothing
                   continue       
                elseif(dataEntry%timeAvgOpt.eq.2.or.&
                     dataEntry%timeAvgOpt.eq.1) then
                   dataEntry%modelOutput(1,t,k) = dataEntry%modelOutput(1,t,k)/&
                        dataEntry%count(t,k)
                else !do nothing
                   continue
                endif
             else
                dataEntry%modelOutput(1,t,k) = LIS_rc%udef
             endif
          enddo
       enddo

       
       if(dataEntry%timeAvgOpt.eq.0) then
          mvar = trim(dataEntry%short_name)//'_inst'
       elseif(dataEntry%timeAvgOpt.eq.1) then
          mvar = trim(dataEntry%short_name)//'_tavg'
       elseif(dataEntry%timeAvgOpt.eq.3) then
          mvar = trim(dataEntry%short_name)//'_acc'
       endif
       units = dataEntry%units
       
       write(ftn) dataEntry%vlevels
       write(ftn) mvar
       write(ftn) units
       
       do k=1,dataEntry%vlevels
          ! accumulated values
          ! time-averaged values and instantaneous values
          if(dataEntry%timeAvgOpt.eq.2) then 
             call writeroutingvar_dist_bin(ftn,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%form)
             call writeroutingvar_dist_bin(ftn,n,&
                  dataEntry%modelOutput(2,:,k),&
                  dataEntry%form)
             ! time-averaged or instantaneous values
          elseif(dataEntry%timeAvgOpt.eq.3) then 
             call writeroutingvar_dist_bin(ftn,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%form)
          elseif(dataEntry%timeAvgOpt.eq.0) then 
             call writeroutingvar_dist_bin(ftn,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%form)             
          elseif(dataEntry%timeAvgOpt.eq.1) then 
             call writeroutingvar_dist_bin(ftn,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%form)
          endif
          
       enddo
    endif
  end subroutine writeRoutingSingleBinaryDistributedVar
  
!BOP
! !ROUTINE: writeGribOutput
! \label{writeGribOutput}
! 
! !INTERFACE: writeGribOutput
  subroutine writeGribOutput(n, group, ftn,ftn_stats, outInterval, &
                              nsoillayers, lyrthk, nsoillayers2, lyrthk2)
! !USES: 

! !ARGUMENTS: 
    implicit none
    integer,   intent(in)   :: n 
    integer,   intent(in)   :: group
    integer,   intent(in)   :: ftn
    integer,   intent(in)   :: ftn_stats
    real,      intent(in)   :: outInterval
    integer,   intent(in)   :: nsoillayers
    real,      intent(in)   :: lyrthk(nsoillayers)
    integer,   intent(in)   :: nsoillayers2
    real,      intent(in), optional   :: lyrthk2(nsoillayers2)

! 
! !DESCRIPTION: 
!  This routine writes an output file in the GRIB-1 or GRIB-2 format based on the 
!  list of selected output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[outInterval]   history output frequency
!    \item[nsoillayers]  Number of soil layers 
!    \item[lyrthk]   Thickness of soil layers
!                    The depths for layers below the land
!                    surface (surface = 112) are expected in centimeters,
!                    per GRIB-1 specifications.
!                    Convert lyrthk to cm in the call to LIS\_writeModelOutput.
!                    The depth below land surface (surface = 106) in GRIB-2
!                    is in meters. However, no need to change lyrthk (cm) 
!                    because a sacle factor can be applied to reflect correct
!                    depth units.
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[writeSingleGrib1Var](\ref{writeSingleGrib1Var}) \newline
!     writes a single variable into a GRIB-1 formatted file. 
!   \item[writeSingleGrib2Var](\ref{writeSingleGrib2Var})\\
!     writes a single variable into a GRIB-2 formatted file. 
!   \end{description}
!EOP
    integer                              :: i
    integer                              :: time_unit
    integer                              :: time_past
    integer                              :: time_curr
    real, dimension(nsoillayers)         :: toplev, botlev
    real, dimension(1)                   :: toplev0, botlev0
    type(LIS_metadataEntry), pointer     :: dataEntry
    real*8                               :: time
    real                                 :: gmt
    integer                              :: yr,mo,da,hr,mn,ss,doy
    real, dimension(nsoillayers2)        :: topleva, botleva
    real, dimension(1)                   :: toplev0a, botlev0a

    ! Setup of GRIB-1 and GRIB-2 Metadata Section
    
    ! toplev is the depth of the top of each soil layer
    ! botlev is the depth of the bottom of each soil layer
    toplev(1) = 0.0
    botlev(1) = lyrthk(1)

    ! determine bounding levels for each soil moisture layer
    do i = 2, nsoillayers
       toplev(i) = toplev(i-1) + lyrthk(i-1)
       botlev(i) = botlev(i-1) + lyrthk(i)
    enddo
    
    ! Set values for non layered fields (Fluxes, Sfc Fields, etc.)
    toplev0 = 0
    botlev0 = 0
    
    ! set second set of layer values
    if(.NOT.PRESENT(lyrthk2)) then
     topleva = 0.0
     botleva = 0.0
     toplev0a = 0.0
     botlev0a = 0.0
    else
     topleva(1) = 0.0
     botleva(1) = lyrthk2(1)
     do i = 2, nsoillayers2
        topleva(i) = topleva(i-1) + lyrthk2(i-1)
        botleva(i) = botleva(i-1) + lyrthk2(i)
     enddo
     toplev0a = 0
     botlev0a = 0
    endif

    yr = LIS_rc%yr
    mo = LIS_rc%mo
    da = LIS_rc%da
    hr = LIS_rc%hr
    mn = LIS_rc%mn
    ss = LIS_rc%ss

    call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,-1.0*outInterval)

    if(outInterval .GT. 0) then
       time_unit = 254     ! seconds
       time_curr = 0
       time_past = outInterval
    endif
    if(outInterval .GE. 60) then
       time_unit = 0      ! minutes       
       time_curr = 0
       time_past = (outInterval / 60)
    endif
    if(outInterval .GE. 3600) then
       time_unit = 1    ! hours
       time_curr = 0
       time_past = (outInterval / 3600)
    endif
    if(outInterval .GE. 86400) then
       time_unit = 2   ! days
       time_curr = 0
       time_past = (outInterval / 86400)
    endif

    !time_past: from LIS_grib1_finalize
    !time_P1 (Negative Time Unit for avg, or 0 for analysis) 
    !According to the in-line comments, time_past must be negative or 0.
    !Here we are setting it to a positive value.  This produces bad output.
    !Setting it to a negative value also produces bad output.
    !So I am resetting it to zero.  This produces output that matches
    !the binary output.
!    time_past=0

    if(group.eq.1) then !LSM output
       dataEntry => LIS_histData(n)%head_lsm_list
    elseif(group.eq.2) then !ROUTING
       dataEntry => LIS_histData(n)%head_routing_list
    elseif(group.eq.3) then !RTM
       dataEntry => LIS_histData(n)%head_rtm_list
    elseif(group.eq.4) then !Irrigation
       dataEntry => LIS_histData(n)%head_irrig_list
    endif

    do while ( associated(dataEntry) )
       if ( dataEntry%vlevels == 1 ) then
         if(LIS_rc%wout.eq."grib1") then 
          call writeSingleGrib1Var(ftn,ftn_stats,n,time_unit,time_past,&
               time_curr,1,toplev0,botlev0,1,toplev0a,botlev0a, &
               dataEntry%index,dataEntry)
         elseif(LIS_rc%wout.eq."grib2") then 
          call writeSingleGrib2Var(ftn,ftn_stats,n,time_unit,time_past,&
               time_curr,1,toplev0,botlev0,1,toplev0a,botlev0a, &
               dataEntry%index,dataEntry)
         endif
       else
         if(LIS_rc%wout.eq."grib1") then 
          call writeSingleGrib1Var(ftn,ftn_stats,n,time_unit,time_past,&
               time_curr,nsoillayers,toplev,botlev, &
               nsoillayers2,topleva,botleva,dataEntry%index,dataEntry)
         elseif(LIS_rc%wout.eq."grib2") then 
          call writeSingleGrib2Var(ftn,ftn_stats,n,time_unit,time_past,&
               time_curr,nsoillayers,toplev,botlev, &
               nsoillayers2,topleva,botleva,dataEntry%index,dataEntry)
         endif
       endif
       dataEntry => dataEntry%next
    enddo
  end subroutine writeGribOutput

!BOP
! !ROUTINE: writeSingleGrib1Var
! \label{writeSingleGrib1Var}
! 
! !INTERFACE: writeSingleGrib1Var
  subroutine writeSingleGrib1Var(ftn,ftn_stats,n,time_unit, time_past, &
       time_curr, nsoillayers,toplev, botlev, nsoillayers2, topleva, botleva, &
       var_index,dataEntry)
! !USES: 
    
    implicit none
! !ARGUMENTS:
    integer,   intent(in)    :: n 
    integer,   intent(in)    :: ftn
    integer,   intent(in)    :: ftn_stats    
    integer,   intent(in)    :: var_index
    integer,   intent(in)    :: time_unit
    integer,   intent(in)    :: time_past
    integer,   intent(in)    :: time_curr
    integer                  :: nsoillayers
    integer                  :: nsoillayers2
    real, dimension(nsoillayers) :: toplev, botlev
    real, dimension(nsoillayers2) :: topleva, botleva
    real, dimension(:),allocatable :: toplevtemp, botlevtemp
    type(LIS_metadataEntry),pointer :: dataEntry
! 
! !DESCRIPTION: 
!  This routine writes a single variable to a GRIB-1 file
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[var\_index] index of the variable
!    \item[dataEntry]
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_writevar\_grib1](\ref{LIS_writevar_grib1})
!     writes a variable into a GRIB-1 formatted file. 
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \end{description}
!EOP
    integer :: i,k,m,t,timeRange
    character*10 :: stepType
    integer :: ierr

    if(dataEntry%selectOpt.eq.1) then 

       if ((var_index.eq.LIS_MOC_SOILTEMP).and.                     &
           (nsoillayers .ne. nsoillayers2)) then
            allocate(toplevtemp(nsoillayers2))
            allocate(botlevtemp(nsoillayers2))
            toplevtemp(1:nsoillayers2) = topleva(1:nsoillayers2)
            botlevtemp(1:nsoillayers2) = botleva(1:nsoillayers2)
       else
            allocate(toplevtemp(nsoillayers))
            allocate(botlevtemp(nsoillayers))
            toplevtemp(1:nsoillayers) = toplev(1:nsoillayers)
            botlevtemp(1:nsoillayers) = botlev(1:nsoillayers)
       endif
#if (defined AFWA_GRIB_CONFIGS)
! Set timerange indicator equal to 133 for AFWA's specifications
! for surface runoff, baseflow, and total precipitation
! to make the LIS-7 output match the LIS-6 style. - dmm
       if ((var_index.eq.LIS_MOC_QS).or.                               &
           (var_index.eq.LIS_MOC_QSB).or.                              &
           (var_index.eq.LIS_MOC_TOTALPRECIP)) then
!          print *, "inside afwa_grib_configs(1): ", timeRange
          timeRange = 133
!          print *, "inside afwa_grib_configs(2): ", timeRange
!          stop
       endif
! Set the baseflow bottom level to 200-cm and the height
! above ground for the winds to 10-m and for T/Q to 2-m
! to make the LIS-7 output match the LIS-6 style. - dmm
       if (var_index.eq.LIS_MOC_QSB) toplevtemp(1:1) = 0.0
       if (var_index.eq.LIS_MOC_QSB) botlevtemp(1:1) = 200.0
       if (var_index.eq.LIS_MOC_WINDFORC) botlevtemp(1:1) = 10.0
       if (var_index.eq.LIS_MOC_TAIRFORC) botlevtemp(1:1) = 2.0
       if (var_index.eq.LIS_MOC_QAIRFORC) botlevtemp(1:1) = 2.0
       if (var_index.eq.LIS_MOC_RHMIN) botlevtemp(1:1) = 2.0
#endif

       do t=1,LIS_rc%ntiles(n)
          m = LIS_domain(n)%tile(t)%sftype
          do k=1,dataEntry%vlevels
             if(dataEntry%count(t,k).gt.0) then 
                if(dataEntry%timeAvgOpt.eq.3) then  !do nothing
                   continue   
                elseif(dataEntry%timeAvgOpt.eq.2.or.dataEntry%timeAvgOpt.eq.1) then 
                   dataEntry%modelOutput(1,t,k) = dataEntry%modelOutput(1,t,k)/&
                        dataEntry%count(t,k)
                else !do nothing
                   continue   
                endif
             else
                dataEntry%modelOutput(1,t,k) = LIS_rc%udef
             endif
          enddo
       enddo

       do k=1,dataEntry%vlevels
          ! accumulated values
          
          if(dataEntry%timeAvgOpt.eq.3) then
             stepType = 'accum'
#if (defined AFWA_GRIB_CONFIGS)
             if ((var_index.ne.LIS_MOC_QS).and.                     &
                  (var_index.ne.LIS_MOC_QSB).and.                    &
                  (var_index.ne.LIS_MOC_TOTALPRECIP)) then
                timeRange = 4
             endif
#else 
             timeRange = 4
#endif
             call LIS_writevar_grib1(ftn, ftn_stats, n,   &
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%short_name,&
                  dataEntry%varId_def, &
                  dataEntry%gribSF, &
                  dataEntry%gribSfc,&
                  dataEntry%gribLvl,&
                  stepType, & 
                  time_unit, time_past,time_curr,timeRange,&
                  dataEntry%form, &
                  toplevtemp(k:k), botlevtemp(k:k), 1)
             
             ! time-averaged values and instantaneous values
          elseif(dataEntry%timeAvgOpt.eq.2) then  
             stepType = 'avg'
#if (defined AFWA_GRIB_CONFIGS)
             if ((var_index.ne.LIS_MOC_QS).and.                     &
                  (var_index.ne.LIS_MOC_QSB).and.                    &
                  (var_index.ne.LIS_MOC_TOTALPRECIP)) then
                timeRange = 7
             endif
#else
             timeRange = 7
#endif                
             call LIS_writevar_grib1(ftn,ftn_stats,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%short_name,&
                  dataEntry%varId_def, &
                  dataEntry%gribSF, &
                  dataEntry%gribSfc,&
                  dataEntry%gribLvl,&
                  stepType, & 
                  time_unit, time_past,time_curr,timeRange,&
                  dataEntry%form,  &
                  toplevtemp(k:k), botlevtemp(k:k), 1)!dataEntry%vlevels)
             
             stepType = 'instant'
! Set timerange indicator equal to 1 for AFWA's specifications
! for instantaneous output
! to make the LIS-7 output match the LIS-6 style. - dmm
#if (defined AFWA_GRIB_CONFIGS)
             timeRange = 1
#else
             timeRange = 0
#endif                
             call LIS_writevar_grib1(ftn,ftn_stats,n,&
                  dataEntry%modelOutput(2,:,k),&
                  dataEntry%short_name,&
                  dataEntry%varId_def, &
                  dataEntry%gribSF, &
                  dataEntry%gribSfc,&
                  dataEntry%gribLvl,&
                  stepType, & 
                  time_unit, 0,0,timeRange,&
                  dataEntry%form,  &
                  toplevtemp(k:k), botlevtemp(k:k), 1)
             
             ! time-averaged values
          elseif(dataEntry%timeAvgOpt.eq.1) then  
             stepType = 'avg'
#if (defined AFWA_GRIB_CONFIGS)
             if ((var_index.ne.LIS_MOC_QS).and.                     &
                  (var_index.ne.LIS_MOC_QSB).and.                    &
                  (var_index.ne.LIS_MOC_TOTALPRECIP)) then
                timeRange = 7
             endif
#else
             timeRange = 7
#endif                
             call LIS_writevar_grib1(ftn,ftn_stats,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%short_name,&
                  dataEntry%varId_def, &
                  dataEntry%gribSF, &
                  dataEntry%gribSfc,&
                  dataEntry%gribLvl,&
                  stepType, & 
                  time_unit, time_past,time_curr,timeRange,&
                  dataEntry%form,  &
                  toplevtemp(k:k), botlevtemp(k:k), 1)

             ! instantaneous values
             else
#if (defined AFWA_GRIB_CONFIGS)
                if(var_index.ne.LIS_MOC_RHMIN) then 
                   stepType = 'instant'
                   call LIS_writevar_grib1(ftn, ftn_stats, n,   &
                        dataEntry%modelOutput(1,:,k),&
                        dataEntry%short_name,&
                        dataEntry%varId_def, &
                        dataEntry%gribSF, &
                        dataEntry%gribSfc,&
                        dataEntry%gribLvl,&
                        stepType, & 
! Set timerange indicator equal to 1 for AFWA's specifications
! for instantaneous output
! to make the LIS-7 output match the LIS-6 style. - dmm
                        time_unit, 0,0,1,&
                        dataEntry%form, &
                        toplevtemp(k:k), botlevtemp(k:k), 1)
                else
! Set timerange indicator equal to 132 for AFWA's specifications
! and the time intervals match with AFWA for variable RHMIN only
! to make the LIS-7 output match the LIS-6 style. - dmm
                   stepType = 'min'
                   call LIS_writevar_grib1(ftn, ftn_stats, n,   &
                        dataEntry%modelOutput(1,:,k),&
                        dataEntry%short_name,&
                        dataEntry%varId_def, &
                        dataEntry%gribSF, &
                        dataEntry%gribSfc,&
                        dataEntry%gribLvl,&
                        stepType, & 
                        time_unit, time_past,time_curr,132,&
                        dataEntry%form, &
                        toplevtemp(k:k), botlevtemp(k:k), 1)
                endif
#else
                stepType = 'instant'
                call LIS_writevar_grib1(ftn, ftn_stats, n,   &
                     dataEntry%modelOutput(1,:,k),&
                     dataEntry%short_name,&
                     dataEntry%varId_def, &
                     dataEntry%gribSF, &
                     dataEntry%gribSfc,&
                     dataEntry%gribLvl,&
                     stepType, & 
                     time_unit, 0,0,0,&
                     dataEntry%form, &
                     toplevtemp(k:k), botlevtemp(k:k), 1)
#endif
             endif

             if ( dataEntry%minMaxOpt.ne.0 ) then
                if ( var_index.ne.LIS_MOC_TAIRFORC ) then
                   write(LIS_logunit,*) '[ERR]  min/max support is not '//&
                   'implemented in a generic way for GRIB-1 output.'
                   write(LIS_logunit,*) '[ERR]  min/max support is not '//&
                   'supported for variable ',var_index
                else
! Hard-code Tmax and Tmin only because there are kpds values
! for only these temperature variables in GRIB-1 Table 2 - dmm
!  http://www.nco.ncep.noaa.gov/pmb/docs/on388/table2.html
                   botlevtemp(1:1) = 2.0
                   stepType = 'min'
                   call LIS_writevar_grib1(ftn, ftn_stats, n,   &
                        dataEntry%minimum(:,k),&
                        trim(dataEntry%short_name)//'Min ',&
                        16, &
                        dataEntry%gribSF,&
                        105,&
                        2,&
                        stepType, & 
#if (defined AFWA_GRIB_CONFIGS)
! Set timerange values to 132 for AFWA's specifications.
                        time_unit, time_past,time_curr,132,&
#else
                        time_unit, time_past,time_curr,7,&
#endif
                        dataEntry%form, &
                        toplevtemp(k:k), botlevtemp(k:k), 1)
                   stepType = 'max'
                   call LIS_writevar_grib1(ftn, ftn_stats, n,   &
                        dataEntry%maximum(:,k),&
                        trim(dataEntry%short_name)//'Max ',&
                        15, &
                        dataEntry%gribSF,&
                        105,&
                        2,&
                        stepType, & 
#if (defined AFWA_GRIB_CONFIGS)
! Set timerange values to 132 for AFWA's specifications.
                        time_unit, time_past,time_curr,132,&
#else
                        time_unit, time_past,time_curr,7,&
#endif
                        dataEntry%form, &
                        toplevtemp(k:k), botlevtemp(k:k), 1)
                endif

             endif
          enddo
          deallocate(toplevtemp)
          deallocate(botlevtemp)
       endif
     end subroutine writeSingleGrib1Var
!BOP
! !ROUTINE: writeSingleGrib2Var
! \label{writeSingleGrib2Var}
! 
! !INTERFACE: writeSingleGrib2Var
  subroutine writeSingleGrib2Var(ftn,ftn_stats,n,time_unit, time_past, &
       time_curr, nsoillayers,toplev, botlev, nsoillayers2,topleva,botleva, &
       var_index,dataEntry)
! !USES: 
    
    implicit none
! !ARGUMENTS:
    integer,   intent(in)    :: n 
    integer,   intent(in)    :: ftn
    integer,   intent(in)    :: ftn_stats    
    integer,   intent(in)    :: var_index
    integer,   intent(in)    :: time_unit
    integer,   intent(in)    :: time_past
    integer,   intent(in)    :: time_curr
    integer                  :: nsoillayers
    integer                  :: nsoillayers2
    real, dimension(nsoillayers) :: toplev, botlev
    real, dimension(nsoillayers2) :: topleva, botleva
    real, dimension(:),allocatable :: toplevtemp, botlevtemp
    type(LIS_metadataEntry),pointer :: dataEntry
! 
! !DESCRIPTION: 
!  This routine writes a single variable to a GRIB-2 file
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[var\_index] index of the variable
!    \item[dataEntry]
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_writevar\_grib2](\ref{LIS_writevar_grib2})
!     writes a variable into a GRIB-2 formatted file. 
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \end{description}
!EOP
    integer :: i,k,m,t,timeRange, pdTemplate !, depscale
    character*10 :: stepType
    integer :: ierr
    integer, dimension(:),allocatable :: depscale

    if(dataEntry%selectOpt.eq.1) then 

       if ((var_index.eq.LIS_MOC_SOILTEMP).and.                     &
           (nsoillayers .ne. nsoillayers2)) then
            allocate(toplevtemp(nsoillayers2))
            allocate(botlevtemp(nsoillayers2))
            allocate(depscale(nsoillayers2))
            toplevtemp(1:nsoillayers2) = topleva(1:nsoillayers2)
            botlevtemp(1:nsoillayers2) = botleva(1:nsoillayers2)
       else
            allocate(toplevtemp(nsoillayers))
            allocate(botlevtemp(nsoillayers))
            allocate(depscale(nsoillayers))
            toplevtemp(1:nsoillayers) = toplev(1:nsoillayers)
            botlevtemp(1:nsoillayers) = botlev(1:nsoillayers)
       endif
       ! depth are in cm (set in lsmMod) here -> check values are <=255
       ! and change to m if > 255 and set the scaling factor
       ! scale varies in each layer, so thin layers are distinguishable
       ! (depths are truncated to integer with meters)
       depscale = 0
       do k=1,dataEntry%vlevels
        if ( botlevtemp(k).gt.255 ) then
             toplevtemp(k) = toplevtemp(k) * 0.01   ! cm -> m
             botlevtemp(k) = botlevtemp(k) * 0.01   ! cm -> m
             depscale(k) = 0    ! m
        else
             depscale(k) = 2    ! cm 
        endif
       end do

       timeRange = time_past

       do t=1,LIS_rc%ntiles(n)
          m = LIS_domain(n)%tile(t)%sftype
          do k=1,dataEntry%vlevels
             if(dataEntry%count(t,k).gt.0) then 
                if(dataEntry%timeAvgOpt.eq.3) then  !do nothing
                   continue   
                elseif(dataEntry%timeAvgOpt.eq.2.or.dataEntry%timeAvgOpt.eq.1) then 
                   dataEntry%modelOutput(1,t,k) = dataEntry%modelOutput(1,t,k)/&
                        dataEntry%count(t,k)
                else !do nothing
                   continue   
                endif
             else
                dataEntry%modelOutput(1,t,k) = LIS_rc%udef
             endif
          enddo
       enddo

       do k=1,dataEntry%vlevels
          ! accumulated values
          
          if(dataEntry%timeAvgOpt.eq.3) then
             stepType = 'accum'
             pdTemplate = 8
             call LIS_writevar_grib2(ftn, ftn_stats, n,   &
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%short_name,&
                  dataEntry%varId_def, &
                  dataEntry%gribSF, &
                  dataEntry%gribSfc,&
                  dataEntry%gribDis,&
                  dataEntry%gribCat,&
                  stepType, & 
                  time_unit, timeRange, pdTemplate, &
                  dataEntry%form, &
                  toplevtemp(k:k), botlevtemp(k:k), 1, depscale(k))
             
             ! time-averaged values and instantaneous values
          elseif(dataEntry%timeAvgOpt.eq.2) then  
             stepType = 'avg'
             pdTemplate = 8
             call LIS_writevar_grib2(ftn,ftn_stats,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%short_name,&
                  dataEntry%varId_def, &
                  dataEntry%gribSF, &
                  dataEntry%gribSfc,&
                  dataEntry%gribDis,&
                  dataEntry%gribCat,&
                  stepType, & 
                  time_unit, timeRange, pdTemplate,&
                  dataEntry%form,  &
                  toplevtemp(k:k), botlevtemp(k:k), 1, depscale(k))!dataEntry%vlevels)
             
             stepType = 'instant'
             pdTemplate = 0
             call LIS_writevar_grib2(ftn,ftn_stats,n,&
                  dataEntry%modelOutput(2,:,k),&
                  dataEntry%short_name,&
                  dataEntry%varId_def, &
                  dataEntry%gribSF, &
                  dataEntry%gribSfc,&
                  dataEntry%gribDis,&
                  dataEntry%gribCat,&
                  stepType, & 
                  time_unit, timeRange, pdTemplate,&
                  dataEntry%form,  &
                  toplevtemp(k:k), botlevtemp(k:k), 1, depscale(k))
             
             ! time-averaged values
          elseif(dataEntry%timeAvgOpt.eq.1) then  
             stepType = 'avg'
             pdTemplate = 8
             call LIS_writevar_grib2(ftn,ftn_stats,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%short_name,&
                  dataEntry%varId_def, &
                  dataEntry%gribSF, &
                  dataEntry%gribSfc,&
                  dataEntry%gribDis,&
                  dataEntry%gribCat,&
                  stepType, & 
                  time_unit, timeRange, pdTemplate, &
                  dataEntry%form,  &
                  toplevtemp(k:k), botlevtemp(k:k), 1, depscale(k))

             ! instantaneous values
             else
                stepType = 'instant'
                pdTemplate = 0
                call LIS_writevar_grib2(ftn, ftn_stats, n,   &
                     dataEntry%modelOutput(1,:,k),&
                     dataEntry%short_name,&
                     dataEntry%varId_def, &
                     dataEntry%gribSF, &
                     dataEntry%gribSfc,&
                     dataEntry%gribDis,&
                     dataEntry%gribCat,&
                     stepType, & 
                     time_unit, timeRange, pdTemplate, &
                     dataEntry%form, &
                     toplevtemp(k:k), botlevtemp(k:k), 1, depscale(k))
             endif

             if ( dataEntry%minMaxOpt.ne.0 ) then
! GRIB2 supports generic min/max for any parameter, however, the variable
! needs to be computed in model or template main routines (not yet 
! implemented). -hkb
                   write(LIS_logunit,*) 'Warning: min/max support is not '//&
                   'implemented in a generic way in LIS for ',var_index
                   stepType = 'min'
                   pdTemplate = 8
                   call LIS_writevar_grib2(ftn,ftn_stats,n,&
                        dataEntry%minimum(:,k),&
                        trim(dataEntry%short_name)//'Min ',&
                        dataEntry%varId_def, &
                        dataEntry%gribSF, &
                        dataEntry%gribSfc,&
                        dataEntry%gribDis,&
                        dataEntry%gribCat,&
                        stepType, & 
                        time_unit, timeRange, pdTemplate, &
                        dataEntry%form,  &
                        toplevtemp(k:k), botlevtemp(k:k), 1, depscale(k))
                   stepType = 'max'
                   pdTemplate = 8
                   call LIS_writevar_grib2(ftn, ftn_stats, n,   &
                        dataEntry%maximum(:,k),&
                        trim(dataEntry%short_name)//'Max ',&
                        dataEntry%varId_def, &
                        dataEntry%gribSF, &
                        dataEntry%gribSfc,&
                        dataEntry%gribDis,&
                        dataEntry%gribCat,&
                        stepType, & 
                        time_unit, timeRange, pdTemplate, &
                        dataEntry%form,  &
                        toplevtemp(k:k), botlevtemp(k:k), 1, depscale(k))

             endif
          enddo
          deallocate(toplevtemp)
          deallocate(botlevtemp)
          deallocate(depscale)
       endif
     end subroutine writeSingleGrib2Var
!BOP
! !ROUTINE: writeNetcdfOutput
! \label{writeNetcdfOutput}
! 
! !INTERFACE: writeNetcdfOutput
  subroutine writeNetcdfOutput(n, group, ftn, ftn_stats, outInterval, &
       nsoillayers, lyrthk, model_name)
! !USES: 

! !ARGUMENTS: 
    integer,   intent(in)   :: n 
    integer,   intent(in)   :: group
    integer,   intent(in)   :: ftn
    integer,   intent(in)   :: ftn_stats
    real,      intent(in)   :: outInterval
    integer,   intent(in)   :: nsoillayers
    real,      intent(in)   :: lyrthk(nsoillayers)
    character*100, intent(in) :: model_name
! 
! !DESCRIPTION: 
!  This routine writes an output file in the NETCDF format based on the 
!  list of selected output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[outInterval]   history output frequency
!    \item[nsoillayers]  Number of soil layers
!    \item[lyrthk]   Thickness of soil layers
!    \item[model\_name] Name of the model that generates the output
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[defineNETCDFheadervar](\ref{defineNETCDFheaderVar})
!     writes the required headers for a single variable
!   \item[writeSingleNETCDFvar](\ref{writeSingleNETCDFvar})
!     writes a single variable into a netcdf formatted file. 
!   \item[LIS\_verify](\ref{LIS_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP

    integer                 :: dimID(4)
    integer                 :: tdimID,xtimeID,ensID
    integer                 :: t,c,r,i,index1
    real, allocatable       :: ensval(:) 
    type(LIS_metadataEntry), pointer :: xlat, xlong
    character*8             :: xtime_begin_date
    character*6             :: xtime_begin_time
    character*50            :: xtime_units
    character*50            :: xtime_timeInc
    integer                 :: iret
! Note that the fix to add lat/lon to the NETCDF output will output
! undefined values for the water points. 
    character(len=8)        :: date
    character(len=10)       :: time
    character(len=5)        :: zone
    integer, dimension(8)   :: values
    type(LIS_metadataEntry), pointer :: dataEntry

#if (defined USE_NETCDF3 || defined USE_NETCDF4)           
    call date_and_time(date,time,zone,values)
    
    allocate(xlat)
    allocate(xlong)

    xlat%short_name = "lat"
    xlat%long_name = "latitude"
    xlat%standard_name = "latitude"
    xlat%units = "degree_north"
    xlat%nunits = 1
    xlat%format = 'F'
    xlat%form = 1
    xlat%vlevels = 1
    xlat%timeAvgOpt = 0 
    xlat%selectOpt = 1
    xlat%minMaxOpt = 0 
    xlat%stdOpt = 0 
    allocate(xlat%modelOutput(1,LIS_rc%ntiles(n),xlat%vlevels))
    allocate(xlat%count(1,xlat%vlevels))
    xlat%count = 1
    allocate(xlat%unittypes(1))
    xlat%unittypes(1) = "degree_north"
    xlat%valid_min = 0.0
    xlat%valid_max = 0.0

    xlong%short_name = "lon"
    xlong%long_name = "longitude"
    xlong%standard_name = "longitude"
    xlong%units = "degree_east"
    xlong%nunits = 1
    xlong%format = 'F'
    xlong%form = 1
    xlong%vlevels = 1
    xlong%timeAvgOpt = 0 
    xlong%selectOpt = 1
    xlong%minMaxOpt = 0 
    xlong%stdOpt = 0 
    allocate(xlong%modelOutput(1,LIS_rc%ntiles(n),xlong%vlevels))
    allocate(xlong%count(1,xlong%vlevels))
    xlong%count = 1
    allocate(xlong%unittypes(1))
    xlong%unittypes(1) = "degree_east"
    xlong%valid_min = 0.0
    xlong%valid_max = 0.0

    if(LIS_masterproc) then 
       if(LIS_rc%wopt.eq."1d tilespace") then 
          call LIS_verify(nf90_def_dim(ftn,'ntiles',LIS_rc%glbntiles(n),&
               dimID(1)),&
               'nf90_def_dim for ntiles failed in LIS_historyMod')

       elseif(LIS_rc%wopt.eq."2d gridspace") then 
          call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
               dimID(1)),&
               'nf90_def_dim for east_west failed in LIS_historyMod')
          call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
               dimID(2)),&
               'nf90_def_dim for north_south failed in LIS_historyMod')

       elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
          call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
               dimID(1)),&
               'nf90_def_dim for east_west failed in LIS_historyMod')
          call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
               dimID(2)),&
               'nf90_def_dim for north_south failed in LIS_historyMod')
          call LIS_verify(nf90_def_dim(ftn,'ensemble',LIS_rc%nensem(n),&
               dimID(3)),&
               'nf90_def_dim for ensemble failed in LIS_historyMod')
       endif

       ! LIS output is always writing output for a single time
       call LIS_verify(nf90_def_dim(ftn,'time',NF90_UNLIMITED,tdimID),&
            'nf90_def_dim for time failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value", LIS_rc%udef),&
            'nf90_put_att for missing_value failed in LIS_historyMod')
       
       call defineNETCDFheaderVar(n,ftn,dimID, xlat,&
            non_model_fields = 1 )       
       call defineNETCDFheaderVar(n,ftn,dimID, xlong, &
            non_model_fields = 2 )              

       ! defining time field
       call LIS_verify(nf90_def_var(ftn,'time',&
            nf90_float,dimids = tdimID, varID=xtimeID),&
            'nf90_def_var for time failed in LIS_historyMod')
       
       write(xtime_units,200) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
            LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
200    format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
            I2.2,':',I2.2)
       write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
            LIS_rc%yr, LIS_rc%mo, LIS_rc%da
       write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
            LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
       write(xtime_timeInc, fmt='(I20)') &
            nint(outInterval)
       ! time field attributes
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "units",trim(xtime_units)),&
            'nf90_put_att for units failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "long_name","time"),&
            'nf90_put_att for long_name failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "time_increment",trim(adjustl(xtime_timeInc))),&
            'nf90_put_att for time_increment failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "begin_date",xtime_begin_date),&
            'nf90_put_att for begin_date failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "begin_time",xtime_begin_time),&
            'nf90_put_att for begin_time failed in LIS_historyMod')

       ! Write ensemble information as a variable:
       if( LIS_rc%wopt.eq."2d ensemble gridspace" ) then
          call LIS_verify(nf90_def_var(ftn,'ensemble',&
               nf90_float, dimids=dimID(3), varID=ensID),&
              'nf90_def_var for ensemble failed in LIS_historyMod')
          ! ensemble var attributes
          call LIS_verify(nf90_put_att(ftn,ensID,&
              "units","ensemble number"),&
              'nf90_put_att for ensemble units failed in LIS_historyMod')
          call LIS_verify(nf90_put_att(ftn,ensID,&
              "long_name","Ensemble numbers"),&
              'nf90_put_att for ensemble long_name failed in LIS_historyMod')
          allocate(ensval(LIS_rc%nensem(n)))
          do i = 1, LIS_rc%nensem(n) 
             ensval(i) = float(i)
          end do
       endif

       ! Pointer to header information
       if(group.eq.1) then     ! LSM output
          dataEntry => LIS_histData(n)%head_lsm_list
       elseif(group.eq.2) then ! ROUTING
          dataEntry => LIS_histData(n)%head_routing_list
       elseif(group.eq.3) then ! RTM
          dataEntry => LIS_histData(n)%head_rtm_list
       elseif(group.eq.4) then ! Irrigation
          dataEntry => LIS_histData(n)%head_irrig_list
       endif

       do while ( associated(dataEntry) )
          call defineNETCDFheaderVar(n,ftn,dimId,&
               dataEntry)
          dataEntry => dataEntry%next
       enddo

       ! Global attributes
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUM_SOIL_LAYERS", &
            nsoillayers),&
            'nf90_put_att for title failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOIL_LAYER_THICKNESSES", &
            lyrthk),&
            'nf90_put_att for title failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
            "LIS land surface model output"),&
            'nf90_put_att for title failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
            trim(LIS_rc%institution)),&
            'nf90_put_att for institution failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"source",&
            trim(model_name)),&
            'nf90_put_att for source failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
            "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
            date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)),&
            'nf90_put_att for history failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
            "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"),&
            'nf90_put_att for references failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"conventions", &
            "CF-1.6"),'nf90_put_att for conventions failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
            "website: http://lis.gsfc.nasa.gov/"),&
            'nf90_put_att for comment failed in LIS_historyMod')

       ! Grid information
       if(LIS_rc%lis_map_proj.eq."latlon") then   ! latlon
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "EQUIDISTANT CYLINDRICAL"))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,9)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,10)))       
       elseif(LIS_rc%lis_map_proj.eq."mercator") then 
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "MERCATOR"))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LIS_rc%gridDesc(n,10)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LIS_rc%gridDesc(n,11)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,8)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,9)))
       elseif(LIS_rc%lis_map_proj.eq."lambert") then ! lambert conformal
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "LAMBERT CONFORMAL"))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LIS_rc%gridDesc(n,10)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
               LIS_rc%gridDesc(n,7)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LIS_rc%gridDesc(n,11)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,8)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,9)))

       elseif(LIS_rc%lis_map_proj.eq."polar") then ! polar stereographic
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "POLAR STEREOGRAPHIC"))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LIS_rc%gridDesc(n,10)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
               LIS_rc%gridDesc(n,7)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LIS_rc%gridDesc(n,11)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,8)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,9)))
       endif       
       call LIS_verify(nf90_enddef(ftn))
    endif

    ! Write variable data (both non-model and model-based):
    do t=1,LIS_rc%ntiles(n)
       c = LIS_domain(n)%tile(t)%col
       r = LIS_domain(n)%tile(t)%row
       index1 = LIS_domain(n)%gindex(c,r)
       xlat%modelOutput(1,t,1) = LIS_domain(n)%grid(index1)%lat
       xlong%modelOutput(1,t,1) = LIS_domain(n)%grid(index1)%lon
    enddo
    call writeSingleNETCDFvar(ftn,ftn_stats,n,xlat,non_model_fields =1)
    call writeSingleNETCDFvar(ftn,ftn_stats,n,xlong, non_model_fields = 2)

    if(LIS_masterproc) then 
       call LIS_verify(nf90_put_var(ftn,xtimeID,0.0),&
            'nf90_put_var for xtimeID failed in LIS_historyMod')
       ! Write ensemble var info:
       if( LIS_rc%wopt.eq."2d ensemble gridspace" ) then
         call LIS_verify(nf90_put_var(ftn,ensID,ensval,(/1/),(/LIS_rc%nensem(n)/)),&
            'nf90_put_var for ensID failed in LIS_historyMod')
       endif
    endif

    if(group.eq.1) then     ! LSM output
       dataEntry => LIS_histData(n)%head_lsm_list
    elseif(group.eq.2) then ! ROUTING
       dataEntry => LIS_histData(n)%head_routing_list
    elseif(group.eq.3) then ! RTM
       dataEntry => LIS_histData(n)%head_rtm_list
    elseif(group.eq.4) then ! Irrigation
       dataEntry => LIS_histData(n)%head_irrig_list
    endif

    do while ( associated(dataEntry) )
       call writeSingleNETCDFvar(ftn,ftn_stats,n,&
                                 dataEntry)
       dataEntry => dataEntry%next
    enddo
    deallocate(xlat%modelOutput)
    deallocate(xlat%count)
    deallocate(xlat%unittypes)
    deallocate(xlat)

    deallocate(xlong%modelOutput)
    deallocate(xlong%count)
    deallocate(xlong%unittypes)
    deallocate(xlong)
    
    if( LIS_masterproc ) then
      if( LIS_rc%wopt.eq."2d ensemble gridspace" ) then
        deallocate(ensval)
      endif
    endif

#endif
  end subroutine writeNetcdfOutput

!BOP
! !ROUTINE: writeRoutingNetcdfOutput
! \label{writeRoutingNetcdfOutput}
! 
! !INTERFACE: writeRoutingNetcdfOutput
  subroutine writeRoutingNetcdfOutput(n, group, ftn, ftn_stats, outInterval, &
       nsoillayers, lyrthk, model_name)
! !USES: 

! !ARGUMENTS: 
    integer,   intent(in)   :: n 
    integer,   intent(in)   :: group
    integer,   intent(in)   :: ftn
    integer,   intent(in)   :: ftn_stats
    real,      intent(in)   :: outInterval
    integer,   intent(in)   :: nsoillayers
    real,      intent(in)   :: lyrthk(nsoillayers)
    character*100, intent(in) :: model_name
! 
! !DESCRIPTION: 
!  This routine writes an output file in the NETCDF format based on the 
!  list of selected output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[outInterval]   history output frequency
!    \item[nsoillayers]  Number of soil layers
!    \item[lyrthk]   Thickness of soil layers
!    \item[model\_name] Name of the model that generates the output
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[defineNETCDFheadervar](\ref{defineNETCDFheaderVar})
!     writes the required headers for a single variable
!   \item[writeSingleNETCDFvar](\ref{writeSingleNETCDFvar})
!     writes a single variable into a netcdf formatted file. 
!   \item[LIS\_verify](\ref{LIS_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP

    integer                 :: dimID(4)
    integer                 :: tdimID,xtimeID,ensID
    integer                 :: t,c,r,m,i,index1
    real, allocatable       :: ensval(:) 
    type(LIS_metadataEntry), pointer :: xlat, xlong
    character*8             :: xtime_begin_date
    character*6             :: xtime_begin_time
    character*50            :: xtime_units
    character*50            :: xtime_timeInc
    integer                 :: iret
! Note that the fix to add lat/lon to the NETCDF output will output
! undefined values for the water points. 
    character(len=8)        :: date
    character(len=10)       :: time
    character(len=5)        :: zone
    integer, dimension(8)   :: values
    type(LIS_metadataEntry), pointer :: dataEntry

#if (defined USE_NETCDF3 || defined USE_NETCDF4)           
    call date_and_time(date,time,zone,values)
    
    allocate(xlat)
    allocate(xlong)

    xlat%short_name = "lat"
    xlat%long_name = "latitude"
    xlat%standard_name = "latitude"
    xlat%units = "degree_north"
    xlat%nunits = 1
    xlat%format = 'F'
    xlat%form = 1
    xlat%vlevels = 1
    xlat%timeAvgOpt = 0 
    xlat%selectOpt = 1
    xlat%minMaxOpt = 0 
    xlat%stdOpt = 0 
    allocate(xlat%modelOutput(1,LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n),&
         xlat%vlevels))
    allocate(xlat%count(1,xlat%vlevels))
    xlat%count = 1
    allocate(xlat%unittypes(1))
    xlat%unittypes(1) = "degree_north"
    xlat%valid_min = 0.0
    xlat%valid_max = 0.0

    xlong%short_name = "lon"
    xlong%long_name = "longitude"
    xlong%standard_name = "longitude"
    xlong%units = "degree_east"
    xlong%nunits = 1
    xlong%format = 'F'
    xlong%form = 1
    xlong%vlevels = 1
    xlong%timeAvgOpt = 0 
    xlong%selectOpt = 1
    xlong%minMaxOpt = 0 
    xlong%stdOpt = 0 
    allocate(xlong%modelOutput(1,LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n),&
         xlong%vlevels))
    allocate(xlong%count(1,xlong%vlevels))
    xlong%count = 1
    allocate(xlong%unittypes(1))
    xlong%unittypes(1) = "degree_east"
    xlong%valid_min = 0.0
    xlong%valid_max = 0.0

    if(LIS_masterproc) then 
       if(LIS_rc%wopt.eq."1d tilespace") then 
          write(LIS_logunit,*) '[ERR] 1d tilespace output for routing models'
          write(LIS_logunit,*) '[ERR] is not supported currently'
          call LIS_endrun()

       elseif(LIS_rc%wopt.eq."2d gridspace") then 
          call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
               dimID(1)),&
               'nf90_def_dim for east_west failed in LIS_historyMod')
          call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
               dimID(2)),&
               'nf90_def_dim for north_south failed in LIS_historyMod')

       elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
          call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
               dimID(1)),&
               'nf90_def_dim for east_west failed in LIS_historyMod')
          call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
               dimID(2)),&
               'nf90_def_dim for north_south failed in LIS_historyMod')
          call LIS_verify(nf90_def_dim(ftn,'ensemble',LIS_rc%nensem(n),&
               dimID(3)),&
               'nf90_def_dim for ensemble failed in LIS_historyMod')
       endif

       ! LIS output is always writing output for a single time
       call LIS_verify(nf90_def_dim(ftn,'time',NF90_UNLIMITED,tdimID),&
            'nf90_def_dim for time failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value", LIS_rc%udef),&
            'nf90_put_att for missing_value failed in LIS_historyMod')
       
       call defineNETCDFheaderVar(n,ftn,dimID, xlat,&
            non_model_fields = 1 )       
       call defineNETCDFheaderVar(n,ftn,dimID, xlong, &
            non_model_fields = 2 )              

       ! defining time field
       call LIS_verify(nf90_def_var(ftn,'time',&
            nf90_float,dimids = tdimID, varID=xtimeID),&
            'nf90_def_var for time failed in LIS_historyMod')
       
       write(xtime_units,200) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
            LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
200    format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
            I2.2,':',I2.2)
       write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
            LIS_rc%yr, LIS_rc%mo, LIS_rc%da
       write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
            LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
       write(xtime_timeInc, fmt='(I20)') &
            nint(outInterval)
       ! time field attributes
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "units",trim(xtime_units)),&
            'nf90_put_att for units failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "long_name","time"),&
            'nf90_put_att for long_name failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "time_increment",trim(adjustl(xtime_timeInc))),&
            'nf90_put_att for time_increment failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "begin_date",xtime_begin_date),&
            'nf90_put_att for begin_date failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "begin_time",xtime_begin_time),&
            'nf90_put_att for begin_time failed in LIS_historyMod')

       ! Write ensemble information as a variable:
       if( LIS_rc%wopt.eq."2d ensemble gridspace" ) then
          call LIS_verify(nf90_def_var(ftn,'ensemble',&
               nf90_float, dimids=dimID(3), varID=ensID),&
              'nf90_def_var for ensemble failed in LIS_historyMod')
          ! ensemble var attributes
          call LIS_verify(nf90_put_att(ftn,ensID,&
              "units","ensemble number"),&
              'nf90_put_att for ensemble units failed in LIS_historyMod')
          call LIS_verify(nf90_put_att(ftn,ensID,&
              "long_name","Ensemble numbers"),&
              'nf90_put_att for ensemble long_name failed in LIS_historyMod')
          allocate(ensval(LIS_rc%nensem(n)))
          do i = 1, LIS_rc%nensem(n) 
             ensval(i) = float(i)
          end do
       endif

       ! Pointer to header information
       dataEntry => LIS_histData(n)%head_routing_list

       do while ( associated(dataEntry) )
          call defineNETCDFheaderVar(n,ftn,dimId,&
               dataEntry)
          dataEntry => dataEntry%next
       enddo

       ! Global attributes
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUM_SOIL_LAYERS", &
            nsoillayers),&
            'nf90_put_att for title failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOIL_LAYER_THICKNESSES", &
            lyrthk),&
            'nf90_put_att for title failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
            "LIS land surface model output"),&
            'nf90_put_att for title failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
            trim(LIS_rc%institution)),&
            'nf90_put_att for institution failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"source",&
            trim(model_name)),&
            'nf90_put_att for source failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
            "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
            date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)),&
            'nf90_put_att for history failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
            "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"),&
            'nf90_put_att for references failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"conventions", &
            "CF-1.6"),'nf90_put_att for conventions failed in LIS_historyMod')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
            "website: http://lis.gsfc.nasa.gov/"),&
            'nf90_put_att for comment failed in LIS_historyMod')

       ! Grid information
       if(LIS_rc%lis_map_proj.eq."latlon") then   ! latlon
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "EQUIDISTANT CYLINDRICAL"))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,9)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,10)))       
       elseif(LIS_rc%lis_map_proj.eq."mercator") then 
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "MERCATOR"))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LIS_rc%gridDesc(n,10)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LIS_rc%gridDesc(n,11)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,8)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,9)))
       elseif(LIS_rc%lis_map_proj.eq."lambert") then ! lambert conformal
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "LAMBERT CONFORMAL"))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LIS_rc%gridDesc(n,10)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
               LIS_rc%gridDesc(n,7)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LIS_rc%gridDesc(n,11)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,8)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,9)))

       elseif(LIS_rc%lis_map_proj.eq."polar") then ! polar stereographic
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "POLAR STEREOGRAPHIC"))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LIS_rc%gridDesc(n,10)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
               LIS_rc%gridDesc(n,7)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LIS_rc%gridDesc(n,11)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,8)))
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,9)))
       endif       
       call LIS_verify(nf90_enddef(ftn))
    endif

    ! Write variable data (both non-model and model-based):
    do i=1,LIS_rc%nroutinggrid(n)
       do m=1,LIS_rc%nensem(n)
          t = m+(i-1)*LIS_rc%nensem(n)

          c = LIS_routing(n)%tile(t)%col
          r = LIS_routing(n)%tile(t)%row

          xlat%modelOutput(1,t,1) = LIS_domain(n)%lat(c+(r-1)*LIS_rc%lnc(n))
          xlong%modelOutput(1,t,1) = LIS_domain(n)%lon(c+(r-1)*LIS_rc%lnc(n))
       enddo
    enddo

    call writeSingleRoutingNETCDFvar(ftn,ftn_stats,n,xlat,&
         non_model_fields = 1)
    call writeSingleRoutingNETCDFvar(ftn,ftn_stats,n,xlong, &
         non_model_fields = 2)

    if(LIS_masterproc) then 
       call LIS_verify(nf90_put_var(ftn,xtimeID,0.0),&
            'nf90_put_var for xtimeID failed in LIS_historyMod')
       ! Write ensemble var info:
       if( LIS_rc%wopt.eq."2d ensemble gridspace" ) then
         call LIS_verify(nf90_put_var(ftn,ensID,ensval,(/1/),(/LIS_rc%nensem(n)/)),&
            'nf90_put_var for ensID failed in LIS_historyMod')
       endif
    endif

    dataEntry => LIS_histData(n)%head_routing_list

    do while ( associated(dataEntry) )
       call writeSingleRoutingNETCDFvar(ftn,ftn_stats,n,&
                                 dataEntry)
       dataEntry => dataEntry%next
    enddo
    deallocate(xlat%modelOutput)
    deallocate(xlat%count)
    deallocate(xlat%unittypes)
    deallocate(xlat)

    deallocate(xlong%modelOutput)
    deallocate(xlong%count)
    deallocate(xlong%unittypes)
    deallocate(xlong)
    
    if( LIS_masterproc ) then
      if( LIS_rc%wopt.eq."2d ensemble gridspace" ) then
        deallocate(ensval)
      endif
    endif

#endif
  end subroutine writeRoutingNetcdfOutput

!BOP
! !ROUTINE: defineNETCDFheaderVar
! \label{defineNETCDFheaderVar}
! 
! !INTERFACE: 
  subroutine defineNETCDFheaderVar(n,ftn,dimID, dataEntry, non_model_fields)
! !USES: 

! !ARGUMENTS:     
    integer                           :: n
    integer                           :: ftn
    type(LIS_metadataEntry), pointer  :: dataEntry
    integer,   optional               :: non_model_fields
    integer                           :: dimID(4)
! 
! !DESCRIPTION: 
!    This routine writes the required NETCDF header for a single variable
! 
!   The arguments are: 
!   \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    NETCDF file unit handle
!   \item[dimID]
!    NETCDF dimension ID corresponding to the variable
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \item[LIS\_verify](\ref{LIS_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP
    integer       :: nmodel_status
    integer       :: data_index
    integer       :: shuffle, deflate, deflate_level
    character*100 :: short_name
    integer       :: fill_value

    ! EMK FIXME...This subroutine should be refactored to specify the correct length of the dimID array to
    ! pass to the netCDF library.  Once set based on output methodology and number of vertical levels, the
    ! logic for specifying inst, tavg, acc, min, max, etc variables should be straightforward and can eliminate
    ! much redundant code below.

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    nmodel_status = 0
    if(present(non_model_fields)) then 
       nmodel_status = non_model_fields
    endif

    data_index = dataEntry%index

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    if(dataEntry%selectOpt.eq.1)then 
       if(LIS_rc%wopt.eq."1d tilespace") then 
          if(dataEntry%vlevels.gt.1) then 
             call LIS_verify(nf90_def_dim(ftn,&
                  trim(dataEntry%short_name)//'_profiles',&
                  dataEntry%vlevels, dimID(2)),&
                  'nf90_def_dim failed (1d tilespace) in LIS_historyMod')
          endif
       elseif(LIS_rc%wopt.eq."2d gridspace") then 
          if(dataEntry%vlevels.gt.1) then 
             call LIS_verify(nf90_def_dim(ftn,&
                  trim(dataEntry%short_name)//'_profiles',&
                  dataEntry%vlevels, dimID(3)),&
                  'nf90_def_dim failed (2d gridspace) in LIS_historyMod')
          endif
       elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
          if(dataEntry%vlevels.gt.1) then 
             call LIS_verify(nf90_def_dim(ftn,&
                  trim(dataEntry%short_name)//'_profiles',&
                  dataEntry%vlevels, dimID(4)),&
                  'nf90_def_dim failed (2d ensemble gridspace) in LIS_historyMod')
          endif
       endif
       if(LIS_rc%wopt.eq."1d tilespace") then              
          if(dataEntry%timeAvgOpt.eq.2) then 
             ! EMK...Added support for Min/Max
             !if(dataEntry%minMaxOpt.gt.1) then 
             !   write(LIS_logunit,*) '[ERR] Enabling Min/Max option is not supported'
             !   write(LIS_logunit,*) '[ERR] when the time averaging option is set to 2'
             !   call LIS_endrun()
             !endif

             if(dataEntry%vlevels.gt.1) then 
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_tavg',&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')

#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
                
                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
#endif
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_inst',&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_opt1),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')

#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_opt1, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_opt1,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
#endif

!EMK...Add max/min for Air Force
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID(1:2), varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                   
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID(1:2), varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")
                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif

!EMK END

             else ! single level
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_tavg',&
                     nf90_float,&
                     dimids = dimID(1:1), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_inst',&
                     nf90_float,&
                     dimids = dimID(1:1), varID=dataEntry%varId_opt1),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_opt1, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_opt1,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif

!EMK...Add max/min for Air Force
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID(1:1), varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                   
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID(1:1), varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")
                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif

!EMK END

             endif
          else
             if(nmodel_status==0) then 
                if(dataEntry%timeAvgOpt.eq.0) then
                   short_name = trim(dataEntry%short_name)//'_inst'
                elseif(dataEntry%timeAvgOpt.eq.1) then
                   short_name = trim(dataEntry%short_name)//'_tavg'
                elseif(dataEntry%timeAvgOpt.eq.3) then
                   short_name = trim(dataEntry%short_name)//'_acc'
                endif
             else
                short_name = trim(dataEntry%short_name)
             endif
             
             if(dataEntry%vlevels.gt.1) then 
                call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
                

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')
#endif

! EMK...Add max/min support for Air Force
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID(1:1), varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif

                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID(1:1), varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif
!EMK END
             else ! 1 vertical level
                call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                     nf90_float,&
                     dimids = dimID(1:1), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID(1:1), varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif

                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID(1:1), varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif
             endif
          endif
       elseif(LIS_rc%wopt.eq."2d gridspace") then 
          if(dataEntry%timeAvgOpt.eq.2) then 
             ! EMK...Added support for Max/Min
             !if(dataEntry%minMaxOpt.gt.1) then 
             !   write(LIS_logunit,*) '[ERR] Enabling Min/Max option is not supported'
             !   write(LIS_logunit,*) '[ERR] when the time averaging option is set to 2'
             !   call LIS_endrun()
             !endif
             if(dataEntry%vlevels.gt.1) then 
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_tavg',&
                     nf90_float,&
                     dimids = dimID(1:3), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_inst',&
                     nf90_float,&
                     dimids = dimID(1:3), varID=dataEntry%varId_opt1),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_opt1, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
                
                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_opt1,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif

!EMK...Add max/min for Air Force
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID(1:3), varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                   
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID(1:3), varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")
                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif

!EMK END
             else ! 1 vertical level
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_tavg',&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#if(defined USE_NETCDF4)

                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_inst',&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_opt1),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_opt1, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_opt1,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                

!EMK...Add max/min for Air Force
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID(1:2), varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                   
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID(1:2), varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")
                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif

!EMK END

             endif
          else
             if(nmodel_status==0) then 
                if(dataEntry%timeAvgOpt.eq.0) then
                   short_name = trim(dataEntry%short_name)//'_inst'
                elseif(dataEntry%timeAvgOpt.eq.1) then
                   short_name = trim(dataEntry%short_name)//'_tavg'
                elseif(dataEntry%timeAvgOpt.eq.3) then
                   short_name = trim(dataEntry%short_name)//'_acc'
                endif
             else
                short_name = trim(dataEntry%short_name)
             endif
             
             if(dataEntry%vlevels.gt.1) then 
                call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                     nf90_float,&
                     dimids = dimID(1:3), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID(1:3), varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(short_name)//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                   
                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID(1:3), varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(short_name)//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif
             else ! 1 vertical level
                ! lat/lon fields will write in 1D  
                if(LIS_rc%nlatlon_dimensions == '1D') then
                   if(nmodel_status.eq.1) then
                      call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                           nf90_float,&
                           dimids = dimID(2), varID=dataEntry%varId_def),&
                           'nf90_def_var for '//trim(short_name)//&
                           'failed in defineNETCDFheadervar')                     
                   elseif(nmodel_status.eq.2) then
                      call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                           nf90_float,&
                           dimids = dimID(1), varID=dataEntry%varId_def),&
                           'nf90_def_var for '//trim(short_name)//&
                           'failed in defineNETCDFheadervar')                     
                   else
                      call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                           nf90_float,&
                           dimids = dimID(1:2), varID=dataEntry%varId_def),&
                           'nf90_def_var for '//trim(short_name)//&
                           'failed in defineNETCDFheadervar')
                   endif
                ! lat/lon fields will write in 2D
                else
                   call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                        nf90_float,&
                        dimids = dimID(1:2), varID=dataEntry%varID_def),&
                        'nf90_def_var for '//trim(short_name)//&
                        'failed in defineNETCDFheadervar')
                endif
               

#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)                

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID(1:2), varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(short_name)//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                   
                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID(1:2), varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(short_name)//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")
                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif
             endif
          endif
       elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 
          if(dataEntry%timeAvgOpt.eq.2) then 
             
             ! EMK...Added support for Max/Min
!             if(dataEntry%minMaxOpt.gt.1) then 
!                write(LIS_logunit,*) '[ERR] Enabling Min/Max option is not supported'
!                write(LIS_logunit,*) '[ERR] when the time averaging option is set to 2'
!                call LIS_endrun()
!             endif

             if(dataEntry%vlevels.gt.1) then 
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_tavg',&
                     nf90_float,&
                     dimids = dimID, varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_inst',&
                     nf90_float,&
                     dimids = dimID, varID=dataEntry%varId_opt1),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_opt1, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
                
                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_opt1,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif

!EMK...Add max/min for Air Force
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID, varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                   
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID, varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")
                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif

!EMK END
             else  ! 1 vertical level
                ! lat/lon field output will write in 1D 
                if(LIS_rc%nlatlon_dimensions == '1D') then
                   if(nmodel_status.eq.1) then
                      call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                           nf90_float,&
                           dimids = dimID(2), varID=dataEntry%varId_def),&
                           'nf90_def_var for '//trim(short_name)//&
                           'failed in defineNETCDFheadervar')                     
                   elseif(nmodel_status.eq.2) then
                      call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                           nf90_float,&
                           dimids = dimID(1), varID=dataEntry%varId_def),&
                           'nf90_def_var for '//trim(short_name)//&
                           'failed in defineNETCDFheadervar')                     
                   else                
                      call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_tavg',&
                           nf90_float,&
                           dimids = dimID(1:3), varID=dataEntry%varId_def),&
                           'nf90_def_var for '//trim(dataEntry%short_name)//&
                           'failed in defineNETCDFheadervar')
                   endif
                ! latlon field output will write in 2D
                else
                     call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                         nf90_float,&
                         dimids = dimID(1:2), varID=dataEntry%varID_def),&
                         'nf90_def_var for '//trim(short_name)//&
                         'failed in defineNETCDFheadervar')
                endif 



#if(defined USE_NETCDF4)

                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                
                call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_inst',&
                     nf90_float,&
                     dimids = dimID(1:3), varID=dataEntry%varId_opt1),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_opt1, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_opt1,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                

!EMK...Add max/min for Air Force
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID(1:3), varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                   
                   call LIS_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID(1:3), varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")
                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(dataEntry%short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif

!EMK END

             endif
          else
             if(nmodel_status==0) then 
                if(dataEntry%timeAvgOpt.eq.0) then
                   short_name = trim(dataEntry%short_name)//'_inst'
                elseif(dataEntry%timeAvgOpt.eq.1) then
                   short_name = trim(dataEntry%short_name)//'_tavg'
                elseif(dataEntry%timeAvgOpt.eq.3) then
                   short_name = trim(dataEntry%short_name)//'_acc'
                endif
             else
                short_name = trim(dataEntry%short_name)
             endif
             
             if(dataEntry%vlevels.gt.1) then 
                call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                     nf90_float,&
                     dimids = dimID, varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif
                if(dataEntry%minMaxOpt.gt.0) then 
                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID, varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(short_name)//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                   
                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID, varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(short_name)//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif
             else
                if(nmodel_status==0) then                 
                   call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                        nf90_float,&
                        dimids = dimID(1:3), varID=dataEntry%varId_def),&
                        'nf90_def_var for '//trim(short_name)//&
                        'failed in defineNETCDFheadervar')                     
                else
                   ! lat/lon field output will write in 1D
                   if(LIS_rc%nlatlon_dimensions == '1D') then
                      if(nmodel_status.eq.1) then
                         call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                              nf90_float,&
                              dimids = dimID(2), varID=dataEntry%varId_def),&
                              'nf90_def_var for '//trim(short_name)//&
                              'failed in defineNETCDFheadervar')                     
                      elseif(nmodel_status.eq.2) then
                         call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                              nf90_float,&
                              dimids = dimID(1), varID=dataEntry%varId_def),&
                              'nf90_def_var for '//trim(short_name)//&
                              'failed in defineNETCDFheadervar')                     
                      else 
                         call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                              nf90_float,&
                              dimids = dimID(1:2), varID=dataEntry%varId_def),&
                              'nf90_def_var for '//trim(short_name)//&
                              'failed in defineNETCDFheadervar')
                      endif
                   ! lat/lon field output will write in 2D
                   else
                      call LIS_verify(nf90_def_var(ftn,trim(short_name),&
                           nf90_float,&
                           dimids = dimID(1:2), varID=dataEntry%varID_def),&
                           'nf90_def_var for '//trim(short_name)//&
                           'failed in defineNETCDFheadervar')
                   endif
                endif
#if(defined USE_NETCDF4)
                call LIS_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)                

                call LIS_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                
                if(dataEntry%minMaxOpt.gt.0) then 
                   ! FIXME...Should dimID always be 1:3 here?
                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_min",&
                        nf90_float,&
                        dimids = dimID(1:3), varID=dataEntry%varId_min),&
                        'nf90_def_var for '//trim(short_name)//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_min, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_min")

                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_min,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_min"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                   
                   call LIS_verify(nf90_def_var(ftn,trim(short_name)//"_max",&
                        nf90_float,&
                        dimids = dimID(1:3), varID=dataEntry%varId_max),&
                        'nf90_def_var for '//trim(short_name)//&
                        'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)                
                   call LIS_verify(nf90_def_var_fill(ftn,&
                        dataEntry%varId_max, &
                        1,fill_value), 'nf90_def_var_fill failed for '//&
                        trim(dataEntry%short_name)//"_max")
                   call LIS_verify(nf90_def_var_deflate(ftn,&
                        dataEntry%varId_max,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(short_name)//"_max"//&
                        'failed in defineNETCDFheadervar')                     
#endif
                endif
             endif
          endif
       endif
       
       call LIS_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "units",trim(dataEntry%units)),&
            'nf90_put_att for units failed in defineNETCDFheaderVar')
       call LIS_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "standard_name",trim(dataEntry%standard_name)),&
            'nf90_put_att for standard_name failed in defineNETCDFheaderVar')
       call LIS_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "long_name",trim(dataEntry%long_name)),&
            'nf90_put_att for long_name failed in defineNETCDFheaderVar')
       call LIS_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "scale_factor",1.0),&
            'nf90_put_att for scale_factor failed in defineNETCDFheaderVar')
       call LIS_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "add_offset",0.0),&
            'nf90_put_att for add_offset failed in defineNETCDFheaderVar')
       call LIS_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "missing_value",LIS_rc%udef),&
            'nf90_put_att for missing_value failed in defineNETCDFheaderVar')
       call LIS_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "_FillValue",LIS_rc%udef),&
            'nf90_put_att for _FillValue failed in defineNETCDFheaderVar')
       call LIS_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "vmin",dataEntry%valid_min),&
            'nf90_put_att for vmin failed in defineNETCDFheaderVar')
       call LIS_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "vmax",dataEntry%valid_max),&
            'nf90_put_att for vmax failed in defineNETCDFheaderVar')

       if(dataEntry%timeAvgOpt.eq.2) then
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_opt1,&
               "units",trim(dataEntry%units)),&
               'nf90_put_att for units failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_opt1,&
               "standard_name",trim(dataEntry%standard_name)),&
               'nf90_put_att for standard_name failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_opt1,&
               "long_name",trim(dataEntry%long_name)),&
               'nf90_put_att for long_name failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_opt1,&
               "scale_factor",1.0),&
               'nf90_put_att for scale_factor failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_opt1,&
               "add_offset",0.0),&
               'nf90_put_att for add_offset failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_opt1,&
               "missing_value",LIS_rc%udef),&
               'nf90_put_att for missing_value failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_opt1,&
               "_FillValue",LIS_rc%udef),&
               'nf90_put_att for _FillValue failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_opt1,&
               "vmin",dataEntry%valid_min),&
               'nf90_put_att for vmin failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_opt1,&
               "vmax",dataEntry%valid_max),&
               'nf90_put_att for vmax failed in defineNETCDFheaderVar')

       endif

       ! EMK...Add metadata for max/min variables
       if(dataEntry%minMaxOpt.gt.0) then

          ! Min metadata
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_min,&
               "units",trim(dataEntry%units)),&
               'nf90_put_att for units failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_min,&
               "standard_name",trim(dataEntry%standard_name)),&
               'nf90_put_att for standard_name failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_min,&
               "long_name",trim(dataEntry%long_name)),&
               'nf90_put_att for long_name failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_min,&
               "scale_factor",1.0),&
               'nf90_put_att for scale_factor failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_min,&
               "add_offset",0.0),&
               'nf90_put_att for add_offset failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_min,&
               "missing_value",LIS_rc%udef),&
               'nf90_put_att for missing_value failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_min,&
               "_FillValue",LIS_rc%udef),&
               'nf90_put_att for _FillValue failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_min,&
               "vmin",dataEntry%valid_min),&
               'nf90_put_att for vmin failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_min,&
               "vmax",dataEntry%valid_max),&
               'nf90_put_att for vmax failed in defineNETCDFheaderVar')

          ! Max metadata
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_max,&
               "units",trim(dataEntry%units)),&
               'nf90_put_att for units failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_max,&
               "standard_name",trim(dataEntry%standard_name)),&
               'nf90_put_att for standard_name failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_max,&
               "long_name",trim(dataEntry%long_name)),&
               'nf90_put_att for long_name failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_max,&
               "scale_factor",1.0),&
               'nf90_put_att for scale_factor failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_max,&
               "add_offset",0.0),&
               'nf90_put_att for add_offset failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_max,&
               "missing_value",LIS_rc%udef),&
               'nf90_put_att for missing_value failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_max,&
               "_FillValue",LIS_rc%udef),&
               'nf90_put_att for _FillValue failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_max,&
               "vmin",dataEntry%valid_min),&
               'nf90_put_att for vmin failed in defineNETCDFheaderVar')
          call LIS_verify(nf90_put_att(ftn,dataEntry%varId_max,&
               "vmax",dataEntry%valid_max),&
               'nf90_put_att for vmax failed in defineNETCDFheaderVar')
       endif

    endif
#endif
  end subroutine defineNETCDFheaderVar

!BOP
! !ROUTINE: writeSingleNETCDFvar
! \label{writeSingleNETCDFvar}
!
! !INTERFACE: 
  subroutine writeSingleNETCDFvar(ftn,ftn_stats,n,dataEntry,&
       non_model_fields)
! !USES: 
    use LIS_coreMod,   only : LIS_rc

    implicit none

    integer,   intent(in)   :: n 
    integer,   intent(in)   :: ftn
    integer,   intent(in)   :: ftn_stats
    type(LIS_metadataEntry), pointer :: dataEntry
    integer,   optional     :: non_model_fields
! 
! !DESCRIPTION: 
!  This routine writes a single variable to a NETCDF file
!  The arguments are: 
!  \begin{description}
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[n] index of the nest
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_writevar\_netcdf](\ref{LIS_writevar_netcdf})
!     writes a variable into a netcdf formatted file. 
!   \end{description}
!EOP    
    integer       :: i,k,m,t
    integer       :: nmodel_status

    nmodel_status = 0
    if(present(non_model_fields)) then 
       nmodel_status = non_model_fields
    endif

    if(dataEntry%selectOpt.eq.1) then

       if(nmodel_status==0) then

          do t=1,LIS_rc%ntiles(n)
             if(nmodel_status==0) then 
                m = 1
             else
                m = LIS_domain(n)%tile(t)%sftype
             endif
             
             do k=1,dataEntry%vlevels
                if(dataEntry%count(t,k).gt.0) then 
                   if(dataEntry%timeAvgOpt.eq.3) then  !do nothing
                      continue   
                   elseif(dataEntry%timeAvgOpt.eq.2.or.&
                        dataEntry%timeAvgOpt.eq.1) then 
                      dataEntry%modelOutput(1,t,k) = &
                           dataEntry%modelOutput(1,t,k)/&
                           dataEntry%count(t,k)
                   else !do nothing
                      continue   
                   endif
                else
                   dataEntry%modelOutput(1,t,k) = LIS_rc%udef
                endif
             enddo
          enddo
       endif

       do k=1,dataEntry%vlevels
          ! accumulated values
          ! time-averaged values and instantaneous values
          if(dataEntry%timeAvgOpt.eq.2) then 
             call LIS_writevar_netcdf(ftn,ftn_stats, n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%varId_def, &
                  trim(dataEntry%short_name)//'('//&
                  trim(dataEntry%units)//')',&
                  dataEntry%form,nmodel_status,dim1=k)
             call LIS_writevar_netcdf(ftn,ftn_stats, n,&
                  dataEntry%modelOutput(2,:,k),&
                  dataEntry%varId_opt1, &
                  trim(dataEntry%short_name)//'('//&
                  trim(dataEntry%units)//')',&
                  dataEntry%form,nmodel_status,dim1=k)
          ! time-averaged values or instantaneous values
          else
             call LIS_writevar_netcdf(ftn,ftn_stats, n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%varId_def, &
                  trim(dataEntry%short_name)//'('//&
                  trim(dataEntry%units)//')',&
                  dataEntry%form,nmodel_status,dim1=k)
          end if ! EMK
          if(dataEntry%minmaxOpt.gt.0) then 
             call LIS_writevar_netcdf(ftn,ftn_stats, n,&
                  dataEntry%minimum(:,k),&
                  dataEntry%varId_min, &
                  trim(dataEntry%short_name)//'_min ('//&
                  trim(dataEntry%units)//')',&
                  dataEntry%form,nmodel_status,dim1=k)
             
             call LIS_writevar_netcdf(ftn,ftn_stats, n,&
                  dataEntry%maximum(:,k),&
                  dataEntry%varId_max, &
                  trim(dataEntry%short_name)//'_max ('//&
                  trim(dataEntry%units)//')',&
                  dataEntry%form,nmodel_status,dim1=k)
             
          endif
       enddo
    endif

  end subroutine writeSingleNETCDFvar

!BOP
! !ROUTINE: writeSingleRoutingNETCDFvar
! \label{writeSingleRoutingNETCDFvar}
!
! !INTERFACE: 
  subroutine writeSingleRoutingNETCDFvar(ftn,ftn_stats,n,dataEntry,&
       non_model_fields)
! !USES: 
    use LIS_coreMod,   only : LIS_rc

    implicit none

    integer,   intent(in)   :: n 
    integer,   intent(in)   :: ftn
    integer,   intent(in)   :: ftn_stats
    type(LIS_metadataEntry), pointer :: dataEntry
    integer,   optional     :: non_model_fields
! 
! !DESCRIPTION: 
!  This routine writes a single variable to a NETCDF file
!  The arguments are: 
!  \begin{description}
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[n] index of the nest
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[writeroutingvar\_netcdf\_real](\ref{writeroutingvar_netcdf_real})
!     writes a variable into a netcdf formatted file. 
!   \end{description}
!EOP    
    integer       :: i,k,t
    integer       :: nmodel_status

    nmodel_status = 0
    if(present(non_model_fields)) then 
       nmodel_status = non_model_fields
    endif

    if(dataEntry%selectOpt.eq.1) then

       if(nmodel_status==0) then

          do t=1,LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)
             do k=1,dataEntry%vlevels
                if(dataEntry%count(t,k).gt.0) then 
                   if(dataEntry%timeAvgOpt.eq.3) then  !do nothing
                      continue   
                   elseif(dataEntry%timeAvgOpt.eq.2.or.&
                        dataEntry%timeAvgOpt.eq.1) then 
                      dataEntry%modelOutput(1,t,k) =&
                           dataEntry%modelOutput(1,t,k)/&
                           dataEntry%count(t,k)
                   else !do nothing
                      continue   
                   endif
                else
                   dataEntry%modelOutput(1,t,k) = LIS_rc%udef
                endif
             enddo
          enddo
       endif

       do k=1,dataEntry%vlevels
          ! accumulated values
          ! time-averaged values and instantaneous values
          if(dataEntry%timeAvgOpt.eq.2) then 
             call writeroutingvar_netcdf_real(ftn,ftn_stats, n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%varId_def, &
                  trim(dataEntry%short_name)//'('//&
                  trim(dataEntry%units)//')',&
                  dataEntry%form,nmodel_status,dim1=k)
             call writeroutingvar_netcdf_real(ftn,ftn_stats, n,&
                  dataEntry%modelOutput(2,:,k),&
                  dataEntry%varId_opt1, &
                  trim(dataEntry%short_name)//'('//&
                  trim(dataEntry%units)//')',&
                  dataEntry%form,nmodel_status,dim1=k)
          ! time-averaged values or instantaneous values
          else
             call writeroutingvar_netcdf_real(ftn,ftn_stats, n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%varId_def, &
                  trim(dataEntry%short_name)//'('//&
                  trim(dataEntry%units)//')',&
                  dataEntry%form,nmodel_status,dim1=k)
          end if ! EMK
          if(dataEntry%minmaxOpt.gt.0) then 
             call writeroutingvar_netcdf_real(ftn,ftn_stats, n,&
                  dataEntry%minimum(:,k),&
                  dataEntry%varId_min, &
                  trim(dataEntry%short_name)//'_min ('//&
                  trim(dataEntry%units)//')',&
                  dataEntry%form,nmodel_status,dim1=k)
             
             call writeroutingvar_netcdf_real(ftn,ftn_stats, n,&
                  dataEntry%maximum(:,k),&
                  dataEntry%varId_max, &
                  trim(dataEntry%short_name)//'_max ('//&
                  trim(dataEntry%units)//')',&
                  dataEntry%form,nmodel_status,dim1=k)
             
          endif
       enddo
    endif

  end subroutine writeSingleRoutingNETCDFvar


!BOP
! !ROUTINE: LIS_writeGlobalHeader_restart
! \label{LIS_writeGlobalHeader_restart}
! 
! !INTERFACE: LIS_writeGlobalHeader_restart
  subroutine LIS_writeGlobalHeader_restart(ftn,n,m, &
       model_name, dimID, dim1,dim2,dim3,dim4,dim5,dim6, dim7, dim8, dim9, dim10,&
       output_format)
! !USES: 

! !ARGUMENTS: 
    integer,   intent(in)     :: n 
    integer,   intent(in)     :: m
    integer,   intent(in)     :: ftn
    character(len=*), intent(in) :: model_name
    integer                   :: dimID(11)
    integer,    optional      :: dim1
    integer,    optional      :: dim2
    integer,    optional      :: dim3
    integer,    optional      :: dim4
    integer,    optional      :: dim5
    integer,    optional      :: dim6
    integer,    optional      :: dim7
    integer,    optional      :: dim8
    integer,    optional      :: dim9
    integer,    optional      :: dim10
    character(len=*), optional :: output_format

! 
! !DESCRIPTION: 
!  This routine writes an output file in the NETCDF format based on the 
!  list of selected output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[model\_name] Name of the model that generates the output
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[defineNETCDFheadervar](\ref{defineNETCDFheaderVar})
!     writes the required headers for a single variable
!   \item[writeSingleNETCDFvar](\ref{writeSingleNETCDFvar})
!     writes a single variable into a netcdf formatted file. 
!   \item[LIS\_verify](\ref{LIS_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP

    integer, dimension(8) :: values
    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    character*20      :: wout

    if(.NOT.PRESENT(output_format)) then 
       wout = LIS_rc%wout
    else
       wout = trim(output_format)
    endif
       
    if(wout.eq."binary") then 
       if(LIS_masterproc) then 
          write(ftn) LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%glbnpatch_red(n,m)
       endif
    elseif(wout.eq."netcdf") then        
#if (defined USE_NETCDF3 || defined USE_NETCDF4)           
       call date_and_time(date,time,zone,values)       
       if(LIS_masterproc) then 
          call LIS_verify(nf90_def_dim(ftn,'ntiles',LIS_rc%glbnpatch_red(n,m),&
               dimID(1)),&
               'nf90_def_dim failed for ntiles in LIS_writeGlobalHeader_restart')
          if(present(dim1)) then
             call LIS_verify(nf90_def_dim(ftn,"dim1",&
                  dim1,dimID(2)),&
                  'nf90_def_dim failed for dim1 in LIS__writeGlobalHeader_restart')
          endif
          if(present(dim2)) then 
             call LIS_verify(nf90_def_dim(ftn,"dim2",&
                  dim2,dimID(3)),&
                  'nf90_def_dim failed for dim2 in LIS__writeGlobalHeader_restart')
          endif
          if(present(dim3)) then
             call LIS_verify(nf90_def_dim(ftn,"dim3",&
                  dim3,dimID(4)),&
                  'nf90_def_dim failed for dim3 in LIS_writeGlobalHeader_restart')
          endif
          if(present(dim4)) then
             call LIS_verify(nf90_def_dim(ftn,"dim4",&
                  dim4,dimID(5)),&
                  'nf90_def_dim failed for dim4 in LIS_writeGlobalHeader_restart')
          endif
          if(present(dim5)) then
             call LIS_verify(nf90_def_dim(ftn,"dim5",&
                  dim5,dimID(6)),&
                  'nf90_def_dim failed for dim5 in LIS_writeGlobalHeader_restart')
          endif
          if(present(dim6)) then
             call LIS_verify(nf90_def_dim(ftn,"dim6",&
                  dim6,dimID(7)),&
                  'nf90_def_dim failed for dim6 in LIS_writeGlobalHeader_restart')
          endif
          if(present(dim7)) then
             call LIS_verify(nf90_def_dim(ftn,"dim7",&
                  dim7,dimID(8)),&
                  'nf90_def_dim failed for dim7 in LIS_writeGlobalHeader_restart')
          endif
          if(present(dim8)) then
             call LIS_verify(nf90_def_dim(ftn,"dim8",&
                  dim8,dimID(9)),&
                  'nf90_def_dim failed for dim8 in LIS_writeGlobalHeader_restart')
          endif
          if(present(dim9)) then
             call LIS_verify(nf90_def_dim(ftn,"dim9",&
                  dim9,dimID(10)),&
                  'nf90_def_dim failed for dim9 in LIS_writeGlobalHeader_restart')
          endif
          if(present(dim10)) then
             call LIS_verify(nf90_def_dim(ftn,"dim10",&
                  dim10,dimID(11)),&
                  'nf90_def_dim failed for dim10 in LIS_writeGlobalHeader_restart')
          endif
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value", &
               LIS_rc%udef),'nf90_put_att failed for missing_value')
          
          
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
               "LIS land surface model restart"),&
               'nf90_put_att failed for title')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
               trim(LIS_rc%institution)),&
               'nf90_put_att failed for institution')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"source",&
               trim(model_name)),&
               'nf90_put_att failed for source')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
               "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
               date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)),&
               'nf90_put_att failed for history')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
               "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"),&
               'nf90_put_att failed for references')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"conventions", &
               "CF-1.6"),'nf90_put_att failed for conventions') !CF version 1.6
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
               "website: http://lis.gsfc.nasa.gov/"),&
               'nf90_put_att failed for comment')

!grid information
          if(LIS_rc%lis_map_proj.eq."latlon") then !latlon
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "EQUIDISTANT CYLINDRICAL"),&
                  'nf90_put_att failed for MAP_PROJECTION')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
                  "SOUTH_WEST_CORNER_LAT", &
                  LIS_rc%gridDesc(n,4)),&
                  'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
                  "SOUTH_WEST_CORNER_LON", &
                  LIS_rc%gridDesc(n,5)),&
                  'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LIS_rc%gridDesc(n,9)),&
                  'nf90_put_att failed for DX')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LIS_rc%gridDesc(n,10)),&
                  'nf90_put_att failed for DY')
          elseif(LIS_rc%lis_map_proj.eq."mercator") then 
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "MERCATOR"),&
                  'nf90_put_att failed for MAP_PROJECTION')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
                  "SOUTH_WEST_CORNER_LAT", &
                  LIS_rc%gridDesc(n,4)),&
                  'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
                  "SOUTH_WEST_CORNER_LON", &
                  LIS_rc%gridDesc(n,5)),&
                  'nf90_put_att failed for SOUTH_WEST_CORNER_LON') 
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LIS_rc%gridDesc(n,10)),&
                  'nf90_put_att failed for TRUELAT1')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LIS_rc%gridDesc(n,11)),&
                  'nf90_put_att failed for STANDARD_LON')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LIS_rc%gridDesc(n,8)),&
                  'nf90_put_att failed for DX')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LIS_rc%gridDesc(n,9)),&
                  'nf90_put_att failed for DY')
          elseif(LIS_rc%lis_map_proj.eq."lambert") then !lambert conformal
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "LAMBERT CONFORMAL"),&
                  'nf90_put_att failed for MAP_PROJECTION')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
                  "SOUTH_WEST_CORNER_LAT", &
                  LIS_rc%gridDesc(n,4)),&
                  'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
                  "SOUTH_WEST_CORNER_LON", &
                  LIS_rc%gridDesc(n,5)),&
                  'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LIS_rc%gridDesc(n,10)),&
                  'nf90_put_att failed for TRUELAT1')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
                  LIS_rc%gridDesc(n,7)),&
                  'nf90_put_att failed for TRUELAT2')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LIS_rc%gridDesc(n,11)),&
                  'nf90_put_att failed for STANDARD_LON')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LIS_rc%gridDesc(n,8)),&
                  'nf90_put_att failed for DX')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LIS_rc%gridDesc(n,9)),&
                  'nf90_put_att failed for DY')
             
          elseif(LIS_rc%lis_map_proj.eq."polar") then ! polar stereographic
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "POLAR STEREOGRAPHIC"),&
                  'nf90_put_att failed for MAP_PROJECTION')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
                  "SOUTH_WEST_CORNER_LAT", &
                  LIS_rc%gridDesc(n,4)),&
                  'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
                  "SOUTH_WEST_CORNER_LON", &
                  LIS_rc%gridDesc(n,5)),&
                  'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LIS_rc%gridDesc(n,10)),&
                  'nf90_put_att failed for TRUELAT1')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
                  LIS_rc%gridDesc(n,7)),&
                  'nf90_put_att failed for ORIENT')
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LIS_rc%gridDesc(n,11)),&
                  'nf90_put_att failed for STANDARD_LON')                  
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LIS_rc%gridDesc(n,8)),&
                  'nf90_put_att failed for DX')                  
             call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LIS_rc%gridDesc(n,9)),&
                  'nf90_put_att failed for DY')                  
          endif
       endif
#endif
    endif
  end subroutine LIS_writeGlobalHeader_restart

!BOP
! !ROUTINE: LIS_writeHeader_restart
! \label{LIS_writeHeader_restart}
! 
! !INTERFACE: 
  subroutine LIS_writeHeader_restart(ftn,n,dimID, vid, standard_name, &
       long_name, units, vlevels, valid_min, valid_max, var_flag)
! !USES: 

! !ARGUMENTS:     
    integer                    :: ftn
    integer                    :: n
    integer                    :: dimID(11)
    integer                    :: vid
    character(len=*)           :: standard_name
    character(len=*)           :: long_name
    character(len=*)           :: units
    integer                    :: vlevels
    real                       :: valid_min
    real                       :: valid_max
    character(len=*), optional :: var_flag
! 
! !DESCRIPTION: 
!    This routine writes the required NETCDF header for a single variable
! 
!   The arguments are: 
!   \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    NETCDF file unit handle
!   \item[dimID]
!    NETCDF dimension ID corresponding to the variable
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \item[LIS\_verify](\ref{LIS_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    integer :: data_index
    integer :: shuffle, deflate, deflate_level
    integer :: dimID_t(2)
    character*50 :: var_flag_tmp
    integer :: fill_value

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    if(present(var_flag)) then 
       var_flag_tmp = var_flag
    else
       var_flag_tmp = ""
    endif

    dimID_t(1) = dimID(1)
    if(LIS_masterproc) then 
       if(var_flag_tmp.eq."dim1") then 
          dimID_t(2) = dimID(2)
       elseif(var_flag_tmp.eq."dim2") then 
          dimID_t(2) = dimID(3)
       elseif(var_flag_tmp.eq."dim3") then 
          dimID_t(2) = dimID(4)
       elseif(var_flag_tmp.eq."dim4") then 
          dimID_t(2) = dimID(5)
       elseif(var_flag_tmp.eq."dim5") then 
          dimID_t(2) = dimID(6)
       elseif(var_flag_tmp.eq."dim6") then 
          dimID_t(2) = dimID(7)
       elseif(var_flag_tmp.eq."dim7") then 
          dimID_t(2) = dimID(8)
       elseif(var_flag_tmp.eq."dim8") then 
          dimID_t(2) = dimID(9)
       elseif(var_flag_tmp.eq."dim9") then 
          dimID_t(2) = dimID(10)
       elseif(var_flag_tmp.eq."dim10") then 
          dimID_t(2) = dimID(11)
       elseif(var_flag_tmp.eq."tbot_lagday") then 
          call LIS_verify(nf90_inq_dimid(ftn, "tbot_lagday", dimID_t(2)),&
               'nf90_inq_dimid for tbot_lagday failed '//&
               'in LIS_writeHeader_restart')
       endif
       if(vlevels.gt.1) then 
          call LIS_verify(nf90_def_var(ftn,trim(standard_name),&
               nf90_float, dimids = dimID_t(1:2), varID=vid),&
               'nf90_def_var(2d) failed in LIS_writeHeader_restart')
#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_fill(ftn,&
               vid, 1,fill_value), 'nf90_def_var_fill failed for '//&
               standard_name)

          call LIS_verify(nf90_def_var_deflate(ftn,&
               vid,shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate(2d) failed in LIS_writeHeader_restart')
#endif
       else
          call LIS_verify(nf90_def_var(ftn,trim(standard_name),&
               nf90_float,dimids = dimID_t(1:1), varID=vid),&
               'nf90_def_var(1d) failed in LIS_writeHeader_restart')

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_fill(ftn,&
               vid, 1,fill_value), 'nf90_def_var_fill failed for '//&
               standard_name)

          call LIS_verify(nf90_def_var_deflate(ftn,&
               vid, shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate(1d) failed in LIS_writeHeader_restart')
#endif
       endif
       call LIS_verify(nf90_put_att(ftn,vid,&
            "units",trim(units)),&
            'nf90_put_att failed for units')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "standard_name",trim(standard_name)))
       call LIS_verify(nf90_put_att(ftn,vid,&
            "long_name",trim(long_name)),&
            'nf90_put_att failed for long_name')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "scale_factor",1.0),&
            'nf90_put_att failed for scale_factor')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "add_offset",0.0),&
            'nf90_put_att failed for add_offset')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "missing_value",LIS_rc%udef),&
            'nf90_put_att failed for missing_value')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "_FillValue",LIS_rc%udef),&
            'nf90_put_att failed for _FillValue')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "vmin",valid_min),&
            'nf90_put_att failed for vmin')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "vmax",valid_max),&
            'nf90_put_att failed for vmax')
    endif
#endif
  end subroutine LIS_writeHeader_restart

!BOP
! !ROUTINE: LIS_closeHeader_restart
! \label{LIS_closeHeader_restart}
! 
! !INTERFACE: 
  subroutine LIS_closeHeader_restart(ftn,n,m,dimID, rstInterval)

    implicit none

! !ARGUMENTS:     
    integer            :: ftn
    integer            :: n
    integer            :: m
    real               :: rstInterval
    integer            :: dimID(4)
! 
! !DESCRIPTION: 
!    This routine closes the required NETCDF header.
! 
!   The arguments are: 
!   \begin{description}
!   \item[ftn]
!    NETCDF file unit handle
!   \item[n]
!    index of the nest
!   \item[m]
!    index of the surface type
!   \item[dimID]
!    NETCDF dimension ID corresponding to the variable
!   \item[rstInterval]
!    the restart interval
!   \end{description}

!EOP
    integer                 :: tdimID,xtimeID
    integer                 :: t,c,r,index1
    type(LIS_metadataEntry) :: xlat,xlong,xtime
    integer                 :: sindex, eindex
    character*8             :: xtime_begin_date
    character*6             :: xtime_begin_time
    character*50            :: xtime_units
    character*50            :: xtime_timeInc
    integer                 :: ftn_stats
    integer                 :: iret
    real,     allocatable   :: gvar(:)
    integer                 :: shuffle, deflate, deflate_level
! Note that the fix to add lat/lon to the NETCDF output will output
! undefined values for the water points. 

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
       
    xlat%short_name = "lat"
    xlat%long_name = "latitude"
    xlat%standard_name = "latitude"
    xlat%units = "degree_north"
    xlat%nunits = 1
    xlat%format = 'F'
    xlat%form = 1
    xlat%vlevels = 1
    xlat%timeAvgOpt = 0 
    xlat%selectOpt = 1
    xlat%minMaxOpt = 0 
    xlat%stdOpt = 0 
    allocate(xlat%modelOutput(1,LIS_rc%npatch(n,m),xlat%vlevels))
    allocate(xlat%count(1,xlat%vlevels))
    xlat%count = 1
    allocate(xlat%unittypes(1))
    xlat%unittypes(1) = "degree_north"
    
    xlong%short_name = "lon"
    xlong%long_name = "longitude"
    xlong%standard_name = "longitude"
    xlong%units = "degree_east"
    xlong%nunits = 1
    xlong%format = 'F'
    xlong%form = 1
    xlong%vlevels = 1
    xlong%timeAvgOpt = 0 
    xlong%selectOpt = 1
    xlong%minMaxOpt = 0 
    xlong%stdOpt = 0 
    allocate(xlong%modelOutput(1,LIS_rc%npatch(n,m),xlong%vlevels))
    allocate(xlong%count(1,xlong%vlevels))
    xlong%count = 1
    allocate(xlong%unittypes(1))
    xlong%unittypes(1) = "degree_east"

!lat header
    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    if(LIS_masterproc) then 
       call LIS_verify(nf90_def_var(ftn,trim(xlat%short_name),&
            nf90_float,&
            dimids = dimID(1:1), varID=xlat%varId_def),&
            'nf90_def_var failed for xlat:short_name')
#if(defined USE_NETCDF4)                
       call LIS_verify(nf90_def_var_deflate(ftn,&
            xlat%varId_def, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for xlat')
#endif
       call LIS_verify(nf90_put_att(ftn,xlat%varId_def,&
            "units",trim(xlat%units)),&
            'nf90_put_att failed for xlat:units')
       call LIS_verify(nf90_put_att(ftn,xlat%varId_def,&
            "standard_name",trim(xlat%standard_name)),&
            'nf90_put_att failed for xlat:standard_name')
       call LIS_verify(nf90_put_att(ftn,xlat%varId_def,&
            "long_name",trim(xlat%long_name)),&
            'nf90_put_att failed for xlat:long_name')
       call LIS_verify(nf90_put_att(ftn,xlat%varId_def,&
            "scale_factor",1.0),&
            'nf90_put_att failed for xlat:scale_factor')
       call LIS_verify(nf90_put_att(ftn,xlat%varId_def,&
            "add_offset",0.0),&
            'nf90_put_att failed for xlat:add_offset')
       call LIS_verify(nf90_put_att(ftn,xlat%varId_def,&
            "missing_value",LIS_rc%udef),&
            'nf90_put_att failed for xlat:missing_value')
       call LIS_verify(nf90_put_att(ftn,xlat%varId_def,&
            "_FillValue",LIS_rc%udef),&
            'nf90_put_att failed for xlat:_FillValue')
       call LIS_verify(nf90_put_att(ftn,xlat%varId_def,&
            "vmin",0.0),&
            'nf90_put_att failed for xlat:vmin')
       call LIS_verify(nf90_put_att(ftn,xlat%varId_def,&
            "vmax",0.0),&
            'nf90_put_att failed for xlat:vmax')
!lon
       call LIS_verify(nf90_def_var(ftn,trim(xlong%short_name),&
            nf90_float,&
            dimids = dimID(1:1), varID=xlong%varId_def),&
            'nf90_def_var failed for xlong:short_name')
#if(defined USE_NETCDF4)                
       call LIS_verify(nf90_def_var_deflate(ftn,&
            xlong%varId_def, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for xlong')
#endif
       call LIS_verify(nf90_put_att(ftn,xlong%varId_def,&
            "units",trim(xlong%units)),&
            'nf90_put_att failed for xlong:units')
       call LIS_verify(nf90_put_att(ftn,xlong%varId_def,&
            "standard_name",trim(xlong%standard_name)),&
            'nf90_put_att failed for xlong:standard_name')
       call LIS_verify(nf90_put_att(ftn,xlong%varId_def,&
            "long_name",trim(xlong%long_name)),&
            'nf90_put_att failed for xlong:long_name')
       call LIS_verify(nf90_put_att(ftn,xlong%varId_def,&
            "scale_factor",1.0),&
            'nf90_put_att failed for xlong:scale_factor')
       call LIS_verify(nf90_put_att(ftn,xlong%varId_def,&
            "add_offset",0.0),&
            'nf90_put_att failed for xlong:add_offset')
       call LIS_verify(nf90_put_att(ftn,xlong%varId_def,&
            "missing_value",LIS_rc%udef),&
            'nf90_put_att failed for xlong:missing_value')
       call LIS_verify(nf90_put_att(ftn,xlong%varId_def,&
            "_FillValue",LIS_rc%udef),&
            'nf90_put_att failed for xlong:_FillValue')
       call LIS_verify(nf90_put_att(ftn,xlong%varId_def,&
            "vmin",0.0),&
            'nf90_put_att failed for xlong:vmin')
       call LIS_verify(nf90_put_att(ftn,xlong%varId_def,&
            "vmax",0.0),&
            'nf90_put_att failed for xlong:vmax')
       !defining time field
       call LIS_verify(nf90_def_dim(ftn,'time',NF90_UNLIMITED,tdimID),&
            'nf90_def_dim failed for time')
       
       call LIS_verify(nf90_def_var(ftn,'time',&
            nf90_float,dimids = tdimID, varID=xtimeID),&
            'nf90_def_var failed for time')
    
       write(xtime_units,200) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
            LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
200 format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
         I2.2,':',I2.2)
       write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
            LIS_rc%yr, LIS_rc%mo, LIS_rc%da
       write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
            LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
       write(xtime_timeInc, fmt='(I20)') &
            nint(rstInterval)
       
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "units",trim(xtime_units)),&
            'nf90_put_att failed for time:units')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "long_name","time"),&
            'nf90_put_att failed for time:long_name')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "time_increment",trim(adjustl(xtime_timeInc))),&
            'nf90_put_att failed for time:time_increment')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "begin_date",xtime_begin_date),&
            'nf90_put_att failed for time:begin_date')
       call LIS_verify(nf90_put_att(ftn,xtimeID,&
            "begin_time",xtime_begin_time),&
            'nf90_put_att failed for time:begin_time')
       
       call LIS_verify(nf90_enddef(ftn),&
            'nf90_enddef failed')       
    endif
    do t=1,LIS_rc%npatch(n,m)
       c = LIS_domain(n)%tile(t)%col
       r = LIS_domain(n)%tile(t)%row
       index1 = LIS_domain(n)%gindex(c,r)
       xlat%modelOutput(1,t,1) = LIS_domain(n)%grid(index1)%lat
       xlong%modelOutput(1,t,1) = LIS_domain(n)%grid(index1)%lon
    enddo

    call LIS_gather_patch_vector_output(n,m,gvar,xlat%modelOutput(1,:,1))
    if(LIS_masterproc) then 
       iret = nf90_put_var(ftn,xlat%varId_def,gvar,(/1,1/),&
            (/LIS_rc%glbnpatch_red(n,m),1/))
       call LIS_verify(iret,'nf90_put_var failed for xlat ')
       deallocate(gvar)
    endif

    call LIS_gather_patch_vector_output(n,m,gvar,xlong%modelOutput(1,:,1))
    if(LIS_masterproc) then 
       iret = nf90_put_var(ftn,xlong%varId_def,gvar,(/1,1/),&
            (/LIS_rc%glbnpatch_red(n,m),1/))
       call LIS_verify(iret,'nf90_put_var failed for xlong ')
       deallocate(gvar)
    endif
    if(LIS_masterproc) then 
       call LIS_verify(nf90_put_var(ftn,xtimeID,0.0),&
            'nf90_put_var failed for time')
    endif

    deallocate(xlat%modelOutput)
    deallocate(xlat%count)
    deallocate(xlat%unittypes)
    
    deallocate(xlong%modelOutput)
    deallocate(xlong%count)
    deallocate(xlong%unittypes)    
#endif
  end subroutine LIS_closeHeader_restart

!BOP
! !ROUTINE: writevar_bin_withstats_real
! \label{writevar_bin_withstats_real}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_bin
  subroutine writevar_bin_withstats_real(ftn,ftn_stats,n,var,mvar,form)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: ftn
    integer, intent(in) :: ftn_stats
    real, intent(in)    :: var(LIS_rc%ntiles(n))
    character (len=*)   :: mvar
    integer, intent(in) :: form
! !DESCRIPTION:
!  Write a real variable to a binary output file with a number of diagnostic 
!  statistics (mean, standardar deviation, min/max during the temporal output
!  interval) written to a text file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [ftn\_stats]
!     unit number of the ASCII text statistics file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item[mvar]
!     short name of the variable being written (will be used in the stats file)
!   \item[form]
!     format to be used in the stats file (1-decimal format, 
!     2-scientific format)
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_gather\_tiled\_vector\_output](\ref{LIS_gather_tiled_vector_output})
!     call to gather the 1d tiled output variables in tile space
!   \item[LIS\_gather\_gridded\_output](\ref{LIS_gather_gridded_output})
!     call to gather the 1d tiled output variable into a 2d gridded array
!   \item[LIS\_gather\_gridded\_vector\_output](\ref{LIS_gather_gridded_vector_output})
!     call to gather the 1d tiled output variable into a 1d gridded array
!   \item[write\_stats](\ref{write_stats})
!     call to compute the diagnostic statistics
!  \end{description}
!
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)

    if(LIS_rc%wopt.eq."1d tilespace") then !tiled output
       call LIS_gather_tiled_vector_output(n,gtmp1,var)
       
       if(LIS_masterproc) then 
          write(ftn) gtmp1
          if(ftn_stats.ne.-1) then
             call write_stats(gtmp1, LIS_rc%glbntiles_red(n), mvar, ftn_stats,form)
          endif
          deallocate(gtmp1)
       endif
       
    elseif(LIS_rc%wopt.eq."2d gridspace") then !2d gridded output
       call LIS_gather_gridded_output(n, gtmp, var)
       
       if ( LIS_masterproc ) then
          write(ftn) gtmp
          if(ftn_stats.ne.-1) then
             call write_stats(gtmp, LIS_rc%gnc(n)*LIS_rc%gnr(n), &
                  mvar, ftn_stats, form)
          endif
          deallocate(gtmp)
       endif
    elseif(LIS_rc%wopt.eq."1d gridspace") then !1d gridded output
       call LIS_gather_gridded_vector_output(n, gtmp1, var)
       
       if ( LIS_masterproc ) then
          write(ftn) gtmp1
          if(ftn_stats.ne.-1) then
             call write_stats(gtmp1, LIS_rc%glbngrid_red(n), &
                  mvar, ftn_stats, form)
          endif
          deallocate(gtmp1)
       endif
    endif
  end subroutine writevar_bin_withstats_real

!BOP
! !ROUTINE: writevar_dist_bin
! \label{writevar_dist_bin}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_bin
  subroutine writevar_dist_bin(ftn,n,var,form)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: ftn
    real, intent(in)    :: var(LIS_rc%ntiles(n))
    integer, intent(in) :: form
! !DESCRIPTION:
!  Write a real variable to a binary output file with a number of diagnostic 
!  statistics (mean, standardar deviation, min/max during the temporal output
!  interval) written to a text file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [ftn\_stats]
!     unit number of the ASCII text statistics file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item[mvar]
!     short name of the variable being written (will be used in the stats file)
!   \item[form]
!     format to be used in the stats file (1-decimal format, 
!     2-scientific format)
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_gather\_tiled\_vector\_output](\ref{LIS_gather_tiled_vector_output})
!     call to gather the 1d tiled output variables in tile space
!   \item[LIS\_gather\_gridded\_output](\ref{LIS_gather_gridded_output})
!     call to gather the 1d tiled output variable into a 2d gridded array
!   \item[LIS\_gather\_gridded\_vector\_output](\ref{LIS_gather_gridded_vector_output})
!     call to gather the 1d tiled output variable into a 1d gridded array
!   \item[write\_stats](\ref{write_stats})
!     call to compute the diagnostic statistics
!  \end{description}
!
!EOP
    REAL          :: gtmp(LIS_rc%lnc(n),LIS_rc%lnr(n))
    integer       :: gid, c,r,i,t,m
    
    if(LIS_rc%wopt.eq."1d tilespace") then !tiled output
       write(LIS_logunit,*) '[ERR] 1d tilespace option is not supported'
       write(LIS_logunit,*) '[ERR] in the distributed binary mode'
       call LIS_endrun()
       
    elseif(LIS_rc%wopt.eq."2d gridspace") then !2d gridded output

       gtmp = LIS_rc%udef
       
       do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)

          gid = LIS_domain(n)%tile(i)%index

          c = LIS_domain(n)%grid(gid)%col
          r = LIS_domain(n)%grid(gid)%row

          gtmp(c,r) = 0.0
          
          do m=1,LIS_rc%nensem(n)
             t = i+m-1
             if ( var(t) == -9999.0 ) then
                gtmp(c,r) = -9999.0
             else
                gtmp(c,r) = gtmp(c,r)+&
                     var(t)*LIS_domain(n)%tile(t)%fgrd*&
                     LIS_domain(n)%tile(t)%pens
             endif
          enddo                    
       enddo       
       
       write(ftn) gtmp

    elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then !2d gridded output

       do m=1, LIS_rc%nensem(n)

          gtmp = -9999.0
          
          do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
             
             gid = LIS_domain(n)%tile(i)%index
             
             c = LIS_domain(n)%grid(gid)%col
             r = LIS_domain(n)%grid(gid)%row
             
             t = i+m-1
             
             if ( var(t) == -9999.0 ) then
                gtmp(c,r) = -9999.0
             else
                gtmp(c,r) = var(t)
             endif
          enddo
          write(ftn) gtmp
       enddo             
       
    elseif(LIS_rc%wopt.eq."1d gridspace") then !1d gridded output

       write(LIS_logunit,*) '[ERR] 1d gridspace option is not supported'
       write(LIS_logunit,*) '[ERR] in the distributed binary mode'
       call LIS_endrun()
       
    endif
  end subroutine writevar_dist_bin

!BOP
! !ROUTINE: writeroutingvar_dist_bin
! \label{writeroutingvar_dist_bin}
! 
! !INTERFACE:
! Private name: call using LIS_writeroutingvar_bin
  subroutine writeroutingvar_dist_bin(ftn,n,var,form)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: ftn
    real, intent(in)    :: var(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n))
    integer, intent(in) :: form
! !DESCRIPTION:
!  Write a real variable to a binary output file with a number of diagnostic 
!  statistics (mean, standardar deviation, min/max during the temporal output
!  interval) written to a text file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [ftn\_stats]
!     unit number of the ASCII text statistics file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item[mvar]
!     short name of the variable being written (will be used in the stats file)
!   \item[form]
!     format to be used in the stats file (1-decimal format, 
!     2-scientific format)
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_gather\_tiled\_vector\_output](\ref{LIS_gather_tiled_vector_output})
!     call to gather the 1d tiled output variables in tile space
!   \item[LIS\_gather\_gridded\_output](\ref{LIS_gather_gridded_output})
!     call to gather the 1d tiled output variable into a 2d gridded array
!   \item[LIS\_gather\_gridded\_vector\_output](\ref{LIS_gather_gridded_vector_output})
!     call to gather the 1d tiled output variable into a 1d gridded array
!   \item[write\_stats](\ref{write_stats})
!     call to compute the diagnostic statistics
!  \end{description}
!
!EOP
    REAL          :: gtmp(LIS_rc%lnc(n),LIS_rc%lnr(n))
    integer       :: gid, c,r,i,t,m
    
    if(LIS_rc%wopt.eq."1d tilespace") then !tiled output
       write(LIS_logunit,*) '[ERR] 1d tilespace option is not supported'
       write(LIS_logunit,*) '[ERR] in the distributed binary mode'
       call LIS_endrun()
       
    elseif(LIS_rc%wopt.eq."2d gridspace") then !2d gridded output

       gtmp = LIS_rc%udef
       
       do i=1,LIS_rc%nroutinggrid(n),LIS_rc%nensem(n)

          gid = LIS_routing(n)%tile(i)%index

          c = LIS_routing(n)%grid(gid)%col
          r = LIS_routing(n)%grid(gid)%row

          gtmp(c,r) = 0.0
          
          do m=1,LIS_rc%nensem(n)
             t = i+m-1
             if ( var(t) == -9999.0 ) then
                gtmp(c,r) = -9999.0
             else
                gtmp(c,r) = gtmp(c,r)+&
                     var(t)*LIS_routing(n)%tile(t)%fgrd*&
                     LIS_routing(n)%tile(t)%pens
             endif
          enddo                    
       enddo       
       
       write(ftn) gtmp

    elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then !2d gridded output

       do m=1, LIS_rc%nensem(n)

          gtmp = -9999.0
          
          do i=1,LIS_rc%nroutinggrid(n),LIS_rc%nensem(n)
             
             gid = LIS_routing(n)%tile(i)%index
             
             c = LIS_routing(n)%grid(gid)%col
             r = LIS_routing(n)%grid(gid)%row
             
             t = i+m-1
             
             if ( var(t) == -9999.0 ) then
                gtmp(c,r) = -9999.0
             else
                gtmp(c,r) = var(t)
             endif
          enddo
          write(ftn) gtmp
       enddo             
       
    elseif(LIS_rc%wopt.eq."1d gridspace") then !1d gridded output

       write(LIS_logunit,*) '[ERR] 1d gridspace option is not supported'
       write(LIS_logunit,*) '[ERR] in the distributed binary mode'
       call LIS_endrun()
       
    endif
  end subroutine writeroutingvar_dist_bin
  
!BOP
! !ROUTINE: writevar_bin_real
! \label{writevar_bin_real}
!
! !INTERFACE:
! Private name: call using LIS_writevar_bin
  subroutine writevar_bin_real(ftn, n, var)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: ftn
    real, intent(in)    :: var(LIS_rc%ntiles(n))
! !DESCRIPTION:
!  Write a real variable to a binary output file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_gather\_tiled\_vector\_output](\ref{LIS_gather_tiled_vector_output})
!     call to gather the 1d tiled output variables in tile space
!   \item[LIS\_gather\_gridded\_output](\ref{LIS_gather_gridded_output})
!     call to gather the 1d tiled output variable into a 2d gridded array
!   \item[LIS\_gather\_gridded\_vector\_output](\ref{LIS_gather_gridded_vector_output})
!     call to gather the 1d tiled output variable into a 1d gridded array
!  \end{description}
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)

    if(LIS_rc%wopt.eq."1d tilespace") then 

       call LIS_gather_tiled_vector_output(n,gtmp1,var)

       if ( LIS_masterproc ) then
          write(ftn) gtmp1
          deallocate(gtmp1)
       endif

    elseif(LIS_rc%wopt.eq."2d gridspace") then 
       call LIS_gather_gridded_output(n, gtmp, var)

       if ( LIS_masterproc ) then
          write(ftn) gtmp
          deallocate(gtmp)
       endif
    elseif(LIS_rc%wopt.eq."1d gridspace") then 
       call LIS_gather_gridded_vector_output(n, gtmp1, var)

       if ( LIS_masterproc ) then
          write(ftn) gtmp1
          deallocate(gtmp1)
       endif
    endif
  end subroutine writevar_bin_real

!BOP
! !ROUTINE: writevar_bin_real_direct
! \label{writevar_bin_real_direct}
!
! !INTERFACE:
! Private name: call using LIS_writevar_bin
  subroutine writevar_bin_real_direct(ftn, n, var, direct)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: ftn
    real,    intent(in) :: var(LIS_rc%lnc(n), LIS_rc%lnr(n))
    integer, intent(in) :: direct 
! !DESCRIPTION:
!  Write a real variable to a binary output file in a direct access format
!
!  The arguments are: 
!  \begin{description}
!   \item [ftn]
!     unit number of the binary output file
!   \item [n]
!     index of the domain or nest.
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item [direct]
!     dummy argument indicating direct access output 
!  \end{description}
!EOP
    real, allocatable :: var1(:) 
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)
    integer :: count1 ,c,r,gid,ntiles,ierr,l
    integer :: gdeltas

    allocate(var1(LIS_rc%ngrid(n)))
    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       allocate(gtmp1(LIS_rc%glbngrid(n)))
    else
       allocate(gtmp1(1))
    endif
    do r=1,LIS_rc%lnr(n)
       do c=1,LIS_rc%lnc(n)
          if(LIS_domain(n)%gindex(c,r).ne.-1) then 
             var1(LIS_domain(n)%gindex(c,r)) = var(c,r)
          endif
       enddo
    enddo
#if (defined SPMD)      
    gdeltas = LIS_gdeltas(n,LIS_localPet)
    call MPI_GATHERV(var1,gdeltas,&
         MPI_REAL,gtmp1,LIS_gdeltas(n,:),LIS_goffsets(n,:),&
         MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
    gtmp1 = var1
#endif
    if(LIS_masterproc) then 
       gtmp = LIS_rc%udef
       count1=1
       do l=1,LIS_npes
          do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                gid = c+(r-1)*LIS_rc%gnc(n)
                ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                if(ntiles.ne.0) then 
                   if(r.ge.LIS_nss_ind(n,l).and.&
                        r.le.LIS_nse_ind(n,l).and.&
                        c.ge.LIS_ews_ind(n,l).and.&
                        c.le.LIS_ewe_ind(n,l))then !points not in halo                      
                      gtmp(c,r) = gtmp1(count1)
                   endif
                   count1 = count1 + 1
                endif
             enddo
          enddo
       enddo
       
       write(ftn,rec=direct) gtmp

!       do r=1,LIS_rc%gnr(n)
!          do c=1,LIS_rc%gnc(n)
!             line = c+(r-1)*LIS_rc%gnc(n)
!             write(ftn,rec=line) gtmp(c,r)
!          enddo
!       enddo
       deallocate(gtmp)
    endif
    deallocate(gtmp1)
    deallocate(var1)
  end subroutine writevar_bin_real_direct

!BOP
! !ROUTINE: writevar_restart_tile_int
! \label{writevar_restart_tile_int}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_restart
  subroutine writevar_restart_tile_int(ftn, n, var, varid, dim, wformat)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    integer             :: var(LIS_rc%ntiles(n))
    integer,          optional :: varid
    integer,          optional :: dim
    character(len=*), optional :: wformat
! !DESCRIPTION:
!  Writes a int variable to a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    integer, allocatable :: gtmp(:)
    integer, allocatable :: gtmp1(:)
    integer :: count1 ,c,r,ntiles,t,gid,stid,tid,l
    integer :: ierr
    character*20 :: wform
    integer      :: tdeltas 

    if(present(wformat)) then 
       wform = wformat
    else
       wform = "binary"
    endif

    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%glbntiles_red(n)))
       allocate(gtmp1(LIS_rc%glbntiles(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)      
    tdeltas = LIS_tdeltas(n,LIS_localPet)
    call MPI_GATHERV(var,tdeltas,&
         MPI_INTEGER,gtmp1,LIS_tdeltas(n,:),LIS_toffsets(n,:),&
         MPI_INTEGER,0,LIS_mpi_comm,ierr)
#else 
    gtmp1 = var
#endif
    if(LIS_masterproc) then
       count1=1
       do l=1,LIS_npes
          do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                gid = c+(r-1)*LIS_rc%gnc(n)
                ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                stid = LIS_domain(n)%str_tind(gid)
                if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                   do t=1,ntiles
                      tid = stid + t-1
                      gtmp(tid) = gtmp1(count1)
                      count1 = count1 + 1
                   enddo
                else
                   count1 = count1 + ntiles
                endif
             enddo
          enddo
       enddo
       if(wform.eq."binary") then 
          write(ftn) gtmp
       elseif(wform.eq."netcdf") then 
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
          if(PRESENT(dim)) then 
             ierr = nf90_put_var(ftn,varid,gtmp,(/1,dim/),&
                  (/LIS_rc%glbntiles_red(n),1/))
             call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')

          else            
             ierr = nf90_put_var(ftn,varid,gtmp,(/1/),&
                  (/LIS_rc%glbntiles_red(n)/))
             call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')
          endif
#endif          
       endif
       deallocate(gtmp)
    endif
    deallocate(gtmp1)

  end subroutine writevar_restart_tile_int

!BOP
! !ROUTINE: writevar_restart_tile_real
! \label{writevar_restart_tile_real}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_restart
  subroutine writevar_restart_tile_real(ftn, n, var, varid, dim, wformat)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    real                :: var(LIS_rc%ntiles(n))
    integer,          optional :: varid
    integer,          optional :: dim
    character(len=*), optional :: wformat
! !DESCRIPTION:
!  Writes a real variable to a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp1(:)
    integer :: count1 ,c,r,ntiles,t,gid,stid,tid,l
    integer :: ierr
    character*20 :: wform
    integer      :: tdeltas 

    if(present(wformat)) then 
       wform = wformat
    else
       wform = "binary"
    endif

    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%glbntiles_red(n)))
       allocate(gtmp1(LIS_rc%glbntiles(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)      
    tdeltas = LIS_tdeltas(n,LIS_localPet)
    call MPI_GATHERV(var,tdeltas,&
         MPI_REAL,gtmp1,LIS_tdeltas(n,:),LIS_toffsets(n,:),&
         MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
    gtmp1 = var
#endif
    if(LIS_masterproc) then
       count1=1
       do l=1,LIS_npes
          do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                gid = c+(r-1)*LIS_rc%gnc(n)
                ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                stid = LIS_domain(n)%str_tind(gid)
                if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                   do t=1,ntiles
                      tid = stid + t-1
                      gtmp(tid) = gtmp1(count1)
                      count1 = count1 + 1
                   enddo
                else
                   count1 = count1 + ntiles
                endif
             enddo
          enddo
       enddo
       if(wform.eq."binary") then 
          write(ftn) gtmp
       elseif(wform.eq."netcdf") then 
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
          if(PRESENT(dim)) then 
             ierr = nf90_put_var(ftn,varid,gtmp,(/1,dim/),&
                  (/LIS_rc%glbntiles_red(n),1/))
             call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')

          else            
             ierr = nf90_put_var(ftn,varid,gtmp,(/1/),&
                  (/LIS_rc%glbntiles_red(n)/))
             call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')
          endif
#endif          
       endif
       deallocate(gtmp)
    endif
    deallocate(gtmp1)

  end subroutine writevar_restart_tile_real

!BOP
! !ROUTINE: writevar_restart_patch_int
! \label{writevar_restart_patch_int}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_restart
  subroutine writevar_restart_patch_int(ftn, n, m, var, varid, dim, wformat)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer             :: var(LIS_rc%npatch(n,m))
    integer,          optional :: varid
    integer,          optional :: dim
    character(len=*), optional :: wformat
! !DESCRIPTION:
!  Writes a real variable to a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    integer, allocatable :: gtmp(:)
    integer, allocatable :: gtmp1(:)
    integer :: count1 ,c,r,npatch,t,gid,stid,tid,l
    integer :: ierr
    integer :: patch_deltas
    character*20 :: wform

    if(present(wformat)) then
       wform = wformat
    else
       wform = "binary"
    endif

    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%glbnpatch_red(n,m)))
       allocate(gtmp1(LIS_rc%glbnpatch(n,m)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)      
    patch_deltas = LIS_patch_deltas(n,m,LIS_localPet)
    call MPI_GATHERV(var,patch_deltas,&
         MPI_INTEGER,gtmp1,LIS_patch_deltas(n,m,:),LIS_patch_offsets(n,m,:),&
         MPI_INTEGER,0,LIS_mpi_comm,ierr)
#else 
    gtmp1 = var
#endif

    if(LIS_masterproc) then
       count1=1
       do l=1,LIS_npes
          do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                gid = c+(r-1)*LIS_rc%gnc(n)
                npatch = LIS_surface(n,m)%npatch_pergrid(gid)
                stid = LIS_surface(n,m)%str_patch_ind(gid)
                if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                   do t=1,npatch
                      tid = stid + t-1
                      gtmp(tid) = gtmp1(count1)
                      count1 = count1 + 1
                   enddo
                else
                   count1 = count1 + npatch
                endif
             enddo
          enddo
       enddo
       if(wform.eq."binary") then 
          write(ftn) gtmp
       elseif(wform.eq."netcdf") then 
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
          if(PRESENT(dim)) then 
             ierr = nf90_put_var(ftn,varid,gtmp,(/1,dim/),&
                  (/LIS_rc%glbnpatch_red(n,m),1/))
             call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')

          else            
             ierr = nf90_put_var(ftn,varid,gtmp,(/1/),&
                  (/LIS_rc%glbnpatch_red(n,m)/))
             call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')
          endif
#endif          
       endif
       deallocate(gtmp)
    endif
    deallocate(gtmp1)
  end subroutine writevar_restart_patch_int

!BOP
! !ROUTINE: writevar_restart_patch_real
! \label{writevar_restart_patch_real}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_restart
  subroutine writevar_restart_patch_real(ftn, n, m, var, varid, dim, wformat)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    integer, intent(in) :: m
    real                :: var(LIS_rc%npatch(n,m))
    integer,          optional :: varid
    integer,          optional :: dim
    character(len=*), optional :: wformat
! !DESCRIPTION:
!  Writes a real variable to a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp1(:)
    integer :: count1 ,c,r,npatch,t,gid,stid,tid,l
    integer :: ierr
    integer      :: patch_deltas
    character*20 :: wform

    if(present(wformat)) then
       wform = wformat
    else
       wform = "binary"
    endif

    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%glbnpatch_red(n,m)))
       allocate(gtmp1(LIS_rc%glbnpatch(n,m)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)      
    patch_deltas = LIS_patch_deltas(n,m,LIS_localPet)
    call MPI_GATHERV(var,patch_deltas,&
         MPI_REAL,gtmp1,LIS_patch_deltas(n,m,:),LIS_patch_offsets(n,m,:),&
         MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
    gtmp1 = var
#endif

    if(LIS_masterproc) then
       count1=1
       do l=1,LIS_npes
          do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                gid = c+(r-1)*LIS_rc%gnc(n)
                npatch = LIS_surface(n,m)%npatch_pergrid(gid)
                stid = LIS_surface(n,m)%str_patch_ind(gid)
                if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                   do t=1,npatch
                      tid = stid + t-1
                      gtmp(tid) = gtmp1(count1)
                      count1 = count1 + 1
                   enddo
                else
                   count1 = count1 + npatch
                endif
             enddo
          enddo
       enddo
       if(wform.eq."binary") then 
          write(ftn) gtmp
       elseif(wform.eq."netcdf") then 
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
          if(PRESENT(dim)) then 
             ierr = nf90_put_var(ftn,varid,gtmp,(/1,dim/),&
                  (/LIS_rc%glbnpatch_red(n,m),1/))
             call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')

          else            
             ierr = nf90_put_var(ftn,varid,gtmp,(/1/),&
                  (/LIS_rc%glbnpatch_red(n,m)/))
             call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')
          endif
#endif          
       endif
       deallocate(gtmp)
    endif
    deallocate(gtmp1)
  end subroutine writevar_restart_patch_real

!BOP
! !ROUTINE: writevar_tile_noensemble_real
! \label{writevar_tile_noensemble_real}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_reduced_tilespace
  subroutine writevar_tile_noensemble_real(ftn, n, var)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    real, intent(in)    :: var(LIS_rc%ntiles(n)/LIS_rc%nensem(n))
! !DESCRIPTION:
!  Writes a real variable to a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp1(:)
    integer :: count1 ,c,r,ntiles,t,gid,stid,tid,l
    integer :: tdeltas
    integer :: ierr
    
    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%glbntiles_red(n)/LIS_rc%nensem(n)))
       allocate(gtmp1(LIS_rc%glbntiles(n)/LIS_rc%nensem(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)
    tdeltas = LIS_tdeltas(n,LIS_localPet)/LIS_rc%nensem(n)
    call MPI_GATHERV(var,tdeltas,&
         MPI_REAL,gtmp1,LIS_tdeltas(n,:)/LIS_rc%nensem(n),&
         LIS_toffsets(n,:)/LIS_rc%nensem(n),MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
    gtmp1 = var
#endif
    if(LIS_masterproc) then 
       count1=1
       do l=1,LIS_npes
          do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                gid = c+(r-1)*LIS_rc%gnc(n)
                ntiles = LIS_domain(n)%ntiles_pergrid(gid)/LIS_rc%nensem(n)
                stid = (LIS_domain(n)%str_tind(gid)-1)/LIS_rc%nensem(n)+1
                if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                   do t=1,ntiles
                      tid = stid + t-1
                      gtmp(tid) = gtmp1(count1)                      
                      count1 = count1 + 1
                   enddo
                else
                   count1 = count1 + ntiles
                endif
             enddo
          enddo
       enddo
       
       write(ftn) gtmp
       deallocate(gtmp)
    endif
    deallocate(gtmp1)
  end subroutine writevar_tile_noensemble_real


!BOP
! !ROUTINE: readvar_tile_noensemble_real
! \label{readvar_tile_noensemble_real}
! 
! !INTERFACE:
! Private name: call LIS_readvar_reduced_tilespace
  subroutine readvar_tile_noensemble_real(ftn, n, var)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: n
    real, intent(inout)   :: var(LIS_rc%ntiles(n)/LIS_rc%nensem(n))
! !DESCRIPTION:
!  Reads a real variable from a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:)
    integer :: count1 ,c,r,ntiles,t,gid,stid,tid

    allocate(gtmp(LIS_rc%glbntiles_red(n)/LIS_rc%nensem(n)))
    read(ftn) gtmp

    count1=1
    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          gid = c+(r-1)*LIS_rc%gnc(n)
          ntiles = LIS_domain(n)%ntiles_pergrid(gid)/LIS_rc%nensem(n)
          stid = (LIS_domain(n)%str_tind(gid)-1)/LIS_rc%nensem(n)+1
          do t=1,ntiles
             tid = stid + t-1
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
    deallocate(gtmp)   
  end subroutine readvar_tile_noensemble_real

#if 0 
!BOP
! !ROUTINE: readvar_restart_tile_int
! \label{readvar_restart_tile_int}
!
! !INTERFACE:
! Private name: call using LIS_readvar_restart
  subroutine readvar_restart_tile_int(ftn, n, var)
! !USES:
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in)    :: ftn
    integer, intent(in)    :: n
    integer, intent(inout) :: var(LIS_rc%ntiles(n))
! !DESCRIPTION:
!  Reads an integer variable from a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    integer, allocatable :: gtmp(:)
    integer :: count1 ,c,r,ntiles,t,gid,stid,tid

    allocate(gtmp(LIS_rc%glbntiles_red(n)))
    read(ftn) gtmp
    
    count1=1
    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          gid = c+(r-1)*LIS_rc%gnc(n)
          ntiles = LIS_domain(n)%ntiles_pergrid(gid)
          stid = LIS_domain(n)%str_tind(gid)
          do t=1,ntiles
             tid = stid + t-1
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
    deallocate(gtmp)

  end subroutine readvar_restart_tile_int
#endif
!BOP
! !ROUTINE: readvar_restart_patch_int
! \label{readvar_restart_patch_int}
! 
! !INTERFACE:
! Private name: call LIS_readvar_restart
  subroutine readvar_restart_patch_int(ftn, n, m,var, varname,dim,vlevels, wformat)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)    :: ftn
    integer, intent(in)    :: n
    integer, intent(in)    :: m
    integer, intent(inout) :: var(LIS_rc%npatch(n,m))
    character(len=*), optional :: varname
    integer,          optional :: dim
    integer,          optional :: vlevels
    character(len=*), optional :: wformat
    
! !DESCRIPTION:
!  Reads a real variable from a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    integer, allocatable :: gtmp(:)
    integer, allocatable :: gtmp_v(:,:)
    integer :: varid
    integer :: status
    integer :: count1 ,c,r,npatch,t,gid,stid,tid
    character*20 :: wform

    if(present(wformat)) then 
       wform = wformat
    else
       wform = "binary"
    endif
    allocate(gtmp(LIS_rc%glbnpatch_red(n,m)))
    if(wform.eq."binary") then 
       read(ftn) gtmp
    elseif(wform.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
       status = nf90_inq_varid(ftn,trim(varname),varid)
       call LIS_verify(status,'Error in nf90_inq_varid in LIS_readvar_restart')

       if(present(dim).and.present(vlevels)) then 
          allocate(gtmp_v(LIS_rc%glbnpatch_red(n,m),vlevels))
          status = nf90_get_var(ftn,varid,gtmp_v)
          call LIS_verify(status,'Error in nf90_get_var in LIS_readvar_restart')
          gtmp = gtmp_v(:,dim)
          deallocate(gtmp_v)
       else
          status = nf90_get_var(ftn,varid,gtmp)
          call LIS_verify(status,'Error in nf90_get_var in LIS_readvar_restart')
       endif
#endif       
    endif
    count1=1
    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          gid = c+(r-1)*LIS_rc%gnc(n)
          npatch = LIS_surface(n,m)%npatch_pergrid(gid)
          stid = LIS_surface(n,m)%str_patch_ind(gid)
          do t=1,npatch
             tid = stid + t-1
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
    deallocate(gtmp)   
  end subroutine readvar_restart_patch_int


!BOP
! !ROUTINE: readvar_restart_patch_real
! \label{readvar_restart_patch_real}
! 
! !INTERFACE:
! Private name: call LIS_readvar_restart
  subroutine readvar_restart_patch_real(ftn, n, m,var, varname,dim,vlevels, wformat)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: n
    integer, intent(in)   :: m
    real, intent(inout)   :: var(LIS_rc%npatch(n,m))
    character(len=*), optional :: varname
    integer,          optional :: dim
    integer,          optional :: vlevels
    character(len=*), optional :: wformat
    
! !DESCRIPTION:
!  Reads a real variable from a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp_v(:,:)
    integer :: varid
    integer :: status
    integer :: count1 ,c,r,npatch,t,gid,stid,tid
    integer :: glbnpatch_size
    character*20 :: wform

    if(present(wformat)) then 
       wform = wformat
    else
       wform = "binary"
    endif
    glbnpatch_size = LIS_rc%glbnpatch_red(n,m)
    allocate(gtmp(glbnpatch_size))
    if(wform.eq."binary") then 
       read(ftn) gtmp
    elseif(wform.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
       status = nf90_inq_varid(ftn,trim(varname),varid)
       call LIS_verify(status,'Error in nf90_inq_varid in LIS_readvar_restart')

       if(present(dim).and.present(vlevels)) then 
          if ( dim > vlevels ) then
             write(LIS_logunit,*) '[ERR] LIS_readvar_restart: ' // &
                'requested level greater than total number of levels'
             call LIS_endrun
          endif
          allocate(gtmp_v(glbnpatch_size,1))
          status = nf90_get_var(ftn,varid,gtmp_v,&
             start=(/1,dim/),count=(/glbnpatch_size,1/))
          call LIS_verify(status,'Error in nf90_get_var in LIS_readvar_restart')
          gtmp = gtmp_v(:,1)
          deallocate(gtmp_v)
       else
          status = nf90_get_var(ftn,varid,gtmp)
          call LIS_verify(status,'Error in nf90_get_var in LIS_readvar_restart')
       endif
#endif       
    endif
    count1=1
    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          gid = c+(r-1)*LIS_rc%gnc(n)
          npatch = LIS_surface(n,m)%npatch_pergrid(gid)
          stid = LIS_surface(n,m)%str_patch_ind(gid)
          do t=1,npatch
             tid = stid + t-1
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
    deallocate(gtmp)   
  end subroutine readvar_restart_patch_real

!BOP
! !ROUTINE: readvar_restart_tile_real
! \label{readvar_restart_tile_real}
! 
! !INTERFACE:
! Private name: call LIS_readvar_restart
  subroutine readvar_restart_tile_real(ftn, n, var, varname,dim,vlevels, wformat)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: n
    real, intent(inout)   :: var(LIS_rc%ntiles(n))
    character(len=*), optional :: varname
    integer,          optional :: dim
    integer,          optional :: vlevels
    character(len=*), optional :: wformat
    
! !DESCRIPTION:
!  Reads a real variable from a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp_v(:,:)
    integer :: varid
    integer :: status
    integer :: count1 ,c,r,ntiles,t,gid,stid,tid
    character*20 :: wform

    if(present(wformat)) then 
       wform = wformat
    else
       wform = "binary"
    endif
    allocate(gtmp(LIS_rc%glbntiles_red(n)))
    if(wform.eq."binary") then 
       read(ftn) gtmp
    elseif(wform.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
       status = nf90_inq_varid(ftn,trim(varname),varid)
       call LIS_verify(status,'Error in nf90_inq_varid in LIS_readvar_restart')

       if(present(dim).and.present(vlevels)) then 
          allocate(gtmp_v(LIS_rc%glbntiles_red(n),vlevels))
          status = nf90_get_var(ftn,varid,gtmp_v)
          call LIS_verify(status,'Error in nf90_get_var in LIS_readvar_restart')
          gtmp = gtmp_v(:,dim)
          deallocate(gtmp_v)
       else
          status = nf90_get_var(ftn,varid,gtmp)
          call LIS_verify(status,'Error in nf90_get_var in LIS_readvar_restart')
       endif
#endif       
    endif
    count1=1
    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          gid = c+(r-1)*LIS_rc%gnc(n)
          ntiles = LIS_domain(n)%ntiles_pergrid(gid)
          stid = LIS_domain(n)%str_tind(gid)
          do t=1,ntiles
             tid = stid + t-1
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
    deallocate(gtmp)   
  end subroutine readvar_restart_tile_real

!BOP
! !ROUTINE: readvar_restart_tile_int
! \label{readvar_restart_tile_int}
! 
! !INTERFACE:
! Private name: call LIS_readvar_restart
  subroutine readvar_restart_tile_int(ftn, n, var, varname,dim,vlevels, wformat)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: n
    integer, intent(inout)   :: var(LIS_rc%ntiles(n))
    character(len=*), optional :: varname
    integer,          optional :: dim
    integer,          optional :: vlevels
    character(len=*), optional :: wformat
    
! !DESCRIPTION:
!  Reads a integer variable from a binary restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    integer, allocatable :: gtmp(:)
    integer, allocatable :: gtmp_v(:,:)
    integer :: varid
    integer :: status
    integer :: count1 ,c,r,ntiles,t,gid,stid,tid
    character*20 :: wform

    if(present(wformat)) then 
       wform = wformat
    else
       wform = "binary"
    endif
    allocate(gtmp(LIS_rc%glbntiles_red(n)))
    if(wform.eq."binary") then 
       read(ftn) gtmp
    elseif(wform.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
       status = nf90_inq_varid(ftn,trim(varname),varid)
       call LIS_verify(status,'Error in nf90_inq_varid in LIS_readvar_restart')

       if(present(dim).and.present(vlevels)) then 
          allocate(gtmp_v(LIS_rc%glbntiles_red(n),vlevels))
          status = nf90_get_var(ftn,varid,gtmp_v)
          call LIS_verify(status,'Error in nf90_get_var in LIS_readvar_restart')
          gtmp = gtmp_v(:,dim)
          deallocate(gtmp_v)
       else
          status = nf90_get_var(ftn,varid,gtmp)
          call LIS_verify(status,'Error in nf90_get_var in LIS_readvar_restart')
       endif
#endif       
    endif
    count1=1
    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          gid = c+(r-1)*LIS_rc%gnc(n)
          ntiles = LIS_domain(n)%ntiles_pergrid(gid)
          stid = LIS_domain(n)%str_tind(gid)
          do t=1,ntiles
             tid = stid + t-1
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
    deallocate(gtmp)   
  end subroutine readvar_restart_tile_int


!BOP
! !ROUTINE: readvar_1dgridded_real
! \label{readvar_1dgridded_real}
! 
! !INTERFACE:
! Private name: call LIS_readvar_gridded
  subroutine readvar_1dgridded_real(ftn, n, var)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: n
    real, intent(inout)   :: var(LIS_rc%ngrid(n))
! !DESCRIPTION:
!  Reads a real variable from a binary 
!  sequential access, gridded file. After reading
!  the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:,:)
    integer :: c,r,gid
    integer :: nc,c1,r1

    allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
    read(ftn) gtmp
    
    nc = (LIS_ewe_halo_ind(n,LIS_localPet+1)-LIS_ews_halo_ind(n,LIS_localPet+1))+1

    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          c1 = c-LIS_ews_halo_ind(n,LIS_localPet+1)+1
          r1 = r-LIS_nss_halo_ind(n,LIS_localPet+1)+1
!          gid = r1+(c1-1)*nc
          gid = LIS_domain(n)%gindex(c1,r1)
          if(gid.ne.-1) then
             var(gid) = gtmp(c,r)
          endif
       enddo
    enddo
    deallocate(gtmp)

  end subroutine readvar_1dgridded_real

!BOP
! !ROUTINE: readvar_1dgridded_fromvector_real
! \label{readvar_1dgridded_fromvector_real}
! 
! !INTERFACE:
! Private name: call LIS_readvar_gridded
  subroutine readvar_1dgridded_fromvector_real(ftn, n, var, oned)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: n
    real, intent(inout)   :: var(LIS_rc%ngrid(n))
    integer, intent(in)   :: oned
! !DESCRIPTION:
!  Reads a real variable from a binary 
!  sequential access, 1d gridded file. After reading
!  the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item [oned]
!     dummy variable to distinguish the interface. 
!  \end{description}
!EOP
    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp2d(:,:)
    integer :: i, gid
    integer :: c,r, nc, cnt
    integer :: c1,c2,r1,r2

    allocate(gtmp(LIS_rc%glbngrid(n)))
    read(ftn) gtmp

    allocate(gtmp2d(LIS_rc%gnc(n), LIS_rc%gnr(n)))

    cnt = 1
    do r=1,LIS_rc%gnr(n)
       do c=1,LIS_rc%gnc(n)
          gid = c+(r-1)*LIS_rc%gnc(n)
          if(LIS_domain(n)%ntiles_pergrid(gid).gt.0) then 
             gtmp2d(c,r) = gtmp(cnt)
             cnt = cnt+1
          endif
       enddo
    enddo

    nc = (LIS_ewe_halo_ind(n,LIS_localPet+1)-LIS_ews_halo_ind(n,LIS_localPet+1))+1

    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          c1 = c-LIS_ews_halo_ind(n,LIS_localPet+1)+1
          r1 = r-LIS_nss_halo_ind(n,LIS_localPet+1)+1
!          gid = r1+(c1-1)*nc
          gid = LIS_domain(n)%gindex(c1,r1)
          if(gid.ne.-1) then
             var(gid) = gtmp2d(c,r)
          endif
       enddo
    enddo

    deallocate(gtmp)
    deallocate(gtmp2d)


  end subroutine readvar_1dgridded_fromvector_real

!BOP
! !ROUTINE: readvar_2dgridded_real
! \label{readvar_2dgridded_real}
! 
! !INTERFACE:
! Private name: call LIS_readvar_gridded
  subroutine readvar_2dgridded_real(ftn, n, var)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: n
    real, intent(inout)   :: var(LIS_rc%lnc(n), LIS_rc%lnr(n))
! !DESCRIPTION:
!  Reads a real variable from a binary 
!  sequential access, gridded file. After reading
!  the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:,:)
    integer :: c,r
    integer :: nc,c1,r1

    allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
    read(ftn) gtmp
    
    nc = (LIS_ewe_halo_ind(n,LIS_localPet+1)-LIS_ews_halo_ind(n,LIS_localPet+1))+1

    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          c1 = c-LIS_ews_halo_ind(n,LIS_localPet+1)+1
          r1 = r-LIS_nss_halo_ind(n,LIS_localPet+1)+1
          var(c1,r1) = gtmp(c,r)
       enddo
    enddo
    deallocate(gtmp)

  end subroutine readvar_2dgridded_real

!BOP
! !ROUTINE: writevar_gridded_real
! \label{writevar_gridded_real}
! 
! !INTERFACE:
! Private name: call LIS_writevar_gridded
  subroutine writevar_gridded_real(ftn, n, var, wopt)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: n
    real, intent(inout)   :: var(LIS_rc%ngrid(n))
    character(len=*), intent(in), optional :: wopt
! !DESCRIPTION:
!  Writes a real variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)
    real, allocatable :: gtmp2(:)
    integer       :: l
    integer       :: c,r,gid,ntiles, count1, ierr
    character*20  :: wopt_temp
    integer       :: gdeltas

    if(present(wopt)) then 
       wopt_temp = wopt
    else
       wopt_temp = LIS_rc%wopt
    endif

    if(wopt_temp.eq."2d gridspace") then 
       if(LIS_masterproc) then 
          allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(gtmp1(LIS_rc%glbngrid(n)))
          gtmp = 0.0
          gtmp1 = 0.0
       else
          allocate(gtmp1(1))
          gtmp1 = 0.0
       endif
#if (defined SPMD)
       gdeltas = LIS_gdeltas(n,LIS_localPet)
       call MPI_GATHERV(var,gdeltas,MPI_REAL,gtmp1,&
            LIS_gdeltas(n,:),LIS_goffsets(n,:),MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
       gtmp1 = var
#endif
       if(LIS_masterproc) then 
          gtmp = LIS_rc%udef
          count1=1
          do l=1,LIS_npes
             do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                   gid = c+(r-1)*LIS_rc%gnc(n)
                   ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                   if(ntiles.ne.0) then       
                      if(r.ge.LIS_nss_ind(n,l).and.&
                           r.le.LIS_nse_ind(n,l).and.&
                           c.ge.LIS_ews_ind(n,l).and.&
                           c.le.LIS_ewe_ind(n,l)) then !points not in halo                   
                         gtmp(c,r) = gtmp1(count1)
                      endif
                      count1 = count1 + 1
                   endif
                enddo
             enddo
          enddo
          write(ftn) gtmp
          deallocate(gtmp)
       endif
       deallocate(gtmp1)
    elseif(wopt_temp.eq."1d gridspace") then 
       if(LIS_masterproc) then 
          allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(gtmp1(LIS_rc%glbngrid(n)))
          allocate(gtmp2(LIS_rc%glbngrid_red(n)))
          gtmp = 0.0
          gtmp1 = 0.0
       else
          allocate(gtmp1(1))
          gtmp1 = 0.0
       endif
#if (defined SPMD)
       gdeltas = LIS_gdeltas(n,LIS_localPet)
       call MPI_GATHERV(var,gdeltas,MPI_REAL,gtmp1,&
            LIS_gdeltas(n,:),LIS_goffsets(n,:),MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
       gtmp1 = var
#endif
       if(LIS_masterproc) then 
          gtmp = LIS_rc%udef
          count1=1
          do l=1,LIS_npes
             do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                   gid = c+(r-1)*LIS_rc%gnc(n)
                   ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                   if(ntiles.ne.0) then       
                      if(r.ge.LIS_nss_ind(n,l).and.&
                           r.le.LIS_nse_ind(n,l).and.&
                           c.ge.LIS_ews_ind(n,l).and.&
                           c.le.LIS_ewe_ind(n,l)) then !points not in halo                   
                         gtmp(c,r) = gtmp1(count1)
                      endif
                      count1 = count1 + 1
                   endif
                enddo
             enddo
          enddo

          count1 = 1
          do r=1,LIS_rc%gnr(n)
             do c=1,LIS_rc%gnc(n)
                gid = c+(r-1)*LIS_rc%gnc(n)
                ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                if(ntiles.ne.0) then    
                   gtmp2(count1) = gtmp(c,r)
                   count1 = count1 + 1
                endif
             enddo
          enddo

          write(ftn) gtmp2
          deallocate(gtmp)
          deallocate(gtmp2)
       endif
       deallocate(gtmp1)
    endif

  end subroutine writevar_gridded_real

!BOP
! !ROUTINE: writevar_gridded_real_withstats
! \label{writevar_gridded_real_withstats}
! 
! !INTERFACE:
! Private name: call LIS_writevar_gridded
  subroutine writevar_gridded_real_withstats(ftn,ftn_stats,n,var,mvar,form)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: ftn_stats
    integer, intent(in)   :: n
    real, intent(inout)   :: var(LIS_rc%ngrid(n))
    character(len=*)      :: mvar
    integer, intent(in)   :: form
! !DESCRIPTION:
!  Writes a real variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)
    real, allocatable :: gtmp2(:)
    integer       :: l
    integer       :: gdeltas
    integer       :: c,r,gid,ntiles, count1, ierr

    if(LIS_rc%wopt.eq."2d gridspace") then 
       if(LIS_masterproc) then 
          allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(gtmp1(LIS_rc%glbngrid(n)))
          gtmp = 0.0
          gtmp1 = 0.0
       else
          allocate(gtmp1(1))
          gtmp1 = 0.0
       endif
#if (defined SPMD)
       gdeltas = LIS_gdeltas(n,LIS_localPet)
       call MPI_GATHERV(var,gdeltas,MPI_REAL,gtmp1,&
            LIS_gdeltas(n,:),LIS_goffsets(n,:),MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
       gtmp1 = var
#endif
       if(LIS_masterproc) then 
          gtmp = LIS_rc%udef
          count1=1
          do l=1,LIS_npes
             do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                   gid = c+(r-1)*LIS_rc%gnc(n)
                   ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                   if(ntiles.ne.0) then       
                      if(r.ge.LIS_nss_ind(n,l).and.&
                           r.le.LIS_nse_ind(n,l).and.&
                           c.ge.LIS_ews_ind(n,l).and.&
                           c.le.LIS_ewe_ind(n,l)) then !points not in halo                   
                         gtmp(c,r) = gtmp1(count1)
                      endif
                      count1 = count1 + 1
                   endif
                enddo
             enddo
          enddo
          write(ftn) gtmp
          call write_stats(gtmp, LIS_rc%gnc(n)*LIS_rc%gnr(n), mvar, ftn_stats, form)
          deallocate(gtmp)
       endif
       deallocate(gtmp1)
    elseif(LIS_rc%wopt.eq."1d gridspace") then 
       if(LIS_masterproc) then 
          allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(gtmp1(LIS_rc%glbngrid(n)))
          allocate(gtmp2(LIS_rc%glbngrid_red(n)))
          gtmp = 0.0
          gtmp1 = 0.0
       else
          allocate(gtmp1(1))
          gtmp1 = 0.0
       endif
#if (defined SPMD)
       gdeltas = LIS_gdeltas(n,LIS_localPet)
       call MPI_GATHERV(var,LIS_gdeltas(n,LIS_localPet),MPI_REAL,gtmp1,&
            LIS_gdeltas(n,:),LIS_goffsets(n,:),MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
       gtmp1 = var
#endif
       if(LIS_masterproc) then 
          gtmp = LIS_rc%udef
          count1=1
          do l=1,LIS_npes
             do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                   gid = c+(r-1)*LIS_rc%gnc(n)
                   ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                   if(ntiles.ne.0) then       
                      if(r.ge.LIS_nss_ind(n,l).and.&
                           r.le.LIS_nse_ind(n,l).and.&
                           c.ge.LIS_ews_ind(n,l).and.&
                           c.le.LIS_ewe_ind(n,l)) then !points not in halo                   
                         gtmp(c,r) = gtmp1(count1)
                      endif
                      count1 = count1 + 1
                   endif
                enddo
             enddo
          enddo

          count1 = 1
          do r=1,LIS_rc%gnr(n)
             do c=1,LIS_rc%gnc(n)
                gid = c+(r-1)*LIS_rc%gnc(n)
                ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                if(ntiles.ne.0) then    
                   gtmp2(count1) = gtmp(c,r)
                   count1 = count1 + 1
                endif
             enddo
          enddo

          write(ftn) gtmp2
          call write_stats(gtmp2, LIS_rc%glbngrid_red(n), mvar, ftn_stats, form)
          deallocate(gtmp)
          deallocate(gtmp2)
       endif
       deallocate(gtmp1)
    endif

  end subroutine writevar_gridded_real_withstats

!BOP
! !ROUTINE: writevar_netcdf_withstats_real
! \label{writevar_netcdf_withstats_real}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_netcdf
  subroutine writevar_netcdf_withstats_real(ftn,ftn_stats, n, var,varid, mvar,&
       form, nmodel_status,dim1)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: ftn
    integer, intent(in) :: ftn_stats
    integer             :: varid
    real, intent(in)    :: var(LIS_rc%ntiles(n))
    character (len=*)   :: mvar
    integer, intent(in) :: form
    integer, intent(in) :: nmodel_status
    integer, intent(in), optional :: dim1
    integer             :: gindex
!
! !DESCRIPTION:
!  Write a real variable to a netcdf output file with some diagnostic 
!  statistics written to a text file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the netcdf output file
!   \item [ftn\_stats]
!     unit number of the ASCII text statistics file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item[mvar]
!     name of the variable being written (will be used in the stats file)
!   \item[form]
!     format to be used in the stats file (1-decimal format, 
!     2-scientific format)
!   \item [flag]
!    option to determine if the variable needs to be written (1-write, 
!    0-do not write)
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[stats](\ref{stats}) \newline
!     call to compute the diagnostic statistics
!  \end{description}
!
!EOP
    integer             :: l, iret
    real :: vmean,vstdev,vmin,vmax

    real, allocatable :: var1(:)
    real, allocatable :: var1_ens(:,:)
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmplat(:),gtmplon(:)
    real, allocatable :: gtmp_ens(:,:,:)
    real, allocatable :: gtmp1(:)
    real, allocatable :: gtmp1_ens(:,:)
    integer :: gdeltas
    integer :: count1 ,c,r,m,gid,ntiles,ierr,i,t

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    ! Write output in 1-D ensemble tile array space:
    if(LIS_rc%wopt.eq."1d tilespace") then !tiled output
       call LIS_gather_tiled_vector_output(n,gtmp1,var)
       
       if(LIS_masterproc) then 
          
          if(PRESENT(dim1)) then 
             iret = nf90_put_var(ftn,varid,gtmp1,(/1,dim1/),&
                  (/LIS_rc%glbntiles_red(n),1/))
             call LIS_verify(iret,'nf90_put_var failed in LIS_historyMod')
          else            
             iret = nf90_put_var(ftn,varid,gtmp1,(/1/),&
                  (/LIS_rc%glbntiles_red(n)/))
             call LIS_verify(iret,'nf90_put_var failed in LIS_historyMod')
          endif
          if(ftn_stats.ne.-1) then
             if ( LIS_rc%sout ) then
                call stats(gtmp1,LIS_rc%udef,LIS_rc%glbntiles_red(n), &
                           vmean,vstdev,vmin,vmax)
                if(form==1) then
                   write(ftn_stats,999) mvar,vmean,vstdev,vmin,vmax
                elseif(form==2) then
                   write(ftn_stats,998) mvar,vmean,vstdev,vmin,vmax
                endif
                flush(ftn_stats)
             endif
          endif
          deallocate(gtmp1)
       endif

    ! Write output in 2d grid space:
    elseif(LIS_rc%wopt.eq."2d gridspace") then
       
       allocate(var1(LIS_rc%ngrid(n)))
       if(LIS_masterproc) then 
          allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(gtmp1(LIS_rc%glbngrid(n)))
          gtmp = 0.0
          gtmp1 = 0.0
       else
          allocate(gtmp1(1))
          gtmp1 = 0.0
       endif
       var1 = 0 
       do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
          c = LIS_domain(n)%tile(i)%index
          do m=1,LIS_rc%nensem(n)
             t = i+m-1
             if ( var(t) == -9999.0 ) then
                var1(c) = -9999.0
             else
                var1(c) = var1(c) + &
                     var(t)*LIS_domain(n)%tile(t)%fgrd*&
                     LIS_domain(n)%tile(t)%pens
             endif
          enddo
       enddo

#if (defined SPMD)      
       gdeltas = LIS_gdeltas(n,LIS_localPet)
       call MPI_GATHERV(var1,gdeltas,&
            MPI_REAL,gtmp1,LIS_gdeltas(n,:),LIS_goffsets(n,:),&
            MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
       gtmp1 = var1
#endif
       deallocate(var1)
       if(LIS_masterproc) then 
          gtmp = LIS_rc%udef
          count1=1
          do l=1,LIS_npes
             do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                   gid = c+(r-1)*LIS_rc%gnc(n)
                   ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                   if(ntiles.ne.0) then                 
                      if(r.ge.LIS_nss_ind(n,l).and.&
                           r.le.LIS_nse_ind(n,l).and.&
                           c.ge.LIS_ews_ind(n,l).and.&
                           c.le.LIS_ewe_ind(n,l))then !points not in halo
                         gtmp(c,r) = gtmp1(count1)
                      endif
                      count1 = count1 + 1
                   endif
                enddo
             enddo
          enddo

          ! The latlon fields are written to 1D
          if(LIS_rc%nlatlon_dimensions == '1D') then
             if(nmodel_status.eq.1) then   ! lat
                allocate(gtmplat(LIS_rc%gnr(n)))
                gtmplat = LIS_rc%udef
             
                do r=1,LIS_rc%gnr(n)
                   do c=1,LIS_rc%gnc(n)
                     gindex = c+(r-1)*LIS_rc%gnc(n) 
                     gtmplat(r) = LIS_domain(n)%glat(gindex)
                   enddo
                enddo
             
                iret = nf90_put_var(ftn,varid,gtmplat,(/1/),&
                     (/LIS_rc%gnr(n)/))
                deallocate(gtmplat) 

             elseif(nmodel_status.eq.2) then !lon
                allocate(gtmplon(LIS_rc%gnc(n)))
                gtmplon = LIS_rc%udef

                do r=1,LIS_rc%gnr(n)
                   do c=1,LIS_rc%gnc(n)
                      gindex = c+(r-1)*LIS_rc%gnc(n)
                      gtmplon(c) = LIS_domain(n)%glon(gindex)
                   enddo
                enddo             
             
                iret = nf90_put_var(ftn,varid,gtmplon,(/1/),&
                     (/LIS_rc%gnc(n)/))
                deallocate(gtmplon)

             else
                if(PRESENT(dim1)) then 
                   iret = nf90_put_var(ftn,varid,gtmp,(/1,1,dim1/),&
                        (/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
                else            
                   iret = nf90_put_var(ftn,varid,gtmp,(/1,1/),&
                        (/LIS_rc%gnc(n),LIS_rc%gnr(n)/))
                endif
             endif

          ! The latlon fields are written to 2D
          else 
             if(PRESENT(dim1)) then
                iret = nf90_put_var(ftn,varid,gtmp,(/1,1,dim1/),&
                   (/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
             else
                iret = nf90_put_var(ftn,varid,gtmp,(/1,1/),&
                   (/LIS_rc%gnc(n),LIS_rc%gnr(n)/))
             endif
          endif


          if(ftn_stats.ne.-1) then
             if ( LIS_rc%sout ) then
                call stats(gtmp,LIS_rc%udef,LIS_rc%gnc(n)*LIS_rc%gnr(n),&
                           vmean,vstdev,vmin,vmax)
                if(form==1) then
                   write(ftn_stats,999) mvar,vmean,vstdev,vmin,vmax
                elseif(form==2) then
                   write(ftn_stats,998) mvar,vmean,vstdev,vmin,vmax
                endif
                flush(ftn_stats)
             endif
          endif
          deallocate(gtmp)
       endif
       deallocate(gtmp1)

    ! Write output in 2D ensemble grid space:
    elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 

       ! Non-model output field status (T=non-model; F=model-based):
       if(nmodel_status.ne.0) then   ! non-model output field status
          allocate(var1(LIS_rc%ngrid(n)))
          if(LIS_masterproc) then 
             allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
             allocate(gtmp1(LIS_rc%glbngrid(n)))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1))
             gtmp1 = 0.0
          endif
          var1 = 0 
          do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
             c = LIS_domain(n)%tile(i)%index
             do m=1,LIS_rc%nensem(n)
                t = i+m-1
                if ( var(t) == -9999.0 ) then
                   var1(c) = -9999.0
                else
                   var1(c) = var1(c) + &
                        var(t)*LIS_domain(n)%tile(t)%fgrd*&
                        LIS_domain(n)%tile(t)%pens
                endif
             enddo
          enddo
          
#if (defined SPMD)      
          gdeltas = LIS_gdeltas(n,LIS_localPet)
          call MPI_GATHERV(var1,gdeltas,&
               MPI_REAL,gtmp1,LIS_gdeltas(n,:),LIS_goffsets(n,:),&
               MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
          gtmp1 = var1
#endif
          deallocate(var1)
          if(LIS_masterproc) then 
             gtmp = LIS_rc%udef
             count1=1
             do l=1,LIS_npes
                do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                   do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                      gid = c+(r-1)*LIS_rc%gnc(n)
                      ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                      if(ntiles.ne.0) then                 
                         if(r.ge.LIS_nss_ind(n,l).and.&
                              r.le.LIS_nse_ind(n,l).and.&
                              c.ge.LIS_ews_ind(n,l).and.&
                              c.le.LIS_ewe_ind(n,l))then !points not in halo
                            gtmp(c,r) = gtmp1(count1)
                         endif
                         count1 = count1 + 1
                      endif
                   enddo
                enddo
             enddo

             ! The lat/lon fields are written to 1D
             if(LIS_rc%nlatlon_dimensions == '1D') then             
                if(nmodel_status.eq.1) then   ! lat
                   allocate(gtmplat(LIS_rc%gnr(n)))
                   gtmplat = LIS_rc%udef
                
                   do r=1,LIS_rc%gnr(n)
                      do c=1,LIS_rc%gnc(n)
                         gindex = c+(r-1)*LIS_rc%gnc(n)
                         gtmplat(r) = LIS_domain(n)%glat(gindex)
                      enddo
                   enddo
                
                   iret = nf90_put_var(ftn,varid,gtmplat,(/1/),&
                        (/LIS_rc%gnr(n)/))
                   deallocate(gtmplat)
                
                elseif(nmodel_status.eq.2) then !lon
                   allocate(gtmplon(LIS_rc%gnc(n)))
                   gtmplon = LIS_rc%udef
                
                   do r=1,LIS_rc%gnr(n)
                      do c=1,LIS_rc%gnc(n)
                         gindex = c+(r-1)*LIS_rc%gnc(n)
                         gtmplon(c) = LIS_domain(n)%glon(gindex)
                      enddo
                   enddo
                
                   iret = nf90_put_var(ftn,varid,gtmplon,(/1/),&
                        (/LIS_rc%gnc(n)/))
                   deallocate(gtmplon) 

                else             
                   if(PRESENT(dim1)) then 
                      iret = nf90_put_var(ftn,varid,gtmp,(/1,1,dim1/),&
                           (/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
                   else            
                      iret = nf90_put_var(ftn,varid,gtmp,(/1,1/),&
                           (/LIS_rc%gnc(n),LIS_rc%gnr(n)/))
                   endif
                endif

             ! The latlon fields are written to 2D
             else 
                if(PRESENT(dim1)) then
                   iret = nf90_put_var(ftn,varid,gtmp,(/1,1,dim1/),&
                        (/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
                else
                   iret = nf90_put_var(ftn,varid,gtmp,(/1,1/),&
                        (/LIS_rc%gnc(n),LIS_rc%gnr(n)/))
                endif
             endif
             
             if(ftn_stats.ne.-1) then
                if ( LIS_rc%sout ) then
                   call stats(gtmp,LIS_rc%udef,LIS_rc%gnc(n)*LIS_rc%gnr(n),&
                              vmean,vstdev,vmin,vmax)
                   if(form==1) then
                      write(ftn_stats,999) mvar,vmean,vstdev,vmin,vmax
                   elseif(form==2) then
                      write(ftn_stats,998) mvar,vmean,vstdev,vmin,vmax
                   endif
                   flush(ftn_stats)
                endif
             endif
             deallocate(gtmp)
          endif
          deallocate(gtmp1)
          
       ! Model-based field output:
       else

          allocate(var1_ens(LIS_rc%ngrid(n), LIS_rc%nensem(n)))
          allocate(var1(LIS_rc%ngrid(n))) ! EMK 
          if(LIS_masterproc) then 
             allocate(gtmp_ens(LIS_rc%gnc(n),LIS_rc%gnr(n), LIS_rc%nensem(n)))
             allocate(gtmp1_ens(LIS_rc%glbngrid(n),LIS_rc%nensem(n)))
             gtmp_ens = 0.0
             gtmp1_ens = 0.0
          else
             allocate(gtmp1_ens(1,LIS_rc%nensem(n)))
             gtmp1_ens = 0.0
          endif

          var1_ens = 0 
          do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
             c = LIS_domain(n)%tile(i)%index
             do m=1,LIS_rc%nensem(n)
                t = i+m-1
                if ( var(t) == -9999.0 ) then
                   var1_ens(c,m) = -9999.0
                else
                   var1_ens(c,m) = &
                        var(t)*LIS_domain(n)%tile(t)%fgrd
                endif
             enddo
          enddo
       
#if (defined SPMD)      
          gdeltas = LIS_gdeltas(n,LIS_localPet)
          do m=1,LIS_rc%nensem(n)
             ! EMK: It is possible that the first dimension of var1_ens is 0 
             ! (no grid points with tiles in the PET).  Unfortunately, slicing
             ! such an array [e.g., var1_ens(:,m)] will cause an array bounds 
             ! error.  So, we add some defensive code here to (a) copy a 
             ! slice to a 1-d array only if the dimension is > 0; and (b) 
             ! always pass the 1d array to MPI_GATHERV.  Note that no memory 
             ! access error will occur in MPI_GATHERV for the zero-grid count 
             ! case as long as gdeltas is also zero.
             if (LIS_rc%ngrid(n) > 0) then
                var1(:) = var1_ens(:,m)
             end if
             call MPI_GATHERV(var1,gdeltas,&
                  MPI_REAL,gtmp1_ens(:,m),LIS_gdeltas(n,:),LIS_goffsets(n,:),&
                  MPI_REAL,0,LIS_mpi_comm,ierr)
          enddo
#else 
          do m=1,LIS_rc%nensem(n)
             gtmp1_ens(:,m) = var1_ens(:,m)
          enddo
#endif

          deallocate(var1)     ! EMK...Avoid memory leak
          deallocate(var1_ens) ! EMK...Avoid memory leak

          if(LIS_masterproc) then
             gtmp_ens = LIS_rc%udef
             do m=1,LIS_rc%nensem(n)
                count1=1
                do l=1,LIS_npes
                   do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                      do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                         gid = c+(r-1)*LIS_rc%gnc(n)
                         ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                         if(ntiles.ne.0) then                 
                            if(r.ge.LIS_nss_ind(n,l).and.&
                                 r.le.LIS_nse_ind(n,l).and.&
                                 c.ge.LIS_ews_ind(n,l).and.&
                                 c.le.LIS_ewe_ind(n,l))then !points not in halo
                               gtmp_ens(c,r,m) = gtmp1_ens(count1,m)
                            endif
                            count1 = count1 + 1
                         endif
                      enddo
                   enddo
                enddo
                if(PRESENT(dim1)) then 
                   iret = nf90_put_var(ftn,varid,gtmp_ens(:,:,m),(/1,1,m,dim1/),&
                        (/LIS_rc%gnc(n),LIS_rc%gnr(n),1,1/))
                else            
                   iret = nf90_put_var(ftn,varid,gtmp_ens(:,:,m),(/1,1,m/),&
                        (/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
                endif
             enddo
             if(ftn_stats.ne.-1) then
                if ( LIS_rc%sout ) then
                   call stats_ens(gtmp_ens,LIS_rc%udef,               &
                        LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n), &
                        vmean,vstdev,vmin,vmax)
                   if(form==1) then
                      write(ftn_stats,999) mvar,vmean,vstdev,vmin,vmax
                   elseif(form==2) then
                      write(ftn_stats,998) mvar,vmean,vstdev,vmin,vmax
                   endif
                   flush(ftn_stats)
                endif
             endif
             deallocate(gtmp_ens)
          endif
          deallocate(gtmp1_ens)
       endif
    endif
#endif
998 FORMAT(1X,A18,4E14.3)
999 FORMAT(1X,A18,4F14.3)
  end subroutine writevar_netcdf_withstats_real

!BOP
! !ROUTINE: writeroutingvar_netcdf_real
! \label{writeroutingvar_netcdf_real}
! 
! !INTERFACE:
  subroutine writeroutingvar_netcdf_real(ftn,ftn_stats, n, var,varid, mvar,&
       form, nmodel_status,dim1)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: ftn
    integer, intent(in) :: ftn_stats
    integer             :: varid
    real, intent(in)    :: var(LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n))
    character (len=*)   :: mvar
    integer, intent(in) :: form
    integer, intent(in) :: nmodel_status
    integer, intent(in), optional :: dim1
!
! !DESCRIPTION:
!  Write a real variable to a netcdf output file with some diagnostic 
!  statistics written to a text file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the netcdf output file
!   \item [ftn\_stats]
!     unit number of the ASCII text statistics file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item[mvar]
!     name of the variable being written (will be used in the stats file)
!   \item[form]
!     format to be used in the stats file (1-decimal format, 
!     2-scientific format)
!   \item [flag]
!    option to determine if the variable needs to be written (1-write, 
!    0-do not write)
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[stats](\ref{stats}) \newline
!     call to compute the diagnostic statistics
!  \end{description}
!
!EOP
    integer             :: l, iret
    real :: vmean,vstdev,vmin,vmax

    real, allocatable :: var1(:)
    real, allocatable :: var1_ens(:,:)
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp_ens(:,:,:)
    real, allocatable :: gtmp1(:)
    real, allocatable :: gtmp1_ens(:,:)
    integer :: gdeltas
    integer :: count1 ,c,r,m,gid,ntiles,ierr,i,t

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    if(LIS_rc%wopt.eq."2d gridspace") then 
       allocate(var1(LIS_rc%nroutinggrid(n)))
       if(LIS_masterproc) then 
          allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(gtmp1(LIS_rc%glbnroutinggrid(n)))
          gtmp = 0.0
          gtmp1 = 0.0
       else
          allocate(gtmp1(1))
          gtmp1 = 0.0
       endif
       var1 = 0 

       do i=1,LIS_rc%nroutinggrid(n)
          do m=1,LIS_rc%nensem(n)
             t = m+(i-1)*LIS_rc%nensem(n)
             if ( var(t) == -9999.0 ) then
                var1(i) = -9999.0
             else
                var1(i) = var1(i) + &
                     var(t)*LIS_routing(n)%tile(t)%fgrd*&
                     LIS_routing(n)%tile(t)%pens
             endif
          enddo
       enddo
#if (defined SPMD)      
       gdeltas = LIS_routing_gdeltas(n,LIS_localPet)
       call MPI_GATHERV(var1,gdeltas,&
            MPI_REAL,gtmp1,LIS_routing_gdeltas(n,:),LIS_routing_goffsets(n,:),&
            MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
       gtmp1 = var1
#endif
       deallocate(var1)
       if(LIS_masterproc) then 
          gtmp = LIS_rc%udef
          count1=1
          do l=1,LIS_npes
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                   gid = r+(c-1)*LIS_rc%gnr(n)
                   ntiles = LIS_routing(n)%ntiles_pergrid(gid)
                   if(ntiles.ne.0) then                 
                      if(r.ge.LIS_nss_ind(n,l).and.&
                           r.le.LIS_nse_ind(n,l).and.&
                           c.ge.LIS_ews_ind(n,l).and.&
                           c.le.LIS_ewe_ind(n,l))then !points not in halo
                         gtmp(c,r) = gtmp1(count1)
                      endif
                      count1 = count1 + 1
                   endif
                enddo
             enddo
          enddo
          
          if(PRESENT(dim1)) then 
             iret = nf90_put_var(ftn,varid,gtmp,(/1,1,dim1/),&
                                 (/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
          else            
             iret = nf90_put_var(ftn,varid,gtmp,(/1,1/),&
                                 (/LIS_rc%gnc(n),LIS_rc%gnr(n)/))
          endif 
          if(ftn_stats.ne.-1) then
             if ( LIS_rc%sout ) then
                call stats(gtmp,LIS_rc%udef,LIS_rc%gnc(n)*LIS_rc%gnr(n),&
                           vmean,vstdev,vmin,vmax)
                if(form==1) then
                   write(ftn_stats,999) mvar,vmean,vstdev,vmin,vmax
                elseif(form==2) then
                   write(ftn_stats,998) mvar,vmean,vstdev,vmin,vmax
                endif
                flush(ftn_stats)
             endif
          endif
          deallocate(gtmp)
       endif
       deallocate(gtmp1)

    ! Write output in 2D ensemble grid space:
    elseif(LIS_rc%wopt.eq."2d ensemble gridspace") then 

       ! Non-model output field status (T=non-model; F=model-based):
       if(nmodel_status.ne.0) then   ! non-model output field status
          allocate(var1(LIS_rc%nroutinggrid(n)))
          if(LIS_masterproc) then 
             allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
             allocate(gtmp1(LIS_rc%glbnroutinggrid(n)))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1))
             gtmp1 = 0.0
          endif
          var1 = 0 
          do i=1,LIS_rc%nroutinggrid(n)
             do m=1,LIS_rc%nensem(n)
                t = m + (i-1)*LIS_rc%nensem(n)
                if ( var(t) == -9999.0 ) then
                   var1(i) = -9999.0
                else
                   var1(i) = var1(i) + &
                        var(t)*LIS_routing(n)%tile(t)%fgrd*&
                        LIS_routing(n)%tile(t)%pens
                endif
             enddo
          enddo
          
#if (defined SPMD)      
          gdeltas = LIS_routing_gdeltas(n,LIS_localPet)
          call MPI_GATHERV(var1,gdeltas,&
               MPI_REAL,gtmp1,LIS_routing_gdeltas(n,:),&
               LIS_routing_goffsets(n,:),&
               MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
          gtmp1 = var1
#endif
          deallocate(var1)
          if(LIS_masterproc) then 
             gtmp = LIS_rc%udef
             count1=1
             do l=1,LIS_npes
                do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                   do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                      gid = r+(c-1)*LIS_rc%gnr(n)
                      ntiles = LIS_routing(n)%ntiles_pergrid(gid)
                      if(ntiles.ne.0) then                 
                         if(r.ge.LIS_nss_ind(n,l).and.&
                              r.le.LIS_nse_ind(n,l).and.&
                              c.ge.LIS_ews_ind(n,l).and.&
                              c.le.LIS_ewe_ind(n,l))then !points not in halo
                            gtmp(c,r) = gtmp1(count1)
                         endif
                         count1 = count1 + 1
                      endif
                   enddo
                enddo
             enddo
             if(PRESENT(dim1)) then 
                iret = nf90_put_var(ftn,varid,gtmp,(/1,1,dim1/),&
                     (/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
             else            
                iret = nf90_put_var(ftn,varid,gtmp,(/1,1/),&
                     (/LIS_rc%gnc(n),LIS_rc%gnr(n)/))
             endif
             if(ftn_stats.ne.-1) then
                if ( LIS_rc%sout ) then
                   call stats(gtmp,LIS_rc%udef,LIS_rc%gnc(n)*LIS_rc%gnr(n),&
                              vmean,vstdev,vmin,vmax)
                   if(form==1) then
                      write(ftn_stats,999) mvar,vmean,vstdev,vmin,vmax
                   elseif(form==2) then
                      write(ftn_stats,998) mvar,vmean,vstdev,vmin,vmax
                   endif
                   flush(ftn_stats)
                endif
             endif
             deallocate(gtmp)
          endif
          deallocate(gtmp1)
          
       ! Model-based field output:
       else

          allocate(var1_ens(LIS_rc%nroutinggrid(n), LIS_rc%nensem(n)))
          allocate(var1(LIS_rc%nroutinggrid(n)))
          if(LIS_masterproc) then 
             allocate(gtmp_ens(LIS_rc%gnc(n),LIS_rc%gnr(n), LIS_rc%nensem(n)))
             allocate(gtmp1_ens(LIS_rc%glbnroutinggrid(n),LIS_rc%nensem(n)))
             gtmp_ens = 0.0
             gtmp1_ens = 0.0
          else
             allocate(gtmp1_ens(1,LIS_rc%nensem(n)))
             gtmp1_ens = 0.0
          endif

          var1_ens = 0 
          do i=1,LIS_rc%nroutinggrid(n)
             do m=1,LIS_rc%nensem(n)
                t = m + (i-1)*LIS_rc%nensem(n)
                if ( var(t) == -9999.0 ) then
                   var1_ens(i,m) = -9999.0
                else
                   var1_ens(i,m) = &
                        var(t)*LIS_domain(n)%tile(t)%fgrd
                endif
             enddo
          enddo
       
#if (defined SPMD)      
          gdeltas = LIS_routing_gdeltas(n,LIS_localPet)
          do m=1,LIS_rc%nensem(n)
             ! EMK: It is possible that the first dimension of var1_ens is 0 
             ! (no grid points with tiles in the PET).  Unfortunately, slicing
             ! such an array [e.g., var1_ens(:,m)] will cause an array bounds 
             ! error.  So, we add some defensive code here to (a) copy a 
             ! slice to a 1-d array only if the dimension is > 0; and (b) 
             ! always pass the 1d array to MPI_GATHERV.  Note that no memory 
             ! access error will occur in MPI_GATHERV for the zero-grid count 
             ! case as long as gdeltas is also zero.
             if (LIS_rc%nroutinggrid(n) > 0) then
                var1(:) = var1_ens(:,m)
             end if
             call MPI_GATHERV(var1,gdeltas,&
                  MPI_REAL,gtmp1_ens(:,m),LIS_routing_gdeltas(n,:),&
                  LIS_routing_goffsets(n,:),&
                  MPI_REAL,0,LIS_mpi_comm,ierr)
          enddo
#else 
          do m=1,LIS_rc%nensem(n)
             gtmp1_ens(:,m) = var1_ens(:,m)
          enddo
#endif

          deallocate(var1)     ! EMK...Avoid memory leak
          deallocate(var1_ens) ! EMK...Avoid memory leak

          if(LIS_masterproc) then
             gtmp_ens = LIS_rc%udef
             do m=1,LIS_rc%nensem(n)
                count1=1
                do l=1,LIS_npes
                   do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                      do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                         gid = r+(c-1)*LIS_rc%gnr(n)
                         ntiles = LIS_routing(n)%ntiles_pergrid(gid)
                         if(ntiles.ne.0) then                 
                            if(r.ge.LIS_nss_ind(n,l).and.&
                                 r.le.LIS_nse_ind(n,l).and.&
                                 c.ge.LIS_ews_ind(n,l).and.&
                                 c.le.LIS_ewe_ind(n,l))then !points not in halo
                               gtmp_ens(c,r,m) = gtmp1_ens(count1,m)
                            endif
                            count1 = count1 + 1
                         endif
                      enddo
                   enddo
                enddo
                if(PRESENT(dim1)) then 
                   iret = nf90_put_var(ftn,varid,gtmp_ens(:,:,m),(/1,1,m,dim1/),&
                        (/LIS_rc%gnc(n),LIS_rc%gnr(n),1,1/))
                else            
                   iret = nf90_put_var(ftn,varid,gtmp_ens(:,:,m),(/1,1,m/),&
                        (/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
                endif
             enddo
             if(ftn_stats.ne.-1) then
                if ( LIS_rc%sout ) then
                   call stats_ens(gtmp_ens,LIS_rc%udef,               &
                        LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n), &
                        vmean,vstdev,vmin,vmax)
                   if(form==1) then
                      write(ftn_stats,999) mvar,vmean,vstdev,vmin,vmax
                   elseif(form==2) then
                      write(ftn_stats,998) mvar,vmean,vstdev,vmin,vmax
                   endif
                   flush(ftn_stats)
                endif
             endif
             deallocate(gtmp_ens)
          endif
          deallocate(gtmp1_ens)
       endif
    endif
#endif
998 FORMAT(1X,A18,4E14.3)
999 FORMAT(1X,A18,4F14.3)
  end subroutine writeroutingvar_netcdf_real

!BOP
! !ROUTINE: writevar_grib1_withstats_real
! \label{writevar_grib1_withstats_real}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_grib1
subroutine writevar_grib1_withstats_real(ftn, ftn_stats, n,   &
                                         var, mvar, &
                                         gribid, gribSF, gribSfc, &
                                         gribLvl, sType, &
                                         time_unit, time_p1, time_p2, &
                                         timeRange, &
                                         form,     &
                                         toplev, botlev, kdim)

! !USES:
  use map_utils
  use LIS_pluginIndices
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in)         :: ftn
  integer, intent(in)         :: ftn_stats
  integer, intent(in)         :: n 
  real, intent(in)            :: var(LIS_rc%ntiles(n))
  character(len=*),intent(in) :: mvar
  integer, intent(in)         :: gribid
  integer, intent(in)         :: gribSF
  integer, intent(in)         :: gribSfc
  integer, intent(in)         :: gribLvl
  character(len=*), intent(in)  :: sType
  integer, intent(in)         :: timeRange
  integer, intent(in)         :: form
  integer, intent(in)         :: time_unit
  integer, intent(in)         :: time_p1
  integer, intent(in)         :: time_p2
  integer                     :: kdim
  real,    intent(in)         :: toplev(kdim), botlev(kdim)
! !DESCRIPTION:
!  Write a real variable to a GRIB-1 output file with some diagnostic 
!  statistics written to a text file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [ftn\_stats]
!     unit number of the ASCII text statistics file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item[mvar]
!     name of the variable being written (will be used in the stats file)
!   \item [flag]
!    option to determine if the variable needs to be written (1-write, 
!    0-do not write)
!   \item[form]
!     format to be used in the stats file (1-decimal format, 
!     2-scientific format)
!   \item [info]
!     Array containing the sizes of the various GRIB-1 sections
!   \item [sect1]
!     Array for storing the contents of Section 1.
!   \item [s1opt]
!     Array for storing optional Section 1 information (not used)
!   \item [sect2]
!     Array for storing the contents of Section 2.
!   \item [sect3]
!     Array for storing the contents of Section 3.
!   \item [toplev]
!     Array of the heights of the top of each layer
!   \item [botlev]
!     Array of the heights of the bottom of each layer
!   \item [kdim]
!     Grid dimension k-direction (number of levels/layers)
!   \item [size\_3]
!     Size of the sect3 array
!  \end{description} 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_gather\_gridded\_output](\ref{LIS_gather_gridded_output})
!     call to gather the 1d tiled output variable into a 2d gridded array
!   \item[write\_stats](\ref{write_stats})
!    call to compute diagnostic statistics of a variable
!  \end{description}
!EOP
  character*255        :: message(20)
  integer              :: igrib
  character*8          :: date
  integer              :: idate,idate1
  real, allocatable        :: gtmp(:,:)
  real, allocatable    :: gtmp1(:)
  integer              :: c,r,count1
  integer              :: gid,ntiles
  real                 :: lat_ur, lon_ur
  real                 :: lat_ll, lon_ll
  integer              :: yr1, mo1,da1,hr1,mn1
  integer              :: iret,decimalPrecision,gribSFtemp

  call LIS_gather_gridded_output(n, gtmp, var)

#if (defined USE_GRIBAPI)   
!  LIS_rc%grib_table        = GRIB_Table_Version
!  LIS_rc%grib_center_id    = GRIB_Center_Id
!  LIS_rc%grib_subcenter_id = GRIB_Subcenter_Id
!  LIS_rc%grib_grid_id      = GRIB_Grid_Id
!  LIS_rc%grib_process_id   = GRIB_Process_Id

  !putting it back in the right order into gtmp1          
  if ( LIS_masterproc ) then
     allocate(gtmp1(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
     do r=1,LIS_rc%gnr(n)
        do c=1,LIS_rc%gnc(n)
           gid = c+(r-1)*LIS_rc%gnc(n)
           gtmp1(gid) = gtmp(c,r)
        enddo
     enddo

      ! Note passing string of defined points only to output
      ! because bitmap in GRIB-1 file will fill in the rest
     
#if (defined USE_ECCODES)   
     call grib_new_from_samples(igrib,"GRIB1",iret)
#else
     call grib_new_from_template(igrib,"GRIB1",iret)
#endif
     call LIS_verify(iret, 'grib_new_from_template failed in LIS_historyMod')
     
     call grib_set(igrib,'table2Version',LIS_rc%grib_table,iret)
     call LIS_verify(iret,'grib_set:table2version failed in LIS_historyMod')
     
     call grib_set(igrib,'generatingProcessIdentifier',LIS_rc%grib_process_id,iret)
     call LIS_verify(iret,'grib_set:generatingProcessIdentifier failed in LIS_historyMod')
     
     call grib_set(igrib,'gridDefinition',LIS_rc%grib_grid_id,iret)
     call LIS_verify(iret,'grib_set:grid ID failed in LIS_historyMod')
     
     call grib_set(igrib,'indicatorOfParameter',gribid, iret)
     call LIS_verify(iret,'grib_set:indicatorOfParameter failed in LIS_historyMod')
     
 !    call grib_set(igrib,'paramId',gribid, iret)
 !    call LIS_verify(iret,'grib_set:paramId failed in LIS_historyMod')
     
     call grib_set(igrib,'indicatorOfTypeOfLevel',gribSfc, iret)
     call LIS_verify(iret,'grib_set:indicatorOfTypeOfLevel failed in LIS_historyMod')
     
     call grib_set(igrib,'level',gribLvl, iret)
     call LIS_verify(iret,'grib_set:level failed in LIS_historyMod')
     
     call grib_set(igrib,'topLevel',toplev(1), iret)
     call LIS_verify(iret,'grib_set:topLevel failed in LIS_historyMod')

     call grib_set(igrib,'bottomLevel',botlev(1), iret)
     call LIS_verify(iret,'grib_set:bottomLevel failed in LIS_historyMod')

     call grib_set(igrib,'stepType',sType, iret)
     call LIS_verify(iret,'grib_set:stepType failed in LIS_historyMod')

     call grib_set(igrib,'stepUnits',time_unit, iret)
     call LIS_verify(iret,'grib_set:stepUnits failed in LIS_historyMod')

     call grib_set(igrib,'startStep',time_p1, iret)
     call LIS_verify(iret,'grib_set:startStep failed in LIS_historyMod')

     call grib_set(igrib,'endStep',time_p2, iret)
     call LIS_verify(iret,'grib_set:endStep failed in LIS_historyMod')
     
     call grib_set(igrib,'timeRangeIndicator',timeRange, iret)
     call LIS_verify(iret,'grib_set:timeRangeIndicator failed in LIS_historyMod')

     call grib_set(igrib,'swapScanningLat',1, iret)
     call LIS_verify(iret,'grib_set:swapScanningLat failed in LIS_historyMod')

     call grib_set(igrib,'Ni',LIS_rc%gnc(n),iret)
     call LIS_verify(iret, 'grib_set:Ni failed in LIS_historyMod')
     
     call grib_set(igrib,'Nj',LIS_rc%gnr(n),iret)
     call LIS_verify(iret, 'grib_set:Ni failed in LIS_historyMod')
     
     call ij_to_latlon(LIS_domain(n)%lisproj,float(LIS_rc%gnc(n)),&
          float(LIS_rc%gnr(n)),lat_ur,lon_ur)      
     call ij_to_latlon(LIS_domain(n)%lisproj,1.0, 1.0, &
          lat_ll,lon_ll)      
     
     call grib_set(igrib, 'latitudeOfFirstGridPointInDegrees',lat_ll,iret)
     call LIS_verify(iret, 'grib_set:latitudeOfFirstGridPointInDegrees failed in LIS_historyMod')
     
     call grib_set(igrib, 'longitudeOfFirstGridPointInDegrees',lon_ll,iret)
     call LIS_verify(iret, 'grib_set:longitudeOfFirstGridPointInDegrees failed in LIS_historyMod')
     
     call grib_set(igrib, 'latitudeOfLastGridPointInDegrees',lat_ur,iret)
     call LIS_verify(iret, 'grib_set:latitudeOfLastGridPointInDegrees failed in LIS_historyMod')
     
     call grib_set(igrib, 'longitudeOfLastGridPointInDegrees',lon_ur,iret)
     call LIS_verify(iret, 'grib_set:longitudeOfLastGridPointInDegrees failed in LIS_historyMod')
     
     call grib_set(igrib, 'missingValue',LIS_rc%udef,iret)
     call LIS_verify(iret, 'grib_set:missingValue failed in LIS_historyMod')
     
! Should not need to fix the "num bits" value for each parameter
! if the "decimalPrecision" (aka, "DecScale") is set properly. - dmm
!     call grib_set(igrib, 'bitsPerValue',12,iret)
!     call LIS_verify(iret, 'grib_set:bitsPerValue failed in LIS_historyMod')

! Set the "decimalPrecision" (aka, "DecScale") based on the
! gribSF (grib scale factor) set in the MODEL OUTPUT TBL. - dmm
     gribSFtemp = gribSF
     decimalPrecision = 0
     do while (gribSFtemp.ge.10)
        decimalPrecision = decimalPrecision + 1
        gribSFtemp = gribSFtemp / 10
     enddo
     call grib_set(igrib, 'decimalPrecision',decimalPrecision,iret)
     call LIS_verify(iret, 'grib_set:decimalPrecision failed in LIS_historyMod')

     call grib_set(igrib, 'bitmapPresent',1,iret)
     call LIS_verify(iret, 'grib_set:bitmapPresent failed in LIS_historyMod')
               
     if (LIS_rc%lis_map_proj.eq."latlon") then 
        call grib_set(igrib,'gridType','regular_ll',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
        
        call grib_set(igrib,'iDirectionIncrementInDegrees',LIS_rc%gridDesc(n,9),iret)
        call LIS_verify(iret,'grib_set:iDirectionIncrementInDegrees failed in LIS_historyMod')
        
        call grib_set(igrib,'jDirectionIncrementInDegrees',LIS_rc%gridDesc(n,10),iret)
        call LIS_verify(iret,'grib_set:jDirectionIncrementInDegrees failed in LIS_historyMod')
        
     elseif (LIS_rc%lis_map_proj.eq. "polar") then ! for polar stereographic
        call grib_set(igrib,'gridType','polar_stereographic',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
        ! IF LAMBERT CONFORMAL PROJECTION  
     elseif (LIS_rc%lis_map_proj .eq."lambert") then
        call grib_set(igrib,'gridType','lambert',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
        call grib_set(igrib,'Nx',LIS_rc%gnc(n),iret)
        call LIS_verify(iret,'grib_set: Nx failed in LIS_historyMod')
        call grib_set(igrib,'Ny',LIS_rc%gnr(n),iret)
        call LIS_verify(iret,'grib_set: Ny failed in LIS_historyMod')
        call grib_set(igrib,'latitudeOfFirstGridPoint',&
             LIS_domain(n)%stlat*1000,iret)
        call LIS_verify(iret,'grib_set: latitudeOfFirstGridPoint failed in LIS_historyMod')
        call grib_set(igrib,'longitudeOfFirstGridPoint',&
             LIS_domain(n)%stlon*1000,iret)
        call LIS_verify(iret,'grib_set: longitudeOfFirstGridPoint failed in LIS_historyMod')
        call grib_set(igrib,'LoV',LIS_domain(n)%truelon*1000,iret)
        call LIS_verify(iret,'grib_set: LoV failed in LIS_historyMod')
        call grib_set(igrib,'DxInMetres',LIS_domain(n)%dx,iret)
        call LIS_verify(iret,'grib_set: DxInMetres failed in LIS_historyMod')
        call grib_set(igrib,'DyInMetres',LIS_domain(n)%dy,iret)
        call LIS_verify(iret,'grib_set: DyInMetres failed in LIS_historyMod')
        call grib_set(igrib,'Latin1',LIS_domain(n)%truelat1*1000,iret)
        call LIS_verify(iret,'grib_set: Latin1 failed in LIS_historyMod')
        call grib_set(igrib,'Latin2',LIS_domain(n)%truelat2*1000,iret)
        call LIS_verify(iret,'grib_set: Latin2 failed in LIS_historyMod')

        ! These are the default values for earthIsOblate, uvRelativeToGrid,
        ! iScansNegatively, jScansPositively, and jPointsAreConsecutive.
        ! These defaults seem reasonable for LIS.  They are included here
        ! to expose them, should they be needed in the future.
        !
        ! Assume that the Earth is a sphere
        call grib_set(igrib,'earthIsOblate',0,iret)
        call LIS_verify(iret,'grib_set: earthIsOblate failed in LIS_historyMod')
        ! u and v components are easterly and northerly
        call grib_set(igrib,'uvRelativeToGrid',0,iret)
        call LIS_verify(iret,'grib_set: uvRelativeToGrid failed in LIS_historyMod')
        call grib_set(igrib,'iScansNegatively',0,iret)
        call LIS_verify(iret,'grib_set: iScansNegatively failed in LIS_historyMod')
        call grib_set(igrib,'jScansPositively',1,iret)
        call LIS_verify(iret,'grib_set: jScansPositively failed in LIS_historyMod')
        call grib_set(igrib,'jPointsAreConsecutive',0,iret)
        call LIS_verify(iret,'grib_set: jPointsAreConsecutive failed in LIS_historyMod')
        ! I am not sure how to set projectionCenterFlag.  Its default
        ! value is 0.
        call grib_set(igrib,'projectionCenterFlag',0,iret)
        call LIS_verify(iret,'grib_set: projectionCenterFlag failed in LIS_historyMod')
        ! I am not sure how to set latitudeOfSouthernPole and
        ! longitudeOfSouthernPole.  They default to 0.
        !call grib_set(igrib,'latitudeOfSouthernPole',0,iret)
        !call grib_set(igrib,'longitudeOfSouthernPole',0,iret)

        ! IF MERCATOR PROJECTION  
     elseif (LIS_rc%lis_map_proj .eq. "mercator") then
        
        call grib_set(igrib,'gridType','mercator',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
        
     elseif (LIS_rc%lis_map_proj .eq. "UTM") then
        
        call grib_set(igrib,'gridType','UTM',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
     elseif (LIS_rc%lis_map_proj .eq. "gaussian") then
        
        call grib_set(igrib,'gridType','regular_gg',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
     else  !Unsupported Map Projection for GRIB output
        
        message(1)='program:  LIS_historyMod'
        message(2)=' subroutine:  writevar_grib1_withstats_real'
        message(3)='  Unsupported map projection for GRIB1 output!'
        call lis_abort(message)
        stop
        
     endif
     
     da1=LIS_rc%da
     mo1=LIS_rc%mo
     yr1=LIS_rc%yr
     
     write(unit=date,fmt='(i4.4,i2.2,i2.2)') yr1,mo1,da1
     read(date,'(I8)') idate
     
     call grib_set(igrib,'dataDate',idate,iret)
     call LIS_verify(iret, 'grib_set:dataDate failed in LIS_historyMod')

     hr1=LIS_rc%hr
     mn1=LIS_rc%mn

     write(unit=date,fmt='(i2.2,i2.2)') hr1,mn1
     read(date,'(I4)') idate1

     call grib_set(igrib,'dataTime',idate1,iret)
     call LIS_verify(iret, 'grib_set:dataTime failed in LIS_historyMod')

     call grib_set(igrib,'values',gtmp1,iret)
     call LIS_verify(iret, 'grib_set:values failed in LIS_historyMod')
     
! Move setting of centre and subCentre to the end of the settings.
! The order these are written is important and will affect output. - dmm
     call grib_set(igrib,'centre',LIS_rc%grib_center_id,iret)
     call LIS_verify(iret,'grib_set:centre failed in LIS_historyMod')
     
     call grib_set(igrib,'subCentre',LIS_rc%grib_subcenter_id,iret)
     call LIS_verify(iret,'grib_set:subCentre failed in LIS_historyMod')
     
     call grib_write(igrib,ftn,iret)
     call LIS_verify(iret, 'grib_write failed in LIS_historyMod')
     
     call grib_release(igrib,iret)
     call LIS_verify(iret,'grib_release failed in LIS_historyMod')

     call write_stats(gtmp, LIS_rc%gnc(n)*LIS_rc%gnr(n), &
          mvar, ftn_stats, form)
     
     deallocate(gtmp)
     deallocate(gtmp1)
  endif
#endif
          
end subroutine writevar_grib1_withstats_real

!BOP
! !ROUTINE: writevar_grib2_withstats_real
! \label{writevar_grib2_withstats_real}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_grib2
subroutine writevar_grib2_withstats_real(ftn, ftn_stats, n,   &
                                         var, mvar, &
                                         gribid, gribSF, gribSfc, &
                                         gribDis, gribCat,  &
                                         sType, &
                                         time_unit, timeRange, pdTemplate, &
                                         form,     &
                                         toplev, botlev, kdim, depscale)

! !USES:
  use map_utils
  use LIS_pluginIndices
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in)         :: ftn
  integer, intent(in)         :: ftn_stats
  integer, intent(in)         :: n 
  real, intent(in)            :: var(LIS_rc%ntiles(n))
  character(len=*),intent(in) :: mvar
  integer, intent(in)         :: gribid
  integer, intent(in)         :: gribSF
  integer, intent(in)         :: gribSfc
  integer, intent(in)         :: gribDis
  integer, intent(in)         :: gribCat
  character(len=*), intent(in)  :: sType
  integer, intent(in)         :: form
  integer, intent(in)         :: time_unit
  integer, intent(in)         :: timeRange
  integer, intent(in)         :: pdTemplate
  integer                     :: kdim
  real,    intent(in)         :: toplev(kdim), botlev(kdim)
  integer, intent(in)         :: depscale
! !DESCRIPTION:
!  Write a real variable to a GRIB-2 output file with some diagnostic 
!  statistics written to a text file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [ftn\_stats]
!     unit number of the ASCII text statistics file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item[mvar]
!     name of the variable being written (will be used in the stats file)
!   \item[gribid]
!    paremter number per discipline per category 4.2.?.?.table in Section 4
!   \item[gribSF]
!    Decimal scale factor in Section 5.
!   \item[gribSfc]
!    Fixed surface types and units in 4.5.table in Section 4. 
!   \item[sType]
!    product description used in 4.10.table in Section 4.
!   \item[pdTemplate]
!    Product definition template number 4.0.table in Section 4.
!    separates 'inst' and 'avg/accum/max/min' description.
!   \item[time\_unit]
!    Indicator of unit of time range in 4.4.table in Section 4.
!   \item[form]
!     format to be used in the stats file (1-decimal format, 
!     2-scientific format)
!   \item [toplev]
!     Array of the heights of the top of each layer
!   \item [botlev]
!     Array of the heights of the bottom of each layer
!   \item [kdim]
!     Grid dimension k-direction (number of levels/layers)
!   \item [depscale]
!     scale for soil layer depth (0=meters, 2=cm)
!  \end{description} 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_gather\_gridded\_output](\ref{LIS_gather_gridded_output})
!     call to gather the 1d tiled output variable into a 2d gridded array
!   \item[write\_stats](\ref{write_stats})
!    call to compute diagnostic statistics of a variable
!  \end{description}
!EOP
  character*255        :: message(20)
  integer              :: igrib
  character*8          :: date
  integer              :: idate,idate1
  real, allocatable        :: gtmp(:,:)
  real, allocatable    :: gtmp1(:)
  integer              :: c,r,count1
  integer              :: gid,ntiles
  real                 :: lat_ur, lon_ur
  real                 :: lat_ll, lon_ll
  integer              :: yr1, mo1,da1,hr1,mn1
  integer              :: iret,decimalPrecision,gribSFtemp

  call LIS_gather_gridded_output(n, gtmp, var)

#if (defined USE_GRIBAPI)   
!  LIS_rc%grib_table        = GRIB_Table_Version
!  LIS_rc%grib_center_id    = GRIB_Center_Id
!  LIS_rc%grib_subcenter_id = GRIB_Subcenter_Id
!  LIS_rc%grib_grid_id      = GRIB_Grid_Id
!  LIS_rc%grib_process_id   = GRIB_Process_Id

  !putting it back in the right order into gtmp1          
  if ( LIS_masterproc ) then
     allocate(gtmp1(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
     do r=1,LIS_rc%gnr(n)
        do c=1,LIS_rc%gnc(n)
           gid = c+(r-1)*LIS_rc%gnc(n)
           gtmp1(gid) = gtmp(c,r)
        enddo
     enddo

      ! Note passing string of defined points only to output
      ! because bitmap in GRIB-2 file will fill in the rest.

      ! Populate message with a sample file initially.
     call grib_new_from_samples(igrib,"GRIB2",iret)
!     call grib_new_from_samples(igrib,"reduced_gg_sfc_jpeg_grib2",iret)
     call LIS_verify(iret, 'grib_new_from_sample failed in LIS_historyMod')

! Section 0: Indicator     

     call grib_set(igrib,'discipline',gribDis,iret)
     call LIS_verify(iret,'grib_set:discipline failed in LIS_historyMod')

! Section 1: Identification

     call grib_set(igrib,'centre',LIS_rc%grib_center_id,iret)
     call LIS_verify(iret,'grib_set:centre failed in LIS_historyMod')
     call grib_set(igrib,'subCentre',LIS_rc%grib_subcenter_id,iret)
     call LIS_verify(iret,'grib_set:subCentre failed in LIS_historyMod')
     call grib_set(igrib,'tablesVersion',LIS_rc%grib_table,iret)
     call LIS_verify(iret,'grib_set:table2version failed in LIS_historyMod')
      ! Hard-coded: significance of reference time 0=Analysis, 1.2.table
     call grib_set(igrib,'significanceOfReferenceTime',0,iret)
     call LIS_verify(iret,'grib_set:significanceOfReferenceTime failed in LIS_historyMod')
      ! Reference time is output file time stamp
     da1=LIS_rc%da
     mo1=LIS_rc%mo
     yr1=LIS_rc%yr
     
     write(unit=date,fmt='(i4.4,i2.2,i2.2)') yr1,mo1,da1
     read(date,'(I8)') idate
     
     call grib_set(igrib,'dataDate',idate,iret)
     call LIS_verify(iret, 'grib_set:dataDate failed in LIS_historyMod')

     hr1=LIS_rc%hr
     mn1=LIS_rc%mn

     write(unit=date,fmt='(i2.2,i2.2)') hr1,mn1
     read(date,'(I4)') idate1

     call grib_set(igrib,'dataTime',idate1,iret)
     call LIS_verify(iret, 'grib_set:dataTime failed in LIS_historyMod')
      ! Hard-coded: production status of data 2=Research products, 1.3.table
     call grib_set(igrib,'productionStatusOfProcessedData',2,iret)
     call LIS_verify(iret,'grib_set:productionStatusOfProcessedData failed in LIS_historyMod')
      ! Hard-coded: type of data 0=Analysis, 1.4.table
     call grib_set(igrib,'typeOfProcessedData',0,iret)
     call LIS_verify(iret,'grib_set:typeOfProcessedData failed in LIS_historyMod')
     call grib_set(igrib,'stepType',sType, iret)
     call LIS_verify(iret,'grib_set:stepType failed in LIS_historyMod')

! Section 2: Local Use Section (Optional) --none for now

! Section 3: Grid

      ! Grid definition template number 0=Lat-long; 3.1.table 
     call grib_set(igrib,'gridDefinitionTemplateNumber',LIS_rc%grib_grid_id,iret)
     call LIS_verify(iret,'grib_set:gridDefinitionTemplateNumber failed in LIS_historyMod')
      ! Hard-coded: shape of the Earth 0=radius = 6,367,470.0 m; 3.2.table 
     call grib_set(igrib,'shapeOfTheEarth',0,iret)
     call LIS_verify(iret,'grib_set:shapeOfTheEarth failed in LIS_historyMod')

     call grib_set(igrib,'swapScanningLat',1, iret)
     call LIS_verify(iret,'grib_set:swapScanningLat failed in LIS_historyMod')

     call grib_set(igrib,'Ni',LIS_rc%gnc(n),iret)
     call LIS_verify(iret, 'grib_set:Ni failed in LIS_historyMod')
     
     call grib_set(igrib,'Nj',LIS_rc%gnr(n),iret)
     call LIS_verify(iret, 'grib_set:Ni failed in LIS_historyMod')
     
     call ij_to_latlon(LIS_domain(n)%lisproj,float(LIS_rc%gnc(n)),&
          float(LIS_rc%gnr(n)),lat_ur,lon_ur)      
     call ij_to_latlon(LIS_domain(n)%lisproj,1.0, 1.0, &
          lat_ll,lon_ll)      
     
     call grib_set(igrib,'latitudeOfFirstGridPointInDegrees',lat_ll,iret)
     call LIS_verify(iret, 'grib_set:latitudeOfFirstGridPointInDegrees failed in LIS_historyMod')
     
     call grib_set(igrib,'longitudeOfFirstGridPointInDegrees',lon_ll,iret)
     call LIS_verify(iret, 'grib_set:longitudeOfFirstGridPointInDegrees failed in LIS_historyMod')
     
     call grib_set(igrib,'latitudeOfLastGridPointInDegrees',lat_ur,iret)
     call LIS_verify(iret, 'grib_set:latitudeOfLastGridPointInDegrees failed in LIS_historyMod')
     
     call grib_set(igrib,'longitudeOfLastGridPointInDegrees',lon_ur,iret)
     call LIS_verify(iret, 'grib_set:longitudeOfLastGridPointInDegrees failed in LIS_historyMod')

     if (LIS_rc%lis_map_proj.eq."latlon") then 
        call grib_set(igrib,'gridType','regular_ll',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
        
        call grib_set(igrib,'iDirectionIncrementInDegrees',LIS_rc%gridDesc(n,9),iret)
        call LIS_verify(iret,'grib_set:iDirectionIncrementInDegrees failed in LIS_historyMod')
        
        call grib_set(igrib,'jDirectionIncrementInDegrees',LIS_rc%gridDesc(n,10),iret)
        call LIS_verify(iret,'grib_set:jDirectionIncrementInDegrees failed in LIS_historyMod')
        
     elseif (LIS_rc%lis_map_proj.eq. "polar") then ! for polar stereographic
        call grib_set(igrib,'gridType','polar_stereographic',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
        ! IF LAMBERT CONFORMAL PROJECTION  
     elseif (LIS_rc%lis_map_proj .eq."lambert") then
        call grib_set(igrib,'gridType','lambert',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
        call grib_set(igrib,'Nx',LIS_rc%gnc(n),iret)
        call LIS_verify(iret,'grib_set: Nx failed in LIS_historyMod')
        call grib_set(igrib,'Ny',LIS_rc%gnr(n),iret)
        call LIS_verify(iret,'grib_set: Ny failed in LIS_historyMod')
        call grib_set(igrib,'latitudeOfFirstGridPoint',&
             LIS_domain(n)%stlat*1000,iret)
        call LIS_verify(iret,'grib_set: latitudeOfFirstGridPoint failed in LIS_historyMod')
        call grib_set(igrib,'longitudeOfFirstGridPoint',&
             LIS_domain(n)%stlon*1000,iret)
        call LIS_verify(iret,'grib_set: longitudeOfFirstGridPoint failed in LIS_historyMod')
        call grib_set(igrib,'LoV',LIS_domain(n)%truelon*1000,iret)
        call LIS_verify(iret,'grib_set: LoV failed in LIS_historyMod')
        call grib_set(igrib,'DxInMetres',LIS_domain(n)%dx,iret)
        call LIS_verify(iret,'grib_set: DxInMetres failed in LIS_historyMod')
        call grib_set(igrib,'DyInMetres',LIS_domain(n)%dy,iret)
        call LIS_verify(iret,'grib_set: DyInMetres failed in LIS_historyMod')
        call grib_set(igrib,'Latin1',LIS_domain(n)%truelat1*1000,iret)
        call LIS_verify(iret,'grib_set: Latin1 failed in LIS_historyMod')
        call grib_set(igrib,'Latin2',LIS_domain(n)%truelat2*1000,iret)
        call LIS_verify(iret,'grib_set: Latin2 failed in LIS_historyMod')

        ! These are the default values for earthIsOblate, uvRelativeToGrid,
        ! iScansNegatively, jScansPositively, and jPointsAreConsecutive.
        ! These defaults seem reasonable for LIS.  They are included here
        ! to expose them, should they be needed in the future.
        !
        ! Assume that the Earth is a sphere
        call grib_set(igrib,'earthIsOblate',0,iret)
        call LIS_verify(iret,'grib_set: earthIsOblate failed in LIS_historyMod')
        ! u and v components are easterly and northerly
        call grib_set(igrib,'uvRelativeToGrid',0,iret)
        call LIS_verify(iret,'grib_set: uvRelativeToGrid failed in LIS_historyMod')
        call grib_set(igrib,'iScansNegatively',0,iret)
        call LIS_verify(iret,'grib_set: iScansNegatively failed in LIS_historyMod')
        call grib_set(igrib,'jScansPositively',1,iret)
        call LIS_verify(iret,'grib_set: jScansPositively failed in LIS_historyMod')
        call grib_set(igrib,'jPointsAreConsecutive',0,iret)
        call LIS_verify(iret,'grib_set: jPointsAreConsecutive failed in LIS_historyMod')
        ! I am not sure how to set projectionCenterFlag.  Its default
        ! value is 0.
        call grib_set(igrib,'projectionCenterFlag',0,iret)
        call LIS_verify(iret,'grib_set: projectionCenterFlag failed in LIS_historyMod')
        ! I am not sure how to set latitudeOfSouthernPole and
        ! longitudeOfSouthernPole.  They default to 0.
        !call grib_set(igrib,'latitudeOfSouthernPole',0,iret)
        !call grib_set(igrib,'longitudeOfSouthernPole',0,iret)

        ! IF MERCATOR PROJECTION  
     elseif (LIS_rc%lis_map_proj .eq. "mercator") then
        
        call grib_set(igrib,'gridType','mercator',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
        
     elseif (LIS_rc%lis_map_proj .eq. "UTM") then
        
        call grib_set(igrib,'gridType','UTM',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
     elseif (LIS_rc%lis_map_proj .eq. "gaussian") then
        
        call grib_set(igrib,'gridType','regular_gg',iret)
        call LIS_verify(iret,'grib_set: gridType failed in LIS_historyMod')
     else  !Unsupported Map Projection for GRIB output
        
        message(1)='program:  LIS_historyMod'
        message(2)=' subroutine:  writevar_grib2_withstats_real'
        message(3)='  Unsupported map projection for GRIB1 output!'
        call lis_abort(message)
        stop
        
     endif
     
! Section 4: Product Definition Section

     call grib_set(igrib,'productDefinitionTemplateNumber',pdTemplate, iret)
     call LIS_verify(iret,'grib_set:productDefinitionTemplateNumber failed in LIS_historyMod')
     call grib_set(igrib,'parameterCategory',gribCat, iret)
     call LIS_verify(iret,'grib_set:parameterCategory failed in LIS_historyMod')
     call grib_set(igrib,'parameterNumber',gribid, iret)
     call LIS_verify(iret,'grib_set:parameterNumber failed in LIS_historyMod')
     !Hard-coded: Type of generating process 0=Analysis, 4.3.table
     call grib_set(igrib,'typeOfGeneratingProcess',0, iret)
     call LIS_verify(iret,'grib_set:typeOfGeneratingProcess failed in LIS_historyMod')
     call grib_set(igrib,'generatingProcessIdentifier',LIS_rc%grib_process_id,iret)
     call LIS_verify(iret,'grib_set:generatingProcessIdentifier failed in LIS_historyMod')
     call grib_set(igrib,'indicatorOfUnitOfTimeRange',time_unit, iret)
     call LIS_verify(iret,'grib_set:indicatorOfUnitOfTimeRange failed in LIS_historyMod')
     if ( sType .eq. 'instant' ) then 
        ! no need to specify statistical information template.4.statistcal.def
     else
      call grib_set(igrib,'typeOfStatisticalProcessing',sType, iret)
      call LIS_verify(iret,'grib_set:typeOfStatisticalProcessing failed in LIS_historyMod')
      ! Hard-coded: Type of time increment between successive fields used in the
      ! statistical processing (4.11.table)
      ! 6=local use for the backward (past) avg/accum/max/min in LIS
      call grib_set(igrib,'typeOfTimeIncrement',6, iret)
      call LIS_verify(iret,'grib_set:typeOfTimeIncrement failed in LIS_historyMod')
      call grib_set(igrib,'indicatorOfUnitForTimeRange',time_unit, iret)
      call LIS_verify(iret,'grib_set:indicatorOfUnitForTimeRange failed in LIS_historyMod')
      call grib_set(igrib,'lengthOfTimeRange',timeRange, iret)
      call LIS_verify(iret,'grib_set:lengthOfTimeRange failed in LIS_historyMod')
      call grib_set(igrib,'indicatorOfUnitForTimeIncrement',time_unit, iret)
      call LIS_verify(iret,'grib_set:indicatorOfUnitForTimeIncrement failed in LIS_historyMod')
      ! Hard-coded: Time increment between successive fields, in units defined
      ! by the previous Octets.  0 means that the statistical processing is the 
      ! result of a continuous processes, not the processing of a number of
      ! descrete samples.
      call grib_set(igrib,'timeIncrement',0, iret)
      call LIS_verify(iret,'grib_set:timeIncrement failed in LIS_historyMod')
     endif
     
     ! First and second surface types need to be set prior to scale and depth
     ! topLevel and bottomLevel or level will be derived by the settings below.
     ! Levels other than Surface and Soil layers are not considered yet--hkb
     call grib_set(igrib,'typeOfFirstFixedSurface',gribSfc, iret)
     call LIS_verify(iret,'grib_set:typeOfFirstFixedSurface failed in LIS_historyMod')
     if ( gribSfc .eq. 106 ) then   ! soil layers
      call grib_set(igrib,'typeOfSecondFixedSurface',gribSfc, iret)
      call LIS_verify(iret,'grib_set:typeOfFirstFixedSurface failed in LIS_historyMod')
     ! Removed Hard-coded scale: depscale is passed from writeSingleGrib2Var
      call grib_set(igrib,'scaleFactorOfFirstFixedSurface',depscale, iret)
      call LIS_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LIS_historyMod')
      call grib_set(igrib,'scaledValueOfFirstFixedSurface',toplev(1), iret)
      call LIS_verify(iret,'grib_set:scaledValueOfFirstFixedSurface failed in LIS_historyMod')
      call grib_set(igrib,'scaleFactorOfSecondFixedSurface',depscale, iret)
      call LIS_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LIS_historyMod')
      call grib_set(igrib,'scaledValueOfSecondFixedSurface',botlev(1), iret)
      call LIS_verify(iret,'grib_set:scaledValueOfSecondFixedSurface failed in LIS_historyMod')
     elseif ( gribSfc .eq. 1 ) then    ! surface
      call grib_set(igrib,'typeOfSecondFixedSurface',255, iret)
      call LIS_verify(iret,'grib_set:typeOfFirstFixedSurface failed in LIS_historyMod')
      call grib_set(igrib,'scaleFactorOfFirstFixedSurface',0, iret)
      call LIS_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LIS_historyMod')
      call grib_set(igrib,'scaledValueOfFirstFixedSurface',toplev(1), iret)
      call LIS_verify(iret,'grib_set:scaledValueOfFirstFixedSurface failed in LIS_historyMod')

      call grib_set(igrib,'scaleFactorOfSecondFixedSurface',255, iret)
      call LIS_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LIS_historyMod')
      call grib_set(igrib,'scaledValueOfSecondFixedSurface',255, iret)
      call LIS_verify(iret,'grib_set:scaledValueOfSecondFixedSurface failed in LIS_historyMod')
     else   ! 114 (snow level) or old 112 ??
       write(LIS_logunit,*) 'Warning: special surface type !! '//&
                     'verify scale/depth for ',gribSfc
      call grib_set(igrib,'typeOfSecondFixedSurface',gribSfc, iret)
      call LIS_verify(iret,'grib_set:typeOfFirstFixedSurface failed in LIS_historyMod')
      call grib_set(igrib,'scaleFactorOfFirstFixedSurface',0, iret)
      call LIS_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LIS_historyMod')
      call grib_set(igrib,'scaledValueOfFirstFixedSurface',toplev(1), iret)
      call LIS_verify(iret,'grib_set:scaledValueOfFirstFixedSurface failed in LIS_historyMod')

      call grib_set(igrib,'scaleFactorOfSecondFixedSurface',0, iret)
      call LIS_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LIS_historyMod')
      call grib_set(igrib,'scaledValueOfSecondFixedSurface',botlev(1), iret)
      call LIS_verify(iret,'grib_set:scaledValueOfSecondFixedSurface failed in LIS_historyMod')
     endif

! Section 5: Data Representation

     call grib_set(igrib,'packingType',LIS_rc%grib_packing_type,iret)
     call LIS_verify(iret, 'grib_set:packingType failed in LIS_historyMod')
     call grib_set(igrib, 'missingValue',LIS_rc%udef,iret)
     call LIS_verify(iret, 'grib_set:missingValue failed in LIS_historyMod')
     
! Should not need to fix the "num bits" value for each parameter
! if the "decimalPrecision" (aka, "DecScale") is set properly. - dmm
!     call grib_set(igrib, 'bitsPerValue',12,iret)
!     call LIS_verify(iret, 'grib_set:bitsPerValue failed in LIS_historyMod')

! Set the "decimalPrecision" (aka, "DecScale") based on the
! gribSF (grib scale factor) set in the MODEL OUTPUT TBL. - dmm
     gribSFtemp = gribSF
     decimalPrecision = 0
     do while (gribSFtemp.ge.10)
        decimalPrecision = decimalPrecision + 1
        gribSFtemp = gribSFtemp / 10
     enddo
     call grib_set(igrib, 'decimalPrecision',decimalPrecision,iret)
     call LIS_verify(iret, 'grib_set:decimalPrecision failed in LIS_historyMod')

! Section 6: Bit-Map

     call grib_set(igrib, 'bitmapPresent',1,iret)
     call LIS_verify(iret, 'grib_set:bitmapPresent failed in LIS_historyMod')
     
     call grib_set(igrib,'values',gtmp1,iret)
     call LIS_verify(iret, 'grib_set:values failed in LIS_historyMod')
     
     call grib_write(igrib,ftn,iret)
     call LIS_verify(iret, 'grib_write failed in LIS_historyMod')
     
     call grib_release(igrib,iret)
     call LIS_verify(iret,'grib_release failed in LIS_historyMod')

     call write_stats(gtmp, LIS_rc%gnc(n)*LIS_rc%gnr(n), &
          mvar, ftn_stats, form)
     
     deallocate(gtmp)
     deallocate(gtmp1)
  endif
#endif
          
end subroutine writevar_grib2_withstats_real

#if 0 
!BOP
! !ROUTINE: grid2tile_ens
! \label{grid2tile_ens}
!
! !INTERFACE:
  subroutine grid2tile_ens(n,m,gvar,tvar)
! !USES:

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n
    integer, intent(in) :: m
    real                :: gvar(LIS_rc%lnc(n),LIS_rc%lnr(n))
    real                :: tvar(LIS_rc%ntiles(n))

! !DESCRIPTION:
!  This routine converts a tile space variable to the corresponding
!  grid space. The aggregation involves weighted average of each tile
!  in a grid cell based on the vegetation distribution. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [tvar]
!     variable dimensioned in the tile space. 
!   \item [gvar]
!     variable after converstion to the grid space
!  \end{description}
!
!EOP
    integer           :: i,t,c,r

    do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
       r = LIS_domain(n)%tile(i)%row
       c = LIS_domain(n)%tile(i)%col
       t = i+m-1
       tvar(t) = gvar(c,r)
    enddo

  end subroutine grid2tile_ens
#endif
!BOP
! !ROUTINE: grid2tile
! \label{grid2tile}
!
! !INTERFACE:
!
! !INTERFACE:
  subroutine grid2tile_ens(n,m,gvar,tvar)
! !USES:

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n
    integer, intent(in) :: m
    real                :: gvar(LIS_rc%lnc(n),LIS_rc%lnr(n))
    real                :: tvar(LIS_rc%ntiles(n))

! !DESCRIPTION:
!  This routine converts a tile space variable to the corresponding
!  grid space. The aggregation involves weighted average of each tile
!  in a grid cell based on the vegetation distribution. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [tvar]
!     variable dimensioned in the tile space. 
!   \item [gvar]
!     variable after converstion to the grid space
!  \end{description}
!
!EOP
    integer           :: i,t,c,r

    do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
       r = LIS_domain(n)%tile(i)%row
       c = LIS_domain(n)%tile(i)%col
       t = i+m-1
       tvar(t) = gvar(c,r)
    enddo

  end subroutine grid2tile_ens

!BOP
! !ROUTINE: grid2tile
! \label{grid2tile}
!
! !INTERFACE:
  subroutine grid2tile(n,gvar,tvar)
! !USES:

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n
    real                :: tvar(LIS_rc%ntiles(n))
    real                :: gvar(LIS_rc%lnc(n),LIS_rc%lnr(n))
! !DESCRIPTION:
!  This routine converts a tile space variable to the corresponding
!  grid space. The aggregation involves weighted average of each tile
!  in a grid cell based on the vegetation distribution. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [tvar]
!     variable dimensioned in the tile space. 
!   \item [gvar]
!     variable after converstion to the grid space
!  \end{description}
!
!EOP
    integer           :: t,c,r

    do t=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(t)%row
       c = LIS_domain(n)%tile(t)%col
       tvar(t) = gvar(c,r)
    enddo
  end subroutine grid2tile


!BOP
! !ROUTINE: tile2grid_local
! \label{tile2grid_local}
!
! !INTERFACE:
  subroutine tile2grid_local(n,gvar,tvar)
! !USES:

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n
    real              :: gvar(LIS_rc%lnc(n),LIS_rc%lnr(n))
    real, intent(in)  :: tvar(LIS_rc%ntiles(n))
! !DESCRIPTION:
!  This routine converts a tile space variable to the corresponding
!  grid space. The aggregation involves weighted average of each tile
!  in a grid cell based on the vegetation distribution. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [tvar]
!     variable dimensioned in the tile space. 
!   \item [gvar]
!     variable after converstion to the grid space
!  \end{description}
!
!EOP
    integer           :: i,c,r,m,t

    gvar = 0.0
    do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
       c = LIS_domain(n)%tile(i)%col
       r = LIS_domain(n)%tile(i)%row
       do m=1,LIS_rc%nensem(n)
          t = i+m-1
          gvar(c,r) = gvar(c,r)+&
               tvar(t)*LIS_domain(n)%tile(t)%fgrd*&
               LIS_domain(n)%tile(t)%pens
          
       enddo
    enddo
    
  end subroutine tile2grid_local

!BOP
! !ROUTINE: tile2grid_local_ens
! \label{tile2grid_local_ens}
!
! !INTERFACE:
  subroutine tile2grid_local_ens(n,m,gvar,tvar)
! !USES:

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n
    integer, intent(in) :: m
    real              :: gvar(LIS_rc%lnc(n),LIS_rc%lnr(n))
    real, intent(in)  :: tvar(LIS_rc%ntiles(n))
! !DESCRIPTION:
!  This routine converts a tile space variable to the corresponding
!  grid space. The aggregation involves weighted average of each tile
!  in a grid cell based on the vegetation distribution. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [tvar]
!     variable dimensioned in the tile space. 
!   \item [gvar]
!     variable after converstion to the grid space
!  \end{description}
!
!EOP
    integer           :: i,c,r,t

    gvar = 0.0
    do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
       c = LIS_domain(n)%tile(i)%col
       r = LIS_domain(n)%tile(i)%row
       t = i+m-1
       gvar(c,r) = tvar(t)          
    enddo
    
  end subroutine tile2grid_local_ens

  subroutine tile2grid_global_ens(n,ensid,gvar,gvar_tile,global)

    implicit none
    
    integer, intent(in) :: n 
    integer, intent(in) :: ensid
    real                :: gvar(LIS_rc%gnc(n), LIS_rc%gnr(n))
    real                :: gvar_tile(LIS_rc%glbntiles_red(n))
    integer             :: global
    
    integer             :: count1
    integer             :: l,r,c,gid,stid,t, tid
    integer             :: ierr

    if(LIS_masterproc) then      
       do r=1,LIS_rc%gnr(n)
          do c=1,LIS_rc%gnc(n)
             gid = c+(r-1)*LIS_rc%gnc(n)
             stid = LIS_domain(n)%str_tind(gid)
             if(LIS_domain(n)%ntiles_pergrid(gid).gt.0) then 
                tid = stid + ensid-1
                gvar(c,r) = gvar_tile(tid) 
             endif
          enddo
       enddo
    endif

  end subroutine tile2grid_global_ens

  subroutine tile2grid_global_noens(n,gvar,gvar_tile,global)

    implicit none
    
    integer, intent(in) :: n 
    real                :: gvar(LIS_rc%gnc(n), LIS_rc%gnr(n))
    real                :: gvar_tile(LIS_rc%glbntiles_red(n))
    integer             :: global
    
    integer             :: count1
    integer             :: l,m,r,c,gid,stid,t, tid
    integer             :: ierr

    gvar = 0.0 
    if(LIS_masterproc) then      
       do r=1,LIS_rc%gnr(n)
          do c=1,LIS_rc%gnc(n)
             gid = c+(r-1)*LIS_rc%gnc(n)
             stid = LIS_domain(n)%str_tind(gid)
             if(LIS_domain(n)%ntiles_pergrid(gid).gt.0) then 
                do m=1,LIS_rc%nensem(n)
                   tid = stid + m-1
                   gvar(c,r) = gvar(c,r) + gvar_tile(tid) 
                enddo
                gvar(c,r) = gvar(c,r)/LIS_rc%nensem(n)
             endif
          enddo
       enddo
    endif

  end subroutine tile2grid_global_noens


!BOP
! !ROUTINE: LIS_patch2tile
! \label{LIS_patch2tile}
!
! !INTERFACE:
  subroutine LIS_patch2tile(n,m,tvar,pvar)
! !USES:

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n
    integer, intent(in) :: m
    real                :: pvar(LIS_rc%npatch(n,m))
    real                :: tvar(LIS_rc%ntiles(n))
! !DESCRIPTION:
!  This routine converts a tile space variable to the corresponding
!  grid space. The aggregation involves weighted average of each tile
!  in a grid cell based on the vegetation distribution. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [tvar]
!     variable dimensioned in the tile space. 
!   \item [gvar]
!     variable after converstion to the grid space
!  \end{description}
!
!EOP
    integer           :: i,t,tid

    do i=1,LIS_rc%npatch(n,m)
       tid = LIS_surface(n,m)%tile(i)%tile_id
       tvar(tid) = pvar(i) 
    enddo
    
  end subroutine LIS_patch2tile

  subroutine grid2patch_local(n,m,gvar,tvar)

    implicit none
    
    integer, intent(in) :: n 
    integer, intent(in) :: m
    real                :: gvar(LIS_rc%lnc(n), LIS_rc%lnr(n))
    real                :: tvar(LIS_rc%npatch(n,m))
    integer             :: t,r,c

    do t=1,LIS_rc%npatch(n,m)
       r = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
       c = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
       tvar(t) = gvar(c,r)
    enddo

  end subroutine grid2patch_local

  subroutine grid2patch_global(n,m,gvar,tvar,dummy)

    implicit none
    
    integer, intent(in) :: n 
    integer, intent(in) :: m
    real                :: gvar(LIS_rc%gnc(n), LIS_rc%gnr(n))
    real                :: tvar(LIS_rc%npatch(n,m))
    real,     allocatable   :: gvar_patch(:)
    
    integer             :: count1
    integer             :: l,r,c,gid,stid,t,npatch, tid
    logical             :: dummy
    integer             :: ierr

!    if(LIS_masterproc) then        
       allocate(gvar_patch(LIS_rc%glbnpatch(n,m)))
!    else
!       allocate(gvar_patch(1))
!    endif

    if(LIS_masterproc) then      
       do r=1,LIS_rc%gnr(n)
          do c=1,LIS_rc%gnc(n)
             gid = c+(r-1)*LIS_rc%gnc(n)             
             npatch = LIS_surface(n,m)%npatch_pergrid(gid)
             stid = LIS_surface(n,m)%str_patch_ind(gid)
             do t=1,npatch
                tid = stid + t-1
                gvar_patch(tid) = gvar(c,r)
             enddo
          enddo
       enddo
    endif
!scatterv will not work since the individual patch space is ordered
!differently when subset from the global array. If we use scatterv
!we need to reorder the arrays.So for now, using the BCAST
!#if (defined SPMD)   
!    call MPI_SCATTERV(gvar_patch, LIS_tdeltas(n,:),&
!         LIS_toffsets(n,:),MPI_REAL,tvar, &
!         LIS_tdeltas(n,LIS_localPet),MPI_REAL,0,LIS_mpi_comm,ierr)
!#else 
!    tvar = gvar_patch
!#endif

#if (defined SPMD)
    call MPI_BCAST(gvar_patch, LIS_rc%glbnpatch(n,m),MPI_REAL, &
         0, LIS_mpi_comm,ierr)
#endif

    count1=1
    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          gid = c+(r-1)*LIS_rc%gnc(n)
          npatch = LIS_surface(n,m)%npatch_pergrid(gid)
          stid = LIS_surface(n,m)%str_patch_ind(gid)
          do t=1,npatch
             tid = stid + t-1
             tvar(count1) = gvar_patch(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo

    deallocate(gvar_patch)

  end subroutine grid2patch_global


  subroutine grid2patch_global_ens(n,m,gvar,tvar,ensmode)

    implicit none
    
    integer, intent(in) :: n 
    integer, intent(in) :: m
    real                :: gvar(LIS_rc%gnc(n), LIS_rc%gnr(n),LIS_rc%nensem(n))
    real                :: tvar(LIS_rc%npatch(n,m))
    real,     allocatable   :: gvar_patch(:)
    integer             :: ensmode
    
    integer             :: count1
    integer             :: l,r,c,gid,stid,t,npatch, tid,kk
    integer             :: ierr

!    if(LIS_masterproc) then        
       allocate(gvar_patch(LIS_rc%glbnpatch(n,m)))
!    else
!       allocate(gvar_patch(1))
!    endif
    if(LIS_masterproc) then      
       do r=1,LIS_rc%gnr(n)
          do c=1,LIS_rc%gnc(n)
!             gid = c+(r-1)*LIS_rc%gnc(n)
!             npatch = LIS_surface(n,m)%npatch_pergrid(gid)
!             stid = LIS_surface(n,m)%str_patch_ind(gid)
!             do t=1,LIS_rc%nensem(n)
!                tid = stid + t-1
!                gvar_patch(tid) = gvar(c,r,t)
!             enddo
             gid = c+(r-1)*LIS_rc%gnc(n)
             npatch = LIS_surface(n,m)%npatch_pergrid(gid)
             stid = LIS_surface(n,m)%str_patch_ind(gid)
             do t=1,npatch/LIS_rc%nensem(n)
                do kk=1,LIS_rc%nensem(n)
                   tid = stid + (t-1)*LIS_rc%nensem(n)+kk-1
                   gvar_patch(tid) = gvar(c,r,kk)
                enddo
             enddo
          enddo
       enddo
    endif

!scatterv will not work since the individual tile space is ordered
!differently when subset from the global array. If we use scatterv
!we need to reorder the arrays.So for now, using the BCAST
!#if (defined SPMD)   
!    call MPI_SCATTERV(gvar_tile, LIS_tdeltas(n,:),&
!         LIS_toffsets(n,:),MPI_REAL,tvar, &
!         LIS_tdeltas(n,LIS_localPet),MPI_REAL,0,LIS_mpi_comm,ierr)
!#else 
!    tvar = gvar_tile
!#endif

#if (defined SPMD)
    call MPI_BCAST(gvar_patch, LIS_rc%glbnpatch(n,m),MPI_REAL, &
         0, LIS_mpi_comm,ierr)
#endif

    count1=1
    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          gid = c+(r-1)*LIS_rc%gnc(n)
          npatch = LIS_surface(n,m)%npatch_pergrid(gid)
          stid = LIS_surface(n,m)%str_patch_ind(gid)
          do t=1,npatch
             tid = stid + t-1
             tvar(count1) = gvar_patch(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo

    deallocate(gvar_patch)

  end subroutine grid2patch_global_ens


!BOP
! !ROUTINE: LIS_writevar_spread
! \label{LIS_writevar_spread}
!
! !INTERFACE:
  subroutine LIS_writevar_spread(ftn, n, m, varid, var,v)
! !USES:
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: varid
    real, intent(in)    :: var(LIS_rc%npatch(n,LIS_rc%lsm_index))
    integer, intent(in) :: v
!
! !DESCRIPTION:
!  Writes the innovations to a binary file in a gridded format. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the observation space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)

    integer :: ierr
    integer :: npatch
    integer :: stid, tid
    integer :: patch_deltas
    real    :: mean_v, std_v, var_v, max_v, min_v
    integer :: i,c,r,l,t,ntiles,gid,count1

    if(LIS_rc%wopt.eq."2d gridspace") then !gridded output 
       if(LIS_masterproc) then 
          allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(gtmp1(LIS_rc%glbnpatch(n,m)))
       else
          allocate(gtmp1(1))
       endif
#if (defined SPMD)      
       patch_deltas = LIS_patch_deltas(n,m,LIS_localPet)
       call MPI_GATHERV(var,patch_deltas,&
            MPI_REAL,gtmp1,LIS_patch_deltas(n,m,:),LIS_patch_offsets(n,m,:),&
            MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
       gtmp1 = var
#endif

       if(LIS_masterproc) then
          count1=1
          do l=1,LIS_npes
             do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                   gid = c+(r-1)*LIS_rc%gnc(n)
                   npatch = LIS_surface(n,m)%npatch_pergrid(gid)
                   stid = LIS_surface(n,m)%str_patch_ind(gid)
                   if(r.ge.LIS_nss_ind(n,l).and.&
                        r.le.LIS_nse_ind(n,l).and.&
                        c.ge.LIS_ews_ind(n,l).and.&
                        c.le.LIS_ewe_ind(n,l))then !points not in halo

#if 0 
                      mean_v = 0 
                      var_v  = 0 

                      do t=1,npatch
                         mean_v = mean_v + gtmp1(count1)
                         var_v  = var_v + gtmp1(count1)*gtmp1(count1)
                         count1 = count1 + 1
                      enddo
                      if(npatch.gt.0) then 
                         mean_v = mean_v/npatch
                         std_v  = (var_v/npatch - mean_v*mean_v)
                         if(std_v.ge.0) then 
                            std_v = sqrt(std_v)
                         else
                            std_v = LIS_rc%udef
                         endif
                      else
                         std_v = LIS_rc%udef
                      endif
                      gtmp(c,r) = std_v
#endif
!#if 0 

                      min_v = 1000000
                      max_v = -1000000 
                      do t=1,npatch
                         if(gtmp1(count1).lt.min_v) min_v = gtmp1(count1)
                         if(gtmp1(count1).gt.max_v) max_v = gtmp1(count1)
                         count1 = count1 + 1
                      enddo                      
                      gtmp(c,r) = max_v - min_v
                      if(min_v.eq. 1000000.or.&
                           max_v.eq.-1000000 ) then 
                         gtmp(c,r) = LIS_rc%udef
                      endif
!#endif                      
                   else
                      count1 = count1 + npatch
                   endif
                enddo
             enddo
          enddo
          
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          ierr = nf90_put_var(ftn,varid,gtmp,(/1,1,v/),&
               (/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
#endif
       end if
    else
       write(LIS_logunit,*) '[ERR] Writing ensemble spread output is supported '
       write(LIS_logunit,*) '[ERR] only in "2d gridspace"'
       call LIS_endrun()
    endif
  end subroutine LIS_writevar_spread


!BOP
! !ROUTINE: LIS_writevar_incr
! \label{LIS_writevar_incr}
!
! !INTERFACE:
  subroutine LIS_writevar_incr(ftn, n, m, varid, var,v)
! !USES:
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: varid
    real, intent(in)    :: var(LIS_rc%npatch(n,LIS_rc%lsm_index))
    integer, intent(in) :: v
!
! !DESCRIPTION:
!  Writes the analysis increments to a gridded output file
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the observation space
!  \end{description}
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)

    integer :: ierr
    integer :: npatch
    integer :: stid, tid
    integer :: patch_deltas
    real    :: mean_v, std_v, var_avg
    integer :: i,c,r,l,t,ntiles,gid,count1

    if(LIS_rc%wopt.eq."2d gridspace") then !gridded output 
       if(LIS_masterproc) then 
          allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(gtmp1(LIS_rc%glbnpatch(n,m)))
       else
          allocate(gtmp1(1))
       endif
#if (defined SPMD)      
       patch_deltas = LIS_patch_deltas(n,m,LIS_localPet)
       call MPI_GATHERV(var,patch_deltas,&
            MPI_REAL,gtmp1,LIS_patch_deltas(n,m,:),LIS_patch_offsets(n,m,:),&
            MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
       gtmp1 = var
#endif

       if(LIS_masterproc) then
          count1=1
          do l=1,LIS_npes
             do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                   gid = c+(r-1)*LIS_rc%gnc(n)
                   npatch = LIS_surface(n,m)%npatch_pergrid(gid)
                   stid = LIS_surface(n,m)%str_patch_ind(gid)
                   if(r.ge.LIS_nss_ind(n,l).and.&
                        r.le.LIS_nse_ind(n,l).and.&
                        c.ge.LIS_ews_ind(n,l).and.&
                        c.le.LIS_ewe_ind(n,l))then !points not in halo

                      var_avg = 0.0
                      
                      do t=1,npatch
                         var_avg = var_avg + gtmp1(count1)
                         count1 = count1 + 1                         
                      enddo                      
                      if(npatch.gt.0) then 
                         gtmp(c,r) = var_avg/npatch
                      else
                         gtmp(c,r) = LIS_rc%udef
                      endif
                   else
                      count1 = count1 + npatch
                   endif
                enddo
             enddo
          enddo
          
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          ierr = nf90_put_var(ftn,varid,gtmp,(/1,1,v/),&
               (/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
#endif
       end if
    else
       write(LIS_logunit,*) '[ERR] writing ensemble incr output is supported '
       write(LIS_logunit,*) '[ERR] only in "2d gridspace"'
       call LIS_endrun()
    endif
  end subroutine LIS_writevar_incr

!BOP
!
! !ROUTINE: stats
! \label{stats}
!
! !DESCRIPTION:
!  Calculates some diagnostic statistics for a given variable for model
!  output. The routine provides statistics such as mean, standard 
!  deviation, minimum and maximum values of a given variable. 
!
! !REVISION HISTORY:
!  Nov 11 1999:  Jon Radakovich; Initial code
! 
! !INTERFACE:
  subroutine stats(var,udef,ntiles,mean,stdev,min,max)
! !ARGUMENTS: 
    integer, intent(in) :: ntiles
    real, intent(in)    :: var(ntiles), udef
    real, intent(out)   :: mean,stdev,min,max
!EOP
    integer :: t, count1
    real :: dev, vsum

    vsum=0.
    mean=0.
    dev=0.
    stdev=0.
    min=100000.
    max=-100000.
    count1 = 0 
    do t=1,ntiles
       if(var(t).ne.udef)then
          count1 = count1 +1
          vsum=vsum+var(t)
          if(var(t).gt.max)max=var(t)
          if(var(t).lt.min)min=var(t)
       endif
    enddo
    if(vsum.eq.0.)then
       max=0.
       min=0.
    endif
    if(count1 .ge.1) then 
       mean=vsum/float(count1)
    else
       mean = 0 
    endif
    count1 = 0 
    do t=1,ntiles
       if(var(t).ne.udef)then
          count1 = count1 + 1
          dev=dev+(var(t)-mean)**2
       endif
    enddo
    if(count1 .gt.1) then 
       stdev=(dev*(float(count1)-1)**(-1))**(0.5)
    else
       stdev = 0
    endif
    return
  end subroutine stats

!BOP
!
! !ROUTINE: stats_ens
! \label{stats_ens}
!
! !DESCRIPTION:
!  Calculates some diagnostic statistics for a given variable for model
!  output. The routine provides statistics such as mean, standard 
!  deviation, minimum and maximum values of a given variable. 
!
! !REVISION HISTORY:
!  Nov 11 1999:  Jon Radakovich; Initial code
! 
! !INTERFACE:
  subroutine stats_ens(var,udef,nc,nr,nensem,mean,stdev,min,max)
! !ARGUMENTS: 
    integer, intent(in) :: nc,nr,nensem
    real, intent(in)    :: var(nc,nr,nensem), udef
    real, intent(out)   :: mean,stdev,min,max
!EOP
    integer :: c,r,m, count1
    real :: dev, vsum

    vsum=0.
    mean=0.
    dev=0.
    stdev=0.
    min=100000.
    max=-100000.
    count1 = 0 
    do m=1,nensem
       do r=1,nr
          do c=1,nc
             if(var(c,r,m).ne.udef)then
                count1 = count1 +1
                vsum=vsum+var(c,r,m)
                if(var(c,r,m).gt.max)max=var(c,r,m)
                if(var(c,r,m).lt.min)min=var(c,r,m)
             endif
          enddo
       enddo
    enddo

    if(vsum.eq.0.)then
       max=0.
       min=0.
    endif
    if(count1 .ge.1) then 
       mean=vsum/float(count1)
    else
       mean = 0 
    endif
    count1 = 0 
    do m=1,nensem
       do r=1,nr
          do c=1,nc
             if(var(c,r,m).ne.udef)then
                count1 = count1 + 1
                dev=dev+(var(c,r,m)-mean)**2
             endif
          enddo
       enddo
    enddo
    if(count1 .gt.1) then 
       stdev=(dev*(float(count1)-1)**(-1))**(0.5)
    else
       stdev = 0
    endif
    return
  end subroutine stats_ens
  
!BOP
!
! !ROUTINE: stats_da
! \label{stats_da}
!
! !DESCRIPTION:
!  Calculates some diagnostic statistics for a given variable for 
!  data assimilation output. The routine provides statistics 
!  such as mean, standard deviation, minimum and maximum values 
!  and spread of a given variable. 
!
! !REVISION HISTORY:
!  aug 22 2006:  Sujay Kumar; Initial code
! 
! !INTERFACE:
  subroutine stats_da(var,udef,ntiles,mean,spread, stdev,min,max)
! !ARGUMENTS: 
    integer, intent(in) :: ntiles
    real, intent(in)    :: var(ntiles), udef
    real, intent(out)   :: mean,stdev,min,max,spread
!EOP
    integer :: t, count1
    real :: dev, vsum

    vsum=0.
    mean=0.
    dev=0.
    stdev=0.
    min=100000.
    max=-100000.
    count1 = 0 
    do t=1,ntiles
       if(var(t).ne.udef)then
          count1 = count1 +1
          vsum=vsum+var(t)
          if(var(t).gt.max)max=var(t)
          if(var(t).lt.min)min=var(t)
       endif
    enddo
    if(vsum.eq.0.)then
       max=0.
       min=0.
    endif
    if(count1 .ge.1) then 
       mean=vsum/float(count1)
    else
       mean = 0 
    endif
    count1 = 0 
    do t=1,ntiles
       if(var(t).ne.udef)then
          count1 = count1 + 1
          dev=dev+(var(t)-mean)**2
       endif
    enddo
    if(count1 .gt.1) then 
       stdev=(dev*(float(count1)-1)**(-1))**(0.5)
    else
       stdev = 0
    endif
    !  spread = max-min
    spread = dev
    return
  end subroutine stats_da

!BOP
!
! !ROUTINE: write_stats
! \label{write_stats}
!
! !REVISION HISTORY:
!  29 Oct 2008:  James Geiger; Initial code
! 
! !INTERFACE: 
subroutine write_stats(var, size, mvar, ftn_stats, form)
! !USES: 

! !ARGUMENTS: 

   implicit none

   real, dimension(*), intent(in) :: var
   integer, intent(in) :: size
   character (len=*)   :: mvar
   integer, intent(in) :: ftn_stats
   integer, intent(in) :: form
!
! !DESCRIPTION:
!  Top level call to write statistics to the STATS file.
!  The arguments are:
!  \begin{description}
!   \item [var]
!     output data to process
!   \item [size]
!     size of array var
!   \item[mvar]
!     name of the variable being written
!   \item [ftn\_stats]
!     unit number of the ASCII text statistics file
!   \item[form]
!     format to be used in the stats file (1-decimal format, 
!     2-scientific format)
!  \end{description}
!EOP


   real :: vmean,vstdev,vmin,vmax

   if ( LIS_rc%sout ) then
      call stats(var, LIS_rc%udef, size, vmean, vstdev, vmin, vmax)

      if ( form == 1 ) then 
         write(ftn_stats,999) mvar,vmean,vstdev,vmin,vmax
      elseif ( form == 2 ) then 
         write(ftn_stats,998) mvar,vmean,vstdev,vmin,vmax
      endif
      flush(ftn_stats)
   endif

998 FORMAT(1X,A18,4E14.3)
999 FORMAT(1X,A18,4F14.3)

end subroutine write_stats

!BOP
!
! !ROUTINE: LIS_gather_tiled_vector_output
! \label{LIS_gather_tiled_vector_output}
!
! !REVISION HISTORY:
!  30 Jan 2009:  Sujay Kumar; Initial code
! 
! !INTERFACE:
subroutine LIS_gather_tiled_vector_output(n, gtmp, var)
! !USES: 

! !ARGUMENTS: 

   implicit none

   integer                       :: n
   real, allocatable             :: gtmp(:)
   real, intent(in)              :: var(LIS_rc%ntiles(n))

! !DESCRIPTION:
! This routine gathers the output data into a tiled 1d array.
!
! This process aggregates the variable into a tile space. 
!
! This process accounts for the halo.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [gtmp]
!     return array for the tiled output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

    real, allocatable :: gtmp1(:)
    integer :: count1 ,c,r,ntiles,t,gid,stid,tid,l
    integer :: tdeltas
    integer :: ierr
    
    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%glbntiles_red(n)))
       allocate(gtmp1(LIS_rc%glbntiles(n)))
    else
       allocate(gtmp(1))
       allocate(gtmp1(1))
    endif
#if (defined SPMD)      
    tdeltas = LIS_tdeltas(n,LIS_localPet)
    call MPI_GATHERV(var,tdeltas,&
         MPI_REAL,gtmp1,LIS_tdeltas(n,:),LIS_toffsets(n,:),MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
    gtmp1 = var
#endif
    if(LIS_masterproc) then 
       count1=1
       do l=1,LIS_npes
          do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                gid = c+(r-1)*LIS_rc%gnc(n)
                ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                stid = LIS_domain(n)%str_tind(gid)
                if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                   do t=1,ntiles
                      tid = stid + t-1
                      gtmp(tid) = gtmp1(count1)
                      count1 = count1 + 1
                   enddo
                else
                   count1 = count1 + ntiles
                endif
             enddo
          enddo
       enddo
    endif
    deallocate(gtmp1)
 end subroutine LIS_gather_tiled_vector_output


!BOP
!
! !ROUTINE: LIS_gather_tiled_vector_withhalo_output
! \label{LIS_gather_tiled_vector_withhalo_output}
!
! !REVISION HISTORY:
!  30 Jan 2009:  Sujay Kumar; Initial code
! 
! !INTERFACE:
subroutine LIS_gather_tiled_vector_withhalo_output(n, gtmp, var)
! !USES: 

! !ARGUMENTS: 

   implicit none

   integer                       :: n
   real, allocatable                :: gtmp(:)
   real, intent(in)             :: var(LIS_rc%ntiles(n))
! !DESCRIPTION:
! This routine gathers the output data into a tiled 1d array.
!
! This process aggregates the variable into a tile space. 
!
! This process accounts for the halo.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [gtmp]
!     return array for the tiled output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

    real, allocatable :: gtmp1(:)
    integer :: count1 ,c,r,ntiles,t,gid,stid,tid,l
    integer :: tdeltas
    integer :: ierr
    
    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%glbntiles(n)))
       allocate(gtmp1(LIS_rc%glbntiles(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)
    tdeltas = LIS_tdeltas(n,LIS_localPet)
    call MPI_GATHERV(var,tdeltas,&
         MPI_REAL,gtmp1,LIS_tdeltas(n,:),LIS_toffsets(n,:),MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
    gtmp1 = var
#endif
    if(LIS_masterproc) then 
       count1=1
       do l=1,LIS_npes
          do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                gid = c+(r-1)*LIS_rc%gnc(n)
                ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                stid = LIS_domain(n)%str_tind(gid)
!                if(r.ge.LIS_nss_ind(n,l).and.&
!                     r.le.LIS_nse_ind(n,l).and.&
!                     c.ge.LIS_ews_ind(n,l).and.&
!                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                do t=1,ntiles
                   tid = stid + t-1
                   gtmp(tid) = gtmp1(count1)
                   count1 = count1 + 1
                enddo
!                else
!                   count1 = count1 + ntiles
!                endif
             enddo
          enddo
       enddo
    endif
    deallocate(gtmp1)
  end subroutine LIS_gather_tiled_vector_withhalo_output

!BOP
!
! !ROUTINE: LIS_gather_patch_vector_output
! \label{LIS_gather_patch_vector_output}
!
! !REVISION HISTORY:
!  30 Jan 2009:  Sujay Kumar; Initial code
! 
! !INTERFACE:
subroutine LIS_gather_patch_vector_output(n, m, gtmp, var)
! !USES: 

! !ARGUMENTS: 
   implicit none

   integer                       :: n
   integer                       :: m
   real, allocatable             :: gtmp(:)
   real, intent(in)              :: var(LIS_rc%npatch(n,m))

! !DESCRIPTION:
! This routine gathers the output data into a patch 1d array.
!
! This process aggregates the variable into a tile space. 
!
! This process accounts for the halo.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [gtmp]
!     return array for the patch output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

    real, allocatable :: gtmp1(:)
    integer :: patch_deltas
    integer :: count1 ,c,r,npatch,t,gid,stid,tid,l
    integer :: ierr
    
    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%glbnpatch_red(n,m)))
       allocate(gtmp1(LIS_rc%glbnpatch(n,m)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)      
    patch_deltas = LIS_patch_deltas(n,m,LIS_localPet)
    call MPI_GATHERV(var,patch_deltas,&
         MPI_REAL,gtmp1,LIS_patch_deltas(n,m,:),LIS_patch_offsets(n,m,:),&
         MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
    gtmp1 = var
#endif
    if(LIS_masterproc) then 
       count1=1
       do l=1,LIS_npes
          do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                gid = c+(r-1)*LIS_rc%gnc(n)
                npatch = LIS_surface(n,m)%npatch_pergrid(gid)
                stid = LIS_surface(n,m)%str_patch_ind(gid)
                if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                   do t=1,npatch
                      tid = stid + t-1
                      gtmp(tid) = gtmp1(count1)
                      count1 = count1 + 1
                   enddo
                else
                   count1 = count1 + npatch
                endif
             enddo
          enddo
       enddo
    endif
    deallocate(gtmp1)
 end subroutine LIS_gather_patch_vector_output

!BOP
!
! !ROUTINE: gather_gridded_output_tile
! \label{gather_gridded_output_tile}
!
! !REVISION HISTORY:
!  29 Oct 2008:  James Geiger; Initial code
! 
! !INTERFACE:
subroutine gather_gridded_output_tile(n, gtmp, var)
! !USES: 

! !ARGUMENTS: 
   implicit none

   integer, intent(in)          :: n
   real, allocatable            :: gtmp(:,:)
   real, intent(in)             :: var(LIS_rc%ntiles(n))

! !DESCRIPTION:
! This routine gathers the output data into a gridded array.
!
! This process aggregates tiles to their grid, accounting for ensemble runs.
!
! This process accounts for the halo.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [gtmp]
!     return array for the gridded output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

   real, allocatable :: var1(:)
   real, allocatable :: gtmp1(:)
   integer :: i,c,r,m,t,l
   integer :: gdeltas
   integer :: count1,gid,ntiles,ierr

   allocate(var1(LIS_rc%ngrid(n)))
   if ( LIS_masterproc ) then 
      allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
      allocate(gtmp1(LIS_rc%glbngrid(n)))
      gtmp = 0.0
      gtmp1 = 0.0
   else
      allocate(gtmp1(1))
      gtmp1 = 0.0
   endif

   var1 = 0 

   do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
      c = LIS_domain(n)%tile(i)%index
      do m=1,LIS_rc%nensem(n)
         t = i+m-1
         if ( var(t) == -9999.0 ) then
            var1(c) = -9999.0
         else
            var1(c) = var1(c) + &
                 var(t)*LIS_domain(n)%tile(t)%fgrd*&
                 LIS_domain(n)%tile(t)%pens
         endif
      enddo
   enddo
#if (defined SPMD)      
   gdeltas = LIS_gdeltas(n,LIS_localPet)
   call MPI_GATHERV(var1,gdeltas,&
        MPI_REAL,gtmp1,LIS_gdeltas(n,:),LIS_goffsets(n,:),&
        MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
   gtmp1 = var1
#endif

   if ( LIS_masterproc ) then 
      gtmp = LIS_rc%udef
      count1=1

      do l=1,LIS_npes
         do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
            do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
               gid = c+(r-1)*LIS_rc%gnc(n)
               ntiles = LIS_domain(n)%ntiles_pergrid(gid)
               if(ntiles.ne.0) then                          
                  if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                     gtmp(c,r) = gtmp1(count1)
                  endif
                  count1 = count1 + 1
               endif
            enddo
         enddo
      enddo
   endif
   deallocate(gtmp1)
   deallocate(var1)

end subroutine gather_gridded_output_tile

!BOP
!
! !ROUTINE: gather_gridded_output_patch
! \label{gather_gridded_output_patch}
!
! !REVISION HISTORY:
!  29 Oct 2008:  James Geiger; Initial code
! 
! !INTERFACE:
subroutine gather_gridded_output_patch(n, mtype, gtmp, var)
! !USES: 

! !ARGUMENTS: 
   implicit none

   integer, intent(in)          :: n
   integer, intent(in)          :: mtype
   real, allocatable                :: gtmp(:,:)
   real, intent(in)             :: var(LIS_rc%npatch(n,mtype))


! !DESCRIPTION:
! This routine gathers the output data into a gridded array.
!
! This process aggregates tiles to their grid, accounting for ensemble runs.
!
! This process accounts for the halo.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [mtype]
!     index of the surface model type
!   \item [gtmp]
!     return array for the gridded output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

   real, allocatable :: var1(:)
   real, allocatable :: gtmp1(:)
   integer :: i,c,r,m,t,l
   integer :: gdeltas
   integer :: count1,gid,ntiles,ierr

   allocate(var1(LIS_rc%ngrid(n)))
   if ( LIS_masterproc ) then 
      allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
      allocate(gtmp1(LIS_rc%glbngrid(n)))
      gtmp = 0.0
      gtmp1 = 0.0
   else
      allocate(gtmp1(1))
      gtmp1 = 0.0
   endif

   var1 = 0 

   do i=1,LIS_rc%npatch(n,mtype),LIS_rc%nensem(n)
      c = LIS_surface(n,mtype)%tile(i)%index
      do m=1,LIS_rc%nensem(n)
         t = i+m-1
         if ( var(t) == -9999.0 ) then 
            var1(c) = -9999.0
         else
            var1(c) = var1(c) + &
            var(t)*LIS_surface(n,mtype)%tile(t)%fgrd*&
            LIS_surface(n,mtype)%tile(t)%pens
         endif
      enddo
   enddo

#if (defined SPMD)  
   gdeltas = LIS_gdeltas(n,LIS_localPet)
   call MPI_GATHERV(var1,gdeltas,&
        MPI_REAL,gtmp1,LIS_gdeltas(n,:),LIS_goffsets(n,:),MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
   gtmp1 = var1
#endif

   if ( LIS_masterproc ) then 
      gtmp = LIS_rc%udef
      count1=1

      do l=1,LIS_npes
         do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
            do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
               gid = c+(r-1)*LIS_rc%gnc(n)
               ntiles = LIS_domain(n)%ntiles_pergrid(gid)
               if(ntiles.ne.0) then                          
                  if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                     gtmp(c,r) = gtmp1(count1)
                  endif
                  count1 = count1 + 1
               endif
            enddo
         enddo
      enddo
   endif
   deallocate(gtmp1)
   deallocate(var1)
 end subroutine gather_gridded_output_patch


!BOP
!
! !ROUTINE: LIS_gather_gridded_vector_output
! \label{LIS_gather_gridded_vector_output}
!
! !REVISION HISTORY:
!  30 Jan 2009: Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine LIS_gather_gridded_vector_output(n, gtmp, var)
! !USES: 

! !ARGUMENTS: 
   implicit none

   integer, intent(in)          :: n
   real, allocatable            :: gtmp(:)
   real, intent(in)             :: var(LIS_rc%ntiles(n))

!
! !DESCRIPTION:
! This routine gathers the output data into a gridded 1d array.
!
! This process aggregates tiles to their grid, accounting for ensemble runs.
!
! This process accounts for the halo.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [gtmp]
!     return array for the gridded output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

   real, allocatable  :: gtmp2d(:,:)
   real, allocatable :: var1(:)
   integer :: i,c,r,m,t,l
   integer :: ntiles, gid, count1
   integer :: ierr
   integer :: gdeltas

   allocate(var1(LIS_rc%ngrid(n)))
   if ( LIS_masterproc ) then 
      allocate(gtmp(LIS_rc%glbngrid(n)))
      allocate(gtmp2d(LIS_rc%gnc(n), LIS_rc%gnr(n)))
      gtmp = 0.0     
   else
      allocate(gtmp(1))
      allocate(gtmp2d(1,1))
      gtmp = 0.0
   endif

   var1 = 0 

   do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
      c = LIS_domain(n)%tile(i)%index
      do m=1,LIS_rc%nensem(n)
         t = i+m-1
         if ( var(t) == -9999.0 ) then
            var1(c) = -9999.0
         else
            var1(c) = var1(c) + &
            var(t)*LIS_domain(n)%tile(t)%fgrd*&
            LIS_domain(n)%tile(t)%pens
         endif
      enddo
   enddo

#if (defined SPMD)      
   gdeltas = LIS_gdeltas(n,LIS_localPet)
   call MPI_GATHERV(var1,gdeltas,&
   MPI_REAL,gtmp,LIS_gdeltas(n,:),LIS_goffsets(n,:),MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
   gtmp= var1
#endif
! Though gtmp now has a size of glbngrid, the points are not ordered correctly. 
   if ( LIS_masterproc ) then 
      gtmp2d = LIS_rc%udef
      count1=1

      do l=1,LIS_npes
         do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
            do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
               gid = c+(r-1)*LIS_rc%gnc(n)
               ntiles = LIS_domain(n)%ntiles_pergrid(gid)
               if(ntiles.ne.0) then                          
                  if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                     gtmp2d(c,r) = gtmp(count1)
                  endif
                  count1 = count1 + 1
               endif
            enddo
         enddo
      enddo
! now reorder the gtmp2d to create a 1d vector
      deallocate(gtmp)
      allocate(gtmp(LIS_rc%glbngrid_red(n)))
      count1 = 1
      do r=1,LIS_rc%gnr(n)
         do c=1,LIS_rc%gnc(n)
            if(LIS_domain(n)%ntiles_pergrid(c+(r-1)*LIS_rc%gnc(n)).gt.0) then 
               gtmp(count1) = gtmp2d(c,r)
               count1 = count1+1
            endif
         enddo
      enddo
      
   endif
   deallocate(gtmp2d)
   deallocate(var1)
  
 end subroutine LIS_gather_gridded_vector_output

!BOP
!
! !ROUTINE: LIS_gather_1dgrid_to_2dgrid
! \label{LIS_gather_1dgrid_to_2dgrid}
!
! !REVISION HISTORY:
!  30 Jan 2009: Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine LIS_gather_1dgrid_to_2dgrid(n, gtmp, var)
! !USES: 

! !ARGUMENTS: 
   implicit none

   integer, intent(in)          :: n
   real, allocatable            :: gtmp(:,:)
   real, intent(in)             :: var(LIS_rc%ngrid(n))

!
! !DESCRIPTION:
! This routine gathers a gridded 1d array into a 2d global array.
!
! This process accounts for the halo.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [gtmp]
!     return array for the gridded output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

   real, allocatable  :: gtmp1d(:)
   integer :: i,c,r,m,t,l
   integer :: ntiles, gid, count1
   integer :: ierr
   integer :: gdeltas

   if ( LIS_masterproc ) then 
      allocate(gtmp1d(LIS_rc%glbngrid(n)))
      allocate(gtmp(LIS_rc%gnc(n), LIS_rc%gnr(n)))
      gtmp = 0.0     
   else
      allocate(gtmp1d(1))
      allocate(gtmp(1,1))
      gtmp = 0.0
   endif

#if (defined SPMD)      
   gdeltas = LIS_gdeltas(n,LIS_localPet)
   call MPI_GATHERV(var,gdeltas,&
        MPI_REAL,gtmp1d,LIS_gdeltas(n,:),LIS_goffsets(n,:),&
        MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
   gtmp1d= var
#endif
! Though gtmp now has a size of glbngrid, the points are not ordered correctly. 
   if ( LIS_masterproc ) then 
      gtmp = LIS_rc%udef
      count1=1

      do l=1,LIS_npes
         do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
            do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
               gid = c+(r-1)*LIS_rc%gnc(n)
               ntiles = LIS_domain(n)%ntiles_pergrid(gid)
               if(ntiles.ne.0) then                          
                  if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                     gtmp(c,r) = gtmp1d(count1)
                  endif
                  count1 = count1 + 1
               endif
            enddo
         enddo
      enddo
   endif
   deallocate(gtmp1d)
  
 end subroutine LIS_gather_1dgrid_to_2dgrid


!BOP
!
! !ROUTINE: LIS_scatter_global_to_local_grid
! \label{LIS_scatter_global_to_local_grid}
!
! !REVISION HISTORY:
!  30 Jan 2009: Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine LIS_scatter_global_to_local_grid(n, gtmp,ltmp)
! !USES: 

! !ARGUMENTS: 
   implicit none

   integer, intent(in)          :: n
   real, intent(in)             :: gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n))
   real                         :: ltmp(LIS_rc%lnc(n),LIS_rc%lnr(n))

!
! !DESCRIPTION:
! This routine gathers a gridded 1d array into a 2d global array.
!
! This process accounts for the halo.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [gtmp]
!     return array for the gridded output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

   real, allocatable :: gtmp1d(:)
   real, allocatable :: gtmp2d(:,:)
   integer :: i,c,r,m,t,l
   integer :: gid
   integer :: ierr


   allocate(gtmp1d(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
   gtmp1d = LIS_rc%udef

   if ( LIS_masterproc ) then 
      do r=1,LIS_rc%gnr(n)
         do c=1,LIS_rc%gnc(n)
            gid = c+(r-1)*LIS_rc%gnc(n)
            gtmp1d(gid) = gtmp(c,r)
         enddo
      enddo
   endif
#if (defined SPMD)      
   call MPI_BCAST(gtmp1d, LIS_rc%gnc(n)*LIS_rc%gnr(n), MPI_REAL,0, &
        LIS_mpi_comm, ierr)
   call LIS_verify(ierr, 'MPI_BCAST failed in LIS_scatter_global_to_local_grid')
#endif
   allocate(gtmp2d(LIS_rc%gnc(n),LIS_rc%gnr(n)))
   do r=1,LIS_rc%gnr(n)
      do c=1,LIS_rc%gnc(n)
         gid = c+(r-1)*LIS_rc%gnc(n)
         gtmp2d(c,r) = gtmp1d(gid)
      enddo
   enddo
   deallocate(gtmp1d)
   ltmp(:,:) = gtmp2d(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1))
   deallocate(gtmp2d)

 end subroutine LIS_scatter_global_to_local_grid


!BOP
! !ROUTINE: LIS_convertVarToLocalSpace
! \label{LIS_convertVarToLocalSpace}
! 
! !INTERFACE:
! 
  subroutine LIS_convertVarToLocalSpace(n,gvar,lvar)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: n
    real                  :: gvar(LIS_rc%glbngrid_red(n))
    real                  :: lvar(LIS_rc%ngrid(n))
! !DESCRIPTION:
!  Reads a real variable from a binary 
!  sequential access, 1d gridded file. After reading
!  the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item [oned]
!     dummy variable to distinguish the interface. 
!  \end{description}
!EOP
    real, allocatable :: gtmp2d(:,:)
    integer :: i, gid
    integer :: c,r, nc, cnt
    integer :: c1,c2,r1,r2

    allocate(gtmp2d(LIS_rc%gnc(n), LIS_rc%gnr(n)))

    cnt = 1
    do r=1,LIS_rc%gnr(n)
       do c=1,LIS_rc%gnc(n)
          gid = c+(r-1)*LIS_rc%gnc(n)
          if(LIS_domain(n)%ntiles_pergrid(gid).gt.0) then 
             gtmp2d(c,r) = gvar(cnt)
             cnt = cnt+1
          endif
       enddo
    enddo

    nc = (LIS_ewe_halo_ind(n,LIS_localPet+1)-LIS_ews_halo_ind(n,LIS_localPet+1))+1

    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          c1 = c-LIS_ews_halo_ind(n,LIS_localPet+1)+1
          r1 = r-LIS_nss_halo_ind(n,LIS_localPet+1)+1
          gid = LIS_domain(n)%gindex(c1,r1)
          if(gid.ne.-1) then
             lvar(gid) = gtmp2d(c,r)
          endif
       enddo
    enddo

    deallocate(gtmp2d)

  end subroutine LIS_convertVarToLocalSpace


!BOP
! !ROUTINE: LIS_gather_2d_local_to_global
! \label{LIS_gather_2d_local_to_global}
!
! !INTERFACE:
subroutine LIS_gather_2d_local_to_global(n, lvar, gvar)
! !USES:

   implicit none
! !ARGUMENTS: 
   integer, intent(in)  :: n
   real,    intent(in)  :: lvar(LIS_rc%lnc_red(n), LIS_rc%lnr_red(n))
   real,    intent(out) :: gvar(LIS_rc%gnc(n), LIS_rc%gnr(n))

! !DESCRIPTION:
! This routine gathers local (sub-domain) gridded 2d arrays into a
! 2d global array.
!
! The local (sub-domain) arrays must *not* have a halo.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [lvar]
!     local (per-process) 2d array
!   \item [gvar]
!     global 2d array
!  \end{description}
!EOP

   real, allocatable :: var1(:)
   real, allocatable :: gtmp1(:)
   integer           :: count1,c,r,ierr,l
   integer           :: deltas

   allocate(var1(LIS_rc%lnc_red(n)*LIS_rc%lnr_red(n)))
   allocate(gtmp1(LIS_rc%gnc(n)*LIS_rc%gnr(n)))

   var1 = reshape(lvar, (/LIS_rc%lnc_red(n)*LIS_rc%lnr_red(n)/))

#if (defined SPMD)      
   deltas = LIS_deltas(n,LIS_localPet)
   call MPI_ALLGATHERV(var1,deltas,MPI_REAL,                             &
                       gtmp1,LIS_deltas(n,:),LIS_offsets(n,:), MPI_REAL, &
                       LIS_mpi_comm,ierr)
#else 
   gtmp1 = var1
#endif

   count1 = 1
   do l = 1,LIS_npes
      do r = LIS_nss_ind(n,l),LIS_nse_ind(n,l)
         do c = LIS_ews_ind(n,l),LIS_ewe_ind(n,l)
            gvar(c,r) = gtmp1(count1)
            count1 = count1 + 1
         enddo
      enddo
   enddo

   deallocate(gtmp1)
   deallocate(var1)
end subroutine LIS_gather_2d_local_to_global

end module LIS_historyMod
