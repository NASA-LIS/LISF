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
#include "LDT_NetCDF_inc.h"
module LDT_historyMod
!BOP
!
! !MODULE: LDT_historyMod
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
! \item{grib1: WMO grib1 format}
! \item{netcdf: NCAR unidata netcdf format}
! \end{itemize}
! This module also provides generic interfaces for writing and reading
! model restart files. LIS uses the following convention on the 
! filename extensions: 
! \begin{itemize}
! \item{.1gd4r: binary output, 1 correspondes to one variable, g represents
! gridded data, d represnts direct access, 4r represents 4 byte real}
! \item{.grb: files in grib1 format}
! \item{.nc: files in netcdf format}
! \item{*rst: restart files. Currently all restart files are written as
! binary format}
! \end{itemize}
! This module also manages any data aggregatiion from computing nodes
! when parallel processing is employed. 
!
! !REVISION HISTORY:
!  10 Feb 2004: Sujay Kumar; Initial Specification
!   4 Jul 2008: Sujay Kumar; Redesigned with generic routines that can handle
!                   model output in different formats for all LSMs
!  13 Jul 2008: Sujay Kumar; Added support for GRIB2
!  01 Sep 2015: KR Arsenault; Adding support for parallel computing in LDT
! 
! !USES: 
  use LDT_paramDataMod
  use LDT_coreMod
  use LDT_logMod

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

#if (defined SPMD)
#if (defined USE_INCLUDE_MPI)
  implicit none
  include 'mpif.h' 
#else
  use mpi
  implicit none
#endif
#endif

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_writevar_bin               ! write a variable into a binary file
  public :: LDT_writevar_restart           ! writes a variable into a restart file
  public :: LDT_readvar_restart            ! reads a variable from a restart file
  public :: LDT_writevar_reduced_tilespace ! writes a variable in reduced tilespace (ngrid)
  public :: LDT_readvar_reduced_tilespace  ! reads a reduced tilespace variable
  public :: LDT_readvar_gridded            ! read a variable from a gridded file
  public :: LDT_writevar_gridded           ! write a variable in grid space to a file
  public :: LDT_writevar_data_header       ! writes header information into the output file
  public :: LDT_close_data_header          ! closes entries for header information
  public :: LDT_tile2grid                  ! convert a tilespace variable to gridspace
  public :: LDT_writeNETCDFdataHeader
  public :: LDT_writeNETCDFdata
  public :: LDT_readLISSingleNetcdfVar
  public :: LDT_writevar_netcdf

!EOP

!BOP 
! 
! !ROUTINE: LDT_writevar_bin
! \label{LDT_writevar_bin}
! 
! !INTERFACE:
  interface LDT_writevar_bin
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure writevar_bin_real
     module procedure writevar_bin_real_direct
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
! !ROUTINE: LDT_writeNETCDFdata
! \label{LDT_writeNETCDFdata}
!
! !INTERFACE:
  interface LDT_writeNETCDFdata

! !PRIVATE MEMBER FUNCTIONS:
     module procedure writeNETCDFdata_default
     module procedure writeNETCDFdata_givendomain
     module procedure writeNETCDFdata_giventype
     module procedure writeNETCDFdata_LISHydro

!
! !DESCRIPTION:
! This interface provides routines for writing data variables 
! in a netcdf file. The aggregation from decomposed tile spaces in
! each processor to their respective grid spaces and the aggregations
! of these grid spaces from each processor are also performed by
! this routine.
!
!EOP
  end interface

!BOP 
! 
! !ROUTINE: LDT_writevar_restart
! \label{LDT_writevar_restart}
!
! !INTERFACE:
  interface LDT_writevar_restart
! !PRIVATE MEMBER FUNCTIONS:

     module procedure writevar_restart_int
     module procedure writevar_restart_real
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
! 
! !ROUTINE: LDT_writevar_reduced_tilespace
! \label{LDT_writevar_reduced_tilespace}
! 
! !INTERFACE: 
  interface LDT_writevar_reduced_tilespace
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
! !ROUTINE: LDT_readvar_reduced_tilespace
! \label{LDT_readvar_reduced_tilespace}
! 
! !INTERFACE: 
  interface LDT_readvar_reduced_tilespace
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
! !ROUTINE: LDT_readvar_restart
! \label{LDT_readvar_restart}
!
! !INTERFACE:
  interface LDT_readvar_restart
! !PRIVATE MEMBER FUNCTIONS:
     module procedure readvar_restart_int
     module procedure readvar_restart_real
!     module procedure readvar_restart_char

! !DESCRIPTION:
! This interface provides routines for reading variables (integer or real)
! from a restart file. The restart files are read by the master processor. 
! The domain decomposition of the "global" tile space on the master processor
! to individual processors are also performed by this routine. 
!EOP
  end interface

!BOP 
! 
! !ROUTINE: LDT_readvar_gridded
! \label{LDT_readvar_gridded}
!
! !INTERFACE:
  interface LDT_readvar_gridded
! !PRIVATE MEMBER FUNCTIONS:

     module procedure readvar_1dgridded_real
     module procedure readvar_2dgridded_real
     module procedure readvar_1dgridded_fromvector_real

! !DESCRIPTION:
! This interface provides routines for reading variables (real)
! from a gridded (output) file into a 1d or 2d local gridded space. 
! The gridded files are read by the master processor. 
! The domain decomposition of the "global" grid space on the master processor
! to individual processors are also performed by this routine. 
!EOP
  end interface

!BOP 
! 
! !ROUTINE: LDT_writevar_gridded
! \label{LDT_writevar_gridded}
!
! !INTERFACE:
  interface LDT_writevar_gridded
! !PRIVATE MEMBER FUNCTIONS:

     module procedure writevar_gridded_real
     module procedure writevar_gridded_real_2d
     module procedure writevar_gridded_real_3d
     module procedure writevar_gridded_real_4d
     module procedure writevar_gridded_integer_1d  !Y.Kwon

! !DESCRIPTION:
! This interface provides routines for writing variables (real)
! in grid space to a file. 
!EOP
  end interface


!BOP 
! 
! !ROUTINE: LDT_readLISSingleNetcdfVar
! \label{LDT_readLISSingleNetcdfVar}
! 
! !INTERFACE:
  interface LDT_readLISSingleNetcdfVar
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure readLISSingleNetcdfVar_LDTgrid
     module procedure readLISSingleNetcdfVar_LDTgrid_withunits
     module procedure readLISSingleNetcdfVar_Inputgrid
! 
! !DESCRIPTION:
!  This interface provides routines for reading a variable 
!  from a LIS netcdf file. The routines allow the specification
!  of an input grid or assumes that the field to be read
!  is in the LDT grid space. 
!
!EOP 
  end interface

!BOP 
! 
! !ROUTINE: LDT_writeNETCDFdataHeader
! \label{LDT_writeNETCDFdataHeader}
! 
! !INTERFACE:
  interface LDT_writeNETCDFdataHeader
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure writeNETCDFdataHeader_LIS
     module procedure writeNETCDFdataHeader_LISHydro
! 
! !DESCRIPTION:
!  This interface provides routines for reading a variable 
!  from a LIS netcdf file. The routines allow the specification
!  of an input grid or assumes that the field to be read
!  is in the LDT grid space. 
!
!EOP 
  end interface

contains

!BOP
! !ROUTINE: LDT_writevar_data_header
! \label{LDT_writevar_data_header}
!
! !INTERFACE: 
  subroutine LDT_writevar_data_header(n, ftn, name, varid, vlevels, entryno)
! !USES:     
    use LDT_coreMod,   only : LDT_rc
    use LDT_logMod,    only : LDT_verify

    implicit none
!
! !DESCRIPTION: 
!   This routine writes the global header information and the 
!   variable-specific header information into the output file.   
!
!EOP    
    integer,         intent(in) :: n
    integer,         intent(in) :: ftn
    character(len=*),intent(in) :: name
    integer                     :: varid
    integer                     :: vlevels
    integer                     :: entryno

    integer                     :: dimID(3)

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    if(LDT_rc%wout_form.eq.3) then 
       if(entryno.eq.1) then 
          call LDT_verify(nf90_def_dim(ftn,'east_west',LDT_rc%gnc(n),dimID(1)))
          call LDT_verify(nf90_def_dim(ftn,'north_south',LDT_rc%gnr(n),dimID(2)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value", -9999.0))          

          !grid information
          if(trim(LDT_rc%lis_map_proj(n)).eq."latlon") then !latlon
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "EQUIDISTANT CYLINDRICAL"))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LAT", &
                  LDT_rc%gridDesc(n,4)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LDT_rc%gridDesc(n,9)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LDT_rc%gridDesc(n,10)))       

          elseif(trim(LDT_rc%lis_map_proj(n)).eq."mercator") then 
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "MERCATOR"))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LAT", &
                  LDT_rc%gridDesc(n,4)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LON", &
                  LDT_rc%gridDesc(n,5)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LDT_rc%gridDesc(n,10)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LDT_rc%gridDesc(n,11)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LDT_rc%gridDesc(n,8)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LDT_rc%gridDesc(n,9)))

          elseif(trim(LDT_rc%lis_map_proj(n)).eq."lambert") then !lambert conformal
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "LAMBERT CONFORMAL"))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LAT", &
                  LDT_rc%gridDesc(n,4)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LON", &
                  LDT_rc%gridDesc(n,5)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LDT_rc%gridDesc(n,10)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
                  LDT_rc%gridDesc(n,7)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LDT_rc%gridDesc(n,11)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LDT_rc%gridDesc(n,8)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LDT_rc%gridDesc(n,9)))
             
          elseif(trim(LDT_rc%lis_map_proj(n)).eq."polar") then ! polar stereographic
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "POLAR STEREOGRAPHIC"))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LAT", &
                  LDT_rc%gridDesc(n,4)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LON", &
                  LDT_rc%gridDesc(n,5)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LDT_rc%gridDesc(n,10)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
                  LDT_rc%gridDesc(n,7)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LDT_rc%gridDesc(n,11)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LDT_rc%gridDesc(n,8)))
             call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LDT_rc%gridDesc(n,9)))
          endif
       endif

! The following are the 3-d fields
       if(trim(name).eq."SoilMoist") then 
          call LDT_verify(nf90_def_dim(ftn,"soilm_profiles", vlevels,dimID(3)))             
       endif
       if(trim(name).eq."SoilTemp") then 
          call LDT_verify(nf90_def_dim(ftn,"soilt_profiles", vlevels,dimID(3)))             
       endif
       if(vlevels.eq.1) then 
          call LDT_verify(nf90_def_var(ftn,trim(name),nf90_float,&
               dimids = dimID(1:2), varID=varId))
       else
          call LDT_verify(nf90_def_var(ftn,trim(name),nf90_float,&
               dimids = dimID, varID=varId))
       endif

    endif
#endif    

  end subroutine LDT_writevar_data_header

  subroutine LDT_close_data_header(ftn)
    use LDT_coreMod, only : LDT_rc
    use LDT_logMod, only : LDT_verify

    integer, intent(in) :: ftn 
    if(LDT_rc%wout_form.eq.3) then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)     
       call LDT_verify(nf90_enddef(ftn))
#endif
    endif
  end subroutine LDT_close_data_header

!BOP
! !ROUTINE: writevar_bin_real
! \label{writevar_bin_real}
!
! !INTERFACE:
! Private name: call using LDT_writevar_bin
  subroutine writevar_bin_real(n,ftn, var)
! !USES:
    use LDT_coreMod, only : LDT_rc, LDT_domain, LDT_localPet, LDT_masterproc

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in) :: ftn
    real                :: var(LDT_rc%ntiles(n))
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
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)

    if(LDT_rc%wopt.eq."1d tilespace") then     ! 1d tiled output
       call gather_tiled_vector_output(n,gtmp1,var)

       if ( LDT_masterproc ) then
          write(ftn) gtmp1
          deallocate(gtmp1)
       endif

    elseif(LDT_rc%wopt.eq."2d gridspace") then !2d gridded output
       call gather_gridded_output(n,gtmp, var)

       if ( LDT_masterproc ) then
          write(ftn) gtmp
          deallocate(gtmp)
       endif
    elseif(LDT_rc%wopt.eq."1d gridspace") then !1d gridded output
       call gather_gridded_vector_output(n,gtmp1, var)

       if ( LDT_masterproc ) then
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
! Private name: call using LDT_writevar_bin
  subroutine writevar_bin_real_direct(n,ftn, var, direct)

! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, &
         LDT_localPet, LDT_masterproc,LDT_gdeltas,LDT_goffsets,&
         LDT_nss_ind,LDT_ews_ind,LDT_nse_ind,LDT_ewe_ind,&
         LDT_nss_halo_ind,LDT_ews_halo_ind,LDT_nse_halo_ind,LDT_ewe_halo_ind,&
         LDT_npes

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in) :: ftn
    real,    intent(in) :: var(LDT_rc%lnc(n), LDT_rc%lnr(n))
    integer, intent(in) :: direct 
! !DESCRIPTION:
!  Write a real variable to a binary output file in a direct access mode
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!   \item [direct]
!     dummy argument for direct access output 
!  \end{description}
!EOP
    real, allocatable :: var1(:)
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)
    integer :: count1 ,c,r,gid,ntiles,ierr,l

    allocate(var1(LDT_rc%ngrid(n)))
    if(LDT_masterproc) then 
       allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
       allocate(gtmp1(LDT_rc%glbngrid(n)))
    else
       allocate(gtmp1(1))
    endif
    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)
          if(LDT_domain(n)%gindex(c,r).ne.-1) then 
             var1(LDT_domain(n)%gindex(c,r)) = var(c,r)
          endif
       enddo
    enddo
#if (defined SPMD)      
    call MPI_GATHERV(var1,LDT_gdeltas(n,LDT_localPet),&
         MPI_REAL,gtmp1,LDT_gdeltas(n,:),LDT_goffsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
    gtmp1 = var1
#endif
    if(LDT_masterproc) then 
       gtmp = LDT_rc%udef
       count1=1
       do l=1,LDT_npes
          do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
             do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                gid = c+(r-1)*LDT_rc%gnc(n)
                ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                if(ntiles.ne.0) then 
                   if(r.ge.LDT_nss_ind(n,l).and.&
                        r.le.LDT_nse_ind(n,l).and.&
                        c.ge.LDT_ews_ind(n,l).and.&
                        c.le.LDT_ewe_ind(n,l))then !points not in halo                      
                      gtmp(c,r) = gtmp1(count1)
                   endif
                   count1 = count1 + 1
                endif
             enddo
          enddo
       enddo
       
       write(ftn,rec=1) gtmp

!       do r=1,LDT_rc%gnr(n)
!          do c=1,LDT_rc%gnc(n)
!             line = c+(r-1)*LDT_rc%gnc(n)
!             write(ftn,rec=line) gtmp(c,r)
!          enddo
!       enddo
       deallocate(gtmp)
    endif
    deallocate(gtmp1)
    deallocate(var1)
  end subroutine writevar_bin_real_direct

!BOP
! !ROUTINE: writevar_restart_int
! \label{writevar_restart_int}
!
! !INTERFACE:
! Private name: call using LDT_writevar_restart
  subroutine writevar_restart_int(n,ftn,  var)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_masterproc, LDT_localPet, LDT_tdeltas,LDT_toffsets,&
         LDT_nss_halo_ind,LDT_nse_halo_ind,LDT_ews_halo_ind,&
         LDT_nss_ind,LDT_nse_ind,LDT_ews_ind,LDT_ewe_ind,&
         LDT_ewe_halo_ind,LDT_npes

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in) :: ftn
    integer, intent(inout) :: var(LDT_rc%ntiles(n))
! !DESCRIPTION:
!  Writes an integer variable to a binary restart file. 
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

    if(LDT_masterproc) then 
       allocate(gtmp(LDT_rc%glbntiles_red(n)))
       allocate(gtmp1(LDT_rc%glbntiles(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)      
    call MPI_GATHERV(var,LDT_tdeltas(n,LDT_localPet),&
         MPI_INTEGER,gtmp1,LDT_tdeltas(n,:),LDT_toffsets(n,:),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else 
    gtmp1 = var
#endif
    if(LDT_masterproc) then 
       count1=1
       do l=1,LDT_npes
          do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
             do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                if(r.ge.LDT_nss_ind(n,l).and.&
                     r.le.LDT_nse_ind(n,l).and.&
                     c.ge.LDT_ews_ind(n,l).and.&
                     c.le.LDT_ewe_ind(n,l))then !points not in halo

                   gid = c+(r-1)*LDT_rc%gnc(n)
                   ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                   stid = LDT_domain(n)%str_tind(gid)
                   do t=1,ntiles
                      tid = stid + t-1
                      gtmp(tid) = gtmp1(count1)
                      count1 = count1 + 1
                   enddo
                endif
             enddo
          enddo
       enddo
       write(ftn) gtmp
       deallocate(gtmp)
    endif
    deallocate(gtmp1)
  end subroutine writevar_restart_int

!BOP
! !ROUTINE: writevar_restart_real
! \label{writevar_restart_real}
! 
! !INTERFACE:
! Private name: call using LDT_writevar_restart
  subroutine writevar_restart_real(n,ftn,  var)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_masterproc,LDT_localPet, LDT_toffsets,LDT_tdeltas,&
         LDT_nss_halo_ind,LDT_nse_halo_ind,LDT_ews_halo_ind,&
         LDT_nss_ind,LDT_nse_ind,LDT_ews_ind,LDT_ewe_ind,&
         LDT_ewe_halo_ind,LDT_npes

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in) :: ftn
    real, intent(in)    :: var(LDT_rc%ntiles(n))
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
    
    if(LDT_masterproc) then 
       allocate(gtmp(LDT_rc%glbntiles_red(n)))
       allocate(gtmp1(LDT_rc%glbntiles(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)      
    call MPI_GATHERV(var,LDT_tdeltas(n,LDT_localPet),&
         MPI_REAL,gtmp1,LDT_tdeltas(n,:),LDT_toffsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
    gtmp1 = var
#endif
    if(LDT_masterproc) then 
       count1=1
       do l=1,LDT_npes
          do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
             do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                gid = c+(r-1)*LDT_rc%gnc(n)
                ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                stid = LDT_domain(n)%str_tind(gid)
                if(r.ge.LDT_nss_ind(n,l).and.&
                     r.le.LDT_nse_ind(n,l).and.&
                     c.ge.LDT_ews_ind(n,l).and.&
                     c.le.LDT_ewe_ind(n,l))then !points not in halo
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
  end subroutine writevar_restart_real

!BOP
! !ROUTINE: writevar_tile_noensemble_real
! \label{writevar_tile_noensemble_real}
! 
! !INTERFACE:
! Private name: call using LDT_writevar_reduced_tilespace
  subroutine writevar_tile_noensemble_real(n,ftn,  var)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_masterproc,LDT_localPet, LDT_toffsets,LDT_tdeltas,&
         LDT_nss_halo_ind,LDT_nse_halo_ind,LDT_ews_halo_ind,&
         LDT_nss_ind,LDT_nse_ind,LDT_ews_ind,LDT_ewe_ind,&
         LDT_ewe_halo_ind,LDT_npes

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in) :: ftn
    real, intent(in)    :: var(LDT_rc%ntiles(n)/LDT_rc%nensem(n))
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
    
    if(LDT_masterproc) then 
       allocate(gtmp(LDT_rc%glbntiles_red(n)/LDT_rc%nensem(n)))
       allocate(gtmp1(LDT_rc%glbntiles(n)/LDT_rc%nensem(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)      
    call MPI_GATHERV(var,LDT_tdeltas(n,LDT_localPet)/LDT_rc%nensem(n),&
         MPI_REAL,gtmp1,LDT_tdeltas(n,:)/LDT_rc%nensem(n),&
         LDT_toffsets(n,:)/LDT_rc%nensem(n),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
    gtmp1 = var
#endif
    if(LDT_masterproc) then 
       count1=1
       do l=1,LDT_npes
          do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
             do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                gid = c+(r-1)*LDT_rc%gnc(n)
                ntiles = LDT_domain(n)%ntiles_pergrid(gid)/LDT_rc%nensem(n)
                stid = (LDT_domain(n)%str_tind(gid)-1)/LDT_rc%nensem(n)+1
                if(r.ge.LDT_nss_ind(n,l).and.&
                     r.le.LDT_nse_ind(n,l).and.&
                     c.ge.LDT_ews_ind(n,l).and.&
                     c.le.LDT_ewe_ind(n,l))then !points not in halo
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
! Private name: call LDT_readvar_reduced_tilespace
  subroutine readvar_tile_noensemble_real(n,ftn,  var)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain,                 & 
                            LDT_nss_halo_ind,LDT_nse_halo_ind, &
                            LDT_ews_halo_ind,LDT_ewe_halo_ind, &
                            LDT_localPet

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in)   :: ftn
    real, intent(inout)   :: var(LDT_rc%ntiles(n)/LDT_rc%nensem(n))
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

    allocate(gtmp(LDT_rc%glbntiles_red(n)/LDT_rc%nensem(n)))
    read(ftn) gtmp

    count1=1
    do r=LDT_nss_halo_ind(n,LDT_localPet+1),LDT_nse_halo_ind(n,LDT_localPet+1)
       do c=LDT_ews_halo_ind(n,LDT_localPet+1),LDT_ewe_halo_ind(n,LDT_localPet+1)
          gid = c+(r-1)*LDT_rc%gnc(n)
          ntiles = LDT_domain(n)%ntiles_pergrid(gid)/LDT_rc%nensem(n)
          stid = (LDT_domain(n)%str_tind(gid)-1)/LDT_rc%nensem(n)+1
          do t=1,ntiles
             tid = stid + t-1
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
    deallocate(gtmp)   
  end subroutine readvar_tile_noensemble_real

!BOP
! !ROUTINE: readvar_restart_int
! \label{readvar_restart_int}
!
! !INTERFACE:
! Private name: call using LDT_readvar_restart
  subroutine readvar_restart_int(n,ftn,  var)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_nss_halo_ind,LDT_nse_halo_ind,LDT_ews_halo_ind,LDT_ewe_halo_ind,LDT_localPet
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in)    :: ftn
    integer, intent(inout) :: var(LDT_rc%ntiles(n))
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

    allocate(gtmp(LDT_rc%glbntiles_red(n)))
    read(ftn) gtmp
    
    count1=1
    do r=LDT_nss_halo_ind(n,LDT_localPet+1),LDT_nse_halo_ind(n,LDT_localPet+1)
       do c=LDT_ews_halo_ind(n,LDT_localPet+1),LDT_ewe_halo_ind(n,LDT_localPet+1)
          gid = c+(r-1)*LDT_rc%gnc(n)
          ntiles = LDT_domain(n)%ntiles_pergrid(gid)
          stid = LDT_domain(n)%str_tind(gid)
          do t=1,ntiles
             tid = stid + t-1
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
    deallocate(gtmp)

  end subroutine readvar_restart_int

!BOP
! !ROUTINE: readvar_restart_real
! \label{readvar_restart_real}
! 
! !INTERFACE:
! Private name: call LDT_readvar_restart
  subroutine readvar_restart_real(n,ftn,  var)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain,                 & 
                            LDT_nss_halo_ind,LDT_nse_halo_ind, &
                            LDT_ews_halo_ind,LDT_ewe_halo_ind, &
                            LDT_localPet

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in)   :: ftn
    real, intent(inout)   :: var(LDT_rc%ntiles(n))
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

    allocate(gtmp(LDT_rc%glbntiles_red(n)))
    read(ftn) gtmp

    count1=1
    do r=LDT_nss_halo_ind(n,LDT_localPet+1),LDT_nse_halo_ind(n,LDT_localPet+1)
       do c=LDT_ews_halo_ind(n,LDT_localPet+1),LDT_ewe_halo_ind(n,LDT_localPet+1)
          gid = c+(r-1)*LDT_rc%gnc(n)
          ntiles = LDT_domain(n)%ntiles_pergrid(gid)
          stid = LDT_domain(n)%str_tind(gid)
          do t=1,ntiles
             tid = stid + t-1
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
    deallocate(gtmp)   
  end subroutine readvar_restart_real

#if 0 
!BOP
! !ROUTINE: readvar_restart_char
! \label{readvar_restart_char}
!
! !INTERFACE:
! Private name: call using LDT_readvar_restart
  subroutine readvar_restart_char(ftn,  var)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain
    use spmdMod,       only : LDT_nss_ind,LDT_nse_ind,LDT_ews_ind,LDT_ewe_ind,LDT_localPet
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in)    :: ftn
    character*1, intent(inout) :: var(LDT_rc%ntiles(n))
! !DESCRIPTION:
!  Reads an character variable from a binary restart file. 
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
    character*1, allocatable :: gtmp(:)
    integer :: count1 ,c,r,ntiles,t,gid,stid,tid

    allocate(gtmp(LDT_rc%glbntiles(n)))
    read(ftn) gtmp
    
    count1=1
    do r=LDT_nss_ind(LDT_localPet+1),LDT_nse_ind(LDT_localPet+1)
       do c=LDT_ews_ind(LDT_localPet+1),LDT_ewe_ind(LDT_localPet+1)
          gid = c+(r-1)*LDT_rc%gnc(n)
          ntiles = LDT_domain(n)%ntiles_pergrid(gid)
          stid = LDT_domain(n)%str_tind(gid)
          do t=1,ntiles
             tid = stid + t-1
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
    deallocate(gtmp)

  end subroutine readvar_restart_char
#endif

!BOP
! !ROUTINE: readvar_1dgridded_real
! \label{readvar_1dgridded_real}
! 
! !INTERFACE:
! Private name: call LDT_readvar_gridded
  subroutine readvar_1dgridded_real(n,ftn,  var)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_localPet, LDT_nss_halo_ind,LDT_nse_halo_ind,&
         LDT_ews_halo_ind,LDT_ewe_halo_ind

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in)   :: ftn
    real, intent(inout)   :: var(LDT_rc%ngrid(n))
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

    allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
    read(ftn) gtmp
    
    nc = (LDT_ewe_halo_ind(n,LDT_localPet+1)-LDT_ews_halo_ind(n,LDT_localPet+1))+1

    do r=LDT_nss_halo_ind(n,LDT_localPet+1),LDT_nse_halo_ind(n,LDT_localPet+1)
       do c=LDT_ews_halo_ind(n,LDT_localPet+1),LDT_ewe_halo_ind(n,LDT_localPet+1)
          c1 = c-LDT_ews_halo_ind(n,LDT_localPet+1)+1
          r1 = r-LDT_nss_halo_ind(n,LDT_localPet+1)+1
!          gid = r1+(c1-1)*nc
          gid = LDT_domain(n)%gindex(c1,r1)
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
! Private name: call LDT_readvar_gridded
  subroutine readvar_1dgridded_fromvector_real(n,ftn,  var, oned)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_localPet,LDT_nss_halo_ind,LDT_nse_halo_ind,&
         LDT_ews_halo_ind,LDT_ewe_halo_ind

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in)   :: ftn
    real, intent(inout)   :: var(LDT_rc%ngrid(n))
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

    allocate(gtmp(LDT_rc%glbngrid(n)))
    read(ftn) gtmp

    allocate(gtmp2d(LDT_rc%gnc(n), LDT_rc%gnr(n)))

    cnt = 1
    do r=1,LDT_rc%gnr(n)
       do c=1,LDT_rc%gnc(n)
          gid = c+(r-1)*LDT_rc%gnc(n)
          if(LDT_domain(n)%ntiles_pergrid(gid).gt.0) then 
             gtmp2d(c,r) = gtmp(cnt)
             cnt = cnt+1
          endif
       enddo
    enddo

    nc = (LDT_ewe_halo_ind(n,LDT_localPet+1)-LDT_ews_halo_ind(n,LDT_localPet+1))+1

    do r=LDT_nss_halo_ind(n,LDT_localPet+1),LDT_nse_halo_ind(n,LDT_localPet+1)
       do c=LDT_ews_halo_ind(n,LDT_localPet+1),LDT_ewe_halo_ind(n,LDT_localPet+1)
          c1 = c-LDT_ews_halo_ind(n,LDT_localPet+1)+1
          r1 = r-LDT_nss_halo_ind(n,LDT_localPet+1)+1
!          gid = r1+(c1-1)*nc
          gid = LDT_domain(n)%gindex(c1,r1)
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
! Private name: call LDT_readvar_gridded
  subroutine readvar_2dgridded_real(n,ftn,  var)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_localPet, LDT_nss_halo_ind,LDT_nse_halo_ind,&
         LDT_ews_halo_ind,LDT_ewe_halo_ind

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in)   :: ftn
    real, intent(inout)   :: var(LDT_rc%lnc(n), LDT_rc%lnr(n))
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

    allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
    read(ftn) gtmp
    
    nc = (LDT_ewe_halo_ind(n,LDT_localPet+1)-LDT_ews_halo_ind(n,LDT_localPet+1))+1

    do r=LDT_nss_halo_ind(n,LDT_localPet+1),LDT_nse_halo_ind(n,LDT_localPet+1)
       do c=LDT_ews_halo_ind(n,LDT_localPet+1),LDT_ewe_halo_ind(n,LDT_localPet+1)
          c1 = c-LDT_ews_halo_ind(n,LDT_localPet+1)+1
          r1 = r-LDT_nss_halo_ind(n,LDT_localPet+1)+1
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
! Private name: call LDT_writevar_gridded
  subroutine writevar_gridded_real(n, ftn, var, varid, dim1, dim2, wopt)
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_localPet, LDT_masterproc, LDT_npes, &
         LDT_nss_ind,LDT_nse_ind,LDT_ews_ind,LDT_ewe_ind, &
         LDT_nss_halo_ind,LDT_nse_halo_ind,LDT_ews_halo_ind,LDT_ewe_halo_ind, &
         LDT_gdeltas, LDT_goffsets

    implicit none
! !ARGUMENTS: 
    integer, intent(in)           :: n 
    integer, intent(in)           :: ftn
    real                          :: var(LDT_rc%ngrid(n))
    integer                       :: varid
    integer, intent(in), optional :: dim1
    integer, intent(in), optional :: dim2
    character(len=*), intent(in), optional :: wopt
!
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
    integer           :: l,t
    integer           :: iret
    integer           :: c,r,gid,ntiles, count1, ierr
    character(20)     :: wopt_temp

    if(.not.(PRESENT(wopt))) then 
       wopt_temp = LDT_rc%wopt
    else
       wopt_temp = wopt
    endif
    
    allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))

    gtmp = LDT_rc%udef
    
    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)
          gid = LDT_domain(n)%gindex(c,r)
          if (gid.ne.-1) then
             gtmp(c,r) = var(gid)
          endif
       enddo
    enddo
    
    if(PRESENT(dim1).and.PRESENT(dim2)) then 
       iret = nf90_put_var(ftn,varid, gtmp, (/1,1,dim2,dim1/),&
            (/LDT_rc%gnc(n),LDT_rc%gnr(n),1,1/))
    elseif(PRESENT(dim1).and.(.not.(present(dim2)))) then 
       iret = nf90_put_var(ftn,varid, gtmp, (/1,1,dim1/),&
            (/LDT_rc%gnc(n),LDT_rc%gnr(n),1/))
    else
       iret = nf90_put_var(ftn,varid, gtmp, (/1,1/),&
            (/LDT_rc%gnc(n),LDT_rc%gnr(n)/))
    endif
    deallocate(gtmp)

#if 0 
    if(LDT_rc%wout_form.eq.1) then 
       if(wopt_temp.eq."2d gridspace") then 
          if(LDT_masterproc) then 
             allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
             allocate(gtmp1(LDT_rc%glbngrid(n)))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1))
             gtmp1 = 0.0
          endif
#if (defined SPMD)      
          call MPI_GATHERV(var,LDT_gdeltas(n,LDT_localPet),MPI_REAL,gtmp1,&
               LDT_gdeltas(n,:),LDT_goffsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
          gtmp1 = var
#endif
          if(LDT_masterproc) then 
             gtmp = LDT_rc%udef
             count1=1
             do l=1,LDT_npes
                do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
                   do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                      gid = c+(r-1)*LDT_rc%gnc(n)
                      ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                      if(ntiles.ne.0) then       
                         if(r.ge.LDT_nss_ind(n,l).and.&
                              r.le.LDT_nse_ind(n,l).and.&
                              c.ge.LDT_ews_ind(n,l).and.&
                              c.le.LDT_ewe_ind(n,l)) then !points not in halo                   
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
       elseif(wopt.eq."1d gridspace") then 
          if(LDT_masterproc) then 
             allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
             allocate(gtmp1(LDT_rc%glbngrid(n)))
             allocate(gtmp2(LDT_rc%glbngrid_red(n)))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1))
             gtmp1 = 0.0
          endif
#if (defined SPMD)      
          call MPI_GATHERV(var,LDT_gdeltas(n,LDT_localPet),MPI_REAL,gtmp1,&
               LDT_gdeltas(n,:),LDT_goffsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
          gtmp1 = var
#endif
          if(LDT_masterproc) then 
             gtmp = LDT_rc%udef
             count1=1
             do l=1,LDT_npes
                do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
                   do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                      gid = c+(r-1)*LDT_rc%gnc(n)
                      ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                      if(ntiles.ne.0) then       
                         if(r.ge.LDT_nss_ind(n,l).and.&
                              r.le.LDT_nse_ind(n,l).and.&
                              c.ge.LDT_ews_ind(n,l).and.&
                              c.le.LDT_ewe_ind(n,l)) then !points not in halo                   
                            gtmp(c,r) = gtmp1(count1)
                         endif
                         count1 = count1 + 1
                      endif
                   enddo
                enddo
             enddo
             
             count1 = 1
             do r=1,LDT_rc%gnr(n)
                do c=1,LDT_rc%gnc(n)
                   gid = c+(r-1)*LDT_rc%gnc(n)
                   ntiles = LDT_domain(n)%ntiles_pergrid(gid)
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
    elseif(LDT_rc%wout_form.eq.3) then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
       if(wopt.eq."2d gridspace") then 
          if(LDT_masterproc) then 
             allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
             allocate(gtmp1(LDT_rc%glbngrid(n)))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1))
             gtmp1 = 0.0
          endif
#if (defined SPMD)      
          call MPI_GATHERV(var,LDT_gdeltas(n,LDT_localPet),MPI_REAL,gtmp1,&
               LDT_gdeltas(n,:),LDT_goffsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
          gtmp1 = var
#endif
          if(LDT_masterproc) then 
             gtmp = LDT_rc%udef
             count1=1
             do l=1,LDT_npes
                do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
                   do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                      gid = c+(r-1)*LDT_rc%gnc(n)
                      ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                      if(ntiles.ne.0) then       
                         if(r.ge.LDT_nss_ind(n,l).and.&
                              r.le.LDT_nse_ind(n,l).and.&
                              c.ge.LDT_ews_ind(n,l).and.&
                              c.le.LDT_ewe_ind(n,l)) then !points not in halo                   
                            gtmp(c,r) = gtmp1(count1)
                         endif
                         count1 = count1 + 1
                      endif
                   enddo
                enddo
             enddo
             if(PRESENT(dim1)) then 
                iret = nf90_put_var(ftn,varid, gtmp, (/1,1,dim1/),&
                     (/LDT_rc%gnc(n),LDT_rc%gnr(n),1/))
             else
                iret = nf90_put_var(ftn,varid, gtmp, (/1,1/),&
                     (/LDT_rc%gnc(n),LDT_rc%gnr(n)/))
             endif
             deallocate(gtmp)
          endif
          deallocate(gtmp1)
       elseif(wopt.eq."1d gridspace") then 
          if(LDT_masterproc) then 
             allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
             allocate(gtmp1(LDT_rc%glbngrid(n)))
             allocate(gtmp2(LDT_rc%glbngrid_red(n)))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1))
             gtmp1 = 0.0
          endif
#if (defined SPMD)      
          call MPI_GATHERV(var,LDT_gdeltas(n,LDT_localPet),MPI_REAL,gtmp1,&
               LDT_gdeltas(n,:),LDT_goffsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
          gtmp1 = var
#endif
          if(LDT_masterproc) then 
             gtmp = LDT_rc%udef
             count1=1
             do l=1,LDT_npes
                do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
                   do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                      gid = c+(r-1)*LDT_rc%gnc(n)
                      ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                      if(ntiles.ne.0) then       
                         if(r.ge.LDT_nss_ind(n,l).and.&
                              r.le.LDT_nse_ind(n,l).and.&
                              c.ge.LDT_ews_ind(n,l).and.&
                              c.le.LDT_ewe_ind(n,l)) then !points not in halo                   
                            gtmp(c,r) = gtmp1(count1)
                         endif
                         count1 = count1 + 1
                      endif
                   enddo
                enddo
             enddo
             
             count1 = 1
             do r=1,LDT_rc%gnr(n)
                do c=1,LDT_rc%gnc(n)
                   gid = c+(r-1)*LDT_rc%gnc(n)
                   ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                   if(ntiles.ne.0) then    
                      gtmp2(count1) = gtmp(c,r)
                      count1 = count1 + 1
                   endif
                enddo
             enddo
             
             iret = nf90_put_var(ftn,varid,gtmp2,(/1/),(/LDT_rc%glbngrid_red(n)/))
             deallocate(gtmp)
             deallocate(gtmp2)
          endif
          deallocate(gtmp1)
       endif
#endif
    endif
#endif
  end subroutine writevar_gridded_real

!BOP
! !ROUTINE: writevar_gridded_real_3d
! \label{writevar_gridded_real_3d}
! 
! !INTERFACE:
! Private name: call LDT_writevar_gridded
  subroutine writevar_gridded_real_3d(n,ftn,  var, l1, l2, varid )
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_localPet, LDT_masterproc, LDT_npes, &
         LDT_nss_ind,LDT_nse_ind,LDT_ews_ind,LDT_ewe_ind, &
         LDT_nss_halo_ind,LDT_nse_halo_ind,LDT_ews_halo_ind,LDT_ewe_halo_ind, &
         LDT_gdeltas, LDT_goffsets

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in)             :: ftn
    integer                         :: l1
    integer                         :: l2
    real                            :: var(LDT_rc%ngrid(n),l1,l2)
    integer                         :: varid

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
    integer       :: iret

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
    
    iret = nf90_put_var(ftn,varid,var,(/1,1,1/),&
         (/LDT_rc%glbngrid(n),l1,l2/))
    call LDT_verify(iret, 'nf90_put_var failed in writevar_gridded_real_3d')

#endif


  end subroutine writevar_gridded_real_3d

!BOP
! !ROUTINE: writevar_gridded_real_4d
! \label{writevar_gridded_real_4d}
! 
! !INTERFACE:
! Private name: call LDT_writevar_gridded
  subroutine writevar_gridded_real_4d(n,ftn,  var, l1, l2, l3, varid )
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_localPet, LDT_masterproc, LDT_npes, &
         LDT_nss_ind,LDT_nse_ind,LDT_ews_ind,LDT_ewe_ind, &
         LDT_nss_halo_ind,LDT_nse_halo_ind,LDT_ews_halo_ind,LDT_ewe_halo_ind, &
         LDT_gdeltas, LDT_goffsets

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in)             :: ftn
    integer                         :: l1
    integer                         :: l2
    integer                         :: l3
    real                            :: var(LDT_rc%ngrid(n),l1,l2,l3)
    integer                         :: varid

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
    integer       :: iret

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
    
    iret = nf90_put_var(ftn,varid,var,(/1,1,1,1/),&
         (/LDT_rc%glbngrid(n),l1,l2,l3/))
    call LDT_verify(iret, 'nf90_put_var failed in writevar_gridded_real_4d')

#endif


  end subroutine writevar_gridded_real_4d

!BOP
! !ROUTINE: writevar_gridded_real_2d
! \label{writevar_gridded_real_2d}
! 
! !INTERFACE:
! Private name: call LDT_writevar_gridded
  subroutine writevar_gridded_real_2d(n,ftn,  var, l1, varid )
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, & 
         LDT_localPet, LDT_masterproc, LDT_npes, &
         LDT_nss_ind,LDT_nse_ind,LDT_ews_ind,LDT_ewe_ind, &
         LDT_nss_halo_ind,LDT_nse_halo_ind,LDT_ews_halo_ind,LDT_ewe_halo_ind, &
         LDT_gdeltas, LDT_goffsets

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in)             :: ftn
    integer                         :: l1
    real                            :: var(LDT_rc%ngrid(n),l1)
    integer                         :: varid

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
    integer       :: iret

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
    
    iret = nf90_put_var(ftn,varid,var,(/1,1/),&
         (/LDT_rc%glbngrid(n),l1/))
    call LDT_verify(iret, 'nf90_put_var failed in writevar_gridded_real_2d')

#endif


  end subroutine writevar_gridded_real_2d

!Y.Kwon
!BOP
! !ROUTINE: writevar_gridded_integer_1d
! \label{writevar_gridded_integer_1d}
! 
! !INTERFACE:
! Private name: call LDT_writevar_gridded
  subroutine writevar_gridded_integer_1d(n,ftn, var, varid )
! !USES:
    use LDT_coreMod, only : LDT_rc,LDT_domain, &
         LDT_localPet, LDT_masterproc, LDT_npes, &
         LDT_nss_ind,LDT_nse_ind,LDT_ews_ind,LDT_ewe_ind, &
         LDT_nss_halo_ind,LDT_nse_halo_ind,LDT_ews_halo_ind,LDT_ewe_halo_ind, &
         LDT_gdeltas, LDT_goffsets

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in)             :: ftn
    integer                         :: var(LDT_rc%ngrid(n))
    integer                         :: varid

! !DESCRIPTION:
!  Writes an integer variable to a binary 
!  sequential access, gridded file as either 1-dimensional gridded field. 
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
    integer       :: iret

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 

    iret = nf90_put_var(ftn,varid,var,(/1/),&
         (/LDT_rc%glbngrid(n)/))
    call LDT_verify(iret, 'nf90_put_var failed in writevar_gridded_integer_1d')

#endif

end subroutine writevar_gridded_integer_1d

!BOP
! !ROUTINE: LDT_tile2grid
! \label{LDT_tile2grid}
!
! !INTERFACE:
  subroutine LDT_tile2grid(n,gvar,tvar)
! !USES:
    use LDT_coreMod,      only : LDT_rc,LDT_domain

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n 
    real              :: gvar(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real, intent(in)  :: tvar(LDT_rc%ntiles(n))
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
    do i=1,LDT_rc%ntiles(n),LDT_rc%nensem(n)
       c = LDT_domain(n)%tile(i)%col
       r = LDT_domain(n)%tile(i)%row
       do m=1,LDT_rc%nensem(n)
          t = i+m-1
          gvar(c,r) = gvar(c,r)+&
               tvar(t)*LDT_domain(n)%tile(t)%fgrd/float(LDT_rc%nensem(n))
       enddo
    enddo
    
  end subroutine LDT_tile2grid


!BOP
!
! !ROUTINE: gather_tiled_vector_output
! \label{gather_tiled_vector_output}
!
! !REVISION HISTORY:
!  30 Jan 2009:  Sujay Kumar; Initial code
! 
! !INTERFACE:
subroutine gather_tiled_vector_output(n, gtmp, var)
! !USES: 
   use LDT_coreMod, only : LDT_rc, LDT_domain, LDT_localPet, LDT_masterproc, &
                           LDT_gdeltas, LDT_goffsets,&
                           LDT_nss_ind,LDT_ews_ind,LDT_nse_ind,&
                           LDT_ewe_ind, LDT_npes, &
                           LDT_nss_halo_ind,LDT_ews_halo_ind,&
                           LDT_nse_halo_ind,LDT_ewe_halo_ind,&
                           LDT_tdeltas, LDT_toffsets

   implicit none

! !ARGUMENTS: 
   integer, intent(in) :: n 
   real, allocatable :: gtmp(:)
   real          :: var(LDT_rc%ntiles(n))
!
! !DESCRIPTION:
! This routine gathers the output data into a tiled 1d array.
!
! This process aggregates the variable into a tile space. 
!
! This process accounts for the halo.
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [gtmp]
!     return array for the tiled output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

   integer :: ierr

   if(LDT_masterproc) then 
      allocate(gtmp(LDT_rc%glbntiles(n)))
      gtmp = 0.0
   else
      allocate(gtmp(1))
      gtmp = 0.0
   endif
   
#if (defined SPMD)      
   call MPI_GATHERV(var,LDT_tdeltas(n,LDT_localPet),&
        MPI_REAL,gtmp,LDT_tdeltas(n,:),LDT_toffsets(n,:),&
        MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
   gtmp = var
#endif
   
 end subroutine gather_tiled_vector_output
!BOP
!
! !ROUTINE: gather_gridded_output
! \label{gather_gridded_output}
!
! !REVISION HISTORY:
!  29 Oct 2008:  James Geiger; Initial code
! 
! !INTERFACE:
subroutine gather_gridded_output(n,gtmp, var)
! !USES: 
   use LDT_coreMod, only : LDT_rc, LDT_domain, LDT_localPet, LDT_masterproc, &
                           LDT_gdeltas, LDT_goffsets,&
                           LDT_nss_ind,LDT_ews_ind,LDT_nse_ind,&
                           LDT_ewe_ind, LDT_npes, &
                           LDT_nss_halo_ind,LDT_ews_halo_ind,&
                           LDT_nse_halo_ind,LDT_ewe_halo_ind,&
                           LDT_tdeltas, LDT_toffsets

   implicit none

! !ARGUMENTS: 
   integer, intent(in) :: n 
   real, allocatable :: gtmp(:,:)
   real          :: var(LDT_rc%ntiles(n))

! !DESCRIPTION:
! This routine gathers the output data into a gridded array.
!
! This process aggregates tiles to their grid, accounting for ensemble runs.
!
! This process accounts for the halo.
!
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
   integer :: count1,gid,ntiles,ierr

   allocate(var1(LDT_rc%ngrid(n)))

   if ( LDT_masterproc ) then 
      allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
      allocate(gtmp1(LDT_rc%glbngrid(n)))
      gtmp = 0.0
      gtmp1 = 0.0
   else
      allocate(gtmp1(1))
      gtmp1 = 0.0
   endif

   var1 = 0 

   do i=1,LDT_rc%ntiles(n),LDT_rc%nensem(n)
      c = LDT_domain(n)%tile(i)%index
      do m=1,LDT_rc%nensem(n)
         t = i+m-1
         if ( var(t) == -9999.0 ) then
            var1(c) = -9999.0
         else
            var1(c) = var1(c) + &
            var(t)*LDT_domain(n)%tile(t)%fgrd/float(LDT_rc%nensem(n))
         endif
      enddo
   enddo

#if (defined SPMD)      
   call MPI_GATHERV(var1,LDT_gdeltas(n,LDT_localPet),&
   MPI_REAL,gtmp1,LDT_gdeltas(n,:),LDT_goffsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
   gtmp1 = var1
#endif

   if ( LDT_masterproc ) then 
      gtmp = LDT_rc%udef
      count1=1

      do l=1,LDT_npes
         do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
            do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
               gid = c+(r-1)*LDT_rc%gnc(n)
               ntiles = LDT_domain(n)%ntiles_pergrid(gid)
               if(ntiles.ne.0) then                          
                  if(r.ge.LDT_nss_ind(n,l).and.&
                     r.le.LDT_nse_ind(n,l).and.&
                     c.ge.LDT_ews_ind(n,l).and.&
                     c.le.LDT_ewe_ind(n,l))then !points not in halo
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
end subroutine gather_gridded_output


!BOP
!
! !ROUTINE: gather_gridded_vector_output
! \label{gather_gridded_vector_output}
!
! !REVISION HISTORY:
!  30 Jan 2009: Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine gather_gridded_vector_output(n, gtmp, var)
! !USES: 
   use LDT_coreMod, only : LDT_rc, LDT_domain, LDT_localPet, LDT_masterproc, &
                           LDT_gdeltas, LDT_goffsets,&
                           LDT_nss_ind,LDT_ews_ind,LDT_nse_ind,&
                           LDT_ewe_ind, LDT_npes, &
                           LDT_nss_halo_ind,LDT_ews_halo_ind,&
                           LDT_nse_halo_ind,LDT_ewe_halo_ind,&
                           LDT_tdeltas, LDT_toffsets

   implicit none

! !ARGUMENTS: 
   integer, intent(in) :: n 
   real, allocatable :: gtmp(:)
   real          :: var(LDT_rc%ntiles(n))

! 
! !DESCRIPTION:
! This routine gathers the output data into a gridded 1d array.
!
! This process aggregates tiles to their grid, accounting for ensemble runs.
!
! This process accounts for the halo.
!    
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

   allocate(var1(LDT_rc%ngrid(n)))
   if ( LDT_masterproc ) then 
      allocate(gtmp(LDT_rc%glbngrid(n)))
      allocate(gtmp2d(LDT_rc%gnc(n), LDT_rc%gnr(n)))
      gtmp = 0.0     
   else
      allocate(gtmp(1))
      allocate(gtmp2d(1,1))
      gtmp = 0.0
   endif

   var1 = 0 

   do i=1,LDT_rc%ntiles(n),LDT_rc%nensem(n)
      c = LDT_domain(n)%tile(i)%index
      do m=1,LDT_rc%nensem(n)
         t = i+m-1
         if ( var(t) == -9999.0 ) then
            var1(c) = -9999.0
         else
            var1(c) = var1(c) + &
            var(t)*LDT_domain(n)%tile(t)%fgrd/float(LDT_rc%nensem(n))
         endif
      enddo
   enddo

#if (defined SPMD)      
   call MPI_GATHERV(var1,LDT_gdeltas(n,LDT_localPet),&
   MPI_REAL,gtmp,LDT_gdeltas(n,:),LDT_goffsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
   gtmp= var1
#endif
! Though gtmp now has a size of glbngrid, the points are not ordered correctly. 
   if ( LDT_masterproc ) then 
      gtmp2d = LDT_rc%udef
      count1=1

      do l=1,LDT_npes
         do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
            do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
               gid = c+(r-1)*LDT_rc%gnc(n)
               ntiles = LDT_domain(n)%ntiles_pergrid(gid)
               if(ntiles.ne.0) then                          
                  if(r.ge.LDT_nss_ind(n,l).and.&
                     r.le.LDT_nse_ind(n,l).and.&
                     c.ge.LDT_ews_ind(n,l).and.&
                     c.le.LDT_ewe_ind(n,l))then !points not in halo
                     gtmp2d(c,r) = gtmp(count1)
                  endif
                  count1 = count1 + 1
               endif
            enddo
         enddo
      enddo
! now reorder the gtmp2d to create a 1d vector
      count1 = 1
      do r=1,LDT_rc%gnr(n)
         do c=1,LDT_rc%gnc(n)
            if(LDT_domain(n)%ntiles_pergrid(c+(r-1)*LDT_rc%gnc(n)).gt.0) then 
               gtmp(count1) = gtmp2d(c,r)
               count1 = count1+1
            endif
         enddo
      enddo
      
   endif
   deallocate(gtmp2d)
   deallocate(var1)
  
 end subroutine gather_gridded_vector_output

!BOP
!
! !ROUTINE: writeNETCDFdataHeader_LIS
! \label{writeNETCDFdataHeader_LIS}
! 
! !INTERFACE: 
 subroutine writeNETCDFdataHeader_LIS(n,ftn,dimID, paramEntry, wtype, fdimID)
! 
! !DESCRIPTION: 
!  This subroutine writes the NetCDF data headers in the 
!  standard preprocessing mode for LIS
!EOP
 
   integer               :: n 
   integer               :: ftn
   integer               :: dimID(3)
   type(LDT_paramEntry)  :: paramEntry
   character(len=*), intent(in), optional :: wtype
   integer, intent(in),optional           :: fdimID(4)  !4D dimid
   integer               :: shuffle, deflate, deflate_level
   integer               :: wtype_temp

! Added by Hiroko to allow handling of double or integer data types
! in netCDF format for CLM-4.5
   if(.not.(PRESENT(wtype))) then
      wtype_temp = nf90_float
   else 
    if ( trim(wtype) .eq. "double" ) then
      wtype_temp = nf90_double
    else
      wtype_temp = nf90_int
    endif
   endif

   shuffle = NETCDF_shuffle
   deflate = NETCDF_deflate
   deflate_level =NETCDF_deflate_level
   
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    if(paramEntry%selectOpt.gt.0) then 
       if(paramEntry%vlevels.gt.1) then 
         if(.not.(PRESENT(fdimID))) then
          call LDT_verify(nf90_def_var(ftn,trim(paramEntry%short_name),&
               wtype_temp, dimids = dimID, varID=paramEntry%vid))
         else
          call LDT_verify(nf90_def_var(ftn,trim(paramEntry%short_name),&
               wtype_temp, dimids = fdimID, varID=paramEntry%vid))
         endif
       else
          call LDT_verify(nf90_def_var(ftn,trim(paramEntry%short_name),&
               wtype_temp, dimids = dimID(1:2), varID=paramEntry%vid), &
               'nf90_def_var failed')
       endif

#if (defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(ftn,&
            paramEntry%vid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed')
#endif
    endif
#endif
  end subroutine WriteNETCDFdataHeader_LIS

!BOP
!
! !ROUTINE: writeNETCDFdataHeader_LISHydro
! \label{writeNETCDFdataHeader_LISHydro}
! 
! !INTERFACE: 
 subroutine writeNETCDFdataHeader_LISHydro(n,ftn,dimID, paramEntry, flag)
   
! !DESCRIPTION: 
!  This subroutine writes the NetCDF data headers in the 
!  preprocessing mode for LISHydro
!EOP

   integer               :: n 
   integer               :: ftn
   integer               :: dimID(4)
   integer               :: shuffle, deflate, deflate_level
   type(LDT_paramEntry)  :: paramEntry
   integer               :: flag
    
   shuffle = NETCDF_shuffle
   deflate = NETCDF_deflate
   deflate_level =NETCDF_deflate_level
   
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    if(paramEntry%selectOpt.gt.0) then 
       if(paramEntry%vlevels.gt.1) then 
          call LDT_verify(nf90_def_var(ftn,trim(paramEntry%short_name),&
               nf90_float, dimids = dimID, varID=paramEntry%vid))
       else if (dimID(3) == 3 ) then
          call LDT_verify(nf90_def_var(ftn,trim(paramEntry%short_name),&
               nf90_float, dimids = dimID(1:3), varID=paramEntry%vid))
       else
          call LDT_verify(nf90_def_var(ftn,trim(paramEntry%short_name),&
               nf90_float, dimids = dimID(1:2), varID=paramEntry%vid))
       endif

#if (defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(ftn,&
            paramEntry%vid, shuffle, deflate, deflate_level))
#endif
    endif
#endif
  end subroutine writeNETCDFdataHeader_LISHydro
  
  subroutine writeNETCDFdata_default( n, ftn, paramEntry )
    
   integer              :: n 
   integer              :: ftn
   type(LDT_paramEntry) :: paramEntry    
   
   integer           :: c,r,k,l,count,ierr
   real, allocatable :: ltmp(:,:)
   real, allocatable :: gtmp1(:,:),gtmp2(:,:,:)

!SVK edit
#if(defined USE_NETCDF3 || defined USE_NETCDF4)   
   if(paramEntry%selectOpt.gt.0) then 
      
      allocate(ltmp(LDT_rc%lnc(n)*LDT_rc%lnr(n),paramEntry%vlevels))
      
      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            ltmp(c+(r-1)*LDT_rc%lnc(n),:) = paramEntry%value(c,r,:)
         enddo
      enddo
      if(LDT_masterproc) then 
         allocate(gtmp1(LDT_rc%gnc(n)*LDT_rc%gnr(n),paramEntry%vlevels))
         allocate(gtmp2(LDT_rc%gnc(n),LDT_rc%gnr(n),paramEntry%vlevels))
      else
         allocate(gtmp1(1,paramEntry%vlevels))
         allocate(gtmp2(1,1,paramEntry%vlevels))
      endif
      
      do k=1,paramEntry%vlevels
#if (defined SPMD) 
         call MPI_GATHERV(ltmp(:,k),LDT_deltas(n,LDT_localPet),MPI_REAL,gtmp1(:,k),&
              LDT_deltas(n,:),LDT_offsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else
         gtmp1(:,k) = ltmp(:,k)
#endif
      enddo
      if(LDT_masterproc) then 

         count =1 
         do l=1,LDT_npes
            do r=LDT_nss_ind(n,l), LDT_nse_ind(n,l)
                do c=LDT_ews_ind(n,l), LDT_ewe_ind(n,l)
                   gtmp2(c,r,:) = gtmp1(count,:)
                   count = count+1
                enddo
             enddo
          enddo
         if(paramEntry%vlevels.gt.1) then 
            call LDT_verify(nf90_put_var(ftn, paramEntry%vid,gtmp2,&
                 (/1,1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n),&
                 paramEntry%vlevels/)), &
                 'error in nf90_put_var in LDT_writeNETCDFdata (writeNETCDFdata_default)')
         else
            call LDT_verify(nf90_put_var(ftn, paramEntry%vid,gtmp2(:,:,1),&
                 (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
                 'nf90_put_var failed in LDT_writeNETCDFdata (writeNETCDFdata_default)')
         endif

#if ( defined USE_NETCDF3 )
         call LDT_verify(nf90_redef(ftn))
#endif
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "standard_name",trim(paramEntry%standard_name)))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "units",trim(paramEntry%units)))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "scale_factor",1.0))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "add_offset",0.0))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "missing_value",LDT_rc%udef))
         !      call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
         !           "_FillValue",LDT_rc%udef))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "vmin",paramEntry%valid_min))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "vmax",paramEntry%valid_max))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "num_bins",paramEntry%num_bins))

#if ( defined USE_NETCDF3 )
         call LDT_verify(nf90_enddef(ftn))
#endif
      endif
      deallocate(ltmp)
      deallocate(gtmp1)
      deallocate(gtmp2)
   endif
#endif

 end subroutine writeNETCDFdata_default


 subroutine writeNETCDFdata_givendomain( n, ftn, &
                 paramEntry, gnc, gnr )

   integer              :: n
   integer              :: ftn
   type(LDT_paramEntry) :: paramEntry
   integer              :: gnc
   integer              :: gnr

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   if(paramEntry%selectOpt.gt.0) then

   ! For parallel processing, only write data once for master processor:
     if(LDT_masterproc) then

       if(paramEntry%vlevels.gt.1) then
          call LDT_verify(nf90_put_var(ftn, paramEntry%vid, paramEntry%value,&
               (/1,1,1/),(/gnc,gnr,&
               paramEntry%vlevels/)),&
               'nf90_put_var failed in LDT_writeNETCDFdata (writeNETCDFdata_givendomain)')
       else
          call LDT_verify(nf90_put_var(ftn, paramEntry%vid, paramEntry%value,&
               (/1,1/),(/gnc,gnr/)),&
               'nf90_put_var failed in LDT_writeNETCDFdata (writeNETCDFdata_givendomain)')
       endif

       call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
            "standard_name",trim(paramEntry%standard_name)))
       call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
            "units",trim(paramEntry%units)))
       call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
            "scale_factor",1.0))
       call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
            "add_offset",0.0))
       call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
            "missing_value",LDT_rc%udef))
 !      call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
 !           "_FillValue",LDT_rc%udef))
       call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
            "vmin",paramEntry%valid_min))
       call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
            "vmax",paramEntry%valid_max))
       call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
            "num_bins",paramEntry%num_bins))

     endif
   endif
#endif

 end subroutine writeNETCDFdata_givendomain

 subroutine writeNETCDFdata_giventype( n, ftn, paramEntry, wtype )
   
   integer              :: n 
   integer              :: ftn
   type(LDT_paramEntry) :: paramEntry    
   character(len=*)     :: wtype
   
   integer           :: c,r,k,l,count,ierr,z
   real*8, allocatable :: dltmp(:,:,:)
   real*8, allocatable :: dgtmp1(:,:,:),dgtmp2(:,:,:,:)
   integer, allocatable :: iltmp(:,:,:)
   integer, allocatable :: igtmp1(:,:,:),igtmp2(:,:,:,:)

!SVK edit
#if(defined USE_NETCDF3 || defined USE_NETCDF4)   
   if(paramEntry%selectOpt.gt.0) then 
      
    if ( trim(wtype) .eq. "double" ) then
      allocate(dltmp(LDT_rc%lnc(n)*LDT_rc%lnr(n),paramEntry%vlevels,paramEntry%zlevels))
      
      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            if ( paramEntry%zlevels.gt.1 ) then
             dltmp(c+(r-1)*LDT_rc%lnc(n),:,:) = paramEntry%dvalue(c,r,:,:)
            else
             dltmp(c+(r-1)*LDT_rc%lnc(n),:,1) = paramEntry%dvalue(c,r,:,1)
            endif
         enddo
      enddo
      if(LDT_masterproc) then 
         allocate(dgtmp1(LDT_rc%gnc(n)*LDT_rc%gnr(n),paramEntry%vlevels,paramEntry%zlevels))
         allocate(dgtmp2(LDT_rc%gnc(n),LDT_rc%gnr(n),paramEntry%vlevels,paramEntry%zlevels))
      else
         allocate(dgtmp1(1,paramEntry%vlevels,paramEntry%zlevels))
         allocate(dgtmp2(1,1,paramEntry%vlevels,paramEntry%zlevels))
      endif
      
     do z=1,paramEntry%zlevels
      do k=1,paramEntry%vlevels
#if (defined SPMD) 
         call MPI_GATHERV(dltmp(:,k,z),LDT_deltas(n,LDT_localPet),MPI_DOUBLE_PRECISION,dgtmp1(:,k,z),&
              LDT_deltas(n,:),LDT_offsets(n,:),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#else
         dgtmp1(:,k,z) = dltmp(:,k,z)
#endif
      enddo
     enddo
      if(LDT_masterproc) then 
!         print*, 'writing att ',trim(paramEntry%short_name), LDT_localPet
         count =1 
         do l=1,LDT_npes
            do r=LDT_nss_ind(n,l), LDT_nse_ind(n,l)
                do c=LDT_ews_ind(n,l), LDT_ewe_ind(n,l)
                   dgtmp2(c,r,:,:) = dgtmp1(count,:,:)
                   count = count+1
                enddo
             enddo
          enddo
         if(paramEntry%vlevels.gt.1) then 
          if(paramEntry%zlevels.gt.1) then 
            call LDT_verify(nf90_put_var(ftn, paramEntry%vid,dgtmp2,&
                 (/1,1,1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n),&
                 paramEntry%vlevels,paramEntry%zlevels/)), &
                 'error in nf90_put_var in LDT_writeNETCDFdata (writeNETCDFdata_giventype)')
          else
            call LDT_verify(nf90_put_var(ftn, paramEntry%vid,dgtmp2(:,:,:,1),&
                 (/1,1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n),&
                 paramEntry%vlevels/)), &
                 'error in nf90_put_var in LDT_writeNETCDFdata (writeNETCDFdata_giventype)')
          endif
         else
            call LDT_verify(nf90_put_var(ftn, paramEntry%vid,dgtmp2(:,:,1,1),&
                 (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
                 'nf90_put_var failed in LDT_writeNETCDFdata (writeNETCDFdata_giventype)')
         endif

#if ( defined USE_NETCDF3 )
         call LDT_verify(nf90_redef(ftn))
#endif
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "standard_name",trim(paramEntry%standard_name)))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "units",trim(paramEntry%units)))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "scale_factor",1.0))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "add_offset",0.0))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "missing_value",LDT_rc%udef))
         !      call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
         !           "_FillValue",LDT_rc%udef))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "vmin",paramEntry%valid_min))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "vmax",paramEntry%valid_max))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "num_bins",paramEntry%num_bins))

#if ( defined USE_NETCDF3 )
         call LDT_verify(nf90_enddef(ftn))
#endif
      endif
      deallocate(dltmp)
      deallocate(dgtmp1)
      deallocate(dgtmp2)

     elseif ( trim(wtype) .eq. "int" ) then
      allocate(iltmp(LDT_rc%lnc(n)*LDT_rc%lnr(n),paramEntry%vlevels,paramEntry%zlevels))
      
      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            if(paramEntry%zlevels.gt.1) then 
             iltmp(c+(r-1)*LDT_rc%lnc(n),:,:) = paramEntry%dvalue(c,r,:,:)
            else
             iltmp(c+(r-1)*LDT_rc%lnc(n),:,1) = paramEntry%dvalue(c,r,:,1)
            endif
         enddo
      enddo
      if(LDT_masterproc) then 
         allocate(igtmp1(LDT_rc%gnc(n)*LDT_rc%gnr(n),paramEntry%vlevels,paramEntry%zlevels))
         allocate(igtmp2(LDT_rc%gnc(n),LDT_rc%gnr(n),paramEntry%vlevels,paramEntry%zlevels))
      else
         allocate(igtmp1(1,paramEntry%vlevels,paramEntry%zlevels))
         allocate(igtmp2(1,1,paramEntry%vlevels,paramEntry%zlevels))
      endif
      
     do z=1,paramEntry%zlevels
      do k=1,paramEntry%vlevels
#if (defined SPMD) 
         call MPI_GATHERV(iltmp(:,k,z),LDT_deltas(n,LDT_localPet),MPI_INTEGER,igtmp1(:,k,z),&
              LDT_deltas(n,:),LDT_offsets(n,:),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else
         igtmp1(:,k,z) = iltmp(:,k,z)
#endif
      enddo
     enddo
      if(LDT_masterproc) then 
!         print*, 'writing att ',trim(paramEntry%short_name), LDT_localPet
         count =1 
         do l=1,LDT_npes
            do r=LDT_nss_ind(n,l), LDT_nse_ind(n,l)
                do c=LDT_ews_ind(n,l), LDT_ewe_ind(n,l)
                   igtmp2(c,r,:,:) = igtmp1(count,:,:)
                   count = count+1
                enddo
             enddo
          enddo
         if(paramEntry%vlevels.gt.1) then 
          if(paramEntry%zlevels.gt.1) then 
            call LDT_verify(nf90_put_var(ftn, paramEntry%vid,igtmp2,&
                 (/1,1,1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n),&
                 paramEntry%vlevels,paramEntry%zlevels/)), &
                 'error in nf90_put_var in LDT_writeNETCDFdata (writeNETCDFdata_giventype)')
          else
            call LDT_verify(nf90_put_var(ftn, paramEntry%vid,igtmp2(:,:,:,1),&
                 (/1,1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n),&
                 paramEntry%vlevels/)), &
                 'error in nf90_put_var in LDT_writeNETCDFdata (writeNETCDFdata_giventype)')
          endif
         else
            call LDT_verify(nf90_put_var(ftn, paramEntry%vid,igtmp2(:,:,1,1),&
                 (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
                 'nf90_put_var failed in LDT_writeNETCDFdata (writeNETCDFdata_giventype)')
         endif


#if ( defined USE_NETCDF3 )
         call LDT_verify(nf90_redef(ftn))
#endif
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "standard_name",trim(paramEntry%standard_name)))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "units",trim(paramEntry%units)))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "scale_factor",1.0))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "add_offset",0.0))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "missing_value",LDT_rc%udef))
         !      call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
         !           "_FillValue",LDT_rc%udef))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "vmin",paramEntry%valid_min))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "vmax",paramEntry%valid_max))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "num_bins",paramEntry%num_bins))

#if ( defined USE_NETCDF3 )
         call LDT_verify(nf90_enddef(ftn))
#endif
      endif
      deallocate(iltmp)
      deallocate(igtmp1)
      deallocate(igtmp2)
     endif ! if wtype
   endif
#endif

 end subroutine writeNETCDFdata_giventype

!BOP
! !ROUTINE: writeNETCDFdata_LISHydro
!  \label{writeNETCDFdata_LISHydro}
! 
! !INTERFACE: 
 subroutine writeNETCDFdata_LISHydro( n, ftn, paramEntry,flag )
! 
! !DESCRIPTION: 
!  This routine writes the data for a given parameter entry 
!  to the NetCDF file, in the preprocessing mode for LISHydro
!
!EOP    
   integer              :: n 
   integer              :: ftn
   type(LDT_paramEntry) :: paramEntry    
   
   integer           :: c,r,k,l,count,ierr
   real, allocatable :: ltmp(:,:)
   real, allocatable :: gtmp1(:,:),gtmp2(:,:,:)
   integer           :: flag
!SVK edit
#if(defined USE_NETCDF3 || defined USE_NETCDF4)   
   if(paramEntry%selectOpt.gt.0) then 
      
      allocate(ltmp(LDT_rc%lnc(n)*LDT_rc%lnr(n),paramEntry%vlevels))
      
      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            ltmp(c+(r-1)*LDT_rc%lnc(n),:) = paramEntry%value(c,r,:)
         enddo
      enddo
      if(LDT_masterproc) then 
         allocate(gtmp1(LDT_rc%gnc(n)*LDT_rc%gnr(n),paramEntry%vlevels))
         allocate(gtmp2(LDT_rc%gnc(n),LDT_rc%gnr(n),paramEntry%vlevels))
      else
         allocate(gtmp1(1,paramEntry%vlevels))
         allocate(gtmp2(1,1,paramEntry%vlevels))
      endif
      
      do k=1,paramEntry%vlevels
#if (defined SPMD) 
         call MPI_GATHERV(ltmp(:,k),LDT_deltas(n,LDT_localPet),MPI_REAL,gtmp1(:,k),&
              LDT_deltas(n,:),LDT_offsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else
         gtmp1(:,k) = ltmp(:,k)
#endif
      enddo
      if(LDT_masterproc) then 
         count =1 
         do l=1,LDT_npes
            do r=LDT_nss_ind(n,l), LDT_nse_ind(n,l)
                do c=LDT_ews_ind(n,l), LDT_ewe_ind(n,l)
                   gtmp2(c,r,:) = gtmp1(count,:)
                   count = count+1
                enddo
             enddo
          enddo
         if(paramEntry%vlevels.gt.1) then 
            call LDT_verify(nf90_put_var(ftn, paramEntry%vid,gtmp2,&
                 (/1,1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n),&
                 paramEntry%vlevels/)), &
                 'error in nf90_put_var in LDT_writeNETCDFdata (writeNETCDFdata_LISHydro)')
         else
            call LDT_verify(nf90_put_var(ftn, paramEntry%vid,gtmp2(:,:,1),&
                 (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
                 'nf90_put_var failed in LDT_writeNETCDFdata (writeNETCDFdata_LISHydro)')
         endif

#if ( defined USE_NETCDF3 )
         call LDT_verify(nf90_redef(ftn))
#endif
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "standard_name",trim(paramEntry%standard_name)))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "units",trim(paramEntry%units)))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "scale_factor",1.0))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "add_offset",0.0))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "vmin",paramEntry%valid_min))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "vmax",paramEntry%valid_max))
         call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "num_bins",paramEntry%num_bins))
           call LDT_verify(nf90_put_att(ftn,paramEntry%vid,&
              "stagger", "M" ))

#if ( defined USE_NETCDF3 )
         call LDT_verify(nf90_enddef(ftn))
#endif
      endif
      deallocate(ltmp)
      deallocate(gtmp1)
      deallocate(gtmp2)
   endif
#endif

 end subroutine writeNETCDFdata_LISHydro

!BOP
! 
! !ROUTINE: readLISSingleNetcdfVar_LDTgrid
! \label{LDT_readLISSingleNetcdfVar_LDTgrid}
! 
! !INTERFACE:
  subroutine readLISSingleNetcdfVar_LDTgrid(n, ftn, vname, vlevel, value2d)
! !USES: 
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    use netcdf
#endif
    use LDT_coreMod, only     : LDT_rc, LDT_domain
    use LDT_logMod,  only     : LDT_logunit, LDT_verify, LDT_endrun

    implicit none
! !ARGUMENTS: 
    integer, intent(in)     :: n 
    integer                 :: ftn
    character(len=*)        :: vname
    integer                 :: vlevel
    real                    :: value2d(LDT_rc%lnc(n),LDT_rc%lnr(n))
! 
! !DESCRIPTION: 
!   This routine reads a single variable from the LIS output file in the 
!   NetCDF format, assumed to be in the LDT grid.
!   Based on the external datamask, the routine also filters
!   each variable. 
!
!EOP
    real, allocatable :: value1d(:)
    integer :: unit_id
    integer :: varid
    integer :: iret
    integer :: i,c,r,gid
    character*50 :: units_tmp

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    if(LDT_rc%lis_wopt.eq.1) then 
       
       write(*,*) '[ERR] Reading in tile space format is not yet supported'
       write(*,*) '[ERR] program stopping ...'
       stop
       
    elseif(LDT_rc%lis_wopt.eq.2) then 
       
       iret = nf90_inq_varid(ftn, trim(vname), varid)
       call LDT_verify(iret, 'Error in nf90_inq_varid: '//trim(vname))
       
       iret = nf90_get_att(ftn,varid, "units", units_tmp)
       call LDT_verify(iret, 'Error in nf90_get_att:units for '//trim(vname))

       iret = nf90_get_var(ftn, varid, value2d,start = (/1,1,vlevel/),&
          count = (/LDT_rc%lnc(n),LDT_rc%lnr(n),1/))
       call LDT_verify(iret,'Error in nf90_get_var: '//trim(vname))
       
       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(LDT_domain(n)%datamask(c,r).ne.1) then 
                value2d(c,r) = LDT_rc%udef
             endif
          enddo
       enddo
       
    elseif(LDT_rc%lis_wopt.eq.3) then 
       
       allocate(value1d(LDT_rc%ngrid(n)))
       
       iret = nf90_inq_varid(ftn, trim(vname), varid)
       call LDT_verify(iret, 'Error in nf90_inq_varid: '//trim(vname))

       iret = nf90_get_att(ftn,varid, "units", units_tmp)
       call LDT_verify(iret, 'Error in nf90_get_att:units for '//trim(vname))
       
       iret = nf90_get_var(ftn, varid, value1d,&
            start=(/1,vlevel/),count=(/LDT_rc%glbngrid(n),1/))
       call LDT_verify(iret,'Error in nf90_get_var: '//trim(vname))
              
       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(LDT_domain(n)%datamask(c,r).ne.1) then 
                value2d(c,r) = LDT_rc%udef
             else
                gid = LDT_domain(n)%gindex(c,r)
                value2d(c,r) = value1d(gid)
             endif
          enddo
       enddo
          
       deallocate(value1d)
       
    endif
    
#endif
  end subroutine readLISSingleNetcdfVar_LDTgrid

!BOP
! 
! !ROUTINE: readLISSingleNetcdfVar_LDTgrid_withunits
! \label{LDT_readLISSingleNetcdfVar_LDTgrid_withunits}
! 
! !INTERFACE:
  subroutine readLISSingleNetcdfVar_LDTgrid_withunits(n, &
       ftn, vname, vlevel, value2d, &
       units)
! !USES: 
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    use netcdf
#endif
    use LDT_coreMod, only     : LDT_rc, LDT_domain
    use LDT_logMod,  only     : LDT_logunit, LDT_verify, LDT_endrun

    implicit none
! !ARGUMENTS: 
    integer, intent(in)     :: n 
    integer                 :: ftn
    character(len=*)        :: vname
    integer                 :: vlevel
    real                    :: value2d(LDT_rc%lnc(n),LDT_rc%lnr(n))
    character(len=*)        :: units
! 
! !DESCRIPTION: 
!   This routine reads a single variable from the LIS output file in the 
!   NetCDF format, assumed to be in the LDT grid.
!   Based on the external datamask, the routine also filters
!   each variable. 
!
!EOP
    real, allocatable :: value1d(:)
    integer :: unit_id
    integer :: varid
    integer :: iret
    integer :: i,c,r,gid
    character*50 :: units_tmp

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    if(LDT_rc%lis_wopt.eq.1) then 
       
       write(*,*) '[ERR] Reading in tile space format is not yet supported'
       write(*,*) '[ERR] program stopping ...'
       stop
       
    elseif(LDT_rc%lis_wopt.eq.2) then 
       
       iret = nf90_inq_varid(ftn, trim(vname), varid)
       call LDT_verify(iret, 'Error in nf90_inq_varid: '//trim(vname))
       
       iret = nf90_get_att(ftn,varid, "units", units_tmp)
       call LDT_verify(iret, 'Error in nf90_get_att:units for '//trim(vname))

       iret = nf90_get_var(ftn, varid, value2d,start = (/1,1,vlevel/),&
          count = (/LDT_rc%lnc(n),LDT_rc%lnr(n),1/))
       call LDT_verify(iret,'Error in nf90_get_var: '//trim(vname))
       
       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(LDT_domain(n)%datamask(c,r).ne.1) then 
                value2d(c,r) = LDT_rc%udef
             endif
          enddo
       enddo
       
    elseif(LDT_rc%lis_wopt.eq.3) then 
       
       allocate(value1d(LDT_rc%ngrid(n)))
       
       iret = nf90_inq_varid(ftn, trim(vname), varid)
       call LDT_verify(iret, 'Error in nf90_inq_varid: '//trim(vname))

       iret = nf90_get_att(ftn,varid, "units", units_tmp)
       call LDT_verify(iret, 'Error in nf90_get_att:units for '//trim(vname))
       
       iret = nf90_get_var(ftn, varid, value1d,&
            start=(/1,vlevel/),count=(/LDT_rc%glbngrid(n),1/))
       call LDT_verify(iret,'Error in nf90_get_var: '//trim(vname))
              
       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(LDT_domain(n)%datamask(c,r).ne.1) then 
                value2d(c,r) = LDT_rc%udef
             else
                gid = LDT_domain(n)%gindex(c,r)
                value2d(c,r) = value1d(gid)
             endif
          enddo
       enddo
          
       deallocate(value1d)
       
    endif
    
    units = units_tmp 

#endif
  end subroutine readLISSingleNetcdfVar_LDTgrid_withunits


!BOP
! 
! !ROUTINE: readLISSingleNetcdfVar_Inputgrid
! \label{LDT_readLISSingleNetcdfVar_Inputgrid}
! 
! !INTERFACE:
  subroutine readLISSingleNetcdfVar_Inputgrid(n, ftn, vname, vlevel, &
       inc,inr,value2d)
! !USES: 
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    use netcdf
#endif
    use LDT_coreMod, only     : LDT_rc, LDT_domain
    use LDT_logMod,  only     : LDT_logunit, LDT_verify, LDT_endrun

    implicit none
! !ARGUMENTS: 
    integer, intent(in)        :: n 
    integer                    :: ftn
    character(len=*)           :: vname
    integer                    :: vlevel
    integer                    :: inc
    integer                    :: inr
    real                       :: value2d(inc, inr)
! 
! !DESCRIPTION: 
!   This routine reads a single variable from the LIS output file in the 
!   NetCDF format, assumed to be in the LDT grid.
!   Based on the external datamask, the routine also filters
!   each variable. 
!
!EOP
    real, allocatable :: value1d(:)
    integer :: unit_id
    integer :: varid
    integer :: iret
    integer :: i,c,r,gid
    character*50 :: units_tmp

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    if(LDT_rc%lis_wopt.eq.1) then 
       
       write(*,*) '[ERR] Reading in tile space format is not yet supported'
       write(*,*) '[ERR] program stopping ...'
       stop
       
    elseif(LDT_rc%lis_wopt.eq.2) then 
       
       iret = nf90_inq_varid(ftn, trim(vname), varid)
       call LDT_verify(iret, 'Error in nf90_inq_varid: '//trim(vname))
       
       iret = nf90_get_att(ftn,varid, "units", units_tmp)
       call LDT_verify(iret, 'Error in nf90_get_att:units for '//trim(vname))

       iret = nf90_get_var(ftn, varid, value2d,start = (/1,1,vlevel/),&
          count = (/inc,inr,1/))
       call LDT_verify(iret,'Error in nf90_get_var: '//trim(vname))
       
       
    elseif(LDT_rc%lis_wopt.eq.3) then 
       
       write(*,*) '[ERR] Reading in 1d grid space format is not yet supported'
       write(*,*) '[ERR] program stopping ...'
       stop
    endif
    
#endif
  end subroutine readLISSingleNetcdfVar_Inputgrid


!BOP
! !ROUTINE: LDT_writevar_netcdf
! \label{LDT_writevar_netcdf}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_netcdf
  subroutine LDT_writevar_netcdf(ftn,n, var,varid,dim1)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: ftn
    integer             :: varid
    real, intent(in)    :: var(LDT_rc%ntiles(n))
    integer, intent(in), optional :: dim1
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
    real :: var1(LDT_rc%ngrid(n))
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)
    integer :: gdeltas
    integer :: count1 ,c,r,m,gid,ntiles,ierr,i,t

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    if(LDT_rc%wopt.eq."1d tilespace") then !tiled output
       call LDT_gather_tiled_vector_output(n,gtmp1,var)
       
       if(LDT_masterproc) then 
          
          if(PRESENT(dim1)) then 
             iret = nf90_put_var(ftn,varid,gtmp1,(/1,dim1/),&
                  (/LDT_rc%glbntiles_red(n),1/))
             call LDT_verify(iret,'nf90_put_var failed in LDT_historyMod')
          else            
             iret = nf90_put_var(ftn,varid,gtmp1,(/1/),&
                  (/LDT_rc%glbntiles_red(n)/))
             call LDT_verify(iret,'nf90_put_var failed in LDT_historyMod')
          endif
          deallocate(gtmp1)
       endif

    elseif(LDT_rc%wopt.eq."2d gridspace") then 
       if(LDT_masterproc) then 
          allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
          allocate(gtmp1(LDT_rc%glbngrid(n)))
          gtmp = 0.0
          gtmp1 = 0.0
       else
          allocate(gtmp1(1))
          gtmp1 = 0.0
       endif
       var1 = 0 
       do i=1,LDT_rc%ntiles(n),LDT_rc%nensem(n)
          c = LDT_domain(n)%tile(i)%index
          do m=1,LDT_rc%nensem(n)
             t = i+m-1
             if ( var(t) == -9999.0 ) then
                var1(c) = -9999.0
             else
                var1(c) = var1(c) + &
                     var(t)*LDT_domain(n)%tile(t)%fgrd*&
                     LDT_domain(n)%tile(t)%pens
             endif
          enddo
       enddo
#if (defined SPMD)      
       gdeltas = LDT_gdeltas(n,LDT_localPet)
       call MPI_GATHERV(var1,gdeltas,&
            MPI_REAL,gtmp1,LDT_gdeltas(n,:),LDT_goffsets(n,:),&
            MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
       gtmp1 = var1
#endif
       if(LDT_masterproc) then 
          gtmp = LDT_rc%udef
          count1=1
          do l=1,LDT_npes
             do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
                do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                   gid = c+(r-1)*LDT_rc%gnc(n)
                   ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                   if(ntiles.ne.0) then                 
                      if(r.ge.LDT_nss_ind(n,l).and.&
                           r.le.LDT_nse_ind(n,l).and.&
                           c.ge.LDT_ews_ind(n,l).and.&
                           c.le.LDT_ewe_ind(n,l))then !points not in halo
                         gtmp(c,r) = gtmp1(count1)
                      endif
                      count1 = count1 + 1
                   endif
                enddo
             enddo
          enddo
          if(PRESENT(dim1)) then 
             iret = nf90_put_var(ftn,varid,gtmp,(/1,1,dim1/),&
                                 (/LDT_rc%gnc(n),LDT_rc%gnr(n),1/))
          else            
             iret = nf90_put_var(ftn,varid,gtmp,(/1,1/),&
                                 (/LDT_rc%gnc(n),LDT_rc%gnr(n)/))
          endif
          deallocate(gtmp)
       endif
       deallocate(gtmp1)
    endif
#endif
998 FORMAT(1X,A18,4E14.3)
999 FORMAT(1X,A18,4F14.3)
  end subroutine LDT_writevar_netcdf


!BOP
!
! !ROUTINE: LDT_gather_tiled_vector_output
! \label{LDT_gather_tiled_vector_output}
!
! !REVISION HISTORY:
!  30 Jan 2009:  Sujay Kumar; Initial code
! 
! !INTERFACE:
subroutine LDT_gather_tiled_vector_output(n, gtmp, var)
! !USES: 

! !ARGUMENTS: 

   implicit none

   integer                       :: n
   real, allocatable                :: gtmp(:)
   real, intent(in)             :: var(LDT_rc%ntiles(n))

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
    
    if(LDT_masterproc) then 
       allocate(gtmp(LDT_rc%glbntiles_red(n)))
       allocate(gtmp1(LDT_rc%glbntiles(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)      
    tdeltas = LDT_tdeltas(n,LDT_localPet)
    call MPI_GATHERV(var,tdeltas,&
         MPI_REAL,gtmp1,LDT_tdeltas(n,:),LDT_toffsets(n,:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
    gtmp1 = var
#endif
    if(LDT_masterproc) then 
       count1=1
       do l=1,LDT_npes
          do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
             do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                gid = c+(r-1)*LDT_rc%gnc(n)
                ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                stid = LDT_domain(n)%str_tind(gid)
                if(r.ge.LDT_nss_ind(n,l).and.&
                     r.le.LDT_nse_ind(n,l).and.&
                     c.ge.LDT_ews_ind(n,l).and.&
                     c.le.LDT_ewe_ind(n,l))then !points not in halo
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
 end subroutine LDT_gather_tiled_vector_output

end module LDT_historyMod
