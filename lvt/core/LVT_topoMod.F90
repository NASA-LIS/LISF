!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
module LVT_topoMod
!BOP
!
! !MODULE: LVT_topoMod
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
!  3  Apr 2012: Sujay Kumar; Switched to the use of LPT based parameter file
!
  use LVT_fileIOMod
  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LVT_topo_init  ! initializes data structures and memory
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LVT_topo
!EOP
  type, public :: topo_type_dec
     real, allocatable :: elevation(:,:,:)
     real, allocatable :: elevfgrd(:,:,:)
     real, allocatable :: slope(:,:,:)
     real, allocatable :: slopefgrd(:,:,:)
     real, allocatable :: aspect(:,:,:)
     real, allocatable :: aspectfgrd(:,:,:)
     real, allocatable :: curvature(:,:,:)
  end type topo_type_dec
  
  type(topo_type_dec) :: LVT_topo(3)
!EOP

contains

!BOP
! 
! !ROUTINE: LVT_topo_init
! \label{LVT_topo_init}
! 
! !INTERFACE:
  subroutine LVT_topo_init(k,paramfile)
! !USES:
    use ESMF
    use LVT_coreMod

    implicit none
    integer           :: k
    character(len=*)  :: paramfile
    
! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! topo datasets
!
!EOP

    integer :: n, i
    integer :: rc
    integer :: ndoms

    ndoms = 0 
    
    if(LVT_LIS_rc(k)%useelevationmap.ne."none") then 
       call read_elevation(k,paramfile)
    endif
    
    if(LVT_LIS_rc(k)%useslopemap.ne."none") then 
       call read_slope(k,paramfile)
    endif
    if(LVT_LIS_rc(k)%useaspectmap.ne."none") then 
       call read_aspect(k,paramfile)
    endif
    
  end subroutine LVT_topo_init


!BOP
!
! !ROUTINE: read_elevation
!  \label{read_elevation}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_elevation(k,paramfile)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LVT_coreMod
    use LVT_logMod
    use LVT_fileIOMod
    
    implicit none

    integer           :: k
    character(len=*)  :: paramfile

! !ARGUMENTS: 

  
! !DESCRIPTION:
!  This subroutine reads the elevation data
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[locallc]
!    landlc for the region of interest
!   \end{description}
!
!EOP      

    integer          :: ios1
    integer          :: ios,nid,ntypesId, elevid,ncId, nrId
    integer          :: elevfgrdId
    integer          :: nc,nr
    integer          :: c,r,t
    real, allocatable    :: elev(:,:,:)
    real, allocatable    :: elevfgrd(:,:,:)
    real, allocatable    :: elevfgrd1(:,:,:)
    logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    
    inquire(file=paramfile, exist=file_exists)
    if(file_exists) then 
       
       write(LVT_logunit,*)'[INFO] Reading elevation map...'
       ios = nf90_open(path=paramfile,&
            mode=NF90_NOWRITE,ncid=nid)
       call LVT_verify(ios,'Error in nf90_open in read_elevation')
       
       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in read_elevation')
       
       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in read_elevation')
       
       ios = nf90_inq_dimid(nid,"elevbins",ntypesId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in read_elevation')
       
       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in read_elevation')
       
       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in read_elevation')
       
       ios = nf90_inquire_dimension(nid,ntypesId, len=LVT_LIS_rc(k)%nelevbands)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in read_elevation')

! EMK
!       allocate(LVT_topo(k)%elevation(LVT_rc%gnc,LVT_rc%gnr,&
!            LVT_LIS_rc(k)%nelevbands))
       allocate(LVT_topo(k)%elevation(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,&
            LVT_LIS_rc(k)%nelevbands))

       allocate(LVT_topo(k)%elevfgrd(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,&
            LVT_LIS_rc(k)%nelevbands))
       
       allocate(elev(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,LVT_LIS_rc(k)%nelevbands))
       
       ios = nf90_inq_varid(nid,'ELEVATION',elevid)
       call LVT_verify(ios,'ELEVATION field not found in the LVT param file')
       
       ios = nf90_get_var(nid,elevid,elev)
       call LVT_verify(ios,'Error in nf90_get_var in read_elevation')
       
       LVT_topo(k)%elevation = elev
       deallocate(elev)
       
       allocate(elevfgrd(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,LVT_LIS_rc(k)%nelevbands))

       ios = nf90_inq_varid(nid,'ELEVFGRD',elevfgrdid)
       call LVT_verify(ios,'ELEVFGRD field not found in the LVT param file')
       
       ios = nf90_get_var(nid,elevfgrdid,elevfgrd)
       call LVT_verify(ios,'Error in nf90_get_var in read_elevation')
       
       ios = nf90_close(nid)
       call LVT_verify(ios,'Error in nf90_close in read_elevation')
       
       LVT_topo(k)%elevfgrd = elevfgrd

       deallocate(elevfgrd)

    else
       write(LVT_logunit,*) '[ERR] elevation map: ',paramfile, &
            ' does not exist'
       call LVT_endrun
    endif
#endif

  end subroutine read_elevation

!BOP
!
! !ROUTINE: read_slope
!  \label{read_slope}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_slope(k,paramfile)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LVT_coreMod
    use LVT_logMod
    use LVT_fileIOMod
    
    implicit none

    integer          :: k
    character(len=*) :: paramfile
    
! !ARGUMENTS: 

  
! !DESCRIPTION:
!  This subroutine reads the slope data
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[locallc]
!    landlc for the region of interest
!   \end{description}
!
!EOP      

    integer          :: ios1
    integer          :: ios,nid,ntypesId, slopeid,ncId, nrId
    integer          :: slopefgrdId
    integer          :: nc,nr
    integer          :: c,r,t
    real, allocatable    :: slope(:,:,:)
    real, allocatable    :: slopefgrd(:,:,:)
    logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    
    inquire(file=paramfile, exist=file_exists)
    if(file_exists) then 
       
       write(LVT_logunit,*)'[INFO] Reading slope map...'
       ios = nf90_open(path=paramfile,&
            mode=NF90_NOWRITE,ncid=nid)
       call LVT_verify(ios,'Error in nf90_open in read_slope')
       
       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in read_slope')
       
       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in read_slope')
       
       ios = nf90_inq_dimid(nid,"slopebins",ntypesId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in read_slope')
       
       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in read_slope')
       
       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in read_slope')
       
       ios = nf90_inquire_dimension(nid,ntypesId, len=LVT_LIS_rc(k)%nslopebands)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in read_slope')

       allocate(LVT_topo(k)%slopefgrd(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,&
            LVT_LIS_rc(k)%nslopebands))
       
       allocate(slope(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,LVT_LIS_rc(k)%nslopebands))
              
       ios = nf90_inq_varid(nid,'SLOPE',slopeid)
       call LVT_verify(ios,'SLOPE field not found in the LVT param file')
       
       ios = nf90_get_var(nid,slopeid,slope)
       call LVT_verify(ios,'Error in nf90_get_var in read_slope')
       
       LVT_topo(k)%slope = slope
       deallocate(slope)

       allocate(slopefgrd(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,LVT_LIS_rc(k)%nslopebands))

       ios = nf90_inq_varid(nid,'SLOPEFGRD',slopefgrdid)
       call LVT_verify(ios,'SLOPEFGRD field not found in the LVT param file')
       
       ios = nf90_get_var(nid,slopefgrdid,slopefgrd)
       call LVT_verify(ios,'Error in nf90_get_var in read_slope')
       
       ios = nf90_close(nid)
       call LVT_verify(ios,'Error in nf90_close in read_slope')
       
       LVT_topo(k)%slopefgrd = slopefgrd
       deallocate(slopefgrd)
       
    else
       write(LVT_logunit,*) '[ERR] slope map: ',paramfile, &
            ' does not exist'
       call LVT_endrun
    endif
#endif

  end subroutine read_slope

!BOP
!
! !ROUTINE: read_aspect
!  \label{read_aspect}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_aspect(k,paramfile)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LVT_coreMod
    use LVT_logMod
    use LVT_fileIOMod
    
    implicit none

    integer          :: k
    character(len=*) :: paramfile
! !ARGUMENTS: 

  
! !DESCRIPTION:
!  This subroutine reads the aspect data
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[locallc]
!    landlc for the region of interest
!   \end{description}
!
!EOP      

    integer          :: ios1
    integer          :: ios,nid,ntypesId, aspectid,ncId, nrId
    integer          :: aspectfgrdId
    integer          :: nc,nr
    integer          :: c,r,t
    real, allocatable    :: aspect(:,:,:)
    real, allocatable    :: aspectfgrd(:,:,:)
    logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    
    inquire(file=paramfile, exist=file_exists)
    if(file_exists) then 
       
       write(LVT_logunit,*)'[INFO] Reading aspect map...'
       ios = nf90_open(path=paramfile,&
            mode=NF90_NOWRITE,ncid=nid)
       call LVT_verify(ios,'Error in nf90_open in read_aspect')
       
       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in read_aspect')
       
       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in read_aspect')
       
       ios = nf90_inq_dimid(nid,"aspectbins",ntypesId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in read_aspect')
       
       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in read_aspect')
       
       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in read_aspect')
       
       ios = nf90_inquire_dimension(nid,ntypesId, len=LVT_LIS_rc(k)%naspectbands)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in read_aspect')

       allocate(LVT_topo(k)%aspect(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,&
            LVT_LIS_rc(k)%naspectbands))
       allocate(LVT_topo(k)%aspectfgrd(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,&
            LVT_LIS_rc(k)%naspectbands))
       
       allocate(aspect(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,LVT_LIS_rc(k)%naspectbands))
       
       ios = nf90_inq_varid(nid,'ASPECT',aspectid)
       call LVT_verify(ios,'ASPECT field not found in the LVT param file')
       
       ios = nf90_get_var(nid,aspectid,aspect)
       call LVT_verify(ios,'Error in nf90_get_var in read_aspect')
       
       LVT_topo(k)%aspect = aspect
       deallocate(aspect)

       allocate(aspectfgrd(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,LVT_LIS_rc(k)%naspectbands))

       ios = nf90_inq_varid(nid,'ASPECTFGRD',aspectfgrdid)
       call LVT_verify(ios,'ASPECTFGRD field not found in the LVT param file')
       
       ios = nf90_get_var(nid,aspectfgrdid,aspectfgrd)
       call LVT_verify(ios,'Error in nf90_get_var in read_aspect')
       
       ios = nf90_close(nid)
       call LVT_verify(ios,'Error in nf90_close in read_aspect')
       
       LVT_topo(k)%aspectfgrd = aspectfgrd

       deallocate(aspectfgrd)

    else
       write(LVT_logunit,*) '[ERR] aspect map: ',paramfile, &
            ' does not exist'
       call LVT_endrun
    endif
#endif

  end subroutine read_aspect


end module LVT_topoMod
