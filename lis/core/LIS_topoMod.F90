!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module LIS_topoMod
!BOP
!
! !MODULE: LIS_topoMod
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
  use LIS_fileIOMod
  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_topo_init  ! initializes data structures and memory
  public :: LIS_diagnosetopography
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_topo
!EOP
  type, public :: topo_type_dec
     real, allocatable :: elevation(:,:,:)
     real, allocatable :: elevfgrd(:,:,:)
     real, allocatable :: slope(:,:,:)
     real, allocatable :: slopefgrd(:,:,:)
     real, allocatable :: aspect(:,:,:)
     real, allocatable :: aspectfgrd(:,:,:)
     real, allocatable :: curvature(:,:,:)
     real, allocatable :: curvfgrd(:,:,:)
  end type topo_type_dec
  
  type(topo_type_dec), allocatable :: LIS_topo(:)
!EOP

contains

!BOP
! 
! !ROUTINE: LIS_topo_init
! \label{LIS_topo_init}
! 
! !INTERFACE:
  subroutine LIS_topo_init()
! !USES:
    use ESMF
    use LIS_coreMod,  only : LIS_rc, LIS_config
    use LIS_fileIOMod, only : LIS_readDomainConfigSpecs
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
    integer :: ndoms

    ndoms = 0 
    
    allocate(LIS_topo(LIS_rc%nnest))

    do n=1,LIS_rc%nnest
       if(LIS_rc%useelevationmap(n).ne."none") then 
          call read_elevation(n)
       endif
       
       if(LIS_rc%useslopemap(n).ne."none") then 
          call read_slope(n)
       endif
       if(LIS_rc%useaspectmap(n).ne."none") then 
          call read_aspect(n)
       endif

       if(LIS_rc%usecurvaturemap(n).ne."none") then 
          call read_curvature(n)
       endif
    enddo
 
  end subroutine LIS_topo_init


!BOP
! 
! !ROUTINE: LIS_diagnoseTopography
! \label{LIS_diagnoseTopography}
! 
! !INTERFACE: 
  subroutine LIS_diagnoseTopography(n)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_histDataMod
! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n 

! !DESCRIPTION: 
!  This routine writes the LIS topography to history writer
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \end{description}
! 
!  The routines called are: 
!  \begin{description}
!  \item LIS\_diagnoseOutputVar \ref{LIS_diagnoseSurfaceOutputVar}
!  \end{description}
!EOP
    real, allocatable :: temp(:)
    integer :: t
    
    allocate(temp(LIS_rc%ntiles(n)))

    if(LIS_rc%useelevationmap(n).ne."none") then 
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then 
             temp(t) = LIS_domain(n)%tile(t)%elev
             call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ELEVATION,vlevel=1,&
                  value=temp(t),unit="m",direction="-")
          endif
       enddo
    endif
    if(LIS_rc%useslopemap(n).ne."none") then 
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then 
             temp(t) = LIS_domain(n)%tile(t)%slope
             call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SLOPE,vlevel=1,&
                  value=temp(t),unit="-",direction="-")
          endif
       enddo
    endif
    deallocate(temp)

  end subroutine LIS_diagnoseTopography


!BOP
!
! !ROUTINE: read_elevation
!  \label{read_elevation}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_elevation(n)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,        only : LIS_rc, LIS_localPet,&
         LIS_ews_ind, LIS_ewe_ind,&
         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
         LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
  
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
    logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       
       write(LIS_logunit,*)'[INFO] Reading elevation map from ' & 
            //trim(LIS_rc%paramfile(n))
           
       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_elevation')
       
       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_elevation')
       
       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_elevation')
       
       ios = nf90_inq_dimid(nid,"elevbins",ntypesId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_elevation')
       
       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_elevation')
       
       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_elevation')
       
       ios = nf90_inquire_dimension(nid,ntypesId, len=LIS_rc%nelevbands)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_elevation')

       allocate(LIS_topo(n)%elevation(LIS_rc%lnc(n),LIS_rc%lnr(n), &
            LIS_rc%nelevbands))
       allocate(LIS_topo(n)%elevfgrd(LIS_rc%lnc(n),LIS_rc%lnr(n), &
            LIS_rc%nelevbands))
       
       ios = nf90_inq_varid(nid,'ELEVATION',elevid)
       call LIS_verify(ios,'ELEVATION field not found in the LIS param file')
       
       ios = nf90_get_var(nid,elevid,LIS_topo(n)%elevation,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nelevbands/))
       call LIS_verify(ios,'Error in nf90_get_var in read_elevation')
       
       ios = nf90_inq_varid(nid,'ELEVFGRD',elevfgrdid)
       call LIS_verify(ios,'ELEVFGRD field not found in the LIS param file')
       
       ios = nf90_get_var(nid,elevfgrdid,LIS_topo(n)%elevfgrd,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nelevbands/))
       call LIS_verify(ios,'Error in nf90_get_var in read_elevation')
       
       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in read_elevation')
       
    else
       write(LIS_logunit,*) '[ERR] elevation map: ',LIS_rc%paramfile(n), &
            ' does not exist'
       write(LIS_logunit,*) '[ERR] program stopping ...'
       call LIS_endrun
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
  subroutine read_slope(n)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,        only : LIS_rc, LIS_localPet,&
         LIS_ews_ind, LIS_ewe_ind,&
         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
         LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
  
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
    real, allocatable    :: slope1(:,:,:)
    real, allocatable    :: slopefgrd1(:,:,:)
    logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       
       write(LIS_logunit,*)'[INFO] Reading slope map from '&
            //trim(LIS_rc%paramfile(n))
       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_slope')
       
       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_slope')
       
       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_slope')
       
       ios = nf90_inq_dimid(nid,"slopebins",ntypesId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_slope')
       
       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_slope')
       
       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_slope')
       
       ios = nf90_inquire_dimension(nid,ntypesId, len=LIS_rc%nslopebands)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_slope')

       allocate(LIS_topo(n)%slope(LIS_rc%lnc(n),LIS_rc%lnr(n), &
            LIS_rc%nslopebands))
       allocate(LIS_topo(n)%slopefgrd(LIS_rc%lnc(n),LIS_rc%lnr(n), &
            LIS_rc%nslopebands))
       
       ios = nf90_inq_varid(nid,'SLOPE',slopeid)
       call LIS_verify(ios,'SLOPE field not found in the LIS param file')
       
       ios = nf90_get_var(nid,slopeid,LIS_topo(n)%slope,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nslopebands/))
       call LIS_verify(ios,'Error in nf90_get_var in read_slope')
      
       ios = nf90_inq_varid(nid,'SLOPEFGRD',slopefgrdid)
       call LIS_verify(ios,'SLOPEFGRD field not found in the LIS param file')
       
       ios = nf90_get_var(nid,slopefgrdid,LIS_topo(n)%slopefgrd,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nslopebands/))
       call LIS_verify(ios,'Error in nf90_get_var in read_slope')
       
       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in read_slope')
       
    else
       write(LIS_logunit,*) '[ERR] slope map: ',LIS_rc%paramfile(n), &
            ' does not exist'
       write(LIS_logunit,*) '[ERR] program stopping ...'
       call LIS_endrun
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
  subroutine read_aspect(n)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,        only : LIS_rc, LIS_localPet,&
         LIS_ews_ind, LIS_ewe_ind,&
         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
         LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
  
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
    real, allocatable    :: aspect1(:,:,:)
    real, allocatable    :: aspectfgrd1(:,:,:)
    logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       
       write(LIS_logunit,*)'[INFO] Reading aspect map from '&
            //trim(LIS_rc%paramfile(n))

       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_aspect')
       
       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_aspect')
       
       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_aspect')
       
       ios = nf90_inq_dimid(nid,"aspectbins",ntypesId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_aspect')
       
       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_aspect')
       
       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_aspect')
       
       ios = nf90_inquire_dimension(nid,ntypesId, len=LIS_rc%naspectbands)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_aspect')

       allocate(LIS_topo(n)%aspect(LIS_rc%lnc(n),LIS_rc%lnr(n), &
            LIS_rc%naspectbands))
       allocate(LIS_topo(n)%aspectfgrd(LIS_rc%lnc(n),LIS_rc%lnr(n), &
            LIS_rc%naspectbands))

       ios = nf90_inq_varid(nid,'ASPECT',aspectid)
       call LIS_verify(ios,'ASPECT field not found in the LIS param file')
       
       ios = nf90_get_var(nid,aspectid,LIS_topo(n)%aspect,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%naspectbands/))
       call LIS_verify(ios,'Error in nf90_get_var in read_aspect')
       
       ios = nf90_inq_varid(nid,'ASPECTFGRD',aspectfgrdid)
       call LIS_verify(ios,'ASPECTFGRD field not found in the LIS param file')
       
       ios = nf90_get_var(nid,aspectfgrdid,LIS_topo(n)%aspectfgrd,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%naspectbands/))
       call LIS_verify(ios,'Error in nf90_get_var in read_aspect')
       
       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in read_aspect')
       
    else
       write(LIS_logunit,*) '[ERR] aspect map: ',LIS_rc%paramfile(n), &
            ' does not exist'
       write(LIS_logunit,*) '[ERR] program stopping ...'
       call LIS_endrun
    endif
#endif

  end subroutine read_aspect

!BOP
!
! !ROUTINE: read_curvature
!  \label{read_curvature}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  Mar 2022: K.R. Arsenault;  Added curvature
!
! !INTERFACE:
  subroutine read_curvature(n)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,        only : LIS_rc, LIS_localPet,&
         LIS_ews_ind, LIS_ewe_ind,&
         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
         LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n

! !DESCRIPTION:
!  This subroutine reads the curvature data
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
    integer          :: ios,nid,ntypesId,ncId, nrId
    integer          :: curvid, curvfgrdId
    integer          :: nc,nr,c,r
    logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then

       write(LIS_logunit,*)'[INFO] Reading curvature map from '&
            //trim(LIS_rc%paramfile(n))

       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_curvature')

       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_curvature')

       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_curvature')

!       ios = nf90_inq_dimid(nid,"curvbins",ntypesId)
!       call LIS_verify(ios,'Error in nf90_inq_dimid in read_curvature')

       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_curvature')

       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_curvature')

!       ios = nf90_inquire_dimension(nid,ntypesId, len=LIS_rc%curvbands)
!       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_curvature')

       allocate(LIS_topo(n)%curvature(LIS_rc%lnc(n),LIS_rc%lnr(n), &
                1))
       allocate(LIS_topo(n)%curvfgrd(LIS_rc%lnc(n),LIS_rc%lnr(n), &
                1))

       ios = nf90_inq_varid(nid,'CURVATURE',curvid)
       call LIS_verify(ios,'CURVATURE field not found in the LIS param file')


       ios = nf90_get_var(nid,curvid,LIS_topo(n)%curvature,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),1/))
       call LIS_verify(ios,'Error in nf90_get_var in read_aspect')

       do r = 1,LIS_rc%lnr(n)
          do c = 1,LIS_rc%lnc(n)
             if( LIS_topo(n)%curvature(c,r,1) .ne. LIS_rc%udef ) then
                LIS_topo(n)%curvfgrd(c,r,1) = 1.0
             else
                LIS_topo(n)%curvfgrd(c,r,1) = 0.0
             endif
          enddo
       enddo

    else
       write(LIS_logunit,*) '[ERR] curvature map in: ',LIS_rc%paramfile(n), &
            ' does not exist.'
       call LIS_endrun
    endif
#endif

  end subroutine read_curvature

end module LIS_topoMod
