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
module LIS_LMLCMod
!BOP
!
! !MODULE: LIS_LMLCMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read various sources of
!  landmask and landcover data. 
! 
!  \subsubsection{Overview}
!  This module provides routines for reading and modifying landmask 
!  and landcover data. 
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
  public :: LIS_LMLC_init          ! initializes data structures and memory
  public :: LIS_diagnoselandmask   ! diagnoses the LIS landmask for 
                                   ! history output
  public :: LIS_diagnoselandcover  ! diagnoses the LIS landcover data for
                                   !  history output

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_LMLC  ! derived datatype that stores the landmask and landcover data
!EOP
  type, public :: lmlc_type_dec
     real, allocatable :: dommask(:,:)
     real, allocatable :: landmask(:,:)
     real, allocatable :: glandmask(:,:)
     real, allocatable :: landcover(:,:,:)
     real, allocatable :: surfacetype(:,:,:)
  end type lmlc_type_dec
  
  type(lmlc_type_dec), allocatable :: LIS_LMLC(:)
contains

!BOP
! 
! !ROUTINE: LIS_LMLC_init
! \label{LIS_LMLC_init}
! 
! !INTERFACE:
  subroutine LIS_LMLC_init()
! !USES:
    use ESMF
    use LIS_coreMod,  only : LIS_rc, LIS_config, LIS_domain
    use LIS_logMod,   only : LIS_verify, LIS_logunit
    use LIS_fileIOMod, only : LIS_readDomainConfigSpecs
! 
! !DESCRIPTION:
!
! Reads the configurable options related to the landmask and 
! landcover datasets (filenames, geographical extent information)
! and allocates memory for data structures for reading 
! landmask and landcover datasets
!
!EOP
    implicit none
    integer :: n, i
    integer :: rc, iret

    allocate(LIS_LMLC(LIS_rc%nnest))

    do n=1,LIS_rc%nnest
       write(LIS_logunit,*) '[INFO] Reading mask from '& 
            //trim(LIS_rc%paramfile(n))

       allocate(LIS_LMLC(n)%dommask(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       call LIS_read_param(n,"DOMAINMASK",LIS_LMLC(n)%dommask,iret)

       !If the DOMAINMASK field doesn't exist, then use the LANDMASK
       if(iret.ne.0) then
          call LIS_read_param(n,"LANDMASK",LIS_LMLC(n)%dommask)
       endif

       allocate(LIS_LMLC(n)%landmask(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       call LIS_read_param(n,"LANDMASK",LIS_LMLC(n)%landmask)

       write(LIS_logunit,*) '[INFO] Saving global land mask...'
       allocate(LIS_LMLC(n)%glandmask(LIS_rc%gnc(n),LIS_rc%gnr(n)))

       call LIS_read_gparam(n,"LANDMASK",LIS_LMLC(n)%glandmask)

       call read_landcover(n)
       call read_surfacetype(n)
    enddo

!------------------------------------------------------------------------
! Fix SURFACE_MAXT out of bounds problems
!------------------------------------------------------------------------
    if(LIS_rc%surface_maxt.gt.LIS_rc%nsurfacetypes) &
         LIS_rc%surface_maxt=LIS_rc%nsurfacetypes
    if(LIS_rc%surface_maxt.lt.1) LIS_rc%surface_maxt=1

    write(LIS_logunit,*) '[INFO] Finished reading landmask and landcover data'
  end subroutine LIS_LMLC_init

!BOP
! 
! !ROUTINE: LIS_diagnoselandmask
! \label{LIS_diagnoselandmask}
! 
! !INTERFACE: 
  subroutine LIS_diagnoselandmask(n)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_histDataMod
! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n 

! !DESCRIPTION: 
!  This routine diagnoses the LIS landmask to be used later
!  for history output 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_diagnoseOutputVar] (\ref{LIS_diagnoseSurfaceOutputVar})  \newline
!   This routine maps landmask data to the history writing routines
!  \end{description}
! 
!EOP
    integer :: t
    real, allocatable    :: temp(:)
    
    allocate(temp(LIS_rc%ntiles(n)))

    temp = LIS_rc%udef
    do t=1,LIS_rc%ntiles(n)
       if(LIS_domain(n)%tile(t)%index.ne.-1) then 
          temp(t) = 1.0
       endif
       call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LANDMASK,vlevel=1,&
            value=temp(t),unit="-",direction="-")
    enddo
    deallocate(temp)

  end subroutine LIS_diagnoselandmask

!BOP
! 
! !ROUTINE: LIS_diagnoseLandcover
! \label{LIS_diagnoseLandcover}
! 
! !INTERFACE: 
  subroutine LIS_diagnoseLandcover(n)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_histDataMod    
! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n 

! !DESCRIPTION: 
!  This routine diagnoses the LIS landcover to be later
!  used for history output
!
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_diagnoseOutputVar] (\ref{LIS_diagnoseSurfaceOutputVar})  \newline
!   This routine maps landcover data to the history writing routines
!  \end{description}
!EOP
    integer :: t
    real, allocatable    :: temp(:)
    
    allocate(temp(LIS_rc%ntiles(n)))

    temp = LIS_rc%udef
    do t=1,LIS_rc%ntiles(n)
       if(LIS_domain(n)%tile(t)%index.ne.-1) then 
          temp(t) = real(LIS_domain(n)%tile(t)%vegt)
       endif
       call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LANDCOVER,vlevel=1,&
            value=temp(t),unit="-",direction="-")
    enddo
    deallocate(temp)
  end subroutine LIS_diagnoseLandcover
!BOP
!
! !ROUTINE: read_landcover
!  \label{read_landcover}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_landcover(n)
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
!  This subroutine reads the landcover data
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
  integer          :: ios,nid,vegId, lcid,ncId, nrId
  integer          :: nc,nr
  integer          :: c,r,t
  real, allocatable    :: lc(:,:,:)
  logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then 

     write(LIS_logunit,*)'[INFO] Reading landcover map from '&
          //trim(LIS_rc%paramfile(n))
     ios = nf90_open(path=LIS_rc%paramfile(n),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in read_landcover')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_landcover')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_landcover')

     ios = nf90_inq_dimid(nid,"sfctypes",vegId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_landcover')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_landcover')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_landcover')

     ios = nf90_inquire_dimension(nid,vegId, len=LIS_rc%nsurfacetypes)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_landcover')

     allocate(LIS_LMLC(n)%landcover(LIS_rc%lnc(n),LIS_rc%lnr(n), &
          LIS_rc%nsurfacetypes))

     ios = nf90_inq_varid(nid,'LANDCOVER',lcid)
     call LIS_verify(ios,'LANDCOVER field not found in the LIS param file')

     ios = nf90_get_var(nid,lcid,LIS_LMLC(n)%landcover,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nsurfacetypes/))
     call LIS_verify(ios,'Error in nf90_get_var in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'BARESOILCLASS', LIS_rc%bareclass)
     call LIS_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'URBANCLASS', LIS_rc%urbanclass)
     call LIS_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'SNOWCLASS', LIS_rc%snowclass)
     call LIS_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'WATERCLASS', LIS_rc%waterclass)
     call LIS_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'WETLANDCLASS', LIS_rc%wetlandclass)
     call LIS_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'GLACIERCLASS', LIS_rc%glacierclass)
     call LIS_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'CROPCLASS', LIS_rc%cropclass)
!<kluge -- jim testing>
! Temporarily disable error check for cropclass.  cropclass was added
! to support NoahMP 4.0.1.  No other lsm uses this right now, so everyone's
! runs will fail until they rerun LDT.
!
! Discuss with everyone whether this should be a required attribute and
! give everyone time to reprocess their domain and parameter files.
!     call LIS_verify(ios,'Error in nf90_get_att in read_landcover')
!</kluge -- jim testing>

     ios = nf90_get_att(nid, NF90_GLOBAL, 'NUMVEGTYPES', LIS_rc%nvegtypes)
     call LIS_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'LANDCOVER_SCHEME', LIS_rc%lcscheme)
     call LIS_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in read_landcover')

  else
     write(LIS_logunit,*) '[ERR] landcover map: ',LIS_rc%paramfile(n), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif
#endif

end subroutine read_landcover

!BOP
!
! !ROUTINE: read_surfacetype
!  \label{read_surfacetype}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_surfacetype(n)
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
!  This subroutine reads the surfacetype data
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
  integer          :: ios,nid,vegId, lcid,ncId, nrId
  integer          :: nc,nr
  integer          :: c,r,t
  real, allocatable    :: lc(:,:,:)
  logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then 

     write(LIS_logunit,*)'[INFO] Reading surfacetype map from '&
          //trim(LIS_rc%paramfile(n))
     ios = nf90_open(path=LIS_rc%paramfile(n),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in read_surfacetype')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_surfacetype')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_surfacetype')

     ios = nf90_inq_dimid(nid,"sfctypes",vegId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_surfacetype')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_surfacetype')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_surfacetype')

     ios = nf90_inquire_dimension(nid,vegId, len=LIS_rc%nsurfacetypes)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_surfacetype')

     allocate(LIS_LMLC(n)%surfacetype(LIS_rc%lnc(n),LIS_rc%lnr(n), LIS_rc%nsurfacetypes))

!     allocate(lc(LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nsurfacetypes))

     ios = nf90_inq_varid(nid,'SURFACETYPE',lcid)
     call LIS_verify(ios,'SURFACETYPE field not found in the LIS param file')

     ios = nf90_get_var(nid,lcid,LIS_LMLC(n)%surfacetype,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nsurfacetypes/))
     call LIS_verify(ios,'Error in nf90_get_var in read_surfacetype')

     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in read_surfacetype')

  else
     write(LIS_logunit,*) '[ERR] surfacetype map: ',LIS_rc%paramfile(n), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif
#endif

end subroutine read_surfacetype

end module LIS_LMLCMod
