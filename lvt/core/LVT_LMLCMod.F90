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
!BOP
! 
! !MODULE: LVT_LMLCMod
! \label(LVT_LMLCMod)
!
! !INTERFACE:
module LVT_LMLCMod
! 
! !USES:   

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  The code in this file implements routines to read various sources of
!  landmask and landcover data. 
! 
!  \subsubsection{Overview}
!  This module provides routines for reading and modifying landmask 
!  and landcover data. 
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  18 Jul 2008: Sujay Kumar; Initial implementation
! 
!EOP
!BOP
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LVT_LMLC_init            ! initializes data structures and memory

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LVT_LMLC  ! derived datatype that stores the landmask and landcover data
!EOP
  type, public :: lmlc_type_dec
     real, allocatable :: dommask(:,:)
     real, allocatable :: landmask(:,:)
     real, allocatable :: landcover(:,:,:)
     real, allocatable :: surfacetype(:,:,:)
  end type lmlc_type_dec
  
  type(lmlc_type_dec) :: LVT_LMLC(3)
contains

!BOP
! 
! !ROUTINE: LVT_LMLC_init
! \label{LVT_LMLC_init}
!
! !INTERFACE:
  subroutine LVT_LMLC_init(k, paramfile)
! 
! !USES:
    use ESMF
    use LVT_coreMod
    use LVT_logMod
    use LVT_fileIOMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! Allocates memory for data structures for reading 
! landmask and landcover datasets
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    integer           :: k
    character(len=*)  :: paramfile
    logical           :: file_exists
    integer           :: ios, iret
    integer           :: nid,ncId,nrId,nc,nr,paramId,dommaskId
    real, allocatable :: param(:,:)
    integer           :: rc

    allocate(LVT_LMLC(k)%dommask(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr))
    allocate(LVT_LMLC(k)%landmask(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr))
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    inquire(file=paramfile, exist=file_exists)
    if(file_exists) then 

       write(LVT_logunit,*)'[INFO] Reading LANDMASK map ...'
       ios = nf90_open(path=(paramfile),&
            mode=NF90_NOWRITE,ncid=nid)
       call LVT_verify(ios,'Error in nf90_open in readparam_real_2d')
       
       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in readparam_real_2d')
       
       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LVT_verify(ios,'Error in nf90_inq_dimid in readparam_real_2d')
       
       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in readparam_real_2d')
       
       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LVT_verify(ios,'Error in nf90_inquire_dimension in readparam_real_2d')
       
       iret = nf90_inq_varid(nid,"DOMAINMASK",dommaskid)

       ios = nf90_inq_varid(nid,"LANDMASK",paramid)
       call LVT_verify(ios,'LANDMASK field not found in the LVT param file')
       
       allocate(param(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr))

       ios = nf90_get_var(nid,paramid,param)
       call LVT_verify(ios,'Error in nf90_get_var in readparam_real_2d')
       
       LVT_LMLC(k)%landmask(:,:) =  param(:,:)

!If the DOMAINMASK field doesn't exist, then use the LANDMASK
       if(iret.eq.0) then 
          ios = nf90_get_var(nid,dommaskid,param)
          call LVT_verify(ios,'Error in nf90_get_var in readparam_real_2d')
       endif
       LVT_LMLC(k)%dommask(:,:) =  param(:,:)

       ios = nf90_close(nid)
       call LVT_verify(ios,'Error in nf90_close in readparam_real_2d')       

       deallocate(param)

       write(LVT_logunit,*)'[INFO] Successfully read LANDMASK map '
    else
       write(LVT_logunit,*) '[ERR] LANDMASK map: ',&
            (paramfile), '[ERR] does not exist'
       call LVT_endrun
    endif

#endif


    call read_landcover(k,paramfile)
    call read_surfacetype(k,paramfile)

    if(LVT_LIS_rc(k)%surface_maxt.gt.LVT_LIS_rc(k)%nsurfacetypes) &
         LVT_LIS_rc(k)%surface_maxt=LVT_LIS_rc(k)%nsurfacetypes
    if(LVT_LIS_rc(k)%surface_maxt.lt.1) LVT_LIS_rc(k)%surface_maxt=1

    write(LVT_logunit,*) '[INFO] Finished reading landmask and landcover data'
  end subroutine LVT_LMLC_init

!BOP
!
! !ROUTINE: read_landcover
!  \label{read_landcover}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_landcover(k,paramfile)
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

  inquire(file=paramfile, exist=file_exists)
  if(file_exists) then 

     write(LVT_logunit,*)'[INFO] Reading landcover map...'
     ios = nf90_open(path=paramfile,&
          mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios,'Error in nf90_open in read_landcover')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_landcover')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_landcover')

     ios = nf90_inq_dimid(nid,"sfctypes",vegId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_landcover')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_landcover')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_landcover')

     ios = nf90_inquire_dimension(nid,vegId, len=LVT_LIS_rc(k)%nsurfacetypes)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_landcover')

     allocate(LVT_LMLC(k)%landcover(&
          LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr, &
          LVT_LIS_rc(k)%nsurfacetypes))

     allocate(lc(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,&
          LVT_LIS_rc(k)%nsurfacetypes))

     ios = nf90_inq_varid(nid,'LANDCOVER',lcid)
     call LVT_verify(ios,'LANDCOVER field not found in the LVT param file')

     ios = nf90_get_var(nid,lcid,lc)
     call LVT_verify(ios,'Error in nf90_get_var in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'BARESOILCLASS', LVT_LIS_rc(k)%bareclass)
     call LVT_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'URBANCLASS', LVT_LIS_rc(k)%urbanclass)
     call LVT_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'SNOWCLASS', LVT_LIS_rc(k)%snowclass)
     call LVT_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'WATERCLASS', LVT_LIS_rc(k)%waterclass)
     call LVT_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'WETLANDCLASS', LVT_LIS_rc(k)%wetlandclass)
     call LVT_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_get_att(nid, NF90_GLOBAL, 'GLACIERCLASS', LVT_LIS_rc(k)%glacierclass)
     call LVT_verify(ios,'Error in nf90_get_att in read_landcover')

     ios = nf90_close(nid)
     call LVT_verify(ios,'Error in nf90_close in read_landcover')

     LVT_LMLC(k)%landcover = lc

     deallocate(lc)
  else
     write(LVT_logunit,*) '[ERR] landcover map: ',paramfile, ' does not exist'
     call LVT_endrun
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
subroutine read_surfacetype(k,paramfile)
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

  inquire(file=paramfile, exist=file_exists)
  if(file_exists) then 

     write(LVT_logunit,*)'[INFO] Reading surfacetype map...'
     ios = nf90_open(path=paramfile,&
          mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios,'Error in nf90_open in read_surfacetype')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_surfacetype')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_surfacetype')

     ios = nf90_inq_dimid(nid,"sfctypes",vegId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_surfacetype')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_surfacetype')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_surfacetype')

     ios = nf90_inquire_dimension(nid,vegId, len=LVT_LIS_rc(k)%nsurfacetypes)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_surfacetype')

!to be removed later: 
     LVT_LIS_rc(k)%nvegtypes = LVT_LIS_rc(k)%nsurfacetypes

     allocate(LVT_LMLC(k)%surfacetype(&
          LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,&
          LVT_LIS_rc(k)%nsurfacetypes))

     allocate(lc(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,&
          LVT_LIS_rc(k)%nsurfacetypes))

     ios = nf90_inq_varid(nid,'SURFACETYPE',lcid)
     call LVT_verify(ios,'SURFACETYPE field not found in the LVT param file')

     ios = nf90_get_var(nid,lcid,lc)
     call LVT_verify(ios,'Error in nf90_get_var in read_surfacetype')

     ios = nf90_close(nid)
     call LVT_verify(ios,'Error in nf90_close in read_surfacetype')

     LVT_LMLC(k)%surfacetype = lc
     
     deallocate(lc)
  else
     write(LVT_logunit,*) '[ERR] surfacetype map: ',paramfile, ' does not exist'
     call LVT_endrun
  endif
#endif

end subroutine read_surfacetype

end module LVT_LMLCMod
