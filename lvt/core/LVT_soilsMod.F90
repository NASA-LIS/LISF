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
module LVT_soilsMod
!BOP
!
! !MODULE: LVT_soilsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read various sources of
!  soil parameter data. 
! 
!  \subsubsection{Overview}
!  This module provides routines for reading and manipulating various 
!  parameters related to soil properties. The following list of 
!  soil parameters are currently supported. 
!  \begin{description}
!   \item[sand, silt, clay fractions]
!   \item[soil color data]
!   \item[soil texture data]
!   \item[soil porosity data]
!   \item[hydraulic conductivity data]
!   \item[thermal conductivity data]
!   \item[b parameter data]
!   \item[quartz data]
!  \end{description}
!
! !REVISION HISTORY:
!
!  21 Oct 2005: Sujay Kumar; Initial implementation
!  3  Apr 2012: Sujay Kumar; Switched to the use of LPT based parameter file
!
  use LVT_fileIOMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LVT_soils_init  ! initializes data structures and read soil data
  public :: LVT_soils_finalize !cleanup allocated structures
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LVT_soils      !data structure containing soils data
!EOP

  type, public :: soils_type_dec 
     real, allocatable :: sand(:,:,:)
     real, allocatable :: clay(:,:,:)
     real, allocatable :: silt(:,:,:)
     real, allocatable :: soilffgrd(:,:,:)
     real, allocatable :: color(:,:)     
     real, allocatable :: texture(:,:,:)
     real, allocatable :: porosity(:,:,:)
     real, allocatable :: psisat(:,:)
     real, allocatable :: ksat(:,:)
     real, allocatable :: bexp(:,:)
     real, allocatable :: quartz(:,:)
  end type soils_type_dec
  
  type(soils_type_dec) :: LVT_soils(3)

contains

!BOP
! 
! !ROUTINE: LVT_soils_init
! \label{LVT_soils_init}
! 
! !INTERFACE:
  subroutine LVT_soils_init(k,paramfile)
! !USES:
    use ESMF
    use LVT_coreMod
    use LVT_logMod
! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! soils datasets
!
!  Reads the soils data based on the choice of options specified 
!  in the lis configuration. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readsoiltexture](\ref{readsoiltexture}) \newline
!    invokes the generic method in the registry to read
!    the soil texture data
!   \item[readsand](\ref{readsand}) \newline
!    invokes the generic method in the registry to read
!    the sand fraction data
!   \item[readclay](\ref{readclay}) \newline
!    invokes the generic method in the registry to read
!    the clay fraction data
!   \item[readsilt](\ref{readsilt}) \newline
!    invokes the generic method in the registry to read
!    the silt fraction data
!   \item[readcolor](\ref{readcolor}) \newline
!    invokes the generic method in the registry to read
!    the soil color data
!   \item[readporosity](\ref{readporosity}) \newline
!    invokes the generic method in the registry to read
!    the porosity data
!   \item[readpsisat](\ref{readpsisat}) \newline
!    invokes the generic method in the registry to read
!    the saturated matric potential data
!   \item[readksat](\ref{readksat}) \newline
!    invokes the generic method in the registry to read
!    the saturated hydraulic conductivity data
!   \item[readbexp](\ref{readbexp}) \newline
!    invokes the generic method in the registry to read
!    the b parameter data
!   \item[readquartz](\ref{readquartz}) \newline
!    invokes the generic method in the registry to read
!    the quartz data
!  \end{description}
!
!EOP
    implicit none
    integer          :: k
    character(len=*) :: paramfile

    integer :: i
    integer :: rc

    if(LVT_LIS_rc(k)%usetexturemap.ne."none".and.&
         LVT_LIS_rc(k)%usesoilfractionmap.ne."none") then 
       
       write(LVT_logunit,*) '[ERR] Please select either the soil texture or the soil '
       write(LVT_logunit,*) '[ERR] fraction dataset. Both should not be enabled '
       write(LVT_logunit,*) '[ERR] simultaneously ...'
       call LVT_endrun()
       
    endif
    
    if(LVT_LIS_rc(k)%usetexturemap.ne."none") then           
       call read_soiltexture(k,paramfile)
    else
       if(LVT_LIS_rc(k)%usesoilfractionmap.ne."none") then 
          call read_soilfraction(k,paramfile)
       endif
    endif
  end subroutine LVT_soils_init

  
!BOP
! 
! !ROUTINE: LVT_soils_finalize
! \label{LVT_soils_finalize}
! 
! !INTERFACE:
  subroutine LVT_soils_finalize()
! !USES:
    use LVT_coreMod
! !DESCRIPTION:
!
! deallocates memory for all datastructures used for reading
! soils datasets. This method is typically called after the 
! information is translated to the LSM model tiles. 
! 
!
!EOP    
    implicit none
    integer :: n,k
    
    do k=1,2
       
       if(LVT_LIS_rc(k)%usetexturemap.ne."none") then !read soil texture
          deallocate(LVT_soils(k)%texture)
       else
          if(LVT_LIS_rc(k)%usesoilfractionmap.ne."none") then 
             deallocate(LVT_soils(k)%sand)
             deallocate(LVT_soils(k)%clay)
             deallocate(LVT_soils(k)%silt)
          endif
          deallocate(LVT_soils(k)%texture)
       endif
    enddo

  end subroutine LVT_soils_finalize
  
  

!BOP
!
! !ROUTINE: read_soiltexture
!  \label{read_soiltexture}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_soiltexture(k,paramfile)
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
!  This subroutine reads the soiltexture data
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
  integer          :: ios,nid,ntypesId, txtid,ncId, nrId
  integer          :: nc,nr
  integer          :: c,r,t
  real, allocatable    :: txt(:,:,:)
  logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=paramfile, exist=file_exists)
  if(file_exists) then 

     write(LVT_logunit,*)'[INFO] Reading soiltexture map...'
     ios = nf90_open(path=paramfile,&
          mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios,'Error in nf90_open in read_soiltexture')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_soiltexture')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_soiltexture')

     ios = nf90_inq_dimid(nid,"soiltypes",ntypesId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_soiltexture')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_soiltexture')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_soiltexture')

     ios = nf90_inquire_dimension(nid,ntypesId, len=LVT_LIS_rc(k)%nsoiltypes)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_soiltexture')

     allocate(LVT_soils(k)%texture(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr, &
          LVT_LIS_rc(k)%nsoiltypes))

     allocate(txt(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,LVT_LIS_rc(k)%nsoiltypes))

     ios = nf90_inq_varid(nid,'TEXTURE',txtid)
     call LVT_verify(ios,'TEXTURE field not found in the LVT param file')

     ios = nf90_get_var(nid,txtid,txt)
     call LVT_verify(ios,'Error in nf90_get_var in read_soiltexture')

     ios = nf90_close(nid)
     call LVT_verify(ios,'Error in nf90_close in read_soiltexture')

     LVT_soils(k)%texture = txt
     deallocate(txt)
  else
     write(LVT_logunit,*) '[ERR] soiltexture map: ',paramfile, &
          ' does not exist'
     call LVT_endrun
  endif
#endif

end subroutine read_soiltexture


!BOP
!
! !ROUTINE: read_soilfraction
!  \label{read_soilfraction}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_soilfraction(k,paramfile)
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
!  This subroutine reads the soilfraction data
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
  integer          :: ios,nid,ntypesId, ncId, nrId
  integer          :: sandId, clayId,soilfgrdId
  integer          :: nc,nr
  integer          :: c,r,t,kk
  real, allocatable    :: sand(:,:,:),clay(:,:,:),soilfgrd(:,:,:)
  logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=paramfile, exist=file_exists)
  if(file_exists) then 

     write(LVT_logunit,*)'[INFO] Reading sandfraction map...'
     ios = nf90_open(path=paramfile,&
          mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios,'Error in nf90_open in read_soilfraction')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_soilfraction')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_soilfraction')

     ios = nf90_inq_dimid(nid,"soilfracbins",ntypesId)
     call LVT_verify(ios,'Error in nf90_inq_dimid in read_soilfraction')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_soilfraction')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_soilfraction')

     ios = nf90_inquire_dimension(nid,ntypesId, len=LVT_LIS_rc(k)%nsoilfbands)
     call LVT_verify(ios,'Error in nf90_inquire_dimension in read_soilfraction')

     allocate(LVT_soils(k)%sand(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr, &
          LVT_LIS_rc(k)%nsoilfbands))
     allocate(LVT_soils(k)%clay(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr, &
          LVT_LIS_rc(k)%nsoilfbands))

     allocate(sand(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,LVT_LIS_rc(k)%nsoilfbands))
     allocate(clay(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,LVT_LIS_rc(k)%nsoilfbands))

     ios = nf90_inq_varid(nid,'SAND',sandid)
     call LVT_verify(ios,'SAND field not found in the LVT param file')

     ios = nf90_get_var(nid,sandid,sand)
     call LVT_verify(ios,'Error in nf90_get_var in read_soilfraction')

     LVT_soils(k)%sand = sand
     deallocate(sand)

     ios = nf90_inq_varid(nid,'CLAY',clayid)
     call LVT_verify(ios,'CLAY field not found in the LVT param file')

     ios = nf90_get_var(nid,clayid,clay)
     call LVT_verify(ios,'Error in nf90_get_var in read_soilfraction')

     LVT_soils(k)%clay = clay

     deallocate(clay)

     allocate(soilfgrd(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,LVT_LIS_rc(k)%nsoilfbands))

     allocate(LVT_soils(k)%soilffgrd(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,&
          LVT_LIS_rc(k)%nsoilfbands))

     ios = nf90_inq_varid(nid,'SOILSFGRD',soilfgrdid)
     call LVT_verify(ios,'SOILSFGRD field not found in the LVT param file')

     ios = nf90_get_var(nid,soilfgrdid,soilfgrd)
     call LVT_verify(ios,'Error in nf90_get_var in read_soilfraction')

     LVT_soils(k)%soilffgrd = soilfgrd
     
     deallocate(soilfgrd)

     ios = nf90_close(nid)
     call LVT_verify(ios,'Error in nf90_close in read_soilfraction')

     allocate(LVT_soils(k)%silt(LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr,&
          LVT_LIS_rc(k)%nsoilfbands))
     
     do kk =1,LVT_LIS_rc(k)%nsoilfbands
        do r=1,LVT_LIS_rc(k)%gnr
           do c=1,LVT_LIS_rc(k)%gnc
              if(LVT_soils(k)%sand(c,r,kk).ne.-9999.0.and.&
                   LVT_soils(k)%clay(c,r,kk).ne.-9999.0) then 
                 LVT_soils(k)%silt(c,r,kk) = 1.0 - &
                      LVT_soils(k)%sand(c,r,kk) - & 
                      LVT_soils(k)%clay(c,r,kk)
              else
                 LVT_soils(k)%silt(c,r,kk) = -9999.0
              endif
           enddo
        enddo
     enddo
     
  else
     write(LVT_logunit,*) '[ERR] sand/clay fraction map: ',paramfile, &
          ' does not exist'
     call LVT_endrun
  endif
#endif

end subroutine read_soilfraction


end module LVT_soilsMod
