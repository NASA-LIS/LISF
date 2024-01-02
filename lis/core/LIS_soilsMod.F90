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
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
module LIS_soilsMod
!BOP
!
! !MODULE: LIS_soilsMod
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
  use LIS_fileIOMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_soils_init  ! initializes data structures and read soil data
  public :: LIS_diagnosesoils ! routine that maps soils data to the history writer
  public :: LIS_soils_finalize !cleanup allocated structures
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_soils      !data structure containing soils data
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
  
  type(soils_type_dec), allocatable :: LIS_soils(:)

contains

!BOP
! 
! !ROUTINE: LIS_soils_init
! \label{LIS_soils_init}
! 
! !INTERFACE:
  subroutine LIS_soils_init()
! !USES:
    use ESMF
    use LIS_coreMod,  only : LIS_rc
    use LIS_logMod,   only : LIS_logunit, LIS_endrun
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
!   \item[read\_soiltexture](\ref{read_soiltexture}) \newline
!    invokes the generic method in the registry to read
!    the soil texture data
!   \item[read\_soilfraction](\ref{read_soilfraction}) \newline
!    invokes the generic method in the registry to read
!    the soil fraction data
!   \item[read\_porosity](\ref{read_porosity}) \newline
!    invokes the generic method in the registry to read
!    the soil porosity data
!   \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
!    reads the specified soil property from the 
!    LIS input parameters file (e.g.; lis\_input.d01.nc)
!  \end{description}
!
!EOP
    implicit none
    integer :: n, i
    integer :: rc
    logical :: SubModelIsCrocus    
    

    TRACE_ENTER("soils_init")
    allocate(LIS_soils(LIS_rc%nnest))

    do n=1,LIS_rc%nnest
       if(LIS_rc%usetexturemap(n).ne."none".and.&
            LIS_rc%usesoilfractionmap(n).ne."none") then 

            write(LIS_logunit,*) '[WARN] Both soil texture and soil fraction dataset are selected.'
            write(LIS_logunit,*) '[WARN] Both should generally not be enabled simultaneously.'
            write(LIS_logunit,*) '[WARN] For now, the soil texture will be used; the soil fraction'
            write(LIS_logunit,*) '[WARN] will be ignored. However, users should double-check'
            write(LIS_logunit,*) '[WARN] the output of "Soiltype:" via the MODEL OUTPUT TBL.'
       endif

       if(LIS_rc%usetexturemap(n).ne."none") then           
          call read_soiltexture(n)
       endif
       if(LIS_rc%usesoilfractionmap(n).ne."none") then 
          call read_soilfraction(n)
       endif
       if(LIS_rc%usesoilcolormap(n).ne."none") then 
          allocate(LIS_soils(n)%color(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          call LIS_read_param(n,"SOILCOLOR",LIS_soils(n)%color)
       endif
       if(LIS_rc%useporositymap(n).ne."none") then 
          call read_porosity(n)
       endif
       if(LIS_rc%usepsisatmap(n).ne."none") then 
          allocate(LIS_soils(n)%psisat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          call LIS_read_param(n,"PSISAT",LIS_soils(n)%psisat)
       endif
       if(LIS_rc%useksatmap(n).ne."none") then 
          allocate(LIS_soils(n)%ksat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          call LIS_read_param(n,"KSAT",LIS_soils(n)%ksat)
       endif
       if(LIS_rc%usebexpmap(n).ne."none") then 
          allocate(LIS_soils(n)%bexp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          call LIS_read_param(n,"BEXP",LIS_soils(n)%bexp)
       endif
       if(LIS_rc%usequartzmap(n).ne."none") then 
          allocate(LIS_soils(n)%quartz(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          call LIS_read_param(n,"QUARTZ",LIS_soils(n)%quartz)
       endif
    enddo
    TRACE_EXIT("soils_init")
  end subroutine LIS_soils_init

!BOP
! 
! !ROUTINE: LIS_diagnosesoils
! \label{LIS_diagnosesoils}
! 
! !INTERFACE: 
  subroutine LIS_diagnosesoils(n)
! !USES: 
    use LIS_coreMod,     only : LIS_rc, LIS_domain
    use LIS_histDataMod, only : LIS_diagnoseSurfaceOutputVar, &
                                LIS_MOC_SOILTYPE,             &
                                LIS_MOC_SANDFRAC,             &
                                LIS_MOC_CLAYFRAC,             &
                                LIS_MOC_SILTFRAC,             &
                                LIS_MOC_SOILCOLOR,            &
                                LIS_MOC_POROSITY
#ifdef ESMF_TRACE
    use ESMF
#endif
! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n 

! !DESCRIPTION: 
!  This routine writes the LIS soils to the LIS 
!  history writer
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[ftn] file unit number to be used \newline
!  \end{description}
! 
!  The routines called are: 
!  \begin{description}
!  \item[LIS\_diagnoseOutputVar] (\ref{LIS_diagnoseSurfaceOutputVar})  \newline
!   This routine maps a variable to the history writing routines
!  \end{description}
!EOP
    integer :: t,k
    real, pointer :: temp(:)

    TRACE_ENTER("soils_diag")
    allocate(temp(LIS_rc%ntiles(n)))

    if(LIS_rc%usetexturemap(n).ne."none") then 
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then 
             temp(t) = real(LIS_domain(n)%tile(t)%soilt)
          endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILTYPE,vlevel=1,&
                                           value=temp(t),unit="-",direction="-")
       enddo
    endif
    if(LIS_rc%usesoilfractionmap(n).ne."none") then 
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then 
             temp(t) = LIS_domain(n)%tile(t)%sand
          endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SANDFRAC,vlevel=1,&
               value=temp(t),unit="-",direction="-")
       enddo
       
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then 
             temp(t) =LIS_domain(n)%tile(t)%clay
          endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CLAYFRAC,vlevel=1,&
               value=temp(t),unit="-",direction="-")
       enddo
       
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then 
             temp(t) = LIS_domain(n)%tile(t)%silt
          endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SILTFRAC,vlevel=1,&
               value=temp(t),unit="-",direction="-")
       enddo
    endif

    if(LIS_rc%usesoilcolormap(n).ne."none") then 
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then 
             temp(t) = LIS_soils(n)%color(LIS_domain(n)%tile(t)%col,&
                  LIS_domain(n)%tile(t)%row)
          endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILCOLOR,vlevel=1,&
                                           value=temp(t),unit="-",direction="-")
       enddo

    endif

    if(LIS_rc%useporositymap(n).ne."none") then 
       do k=1,3
          temp = LIS_rc%udef
          do t=1,LIS_rc%ntiles(n)
             if(LIS_domain(n)%tile(t)%index.ne.-1) then 
                temp(t) = LIS_soils(n)%porosity(LIS_domain(n)%tile(t)%col,&
                     LIS_domain(n)%tile(t)%row,k)
             endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POROSITY,vlevel=k,&
                                           value=temp(t),unit="-",direction="-")
          enddo
       enddo
    endif
    deallocate(temp)
    TRACE_EXIT("soils_diag")

  end subroutine LIS_diagnosesoils
  
!BOP
! 
! !ROUTINE: LIS_soils_finalize
! \label{LIS_soils_finalize}
! 
! !INTERFACE:
  subroutine LIS_soils_finalize()
! !USES:
    use LIS_coreMod, only : LIS_rc
! !DESCRIPTION:
!
! deallocates memory for all datastructures used for reading
! soils datasets. This method is typically called after the 
! information is translated to the LSM model tiles. 
! 
!
!EOP    
    implicit none
    integer :: n
    do n=1,LIS_rc%nnest
       if(LIS_rc%usetexturemap(n).ne."none") then !read soil texture
          deallocate(LIS_soils(n)%texture)
       endif
       if(LIS_rc%usesoilfractionmap(n).ne."none") then 
          deallocate(LIS_soils(n)%sand)
          deallocate(LIS_soils(n)%clay)
          deallocate(LIS_soils(n)%silt)
       endif
       if(LIS_rc%usesoilcolormap(n).ne."none") then 
          deallocate(LIS_soils(n)%color)
       endif
       if(LIS_rc%useporositymap(n).ne."none") then 
          deallocate(LIS_soils(n)%porosity)
       endif
       if(LIS_rc%usepsisatmap(n).ne."none") then 
          deallocate(LIS_soils(n)%psisat)
       endif
       if(LIS_rc%useksatmap(n).ne."none") then 
          deallocate(LIS_soils(n)%ksat)
       endif
       if(LIS_rc%usebexpmap(n).ne."none") then 
          deallocate(LIS_soils(n)%bexp)
       endif
       if(LIS_rc%usequartzmap(n).ne."none") then 
          deallocate(LIS_soils(n)%quartz)
       endif
    enddo
  end subroutine LIS_soils_finalize


!BOP
!
! !ROUTINE: read_porosity
!  \label{read_porosity}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_porosity(n)
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
!  This subroutine reads the porosity data 
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[localmask]
!    porosity for the region of interest
!   \end{description}
!
!EOP      

  integer :: ios1
  integer :: ios,nid,porosityid,ncId, nrId
  integer :: nc,nr,t
  real, allocatable :: porosity(:,:,:)
  logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then 

     write(LIS_logunit,*)'[ERR] Reading soil porosity map from '&
          //trim(LIS_rc%paramfile(n))

     ios = nf90_open(path=LIS_rc%paramfile(n),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in read_porosity')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_porosity')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_porosity')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_porosity')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_porosity')

     ios = nf90_inq_varid(nid,'POROSITY',porosityid)
     call LIS_verify(ios,'POROSITY field not found in the LIS param file')

     allocate(LIS_soils(n)%porosity(LIS_rc%lnc(n),LIS_rc%lnr(n),3))

     ios = nf90_get_var(nid,porosityid,LIS_soils(n)%porosity,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),3/))
     call LIS_verify(ios,'Error in nf90_get_var in read_porosity')
     
     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in read_porosity')

  else
     write(LIS_logunit,*) '[ERR] porosity map: ',LIS_rc%paramfile(n), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif
#endif
end subroutine read_porosity

!BOP
!
! !ROUTINE: read_soiltexture
!  \label{read_soiltexture}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_soiltexture(n)
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

  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then 

     write(LIS_logunit,*)'[INFO] Reading soiltexture map from '&
          //trim(LIS_rc%paramfile(n))
     ios = nf90_open(path=LIS_rc%paramfile(n),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in read_soiltexture')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_soiltexture')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_soiltexture')

     ios = nf90_inq_dimid(nid,"soiltypes",ntypesId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_soiltexture')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_soiltexture')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_soiltexture')

     ios = nf90_inquire_dimension(nid,ntypesId, len=LIS_rc%nsoiltypes)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_soiltexture')

     allocate(LIS_soils(n)%texture(LIS_rc%lnc(n),LIS_rc%lnr(n), &
          LIS_rc%nsoiltypes))

     ios = nf90_inq_varid(nid,'TEXTURE',txtid)
     call LIS_verify(ios,'TEXTURE field not found in the LIS param file')

     ios = nf90_get_var(nid,txtid,LIS_soils(n)%texture, &
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nsoiltypes/))
     call LIS_verify(ios,'Error in nf90_get_var in read_soiltexture')

     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in read_soiltexture')

  else
     write(LIS_logunit,*) '[ERR] soiltexture map: ',LIS_rc%paramfile(n), &
          ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
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
subroutine read_soilfraction(n)
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
  integer          :: c,r,t,k
  real, allocatable    :: sand(:,:,:),clay(:,:,:),soilfgrd(:,:,:)
  logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then 

     write(LIS_logunit,*)'[INFO] Reading sandfraction map from '&
          //trim(LIS_rc%paramfile(n))

     ios = nf90_open(path=LIS_rc%paramfile(n),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in read_soilfraction')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_soilfraction')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_soilfraction')

     ios = nf90_inq_dimid(nid,"soilfracbins",ntypesId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_soilfraction')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_soilfraction')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_soilfraction')

     ios = nf90_inquire_dimension(nid,ntypesId, len=LIS_rc%nsoilfbands)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_soilfraction')

     allocate(LIS_soils(n)%sand(LIS_rc%lnc(n),LIS_rc%lnr(n), &
          LIS_rc%nsoilfbands))
     allocate(LIS_soils(n)%clay(LIS_rc%lnc(n),LIS_rc%lnr(n), &
          LIS_rc%nsoilfbands))

     ios = nf90_inq_varid(nid,'SAND',sandid)
     call LIS_verify(ios,'SAND field not found in the LIS param file')

     ios = nf90_get_var(nid,sandid,LIS_soils(n)%sand,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nsoilfbands/))
     call LIS_verify(ios,'Error in nf90_get_var in read_soilfraction')

     ios = nf90_inq_varid(nid,'CLAY',clayid)
     call LIS_verify(ios,'CLAY field not found in the LIS param file')

     ios = nf90_get_var(nid,clayid,LIS_soils(n)%clay,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nsoilfbands/))
     call LIS_verify(ios,'Error in nf90_get_var in read_soilfraction')

     allocate(LIS_soils(n)%soilffgrd(LIS_rc%lnc(n),LIS_rc%lnr(n),&
          LIS_rc%nsoilfbands))

     ios = nf90_inq_varid(nid,'SOILSFGRD',soilfgrdid)
     call LIS_verify(ios,'SOILSFGRD field not found in the LIS param file')

     ios = nf90_get_var(nid,soilfgrdid,LIS_soils(n)%soilffgrd,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1),1/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nsoilfbands/))
     call LIS_verify(ios,'Error in nf90_get_var in read_soilfraction')

     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in read_soilfraction')

     allocate(LIS_soils(n)%silt(LIS_rc%lnc(n),LIS_rc%lnr(n),&
          LIS_rc%nsoilfbands))
     
     do k =1,LIS_rc%nsoilfbands
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_soils(n)%sand(c,r,k).ne.-9999.0.and.&
                   LIS_soils(n)%clay(c,r,k).ne.-9999.0) then 
                 LIS_soils(n)%silt(c,r,k) = 1.0 - &
                      LIS_soils(n)%sand(c,r,k) - & 
                      LIS_soils(n)%clay(c,r,k)
              else
                 LIS_soils(n)%silt(c,r,k) = -9999.0
              endif
           enddo
        enddo
     enddo
     
  else
     write(LIS_logunit,*) '[ERR] sand/clay fraction map: ',LIS_rc%paramfile(n), &
          ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif
#endif

end subroutine read_soilfraction


end module LIS_soilsMod
