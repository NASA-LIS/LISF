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
module LDT_SurfaceTypeMod
!BOP
!
! !MODULE: LDT_SurfaceTypeMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read various sources of
!   surface type data.
! 
!  \subsubsection{Overview}
!  This module provides routines for reading and assigning surface type
!   data.
!
! !REVISION HISTORY:
!
!  18 Jul 2008: Sujay Kumar; Initial implementation
!  28 Apr 2013: KR Arsenault; Expanded surface type specific code
!
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_domainMod, only: isSurfaceTypeSelected

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_surfacetype_init         ! initializes data structures and memory
  public :: LDT_assign_lcsfctype
  public :: LDT_assign_lakesfctype
  public :: LDT_assign_openwatersfctype
  public :: LDT_assign_glaciersfctype
  public :: LDT_surfacetype_writeHeader
  public :: LDT_surfacetype_writeData
  
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
!EOP

!BOP 
! 
! !ROUTINE: LDT_surfacetype_writeHeader 
! \label{LDT_surfacetype_writeHeader}
! 
! !INTERFACE:
  interface LDT_surfacetype_writeHeader
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure surfacetype_writeHeader_LIS
     module procedure surfacetype_writeHeader_LISHydro
! 
! !DESCRIPTION:
! This interface provides routines for writing NETCDF header both 
! in LIS preprocessing requirements as well as LISHydro(WRFHydro) 
! preprocessing requiremetns. A dummy argument call "flagX" was added 
! to overload the LISHydro procedue.
!EOP 
  end interface
contains

!BOP
! 
! !ROUTINE: LDT_surfacetype_init
! \label{LDT_surfacetype_init}
! 
! !INTERFACE:
  subroutine LDT_surfacetype_init()

! !USES:
   use ESMF
   use LDT_coreMod,  only : LDT_rc, LDT_config, LDT_domain
   use LDT_logMod,   only : LDT_verify, LDT_logunit
!   use LDT_fileIOMod,only : LDT_readDomainConfigSpecs
!   use LDT_paramOptCheckMod, only: LDT_maskSfctypeChecks
   use LDT_domainMod, only: isSurfaceTypeSelected
! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! surface type datasets
!
!EOP
   implicit none
   integer :: rc
   integer :: n
   integer :: i
   integer :: num_lctypes
   integer :: num_laketypes
   integer :: num_glaciertypes
   integer :: num_wetlandtypes
   integer :: num_openwatertypes

   character(20) :: lsmname,lakename,wetlandname,glaciername
   character(20) :: openwatername
   character(2)  :: char_numtypes

! ____________________________________________________________

  LDT_rc%inc_water_pts = .false.

  num_lctypes = 0
  num_laketypes = 0
  num_glaciertypes = 0
  num_wetlandtypes = 0
  num_openwatertypes = 0
  lsmname  = ""
  lakename = ""
  glaciername = ""
  wetlandname = ""
  openwatername = ""

  write(LDT_logunit,*) "[INFO]  - - - - - - -  SURFACE TYPES  - - - - - - - - "
  write(LDT_logunit,*) "[INFO] Initializing and Summing Surface Types " 
  write(LDT_logunit,*) "[INFO] ... Surface types selected: " 
  write(char_numtypes,'(i2.0)') LDT_rc%nsf_model_types
  write(LDT_logunit,'(a6,'//char_numtypes//'(a10))') " ... ",&
       (LDT_rc%sf_model_type_name_select(i),i=1,LDT_rc%nsf_model_types)


!- Allocate surface type file arrays:
   do n = 1, LDT_rc%nnest

   !- Landcover types:
      if( isSurfaceTypeSelected(1) ) then
        num_lctypes = LDT_LSMparam_struc(n)%landcover%num_bins
      endif

   !- Lake tiles:
      if( isSurfaceTypeSelected(2) ) then
        if( LDT_LSMparam_struc(n)%lakecover%selectOpt > 0 ) then
           num_laketypes = LDT_LSMparam_struc(n)%lakecover%num_bins 
        endif
        if( .NOT. isSurfaceTypeSelected(5) ) then 
           num_openwatertypes = 1
        endif
        LDT_rc%inc_water_pts = .true.
      endif

   !- Glacier tiles:
      if( isSurfaceTypeSelected(3) ) then
         num_glaciertypes = LDT_LSMparam_struc(n)%glaciermask%num_bins 
      end if

   !- Wetland tiles:
      if( isSurfaceTypeSelected(4) ) then
!         num_wetlandtypes = LDT_LSMparam_struc(n)%wetland%num_bins 
      end if

   !- Openwater type:
      if( isSurfaceTypeSelected(5) ) then
!         num_openwatertypes = LDT_LSMparam_struc(n)%openwater%num_bins 
         num_openwatertypes = 1
         LDT_rc%inc_water_pts = .true.
      end if

   !- Total surface types:
      LDT_LSMparam_struc(n)%sfctype%num_bins = num_lctypes      + &
                                               num_laketypes    + &
                                               num_glaciertypes + &
                                               num_wetlandtypes + &
                                               num_openwatertypes

!      print *, 'sfc types: ',LDT_LSMparam_struc(n)%sfctype%num_bins

   !- Set surface type source name (for output field):
!      write(*,*) (LDT_rc%sf_model_type_name_select(i),i=1, LDT_rc%nsf_model_types )

      if( isSurfaceTypeSelected(1) ) lsmname = trim(LDT_rc%lsm)
      if( isSurfaceTypeSelected(2) ) lakename = "+"//trim(LDT_rc%lakemodel)
!      if( isSurfaceTypeSelected(3) ) glaciername = "+"//trim(LDT_rc%glaciermodel)
!      if( isSurfaceTypeSelected(4) ) wetlandname = "+"//trim(LDT_rc%wetlandmodel)
      if( isSurfaceTypeSelected(5) ) openwatername = "+Openwater"

      LDT_LSMparam_struc(n)%sfctype%source = trim(lsmname)//&
                                             trim(lakename)//&
                                             trim(glaciername)//&
                                             trim(wetlandname)//&
                                             trim(openwatername)

      LDT_LSMparam_struc(n)%sfctype%short_name    = "SURFACETYPE"
      LDT_LSMparam_struc(n)%sfctype%standard_name = "Surface type"
      LDT_LSMparam_struc(n)%sfctype%vlevels       =  &
                   LDT_LSMparam_struc(n)%sfctype%num_bins 
      LDT_LSMparam_struc(n)%sfctype%selectOpt     = 1
      LDT_LSMparam_struc(n)%sfctype%num_times     = 1 
      LDT_LSMparam_struc(n)%sfctype%units         = "-"

   !- Allocate surface type file arrays:
      allocate(LDT_LSMparam_struc(n)%sfctype%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%sfctype%num_bins))

      LDT_LSMparam_struc(n)%sfctype%value = 0.   ! Initialize as "openwater" or "undefined"

   enddo

   write(LDT_logunit,*) "[INFO] Finished initializing surface types"

 end subroutine LDT_surfacetype_init


! -----
!BOP
!
! !ROUTINE: LDT_assign_lcsfctype
!  \label{LDT_assign_lcsfctype}
!
! !REVISION HISTORY:
!  20 Apr  2013: KR Arsenault;  Separated surface type field
!                                assignment into own subroutine
! !INTERFACE:
 subroutine LDT_assign_lcsfctype(n, total_types, nt,  &
                           sfctype, localmask, lcfgrd )
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_endrun
  use LDT_domainMod,only : isSurfaceTypeSelected

  implicit none

! !ARGUMENTS: 
  integer, intent(in)     :: n             ! Nest
  integer, intent(in)     :: total_types   ! Total surface types
  integer, intent(in)     :: nt            ! Landcover types
  real,    intent(in)     :: lcfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),nt)
  real,    intent(inout)  :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real,    intent(out)    :: sfctype(LDT_rc%lnc(n),LDT_rc%lnr(n),total_types)

! !DESCRIPTION:
!  This subroutine assigns surface types based on land/water mask 
!   and surface types selected.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[nt]
!    number of landcover classes or types
!   \item[lcfgrd]
!    landcover fractions of gridcell
!   \item[localmask]
!    landmask for the region of interest
!   \item[sfctype]
!    surface type for the region of interest
!   \end{description}
!EOP      
  integer :: c, r, t, l, i
! ___________________________

  write(LDT_logunit,fmt=*) '[INFO] SURFACETYPE -- Assigning landcover surface types'

   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         do t = 1, nt  ! landcover loop
            if( lcfgrd(c,r,t) > 0. ) then
               sfctype(c,r,t) = float(LDT_rc%lsm_index)
            endif
         end do
      end do
   end do

  end subroutine LDT_assign_lcsfctype
! ___________________________________________________________


!BOP
!
! !ROUTINE: LDT_assign_lakesfctype
!  \label{LDT_assign_lakesfctype}
!
! !REVISION HISTORY:
!  20 Apr  2013: KR Arsenault;  Separated surface type field
!                                assignment into own subroutine
! !INTERFACE:
 subroutine LDT_assign_lakesfctype(n, lc_types, total_types, sfctype, &
                       localmask, lcfgrd, lkfgrd )

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_endrun
  use LDT_domainMod,only : isSurfaceTypeSelected

  implicit none

! !ARGUMENTS: 
  integer,  intent(in)    :: n             ! Nest
  integer,  intent(in)    :: lc_types      ! Total landcover types 
  integer,  intent(in)    :: total_types   ! Total surface types
  real,     intent(in)    :: lkfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),total_types)
  real,     intent(inout) :: lcfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),total_types)
  real,     intent(inout) :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real,     intent(out)   :: sfctype(LDT_rc%lnc(n),LDT_rc%lnr(n),total_types)

! !DESCRIPTION:
!  This subroutine assigns surface types based on land/water mask 
!   and surface types selected.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[total_types]
!    number of total surface types
!   \item[lcfgrd]
!    landcover fractions of gridcell
!   \item[lkfgrd]
!    lake fractions of gridcell
!   \item[localmask]
!    landmask for the region of interest
!   \item[sfctype]
!    surface type for the region of interest
!   \end{description}
!EOP      
  integer :: c, r, t, l, i
  integer :: lake_tile

!  real    :: maxlc
!  real    :: swaptemp
! ___________________________

  write(LDT_logunit,fmt=*) "[INFO] SURFACETYPE -- Assigning lake surface types"

   lake_tile = lc_types + 1

   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)

      !- Landcover water fraction present:
         if( lcfgrd(c,r,LDT_rc%waterclass) > 0. ) then

         !- If Lake is present: 
         ! Note:  Current set-up here works only for 1-lake tile ...
         !        TO BE UPDATED FOR MULTI-LAKE TILES ...
           if( lkfgrd(c,r,lake_tile) > 0. ) then

            ! Set lake surface type tile to lake index = 2
              sfctype(c,r,lake_tile) = float(LDT_rc%lake_index)
              lcfgrd(c,r,lake_tile)  = lcfgrd(c,r,LDT_rc%waterclass) 
            ! Set original landcover water fraction at same point to 0.
              sfctype(c,r,LDT_rc%waterclass) = 0.
              lcfgrd(c,r,LDT_rc%waterclass)  = 0   ! 1-lake tile only

! -- OTHER OPTION: Force lake to be dominant fraction in landcover array ...
!                 maxlc = maxval(lcfgrd(c,r,:))
!                 do i = 1, total_types
!                    if( lcfgrd(c,r,i) == maxlc ) then
!                       swaptemp = lcfgrd(c,r,i) 
!                       lcfgrd(c,r,i) = lcfgrd(c,r,l)
!                       lcfgrd(c,r,l) = swaptemp
!                    endif
!                 end do
! -- 
        !- If Lake is NOT present:
           elseif( lkfgrd(c,r,lake_tile) == 0. ) then

            ! Set openwater surface type tile to index = 5
              sfctype(c,r,total_types) = float(LDT_rc%openwater_index)
              lcfgrd(c,r,total_types)  = lcfgrd(c,r,LDT_rc%waterclass) 

            ! Set original landcover water fraction at same point to 0.
              sfctype(c,r,LDT_rc%waterclass) = 0.
              lcfgrd(c,r,LDT_rc%waterclass)  = 0   

           end if

         else

         end if     ! end land cover water fraction check
      end do
   end do

   write(LDT_logunit,*) "[INFO] Setting lake fractions to landcover water fraction ..."

  end subroutine LDT_assign_lakesfctype


!BOP
!
! !ROUTINE: LDT_assign_glaciersfctype
!  \label{LDT_assign_glaciersfctype}
!
! !REVISION HISTORY:
!  31 Mar  2018: Sujay Kumar;  Separated surface type field
!                                assignment into own subroutine
! !INTERFACE:
 subroutine LDT_assign_glaciersfctype(n, lc_types, total_types, sfctype, &
                       localmask, lcfgrd, glfgrd )

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_endrun
  use LDT_domainMod,only : isSurfaceTypeSelected

  implicit none

! !ARGUMENTS: 
  integer    :: n             ! Nest
  integer    :: lc_types      ! Total landcover types 
  integer    :: total_types   ! Total surface types
  real       :: sfctype(LDT_rc%lnc(n),LDT_rc%lnr(n),total_types)
  real       :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real       :: lcfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),total_types)
  real       :: glfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),total_types)



! !DESCRIPTION:
!  This subroutine assigns surface types based on land/water mask 
!   and surface types selected.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[total_types]
!    number of total surface types
!   \item[lcfgrd]
!    landcover fractions of gridcell
!   \item[glfgrd]
!    glacier fractions of gridcell
!   \item[localmask]
!    landmask for the region of interest
!   \item[sfctype]
!    surface type for the region of interest
!   \end{description}
!EOP      
  integer :: c, r, t, l, i
  integer :: glacier_tile

!  real    :: maxlc
!  real    :: swaptemp
! ___________________________

  write(LDT_logunit,fmt=*) "[INFO] SURFACETYPE -- Assigning glacier surface types"
  
  glacier_tile = lc_types + 1
  
  if(LDT_LSMparam_struc(n)%glaciermask%source.eq."GLIMS") then 
     
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           
           !- If Glacier is present: 
           ! Note:  Current set-up here works only for 1-glacier tile ...
           !        TO BE UPDATED FOR MULTI-GLACIER TILES ...
           if( glfgrd(c,r,glacier_tile) > 0. ) then
              
#if 0 
              ! Set glacier surface type tile to glacier index = 3
              sfctype(c,r,glacier_tile) = float(LDT_rc%glacier_index)
              lcfgrd(c,r,glacier_tile)  = lcfgrd(c,r,LDT_rc%glacierclass) 
              ! Set original landcover glacier fraction at same point to 0.
              sfctype(c,r,LDT_rc%glacierclass) = 0.
              lcfgrd(c,r,LDT_rc%glacierclass)  = 0   
#endif
              ! Set glacier surface type tile to glacier index = 3
              sfctype(c,r,glacier_tile) = float(LDT_rc%glacier_index)
              lcfgrd(c,r,:) = 0.0
              lcfgrd(c,r,glacier_tile)  = 1.0
              ! Set original landcover glacier fraction at same point to 0.
              sfctype(c,r,LDT_rc%glacierclass) = 0.
              
           endif
        enddo
     enddo

     !remap the remaining glacier categories in the landcover data
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if(lcfgrd(c,r,LDT_rc%glacierclass).gt.&
                LDT_rc%gridcell_glacier_frac(n)) then 
              
              if(LDT_LSMparam_struc(n)%landcover%source.eq."MODIS_Native") then 
                 sfctype(c,r,glacier_tile) = float(LDT_rc%glacier_index)
                 lcfgrd(c,r,:) = 0.0
                 lcfgrd(c,r,glacier_tile) = 1.0              
                 sfctype(c,r,LDT_rc%glacierclass) = 0
              else
                 sfctype(c,r,:) = 0.0
                 sfctype(c,r,glacier_tile) = float(LDT_rc%glacier_index)
                 lcfgrd(c,r,:) = 0.0
                 lcfgrd(c,r,glacier_tile) = 1.0    
              endif
           endif
        enddo
     enddo

  else
     !remap the existing glacier categories in the landcover data
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if(lcfgrd(c,r,LDT_rc%glacierclass).gt.&
                LDT_rc%gridcell_glacier_frac(n)) then 
              sfctype(c,r,glacier_tile) = float(LDT_rc%glacier_index)
              lcfgrd(c,r,:) = 0.0
              lcfgrd(c,r,glacier_tile) = 1.0
              sfctype(c,r,LDT_rc%glacierclass) = 0

           else
               !do we need to renormalize?
              sfctype(c,r,glacier_tile) = 0
              lcfgrd(c,r,glacier_tile) = 0
              
              sfctype(c,r,LDT_rc%glacierclass) = 0
              lcfgrd(c,r,LDT_rc%glacierclass) = 0
           endif
        enddo
     enddo
     
  endif

  write(LDT_logunit,*) "[INFO] Setting glacier fractions to landcover water fraction ..."

  end subroutine LDT_assign_glaciersfctype

! ___________________________________________________________

!BOP
!
! !ROUTINE: LDT_assign_openwatersfctype
!  \label{LDT_assign_openwatersfctype}
!
! !REVISION HISTORY:
!  20 Apr  2013: KR Arsenault;  Separated surface type field
!                                assignment into own subroutine
! !INTERFACE:
 subroutine LDT_assign_openwatersfctype(n, total_types, &
                sfctype, localmask, lcfgrd, openwaterfgrd, lake_presence )
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_endrun
  use LDT_domainMod,only : isSurfaceTypeSelected

  implicit none

  integer, intent(in)     :: n             ! Nest
  integer, intent(in)     :: total_types   ! Total surface types
  real,    intent(inout)  :: lcfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),total_types)
  real,    intent(inout)  :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real,    intent(inout)  :: sfctype(LDT_rc%lnc(n),LDT_rc%lnr(n),total_types)
!  real,    intent(out)    :: openwaterfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),total_types)
  real,    intent(out)    :: openwaterfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
  real,intent(in),optional :: lake_presence(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine assigns surface types based on land/water mask 
!   and surface types selected.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[total_types]
!    number of total surface types
!   \item[lcfgrd]
!    landcover fractions of gridcell
!   \item[openwaterfgrd]
!    openwater fractions of gridcell
!   \item[localmask]
!    landmask for the region of interest
!   \item[sfctype]
!    surface type for the region of interest
!   \end{description}
!EOP      

   integer :: c, r, t, l, i

   if( LDT_rc%waterclass == 0. ) then
     write(LDT_logunit,*) "[ERR] No water class for associated landcover ..."
     write(LDT_logunit,*) " Stopping for now until solution entered."
     call LDT_endrun
   endif
   write(LDT_logunit,fmt=*) '[INFO] SURFACETYPE -- Assigning openwater surface types'

   openwaterfgrd = 0.
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)

      !- Landcover water fraction present:
         if( lcfgrd(c,r,LDT_rc%waterclass) > 0. ) then

           sfctype(c,r,total_types) = LDT_rc%openwater_index
           sfctype(c,r,LDT_rc%waterclass) = 0.

           lcfgrd(c,r,total_types)  = lcfgrd(c,r,LDT_rc%waterclass)
           lcfgrd(c,r,LDT_rc%waterclass) = 0.

           if( localmask(c,r) == 0 ) then
              localmask(c,r) = 1.
           endif

         elseif( lcfgrd(c,r,LDT_rc%waterclass) == 0. ) then
           lcfgrd(c,r,total_types)  = 0.

         end if     ! end land cover water fraction check

      end do
   end do
   openwaterfgrd(:,:,1) = lcfgrd(:,:,total_types)  ! Temporary assignment

   write(LDT_logunit,*) "[INFO] Setting openwater fractions to surface type array ..."

 end subroutine LDT_assign_openwatersfctype

! ___________________________________________________________


  subroutine surfacetype_writeHeader_LIS(n,ftn,dimID)
    
    use LDT_coreMod, only : LDT_rc

    integer      :: n 
    integer      :: ftn
    integer      :: dimID(3)
    integer      :: tdimID(3)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)
   
    call LDT_verify(nf90_def_dim(ftn,'sfctypes',&
         LDT_LSMparam_struc(n)%sfctype%vlevels,tdimID(3)))

    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LDT_LSMparam_struc(n)%sfctype)

    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LDT_LSMparam_struc(n)%landcover)

! - Enter surface model types, names:
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SFCMODELS", &
         LDT_LSMparam_struc(n)%sfctype%source))

#endif
  end subroutine Surfacetype_writeHeader_LIS

  subroutine surfacetype_writeHeader_LISHydro(n,ftn,dimID,flag)
    
    use LDT_coreMod, only : LDT_rc

    integer      :: n 
    integer      :: ftn
    integer      :: dimID(4)
    integer      :: tdimID(4)
    integer      :: flagn
    integer      :: flag

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)
    tdimID(4) = dimID(4)
   
    call LDT_verify(nf90_def_dim(ftn,'land_cat',&
         LDT_LSMparam_struc(n)%sfctype%vlevels,tdimID(3)))

    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LDT_LSMparam_struc(n)%sfctype,flagn)

    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         LDT_LSMparam_struc(n)%landcover,flagn)

! - Enter surface model types, names:
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SFCMODELS", &
         LDT_LSMparam_struc(n)%sfctype%source))

#endif
  end subroutine surfacetype_writeHeader_LISHydro


  subroutine LDT_surfacetype_writeData(n,ftn)

    use LDT_coreMod
#if ( defined SPMD )
    use mpi
#endif

    integer  :: n 
    integer  :: ftn
    integer  :: ierr

    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%sfctype)

!    print*, 'processed surface type ', LDT_localPet
#if ( defined SPMD )
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

    call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%landcover)

  end subroutine LDT_surfacetype_writeData

end module LDT_SurfaceTypeMod
