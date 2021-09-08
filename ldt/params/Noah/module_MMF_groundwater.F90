!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: MMF_groundwater
! \label{MMF_groundwater}
!
! !REVISION HISTORY:
!  27 Aug 2021: Sarith Mahanama
!     with inputs from David Mocko, Zhuo Wang ; Initial Specification
!
! !INTERFACE:

#include "LDT_misc.h"
module  MMF_groundwater
  
  use LDT_logMod,  only : LDT_logunit, LDT_endrun
  use LDT_coreMod
  use map_utils,   only : ij_to_latlon
  use LDT_gridmappingMod
  use LDT_paramMaskCheckMod
  
  implicit none
  private

  type, public :: MMF_BCsReader
   
     real, dimension(20)                 :: param_gridDesc, subparam_gridDesc
     integer                             :: glpnc, glpnr, subpnc, subpnr   
     integer,allocatable,dimension (:,:) :: lat_line, lon_line, local_mask
     
     
   contains

     procedure, public :: mi => mmf_init
     procedure, public :: mr => mmf_data_reader
     
  end type MMF_BCsReader
  
  integer,parameter                    :: NY_MMF = 21600
  integer,parameter                    :: NX_MMF = 43200

contains 
  
  SUBROUTINE mmf_init (MBR, nest, project)
    
    implicit none
    
    class (MMF_BCsReader), intent(inout) :: MBR
    integer, intent (in)                 :: nest
    character(*), intent (in)            :: project
    real                                 :: IN_xres, IN_yres
    integer :: c,r,gr,gc, glpnc, glpnr
    real    :: rlon(LDT_rc%lnc(nest),LDT_rc%lnr(nest)),rlat(LDT_rc%lnc(nest),LDT_rc%lnr(nest)) 
    real    :: param_grid(20)
      
    IN_xres    = 360./REAL(NX_MMF)
    IN_yres    = 180./REAL(NY_MMF)
    MBR%param_gridDesc(1)  = 0.          ! Latlon
    MBR%param_gridDesc(2)  = real(NX_MMF)
    MBR%param_gridDesc(3)  = real(NY_MMF)
    MBR%param_gridDesc(4)  = -90.0  + (IN_yres/2) ! LL lat 
    MBR%param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon 
    MBR%param_gridDesc(6)  = 128
    MBR%param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat
    MBR%param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon
    MBR%param_gridDesc(9)  = IN_yres     
    MBR%param_gridDesc(10) = IN_xres     
    MBR%param_gridDesc(20) = 64
    
    param_grid(:) = LDT_rc%mask_gridDesc(nest,:)
    glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
    glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1
    allocate (MBR%local_mask (LDT_rc%lnc(nest),LDT_rc%lnr(nest)))
    
    MBR%local_mask = LDT_rc%udef
    
    do r = 1, LDT_rc%lnr(nest)
       do c = 1, LDT_rc%lnc(nest)
          call ij_to_latlon(LDT_domain(nest)%ldtproj,float(c),float(r),&
               rlat(c,r),rlon(c,r))
          gr = nint((rlat(c,r)-param_grid(4))/param_grid(10)) + 1
          gc = nint((rlon(c,r)-param_grid(5))/param_grid( 9)) + 1
          if(LDT_rc%global_mask(gc,gr) > 0. ) MBR%local_mask(c,r) = 1
       end do
    end do
    
    ! ------------------------------------------------------------
    !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN BCS DATA
    ! ------------------------------------------------------------
    
    !- Map Parameter Grid Info to LIS Target Grid/Projection Info --
    
    MBR%subparam_gridDesc = 0.

    call LDT_RunDomainPts( nest, project, MBR%param_gridDesc(:),        &
         MBR%glpnc, MBR%glpnr, MBR%subpnc, MBR%subpnr, MBR%subparam_gridDesc, &
         MBR%lat_line, MBR%lon_line)

  END SUBROUTINE mmf_init
  
  ! ----------------------------------------------------------------
  
  SUBROUTINE mmf_data_reader (MBR, nest, varname, datadir, lisout, MMF_fillopts)
    
    implicit none

    class (MMF_BCsReader), intent(inout)    :: MBR
    character(*), intent (in)               :: varname, datadir
    real, dimension (:,:,:), intent (inout) :: lisout
    integer, intent (in)                    :: nest
    type(LDT_fillopts)                      :: MMF_fillopts
    
    ! Adapted from Zhuo Wang's
    ! /discover/nobackup/projects/lis_aist17/zwang9/Geogrid/geogrid2netcdf.f90
    
    integer, parameter :: iSigned = 1
    integer, parameter :: endian = 0
    integer, parameter :: wordsize = 2
    integer, parameter :: missing = -9999
    integer, parameter :: ntsteps =1
    integer, parameter :: nlat_tile = 1200
    integer, parameter :: nlon_tile = 1200
    character(len=1024):: geogridFile
    character(len=80)  :: tileName
    integer            :: num_tile_lon,num_tile_lat
    integer            :: rc,status,tx,ty,txs,txe,tys,tye    
    real,allocatable   :: tarray(:,:,:), garray(:,:,:)
    real               :: SF
    
    ! c function from Michael G. Duda, NCAR/MMM
    integer,  external :: read_geogrid

    ! set scale factor

    SF = 1.
    select case (trim(varname))

    case ('T')
       ! transmissivity
       SF = 0.01
       
    case ('R')
       ! recharge
       SF = 0.03
       
    case ('E')
       ! river bed elevation
       SF = 0.3
       
    case ('W')
       ! water table depth
       SF = -0.03
       
    case default  
       write(LDT_logunit,*) '[ERROR] Unknown MMF data field : ', trim (varname)
       call LDT_endrun
    end select

    ! Reading GEOGRID data
    ! --------------------
    
    num_tile_lon = NX_MMF / nlon_tile
    num_tile_lat = NY_MMF / nlat_tile

    allocate(tarray(nlon_tile,nlat_tile,ntsteps))
    allocate(garray(NX_MMF, NY_MMF, ntsteps))
    garray = LDT_rc%udef  ! global array (43200,21600) is constructed by assebling 36x18 # of tarrays (1200,1200) 

    TILE_COLS: do tx = 1, num_tile_lon
       
       txs = 1 + (tx-1)*nlon_tile
       txe = tx * nlon_tile

       TILE_ROWS: do ty = 1, num_tile_lat
          
          tarray = missing ! GEOGRID tile array 1200,1200
          tys = 1 + (ty-1)*nlat_tile
          tye = ty * nlat_tile

          write(tileName, fmt='(4(a1,i5.5))') '/',txs, '-', txe, '.', tys, '-', tye
          geogridFile = trim(datadir)//trim(tileName)
          
          rc = read_geogrid(trim(geogridFile),len(trim(geogridFile)),tarray, &
               nlon_tile,nlat_tile,ntsteps,isigned,endian,1.,wordsize,status)

          if (rc == 1 .or. status == 1) then
             write(LDT_logunit,*) '[ERROR] reading GEOGRID file : ',trim(geogridFile)
             call LDT_endrun
           end if
             
          where(tarray == missing)
             tarray = LDT_rc%udef
          elsewhere
             tarray = SF * tarray
          end where
          
          garray(txs:txe,tys:tye,1) = tarray (1:nlon_tile, 1:nlat_tile,1)
          
       end do TILE_ROWS
    end do TILE_COLS

    ! regrid garray and construct lisout array on the LIS grid
    
    call regrid_to_lisgrid (nest, garray, lisout)

    ! NOTE : Use MMF_fillopts%filltype and MMF_fillopts%fillradius
    ! to fill  where lisout values are missing compared
    ! to land/water mask:
    
    deallocate (tarray, garray)
    
  contains
      
    ! ---------------------------------------------------------------------
    
    SUBROUTINE regrid_to_lisgrid (nest, data_in, var_subset)
      
      ! adapted from params/irrigation/read_GRIPC_irrigfrac.F90 
      
      implicit none
      integer, intent (in)                  :: nest
      real, dimension(:,:,:), intent   (in) :: data_in
      real, dimension(:,:,:), intent(inout) :: var_subset
      
      ! arrays read from LL
      integer,parameter :: k10 = selected_int_kind(10)
      integer   :: mi                     ! Total number of input param grid array points
      integer   :: mo                     ! Total number of output LIS grid array points
      integer   :: nc, nr, i, j
      integer, allocatable  :: n11(:)     ! Array that maps the location of each input grid
      !   point in the output grid. 
      real,    allocatable  :: gi(:)      ! Input parameter 1d grid
      logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)
      real, allocatable, dimension (:,:)    :: var_in
      real, allocatable, dimension (:)      :: go2     ! Output lis 1d grid
      logical*1, allocatable, dimension (:) :: lo2  ! Output logical mask (to match go)
      
      !mi = INT(dble(MBR%NX)*dble(MBR%NY), 8)
      mi = MBR%subpnc*MBR%subpnr
      mo = LDT_rc%lnc(nest)*LDT_rc%lnr(nest)
      
      allocate (var_in(MBR%subpnc,MBR%subpnr))                
      allocate ( li(mi), gi (mi), n11(mi))
      allocate (go2(mo), lo2(mo))
      
      !- Create mapping between parameter domain and LIS grid domain:
      call upscaleByAveraging_input( MBR%subparam_gridDesc, &
           LDT_rc%gridDesc(nest,:), mi, mo, n11)
      
      Data_Fields: do j = 1, size (data_in, 3)
         
         var_in = LDT_rc%udef
         do nr = 1, MBR%subpnr
            do nc = 1, MBR%subpnc
               var_in (nc,nr) = data_in (MBR%lon_line(nc,1),MBR%lat_line (1,nr),j)
            end do
         end do
         
         ! -------------------------------------------------------------------
         !     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
         ! -------------------------------------------------------------------
         
         gi  = LDT_rc%udef
         li  = .false.
         
         !- Assign 2-D array to 1-D for aggregation routines:
         i = 0
         do nr = 1, MBR%subpnr
            do nc = 1, MBR%subpnc
               i = i + 1            
               if( var_in(nc,nr) .NE. LDT_rc%udef) then  
                  gi(i) = var_in(nc,nr)  
                  li(i) = .true.
               endif
            end do
         enddo
         
         !- Spatial average within each coarse gridcell:
         call upscaleByAveraging ( mi, mo, LDT_rc%udef, n11, li, gi, &
              lo2(:), go2 (:))
         i = 0
         do nr = 1, LDT_rc%lnr(nest)
            do nc = 1, LDT_rc%lnc(nest)
               i = i + 1
               var_subset (nc,nr,j) = go2(i)
            enddo
         enddo
         
      end do Data_Fields
      deallocate (var_in, li, gi, n11, go2, lo2)
      
    END SUBROUTINE regrid_to_lisgrid
    
  END SUBROUTINE mmf_data_reader

  ! ----------------------------------------------------------------

end module MMF_groundwater

  
