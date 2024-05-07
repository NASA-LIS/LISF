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
#include "vic412_usr_def.h"

module LDT_vic412rstMod
! !MODULE: LDT_vic412rstMod
!
! !DESCRIPTION:
!   The code in this file provides interfaces to manage
!   the processing of climatological restart files specifically
!   for VIC4.1.2.  VIC needs additional layer of restart data structure
!   to unpack each variables and accumulate.
!   Routines are based on the VIC C codes in LIS7.
!   Pre-processor definitions need to be specified in vic412_usr_def.h
!   according to the vic simulation specification.
!
! !REVISION HISTORY:
!  18 Oct 2017    H. Beaudoing:  Initial Specification

 use VIC_parmsMod

 implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_vic412rstInit    !allocates memory for required structures
  public :: LDT_vic412rstDiagnose
  public :: LDT_vic412rstAvePack
  public :: LDT_vic412rstFinalize
  public :: count_model_state_412

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: rstdec

  type, public :: vic_state_dec
   integer :: tmp_cellnum
   integer :: tmp_Nveg
   integer :: tmp_Nband
   integer, allocatable, dimension(:,:) :: tveg, tband  !Nveg,Nband
   real*8, allocatable,dimension(:) :: dz_node   !Nnode
   real*8, allocatable,dimension(:) :: Zsum_node  !Nnode
   real*8, allocatable,dimension(:) :: mu  !Nveg
   logical*1, allocatable,dimension(:) :: init_STILL_STORM  !Nveg: charac*1
   integer, allocatable,dimension(:) :: init_DRY_TIME  !Nveg
#if EXCESS_ICE
   real*8, allocatable,dimension(:) :: depth !Nlayer
   real*8, allocatable,dimension(:) :: effective_porosity !Nlayer
   real*8  :: dp
#endif
   real*8, allocatable,dimension(:,:,:,:) :: moist  !NdistxNvegxNbandxNlayer
#if SPATIAL_FROST
   real*8, allocatable,dimension(:,:,:,:,:) :: ice  !NdistxNvegxNbandxNlayerxFROST_SUBAREAS];
#else
   real*8, allocatable,dimension(:,:,:,:) :: ice  !NdistxNvegxNbandxNlayer
#endif
   real*8, allocatable,dimension(:,:,:) :: Wdew  !NdistxNvegxNband
! NvegxNband
   integer, allocatable,dimension(:,:) :: last_snow  
   logical*1, allocatable,dimension(:,:) :: MELTING  !character*1
   real*8, allocatable,dimension(:,:) :: coverage
   real*8, allocatable,dimension(:,:) :: swq
   real*8, allocatable,dimension(:,:) :: surf_temp
   real*8, allocatable,dimension(:,:) :: surf_water
   real*8, allocatable,dimension(:,:) :: pack_temp, pack_water, density
   real*8, allocatable,dimension(:,:) :: coldcontent, snow_canopy
   real*8, allocatable,dimension(:,:,:) :: T1 ! NvegxNbandxNnode
#if LAKES
   real*8, allocatable,dimension(:) :: lake_moist  ![MAX_LAYERS];
#if SPATIAL_FROST
   real*8, allocatable,dimension(:,:) :: lake_ice  ![MAX_LAYERS][FROST_SUBAREAS];
#else
   real*8, allocatable,dimension(:) :: lake_ice  ![MAX_LAYERS];
#endif
   integer :: lake_last_snow  
   logical*1 :: lake_MELTING  ! character(len=1)
   real*8 :: lake_coverage, lake_swq, lake_surf_temp,lake_surf_water
   real*8 :: lake_pack_temp, lake_pack_water, lake_density
   real*8 :: lake_coldcontent, lake_snow_canopy
   real*8, allocatable,dimension(:) :: lake_T1 ! Nnode
   integer :: lake_activenod
   real*8  :: lake_dz,lake_surfdz,lake_ldepth
   real*8, allocatable, dimension(:) :: lake_surface,lake_temp  !node
   real*8  :: lake_tempavg,lake_areai,lake_new_ice_area,lake_ice_water_eq
   real*8  :: lake_hice,lake_tempi,lake_swe,lake_SAlbedo,lake_sdepth
   real*8  :: lake_sarea,lake_volume
   real*8  :: lake_depth
#endif
  integer, allocatable,dimension(:,:) :: last_snow_cnt  ! count  
  end type vic_state_dec
  type, public :: modeltiledec
     type(vic_state_dec), allocatable :: vic(:)
  end type modeltiledec
  type (modeltiledec), allocatable :: rstdec(:)

  contains 

!BOP
! !ROUTINE: LDT_vic412rstInit
! label{LDT_vic412rstInit}
!
! !INTERFACE:
    subroutine LDT_vic412rstInit(n,nfiles,npatch,state_chunk_size)
!
! !DESCRIPTION:
!
!  This routine performs initialization in the restart processing mode
!  for VIC412.

  integer, intent(in) :: n
  integer, intent(in) :: nfiles
  integer, intent(in) :: npatch
  integer, intent(in) :: state_chunk_size != 1610
! moved to ldt.config
!  integer, parameter :: Nnode = 12
!  integer, parameter :: Nlayer = 3
!  integer, parameter :: Ndist = 1
!  integer, parameter :: Nveg = 1   !13
!  integer, parameter :: Nbands = 25
  integer             :: nf, t

!-------------------------------------------------------------------
!-- allocate vic restart array variables to be accumulated
    allocate(rstdec(nfiles))
    do nf=1,nfiles
     allocate(rstdec(nf)%vic(npatch))
    end do
    do nf=1,nfiles
     do t=1, npatch
      allocate(rstdec(nf)%vic(t)%dz_node(VIC_struc(n)%Nnode))
      allocate(rstdec(nf)%vic(t)%Zsum_node(VIC_struc(n)%Nnode))
      allocate(rstdec(nf)%vic(t)%mu(VIC_struc(n)%Nveg+1))
      allocate(rstdec(nf)%vic(t)%tveg(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%tband(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%init_STILL_STORM(VIC_struc(n)%Nveg+1))
      allocate(rstdec(nf)%vic(t)%init_DRY_TIME(VIC_struc(n)%Nveg+1))
#if EXCESS_ICE
      allocate(rstdec(nf)%vic(t)%depth(VIC_struc(n)%Nlayer))
      allocate(rstdec(nf)%vic(t)%effective_porosity(VIC_struc(n)%Nlayer))
#endif
      allocate(rstdec(nf)%vic(t)%moist(VIC_struc(n)%Ndist,VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands,VIC_struc(n)%Nlayer)); rstdec(nf)%vic(t)%moist = 0.0
#if SPATIAL_FROST
      allocate(rstdec(nf)%vic(t)%ice(VIC_struc(n)%Ndist,VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands,VIC_struc(n)%Nlayer,FROST_SUBAREAS))
#else
      allocate(rstdec(nf)%vic(t)%ice(VIC_struc(n)%Ndist,VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands,VIC_struc(n)%Nlayer))
#endif
      allocate(rstdec(nf)%vic(t)%Wdew(VIC_struc(n)%Ndist,VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%last_snow(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%last_snow_cnt(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%MELTING(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%coverage(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%swq(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%surf_temp(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%surf_water(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%pack_temp(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%pack_water(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%density(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%coldcontent(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%snow_canopy(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands))
      allocate(rstdec(nf)%vic(t)%T1(VIC_struc(n)%Nveg+1,VIC_struc(n)%Nbands,VIC_struc(n)%Nnode))
#if LAKES
      allocate(rstdec(nf)%vic(t)%lake_moist(VIC_struc(n)%Nlayer))
#if SPATIAL_FROST
      allocate(rstdec(nf)%vic(t)%lake_ice(VIC_struc(n)%Nlayer,FROST_SUBAREAS))
#else
      allocate(rstdec(nf)%vic(t)%lake_ice(VIC_struc(n)%Nlayer))
#endif
      allocate(rstdec(nf)%vic(t)%lake_T1(VIC_struc(n)%Nnode))
      allocate(rstdec(nf)%vic(t)%lake_surface(VIC_struc(n)%Nnode))
      allocate(rstdec(nf)%vic(t)%lake_temp(VIC_struc(n)%Nnode))
#endif
!-- and initialize summing arrays
      rstdec(nf)%vic(t)%dz_node = 0.0
      rstdec(nf)%vic(t)%Zsum_node = 0.0
      rstdec(nf)%vic(t)%mu = 0.0
      rstdec(nf)%vic(t)%tveg = 0
      rstdec(nf)%vic(t)%tband = 0
      rstdec(nf)%vic(t)%init_STILL_STORM = .false.
      rstdec(nf)%vic(t)%init_DRY_TIME = 0
#if EXCESS_ICE
      rstdec(nf)%vic(t)%depth = 0.0
      rstdec(nf)%vic(t)%effective_porosity = 0.0
#endif
      rstdec(nf)%vic(t)%moist = 0.0
#if SPATIAL_FROST
      rstdec(nf)%vic(t)%ice = 0.0
#else
      rstdec(nf)%vic(t)%ice = 0.0
#endif
      rstdec(nf)%vic(t)%last_snow = 0
      rstdec(nf)%vic(t)%Wdew = 0.0
      rstdec(nf)%vic(t)%last_snow_cnt = 0
      rstdec(nf)%vic(t)%MELTING = .false.
      rstdec(nf)%vic(t)%coverage = 0.0
      rstdec(nf)%vic(t)%swq = 0.0
      rstdec(nf)%vic(t)%surf_temp = 0.0
      rstdec(nf)%vic(t)%surf_water = 0.0
      rstdec(nf)%vic(t)%pack_temp = 0.0
      rstdec(nf)%vic(t)%pack_water = 0.0
      rstdec(nf)%vic(t)%density = 0.0
      rstdec(nf)%vic(t)%coldcontent = 0.0
      rstdec(nf)%vic(t)%snow_canopy = 0.0
      rstdec(nf)%vic(t)%T1 = 0.0
#if LAKES
      rstdec(nf)%vic(t)%lake_moist = 0.0
#if SPATIAL_FROST
      rstdec(nf)%vic(t)%lake_ice = 0.0
#else
      rstdec(nf)%vic(t)%lake_ice = 0.0
#endif
      rstdec(nf)%vic(t)%lake_T1 = 0.0
      rstdec(nf)%vic(t)%lake_surface = 0.0
      rstdec(nf)%vic(t)%lake_temp = 0.0
#endif
    end do   ! t
   end do   ! nf

    end subroutine LDT_vic412rstInit

    subroutine LDT_vic412rstDiagnose(n,tindex,npatch,state_chunk_size,var)
! The routine is based on unpack_model_state.c
! Accumulate VIC variables in rstdec%vic%variable structure.
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: tindex  ! month index
  integer, intent(in) :: npatch  ! tile 
  integer, intent(in) :: state_chunk_size  
  real,dimension(npatch,state_chunk_size),intent(in) :: var
  real,dimension(state_chunk_size) :: states
  integer :: t, l
  integer :: count
  integer :: iveg, iband
  integer :: nidx,lidx,veg,band,dist,frost_area,node
  integer :: activenod
  integer :: ival
  real*8  :: dval, tmpval
  real    :: fval
  character(len=1) :: cval
  logical*1 :: lval
!------------------------------------------------------ 
! Unpack state_chunks and accumulate
!------------------------------------------------------ 
 do t=1, npatch
  states = var(t,:)
  count = 1
  call unpack_state_int(states(count),count,ival)
  rstdec(tindex)%vic(t)%tmp_cellnum = ival
  call unpack_state_int(states(count),count,ival)
  rstdec(tindex)%vic(t)%tmp_Nveg = ival
  call unpack_state_int(states(count),count,ival)
  rstdec(tindex)%vic(t)%tmp_Nband = ival
!    /* Read soil thermal node deltas */
    do nidx=1, VIC_struc(n)%Nnode
       call unpack_state_double(states(count), count, dval)
       rstdec(tindex)%vic(t)%dz_node(nidx) =  &
                     rstdec(tindex)%vic(t)%dz_node(nidx) + dval
    end do ! Nnode
    if ( VIC_struc(n)%Nnode == 1 ) rstdec(tindex)%vic(t)%dz_node(1) = 0
!    /* Read soil thermal node depths */
    do nidx=1, VIC_struc(n)%Nnode
       call unpack_state_double(states(count), count, dval)
       rstdec(tindex)%vic(t)%Zsum_node(nidx) = &
                     rstdec(tindex)%vic(t)%Zsum_node(nidx) + dval
    end do ! Nnode
    if ( VIC_struc(n)%Nnode == 1 ) rstdec(tindex)%vic(t)%Zsum_node(1) = 0
#if EXCESS_ICE
!    /* Read soil depth */
    do lidx=1, VIC_struc(n)%Nlayer
       call unpack_state_double(states(count), count, dval)
       rstdec(tindex)%vic(t)%depth(lidx) = &
                      rstdec(tindex)%vic(t)%depth(lidx) + dval
    end do ! Nlayer
!    /* Read effective porosity */
    do lidx=1, VIC_struc(n)%Nlayer
       call unpack_state_double(states(count), count, dval)
       rstdec(tindex)%vic(t)%effective_porosity(lidx) =  &
                      rstdec(tindex)%vic(t)%effective_porosity(lidx) + dval
    end do ! Nlayer
!    /* Reading damping depth */
    call unpack_state_double(states(count), count, dval);
    rstdec(tindex)%vic(t)%dp = rstdec(tindex)%vic(t)%dp + dval
#endif
!    /* Input for all vegetation types */
    do veg = 1, VIC_struc(n)%Nveg+1
!        // read distributed precipitation variables
        call unpack_state_double(states(count), count, dval)
        rstdec(tindex)%vic(t)%mu(veg) = rstdec(tindex)%vic(t)%mu(veg) + dval
        call unpack_state_char(states(count), count, lval)
        rstdec(tindex)%vic(t)%init_STILL_STORM(veg) = lval
        call unpack_state_int(states(count), count, ival)
        rstdec(tindex)%vic(t)%init_DRY_TIME(veg) = rstdec(tindex)%vic(t)%init_DRY_TIME(veg) + ival

!        /* Input for all snow bands */
        do band = 1, VIC_struc(n)%Nbands
!            /* Read cell identification information */
            call unpack_state_int(states(count), count, iveg)
            rstdec(tindex)%vic(t)%tveg(veg,band) = iveg
            call unpack_state_int(states(count), count, iband)
            rstdec(tindex)%vic(t)%tband(veg,band) = iband

            if ( iveg /= (veg-1) .or. iband /= (band-1) ) then
                print*,"The vegetation and snow band indices in the model &
                        state file (veg = ",iveg,", band = ",iband,") do not &
                        match those currently requested (veg = ",(veg-1), &
                        ", band = ",(band-1),").  Model state file must be &
                        stored with variables for all vegetation indexed by &
                        variables for all snow bands."
                !nrerror(ErrStr);
            endif  

!            // Read both wet and dry fractions if using distributed precipitation
            do dist = 1, VIC_struc(n)%Ndist
!                /* Read total soil moisture */
                do lidx = 1, VIC_struc(n)%Nlayer
                    call unpack_state_double(states(count), count, dval)
                    rstdec(tindex)%vic(t)%moist(dist,veg,band,lidx) = &
                      rstdec(tindex)%vic(t)%moist(dist,veg,band,lidx) + dval
                end do ! lidx 

!                /* Read average ice content */
                do lidx = 1, VIC_struc(n)%Nlayer
#if SPATIAL_FROST
                    do frost_area = 1, FROST_SUBAREAS
                        call unpack_state_double(states(count), count, dval)
                        rstdec(tindex)%vic(t)%ice(dist,veg,band,lidx,frost_area) = rstdec(tindex)%vic(t)%ice(dist,veg,band,lidx,forst_area) + dval
                    end do
#else
                    call unpack_state_double(states(count), count, dval)
                    rstdec(tindex)%vic(t)%ice(dist,veg,band,lidx) = rstdec(tindex)%vic(t)%ice(dist,veg,band,lidx) + dval
#endif
                end do ! lidx

!                /* Read dew storage */
                if ( veg < VIC_struc(n)%Nveg+1 ) then
                    call unpack_state_double(states(count), count, dval) 
                    rstdec(tindex)%vic(t)%Wdew(dist,veg,band) =  &
                       rstdec(tindex)%vic(t)%Wdew(dist,veg,band) + dval
                endif
            end do ! dist 

!            /* Read snow data */
            call unpack_state_int(states(count), count, ival)
            if ( ival .ne. -99999 ) then
             rstdec(tindex)%vic(t)%last_snow(veg,band) = &
                rstdec(tindex)%vic(t)%last_snow(veg,band) + ival
             rstdec(tindex)%vic(t)%last_snow_cnt(veg,band) =  &
                rstdec(tindex)%vic(t)%last_snow_cnt(veg,band)+1
            endif
            call unpack_state_char(states(count), count, lval)
            rstdec(tindex)%vic(t)%MELTING(veg,band) = lval 
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%coverage(veg,band) =  &
               rstdec(tindex)%vic(t)%coverage(veg,band) + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%swq(veg,band) = &
               rstdec(tindex)%vic(t)%swq(veg,band) + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%surf_temp(veg,band) = &
               rstdec(tindex)%vic(t)%surf_temp(veg,band) + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%surf_water(veg,band) = &
               rstdec(tindex)%vic(t)%surf_water(veg,band) + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%pack_temp(veg,band) =  &
               rstdec(tindex)%vic(t)%pack_temp(veg,band) + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%pack_water(veg,band) = &
               rstdec(tindex)%vic(t)%pack_water(veg,band) + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%density(veg,band) =  &
               rstdec(tindex)%vic(t)%density(veg,band) + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%coldcontent(veg,band) = &
               rstdec(tindex)%vic(t)%coldcontent(veg,band) + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%snow_canopy(veg,band) =  &
               rstdec(tindex)%vic(t)%snow_canopy(veg,band) + dval

!            /* Read soil thermal node temperatures */
            do nidx = 1, VIC_struc(n)%Nnode
                call unpack_state_double(states(count), count, dval)
                rstdec(tindex)%vic(t)%T1(veg,band,nidx) =  &
                   rstdec(tindex)%vic(t)%T1(veg,band,nidx) + dval
            end do  !nidx 
        end do  !band

#if LAKES
!        if ( LAKES ) then
!            // Read both wet and dry fractions if using distributed precipitation
            do dist = 1, VIC_struc(n)%Ndist
!                /* Read total soil moisture */
                do lidx = 1, VIC_struc(n)%Nlayer
                    call unpack_state_double(states(count), count, dval)
                    rstdec(tindex)%vic(t)%lake_moist(lidx) = rstdec(tindex)%vic(t)%lake_moist(lidx) + dval
                end do

!                /* Read average ice content */
                do lidx = 1, VIC_struc(n)%Nlayer
#if SPATIAL_FROST
                    do frost_area = 1, FROST_SUBAREAS
                        call unpack_state_double(states(count), count, dval)
                        rstdec(tindex)%vic(t)%lake_ice(lidx,frost_area) = &
                           rstdec(tindex)%vic(t)%lake_ice(lidx,frost_area) + dval
                    end do 
#else
                    call unpack_state_double(states(count), count, dval)
                    rstdec(tindex)%vic(t)%lake_ice(lidx) = rstdec(tindex)%vic(t)%lake_ice(lidx) + dval
#endif
                end do

            end do ! dist

!            /* Read snow data */
            call unpack_state_int(states(count), count, ival)
            rstdec(tindex)%vic(t)%lake_last_snow = rstdec(tindex)%vic(t)%lake_last_snow + ival
            call unpack_state_char(states(count), count, lval)
            rstdec(tindex)%vic(t)%lake_MELTING = lval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_coverage = rstdec(tindex)%vic(t)%lake_coverage + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_swq = rstdec(tindex)%vic(t)%lake_swq + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_surf_temp = rstdec(tindex)%vic(t)%lake_surf_temp + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_surf_water = rstdec(tindex)%vic(t)%lake_surf_water + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_pack_temp = rstdec(tindex)%vic(t)%lake_pack_temp + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_pack_water = rstdec(tindex)%vic(t)%lake_pack_water + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_density = rstdec(tindex)%vic(t)%lake_density + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_coldcontent = rstdec(tindex)%vic(t)%lake_coldcontent + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_snow_canopy = rstdec(tindex)%vic(t)%lake_snow_canopy + dval

!            /* Read soil thermal node temperatures */
            do nidx = 1, VIC_struc(n)%Nnode
                call unpack_state_double(states(count), count, dval)
                rstdec(tindex)%vic(t)%lake_T1(nidx) = rstdec(tindex)%vic(t)%lake_T1(nidx) + dval
            end do

!            /* Read lake-specific variables */
            call unpack_state_int(states(count), count, ival)
            activenod = ival
            rstdec(tindex)%vic(t)%lake_activenod = rstdec(tindex)%vic(t)%lake_activenod + ival
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_dz = rstdec(tindex)%vic(t)%lake_dz + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_surfdz = rstdec(tindex)%vic(t)%lake_surfdz + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_ldepth = rstdec(tindex)%vic(t)%lake_ldepth + dval

            do node = 1, activenod
                call unpack_state_double(states(count), count, dval)
                rstdec(tindex)%vic(t)%lake_surface(node) = rstdec(tindex)%vic(t)%lake_surface(node) + dval
            end do

            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_sarea = rstdec(tindex)%vic(t)%lake_sarea + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_volume = rstdec(tindex)%vic(t)%lake_volume + dval

            do node = 1, activenod
                call unpack_state_double(states(count), count, dval)
                rstdec(tindex)%vic(t)%lake_temp(node) = rstdec(tindex)%vic(t)%lake_temp(node) + dval
            end do 

            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_tempavg = rstdec(tindex)%vic(t)%lake_tempavg + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_areai = rstdec(tindex)%vic(t)%lake_areai + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_new_ice_area = rstdec(tindex)%vic(t)%lake_new_ice_area + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_ice_water_eq = rstdec(tindex)%vic(t)%lake_ice_water_eq + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_hice = rstdec(tindex)%vic(t)%lake_hice + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_tempi = rstdec(tindex)%vic(t)%lake_tempi + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_swe = rstdec(tindex)%vic(t)%lake_swe + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_surf_temp = rstdec(tindex)%vic(t)%lake_surf_temp + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_pack_temp = rstdec(tindex)%vic(t)%lake_pack_temp + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_coldcontent = rstdec(tindex)%vic(t)%lake_coldcontent + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_surf_water = rstdec(tindex)%vic(t)%lake_surf_water + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_pack_water = rstdec(tindex)%vic(t)%lake_pack_water + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_SAlbedo = rstdec(tindex)%vic(t)%lake_SAlbedo + dval
            call unpack_state_double(states(count), count, dval)
            rstdec(tindex)%vic(t)%lake_sdepth = rstdec(tindex)%vic(t)%lake_sdepth + dval
!        endif  ! LAKE 
#endif
    end do ! Nveg 

 end do  ! t

    end subroutine LDT_vic412rstDiagnose

    subroutine LDT_vic412rstAvePack(n,tindex,npatch,state_chunk_size,noutVar)
! Compute average of vic variables and pack into a 2D variable.
! The packing routine is based on pack_model_state.c
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: tindex  ! month index
  integer, intent(in) :: npatch  ! tile 
  integer, intent(in) :: state_chunk_size  
  integer, dimension(npatch), intent(in) :: noutVar
  real,dimension(state_chunk_size) :: states
  integer :: t, l
  integer :: count
  integer :: iveg, iband
  integer :: nidx,lidx,veg,band,dist,frost_area,node
  integer :: activenod
  integer :: ival
  real*8  :: dval, tmpval
  real    :: fval
  character(len=1) :: cval
  logical*1 :: lval

!------------------------------------------------------ 
! Take average
!------------------------------------------------------ 
 do t=1, npatch
    do nidx=1, VIC_struc(n)%Nnode
       rstdec(tindex)%vic(t)%dz_node(nidx) = rstdec(tindex)%vic(t)%dz_node(nidx)/noutVar(t)
    end do ! Nnode
    do nidx=1, VIC_struc(n)%Nnode
       rstdec(tindex)%vic(t)%Zsum_node(nidx) = rstdec(tindex)%vic(t)%Zsum_node(nidx)/noutVar(t)
    end do ! Nnode
#if EXCESS_ICE
    do lidx=1, VIC_struc(n)%Nlayer
       rstdec(tindex)%vic(t)%depth(lidx) = rstdec(tindex)%vic(t)%depth(lidx)/noutVar(t)
    end do ! Nlayer
    do lidx=1, VIC_struc(n)%Nlayer
       rstdec(tindex)%vic(t)%effective_porosity(lidx) = rstdec(tindex)%vic(t)%effective_porosity(lidx)/noutVar(t)
    end do ! Nlayer
    rstdec(tindex)%vic(t)%dp = rstdec(tindex)%vic(t)%dp/noutVar(t)
#endif
    do veg = 1, VIC_struc(n)%Nveg+1
        rstdec(tindex)%vic(t)%mu(veg) = rstdec(tindex)%vic(t)%mu(veg)/noutVar(t)
        rstdec(tindex)%vic(t)%init_DRY_TIME(veg) = rstdec(tindex)%vic(t)%init_DRY_TIME(veg)/noutVar(t)

        do band = 1, VIC_struc(n)%Nbands
            do dist = 1, VIC_struc(n)%Ndist
                do lidx = 1, VIC_struc(n)%Nlayer
                    rstdec(tindex)%vic(t)%moist(dist,veg,band,lidx) = rstdec(tindex)%vic(t)%moist(dist,veg,band,lidx)/noutVar(t)
                end do ! lidx 

!                /* Read average ice content */
                do lidx = 1, VIC_struc(n)%Nlayer
#if SPATIAL_FROST
                    do frost_area = 1, FROST_SUBAREAS
                        rstdec(tindex)%vic(t)%ice(dist,veg,band,lidx,frost_area) = rstdec(tindex)%vic(t)%ice(dist,veg,band,lidx,forst_area)/noutVar(t)
                    end do
#else
                    rstdec(tindex)%vic(t)%ice(dist,veg,band,lidx) = rstdec(tindex)%vic(t)%ice(dist,veg,band,lidx)/noutVar(t)
#endif
                end do ! lidx

!                /* Read dew storage */
                if ( veg < VIC_struc(n)%Nveg+1 ) then
                    rstdec(tindex)%vic(t)%Wdew(dist,veg,band) = &
                                    rstdec(tindex)%vic(t)%Wdew(dist,veg,band)/noutVar(t)
                endif
            end do ! dist 

!            rstdec(tindex)%vic(t)%last_snow(veg,band) = rstdec(tindex)%vic(t)%last_snow(veg,band)/noutVar(t)
            if ( rstdec(tindex)%vic(t)%last_snow_cnt(veg,band) .gt. 0 ) then
            rstdec(tindex)%vic(t)%last_snow(veg,band) = &
                           rstdec(tindex)%vic(t)%last_snow(veg,band) / &
                           rstdec(tindex)%vic(t)%last_snow_cnt(veg,band)
            else
            rstdec(tindex)%vic(t)%last_snow(veg,band) = -99999
            endif
            rstdec(tindex)%vic(t)%coverage(veg,band) = rstdec(tindex)%vic(t)%coverage(veg,band)/noutVar(t)
            rstdec(tindex)%vic(t)%swq(veg,band) = rstdec(tindex)%vic(t)%swq(veg,band)/noutVar(t)
            rstdec(tindex)%vic(t)%surf_temp(veg,band) = rstdec(tindex)%vic(t)%surf_temp(veg,band)/noutVar(t)
            rstdec(tindex)%vic(t)%surf_water(veg,band) = rstdec(tindex)%vic(t)%surf_water(veg,band)/noutVar(t)
            rstdec(tindex)%vic(t)%pack_temp(veg,band) = rstdec(tindex)%vic(t)%pack_temp(veg,band)/noutVar(t)
            rstdec(tindex)%vic(t)%pack_water(veg,band) = rstdec(tindex)%vic(t)%pack_water(veg,band)/noutVar(t)
            rstdec(tindex)%vic(t)%density(veg,band) = rstdec(tindex)%vic(t)%density(veg,band)/noutVar(t)
            rstdec(tindex)%vic(t)%coldcontent(veg,band) = rstdec(tindex)%vic(t)%coldcontent(veg,band)/noutVar(t)
            rstdec(tindex)%vic(t)%snow_canopy(veg,band) = rstdec(tindex)%vic(t)%snow_canopy(veg,band)/noutVar(t)

!            /* Read soil thermal node temperatures */
            do nidx = 1, VIC_struc(n)%Nnode
                rstdec(tindex)%vic(t)%T1(veg,band,nidx) = rstdec(tindex)%vic(t)%T1(veg,band,nidx)/noutVar(t)
            end do  !nidx 
        end do  !band

#if LAKES
!        if ( LAKES ) then
!            // Read both wet and dry fractions if using distributed precipitation
            do dist = 1, VIC_struc(n)%Ndist
                do lidx = 1, VIC_struc(n)%Nlayer
                    rstdec(tindex)%vic(t)%lake_moist(lidx) = rstdec(tindex)%vic(t)%lake_moist(lidx)/noutVar(t)
                end do

                do lidx = 1, VIC_struc(n)%Nlayer
#if SPATIAL_FROST
                    do frost_area = 1, FROST_SUBAREAS
                        rstdec(tindex)%vic(t)%lake_ice(lidx,frost_area) = &
                           rstdec(tindex)%vic(t)%lake_ice(lidx,frost_area)/noutVar(t)
                    end do 
#else
                    rstdec(tindex)%vic(t)%lake_ice(lidx) = rstdec(tindex)%vic(t)%lake_ice(lidx)/noutVar(t)
#endif
                end do

            end do ! dist

            rstdec(tindex)%vic(t)%lake_last_snow = rstdec(tindex)%vic(t)%lake_last_snow/noutVar(t)
            rstdec(tindex)%vic(t)%lake_coverage = rstdec(tindex)%vic(t)%lake_coverage/noutVar(t)
            rstdec(tindex)%vic(t)%lake_swq = rstdec(tindex)%vic(t)%lake_swq/noutVar(t)
            rstdec(tindex)%vic(t)%lake_surf_temp = rstdec(tindex)%vic(t)%lake_surf_temp/noutVar(t)
            rstdec(tindex)%vic(t)%lake_surf_water = rstdec(tindex)%vic(t)%lake_surf_water/noutVar(t)
            rstdec(tindex)%vic(t)%lake_pack_temp = rstdec(tindex)%vic(t)%lake_pack_temp/noutVar(t)
            rstdec(tindex)%vic(t)%lake_pack_water = rstdec(tindex)%vic(t)%lake_pack_water/noutVar(t)
            rstdec(tindex)%vic(t)%lake_density = rstdec(tindex)%vic(t)%lake_density/noutVar(t)
            rstdec(tindex)%vic(t)%lake_coldcontent = rstdec(tindex)%vic(t)%lake_coldcontent/noutVar(t)
            rstdec(tindex)%vic(t)%lake_snow_canopy = rstdec(tindex)%vic(t)%lake_snow_canopy/noutVar(t)

            do nidx = 1, VIC_struc(n)%Nnode
                rstdec(tindex)%vic(t)%lake_T1(nidx) = rstdec(tindex)%vic(t)%lake_T1(nidx)/noutVar(t)
            end do

            rstdec(tindex)%vic(t)%lake_activenod = rstdec(tindex)%vic(t)%lake_activenod/noutVar(t)
            rstdec(tindex)%vic(t)%lake_dz = rstdec(tindex)%vic(t)%lake_dz/noutVar(t)
            rstdec(tindex)%vic(t)%lake_surfdz = rstdec(tindex)%vic(t)%lake_surfdz/noutVar(t)
            rstdec(tindex)%vic(t)%lake_ldepth = rstdec(tindex)%vic(t)%lake_ldepth/noutVar(t)

            do node = 1,  rstdec(tindex)%vic(t)%lake_activenod
                rstdec(tindex)%vic(t)%lake_surface(node) = rstdec(tindex)%vic(t)%lake_surface(node)/noutVar(t)
            end do

            rstdec(tindex)%vic(t)%lake_sarea = rstdec(tindex)%vic(t)%lake_sarea/noutVar(t)
            rstdec(tindex)%vic(t)%lake_volume = rstdec(tindex)%vic(t)%lake_volume/noutVar(t)

            do node = 1,  rstdec(tindex)%vic(t)%lake_activenod
                rstdec(tindex)%vic(t)%lake_temp(node) = rstdec(tindex)%vic(t)%lake_temp(node)/noutVar(t)
            end do 

            rstdec(tindex)%vic(t)%lake_tempavg = rstdec(tindex)%vic(t)%lake_tempavg/noutVar(t)
            rstdec(tindex)%vic(t)%lake_areai = rstdec(tindex)%vic(t)%lake_areai/noutVar(t)
            rstdec(tindex)%vic(t)%lake_new_ice_area = rstdec(tindex)%vic(t)%lake_new_ice_area/noutVar(t)
            rstdec(tindex)%vic(t)%lake_ice_water_eq = rstdec(tindex)%vic(t)%lake_ice_water_eq/noutVar(t)
            rstdec(tindex)%vic(t)%lake_hice = rstdec(tindex)%vic(t)%lake_hice/noutVar(t)
            rstdec(tindex)%vic(t)%lake_tempi = rstdec(tindex)%vic(t)%lake_tempi/noutVar(t)
            rstdec(tindex)%vic(t)%lake_swe = rstdec(tindex)%vic(t)%lake_swe/noutVar(t)
            rstdec(tindex)%vic(t)%lake_surf_temp = rstdec(tindex)%vic(t)%lake_surf_temp/noutVar(t)
            rstdec(tindex)%vic(t)%lake_pack_temp = rstdec(tindex)%vic(t)%lake_pack_temp/noutVar(t)
            rstdec(tindex)%vic(t)%lake_coldcontent = rstdec(tindex)%vic(t)%lake_coldcontent/noutVar(t)
            rstdec(tindex)%vic(t)%lake_surf_water = rstdec(tindex)%vic(t)%lake_surf_water/noutVar(t)
            rstdec(tindex)%vic(t)%lake_pack_water = rstdec(tindex)%vic(t)%lake_pack_water/noutVar(t)
            rstdec(tindex)%vic(t)%lake_SAlbedo = rstdec(tindex)%vic(t)%lake_SAlbedo/noutVar(t)
            rstdec(tindex)%vic(t)%lake_sdepth = rstdec(tindex)%vic(t)%lake_sdepth/noutVar(t)
!        endif  ! LAKE 
#endif
    end do ! Nveg 
 end do  ! t
!------------------------------------------------------ 
! Pack state_chunks 
!------------------------------------------------------ 
 do t=1, npatch
  states = 0.0
  count = 1
!    /* write cell information */
    call pack_data(states(count), count, 1.0*rstdec(tindex)%vic(t)%tmp_cellnum)
    call pack_data(states(count), count, 1.0*rstdec(tindex)%vic(t)%tmp_Nveg)
    call pack_data(states(count), count, 1.0*rstdec(tindex)%vic(t)%tmp_Nband)

!    /* Write soil thermal node deltas */
    do nidx = 1, VIC_struc(n)%Nnode
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%dz_node(nidx)))
    end do  

!    /* Write soil thermal node depths */
    do nidx = 1, VIC_struc(n)%Nnode
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%Zsum_node(nidx)))
    end do 

!    /* Write dynamic soil properties */
#if EXCESS_ICE
!    /* Write soil depth */
    do lidx = 1, VIC_struc(n)%Nlayer
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%depth(lidx)))
    end do 

!    /* Write effective porosity */
    do lidx = 1, VIC_struc(n)%Nlayer
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%effective_porosity(lidx)))
    end do

!    /* Write damping depth */
    call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%dp))
#endif

!    /* Output for all vegetation types */
    do veg = 1, VIC_struc(n)%Nveg+1
!        // Store distributed precipitation fraction
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%mu(veg)))

!        // Store distributed precipitation variables
!        //pack_data(states(count), count, 1.0*STILL_STORM[veg]);
        if(rstdec(tindex)%vic(t)%init_STILL_STORM(veg))  then    ! .true.
            call pack_data(states(count), count, 1.0)
        else
            call pack_data(states(count), count, 0.0)
        endif

        call pack_data(states(count), count, 1.0*rstdec(tindex)%vic(t)%init_DRY_TIME(veg))

!        /* Output for all snow bands */
        do band = 1, VIC_struc(n)%Nbands
!            /* Write cell identification information */
            call pack_data(states(count), count, 1.0*rstdec(tindex)%vic(t)%tveg(veg,band))
            call pack_data(states(count), count, 1.0*rstdec(tindex)%vic(t)%tband(veg,band))

            do dist = 1, VIC_struc(n)%Ndist
!                // Store both wet and dry fractions if using distributed precipitation
!                /* Write total soil moisture */
                do lidx = 1, VIC_struc(n)%Nlayer
                    tmpval = rstdec(tindex)%vic(t)%moist(dist,veg,band,lidx)
                    call pack_data(states(count), count, real(tmpval))
                end do 

!                /* Write average ice content */
                do lidx = 1, VIC_struc(n)%Nlayer
#if SPATIAL_FROST
                    do frost_area = 1, FROST_SUBAREAS
                        tmpval = rstdec(tindex)%vic(t)%ice(dist,veg,band,lidx,frost_area)
                        call pack_data(states(count), count, real(tmpval))
                    end do 
#else
                    tmpval = rstdec(tindex)%vic(t)%ice(dist,veg,band,lidx)
                    call pack_data(states(count), count, real(tmpval))
#endif
                end do  ! lidx

!                /* Write dew storage */
                if ( veg < VIC_struc(n)%Nveg+1 ) then
                    tmpval = rstdec(tindex)%vic(t)%Wdew(dist,veg,band)
                    call pack_data(states(count), count, real(tmpval))
                endif 
            end do  ! dist

!            /* Write snow data */
            call pack_data(states(count), count, 1.0*rstdec(tindex)%vic(t)%last_snow(veg,band))
            if(rstdec(tindex)%vic(t)%MELTING(veg,band))  then    ! .true.
             call pack_data(states(count), count, 1.0)
            else
             call pack_data(states(count), count, 0.0)
            endif
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%coverage(veg,band)))
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%swq(veg,band)))
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%surf_temp(veg,band)))
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%surf_water(veg,band)))
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%pack_temp(veg,band)))
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%pack_water(veg,band)))
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%density(veg,band)))
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%coldcontent(veg,band)))
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%snow_canopy(veg,band)))

!            /* Write soil thermal node temperatures */
            do nidx = 1, VIC_struc(n)%Nnode
                call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%T1(veg,band,nidx)))
            end do 

        end do  ! band 
    end do  ! veg 

#if LAKES
        do dist = 1, VIC_struc(n)%Ndist
!            // Store both wet and dry fractions if using distributed precipitation
!            /* Write total soil moisture */
            do lidx = 1, VIC_struc(n)%Nlayer
                call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_moist(lidx)))
            end do 

!            /* Write average ice content */
            do lidx = 1, VIC_struc(n)%Nlayer
#if SPATIAL_FROST
                do frost_area = 1, FROST_SUBAREAS
                    call pack_data(states(count), count, rstdec(tindex)%vic(t)%lake_ice(lidx,frost_area))

                end do
#else
                call pack_data(states(count), count, rstdec(tindex)%vic(t)%lake_ice(lidx))
#endif
            end do   ! lidx 
        end do   ! dist
!        /* Write snow data */
        call pack_data(states(count), count, 1.0*rstdec(tindex)%vic(t)%lake_last_snow)
        read(rstdec(tindex)%vic(t)%lake_MELTING,*) fval
        call pack_data(states(count), count, fval)
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_coverage))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_swq))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_surf_temp))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_surf_water))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_pack_temp))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_pack_water))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_density))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_coldcontent))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_snow_canopy))

!        /* Write soil thermal node temperatures */
        do nidx = 1, VIC_struc(n)%Nnode
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_T1(nidx)))
        end do 

!        /* Write lake-specific variables */
        call pack_data(states(count), count, 1.0*rstdec(tindex)%vic(t)%lake_activenod)
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_dz))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_surfdz))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_ldepth))

        do node = 1, rstdec(tindex)%vic(t)%lake_activenod
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_surface(node)))
        end do  

        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_sarea))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_volume))

        do node = 1, rstdec(tindex)%vic(t)%lake_activenod
            call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_temp(node)))
        end do 

        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_tempavg))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_areai))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_new_ice_area))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_ice_water_eq))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_hice))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_tempi))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_swe))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_surf_temp))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_pack_temp))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_coldcontent))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_surf_water))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_pack_water))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_SAlbedo))
        call pack_data(states(count), count, real(rstdec(tindex)%vic(t)%lake_sdepth))

#endif

 VIC_struc(n)%state_chunk(t,:) = states(:)
 end do  ! t

    end subroutine LDT_vic412rstAvePack

    subroutine LDT_vic412rstFinalize(nfiles,npatch)

  implicit none
  integer, intent(in) :: nfiles
  integer, intent(in) :: npatch
  integer :: nf, t

!---deallocate
    do nf=1, nfiles
     do t=1, npatch
     deallocate(rstdec(nf)%vic(t)%dz_node,rstdec(nf)%vic(t)%Zsum_node)
     deallocate(rstdec(nf)%vic(t)%tveg,rstdec(nf)%vic(t)%tband)
     deallocate(rstdec(nf)%vic(t)%mu,rstdec(nf)%vic(t)%init_STILL_STORM,rstdec(nf)%vic(t)%init_DRY_TIME)
#if EXCESS_ICE
     deallocate(rstdec(nf)%vic(t)%depth,rstdec(nf)%vic(t)%effective_porosity)
#endif
     deallocate(rstdec(nf)%vic(t)%moist)
     deallocate(rstdec(nf)%vic(t)%ice)
     deallocate(rstdec(nf)%vic(t)%Wdew)
     deallocate(rstdec(nf)%vic(t)%last_snow,rstdec(nf)%vic(t)%last_snow_cnt,rstdec(nf)%vic(t)%MELTING)
     deallocate(rstdec(nf)%vic(t)%coverage,rstdec(nf)%vic(t)%swq,rstdec(nf)%vic(t)%surf_temp)
     deallocate(rstdec(nf)%vic(t)%surf_water,rstdec(nf)%vic(t)%pack_temp,rstdec(nf)%vic(t)%pack_water)
     deallocate(rstdec(nf)%vic(t)%density,rstdec(nf)%vic(t)%coldcontent,rstdec(nf)%vic(t)%snow_canopy)
     deallocate(rstdec(nf)%vic(t)%T1)
#if LAKES
     deallocate(rstdec(nf)%vic(t)%lake_moist)
     deallocate(rstdec(nf)%vic(t)%lake_ice)
     deallocate(rstdec(nf)%vic(t)%lake_T1)
     deallocate(rstdec(nf)%vic(t)%lake_surface)
     deallocate(rstdec(nf)%vic(t)%lake_temp)
#endif
     end do  ! t
    end do
    do nf=1, nfiles
     deallocate(rstdec(nf)%vic)
    end do
    deallocate(rstdec)

    end subroutine LDT_vic412rstFinalize

 subroutine unpack_state_int (statev,count,outv)
 implicit none
  real, intent(in)       :: statev
  integer, intent(inout) :: count
  integer, intent(out)   :: outv
   outv = int(statev)
   count = count + 1
   return
 end subroutine unpack_state_int
 subroutine unpack_state_char (statev,count,outv)
 implicit none
  real, intent(in)       :: statev
  integer, intent(inout) :: count
  logical*1, intent(out)   :: outv
   if ( statev == 1.0 ) then
    outv = .true.
   else
    outv = .false.
   endif
   count = count + 1
   return
 end subroutine unpack_state_char
 subroutine unpack_state_double (statev,count,outv)
 implicit none
  real, intent(in)       :: statev
  integer, intent(inout) :: count
  real*8, intent(out)   :: outv
   outv = dble(statev)
   count = count + 1
   return
 end subroutine unpack_state_double

 subroutine pack_data (ostates,count,outv)
 implicit none
  real, intent(inout):: ostates
  integer, intent(inout) :: count
  real, intent(in)   :: outv
   ostates = outv
   count = count + 1
   return
 end subroutine pack_data

 subroutine count_model_state_412(n,npatch,state_chunk_size)
!Find "state_chunk_size" that is a dimension of vic restart variable
!for binary format restart files.
!This routine is based on vic412_count_model_state.c
 implicit none
  integer, intent(in)  :: n
  integer, intent(in)  :: npatch
  integer, intent(out) :: state_chunk_size

  integer              :: count
  integer :: nidx,lidx,veg,band,dist,frost_area,node
  integer :: activenod

    count = 0   ! start with zero, to be consistent with counter in C

!    /* write cell information */
    count = count + 1 ! count_data_412(states, count, 1.0*cellnum);
    count = count + 1 ! count_data_412(states, count, 1.0*Nveg);
    count = count + 1 ! count_data_412(states, count, 1.0*Nbands);

!    /* Write soil thermal node deltas */
    do nidx = 1, VIC_struc(n)%Nnode
        count = count + 1  ! count_data_412(states, count, 1.0*soil_con->dz_node[nidx]);
    end do
    do nidx = 1, VIC_struc(n)%Nnode
        count = count + 1  ! count_data_412(states, count, 1.0*soil_con->Zsum_node[nidx]);
    end do
!    /* Write dynamic soil properties */
#if EXCESS_ICE
!    /* Write soil depth */
    do lidx = 1, VIC_struc(n)%Nlayer
        count = count + 1  !count_data_412(states, count,1.0*soil_con->depth[lidx]);
    end do
!    /* Write effective porosity */
    do lidx = 1, VIC_struc(n)%Nlayer
        count = count + 1  !count_data_412(states, count,1.0*soil_con->effective_porosity[lidx]);
    end do
!    /* Write damping depth */
    count = count + 1  !count_data_412(states, count, 1.0*soil_con->dp);
#endif

!    /* Output for all vegetation types */
    do veg = 1, VIC_struc(n)%Nveg+1
!        // Store distributed precipitation fraction
        count = count + 1  !count_data_412(states, count, 1.0*prcp->mu[veg]);

!        // Store distributed precipitation variables
        count = count + 1  ! count_data_412(states, count, 1.0*STILL_STORM[veg]);

!        /*
!        if(STILL_STORM[veg]==TRUE)
!            *count += 1; // count_data_412(states, count, 1.0);
!        else
!            *count += 1; // count_data_412(states, count, 0.0);
!        */
        count = count + 1

        count = count + 1 ! count_data_412(states, count, 1.0*DRY_TIME[veg]);

!        /* Output for all snow bands */
        do band = 1, VIC_struc(n)%Nbands
!            /* Write cell identification information */
            count = count + 1 ! count_data_412(states, count, 1.0*veg);
            count = count + 1 ! count_data_412(states, count, 1.0*band);

            do dist = 1, VIC_struc(n)%Ndist
!                // Store both wet and dry fractions if using distributed precipitation

!                /* Write total soil moisture */
                do lidx = 1, VIC_struc(n)%Nlayer
!                  tmpval = cell[dist][veg][band].layer[lidx].moist;
                   count = count + 1 !count_data_412(states, count, 1.0*tmpval);
                end do 

!                /* Write average ice content */
                do lidx = 1, VIC_struc(n)%Nlayer
#if SPATIAL_FROST
                    do frost_area = 1, FROST_SUBAREAS
!                    tmpval =cell[dist][veg][band].layer[lidx].ice[frost_area];
                        count = count + 1 !count_data_412(states, count,1.0*tmpval);
                    end do 
#else
!                   tmpval = cell[dist][veg][band].layer[lidx].ice;
                    count = count + 1 !count_data_412(states, count, 1.0*tmpval);
#endif
                end do ! lidx 

!                /* Write dew storage */
                if ( veg < VIC_struc(n)%Nveg+1 ) then
!                    tmpval = veg_var[dist][veg][band].Wdew;
                    count = count + 1 ! count_data_412(states, count, 1.0*tmpval);
                endif 
            end do  ! dist

!            /* Write snow data */
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].last_snow);
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].MELTING);
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].coverage);
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].swq);
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].surf_temp);
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].surf_water);
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].pack_temp);
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].pack_water);
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].density);
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].coldcontent);
            count = count + 1 ! count_data_412(states, count,1.0*snow[veg][band].snow_canopy);

!            /* Write soil thermal node temperatures */
            do nidx = 1, VIC_struc(n)%Nnode
                count = count + 1 ! count_data_412(states, count,1.0*energy[veg][band].T[nidx]);
            end do

        end do  ! Nband 
    end do ! Nveg 

#if LAKES
        do dist = 1, VIC_struc(n)%Ndist
!            // Store both wet and dry fractions if using distributed precipitation
!            /* Write total soil moisture */
            do lidx = 1, VIC_struc(n)%Nlayer
                count = count + 1 ! count_data_412(states, count,1.0*lake_var.soil.layer[lidx].moist);
            end do 
!            /* Write average ice content */
            do lidx = 1, VIC_struc(n)%Nlayer
#if SPATIAL_FROST
                do frost_area = 1, FROST_SUBAREAS
                    count = count + 1 ! count_data_412(states, count,1.0*lake_var.soil.layer[lidx].ice[frost_area]);
                end do 
#else
                count = count + 1 ! count_data_412(states, count,1.0*lake_var.soil.layer[lidx].ice);
#endif
            end do 
        end do  ! Ndist 
!        /* Write snow data */
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.snow.last_snow);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.snow.MELTING);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.snow.coverage);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.snow.swq);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.snow.surf_temp);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.snow.surf_water);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.snow.pack_temp);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.snow.pack_water);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.snow.density);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.snow.coldcontent);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.snow.snow_canopy);

!        /* Write soil thermal node temperatures */
        do nidx = 1, VIC_struc(n)%Nnode
            count = count + 1 ! count_data_412(states, count,1.0*lake_var.energy.T[nidx]);
        end do 

!        /* Write lake-specific variables */
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.activenod);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.dz);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.surfdz);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.ldepth);

!hkb lake_var.activenod needs to be read in if LAKES option enabled!!!!
        do node = 1, activenod+1
            count = count + 1 ! count_data_412(states, count,1.0*lake_var.surface[node]);
        end do 

        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.sarea);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.volume);

        do node = 1, activenod
            count = count + 1 ! count_data_412(states, count,1.0*lake_var.temp[node]);
        end do 

        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.tempavg);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.areai);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.new_ice_area);
        count = count + 1 ! count_data_412(states, count,1.0*lake_var.ice_water_eq);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.hice);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.tempi);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.swe);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.surf_temp);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.pack_temp);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.coldcontent);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.surf_water);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.pack_water);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.SAlbedo);
        count = count + 1 ! count_data_412(states, count, 1.0*lake_var.sdepth);
#endif
 state_chunk_size = count
 return 
 end subroutine count_model_state_412

end module LDT_vic412rstMod
