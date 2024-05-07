!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!*****************************************************************************************
!*****************************************************************************************
!**
!**  NAME: USAFSI_ARRAYS
!**
!**  PURPOSE: HOLDS USAFSI MODEL SHARED ARRAYS
!**
!**  UPDATES
!**  =======
!**  21 JUL 04  INITIAL VERSION...........................................MR EYLANDER/DNXM
!**  14 APR 05  SEPARATED SNOW AND ICE INTO SEPARATE ARRAYS.
!**             REMOVED WATERPTS (USED FOR GRIBBING)......................MR LEWISTON/DNXM
!**  08 MAY 09  ADDED AMSR-E SNOW. REMOVED SEA ICE....................MR LEWISTON/2WXG/WEA
!**  01 JUL 09  ADDED AMSR-E STATIC FILES.............................MR LEWISTON/2WXG/WEA
!**  05 JAN 10  CHANGED SSMIS FROM SDRS TO EDRS.......................MR LEWISTON/16WS/WXE
!**  15 SEP 10  DELETED SATELLITE SNOW FLAGS..........................MR LEWISTON/16WS/WXE
!**  08 NOV 11  REMOVED AMSR-E ARRAYS AND PORTED TO LINUX.............MR LEWISTON/16WS/WXE
!**  11 SEP 12  REMOVED SNOW_DENSITY AND SNOW WATER EQUIVALENT (SWE)..MR LEWISTON/16WS/WXE
!**  07 NOV 12  ADDED FRACTIONAL SNOW.................................MR LEWISTON/16WS/WXE
!**  10 OCT 13  ADDED ICEAGE12Z AND SNOAGE12Z.........................MR LEWISTON/16WS/WXE
!**  15 FEB 17  ADDED VIIRS DATA........................................MR PUSKAR/16WS/WXE
!**  22 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
!**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
!**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
!** 
!*****************************************************************************************
!*****************************************************************************************

#include "LDT_misc.h"

module USAFSI_arraysMod

   ! Defaults
   implicit none
   private

   type USAFSI_arrays_t
      integer*1,  allocatable    :: snow_poss        ( : , : )    ! SNOW POSSIBLE MASK (0=NO SNOW; 1=SNOW)
      integer,    allocatable    :: iceage           ( : , : )    ! ICE AGE (DAYS)
      integer,    allocatable    :: iceage12z        ( : , : )    ! ICE AGE (DAYS) PROM PREVIOUS 12Z
      integer,    allocatable    :: icecon           ( : , : )    ! ICE CONCENTRATIONS (PERCENT, 0-100)
      integer,    allocatable    :: icemask          ( : , : )    ! ICE MASK (0 = NO ICE; 1 = ICE)
      integer,    allocatable    :: oldcon           ( : , : )    ! PREVIOUS DAY'S ICE CONCENTRATIONS
      integer,    allocatable    :: oldmask          ( : , : )    ! PREVIOUS DAY'S ICE FLAGS
      integer,    allocatable    :: ssmis_icecon     ( : , : )    ! SSMIS ICE CONCENTRATIONS (PERCENT, 0-100)
      integer,    allocatable    :: snoage           ( : , : )    ! SNOW AGE (DAYS)
      integer,    allocatable    :: snoage12z        ( : , : )    ! SNOW AGE (DAYS) FROM PREVIOUS 12Z
      integer,    allocatable    :: viirsmap         ( : , : )    ! VIIRS SNOW COVERED AREA MAP   
      real,       allocatable    :: climo            ( : , : )    ! THIS MONTH'S SNOW CLIMATOLOGY (METERS)
      real,       allocatable    :: elevat           ( : , : )    ! TERRAIN ELEVATION (METERS)
      real,       allocatable    :: olddep           ( : , : )    ! PREVIOUS DAY'S SNOW ANALYSIS (METERS)
      real,       allocatable    :: ptlat            ( : , : )    ! GRID POINT LATITUDES
      real,       allocatable    :: ptlon            ( : , : )    ! GRID POINT LONGITUDES
      real,       allocatable    :: snoanl           ( : , : )    ! CURRENT SNOW DEPTH ANALYSIS (METERS)
      real,       allocatable    :: snofrac          ( : , : )    ! FRACTIONAL SNOW DATA ON USAFSI GRID
      real,       allocatable    :: ssmis_depth      ( : , : )    ! SNOW DEPTH FROM SSMIS EDRS
      real,       allocatable    :: sst              ( : , : )    ! NAVY SEA SURFACE TEMPERATURES (KELVIN)
      real, allocatable :: gofs_icecon(:,:)
   end type USAFSI_arrays_t

   type(USAFSI_arrays_t), public :: USAFSI_arrays

end module USAFSI_arraysMod
