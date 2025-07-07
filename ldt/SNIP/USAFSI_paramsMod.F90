!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!*******************************************************************************
!*******************************************************************************
!**
!**  NAME: PARAMETERS
!**
!**  PURPOSE: HOLDS PARAMETER STATEMENTS FOR USAFSI
!**
!**  FILES ACCESSED: NONE
!**
!**  UPDATES
!**  =======
!**  26 APR 95  INITIAL UNISYS VERSION...........................SSGT CONRY/SYSM
!**  08 FEB 96  ADDED ELVCHK, ELVCOR, AND NUMCHL.................SSGT CONRY/SYSM
!**  23 OCT 97  ADDED CLOBMX.  MOVED CLOSOB, ELVCHK, ELVCOR, AND ROI TO
!**             TUNES COMMON BLOCK...............................SSGT CONRY/DNXM
!**  05 MAY 98  CHANGED NUMCHL TO 7.................................DR KOPP/DNXM
!**  26 FEB 01  PORTED FROM UNISYS MAINFRAME TO UNIX............SSGT MILLER/DNXM
!**  03 MAR 03  REMOVED VARIABLE CLOBMX, NOT NEEDED................MR GAYNO/DNXM
!**  17 MAY 04  CONVERTED TO FORTRAN 90 MODULE AND ADDED NT2ICE
!**             PARAMETERS......................................MR EYLANDER/DNXM
!**  15 APR 05  CHANGED ICEPNT & MAXAGE TO INT AND ICEPNT TO 1.
!**             ADDED GRIDSIZE, MINSNO, NOCOSNO, & SNOTHRESH....MR LEWISTON/DNXM
!**  09 MAY 05  CHANGED MISVAL FROM -99999998 TO -9999..........MR LEWISTON/DNXM
!**  20 MAY 05  ADDED MISANL....................................MR LEWISTON/DNXM
!**  27 AUG 08  CHANGED MISVAL BACK TO -99999998............MR LEWISTON/2WXG/WEA
!**  11 JUN 09  ADDED 1/4 DEGREE EQUIDISTANT CYLIND PARAMS..MR LEWISTON/2WXG/WEA
!**  09 SEP 09  REMOVED SEA ICE AND SDR PARAMETERS..........MR LEWISTON/2WXG/WEA
!**  13 MAY 10  REVERSED LATITIUDE ORIENTATION..............MR LEWISTON/16WS/WXE
!**  16 SEP 10  REMOVED NT2ICE PARAMETERS...................MR LEWISTON/16WS/WXE
!**  29 DEC 10  ADDED MESHNAME..............................MR LEWISTON/16WS/WXE
!**  10 FEB 11  ADDED MESHDEG FOR OBS SPREADING.............MR LEWISTON/16WS/WXE
!**  12 SEP 11  REMOVED MAXSWE..............................MR LEWISTON/16WS/WXE
!**  14 DEC 12  REMOVED RECLEN_LIS..........................MR LEWISTON/16WS/WXE
!**  22 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
!**  09 Mar 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
!**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
!**
!*******************************************************************************
!*******************************************************************************

#include "LDT_misc.h"

module USAFSI_paramsMod

   ! Defaults
   implicit none
   public

   character*12, parameter     :: program_name = 'USAFSI'        ! NAME OF MAIN PROGRAM
   
   integer, parameter    :: msglns      = 20                     ! MAXIMUM NUMBER OF LINES IN ERROR MESSAGE

   ! Parameters for reading legacy 0.25 deg SNODEP product
   character*8,  parameter     :: meshname = '_0p25deg'          ! GRID AND MESH FOR DATA FILE NAMES
   integer, parameter    :: mesh        = 4                      ! EQUIDISTANT CYL GRID MESH (4 = 1/4 DEGREE)
   integer, parameter    :: igrid       = 360 * mesh             ! SIZE OF GRID IN THE I-DIRECTION
   integer, parameter    :: jgrid       = 180 * mesh             ! SIZE OF GRID IN THE J-DIRECTION   
   real,    parameter    :: begin_lat   = -89.875                ! LATITUDE OF LOWER LEFT CORNER OF GRID
   real,    parameter    :: begin_lon   = -179.875               ! LONGITUDE OF LOWER LEFT CORNER OF GRID
   
   ! Parameter for reading legacy 0.25 deg degribbed LIS product
   integer, parameter    :: jgrid_lis   = 150 * mesh             ! SIZE OF LIS GRID IN J-DIRECTION
   
   ! Parameters for reading legacy 0.25 deg US Navy SST product
   integer, parameter    :: sst_igrid   = 1440                   ! SIZE OF SST GRID IN THE I-DIRECTION
   integer, parameter    :: sst_jgrid   = 721                    ! SIZE OF SST GRID IN THE J-DIRECTION
   
   ! PARAMETERS NEEDED BY SNOW DEPTH ALGORITHM.
   integer, parameter    :: arctlat     = 6650                   ! LATITUDE OF ARCTIC CIRCLE (*100)
   real,    parameter    :: misanl      = -1.0                   ! MISSING ANALYSIS FLAG   
   integer, parameter    :: icepnt      = 1                      ! ICE FLAG (0 = NO ICE, 1 = ICE)
   integer, parameter    :: maxage      = 365                    ! MAXIMUM SNOW AGE IN DAYS
   integer, parameter    :: missing     = -1                     ! MISSING VALUE FLAG
   integer, parameter    :: misval      = -99999998              ! CDMS MISSING VALUE FLAG   
   real,    parameter    :: snothresh   = 0.02                   ! THRESHOLD VALUE FOR TRACKING SNOW DEPTH

   integer, parameter :: icedef = 95 ! Bogus value for ice concentration
end module USAFSI_paramsMod
