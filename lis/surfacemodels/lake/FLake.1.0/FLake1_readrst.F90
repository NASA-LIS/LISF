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
!BOP
!
! !ROUTINE: FLake1_readrst
! \label{FLake1_readrst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   6/4/13: Shugong Wang; initial implementation for LIS 7 and FLake1
!
! !INTERFACE:
subroutine FLake1_readrst()
! !USES:
  use LIS_coreMod, only    : LIS_rc, LIS_masterproc
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
       LIS_getNextUnitNumber,   &
       LIS_releaseUnitNumber,   &
       LIS_verify                
  use FLake1_Mod
  
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for FLake1.  This
!  includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the FLake1
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    T_snow                     - FLake1 temperature at the air-snow interface [K]
!    T_ice                      - FLake1 temperature at the snow-ice interface [K]
!    T_mnw                      - FLake1 mean temperature of the water column [K]
!    T_wML                      - FLake1 temperature of mixed layer [K]
!    T_bot                      - FLake1 temperature at the water-bottom sediment interface [K]
!    T_b1                       - FLake1 temperature at the bottom of the upper layer of the sediments [K]
!    C_T                        - FLake1 thermocline shape factor [-]
!    H_snow                     - FLake1 snow thickness [m]
!    H_ice                      - FLake1 ice thickness [m]
!    H_ML                       - FLake1 thickness of mixed layer [m]
!    H_B1                       - FLake1 thickness of the upper layer of bottom sediments [m]
!    T_sfc                      - FLake1 surface temperature [K]
!    albedo_water               - FLake1 water surface albedo with resect to solar radiation [-]
!    albedo_ice                 - FLake1 ice surface albedo with respect to the solar radiation [-]
!    albedo_snow                - FLake1 snow surface albedo with respect to the solar radiation [-]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!      reads a variable from the restart file
!   \item[FLake1\_coldstart](\ref{FLake1_coldstart}) \newline
!      initializes the FLake1 state variables
! \end{description}
!EOP
 
  implicit none
  
  integer           :: t, l
  integer           :: nc, nr, npatch
  integer           :: n
  integer           :: ftn
  integer           :: status
  real, allocatable :: tmptilen(:)
  logical           :: file_exists
  character*20      :: wformat
  
  do n=1, LIS_rc%nnest
     wformat = trim(FLAKE1_struc(n)%rformat)
     ! coldstart
     if(LIS_rc%startcode .eq. "coldstart") then  
        call FLake1_coldstart(LIS_rc%lake_index)
        ! restart
     elseif(LIS_rc%startcode .eq. "restart") then
        allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lake_index)))
        ! check the existance of restart file
        inquire(file=FLAKE1_struc(n)%rfile, exist=file_exists)
        If (.not. file_exists) then 
           write(LIS_logunit,*) "FLake1 restart file ", FLAKE1_struc(n)%rfile," does not exist "
           write(LIS_logunit,*) "Program stopping ..."
           call LIS_endrun
        endif
        write(LIS_logunit,*) "FLake1 restart file used: ", FLAKE1_struc(n)%rfile
        
        ! open restart file
        if(wformat .eq. "binary") then
           ftn = LIS_getNextUnitNumber()
           open(ftn, file=FLAKE1_struc(n)%rfile, form="unformatted")
           read(ftn) nc, nr, npatch  !time, veg class, no. tiles
           
           ! check for grid space conflict
           if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
              write(LIS_logunit,*) FLAKE1_struc(n)%rfile, "grid space mismatch - FLake1 halted"
              call LIS_endrun
           endif
           
           if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lake_index)) then
              write(LIS_logunit,*) "restart tile space mismatch, halting..."
              call LIS_endrun
           endif
        elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_open(path=FLAKE1_struc(n)%rfile, &
                mode=NF90_NOWRITE, ncid=ftn)
           call LIS_verify(status, "Error opening file "//FLAKE1_struc(n)%rfile)
#endif
        endif
        
        ! read: temperature at the air-snow interface
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_snow, &
             varname="T_SNOW", wformat=wformat)
        
        ! read: temperature at the snow-ice interface
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_ice, &
             varname="T_ICE", wformat=wformat)
        
        ! read: mean temperature of the water column
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_mnw, &
             varname="T_MNW", wformat=wformat)
        
        ! read: temperature of mixed layer
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_wML, &
             varname="T_WML", wformat=wformat)
        
        ! read: temperature at the water-bottom sediment interface
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_bot, &
             varname="T_BOT", wformat=wformat)
        
        ! read: temperature at the bottom of the upper layer of the sediments
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_b1, &
             varname="T_B1", wformat=wformat)
        
        ! read: thermocline shape factor
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%C_T, &
             varname="C_T", wformat=wformat)
        
        ! read: snow thickness
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%H_snow, &
             varname="H_SNOW", wformat=wformat)
        
        ! read: ice thickness
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%H_ice, &
             varname="H_ICE", wformat=wformat)
        
        ! read: thickness of mixed layer
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%H_ML, &
             varname="H_ML", wformat=wformat)
        
        ! read: thickness of the upper layer of bottom sediments
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%H_B1, &
             varname="H_B1", wformat=wformat)
        
        ! read: surface temperature
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_sfc, &
             varname="T_SFC", wformat=wformat)
        
        ! read: water surface albedo with resect to solar radiation
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%albedo_water, &
             varname="ALBEDO_WATER", wformat=wformat)
        
        ! read: ice surface albedo with respect to the solar radiation
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%albedo_ice, &
             varname="ALBEDO_ICE", wformat=wformat)
        
        ! read: snow surface albedo with respect to the solar radiation
        call LIS_readvar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%albedo_snow, &
             varname="ALBEDO_SNOW", wformat=wformat)
        
        ! close restart file
        if(wformat .eq. "binary") then
           call LIS_releaseUnitNumber(ftn)
        elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status, "Error in nf90_close in FLake1_readrst")
#endif
        endif
        deallocate(tmptilen)
     endif
  enddo
end subroutine FLake1_readrst
        
