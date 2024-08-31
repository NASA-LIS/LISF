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
! !ROUTINE: FLake1_writerst
! \label{FLake1_writerst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   6/4/13: Shugong Wang; initial implementation for LIS 7 and FLake1
!
! !INTERFACE:
subroutine FLake1_writerst(n)
! !USES:
  use LIS_coreMod, only    : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber , LIS_verify
  use LIS_fileIOMod, only  : LIS_create_output_directory, &
       LIS_create_restart_filename
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use FLake1_Mod
  
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for FLake1.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[FLake1\_dump\_restart](\ref{FLake1_dump_restart}) \newline
!   writes the FLake1 variables into the restart file
! \end{description}
!EOP

  character(len=LIS_CONST_PATH_LEN) :: filen
  character*20  :: wformat
  logical       :: alarmCheck
  integer       :: ftn
  integer       :: status
  character*3   :: fnest

  write(fnest,'(i3.3)') n
  ! set restart alarm
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "FLAKE1 restart alarm"//trim(fnest))

  ! set restart file format (read from LIS configration file_
  wformat = trim(FLAKE1_struc(n)%rformat)

  if(alarmCheck .or. (LIS_rc%endtime ==1)) then
     If (LIS_masterproc) Then
        call LIS_create_output_directory("SURFACEMODEL")
        call LIS_create_restart_filename(n, filen, "SURFACEMODEL", &
             "FLAKE1",&
             wformat=wformat)
        if(wformat .eq. "binary") then
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=filen,status="unknown", form="unformatted")
        elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
           status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
           call LIS_verify(status,"Error in nf90_open in FLake1_writerst")
#endif
#if (defined USE_NETCDF3)
           status = nf90_create(Path = filen, cmode = nf90_clobber, &
                ncid = ftn)
           call LIS_verify(status, "Error in nf90_open in FLake1_writerst")
#endif
        endif
     endif

     call FLake1_dump_restart(n, ftn, wformat)

     if (LIS_masterproc) then
        if(wformat .eq. "binary") then
           call LIS_releaseUnitNumber(ftn)
        elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status, "Error in nf90_close in FLake1_writerst")
#endif
        endif
        write(LIS_logunit, *) "FLake1 archive restart written: ", trim(filen)
     endif
  endif
end subroutine FLake1_writerst

!BOP
!
! !ROUTINE: FLake1_dump_restart
! \label{FLake1_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  6/4/13: Shugong Wang, initial implementation for LIS 7 and FLake1
! !INTERFACE:
subroutine FLake1_dump_restart(n, ftn, wformat)

! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_logMod, only  : LIS_logunit
  use LIS_historyMod
  use FLake1_Mod
  
  implicit none
  
  integer, intent(in) :: ftn
  integer, intent(in) :: n
  character(len=*), intent(in) :: wformat
!
! !DESCRIPTION:
!  This routine gathers the necessary restart variables and performs
!  the actual write statements to create the restart files.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    unit number for the restart file
!   \item[wformat]
!    restart file format (binary/netcdf)
!  \end{description}
!
!
!  The following is the list of variables written in the FLake1
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
! The routines invoked are:
! \begin{description}
!   \item[LIS\_writeGlobalHeader\_restart](\ref{LIS_writeGlobalHeader_restart}) \newline
!      writes the global header information
!   \item[LIS\_writeHeader\_restart](\ref{LIS_writeHeader_restart}) \newline
!      writes the header information for a variable
!   \item[LIS\_closeHeader\_restart](\ref{LIS_closeHeader_restart}) \newline
!      close the header
!   \item[LIS\_writevar\_restart](\ref{LIS_writevar_restart}) \newline
!      writes a variable to the restart file
! \end{description}
! 
!EOP 
  
  integer :: l, t 
  real    :: tmptilen(LIS_rc%npatch(n, LIS_rc%lake_index))
  integer :: dimID(11)
  integer :: T_snow_ID
  integer :: T_ice_ID
  integer :: T_mnw_ID
  integer :: T_wML_ID
  integer :: T_bot_ID
  integer :: T_b1_ID
  integer :: C_T_ID
  integer :: H_snow_ID
  integer :: H_ice_ID
  integer :: H_ML_ID
  integer :: H_B1_ID
  integer :: T_sfc_ID
  integer :: albedo_water_ID
  integer :: albedo_ice_ID
  integer :: albedo_snow_ID

  ! write the header of the restart file
  call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lake_index, &
       "FLAKE1", dim1=4, dim2=4, dimID=dimID, &
       output_format = trim(wformat))

  ! write the header for state variable T_snow
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, T_snow_ID, "T_SNOW", &
       "temperature at the air-snow interface", &
       "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable T_ice
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, T_ice_ID, "T_ICE", &
       "temperature at the snow-ice interface", &
       "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable T_mnw
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, T_mnw_ID, "T_MNW", &
       "mean temperature of the water column", &
       "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable T_wML
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, T_wML_ID, "T_WML", &
       "temperature of mixed layer", &
       "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable T_bot
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, T_bot_ID, "T_BOT", &
       "temperature at the water-bottom sediment interface", &
       "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable T_b1
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, T_b1_ID, "T_B1", &
       "temperature at the bottom of the upper layer of the sediments", &
       "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable C_T
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, C_T_ID, "C_T", &
       "thermocline shape factor", &
       "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable H_snow
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, H_snow_ID, "H_SNOW", &
       "snow thickness", &
       "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable H_ice
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, H_ice_ID, "H_ICE", &
       "ice thickness", &
       "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable H_ML
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, H_ML_ID, "H_ML", &
       "thickness of mixed layer", &
       "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable H_B1
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, H_B1_ID, "H_B1", &
       "thickness of the upper layer of bottom sediments", &
       "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable T_sfc
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, T_sfc_ID, "T_SFC", &
       "surface temperature", &
       "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable albedo_water
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, albedo_water_ID, "ALBEDO_WATER", &
       "water surface albedo with resect to solar radiation", &
       "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable albedo_ice
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, albedo_ice_ID, "ALBEDO_ICE", &
       "ice surface albedo with respect to the solar radiation", &
       "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! write the header for state variable albedo_snow
  !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
  call LIS_writeHeader_restart(ftn, n, dimID, albedo_snow_ID, "ALBEDO_SNOW", &
       "snow surface albedo with respect to the solar radiation", &
       "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

  ! close header of restart file
  call LIS_closeHeader_restart(ftn, n, LIS_rc%lake_index, dimID, FLAKE1_struc(n)%rstInterval)

  ! write state variables into restart file
  ! temperature at the air-snow interface
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_snow, &
       varid=T_snow_ID, dim=1, wformat=wformat)

  ! temperature at the snow-ice interface
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_ice, &
       varid=T_ice_ID, dim=1, wformat=wformat)

  ! mean temperature of the water column
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_mnw, &
       varid=T_mnw_ID, dim=1, wformat=wformat)

  ! temperature of mixed layer
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_wML, &
       varid=T_wML_ID, dim=1, wformat=wformat)

  ! temperature at the water-bottom sediment interface
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_bot, &
       varid=T_bot_ID, dim=1, wformat=wformat)

  ! temperature at the bottom of the upper layer of the sediments
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_b1, &
       varid=T_b1_ID, dim=1, wformat=wformat)

  ! thermocline shape factor
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%C_T, &
       varid=C_T_ID, dim=1, wformat=wformat)

  ! snow thickness
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%H_snow, &
       varid=H_snow_ID, dim=1, wformat=wformat)

  ! ice thickness
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%H_ice, &
       varid=H_ice_ID, dim=1, wformat=wformat)

  ! thickness of mixed layer
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%H_ML, &
       varid=H_ML_ID, dim=1, wformat=wformat)

  ! thickness of the upper layer of bottom sediments
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%H_B1, &
       varid=H_B1_ID, dim=1, wformat=wformat)

  ! surface temperature
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%T_sfc, &
       varid=T_sfc_ID, dim=1, wformat=wformat)

  ! water surface albedo with resect to solar radiation
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%albedo_water, &
       varid=albedo_water_ID, dim=1, wformat=wformat)

  ! ice surface albedo with respect to the solar radiation
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%albedo_ice, &
       varid=albedo_ice_ID, dim=1, wformat=wformat)

  ! snow surface albedo with respect to the solar radiation
  call LIS_writevar_restart(ftn, n, LIS_rc%lake_index, FLAKE1_struc(n)%flake1%albedo_snow, &
       varid=albedo_snow_ID, dim=1, wformat=wformat)

end subroutine FLake1_dump_restart

