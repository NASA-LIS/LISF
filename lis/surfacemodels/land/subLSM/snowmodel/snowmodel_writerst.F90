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
! !ROUTINE: snowmodel_writerst
! \label{snowmodel_writerst}
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel 
!  05 Aug 2022: Kristi Arsenault; Updated SnowModel state restart writer
!
! !INTERFACE:
subroutine snowmodel_writerst(n)
! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod,     only : LIS_logunit, LIS_getNextUnitNumber, &
                             LIS_releaseUnitNumber, LIS_verify
  use LIS_fileIOMod,  only : LIS_create_output_directory, &
                             LIS_create_restart_filename
  use snowmodel_lsmMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for SnowModel.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[snowmodel\_dump\_restart](\ref{snowmodel_dump_restart}) \newline
!   writes the SnowModel variables into the restart file
! \end{description}
!EOP
  character(len=LIS_CONST_PATH_LEN) :: filen
  logical       :: alarmCheck
  integer       :: ftn
  integer       :: status
  character*20  :: wformat
  character*3   :: fnest

  write(fnest,'(i3.3)') n
  ! set restart alarm
  alarmCheck = LIS_isAlarmRinging(LIS_rc,"SnowModel restart alarm "//trim(fnest))

  write(LIS_logunit,*) '[INFO] Call to the SnowModel write restart routine ...'

  ! set restart file format (read from LIS configration file_
  wformat = trim(snowmodel_struc(n)%rformat)

  if(alarmCheck .or. (LIS_rc%endtime ==1)) then
     if (LIS_masterproc) Then
        call LIS_create_output_directory("SURFACEMODEL")
        call LIS_create_restart_filename(n, filen, "SURFACEMODEL", &
                                         "SNOWMODEL",wformat=wformat)
        if(wformat .eq. "binary") then
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=filen,status="unknown", form="unformatted")
        elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
           status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
           call LIS_verify(status, &
                "Error in nf90_open in snowmodel_writerst")
#endif
#if (defined USE_NETCDF3)
           status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
           call LIS_verify(status, &
                "Error in nf90_open in snowmodel_writerst")
#endif
        endif
     endif

     ! Call routine to write out the model states to the file:
     call SnowModel_dump_restart(n, ftn, wformat)

     if (LIS_masterproc) then
        if(wformat .eq. "binary") then
           call LIS_releaseUnitNumber(ftn)
        elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status, &
                "Error in nf90_close in snowmodel_writerst")
#endif
        endif
        write(LIS_logunit,*)&
             "[INFO] SnowModel archive restart written: ",trim(filen)
     endif
  endif
end subroutine snowmodel_writerst

!BOP
!
! !ROUTINE: SnowModel_dump_restart
! \label{SnowModel_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!  05 Aug 2022: Kristi Arsenault; Updated SnowModel state restart writer
!
! !INTERFACE:
subroutine SnowModel_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use snowmodel_lsmMod

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
!  The following is the list of variables written in the SnowModel
!  restart file:
!
!  \begin{verbatim}
!    nc, nr, ntiles        - grid and tile space dimensions
!    snow_d                - Snow depth (m), density-adjusted and used in subgrid-scale snow
!    snow_depth            - Snow depth (m)
!    canopy_int            - Canopy interception store (m)
!    soft_snow_d           - Soft snow layer depth (m)
!    ro_snow_grid          - Snow grid density (kg/m3)
!    swe_depth             - Snow water equivalent depth (m)
!    ro_soft_snow_old      - Density of former soft snow layer (kg/m3)
!    snow_d_init           - Initial snow depth (m)
!    swe_depth_old         - Former SWE depth (m) step
!    canopy_int_old        - Former canopy interception store (m)
!    topo                  - Snow-depth changing grid topography level (m)
!    sum_sprec             - Summed snowfall (m)
!  \end{verbatim}
!
! The routines invoked are:
! \begin{description}
!   \item[LIS\_writeGlobalHeader\_restart](\ref{LIS_writeGlobalHeader_restart})\\
!      writes the global header information
!   \item[LIS\_writeHeader\_restart](\ref{LIS_writeHeader_restart})\\
!      writes the header information for a variable
!   \item[LIS\_closeHeader\_restart](\ref{LIS_closeHeader_restart})\\
!      close the header
!   \item[LIS\_writevar\_restart](\ref{LIS_writevar_restart})\\
!      writes a variable to the restart file
! \end{description}
! 
!EOP 

    integer :: l, t
    real    :: tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index))
    integer :: dimID(11)

    integer :: snow_d_ID          ! 1
    integer :: snow_depth_ID      ! 2
    integer :: canopy_int_ID      ! 3
    integer :: softsnow_d_ID      ! 4
    integer :: ro_snowgrid_ID     ! 5
    integer :: swe_depth_ID       ! 6
    integer :: ro_softsnow_old_ID ! 7
    integer :: snow_d_init_ID     ! 8
    integer :: swe_depth_old_ID   ! 9
    integer :: canopy_int_old_ID  ! 10
    integer :: topo_ID            ! 11
    integer :: sum_sprec_ID       ! 12


!- Write the header of the restart file
   call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, &
                                       "SNOWMODEL", &
                                       dim1=snowmodel_struc(n)%nsnow, &   ! Set as one for now
                                       dimID=dimID, &
                                       output_format = trim(wformat))

   !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max

!-  SnowModel states, based on preprocess_code: HRESTART_SAVE routine

   ! write header for Snow depth (m), density-adjusted and used in subgrid-scale snow (1)
   call LIS_writeHeader_restart(ftn, n, dimID, snow_d_ID, "SNOWD", &
                                "Snow depth, density-adjusted and used in SnowTran subgrid-scale snow", &
                                "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- snow depth (2)
   call LIS_writeHeader_restart(ftn, n, dimID, snow_depth_ID, "SNOWDEPTH", &
                                "Snow depth", &
                                "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- canopy interception store (3)
   call LIS_writeHeader_restart(ftn, n, dimID, canopy_int_ID, "CANOPYINT", &
                                "Canopy interception storage term", &
                                "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- soft snow depth layer (4)
   call LIS_writeHeader_restart(ftn, n, dimID, softsnow_d_ID, "SOFTSNOWD", &
                                "Soft snow layer depth", &
                                "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- snow grid density (5)
   call LIS_writeHeader_restart(ftn, n, dimID, ro_snowgrid_ID, "ROSNOWGRID", &
                                "Snow grid density", &
                                "kg/m3", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- swe depth (6)
   call LIS_writeHeader_restart(ftn, n, dimID, swe_depth_ID, "SWEDEPTH", &
                                "Snow water equivalent depth", &
                                "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- density of former soft snow layer (7)
   call LIS_writeHeader_restart(ftn, n, dimID, ro_softsnow_old_ID, "ROSOFTSNOWOLD", &
                                "Density of former soft snow layer", &
                                "kg/m3", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- Initial snow depth, density-layer and subgrid scale snow (8)
   call LIS_writeHeader_restart(ftn, n, dimID, snow_d_init_ID, "SNOWDINIT", &
                                "Initial snow depth, density-layer and subgrid scale snow", &
                                "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- former SWE depth step (9)
   call LIS_writeHeader_restart(ftn, n, dimID, swe_depth_old_ID, "SWEDEPTHOLD", &
                                "Former SWE depth step", &
                                "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- Former canopy interception store (10)
   call LIS_writeHeader_restart(ftn, n, dimID, canopy_int_old_ID, "CANOPYINTOLD", &
                                "Former canopy interception store", &
                                "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- Snow-depth changing grid topography level (11) 
   call LIS_writeHeader_restart(ftn, n, dimID, topo_ID, "TOPO", &
                                "Snow-depth changing grid topography level", &
                                "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write the header for state variable -- summed snowfall (12) 
   call LIS_writeHeader_restart(ftn, n, dimID, sum_sprec_ID, "SUMSPREC", &
                                "Sum of (frozen) snow precipitation", &
                                "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! close header of restart file
   call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, snowmodel_struc(n)%rstInterval)


!- Write state variables into restart file:

   ! (1) Snow depth (m), density-adjusted and used in subgrid-scale snow
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%snow_d, &
                             varid=snow_d_ID, dim=1, wformat=wformat)

   ! (2) Snow depth (m)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%snow_depth, &
                             varid=snow_depth_ID, dim=1, wformat=wformat)

   ! (3) Canopy interception store (m)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%canopy_int, &
                             varid=canopy_int_ID, dim=1, wformat=wformat)

   ! (4) Soft snow depth layer (m)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%soft_snow_d, &
                             varid=softsnow_d_ID, dim=1, wformat=wformat)

   ! (5) Snowdepth grid density (kg/m3)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%ro_snow_grid, &
                             varid=ro_snowgrid_ID, dim=1, wformat=wformat)

   ! (6) SWE depth (m)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%swe_depth, &
                             varid=swe_depth_ID, dim=1, wformat=wformat)

   ! (7) Former soft snow density step (kg/m3)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%ro_soft_snow_old, &
                             varid=ro_softsnow_old_ID, dim=1, wformat=wformat)

   ! (8) Snowdepth init (m)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%snow_d_init, &
                             varid=snow_d_init_ID, dim=1, wformat=wformat)

   ! (9) Former SWE depth (m)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%swe_depth_old, &
                             varid=swe_depth_old_ID, dim=1, wformat=wformat)

   ! (10) Former canopy interception storage (m)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%canopy_int_old, &
                             varid=canopy_int_old_ID, dim=1, wformat=wformat)

   ! (11) Topographic grid (with snowdepth included) (m)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%topo, &
                             varid=topo_ID, dim=1, wformat=wformat)

   ! (12) Sum of snowfall (m)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             snowmodel_struc(n)%sm%sum_sprec, &
                             varid=sum_sprec_ID, dim=1, wformat=wformat)


end subroutine SnowModel_dump_restart



