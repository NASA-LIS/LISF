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
! !ROUTINE: snowmodel_readrst
! \label{snowmodel_readrst}
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel 
!  05 Aug 2022: Kristi Arsenault; Updated SnowModel state restart reader
!
! !INTERFACE:
subroutine snowmodel_readrst()
! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc, LIS_surface
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, &
                             LIS_getNextUnitNumber,   &
                             LIS_releaseUnitNumber,   &
                             LIS_verify
  use snowmodel_lsmMod
  use snowmodel_inc
  use snowmodel_vars
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
!
! !DESCRIPTION:
!  This program reads restart files for SnowModel.  This
!  includes all relevant water/energy storages and tile information. 
!  The following is the list of variables specified in the SnowModel 
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
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!  reads a variable from the restart file
! \item[snowmodel\_coldstart](\ref{snowmodel_coldstart}) \newline
!   initializes the SnowModel state variables
! \end{description}
!EOP

  implicit none

  integer           :: t, l, col, row
  integer           :: nc, nr, npatch
  integer           :: n
  integer           :: ftn
  integer           :: status
  logical           :: file_exists
  character*20      :: wformat
! ___________________________________

   write(LIS_logunit,*) '[INFO] Call to the SnowModel Read restart routine ...'

   do n=1, LIS_rc%nnest
      wformat = trim(snowmodel_struc(n)%rformat)

      ! Coldstart
      if(LIS_rc%startcode .eq. "coldstart") then
         call SnowModel_coldstart(LIS_rc%lsm_index)
      ! Restart
      elseif(LIS_rc%startcode .eq. "restart") then

         ! Check if MicroMet option is set to "SnwoModel";
         !  If so, alert user that this option currently
         !  is not supported with reading in LIS restart files:
         if( snowmodel_struc(n)%sm_micromet_opt == "SnowModel" ) then
            write(LIS_logunit,*) "[ERR] SnowModel restart file read is not currently"
            write(LIS_logunit,*) "  supported with 'SnowModel' forcing-format option."
            call LIS_endrun
         endif

         ! check the existance of restart file
         inquire(file=snowmodel_struc(n)%rfile, exist=file_exists)
         if(.not. file_exists) then
            write(LIS_logunit,*) "[ERR] SnowModel restart file ", &
                                 snowmodel_struc(n)%rfile," does not exist "
            write(LIS_logunit,*) " Program stopping ..."
            call LIS_endrun
         endif
         write(LIS_logunit,*) "[INFO] SnowModel restart file used: ",&
                               trim(snowmodel_struc(n)%rfile)

         ! open restart file
         if(wformat .eq. "binary") then
            ftn = LIS_getNextUnitNumber()
            open(ftn, file=snowmodel_struc(n)%rfile, form="unformatted")
            read(ftn) nc, nr, npatch  ! time, veg class, no. tiles

            ! check for grid space conflict
            if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
               write(LIS_logunit,*) "[ERR] "//trim(snowmodel_struc(n)%rfile), &
                             ":: grid space mismatch - SnowModel run halted"
               call LIS_endrun
            endif

            if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
               write(LIS_logunit,*) "[ERR] SnowModel restart tile space mismatch"
               call LIS_endrun
            endif

         elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
            status = nf90_open(path=snowmodel_struc(n)%rfile, &
                               mode=NF90_NOWRITE, ncid=ftn)
            call LIS_verify(status, "Error opening file "//snowmodel_struc(n)%rfile)
#endif
         endif

         ! read(1): Snow depth (m), density-adjusted and used in subgrid-scale snow
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%snow_d, &
                                  varname="SNOWD", &
                                  wformat=wformat)

         ! read(2): Snow depth (m)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%snow_depth, &
                                  varname="SNOWDEPTH", &
                                  wformat=wformat)

         ! read(3): Canopy interception store (m)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%canopy_int, &
                                  varname="CANOPYINT", &
                                  wformat=wformat)

         ! read(4): Soft snow layer depth (m)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%soft_snow_d, &
                                  varname="SOFTSNOWD", &
                                  wformat=wformat)

         ! read(5): Snow grid density (kg/m3)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%ro_snow_grid, &
                                  varname="ROSNOWGRID", &
                                  wformat=wformat)

         ! read(6): SWE depth (m)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%swe_depth, &
                                  varname="SWEDEPTH", &
                                  wformat=wformat)

         ! read(7): Density of former soft snow layer (kg/m3)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%ro_soft_snow_old, &
                                  varname="ROSOFTSNOWOLD", &
                                  wformat=wformat)

         ! read(8): Initial snow depth, density-layer and subgrid scale snow (m)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%snow_d_init, &
                                  varname="SNOWDINIT", &
                                  wformat=wformat)

         ! read(9): Former SWE depth, step (m)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%swe_depth_old, &
                                  varname="SWEDEPTHOLD", &
                                  wformat=wformat)

         ! read(10): Former canopy interception store (m)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%canopy_int_old, &
                                  varname="CANOPYINTOLD", &
                                  wformat=wformat)

         ! read(11): Snow-depth changing grid topography level (m)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%topo, &
                                  varname="TOPO", &
                                  wformat=wformat)

         ! read(12): Sum of snowfall, over time (m)
         call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  snowmodel_struc(n)%sm%sum_sprec, &
                                  varname="SUMSPREC", &
                                  wformat=wformat)

         ! Close the restart file
         if(wformat .eq. "binary") then
            call LIS_releaseUnitNumber(ftn)
         elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
            status = nf90_close(ftn)
            call LIS_verify(status, "Error in nf90_close in SnowModel_readrst")
#endif
         endif

         ! Assign the read-in states to the SnowModel local states:
         do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
            col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
            row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row

            snow_d(col,row)        = snowmodel_struc(n)%sm(t)%snow_d 
            snow_depth(col,row)    = snowmodel_struc(n)%sm(t)%snow_depth
            canopy_int(col,row)    = snowmodel_struc(n)%sm(t)%canopy_int 
            soft_snow_d(col,row)   = snowmodel_struc(n)%sm(t)%soft_snow_d 
            ro_snow_grid(col,row)  = snowmodel_struc(n)%sm(t)%ro_snow_grid 
            swe_depth(col,row)     = snowmodel_struc(n)%sm(t)%swe_depth
            ro_soft_snow_old(col,row) = snowmodel_struc(n)%sm(t)%ro_soft_snow_old 
            snow_d_init(col,row)   = snowmodel_struc(n)%sm(t)%snow_d_init 
            swe_depth_old(col,row) = snowmodel_struc(n)%sm(t)%swe_depth_old 
            canopy_int_old(col,row)= snowmodel_struc(n)%sm(t)%canopy_int_old 
            topo(col,row)          = snowmodel_struc(n)%sm(t)%topo 
            sum_sprec(col,row)     = snowmodel_struc(n)%sm(t)%sum_sprec 
         enddo

      endif
   enddo

end subroutine snowmodel_readrst

