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
!BOP
! !ROUTINE: read_era5_elev
!  \label{read_era5_elev}
!
! !REVISION HISTORY:
!
!  23 dec 2019: Sujay Kumar; Added ERA5 terrain height reader
!
! !INTERFACE:
subroutine read_era5_elev( n, findex, era5elev, elevdiff )

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_metforcingMod, only : LDT_forc
  use era5_forcingMod, only : era5_struc
  use LDT_logMod,        only : LDT_logunit, LDT_verify, &
                                LDT_endrun
  use LDT_fileIOMod,     only : LDT_transform_paramgrid
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!- Terrain height will be set to run domain:
  real, intent(inout) :: era5elev(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
  real, intent(inout) :: elevdiff(LDT_rc%met_nc(findex), LDT_rc%met_nr(findex))

! !DESCRIPTION:
!
!  Opens, reads, and interpolates ERA5 model elevation to the LDT
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the forcing dataset selected
!  \end{description}
! 
!  The routines invoked are: 
!   \begin{description}
!   \end{description}
!EOP
   logical :: file_exists
   integer :: ftn_const
   integer :: i,c,r,k,iret

   integer :: phisId
   real    :: read_phis(era5_struc(n)%ncold, era5_struc(n)%nrold,1)

   ! Grid transform fields:
   integer   :: inpts, outpts
   real      :: elev1d(era5_struc(n)%ncold*era5_struc(n)%nrold)
   logical*1 :: lb(era5_struc(n)%ncold*era5_struc(n)%nrold)
   real      :: elev_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   logical*1 :: lb_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))

! _____________________________________________________________________________

#if (defined USE_NETCDF3) 
  write(LDT_logunit,*) "ERR: ERA5 terrain height reader requires NetCDF4"
  call LDT_endrun()
#endif

   era5elev = LDT_rc%udef
   elevdiff = LDT_rc%udef

   ! Check if ERA5 grid is selected but downscaling options, 
   !  like bilinear, is incorrectly selected.
   if( era5_struc(n)%gridDesc(9)  == LDT_rc%gridDesc(n,9) .and. &
       era5_struc(n)%gridDesc(10) == LDT_rc%gridDesc(n,10).and. &
       LDT_rc%gridDesc(n,1) == 0 .and. &
       LDT_rc%met_gridtransform_parms(findex) .ne. "neighbor" ) then
      write(LDT_logunit,*) "[WARN] The ERA5 grid was selected for the"
      write(LDT_logunit,*) "  LDT run domain; however, 'bilinear', 'budget-bilinear',"
      write(LDT_logunit,*) "  or some other unknown option was selected to spatially"
      write(LDT_logunit,*) "  downscale the grid, which will cause errors during runtime."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun()
   endif

   inquire(file = trim(era5_struc(n)%era5hgt_file), exist=file_exists)
   if(.not. file_exists) then
      write(LDT_logunit,*) "[ERR] The ERA5 terrain height file ",&
            trim(era5_struc(n)%era5hgt_file)," is not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

! -------------------------------------------------------------------
! Open and Read-in Forcing Terrain Hght File - Bring to LIS run domain
! -------------------------------------------------------------------
   write(LDT_logunit,*) "[INFO] Reading the ERA5 terrain height file: ", &
        trim(era5_struc(n)%era5hgt_file)

#if (defined USE_NETCDF4) 

   ! Open the MERRA-2 "constants" file to read in the geopotential
   !  height field (PHIS):
   call LDT_verify(nf90_open(path=trim(era5_struc(n)%era5hgt_file), &
            mode=NF90_NOWRITE, ncid=ftn_const), &
           'nf90_open failed in read_era5_elev')

   call LDT_verify(nf90_inq_varid(ftn_const,'elev',phisId), &
           'nf90_inq_varid failed for elev in read_era5_elev')

   ! Reading in PHIS field:
   call LDT_verify(nf90_get_var(ftn_const, phisId, read_phis), &
           'nf90_get_var failed for elev in read_era5_elev')


   ! Initialize arrays for grid transformation:
   inpts  = LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex)
   outpts = LDT_rc%lnr(n)*LDT_rc%lnc(n)
   lb     = .true.
   lb_regrid = .true.
   elev1d = -9999.0

   ! Convert 2D to 1D array for the interplation call:
   do r = 1, LDT_rc%met_nr(findex)
      do c = 1, LDT_rc%met_nc(findex)
         k= c+(r-1)*LDT_rc%met_nc(findex)
         elev1d(k) = read_phis(c,r,1)
         if ( elev1d(k) == -9999.0 ) then
           elev1d(k)  = LDT_rc%udef
           lb(k) = .false.
         endif
      enddo
   enddo

   ! Interp elevation field to output field:
   call LDT_transform_paramgrid(n, LDT_rc%met_gridtransform_parms(findex), &
            LDT_rc%met_gridDesc(findex,:), inpts, 1, elev1d, lb, &
            outpts, elev_regrid, lb_regrid )

   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         i = i + 1
         era5elev(c,r,1) = elev_regrid(i) 
      end do
   end do


#endif
end subroutine read_era5_elev
