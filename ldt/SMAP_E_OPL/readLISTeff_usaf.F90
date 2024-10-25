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
#include "LDT_NetCDF_inc.h"

!BOP
!
! !ROUTINE: readLIS_Teff_usaf
! \label{readLISTeff}
!
! !REVISION HISTORY:
!  14 FEB 2023: Eric Kemp, Initial Specification
!  1 July 2023: Mahdi Navari, Modified the code to read the unperturbed
!        ensemble. This is a temporary fix to overcome the bias in
!        soil temperature from a bug in SMAP soil moisture DA
!        PR#1385
!
! !INTERFACE:
subroutine readLIS_Teff_usaf(n, yyyymmdd, hh, Orbit, teff, rc)
! !USES:
  use LDT_coreMod
  use LDT_domainMod
  use LDT_logMod
  use LDT_smap_e_oplMod

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(8), intent(in) :: yyyymmdd
  character(2), intent(in) :: hh
  character(1), intent(in) :: Orbit
  real, intent(out) :: teff(LDT_rc%lnc(n),LDT_rc%lnr(n))
  integer, intent(out) :: rc

!EOP
  integer :: c,r
  real    :: tsoil(LDT_rc%lnc(n),LDT_rc%lnr(n),4)
  character(255) :: fname
  logical :: file_exists
  integer :: rc1
  integer :: gid
  integer :: nens
  integer, allocatable :: str_tind(:)
  integer, allocatable :: ntiles_pergrid(:)
  real, parameter :: kk = 1.007
  real, parameter :: cc_6am = 0.246
  real, parameter :: cc_6pm = 1.000

  external :: create_LISsoilT_filename_usaf
  external :: read_LIStsoil_data_usaf

  ! Initializations
  teff = LDT_rc%udef
  tsoil = LDT_rc%udef
  rc = 1 ! Assume error by default, update below

  ! Set up basic info on Air Force product
  nens = SMAPeOPL%num_ens
  allocate(str_tind(LDT_rc%gnc(n) * LDT_rc%gnr(n)))
  if (SMAPeOPL%ntiles_pergrid .ne. 1) then
     write(LDT_logunit,*) &
          '[ERR] Current SMAP_E_OPL code assumes ntiles_pergrid = 1'
     write(LDT_logunit,*) &
          '[ERR] Actual value is ', SMAPeOPL%ntiles_pergrid
     call LDT_endrun()
  end if
  do gid = 1, (LDT_rc%gnc(n) * LDT_rc%gnr(n))
     str_tind(gid) = ((gid - 1) * nens) + 1
  end do
  allocate(ntiles_pergrid(LDT_rc%gnc(n) * LDT_rc%gnr(n)))
  ntiles_pergrid = SMAPeOPL%ntiles_pergrid

  call create_LISsoilT_filename_usaf(SMAPeOPL%LISdir, &
       yyyymmdd, hh, fname)

  inquire(file=trim(fname), exist=file_exists)
  if (file_exists) then
     write(LDT_logunit,*) '[INFO] Reading ', trim(fname)
     call read_LIStsoil_data_usaf(n, SMAPeOPL%num_tiles, str_tind, &
          ntiles_pergrid, nens, &
          fname, tsoil, rc1)
     if (rc1 .ne. 0) then
        write(LDT_logunit,*) '[ERR] Cannot read from ', trim(fname)
        return
     endif
     write(LDT_logunit,*) '[INFO] Finished reading ', trim(fname)

     ! Calculate effective temperature
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if (Orbit == "D") then
              if (tsoil(c,r,1) > 273.15 .and. tsoil(c,r,2) > 273.15) then
                 teff(c,r) = kk * &
                      (cc_6am * tsoil(c,r,1) + (1 - cc_6am) * tsoil(c,r,2))
              endif
           elseif (Orbit == "A") then
              if (tsoil(c,r,1) > 273.15 .and. tsoil(c,r,2) > 273.15) then
                 teff(c,r) = kk * &
                      (cc_6pm * tsoil(c,r,1) + (1 - cc_6pm) * tsoil(c,r,2))
              endif
           endif
        end do
     end do
     rc = 0
  end if

  if (allocated(str_tind)) deallocate(str_tind)
  if (allocated(ntiles_pergrid)) deallocate(ntiles_pergrid)

end subroutine readLIS_Teff_usaf

!BOP
!
! !ROUTINE: read_LIStsoil_data_usaf
! \label{read_LIStsoil_data_usaf}
!
! !INTERFACE:
subroutine read_LIStsoil_data_usaf(n, ntiles, str_tind, ntiles_pergrid, nens, &
     fname, tsoil, rc)
!
! !USES:
  use LDT_logMod
  use LDT_coreMod
  use LDT_smap_e_oplMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

!
! !INPUT PARAMETERS:
!
  integer, intent(in) :: n
  integer, intent(in) :: ntiles
  integer, intent(in) :: str_tind(LDT_rc%gnc(n) * LDT_rc%gnc(n))
  integer, intent(in) :: ntiles_pergrid(LDT_rc%gnc(n) * LDT_rc%gnc(n))
  integer, intent(in) :: nens
  character(*), intent(in) :: fname
  real, intent(inout) :: tsoil(LDT_rc%lnc(n),LDT_rc%lnr(n),4)
  integer, intent(out) :: rc
!EOP

  integer :: ios, ncid
  integer :: ntiles_dimid, ntiles_local
  integer :: SoilTemp_profiles_dimid, SoilTemp_profiles
  integer :: SoilTemp_inst_id
  real, allocatable :: SoilTemp_inst_tiles(:,:)
  real :: SoilTemp_inst_ensmean_1layer(LDT_rc%gnc(n), LDT_rc%gnr(n))
  integer :: k, c, r

  external :: calc_gridded_ensmean_1layer

  rc = 1 ! Initialize as error, reset near bottom

  ! Sanity checks
  if (LDT_rc%gnc(n) .ne. LDT_rc%lnc(n)) then
     write(LDT_logunit,*)'[ERR] Mismatched dimensions!'
     write(LDT_logunit,*)'[ERR] LDT_rc%gnc(n) = ', LDT_rc%gnc(n)
     write(LDT_logunit,*)'[ERR] LDT_rc%lnc(n) = ', LDT_rc%lnc(n)
     call LDT_endrun()
  end if
  if (LDT_rc%gnr(n) .ne. LDT_rc%lnr(n)) then
     write(LDT_logunit,*)'[ERR] Mismatched dimensions!'
     write(LDT_logunit,*)'[ERR] LDT_rc%gnr(n) = ', LDT_rc%gnr(n)
     write(LDT_logunit,*)'[ERR] LDT_rc%lnr(n) = ', LDT_rc%lnr(n)
     call LDT_endrun()
  end if

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname), mode=NF90_NOWRITE ,ncid=ncid)
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Error opening file' // trim(fname)
     return
  end if

  ios = nf90_inq_dimid(ncid, "ntiles", ntiles_dimid)
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Cannot find ntiles in ' // trim(fname)
     ios = nf90_close(ncid=ncid)
     return
  end if

  ios = nf90_inquire_dimension(ncid, ntiles_dimid, len=ntiles_local)
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Cannot get dimension ntiles in ' &
          // trim(fname)
     ios = nf90_close(ncid=ncid)
     return
  end if

  ! if (ntiles .ne. LDT_rc%glbntiles_red(n)) then
  !    write(LDT_logunit,*)'[ERR] Dimension mismatch!'
  !    write(LDT_logunit,*)'[ERR] ntiles = ', ntiles
  !    write(LDT_logunit,*)'[ERR] LDT_rc%glbntiles_red(n) = ', &
  !         LDT_rc%glbntiles_red(n)
  !    call LDT_endrun()
  ! end if
  if (ntiles_local .ne. ntiles) then
     write(LDT_logunit,*)'[ERR] Dimension mismatch!'
     write(LDT_logunit,*)'[ERR] ntiles = ', ntiles_local
     write(LDT_logunit,*)'[ERR] Expected ', ntiles
     call LDT_endrun()
  end if

  ios = nf90_inq_dimid(ncid, "SoilTemp_profiles", SoilTemp_profiles_dimid)
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Cannot find SoilTemp_profiles in ' &
          // trim(fname)
     ios = nf90_close(ncid=ncid)
     return
  end if

  ios = nf90_inquire_dimension(ncid, SoilTemp_profiles_dimid, &
       len=SoilTemp_profiles)
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Cannot get dimension SoilTemp_profiles in ' &
          // trim(fname)
     ios = nf90_close(ncid=ncid)
     return
  end if
  if (SoilTemp_profiles .ne. 4) then
     write(LDT_logunit,*)'[ERR] Dimension mismatch!'
     write(LDT_logunit,*)'[ERR] SoilTemp_profiles should be 4, but is ', &
          SoilTemp_profiles
     call LDT_endrun()
  end if

  ios = nf90_inq_varid(ncid, 'SoilTemp_inst', SoilTemp_inst_id)
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Cannot find SoilTemp_inst in ' // trim(fname)
     ios = nf90_close(ncid=ncid)
     return
  end if

  allocate(SoilTemp_inst_tiles(ntiles, SoilTemp_profiles))
  SoilTemp_inst_tiles = 0

  ios = nf90_get_var(ncid, SoilTemp_inst_id, SoilTemp_inst_tiles, &
       start=(/1, 1/), &
       count=(/ntiles, SoilTemp_profiles/))
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Cannot read SoilTemp_inst in ' // trim(fname)
     deallocate(SoilTemp_inst_tiles)
     ios = nf90_close(ncid=ncid)
     return
  end if

  ! Calculate ensemble mean in 2d grid space, for each soil layer
  do k = 1, SoilTemp_profiles
     !call calc_gridded_ensmean_1layer(n, ntiles, str_tind, ntiles_pergrid, &
     !     nens, &
     !     SoilTemp_inst_tiles(:,k), &
     !     SoilTemp_inst_ensmean_1layer)

     call calc_gridded_lastens_1layer(n, ntiles, str_tind, ntiles_pergrid, &
          nens, &
          SoilTemp_inst_tiles(:,k), &
          SoilTemp_inst_ensmean_1layer) ! Note to minimize the code chnage: SoilTemp_inst_ensmean_1layer is actually
                                        ! SoilTemp_inst_lastens_1layer (ens #12)

  do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if (SoilTemp_inst_ensmean_1layer(c,r) > 0) then
              tsoil(c,r,k) = SoilTemp_inst_ensmean_1layer(c,r)
           end if
        end do
     end do
  end do

  ! Clean up
  deallocate(SoilTemp_inst_tiles)

  ios = nf90_close(ncid=ncid)
  if (ios .ne. 0) then
     write(LDT_logunit,*) '[WARN] Error closing ' // trim(fname)
     return
  end if

  rc = 0 ! No error detected

#endif

end subroutine read_LIStsoil_data_usaf

! Subroutine for calculating 2d gridded ensemble mean for a single soil layer,
! from tiled data.
subroutine calc_gridded_ensmean_1layer(n, ntiles, str_tind, ntiles_pergrid, &
     nens, gvar_tile, gvar)

  ! Imports
  use LDT_coreMod, only: LDT_rc, LDT_domain, LDT_masterproc
  use LDT_logMod, only: LDT_logunit

  ! Defaults
  implicit none

  ! Arguments
  integer, intent(in) :: n
  integer, intent(in) :: ntiles
  integer, intent(in) :: str_tind(LDT_rc%gnc(n) * LDT_rc%gnr(n))
  integer, intent(in) :: ntiles_pergrid(LDT_rc%gnc(n) * LDT_rc%gnr(n))
  integer, intent(in) :: nens
  real, intent(in) :: gvar_tile(ntiles)
  real, intent(out) :: gvar(LDT_rc%gnc(n), LDT_rc%gnr(n))

  ! Locals
  integer :: m, r, c, gid, stid, tid

  gvar = 0
  if (LDT_masterproc) then
     do r = 1, LDT_rc%gnr(n)
        do c = 1, LDT_rc%gnc(n)
           gid = c + ((r-1) * LDT_rc%gnc(n))
           stid = str_tind(gid)
           if (ntiles_pergrid(gid) > 0) then
              do m = 1, nens
                 tid = stid + m - 1
                 gvar(c,r) = gvar(c,r) + gvar_tile(tid)
              enddo
              gvar(c,r) = gvar(c,r) / nens
           end if
        end do
     end do
  end if
end subroutine calc_gridded_ensmean_1layer

! Subroutine for extracting last ensemble member for a single soil layer,
! from tiled data.
 subroutine calc_gridded_lastens_1layer(n, ntiles, str_tind, ntiles_pergrid, &
            nens, gvar_tile, gvar)

   ! Imports
   use LDT_coreMod, only: LDT_rc, LDT_domain, LDT_masterproc

   ! Defaults
   implicit none

   ! Arguments
   integer, intent(in) :: n
   integer, intent(in) :: ntiles
   !real, intent(in) :: gvar_tile(LDT_rc%glbntiles_red(n))
   integer, intent(in) :: str_tind(LDT_rc%gnc(n) * LDT_rc%gnr(n))
   integer, intent(in) :: ntiles_pergrid(LDT_rc%gnc(n) * LDT_rc%gnr(n))
   real, intent(out) :: gvar(LDT_rc%gnc(n), LDT_rc%gnr(n))
   real, intent(in)  :: gvar_tile(ntiles)
   integer, intent(in) :: nens

   ! Locals
   integer :: m, r, c, gid, stid, tid

   if (LDT_masterproc) then
      gvar = 0
      do r = 1, LDT_rc%gnr(n)
         do c = 1, LDT_rc%gnc(n)
            gid = c + ((r-1) * LDT_rc%gnc(n))
            stid = str_tind(gid)
            if (ntiles_pergrid(gid) > 0) then
                    m = nens
               tid = stid + m - 1
               gvar(c,r) = gvar_tile(tid)
            end if
         end do
      end do
   end if
 end subroutine calc_gridded_lastens_1layer

!BOP
! !ROUTINE: create_LISsoilT_filename_usaf
! \label{create_LISsoilT_filename_usaf}
!
! !INTERFACE:
subroutine create_LISsoilT_filename_usaf(LISdir, yyyymmdd, hh, filename)
 !USES:

  implicit none
 !ARGUMENTS:
  character(*), intent(in) :: LISdir
  character(8), intent(in) :: yyyymmdd
  character(2), intent(in) :: hh
  character(*), intent(out) :: filename
!EOP

  filename = trim(LISdir) &
       // '/PS.AFWA' &
       // '_SC.U' &
       // '_DI.C' &
       // '_DC.ANLYS' &
       // '_GP.LIS' &
       // '_GR.C0P09DEG' &
       // '_AR.GLOBAL' &
       // '_PA.03-HR-SUM' &
       // '_DD.' // yyyymmdd &
       // '_DT.' // hh // '00' &
       // '_DF.nc'

end subroutine create_LISsoilT_filename_usaf
