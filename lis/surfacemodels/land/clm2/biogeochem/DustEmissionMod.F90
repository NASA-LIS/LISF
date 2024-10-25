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

module DustEmissionMod

 use LIS_precisionMod, only : r8
 use clm2type
 use spmdMod,   only : LIS_masterproc

implicit none

 private
 public Dust
 public Dustini

!=======================================================================
CONTAINS
!=======================================================================

subroutine Dust(clm)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Dust mobilization (adapted from C. Zender's dust model)
! 
! Method: 
! This code simulates dust mobilization from the surface
! into the lowest atmospheric layer
!
! Written to operate at the patch level
!
! On output flx_mss_vrt_dst(ndst) is the surface dust emission 
! (due to wind mobilization) (kg/m**2/s) [ + = to atm]
! 
! Author: Sam Levis
! 
!-----------------------------------------------------------------------
! $Id: DustEmissionMod.F90,v 1.6 2004/11/24 22:56:13 jim Exp $
!-----------------------------------------------------------------------

  use clm2_varpar, only : dst_src_nbr, ndst
  use clm2_varcon, only : grav, istsoil
  implicit none

! ----------------------- arguments -------------------------------
  type (clm1d), intent(inout) :: clm       !CLM 1-D Module
! -----------------------------------------------------------------

! ----------------------- local variables -------------------------
  integer m,n              !indices
  integer mbl_nbr          ![flg] Number of mobilization candidates

  real(r8) dum             !dummy
  real(r8) sqrt2lngsdi     ![frc] Factor in erf argument
  real(r8) lndmaxjovrdmdni ![frc] Factor in erf argument
  real(r8) lndminjovrdmdni ![frc] Factor in erf argument
  real(r8) ryn_nbr_frc_thr_prx_opt ![frc] Threshold friction Reynolds number approximation for optimal size
  real(r8) ryn_nbr_frc_thr_opt_fnc ![frc] Threshold friction Reynolds factor for saltation calculation
  real(r8) icf_fct         !Interpartical cohesive forces factor for saltation calculation
  real(r8) dns_fct         !Density ratio factor for saltation calculation
  real(r8) wnd_frc_rat     ![frc] Wind friction threshold over wind friction
  real(r8) wnd_frc_slt_dlt ![m s-1] Friction velocity increase from saltatn
  real(r8) wnd_rfr_dlt     ![m s-1] Reference windspeed excess over threshld
  real(r8) ovr_src_snk_frc
  real(r8) dst_slt_flx_rat_ttl
  real(r8) flx_mss_hrz_slt_ttl
  real(r8) frc_thr_wet_fct
  real(r8) frc_thr_rgh_fct
  real(r8) wnd_frc_thr_slt
  real(r8) wnd_rfr_thr_slt
  real(r8) wnd_frc_slt

! variables determined only upon firstcall = true

  logical :: firstcall = .true.
  real(r8), save :: tmp1   !Factor in saltation computation
  real(r8), save :: ovr_src_snk_mss(dst_src_nbr,ndst)
  real(r8), save :: lnd_frc_mbl
  real(r8), save :: flx_mss_vrt_dst_ttl

! constants

  real(r8) :: dmt_vma_src(dst_src_nbr) =    &     ![m] Mass median diameter
       (/ 0.0111e-6, 2.524e-6, 42.10e-6 /)
  real(r8) :: gsd_anl_src(dst_src_nbr) =    &     ![frc] Geometric std deviation
       (/ 1.89, 2.0, 2.13 /)
  real(r8) :: mss_frc_src(dst_src_nbr) =    &     ![frc] Mass fraction 
       (/ 0.036, 0.957, 0.007 /)
  real(r8) :: dmt_grd(5) =                  &     ![m] Particle diameter grid
       (/ 0.1e-6, 1.0e-6, 2.5e-6, 5.0e-6, 10.0e-6 /)
  real(r8), parameter :: cst_slt = 2.61           ![frc] Saltation constant
  real(r8), parameter :: dmt_slt_opt = 75.0e-6    ![m] Optim diam for saltation
  real(r8), parameter :: dns_slt = 2650.0         ![kg m-3] Density of optimal saltation particles
  real(r8), parameter :: flx_mss_fdg_fct = 5.0e-4 ![frc] Empir. mass flx tuning eflx_lh_vegt
  real(r8), parameter :: vai_mbl_thr = 0.1        ![m2 m-2] VAI threshold quenching dust mobilization

! declare this intrinsic function

#if (defined AIX) 
#define ERF erf
#else
#define ERF derf
  real(r8) derf
#endif

! -----------------------------------------------------------------
#if 0 
  if (firstcall) then

! the following comes from (1) szdstlgn.F subroutine ovr_src_snk_frc_get
!                      and (2) dstszdst.F subroutine dst_szdst_ini
! purpose(1): given one set (the "source") of lognormal distributions,
!             and one set of bin boundaries (the "sink"), compute and return
!             the overlap factors between the source and sink distributions
! purpose(2): set important statistics of size distributions

! Sanity check: erf() in SGI /usr/lib64/mips4/libftn.so is bogus

     dum = 1.0
     if (abs(0.8427-ERF(dum))/0.8427>0.001) then
        write (6,'(a,f12.10)') 'erf(1.0) = ',ERF(dum)
        write (6,*) 'dst: ovr_src_snk_frc_get() reports Error function error'
        call endrun
     endif
     dum = 0.0
     if (ERF(dum) /= 0.0) then
        write (6,'(a,f12.10)') 'erf(0.0) = ',ERF(dum)
        write (6,*) 'dst: ovr_src_snk_frc_get() reports Error function error'
        call endrun
     endif

     do m = 1, dst_src_nbr
        sqrt2lngsdi = sqrt(2.0) * log(gsd_anl_src(m))
        do n = 1, ndst
           lndmaxjovrdmdni = log(dmt_grd(n+1)/dmt_vma_src(m))
           lndminjovrdmdni = log(dmt_grd(n  )/dmt_vma_src(m))
           ovr_src_snk_frc = 0.5 * (ERF(lndmaxjovrdmdni/sqrt2lngsdi) - &
                                    ERF(lndminjovrdmdni/sqrt2lngsdi))
           ovr_src_snk_mss(m,n) = ovr_src_snk_frc * mss_frc_src(m)
        enddo
     enddo

! got mbl_bsn_fct in Dustini

! the following code from subr. wnd_frc_thr_slt_get was placed here
! because tmp1 needs to be defined just once

     ryn_nbr_frc_thr_prx_opt = 0.38 + 1331.0 * (100.0*dmt_slt_opt)**1.56

     if (ryn_nbr_frc_thr_prx_opt < 0.03) then
        write (6,*) 'dstmbl: ryn_nbr_frc_thr_prx_opt < 0.03'
        call endrun
     else if (ryn_nbr_frc_thr_prx_opt < 10.0) then
        ryn_nbr_frc_thr_opt_fnc = -1.0 + 1.928 * (ryn_nbr_frc_thr_prx_opt**0.0922)
        ryn_nbr_frc_thr_opt_fnc = 0.1291 * 0.1291 / ryn_nbr_frc_thr_opt_fnc
     else
        ryn_nbr_frc_thr_opt_fnc = 1.0 - 0.0858 * exp(-0.0617*(ryn_nbr_frc_thr_prx_opt-10.0))
        ryn_nbr_frc_thr_opt_fnc = 0.120 * 0.120 * ryn_nbr_frc_thr_opt_fnc * ryn_nbr_frc_thr_opt_fnc
     endif

     icf_fct = 1.0 + 6.0e-07 / (dns_slt * grav * (dmt_slt_opt**2.5))
     dns_fct = dns_slt * grav * dmt_slt_opt
     tmp1 = sqrt(icf_fct * dns_fct * ryn_nbr_frc_thr_opt_fnc)

     firstcall = .false.
  end if

! the following code from subr. lnd_frc_mbl_get was adapted for lsm use
! purpose: return fraction of each gridcell suitable for dust mobilization

! reset this integer every time

  mbl_nbr = 0

! the "bare ground" fraction of the current sub-gridscale cell decreases
! linearly from 1 to 0 as VAI(=tlai+clm%tsai) increases from 0 to vai_mbl_thr

  if (clm%itypwat == istsoil) then
     lnd_frc_mbl = 1.0 - min(clm%tlai+clm%tsai,vai_mbl_thr) / vai_mbl_thr
     lnd_frc_mbl = lnd_frc_mbl * (1.0 - clm%frac_sno)
     if (lnd_frc_mbl > 0.0) then
        mbl_nbr = mbl_nbr + 1
     endif
  else !if ice sheet, wetland, or lake, no dust allowed
     lnd_frc_mbl = 0.0
  end if

  if (lnd_frc_mbl>1.0 .or. lnd_frc_mbl<0.0) then
     write (6,'(a,i3,a,f6.4)')'Error dstmbl: lnd_frc_mbl(',clm%kpatch,') = ',lnd_frc_mbl
     call endrun
  end if

! reset output variable before next if-statement to avoid output = inf

  do n = 1, ndst
     clm%flx_mss_vrt_dst(n) = 0.0
  end do

! return if nowhere on earth lnd_frc_mbl > 0.0

  if (mbl_nbr==0) RETURN

! the following comes from subr. frc_thr_rgh_fct_get
! purpose: compute factor by which surface roughness increases threshold
!          friction velocity (currently a constant)

  frc_thr_rgh_fct = 1.0

! the following comes from subr. frc_thr_wet_fct_get
! purpose: compute factor by which soil moisture increases threshold friction
!          velocity

! adjust threshold velocity for inhibition by moisture

  if (clm%h2osoi_vol(1)>clm%vwc_thr .and. lnd_frc_mbl > 0.0) then
     frc_thr_wet_fct = sqrt( 1.0 + 1.21 * (100.0*(clm%h2osoi_vol(1)-clm%vwc_thr))**0.68 )
  else
     frc_thr_wet_fct = 1.0
  end if

! the following lines come from subr. dst_mbl
! purpose: adjust threshold friction velocity to acct for moisture and
!          roughness. The ratio tmp1 / sqrt(forc_rho) comes from
!          subr. wnd_frc_thr_slt_get which computes dry threshold
!          friction velocity for saltation

  wnd_frc_thr_slt = tmp1 / sqrt(clm%forc_rho) * frc_thr_wet_fct * frc_thr_rgh_fct

! reset these variables which will be updated in the following if-block

!  wnd_frc_slt = clm%fv
! to make sure that they are initialized
  wnd_frc_slt = 0
  flx_mss_hrz_slt_ttl = 0.0
  flx_mss_vrt_dst_ttl = 0.0


  if (lnd_frc_mbl > 0.0) then


! the following line comes from subr. dst_mbl
! purpose: threshold saltation wind speed

!     wnd_rfr_thr_slt = clm%u10 * wnd_frc_thr_slt / clm%fv
     wnd_rfr_thr_slt = 0

! the following if-block comes from subr. wnd_frc_slt_get 
! purpose: compute the saltating friction velocity
! theory: saltation roughens the boundary layer, AKA "Owen's effect"

!     if (clm%u10 >= wnd_rfr_thr_slt) then
!        wnd_rfr_dlt = clm%u10 - wnd_rfr_thr_slt
!        wnd_frc_slt_dlt = 0.003 * wnd_rfr_dlt * wnd_rfr_dlt
!        wnd_frc_slt = clm%fv + wnd_frc_slt_dlt
!     endif

! the following comes from subr. flx_mss_hrz_slt_ttl_Whi79_get
! purpose: compute vertically integrated streamwise mass flux of particles

     if (wnd_frc_slt > wnd_frc_thr_slt) then
        wnd_frc_rat = wnd_frc_thr_slt / wnd_frc_slt
        flx_mss_hrz_slt_ttl = cst_slt * clm%forc_rho * (wnd_frc_slt**3.0) * &
             (1.0 - wnd_frc_rat) * (1.0 + wnd_frc_rat) * (1.0 + wnd_frc_rat) / grav

! the following loop originates from subr. dst_mbl
! purpose: apply land sfc and veg limitations and global tuning factor

        flx_mss_hrz_slt_ttl = flx_mss_hrz_slt_ttl * lnd_frc_mbl * clm%mbl_bsn_fct * &
             flx_mss_fdg_fct
     endif

! the following comes from subr. flx_mss_vrt_dst_ttl_MaB95_get
! purpose: diagnose total vertical mass flux of dust from vertically
!          integrated streamwise mass flux

     dst_slt_flx_rat_ttl = 100.0 * exp( log(10.0) * (13.4*clm%mss_frc_cly_vld - 6.0) )
     flx_mss_vrt_dst_ttl = flx_mss_hrz_slt_ttl * dst_slt_flx_rat_ttl
     
  endif !lnd_frc_mbl > 0.0

! the following comes from subr. flx_mss_vrt_dst_prt in C. Zender's code
! purpose: partition total vertical mass flux of dust into transport bins

  do n = 1, ndst
     if (lnd_frc_mbl > 0.0) then
        do m = 1, dst_src_nbr
           clm%flx_mss_vrt_dst(n) = clm%flx_mss_vrt_dst(n) + &
                ovr_src_snk_mss(m,n) * flx_mss_vrt_dst_ttl
        end do
     end if
  end do
#endif
  return
end subroutine Dust

!=======================================================================

  subroutine Dustini

!-----------------------------------------------------------------------
!
! Purpose: Compute source efficiency factor from topography
! Source: Paul Ginoux
! Modifications by C. Zender and later by S. Levis
!
!-----------------------------------------------------------------------
! $Id: DustEmissionMod.F90,v 1.6 2004/11/24 22:56:13 jim Exp $
!-----------------------------------------------------------------------

    use clm2_lsmMod
    use clm_varmap, only : begpatch, endpatch          
#if (defined MAKE_DUST_BASIN)
    use clm2_varsur
    use clm2_varcon, only : re
    use clm2_areaMod
    include 'netcdf.inc'

! ------------------------ local variables ---------------------------
  character(len=80) felev     !elevation file

! NETCDF related

  integer status              !netCDF error status
  integer ncid                !netCDF dataset id
  integer ngat                !netCDF global attributes
  integer nvar                !netCDF number of variables
  integer tim_id              !netCDF id for time dimension
  integer ndim                !number of dimensions in the data file
  integer varid               !netCDF variable id
  integer dimlen(10)          !dimension lengths
  character(len=8) dimnam(10) !dimension names of fields

  integer longxy_id           !id for this variable
  integer latixy_id           !id for this variable
  integer edgew_id            !id for this variable
  integer edgee_id            !id for this variable
  integer edges_id            !id for this variable
  integer edgen_id            !id for this variable

  integer, parameter :: toplon = 2160
  integer, parameter :: toplat = 1080

  real(r8) lmask(toplon,toplat)
  real(r8) mask_i(toplon,toplat)
  real(r8) elev(toplon,toplat)
  real(r8) latxy(toplon,toplat)
  real(r8) lonxy(toplon,toplat)
  real(r8) bsn_fct_2d(toplon,toplat)
  real(r8) bsn_fct_2d_o(lsmlon,lsmlat)
  real(r8) latst(toplat+1)
  real(r8) lonwt(toplon+1,toplat)
  real(r8) edgew, edgee
  real(r8) edges, edgen       !edges read in from riverinputs
  real(r8) col                !read in and compared against
  real(r8) row                !expected values

  integer n,ii,ji,io,jo       !indices
  integer numlon_t(toplat)
  real(r8) area_t(toplon,toplat)

  integer bsn_pnt_nbr       ! [nbr] Number of gridpoints in basin
  integer itr_idx           ! [idx] Iteration index
  integer lat_bsn_idx       ! [idx] Counting index
  integer lat_bsn_idx_bnd   ! [idx] Bounded counting index
  integer lon_bsn_idx_bnd   ! [idx] Bounded counting index
  integer bsn_lat_nbr_hlf   ! [nbr] Half number of latitude points in basin
  integer bsn_lon_nbr_hlf   ! [nbr] Half number of longitude points in basin
  integer lat_idx           ! [idx] Counting index
  integer lat_ngh_idx       ! [idx] Latitude neighbor index
  integer lat_ngh_idx_bnd   ! [idx] Bounded latitude neighbor index
  integer lon_bsn_idx       ! [idx] Counting index
  integer lon_idx           ! [idx] Counting index
  integer lon_ngh_idx       ! [idx] Longitude neighbor index
  integer lon_ngh_idx_bnd   ! [idx] Bounded longitude neighbor index
  integer ocn_ngh_nbr       ! [nbr] Number of neighbors which are ocean
  real(r8) pi               ! 3.14...
  real(r8) earth_crc        ! [m] Earth's circumference
  real(r8) hgt_dlt_bsn      ! [m] Elevation range spanned by basin
  real(r8) hgt_dlt_bsn_bnd  ! [m] Bounded elevation range spanned by basin
  real(r8) hgt_rng_lcl_nrm  ! [m] Basin peak minus local mean elevation
!                             normalized by basin elevation range
  real(r8) hgt_max          ! [m] Highest elevation in basin
  real(r8) hgt_min          ! [m] Lowest elevation in basin
  real(r8) hgt_sfc_bnd(toplon,toplat) ! [m] Surface height bounded
  real(r8) hgt_sfc_bsn_avg  ! [m] Mean elevation in basin
  real(r8) lat_dlt_dgr      ! [dgr] Cell latitude size
  real(r8) lat_dlt_m        ! [m] Cell latitude size
  real(r8) lon_dlt_dgr      ! [dgr] Cell longitude size
  real(r8) lon_dlt_m        ! [m] Cell longitude size
  real(r8) wgt_hgt_rng      ! [frc] Elevation-range-dependent weight
  real(r8) wgt_ocn_ngh      ! [frc] Coastline-dependent weight

  real(r8):: mask_o(lsmlon,lsmlat)         !output grid: mask (0, 1)
  integer :: mxovr                         !max # of overlapping input cells
  integer :: novr_i2o(lsmlon,lsmlat)       !number of overlapping input cells
  integer , allocatable :: iovr_i2o(:,:,:) !lon index of overlap input cell
  integer , allocatable :: jovr_i2o(:,:,:) !lat index of overlap input cell
  real(r8), allocatable :: wovr_i2o(:,:,:) !weight    of overlap input cell

! Parameters
  integer,parameter::itr_nbr_max=1           ![nbr] Maximum number of iterations
  real(r8),parameter::hgt_dlt_bsn_max=8000.0 ![m] Max elevn diff in a basin
  real(r8),parameter::hgt_dlt_bsn_min=1.0    ![m] Min elevn diff within a basin
  real(r8),parameter::lat_mbl_max_dgr=60.0   ![dgr] Lat of max dust mobilizn

! Ginoux uses basins that are 13 degrees on a side at 1 degree resolution
  real(r8),parameter::bsn_sz_lon=1000.0e3    ![m] Zonal size of basin
  real(r8),parameter::bsn_sz_lat=1000.0e3    ![m] Meridional size of basin

  real(r8) buf1d(begpatch:endpatch)          !temporary array 
#else
  integer :: k                               !index
#endif
! --------------------------------------------------------------------
#if 0 
#if (!defined MAKE_DUST_BASIN)

  do k = begpatch, endpatch
     clm(k)%mbl_bsn_fct = 1.0
  end do

#else

  felev = '/ptmp/slevis/lsminput/new/navytopo/mksrf_elev.nc'
! felev = '/data/slevis/_aux0_/slevis/lsmv2/new/input/raw/mksrf_elev.nc'

! read navy elevation data at 10 x 10 minutes (2160 columns x 1080 rows)
! Open netCDF data file

  call wrap_open (felev, nf_nowrite, ncid)
  write(*,*)'opened navy elevation data in Dustini'

! number of dimensions, variables, global attributes, and
! id of unlimited (time) dimension

  status = nf_inq (ncid, ndim, nvar, ngat, tim_id)
  if (status /= nf_noerr) then
     write (6,*) ' Dustini netCDF error = ',nf_strerror(status)
     call endrun 
  end if

! get lengths of lat,lon dimensions

  col = 0
  row = 0

  do ii = 1, ndim
     call wrap_inq_dim (ncid, ii, dimnam(ii), dimlen(ii))
     if (dimlen(ii) == 2160) col = 2160
     if (dimlen(ii) == 1080) row = 1080
  end do

  if (col /= 2160) then
     write (6,*) 'Dustini error: col = ',col, &
          ' in data file not equal to 2160'
     call endrun
  end if

  if (row /= 1080) then
     write (6,*) 'Dustini error: row = ',row, &
          ' in data file not equal to 1080'
     call endrun
  end if

! Extract edges and lat/lon fields

  status = nf_inq_varid(ncid,'edgew',edgew_id)
  if (status /= nf_noerr) then
     call wrap_inq_varid(ncid,'EDGEW',edgew_id)
  end if
  call wrap_get_var_realx(ncid, edgew_id, edgew)

  status = nf_inq_varid(ncid,'edgee',edgee_id)
  if (status /= nf_noerr) then
     call wrap_inq_varid(ncid,'EDGEE',edgee_id)
  end if
  call wrap_get_var_realx(ncid, edgee_id, edgee)

  status = nf_inq_varid(ncid,'edges',edges_id)
  if (status /= nf_noerr) then
     call wrap_inq_varid(ncid,'EDGES',edges_id)
  end if
  call wrap_get_var_realx(ncid, edges_id, edges)

  status = nf_inq_varid(ncid,'edgen',edgen_id)
  if (status /= nf_noerr) then
     call wrap_inq_varid(ncid,'EDGEN',edgen_id)
  end if
  call wrap_get_var_realx(ncid, edgen_id, edgen)

  call wrap_inq_varid(ncid,'LONGXY',longxy_id)
  call wrap_get_var_realx(ncid, longxy_id, lonxy)

  call wrap_inq_varid(ncid,'LATIXY',latixy_id)
  call wrap_get_var_realx(ncid, latixy_id, latxy)

  call wrap_inq_varid (ncid, 'ELEVATION', varid)
  call wrap_get_var_realx(ncid,varid,elev)

  call wrap_inq_varid (ncid, 'LANDMASK', varid)
  call wrap_get_var_realx(ncid,varid,lmask)

! Close netCDF file

  call wrap_close (ncid)
  if ( LIS_masterproc )then
     write (6,*) '---------------------------------------'
     write (6,*) 'Dustini: closing navy elevation data'
     write (6,*) '---------------------------------------'
     write (6,*)
  end if

! the following code is taken from mkglacier
! -----------------------------------------------------------------
! Map data from input grid to LSM grid. Get:
!    o mxovr    - maxium number of overlapping input cells on LSM grid
!    o novr_i2o - number of input grid cells that overlap each LSM grid cell
!    o iovr_i2o - longitude index of overlapping input grid cell
!    o jovr_i2o - latitude  index of overlapping input grid cell
!    o wovr_i2o - fraction of LSM grid cell overlapped by input grid cell
! -----------------------------------------------------------------

  do ji = 1,toplat
     numlon_t(ji) = 0
     do ii = 1,toplon
        if (lonxy(ii,ji) /= 1.e36) numlon_t(ji) = numlon_t(ji) + 1
     enddo
  enddo

  call celledge(toplat, toplon, numlon_t, lonxy, &
                latxy , edgen , edgee   , edges, &
                edgew , latst , lonwt   )

  call cellarea (toplat, toplon, numlon_t, latst, &
                 lonwt , edgen , edgee   , edges, &
                 edgew , area_t)

  call mkmxovr (toplon, toplat  , numlon_t, lonwt, latst, &
                lsmlon, lsmlat  , numlon  , lonw , lats , &
                mxovr , novr_i2o)

  if ( LIS_masterproc ) write (6,*)'mxovr= ',mxovr
  allocate(iovr_i2o(lsmlon,lsmlat,mxovr))
  allocate(jovr_i2o(lsmlon,lsmlat,mxovr))
  allocate(wovr_i2o(lsmlon,lsmlat,mxovr))

  do ji = 1, toplat
     do ii = 1, numlon_t(ji)
        mask_i(ii,ji) = 1.
     end do
  end do

  do jo = 1, lsmlat
     do io = 1, numlon(jo)
        mask_o(io,jo) = 1.
     end do
  end do

  call areaini ( toplon, toplat, numlon_t, lonwt, latst, area_t, mask_i, &
                 lsmlon, lsmlat, numlon  , lonw , lats , area  , mask_o, &
                 mxovr , novr_i2o, iovr_i2o, jovr_i2o, wovr_i2o )

  do ji = 1, toplat
     do ii = 1, numlon_t(ji)
        mask_i(ii,ji) = lmask(ii,ji)
     end do
  end do

  do jo = 1, lsmlat
     do io = 1, numlon(jo)
        mask_o(io,jo) = 0.
        do n = 1, novr_i2o(io,jo) !overlap cell index
           ii = iovr_i2o(io,jo,n) !lon index (input grid) of overlap cell
           ji = jovr_i2o(io,jo,n) !lat index (input grid) of overlap cell
           mask_o(io,jo) = mask_o(io,jo) + mask_i(ii,ji) * wovr_i2o(io,jo,n)
        end do
     end do
  end do

  call areaini ( toplon, toplat, numlon_t, lonwt, latst, area_t, mask_i, &
                 lsmlon, lsmlat, numlon  , lonw , lats , area  , mask_o, &
                 mxovr , novr_i2o, iovr_i2o, jovr_i2o, wovr_i2o )

! ------------------------------------------------------------------------
! from Charlie's subroutine dst_src_bsn
! ------------------------------------------------------------------------

! Initialize scalars

  pi=SHR_CONST_PI           ! 3.14...
  earth_crc=2.0*pi*re*1000. ! [m] Earth's circumference

  do lat_idx=1,toplat
     do lon_idx=1,toplon

! Initialize output array

        bsn_fct_2d(lon_idx,lat_idx)=0.0

! Assign minimal surface height of 0.0 m

        hgt_sfc_bnd(lon_idx,lat_idx)=max(0._r8, elev(lon_idx,lat_idx))

     end do
  end do

! Iterate to enlarge surrounding domain
! This reinforces weight of large basins compared to small depressions

  do itr_idx=1,itr_nbr_max

! Compute properties of local "basin" for every latitude slice

     do lat_idx=1,toplat

! Number of latitude points in basin depends on latitude if grid is irregular

        lat_dlt_dgr=abs((latst(lat_idx+1)-latst(lat_idx))) ! [dgr]
        if (lat_dlt_dgr > 100.) then
           write (6,*) 'ERROR: lat_dlt_dgr>100. in dst_src_bsn()'
           call endrun
        end if
        lat_dlt_m=lat_dlt_dgr*earth_crc/180.0 ! [m/latitude point]
        bsn_lat_nbr_hlf=0.5*(itr_idx+bsn_sz_lat/lat_dlt_m) ! [nbr] Half number of latitude points in basin
        if (2*bsn_lat_nbr_hlf+1.gt.toplat) then
           write (6,*) 'ERROR: 2*bsn_lat_nbr_hlf+1. gt.toplat in dst_src_bsn()'
           call endrun
        end if

        do lon_idx=1,toplon

! Number of longitude points in basin depends on latitude for irregular grids
! and may even depend on longitude for very irregular (eg reduced) grids

           lon_dlt_dgr=abs((lonwt(lon_idx+1,1)-lonwt(lon_idx,1))) ! [dgr]
           if (lon_dlt_dgr.gt.100.0) then
              write (6,*) 'ERROR: lon_dlt_dgr.gt.100.0 in dst_src_bsn()'
              call endrun
           end if
           lon_dlt_m=lon_dlt_dgr*earth_crc/360.0 ! [m/longitude point]
           bsn_lon_nbr_hlf=0.5*(itr_idx+bsn_sz_lon/lon_dlt_m) ! [nbr] Half number of longitude points in basin
           if (2*bsn_lon_nbr_hlf+1.gt.toplon) then
              write (6,*) 'ERROR: 2*bsn_lon_nbr_hlf+1.gt.toplon in dst_src_bsn()'
              call endrun
           end if

! initialize description of current basin

           bsn_pnt_nbr=0       ! [nbr] Number of gridpoints in basin
           hgt_min=1.0e36      ! [m] Lowest elevation in basin
           hgt_max=-1.0e36     ! [m] Highest elevation in basin
           hgt_sfc_bsn_avg=0.0 ! [m] Mean elevation in basin
           ocn_ngh_nbr=0       ! [nbr] Number of neighbors which are ocean

! Count neighboring land points

           do lat_ngh_idx=lat_idx-1,lat_idx+1
              do lon_ngh_idx=lon_idx-1,lon_idx+1
                 lon_ngh_idx_bnd=min(max(1,lon_ngh_idx),toplon)
                 lat_ngh_idx_bnd=min(max(1,lat_ngh_idx),toplat)
                 if (lmask(lon_ngh_idx_bnd,lat_ngh_idx_bnd) == 0) ocn_ngh_nbr=ocn_ngh_nbr+1
              end do        ! end loop over neighbor lons
           end do           ! end loop over neighbor lats

! Loop over grid points contributing to current basin

           do lat_bsn_idx=lat_idx-bsn_lat_nbr_hlf,lat_idx+bsn_lat_nbr_hlf

! Bound [lon,lat]_bsn_idx indices to set of valid array indices
! Latitude coordinate does not wrap

              lat_bsn_idx_bnd=lat_bsn_idx ! Normal point

! Jump to next iteration of lat_bsn_idx loop

              if (lat_bsn_idx.lt.1.or.lat_bsn_idx.gt.toplat) cycle
              do lon_bsn_idx=lon_idx-bsn_lon_nbr_hlf,lon_idx+bsn_lon_nbr_hlf

! Longitude coordinate wraps

                 lon_bsn_idx_bnd=lon_bsn_idx
                 if (lon_bsn_idx.lt.1) then
                    lon_bsn_idx_bnd=toplon-abs(lon_bsn_idx)
                 else if (lon_bsn_idx.gt.toplon) then
                    lon_bsn_idx_bnd=lon_bsn_idx-toplon
                 endif      ! endif

                 bsn_pnt_nbr=bsn_pnt_nbr+1
                 hgt_sfc_bsn_avg=hgt_sfc_bsn_avg+hgt_sfc_bnd(lon_bsn_idx_bnd,lat_bsn_idx_bnd)

! Elevation extrema in current basin

                 hgt_max=max(hgt_max,hgt_sfc_bnd(lon_bsn_idx_bnd,lat_bsn_idx_bnd))
                 hgt_min=min(hgt_min,hgt_sfc_bnd(lon_bsn_idx_bnd,lat_bsn_idx_bnd))

              enddo         ! end loop over lon_bsn_idx
           enddo            ! end loop over lat_bsn_idx

! Sanity checks

           if (bsn_pnt_nbr > (2*bsn_lat_nbr_hlf+1)*(2*bsn_lon_nbr_hlf+1)) then
              write (6,*) 'ERROR: bsn_pnt_nbr too large in dst_src_bsn'
              call endrun
           end if
           if (bsn_pnt_nbr > 0) then
              hgt_sfc_bsn_avg=hgt_sfc_bsn_avg/bsn_pnt_nbr
           else
              hgt_sfc_bsn_avg=0.0
           end if

           hgt_dlt_bsn=hgt_max-hgt_min
           if (hgt_dlt_bsn < 0.0) then
              write (6,*) 'ERROR: hgt_dlt_bsn<0. in dst_src_bsn()'
              call endrun
           end if

! Approximately 30 basins (at T42) have hgt_dlt_bsn = 0.0 m
! Enforce minimum hgt_dlt_bsn = 1.0 m in these regions to avoid
! divide-by-zero below in wgt_bsn

           hgt_dlt_bsn_bnd=max(hgt_dlt_bsn_min,hgt_dlt_bsn)

! Potential sources are points above MSL with some land

           if (lmask(lon_idx,lat_idx) == 1) then

! Basin characteristics
! Find elevn diff bet peak of basin and local mean elevn and normalize this by
! (bounded) height difference anywhere in basin
! Combined with following step, this ensures deeper portions of basin have more
! mobilization potential

              hgt_rng_lcl_nrm=(hgt_max-hgt_sfc_bnd(lon_idx,lat_idx))/hgt_dlt_bsn_bnd

! 0.0 < wgt_hgt_rng <= 1.0:
! wgt_hgt_rng is 0.0 when local point is highest point in basin
! wgt_hgt_rng is 1.0 when local point is lowest point in basin
! Fifth power dependence comes from Ginoux's tuning

              wgt_hgt_rng=hgt_rng_lcl_nrm**5

! Coastal points are usually lower than surrounding inland points
! However, they are not really basins, just victims of geography
! Taper mobilization efficiency at coastlines by dividing by two for
! every neighboring non-land point
! Gridpoint islands should not be basins either

              wgt_ocn_ngh=2.0**(-ocn_ngh_nbr) ! [frc] Coastline-dep weight

! Sanity check

              if (wgt_ocn_ngh < 1.0 .and. ocn_ngh_nbr == 0) then
                 write (6,*) 'Dustini:: Error here (wgt_ocn_ngh < 1.0 .and. ocn_ngh_nbr == 0)!'
                 call endrun
              end if        ! endif insane

! Increment bsn_fct_2d (rather than straight assign) to take advantage of
! iteration capability

              bsn_fct_2d(lon_idx,lat_idx)=bsn_fct_2d(lon_idx,lat_idx)+wgt_hgt_rng*wgt_ocn_ngh/itr_nbr_max

           endif            ! end if point is potential source

        enddo               ! end loop over lon
     enddo                  ! end loop over lat

  enddo                     ! end loop over itr

! Sanity check
  do lat_idx=1,toplat
     do lon_idx=1,toplon
        if (bsn_fct_2d(lon_idx,lat_idx) < 0. .or. bsn_fct_2d(lon_idx,lat_idx) > 1.) then
           write (6,'(a)') 'dst_src_bsn(): ERROR here too!'
           call endrun
        end if              ! end if error
     enddo                  ! end loop over lon
  enddo                     ! end loop over lat

! next must obtain grid avg of bsn_fct_2d
! -----------------------------------------------------------------
! Process each cell on LSM grid (code from mkglacier)
! -----------------------------------------------------------------

  do jo = 1, lsmlat
     do io = 1, numlon(jo)

! Make area average

        bsn_fct_2d_o(io,jo) = 0.
        do n = 1, novr_i2o(io,jo)  !overlap cell index
           ii = iovr_i2o(io,jo,n)  !lon index (input grid) of overlap cell
           ji = jovr_i2o(io,jo,n)  !lat index (input grid) of overlap cell
           bsn_fct_2d_o(io,jo) = bsn_fct_2d_o(io,jo) + bsn_fct_2d(ii,ji) * wovr_i2o(io,jo,n)
        end do

! Corrections: set oceans to zero

        if (landmask(io,jo) == 0) then
           bsn_fct_2d_o(io,jo) = 0.
        end if

! Error checks

        if (bsn_fct_2d_o(io,jo) > 1.) then
           write (6,*) 'Dustini error: bsn_fct_2d_o > 1. at io, jo = ', io, jo
           call endrun
        else if (bsn_fct_2d_o(io,jo) < 0.) then
           write (6,*) 'Dustini error: bsn_fct_2d_o < 0. at io, jo = ', io, jo
           call endrun
        end if

     end do
  end do

  call xy2v(lsmlon, lsmlat, bsn_fct_2d_o, begpatch, endpatch, buf1d)
  clm(begpatch:endpatch)%mbl_bsn_fct = buf1d(begpatch:endpatch)

  deallocate(iovr_i2o)
  deallocate(jovr_i2o)
  deallocate(wovr_i2o)

#endif
#endif
  return
end subroutine Dustini

!=======================================================================

end module DustEmissionMod
