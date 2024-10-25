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

module atm_lndMod

#if (defined COUP_CAM)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Atm - Land interface module
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: atm_lndMod.F90,v 1.6 2004/11/24 22:57:09 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use pmgrid, only: plon, plond, plat
  use tracers, only: pcnst, pnats
  use rgrid, only: nlon
  use ppgrid, only: pcols, begchunk, endchunk
  use phys_grid
  use comsrf, only : oro, ts, asdir, aldir, asdif, aldif, snowh, lwup, &
                     wsx, wsy, lhf, shf, cflx, tref, &
                     zbot, ubot, vbot, tbot, thbot, qbot, pbot, flwds, &
                     precsc, precsl, precc, precl, soll, sols, solld, solsd 
  use history, only : caseid, ctitle, inithist, nhtfrq, mfilt
  use clm2_shr_const_mod, only: SHR_CONST_PI
  implicit none
  
  private              ! By default make data private
  integer :: landmask(plon,plat) !2d land mask
  integer, allocatable, dimension(:,:) :: landmask_chunk
  
  integer , private, parameter :: nsend_atm = 16
  real(r8), private :: send2d(plon,nsend_atm, plat) !output to clm
  real(r8), allocatable, dimension(:,:,:) :: send2d_chunk
  
  integer , private, parameter :: nrecv_atm = 13
  real(r8), private :: recv2d(plon,nrecv_atm, plat) !input from clm
  real(r8), allocatable, dimension(:,:,:) :: recv2d_chunk

  public atmlnd_ini, atmlnd_drv, is_lsm   ! Public interfaces


!===============================================================================
CONTAINS
!===============================================================================

  subroutine atmlnd_ini()

    use initializeMod, only : initialize           !initialization of clm 
    use lnd_atmMod, only : lnd_to_atm_mapping_ini  !mapping from atm grid space <-> clm tile space
    use error_messages, only: alloc_err
#if ( defined SPMD )
    use mpishorthand
#endif
    use commap    
    use LIS_timeMgrMod, only: get_nstep

#include <comsol.h>
#include <comctl.h>
#include <commss.h>

!-----------------------------------------------------------------------
! Initialize land surface model and obtain relevant atmospheric model 
! arrays back from (i.e. albedos, surface temperature and snow cover over land)
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer  :: i,lat,n,lchnk,ncols !indices
    integer  :: istat               !error return
    integer  :: nstep               !current timestep number
    integer  :: lats(pcols)         !chunk latitudes
    integer  :: lons(pcols)         !chunk longitudes
    real(r8) :: oro_glob(plond,plat)!global oro field
    real(r8) :: landfrac(plon,plat) !2d fractional land
    real(r8) :: latixy(plon,plat)   !2d latitude  grid (degrees)
    real(r8) :: longxy(plon,plat)   !2d longitude grid (degrees)
    real(r8) :: pi                  !3.14159...
!-----------------------------------------------------------------------

! Time management variables.

    nstep = get_nstep()

! Allocate land model chunk data structures

   allocate( landmask_chunk(pcols,begchunk:endchunk), stat=istat )
   call alloc_err( istat, 'atmlnd_ini', 'landmask_chunk', &
                   pcols*(endchunk-begchunk+1) )
   allocate( send2d_chunk(pcols,nsend_atm,begchunk:endchunk), stat=istat )
   call alloc_err( istat, 'atmlnd_ini', 'send2d_chunk', &
                   pcols*nsend_atm*(endchunk-begchunk+1) )
   allocate( recv2d_chunk(pcols,nrecv_atm,begchunk:endchunk), stat=istat )
   call alloc_err( istat, 'atmlnd_ini', 'recv2d_chunk', &
                   pcols*nrecv_atm*(endchunk-begchunk+1) )

! Initialize land model

    call gather_chunk_to_field(1,1,1,plond,oro,oro_glob)
#if (defined SPMD)
    call mpibcast (oro_glob, size(oro_glob), mpir8, 0, mpicom)
#endif

    pi = SHR_CONST_PI
    longxy(:,:) = 1.e36
    do lat = 1,plat
       do i = 1,nlon(lat)
          longxy(i,lat) = (i-1)*360.0/nlon(lat)
          latixy(i,lat) = (180./pi)*clat(lat)
          if (nint(oro_glob(i,lat)) == 1) then
             landmask(i,lat) = 1
             landfrac(i,lat) = 1.0
          else
             landmask(i,lat) = 0
             landfrac(i,lat) = 0.
          endif
       end do
    end do

    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       call get_lat_all_p(lchnk,pcols,lats)
       call get_lon_all_p(lchnk,pcols,lons)
       do i=1,ncols
          landmask_chunk(i,lchnk) = landmask(lons(i),lats(i))
       enddo
    enddo

    call  initialize(eccen    , obliqr   , lambm0  , mvelpp  , caseid  , &
                     ctitle   , nsrest   , nstep   , iradsw  , inithist, &
                     nhtfrq(1), mfilt(1) , longxy  , latixy  , nlon    , &
                     landmask , landfrac , irt)

! For initial run only - get 2d data back from land model (Note that 
! in SPMD case, only LIS_masterproc contains valid recv2d data) and 
! split 2d data into appropriate arrays contained in module comsrf. 

    if (nstep == 0) then

       call lnd_to_atm_mapping_ini(recv2d)
       call scatter_field_to_chunk(1,nrecv_atm,1,plon,recv2d,recv2d_chunk)

       do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          do i=1,ncols
             if (landmask_chunk(i,lchnk) == 1) then
                ts(i,lchnk)    = recv2d_chunk(i, 1,lchnk) 
                asdir(i,lchnk) = recv2d_chunk(i, 2,lchnk) 
                aldir(i,lchnk) = recv2d_chunk(i, 3,lchnk) 
                asdif(i,lchnk) = recv2d_chunk(i, 4,lchnk) 
                aldif(i,lchnk) = recv2d_chunk(i, 5,lchnk) 
                snowh(i,lchnk) = recv2d_chunk(i, 6,lchnk) 
                lwup(i,lchnk)  = recv2d_chunk(i,11,lchnk) 
             endif
          end do
       end do

    endif
    
    return
  end subroutine atmlnd_ini

!===============================================================================

  subroutine atmlnd_drv (nstep, iradsw, eccen, obliqr, lambm0, mvelpp)

!-----------------------------------------------------------------------
! Pack data to be sent to land model into a single array. 
! Send data to land model and call land model driver. 
! Receive data back from land model in a single array.
! Unpack this data into component arrays. 
! NOTE: component arrays are contained in module comsrf.
!-----------------------------------------------------------------------

#if ( defined SPMD )
    use mpishorthand
#endif
    use lnd_atmMod  !mapping from atm grid space <-> clm tile space

!---------------------------Arguments----------------------------------- 
    integer , intent(in) :: nstep    !Current time index
    integer , intent(in) :: iradsw   !Iteration frequency for shortwave radiation
    real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
    real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
    real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox longitude of perihelion + pi (radians)
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer :: i,lat,m,n,lchnk,ncols !indices
    logical doalb          !true if surface albedo calculation time step
!-----------------------------------------------------------------------

! -----------------------------------------------------------------
! Determine doalb
! [doalb] is a logical variable that is true when the next time
! step is a radiation time step. This allows for the fact that
! an atmospheric model may not do the radiative calculations 
! every time step. For example:
!      nstep dorad doalb
!        1     F     F
!        2     F     T
!        3     T     F
!        4     F     F
!        5     F     T
!        6     T     F
! The following expression for doalb is for example only (it is 
! specific to the NCAR CAM). This variable must be calculated
! appropriately for the host atmospheric model
! -----------------------------------------------------------------

    doalb = iradsw==1 .or. (mod(nstep,iradsw)==0 .and. nstep+1/=1)

! Condense the 2d atmospheric data needed by the land surface model into 
! one array. Note that precc and precl precipitation rates are in units 
! of m s-1ec. They are turned into fluxes by multiplying by 1000 kg/m^3.

!$OMP PARALLEL DO PRIVATE(lchnk,ncols,i)
    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       do i=1,ncols
          send2d_chunk(i, 1,lchnk)  =  zbot(i,lchnk)  ! Atmospheric state variable m
          send2d_chunk(i, 2,lchnk)  =  ubot(i,lchnk)  ! Atmospheric state variable m s-1
          send2d_chunk(i, 3,lchnk)  =  vbot(i,lchnk)  ! Atmospheric state variable m s-1
          send2d_chunk(i, 4,lchnk)  =  thbot(i,lchnk) ! Atmospheric state variable K
          send2d_chunk(i, 5,lchnk)  =  qbot(i,lchnk)  ! Atmospheric state variable kg kg-1
          send2d_chunk(i, 6,lchnk)  =  pbot(i,lchnk)  ! Atmospheric state variable Pa
          send2d_chunk(i, 7,lchnk)  =  tbot(i,lchnk)  ! Atmospheric state variable K
          send2d_chunk(i, 8,lchnk)  =  flwds(i,lchnk) ! Atmospheric flux W/m^2
          send2d_chunk(i, 9,lchnk)  =  precsc(i,lchnk)*1000.                  !convert from m s-1ec to mm/sec
          send2d_chunk(i,10,lchnk)  =  precsl(i,lchnk)*1000.                  !convert from m s-1ec to mm/sec
          send2d_chunk(i,11,lchnk)  =  (precc(i,lchnk) - precsc(i,lchnk))*1000. !convert from m s-1ec to mm/sec
          send2d_chunk(i,12,lchnk)  =  (precl(i,lchnk) - precsl(i,lchnk))*1000. !convert from m s-1ec to mm/sec
          send2d_chunk(i,13,lchnk)  =  soll(i,lchnk)  ! Atmospheric flux W/m^2
          send2d_chunk(i,14,lchnk)  =  sols(i,lchnk)  ! Atmospheric flux W/m^2
          send2d_chunk(i,15,lchnk)  =  solld(i,lchnk) ! Atmospheric flux W/m^2
          send2d_chunk(i,16,lchnk)  =  solsd(i,lchnk) ! Atmospheric flux W/m^2
       end do
    end do

    call gather_chunk_to_field(1,nsend_atm,1,plon,send2d_chunk,send2d)

! Convert two dimensional atm input data to one dimensional land model data
 
    call atm_to_lnd_mapping(send2d)

! Call land model driver

    call driver (doalb, eccen, obliqr, lambm0, mvelpp)

! Convert one dimensional land model output data to two dimensional atm data 

    call lnd_to_atm_mapping(recv2d) 
    call scatter_field_to_chunk(1,nrecv_atm,1,plon,recv2d,recv2d_chunk)

! Split 2d recv array into component arrays (in module comsrf)

!$OMP PARALLEL DO PRIVATE(lchnk,ncols,i)
    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       do i=1,ncols
          if (landmask_chunk(i,lchnk) == 1) then
             ts(i,lchnk)     =  recv2d_chunk(i, 1,lchnk) 
             asdir(i,lchnk)  =  recv2d_chunk(i, 2,lchnk) 
             aldir(i,lchnk)  =  recv2d_chunk(i, 3,lchnk) 
             asdif(i,lchnk)  =  recv2d_chunk(i, 4,lchnk) 
             aldif(i,lchnk)  =  recv2d_chunk(i, 5,lchnk) 
             snowh(i,lchnk)  =  recv2d_chunk(i, 6,lchnk) 
             wsx(i,lchnk)    =  recv2d_chunk(i, 7,lchnk) 
             wsy(i,lchnk)    =  recv2d_chunk(i, 8,lchnk) 
             lhf(i,lchnk)    =  recv2d_chunk(i, 9,lchnk) 
             shf(i,lchnk)    =  recv2d_chunk(i,10,lchnk) 
             lwup(i,lchnk)   =  recv2d_chunk(i,11,lchnk) 
             cflx(i,1,lchnk) =  recv2d_chunk(i,12,lchnk) 
             tref(i,lchnk)   =  recv2d_chunk(i,13,lchnk) 
          endif
       end do
    end do
    
! Reset all other consitutent surfaces fluxes to zero over land

    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       do i=1,ncols
          if (landmask_chunk(i,lchnk) == 1) then
             do m = 2,pcnst+pnats
                cflx(i,m,lchnk) = 0.
             end do
          endif
       end do
    end do
    
    return
  end subroutine atmlnd_drv

   logical function is_lsm ( )
     is_lsm = .false.
     return
   end function is_lsm


!===============================================================================

#endif        

end module atm_lndMod

