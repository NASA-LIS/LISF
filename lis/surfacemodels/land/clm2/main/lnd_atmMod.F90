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

module lnd_atmMod

#if (defined COUP_CAM)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Does atm to land and land to atm mapping
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: lnd_atmMod.F90,v 1.6 2004/11/24 22:57:18 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use infnan
  use clm2_varpar             !parameters
  use clm_varmap             !mapping variables 
  use spmdMod
#if (defined SPMD)
  use mpishorthand
#endif
  implicit none

  integer , private, parameter   :: nrecv_lnd = 16
  real(r8), private, allocatable :: recv1d(:,:) 
  
  integer , private, parameter   :: nsend_lnd = 13
  real(r8), private, allocatable :: send1d(:,:)  
  
  SAVE

!===============================================================================
CONTAINS
!===============================================================================

  subroutine atm_to_lnd_mapping (recv2d)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Receive data from the atm
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: lnd_atmMod.F90,v 1.6 2004/11/24 22:57:18 jim Exp $
!-----------------------------------------------------------------------

    use clm2_lsmMod
    use clm2_varcon             !physical constants
    use clm2_varsur             !surface variables

! --------------------------- arguments------ ---------------------
    real(r8), intent(in) :: recv2d(lsmlon,nrecv_lnd,lsmlat) !input from atm
! -----------------------------------------------------------------

! --------------------------- Local variables ---------------------
    integer  :: i,j,k,n             !indices 
    real(r8) :: forc_rainc          !rainxy Atm flux mm s-1   
    real(r8) :: forc_rainl          !rainxy Atm flux mm s-1   
    real(r8) :: forc_snowc          !snowfxy Atm flux  mm s-1 
    real(r8) :: forc_snowl          !snowfxl Atm flux  mm s-1 
#if (defined SPMD)
    integer  :: ier                 !error code 
    integer  :: numsendv(0:npes-1)  !vector of items to be sent
    integer  :: displsv(0:npes-1)   !displacement vector
    integer  :: numrecv             !number of items to be received
#endif
! -----------------------------------------------------------------

! Map received fields on [lsmlon]x[lsmlat] grid to subgrid vectors 
! of [begpatch:endpatch] points

    if (.not. allocated(recv1d)) then
       allocate (recv1d(nrecv_lnd,numpatch)) 
       recv1d(:,:) = inf
    endif
    if (LIS_masterproc) then
       do k = 1,numpatch
          i = patchvec%ixy(k) 
          j = patchvec%jxy(k) 
          do n = 1,nrecv_lnd
             recv1d(n,k) = recv2d(i,n,j)
          end do
       end do
    end if

#if (defined SPMD)
    call compute_mpigs_patch(nrecv_lnd, numrecv, numsendv, displsv)
    call mpi_scatterv (recv1d            , numsendv, displsv, mpir8, &
                       recv1d(1,begpatch), numrecv , mpir8, 0, mpicom, ier)
#endif

! Split data from atm into component arrays and also determine
! derived quantities. Note that atm precipitation is input in 
! units of m s-1ec and must be converted to units of mm/s.

    do k = begpatch, endpatch
       clm(k)%forc_hgt      = recv1d( 1,k)       !zgcmxy  Atm state m
       clm(k)%forc_u        = recv1d( 2,k)       !forc_uxy  Atm state m s-1
       clm(k)%forc_v        = recv1d( 3,k)       !forc_vxy  Atm state m s-1
       clm(k)%forc_th       = recv1d( 4,k)       !forc_thxy Atm state K
       clm(k)%forc_q        = recv1d( 5,k)       !forc_qxy  Atm state kg kg-1
       clm(k)%forc_pbot     = recv1d( 6,k)       !ptcmxy  Atm state Pa
       clm(k)%forc_t        = recv1d( 7,k)       !forc_txy  Atm state K
       clm(k)%forc_lwrad    = recv1d( 8,k)       !flwdsxy Atm flux  W/m^2
       forc_snowc           = recv1d( 9,k)       !mm s-1
       forc_snowl           = recv1d(10,k)       !mm s-1
       forc_rainc           = recv1d(11,k)       !mm s-1 
       forc_rainl           = recv1d(12,k)       !mm s-1 
#if defined(PERGRO)
!
! For error-growth only allow rain not snowfall
!
       forc_rainc           = forc_rainc + forc_snowc
       forc_rainl           = forc_rainl + forc_snowl
       forc_snowc           = 0.0_r4
       forc_snowl           = 0.0_r4
#endif
       clm(k)%forc_solad(2) = recv1d(13,k)       !forc_sollxy  Atm flux  W/m^2
       clm(k)%forc_solad(1) = recv1d(14,k)       !forc_solsxy  Atm flux  W/m^2 
       clm(k)%forc_solai(2) = recv1d(15,k)       !forc_solldxy Atm flux  W/m^2
       clm(k)%forc_solai(1) = recv1d(16,k)       !forc_solsdxy Atm flux  W/m^2

       ! determine derived quantities

       clm(k)%forc_hgt_u = clm(k)%forc_hgt       !observational height of wind [m] 
       clm(k)%forc_hgt_t = clm(k)%forc_hgt       !observational height of temperature [m]  
       clm(k)%forc_hgt_q = clm(k)%forc_hgt       !observational height of humidity [m]      
       clm(k)%forc_vp    = clm(k)%forc_q*clm(k)%forc_pbot / (0.622+0.378*clm(k)%forc_q)   
       clm(k)%forc_rho   = (clm(k)%forc_pbot-0.378*clm(k)%forc_vp) / (rair*clm(k)%forc_t) 
       clm(k)%forc_co2   = pco2*clm(k)%forc_pbot                                          
       clm(k)%forc_o2    = po2*clm(k)%forc_pbot                                           

       ! Determine precipitation needed by clm

       clm(k)%forc_rain = forc_rainc + forc_rainl
       clm(k)%forc_snow = forc_snowc + forc_snowl

       if ( clm(k)%forc_snow > 0.0_r4  .and. clm(k)%forc_rain > 0.0_r4 ) then
          write(6,*) 'kpatch= ',k,' snow= ',clm(k)%forc_snow,' rain= ',clm(k)%forc_rain, &
               ' CLM cannot currently handle both non-zero rain and snow'
          call endrun
       elseif (clm(k)%forc_rain > 0.) then
          clm(k)%itypprc = 1
       elseif (clm(k)%forc_snow > 0.) then
          clm(k)%itypprc = 2
       else
          clm(k)%itypprc = 0
       endif

    end do

    return
  end subroutine atm_to_lnd_mapping

!===============================================================================

  subroutine lnd_to_atm_mapping_ini (send2d)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Send initial land model data back to the atm model
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: lnd_atmMod.F90,v 1.6 2004/11/24 22:57:18 jim Exp $
!-----------------------------------------------------------------------

    use clm2_lsmMod
    use clm2_varcon, only : sb
    use clm2_varsur, only : landmask

! --------------------------- Arguments------ ---------------------
    real(r8), intent(out) :: send2d(lsmlon,nsend_lnd,lsmlat) !output to atm
! -----------------------------------------------------------------

! --------------------------- Local variables ---------------------
    integer :: i,j,k,n               !loop indices
    integer :: ilen                  !temporary       
    real(r8):: wt                    !remapping weight
#if (defined SPMD)
    integer :: ier                   !error code
    integer :: numrecvv(0:npes-1)    !vector of items to be received  
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numsend               !number of items to be sent
#endif
! -----------------------------------------------------------------

! Determine vector of fields that will be sent to the atm

    if (.not. allocated(send1d)) then
       allocate (send1d(nsend_lnd,numpatch))
       send1d(:,:) = inf
    endif

    do k= begpatch, endpatch
       send1d( 1,k) = clm(k)%t_grnd       !tsxy
       send1d( 2,k) = clm(k)%albd(1)      !asdir
       send1d( 3,k) = clm(k)%albd(2)      !aldir
       send1d( 4,k) = clm(k)%albi(1)      !asdif
       send1d( 5,k) = clm(k)%albi(2)      !aldif
       send1d( 6,k) = clm(k)%h2osno/1000. !snow (convert mm->m)
       send1d( 7,k) = 1.e36
       send1d( 8,k) = 1.e36
       send1d( 9,k) = 1.e36
       send1d(10,k) = 1.e36
       send1d(11,k) = sb*(clm(k)%t_grnd**4)   !lwup
       send1d(12,k) = 1.e36
       send1d(13,k) = 1.e36
    end do
#if (defined SPMD)
    call compute_mpigs_patch(nsend_lnd, numsend, numrecvv, displsv)
    call mpi_gatherv (send1d(1,begpatch), numsend , mpir8, &
                      send1d            , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif

! Map fields from subgrid vector with length [numpatch] to [lsmlon]x[lsmlat] grid.
! NOTE: snow is sent as zero over non-land to be consistent with csm cpl code. 
! NOTE: do not set values over lon-land because that can cause problems with the
! atm values for sea ice temperatures. 

    if (LIS_masterproc ) then
       do n = 1,nsend_lnd
          where(landmask == 1) 
             send2d(:,n,:) = 0. 
          endwhere
       end do
       do k = 1,numpatch
          if (patchvec%wtxy(k) /= 0.) then
             i  = patchvec%ixy(k)    
             j  = patchvec%jxy(k)    
             wt = patchvec%wtxy(k) 
             do n = 1,nsend_lnd
                send2d(i,n,j) = send2d(i,n,j) + send1d(n,k)*wt
             end do
          end if
       end do
    endif
    
    return
  end subroutine lnd_to_atm_mapping_ini

!===============================================================================

  subroutine lnd_to_atm_mapping(send2d)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Send land model data back to the atm
! 
! Method: 
! 
! Author:
! 
!-----------------------------------------------------------------------

    use clm2_lsmMod
    use clm2_varsur, only : landmask

! --------------------------- Arguments------ ---------------------
    real(r8), intent(out) :: send2d(lsmlon,nsend_lnd,lsmlat) !output to atm
! -----------------------------------------------------------------

! --------------------------- Local variables ---------------------
    integer  :: i,j,k,l,m,n         !loop indices
    real(r8) :: wt                  !remapping weight
#if (defined SPMD)
    integer  :: ier                 !error code
    integer  :: numrecvv(0:npes-1)  !vector of items to be received  
    integer  :: displsv(0:npes-1)   !displacement vector
    integer  :: numsend             !number of items to be sent
#endif
! -----------------------------------------------------------------

! Determine vector of fields that will be sent to the atm

    if (.not. allocated(send1d)) then
       allocate (send1d(nsend_lnd,numpatch)) ; send1d(:,:) = inf
    endif
    do k= begpatch, endpatch
       send1d( 1,k) = clm(k)%t_rad                           !tsxy 
       send1d( 2,k) = clm(k)%albd(1)                         !asdir
       send1d( 3,k) = clm(k)%albd(2)                         !aldir
       send1d( 4,k) = clm(k)%albi(1)                         !asdif
       send1d( 5,k) = clm(k)%albi(2)                         !aldif
       send1d( 6,k) = clm(k)%h2osno/1000.                    !snow (convert mm->m)
       send1d( 7,k) = clm(k)%taux                            !taux 
       send1d( 8,k) = clm(k)%tauy                            !tauy
       send1d( 9,k) = clm(k)%eflx_lh_tot                     !lhflx 
       send1d(10,k) = clm(k)%eflx_sh_tot                     !shflx 
       send1d(11,k) = clm(k)%eflx_lwrad_out                  !lwup
       send1d(12,k) = clm(k)%qflx_evap_tot                   !qflx 
       send1d(13,k) = clm(k)%t_ref2m                         !tref
    end do
#if (defined SPMD)
    call compute_mpigs_patch(nsend_lnd, numsend, numrecvv, displsv)
    call mpi_gatherv (send1d(1,begpatch), numsend , mpir8, &
                      send1d            , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif

! Map fields from subgrid vector with length [numpatch] to [lsmlon]x[lsmlat] grid.
! NOTE: use only points with wt > 0 so SPMD code will not use uninitialized 
! stack memory values for arrays like taux. 
! NOTE: do not set values over lon-land because that can cause problems with the
! atm values for sea ice temperatures. 

    if (LIS_masterproc ) then
       do n = 1, nsend_lnd
          where(landmask == 1) 
             send2d(:,n,:) = 0.
          endwhere
       end do
       do k = 1,numpatch
          if (patchvec%wtxy(k) /= 0.) then
             i  = patchvec%ixy(k)    
             j  = patchvec%jxy(k)    
             wt = patchvec%wtxy(k) 
             do n = 1,nsend_lnd
                send2d(i,n,j) = send2d(i,n,j) + send1d(n,k)*wt
             end do
          end if
       end do
    endif

    return
  end subroutine lnd_to_atm_mapping

!===============================================================================

#endif

end module lnd_atmMod





