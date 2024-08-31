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
! !ROUTINE: iniTimeConst
! \label{iniTimeConst} 
! 
! !INTERFACE: 
subroutine iniTimeConst (n)
! !USES:
  use LIS_precisionMod
  use clm2_lsmMod
  use clm2_varpar   , only : nlevsoi, nlevlak, numrad
  use clm2_varcon   , only : istsoil, istice, istdlak, istslak, istwet, spval
  use pft_varcon   , only : ncorn, nwheat, roota_par, rootb_par,  &
                            z0mr, displar, dleaf, rhol, rhos, taul, taus, xl, &
                            qe25, vcmx25, mp, c3psn
  use clm2_varsur   , only : zlak, dzlak, zsoi, dzsoi, zisoi
  use clm2_shr_const_mod, only : SHR_CONST_PI
  use LIS_coreMod,   only : LIS_rc, LIS_surface,LIS_domain
  use LIS_soilsMod, only  : LIS_soils
  use LIS_logMod, only    : LIS_logunit, LIS_endrun

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
!
! !DESCRIPTION: 
!
! This routine initializes the time invariant CLM variables
! 
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest
! \end{description}
!
!EOP

! ------------------------ arguments --------------------------------- 
!  integer  :: maxpatch=17,maxpatch_pft = 13, npatch_urban = 14, & 
!       npatch_lake = 15, npatch_wet=16, npatch_gla = 17

  integer  :: maxpatch,maxpatch_pft, npatch_urban, & 
       npatch_lake, npatch_wet, npatch_gla
  
! --------------------------------------------------------------------
  
! ------------------------ local variables ---------------------------
  integer  :: j,k,l,ib !indices
  real(r8) :: pi              !3.159...
  integer  :: ivt             !vegetation type index
  real(r8) :: bd              !bulk density of dry soil material [kg/m^3]
  real(r8) :: tkm             !mineral conductivity
  real(r8) :: xksat           !maximum hydraulic conductivity of soil [mm s-1]
  real(r8) :: scalez = 0.025  !Soil layer thickness discretization (m)     
  real(r8) :: hkdepth = 0.5   !Length scale for Ksat decrease (m)
  real(r8) :: sand1(1:nlevsoi,1:LIS_rc%npatch(n,LIS_rc%lsm_index)) !temporary sand
  real(r8) :: clay1(1:nlevsoi,1:LIS_rc%npatch(n,LIS_rc%lsm_index)) !temporary clay
! --------------------------------------------------------------------
! --------------------------------------------------------------------
! Initialize local variables for error checking
! --------------------------------------------------------------------
  maxpatch_pft = LIS_rc%nvegtypes
  npatch_urban = LIS_rc%urbanclass
  npatch_lake  = LIS_rc%waterclass
  npatch_wet   = LIS_rc%wetlandclass
  npatch_gla   = LIS_rc%glacierclass
  maxpatch     = max(LIS_rc%urbanclass, LIS_rc%waterclass, &
       LIS_rc%wetlandclass, LIS_rc%glacierclass, LIS_rc%nvegtypes)

  sand1(:,:) = inf
  clay1(:,:) = inf

! --------------------------------------------------------------------
! itypveg, isoicol, itypwat, sand, and clay from 2-d surface type
! LIS modification: Read in soil files and equate lat,lon 
! LIS grid arrays with clm grid
! Reynolds soil parameterization scheme
! Read in soil data maps to gridded arrays
! --------------------------------------------------------------------
  
  write(LIS_logunit,*) 'MSG: iniTimeConst -- Reading soil, clay and color files'
  
  pi = SHR_CONST_PI
   
  write(LIS_logunit,*) 'MSG: iniTimeConst -- Read soil, clay and color files'
  do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%fgrd > 0.0) then
     !valid subgrid patch
        clm2_struc(n)%clm(k)%kpatch = k
        clm2_struc(n)%clm(k)%itypveg = &
             LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt
        clm2_struc(n)%clm(k)%isoicol = int(lis_soils(n)%color(&
             LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col,&
             LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row))

!           clm2_struc(n)%clm(k)%isoicol = color(LIS_domain(n)%tile(k)%col,LIS_domain(n)%tile(k)%row) 
!           clm2_struc(n)%clm(k)%isoicol = 2
!        if(LIS_rc%inc_water_pts) then
           if(clm2_struc(n)%clm(k)%isoicol.eq.LIS_rc%udef) then 
              clm2_struc(n)%clm(k)%isoicol = 2
!           endif
        endif

        if (clm2_struc(n)%clm(k)%itypveg == npatch_urban) then !urban, from pcturb
           clm2_struc(n)%clm(k)%itypwat = istsoil
           clm2_struc(n)%clm(k)%itypveg = 0 !noveg
!           write(LIS_logunit,*) 'ERR: iniTimeConst -- not supposed to be in urban block'
        else if (clm2_struc(n)%clm(k)%itypveg == LIS_rc%bareclass) then
           clm2_struc(n)%clm(k)%itypwat = istsoil
           clm2_struc(n)%clm(k)%itypveg = 0 !noveg
           
        else if (clm2_struc(n)%clm(k)%itypveg == npatch_lake) then  !deep lake, from pctlak
!           write(LIS_logunit,*) 'Lake point ',k
           clm2_struc(n)%clm(k)%itypwat = istdlak
           clm2_struc(n)%clm(k)%itypveg = 0 !noveg
           do l = 1,nlevsoi
              sand1(l,k) = 0._r8
              clay1(l,k) = 0._r8
           end do
        else if (clm2_struc(n)%clm(k)%itypveg == npatch_wet) then   !wetland, from pctwet
!           write(LIS_logunit,*) 'wetland ',k
           clm2_struc(n)%clm(k)%itypwat = istwet
           clm2_struc(n)%clm(k)%itypveg = 0 !noveg
           do l = 1,nlevsoi
              sand1(l,k) = 0._r8 
              clay1(l,k) = 0._r8
           end do
        else if (clm2_struc(n)%clm(k)%itypveg == npatch_gla) then   !glacier, from pctgla
           clm2_struc(n)%clm(k)%itypwat = istice
           clm2_struc(n)%clm(k)%itypveg = 0 !noveg
!           write(LIS_logunit,*) 'glacier point ',k
        else                             !soil 
           clm2_struc(n)%clm(k)%itypwat = istsoil
           do l = 1, nlevsoi
              sand1(l,k)=LIS_surface(n,LIS_rc%lsm_index)%tile(k)%sand
              clay1(l,k)=LIS_surface(n,LIS_rc%lsm_index)%tile(k)%clay

              if ( clay1(l,k) .lt. 0 ) then
                 write(LIS_logunit,*) 'missing clay',k,      &
                      LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col, &
                      LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row, &
                      clay1(l,k)
              endif
              if ( sand1(l,k) .lt. 0 ) then
                 write(LIS_logunit,*) 'missing sand..',k,    &
                      LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col, &
                      LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row, &
                      sand1(l,k)
              endif
           enddo
        end if
     end if
  enddo
! --------------------------------------------------------------------
! tag lake points
! --------------------------------------------------------------------

  do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (clm2_struc(n)%clm(k)%itypwat==istdlak .or. clm2_struc(n)%clm(k)%itypwat==istslak) then
        clm2_struc(n)%clm(k)%lakpoi = .true.
     else
        clm2_struc(n)%clm(k)%lakpoi = .false.
     end if
  end do

! --------------------------------------------------------------------
! latitudes and longitudes
! --------------------------------------------------------------------
  if(LIS_rc%gridDesc(n,9) .ne. 0.01) then 
     pi = SHR_CONST_PI
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        clm2_struc(n)%clm(k)%lat = LIS_domain(n)%grid(&
             LIS_surface(n,LIS_rc%lsm_index)%tile(k)%index)%lat *pi /180
        clm2_struc(n)%clm(k)%lon = LIS_domain(n)%grid(&
             LIS_surface(n,LIS_rc%lsm_index)%tile(k)%index)%lon *pi /180
     enddo
  endif
! --------------------------------------------------------------------
! Define layer structure for soil and lakes 
! Vertical profile of snow is initialized in routine iniTimeVar
! --------------------------------------------------------------------

  if (nlevlak /= nlevsoi) then
     write(LIS_logunit,*)'number of soil levels and number of lake levels must be the same'
     write(LIS_logunit,*)'nlevsoi= ',nlevsoi,' nlevlak= ',nlevlak
     call LIS_endrun
  endif
  
  dzlak(1) = 1.               
  dzlak(2) = 2.               
  dzlak(3) = 3.               
  dzlak(4) = 4.               
  dzlak(5) = 5.
  dzlak(6) = 7.
  dzlak(7) = 7.
  dzlak(8) = 7.
  dzlak(9) = 7.
  dzlak(10)= 7.

  zlak(1) =  0.5
  zlak(2) =  1.5
  zlak(3) =  4.5
  zlak(4) =  8.0
  zlak(5) = 12.5
  zlak(6) = 18.5
  zlak(7) = 25.5
  zlak(8) = 32.5
  zlak(9) = 39.5
  zlak(10)= 46.5

  do j = 1, nlevsoi
     zsoi(j) = scalez*(exp(0.5*(j-0.5))-1.)    !node depths
  enddo
  
  dzsoi(1) = 0.5*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
  do j = 2,nlevsoi-1
     dzsoi(j)= 0.5*(zsoi(j+1)-zsoi(j-1)) 
  enddo
  dzsoi(nlevsoi) = zsoi(nlevsoi)-zsoi(nlevsoi-1)
  
  zisoi(0) = 0.
  do j = 1, nlevsoi-1
     zisoi(j) = 0.5*(zsoi(j)+zsoi(j+1))         !interface depths
  enddo
  zisoi(nlevsoi) = zsoi(nlevsoi) + 0.5*dzsoi(nlevsoi)

!  write(LIS_logunit,*) 'zsoi ',zsoi
!  write(LIS_logunit,*) 'dzsoi ',dzsoi
!  write(LIS_logunit,*) 'zisoi ',zisoi
!  stop

  do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (clm2_struc(n)%clm(k)%itypwat == istdlak) then        !assume all lakes are deep lakes
!        clm2_struc(n)%clm(k)%z(1:nlevlak)  = zlak(1:nlevlak)
!        clm2_struc(n)%clm(k)%dz(1:nlevlak) = dzlak(1:nlevlak)
        clm2_struc(n)%clm(k)%z(1:nlevsoi)  = zsoi(1:nlevsoi)
        clm2_struc(n)%clm(k)%dz(1:nlevsoi) = dzsoi(1:nlevsoi)
        clm2_struc(n)%clm(k)%zi(0:nlevsoi) = zisoi(0:nlevsoi)
     else if (clm2_struc(n)%clm(k)%itypwat == istslak) then   !shallow lake (not used)
        clm2_struc(n)%clm(k)%dz(1:nlevlak) = NaN
        clm2_struc(n)%clm(k)%z(1:nlevlak)  = NaN
     else                                       !soil, ice, wetland
        clm2_struc(n)%clm(k)%z(1:nlevsoi)  = zsoi(1:nlevsoi)
        clm2_struc(n)%clm(k)%dz(1:nlevsoi) = dzsoi(1:nlevsoi)
        clm2_struc(n)%clm(k)%zi(0:nlevsoi) = zisoi(0:nlevsoi) 
     endif
  end do
  
! --------------------------------------------------------------------
! Initialize root fraction (computing from surface, d is depth in meter):
! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with beta & d_obs
! given in Zeng et al. (1998).
! --------------------------------------------------------------------
!  do k = begpatch, endpatch
     do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        if (.not. clm2_struc(n)%clm(k)%lakpoi) then
           ivt = clm2_struc(n)%clm(k)%itypveg
           do j = 1, nlevsoi-1
              clm2_struc(n)%clm(k)%rootfr(j) = .5*( exp(-roota_par(ivt)*clm2_struc(n)%clm(k)%zi(j-1))  &
                   + exp(-rootb_par(ivt)*clm2_struc(n)%clm(k)%zi(j-1))  &
                   - exp(-roota_par(ivt)*clm2_struc(n)%clm(k)%zi(j  ))  &
                   - exp(-rootb_par(ivt)*clm2_struc(n)%clm(k)%zi(j  )) )
           end do
           clm2_struc(n)%clm(k)%rootfr(nlevsoi) = .5*( exp(-roota_par(ivt)*&
                clm2_struc(n)%clm(k)%zi(nlevsoi-1))  &
                + exp(-rootb_par(ivt)*clm2_struc(n)%clm(k)%zi(nlevsoi-1)) )
           
        else
           clm2_struc(n)%clm(k)%rootfr(1:nlevsoi) = spval
        end if
     end do

! --------------------------------------------------------------------
! Initialize soil thermal and hydraulic properties 
! --------------------------------------------------------------------
     do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        if (clm2_struc(n)%clm(k)%itypwat == istsoil) then ! .and. & 
!          clm2_struc(n)%clm(k)%itypveg.ne.13) then                !soil
           do j = 1, nlevsoi
              clm2_struc(n)%clm(k)%bsw(j)    = 2.91 + 0.159*clay1(j,k)*100.0
              clm2_struc(n)%clm(k)%watsat(j) = 0.489 - 0.00126*sand1(j,k)*100.0 

              xksat            = 0.0070556 *( 10.**(-0.884+&
                   0.0153*sand1(j,k)*100.0) ) ! mm s-1

              clm2_struc(n)%clm(k)%hksat(j)  = xksat * exp(-clm2_struc(n)%clm(k)%zi(j)/hkdepth)
              clm2_struc(n)%clm(k)%sucsat(j) = 10. * ( 10.**(1.88-0.0131*&
                   sand1(j,k)*100.0) )
              if(sand1(j,k).eq.0.and.clay1(j,k).eq.0) then
                 sand1(j,k) = 0.00001
                 clay1(j,k) = 0.00001
              endif
              tkm              = (8.80*sand1(j,k)*100.0+2.92*&
                   clay1(j,k)*100.0)/(sand1(j,k)*100.0+clay1(j,k)*100.0) ! W/(m K)
              bd               = (1.-clm2_struc(n)%clm(k)%watsat(j))*2.7e3
              clm2_struc(n)%clm(k)%tkmg(j)   = tkm ** (1.- clm2_struc(n)%clm(k)%watsat(j))           
              clm2_struc(n)%clm(k)%tksatu(j) = clm2_struc(n)%clm(k)%tkmg(j)*0.57**clm2_struc(n)%clm(k)%watsat(j)   
              clm2_struc(n)%clm(k)%tkdry(j)  = (0.135*bd + 64.7) / (2.7e3 - 0.947*bd)  
              clm2_struc(n)%clm(k)%csol(j)   = (2.128*sand1(j,k)*100.0+2.385*&
                   clay1(j,k)*100.0)/ (sand1(j,k)*100.0+&
                   clay1(j,k)*100.0)*1.e6  ! J/(m3 K)
           end do
        elseif(clm2_struc(n)%clm(k)%itypveg.eq.13) then
           do j = 1, nlevsoi
              clm2_struc(n)%clm(k)%csol(j) = 1.940E6
              clm2_struc(n)%clm(k)%tkdry(j) = 1.28
              clm2_struc(n)%clm(k)%tkmg(j) = 1.28
              clm2_struc(n)%clm(k)%tksatu(j)  =1.28
              clm2_struc(n)%clm(k)%hksat(j) = 1e-10
              clm2_struc(n)%clm(k)%watsat(j) = 0.001
           enddo
        else                                        
           do j = 1, nlevsoi
              clm2_struc(n)%clm(k)%bsw(j)    =  spval
              clm2_struc(n)%clm(k)%watsat(j) =  spval
              clm2_struc(n)%clm(k)%hksat(j)  =  spval
              clm2_struc(n)%clm(k)%sucsat(j) =  spval
              clm2_struc(n)%clm(k)%tkmg(j)   =  spval
              clm2_struc(n)%clm(k)%tksatu(j) =  spval
              clm2_struc(n)%clm(k)%tkdry(j)  =  spval
              clm2_struc(n)%clm(k)%csol(j)   =  spval
           end do
        end if
     end do
     
! --------------------------------------------------------------------
! Initialize clm derived type components from pft_varcon to avoid
! indirect addressing and be compatible with offline CLM code
! --------------------------------------------------------------------

!  do k = begpatch, endpatch
  do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ivt = clm2_struc(n)%clm(k)%itypveg
     clm2_struc(n)%clm(k)%z0mr    = z0mr(ivt)
     clm2_struc(n)%clm(k)%displar = displar(ivt)
     clm2_struc(n)%clm(k)%dleaf  = dleaf(ivt)
     clm2_struc(n)%clm(k)%xl     = xl(ivt)     
     do ib = 1,numrad 
        clm2_struc(n)%clm(k)%rhol(ib) = rhol(ivt,ib)
        clm2_struc(n)%clm(k)%rhos(ib) = rhos(ivt,ib)
        clm2_struc(n)%clm(k)%taul(ib) = taul(ivt,ib)
        clm2_struc(n)%clm(k)%taus(ib) = taus(ivt,ib)
     end do
     clm2_struc(n)%clm(k)%qe25   = qe25(ivt)      ! quantum efficiency at 25c (umol co2 / umol photon)
     clm2_struc(n)%clm(k)%vcmx25 = vcmx25(ivt)    ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
     clm2_struc(n)%clm(k)%mp     = mp(ivt)        ! slope for conductance-to-photosynthesis relationship
     clm2_struc(n)%clm(k)%c3psn  = c3psn(ivt)     ! photosynthetic pathway: 0. = c4, 1. = c3
  end do

  do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     clm2_struc(n)%clm(k)%dtime = LIS_rc%nts(n)
  end do
  
  write(LIS_logunit,*)
  write(LIS_logunit,30)
  do j = 1,nlevlak
     write(LIS_logunit,40)zlak(j),dzlak(j)
  end do
  write(LIS_logunit,*)
  write(LIS_logunit,35)
  do j = 1,nlevsoi
     write(LIS_logunit,45)zsoi(j),dzsoi(j),zisoi(j)
  end do
  write(LIS_logunit,50)
  write(LIS_logunit,*)
     
30 format(' ',' lake levels ',' lake thickness(m)')
35 format(' ',' soil levels ',' soil thickness(m)',' soil interfaces(m)')
40 format(' ',2(f7.3,8x))
45 format(' ',3(f7.3,8x))
50 format(' ','Note: top level soil interface is set to 0')

  return

end subroutine iniTimeConst 
