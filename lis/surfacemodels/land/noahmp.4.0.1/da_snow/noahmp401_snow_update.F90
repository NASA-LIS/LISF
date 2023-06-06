!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_snow_update
! \label{noahmp401_snow_update}
!
! !REVISION HISTORY:
!  13 Aug 2017: Sujay Kumar; Initial specification
!  14 Dec 2018: Yeosang Yoon; Modified code for NoahMP 4.0.1 and SNODEP
!  15 May 2019: Yeosang Yoon; Modified for NoahMP 4.0.1 and LDTSI
!  13 Dec 2019: Eric Kemp; Replaced LDTSI with SNOW
!  05 Jun 2023: Justin Pflug; fixes for SnowModel-defined snow updates

!
! !INTERFACE
subroutine noahmp401_snow_update(n, t, dsneqv, dsnowh)

  use LIS_coreMod
  use NoahMP401_lsmMod
  use module_sf_noahmplsm_401
  use noahmp_tables_401

  implicit none
! 
! !DESCRIPTION: 
!  This subroutine updates relevant snow prognostics based
!  on the update to the total SWE (dsneqv) and total
!  snow depth (dsnowh). The updated variables include
!  number of snow layers, snice, snliq, snow temperature 
!  and snow thickness. 
! 
! !ARGUMENTS:
  integer, intent(in)  :: n
  integer, intent(in)  :: t
  real                 :: dsneqv !mm
  real                 :: dsnowh !m
!EOP

  real, allocatable, dimension(:) :: zsoil
  real, allocatable, dimension(:) :: ficeold
  real, allocatable, dimension(:) :: snice
  real, allocatable, dimension(:) :: snliq
  real, allocatable, dimension(:) :: stc
  real, allocatable, dimension(:) :: supercool
  real, allocatable, dimension(:) :: mice
  real, allocatable, dimension(:) :: mliq
  real, allocatable, dimension(:) :: dzsnso
  real, allocatable, dimension(:) :: zsnso

  integer, allocatable, dimension(:) :: imelt  !phase change index
  real,    allocatable, dimension(:) :: sice

  integer :: snl_idx,i,j,iz
  integer :: iloc, jloc        ! needed, but not use
  real    :: smp,sneqv,snowh
  real    :: snoden
  real    :: sneqv1,snowh1
  real    :: ponding1,ponding2
  integer :: newnode
  integer :: isnow, nsoil, nsnow, soiltype(4), isoil

! local
  real    :: SNOFLOW, BDSNOW
  type (noahmp_parameters)  :: parameters

  isnow = noahmp401_struc(n)%noahmp401(t)%isnow
  nsoil = noahmp401_struc(n)%nsoil
  nsnow = noahmp401_struc(n)%nsnow

  allocate(ficeold(-nsnow+1:0))
  allocate(snice(-nsnow+1:0))
  allocate(snliq(-nsnow+1:0))
  allocate(stc(-nsnow+1:nsoil))
  allocate(imelt(-nsnow+1:nsoil))
  allocate(supercool(-nsnow+1:nsoil))
  allocate(mice(-nsnow+1:nsoil))
  allocate(mliq(-nsnow+1:nsoil))
  allocate(dzsnso(-nsnow+1:nsoil))
  allocate(zsnso(-nsnow+1:nsoil))
  allocate(sice(nsoil))

  imelt = 0

  !set empty snow layers to zero
  do iz = -nsnow+1, isnow
     snice(iz) = 0.
     snliq(iz) = 0.
     stc(iz)   = 0.
     dzsnso(iz)= 0.
     zsnso(iz) = 0.
  enddo

  ! initialize the variables
  soiltype = noahmp401_struc(n)%noahmp401(t)%soiltype
  do isoil = 1, size(soiltype)
    parameters%BEXP(isoil)   = BEXP_TABLE   (SOILTYPE(isoil))
    parameters%PSISAT(isoil) = PSISAT_TABLE (SOILTYPE(isoil))
    parameters%SMCMAX(isoil) = SMCMAX_TABLE (SOILTYPE(isoil))
  end do

  sneqv = noahmp401_struc(n)%noahmp401(t)%sneqv
  snowh = noahmp401_struc(n)%noahmp401(t)%snowh

  zsnso(-nsnow+1:nsoil) = noahmp401_struc(n)%noahmp401(t)%zss(1:nsnow+nsoil) 

! snow/soil layer thickness (m)
  do iz = isnow+1, nsoil
     if(iz == isnow+1) then
        dzsnso(iz) = - zsnso(iz)
     else
        dzsnso(iz) = zsnso(iz-1) - zsnso(iz)
     end if
  end do 

  ! set ZSOIL 
  allocate(zsoil(nsoil))
  ! zsoil is negative.
  zsoil(1) = -noahmp401_struc(n)%sldpth(1)
  do i = 2, nsoil
     zsoil(i) = zsoil(i-1) - NOAHMP401_struc(n)%sldpth(i)
  enddo


  ! state variables 
  snice(-nsnow+1:0) = &
       noahmp401_struc(n)%noahmp401(t)%snowice(1:nsnow)
  snliq(-nsnow+1:0) = &
       noahmp401_struc(n)%noahmp401(t)%snowliq(1:nsnow) 
  stc(-nsnow+1:0) = &
       noahmp401_struc(n)%noahmp401(t)%tsno(1:nsnow)
  ! soil temperature
  stc(1:nsoil) = &
       noahmp401_struc(n)%noahmp401(t)%tslb(1:nsoil)    

  ! NMP snow density calculation
  if(snowh.gt.0) then
     snoden = sneqv/(snowh*1000)
  else
     snoden = 0.55
  endif

  ! allow snow update even in cases where changes opp. directions
  ! alter snow depth change to be in direction of SWE change
  if((dsneqv.gt.0.and.dsnowh.le.0).or.&
          (dsneqv.lt.0.and.dsnowh.ge.0)) then
     dsnowh = (dsneqv/snoden)/1000
  ! set snow depth change to zero in instance where no SWE change
  elseif(dsneqv.eq.0.and.dsnowh.ne.0) then
     dsnowh = 0.
  endif

  ! from snowfall routine
  ! creating a new layer
  if(isnow == 0.and.(dsneqv.gt.0.and.dsnowh.gt.0))  then
     snowh = snowh + dsnowh
     sneqv = sneqv + dsneqv
  end if
  
  newnode = 0 

  if(isnow.eq.0.and.dsneqv.gt.0.and.dsnowh.gt.0) then
     if(snowh.ge.0.025) then
        isnow    = -1
        newnode  =  1
     endif
     dzsnso(0)= snowh
     stc(0)   = min(273.16, noahmp401_struc(n)%noahmp401(t)%sfctmp) 
     snice(0) = sneqv
     snliq(0) = 0.
  end if
  
  ! snow with layers
  if(isnow <  0 .and. newnode == 0 .and. &
       (dsneqv.gt.0.and.dsnowh.gt.0)) then
     snice(isnow+1)  = snice(isnow+1)   + dsneqv
     dzsnso(isnow+1) = dzsnso(isnow+1)  + dsnowh
  endif
  
  if(dsneqv.lt.0.and.dsnowh.lt.0) then
     snowh1 = snowh + dsnowh
     sneqv1 = sneqv + dsneqv
! if dsnowh adjusted since dsneqv and dsnowh in opp. directions
! can cause one or other snowh1 or sneqv1 to be negative
     if(sneqv1.gt.0.and.snowh1.le.0) then
        snowh = ((sneqv1/snoden)/1000)-dsnowh
        snowh1 = snowh + dsnowh
! if SWE disappears, also make sure snow depth disappears
     elseif(sneqv.le.0) then
        sneqv = -dsneqv
        sneqv1 = sneqv + dsneqv
        snowh = -dsnowh
        snowh1 = snowh + dsnowh
     endif
! make sure snow layers currently exist in decrease case
     if(dzsnso(0).eq.0) then
        if(snowh.ge.0.025) then
           isnow = -1
        else
           isnow = 0
        endif
        dzsnso(0)= snowh
        stc(0)   = min(273.16, noahmp401_struc(n)%noahmp401(t)%sfctmp)
        snice(0) = sneqv
        snliq(0) = 0.
     endif
     if(snowh1.ge.0.and.sneqv1.ge.0) then         
        snowh = snowh + dsnowh
        sneqv = sneqv + dsneqv
! snow can no longer fill layer 1
        if(snowh.le.dzsnso(0)) then 
           isnow = 0
           dzsnso(-nsnow+1:(isnow-1)) = 0 
           dzsnso(isnow) = snowh
           snice(-nsnow+1:(isnow-1)) = 0
           snice(isnow) = sneqv
           snliq(-nsnow+1:isnow) = 0
! snow can no longer fill layer 1 and 2
        elseif(snowh.le.(dzsnso(0)+dzsnso(-1))) then 
           isnow = -2
           dzsnso(-nsnow+1:isnow) = 0 
           dzsnso(isnow+1) = snowh -dzsnso(0)
           ! scale swe in layers by ratio of depth to pack
           do snl_idx=-nsnow+1,0
              snice(snl_idx) = sneqv*(dzsnso(snl_idx)/snowh)
           enddo
           snliq(-nsnow+1:isnow) = 0
! all other cases
        elseif(snowh.le.(dzsnso(0)+dzsnso(-1)+dzsnso(-2))) then 
           isnow = -3
           dzsnso(isnow+1) = snowh -dzsnso(-1) -dzsnso(0)
           ! scale swe in layers by ratio of depth to pack
           do snl_idx=-nsnow+1,0
              snice(snl_idx) = sneqv*(dzsnso(snl_idx)/snowh)
           enddo
           snliq(-nsnow+1:isnow) = 0
        endif           
     endif
  endif

  ! ice fraction at the last timestep, add check for both snice and snliq are 0.0
  do snl_idx=isnow+1,0
    if(snice(snl_idx)+snliq(snl_idx)>0.0) then
      ficeold(snl_idx)  = snice(snl_idx) / (snice(snl_idx)+snliq(snl_idx))
    else 
      ficeold(snl_idx)  = 0.0
    endif
  enddo

  sice(:) = max(0.0, noahmp401_struc(n)%noahmp401(t)%smc(:)&
       - noahmp401_struc(n)%noahmp401(t)%sh2o(:))   

  !imelt
  do j = -nsnow+1, nsoil
     supercool(j) = 0.0
  end do
  
  do j = isnow+1,0       ! all layers
     mice(j) = snice(j)
     mliq(j) = snliq(j)
  end do
  
  do j = 1, nsoil         ! soil
     mliq(j) =  noahmp401_struc(n)%noahmp401(t)%sh2o(j) * dzsnso(j) * 1000.
     mice(j) = (noahmp401_struc(n)%noahmp401(t)%smc(j) - &
          noahmp401_struc(n)%noahmp401(t)%sh2o(j))  * dzsnso(j) * 1000.
  end do
  
  do j = isnow+1,nsoil    ! all layers
     imelt(j)    = 0
  enddo
  
  do j = 1,nsoil
     if(stc(j) < tfrz) then
        smp = hfus*(tfrz-stc(j))/(grav*stc(j))             !(m)
        supercool(j) = parameters%smcmax(j)*(smp/parameters%psisat(j))**(-1./parameters%bexp(j))
        supercool(j) = supercool(j)*dzsnso(j)*1000.        !(mm)
     end if
  enddo

  do j = isnow+1,nsoil
     if (mice(j) > 0. .and. stc(j) >= tfrz) then  !melting 
        imelt(j) = 1
     endif
     if (mliq(j) > supercool(j) .and. stc(j) < tfrz) then
        imelt(j) = 2
     endif
     
     ! if snow exists, but its thickness is not enough to create a layer
     if (isnow == 0 &
          .and. sneqv > 0. .and. j == 1) then
        if (stc(j) >= tfrz) then
           imelt(j) = 1
        endif
     endif
  enddo

  ! from snowwater 
  snoflow = 0.0
  ponding1 = 0.0
  ponding2 = 0.0  

  if(isnow < 0) &     ! when multi-layer
       call  compact (parameters, nsnow, nsoil, noahmp401_struc(n)%ts,     & !in
                     stc, snice, snliq, zsoil, imelt, ficeold, iloc, jloc, & !in
                     isnow, dzsnso ,zsnso)                                   !inout
  if(isnow < 0) &
       call  combine (parameters, nsnow, nsoil ,iloc, jloc,          & !in
                      isnow, noahmp401_struc(n)%noahmp401(t)%sh2o,   & !inout
                      stc, snice, snliq, dzsnso, sice, snowh, sneqv, & !inout
                      ponding1, ponding2)                              !out
  if(isnow < 0) &        
       call divide (parameters, nsnow, nsoil,      & !in
                   isnow, stc, snice, snliq, dzsnso) !inout

  !set empty snow layers to zero
  do iz = -nsnow+1, isnow
     snice(iz) = 0.
     snliq(iz) = 0.
     stc(iz)   = 0.
     dzsnso(iz)= 0.
     zsnso(iz) = 0.
  enddo
  
  !to obtain equilibrium state of snow in glacier region
  if(sneqv > 2000.) then   ! 2000 mm -> maximum water depth
     bdsnow      = snice(0) / dzsnso(0)
     snoflow     = (sneqv - 2000.)
     snice(0)    = snice(0)  - snoflow
     dzsnso(0)   = dzsnso(0) - snoflow/bdsnow
     !snoflow     = snoflow / dt
  end if

  ! sum up snow mass for layered snow
  if(isnow < 0) then  ! mb: only do for multi-layer
     sneqv = 0.
     snowh = 0.      ! yeosang yoon
     do iz = isnow+1,0
        sneqv = sneqv + snice(iz) + snliq(iz)
        snowh = snowh + dzsnso(iz)             ! yeosang yoon
     enddo
  end if

  ! no snow layer case, limit snow density to 1000
   if (isnow == 0 .and. sneqv > 0. .and. snowh > 0.) then
        bdsnow = sneqv/snowh
        if (bdsnow >= denh2o) then
            snowh  = snowh*(bdsnow/1000.) ! change unit, sneqv=[mm] snowh=[m]
        end if
   end if

  ! reset zsnso and layer thinkness dzsnso
  do iz = isnow+1, 0
     dzsnso(iz) = -dzsnso(iz)
  end do
  
  dzsnso(1) = zsoil(1)
  do iz = 2,nsoil
     dzsnso(iz) = (zsoil(iz) - zsoil(iz-1))
  end do
  
  zsnso(isnow+1) = dzsnso(isnow+1)
  do iz = isnow+2 ,nsoil
     zsnso(iz) = zsnso(iz-1) + dzsnso(iz)
  enddo

  do iz = isnow+1 ,nsoil
     dzsnso(iz) = -dzsnso(iz)
  end do

  ! update state vars
  noahmp401_struc(n)%noahmp401(t)%isnow = isnow
  noahmp401_struc(n)%noahmp401(t)%sneqv = sneqv
  noahmp401_struc(n)%noahmp401(t)%snowh = snowh 

  noahmp401_struc(n)%noahmp401(t)%zss(1:nsnow+&
       nsoil) = zsnso(-nsnow+1:nsoil)
  noahmp401_struc(n)%noahmp401(t)%snowice(1:nsnow) = & 
       snice(-nsnow+1:0) 
  noahmp401_struc(n)%noahmp401(t)%snowliq(1:nsnow)  = &        
       snliq(-nsnow+1:0) 
  noahmp401_struc(n)%noahmp401(t)%tsno(1:nsnow) = stc(-nsnow+1:0)
  noahmp401_struc(n)%noahmp401(t)%tslb(1:nsoil) = stc(1:nsoil)

  deallocate(ficeold)
  deallocate(snice)
  deallocate(snliq)
  deallocate(stc)
  deallocate(imelt)
  deallocate(supercool)
  deallocate(mice)
  deallocate(mliq)
  deallocate(dzsnso)
  deallocate(zsnso)
  deallocate(sice)

end subroutine noahmp401_snow_update
