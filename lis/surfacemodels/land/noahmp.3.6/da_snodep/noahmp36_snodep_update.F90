!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp36_snodep_update
! \label{noahmp36_snodep_update}
!
! !REVISION HISTORY:
!  13 Aug 2017: Sujay Kumar; Initial specification
!
! !INTERFACE:
subroutine noahmp36_snodep_update(n, t, dsneqv, dsnowh)

  use LIS_coreMod
  use NoahMP36_lsmMod
  use module_sf_noahlsm_36 
  use NOAHMP_ROUTINES_36
  
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
  integer :: iloc, jloc
!  real    :: smcmax,psisat,bexp
  real    :: smp,sneqv,snowh
  real    :: sneqv1,snowh1
  real    :: ponding1,ponding2
  integer :: newnode
  integer :: isnow, nsoil, nsnow,soiltyp

! local
  real    :: SNOFLOW, BDSNOW

  isnow = noahmp36_struc(n)%noahmp36(t)%isnow
  nsoil = noahmp36_struc(n)%nsoil
  nsnow = noahmp36_struc(n)%nsnow

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

!  soiltyp = noahmp36_struc(n)%noahmp36(t)%soiltype        
!  smcmax = maxsmc (soiltyp)
!  bexp   = bb (soiltyp)
!  psisat = satpsi (soiltyp)
  
  sneqv = noahmp36_struc(n)%noahmp36(t)%sneqv
  snowh = noahmp36_struc(n)%noahmp36(t)%snowh

  zsnso(-nsnow+1:nsoil) = noahmp36_struc(n)%noahmp36(t)%zss(1:nsnow+nsoil) 

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
  zsoil(1) = -NOAHMP36_struc(n)%sldpth(1)
  do i = 2, nsoil
     zsoil(i) = zsoil(i-1) - NOAHMP36_struc(n)%sldpth(i)
  enddo


  ! state variables 
  snice(-nsnow+1:0) = &
       NOAHMP36_struc(n)%noahmp36(t)%snowice(1:nsnow)
  snliq(-nsnow+1:0) = &
       NOAHMP36_struc(n)%noahmp36(t)%snowliq(1:nsnow) 
  stc(-nsnow+1:nsoil) = &
       NOAHMP36_struc(n)%noahmp36(t)%sstc(1:nsnow+&
       nsoil) 
! from snowfall routine
  ! creating a new layer
 
  IF(ISNOW == 0.and.(dsneqv.gt.0.and.dsnowh.gt.0))  THEN
     SNOWH = SNOWH + dsnowh
     SNEQV = SNEQV + dsneqv
  END IF
  
  NEWNODE = 0 

  IF(ISNOW == 0 .AND. SNOWH >= 0.025.and.&
       (dsneqv.gt.0.and.dsnowh.gt.0))  THEN !MB: change limit
     !    IF(ISNOW == 0  .AND. QSNOW>0. .AND. SNOWH >= 0.05) THEN
     ISNOW    = -1
     NEWNODE  =  1
     DZSNSO(0)= SNOWH
     SNOWH    = 0.
     STC(0)   = MIN(273.16, NOAHMP36_struc(n)%noahmp36(t)%sfctmp)   ! temporary setup
     SNICE(0) = SNEQV
     SNLIQ(0) = 0.
  END IF
  
  ! snow with layers

  IF(ISNOW <  0 .AND. NEWNODE == 0 .and. &
       (dsneqv.gt.0.and.dsnowh.gt.0)) then
     SNICE(ISNOW+1)  = SNICE(ISNOW+1)   + dsneqv
     DZSNSO(ISNOW+1) = DZSNSO(ISNOW+1)  + dsnowh
  ENDIF
  
  if(dsneqv.lt.0.and.dsnowh.lt.0) then
     snowh1 = snowh + dsnowh
     sneqv1 = sneqv + dsneqv
     if(snowh1.ge.0.and.sneqv1.ge.0) then         
        SNOWH = SNOWH + dsnowh
        SNEQV = SNEQV + dsneqv
! Update dzsnso
! how do you determine the thickness of a layer?
        if(snowh.le.dzsnso(0)) then 
           isnow = 0
           dzsnso(-nsnow+1:(isnow-1)) = 0 
           dzsnso(isnow) = snowh
        elseif(snowh.le.(dzsnso(0)+dzsnso(-1))) then 
           isnow = -1
           dzsnso(-nsnow+1:(isnow-1)) = 0 
           dzsnso(isnow) = snowh -dzsnso(isnow+1)
        elseif(snowh.le.(dzsnso(0)+dzsnso(-1)+dzsnso(-2))) then 
           isnow = -2
           dzsnso(-nsnow+1:(isnow-2)) = 0 
           dzsnso(isnow) = snowh -dzsnso(isnow+2)
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

  sice(:) = max(0.0, NOAHMP36_struc(n)%noahmp36(t)%smc(:)&
       - NOAHMP36_struc(n)%noahmp36(t)%sh2o(:))   

!imelt
  do j = -nsnow+1, nsoil
     supercool(j) = 0.0
  end do
  
  do j = isnow+1,0       ! all layers
     mice(j) = snice(j)
     mliq(j) = snliq(j)
  end do
  
  do j = 1, nsoil               ! soil
     mliq(j) =  NOAHMP36_struc(n)%noahmp36(t)%sh2o(j) * dzsnso(j) * 1000.
     mice(j) = (NOAHMP36_struc(n)%noahmp36(t)%smc(j) - &
          NOAHMP36_struc(n)%noahmp36(t)%sh2o(j))  * dzsnso(j) * 1000.
  end do
  
  do j = isnow+1,nsoil       ! all layers
     imelt(j)    = 0
  enddo
  
  do j = 1,nsoil
!     if (opt_frz == 1) then
! Assuming the use of option 1 for now
        if(stc(j) < tfrz) then
           smp = hfus*(tfrz-stc(j))/&
                (grav*stc(j))             !(m)
           supercool(j) = smcmax*(smp/psisat)**(-1./bexp)
           supercool(j) = supercool(j)*dzsnso(j)*1000.        !(mm)
        end if
!     end if
!     if (opt_frz == 2) then
!        call frh2o (supercool(j),&
!             NOAHMP36_struc(n)%noahmp36(t)%sstc(j),&
!             NOAHMP36_struc(n)%noahmp36(t)%smc(j),&
!             NOAHMP36_struc(n)%noahmp36(t)%sh2o(j))
!        supercool(j) = supercool(j)*dzsnso(j)*1000.        !(mm)
!     end if
  enddo

  do j = isnow+1,nsoil
     if (mice(j) > 0. .and. stc(j) >= tfrz) then  !melting 
        imelt(j) = 1
     endif
     if (mliq(j) > supercool(j) .and. stc(j) < tfrz) then
        imelt(j) = 2
     endif
     
     ! If snow exists, but its thickness is not enough to create a layer
     if (isnow == 0 &
          .and. sneqv > 0. .and. j == 1) then
        if (stc(j) >= tfrz) then
           imelt(j) = 1
        endif
     endif
  enddo


  if(isnow < 0) &     
       call  compact (nsnow  ,&
       nsoil  ,&
       noahmp36_struc(n)%dt     ,&
       stc    , &
       snice  , & 
       snliq  , &
       zsoil  , & 
       imelt  , &
       ficeold, & 
       iloc   , jloc ,& !in
       isnow  ,&
       dzsnso ,zsnso  )                   !inoutl

  if(isnow < 0) &
       call  combine (nsnow  ,&
       nsoil  ,iloc   ,jloc   ,         & !in
       isnow  ,&
       noahmp36_struc(n)%noahmp36(t)%sh2o   ,&
       stc    ,snice  ,snliq  , & !inout
       dzsnso ,sice   ,snowh  ,sneqv  ,         & !inout
       ponding1       ,ponding2)                  !out


  if(isnow < 0) &        
       call divide (nsnow  ,&
       nsoil  ,                         & !in
       isnow  ,&
       stc    ,&
       snice  ,snliq  ,dzsnso )   !inout
  
#if 0 

  call  snowh2o (nsnow  ,&
       nsoil  ,&
       noahmp36_struc(n)%dt     ,qsnfro ,qsnsub , & !in 
       qrain  ,iloc   ,jloc   ,                 & !in
       isnow  ,&
       dzsnso ,snowh  ,sneqv  ,snice  , & !inout
       snliq  ,noahmp36_struc(n)%noahmp36(t)%sh2o   ,&
       sice   ,noahmp36_struc(n)%noahmp36(t)%stc    ,         & !inout
       qsnbot ,ponding1       ,ponding2)           !out
  
#endif

!set empty snow layers to zero

  do iz = -nsnow+1, isnow
     snice(iz) = 0.
     snliq(iz) = 0.
     stc(iz)   = 0.
     dzsnso(iz)= 0.
     zsnso(iz) = 0.
  enddo
  
!to obtain equilibrium state of snow in glacier region
  IF(SNEQV > 2000.) THEN   ! 2000 mm -> maximum water depth
     BDSNOW      = SNICE(0) / DZSNSO(0)
     SNOFLOW     = (SNEQV - 2000.)
     SNICE(0)    = SNICE(0)  - SNOFLOW
     DZSNSO(0)   = DZSNSO(0) - SNOFLOW/BDSNOW
!     SNOFLOW     = SNOFLOW / DT
  END IF

! sum up snow mass for layered snow
  IF(ISNOW < 0) THEN  ! MB: only do for multi-layer
     SNEQV = 0.
     SNOWH = 0.      ! Yeosang Yoon
     DO IZ = ISNOW+1,0
        SNEQV = SNEQV + SNICE(IZ) + SNLIQ(IZ)
        SNOWH = SNOWH + DZSNSO(IZ)             ! Yeosang Yoon
     ENDDO
  END IF

! Yeosag Yoon, no snow layer case, limit snow density to 1000
   IF (ISNOW == 0 .AND. SNEQV > 0. .AND. SNOWH > 0.) THEN
        BDSNOW = SNEQV/SNOWH
        IF (BDSNOW >= DENH2O) THEN
            SNOWH  = SNOWH*(BDSNOW/1000.) ! change unit, SNEQV=[mm] SNOWH=[m]
        END IF
   END IF

! Reset ZSNSO and layer thinkness DZSNSO

  DO IZ = ISNOW+1, 0
     DZSNSO(IZ) = -DZSNSO(IZ)
  END DO
  
  DZSNSO(1) = ZSOIL(1)
  DO IZ = 2,NSOIL
     DZSNSO(IZ) = (ZSOIL(IZ) - ZSOIL(IZ-1))
  END DO
  
  ZSNSO(ISNOW+1) = DZSNSO(ISNOW+1)
  DO IZ = ISNOW+2 ,NSOIL
     ZSNSO(IZ) = ZSNSO(IZ-1) + DZSNSO(IZ)
  ENDDO

  DO IZ = ISNOW+1 ,NSOIL
     DZSNSO(IZ) = -DZSNSO(IZ)
  END DO

  noahmp36_struc(n)%noahmp36(t)%isnow = isnow
  noahmp36_struc(n)%noahmp36(t)%sneqv = sneqv
  noahmp36_struc(n)%noahmp36(t)%snowh = snowh 

  NOAHMP36_struc(n)%noahmp36(t)%zss(1:nsnow+&
       nsoil) = & 
       zsnso(-nsnow+1:nsoil)

  NOAHMP36_struc(n)%noahmp36(t)%snowice(1:nsnow) = & 
       snice(-nsnow+1:0) 
  NOAHMP36_struc(n)%noahmp36(t)%snowliq(1:nsnow)  = &        
       snliq(-nsnow+1:0) 
  NOAHMP36_struc(n)%noahmp36(t)%sstc(1:nsnow+&
       nsoil)   = & 
       stc(-nsnow+1:nsoil) 

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

end subroutine noahmp36_snodep_update
