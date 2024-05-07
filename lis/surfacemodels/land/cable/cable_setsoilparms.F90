!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: cable_setsoilparms
! \label{cable_setsoilparms}
!
! !REVISION HISTORY:
!  25 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!  13 Sep 2011: Claire Carouge (ccc), CABLE LSM improvements
!
! !INTERFACE:
subroutine cable_setsoilparms
! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_surface
  use LIS_soilsMod,       only : LIS_soils
  use LIS_logMod,         only : LIS_logunit
  use cable_lsmMod
!
! !DESCRIPTION:
!  This subroutine retrieves CABLE soil parameters.
!  The current implementation uses a table-based lookup based
!  on soil texture classes to initialize the following parameters.
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_mapSoilType](\ref{LIS_mapSoilType}) \newline
!    Method to derive the soil texture type from the
!    sand, silt, and clay fractions
!  \end{description}
!EOP
  implicit none

  character(len=70), dimension(:), allocatable :: soil_desc
  integer :: n,t,a,nsoilt
  integer :: ibp2

  real, dimension(:), allocatable :: soilin_silt
  real, dimension(:), allocatable :: soilin_clay
  real, dimension(:), allocatable :: soilin_sand
  real, dimension(:), allocatable :: soilin_swilt
  real, dimension(:), allocatable :: soilin_sfc
  real, dimension(:), allocatable :: soilin_ssat
  real, dimension(:), allocatable :: soilin_bch
  real, dimension(:), allocatable :: soilin_hyds
  real, dimension(:), allocatable :: soilin_sucs
  real, dimension(:), allocatable :: soilin_rhosoil
  real, dimension(:), allocatable :: soilin_css

  do n=1,LIS_rc%nnest

     write(LIS_logunit,*)                                          &
          'MSG: cable_setsoilparms -- reading soil files'

     if (LIS_rc%usetexturemap(n).eq."none") then
        ! As of now, CABLE only uses the Zobler lookup table.  If another
        ! lookup table is added in the future (such as STATSGO), change the
        ! "1" in the call to the following subroutine to a flag variable. - dmm
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           call LIS_mapSoilType(1,&
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%sand, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%clay, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%silt,&
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%soilt)
        enddo
     endif

     !-----------------------------------------------------------------------
     ! Set CABLE soil type at tile from the LIS domain
     !-----------------------------------------------------------------------
     if (cable_struc(n)%fixedsoiltype.ne.0) then
        write(LIS_logunit,*) 'Fixing CABLE soil index to type: ',  &
             cable_struc(n)%fixedsoiltype
     endif

     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        if (cable_struc(n)%fixedsoiltype.eq.0) then
           cable_struc(n)%cable(t)%soiltype =                      &
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%soilt
           ! Kluge to not let soil texture 10 through - dmm.
           !               if (cable_struc(n)%cable(t)%soiltype.eq.8) &
           !                   cable_struc(n)%cable(t)%soiltype = 2
           if (cable_struc(n)%cable(t)%soiltype.eq.10) &
                cable_struc(n)%cable(t)%soiltype = 2
           !               if (cable_struc(n)%cable(t)%soiltype.ge.8) &
           !                 print *,'Soil type = ',cable_struc(n)%cable(t)%soiltype
        else
           cable_struc(n)%cable(t)%soiltype =                      &
                cable_struc(n)%fixedsoiltype
        endif

        ! Set undefined points to Zobler soil type index "6"
        if ((cable_struc(n)%cable(t)%soiltype.eq.-9999).or.        &
             (cable_struc(n)%cable(t)%soiltype.eq.0)) then
           cable_struc(n)%cable(t)%soiltype = 6 ! clay loam
        endif
     enddo

     !-----------------------------------------------------------------------
     ! Read in the CABLE Soil Parameter File
     !-----------------------------------------------------------------------
     write(LIS_logunit,*) 'Reading CABLE soil parameter file: ',   &
          trim(cable_struc(n)%sfile)
     open(unit=18,file=cable_struc(n)%sfile,status='old',          &
          access='sequential')

     read(18,*)
     read(18,*)
     read(18,*) nsoilt      ! Number of soil types
     read(18,*)
     read(18,*)
     ! BP added the ALLOCATE statements here. (Dec 2007)
     allocate(soil_desc(nsoilt))
     allocate(soilin_silt(nsoilt),soilin_clay(nsoilt))
     allocate(soilin_sand(nsoilt),soilin_swilt(nsoilt))
     allocate(soilin_sfc(nsoilt),soilin_ssat(nsoilt))
     allocate(soilin_bch(nsoilt),soilin_hyds(nsoilt))
     allocate(soilin_sucs(nsoilt),soilin_rhosoil(nsoilt))
     allocate(soilin_css(nsoilt))
     ! Read description of each soil type
     do a = 1,nsoilt
        read(18,'(8X,A70)') soil_desc(a)
     enddo
     read(18,*)
     read(18,*)
     read(18,*) soilin_silt
     read(18,*) soilin_clay
     read(18,*) soilin_sand
     read(18,*) soilin_swilt
     read(18,*) soilin_sfc
     read(18,*) soilin_ssat
     read(18,*) soilin_bch
     read(18,*) soilin_hyds
     read(18,*) soilin_sucs
     read(18,*) soilin_rhosoil
     read(18,*) soilin_css
     close(18)

     cable_struc(n)%nsoilt = nsoilt

     !-----------------------------------------------------------------------
     ! Assign soil parameters to each tile based on
     ! the type of soil class present in that tile.
     !-----------------------------------------------------------------------
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        ! Parameters that are not spatially dependent
        cable_struc(n)%cable(t)%albsoil= cable_struc(n)%fixedalbsoil
        ! Add %zse to subroutine_readcrd instead of hard-coding it. - dmm
        cable_struc(n)%cable(t)%zse =                              &
             (/.022, .058, .154, .409, 1.085, 2.872/) ! layer thickness nov03
        cable_struc(n)%cable(t)%zshh(1) = 0.5*cable_struc(n)%cable(t)%zse(1)
        cable_struc(n)%cable(t)%zshh(7) = 0.5*cable_struc(n)%cable(t)%zse(6)
        cable_struc(n)%cable(t)%zshh(2:6) = 0.5 * (cable_struc(n)%cable(t)%zse(1:5) + cable_struc(n)%cable(t)%zse(2:6))

        ! Parameters that vary spatially by soil type
        cable_struc(n)%cable(t)%silt =                             &
             soilin_silt(cable_struc(n)%cable(t)%soiltype)
        cable_struc(n)%cable(t)%clay =                             &
             soilin_clay(cable_struc(n)%cable(t)%soiltype)
        cable_struc(n)%cable(t)%sand =                             &
             soilin_sand(cable_struc(n)%cable(t)%soiltype)
        cable_struc(n)%cable(t)%swilt =                            &
             soilin_swilt(cable_struc(n)%cable(t)%soiltype)
        cable_struc(n)%cable(t)%sfc =                              &
             soilin_sfc(cable_struc(n)%cable(t)%soiltype)
        cable_struc(n)%cable(t)%rhosoil =                          &
             soilin_rhosoil(cable_struc(n)%cable(t)%soiltype)
        cable_struc(n)%cable(t)%css =                              &
             soilin_css(cable_struc(n)%cable(t)%soiltype)
     enddo

     ! Set tile soil porosity
     if (LIS_rc%useporositymap(n).eq."none") then ! default, from look-up table
        do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           cable_struc(n)%cable(t)%ssat =                          &
                soilin_ssat(cable_struc(n)%cable(t)%soiltype)
        enddo
     else
        do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           cable_struc(n)%cable(t)%ssat =                          &
                LIS_soils(n)%porosity(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,   &
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row,1)
        enddo
     endif

     ! Set tile soil b-parameter
     if (LIS_rc%usebexpmap(n).eq."none") then ! default, from look-up table
        do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           cable_struc(n)%cable(t)%bch =                           &
                soilin_bch(cable_struc(n)%cable(t)%soiltype)
        enddo
     else
        do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           cable_struc(n)%cable(t)%bch =                           &
                LIS_soils(n)%bexp(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,       &
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
        enddo
     endif

     ! Set tile soil K-sat
     if (LIS_rc%useksatmap(n).eq."none") then ! default, from look-up table
        do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           cable_struc(n)%cable(t)%hyds =                          &
                soilin_hyds(cable_struc(n)%cable(t)%soiltype)
        enddo
     else
        do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           cable_struc(n)%cable(t)%hyds =                          &
                LIS_soils(n)%ksat(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,       &
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
        enddo
     endif

     ! Set tile soil Psi-sat
     if (LIS_rc%usepsisatmap(n).eq."none") then ! default, from look-up table
        do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           cable_struc(n)%cable(t)%sucs =                          &
                soilin_sucs(cable_struc(n)%cable(t)%soiltype)
        enddo
     else
        do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           cable_struc(n)%cable(t)%sucs =                          &
                LIS_soils(n)%psisat(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,     &
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
        enddo
     endif

     ! Calculate pwb_min
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        ibp2 = nint(cable_struc(n)%cable(t)%bch)+2
        cable_struc(n)%cable(t)%pwb_min = ( cable_struc(n)%cable(t)%swilt / &
             cable_struc(n)%cable(t)%ssat )  &
             **ibp2
     enddo

     deallocate(soil_desc)
     deallocate(soilin_silt,soilin_clay,soilin_sand)
     deallocate(soilin_swilt,soilin_sfc,soilin_ssat)
     deallocate(soilin_bch,soilin_hyds,soilin_sucs)
     deallocate(soilin_rhosoil,soilin_css)
  enddo

end subroutine cable_setsoilparms
