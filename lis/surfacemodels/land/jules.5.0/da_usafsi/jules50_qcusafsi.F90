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
! !ROUTINE: jules50_qcusafsi
! \label{jules50_qcusafsi}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 21 Jul 2011: James Geiger; Modified for Noah 3.2
! 30 Jan 2015: Yuqiong Liu; added additional QC
! 05 Nov 2018: Yeosang Yoon; Modified for Jules 5.0 and SNODEP
! 08 Jul 2019: Yeosang Yoon; Modified for Jules.5.0 and LDT-SI data
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with USAFSI
! 30 Dec 2019: Yeosang Yoon; updated QC
! 07 Jan 2021: David Mocko; Use snow density from JULES state;
!              Only check SWE/snow maximum when doing update;
!              Turn off DA in high snow depth regions;
!                  (from Eric Kemp's SNODEP qc change)
! 19 Jan 2020: David Mocko; Further tweaks to checking high snow
!              depth to turn off DA (now using attributes file);
!              Added consistency checks to setting maximum SWE
!              and snow depth after DA against the snow density;
!              Changes recommended after discussion with:
!                   Yeosang Yoon, Yonghwan Kwon, Eric Kemp
!
! !INTERFACE:
subroutine jules50_qcusafsi(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod
  use jules50_lsmMod
  use LIS_logMod
  use c_densty,       only: rho_water
  use jules_snow_mod, only: rho_snow_const, l_snowdep_surf

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  QC's the related state prognostic variable objects for
!  USAFSI data assimilation
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  integer                :: t, gid, pft
  integer                :: status
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)

  real                   :: swemax,snodmax
  real                   :: swemin,snodmin

  real                   :: sndens
  logical                :: update_flag(LIS_rc%ngrid(n))

  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(sweField,"Max Value",swemax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sweField,"Min Value",swemin,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(snodField,"Max Value",snodmax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(snodField,"Min Value",snodmin,rc=status)
  call LIS_verify(status)

  update_flag    = .true.

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     if ((snod(t).lt.snodmin).or.(swe(t).lt.swemin)) then
        update_flag(gid) = .false.
     endif

     if ((snod(t).gt.snodmax).or.(swe(t).gt.swemax)) then
        update_flag(gid) = .false.
     endif

  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

!Use the model's snow density from the previous timestep
     pft = jules50_struc(n)%jules50(t)%pft
     sndens = jules50_struc(n)%jules50(t)%rho_snow_grnd(pft)

     if (l_snowdep_surf) then
       sndens = max(rho_snow_const,sndens)
       sndens = min(rho_water     ,sndens)
     else
       sndens = max(   1.0,sndens)
       sndens = min(1000.0,sndens)
     endif

!Update SWE and snow depth
     if (update_flag(gid)) then
        snod(t) = snod(t)
        swe(t)  = snod(t)*sndens

        if (swe(t).gt.swemax) then
           swe(t) = swemax
           snod(t) = swe(t)/sndens
        endif
        if (snod(t).gt.snodmax) then
           snod(t) = snodmax
           swe(t)  = snod(t)*sndens
        endif

        if (swe(t).lt.swemin) then
           swe(t) = swemin
           snod(t) = swe(t)/sndens
        endif
        if (snod(t).lt.snodmin) then
           snod(t) = snodmin
           swe(t)  = snod(t)*sndens
        endif
        
!If the update is unphysical, do not update
     else
        snod(t) = jules50_struc(n)%jules50(t)%snowdepth(pft)
        swe(t)  = jules50_struc(n)%jules50(t)%snow_mass_ij
     endif

  enddo

end subroutine jules50_qcusafsi

