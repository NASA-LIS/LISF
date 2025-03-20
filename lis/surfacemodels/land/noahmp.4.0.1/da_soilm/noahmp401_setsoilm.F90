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
! !ROUTINE: noahmp401_setsoilm
!  \label{noahmp401_setsoilm}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li; created for Noah-MP4.0.1
! 20 Mar 2025: Eric Kemp; Major rewrite -- uses top soil layer only,
!   overhauled bounds checks and ensemble spread adjusted.  Based on
!   code from Sujay Kumar and Wanshu Nie used for, e.g., HydroGlobe.
!
! !INTERFACE:
subroutine noahmp401_setsoilm(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only: LIS_rc, LIS_domain, LIS_surface
  use LIS_logMod, only: LIS_verify
  use noahMP401_lsmMod, only: noahmp401_struc
  use noahmp_tables_401, only: smcmax_table, smcwlt_table, isice_table

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State), intent(in) :: LSM_State
!
! !DESCRIPTION:
!
!  This routine assigns the soil moisture prognostic variables to NoahMP's
!  model space.
!
!EOP

  real                   :: MIN_THRESHOLD
  real                   :: MAX_THRESHOLD
  real                   :: sm_threshold
  type(ESMF_Field)       :: sm1Field
  real, pointer          :: soilm1(:)
  real                   :: delta
  logical                :: ens_diff(LIS_rc%ngrid(n))
  logical                :: bounds_violation(LIS_rc%ngrid(n))
  integer                :: i, c, r, t, m, gid
  integer                :: SOILTYP           ! soil type index [-]
  integer                :: vegtype
  integer                :: status
  logical                :: update_flag(LIS_rc%ngrid(n))
  logical                :: rc
  external :: noahmp401_sm_reorderEnsForOutliers

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in noahmp401_setsoilm")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmp401_setsoilm")

  update_flag = .true.
  bounds_violation = .false.
  ens_diff = .false.

  ! Identify grid points where ice exists (either glacier point, or
  ! some frozen soil moisture exists).  Note that grid point is flagged
  ! if *any* ensemble member indicates ice.
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     gid = LIS_domain(n)%gindex(c,r)
     if (.not. update_flag(gid)) cycle ! Skip if ice already detected
     vegtype = NOAHMP401_struc(n)%noahmp401(t)%vegetype
     if (vegtype .eq. isice_table) then ! Glacier
        update_flag(gid) = .false.
        cycle
     endif
     ! Check if any frozen soil moisture is detected (not just glacier)
     delta = noahmp401_struc(n)%noahmp401(t)%smc(1) - &
          noahmp401_struc(n)%noahmp401(t)%sh2o(1)
     if (abs(delta).ne.0) then
        update_flag(gid) = .false.
        cycle
     end if
  enddo

  ! Toss DA analysis if ice detected.
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     gid = LIS_domain(n)%gindex(c,r)
     if (.not. update_flag(gid)) then
        soilm1(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
     endif
  enddo

  ! Find grid points with out-of-bounds DA values, and where DA values
  ! differ from control member
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     gid = LIS_domain(n)%gindex(c,r)
     SOILTYP = NOAHMP401_struc(n)%noahmp401(t)%soiltype
     MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP)
     MIN_THRESHOLD = SMCWLT_TABLE(SOILTYP)
     sm_threshold = MAX_THRESHOLD - 0.02
     if (soilm1(t).lt.MIN_THRESHOLD.or.&
          soilm1(t).gt.MAX_THRESHOLD) then
        bounds_violation(gid) = .true.
     endif
     if (soilm1(t).ne.soilm1(i*LIS_rc%nensem(n))) then
        ens_diff(gid) = .true.
     endif
  enddo

  ! Attempt to rescale DA ensembles at non-ice points, if bound violation
  ! is detected and ensemble diff exists.  If not possible, toss the
  ! DA values.
  do i=1,LIS_rc%ngrid(n)
     rc = .true.
     if (bounds_violation(i) .and. ens_diff(i) .and. update_flag(i)) then
        call noahmp401_sm_reorderEnsForOutliers(LIS_rc%nensem(n),&
             soilm1((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_THRESHOLD, MAX_THRESHOLD, rc)
     endif
     if (.not.rc) then
        do m=1,LIS_rc%nensem(n)
           t = (i-1)*LIS_rc%nensem(n)+m
           soilm1(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
        enddo
     endif
  enddo

  ! Apply the DA updates
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
     if (update_flag(gid)) then
        delta = soilm1(t) - NOAHMP401_struc(n)%noahmp401(t)%smc(1)
        NOAHMP401_struc(n)%noahmp401(t)%smc(1) = soilm1(t)
        NOAHMP401_struc(n)%noahmp401(t)%sh2o(1) = &
             NOAHMP401_struc(n)%noahmp401(t)%sh2o(1) + delta
     endif
  enddo

end subroutine noahmp401_setsoilm

subroutine noahmp401_sm_reorderEnsForOutliers(nensem, statevec, &
     minvalue, maxvalue, status)

  implicit none

  integer, intent(in)  :: nensem
  real, intent(inout)  :: statevec(nensem)
  real, intent(in)     :: minvalue,maxvalue
  logical, intent(inout) :: status

  real                 :: minvT, maxvT, minvG, maxvG
  integer              :: k
  real                 :: spread_total, spread_good, spread_ratio

  !Ensemble spread (total and with 'good' ensemble members
  minvT = 1E10
  maxvT = -1E10
  minvG = 1E10
  maxvG = -1E10
  status = .true.

  do k=1, nensem
     if (statevec(k).lt.minvT) then
        minvT = statevec(k)
     endif
     if (statevec(k).gt.maxvT) then
        maxvT = statevec(k)
     endif

     if (statevec(k).gt.minvalue.and.statevec(k).lt.maxvalue) then
        if (statevec(k).lt.minvG) then
           minvG = statevec(k)
        endif
        if (statevec(k).gt.maxvG) then
           maxvG = statevec(k)
        endif
     endif
  enddo

  if (minvG.eq.1E10.and.maxvG.eq.-1E10) then
     !all members are unphysical.
     statevec = minvalue
     status = .false.
  else
     spread_total = (maxvT - minvT)
     spread_good  = (maxvG - minvG)
     spread_ratio = spread_good/spread_total
     !rescale the ensemble
     do k=1,nensem-1
        statevec(k) = statevec(nensem) + &
             (statevec(k) - statevec(nensem))*spread_ratio
     enddo
  endif

end subroutine noahmp401_sm_reorderEnsForOutliers
