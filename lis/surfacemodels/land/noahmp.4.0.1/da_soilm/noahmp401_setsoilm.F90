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
! 25 Mar 2025: Eric Kemp; Changed to use top soil layer only, overhauled
!   bounds checks and ensemble spread adjustment.  Based on code from
!   Sujay Kumar and Wanshu Nie used for, e.g., HydroGlobe.
!
! !INTERFACE:
subroutine noahmp401_setsoilm(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only: LIS_rc, LIS_domain, LIS_surface
  use LIS_logMod, only: LIS_verify
  use noahMP401_lsmMod, only: noahmp401_struc
  use noahmp_tables_401, only : smcmax_table, smcwlt_table, isice_table

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

  real                   :: min_threshold
  real                   :: max_threshold
  type(ESMF_Field)       :: sm1Field
  real, pointer          :: soilm1(:)
  real                   :: delta
  logical                :: nonzero_spread(LIS_rc%ngrid(n))
  logical                :: bounds_violation(LIS_rc%ngrid(n))
  integer                :: c, r, t, gid
  integer                :: t_unpert, t_first_mem
  integer                :: soiltyp           ! soil type index [-]
  integer                :: status
  logical                :: update_flag(LIS_rc%ngrid(n))
  logical                :: rc
  integer                :: vegtype
  external               :: noahmp401_sm_reorderEnsForOutliers

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in noahmp401_setsoilm")
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmp401_setsoilm")

  update_flag = .true.
  bounds_violation = .false.
  nonzero_spread = .false.

  ! Identify grid points where any tile has ice (either glacier point, or
  ! some frozen soil moisture exists).  For other grid points, identify
  ! if a bounds violation exists and/or if non-zero spread exists across
  ! the tiles (w/r/t unperturbed member).
  do t_first_mem = 1, LIS_rc%npatch(n,LIS_rc%lsm_index), LIS_rc%nensem(n)
     c = LIS_surface(n,LIS_rc%lsm_index)%tile(t_first_mem)%col
     r = LIS_surface(n,LIS_rc%lsm_index)%tile(t_first_mem)%row
     gid = LIS_domain(n)%gindex(c,r)
     t_unpert = t_first_mem + LIS_rc%nensem(n) - 1
     do t = t_first_mem, t_unpert
        vegtype = noahmp401_struc(n)%noahmp401(t)%vegetype
        if (vegtype .eq. isice_table) then ! Glacier point
           update_flag(gid) = .false.
           exit
        endif
        ! Check if any frozen soil moisture is detected (not just glacier)
        delta = noahmp401_struc(n)%noahmp401(t)%smc(1) - &
             noahmp401_struc(n)%noahmp401(t)%sh2o(1)
        if (delta > 0) then
           update_flag(gid) = .false.
           exit
        end if
        ! For remainder, check if DA soil moisture is out of bounds.
        soiltyp = noahmp401_struc(n)%noahmp401(t)%soiltype
        max_threshold = smcmax_table(soiltyp)
        min_threshold = smcwlt_table(soiltyp)
        if (soilm1(t).lt.min_threshold .or.&
             soilm1(t).gt.max_threshold) then
           bounds_violation(gid) = .true.
        endif
        ! For remainder, check if DA soil moisture tiles have non-zero
        ! spread.
        if (soilm1(t) .ne. soilm1(t_unpert)) then
           nonzero_spread(gid) = .true.
        endif
     end do ! t loop
  enddo ! t_first_mem loop

  ! For non-ice grid points, attempt to rescale DA ensembles if
  ! bound violation is detected and non-zero spread exists.  If problem
  ! occurs, toss the DA values.
  do t_first_mem = 1, LIS_rc%npatch(n,LIS_rc%lsm_index), LIS_rc%nensem(n)
     gid = LIS_domain(n)%gindex( &
          LIS_surface(n,LIS_rc%lsm_index)%tile(t_first_mem)%col, &
          LIS_surface(n,LIS_rc%lsm_index)%tile(t_first_mem)%row)
     rc = .true.
     if (update_flag(gid) .and. bounds_violation(gid) .and. &
          nonzero_spread(gid)) then
        t_unpert = t_first_mem + LIS_rc%nensem(n) - 1
        call noahmp401_sm_reorderEnsForOutliers( &
             LIS_rc%nensem(n), &
             soilm1(t_first_mem:t_unpert), &
             MIN_THRESHOLD, MAX_THRESHOLD, rc)
     endif
     if (.not.rc) then ! Problem occurred when rescaling
        update_flag(gid) = .false.
     endif
  enddo

  ! Apply the (possibly rescaled) DA updates to non-ice points.
  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex( &
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col, &
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
  real, intent(in)     :: minvalue, maxvalue
  logical, intent(inout) :: status

  real                 :: minvT, maxvT, minvG, maxvG
  logical              :: found_good
  integer              :: k
  real                 :: spread_total, spread_good, spread_ratio

  !Ensemble spread (total and with 'good' ensemble members
  minvT = 1E10
  maxvT = -1E10
  minvG = 1E10
  maxvG = -1E10
  status = .true.
  found_good = .false.

  do k=1,nensem

     if(statevec(k).lt.minvT) then
        minvT = statevec(k)
     endif
     if(statevec(k).gt.maxvT) then
        maxvT = statevec(k)
     endif

     if(statevec(k).gt.minvalue.and.statevec(k).lt.maxvalue) then
        found_good = .true.
        if(statevec(k).lt.minvG) then
           minvG = statevec(k)
        endif
        if(statevec(k).gt.maxvG) then
           maxvG = statevec(k)
        endif
     endif
  enddo

  if (found_good) then
     spread_total = (maxvT - minvT)
     spread_good  = (maxvG - minvG)
     spread_ratio = spread_good/spread_total
     do k=1, nensem-1
        statevec(k) = statevec(nensem) + &
             (statevec(k) - statevec(nensem))*spread_ratio
     enddo
  else
     !all members are unphysical.
     statevec = minvalue
     status = .false.
  endif

end subroutine noahmp401_sm_reorderEnsForOutliers
