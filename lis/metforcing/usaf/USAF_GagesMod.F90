!
! MODULE: USAF_GagesMod
!
! DESCRIPTION: Contains structure and related routines for storing and
! updating new database of rain gauge reports from USAF. Includes
! routines to reconcile reports by intercomparison.
!
! AUTHOR: Eric Kemp, SSAI/NASA GSFC
!
! REFERENCES:
! Durre, I, M J Menne, B E Gleason, T G Houston, and R S Vose, 2010:
!   Comprehensive automated quality assurance of daily surface observations.
!   J Appl Meteor Climatol, 49, 1615-1633,
!   https://doi.org/10.1175/2010JAMC2375.1
! Qi, Y, S Martinaitis, J Zhang, and S Cocks, 2016: A real-time automated
!   quality control of hourly rain gauge data based on multiple sensors in
!   MRMS system. J Hydrometeor, 17, 1675-1691,
!   https://doi.org/10.1175/JHM-D-15-0188.1
! WMO, 2019: Manual on Codes, International Codes, Volume I.1, Part A --
!   Alphanumeric Codes. WMO No. 306,
!   https://library.wmo.int/doc_num.php?explnum_id=10235
! WMO, 2021: Manual on Codes, International Codes, Volume I.2, Part B --
!   Binary Codes. WMO No. 306,
!   https://library.wmo.int/doc_num.php?explnum_id=10722
! WMO, 2018: Manual on Codes, Regional Codes and National Coding Practices,
!   Volume II. WMO No. 306,
!   https://library.wmo.int/doc_num.php?explnum_id=5730
! WMO, 2012: Weather Reporting, Volume A -- Observing Stations.  WMO No. 9,
!   https://library.wmo.int/doc_num.php?explnum_id=9896
!
! REVISION HISTORY:
! 30 Aug 2021...First developmental version committed to git.
! 03 Aug 2023...Added 1-hr and 2-hr accumulations, date/times of each
!   observation, and ob-specific past weather durations.
! 23 Aug 2023...Added both WMO and FIPS country codes.
! 07 Sep 2023...Fixed WMO block check (values in preobs files are typically
!   <= 100, so these are multipied by 1000 or 10000 to check for country
!   ranges).
!------------------------------------------------------------------------------

module USAF_GagesMod

  ! Defaults
  implicit none
  private

  type, public :: USAF_Gages_t
     private
     logical :: cdms_flag
     character(10) :: date10
     integer :: nobs
     character(14), allocatable :: YYYYMMDDhhmmss(:)
     character(32), allocatable :: networks(:)
     character(32), allocatable :: platforms(:)
     character(2), allocatable :: wmocode_id(:)
     character(2), allocatable :: fipscode_id(:)
     integer, allocatable :: wmonumbers(:)
     integer, allocatable :: bsn(:)
     ! NOTE: Latitudes and longitudes are in hundreths of degrees
     ! (100 = 1.00 deg)
     integer, allocatable :: lats(:)
     integer, allocatable :: lons(:)
     ! NOTE:  Accumulations are in tenths of mm (10 = 1 mm)
     integer, allocatable :: amts24(:)
     integer, allocatable :: amts21(:) ! Non-standard
     integer, allocatable :: amts18(:)
     integer, allocatable :: amts15(:)
     integer, allocatable :: amts12(:)
     ! "Original", non-filled or reconciled 12hr amount, needed for
     ! some 12Z South America sanity checks.  These are not written
     ! to the presav2 file.
     integer, allocatable :: amts12_orig(:)
     integer, allocatable :: amts09(:)
     integer, allocatable :: amts06(:)
     ! "Original", non-filled or reconciled 06hr amount, needed for
     ! some 12Z South America sanity checks.  These are not written
     ! to the presav2 file.
     integer, allocatable :: amts06_orig(:)
     integer, allocatable :: amts03(:)
     integer, allocatable :: amts02(:)
     integer, allocatable :: amts01(:)
     integer, allocatable :: amts00(:) ! Non-standard accumulations
     integer, allocatable :: durations(:) ! Durations of non-std accums
     integer, allocatable :: preswx(:)
     integer, allocatable :: pastwx(:)
     integer, allocatable :: pastwx_durations(:)
     character(32), allocatable :: unique_networks(:)
     integer, allocatable :: firsts(:) ! Starting indices for each network
     integer, allocatable :: lasts(:)  ! Ending indices for each network
     integer :: num_unique_networks
   contains
     procedure :: new => USAF_gages_new
     procedure :: read_data => USAF_gages_read_data
     procedure :: delete => USAF_gages_delete
     procedure :: check_gross_errors => &
          USAF_gages_check_gross_errors
     procedure :: use_misc_precip => USAF_gages_use_misc_precip
     procedure :: reconcile_self => USAF_gages_reconcile_self
     procedure :: reconcile_gages01 => USAF_gages_reconcile_gages01
     procedure :: reconcile_gages02 => USAF_gages_reconcile_gages02
     procedure :: reconcile_gages03 => USAF_gages_reconcile_gages03
     procedure :: reconcile_gages06 => USAF_gages_reconcile_gages06
     procedure :: reconcile_gages09 => USAF_gages_reconcile_gages09
     procedure :: reconcile_gages12 => USAF_gages_reconcile_gages12
     procedure :: reconcile_gages15 => USAF_gages_reconcile_gages15
     procedure :: reconcile_gages18 => USAF_gages_reconcile_gages18
     procedure :: reconcile_gages21 => USAF_gages_reconcile_gages21
     procedure :: correct_region3_12Z => USAF_gages_correct_region3_12Z
     procedure :: use_preswx_pastwx => USAF_gages_use_preswx_pastwx
     procedure :: fill_gaps => USAF_gages_fill_gaps
     procedure :: write_data => USAF_gages_write_data
     procedure :: set_pastwx_durations => USAF_gages_set_pastwx_durations
     procedure :: copy_to_usaf_obsdata => USAF_copy_to_usaf_obsdata
  end type USAF_Gages_t

  ! Local parameters
  integer, parameter :: MAX_UNIQUE_NETWORKS = 25
  integer, parameter :: MISSING = -99999999
  character(4), parameter :: MISSING_NAME = "NULL"

  ! Old CDMS replacements for WMO numbers
  ! The following two are guesses
  integer, parameter :: CDMS_SWEDEN_LOWLIMIT =    020010
  integer, parameter :: CDMS_SWEDEN_HIGHLIMIT =   026999
  ! The following two are guesses
  integer, parameter :: CDMS_DENMARK_LOWLIMIT =   060010
  integer, parameter :: CDMS_DENMARK_HIGHLIMIT =  061999
  integer, parameter :: CDMS_RUSSIA_LOWLIMIT =    200000
  integer, parameter :: CDMS_RUSSIA_HIGHLIMIT =   390009
  integer, parameter :: CDMS_INDIA_LOWLIMIT =     420010
  integer, parameter :: CDMS_INDIA_HIGHLIMIT =    433999
  integer, parameter :: CDMS_SRILANKA_LOWLIMIT =  434000
  integer, parameter :: CDMS_SRILANKA_HIGHLIMIT = 434979
  integer, parameter :: CDMS_S_AMER_LOWLIMIT =    800000
  integer, parameter :: CDMS_S_AMER_HIGHLIMIT =   889999

  ! JMOBS WMO Numbers.
  integer, parameter, public :: JMOBS_SWEDEN_LOWLIMIT =    02001
  integer, parameter, public :: JMOBS_SWEDEN_HIGHLIMIT =   02699
  integer, parameter, public :: JMOBS_DENMARK_LOWLIMIT =   06001
  integer, parameter, public :: JMOBS_DENMARK_HIGHLIMIT =  06199
  integer, parameter, public :: JMOBS_RUSSIA_LOWLIMIT =    20000
  integer, parameter, public :: JMOBS_RUSSIA_HIGHLIMIT =   39000
  integer, parameter, public :: JMOBS_INDIA_LOWLIMIT =     42000
  integer, parameter, public :: JMOBS_INDIA_HIGHLIMIT =    43399
  integer, parameter, public :: JMOBS_SRILANKA_LOWLIMIT =  43400
  integer, parameter, public :: JMOBS_SRILANKA_HIGHLIMIT = 43497
  integer, parameter, public :: JMOBS_S_AMER_LOWLIMIT =    80000
  integer, parameter, public :: JMOBS_S_AMER_HIGHLIMIT =   88999

  ! Maximum allowed accumulations
  ! Values below are in tenths of mm (10 = 1 mm)
  ! Qi et al. (2016) use 50.8 mm (2 in) for 1 hour accum, based on
  ! subjective evaluation in CONUS regions with poor radar coverage;
  ! they judged approximately 97% of observations above this amount were
  ! erroneous.
  integer, parameter :: MAX_PCP_1HR =    508 !  2 inches
  integer, parameter :: MAX_PCP_2HR =   1016 !  4 inches
  integer, parameter :: MAX_PCP_3HR =   1524 !  6 inches
  integer, parameter :: MAX_PCP_6HR =   3048 ! 12 inches
  integer, parameter :: MAX_PCP_9HR =   4572 ! 18 inches
  integer, parameter :: MAX_PCP_12HR =  6096 ! 24 inches
  integer, parameter :: MAX_PCP_15HR =  7620 ! 30 inches
  integer, parameter :: MAX_PCP_18HR =  9144 ! 36 inches
  integer, parameter :: MAX_PCP_21HR = 10668 ! 42 inches
  integer, parameter :: MAX_PCP_24HR = 12192 ! 48 inches

  ! Durre et al (2010) use 1828.8 mm (72 in) for 24-hr limit, based on WMO
  ! world record for 24-hr accumulation.
  ! integer, parameter :: MAX_PCP_3HR  =  2286 !  9 inches
  ! integer, parameter :: MAX_PCP_6HR  =  4572 ! 18 inches
  ! integer, parameter :: MAX_PCP_9HR  =  6858 ! 27 inches
  ! integer, parameter :: MAX_PCP_12HR =  9144 ! 36 inches
  ! integer, parameter :: MAX_PCP_15HR = 11430 ! 45 inches
  ! integer, parameter :: MAX_PCP_18HR = 13716 ! 54 inches
  ! integer, parameter :: MAX_PCP_21HR = 16002 ! 63 inches
  ! integer, parameter :: MAX_PCP_24HR = 18288 ! 72 inches

  ! Note:  WMO world records. See https://wmo.asu.edu
  ! 1-hr rainfall:  305 mm (12 inch), Holt MO, USA, June 1947
  ! 12-hr rainfall: 1144 mm (45 inch), Foc-Foc, La Reunion, Jan 1966
  ! 24-hr rainfall: 1825 mm (71.8 inch), Foc-Foc, La Reunion, Jan 1966

contains

  ! Constructor for USAF_gages_t type
  subroutine USAF_gages_new(this, date10, &
       nobs, YYYYMMDDhhmmss, networks, platforms, &
       wmocode_id, fipscode_id, &
       wmonumbers, bsn, lats, lons, &
       amts24, amts21, amts18, amts15, amts12, amts09, amts06, amts03, &
       amts02, amts01, &
       amts00, durations, preswx, pastwx)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(out) :: this
    character(10), intent(in) :: date10
    integer, intent(in) :: nobs
    character(14), intent(in) :: YYYYMMDDhhmmss(nobs)
    character(32), intent(in) :: networks(nobs)
    character(32), intent(in) :: platforms(nobs)
    character(2), intent(in) :: wmocode_id(nobs)
    character(2), intent(in) :: fipscode_id(nobs)
    integer, intent(in) :: wmonumbers(nobs)
    integer, intent(in) :: bsn(nobs)
    integer, intent(in) :: lats(nobs)
    integer, intent(in) :: lons(nobs)
    integer, intent(in) :: amts24(nobs)
    integer, intent(in) :: amts21(nobs)
    integer, intent(in) :: amts18(nobs)
    integer, intent(in) :: amts15(nobs)
    integer, intent(in) :: amts12(nobs)
    integer, intent(in) :: amts09(nobs)
    integer, intent(in) :: amts06(nobs)
    integer, intent(in) :: amts03(nobs)
    integer, intent(in) :: amts02(nobs)
    integer, intent(in) :: amts01(nobs)
    integer, intent(in) :: amts00(nobs)
    integer, intent(in) :: durations(nobs)
    integer, intent(in) :: preswx(nobs)
    integer, intent(in) :: pastwx(nobs)

    this%date10 = date10

    this%nobs = nobs

    allocate(this%YYYYMMDDhhmmss(nobs))
    this%YYYYMMDDhhmmss = YYYYMMDDhhmmss

    allocate(this%networks(nobs))
    this%networks = networks

    allocate(this%platforms(nobs))
    this%platforms = platforms

    allocate(this%wmocode_id(nobs))
    this%wmocode_id = wmocode_id

    allocate(this%fipscode_id(nobs))
    this%fipscode_id = fipscode_id

    allocate(this%wmonumbers(nobs))
    this%wmonumbers = wmonumbers

    allocate(this%bsn(nobs))
    this%bsn = bsn

    allocate(this%lats(nobs))
    this%lats = lats

    allocate(this%lons(nobs))
    this%lons = lons

    allocate(this%amts24(nobs))
    this%amts24 = amts24

    allocate(this%amts21(nobs))
    this%amts21 = amts21

    allocate(this%amts18(nobs))
    this%amts18 = amts18

    allocate(this%amts15(nobs))
    this%amts15 = amts15

    allocate(this%amts12(nobs))
    this%amts12 = amts12

    allocate(this%amts12_orig(nobs))
    this%amts12_orig = amts12

    allocate(this%amts09(nobs))
    this%amts09 = amts09

    allocate(this%amts06(nobs))
    this%amts06 = amts06

    allocate(this%amts06_orig(nobs))
    this%amts06_orig = amts06

    allocate(this%amts03(nobs))
    this%amts03 = amts03

    allocate(this%amts02(nobs))
    this%amts02 = amts02

    allocate(this%amts01(nobs))
    this%amts01 = amts01

    allocate(this%amts00(nobs))
    this%amts00 = amts00

    allocate(this%durations(nobs))
    this%durations = durations

    allocate(this%preswx(nobs))
    this%preswx = preswx

    allocate(this%pastwx(nobs))
    this%pastwx = pastwx

    !this%pastwx_duration = set_pastwx_duration(date10)
    allocate(this%pastwx_durations(nobs))
    call this%set_pastwx_durations()

    ! Make sure lats/lons are valid.
    call check_latlons(this)

    ! Handle arrays for unique networks in private method.  Also set &
    ! cdms_flag.
    call set_unique_networks(this)

  end subroutine USAF_gages_new

  ! Destructor
  subroutine USAF_gages_delete(this)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this

    this%date10 = "NULL"

    this%nobs = 0
    if (allocated(this%YYYYMMDDhhmmss)) deallocate(this%YYYYMMDDhhmmss)
    if (allocated(this%networks)) deallocate(this%networks)
    if (allocated(this%platforms)) deallocate(this%platforms)
    if (allocated(this%wmocode_id)) deallocate(this%wmocode_id)
    if (allocated(this%fipscode_id)) deallocate(this%fipscode_id)
    if (allocated(this%wmonumbers)) deallocate(this%wmonumbers)
    if (allocated(this%bsn)) deallocate(this%bsn)
    if (allocated(this%lats)) deallocate(this%lats)
    if (allocated(this%lons)) deallocate(this%lons)
    if (allocated(this%amts24)) deallocate(this%amts24)
    if (allocated(this%amts21)) deallocate(this%amts21)
    if (allocated(this%amts18)) deallocate(this%amts18)
    if (allocated(this%amts15)) deallocate(this%amts15)
    if (allocated(this%amts12)) deallocate(this%amts12)
    if (allocated(this%amts12_orig)) deallocate(this%amts12_orig)
    if (allocated(this%amts09)) deallocate(this%amts09)
    if (allocated(this%amts06)) deallocate(this%amts06)
    if (allocated(this%amts06_orig)) deallocate(this%amts06_orig)
    if (allocated(this%amts03)) deallocate(this%amts03)
    if (allocated(this%amts02)) deallocate(this%amts02)
    if (allocated(this%amts01)) deallocate(this%amts01)
    if (allocated(this%amts00)) deallocate(this%amts00)
    if (allocated(this%durations)) deallocate(this%durations)
    if (allocated(this%preswx)) deallocate(this%preswx)
    if (allocated(this%pastwx)) deallocate(this%pastwx)
    if (allocated(this%unique_networks)) deallocate(this%unique_networks)
    if (allocated(this%firsts)) deallocate(this%firsts)
    if (allocated(this%lasts)) deallocate(this%lasts)
    if (allocated(this%pastwx_durations)) &
         deallocate(this%pastwx_durations)
    this%num_unique_networks = 0
  end subroutine USAF_gages_delete

  ! Use the miscellaneous precip amount.
  subroutine USAF_gages_use_misc_precip(this)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this

    ! Locals
    integer :: russia_lowlimit
    integer :: russia_highlimit
    integer :: india_lowlimit
    integer :: india_highlimit
    integer :: srilanka_lowlimit
    integer :: srilanka_highlimit
    integer :: wmo_block
    character(10) :: date10
    integer :: nobs
    integer :: i

    if (this%cdms_flag) then
       india_lowlimit = CDMS_INDIA_LOWLIMIT
       india_highlimit = CDMS_INDIA_HIGHLIMIT
       russia_lowlimit = CDMS_RUSSIA_LOWLIMIT
       russia_highlimit = CDMS_RUSSIA_HIGHLIMIT
       srilanka_lowlimit = CDMS_SRILANKA_LOWLIMIT
       srilanka_highlimit = CDMS_SRILANKA_HIGHLIMIT
    else
       india_lowlimit = JMOBS_INDIA_LOWLIMIT
       india_highlimit = JMOBS_INDIA_HIGHLIMIT
       russia_lowlimit = JMOBS_RUSSIA_LOWLIMIT
       russia_highlimit = JMOBS_RUSSIA_HIGHLIMIT
       srilanka_lowlimit = JMOBS_SRILANKA_LOWLIMIT
       srilanka_highlimit = JMOBS_SRILANKA_HIGHLIMIT
    end if

    date10 = this%date10

    nobs = this%nobs
    do i = 1, nobs

       ! WMO block number in preobs files usually range from 1 to 100.
       ! Multiply so they are directly comparable to the
       ! appropriate WMO block ID ranges.
       wmo_block = this%wmonumbers(i)
       if (wmo_block < 1000) then
          if (this%cdms_flag) then
             wmo_block = wmo_block*10000
          else
             wmo_block = wmo_block*1000
          end if
       end if

       ! India rule...Obs are usually relative to 03Z with an invalid
       ! duration flag. So we may need to copy values from miscellaneous
       ! (amts00) to a more standard duration.
       if (this%durations(i) .eq. 0 .or. &
            this%durations(i) .eq. MISSING) then

          if ((wmo_block .ge. india_lowlimit .and. &
               wmo_block .le. india_highlimit) .or. &
               (this%bsn(i) .ge. india_lowlimit .and. &
                this%bsn(i) .le. india_highlimit)) then

             if (this%amts00(i) .ne. MISSING) then

                if (date10(9:10) == "03" .and. &
                     this%amts24(i) .eq. MISSING) then
                   this%amts24(i) = this%amts00(i)
                else if (date10(9:10) == "06" .and. &
                     this%amts03(i) .eq. MISSING) then
                   this%amts03(i) = this%amts00(i)
                else if (date10(9:10) == "09" .and. &
                     this%amts06(i) .eq. MISSING) then
                   this%amts06(i) = this%amts00(i)
                else if (date10(9:10) == "12" .and. &
                     this%amts09(i) .eq. MISSING) then
                   this%amts09(i) = this%amts00(i)
                else if (date10(9:10) == "15" .and. &
                     this%amts12(i) .eq. MISSING) then
                   this%amts12(i) = this%amts00(i)
                else if (date10(9:10) == "18" .and. &
                     this%amts15(i) .eq. MISSING) then
                   this%amts15(i) = this%amts00(i)
                else if (date10(9:10) == "21" .and. &
                     this%amts18(i) .eq. MISSING) then
                   this%amts18(i) = this%amts00(i)
                else if (date10(9:10) == "00" .and. &
                     this%amts21(i) .eq. MISSING) then
                   this%amts21(i) = this%amts00(i)
                end if
             end if
          end if
       end if

       ! Sri Lanka rule.
       if ((wmo_block .ge. srilanka_lowlimit .and. &
            wmo_block .le. srilanka_highlimit) .or. &
            (this%bsn(i) .ge. srilanka_lowlimit .and. &
            this%bsn(i) .le. srilanka_highlimit) ) then

          ! At present, it doesn't appear a special Sri Lanka rule is
          ! required.  We can process 3, 9, and 15 hour reports.
          continue
       end if

       ! Russia rule (actually, post-Soviet rule).  In the past, USAF
       ! personnel noted many "Russian" obs were only reported at 00Z and
       ! 12Z, and did not use a valid duration flag. These are decoded in
       ! the miscellaneous accumulation (amts00), and the duration is set
       ! to 0.  In AGRMET they were treated as 12-hr accumulations.  More
       ! recently, a spot check in July 2021 showed similar invalid
       ! durations in Kazakhstan, Kyrgyzstan, Georgia, Tajiskistan, and
       ! Turkenistan at 06Z and 18Z. Since WMO Region VI (Europe) expects
       ! 12-hr accumulations at 06Z and 18Z, we assume these rogue
       ! reports are 12-hr accumulations.
       if (this%durations(i) .eq. 0 .or. &
            this%durations(i) .eq. MISSING) then

          if ((wmo_block .ge. russia_lowlimit .and. &
               wmo_block .le. russia_highlimit) .or. &
               (this%bsn(i) .ge. russia_lowlimit .and. &
                this%bsn(i) .le. russia_highlimit) ) then

             if (date10(9:10) .eq. "00" .or. date10(9:10) .eq. "06" .or. &
                  date10(9:10) .eq. "12" .or. date10(9:10) .eq. "18") then

                if (this%amts00(i) .ne. MISSING) then
                   if (this%amts12(i) .eq. MISSING) then
                      this%amts12(i) = this%amts00(i)
                   end if
                end if
             end if
          end if
       end if ! Russia rule

       ! Copy from miscellaneous (amts00) to a more standard duration,
       ! according to the reported duration.
       if (this%amts00(i) .ne. MISSING) then
          if (this%amts24(i) .eq. MISSING .and. &
               this%durations(i) .eq. 24) then
             this%amts24(i) = this%amts00(i)
          else if (this%amts21(i) .eq. MISSING .and. &
               this%durations(i) .eq. 21) then
             this%amts21(i) = this%amts00(i)
          else if (this%amts18(i) .eq. MISSING .and. &
               this%durations(i) .eq. 18) then
             this%amts18(i) = this%amts00(i)
          else if (this%amts15(i) .eq. MISSING .and. &
               this%durations(i) .eq. 15) then
             this%amts15(i) = this%amts00(i)
          else if (this%amts12(i) .eq. MISSING .and. &
               this%durations(i) .eq. 12) then
             this%amts12(i) = this%amts00(i)
          else if (this%amts09(i) .eq. MISSING .and. &
               this%durations(i) .eq. 9) then
             this%amts09(i) = this%amts00(i)
          else if (this%amts06(i) .eq. MISSING .and. &
               this%durations(i) .eq. 6) then
             this%amts06(i) = this%amts00(i)
          else if (this%amts03(i) .eq. MISSING .and. &
               this%durations(i) .eq. 3) then
             this%amts03(i) = this%amts00(i)
          ! EMK...Disabled saving 2-hr and 1-hr amounts from misc.
          ! Old preobs files have a mix of different hours, so we
          ! can't expect this to work.  We will revisit with the
          ! new preobs files, likely requiring a new argument to
          ! this subroutine to activate.
          ! else if (this%amts02(i) .eq. MISSING .and. &
          !      this%durations(i) .eq. 2) then
          !    this%amts02(i) = this%amts00(i)
          ! else if (this%amts01(i) .eq. MISSING .and. &
          !      this%durations(i) .eq. 1) then
          !    this%amts01(i) = this%amts00(i)
          end if
       end if

    end do ! i

    ! Make sure shorter durations don't exceed larger durations.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_use_misc_precip

  ! Reconciles different accumulations from the same report.
  subroutine USAF_gages_reconcile_self(this)

    use LIS_logMod, only: LIS_logunit

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this

    ! Locals
    integer :: nobs
    integer :: i
    integer :: amts_tmp(10) ! 1, 2, 3, 6, 9, 12, 15, 18, 21, 24
    integer :: j, jj

    nobs = this%nobs
    do i = 1, nobs

       ! Zero out smaller-duration accumulation if larger one is zero.
       if (this%amts24(i) .eq. 0) then
          this%amts21(i) = 0
       end if
       if (this%amts21(i) .eq. 0) then
          this%amts18(i) = 0
       end if
       if (this%amts18(i) .eq. 0) then
          this%amts15(i) = 0
       end if
       if (this%amts15(i) .eq. 0) then
          this%amts12(i) = 0
       end if
       if (this%amts12(i) .eq. 0) then
          this%amts09(i) = 0
       end if
       if (this%amts09(i) .eq. 0) then
          this%amts06(i) = 0
       end if
       if (this%amts06(i) .eq. 0) then
          this%amts03(i) = 0
       end if
       if (this%amts03(i) .eq. 0) then
          this%amts02(i) = 0
       end if
       if (this%amts02(i) .eq. 0) then
          this%amts01(i) = 0
       end if

       ! Sometimes smaller-duration accumulations exceed larger-durations,
       ! due to precision differences in WMO reporting standard.  We
       ! correct that here.
       amts_tmp(1) = this%amts01(i)
       amts_tmp(2) = this%amts02(i)
       amts_tmp(3) = this%amts03(i)
       amts_tmp(4) = this%amts06(i)
       amts_tmp(5) = this%amts09(i)
       amts_tmp(6) = this%amts12(i)
       amts_tmp(7) = this%amts15(i)
       amts_tmp(8) = this%amts18(i)
       amts_tmp(9) = this%amts21(i)
       amts_tmp(10) = this%amts24(i)

       do j = 10, 2, -1 ! Longer duration limit, 24-hr to 02-hr
          if (amts_tmp(j) .ne. MISSING) then
             do jj = j-1, 1, -1 ! Shorter duration limit, 21-hr to 01-hr
                if (amts_tmp(jj) .ne. MISSING) then
                   amts_tmp(jj) = min(amts_tmp(jj), amts_tmp(j))
                end if
             end do ! jj
          end if
       end do ! j

       this%amts01(i) = amts_tmp(1)
       this%amts02(i) = amts_tmp(2)
       this%amts03(i) = amts_tmp(3)
       this%amts06(i) = amts_tmp(4)
       this%amts09(i) = amts_tmp(5)
       this%amts12(i) = amts_tmp(6)
       this%amts15(i) = amts_tmp(7)
       this%amts18(i) = amts_tmp(8)
       this%amts21(i) = amts_tmp(9)
       this%amts24(i) = amts_tmp(10)

    end do ! i

  end subroutine USAF_gages_reconcile_self

  ! Fill gaps in precip record if bookend accumulations are identical.
  subroutine USAF_gages_fill_gaps(this)

    use LIS_logmod, only:  LIS_logunit

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this

    ! Locals
    integer :: nobs
    integer :: i
    integer :: amts_tmp(10) ! 1, 2, 3, 6, 9, 12, 15, 18, 21, 24
    integer :: j, jj, jjj
    logical :: is_data_gap

    nobs = this%nobs
    do i = 1, nobs

       amts_tmp(1) = this%amts01(i)
       amts_tmp(2) = this%amts02(i)
       amts_tmp(3) = this%amts03(i)
       amts_tmp(4) = this%amts06(i)
       amts_tmp(5) = this%amts09(i)
       amts_tmp(6) = this%amts12(i)
       amts_tmp(7) = this%amts15(i)
       amts_tmp(8) = this%amts18(i)
       amts_tmp(9) = this%amts21(i)
       amts_tmp(10) = this%amts24(i)

       ! Fill in gaps if shorter and longer accumulations are identical.
       do j = 1, 8 ! Shorter duration limit, from 01-hr to 18-hr
          do jj = j+2, 10 ! Longer duration limit, from 03-hr to 24-hr
             if (amts_tmp(j) .ne. MISSING .and. &
                  amts_tmp(jj) .ne. MISSING .and. &
                  amts_tmp(j) == amts_tmp(jj)) then
                ! We have two bookends with identical values.  See
                ! if all the durations in-between are missing.
                is_data_gap = .true.
                do jjj = j+1, jj-1
                   if (amts_tmp(jjj) .ne. MISSING) then
                      is_data_gap = .false.
                      exit
                   end if
                end do ! jjj
                ! If all the durations in between are missing, just
                ! copy the value (we can assume no additional precip
                ! occurred during these periods, at least from this
                ! particular report).
                if (is_data_gap) then
                   do jjj = j+1, jj-1
                      amts_tmp(jjj) = amts_tmp(j)
                   end do ! jjj
                end if
             end if
          end do ! jj
       end do ! j

       this%amts01(i) = amts_tmp(1)
       this%amts02(i) = amts_tmp(2)
       this%amts03(i) = amts_tmp(3)
       this%amts06(i) = amts_tmp(4)
       this%amts09(i) = amts_tmp(5)
       this%amts12(i) = amts_tmp(6)
       this%amts15(i) = amts_tmp(7)
       this%amts18(i) = amts_tmp(8)
       this%amts21(i) = amts_tmp(9)
       this%amts24(i) = amts_tmp(10)

    end do ! i

  end subroutine USAF_gages_fill_gaps

  ! Reconcile current obs with obs from 1 hour ago
  subroutine USAF_gages_reconcile_gages01(this, gages01)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this
    class(USAF_gages_t), intent(in) :: gages01

    ! Locals
    integer :: nobs
    integer tmp
    integer :: i, i1

    nobs = this%nobs

    do i = 1, nobs

       ! Find current station in prior report list
       i1 = search_pcpobs(gages01, this%networks(i), this%platforms(i), &
            this%wmocode_id(i), this%fipscode_id(i))

       if (i1 .eq. MISSING) cycle

       ! Use current 1-hr accumulation, if available.
       if (this%amts01(i) .ne. MISSING) then

          ! Update 3-hr accumulation if missing
          if (this%amts03(i) .eq. MISSING .and. &
               gages01%amts02(i1) .ne. MISSING ) then
             tmp = this%amts01(i) + gages01%amts02(i1)
             if (tmp >= 0 .and. tmp <= MAX_PCP_3HR) this%amts03(i) = tmp
          end if

          ! Update 2-hr accumulation if missing
          if (this%amts02(i) .eq. MISSING .and. &
               gages01%amts01(i1) .ne. MISSING ) then
             tmp = this%amts01(i) + gages01%amts01(i1)
             if (tmp >= 0 .and. tmp <= MAX_PCP_2HR)  this%amts02(i) = tmp
          end if

       end if

       ! Update the 01-hr accumulation
       if (this%amts01(i) .eq. MISSING) then
          if (this%amts02(i) .ne. MISSING .and. &
               gages01%amts01(i1) .ne. MISSING) then
             tmp = this%amts02(i) - gages01%amts01(i1)
             if (tmp >= 0 .and. tmp <= MAX_PCP_1HR) then
                this%amts01(i) = tmp
                cycle
             end if
          end if
          if (this%amts03(i) .ne. MISSING .and. &
               gages01%amts02(i1) .ne. MISSING) then
             tmp = this%amts03(i) - gages01%amts02(i1)
             if (tmp >= 0 .and. tmp <= MAX_PCP_1HR) then
                this%amts01(i) = tmp
                cycle
             end if
          end if
       end if

    end do

    ! Reconcile current obs again to ensure smaller-duration accumulations
    ! don't exceed longer-duration.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_reconcile_gages01

  ! Reconcile current obs with obs from 2 hours ago
  subroutine USAF_gages_reconcile_gages02(this, gages02)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this
    class(USAF_gages_t), intent(in) :: gages02

    ! Locals
    integer :: nobs
    integer :: tmp
    integer :: i, i2

    nobs = this%nobs

    do i = 1, nobs

       ! Find current station in prior report list
       i2 = search_pcpobs(gages02, this%networks(i), this%platforms(i), &
            this%wmocode_id(i), this%fipscode_id(i))

       if (i2 .eq. MISSING) cycle

       ! Use current 2-hr accumulation, if available.
       if (this%amts02(i) .ne. MISSING) then

          ! Update 3-hr accumulation if missing
          if (this%amts03(i) .eq. MISSING .and. &
               gages02%amts01(i2) .ne. MISSING ) then
             tmp = this%amts02(i) + gages02%amts01(i2)
             if (tmp >= 0 .and. tmp <= MAX_PCP_3HR) this%amts03(i) = tmp
          end if
       end if

       ! Update the 2-hr accumulation
       if (this%amts02(i) .eq. MISSING) then
          if (this%amts03(i) .ne. MISSING .and. &
               gages02%amts01(i2) .ne. MISSING) then
             tmp = this%amts03(i) - gages02%amts01(i2)
             if (tmp >= 0 .and. tmp <= MAX_PCP_2HR) this%amts02(i) = tmp
          end if
       end if

    end do

    ! Reconcile current obs again to ensure smaller-duration accumulations
    ! don't exceed longer-duration.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_reconcile_gages02

  ! Reconcile current obs with obs from 3 hours ago
  subroutine USAF_gages_reconcile_gages03(this, gages03)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this
    class(USAF_gages_t), intent(in) :: gages03

    ! Locals
    integer :: nobs
    integer :: tmp
    integer :: i, i3

    nobs = this%nobs
    do i = 1, nobs

       ! Find current station in prior report list
       i3 = search_pcpobs(gages03, this%networks(i), this%platforms(i), &
            this%wmocode_id(i), this%fipscode_id(i))

       if (i3 .eq. MISSING) cycle

       ! Use current 3-hr accumulation, if available.
       if (this%amts03(i) .ne. MISSING) then

          ! Update 24-hr accumulation if missing
          if (this%amts24(i) .eq. MISSING .and. &
               gages03%amts21(i3) .ne. MISSING ) then
             tmp = this%amts03(i) + gages03%amts21(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_24HR) this%amts24(i) = tmp
          end if

          ! Update 21-hr accumulation if missing
          if (this%amts21(i) .eq. MISSING .and. &
               gages03%amts18(i3) .ne. MISSING) then
             tmp = this%amts03(i) + gages03%amts18(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_21HR) this%amts21(i) = tmp
          end if

          ! Update 18-hr accumulation if missing
          if (this%amts18(i) .eq. MISSING .and. &
               gages03%amts15(i3) .ne. MISSING) then
             tmp = this%amts03(i) + gages03%amts15(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_18HR) this%amts18(i) = tmp
          end if

          ! Update 15-hr accumulation if missing
          if (this%amts15(i) .eq. MISSING .and. &
               gages03%amts12(i3) .ne. MISSING) then
             tmp = this%amts03(i) + gages03%amts12(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_15HR) this%amts15(i) = tmp
          end if

          ! Update 12-hr accumulation if missing
          if (this%amts12(i) .eq. MISSING .and. &
               gages03%amts09(i3) .ne. MISSING) then
             tmp = this%amts03(i) + gages03%amts09(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_12HR) this%amts12(i) = tmp
          end if

          ! Update 09-hr accumulation if missing
          if (this%amts09(i) .eq. MISSING .and. &
               gages03%amts06(i3) .ne. MISSING) then
             tmp = this%amts03(i) + gages03%amts06(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_9HR) this%amts09(i) = tmp
          end if

          ! Update 06-hr accumulation if missing
          if (this%amts06(i) .eq. MISSING .and. &
               gages03%amts03(i3) .ne. MISSING) then
             tmp = this%amts03(i) + gages03%amts03(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_6HR) this%amts06(i) = tmp
          end if
       end if

       ! Multiple attempts to update 03-hr accumulation.
       if (this%amts03(i) .eq. MISSING) then
          if (this%amts06(i) .ne. MISSING .and. &
               gages03%amts03(i3) .ne. MISSING) then
             tmp = this%amts06(i) - gages03%amts03(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_3HR) then
                this%amts03(i) = tmp
                cycle
             end if
          end if

          if (this%amts09(i) .ne. MISSING .and. &
               gages03%amts06(i3) .ne. MISSING) then
             tmp = this%amts09(i) - gages03%amts06(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_3HR) then
                this%amts03(i) = tmp
                cycle
             end if
          end if

          if (this%amts12(i) .ne. MISSING .and. &
               gages03%amts09(i3) .ne. MISSING) then
             tmp = this%amts12(i) - gages03%amts09(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_3HR) then
                this%amts03(i) = tmp
                cycle
             end if
          end if

          if (this%amts15(i) .ne. MISSING .and. &
               gages03%amts12(i3) .ne. MISSING) then
             tmp = this%amts15(i) - gages03%amts12(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_3HR) then
                this%amts03(i) = tmp
                cycle
             end if
          end if

          if (this%amts18(i) .ne. MISSING .and. &
               gages03%amts15(i3) .ne. MISSING) then
             tmp = this%amts18(i) - gages03%amts15(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_3HR) then
                this%amts03(i) = tmp
                cycle
             end if
          end if

          if (this%amts21(i) .ne. MISSING .and. &
               gages03%amts18(i3) .ne. MISSING) then
             tmp = this%amts21(i) - gages03%amts18(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_3HR) then
                this%amts03(i) = tmp
                cycle
             end if
          end if

          if (this%amts24(i) .ne. MISSING .and. &
               gages03%amts21(i3) .ne. MISSING ) then
             tmp = this%amts24(i) - gages03%amts21(i3)
             if (tmp >= 0 .and. tmp <= MAX_PCP_3HR) then
                this%amts03(i) = tmp
                cycle
             end if
          end if
       end if

    end do

    ! Reconcile current obs again to ensure smaller-duration accumulations
    ! don't exceed longer-duration.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_reconcile_gages03

  ! Reconcile current obs with obs from 6 hours ago.
  subroutine USAF_gages_reconcile_gages06(this, gages06)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this
    class(USAF_gages_t), intent(in) :: gages06

    ! Locals
    integer :: nobs
    integer :: tmp
    integer :: i, i6

    nobs = this%nobs

    do i = 1, nobs

       ! Find current station in prior report list
       i6 = search_pcpobs(gages06, this%networks(i), this%platforms(i), &
            this%wmocode_id(i), this%fipscode_id(i))

       if (i6 .eq. MISSING) cycle

       ! Use current 6-hr accumulation, if available.
       if (this%amts06(i) .ne. MISSING) then

          ! Update 24-hr accumulation if missing
          if (this%amts24(i) .eq. MISSING .and. &
               gages06%amts18(i6) .ne. MISSING) then
             tmp = this%amts06(i) + gages06%amts18(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_24HR) this%amts24(i) = tmp
          end if

          ! Update 21-hr accumulation if missing
          if (this%amts21(i) .eq. MISSING .and. &
               gages06%amts15(i6) .ne. MISSING) then
             tmp = this%amts06(i) + gages06%amts15(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_21HR) this%amts21(i) = tmp
          end if

          ! Update 18-hr accumulation if missing
          if (this%amts18(i) .eq. MISSING .and. &
               gages06%amts12(i6) .ne. MISSING) then
             tmp = this%amts06(i) + gages06%amts12(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_18HR) this%amts18(i) = tmp
          end if

          ! Update 15-hr accumulation if missing
          if (this%amts15(i) .eq. MISSING .and. &
               gages06%amts09(i6) .ne. MISSING) then
             tmp = this%amts06(i) + gages06%amts09(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_15HR) this%amts15(i) = tmp
          end if

          ! Update 12-hr accumulation if missing
          if (this%amts12(i) .eq. MISSING .and. &
               gages06%amts06(i6) .ne. MISSING) then
             tmp = this%amts06(i) + gages06%amts06(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_12HR) this%amts12(i) = tmp
          end if

          ! Update 09-hr accumulation if missing
          if (this%amts09(i) .eq. MISSING .and. &
               gages06%amts03(i6) .ne. MISSING) then
             tmp = this%amts06(i) + gages06%amts03(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_9HR) this%amts09(i) = tmp
          end if

       end if

       ! Multiple attempts to update 06-hr accumulation.
       if (this%amts06(i) .eq. MISSING) then
          if ( (this%amts12(i) .ne. MISSING) .and. &
               (gages06%amts06(i6) .ne. MISSING) ) then
             tmp = this%amts12(i) - gages06%amts06(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_6HR) then
                this%amts06(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts15(i) .ne. MISSING) .and. &
               (gages06%amts09(i6) .ne. MISSING) ) then
             tmp = this%amts15(i) - gages06%amts09(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_6HR) then
                this%amts06(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts18(i) .ne. MISSING) .and. &
               (gages06%amts12(i6) .ne. MISSING) ) then
             tmp = this%amts18(i) - gages06%amts12(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_6HR) then
                this%amts06(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts21(i) .ne. MISSING) .and. &
               (gages06%amts15(i6) .ne. MISSING) ) then
             tmp = this%amts21(i) - gages06%amts15(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_6HR) then
                this%amts06(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts24(i) .ne. MISSING) .and. &
               (gages06%amts18(i6) .ne. MISSING) ) then
             tmp = this%amts24(i) - gages06%amts18(i6)
             if (tmp >= 0 .and. tmp <= MAX_PCP_6HR) then
                this%amts06(i) = tmp
                cycle
             end if
          end if
       end if

    end do

    ! Reconcile current obs again to ensure smaller-duration
    ! accumulations don't exceed longer-duration.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_reconcile_gages06

  ! Reconcile current obs with obs from 9 hours ago.
  subroutine USAF_gages_reconcile_gages09(this, gages09)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this
    class(USAF_gages_t), intent(in) :: gages09

    ! Locals
    integer :: nobs
    integer :: tmp
    integer :: i, i9

    nobs = this%nobs

    do i = 1, nobs

       ! Find current station in prior report list
       i9 = search_pcpobs(gages09, this%networks(i), this%platforms(i), &
            this%wmocode_id(i), this%fipscode_id(i))

       if (i9 .eq. MISSING) cycle

       ! Use current 9-hr accumulation, if available.
       if (this%amts09(i) .ne. MISSING) then

          ! Update 24-hr accumulation if missing
          if (this%amts24(i) .eq. MISSING .and. &
               gages09%amts15(i9) .ne. MISSING ) then
             tmp = this%amts09(i) + gages09%amts15(i9)
             if (tmp >= 0 .and. tmp <= MAX_PCP_24HR) this%amts24(i) = tmp
          end if

          ! Update 21-hr accumulation if missing
          if (this%amts21(i) .eq. MISSING .and. &
               gages09%amts12(i9) .ne. MISSING ) then
             tmp = this%amts09(i) + gages09%amts12(i9)
             if (tmp >= 0 .and. tmp <= MAX_PCP_21HR) this%amts21(i) = tmp
          end if

          ! Update 18-hr accumulation if missing
          if (this%amts18(i) .eq. MISSING .and. &
               gages09%amts09(i9) .ne. MISSING) then
             tmp = this%amts09(i) + gages09%amts09(i9)
             if (tmp >= 0 .and. tmp <= MAX_PCP_18HR) this%amts18(i) = tmp
          end if

          ! Update 15-hr accumulation if missing
          if (this%amts15(i) .eq. MISSING .and. &
               gages09%amts06(i9) .ne. MISSING) then
             tmp = this%amts09(i) + gages09%amts06(i9)
             if (tmp >= 0 .and. tmp <= MAX_PCP_15HR) this%amts15(i) = tmp
          end if

          ! Update 12-hr accumulation if missing
          if (this%amts12(i) .eq. MISSING .and. &
               gages09%amts03(i9) .ne. MISSING) then
             tmp = this%amts09(i) + gages09%amts03(i9)
             if (tmp >= 0 .and. tmp <= MAX_PCP_12HR) this%amts12(i) = tmp
          end if

       end if

       ! Multiple attempts to update 09-hr accumulation.
       if (this%amts09(i) .eq. MISSING) then
          if ( (this%amts12(i) .ne. MISSING) .and. &
               (gages09%amts03(i9) .ne. MISSING) ) then
             tmp = this%amts12(i) - gages09%amts03(i9)
             if (tmp >= 0 .and. tmp <= MAX_PCP_9HR) then
                this%amts09(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts15(i) .ne. MISSING) .and. &
               (gages09%amts06(i9) .ne. MISSING) ) then
             this%amts09(i) = this%amts15(i) - gages09%amts06(i9)
             this%amts09(i) = max(this%amts09(i), 0)
             tmp = this%amts15(i) - gages09%amts06(i9)
             if (tmp >= 0 .and. tmp <= MAX_PCP_9HR) then
                this%amts09(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts18(i) .ne. MISSING) .and. &
               (gages09%amts09(i9) .ne. MISSING) ) then
             tmp = this%amts18(i) - gages09%amts09(i9)
             if (tmp >= 0 .and. tmp <= MAX_PCP_9HR) then
                this%amts09(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts21(i) .ne. MISSING) .and. &
               (gages09%amts12(i9) .ne. MISSING) ) then
             tmp = this%amts21(i) - gages09%amts12(i9)
             if (tmp >= 0 .and. tmp <= MAX_PCP_9HR) then
                this%amts09(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts24(i) .ne. MISSING) .and. &
               (gages09%amts15(i9) .ne. MISSING) ) then
             tmp = this%amts24(i) - gages09%amts15(i9)
             if (tmp >= 0 .and. tmp <= MAX_PCP_9HR) then
                this%amts09(i) = tmp
                cycle
             end if
          end if
       end if
    end do

    ! Reconcile current obs again to ensure smaller-duration
    ! accumulations don't exceed longer-duration.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_reconcile_gages09

  ! Reconcile current obs with obs from 12 hours ago.
  subroutine USAF_gages_reconcile_gages12(this, gages12)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this
    class(USAF_gages_t), intent(in) :: gages12

    ! Locals
    integer :: nobs
    integer :: tmp
    integer :: i, i12

    nobs = this%nobs

    do i = 1, nobs

       ! Find current station in prior report list
       i12 = search_pcpobs(gages12, this%networks(i), &
            this%platforms(i), this%wmocode_id(i), this%fipscode_id(i))

       if (i12 .eq. MISSING) cycle

       ! Use current 12-hr accumulation, if available
       if (this%amts12(i) .ne. MISSING) then

          ! Update 24-hr accumulation if missing
          if (this%amts24(i) .eq. MISSING .and. &
               gages12%amts12(i12) .ne. MISSING) then
             tmp = this%amts12(i) + gages12%amts12(i12)
             if (tmp >= 0 .and. tmp <= MAX_PCP_24HR) this%amts24(i) = tmp
          end if

          ! Update 21-hr accumulation if missing
          if (this%amts21(i) .eq. MISSING .and. &
               gages12%amts09(i12) .ne. MISSING) then
             tmp = this%amts12(i) + gages12%amts09(i12)
             if (tmp >= 0 .and. tmp <= MAX_PCP_21HR) this%amts21(i) = tmp
          end if

          ! Update 18-hr accumulation if missing
          if (this%amts18(i) .eq. MISSING .and. &
               gages12%amts06(i12) .ne. MISSING) then
             tmp = this%amts12(i) + gages12%amts06(i12)
             if (tmp >= 0 .and. tmp <= MAX_PCP_18HR) this%amts18(i) = tmp
          end if

          ! Update 15-hr accumulation if missing
          if (this%amts15(i) .eq. MISSING .and. &
               gages12%amts03(i12) .ne. MISSING) then
             tmp = this%amts12(i) + gages12%amts03(i12)
             if (tmp >= 0 .and. tmp <= MAX_PCP_15HR) this%amts15(i) = tmp
          end if
       end if

       ! Multiple attempts to update 12-hr accumulation.
       if (this%amts12(i) .eq. MISSING) then
          if ( (this%amts15(i) .ne. MISSING) .and. &
               (gages12%amts03(i12) .ne. MISSING) ) then
             tmp = this%amts15(i) - gages12%amts03(i12)
             if (tmp >= 0 .and. tmp <= MAX_PCP_12HR) then
                this%amts12(i) = tmp
                cycle
             end if
          end if

          if( (this%amts18(i) .ne. MISSING) .and. &
               (gages12%amts06(i12) .ne. MISSING) ) then
             tmp = this%amts18(i) - gages12%amts06(i12)
             if (tmp >= 0 .and. tmp <= MAX_PCP_12HR) then
                this%amts12(i) = tmp
                cycle
             end if
          end if

          if( (this%amts21(i) .ne. MISSING) .and. &
               (gages12%amts09(i12) .ne. MISSING) ) then
             tmp = this%amts21(i) - gages12%amts09(i12)
             if (tmp >= 0 .and. tmp <= MAX_PCP_12HR) then
                this%amts12(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts24(i) .ne. MISSING) .and. &
               (gages12%amts12(i12) .ne. MISSING) ) then
             tmp = this%amts24(i) - gages12%amts12(i12)
             if (tmp >= 0 .and. tmp <= MAX_PCP_12HR) then
                this%amts12(i) = tmp
                cycle
             end if
          end if

       end if

    end do

    ! Reconcile current obs again to ensure smaller-duration
    ! accumulations don't exceed longer-duration.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_reconcile_gages12

  ! Reconcile current obs with obs from 15 hours ago.
  subroutine USAF_gages_reconcile_gages15(this, gages15)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this
    class(USAF_gages_t), intent(in) :: gages15

    ! Locals
    integer :: nobs
    integer :: tmp
    integer :: i, i15

    nobs = this%nobs

    do i = 1, nobs

       ! Find current station in prior report list
       i15 = search_pcpobs(gages15, this%networks(i), &
            this%platforms(i), this%wmocode_id(i), this%fipscode_id(i))

       if (i15 .eq. MISSING) cycle

       ! Use current 15-hr accumulation, if available
       if (this%amts15(i) .ne. MISSING) then

          ! Update 24-hr accumulation if missing
          if (this%amts24(i) .eq. MISSING .and. &
               gages15%amts09(i15) .ne. MISSING) then
             tmp = this%amts15(i) + gages15%amts09(i15)
             if (tmp >= 0 .and. tmp <= MAX_PCP_24HR) this%amts24(i) = tmp
          end if

          ! Update 21-hr accumulation if missing
          if (this%amts21(i) .eq. MISSING .and. &
               gages15%amts06(i15) .ne. MISSING) then
             tmp = this%amts15(i) + gages15%amts06(i15)
             if (tmp >= 0 .and. tmp <= MAX_PCP_21HR) this%amts21(i) = tmp
          end if

          ! Update 18-hr accumulation if missing
          if (this%amts18(i) .eq. MISSING .and. &
               gages15%amts03(i15) .ne. MISSING) then
             tmp = this%amts15(i) + gages15%amts03(i15)
             if (tmp >= 0 .and. tmp <= MAX_PCP_18HR) this%amts18(i) = tmp
          end if

       end if

       ! Multiple attempts to update 15-hr accumulation.
       if (this%amts15(i) .eq. MISSING) then
          if ( (this%amts18(i) .ne. MISSING) .and. &
               (gages15%amts03(i15) .ne. MISSING) ) then
             tmp = this%amts18(i) - gages15%amts03(i15)
             if (tmp >= 0 .and. tmp <= MAX_PCP_15HR) then
                this%amts15(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts21(i) .ne. MISSING) .and. &
               (gages15%amts06(i15) .ne. MISSING) ) then
             tmp = this%amts21(i) - gages15%amts06(i15)
             if (tmp >= 0 .and. tmp <= MAX_PCP_15HR) then
                this%amts15(i) = tmp
                cycle
             end if
          end if

          if ( (this%amts24(i) .ne. MISSING) .and. &
               (gages15%amts09(i15) .ne. MISSING) ) then
             tmp = this%amts24(i) - gages15%amts09(i15)
             if (tmp >= 0 .and. tmp <= MAX_PCP_15HR) then
                this%amts15(i) = tmp
                cycle
             end if
          end if

       end if

    end do

    ! Reconcile current obs again to ensure smaller-duration &
    ! accumulations don't exceed longer-duration.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_reconcile_gages15

  ! Reconcile current obs with obs from 18 hours ago.
  subroutine USAF_gages_reconcile_gages18(this, gages18)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this
    class(USAF_gages_t), intent(in) :: gages18

    ! Locals
    integer :: nobs
    integer :: tmp
    integer :: i, i18

    nobs = this%nobs

    do i = 1, nobs

       ! Find current station in prior report list
       i18 = search_pcpobs(gages18, this%networks(i), &
            this%platforms(i), this%wmocode_id(i), this%fipscode_id(i))

       if (i18 .eq. MISSING) cycle

       ! Use current 18-hr accumulation, if available
       if (this%amts18(i) .ne. MISSING) then

          ! Update 24-hr accumulation if missing
          if (this%amts24(i) .eq. MISSING .and. &
               gages18%amts06(i18) .ne. MISSING) then
             tmp = this%amts18(i) + gages18%amts06(i18)
             if (tmp >= 0 .and. tmp <= MAX_PCP_24HR) this%amts24(i) = tmp
          end if

          ! Update 21-hr accumulation if missing
          if (this%amts21(i) .eq. MISSING .and. &
               gages18%amts03(i18) .ne. MISSING) then
             tmp = this%amts18(i) + gages18%amts03(i18)
             if (tmp >= 0 .and. tmp <= MAX_PCP_21HR) this%amts21(i) = tmp
          end if

       end if

       ! Multiple attempts to update 18-hr accumulation.
       if (this%amts18(i) .eq. MISSING) then
          if (this%amts21(i) .ne. MISSING .and. &
               gages18%amts03(i18) .ne. MISSING) then
             tmp = this%amts21(i) - gages18%amts03(i18)
             if (tmp >= 0 .and. tmp <= MAX_PCP_18HR) then
                this%amts18(i) = tmp
                cycle
             end if
          end if

          if (this%amts24(i) .ne. MISSING .and. &
               gages18%amts06(i18) .ne. MISSING) then
             tmp = this%amts24(i) - gages18%amts06(i18)
             if (tmp >= 0 .and. tmp <= MAX_PCP_18HR) then
                this%amts18(i) = tmp
                cycle
             end if
          end if
       end if

    end do

    ! Reconcile current obs again to ensure smaller-duration
    ! accumulations don't exceed longer-duration.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_reconcile_gages18

  ! Reconcile current obs with obs from 21 hours ago.
  subroutine USAF_gages_reconcile_gages21(this, gages21)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this
    class(USAF_gages_t), intent(in) :: gages21

    ! Locals
    integer :: nobs
    integer :: tmp
    integer :: i, i21

    nobs = this%nobs

    do i = 1, nobs

       ! Find current station in prior report list
       i21 = search_pcpobs(gages21, this%networks(i), &
            this%platforms(i), this%wmocode_id(i), this%fipscode_id(i))

       if (i21 .eq. MISSING) cycle

       ! Use current 21-hr accumulation, if available
       if (this%amts21(i) .ne. MISSING) then

          ! Update 24-hr accumulation if missing
          if (this%amts24(i) .eq. MISSING .and. &
               gages21%amts03(i21) .ne. MISSING) then
             tmp = this%amts21(i) + gages21%amts03(i21)
             if (tmp >= 0 .and. tmp <= MAX_PCP_24HR) this%amts24(i) = tmp
          end if

       end if

       ! Update 21-hr accumulation
       if (this%amts21(i) .eq. MISSING) then
          if (this%amts24(i) .ne. MISSING .and. &
               gages21%amts03(i21) .ne. MISSING) then
             tmp = this%amts24(i) - gages21%amts03(i21)
             if (tmp >= 0 .and. tmp <= MAX_PCP_21HR) this%amts21(i) = tmp
          end if
       end if

    end do

    ! Reconcile current obs again to ensure smaller-duration
    ! accumulations don't exceed longer-duration.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_reconcile_gages21

  ! Correct for missing "zero-precip" 12-hr reports at 12Z for Region III
  ! (South America).  Based on logic in AGRMET_processobs.
  subroutine USAF_gages_correct_region3_12Z(this)

    use LIS_logmod, only:  LIS_logunit
    
    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this

    ! Locals
    integer :: s_amer_lowlimit
    integer :: s_amer_highlimit
    integer :: wmo_block
    character(10) :: date10
    integer :: nobs
    integer :: i

    ! Only work with 12Z observations.
    date10 = this%date10
    if (date10(9:10) .ne. "12") return

    if (this%cdms_flag) then
       s_amer_lowlimit = CDMS_S_AMER_LOWLIMIT
       s_amer_highlimit = CDMS_S_AMER_HIGHLIMIT
    else
       s_amer_lowlimit = JMOBS_S_AMER_LOWLIMIT
       s_amer_highlimit = JMOBS_S_AMER_HIGHLIMIT
    end if

    nobs = this%nobs

    do i = 1, nobs

       ! WMO block number in preobs files usually range from 1 to 100.
       ! Multiply so they are directly comparable to the
       ! appropriate WMO block ID ranges.
       wmo_block = this%wmonumbers(i)
       if (wmo_block < 1000) then
          if (this%cdms_flag) then
             wmo_block = wmo_block*10000
          else
             wmo_block = wmo_block*1000
          end if
       end if

       ! In the past, USAF personnel noticed a lack of "zero-precip"
       ! 12-hr accumulations at 12Z over South America, which they
       ! attributed to two causes:
       ! (1) Many stations fail to report at 06Z; and
       ! (2) Many 12Z reports indicate zero precip without indicating
       !     duration.  (Speculation: SYNOP reports excluded a 6RRRt_r
       !     precip group in the ob, and instead mark the "I_R"
       !     indicator in Section 0 as "3" (meaning precipitation amount
       !     is zero) without indicating the duration.)
       ! Apparently this is decoded by default as "no precip over last
       ! six hours".
       ! To reconcile, South American obs at 12Z with missing 12-hr
       ! accumulations but with "zero" 6-hr accumulations will also
       ! be treated as "zero" 12-hr accumulations.
       ! EXCEPTION:  If a rare 9-hour non-zero report exists, don't
       ! change the 12-hr.
       if (this%amts12(i) .ne. MISSING) cycle
       if (this%amts12_orig(i) .ne. MISSING) cycle
       if (this%amts09(i) .gt. 0) cycle

       if ((wmo_block .ge. s_amer_lowlimit .and. &
            wmo_block .le. s_amer_highlimit) .or. &
            (this%bsn(i) .ge. s_amer_lowlimit .and. &
            this%bsn(i) .le. s_amer_highlimit)) then

          ! We want to make sure the 6-hr report is "original", and
          ! not just filled in or derived from comparing other gage
          ! reports or present/past weather.  That is, the station
          ! really made a "zero precip" ob at 12Z, and nothing is
          ! available from 06Z.
          if (this%amts06(i) .eq. 0 .and. &
               this%amts06_orig(i) .eq. 0) then
             this%amts12(i) = 0
             this%amts09(i) = 0
             this%amts06(i) = 0
             this%amts03(i) = 0
             this%amts02(i) = 0
             this%amts01(i) = 0
          end if
       end if

    end do

    ! Reconcile current obs again to ensure smaller-duration
    ! accumulations don't exceed longer-duration.
    call this%reconcile_self()
    call this%fill_gaps()

  end subroutine USAF_gages_correct_region3_12Z

  ! Write gage data to output file.
  subroutine USAF_gages_write_data(this, filename)

    ! Imports
    use LIS_logMod, only: LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_logunit

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(in) :: this
    character(*), intent(in) :: filename

    ! Locals
    integer :: istat
    integer :: nobs, nobs_good
    integer :: i
    integer :: iunit

    nobs = this%nobs

    ! First, count number of good obs in the structure, after quality
    ! control.
    nobs_good = 0
    do i = 1, nobs
       if (this%amts24(i) .eq. MISSING .and. &
            this%amts21(i) .eq. MISSING .and. &
            this%amts18(i) .eq. MISSING .and. &
            this%amts15(i) .eq. MISSING .and. &
            this%amts12(i) .eq. MISSING .and. &
            this%amts09(i) .eq. MISSING .and. &
            this%amts06(i) .eq. MISSING .and. &
            this%amts03(i) .eq. MISSING .and. &
            this%amts02(i) .eq. MISSING .and. &
            this%amts01(i) .eq. MISSING) cycle
       nobs_good = nobs_good + 1
    end do

    iunit = LIS_getNextUnitNumber()
    open(iunit, file=trim(filename), iostat=istat, err=300)
    write(iunit, *, iostat=istat, err=300) nobs_good
    do i = 1, nobs

       if (this%amts24(i) .eq. MISSING .and. &
            this%amts21(i) .eq. MISSING .and. &
            this%amts18(i) .eq. MISSING .and. &
            this%amts15(i) .eq. MISSING .and. &
            this%amts12(i) .eq. MISSING .and. &
            this%amts09(i) .eq. MISSING .and. &
            this%amts06(i) .eq. MISSING .and. &
            this%amts03(i) .eq. MISSING .and. &
            this%amts02(i) .eq. MISSING .and. &
            this%amts01(i) .eq. MISSING) cycle
       write(iunit, 6000, iostat=istat, err=300) &
            this%YYYYMMDDhhmmss(i), &
            this%networks(i), this%platforms(i), &
            this%wmocode_id(i), this%fipscode_id(i), &
            this%wmonumbers(i), this%bsn(i), &
            this%lats(i), this%lons(i), &
            this%amts24(i), this%amts21(i), this%amts18(i), &
            this%amts15(i), this%amts12(i), this%amts09(i), &
            this%amts06(i), this%amts03(i), &
            this%amts02(i), this%amts01(i), &
            this%amts00(i), this%durations(i), &
            this%preswx(i), this%pastwx(i)
6000   format (a14, 1x, &
            a32, 1x, a32, 1x, &
            a2, 1x, a2, 1x, &
            i6, 1x, i6, 1x, &
            i9, 1x, i9, 1x, &
            i9, 1x, i9, 1x, i9, 1x, &
            i9, 1x, i9, 1x, i9, 1x, &
            i9, 1x, i9, 1x, &
            i9, 1x, i9, 1x, &
            i9, 1x, i9, 1x &
            i9, 1x, i9)

    end do

    300 continue
    close(iunit)
    call LIS_releaseUnitNumber(iunit)

  end subroutine USAF_gages_write_data

  ! Read gage data from file.  Acts as an alternative constructor.
  subroutine USAF_gages_read_data(this, filename, date10, alert_number)

    ! Imports
    use LIS_coreMod, only: LIS_masterproc
    use LIS_logMod, only: LIS_getNextUnitNumber, LIS_releaseUnitNumber, &
         LIS_alert, LIS_logunit

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(out) :: this
    character(*), intent(in) :: filename
    character(10), intent(in) :: date10
    integer, intent(inout) :: alert_number

    ! Locals
    integer :: nobs
    character(14), allocatable :: YYYYMMDDhhmmss(:)
    character(32), allocatable :: networks(:)
    character(32), allocatable :: platforms(:)
    character(2), allocatable :: wmocode_id(:)
    character(2), allocatable :: fipscode_id(:)
    integer, allocatable :: wmonumbers(:)
    integer, allocatable :: bsn(:)
    integer, allocatable :: lats(:)
    integer, allocatable :: lons(:)
    integer, allocatable :: amts24(:)
    integer, allocatable :: amts21(:)
    integer, allocatable :: amts18(:)
    integer, allocatable :: amts15(:)
    integer, allocatable :: amts12(:)
    integer, allocatable :: amts09(:)
    integer, allocatable :: amts06(:)
    integer, allocatable :: amts03(:)
    integer, allocatable :: amts02(:)
    integer, allocatable :: amts01(:)
    integer, allocatable :: amts00(:)
    integer, allocatable :: durations(:)
    integer, allocatable :: preswx(:)
    integer, allocatable :: pastwx(:)
    integer :: istat
    logical :: found
    integer :: i
    integer :: iunit
    character(255) :: message(20)

    message = ''
    call this%delete() ! Make sure structure is empty

    inquire(file=trim(filename), exist=found)
    if (.not. found) then
       write(LIS_logunit,*)'[WARN] Cannot find ', trim(filename)
       message(1) = '[WARN] Program:  LIS'
       message(2) = '  Routine: USAF_read_data'
       message(3) = '  Cannot find presav2 file ' // &
            trim(filename)
       message(4) = ' Observation count will be reduced'
       if (LIS_masterproc) then
          call LIS_alert('LIS.USAF_read_data', &
               alert_number, message)
          alert_number = alert_number + 1
       end if
       return
    end if
    iunit = LIS_getNextUnitNumber()

    open(iunit, file=trim(filename), iostat=istat)
    if (istat .ne. 0) then
       write(LIS_logunit,*)'[WARN] Problem opening ', trim(filename)
       message(1) = '[WARN] Program:  LIS'
       message(2) = '  Routine: USAF_gages_read_data'
       message(3) = '  Cannot open file ' // trim(filename)
       if (LIS_masterproc) then
          call LIS_alert('LIS.USAF_gages_read_data', &
               alert_number, message)
          alert_number = alert_number + 1
       end if
       return
    end if

    write(LIS_logunit,*) '[INFO] Reading ', trim(filename)
    read(iunit, *, iostat=istat) nobs
    if (istat .ne. 0) then
       write(LIS_logunit,*)'[WARN] Problem reading ', trim(filename)
       message(1) = '[WARN] Program:  LIS'
       message(2) = '  Routine: USAF_gages_read_data'
       message(3) = '  Problem reading file ' // trim(filename)
       if (LIS_masterproc) then
          call LIS_alert('LIS.USAF_gages_read_data', &
               alert_number, message)
          alert_number = alert_number + 1
       end if
       close(iunit)
       call LIS_releaseUnitNumber(iunit)
       return
    end if

    if (nobs .le. 0) then
       write(LIS_logunit,*)'[WARN] No precip obs found in ', &
            trim(filename)
       message(1) = '[WARN] Program:  LIS'
       message(2) = '  Routine: USAF_gages_read_data'
       message(3) = '  No precip obs found in file ' // trim(filename)
       if (LIS_masterproc) then
          call LIS_alert('LIS.USAF_gages_read_data', &
               alert_number, message)
          alert_number = alert_number + 1
       end if
       close(iunit)
       call LIS_releaseUnitNumber(iunit)
       return
    end if

    allocate(YYYYMMDDhhmmss(nobs))
    allocate(networks(nobs))
    allocate(platforms(nobs))
    allocate(wmocode_id(nobs))
    allocate(fipscode_id(nobs))
    allocate(wmonumbers(nobs))
    allocate(bsn(nobs))
    allocate(lats(nobs))
    allocate(lons(nobs))
    allocate(amts24(nobs))
    allocate(amts21(nobs))
    allocate(amts18(nobs))
    allocate(amts15(nobs))
    allocate(amts12(nobs))
    allocate(amts09(nobs))
    allocate(amts06(nobs))
    allocate(amts03(nobs))
    allocate(amts02(nobs))
    allocate(amts01(nobs))
    allocate(amts00(nobs))
    allocate(durations(nobs))
    allocate(preswx(nobs))
    allocate(pastwx(nobs))

    do i = 1, nobs
       read(iunit, 6000, iostat=istat, err=300, end=300) &
            YYYYMMDDhhmmss(i), &
            networks(i), platforms(i), &
            wmocode_id(i), fipscode_id(i), &
            wmonumbers(i), bsn(i), &
            lats(i), lons(i), &
            amts24(i), amts21(i), amts18(i), &
            amts15(i), amts12(i), amts09(i), &
            amts06(i), amts03(i), &
            amts02(i), amts01(i), &
            amts00(i), durations(i), &
            preswx(i), pastwx(i)
6000   format (a14, 1x, &
            a32, 1x, a32, 1x, &
            a2, 1x, a2, 1x, &
            i6, 1x, i6, 1x, &
            i9, 1x, i9, 1x, &
            i9, 1x, i9, 1x, i9, 1x, &
            i9, 1x, i9, 1x, i9, 1x, &
            i9, 1x, i9, 1x, &
            i9, 1x, i9, 1x, &
            i9, 1x, i9, 1x, &
            i9, 1x, i9)
    end do

    close(iunit)
    call LIS_releaseUnitNumber(iunit)
300 continue

    ! If read was successful, copy to USAF_gages_t structure.
    if (istat .eq. 0) then
       call this%new(date10, nobs, YYYYMMDDhhmmss, &
            networks, platforms, &
            wmocode_id, fipscode_id, &
            wmonumbers, bsn, lats, lons, &
            amts24, amts21, amts18, amts15, amts12, amts09, amts06, &
            amts03, amts02, amts01, &
            amts00, durations, preswx, pastwx)
    else
       write(LIS_logunit,*)'[WARN] Problem reading from ', &
            trim(filename)
       message(1) = '[WARN] Program:  LIS'
       message(2) = '  Routine: USAF_gages_read_data'
       message(3) = '  Problem reading from file ' // trim(filename)
       if (LIS_masterproc) then
          call LIS_alert('LIS.USAF_gages_read_data', &
               alert_number, message)
          alert_number = alert_number + 1
       end if
    end if

    ! Clean up
    deallocate(YYYYMMDDhhmmss)
    deallocate(networks)
    deallocate(platforms)
    deallocate(wmocode_id)
    deallocate(fipscode_id)
    deallocate(wmonumbers)
    deallocate(bsn)
    deallocate(lats)
    deallocate(lons)
    deallocate(amts24)
    deallocate(amts21)
    deallocate(amts18)
    deallocate(amts15)
    deallocate(amts12)
    deallocate(amts09)
    deallocate(amts06)
    deallocate(amts03)
    deallocate(amts02)
    deallocate(amts01)
    deallocate(amts00)
    deallocate(durations)
    deallocate(preswx)
    deallocate(pastwx)

    ! Handle arrays for unique networks in private method.  Also set
    ! cdms_flag.
    if (istat .eq. 0) then
       call set_unique_networks(this)
    end if

    !write(LIS_logunit,*)'[INFO] Read in ', nobs, ' gage reports'
  end subroutine USAF_gages_read_data

  ! Search USAF_gages_t type for presumed erroneous accumulations.
  subroutine USAF_gages_check_gross_errors(this)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this

    ! Locals
    character(10) :: date10
    integer :: india_lowlimit
    integer :: india_highlimit
    integer :: russia_lowlimit
    integer :: russia_highlimit
    integer :: threshold
    integer :: threshold_india
    integer :: nobs
    integer :: i

    if (this%cdms_flag) then
       india_lowlimit = CDMS_INDIA_LOWLIMIT
       india_highlimit = CDMS_INDIA_HIGHLIMIT
       russia_lowlimit = CDMS_RUSSIA_LOWLIMIT
       russia_highlimit = CDMS_RUSSIA_HIGHLIMIT
    else
       india_lowlimit = JMOBS_INDIA_LOWLIMIT
       india_highlimit = JMOBS_INDIA_HIGHLIMIT
       russia_lowlimit = JMOBS_RUSSIA_LOWLIMIT
       russia_highlimit = JMOBS_RUSSIA_HIGHLIMIT
    end if
    date10 = this%date10

    nobs = this%nobs
    do i = 1, nobs

       ! Check for negatives
       if (this%amts24(i) < 0) this%amts24(i) = MISSING
       if (this%amts21(i) < 0) this%amts21(i) = MISSING
       if (this%amts18(i) < 0) this%amts18(i) = MISSING
       if (this%amts15(i) < 0) this%amts15(i) = MISSING
       if (this%amts12(i) < 0) this%amts12(i) = MISSING
       if (this%amts09(i) < 0) this%amts09(i) = MISSING
       if (this%amts06(i) < 0) this%amts06(i) = MISSING
       if (this%amts03(i) < 0) this%amts03(i) = MISSING
       if (this%amts02(i) < 0) this%amts02(i) = MISSING
       if (this%amts01(i) < 0) this%amts01(i) = MISSING
       if (this%amts00(i) < 0) this%amts00(i) = MISSING

       ! Check for large errors
       if (this%amts24(i) > MAX_PCP_24HR) this%amts24(i) = MISSING
       if (this%amts21(i) > MAX_PCP_21HR) this%amts21(i) = MISSING
       if (this%amts18(i) > MAX_PCP_18HR) this%amts18(i) = MISSING
       if (this%amts15(i) > MAX_PCP_15HR) this%amts15(i) = MISSING
       if (this%amts12(i) > MAX_PCP_12HR) this%amts12(i) = MISSING
       if (this%amts09(i) > MAX_PCP_9HR) this%amts09(i) = MISSING
       if (this%amts06(i) > MAX_PCP_6HR) this%amts06(i) = MISSING
       if (this%amts03(i) > MAX_PCP_3HR) this%amts03(i) = MISSING
       if (this%amts02(i) > MAX_PCP_2HR) this%amts02(i) = MISSING
       if (this%amts01(i) > MAX_PCP_1HR) this%amts01(i) = MISSING

       ! Checking miscellaneous accumulations is more complicated.
       ! First, see if we have a WMO established duration.
       threshold = set_precip_duration_threshold(this%durations(i))
       if (threshold .ne. MISSING) then
          if (this%amts00(i) > threshold) then
             this%amts00(i) = MISSING
             this%durations(i) = MISSING
          end if
          cycle
       end if

       ! Try the "Russia" rule.  Three sets of tests before invoking.
       if (this%durations(i) .eq. 0 .or. &
            this%durations(i) .eq. MISSING) then

          if ((this%wmonumbers(i) .ge. russia_lowlimit .and. &
               this%wmonumbers(i) .le. russia_highlimit) .or. &
               (this%bsn(i) .ge. russia_lowlimit .and. &
                this%bsn(i) .le. russia_highlimit) ) then

             if (date10(9:10) .eq. "00" .or. &
                  date10(9:10) .eq. "06" .or. &
                  date10(9:10) .eq. "12" .or. &
                  date10(9:10) .eq. "18") then

                ! Russia rule (actually, post-Soviet rule).  In the
                ! past, USAF personnel noted many "Russian" obs were
                ! only reported at 00Z and 12Z, and did not use a valid
                ! duration flag. These are decoded in the miscellaneous
                ! accumulation (amts00). In AGRMET they were treated as
                ! 12-hr accumulations.  More recently, a spot check in
                ! July 2021 showed similar invalid durations in
                ! Kazakhstan, Kyrgyzstan, Georgia, Tajiskistan, and
                ! Turkenistan at 06Z and 18Z.  Since WMO Region VI
                ! (Europe) expects 12-hr accumulations at 06Z and 18Z, we
                ! assume these rogue reports are also 12-hr
                ! accumulations.

                if (this%amts00(i) > MAX_PCP_12HR) then
                   this%amts00(i) = MISSING
                   this%durations(i) = MISSING
                end if
                cycle
             end if
          end if
       end if ! Russia rule

       ! Try the India rule.  Two sets of tests.
       if (this%durations(i) .eq. 0 .or. &
            this%durations(i) .eq. MISSING) then

          if ((this%wmonumbers(i) .ge. india_lowlimit .and. &
               this%wmonumbers(i) .le. india_highlimit) .or. &
               (this%bsn(i) .ge. india_lowlimit .or. &
               this%bsn(i) .le. india_highlimit)) then

             ! If these are from India, we assume duration is from 03Z.
             ! But, we will only accept reports at three-hour intervals.
             threshold_india = set_india_precip_threshold(date10(9:10))
             if (threshold_india .ne. MISSING) then
                if (this%amts00(i) > threshold_india) then
                   this%amts00(i) = MISSING
                   this%durations(i) = MISSING
                end if
                cycle
             end if
          end if
       end if ! India rule

       ! If all else fails, assume 24-hr accumulation
       if (this%amts00(i) > MAX_PCP_24HR) then
          this%amts00(i) = MISSING
          this%durations(i) = MISSING
       end if

    end do
  end subroutine USAF_gages_check_gross_errors

  ! Check the observations for bad lat/lon values.  If found, reset the
  ! data to missing.
  subroutine check_latlons(this)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this

    ! Locals
    logical :: bad_latlon
    integer :: nobs
    integer :: i

    nobs = this%nobs
    do i = 1, nobs
       bad_latlon = .false.

       ! Sometimes a station report is listed with a lat/lon of 0,
       ! putting it in the Atlantic Ocean.  We check and flag.
       if (this%lats(i) == 0 .and. this%lons(i) == 0) bad_latlon = .true.

       ! If we don't have a good lat/lon, the information from the
       ! station is useless.  So remove it.
       if (bad_latlon) then
          this%amts24(i) = MISSING
          this%amts21(i) = MISSING
          this%amts18(i) = MISSING
          this%amts15(i) = MISSING
          this%amts12(i) = MISSING
          this%amts09(i) = MISSING
          this%amts06(i) = MISSING
          this%amts03(i) = MISSING
          this%amts02(i) = MISSING
          this%amts01(i) = MISSING
          this%amts00(i) = MISSING
          this%durations(i) = MISSING
          this%preswx(i) = MISSING
          this%pastwx(i) = MISSING
       end if
    end do
  end subroutine check_latlons

  ! Search USAF_gages_t type for unique networks, and establish index
  ! bounds for each network.  Assumes data are sorted.  Borrows logic
  ! from AGRMET_processobs.
  subroutine set_unique_networks(this)

    ! Imports
    use LIS_logMod, only: LIS_logunit, LIS_endrun

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this

    ! Locals
    character(32) :: prior_net
    integer :: net_count
    integer :: nobs
    integer :: i

    if (.not. allocated(this%unique_networks)) then
       allocate(this%unique_networks(MAX_UNIQUE_NETWORKS))
    end if
    this%unique_networks = MISSING_NAME
    if (.not. allocated(this%firsts)) then
       allocate(this%firsts(MAX_UNIQUE_NETWORKS))
    end if
    this%firsts = MISSING
    if (.not. allocated(this%lasts)) then
       allocate(this%lasts(MAX_UNIQUE_NETWORKS))
    end if
    this%lasts = MISSING

    prior_net = MISSING_NAME

    nobs = this%nobs
    do i = 1, nobs
       if (prior_net .eq. MISSING_NAME) then
          this%firsts(1) = i
          this%unique_networks(i) = this%networks(i)
          net_count = 1
          prior_net = this%networks(i)
       else if (this%networks(i) .ne. prior_net) then
          this%lasts(net_count) = i - 1
          net_count = net_count + 1
          if (net_count .gt. MAX_UNIQUE_NETWORKS) then
             write(LIS_logunit,*) '[ERR] Too many unique networks!'
             call LIS_endrun
          end if
          this%firsts(net_count) = i
          this%unique_networks(net_count) = this%networks(i)
          prior_net = this%networks(i)
       end if
    end do
    this%lasts(net_count) = this%nobs

    this%num_unique_networks = net_count

    ! Mark if these obs are all legacy CDMS output
    this%cdms_flag = .false.
    if (net_count .eq. 1) then
       if (this%unique_networks(net_count) .eq. "CDMS") then
          this%cdms_flag = .true.
       end if
    end if
  end subroutine set_unique_networks

  ! Search USAF_gages_t type for given network and platform, and
  ! return index.  Based on AGRMET_pcpobs_search
  function search_pcpobs(this, network, plat_id, &
       wmocode_id, fipscode_id) result(index)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(in) :: this
    character(32), intent(in) :: network
    character(32), intent(in) :: plat_id
    character(2), intent(in) :: wmocode_id
    character(2), intent(in) :: fipscode_id

    ! Result
    integer :: index

    ! Locals
    integer :: first
    integer :: middle
    integer :: last
    integer :: net_index
    logical :: found

    index = MISSING
    first = 1
    last = this%nobs
    found = .false.

    net_index = 1

    ! Select the appropriate 'first' and 'last' search indices based on
    ! requested network
    do while ((net_index .le. this%num_unique_networks) .and. &
         (.not. found))
       if (network .eq. this%unique_networks(net_index)) then
          first = this%firsts(net_index)
          last = this%lasts(net_index)
          found = .true.
       end if
       net_index = net_index + 1
    end do
    if (.not. found) return

    found = .false. ! Reuse for finding station ID below
    do
       if ((first .gt. last) .or. found) return
       ! Use binary search
       middle = (first + last) / 2
       if (plat_id .lt. this%platforms(middle)) then
          last = middle - 1
       elseif (plat_id .gt. this%platforms(middle)) then
          first = middle + 1
       else
          ! Make sure country code matches
          if (this%wmocode_id(middle) .eq. wmocode_id .and. &
               this%fipscode_id(middle) .eq. fipscode_id) then
             found = .true.
             index = middle
          end if
          exit
       end if
    end do

  end function search_pcpobs

  ! Infer zero precip for certain durations if the accumulations are
  ! missing and if present and past weather codes support no precip.
  subroutine USAF_gages_use_preswx_pastwx(this)

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_gages_t), intent(inout) :: this

    ! Locals
    integer :: sweden_lowlimit
    integer :: sweden_highlimit
    integer :: denmark_lowlimit
    integer :: denmark_highlimit
    integer :: pastwx_duration
    integer :: nobs
    integer :: i

    if (this%cdms_flag) then
       sweden_lowlimit = CDMS_SWEDEN_LOWLIMIT
       sweden_highlimit = CDMS_SWEDEN_HIGHLIMIT
       denmark_lowlimit = CDMS_DENMARK_LOWLIMIT
       denmark_highlimit = CDMS_DENMARK_HIGHLIMIT
    else
       sweden_lowlimit = JMOBS_SWEDEN_LOWLIMIT
       sweden_highlimit = JMOBS_SWEDEN_HIGHLIMIT
       denmark_lowlimit = JMOBS_DENMARK_LOWLIMIT
       denmark_highlimit = JMOBS_DENMARK_HIGHLIMIT
    end if

    nobs = this%nobs
    do i = 1, nobs

       ! If we don't have a valid past weather duration, we can't know
       ! which precip accumulation to adjust.
       pastwx_duration = this%pastwx_durations(i)
       if (pastwx_duration .eq. MISSING) cycle

       if (this%preswx(i) .eq. MISSING) cycle
       if (this%pastwx(i) .eq. MISSING) cycle

       ! We will only use past and present weather if the appropriate
       ! accumulations are missing.
       if (pastwx_duration .eq. 6) then
          if (this%amts06(i) .ne. MISSING) cycle
          if (this%amts03(i) .ne. MISSING) cycle
          if (this%amts02(i) .ne. MISSING) cycle
          if (this%amts01(i) .ne. MISSING) cycle
       end if

       if (pastwx_duration .eq. 3) then
          if (this%amts03(i) .ne. MISSING) cycle
          if (this%amts02(i) .ne. MISSING) cycle
          if (this%amts01(i) .ne. MISSING) cycle
       end if

       ! Exclude Danish stations, since reports at non-standard
       ! international hours have past weather of only 1 hour.
       if ((this%wmonumbers(i) .ge. denmark_lowlimit .and. &
            this%wmonumbers(i) .le. denmark_highlimit) .or. &
            (this%bsn(i) .ge. denmark_lowlimit .and. &
            this%bsn(i) .le. denmark_highlimit)) then
          cycle
       end if

       ! Exclude Swedish stations, since durations for past weather
       ! vary from one to six hours depending on when the report was
       ! made.
       if ((this%wmonumbers(i) .ge. sweden_lowlimit .and. &
            this%wmonumbers(i) .le. sweden_highlimit) .or. &
            (this%bsn(i) .ge. sweden_lowlimit .and. &
            this%bsn(i) .le. sweden_highlimit)) then
          cycle
       end if

       ! Exclude Antarctic stations. This is because (1) Australian
       ! Antarctic stations are known to have varying past weather
       ! durations due to staffing; (2) isolating the Australian
       ! stations is difficult due to WMO station ID assignments; and
       ! (3) Antarctic stations will probably be subfreezing and have
       ! their precip accumulations rejected anyway.  For simplicity, we
       ! skip any report south of 60S.
       if (this%lons(i) .le. -6000) cycle

       ! WMO SYNOP reports use different code tables depending on if the
       ! station is manned or automatic.  Unfortunately, the indicator
       ! for manned vs automatic is not included in the USAF decoded
       ! files, so we must be conservative and check for precipitation
       ! codes from either scenario. The decision criteria are kept
       ! separate in case the USAF file format is changed in the future
       ! to include this information.  In addition, BUFR has it's own
       ! tables, and it is not clear from the USAF decoded files if BUFR
       ! or SYNOP was decoded.

       ! First, check present weather for manned station (WMO Code Table
       ! 4677; also WMO BUFR Code Table 0 20 003).
       select case (this%preswx(i))
       case (50:99) ! Precipitation at station at time of observation
          cycle
       case (20:27, 29) ! Precip at station in last hour, but not at
          cycle         ! time of observation
       end select

       ! Next, check present weather for automatic station (WMO Code
       ! Table 4680; see also WMO BUFR Code Table 0 20 003)
       select case (this%preswx(i))
       case (40:48, 50:58, 60:68, 70:78, 89, 90:96)
          cycle ! Precipitation at station at time of observation
       case (21:26)
          cycle ! Precip at station in last hour, but not at time of
                ! observation
       case (140:148, 150:158, 160:168, 170:178, 180:188, 190:196)
          cycle ! BUFR codes for precip at automated station at time of
                ! observation
       case (121:126)
          cycle ! BUFR codes for precip at station in last hour, but
                ! not at time of observation
       case (250:257, 259, 260:267, 270:279, 280:291)
          cycle ! More descriptive BUFR codes for precip at station.
       case (510:511)
          cycle ! BUFR codes for missing data
       end select

       ! Next, check past weather for manned station (WMO Code Table
       ! 4561)
       select case (this%pastwx(i))
       case (5:9)
          cycle
       end select

       ! Next, check past weather for automatic station (WMO Code Table
       ! 4531, and WMO BUFR Code Table 0 20 004 / 0 20 005)
       select case (this%pastwx(i))
       case (4:9) ! WMO Code Table 4531
          cycle
       case (14:19) ! BUFR Code Tables
          cycle
       case (31) ! BUFR code for missing
          cycle
       end select

       ! At this point, we can infer zero precipitation for the
       ! appropriate duration
       if (pastwx_duration .eq. 6) then
          this%amts06(i) = 0
          this%amts03(i) = 0
          this%amts02(i) = 0
          this%amts01(i) = 0
       else if (pastwx_duration .eq. 3) then
          this%amts03(i) = 0
          this%amts02(i) = 0
          this%amts01(i) = 0
       end if
    end do

    call this%reconcile_self()
    call this%fill_gaps()
  end subroutine USAF_gages_use_preswx_pastwx

  ! Method for setting pastwx duration for each report.
  subroutine USAF_gages_set_pastwx_durations(this)
    implicit none
    class(USAF_gages_t), intent(inout) :: this
    integer :: nobs
    character(10) :: date10
    integer :: i
    nobs = this%nobs
    do i = 1, nobs
       if (this%YYYYMMDDhhmmss(i) .eq. "NULL") then
          date10 = this%date10
       else
          date10 = this%YYYYMMDDhhmmss(i)(1:10)
       end if
       this%pastwx_durations(i) = set_pastwx_duration(date10)
    end do
  end subroutine USAF_gages_set_pastwx_durations

  ! Function for setting duration of past weather report.  This is based
  ! on WMO SYNOP definition.  (Two and one hour durations are also
  ! allowed by the WMO, but these require SYNOP reports every two or one
  ! hours; and it is not clear from the USAF decoded files when that
  ! occurs.)
  function set_pastwx_duration(date10) result (duration)
    implicit none
    character(10), intent(in) :: date10
    integer :: duration
    select case (date10(9:10))
    case ('00', '06', '12', '18')
       duration = 6
    case ('03', '09', '15', '21')
       duration = 3
    case default
       duration = MISSING
    end select
  end function set_pastwx_duration

  ! Return appropriate precipitation threshold based on duration.
  ! NOTE:  We currently skip 1 hour and 2-hour accumulations.
  function set_precip_duration_threshold(duration) result (threshold)

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: duration

    ! Return variable
    integer :: threshold

    ! Select appropriate threshold
    select case (duration)
    case (1)
       threshold = MAX_PCP_1HR
    case (2)
       threshold = MAX_PCP_2HR
    case (3)
       threshold = MAX_PCP_3HR
    case (6)
       threshold = MAX_PCP_6HR
    case (9)
       threshold = MAX_PCP_9HR
    case (12)
       threshold = MAX_PCP_12HR
    case (15)
       threshold = MAX_PCP_15HR
    case (18)
       threshold = MAX_PCP_18HR
    case (21) ! Not a WMO established duration, but we'll check anyway.
       threshold = MAX_PCP_21HR
    case (24)
       threshold = MAX_PCP_24HR
    case default
       threshold = MISSING
    end select

  end function set_precip_duration_threshold

  ! Return appropriate precipitation threshold based on report time in
  ! India
  function set_india_precip_threshold(utc_hour) result (threshold)

    ! Defaults
    implicit none

    ! Arguments
    character(2), intent(in) :: utc_hour

    ! Return variable
    integer :: threshold

    ! Select appropriate threshold.  India obs are assumed to accumulate
    ! from 03Z
    select case (utc_hour)
    case ("04")
       threshold = MAX_PCP_1HR
    case ("05")
       threshold = MAX_PCP_2HR
    case ("06")
       threshold = MAX_PCP_3HR
    case ("09")
       threshold = MAX_PCP_6HR
    case ("12")
       threshold = MAX_PCP_9HR
    case ("15")
       threshold = MAX_PCP_12HR
    case ("18")
       threshold = MAX_PCP_15HR
    case ("21")
       threshold = MAX_PCP_18HR
    case ("00")
       threshold = MAX_PCP_21HR
    case ("03")
       threshold = MAX_PCP_24HR
    case default
       threshold = MISSING
    end select

  end function set_india_precip_threshold

  ! Method for copying appropriate data to a USAF_ObsData structure
  ! for use in data assimilation.
  subroutine USAF_copy_to_usaf_obsdata(this, hr, gage_sigma_o_sqr, &
       precipObs)

    ! Imports
    use LIS_logMod, only: LIS_logunit, LIS_endrun
    use USAF_bratsethMod, only: USAF_ObsData, USAF_assignObsData

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_Gages_t), intent(in) :: this
    integer, intent(in) :: hr
    real, intent(in) :: gage_sigma_o_sqr
    type(USAF_ObsData), intent(inout) :: precipObs

    ! Locals
    integer :: i
    integer :: num_obs_copied

    ! Sanity checks
    if (this%nobs == 0) return
    if (hr .ne. 6 .and. hr .ne. 12) then
       write(LIS_logunit,*) &
            '[ERR] Invalid hour passed to USAF_copy_to_usaf_obsdata'
       write(LIS_logunit,*) &
            '[ERR] Should be 6 or 12, received ', hr
       call LIS_endrun()
    end if

    num_obs_copied = 0
    if (hr == 6) then
       do i = 1, this%nobs
          if (this%amts06(i) < 0) cycle
          call USAF_assignObsData(precipObs, &
               this%networks(i), &
               this%platforms(i), &
               real(this%amts06(i)) * 0.1, &
               real(this%lats(i)) * 0.01, &
               real(this%lons(i)) * 0.01, &
               gage_sigma_o_sqr, 0.)
          num_obs_copied = num_obs_copied + 1
       end do
       write(LIS_logunit,*)'[INFO] Copied ', num_obs_copied, &
            ' 6-hr gage reports'
    else if (hr == 12) then
       do i = 1, this%nobs
          if (this%amts12(i) < 0) cycle
          call USAF_assignObsData(precipObs, &
               this%networks(i), &
               this%platforms(i), &
               real(this%amts12(i)) * 0.1, &
               real(this%lats(i)) * 0.01, &
               real(this%lons(i)) * 0.01, &
               gage_sigma_o_sqr, 0.)
          num_obs_copied = num_obs_copied + 1
       end do
       write(LIS_logunit,*)'[INFO] Copied ', num_obs_copied, &
            ' 12-hr gage reports'
    end if


  end subroutine USAF_copy_to_usaf_obsdata
end module USAF_GagesMod
