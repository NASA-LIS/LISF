! ROUTINE: USAF_getpcpcobs
!
! REVISION HISTORY:
! 17 Oct 2023 Initial version.  Eric Kemp/SSAI/NASA.
!
! DESCRIPTION:
! New routine to retrieve precip observations and process with new
! USAF_Gages library.  Supports legacy hemispheric and new global
! preobs file formats, and new presav2 file format.  Borrows heavily
! from older AGRMET_getpcpobs subroutine.
!--------------------------------------------------------------------------

subroutine USAF_getpcpobs(n, j6hr, month, use_twelve, pcp_src, &
     use_expanded_station_ids, alert_number, precip6, precip12)

  ! Imports
  use AGRMET_forcingMod, only: agrmet_struc
  use ESMF
  use LIS_coreMod, only: LIS_rc
  use LIS_logMod, only: LIS_logunit
  use LIS_timeMgrMod, only: LIS_tick, LIS_julhr_date
  use USAF_bratsethMod, only: USAF_ObsData, USAF_setbratsethprecipstats
  use USAF_GagesMod, only : USAF_Gages_t
  use USAF_PreobsReaderMod, only: USAF_read_preobs

  ! Defaults
  implicit none

  ! Arguments
  integer, intent(in) :: n
  integer, intent(in) :: j6hr
  integer, intent(in) :: month
  logical, intent(in) :: use_twelve
  character*6, intent(in) :: pcp_src(4)
  integer, intent(in) :: use_expanded_station_ids
  integer,  intent(inout):: alert_number
  type(USAF_ObsData), intent(inout) :: precip6
  type(USAF_ObsData), intent(inout) :: precip12

  ! Locals
  integer :: j3hr
  integer :: k
  integer :: yr, mo, da, hr
  character*255 :: preobsdir
  character*8 :: yyyymmdd
  character*10 :: yyyymmddhh
  character*255 :: presav_filename
  logical :: file_exists
  type(USAF_Gages_t) :: obscur

  if (use_twelve) then
     k = 3
  else
     k = 1
  end if

  do j3hr = j6hr+3, j6hr+6, 3

     ! Set Bratseth error statistics based on source of background field.
     call USAF_setBratsethPrecipStats(pcp_src(k), n)

     call LIS_julhr_date(j3hr, yr, mo, da, hr)
     if (agrmet_struc(n)%use_timestamp == 1) then
        write(unit=yyyymmdd, fmt='(i4.4, i2.2, i2.2)') &
             yr, mo, da
        preobsdir = trim(agrmet_struc(n)%agrmetdir) // '/' &
             // yyyymmdd // '/' &
             // trim(agrmet_struc(n)%cdmsdir) // '/'
     else
        preobsdir = trim(agrmet_struc(n)%agrmetdir) // '/' &
             // trim(agrmet_struc(n)%cdmsdir) // '/'
     end if

     ! Read appropriate preobs file(s), intercompare with older presav2
     ! files, and create new presav2 file for current date/time.
     call USAF_read_preobs(preobsdir, &
          trim(agrmet_struc(n)%analysisdir), &
          agrmet_struc(n)%use_timestamp, yr, mo, da, hr, &
          use_expanded_station_ids, alert_number)

     ! If this is a synoptic time, read the presav2 file back in and
     ! populate the appropriate USAF_ObsData object.
     if ( mod(j3hr, 6) == 0 ) then
        write(presav_filename,'(A,A,i4.4,i2.2,i2.2,i2.2)') &
             trim(agrmet_struc(n)%analysisdir), '/presav2.03hr.', &
             yr, mo, da, hr
        inquire(file=presav_filename, exist=file_exists)
        if (file_exists) then
           write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') &
                yr, mo, da, hr
           call obscur%read_data(presav_filename, yyyymmddhh, &
                alert_number)
           if (use_twelve) then
              call obscur%copy_to_usaf_obsdata(12, &
                   agrmet_struc(n)%bratseth_precip_gauge_sigma_o_sqr, &
                   precip12)
           else
              call obscur%copy_to_usaf_obsdata(6, &
                   agrmet_struc(n)%bratseth_precip_gauge_sigma_o_sqr, &
                   precip6)
           end if
        end if
     end if
     k = k + 1
  end do

end subroutine USAF_getpcpobs
