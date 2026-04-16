!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: LDT_wsf_oplMod
!
! DESCRIPTION: IMPROVED hour filtering with integrity-based file prioritization
!              Prioritizes complete files (i0) over incomplete (i1) and uses
!              file size as tiebreaker
!
!-------------------------------------------------------------------------

module LDT_wsf_oplMod

  use LDT_constantsMod, only: LDT_CONST_PATH_LEN

  implicit none
  private

  public :: LDT_wsf_oplInit
  public :: LDT_wsf_oplRun
  public :: wsf_file_info

  type, public :: wsf_opl_dec
    character(len=LDT_CONST_PATH_LEN) :: WSFdir
    character(len=LDT_CONST_PATH_LEN) :: WSFoutdir
    character*10  :: date_curr

    real*4, allocatable :: ARFS_TB_10V(:,:)
    real*4, allocatable :: ARFS_TB_10H(:,:)
    real*4, allocatable :: ARFS_TB_18V(:,:)
    real*4, allocatable :: ARFS_TB_18H(:,:)
    real*4, allocatable :: ARFS_TB_23V(:,:)
    real*4, allocatable :: ARFS_TB_23H(:,:)
    real*4, allocatable :: ARFS_TB_36V(:,:)
    real*4, allocatable :: ARFS_TB_36H(:,:)
    real*4, allocatable :: ARFS_TB_89V(:,:)
    real*4, allocatable :: ARFS_TB_89H(:,:)
    real*4, allocatable :: ARFS_LAND_FRAC(:,:)
    integer*1, allocatable :: ARFS_QUALITY_FLAG(:,:)
    integer*4, allocatable :: ARFS_SNOW(:,:)
    integer*4, allocatable :: ARFS_PRECIP(:,:)
    integer*4, allocatable :: ARFS_SAMPLE_V(:,:)
    integer*4, allocatable :: ARFS_SAMPLE_H(:,:)
    integer       :: filter_snow_precip  ! NEW 1=filter (default, SM retrieval), 0=keep all footprints (snow depth)
  end type wsf_opl_dec

  type(wsf_opl_dec), public :: WSFopl

  type :: wsf_file_info
    character(len=LDT_CONST_PATH_LEN) :: filename
    character*8   :: date_str
    character*6   :: start_time
    character*6   :: end_time
    character*2   :: copy_num
    integer       :: copy_int
    character*1   :: integrity     ! NEW: 'i0' = complete, 'i1' = incomplete
    integer       :: integrity_int  ! NEW: 0 = complete, 1 = incomplete
    integer*8     :: file_size     ! NEW: file size in bytes
    character*2   :: hour_str
    integer       :: start_hour
    integer       :: end_hour
    integer       :: pass_type     ! NEW: 1=asc, -1=desc, 0=unknown
  end type wsf_file_info

contains

  subroutine LDT_wsf_oplInit()
    use ESMF
    use LDT_coreMod, only: LDT_config
    use LDT_logMod, only: LDT_logunit, LDT_verify
    
    implicit none
    integer :: rc
    character(len=255) :: cfg_entry
    
    write(LDT_logunit,*) '[INFO] ========================================='
    write(LDT_logunit,*) '[INFO] Initializing WSF low resolution resampling'
    write(LDT_logunit,*) '[INFO] WITH IMPROVED integrity-based filtering'
    write(LDT_logunit,*) '[INFO] ========================================='
    
    cfg_entry = "WSF valid date (YYYYMMDDHH):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, WSFopl%date_curr, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    
    cfg_entry = "WSF input directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, WSFopl%WSFdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    
    cfg_entry = "WSF output directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, WSFopl%WSFoutdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    
    ! Read snow/precip filter option (default = 1 = filter ON for SM retrieval)
    cfg_entry = "WSF filter snow and precip footprints:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    if (rc == 0) then
      call ESMF_ConfigGetAttribute(LDT_config, WSFopl%filter_snow_precip, rc=rc)
      if (rc /= 0) then
        WSFopl%filter_snow_precip = 1  ! Default: filter ON
      endif
    else
      WSFopl%filter_snow_precip = 1    ! Default: filter ON
    endif

    write(LDT_logunit,*) '[INFO] WSF valid date: ', trim(WSFopl%date_curr)
    write(LDT_logunit,*) '[INFO] WSF input directory: ', trim(WSFopl%WSFdir)
    write(LDT_logunit,*) '[INFO] WSF output directory: ', trim(WSFopl%WSFoutdir)
    if (WSFopl%filter_snow_precip == 1) then
      write(LDT_logunit,*) '[INFO] WSF snow/precip filtering: ON (SM retrieval mode)'
    else
      write(LDT_logunit,*) '[INFO] WSF snow/precip filtering: OFF (all footprints kept)'
    endif

  end subroutine LDT_wsf_oplInit

  subroutine LDT_wsf_oplRun(n)
    use LDT_coreMod, only: LDT_rc
    use LDT_logMod
    use TOOLSUBS_WSF

    
    implicit none
    integer, intent(in) :: n
    
    type(wsf_file_info), allocatable :: wsf_files(:)
    type(wsf_file_info), allocatable :: filtered_files(:)
    type(wsf_file_info), allocatable :: hour_group(:)
    character(len=LDT_CONST_PATH_LEN) :: fname
    character*10 :: tmp
    character*2 :: target_hour
    character*8 :: yyyymmdd
    integer :: ftn, ierr, fi, i, j
    integer :: n_filtered, n_hour_group
    integer :: target_hour_int
    logical :: file_exists
    integer :: nscans, nfovs, nchans
    real*4, allocatable :: tb_lowres(:,:,:), lat_in(:,:), lon_in(:,:)
    real*4, allocatable :: land_frac_low(:,:), earth_inc_angle(:,:,:)
    integer*1, allocatable :: quality_flag_in(:,:)
    integer*4, allocatable :: snow_in(:,:), precip_in(:,:)
    real*4, allocatable :: chan_frequencies(:)
    character*1, allocatable :: chan_polarizations(:)
    integer :: n_asc, n_desc
    type(wsf_file_info), allocatable :: filtered_asc(:), filtered_desc(:)

        
    external :: WSF_ARFS_RESAMPLE_HOURLY
    
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] Starting WSF Low Resolution Resampling'
    write(LDT_logunit,*) '[INFO] ========================================'
    
    ! Extract date and target hour
    yyyymmdd = WSFopl%date_curr(1:8)
    target_hour = WSFopl%date_curr(9:10)
    read(target_hour, '(I2)') target_hour_int
    
    write(LDT_logunit,*) '[INFO] Target date: ', yyyymmdd
    write(LDT_logunit,*) '[INFO] Target hour: ', target_hour, 'H (', target_hour_int, ')'
    
    ! Search for files
    tmp = trim(WSFopl%date_curr)    ! date_curr is already YYYYMMDDHH
    call search_WSF_files(WSFopl%WSFdir, WSFopl%date_curr, tmp)   ! tmp already = date_curr
    
    allocate(wsf_files(1000))
    
    ftn = LDT_getNextUnitNumber()
    open(ftn, file='WSF_filelist_'//trim(tmp)//'.dat', status='old', iostat=ierr)
    
    if (ierr /= 0) then
      write(LDT_logunit,*) '[ERR] Cannot open WSF_filelist_'//trim(tmp)//'.dat'
      return
    endif
    
    ! Read and parse all WSF filenames
    fi = 0
    do while (ierr == 0)
      read(ftn, '(a)', iostat=ierr) fname
      if (ierr /= 0) exit
      
      if (len_trim(fname) == 0) cycle
      if (index(fname, 'No such file') > 0 .or. &
          index(fname, 'cannot access') > 0) then
        cycle
      endif
      
      inquire(file=trim(fname), exist=file_exists)
      
      if (.not. file_exists) then
        write(LDT_logunit,*) '[WARN] File does not exist: ', trim(fname)
        cycle
      endif

      fi = fi + 1
      if (fi <= 1000) then
        call parse_wsf_filename_improved(fname, wsf_files(fi))
      endif
      
    end do
    call LDT_releaseUnitNumber(ftn)
    
    write(LDT_logunit,*) '[INFO] Found ', fi, ' total WSF files'
    
    if (fi == 0) then
      write(LDT_logunit,*) '[WARN] No WSF files found for date: ', WSFopl%date_curr
      deallocate(wsf_files)
      return
    endif
    
    ! Filter duplicates with improved logic
    allocate(filtered_files(fi))
    call filter_duplicate_files_improved(wsf_files(1:fi), fi, filtered_files, n_filtered)
    
    write(LDT_logunit,*) '[INFO] After duplicate filtering: ', n_filtered, ' files'
    
    ! Filter by target hour (files that overlap with target hour)
    ! Filter by target hour
    allocate(hour_group(n_filtered))
    n_hour_group = 0
    
    do i = 1, n_filtered
      if (file_overlaps_hour(filtered_files(i), target_hour_int)) then
        n_hour_group = n_hour_group + 1
        hour_group(n_hour_group) = filtered_files(i)
      endif
    end do
    
    write(LDT_logunit,*) '[INFO] Files overlapping hour ', target_hour, ': ', n_hour_group
    
    ! Detect pass type only for hour-relevant files
    n_asc = 0
    n_desc = 0
    
    do i = 1, n_hour_group
      call get_wsf_data_with_flags(hour_group(i)%filename, &
          tb_lowres, lat_in, lon_in, land_frac_low, quality_flag_in, &
          earth_inc_angle, snow_in, precip_in, &
          nscans, nfovs, nchans, chan_frequencies, chan_polarizations, ierr)
      
      if (ierr == 0) then
        hour_group(i)%pass_type = detect_pass_type(lat_in, nscans, nfovs)
        deallocate(tb_lowres, lat_in, lon_in, land_frac_low, quality_flag_in, &
                   earth_inc_angle, snow_in, precip_in, &
                   chan_frequencies, chan_polarizations)
      else
        hour_group(i)%pass_type = 0
        if (allocated(tb_lowres))          deallocate(tb_lowres)
        if (allocated(lat_in))             deallocate(lat_in)
        if (allocated(lon_in))             deallocate(lon_in)
        if (allocated(land_frac_low))      deallocate(land_frac_low)
        if (allocated(quality_flag_in))    deallocate(quality_flag_in)
        if (allocated(earth_inc_angle))    deallocate(earth_inc_angle)
        if (allocated(snow_in))            deallocate(snow_in)
        if (allocated(precip_in))          deallocate(precip_in)
        if (allocated(chan_frequencies))   deallocate(chan_frequencies)
        if (allocated(chan_polarizations)) deallocate(chan_polarizations)
      endif
    end do
    
    ! Separate into ASC and DESC groups
    ! Pack ASC at front, DESC at back of hour_group
    n_asc = 0
    n_desc = 0
    do i = 1, n_hour_group
      if (hour_group(i)%pass_type == 1) then
        n_asc = n_asc + 1
        write(LDT_logunit,*) '[INFO] ASC file: ', trim(hour_group(i)%filename)
      else if (hour_group(i)%pass_type == -1) then
        n_desc = n_desc + 1
        write(LDT_logunit,*) '[INFO] DESC file: ', trim(hour_group(i)%filename)
      else
        write(LDT_logunit,*) '[WARN] Unknown pass type, skipping: ', trim(hour_group(i)%filename)
      endif
    end do
    
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] Processing hour: ', target_hour, 'H'
    write(LDT_logunit,*) '[INFO] Ascending: ', n_asc, ', Descending: ', n_desc
    write(LDT_logunit,*) '[INFO] ========================================'
    
    if (n_asc > 0) then
      allocate(filtered_asc(n_asc))
      j = 0
      do i = 1, n_hour_group
        if (hour_group(i)%pass_type == 1) then
          j = j + 1
          filtered_asc(j) = hour_group(i)
        endif
      end do
      write(LDT_logunit,*) '[INFO] Processing ASCENDING: ', n_asc, ' files'
      call WSF_ARFS_RESAMPLE_HOURLY(filtered_asc, n_asc, &
                                    WSFopl%WSFoutdir, yyyymmdd, target_hour, n, 'ASC', &
                                    WSFopl%filter_snow_precip)
      deallocate(filtered_asc)
    endif
    
    if (n_desc > 0) then
      allocate(filtered_desc(n_desc))
      j = 0
      do i = 1, n_hour_group
        if (hour_group(i)%pass_type == -1) then
          j = j + 1
          filtered_desc(j) = hour_group(i)
        endif
      end do
      write(LDT_logunit,*) '[INFO] Processing DESCENDING: ', n_desc, ' files'
      call WSF_ARFS_RESAMPLE_HOURLY(filtered_desc, n_desc, &
                                    WSFopl%WSFoutdir, yyyymmdd, target_hour, n, 'DESC', &
                                    WSFopl%filter_snow_precip)
      deallocate(filtered_desc)
    endif
    
    if (n_asc == 0 .and. n_desc == 0) then
      write(LDT_logunit,*) '[WARN] No valid files for hour ', target_hour
    endif
    
    deallocate(wsf_files, filtered_files, hour_group)
    
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] WSF resampling complete'
    write(LDT_logunit,*) '[INFO] ========================================'
    
  end subroutine LDT_wsf_oplRun

  function file_overlaps_hour(file_info, target_hour) result(overlaps)
    ! Check if file time range overlaps with target hour
    ! File format: t002900_e004059 means 00:29:00 to 00:40:59
    
    implicit none
    type(wsf_file_info), intent(in) :: file_info
    integer, intent(in) :: target_hour
    logical :: overlaps
    
    integer :: start_hour, end_hour
    
    ! Extract hours from time strings (HHMMSS format)
    read(file_info%start_time(1:2), '(I2)') start_hour
    read(file_info%end_time(1:2), '(I2)') end_hour
    
    ! Check if target hour overlaps with [start_hour, end_hour]
    ! Handle wrap-around case (e.g., 23:xx to 00:xx)
    if (end_hour < start_hour) then
      ! Wraps around midnight
      overlaps = (target_hour >= start_hour .or. target_hour <= end_hour)
    else
      ! Normal case
      overlaps = (target_hour >= start_hour .and. target_hour <= end_hour)
    endif
    
  end function file_overlaps_hour

  subroutine parse_wsf_filename_improved(filename, file_info)
    use LDT_logMod, only: LDT_logunit
    
    implicit none
    character(len=*), intent(in) :: filename
    type(wsf_file_info), intent(out) :: file_info

    character(len=LDT_CONST_PATH_LEN) :: basename
    integer :: pos, pos_d, pos_t, pos_e, pos_i
    integer*8 :: file_size_bytes
    
    file_info%filename = filename
    
    ! Get file size
    inquire(file=trim(filename), size=file_size_bytes)
    file_info%file_size = file_size_bytes
    
    ! Extract basename
    pos = index(filename, '/', back=.true.) + 1
    if (pos == 1) pos = 1
    basename = filename(pos:)
    
    ! Find _d (date), _t (start time), _e (end time) by searching
    pos_d = index(basename, '_d')
    pos_t = index(basename, '_t')
    pos_e = index(basename, '_e')
    
    if (pos_d == 0 .or. pos_t == 0 .or. pos_e == 0) then
      write(LDT_logunit,*) '[WARN] Cannot parse _d/_t/_e from: ', trim(basename)
      file_info%date_str = '00000000'
      file_info%start_time = '000000'
      file_info%end_time = '000000'
      file_info%hour_str = '00'
      file_info%start_hour = -1
      file_info%end_hour = -1
      file_info%copy_num = '00'
      file_info%copy_int = 0
      file_info%integrity = '1'
      file_info%integrity_int = 1
      file_info%pass_type = 0
      return
    endif
    
    ! _d + 8 digits = date
    file_info%date_str = basename(pos_d+2:pos_d+9)
    
    ! _t + 6 digits = start time
    file_info%start_time = basename(pos_t+2:pos_t+7)
    file_info%hour_str   = basename(pos_t+2:pos_t+3)
    read(file_info%hour_str, '(I2)') file_info%start_hour
    
    ! _e + 6 digits = end time
    file_info%end_time = basename(pos_e+2:pos_e+7)
    read(basename(pos_e+2:pos_e+3), '(I2)') file_info%end_hour
    
    ! _c + 2 digits = copy number
    pos = index(basename, '_c')
    if (pos > 0) then
      file_info%copy_num = basename(pos+2:pos+3)
      read(file_info%copy_num, '(I2)') file_info%copy_int
    else
      file_info%copy_num = '00'
      file_info%copy_int = 0
    endif
    
    ! _i + 1 digit = integrity (must come after _c)
    pos_i = index(basename, '_i')
    if (pos_i > 0 .and. pos_i > pos) then
      file_info%integrity = basename(pos_i+2:pos_i+2)
      read(file_info%integrity, '(I1)') file_info%integrity_int
    else
      file_info%integrity = '1'
      file_info%integrity_int = 1
    endif
    
    file_info%pass_type = 0
    
  end subroutine parse_wsf_filename_improved
    
  subroutine filter_duplicate_files_improved(all_files, n_files, filtered, n_filtered)
    use LDT_logMod, only: LDT_logunit
    
    implicit none
    integer, intent(in) :: n_files
    type(wsf_file_info), intent(in) :: all_files(n_files)
    type(wsf_file_info), intent(out) :: filtered(n_files)
    integer, intent(out) :: n_filtered
    
    character*50 :: unique_key
    logical :: is_duplicate, should_replace
    integer :: i, j
    
    n_filtered = 0
    
    do i = 1, n_files
      unique_key = trim(all_files(i)%date_str)//'_'// &
                  trim(all_files(i)%start_time)//'_'// &
                  trim(all_files(i)%end_time)
      
      is_duplicate = .false.
      do j = 1, n_filtered
        if (trim(filtered(j)%date_str)//'_'// &
            trim(filtered(j)%start_time)//'_'// &
            trim(filtered(j)%end_time) == unique_key) then
          is_duplicate = .true.
          
          ! Improved selection logic
          should_replace = .false.
          
          ! Rule 1: If c01 has i0 (complete), keep it
          if (filtered(j)%copy_int == 1 .and. filtered(j)%integrity_int == 0) then
            should_replace = .false.
            write(LDT_logunit,*) '[INFO] Keeping c01_i0 for ', trim(unique_key)
            
          ! Rule 2: Prefer complete files (i0) over incomplete (i1)
          else if (all_files(i)%integrity_int < filtered(j)%integrity_int) then
            should_replace = .true.
            write(LDT_logunit,*) '[INFO] Replacing i', filtered(j)%integrity_int, &
                                ' with i', all_files(i)%integrity_int, &
                                ' for time ', trim(unique_key)
                                
          ! Rule 3: If both have same integrity, prefer lower copy number for complete files
          else if (all_files(i)%integrity_int == 0 .and. &
                  filtered(j)%integrity_int == 0) then
            ! Both complete - prefer lower copy number (c01 over c02)
            if (all_files(i)%copy_int < filtered(j)%copy_int) then
              should_replace = .true.
              write(LDT_logunit,*) '[INFO] Replacing c', filtered(j)%copy_num, &
                                  ' with c', all_files(i)%copy_num, &
                                  ' (both i0) for ', trim(unique_key)
            endif
            
          ! Rule 4: If all are incomplete (i1), select largest file size
          else if (all_files(i)%integrity_int == 1 .and. &
                  filtered(j)%integrity_int == 1) then
            if (all_files(i)%file_size > filtered(j)%file_size) then
              should_replace = .true.
              write(LDT_logunit,*) '[INFO] Replacing with larger file size: ', &
                                  all_files(i)%file_size, ' > ', filtered(j)%file_size, &
                                  ' bytes for ', trim(unique_key)
            endif
          endif
          
          if (should_replace) then
            filtered(j) = all_files(i)
          endif
          exit
        endif
      end do
      
      if (.not. is_duplicate) then
        n_filtered = n_filtered + 1
        filtered(n_filtered) = all_files(i)
        write(LDT_logunit,*) '[INFO] Adding file: c', all_files(i)%copy_num, &
                            '_i', all_files(i)%integrity_int, &
                            ', size=', all_files(i)%file_size, ' bytes'
      endif
    end do
    
    ! Final report on selected files
    write(LDT_logunit,*) '[INFO] ========================================='
    write(LDT_logunit,*) '[INFO] Final file selection summary:'
    write(LDT_logunit,*) '[INFO] Total unique time slots: ', n_filtered
    
    do i = 1, n_filtered
      write(LDT_logunit,*) '[INFO] Selected: ', &
                          filtered(i)%start_time, '-', filtered(i)%end_time, &
                          ' c', filtered(i)%copy_num, &
                          ' i', filtered(i)%integrity_int, &
                          ' size=', filtered(i)%file_size
    end do
    write(LDT_logunit,*) '[INFO] ========================================='
    
  end subroutine filter_duplicate_files_improved

  subroutine search_WSF_files(ndir, date_curr, suffix)
    use LDT_logMod, only: LDT_logunit
    
    implicit none
    character (len=*) :: ndir
    character (len=*) :: date_curr
    character(len=*)  :: suffix

    character*8       :: yyyymmdd
    character*10      :: tmp
    character(len=LDT_CONST_PATH_LEN) :: list_files
    character(len=LDT_CONST_PATH_LEN) :: search_pattern

    external :: system

    yyyymmdd = date_curr(1:8)
    
    tmp = trim(suffix)
    
    search_pattern = trim(ndir)//'/*WSFM_01_d'//trim(yyyymmdd)//'*_res_sdr.nc'
    
    list_files = 'ls '//trim(search_pattern)// &
                 ' > WSF_filelist_'//trim(tmp)//'.dat 2>&1'
    
    write(LDT_logunit,*) '[INFO] ========================================='
    write(LDT_logunit,*) '[INFO] Searching for WSF files'
    write(LDT_logunit,*) '[INFO] Date: ', trim(yyyymmdd)
    write(LDT_logunit,*) '[INFO] Search pattern: ', trim(search_pattern)
    write(LDT_logunit,*) '[INFO] ========================================='
    
    call system(trim(list_files))

  end subroutine search_WSF_files

end module LDT_wsf_oplMod
