!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: LDT_wsf_oplMod
!
! DESCRIPTION: LDT module for WSF low resolution resampling
!              MODIFIED to handle:
!              - Multiple files per hour with duplicate filtering (highest copy number)
!              - Hourly grouping and stitching
!              - Mean averaging for overlapping observations
!
!-------------------------------------------------------------------------

module LDT_wsf_oplMod

  implicit none
  private

  public :: LDT_wsf_oplInit
  public :: LDT_wsf_oplRun

  type, public :: wsf_opl_dec
    character*100 :: WSFdir          ! Input directory
    character*100 :: WSFoutdir       ! Output directory
    character*10  :: date_curr       ! Current date (YYYYMMDDHH)
    integer       :: WSFfilelistSuffixNumber  ! File list suffix number
    
    ! Output arrays in AMSR_OPL style - ONLY V and H polarizations
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
    
    ! Auxiliary data
    real*4, allocatable :: ARFS_LAND_FRAC(:,:)
    integer*1, allocatable :: ARFS_QUALITY_FLAG(:,:)
    integer*4, allocatable :: ARFS_SNOW(:,:)
    integer*4, allocatable :: ARFS_PRECIP(:,:)
    integer*4, allocatable :: ARFS_SAMPLE_V(:,:)
    integer*4, allocatable :: ARFS_SAMPLE_H(:,:)
    
  end type wsf_opl_dec

  type(wsf_opl_dec), public :: WSFopl

  ! File information structure for grouping
  type :: wsf_file_info
    character*255 :: filename
    character*8   :: date_str     ! YYYYMMDD
    character*6   :: start_time   ! HHMMSS
    character*6   :: end_time     ! HHMMSS
    character*2   :: copy_num     ! c00, c01, c02, etc.
    integer       :: copy_int     ! Integer copy number for sorting
    character*2   :: hour_str     ! HH for grouping
  end type wsf_file_info

contains

  subroutine LDT_wsf_oplInit()
    ! Reads WSF-specific entries from ldt.config
    
    use ESMF
    use LDT_coreMod, only: LDT_config
    use LDT_logMod, only: LDT_logunit, LDT_verify
    
    implicit none
    integer :: rc
    character(len=255) :: cfg_entry
    
    write(LDT_logunit,*) '[INFO] ========================================='
    write(LDT_logunit,*) '[INFO] Initializing WSF low resolution resampling'
    write(LDT_logunit,*) '[INFO] Output structure: AMSR_OPL-style (10 channels)'
    write(LDT_logunit,*) '[INFO] WITH hourly grouping and duplicate filtering'
    write(LDT_logunit,*) '[INFO] WITH Snow/Precip Filtering'
    write(LDT_logunit,*) '[INFO] ========================================='
    
    ! Read current date
    cfg_entry = "WSF valid date (YYYYMMDDHH):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, WSFopl%date_curr, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    
    ! Read input directory
    cfg_entry = "WSF input directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, WSFopl%WSFdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    
    ! Read output directory
    cfg_entry = "WSF output directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, WSFopl%WSFoutdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    
    ! Read file list suffix number
    cfg_entry = "WSF filelist suffix number:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, WSFopl%WSFfilelistSuffixNumber, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    
    write(LDT_logunit,*) '[INFO] WSF valid date: ', trim(WSFopl%date_curr)
    write(LDT_logunit,*) '[INFO] WSF input directory: ', trim(WSFopl%WSFdir)
    write(LDT_logunit,*) '[INFO] WSF output directory: ', trim(WSFopl%WSFoutdir)
    write(LDT_logunit,*) '[INFO] WSF filelist suffix: ', WSFopl%WSFfilelistSuffixNumber
    write(LDT_logunit,*) '[INFO] Using low resolution (30 km) data only'
    write(LDT_logunit,*) '[INFO] Duplicate handling: highest copy number selected'
    write(LDT_logunit,*) '[INFO] Hourly grouping: ENABLED'
    
  end subroutine LDT_wsf_oplInit

  subroutine LDT_wsf_oplRun(n)
    ! Main WSF resampling driver with hourly grouping and duplicate filtering
    
    use LDT_coreMod, only: LDT_rc
    use LDT_logMod
    
    implicit none
    integer, intent(in) :: n
    
    type(wsf_file_info), allocatable :: wsf_files(:)
    type(wsf_file_info), allocatable :: filtered_files(:)
    type(wsf_file_info), allocatable :: hour_group(:)
    character(len=255) :: fname
    character*2 :: tmp, current_hour
    character*8 :: yyyymmdd
    integer :: ftn, ierr, fi, i, j, k
    integer :: n_filtered, n_hour_group
    logical :: file_exists
    
    external :: WSF_ARFS_RESAMPLE_HOURLY
    
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] Starting WSF Low Resolution Resampling'
    write(LDT_logunit,*) '[INFO] With hourly grouping and duplicate filtering'
    write(LDT_logunit,*) '[INFO] ========================================'
    
    ! Extract date components
    yyyymmdd = WSFopl%date_curr(1:8)
    
    ! Search for files
    write(tmp,'(I2.2)') WSFopl%WSFfilelistSuffixNumber
    call search_WSF_files(WSFopl%WSFdir, WSFopl%date_curr, WSFopl%WSFfilelistSuffixNumber)
    
    ! Read and parse all filenames
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
      
      ! Skip empty lines and errors
      if (len_trim(fname) == 0) cycle
      if (index(fname, 'No such file') > 0 .or. &
          index(fname, 'cannot access') > 0) then
        cycle
      endif
      
      ! Check if file exists
      inquire(file=trim(fname), exist=file_exists)
      if (.not. file_exists) then
        write(LDT_logunit,*) '[WARN] File does not exist: ', trim(fname)
        cycle
      endif
      
      fi = fi + 1
      if (fi <= 1000) then
        ! Parse filename components
        call parse_wsf_filename(fname, wsf_files(fi))
      endif
    end do
    call LDT_releaseUnitNumber(ftn)
    
    write(LDT_logunit,*) '[INFO] Found ', fi, ' total WSF files'
    
    if (fi == 0) then
      write(LDT_logunit,*) '[WARN] No WSF files found for date: ', WSFopl%date_curr
      deallocate(wsf_files)
      return
    endif
    
    ! Filter duplicates (keep highest copy number)
    allocate(filtered_files(fi))
    call filter_duplicate_files(wsf_files(1:fi), fi, filtered_files, n_filtered)
    
    write(LDT_logunit,*) '[INFO] After duplicate filtering: ', n_filtered, ' files'
    
    ! Process files grouped by hour
    allocate(hour_group(100))  ! Max 100 files per hour
    
    ! Sort filtered files by hour
    call sort_files_by_hour(filtered_files(1:n_filtered), n_filtered)
    
    ! Process each hour group
    i = 1
    do while (i <= n_filtered)
      current_hour = filtered_files(i)%hour_str
      n_hour_group = 0
      
      ! Collect all files for this hour
      do j = i, n_filtered
        if (filtered_files(j)%hour_str == current_hour) then
          n_hour_group = n_hour_group + 1
          hour_group(n_hour_group) = filtered_files(j)
        else
          exit
        endif
      end do
      
      write(LDT_logunit,*) '[INFO] ========================================'
      write(LDT_logunit,*) '[INFO] Processing hour: ', current_hour, 'H'
      write(LDT_logunit,*) '[INFO] Number of files: ', n_hour_group
      write(LDT_logunit,*) '[INFO] ========================================'
      
      ! Process this hour group (stitching multiple files)
      call WSF_ARFS_RESAMPLE_HOURLY(hour_group(1:n_hour_group), n_hour_group, &
                                    WSFopl%WSFoutdir, yyyymmdd, current_hour, n)
      
      i = j
    end do
    
    deallocate(wsf_files)
    deallocate(filtered_files)
    deallocate(hour_group)
    
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] WSF resampling complete'
    write(LDT_logunit,*) '[INFO] All hourly groups processed'
    write(LDT_logunit,*) '[INFO] ========================================'
    
  end subroutine LDT_wsf_oplRun

  subroutine parse_wsf_filename(filename, file_info)
    ! Parse WSF filename to extract components
    ! Format: WSFM_01_d20250601_t002900_e004059_gREEF-A_r05898_c02_i0_v0522_res_sdr.nc
    
    use LDT_logMod, only: LDT_logunit
    
    implicit none
    character(len=*), intent(in) :: filename
    type(wsf_file_info), intent(out) :: file_info
    
    character*255 :: basename
    integer :: pos, i
    
    ! Store full filename
    file_info%filename = filename
    
    ! Extract basename
    pos = index(filename, '/', back=.true.) + 1
    if (pos == 1) pos = 1
    basename = filename(pos:)
    
    ! Extract date (d20250601 -> positions 9-16)
    file_info%date_str = basename(9:16)
    
    ! Extract start time (t002900 -> positions 18-23)
    file_info%start_time = basename(18:23)
    
    ! Extract hour from start time
    file_info%hour_str = basename(18:19)
    
    ! Extract end time (e004059 -> positions 25-30)
    file_info%end_time = basename(25:30)
    
    ! Extract copy number (find _c## pattern)
    pos = index(basename, '_c')
    if (pos > 0) then
      file_info%copy_num = basename(pos+2:pos+3)
      read(file_info%copy_num, '(I2)') file_info%copy_int
    else
      file_info%copy_num = '00'
      file_info%copy_int = 0
    endif
    
  end subroutine parse_wsf_filename

  subroutine filter_duplicate_files(all_files, n_files, filtered, n_filtered)
    ! Filter duplicate files, keeping only highest copy number
    
    use LDT_logMod, only: LDT_logunit
    
    implicit none
    integer, intent(in) :: n_files
    type(wsf_file_info), intent(in) :: all_files(n_files)
    type(wsf_file_info), intent(out) :: filtered(n_files)
    integer, intent(out) :: n_filtered
    
    character*50 :: unique_key
    logical :: is_duplicate
    integer :: i, j, highest_copy
    
    n_filtered = 0
    
    do i = 1, n_files
      ! Create unique key (without copy number)
      unique_key = trim(all_files(i)%date_str)//'_'// &
                  trim(all_files(i)%start_time)//'_'// &
                  trim(all_files(i)%end_time)
      
      ! Check if we already have this time period
      is_duplicate = .false.
      do j = 1, n_filtered
        if (trim(filtered(j)%date_str)//'_'// &
            trim(filtered(j)%start_time)//'_'// &
            trim(filtered(j)%end_time) == unique_key) then
          is_duplicate = .true.
          ! Keep the one with higher copy number
          if (all_files(i)%copy_int > filtered(j)%copy_int) then
            write(LDT_logunit,*) '[INFO] Replacing c', filtered(j)%copy_num, &
                                ' with c', all_files(i)%copy_num, &
                                ' for time ', trim(unique_key)
            filtered(j) = all_files(i)
          endif
          exit
        endif
      end do
      
      if (.not. is_duplicate) then
        n_filtered = n_filtered + 1
        filtered(n_filtered) = all_files(i)
      endif
    end do
    
    write(LDT_logunit,*) '[INFO] Filtered from ', n_files, ' to ', n_filtered, ' files'
    
  end subroutine filter_duplicate_files

  subroutine sort_files_by_hour(files, n_files)
    ! Sort files by hour for grouping
    
    implicit none
    integer, intent(in) :: n_files
    type(wsf_file_info), intent(inout) :: files(n_files)
    
    type(wsf_file_info) :: temp
    integer :: i, j
    
    ! Simple bubble sort by hour and then by start time
    do i = 1, n_files-1
      do j = i+1, n_files
        if (files(i)%hour_str > files(j)%hour_str .or. &
            (files(i)%hour_str == files(j)%hour_str .and. &
             files(i)%start_time > files(j)%start_time)) then
          temp = files(i)
          files(i) = files(j)
          files(j) = temp
        endif
      end do
    end do
    
  end subroutine sort_files_by_hour

  subroutine search_WSF_files(ndir, date_curr, suffix)
    ! Searches for WSF NetCDF files matching the current date
    
    use LDT_logMod, only: LDT_logunit
    
    implicit none
    character (len=*) :: ndir
    character (len=*) :: date_curr
    integer           :: suffix

    ! Local variables
    character*8       :: yyyymmdd
    character*2       :: tmp
    character*255     :: list_files
    character*255     :: search_pattern

    yyyymmdd = date_curr(1:8)
    
    write (tmp,'(I2.2)') suffix
    
    ! Build search pattern - get ALL files for the date
    search_pattern = trim(ndir)//'/WSFM_01_d'//trim(yyyymmdd)//'*_res_sdr.nc'
    
    ! Search for WSF files with the date pattern
    list_files = 'ls '//trim(search_pattern)// &
                 ' > WSF_filelist_'//trim(tmp)//'.dat 2>&1'
    
    write(LDT_logunit,*) '[INFO] ========================================='
    write(LDT_logunit,*) '[INFO] Searching for WSF files'
    write(LDT_logunit,*) '[INFO] Date (YYYYMMDD): ', trim(yyyymmdd)
    write(LDT_logunit,*) '[INFO] Search pattern: ', trim(search_pattern)
    write(LDT_logunit,*) '[INFO] Output file: WSF_filelist_'//trim(tmp)//'.dat'
    write(LDT_logunit,*) '[INFO] ========================================='
    
    call system(trim(list_files))

  end subroutine search_WSF_files

end module LDT_wsf_oplMod