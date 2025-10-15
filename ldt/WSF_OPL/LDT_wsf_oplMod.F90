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
!              Modified to output in AMSR_OPL-style structure with 
!              only 10 channels (no Stokes parameters)
!              Includes file search functionality similar to AMSR_OPL
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
    ! Exactly matching AMSR_OPL structure with 10 channels
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
    write(LDT_logunit,*) '[INFO] Quality filtering: bits 0-4 must be 0'
    write(LDT_logunit,*) '[INFO] Land fraction filter: > 0.8'
    write(LDT_logunit,*) '[INFO] Snow/Precip filtering: ENABLED (matching AMSR_OPL)'
    
  end subroutine LDT_wsf_oplInit

  subroutine LDT_wsf_oplRun(n)
    ! Main WSF resampling driver
    
    use LDT_coreMod, only: LDT_rc
    use LDT_logMod
    
    implicit none
    integer, intent(in) :: n
    ! Make sure it's declared like this:
    character(len=255), allocatable :: wsf_filelist(:)   ! Modern syntax
    character*255 :: fname  ! INCREASED from 100 to 255 for long paths
    integer :: ftn, ierr, fi, i
    character*2 :: tmp
    character*8 :: yyyymmdd
    logical :: file_exists
    
    external :: WSF_ARFS_RESAMPLE
    
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] Starting WSF Low Resolution Resampling'
    write(LDT_logunit,*) '[INFO] Output format: AMSR_OPL-style structure'
    write(LDT_logunit,*) '[INFO] Snow/Precip filtering: ENABLED'
    write(LDT_logunit,*) '[INFO] ========================================'
    
    ! Extract date components
    yyyymmdd = WSFopl%date_curr(1:8)
    
    ! Search for WSF files matching the date
    call search_WSF_files(WSFopl%WSFdir, WSFopl%date_curr, &
                         WSFopl%WSFfilelistSuffixNumber)
    
    ! Allocate output arrays with dimensions matching LDT grid
    write(LDT_logunit,*) '[INFO] Allocating output arrays...'
    write(LDT_logunit,*) '[INFO] Grid dimensions: ', LDT_rc%lnc(n), ' x ', LDT_rc%lnr(n)
    
    ! Allocate all 10 channel arrays (V and H only, no Stokes)
    allocate(WSFopl%ARFS_TB_10V(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_TB_10H(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_TB_18V(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_TB_18H(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_TB_23V(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_TB_23H(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_TB_36V(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_TB_36H(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_TB_89V(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_TB_89H(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    
    ! Allocate auxiliary arrays
    allocate(WSFopl%ARFS_LAND_FRAC(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_QUALITY_FLAG(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_SNOW(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_PRECIP(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_SAMPLE_V(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_SAMPLE_H(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    
    ! Initialize arrays to missing values
    WSFopl%ARFS_TB_10V = -9999.0
    WSFopl%ARFS_TB_10H = -9999.0
    WSFopl%ARFS_TB_18V = -9999.0
    WSFopl%ARFS_TB_18H = -9999.0
    WSFopl%ARFS_TB_23V = -9999.0
    WSFopl%ARFS_TB_23H = -9999.0
    WSFopl%ARFS_TB_36V = -9999.0
    WSFopl%ARFS_TB_36H = -9999.0
    WSFopl%ARFS_TB_89V = -9999.0
    WSFopl%ARFS_TB_89H = -9999.0
    WSFopl%ARFS_LAND_FRAC = -9999.0
    WSFopl%ARFS_QUALITY_FLAG = 0
    WSFopl%ARFS_SNOW = 0
    WSFopl%ARFS_PRECIP = 0
    WSFopl%ARFS_SAMPLE_V = 0
    WSFopl%ARFS_SAMPLE_H = 0
    
    ! Read file list generated by search_WSF_files
    allocate(wsf_filelist(1000))
    write (tmp,'(I2.2)') WSFopl%WSFfilelistSuffixNumber
    ftn = LDT_getNextUnitNumber()
    open(unit=ftn, file='./WSF_filelist_'//trim(tmp)//'.dat', &
         status='old', iostat=ierr)
    
    if (ierr /= 0) then
      write(LDT_logunit,*)'[ERR] Cannot open WSF_filelist_'//trim(tmp)//'.dat'
      write(LDT_logunit,*)'[ERR] File search may have failed'
      call LDT_releaseUnitNumber(ftn)
      return
    end if
    
    ! Read all WSF filenames into array
    fi = 0
    do while (ierr == 0)
      read(ftn, '(a)', iostat=ierr) fname
      if (ierr /= 0) exit
      
      ! Skip empty lines
      if (len_trim(fname) == 0) cycle
      
      ! Skip error messages from ls command
      if (index(fname, 'No such file') > 0 .or. &
          index(fname, 'cannot access') > 0) then
        write(LDT_logunit,*) '[WARN] ls error: ', trim(fname)
        cycle
      endif
      
      fi = fi + 1
      if (fi <= 1000) then
        wsf_filelist(fi) = fname
        write(LDT_logunit,*) '[INFO] File ', fi, ': ', trim(fname)
        
        ! Check if file actually exists and is readable
        inquire(file=trim(fname), exist=file_exists)
        if (.not. file_exists) then
          write(LDT_logunit,*) '[ERR] File does not exist: ', trim(fname)
          fi = fi - 1  ! Don't count this file
        endif
      endif
    end do
    call LDT_releaseUnitNumber(ftn)
    
    write(LDT_logunit,*) '[INFO] Found ', fi, ' valid WSF files to process'
    
    ! NEW CODE - Process one file at a time (like AMSR_OPL)
    if (fi >= 1) then
      write(LDT_logunit,*) '[INFO] Processing ', fi, ' WSF files individually'
      
      ! Loop through files and process one at a time
      do i = 1, fi
        write(LDT_logunit,*) '[INFO] Resampling file ', i, ' of ', fi
        write(LDT_logunit,*) '[INFO] File: ', trim(wsf_filelist(i))
        
        ! Call WSF_ARFS_RESAMPLE with SINGLE filename
        call WSF_ARFS_RESAMPLE(wsf_filelist(i), WSFopl%WSFoutdir, n)
        
        write(LDT_logunit,*) '[INFO] Finished processing file ', i
      end do
      
      write(LDT_logunit,*) '[INFO] All ', fi, ' files processed successfully'
    else
      write(LDT_logunit,*) '[WARN] No WSF files found for date: ', WSFopl%date_curr
    endif
    
    deallocate(wsf_filelist)
    
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] WSF resampling complete'
    write(LDT_logunit,*) '[INFO] Snow/Precip filtering applied successfully'
    write(LDT_logunit,*) '[INFO] ========================================'
    
  end subroutine LDT_wsf_oplRun

  subroutine search_WSF_files(ndir, date_curr, suffix)
    ! Searches for WSF NetCDF files matching the current date
    ! WSF filename format: WSFM_01_dYYYYMMDD_tHHMMSS_eHHMMSS_gCOOK-B_r05899_c02_i0_v0522_res_sdr.nc
    ! Example: WSFM_01_d20250601_t000900_e002059_gCOOK-B_r05899_c02_i0_v0522_res_sdr.nc
    
    use LDT_logMod, only: LDT_logunit
    
    implicit none
    ! ARGUMENTS:
    character (len=*) :: ndir
    character (len=*) :: date_curr
    integer           :: suffix

    ! Local variables
    character*8       :: yyyymmdd
    character*2       :: tmp
    character*255     :: list_files      ! INCREASED from 200 to 255
    character*255     :: search_pattern  ! INCREASED from 200 to 255

    yyyymmdd = date_curr(1:8)
    
    write (tmp,'(I2.2)') suffix
    
    ! Build search pattern
    ! Note: Using wildcard W*M to match both WSFM and WSFM naming conventions
    search_pattern = trim(ndir)//'/WSFM_01_d'//trim(yyyymmdd)//'*_res_sdr.nc'
    
    ! Search for WSF files with the date pattern
    list_files = 'ls '//trim(search_pattern)// &
                 ' > WSF_filelist_'//trim(tmp)//'.dat 2>&1'
    
    write(LDT_logunit,*) '[INFO] ========================================='
    write(LDT_logunit,*) '[INFO] Searching for WSF files'
    write(LDT_logunit,*) '[INFO] Date (YYYYMMDD): ', trim(yyyymmdd)
    write(LDT_logunit,*) '[INFO] Search pattern: ', trim(search_pattern)
    write(LDT_logunit,*) '[INFO] Output file: WSF_filelist_'//trim(tmp)//'.dat'
    write(LDT_logunit,*) '[INFO] Executing: ', trim(list_files)
    write(LDT_logunit,*) '[INFO] ========================================='
    
    call system(trim(list_files))

  end subroutine search_WSF_files

end module LDT_wsf_oplMod