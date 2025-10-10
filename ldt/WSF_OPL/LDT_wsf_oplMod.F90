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
    write(LDT_logunit,*) '[INFO] ========================================='
    
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
    
    write(LDT_logunit,*) '[INFO] WSF input directory: ', trim(WSFopl%WSFdir)
    write(LDT_logunit,*) '[INFO] WSF output directory: ', trim(WSFopl%WSFoutdir)
    write(LDT_logunit,*) '[INFO] Using low resolution (30 km) data only'
    write(LDT_logunit,*) '[INFO] Quality filtering: bits 0-4 must be 0'
    write(LDT_logunit,*) '[INFO] Land fraction filter: > 0.8'
    
  end subroutine LDT_wsf_oplInit

  subroutine LDT_wsf_oplRun(n)
    ! Main WSF resampling driver
    
    use LDT_coreMod, only: LDT_rc
    use LDT_logMod
    
    implicit none
    integer, intent(in) :: n
    character*255 :: wsf_filename
    integer :: ftn, ierr, fi
    character*255, allocatable :: wsf_filelist(:)
    
    external :: WSF_ARFS_RESAMPLE
    
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] Starting WSF Low Resolution Resampling'
    write(LDT_logunit,*) '[INFO] Output format: AMSR_OPL-style structure'
    write(LDT_logunit,*) '[INFO] ========================================'
    
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
    
    allocate(WSFopl%ARFS_LAND_FRAC(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    allocate(WSFopl%ARFS_QUALITY_FLAG(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    
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
    
    ! Read file list
    allocate(wsf_filelist(1000))
    ftn = LDT_getNextUnitNumber()
    open(unit=ftn, file='./WSF_filelist.dat', status='old', iostat=ierr)
    
    if (ierr /= 0) then
      write(LDT_logunit,*)'[ERR] Cannot open WSF_filelist.dat'
      write(LDT_logunit,*)'[ERR] Please create this file with:'
      write(LDT_logunit,*)'[ERR]   ls /path/to/wsf/*.nc > WSF_filelist.dat'
      call LDT_releaseUnitNumber(ftn)
      return
    end if
    
    ! Read all filenames
    fi = 0
    do while (ierr == 0)
      read(ftn, '(a)', iostat=ierr) wsf_filename
      if (ierr /= 0) exit
      fi = fi + 1
      wsf_filelist(fi) = wsf_filename
    end do
    call LDT_releaseUnitNumber(ftn)
    
    write(LDT_logunit,*) '[INFO] Found ', fi, ' WSF files to process'
    
    ! Process all files
    call WSF_ARFS_RESAMPLE(wsf_filelist, fi, WSFopl%WSFoutdir, n)
    
    deallocate(wsf_filelist)
    
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] WSF resampling complete'
    write(LDT_logunit,*) '[INFO] ========================================'
    
  end subroutine LDT_wsf_oplRun

end module LDT_wsf_oplMod