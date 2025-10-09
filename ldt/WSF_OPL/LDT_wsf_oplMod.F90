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
    
    write(LDT_logunit,*) '[INFO] Initializing WSF low resolution resampling'
    
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
    
  end subroutine LDT_wsf_oplInit

  subroutine LDT_wsf_oplRun(n)
    ! Main WSF resampling driver
    
    use LDT_logMod
    
    implicit none
    integer, intent(in) :: n
    character*255 :: wsf_filename
    integer :: ftn, ierr
    
    external :: WSF_ARFS_RESAMPLE
    
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] Starting WSF Low Resolution Resampling'
    write(LDT_logunit,*) '[INFO] ========================================'
    
    ! Read file list
    ftn = LDT_getNextUnitNumber()
    open(unit=ftn, file='./WSF_filelist.dat', status='old', iostat=ierr)
    
    if (ierr /= 0) then
      write(LDT_logunit,*)'[ERR] Cannot open WSF_filelist.dat'
      write(LDT_logunit,*)'[ERR] Please create this file with:'
      write(LDT_logunit,*)'[ERR]   ls /path/to/wsf/*.nc > WSF_filelist.dat'
      call LDT_releaseUnitNumber(ftn)
      return
    end if
    
    ! Process each file
    do while (ierr == 0)
      read(ftn, '(a)', iostat=ierr) wsf_filename
      if (ierr /= 0) exit
      
      write(LDT_logunit,*) ''
      write(LDT_logunit,*) '[INFO] ========================================'
      write(LDT_logunit,*) '[INFO] Processing: ', trim(wsf_filename)
      
      call WSF_ARFS_RESAMPLE(wsf_filename, WSFopl%WSFoutdir)
      
      write(LDT_logunit,*) '[INFO] Completed: ', trim(wsf_filename)
      write(LDT_logunit,*) '[INFO] ========================================'
    end do
    
    call LDT_releaseUnitNumber(ftn)
    
    write(LDT_logunit,*) ''
    write(LDT_logunit,*) '[INFO] ========================================'
    write(LDT_logunit,*) '[INFO] All WSF Files Processed Successfully'
    write(LDT_logunit,*) '[INFO] ========================================'
    
  end subroutine LDT_wsf_oplRun

end module LDT_wsf_oplMod