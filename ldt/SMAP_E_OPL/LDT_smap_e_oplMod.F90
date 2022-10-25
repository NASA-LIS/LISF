!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: LDT_smap_e_oplMod
! 
! !DESCRIPTION: 
! This module handles the run mode plugin for the 
! Operation Enhanced (9-km) SMAP soil moisture
!
! !REVISION HISTORY: 
!  14 Dec 2021: Yonghwan Kwon, Initial Specification
!
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module LDT_smap_e_oplMod

  ! Defaults
  implicit none
  private

  ! Public methods
  public :: LDT_smap_e_oplInit
  public :: LDT_smap_e_oplRun

  ! Public type
  type, public :: smap_e_opl_dec

    character*100        :: L1Bdir, L1Bresampledir, L1Bresampledir_02, SMoutdir 
    character*100        :: LISdir, LISsnowdir
    character*100        :: TAUdir, OMEGAfile, BDfile, &
                            CLAYfile, Hfile, LCfile
    character*100        :: dailystats_ref, dailystats_lis
    character*10         :: date_curr
    integer              :: L1BresampWriteOpt, L1Btype
    integer              :: Teffscale
    integer              :: ntimes,ngrid
    real, allocatable    :: mu_6am_ref(:), mu_6pm_ref(:) !(ngrid)
    real, allocatable    :: sigma_6am_ref(:), sigma_6pm_ref(:) !(ngrid)
    real, allocatable    :: mu_6am_lis(:), mu_6pm_lis(:) !(ngrid)
    real, allocatable    :: sigma_6am_lis(:), sigma_6pm_lis(:) !(ngrid)
    integer, allocatable :: grid_col(:), grid_row(:) !(ngrid)
    real, allocatable    :: ARFS_TBV_COR(:,:)
    real                 :: SD_thold

  end type smap_e_opl_dec

  type(smap_e_opl_dec), public :: SMAPeOPL  

contains

  subroutine LDT_smap_e_oplInit()
  ! Reads Operational Enhanced SMAP-specific entries from ldt.config

    ! Imports
    use ESMF
    use LDT_coreMod, only: LDT_config
    use LDT_logMod, only: LDT_verify

    ! Defaults
    implicit none

    ! Local variables
    character(len=255) :: cfg_entry
    integer :: rc

    ! *** Get former environment variables ***

    ! Get L1Bdir
    cfg_entry = "SMAP_E_OPL valid date (YYYYMMDDHH):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%date_curr, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL soil moisture output directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%SMoutdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL L1B data directory:" 
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%L1Bdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL L1B data type:"    !1: NRT; 2: Historical
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%L1Btype, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL write L1B resampled output:"    !0: off; 1: on
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%L1BresampWriteOpt, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    if(SMAPeOPL%L1BresampWriteOpt.eq.1) then
       cfg_entry = "SMAP_E_OPL L1B resampled output directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%L1Bresampledir, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    endif

    cfg_entry = "SMAP_E_OPL LIS soil temperature directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%LISdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL apply soil temperature bias correction:"  !0: off; 1: on
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%Teffscale, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    if(SMAPeOPL%Teffscale.eq.1) then
       cfg_entry = "SMAP_E_OPL reference Teff daily statistics file:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%dailystats_ref, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")

       cfg_entry = "SMAP_E_OPL LIS Teff daily statistics file:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%dailystats_lis, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    endif

    cfg_entry = "SMAP_E_OPL LIS snow directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%LISsnowdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL snow depth threshold:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%SD_thold, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL TAU directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%TAUdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL OMEGA file:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%OMEGAfile, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL soil bulk density file:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%BDfile, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL soil clay fraction file:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%CLAYfile, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL roughness file:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%Hfile, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "SMAP_E_OPL landcover file:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, SMAPeOPL%LCfile, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

  end subroutine LDT_smap_e_oplInit


  subroutine LDT_smap_e_oplRun(n)
! This calls the actual SMAP_E_OPL driver

! !USES:
    use LDT_logMod
    use LDT_coreMod

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
!EOP

    integer, external       :: LDT_create_subdirs
    integer                 :: i, fi
    integer                 :: ftn, ierr
    character*100           :: fname
    character*100           :: smap_L1B_filename(10)
    character*8             :: yyyymmdd, yyyymmdd_01, yyyymmdd_02
    character*6             :: hhmmss(10)
    character*4             :: yyyy, yyyy_01, yyyy_02
    character*2             :: hh, mm, dd
    character*2             :: hh_01, mm_01, dd_01 
    character*2             :: hh_02, mm_02, dd_02
    character*1             :: Orbit
    integer                 :: yr, mo, da, hr
    integer                 :: yr_pre, mo_pre, da_pre
    integer                 :: yr_02, mo_02, da_02, hr_02
    logical                 :: dir_exists, read_L1Bdata
    real                    :: teff_01(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real                    :: teff_02(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real                    :: SnowDepth(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real                    :: TIMEsec(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real                    :: UTChr(LDT_rc%lnc(n),LDT_rc%lnr(n))
    integer                 :: L1B_dir_len
    integer                 :: doy_pre, doy_curr

  ! Resample SMAP L1B to L1C
    call search_SMAPL1B_files(SMAPeOPL%L1Bdir,SMAPeOPL%date_curr,&
                              SMAPeOPL%L1Btype)

    yyyymmdd = SMAPeOPL%date_curr(1:8)
    yyyy     = SMAPeOPL%date_curr(1:4)
    mm       = SMAPeOPL%date_curr(5:6)
    dd       = SMAPeOPL%date_curr(7:8)
    hh       = SMAPeOPL%date_curr(9:10)

    if(SMAPeOPL%L1BresampWriteOpt.eq.1) then
       SMAPeOPL%L1Bresampledir_02 = trim(SMAPeOPL%L1Bresampledir)//'/'//&
                                    trim(yyyymmdd)//'/'//trim(hh)

       ierr = LDT_create_subdirs(len_trim(SMAPeOPL%L1Bresampledir_02), &
          trim(SMAPeOPL%L1Bresampledir_02))
    endif

    ftn = LDT_getNextUnitNumber()
    open(unit=ftn, file='./SMAP_L1B_filelist.dat',&
         status='old', iostat=ierr)

    fi = 0
    do while (ierr .eq. 0)
       read (ftn, '(a)', iostat=ierr) fname
       if (ierr .ne. 0) then
          exit
       endif
       fi = fi + 1
       smap_L1B_filename(fi) = fname
    enddo
    call LDT_releaseUnitNumber(ftn)

    L1B_dir_len = len_trim(SMAPeOPL%L1Bdir)
    read_L1Bdata = .false.
    if(fi.ge.1) then
       do i=1,fi
          hhmmss(i) = trim(smap_L1B_filename(i)(L1B_dir_len+35:L1B_dir_len+40))
          hhmmss(i+1) = trim(smap_L1B_filename(i+1)(L1B_dir_len+35:L1B_dir_len+40))

          ! use latest version (i.e., highest version number of N**** and/or 00*)
          if(i == fi) then
             write (LDT_logunit,*) '[INFO] Resampling ', trim(smap_L1B_filename(i))
             allocate(SMAPeOPL%ARFS_TBV_COR(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             call SMAPL1BRESAMPLE(smap_L1B_filename(i),SMAPeOPL%L1Bdir,Orbit,TIMEsec)
             write (LDT_logunit,*) '[INFO] Finished resampling ', trim(smap_L1B_filename(i))
             read_L1Bdata = .true. 
          elseif(hhmmss(i) /= hhmmss(i+1)) then
             write (LDT_logunit,*) '[INFO] Resampling ', trim(smap_L1B_filename(i))
             allocate(SMAPeOPL%ARFS_TBV_COR(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             call SMAPL1BRESAMPLE(smap_L1B_filename(i),SMAPeOPL%L1Bdir,Orbit,TIMEsec)
             write (LDT_logunit,*) '[INFO] Finished resampling ', trim(smap_L1B_filename(i))
             read_L1Bdata = .true.
          endif

          if(read_L1Bdata) then
  ! Get effective soil temperature (Teff) from LIS outputs

             ! use LIS outputs from previous day
             read(yyyy,*,iostat=ierr)  yr
             read(mm,*,iostat=ierr)    mo
             read(dd,*,iostat=ierr)    da                                        
             read(hh,*,iostat=ierr)    hr

             yr_pre = yr
             mo_pre = mo
             da_pre = da - 1
             if(da_pre.eq.0) then
                mo_pre = mo - 1

                if(mo_pre.eq.0) then
                   yr_pre = yr - 1
                   mo_pre = 12
                   da_pre = 31
                else
                   if(mo_pre.eq.1.or.&
                    mo_pre.eq.3.or.&
                    mo_pre.eq.5.or.&
                    mo_pre.eq.7.or.&
                    mo_pre.eq.8.or.&
                    mo_pre.eq.10.or.&
                    mo_pre.eq.12) then
                      da_pre = 31
                   elseif(mo_pre.eq.2) then
                      if(mod(yr,4).eq.0) then
                         da_pre = 29
                      else
                         da_pre = 28
                      endif
                   else
                      da_pre = 30
                   endif
                endif
             endif

             write(unit=yyyy_01, fmt='(i4.4)') yr_pre
             write(unit=mm_01, fmt='(i2.2)') mo_pre
             write(unit=dd_01, fmt='(i2.2)') da_pre
             yyyymmdd_01 = trim(yyyy_01)//trim(mm_01)//trim(dd_01)
             hh_01 = hh

             call readLIS_Teff(n,yyyymmdd_01,hh_01,Orbit,teff_01)

             yr_02 = yr_pre
             mo_02 = mo_pre
             da_02 = da_pre
             hr_02 = hr + 1

             if(hr_02.eq.24) then
                hr_02 = 0
                da_02 = da_pre + 1
                
                if(mo_pre.eq.1.or.&
                 mo_pre.eq.3.or.&
                 mo_pre.eq.5.or.&
                 mo_pre.eq.7.or.&
                 mo_pre.eq.8.or.&
                 mo_pre.eq.10.or.&
                 mo_pre.eq.12) then
                   if(da_02.gt.31) then
                      da_02 = 1
                      mo_02 = mo_pre + 1
                   endif
                elseif(mo_pre.eq.2) then
                   if(mod(yr_02,4).eq.0) then
                      if(da_02.gt.29) then
                         da_02 = 1
                         mo_02 = mo_pre + 1
                      endif
                   else
                      if(da_02.gt.28) then
                         da_02 = 1
                         mo_02 = mo_pre + 1
                      endif
                  endif
                else
                   if(da_02.gt.30) then
                      da_02 = 1
                      mo_02 = mo_pre + 1
                   endif
                endif

                if(mo_02.gt.12) then
                   mo_02 = 1
                   yr_02 = yr_pre + 1
                endif
             endif

             write(unit=yyyy_02, fmt='(i4.4)') yr_02
             write(unit=mm_02, fmt='(i2.2)') mo_02
             write(unit=dd_02, fmt='(i2.2)') da_02
             write(unit=hh_02, fmt='(i2.2)') hr_02
             yyyymmdd_02 = trim(yyyy_02)//trim(mm_02)//trim(dd_02)

             call readLIS_Teff(n,yyyymmdd_02,hh_02,Orbit,teff_02)

  ! Scale LIS teff to GEOS teff climatology
             ! get DOY
             call get_doy(mo_pre,da_pre,doy_pre)

             if(SMAPeOPL%Teffscale.eq.1) then
                ! get getattributes
                call getattributes(SMAPeOPL%dailystats_ref,&
                                   SMAPeOPL%ntimes,SMAPeOPL%ngrid)

                ! read 6-yr daily mean and std dev
                allocate(SMAPeOPL%mu_6am_ref(SMAPeOPL%ngrid))
                allocate(SMAPeOPL%mu_6pm_ref(SMAPeOPL%ngrid))
                allocate(SMAPeOPL%sigma_6am_ref(SMAPeOPL%ngrid))
                allocate(SMAPeOPL%sigma_6pm_ref(SMAPeOPL%ngrid))
                allocate(SMAPeOPL%mu_6am_lis(SMAPeOPL%ngrid))
                allocate(SMAPeOPL%mu_6pm_lis(SMAPeOPL%ngrid))
                allocate(SMAPeOPL%sigma_6am_lis(SMAPeOPL%ngrid))
                allocate(SMAPeOPL%sigma_6pm_lis(SMAPeOPL%ngrid))
                allocate(SMAPeOPL%grid_col(SMAPeOPL%ngrid))
                allocate(SMAPeOPL%grid_row(SMAPeOPL%ngrid))

                call read_DailyTeffStats(doy_pre)

                ! scale
                write (LDT_logunit,*) '[INFO] Scaling LIS effective soil temperature'
                call scale_teff(n,Orbit,teff_01,teff_02)
                write (LDT_logunit,*) '[INFO] Finished scaling LIS effective soil temperature'

                deallocate(SMAPeOPL%mu_6am_ref)
                deallocate(SMAPeOPL%mu_6pm_ref)
                deallocate(SMAPeOPL%sigma_6am_ref)
                deallocate(SMAPeOPL%sigma_6pm_ref)
                deallocate(SMAPeOPL%mu_6am_lis)
                deallocate(SMAPeOPL%mu_6pm_lis)
                deallocate(SMAPeOPL%sigma_6am_lis)
                deallocate(SMAPeOPL%sigma_6pm_lis)
                deallocate(SMAPeOPL%grid_col)
                deallocate(SMAPeOPL%grid_row)
             endif

             read_L1Bdata = .false.

  ! Get snow information from LIS outputs
             call readLIS_snow(n,yyyymmdd,hh,SnowDepth)

  ! Retrieve SMAP soil moisture
             ! get DOY
             call get_doy(mo,da,doy_curr)

             ! get UTC
             call get_UTC(n,TIMEsec,UTChr)                         

             ! retrieve
             ierr = LDT_create_subdirs(len_trim(SMAPeOPL%SMoutdir), &
                trim(SMAPeOPL%SMoutdir))
             call ARFSSMRETRIEVAL(smap_L1B_filename(i),teff_01,teff_02,&
                                  SnowDepth,doy_curr,UTChr)

             deallocate(SMAPeOPL%ARFS_TBV_COR)
          endif
       enddo
    endif

  end subroutine LDT_smap_e_oplRun


  subroutine search_SMAPL1B_files(ndir,date_curr,L1Btype)

    implicit none
! !ARGUMENTS:
    character (len=*) :: ndir
    character (len=*) :: date_curr
    integer           :: L1Btype

! !Local variables
    character*8       :: yyyymmdd
    character*2       :: hh
    character*200     :: list_files

    yyyymmdd = date_curr(1:8)
    hh       = date_curr(9:10)

    if(L1Btype.eq.1) then   !NRT
       list_files = 'ls '//trim(ndir)//'/SMAP_L1B_TB_NRT_*'//&
                    trim(yyyymmdd)//'T'//trim(hh) &
                    //"*.h5 > SMAP_L1B_filelist.dat"
    elseif(L1Btype.eq.2) then   !Historical
       list_files = 'ls '//trim(ndir)//'/SMAP_L1B_TB_*'//&
                    trim(yyyymmdd)//'T'//trim(hh) &
                    //"*.h5 > SMAP_L1B_filelist.dat"
    endif

    call system(trim(list_files))

  end subroutine search_SMAPL1B_files

end module LDT_smap_e_oplMod
