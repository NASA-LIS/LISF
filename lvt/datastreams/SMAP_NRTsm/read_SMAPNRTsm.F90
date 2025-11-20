!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.6
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
!
! !ROUTINE: read_SMAPNRTsm
! \label{read_SMAPNRTsm}
!
! !INTERFACE:
subroutine read_SMAPNRTsm(source)
  !
  ! !USES:
  use ESMF
  use LVT_constantsMod, only: LVT_CONST_PATH_LEN
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SMAPNRTsm_Mod, only : SMAPNRT_smobs

  implicit none
  !
  ! !INPUT PARAMETERS:
  integer,   intent(in)       :: source
  !
  ! !OUTPUT PARAMETERS:
  !
  ! !DESCRIPTION:
  !
  ! This subroutine provides the data reader for the SMAP L2
  ! NRT soil moisture retrieval product.
  !
  ! !FILES USED:
  !
  ! !REVISION HISTORY:
  !  02 July 2025: Mahdi Navari, Initial Specification
  !  16 July 2025: Mahdi Navari, qa added
  !EOP

  logical           :: alarmcheck
  real              :: smc(LVT_rc%lnc, LVT_rc%lnr)
  character*8       :: yyyymmdd
  character*4       :: yyyy
  character*2       :: mm,dd,hh
  character(LVT_CONST_PATH_LEN) :: list_files
  integer           :: ftn,ierr
  character(LVT_CONST_PATH_LEN) :: fname
  integer           :: mn_ind
  integer           :: mn, ss
  real              :: gmt
  integer           :: c,r
  integer           :: doy
  real*8            :: timenow
  real              :: smobs(LVT_rc%lnc,LVT_rc%lnr)

  external :: read_SMAPNRTL2sm_data

  smc = LVT_rc%udef
  smobs = LVT_rc%udef

  timenow = float(LVT_rc%dhr(source))*3600 +&
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)

  alarmcheck = (mod(timenow, 3600.0).eq.0)

  if (SMAPNRT_smobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then

     LVT_rc%resetFlag(source) = .false.

     SMAPNRT_smobs(source)%startflag = .false.

     write(yyyymmdd,'(i4.4,2i2.2)') LVT_rc%dyr(source), &
          LVT_rc%dmo(source), &
          LVT_rc%dda(source)
     write(yyyy,'(i4.4)') LVT_rc%dyr(source)
     write(mm,'(i2.2)') LVT_rc%dmo(source)
     write(dd,'(i2.2)') LVT_rc%dda(source)
     write(hh,'(i2.2)') LVT_rc%dhr(source)

     ! SMAP_L2_SM_P_NRT_45372_A_20230730T212058_N17701_001.h5
     list_files = 'ls ' // trim(SMAPNRT_smobs(source)%odir) // &
          '/SMAP_L2_*' // trim(yyyymmdd) // 'T' // trim(hh) &
          // "*.h5 > SMAPsm/SMAP_filelist.dat"

     call system(trim(list_files))

     ftn = LVT_getNextUnitNumber()
     open(ftn, file="./SMAPsm/SMAP_filelist.dat", &
          status='old', iostat=ierr)

     ! if multiple files for the same time and orbits are present, the
     ! latest one will overwrite older ones, though multiple (redundant)
     ! reads occur.
     ! This assumes that the 'ls command' will list the files in that
     ! order.
     do while (ierr.eq.0)
        read(ftn,'(a)',iostat=ierr) fname
        if (ierr.ne.0) then
           exit
        endif
        mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh)) + 11
        read(fname(mn_ind:mn_ind+1),'(i2.2)') mn
        ss = 0
        call LVT_tick(timenow,doy,gmt,LVT_rc%dyr(source), &
             LVT_rc%dmo(source), &
             LVT_rc%dda(source), &
             LVT_rc%dhr(source), mn, ss, 0)

        write(LVT_logunit,*) '[INFO] reading ',trim(fname)

        call read_SMAPNRTL2sm_data(source, fname, smobs, timenow)

     enddo
     call LVT_releaseUnitNumber(ftn)

     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           smc(c,r) = smobs(c,r)
        enddo
     enddo

  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source,&
       smc, vlevel=1, units="m3/m3")

end subroutine read_SMAPNRTsm

!BOP
!
! !ROUTINE: read_SMAPNRTL2sm_data
! \label{read_SMAPNRTL2sm_data}
!
! !INTERFACE:
subroutine read_SMAPNRTL2sm_data(source, fname, smobs_inp, time)
  !
  ! !USES:

#if (defined USE_HDF5)
  use hdf5
#endif
  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SMAPNRTsm_Mod, only : SMAPNRT_smobs

  implicit none
  !
  ! !INPUT PARAMETERS:
  !
  integer                  :: source
  character (len=*)        :: fname
  real                     :: smobs_inp(LVT_rc%lnc,LVT_rc%lnr)
  real*8                   :: time

  ! !OUTPUT PARAMETERS:
  !
  !
  ! !DESCRIPTION:
  !
  !
  !EOP

#if (defined USE_HDF5)

  character*100, parameter :: sm_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100, parameter :: sm_field_name = "soil_moisture"
  character*100, parameter :: sm_qa_name = "retrieval_qual_flag"

  integer(hsize_t), dimension(1) :: dims
  integer(hsize_t), dimension(1) :: maxdims
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: dspace_id
  integer(hid_t)                 :: row_id, col_id
  integer(hid_t)                 :: sm_gr_id, sm_field_id, sm_qa_id
  real,             allocatable  :: sm_field(:)
  integer,          allocatable  :: sm_qa(:)
  integer,          allocatable  :: ease_row(:)
  integer,          allocatable  :: ease_col(:)
  integer                        :: c, r, t
  logical*1                      :: &
       sm_data_b(SMAPNRT_smobs(source)%nc*SMAPNRT_smobs(source)%nr)
  logical*1                      :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: &
       sm_data(SMAPNRT_smobs(source)%nc*SMAPNRT_smobs(source)%nr)
  real                           :: smobs_ip(LVT_rc%lnc*LVT_rc%lnr)
  integer                        :: status, iret

  external :: neighbor_interp

  call h5open_f(status)
  call LVT_verify(status, 'Error opening HDF fortran interface')

  call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, status)
  call LVT_verify(status, 'Error opening SMAP L2 file ')

  call h5gopen_f(file_id, sm_gr_name, sm_gr_id, status)
  call LVT_verify(status, 'Error opening SM group in SMAP L2 file')

  call h5dopen_f(sm_gr_id, sm_field_name, sm_field_id, status)
  call LVT_verify(status, 'Error opening SM field in SMAP L2 file')

  call h5dopen_f(sm_gr_id, "EASE_row_index", row_id, status)
  call LVT_verify(status, 'Error opening row index field in SMAP L2 file')

  call h5dopen_f(sm_gr_id, "EASE_column_index", col_id, status)
  call LVT_verify(status, 'Error opening column index field in SMAP L2 file')

  call h5dopen_f(sm_gr_id, sm_qa_name, sm_qa_id, status)
  call LVT_verify(status, 'Error opening QA field in SMAP L2 file')

  call h5dget_space_f(sm_field_id, dspace_id, status)
  call LVT_verify(status, 'Error in h5dget_space_f: readSMAP L2Obs')

  ! Size of the arrays
  ! This routine returns -1 on failure, rank on success.
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, status)
  if (status.eq.-1) then
     call LVT_verify(status, &
          'Error in h5sget_simple_extent_dims_f: readSMAP L2Obs')
  endif

  allocate(sm_field(maxdims(1)))
  allocate(sm_qa(maxdims(1)))
  allocate(ease_row(maxdims(1)))
  allocate(ease_col(maxdims(1)))

  call h5dread_f(row_id, H5T_NATIVE_INTEGER, ease_row, dims, status)
  call LVT_verify(status, 'Error extracting row index from SMAP L2 file')

  call h5dread_f(col_id, H5T_NATIVE_INTEGER, ease_col, dims, status)
  call LVT_verify(status, 'Error extracting col index from SMAP L2 file')

  call h5dread_f(sm_field_id, H5T_NATIVE_REAL, sm_field, dims, status)
  call LVT_verify(status, 'Error extracting SM field from SMAP L2 file')

  call h5dread_f(sm_qa_id, H5T_NATIVE_INTEGER, sm_qa, dims, status)
  call LVT_verify(status, 'Error extracting SM field from SMAP L2 file')

  call h5dclose_f(sm_qa_id, status)
  call LVT_verify(status, 'Error in H5DCLOSE call')

  call h5dclose_f(row_id, status)
  call LVT_verify(status, 'Error in H5DCLOSE call')

  call h5dclose_f(col_id, status)
  call LVT_verify(status, 'Error in H5DCLOSE call')

  call h5dclose_f(sm_field_id, status)
  call LVT_verify(status, 'Error in H5DCLOSE call')

  call h5gclose_f(sm_gr_id, status)
  call LVT_verify(status, 'Error in H5GCLOSE call')

  call h5fclose_f(file_id, status)
  call LVT_verify(status, 'Error in H5FCLOSE call')

  call h5close_f(status)
  call LVT_verify(status, 'Error in H5CLOSE call')

  sm_data = LVT_rc%udef
  sm_data_b = .false.

  !grid the data in EASE projection
  do t=1,maxdims(1)
     if (ibits(sm_qa(t),0,1).eq.0) then
        if (ease_col(t).gt.0.and.ease_row(t).gt.0) then
           sm_data(ease_col(t) + &
                (ease_row(t)-1)*SMAPNRT_smobs(source)%nc) = sm_field(t)
           if (sm_field(t).ne.-9999.0) then
              sm_data_b(ease_col(t) + &
                   (ease_row(t)-1)*SMAPNRT_smobs(source)%nc) = .true.
           endif
        endif
     endif
  enddo

  t = 1

  !-----------------------------------------------------------------------
  ! Interpolate to the LVT running domain
  !-----------------------------------------------------------------------
  call neighbor_interp(LVT_rc%gridDesc, sm_data_b, sm_data, &
       smobs_b_ip, smobs_ip, &
       SMAPNRT_smobs(source)%nc*SMAPNRT_smobs(source)%nr, &
       LVT_rc%lnc*LVT_rc%lnr, &
       SMAPNRT_smobs(source)%rlat2, SMAPNRT_smobs(source)%rlon2, &
       SMAPNRT_smobs(source)%n112, LVT_rc%udef, iret)

  deallocate(sm_field)
  deallocate(sm_qa)
  deallocate(ease_row)
  deallocate(ease_col)

  !overwrite the input data
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if (smobs_ip(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then
           smobs_inp(c,r) = &
                smobs_ip(c+(r-1)*LVT_rc%lnc)

           SMAPNRT_smobs(source)%smtime(c,r) = time
        endif
     enddo
  enddo

#endif

end subroutine read_SMAPNRTL2sm_data

