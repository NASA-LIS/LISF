!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_SMAPsmobs
! \label{read_SMAPsmobs}
!
! !REVISION HISTORY:
!  21 Sep 2018   Sujay Kumar;   Initial Specification
!
! !INTERFACE: 
subroutine read_SMAPsmobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_fileIOMod
  use SMAPsm_obsMod
  use map_utils

#if (defined USE_NETCDF3 || defined USE_NETCDF4)  
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_Space] Objective Space
!  \end{description}
!
!EOP
  integer                  :: n
  real,    pointer         :: smc(:)
  type(ESMF_Field)         :: smcField
  character*100            :: obsdir
  integer                  :: status
  integer                  :: c,r
  real                     :: lon
  real                     :: lhour
  real                     :: gmt
  integer                  :: zone
  logical                  :: data_update
  logical                  :: alarmCheck
  logical                  :: file_exists
  integer                  :: grid_index
  character*100            :: fname
  real                     :: dt
  real, allocatable        :: smobs_D(:)
  real, allocatable        :: smobs_A(:)

  n = 1

  call ESMF_AttributeGet(Obj_Space,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SMAP data read alarm")
  if(alarmCheck.or.SMAPsm_obs_struc(n)%startMode) then 

     SMAPsm_obs_struc(n)%startMode = .false. 

     if((SMAPsm_obs_struc(n)%data_designation.eq."SPL3SMP_E").or.&
          (SMAPsm_obs_struc(n)%data_designation.eq."SPL3SMP")) then 

        allocate(smobs_A(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
        allocate(smobs_D(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

        call create_SMAPsmobs_filename(&
             obsdir, &
             SMAPsm_obs_struc(n)%data_designation,&
             LIS_rc%yr, LIS_rc%mo, LIS_rc%da, fname)
        inquire(file=fname,exist=file_exists) 
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
           call read_SMAP_E_data(n,'D',fname,smobs_D)
        else
           write(LIS_logunit,*) '[WARN] Missing SMAP file: ',trim(fname)
        endif
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
           call read_SMAP_E_data(n,'A',fname,smobs_A)
        else
           write(LIS_logunit,*) '[WARN] Missing SMAP file: ',trim(fname)
        endif
        SMAPsm_obs_struc(n)%smobs = LIS_rc%udef
        SMAPsm_obs_struc(n)%smtime = -1

!------------------------------------------------------------------------- 
!   Ascending pass assumed to be at 6pm localtime and the descending 
!   pass is assumed to be at 6am local time
!-------------------------------------------------------------------------
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              grid_index = LIS_domain(n)%gindex(c,r)
              if(grid_index.ne.-1) then 
                 
                 if(smobs_D(c+(r-1)*LIS_rc%lnc(n)).ne.-9999.0) then   
                    SMAPsm_obs_struc(n)%smobs(c,r) = &
                         smobs_D(c+(r-1)*LIS_rc%lnc(n))                 
                    lon = LIS_domain(n)%lon(c+(r-1)*LIS_rc%lnc(n))
                    lhour = 6.0
                    call LIS_localtime2gmt (gmt,lon,lhour,zone)
                    SMAPsm_obs_struc(n)%smtime(c,r) = gmt

                 endif
              end if
!-------------------------------------------------------------------------  
! The ascending data is used only over locations where descending data
! doesn't exist. 
!-------------------------------------------------------------------------
              if(smobs_A(c+(r-1)*LIS_rc%lnc(n)).ne.-9999.0.and.&
                   SMAPsm_obs_struc(n)%smobs(c,r).eq.-9999.0) then   
                 SMAPsm_obs_struc(n)%smobs(c,r) = &
                      smobs_A(c+(r-1)*LIS_rc%lnc(n))                 
                 lon = LIS_domain(n)%lon(c+(r-1)*LIS_rc%lnc(n))
                 lhour = 18.0
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 SMAPsm_obs_struc(n)%smtime(c,r) = gmt
              endif
           enddo
        enddo

        deallocate(smobs_A)
        deallocate(smobs_D)

     else
        
     endif
  end if
 
  call ESMF_StateGet(Obj_Space,"SMAP_sm",smcField,&
       rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=smc,rc=status)
  call LIS_verify(status)

  smc = LIS_rc%udef
  
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           grid_index =LIS_domain(n)%gindex(c,r)

           dt = (LIS_rc%gmt - SMAPsm_obs_struc(n)%smtime(c,r))*3600.0
           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
              smc(grid_index) = & 
                   SMAPsm_obs_struc(n)%smobs(c,r)

           endif
        endif
     enddo
  enddo

  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status)

end subroutine read_SMAPsmobs

!BOP
! 
! !ROUTINE: read_SMAP_E_data
! \label{read_SMAP_E_data}
!
! !INTERFACE:
subroutine read_SMAP_E_data(n, pass, fname, smobs_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use SMAPsm_obsMod,  only : SMAPsm_obs_struc
#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: pass
  character (len=*)             :: fname
  real                          :: smobs_ip(LIS_rc%lnc(n)*LIS_rc%lnr(n))


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the SMOS NESDIS binary file and applies the data
!  quality flags to filter the data. !\normalsize
!
!  tb_time_seconds
!  Arithmetic average of the same parameters found in the 
!  fore- and aft-looking groups in the input SPL1CTB granule. 
!  The resulting parameter thus describes the average of UTC 
!  acquisition times of SPL1BTB observations whose boresights 
!  fall within a 36 km EASE-Grid 2.0 cell. The result is then 
!  expressed in J2000 seconds (the number of seconds since 
!  11:58:55.816 on January 1, 2000 UT).
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTNASASMAP AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDF5)

  character*100,    parameter    :: sm_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: sm_field_name_D = "soil_moisture"
  character*100,    parameter    :: sm_qa_name_D = "retrieval_qual_flag"
  character*100,    parameter    :: sm_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: sm_field_name_A = "soil_moisture_pm"
  character*100,    parameter    :: sm_qa_name_A = "retrieval_qual_flag_pm"

  integer(hsize_t), allocatable  :: dims(:)
  integer(hsize_t), dimension(2) :: dimsm
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hid_t)                 :: memspace
  integer(hid_t)                 :: dataspace
  integer                        :: memrank = 2
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer(hsize_t), dimension(2) :: offset_file
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: sm_gr_id_D,sm_field_id_D,sm_qa_id_D
  integer(hid_t)                 :: sm_gr_id_A,sm_field_id_A,sm_qa_id_A
  real,             allocatable  :: sm_field(:,:)
  integer,          allocatable  :: sm_qa(:,:)
  integer                        :: c,r,t
  logical*1                      :: sm_data_b(SMAPsm_obs_struc(n)%nc*SMAPsm_obs_struc(n)%nr)
  logical*1                      :: smobs_b_ip(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                           :: sm_data(SMAPsm_obs_struc(n)%nc*SMAPsm_obs_struc(n)%nr)
  integer                        :: status,ios

  dimsm      = (/SMAPsm_obs_struc(n)%nc, SMAPsm_obs_struc(n)%nr/)
  count_file = (/SMAPsm_obs_struc(n)%nc, SMAPsm_obs_struc(n)%nr/)
  count_mem  = (/SMAPsm_obs_struc(n)%nc, SMAPsm_obs_struc(n)%nr/)
  
  allocate(sm_field(SMAPsm_obs_struc(n)%nc, SMAPsm_obs_struc(n)%nr))
  allocate(sm_qa(SMAPsm_obs_struc(n)%nc, SMAPsm_obs_struc(n)%nr))
  allocate(dims(2))

  dims(1) = SMAPsm_obs_struc(n)%nc
  dims(2) = SMAPsm_obs_struc(n)%nr

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LIS_verify(status, 'Error opening NASASMAP file ')
  
  if(pass.eq.'D') then 
     call h5gopen_f(file_id,sm_gr_name_D,sm_gr_id_D, status)
     call LIS_verify(status, 'Error opening SM group in NASASMAP file')
     
     call h5dopen_f(sm_gr_id_D,sm_field_name_D,sm_field_id_D, status)
     call LIS_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id_D, dataspace, status)
     call LIS_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LIS_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LIS_verify(status, 'Error in h5create_simple_f; read_NASASMAPsm')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_NASASMAPsm')
     
     call h5dread_f(sm_field_id_D, H5T_NATIVE_REAL,sm_field,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM field from NASASMAPfile')

     call h5dclose_f(sm_field_id_D,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

     call h5dopen_f(sm_gr_id_D,sm_qa_name_D,sm_qa_id_D, status)
     call LIS_verify(status, 'Error opening SM QA field in NASASMAP file')
     
     call h5dread_f(sm_qa_id_D, H5T_NATIVE_INTEGER,sm_qa,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM QA field from NASASMAPfile')
     
     call h5dclose_f(sm_qa_id_D,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(sm_gr_id_D,status)
     call LIS_verify(status,'Error in H5GCLOSE call')

  else
     call h5gopen_f(file_id,sm_gr_name_A,sm_gr_id_A, status)
     call LIS_verify(status, 'Error opening SM group in NASASMAP file')
     
     call h5dopen_f(sm_gr_id_A,sm_field_name_A,sm_field_id_A, status)
     call LIS_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id_A, dataspace, status)
     call LIS_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LIS_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LIS_verify(status, 'Error in h5create_simple_f; read_NASASMAPsm')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_NASASMAPsm')
     
     call h5dread_f(sm_field_id_A, H5T_NATIVE_REAL,sm_field,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM field from NASASMAPfile')

     call h5dclose_f(sm_field_id_A,status)
     call LIS_verify(status,'Error in H5DCLOSE call')
     
     call h5dopen_f(sm_gr_id_A,sm_qa_name_A,sm_qa_id_A, status)
     call LIS_verify(status, 'Error opening SM QA field in NASASMAP file')
     
     call h5dread_f(sm_qa_id_A, H5T_NATIVE_INTEGER,sm_qa,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM QA field from NASASMAPfile')

     call h5dclose_f(sm_qa_id_A,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(sm_gr_id_A,status)
     call LIS_verify(status,'Error in H5GCLOSE call')
     
  endif
  
  call h5fclose_f(file_id,status)
  call LIS_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LIS_verify(status,'Error in H5CLOSE call')
  
  sm_data_b = .false. 
  t = 1

! The retrieval_quality_field variable's binary representation consists of bits
! that indicate whether retrieval is performed or not at a given grid cell. 
! When retrieval is performed, it contains additional bits to further 
! indicate the exit status and quality of the retrieval. The first bit 
! indicates the recommended quality (0-means retrieval has recommended quality
!

  do r=1,SMAPsm_obs_struc(n)%nr
     do c=1,SMAPsm_obs_struc(n)%nc        
        sm_data(t) = sm_field(c,r)
        if(sm_data(t).ne.-9999.0) then 
!           if(SMAPsm_obs_struc(n)%qcFlag.eq.1) then 
              if(ibits(sm_qa(c,r),0,1).eq.0) then 
                 sm_data_b(t) = .true.
              else
                 sm_data(t) = -9999.0
              endif
!           else
!              sm_data_b(t) = .true.
!           endif
           endif
        t = t+1
     enddo
  enddo

  deallocate(sm_qa)

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%gridDesc(n,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       SMAPsm_obs_struc(n)%nc*SMAPsm_obs_struc(n)%nr, &
       LIS_rc%lnc(n)*LIS_rc%lnr(n), &
       LIS_domain(n)%lat, LIS_domain(n)%lon,&
       SMAPsm_obs_struc(n)%n11, LIS_rc%udef, ios)

  deallocate(sm_field)
  deallocate(dims)

#endif

end subroutine read_SMAP_E_data

!BOP
!
! !ROUTINE: create_SMAPsmobs_filename
! \label(create_SMAPsmobs_filename)
!
! !INTERFACE:
subroutine create_SMAPsmobs_filename(&
     obsdir, &
     designation,&
     yr, mo, da, fname)
!
! !USES:
  implicit none
!
! !INPUT PARAMETERS:
  character(len=*), intent(in) :: obsdir
  character(len=*), intent(in) :: designation
  integer  ,        intent(in) :: yr
  integer  ,        intent(in) :: mo
  integer  ,        intent(in) :: da
  character(len=*), intent(out) :: fname
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  if(designation.eq."SPL3SMAP") then 
     fname = trim(obsdir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_AP_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '_R12170_001.h5'
  elseif(designation.eq."SPL3SMP") then 
     fname = trim(obsdir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
         '.h5'
  elseif(designation.eq."SPL3SMP_E") then 
     fname = trim(obsdir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_E_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '.h5'
  endif
end subroutine create_SMAPsmobs_filename


