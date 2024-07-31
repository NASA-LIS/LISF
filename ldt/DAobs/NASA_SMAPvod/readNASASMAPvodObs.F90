!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readNASASMAPvodObs
! \label{readNASASMAPvodObs}
! 
! !REVISION HISTORY: 
!  26 Mar 2019: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readNASASMAPvodObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use NASASMAPvod_obsMod
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the NASA
! SMAP vegetation optical depth retrieval product. 
!
!EOP

  real*8            :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  integer           :: ftn
  integer           :: ierr
  integer           :: mn_ind
  integer           :: yr, mo, da, hr, mn, ss
  integer           :: doy
  real              :: gmt
  character*8       :: yyyymmdd
  character*4       :: yyyy
  character*2       :: mm,dd,hh
  character*200     :: list_files
  character(len=LDT_CONST_PATH_LEN)     :: fname
  character(len=LDT_CONST_PATH_LEN)     :: smap_filename(10)
  real              :: vod_out(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real              :: vodobs(LDT_rc%lnc(n),LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  vodobs= LDT_rc%udef

  if(NASASMAPvodobs(n)%data_designation.eq."SPL2SMP".or.&
       NASASMAPvodobs(n)%data_designation.eq."SPL2SMP_E") then 

     if(LDT_rc%ts.gt.3600) then 
        write(LDT_logunit,*)'[ERR] Please set the LDT timestep to 1hr or less'
        write(LDT_logunit,*)'[ERR] This is required for SMAPL2 data processing'
        call LDT_endrun()
     endif

     write(yyyymmdd,'(i4.4,2i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da
     write(yyyy,'(i4.4)') LDT_rc%yr
     write(mm,'(i2.2)') LDT_rc%mo
     write(dd,'(i2.2)') LDT_rc%da
     write(hh,'(i2.2)') LDT_rc%hr
     
     list_files = 'ls '//trim(NASASMAPvodobs(n)%odir)//&
          '/'//trim(yyyy)//'.'//trim(mm)//'.'//dd//&
          '/SMAP_L2_*'//trim(yyyymmdd)//'T'//trim(hh)&
          //"*.h5 > SMAP_filelist.dat"
     
     call system(trim(list_files))

     i = 1
     ftn = LDT_getNextUnitNumber()
     open(ftn,file="./SMAP_filelist.dat",&
          status='old',iostat=ierr)
     
     do while(ierr.eq.0) 
        read(ftn,'(a)',iostat=ierr) fname
        if(ierr.ne.0) then 
           exit
        endif
        mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh))
        
        mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh))+11        
        read(fname(mn_ind:mn_ind+1),'(i2.2)') mn
        ss=0
        call LDT_tick(timenow,doy,gmt,LDT_rc%yr, LDT_rc%mo, LDT_rc%da, &
                LDT_rc%hr, mn, ss, 0.0)
        
        smap_filename(i) = fname
        
        write(LDT_logunit,*) '[INFO] reading ',trim(smap_filename(i))
        
        call read_SMAPL2vod_data(n,smap_filename(i),vodobs,timenow)
        
        i = i+1
     enddo
     call LDT_releaseUnitNumber(ftn)
     
  else
     call create_NASASMAPvod_filename(NASASMAPvodobs(n)%odir, &
          NASASMAPvodobs(n)%data_designation,&
          LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)
     
     inquire(file=trim(fname),exist=file_exists)
     if(file_exists) then
        
        write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
        call read_SMAP_vod_data(n, fname, vod_out)
        write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)
        
        do r=1,LDT_rc%lnr(n)
           do c=1,LDT_rc%lnc(n)
              if(vod_out(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
                 vodobs(c,r) = vod_out(c+(r-1)*LDT_rc%lnc(n))
              endif
           enddo
        enddo
     endif
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%vod_obs,&
       vodobs,vlevel=1)

end subroutine readNASASMAPvodObs


!BOP
! 
! !ROUTINE: read_SMAPL2vod_data
! \label{read_SMAPL2vod_data}
!
! !INTERFACE:
subroutine read_SMAPL2vod_data(n, fname, vodobs_inp, time)
! 
! !USES:   

  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use NASASMAPvod_obsMod
#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n
  character (len=*)        :: fname
  real                     :: vodobs_inp(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real*8                   :: time

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!
!
!EOP

#if (defined USE_HDF5)

  character*100,    parameter    :: vod_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,    parameter    :: vod_field_name = "vegetation_opacity_option3"
  character*100,    parameter    :: vod_qa_name = "retrieval_qual_flag"

  integer(hsize_t), dimension(1) :: dims
  integer(hsize_t), dimension(1) :: maxdims
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: dspace_id
  integer(hid_t)                 :: row_id, col_id
  integer(hid_t)                 :: vod_gr_id,vod_field_id, vod_qa_id
  integer(hid_t)                 :: vod_gr_id_A,vod_field_id_A
  real,             allocatable  :: vod_field(:)
  integer,          allocatable  :: vod_qa(:)
  integer,          allocatable  :: ease_row(:)
  integer,          allocatable  :: ease_col(:)
  integer                        :: c,r,t
  logical*1                      :: vod_data_b(NASASMAPvodobs(n)%nc*NASASMAPvodobs(n)%nr)
  logical*1                      :: vodobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                           :: vod_data(NASASMAPvodobs(n)%nc*NASASMAPvodobs(n)%nr)
  real                           :: vodobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  integer                        :: status,ios,iret

  call h5open_f(status)
  call LDT_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LDT_verify(status, 'Error opening SMAP L2 file ')
  
  call h5gopen_f(file_id,vod_gr_name,vod_gr_id, status)
  call LDT_verify(status, 'Error opening SM group in SMAP L2 file')
  
  call h5dopen_f(vod_gr_id,vod_field_name,vod_field_id, status)
  call LDT_verify(status, 'Error opening SM field in SMAP L2 file')

  call h5dopen_f(vod_gr_id,"EASE_row_index",row_id, status)
  call LDT_verify(status, 'Error opening row index field in SMAP L2 file')

  call h5dopen_f(vod_gr_id,"EASE_column_index",col_id, status)
  call LDT_verify(status, 'Error opening column index field in SMAP L2 file')

  call h5dopen_f(vod_gr_id, vod_qa_name,vod_qa_id, status)
  call LDT_verify(status, 'Error opening QA field in SMAP L2 file')
  
  call h5dget_space_f(vod_field_id, dspace_id, status)
  call LDT_verify(status, 'Error in h5dget_space_f: readSMAP L2Obs')
  
! Size of the arrays
! This routine returns -1 on failure, rank on success. 
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, status) 
  if(status.eq.-1) then 
     call LDT_verify(status, 'Error in h5sget_simple_extent_dims_f: readSMAP L2Obs')
  endif
  
  allocate(vod_field(maxdims(1)))
  allocate(vod_qa(maxdims(1)))
  allocate(ease_row(maxdims(1)))
  allocate(ease_col(maxdims(1)))

  call h5dread_f(row_id, H5T_NATIVE_INTEGER,ease_row,dims,status)
  call LDT_verify(status, 'Error extracting row index from SMAP L2 file')

  call h5dread_f(col_id, H5T_NATIVE_INTEGER,ease_col,dims,status)
  call LDT_verify(status, 'Error extracting col index from SMAP L2 file')
  
  call h5dread_f(vod_field_id, H5T_NATIVE_REAL,vod_field,dims,status)
  call LDT_verify(status, 'Error extracting SM field from SMAP L2 file')

  call h5dread_f(vod_qa_id, H5T_NATIVE_INTEGER,vod_qa,dims,status)
  call LDT_verify(status, 'Error extracting SM field from SMAP L2 file')
  
  call h5dclose_f(vod_qa_id,status)
  call LDT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(row_id,status)
  call LDT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(col_id,status)
  call LDT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(vod_field_id,status)
  call LDT_verify(status,'Error in H5DCLOSE call')
  
  call h5gclose_f(vod_gr_id,status)
  call LDT_verify(status,'Error in H5GCLOSE call')
    
  call h5fclose_f(file_id,status)
  call LDT_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LDT_verify(status,'Error in H5CLOSE call')

  vod_data = LDT_rc%udef
  vod_data_b = .false. 

!grid the data in EASE projection
  do t=1,maxdims(1)
     if(ibits(vod_qa(t),0,1).eq.0) then 
        vod_data(ease_col(t) + &
             (ease_row(t)-1)*NASASMAPvodobs(n)%nc) = vod_field(t) 
        if(vod_field(t).ne.-9999.0) then 
           vod_data_b(ease_col(t) + &
                (ease_row(t)-1)*NASASMAPvodobs(n)%nc) = .true. 
        endif
     endif
  enddo
  
  t = 1
!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LDT_rc%gridDesc, vod_data_b, vod_data, &
       vodobs_b_ip, vodobs_ip, &
       NASASMAPvodobs(n)%nc*NASASMAPvodobs(n)%nr,&
       LDT_rc%lnc(n)*LDT_rc%lnr(n),&
       LDT_domain(n)%lat, LDT_domain(n)%lon,&
       NASASMAPvodobs(n)%n11,LDT_rc%udef, iret)

  deallocate(vod_field)
  deallocate(vod_qa)
  deallocate(ease_row)
  deallocate(ease_col)

!overwrite the input data 
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(vodobs_ip(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
           vodobs_inp(c,r) = & 
                vodobs_ip(c+(r-1)*LDT_rc%lnc(n))

!           NASASMAPvodobs(n)%vodtime(c,r) = & 
!                time
        endif
     enddo
  enddo

#endif

end subroutine read_SMAPL2vod_data


!BOP
! 
! !ROUTINE: read_SMAP_vod_data
! \label(read_SMAP_vod_data)
!
! !INTERFACE:
subroutine read_SMAP_vod_data(n, fname, vodobs_ip)
! 
! !USES:   
  use LDT_coreMod
  use LDT_logMod
  use map_utils
  use NASASMAPvod_obsMod, only : NASASMAPvodobs
  use LDT_paramDataMod, only : LDT_LSMparam_struc
#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname
  real                          :: vodobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SMAP soil moisture file
!  \item[vodobs\_ip]   VOD data processed to the LDT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
#if (defined USE_HDF5)
  character*100,    parameter    :: vod_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,    parameter    :: vod_field_name = "vegetation_opacity"

  character*100,    parameter    :: vod_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: vod_field_name_D = "vegetation_opacity"
  character*100,    parameter    :: vod_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: vod_field_name_A = "vegetation_opacity_pm"

  integer(hsize_t), allocatable  :: dims(:)
  integer(hsize_t), dimension(2) :: dimsm
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hid_t)                 :: memspace
  integer(hid_t)                 :: dataspace
  integer                        :: memrank = 2
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer(hsize_t), dimension(2) :: offset_file = (/0,0/)
  integer(hid_t)                 :: file_id, vod_gr_id,vod_field_id
  integer(hid_t)                 :: vod_gr_id_D,vod_field_id_D
  integer(hid_t)                 :: vod_gr_id_A,vod_field_id_A
  real,             allocatable  :: vod_field(:,:)
  real,             allocatable  :: vod_field_D(:,:)
  real,             allocatable  :: vod_field_A(:,:)
  integer                        :: c,r,t
  logical*1                      :: vod_data_b(NASASMAPvodobs(n)%nc*NASASMAPvodobs(n)%nr)
  logical*1                      :: vodobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                           :: vod_data(NASASMAPvodobs(n)%nc*NASASMAPvodobs(n)%nr)
  integer                        :: search_rad
  integer                        :: ix, jx, c_s, c_e, r_s,r_e
  integer                        :: status

  dimsm      = (/NASASMAPvodobs(n)%nc, NASASMAPvodobs(n)%nr/)
  count_file = (/NASASMAPvodobs(n)%nc, NASASMAPvodobs(n)%nr/)
  count_mem  = (/NASASMAPvodobs(n)%nc, NASASMAPvodobs(n)%nr/)
  
  allocate(vod_field(NASASMAPvodobs(n)%nc, NASASMAPvodobs(n)%nr))
  allocate(vod_field_D(NASASMAPvodobs(n)%nc, NASASMAPvodobs(n)%nr))
  allocate(vod_field_A(NASASMAPvodobs(n)%nc, NASASMAPvodobs(n)%nr))
  allocate(dims(2))

  dims(1) = NASASMAPvodobs(n)%nc
  dims(2) = NASASMAPvodobs(n)%nr

!  if(NASASMAPvodobs(n)%data_designation.eq."SPL3SMP_E") then 

! MN: The structure of the data in the SPL3SMP R14  onward 
! is similar to the SPL3SMP_E
!  if ( (NASASMAPvodobs(n)%data_designation.eq."SPL3SMP_E") .or. &
!       (NASASMAPvodobs(n)%data_designation.eq."SPL3SMP_R14") ) then 

     call h5open_f(status)
     call LDT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LDT_verify(status, 'Error opening NASASMAP file ')     

!Read the AM (descending) data     
     call h5gopen_f(file_id,vod_gr_name_D,vod_gr_id_D, status)
     call LDT_verify(status, 'Error opening SM group (AM) in NASASMAP file')

     call h5dopen_f(vod_gr_id_D,vod_field_name_D,vod_field_id_D, status)
     call LDT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(vod_field_id_D, dataspace, status)
     call LDT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LDT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LDT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LDT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vod_field_id_D, H5T_NATIVE_REAL,vod_field_D,dims,status, &
          memspace, dataspace)
     call LDT_verify(status, 'Error extracting SM (AM) field from NASASMAPfile')

!Read the PM (ascending) data     
     call h5gopen_f(file_id,vod_gr_name_A,vod_gr_id_A, status)
     call LDT_verify(status, 'Error opening SM group (PM) in NASASMAP file')

     call h5dopen_f(vod_gr_id_A,vod_field_name_A,vod_field_id_A, status)
     call LDT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(vod_field_id_A, dataspace, status)
     call LDT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LDT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LDT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LDT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vod_field_id_A, H5T_NATIVE_REAL,vod_field_A,dims,status, &
          memspace, dataspace)
     call LDT_verify(status, 'Error extracting SM (PM) field from NASASMAPfile')


     call h5dclose_f(vod_field_id_D,status)
     call LDT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(vod_field_id_A,status)
     call LDT_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(vod_gr_id_D,status)
     call LDT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(vod_gr_id_A,status)
     call LDT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LDT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LDT_verify(status,'Error in H5CLOSE call')

     vod_field = LDT_rc%udef
     !blend the AM and PM overpasses
     do r=1,NASASMAPvodobs(n)%nr
        do c=1,NASASMAPvodobs(n)%nc
           if(vod_field_D(c,r).ne.LDT_rc%udef) then
              vod_field(c,r) = vod_field_D(c,r)
           endif
        enddo
     enddo

     do r=1,NASASMAPvodobs(n)%nr
        do c=1,NASASMAPvodobs(n)%nc
           if(vod_field_A(c,r).ne.LDT_rc%udef) then
              if(vod_field(c,r).eq.LDT_rc%udef) then 
                 vod_field(c,r) = vod_field_A(c,r)
              endif
           endif
        enddo
     enddo
!  else
!for older versions R13 or older
#if 0 
     call h5open_f(status)
     call LDT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LDT_verify(status, 'Error opening NASASMAP file ')
     
     call h5gopen_f(file_id,vod_gr_name,vod_gr_id, status)
     call LDT_verify(status, 'Error opening SM group in NASASMAP file')


     call h5dopen_f(vod_gr_id,vod_field_name,vod_field_id, status)
     call LDT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(vod_field_id, dataspace, status)
     call LDT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LDT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LDT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LDT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vod_field_id, H5T_NATIVE_REAL,vod_field,dims,status, &
          memspace, dataspace)
     call LDT_verify(status, 'Error extracting SM field from NASASMAPfile')
     
     call h5dclose_f(vod_field_id,status)
     call LDT_verify(status,'Error in H5DCLOSE call')
     
     call h5gclose_f(vod_gr_id,status)
     call LDT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LDT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LDT_verify(status,'Error in H5CLOSE call')
#endif     
!  endif

  vod_data_b = .false. 
  t = 1

  do r=1,NASASMAPvodobs(n)%nr
     do c=1,NASASMAPvodobs(n)%nc        
!may need this for global applications? TBD -SVK
!        vod_data(t) = vod_field(NASASMAPvodobs(n)%nc-c+1,r)
        vod_data(t) = vod_field(c,r)
        if(vod_data(t).ne.-9999.0) then 
           vod_data_b(t) = .true.
        endif
        t = t+1
     enddo
  enddo
  
!  open(100,file='test_inp.bin',form='unformatted')
!  write(100) vod_data
!  close(100)

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LDT_rc%gridDesc(n,:),&
       vod_data_b, vod_data, vodobs_b_ip, vodobs_ip, &
       NASASMAPvodobs(n)%nc*NASASMAPvodobs(n)%nr, &
       LDT_rc%lnc(n)*LDT_rc%lnr(n), &
       LDT_domain(n)%lat, LDT_domain(n)%lon,&
       NASASMAPvodobs(n)%n11,&
       LDT_rc%udef, status)

!  open(100,file='test_out.bin',form='unformatted')
!  write(100) vodobs_ip
!  close(100)
!  stop

!  print*, LDT_rc%lnc(n), LDT_rc%lnr(n)
  deallocate(vod_field)
  deallocate(dims)
  
#endif

end subroutine read_SMAP_vod_data

!BOP
! !ROUTINE: create_NASASMAPvod_filename
! \label{create_NASASMAPvod_filename}
! 
! !INTERFACE: 
subroutine create_NASASMAPvod_filename(ndir, designation,yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: designation
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the timestamped NASASMAP filename 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the NASASMAP soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated NASASMAP filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  if(designation.eq."SPL3SMAP") then 
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_AP_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '_R12170_001.h5'
! MN:   
! The SMAP file names contain components that change in a way that
! is difficult to programatically generate. So after downloading
! a SMAP data file, a symbolic link was created to it which make it 
! easier to generate file name.
!   For example:
!   SMAP_L3_SM_P_20170902.h5 -> SMAP_L3_SM_P_20170902_R15152_001.h5
  elseif(designation.eq."SPL3SMP") then 
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
         '.h5'
! For example:
! SMAP_L3_SM_P_E_20180811.h5 -> SMAP_L3_SM_P_E_20180811_R16010_001.h5
  elseif(designation.eq."SPL3SMP_E") then 
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_E_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '.h5'
  endif
  
end subroutine create_NASASMAPvod_filename
