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
! !ROUTINE: readNASASMAPsmObs
! \label{readNASASMAPsmObs}
!
! !REVISION HISTORY:
!  21 July 2010: Sujay Kumar, Initial Specification
!  12 Feb 2018: Mahdi Navari, openwater proximity detection was added
!                         edited to read New version of the SPL3SMP_R14 (file structure
!                          differs from the previous versions)
!  31 Aug 2018: Mahdi Navari, Edited to read SPL3SMP.005 & SPL3SMP_E.002
!  04 Jun 2019: Sujay Kumar, Updated to support SMAP L2 retrievals
!  15 Aug 2019: Mahdi Navari, There are several version of SMAP sm data available in each directory
!                  with different Release number and different CRID Version Number. The reader was
!                  modified to read the latest version of data (the reader no longer reads the symbolic
!                  link to the SMAP sm data)
!
! !INTERFACE:
subroutine readNASASMAPsmObs(n)
! !USES:
   use ESMF
   use LDT_coreMod
   use LDT_logMod
   use LDT_constantsMod, only : LDT_CONST_PATH_LEN
   use LDT_timeMgrMod
   use LDT_DAobsDataMod
   use NASASMAPsm_obsMod
   use map_utils

   implicit none
! !ARGUMENTS:
   integer, intent(in) :: n
!
! !DESCRIPTION:
!
! This subroutine provides the data reader for the ESACCI
! soil moisture retrieval product.
!
!EOP

   real*8            :: timenow
   logical           :: alarmCheck
   logical           :: file_exists
   integer           :: c, r, i, j
   character(len=LDT_CONST_PATH_LEN)     :: fname
   integer           :: mn_ind
   integer           :: yr, mo, da, hr, mn, ss
   integer           :: doy
   integer           :: ftn
   integer           :: ierr
   real              :: gmt
   character*8       :: yyyymmdd
   character*4       :: yyyy
   character*2       :: mm, dd, hh
   character*200     :: list_files
   character(len=LDT_CONST_PATH_LEN)     :: smap_filename(10)
   real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   character(len=3) :: CRID
!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations.
!-----------------------------------------------------------------------
   NASASMAPsmobs(n)%smobs = LDT_rc%udef
   smobs = LDT_rc%udef

   if (NASASMAPsmobs(n)%data_designation .eq. "SPL2SMP" .or. &
       NASASMAPsmobs(n)%data_designation .eq. "SPL2SMP_E") then

      if (LDT_rc%ts .gt. 3600) then
         write (LDT_logunit, *) '[ERR] Please set the LDT timestep to 1hr or less'
         write (LDT_logunit, *) '[ERR] This is required for SMAPL2 data processing'
         call LDT_endrun()
      endif

      write (yyyymmdd, '(i4.4,2i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da
      write (yyyy, '(i4.4)') LDT_rc%yr
      write (mm, '(i2.2)') LDT_rc%mo
      write (dd, '(i2.2)') LDT_rc%da
      write (hh, '(i2.2)') LDT_rc%hr

      list_files = 'ls '//trim(NASASMAPsmobs(n)%odir)// &
                   '/'//trim(yyyy)//'.'//trim(mm)//'.'//dd// &
                   '/SMAP_L2_*'//trim(yyyymmdd)//'T'//trim(hh) &
                   //"*.h5 > SMAP_filelist.dat"

      call system(trim(list_files))

      i = 1
      ftn = LDT_getNextUnitNumber()
      open (ftn, file="./SMAP_filelist.dat", &
            status='old', iostat=ierr)

      do while (ierr .eq. 0)
         read (ftn, '(a)', iostat=ierr) fname
         if (ierr .ne. 0) then
            exit
         endif
         mn_ind = index(fname, trim(yyyymmdd)//'T'//trim(hh))

         mn_ind = index(fname, trim(yyyymmdd)//'T'//trim(hh)) + 11
         read (fname(mn_ind:mn_ind + 1), '(i2.2)') mn
         ss = 0
         call LDT_tick(timenow, doy, gmt, LDT_rc%yr, LDT_rc%mo, LDT_rc%da, &
                       LDT_rc%hr, mn, ss, 0.0)

         smap_filename(i) = fname

         write (LDT_logunit, *) '[INFO] reading ', trim(smap_filename(i))

         call read_SMAPL2sm_data(n, smap_filename(i), &
                                 NASASMAPsmobs(n)%smobs, timenow)

         i = i + 1
      enddo
      call LDT_releaseUnitNumber(ftn)

   elseif (NASASMAPsmobs(n)%data_designation .eq. "SPL3SMP_E") then

!----------------------------------------------------------------------------------------------------------------
! create filename for 9 km product
!----------------------------------------------------------------------------------------------------------------
      write (yyyy, '(i4.4)') LDT_rc%yr
      write (mm, '(i2.2)') LDT_rc%mo
      write (dd, '(i2.2)') LDT_rc%da
      write (CRID, '(a)') NASASMAPsmobs(n)%release_number

      list_files = 'ls '//trim(NASASMAPsmobs(n)%odir)//'/'//trim(yyyy)//'.'//trim(mm)//'.'// &
                   trim(dd)//'/SMAP_L3_SM_P_E_' &
                   //trim(yyyy)//trim(mm)//trim(dd)//'_'// &
                         trim(CRID)//'*.h5> SMAP_filelist'// &
                         '.dat'

      call system(trim(list_files))
      i = 1
      ftn = LDT_getNextUnitNumber()
      open (ftn, file="./SMAP_filelist.dat", &
            status='old', iostat=ierr)

! if multiple files for the same time and orbits are present, the latest
! one will overwrite older ones, though multiple (redundant) reads occur.
! This assumes that the 'ls command' will list the files in that order.

      do while (ierr .eq. 0)
         read (ftn, '(a)', iostat=ierr) fname
         if (ierr .ne. 0) then
            exit
         endif
         smap_filename(i) = fname
         write (LDT_logunit, *) '[INFO] reading ', trim(smap_filename(i))
         call read_NASASMAP_data(n, smap_filename(i), smobs)
         i = i + 1
      enddo
      call LDT_releaseUnitNumber(ftn)

!        write(LDT_logunit,*) '[INFO] Reading ..',trim(fname)
!        call read_NASASMAP_data(n, fname, smobs)
!        write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            if (smobs(c + (r - 1)*LDT_rc%lnc(n)) .ne. -9999.0) then
               NASASMAPsmobs(n)%smobs(c, r) = smobs(c + (r - 1)*LDT_rc%lnc(n))
            endif
         enddo
      enddo
      !endif
   elseif (NASASMAPsmobs(n)%data_designation .eq. "SPL3SMP") then
!----------------------------------------------------------------------------------------------------------------
! create filename for 36 km product
!----------------------------------------------------------------------------------------------------------------
      write (yyyy, '(i4.4)') LDT_rc%yr
      write (mm, '(i2.2)') LDT_rc%mo
      write (dd, '(i2.2)') LDT_rc%da
      write (CRID, '(a)') NASASMAPsmobs(n)%release_number

      list_files = 'ls '//trim(NASASMAPsmobs(n)%odir)//'/'//trim(yyyy)//'.'//trim(mm)//'.'// &
                   trim(dd)//'/SMAP_L3_SM_P_' &
                   //trim(yyyy)//trim(mm)//trim(dd)//'_'// &
                         trim(CRID)//'*.h5> SMAP_filelist'// &
                         '.dat'

      call system(trim(list_files))
      i = 1
      ftn = LDT_getNextUnitNumber()
      open (ftn, file="./SMAP_filelist.dat", &
            status='old', iostat=ierr)

! if multiple files for the same time and orbits are present, the latest
! one will overwrite older ones, though multiple (redundant) reads occur.
! This assumes that the 'ls command' will list the files in that order.

      do while (ierr .eq. 0)
         read (ftn, '(a)', iostat=ierr) fname
         if (ierr .ne. 0) then
            exit
         endif
         smap_filename(i) = fname
         write (LDT_logunit, *) '[INFO] reading ', trim(smap_filename(i))
         call read_NASASMAP_data(n, smap_filename(i), smobs)
         i = i + 1
      enddo

      call LDT_releaseUnitNumber(ftn)
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            if (smobs(c + (r - 1)*LDT_rc%lnc(n)) .ne. -9999.0) then
               NASASMAPsmobs(n)%smobs(c, r) = smobs(c + (r - 1)*LDT_rc%lnc(n))
            endif
         enddo
      enddo
   endif ! sensor

   call LDT_logSingleDAobs(n, LDT_DAobsData(n)%soilmoist_obs, &
                           NASASMAPsmobs(n)%smobs, vlevel=1)

end subroutine readNASASMAPsmObs

!BOP
! 
! !ROUTINE: read_SMAPL2sm_data
! \label{read_SMAPL2sm_data}
!
! !INTERFACE:
subroutine read_SMAPL2sm_data(n, fname, smobs_inp, time)
! 
! !USES:   

  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use NASASMAPsm_obsMod
#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n
  character (len=*)        :: fname
  real                     :: smobs_inp(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real*8                   :: time

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!
!
!EOP

#if (defined USE_HDF5)

  character*100,    parameter    :: sm_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,    parameter    :: sm_field_name = "soil_moisture"
  character*100,    parameter    :: sm_qa_name = "retrieval_qual_flag"

  integer(hsize_t), dimension(1) :: dims
  integer(hsize_t), dimension(1) :: maxdims
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: dspace_id
  integer(hid_t)                 :: row_id, col_id
  integer(hid_t)                 :: sm_gr_id,sm_field_id, sm_qa_id
  integer(hid_t)                 :: sm_gr_id_A,sm_field_id_A
  real,             allocatable  :: sm_field(:)
  integer,          allocatable  :: sm_qa(:)
  integer,          allocatable  :: ease_row(:)
  integer,          allocatable  :: ease_col(:)
  integer                        :: c,r,t
  logical*1                      :: sm_data_b(NASASMAPsmobs(n)%nc*NASASMAPsmobs(n)%nr)
  logical*1                      :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                           :: sm_data(NASASMAPsmobs(n)%nc*NASASMAPsmobs(n)%nr)
  real                           :: smobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  integer                        :: status,ios,iret

  call h5open_f(status)
  call LDT_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LDT_verify(status, 'Error opening SMAP L2 file ')
  
  call h5gopen_f(file_id,sm_gr_name,sm_gr_id, status)
  call LDT_verify(status, 'Error opening SM group in SMAP L2 file')
  
  call h5dopen_f(sm_gr_id,sm_field_name,sm_field_id, status)
  call LDT_verify(status, 'Error opening SM field in SMAP L2 file')

  call h5dopen_f(sm_gr_id,"EASE_row_index",row_id, status)
  call LDT_verify(status, 'Error opening row index field in SMAP L2 file')

  call h5dopen_f(sm_gr_id,"EASE_column_index",col_id, status)
  call LDT_verify(status, 'Error opening column index field in SMAP L2 file')

!  call h5dopen_f(sm_gr_id, sm_qa_name,sm_qa_id, status)
!  call LDT_verify(status, 'Error opening QA field in SMAP L2 file')
  
  call h5dget_space_f(sm_field_id, dspace_id, status)
  call LDT_verify(status, 'Error in h5dget_space_f: readSMAP L2Obs')
  
! Size of the arrays
! This routine returns -1 on failure, rank on success. 
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, status) 
  if(status.eq.-1) then 
     call LDT_verify(status, 'Error in h5sget_simple_extent_dims_f: readSMAP L2Obs')
  endif
  
  allocate(sm_field(maxdims(1)))
!  allocate(sm_qa(maxdims(1)))
  allocate(ease_row(maxdims(1)))
  allocate(ease_col(maxdims(1)))

  call h5dread_f(row_id, H5T_NATIVE_INTEGER,ease_row,dims,status)
  call LDT_verify(status, 'Error extracting row index from SMAP L2 file')

  call h5dread_f(col_id, H5T_NATIVE_INTEGER,ease_col,dims,status)
  call LDT_verify(status, 'Error extracting col index from SMAP L2 file')
  
  call h5dread_f(sm_field_id, H5T_NATIVE_REAL,sm_field,dims,status)
  call LDT_verify(status, 'Error extracting SM field from SMAP L2 file')

!  call h5dread_f(sm_qa_id, H5T_NATIVE_INTEGER,sm_qa,dims,status)
!  call LDT_verify(status, 'Error extracting SM field from SMAP L2 file')
  
!  call h5dclose_f(sm_qa_id,status)
!  call LDT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(row_id,status)
  call LDT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(col_id,status)
  call LDT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(sm_field_id,status)
  call LDT_verify(status,'Error in H5DCLOSE call')
  
  call h5gclose_f(sm_gr_id,status)
  call LDT_verify(status,'Error in H5GCLOSE call')
    
  call h5fclose_f(file_id,status)
  call LDT_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LDT_verify(status,'Error in H5CLOSE call')

  sm_data = LDT_rc%udef
  sm_data_b = .false. 

!grid the data in EASE projection
  do t=1,maxdims(1)
     if(ease_col(t).gt.0.and.ease_row(t).gt.0) then 
        sm_data(ease_col(t) + &
             (ease_row(t)-1)*NASASMAPsmobs(n)%nc) = sm_field(t) 
        if(sm_field(t).ne.-9999.0) then 
           sm_data_b(ease_col(t) + &
                (ease_row(t)-1)*NASASMAPsmobs(n)%nc) = .true. 
        endif
     endif
  enddo
  
  t = 1
!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LDT_rc%gridDesc, sm_data_b, sm_data, &
       smobs_b_ip, smobs_ip, &
       NASASMAPsmobs(n)%nc*NASASMAPsmobs(n)%nr,&
       LDT_rc%lnc(n)*LDT_rc%lnr(n),&
       LDT_domain(n)%lat, LDT_domain(n)%lon,&
       NASASMAPsmobs(n)%n11,LDT_rc%udef, iret)

  deallocate(sm_field)
!  deallocate(sm_qa)
  deallocate(ease_row)
  deallocate(ease_col)

!overwrite the input data 
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(smobs_ip(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
           smobs_inp(c,r) = & 
                smobs_ip(c+(r-1)*LDT_rc%lnc(n))

!           NASASMAPsmobs(n)%smtime(c,r) = & 
!                time
        endif
     enddo
  enddo

#endif

end subroutine read_SMAPL2sm_data


!BOP
! 
! !ROUTINE: read_NASASMAP_data
! \label(read_NASASMAP_data)
!
! !INTERFACE:
subroutine read_NASASMAP_data(n, fname, smobs_ip)
! 
! !USES:   
  use LDT_coreMod
  use LDT_logMod
  use map_utils
  use NASASMAPsm_obsMod, only : NASASMAPsmobs
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
  real                          :: smobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SMOS NESDIS file
!  \item[smobs\_ip]    soil moisture data processed to the LDT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
#if (defined USE_HDF5)
  character*100,    parameter    :: sm_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,    parameter    :: sm_field_name = "soil_moisture"

  character*100,    parameter    :: sm_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: sm_field_name_D = "soil_moisture"
  character*100,    parameter    :: sm_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: sm_field_name_A = "soil_moisture_pm"

  integer(hsize_t), allocatable  :: dims(:)
  integer(hsize_t), dimension(2) :: dimsm
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hid_t)                 :: memspace
  integer(hid_t)                 :: dataspace
  integer                        :: memrank = 2
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer(hsize_t), dimension(2) :: offset_file = (/0,0/)
  integer(hid_t)                 :: file_id, sm_gr_id,sm_field_id
  integer(hid_t)                 :: sm_gr_id_D,sm_field_id_D
  integer(hid_t)                 :: sm_gr_id_A,sm_field_id_A
  real,             allocatable  :: sm_field(:,:)
  real,             allocatable  :: sm_field_D(:,:)
  real,             allocatable  :: sm_field_A(:,:)
  integer                        :: c,r,t
  logical*1                      :: sm_data_b(NASASMAPsmobs(n)%nc*NASASMAPsmobs(n)%nr)
  logical*1                      :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                           :: sm_data(NASASMAPsmobs(n)%nc*NASASMAPsmobs(n)%nr)
  integer                        :: search_rad
  integer                        :: ix, jx, c_s, c_e, r_s,r_e
  integer                        :: status

  dimsm      = (/NASASMAPsmobs(n)%nc, NASASMAPsmobs(n)%nr/)
  count_file = (/NASASMAPsmobs(n)%nc, NASASMAPsmobs(n)%nr/)
  count_mem  = (/NASASMAPsmobs(n)%nc, NASASMAPsmobs(n)%nr/)
  
  allocate(sm_field(NASASMAPsmobs(n)%nc, NASASMAPsmobs(n)%nr))
  allocate(sm_field_D(NASASMAPsmobs(n)%nc, NASASMAPsmobs(n)%nr))
  allocate(sm_field_A(NASASMAPsmobs(n)%nc, NASASMAPsmobs(n)%nr))
  allocate(dims(2))

  dims(1) = NASASMAPsmobs(n)%nc
  dims(2) = NASASMAPsmobs(n)%nr

!  if(NASASMAPsmobs(n)%data_designation.eq."SPL3SMP_E") then 

! MN: The structure of the data in the SPL3SMP R14  onward 
! is similar to the SPL3SMP_E
!  if ( (NASASMAPsmobs(n)%data_designation.eq."SPL3SMP_E") .or. &
!       (NASASMAPsmobs(n)%data_designation.eq."SPL3SMP_R14") ) then 

     call h5open_f(status)
     call LDT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LDT_verify(status, 'Error opening NASASMAP file ')     

!Read the AM (descending) data     
     call h5gopen_f(file_id,sm_gr_name_D,sm_gr_id_D, status)
     call LDT_verify(status, 'Error opening SM group (AM) in NASASMAP file')

     call h5dopen_f(sm_gr_id_D,sm_field_name_D,sm_field_id_D, status)
     call LDT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id_D, dataspace, status)
     call LDT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LDT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LDT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LDT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(sm_field_id_D, H5T_NATIVE_REAL,sm_field_D,dims,status, &
          memspace, dataspace)
     call LDT_verify(status, 'Error extracting SM (AM) field from NASASMAPfile')

!Read the PM (ascending) data     
     call h5gopen_f(file_id,sm_gr_name_A,sm_gr_id_A, status)
     call LDT_verify(status, 'Error opening SM group (PM) in NASASMAP file')

     call h5dopen_f(sm_gr_id_A,sm_field_name_A,sm_field_id_A, status)
     call LDT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id_A, dataspace, status)
     call LDT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LDT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LDT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LDT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(sm_field_id_A, H5T_NATIVE_REAL,sm_field_A,dims,status, &
          memspace, dataspace)
     call LDT_verify(status, 'Error extracting SM (PM) field from NASASMAPfile')


     call h5dclose_f(sm_field_id_D,status)
     call LDT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(sm_field_id_A,status)
     call LDT_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(sm_gr_id_D,status)
     call LDT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(sm_gr_id_A,status)
     call LDT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LDT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LDT_verify(status,'Error in H5CLOSE call')

     sm_field = LDT_rc%udef
     !blend the AM and PM overpasses
     do r=1,NASASMAPsmobs(n)%nr
        do c=1,NASASMAPsmobs(n)%nc
           if(sm_field_D(c,r).ne.LDT_rc%udef) then
              sm_field(c,r) = sm_field_D(c,r)
           endif
        enddo
     enddo

     do r=1,NASASMAPsmobs(n)%nr
        do c=1,NASASMAPsmobs(n)%nc
           if(sm_field_A(c,r).ne.LDT_rc%udef) then
              if(sm_field(c,r).eq.LDT_rc%udef) then 
                 sm_field(c,r) = sm_field_A(c,r)
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
     
     call h5gopen_f(file_id,sm_gr_name,sm_gr_id, status)
     call LDT_verify(status, 'Error opening SM group in NASASMAP file')


     call h5dopen_f(sm_gr_id,sm_field_name,sm_field_id, status)
     call LDT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id, dataspace, status)
     call LDT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LDT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LDT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LDT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(sm_field_id, H5T_NATIVE_REAL,sm_field,dims,status, &
          memspace, dataspace)
     call LDT_verify(status, 'Error extracting SM field from NASASMAPfile')
     
     call h5dclose_f(sm_field_id,status)
     call LDT_verify(status,'Error in H5DCLOSE call')
     
     call h5gclose_f(sm_gr_id,status)
     call LDT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LDT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LDT_verify(status,'Error in H5CLOSE call')
#endif     
!  endif

  sm_data_b = .false. 
  t = 1

  do r=1,NASASMAPsmobs(n)%nr
     do c=1,NASASMAPsmobs(n)%nc        
!may need this for global applications? TBD -SVK
!        sm_data(t) = sm_field(NASASMAPsmobs(n)%nc-c+1,r)
        sm_data(t) = sm_field(c,r)
!        sm_data(t) = sm_field(c,NASASMAPsmobs(n)%nr-r+1)
        if(sm_data(t).ne.-9999.0) then 
           sm_data_b(t) = .true.
        endif
        t = t+1
     enddo
  enddo
  
!  open(100,file='test_inp.bin',form='unformatted')
!  write(100) sm_data
!  close(100)

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LDT_rc%gridDesc(n,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       NASASMAPsmobs(n)%nc*NASASMAPsmobs(n)%nr, &
       LDT_rc%lnc(n)*LDT_rc%lnr(n), &
       LDT_domain(n)%lat, LDT_domain(n)%lon,&
       NASASMAPsmobs(n)%n11,&
       LDT_rc%udef, status)

!  print*, LDT_rc%lnc(n), LDT_rc%lnr(n)
  deallocate(sm_field)
  deallocate(dims)
  
  !------------------------------------------------------------------------
  !  Remove pixel close to open water
  !------------------------------------------------------------------------
  do r = 1, LDT_rc%lnr(n)
     do c = 1, LDT_rc%lnc(n)
        if (smobs_ip(c+(r-1)*LDT_rc%lnc(n)) .ne. LDT_rc%udef) then 
           search_rad = nint(NASASMAPsmobs(n)%search_radius)
           
           c_s = max(1,c-search_rad)
           c_e = min(LDT_rc%lnc(n),c+search_rad)
           
           r_s = max(1,r-search_rad)
           r_e = min(LDT_rc%lnr(n),r+search_rad)
           
           do ix=c_s,c_e
              do jx=r_s,r_e
                 if( LDT_LSMparam_struc(n)%landmask%value(ix,jx,1) .eq. 0 ) then
                    smobs_ip(c+(r-1)*LDT_rc%lnc(n)) = LDT_rc%udef
                 end if
              enddo
           enddo
        endif
     enddo
  enddo

#endif

end subroutine read_NASASMAP_data


# if 0
!BOP
! !ROUTINE: create_NASASMAPsm_filename
! \label{create_NASASMAPsm_filename}
! 
! !INTERFACE: 
subroutine create_NASASMAPsm_filename(ndir, designation,yr, mo,da, filename)
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
  
end subroutine create_NASASMAPsm_filename
#endif
