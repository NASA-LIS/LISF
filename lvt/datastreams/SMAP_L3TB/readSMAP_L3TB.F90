!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
!
! !ROUTINE: readSMAP_L3TB
! \label{readSMAP_L3TB}
!
! !INTERFACE:
subroutine readSMAP_L3TB(source)
!
! !USES:
   use ESMF
   use LVT_coreMod, only: LVT_rc
   use LVT_histDataMod
   use LVT_logMod !,       only : LVT_logunit
   use SMAP_L3TBMod, only: SMAP_L3TB

   implicit none
!
! !INPUT PARAMETERS:
   integer, intent(in)       :: source
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! This subroutine provides the data reader for the standard
! NASA soil moisture retrieval product.
!
! !FILES USED:
!
! !REVISION HISTORY:
!  19 Nov 2018: Mahdi Navari , Sujay Kumar, Initial Specification
!  March 2019: Peter Shellito  Fix polarization field typos
! 9  July 2019: Mahdi Navari fixed bug in file name generator
! 11  July 2019: Mahdi Navari, There are several version of SMAP sm data available in each directory
!                  with different Release number and different CRID Version Number. The reader was
!                  modified to read the latest version of data (the reader no longer reads the symbolic
!                  link to the SMAP sm data)
!EOP

   logical           :: alarmcheck, file_exists, readflag
   integer           :: iret
   character*200     :: fname
   real              :: Tb(LVT_rc%lnc, LVT_rc%lnr, 4)
   integer           :: fnd
   real              :: timenow
   character*4       :: yyyy
   character*2       :: mm, dd, hh
   integer               :: yr, mo, da, hr, mn, ss
   integer               :: doy
   character*200      :: list_files
   integer               :: ftn, ierr
   character(len=3) :: CRID

   Tb = LVT_rc%udef

   timenow = float(LVT_rc%dhr(source))*3600 + &
             60*LVT_rc%dmn(source) + LVT_rc%dss(source)
   alarmcheck = (mod(timenow, 86400.0) .eq. 0)
   if (SMAP_L3TB(source)%startflag .or. alarmCheck .or. &
       LVT_rc%resetFlag(source)) then
      LVT_rc%resetFlag(source) = .false.
      SMAP_L3TB(source)%startflag = .false.

      if (SMAP_L3TB(source)%data_designation .eq. "SPL3SMP_E") then
!----------------------------------------------------------------------------------------------------------------
! create filename for 9 km product
!----------------------------------------------------------------------------------------------------------------
         write (yyyy, '(i4.4)') LVT_rc%dyr(source)
         write (mm, '(i2.2)') LVT_rc%dmo(source)
         write (dd, '(i2.2)') LVT_rc%dda(source)
         write (CRID, '(a)') SMAP_L3TB(source)%release_number

         list_files = 'ls '//trim(SMAP_L3TB(source)%odir)//'/'//trim(yyyy)//'.'//trim(mm)//'.'// &
                      trim(dd)//'/SMAP_L3_SM_P_E_' &
                      //trim(yyyy)//trim(mm)//trim(dd)//'_'// &
                         trim(CRID)//'*.h5> SMAP_filelist'// &
                         '.dat'

         call system(trim(list_files))
         ftn = LVT_getNextUnitNumber()
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
            write (LVT_logunit, *) '[INFO] Reading SMAP file ', trim(fname)
            call read_SMAP_L3Tb(source, fname, Tb)
         enddo
         call LVT_releaseUnitNumber(ftn)

      elseif (SMAP_L3TB(source)%data_designation .eq. "SPL3SMP") then
!----------------------------------------------------------------------------------------------------------------
! create filename for 36 km product
!----------------------------------------------------------------------------------------------------------------
         write (yyyy, '(i4.4)') LVT_rc%dyr(source)
         write (mm, '(i2.2)') LVT_rc%dmo(source)
         write (dd, '(i2.2)') LVT_rc%dda(source)
         write (CRID, '(a)') SMAP_L3TB(source)%release_number

         list_files = 'ls '//trim(SMAP_L3TB(source)%odir)//'/'&
              //trim(yyyy)//'.'//trim(mm)//'.'// &
              trim(dd)//'/SMAP_L3_SM_P_' &
              //trim(yyyy)//trim(mm)//trim(dd)//'_'// &
              trim(CRID)//'*.h5> SMAP_filelist'// &
              '.dat'
         
         call system(trim(list_files))
         ftn = LVT_getNextUnitNumber()
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
            write (LVT_logunit, *) '[INFO] Reading SMAP file ', trim(fname)
            call read_SMAP_L3Tb(source, fname, Tb)
         enddo

         call LVT_releaseUnitNumber(ftn)
      endif ! sensor
   endif ! alaram

#if 0
! MN : CHECK
   call LVT_logSingleDataStreamVar(LVT_MOC_L3TB, source, &
                                   Tb(:, :, 1), vlevel=1, units="K")
   call LVT_logSingleDataStreamVar(LVT_MOC_L3TB, source, &
                                   Tb(:, :, 2), vlevel=1, units="K")
   call LVT_logSingleDataStreamVar(LVT_MOC_L3TB, source, &
                                   Tb(:, :, 3), vlevel=1, units="K")
   call LVT_logSingleDataStreamVar(LVT_MOC_L3TB, source, &
                                   Tb(:, :, 4), vlevel=1, units="K")
#endif

   if (LVT_MOC_L3TBv_D(source) .ge. 1) then
      call LVT_logSingleDataStreamVar(LVT_MOC_L3TBv_D, source, &
                                      Tb(:, :, 1), vlevel=1, units="K")
   endif
   if (LVT_MOC_L3TBv_D(source) .ge. 1) then
      call LVT_logSingleDataStreamVar(LVT_MOC_L3TBv_A, source, &
                                      Tb(:, :, 2), vlevel=1, units="K")
   endif
   if (LVT_MOC_L3TBv_D(source) .ge. 1) then
      call LVT_logSingleDataStreamVar(LVT_MOC_L3TBh_D, source, &
                                      Tb(:, :, 3), vlevel=1, units="K")
   endif
   if (LVT_MOC_L3TBv_D(source) .ge. 1) then
      call LVT_logSingleDataStreamVar(LVT_MOC_L3TBh_A, source, &
                                      Tb(:, :, 4), vlevel=1, units="K")
   endif

!  open(100,file='test.bin',form='unformatted')
!  write(100) Tbv
!  close(100)
!  stop
end subroutine readSMAP_L3TB

!BOP
! 
! !ROUTINE: read_SMAP_L3Tb
! \label{read_SMAP_L3Tb}
!
! !INTERFACE: 
subroutine read_SMAP_L3Tb(source, fname, L3TB)
! 
! !USES: 
  use LVT_coreMod,         only : LVT_rc,LVT_domain
  use LVT_logMod
  use SMAP_L3TBMod, only : SMAP_L3TB

#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none

  integer                        :: source 
  character(len=*)               :: fname
  real                           :: L3TB(LVT_rc%lnc,LVT_rc%lnr,4)
  
!
! !INPUT PARAMETERS: 
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
!BOP

#if (defined USE_HDF5)

  character*100,   parameter    :: Tb_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,   parameter    :: Tbv_field_name = "soil_moisture"

  character*100,    parameter    :: Tb_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: Tbv_field_name_D = "tb_v_corrected" 
  character*100,    parameter    :: Tbh_field_name_D = "tb_h_corrected"
  character*100,    parameter    :: Tb_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: Tbv_field_name_A = "tb_v_corrected_pm"
  character*100,    parameter    :: Tbh_field_name_A = "tb_h_corrected_pm"

  integer(hid_t)                :: file_id, Tbv_gr_id,Tbv_field_id,Tbh_gr_id,Tbh_field_id
  integer(hid_t)                :: Tbv_gr_id_D,Tbv_field_id_D, Tbh_gr_id_D,Tbh_field_id_D
  integer(hid_t)                :: Tbv_gr_id_A,Tbv_field_id_A, Tbh_gr_id_A,Tbh_field_id_A
  integer(hid_t)                :: dataspace
  integer(hid_t)                :: memspace
  integer                       :: memrank = 2
  integer(hsize_t), allocatable :: dims(:)
  integer(hsize_t), dimension(2) :: dimTb
  integer(hsize_t), dimension(2) :: offset_file
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  !real,             allocatable  :: Tbv_field(:,:)
  real,             allocatable  :: Tbv_field_D(:,:)
  real,             allocatable  :: Tbv_field_A(:,:)
  real,             allocatable  :: Tbh_field_D(:,:)
  real,             allocatable  :: Tbh_field_A(:,:)
  real                           :: Tb1d(SMAP_L3TB(source)%nc*SMAP_L3TB(source)%nr)
  real                           :: Tbv_D_ip(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: Tbv_A_ip(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: Tbh_D_ip(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: Tbh_A_ip(LVT_rc%lnc*LVT_rc%lnr)
  logical*1                      :: li(SMAP_L3TB(source)%nc*SMAP_L3TB(source)%nr)
  logical*1                      :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: udef
  integer                        :: t, c,r,c1,r1
  integer                        :: iret,status
 
  dimTb      = (/SMAP_L3TB(source)%nc, SMAP_L3TB(source)%nr/)
  count_file = (/SMAP_L3TB(source)%nc, SMAP_L3TB(source)%nr/)
  count_mem  = (/SMAP_L3TB(source)%nc, SMAP_L3TB(source)%nr/)
  
  !allocate(Tbv_field(SMAP_L3TB(source)%nc, SMAP_L3TB(source)%nr))
  allocate(Tbv_field_D(SMAP_L3TB(source)%nc, SMAP_L3TB(source)%nr))
  allocate(Tbv_field_A(SMAP_L3TB(source)%nc, SMAP_L3TB(source)%nr))
  allocate(Tbh_field_D(SMAP_L3TB(source)%nc, SMAP_L3TB(source)%nr))
  allocate(Tbh_field_A(SMAP_L3TB(source)%nc, SMAP_L3TB(source)%nr))
  allocate(dims(2))

  dims(1) = SMAP_L3TB(source)%nc
  dims(2) = SMAP_L3TB(source)%nr

  if ((SMAP_L3TB(source)%data_designation.eq."SPL3SMP_E") .or. & 
      (SMAP_L3TB(source)%data_designation.eq."SPL3SMP") ) then 
     call h5open_f(status)
     call LVT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LVT_verify(status, 'Error opening NASASMAP file ')     

!Read  the AM (descending) V-pol  data     
     call h5gopen_f(file_id,Tb_gr_name_D,Tbv_gr_id_D, status)
     call LVT_verify(status, 'Error opening SM group (AM) in NASASMAP file')

     call h5dopen_f(Tbv_gr_id_D,Tbv_field_name_D,Tbv_field_id_D, status)
     call LVT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(Tbv_field_id_D, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimTb, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(Tbv_field_id_D, H5T_NATIVE_REAL,Tbv_field_D,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting SM (AM) field from NASASMAPfile')

!Read the PM (ascending) V-pol data     
     call h5gopen_f(file_id,Tb_gr_name_A,Tbv_gr_id_A, status)
     call LVT_verify(status, 'Error opening SM group (PM) in NASASMAP file')

     call h5dopen_f(Tbv_gr_id_A,Tbv_field_name_A,Tbv_field_id_A, status)
     call LVT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(Tbv_field_id_A, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimTb, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(Tbv_field_id_A, H5T_NATIVE_REAL,Tbv_field_A,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting SM (AM) field from NASASMAPfile')


!Read  the AM (descending) H-pol  data     
     call h5gopen_f(file_id,Tb_gr_name_D,Tbh_gr_id_D, status)
     call LVT_verify(status, 'Error opening SM group (AM) in NASASMAP file')

     call h5dopen_f(Tbh_gr_id_D,Tbh_field_name_D,Tbh_field_id_D, status)
     call LVT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(Tbh_field_id_D, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimTb, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(Tbh_field_id_D, H5T_NATIVE_REAL,Tbh_field_D,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting SM (AM) field from NASASMAPfile')

!Read the PM (ascending) H-pol data     
     call h5gopen_f(file_id,Tb_gr_name_A,Tbh_gr_id_A, status)
     call LVT_verify(status, 'Error opening SM group (PM) in NASASMAP file')

     call h5dopen_f(Tbh_gr_id_A,Tbh_field_name_A,Tbh_field_id_A, status)
     call LVT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(Tbh_field_id_A, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimTb, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(Tbh_field_id_A, H5T_NATIVE_REAL,Tbh_field_A,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting SM (AM) field from NASASMAPfile')


     call h5dclose_f(Tbv_field_id_D,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(Tbv_field_id_A,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(Tbv_gr_id_D,status)
     call LVT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(Tbv_gr_id_A,status)
     call LVT_verify(status,'Error in H5GCLOSE call')
  

     call h5dclose_f(Tbh_field_id_D,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(Tbh_field_id_A,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(Tbh_gr_id_D,status)
     call LVT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(Tbh_gr_id_A,status)
     call LVT_verify(status,'Error in H5GCLOSE call')

   
     call h5fclose_f(file_id,status)
     call LVT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LVT_verify(status,'Error in H5CLOSE call')

! no belnding 
#if 0
     Tbv_field = LVT_rc%udef
     !blend the AM and PM overpasses
     do r=1,SMAP_L3TB(source)%nr
        do c=1,SMAP_L3TB(source)%nc
           if(Tbv_field_D(c,r).ne.LVT_rc%udef) then
              Tbv_field(c,r) = Tbv_field_D(c,r)
           endif
        enddo
     enddo
     do r=1,SMAP_L3TB(source)%nr
        do c=1,SMAP_L3TB(source)%nc
           if(Tbv_field_A(c,r).ne.LVT_rc%udef) then
              if(Tbv_field(c,r).eq.LVT_rc%udef) then 
                 Tbv_field(c,r) = Tbv_field_A(c,r)
              endif
           endif
        enddo
     enddo
#endif

  else

     write(LVT_logunit,*) '[MSG] This reader read SMAP SPL3SMP_E and '//&
        'SPL3SMP V4 onward '
#if 0 
     call h5open_f(status)
     call LVT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LVT_verify(status, 'Error opening NASASMAP file ')
     
     call h5gopen_f(file_id,Tb_gr_name_D,Tbv_gr_id_D, status)
     call LVT_verify(status, 'Error opening SM group in NASASMAP file')

     call h5dopen_f(Tbv_gr_id_D,Tbv_field_name_D,Tbv_field_id_D, status)
     call LVT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(Tbv_field_id_D, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimTb, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(Tbv_field_id_D, H5T_NATIVE_REAL,Tbv_field_D,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting SM field from NASASMAPfile')

     call h5gopen_f(file_id,Tb_gr_name_A,Tbv_gr_id_A, status)
     call LVT_verify(status, 'Error opening SM group in NASASMAP file')

     call h5dopen_f(Tbv_gr_id_A,Tbv_field_name_A,Tbv_field_id_A, status)
     call LVT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(Tbv_field_id_A, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimTb, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(Tbv_field_id_A, H5T_NATIVE_REAL,Tbv_field_A,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting SM field from NASASMAPfile')
     
     call h5dclose_f(Tbv_field_id_D,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(Tbv_field_id_A,status)
     call LVT_verify(status,'Error in H5DCLOSE call')
     
     call h5gclose_f(Tbv_gr_id_D,status)
     call LVT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(Tbv_gr_id_A,status)
     call LVT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LVT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LVT_verify(status,'Error in H5CLOSE call')
#endif


! no belnding 
!#if 0     
!     Tbv_field = LVT_rc%udef
!     !blend the AM and PM overpasses
!     do r=1,SMAP_L3TB(source)%nr
!        do c=1,SMAP_L3TB(source)%nc
!           if(Tbv_field_D(c,r).ne.LVT_rc%udef) then
!              Tbv_field(c,r) = Tbv_field_D(c,r)
!           endif
!        enddo
!     enddo
!     do r=1,SMAP_L3TB(source)%nr
!        do c=1,SMAP_L3TB(source)%nc
!           if(Tbv_field_A(c,r).ne.LVT_rc%udef) then
!              if(Tbv_field(c,r).eq.LVT_rc%udef) then 
 !                Tbv_field(c,r) = Tbv_field_A(c,r)
!              endif
!           endif
!        enddo
!     enddo
! #endif
endif

! V-pol interpolation 
  li = .false. 
  t = 1

  do r=1,SMAP_L3TB(source)%nr
     do c=1,SMAP_L3TB(source)%nc        
        Tb1d(t) = Tbv_field_D(c,r)
        if(Tb1d(t).ne.-9999.0) then 
           li(t) = .true.
        endif
        t = t+1
     enddo
  enddo
  
  udef = -9999.0

  call neighbor_interp(LVT_rc%gridDesc, li, Tb1d, lo, Tbv_D_ip, &
       SMAP_L3TB(source)%nc*SMAP_L3TB(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMAP_L3TB(source)%rlat2, SMAP_L3TB(source)%rlon2,&
       SMAP_L3TB(source)%n112,udef, iret)



  li = .false. 
  t = 1

  do r=1,SMAP_L3TB(source)%nr
     do c=1,SMAP_L3TB(source)%nc        
        Tb1d(t) = Tbv_field_A(c,r)
        if(Tb1d(t).ne.-9999.0) then 
           li(t) = .true.
        endif
        t = t+1
     enddo
  enddo
  
  udef = -9999.0

  call neighbor_interp(LVT_rc%gridDesc, li, Tb1d, lo, Tbv_A_ip, &
       SMAP_L3TB(source)%nc*SMAP_L3TB(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMAP_L3TB(source)%rlat2, SMAP_L3TB(source)%rlon2,&
       SMAP_L3TB(source)%n112,udef, iret)


! H-pol interpolation 
  li = .false. 
  t = 1

  do r=1,SMAP_L3TB(source)%nr
     do c=1,SMAP_L3TB(source)%nc        
        Tb1d(t) = Tbh_field_D(c,r)
        if(Tb1d(t).ne.-9999.0) then 
           li(t) = .true.
        endif
        t = t+1
     enddo
  enddo
  
  udef = -9999.0

  call neighbor_interp(LVT_rc%gridDesc, li, Tb1d, lo, Tbh_D_ip, &
       SMAP_L3TB(source)%nc*SMAP_L3TB(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMAP_L3TB(source)%rlat2, SMAP_L3TB(source)%rlon2,&
       SMAP_L3TB(source)%n112,udef, iret)



  li = .false. 
  t = 1

  do r=1,SMAP_L3TB(source)%nr
     do c=1,SMAP_L3TB(source)%nc        
        Tb1d(t) = Tbh_field_A(c,r)
        if(Tb1d(t).ne.-9999.0) then 
           li(t) = .true.
        endif
        t = t+1
     enddo
  enddo
  
  udef = -9999.0

  call neighbor_interp(LVT_rc%gridDesc, li, Tb1d, lo, Tbh_A_ip, &
       SMAP_L3TB(source)%nc*SMAP_L3TB(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMAP_L3TB(source)%rlat2, SMAP_L3TB(source)%rlon2,&
       SMAP_L3TB(source)%n112,udef, iret)




  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        L3TB(c,r,1) = Tbv_D_ip(c+(r-1)*LVT_rc%lnc)
        L3TB(c,r,2) = Tbv_A_ip(c+(r-1)*LVT_rc%lnc)
        L3TB(c,r,3) = Tbh_D_ip(c+(r-1)*LVT_rc%lnc)
        L3TB(c,r,4) = Tbh_A_ip(c+(r-1)*LVT_rc%lnc)
     enddo
  enddo

  !deallocate(Tbv_field)
  deallocate(dims)
  deallocate(Tbv_field_D)
  deallocate(Tbv_field_A)
  deallocate(Tbh_field_D)
  deallocate(Tbh_field_A)
#endif

end subroutine read_SMAP_L3Tb

#if 0
!BOP
! 
! !ROUTINE: SMAP_L3TB_filename
! \label{SMAP_L3TB_filename}
!
! !INTERFACE: 
subroutine SMAP_L3TB_filename(source, name, designation, ndir, yr, mo,da)
! 
! !USES:   
  use LVT_coreMod,only : LVT_rc
  use LVT_logMod, only : LVT_logunit

  implicit none
!
! !ARGUMENTS: 
  integer            :: source
  character(len=*)   :: name
  character(len=*)   :: designation
  integer            :: yr, mo, da, hr,mn
  character (len=*)  :: ndir
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine creates the NASA AMSRE filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the NASA AMSRE soil moisture filename
!  \item[ndir] name of the NASA AMSRE soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') LVT_rc%dyr(source)
  write(unit=fmo, fmt='(i2.2)') LVT_rc%dmo(source)
  write(unit=fda, fmt='(i2.2)') LVT_rc%dda(source)
  
  if(designation.eq."SPL3SMAP") then 
     name = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
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
     name = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '.h5'

! For example:
! SMAP_L3_SM_P_E_20180811.h5 -> SMAP_L3_SM_P_E_20180811_R16010_001.h5
  elseif(designation.eq."SPL3SMP_E") then 
     name = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_E_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '.h5'
  endif

end subroutine SMAP_L3TB_filename
#endif
