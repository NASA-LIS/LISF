!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
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
  use LDT_logMod
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

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character*100     :: fname
  real              :: vodobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  NASASMAPvodobs(n)%vodobs = LDT_rc%udef
  vodobs= LDT_rc%udef

  call create_NASASMAPvod_filename(NASASMAPvodobs(n)%odir, &
       NASASMAPvodobs(n)%data_designation,&
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)
  
  inquire(file=trim(fname),exist=file_exists)
  if(file_exists) then
     
     write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
     call read_SMAP_vod_data(n, fname, vodobs)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(vodobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
              NASASMAPvodobs(n)%vodobs(c,r) = vodobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%vod_obs,&
       NASASMAPvodobs(n)%vodobs,vlevel=1)

end subroutine readNASASMAPvodObs


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
