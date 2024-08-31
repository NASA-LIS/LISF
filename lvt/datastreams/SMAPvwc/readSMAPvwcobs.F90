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
! !ROUTINE: readSMAPvwcobs
! \label{readSMAPvwcobs}
!
! !INTERFACE:
subroutine readSMAPvwcobs(source)
!
! !USES:
   use ESMF
   use LVT_coreMod, only: LVT_rc
   use LVT_histDataMod
   use LVT_logMod  !,       only : LVT_logunit
   use SMAP_vwcobsMod, only: SMAP_vwcobs

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
!  21 July 2010: Sujay Kumar, Initial Specification
!  28 Aug 2018: Mahdi Navari, Edited to read Vegetation water
!               content from SPL3SMP.005 & SPL3SMP_E.002
! 11 July 2019: Mahdi Navari, There are several version of SMAP sm data available in each directory
!                  with different Release number and different CRID Version Number. The reader was
!                  modified to read the latest version of data (the reader no longer reads the symbolic
!                  link to the SMAP sm data)
!
!EOP

   logical           :: alarmcheck, file_exists, readflag
   integer           :: iret
   character*200     :: fname
   real              :: vwc(LVT_rc%lnc, LVT_rc%lnr)
   integer           :: fnd
   real              :: timenow
   character*4       :: yyyy
   character*2       :: mm, dd, hh
   integer               :: yr, mo, da, hr, mn, ss
   integer               :: doy
   character*200      :: list_files
   integer               :: ftn, ierr
   character(len=3) :: CRID

   vwc = LVT_rc%udef

   timenow = float(LVT_rc%dhr(source))*3600 + &
             60*LVT_rc%dmn(source) + LVT_rc%dss(source)
   alarmcheck = (mod(timenow, 86400.0) .eq. 0)
   if (SMAP_vwcobs(source)%startflag .or. alarmCheck .or. &
       LVT_rc%resetFlag(source)) then
      LVT_rc%resetFlag(source) = .false.
      SMAP_vwcobs(source)%startflag = .false.

      if (SMAP_vwcobs(source)%data_designation .eq. "SPL3SMP_E") then
!----------------------------------------------------------------------------------------------------------------
! create filename for 9 km product
!----------------------------------------------------------------------------------------------------------------
         write (yyyy, '(i4.4)') LVT_rc%dyr(source)
         write (mm, '(i2.2)') LVT_rc%dmo(source)
         write (dd, '(i2.2)') LVT_rc%dda(source)
         write (CRID, '(a)') SMAP_vwcobs(source)%release_number

         list_files = 'ls '//trim(SMAP_vwcobs(source)%odir)//'/'&
              //trim(yyyy)//'.'//trim(mm)//'.'// &
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
            call read_SMAPsm(source, fname, vwc)
         enddo
         call LVT_releaseUnitNumber(ftn)

      elseif (SMAP_vwcobs(source)%data_designation .eq. "SPL3SMP") then
!----------------------------------------------------------------------------------------------------------------
! create filename for 36 km product
!----------------------------------------------------------------------------------------------------------------
         write (yyyy, '(i4.4)') LVT_rc%dyr(source)
         write (mm, '(i2.2)') LVT_rc%dmo(source)
         write (dd, '(i2.2)') LVT_rc%dda(source)
         write (CRID, '(a)') SMAP_vwcobs(source)%release_number

         list_files = 'ls '//trim(SMAP_vwcobs(source)%odir)//'/'&
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
            call read_SMAPsm(source, fname, vwc)
         enddo
         call LVT_releaseUnitNumber(ftn)

      endif   ! sensor
   endif ! alaram

   call LVT_logSingleDataStreamVar(LVT_MOC_VEGWATERCONTENT, source, &
                                   vwc, vlevel=1, units="kg/m2")

!  open(100,file='test.bin',form='unformatted')
!  write(100) vwc
!  close(100)
!  stop
end subroutine readSMAPvwcobs
!BOP
! 
! !ROUTINE: read_SMAPsm
! \label{read_SMAPsm}
!
! !INTERFACE: 
subroutine read_SMAPvwc(source, fname, vwcobs)
! 
! !USES: 
  use LVT_coreMod,         only : LVT_rc,LVT_domain
  use LVT_logMod
  use SMAP_vwcobsMod, only : SMAP_vwcobs

#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none

  integer                       :: source
  character(len=*)              :: fname
  real                          :: vwcobs(LVT_rc%lnc,LVT_rc%lnr)
  
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

  character*100,   parameter    :: vwc_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,   parameter    :: vwc_field_name = "soil_moisture"

  character*100,    parameter    :: vwc_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: vwc_field_name_D = "vegetation_water_content"
  character*100,    parameter    :: vwc_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: vwc_field_name_A = "vegetation_water_content"

  integer(hid_t)                :: file_id, vwc_gr_id,vwc_field_id
  integer(hid_t)                :: vwc_gr_id_D,vwc_field_id_D
  integer(hid_t)                :: vwc_gr_id_A,vwc_field_id_A
  integer(hid_t)                :: dataspace
  integer(hid_t)                :: memspace
  integer                       :: memrank = 2
  integer(hsize_t), allocatable :: dims(:)
  integer(hsize_t), dimension(2) :: dimvwc
  integer(hsize_t), dimension(2) :: offset_file
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  real,             allocatable  :: vwc_field(:,:)
  real,             allocatable  :: vwc_field_D(:,:)
  real,             allocatable  :: vwc_field_A(:,:)
  real                           :: vwc1d(SMAP_vwcobs(source)%nc*SMAP_vwcobs(source)%nr)
  real                           :: vwc_ip(LVT_rc%lnc*LVT_rc%lnr)
  logical*1                      :: li(SMAP_vwcobs(source)%nc*SMAP_vwcobs(source)%nr)
  logical*1                      :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: udef
  integer                        :: t, c,r,c1,r1
  integer                        :: iret,status
 
  dimvwc      = (/SMAP_vwcobs(source)%nc, SMAP_vwcobs(source)%nr/)
  count_file = (/SMAP_vwcobs(source)%nc, SMAP_vwcobs(source)%nr/)
  count_mem  = (/SMAP_vwcobs(source)%nc, SMAP_vwcobs(source)%nr/)
  
  allocate(vwc_field(SMAP_vwcobs(source)%nc, SMAP_vwcobs(source)%nr))
  allocate(vwc_field_D(SMAP_vwcobs(source)%nc, SMAP_vwcobs(source)%nr))
  allocate(vwc_field_A(SMAP_vwcobs(source)%nc, SMAP_vwcobs(source)%nr))
  allocate(dims(2))

  dims(1) = SMAP_vwcobs(source)%nc
  dims(2) = SMAP_vwcobs(source)%nr

  if ((SMAP_vwcobs(source)%data_designation.eq."SPL3SMP_E") .or. & 
      (SMAP_vwcobs(source)%data_designation.eq."SPL3SMP") )then 
     call h5open_f(status)
     call LVT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LVT_verify(status, 'Error opening NASASMAP file ')     

!Read the AM (descending) data     
     call h5gopen_f(file_id,vwc_gr_name_D,vwc_gr_id_D, status)
     call LVT_verify(status, 'Error opening SM group (AM) in NASASMAP file')

     call h5dopen_f(vwc_gr_id_D,vwc_field_name_D,vwc_field_id_D, status)
     call LVT_verify(status, 'Error opening Veg water content field in NASASMAP file')
     
     call h5dget_space_f(vwc_field_id_D, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimvwc, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vwc_field_id_D, H5T_NATIVE_REAL,vwc_field_D,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting Veg water content (AM) field from NASASMAPfile')

!Read the PM (ascending) data     
     call h5gopen_f(file_id,vwc_gr_name_A,vwc_gr_id_A, status)
     call LVT_verify(status, 'Error opening SM group (PM) in NASASMAP file')

     call h5dopen_f(vwc_gr_id_A,vwc_field_name_A,vwc_field_id_A, status)
     call LVT_verify(status, 'Error opening Veg water content field in NASASMAP file')
     
     call h5dget_space_f(vwc_field_id_A, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimvwc, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vwc_field_id_A, H5T_NATIVE_REAL,vwc_field_A,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting Veg water content (AM) field from NASASMAPfile')


     call h5dclose_f(vwc_field_id_D,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(vwc_field_id_A,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(vwc_gr_id_D,status)
     call LVT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(vwc_gr_id_A,status)
     call LVT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LVT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LVT_verify(status,'Error in H5CLOSE call')

     vwc_field = LVT_rc%udef
     !blend the AM and PM overpasses
     do r=1,SMAP_vwcobs(source)%nr
        do c=1,SMAP_vwcobs(source)%nc
           if(vwc_field_D(c,r).ne.LVT_rc%udef) then
              vwc_field(c,r) = vwc_field_D(c,r)
           endif
        enddo
     enddo
     do r=1,SMAP_vwcobs(source)%nr
        do c=1,SMAP_vwcobs(source)%nc
           if(vwc_field_A(c,r).ne.LVT_rc%udef) then
              if(vwc_field(c,r).eq.LVT_rc%udef) then 
                 vwc_field(c,r) = vwc_field_A(c,r)
              endif
           endif
        enddo
     enddo
  else
     call h5open_f(status)
     call LVT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LVT_verify(status, 'Error opening NASASMAP file ')
     
     call h5gopen_f(file_id,vwc_gr_name_D,vwc_gr_id_D, status)
     call LVT_verify(status, 'Error opening SM group in NASASMAP file')

     call h5dopen_f(vwc_gr_id_D,vwc_field_name_D,vwc_field_id_D, status)
     call LVT_verify(status, 'Error opening Veg water content field in NASASMAP file')
     
     call h5dget_space_f(vwc_field_id_D, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimvwc, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vwc_field_id_D, H5T_NATIVE_REAL,vwc_field_D,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting Veg water content field from NASASMAPfile')

     call h5gopen_f(file_id,vwc_gr_name_A,vwc_gr_id_A, status)
     call LVT_verify(status, 'Error opening SM group in NASASMAP file')

     call h5dopen_f(vwc_gr_id_A,vwc_field_name_A,vwc_field_id_A, status)
     call LVT_verify(status, 'Error opening Veg water content field in NASASMAP file')
     
     call h5dget_space_f(vwc_field_id_A, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimvwc, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vwc_field_id_A, H5T_NATIVE_REAL,vwc_field_A,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting Veg water content field from NASASMAPfile')
     
     call h5dclose_f(vwc_field_id_D,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(vwc_field_id_A,status)
     call LVT_verify(status,'Error in H5DCLOSE call')
     
     call h5gclose_f(vwc_gr_id_D,status)
     call LVT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(vwc_gr_id_A,status)
     call LVT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LVT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LVT_verify(status,'Error in H5CLOSE call')
     
     vwc_field = LVT_rc%udef
     !blend the AM and PM overpasses
     do r=1,SMAP_vwcobs(source)%nr
        do c=1,SMAP_vwcobs(source)%nc
           if(vwc_field_D(c,r).ne.LVT_rc%udef) then
              vwc_field(c,r) = vwc_field_D(c,r)
           endif
        enddo
     enddo
     do r=1,SMAP_vwcobs(source)%nr
        do c=1,SMAP_vwcobs(source)%nc
           if(vwc_field_A(c,r).ne.LVT_rc%udef) then
              if(vwc_field(c,r).eq.LVT_rc%udef) then 
                 vwc_field(c,r) = vwc_field_A(c,r)
              endif
           endif
        enddo
     enddo
  endif

  li = .false. 
  t = 1

  do r=1,SMAP_vwcobs(source)%nr
     do c=1,SMAP_vwcobs(source)%nc        
        vwc1d(t) = vwc_field(c,r)
        if(vwc1d(t).ne.-9999.0) then 
           li(t) = .true.
        endif
        t = t+1
     enddo
  enddo
  
  udef = -9999.0

  call neighbor_interp(LVT_rc%gridDesc, li, vwc1d, lo, vwc_ip, &
       SMAP_vwcobs(source)%nc*SMAP_vwcobs(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMAP_vwcobs(source)%rlat2, SMAP_vwcobs(source)%rlon2,&
       SMAP_vwcobs(source)%n112,udef, iret)

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        vwcobs(c,r) = vwc_ip(c+(r-1)*LVT_rc%lnc)
     enddo
  enddo

  deallocate(vwc_field)
  deallocate(dims)

#endif

end subroutine read_SMAPvwc

#if 0

!BOP
! 
! !ROUTINE: SMAP_vwc_filename
! \label{SMAP_vwc_filename}
!
! !INTERFACE: 
subroutine SMAP_vwc_filename(source, name, designation, ndir, yr, mo,da)
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
!  This subroutine creates the NASA SMAP filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the NASA SMAP filename
!  \item[ndir] name of the NASA SMAP directory
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
! is difficult to programatically generate.  So after downloading
! a SMAP data file, a symbolic link was created to it which make it 
! easier to generate file name.
!   For example:
!   SMAP_L3_SM_P_20170902.h5 -> SMAP_L3_SM_P_20170902_R15152_001.h5
  elseif(designation.eq."SPL3SMP") then 
     name = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '.h5'
!          '_R14010_001.h5'
! For example:
! SMAP_L3_SM_P_E_20180811.h5 -> SMAP_L3_SM_P_E_20180811_R16010_001.h5
  elseif(designation.eq."SPL3SMP_E") then 
     name = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_E_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '.h5'
  endif

end subroutine SMAP_vwc_filename
#endif
