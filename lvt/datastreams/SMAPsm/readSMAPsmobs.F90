!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readSMAPsmobs
! \label{readSMAPsmobs}
!
! !INTERFACE: 
subroutine readSMAPsmobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit
  use SMAP_smobsMod, only : SMAP_smobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in)       :: source
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
!  17 Aug 2018: Mahdi Navari, Edited to read SPL3SMP.005 & SPL3SMP_E.002 
!EOP

  logical           :: alarmcheck, file_exists, readflag
  integer           :: iret
  character*200     :: name
  real              :: smc(LVT_rc%lnc, LVT_rc%lnr)
  integer           :: fnd 
  real              :: timenow

  smc = LVT_rc%udef

  timenow = float(LVT_rc%dhr(source))*3600 +&
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  if(SMAP_smobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 

     SMAP_smobs(source)%startflag = .false. 
     call SMAP_sm_filename(source,name,&
          SMAP_smobs(source)%data_designation, & 
          SMAP_smobs(source)%odir, & 
        LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source))
                 
     inquire(file=name, exist=file_exists) 

     if(file_exists) then 
        readflag = .true. 
     else
        readflag = .false. 
     endif
     
     if(readflag) then 
        write(LVT_logunit,*) '[INFO] Reading SMAP file ',name
        call read_SMAPsm(source, name, smc)
        
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source,&
       smc,vlevel=1,units="m3/m3")
 
!  open(100,file='test.bin',form='unformatted')
!  write(100) smc
!  close(100)
!  stop
end subroutine readSMAPsmobs

!BOP
! 
! !ROUTINE: read_SMAPsm
! \label{read_SMAPsm}
!
! !INTERFACE: 
subroutine read_SMAPsm(source, fname, smobs)
! 
! !USES: 
  use LVT_coreMod,         only : LVT_rc,LVT_domain
  use LVT_logMod
  use SMAP_smobsMod, only : SMAP_smobs

#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none

  real                           :: smobs(LVT_rc%lnc,LVT_rc%lnr)
  
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

  integer                       :: source 
  character(len=*)              :: fname
  character*100,   parameter    :: sm_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,   parameter    :: sm_field_name = "soil_moisture"

  character*100,    parameter    :: sm_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: sm_field_name_D = "soil_moisture"
  character*100,    parameter    :: sm_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: sm_field_name_A = "soil_moisture_pm"

  integer(hid_t)                :: file_id, sm_gr_id,sm_field_id
  integer(hid_t)                :: sm_gr_id_D,sm_field_id_D
  integer(hid_t)                :: sm_gr_id_A,sm_field_id_A
  integer(hid_t)                :: dataspace
  integer(hid_t)                :: memspace
  integer                       :: memrank = 2
  integer(hsize_t), allocatable :: dims(:)
  integer(hsize_t), dimension(2) :: dimsm
  integer(hsize_t), dimension(2) :: offset_file = (/0,0/)
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  real,             allocatable  :: sm_field(:,:)
  real,             allocatable  :: sm_field_D(:,:)
  real,             allocatable  :: sm_field_A(:,:)
  real                           :: sm1d(SMAP_smobs(source)%nc*SMAP_smobs(source)%nr)
  real                           :: sm_ip(LVT_rc%lnc*LVT_rc%lnr)
  logical*1                      :: li(SMAP_smobs(source)%nc*SMAP_smobs(source)%nr)
  logical*1                      :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: udef
  integer                        :: t, c,r,c1,r1
  integer                        :: iret,status
 
  dimsm      = (/SMAP_smobs(source)%nc, SMAP_smobs(source)%nr/)
  count_file = (/SMAP_smobs(source)%nc, SMAP_smobs(source)%nr/)
  count_mem  = (/SMAP_smobs(source)%nc, SMAP_smobs(source)%nr/)
  
  allocate(sm_field(SMAP_smobs(source)%nc, SMAP_smobs(source)%nr))
  allocate(sm_field_D(SMAP_smobs(source)%nc, SMAP_smobs(source)%nr))
  allocate(sm_field_A(SMAP_smobs(source)%nc, SMAP_smobs(source)%nr))
  allocate(dims(2))

  dims(1) = SMAP_smobs(source)%nc
  dims(2) = SMAP_smobs(source)%nr

  if ((SMAP_smobs(source)%data_designation.eq."SPL3SMP_E") .or. & 
      (SMAP_smobs(source)%data_designation.eq."SPL3SMP") ) then 
     call h5open_f(status)
     call LVT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LVT_verify(status, 'Error opening NASASMAP file ')     

!Read the AM (descending) data     
     call h5gopen_f(file_id,sm_gr_name_D,sm_gr_id_D, status)
     call LVT_verify(status, 'Error opening SM group (AM) in NASASMAP file')

     call h5dopen_f(sm_gr_id_D,sm_field_name_D,sm_field_id_D, status)
     call LVT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id_D, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(sm_field_id_D, H5T_NATIVE_REAL,sm_field_D,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting SM (AM) field from NASASMAPfile')

!Read the PM (ascending) data     
     call h5gopen_f(file_id,sm_gr_name_A,sm_gr_id_A, status)
     call LVT_verify(status, 'Error opening SM group (PM) in NASASMAP file')

     call h5dopen_f(sm_gr_id_A,sm_field_name_A,sm_field_id_A, status)
     call LVT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id_A, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(sm_field_id_A, H5T_NATIVE_REAL,sm_field_A,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting SM (AM) field from NASASMAPfile')


     call h5dclose_f(sm_field_id_D,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(sm_field_id_A,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(sm_gr_id_D,status)
     call LVT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(sm_gr_id_A,status)
     call LVT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LVT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LVT_verify(status,'Error in H5CLOSE call')

     sm_field = LVT_rc%udef
     !blend the AM and PM overpasses
     do r=1,SMAP_smobs(source)%nr
        do c=1,SMAP_smobs(source)%nc
           if(sm_field_D(c,r).ne.LVT_rc%udef) then
              sm_field(c,r) = sm_field_D(c,r)
           endif
        enddo
     enddo
     do r=1,SMAP_smobs(source)%nr
        do c=1,SMAP_smobs(source)%nc
           if(sm_field_A(c,r).ne.LVT_rc%udef) then
              if(sm_field(c,r).eq.LVT_rc%udef) then 
                 sm_field(c,r) = sm_field_A(c,r)
              endif
           endif
        enddo
     enddo
  else
     call h5open_f(status)
     call LVT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LVT_verify(status, 'Error opening NASASMAP file ')
     
     call h5gopen_f(file_id,sm_gr_name_D,sm_gr_id_D, status)
     call LVT_verify(status, 'Error opening SM group in NASASMAP file')

     call h5dopen_f(sm_gr_id_D,sm_field_name_D,sm_field_id_D, status)
     call LVT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id_D, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(sm_field_id_D, H5T_NATIVE_REAL,sm_field_D,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting SM field from NASASMAPfile')

     call h5gopen_f(file_id,sm_gr_name_A,sm_gr_id_A, status)
     call LVT_verify(status, 'Error opening SM group in NASASMAP file')

     call h5dopen_f(sm_gr_id_A,sm_field_name_A,sm_field_id_A, status)
     call LVT_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id_A, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(sm_field_id_A, H5T_NATIVE_REAL,sm_field_A,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting SM field from NASASMAPfile')
     
     call h5dclose_f(sm_field_id_D,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(sm_field_id_A,status)
     call LVT_verify(status,'Error in H5DCLOSE call')
     
     call h5gclose_f(sm_gr_id_D,status)
     call LVT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(sm_gr_id_A,status)
     call LVT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LVT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LVT_verify(status,'Error in H5CLOSE call')
     
     sm_field = LVT_rc%udef
     !blend the AM and PM overpasses
     do r=1,SMAP_smobs(source)%nr
        do c=1,SMAP_smobs(source)%nc
           if(sm_field_D(c,r).ne.LVT_rc%udef) then
              sm_field(c,r) = sm_field_D(c,r)
           endif
        enddo
     enddo
     do r=1,SMAP_smobs(source)%nr
        do c=1,SMAP_smobs(source)%nc
           if(sm_field_A(c,r).ne.LVT_rc%udef) then
              if(sm_field(c,r).eq.LVT_rc%udef) then 
                 sm_field(c,r) = sm_field_A(c,r)
              endif
           endif
        enddo
     enddo
  endif

  li = .false. 
  t = 1

  do r=1,SMAP_smobs(source)%nr
     do c=1,SMAP_smobs(source)%nc        
        sm1d(t) = sm_field(c,r)
        if(sm1d(t).ne.-9999.0) then 
           li(t) = .true.
        endif
        t = t+1
     enddo
  enddo
  
  udef = -9999.0

  call neighbor_interp(LVT_rc%gridDesc, li, sm1d, lo, sm_ip, &
       SMAP_smobs(source)%nc*SMAP_smobs(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMAP_smobs(source)%rlat2, SMAP_smobs(source)%rlon2,&
       SMAP_smobs(source)%n112,udef, iret)

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        smobs(c,r) = sm_ip(c+(r-1)*LVT_rc%lnc)
     enddo
  enddo

  deallocate(sm_field)
  deallocate(dims)

#endif

end subroutine read_SMAPsm

!BOP
! 
! !ROUTINE: SMAP_sm_filename
! \label{SMAP_sm_filename}
!
! !INTERFACE: 
subroutine SMAP_sm_filename(source, name, designation, ndir, yr, mo,da)
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
!          '_R14010_001.h5'
! For example:
! SMAP_L3_SM_P_E_20180811.h5 -> SMAP_L3_SM_P_E_20180811_R16010_001.h5
  elseif(designation.eq."SPL3SMP_E") then 
     name = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_E_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '.h5'
  endif

end subroutine SMAP_sm_filename
