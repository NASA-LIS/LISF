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
  use LVT_logMod
  use LVT_timeMgrMod
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
!  19 Jun 2019: Sujay Kumar; Added support for SMAP L2 data
!EOP

  logical           :: alarmcheck, file_exists, readflag
  integer           :: iret
  character*200     :: name
  real              :: smc(LVT_rc%lnc, LVT_rc%lnr)
  integer           :: fnd 
  character*8       :: yyyymmdd
  character*4       :: yyyy
  character*2       :: mm,dd,hh
  character*200     :: list_files
  integer           :: i,ftn,ierr
  character*200     :: smap_filename(10)
  character*200     :: fname
  integer           :: mn_ind
  integer           :: yr, mo, da, hr, mn, ss
  real              :: gmt
  integer           :: c,r
  integer           :: doy
  real*8            :: timenow
  real              :: smobs(LVT_rc%lnc,LVT_rc%lnr)

  smc = LVT_rc%udef
  smobs = LVT_rc%udef

  timenow = float(LVT_rc%dhr(source))*3600 +&
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)

  if((SMAP_smobs(source)%data_designation.eq."SPL2SMP_E").or.&
       (SMAP_smobs(source)%data_designation.eq."SPL2SMP")) then     
     alarmcheck = (mod(timenow, 3600.0).eq.0)
  else
     alarmcheck = (mod(timenow, 86400.0).eq.0)
  endif

  if(SMAP_smobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 

     SMAP_smobs(source)%startflag = .false. 

     if((SMAP_smobs(source)%data_designation.eq."SPL2SMP_E").or.&
          (SMAP_smobs(source)%data_designation.eq."SPL2SMP")) then  

        write(yyyymmdd,'(i4.4,2i2.2)') LVT_rc%yr, LVT_rc%mo, LVT_rc%da
        write(yyyy,'(i4.4)') LVT_rc%yr
        write(mm,'(i2.2)') LVT_rc%mo
        write(dd,'(i2.2)') LVT_rc%da
        write(hh,'(i2.2)') LVT_rc%hr
        
        list_files = 'ls '//trim(SMAP_smobs(source)%odir)//&
             '/'//trim(yyyy)//'.'//trim(mm)//'.'//dd//&
             '/SMAP_L2_*'//trim(yyyymmdd)//'T'//trim(hh)&
             //"*.h5 > SMAP_filelist.dat"
        
        call system(trim(list_files))

        i =1
        ftn = LVT_getNextUnitNumber()
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
           call LVT_tick(timenow,doy,gmt,LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
                LVT_rc%hr, mn, ss, 0)
        
           smap_filename(i) = fname
           
           write(LVT_logunit,*) '[INFO] reading ',trim(smap_filename(i))
           
           call read_SMAPL2sm_data(source,smap_filename(i),smobs,timenow)

           i = i+1
        enddo
        call LVT_releaseUnitNumber(ftn)
        
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              smc(c,r) = smobs(c,r)
           enddo
        enddo

     else
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
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source,&
       smc,vlevel=1,units="m3/m3")

!  open(100,file='test_out.bin',form='unformatted')
!  write(100) smc
!  close(100)
!  stop
 
end subroutine readSMAPsmobs

!BOP
! 
! !ROUTINE: read_SMAPL2sm_data
! \label{read_SMAPL2sm_data}
!
! !INTERFACE:
subroutine read_SMAPL2sm_data(source, fname, smobs_inp, time)
! 
! !USES:   

  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SMAP_smobsMod, only : SMAP_smobs
#if (defined USE_HDF5) 
  use hdf5
#endif

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
  logical*1                      :: sm_data_b(SMAP_smobs(source)%nc*SMAP_smobs(source)%nr)
  logical*1                      :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: sm_data(SMAP_smobs(source)%nc*SMAP_smobs(source)%nr)
  real                           :: smobs_ip(LVT_rc%lnc*LVT_rc%lnr)

  integer                        :: status,ios,iret

  smobs_inp = LVT_rc%udef

  call h5open_f(status)
  call LVT_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LVT_verify(status, 'Error opening SMAP L2 file ')
  
  call h5gopen_f(file_id,sm_gr_name,sm_gr_id, status)
  call LVT_verify(status, 'Error opening SM group in SMAP L2 file')
  
  call h5dopen_f(sm_gr_id,sm_field_name,sm_field_id, status)
  call LVT_verify(status, 'Error opening SM field in SMAP L2 file')

  call h5dopen_f(sm_gr_id,"EASE_row_index",row_id, status)
  call LVT_verify(status, 'Error opening row index field in SMAP L2 file')

  call h5dopen_f(sm_gr_id,"EASE_column_index",col_id, status)
  call LVT_verify(status, 'Error opening column index field in SMAP L2 file')

  call h5dopen_f(sm_gr_id, sm_qa_name,sm_qa_id, status)
  call LVT_verify(status, 'Error opening QA field in SMAP L2 file')
  
  call h5dget_space_f(sm_field_id, dspace_id, status)
  call LVT_verify(status, 'Error in h5dget_space_f: readSMAP L2Obs')
  
! Size of the arrays
! This routine returns -1 on failure, rank on success. 
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, status) 
  if(status.eq.-1) then 
     call LVT_verify(status, 'Error in h5sget_simple_extent_dims_f: readSMAP L2Obs')
  endif
  
  allocate(sm_field(maxdims(1)))
  allocate(sm_qa(maxdims(1)))
  allocate(ease_row(maxdims(1)))
  allocate(ease_col(maxdims(1)))

  call h5dread_f(row_id, H5T_NATIVE_INTEGER,ease_row,dims,status)
  call LVT_verify(status, 'Error extracting row index from SMAP L2 file')

  call h5dread_f(col_id, H5T_NATIVE_INTEGER,ease_col,dims,status)
  call LVT_verify(status, 'Error extracting col index from SMAP L2 file')
  
  call h5dread_f(sm_field_id, H5T_NATIVE_REAL,sm_field,dims,status)
  call LVT_verify(status, 'Error extracting SM field from SMAP L2 file')

  call h5dread_f(sm_qa_id, H5T_NATIVE_INTEGER,sm_qa,dims,status)
  call LVT_verify(status, 'Error extracting SM field from SMAP L2 file')
  
  call h5dclose_f(sm_qa_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(row_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(col_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(sm_field_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')
  
  call h5gclose_f(sm_gr_id,status)
  call LVT_verify(status,'Error in H5GCLOSE call')
    
  call h5fclose_f(file_id,status)
  call LVT_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LVT_verify(status,'Error in H5CLOSE call')

  sm_data = LVT_rc%udef
  sm_data_b = .false. 

!grid the data in EASE projection
  do t=1,maxdims(1)
!     if(ibits(sm_qa(t),0,1).eq.0) then 
     if(ease_col(t).gt.0.and.ease_row(t).gt.0) then 
        sm_data(ease_col(t) + &
             (ease_row(t)-1)*SMAP_smobs(source)%nc) = sm_field(t) 
        if(sm_field(t).ne.-9999.0) then 
           sm_data_b(ease_col(t) + &
                (ease_row(t)-1)*SMAP_smobs(source)%nc) = .true. 
        endif
     endif
  enddo
  
  t = 1

!  open(100,file='test_inp.bin',form='unformatted')
!  write(100) sm_data
!  close(100)
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LVT_rc%gridDesc, sm_data_b, sm_data, &
       smobs_b_ip, smobs_ip, &
       SMAP_smobs(source)%nc*SMAP_smobs(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMAP_smobs(source)%rlat2, SMAP_smobs(source)%rlon2,&
       SMAP_smobs(source)%n112,LVT_rc%udef, iret)

!  open(100,file='test_out.bin',form='unformatted')
!  write(100) smobs_ip
!  close(100)
!  stop

  deallocate(sm_field)
  deallocate(sm_qa)
  deallocate(ease_row)
  deallocate(ease_col)

!overwrite the input data 
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(smobs_ip(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
           smobs_inp(c,r) = & 
                smobs_ip(c+(r-1)*LVT_rc%lnc)

           SMAP_smobs(source)%smtime(c,r) = & 
                time
        endif
     enddo
  enddo

#endif

end subroutine read_SMAPL2sm_data

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
