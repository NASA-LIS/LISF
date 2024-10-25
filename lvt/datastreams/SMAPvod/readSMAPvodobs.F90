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
! !ROUTINE: readSMAPvodobs
! \label{readSMAPvodobs}
!
! !INTERFACE: 
subroutine readSMAPvodobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SMAP_vodobsMod, only : SMAP_vodobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in)       :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the 
! vegetation optical depth retrievals from the 
! NASA soil moisture retrieval product.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  10 Apr 2019: Sujay Kumar, Initial Specification
!
!EOP

  logical           :: alarmcheck, file_exists, readflag
  integer           :: iret
  character*200     :: name
  real              :: vod(LVT_rc%lnc, LVT_rc%lnr)
  integer           :: fnd 
  character*8            :: yyyymmdd
  character*4            :: yyyy
  character*2            :: mm,dd,hh
  character*200          :: list_files
  integer                :: i,ftn,ierr
  character*200          :: smap_filename(10)
  character*200          :: fname
  integer                :: mn_ind
  integer                :: yr, mo, da, hr, mn, ss
  real                   :: gmt
  integer                :: c,r
  integer                :: doy
  real*8                 :: timenow
  real                   :: vodobs(LVT_rc%lnc,LVT_rc%lnr)

  vod = LVT_rc%udef
  vodobs = LVT_rc%udef

  timenow = float(LVT_rc%dhr(source))*3600 +&
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  if((SMAP_vodobs(source)%data_designation.eq."SPL2SMP_E").or.&
       (SMAP_vodobs(source)%data_designation.eq."SPL2SMP")) then 
     alarmcheck = (mod(timenow, 3600.0).eq.0)
  else
     alarmcheck = (mod(timenow, 86400.0).eq.0)
  endif

  if(SMAP_vodobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 

     SMAP_vodobs(source)%startflag = .false. 
     
     if((SMAP_vodobs(source)%data_designation.eq."SPL2SMP_E").or.&
          (SMAP_vodobs(source)%data_designation.eq."SPL2SMP")) then 
        
        write(yyyymmdd,'(i4.4,2i2.2)') LVT_rc%dyr(source), &
             LVT_rc%dmo(source), LVT_rc%dda(source)
        write(yyyy,'(i4.4)') LVT_rc%dyr(source)
        write(mm,'(i2.2)') LVT_rc%dmo(source)
        write(dd,'(i2.2)') LVT_rc%dda(source)
        write(hh,'(i2.2)') LVT_rc%dhr(source)

        list_files = 'ls '//trim(SMAP_vodobs(source)%odir)//&
             '/'//trim(yyyy)//'.'//trim(mm)//'.'//dd//&
             '/SMAP_L2_*'//trim(yyyymmdd)//'T'//trim(hh)&
             //"*.h5 > SMAPvod/SMAP_filelist.dat"
        
        call system(trim(list_files))

        i =1
        ftn = LVT_getNextUnitNumber()
        open(ftn,file="./SMAPvod/SMAP_filelist.dat",&
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
           call LVT_tick(timenow,doy,gmt,LVT_rc%dyr(source), &
                LVT_rc%dmo(source), LVT_rc%dda(source), &
                LVT_rc%dhr(source), mn, ss, 0)
        
           smap_filename(i) = fname
           
           write(LVT_logunit,*) '[INFO] reading ',trim(smap_filename(i))
           
           call read_SMAPL2vod_data(source,smap_filename(i),vodobs,timenow)

           i = i+1
        enddo
        call LVT_releaseUnitNumber(ftn)
        
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              vod(c,r) = vodobs(c,r)
           enddo
        enddo

     else
        call SMAP_vod_filename(source,name,&
             SMAP_vodobs(source)%data_designation, & 
             SMAP_vodobs(source)%odir, & 
             LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source))
        
        inquire(file=name, exist=file_exists) 
        
        if(file_exists) then 
           readflag = .true. 
        else
           readflag = .false. 
        endif
        
        if(readflag) then 
           write(LVT_logunit,*) '[INFO] Reading SMAP file ',name
           call read_SMAPvod(source, name, vod)
           
        endif
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_VOD, source,&
       vod,vlevel=1,units="-")
 
!  open(100,file='test_out.bin',form='unformatted')
!  write(100) vod
!  close(100)
!  stop

end subroutine readSMAPvodobs

!BOP
! 
! !ROUTINE: read_SMAPL2vod_data
! \label{read_SMAPL2vod_data}
!
! !INTERFACE:
subroutine read_SMAPL2vod_data(source, fname, vodobs_inp, time)
! 
! !USES:   

  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SMAP_vodobsMod, only : SMAP_vodobs
#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: source
  character (len=*)        :: fname
  real                     :: vodobs_inp(LVT_rc%lnc,LVT_rc%lnr)
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
  logical*1                      :: vod_data_b(SMAP_vodobs(source)%nc*SMAP_vodobs(source)%nr)
  logical*1                      :: vodobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: vod_data(SMAP_vodobs(source)%nc*SMAP_vodobs(source)%nr)
  real                           :: vodobs_ip(LVT_rc%lnc*LVT_rc%lnr)

  integer                        :: status,ios,iret

  call h5open_f(status)
  call LVT_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LVT_verify(status, 'Error opening SMAP L2 file ')
  
  call h5gopen_f(file_id,vod_gr_name,vod_gr_id, status)
  call LVT_verify(status, 'Error opening SM group in SMAP L2 file')
  
  call h5dopen_f(vod_gr_id,vod_field_name,vod_field_id, status)
  call LVT_verify(status, 'Error opening SM field in SMAP L2 file')

  call h5dopen_f(vod_gr_id,"EASE_row_index",row_id, status)
  call LVT_verify(status, 'Error opening row index field in SMAP L2 file')

  call h5dopen_f(vod_gr_id,"EASE_column_index",col_id, status)
  call LVT_verify(status, 'Error opening column index field in SMAP L2 file')

  call h5dopen_f(vod_gr_id, vod_qa_name,vod_qa_id, status)
  call LVT_verify(status, 'Error opening QA field in SMAP L2 file')
  
  call h5dget_space_f(vod_field_id, dspace_id, status)
  call LVT_verify(status, 'Error in h5dget_space_f: readSMAP L2Obs')
  
! Size of the arrays
! This routine returns -1 on failure, rank on success. 
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, status) 
  if(status.eq.-1) then 
     call LVT_verify(status, 'Error in h5sget_simple_extent_dims_f: readSMAP L2Obs')
  endif
  
  allocate(vod_field(maxdims(1)))
  allocate(vod_qa(maxdims(1)))
  allocate(ease_row(maxdims(1)))
  allocate(ease_col(maxdims(1)))

  call h5dread_f(row_id, H5T_NATIVE_INTEGER,ease_row,dims,status)
  call LVT_verify(status, 'Error extracting row index from SMAP L2 file')

  call h5dread_f(col_id, H5T_NATIVE_INTEGER,ease_col,dims,status)
  call LVT_verify(status, 'Error extracting col index from SMAP L2 file')
  
  call h5dread_f(vod_field_id, H5T_NATIVE_REAL,vod_field,dims,status)
  call LVT_verify(status, 'Error extracting SM field from SMAP L2 file')

  call h5dread_f(vod_qa_id, H5T_NATIVE_INTEGER,vod_qa,dims,status)
  call LVT_verify(status, 'Error extracting SM field from SMAP L2 file')
  
  call h5dclose_f(vod_qa_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(row_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(col_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(vod_field_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')
  
  call h5gclose_f(vod_gr_id,status)
  call LVT_verify(status,'Error in H5GCLOSE call')
    
  call h5fclose_f(file_id,status)
  call LVT_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LVT_verify(status,'Error in H5CLOSE call')

  vod_data = LVT_rc%udef
  vod_data_b = .false. 

!grid the data in EASE projection
  do t=1,maxdims(1)
!     if(ibits(vod_qa(t),0,1).eq.0) then 
     if(ease_col(t).gt.0.and.ease_row(t).gt.0) then 
        vod_data(ease_col(t) + &
             (ease_row(t)-1)*SMAP_vodobs(source)%nc) = vod_field(t) 
        if(vod_field(t).ne.-9999.0) then 
           vod_data_b(ease_col(t) + &
                (ease_row(t)-1)*SMAP_vodobs(source)%nc) = .true. 
        endif
     endif
  enddo
  
  t = 1

!  open(100,file='test_inp.bin',form='unformatted')
!  write(100) vod_data
!  close(100)

!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
!  if(LVT_isAtAfinerResolution(SMAP_vodobs(source)%gridDesci(10))) then  
     call neighbor_interp(LVT_rc%gridDesc, vod_data_b, vod_data, &
          vodobs_b_ip, vodobs_ip, &
          SMAP_vodobs(source)%nc*SMAP_vodobs(source)%nr,&
          LVT_rc%lnc*LVT_rc%lnr,&
          SMAP_vodobs(source)%rlat2, SMAP_vodobs(source)%rlon2,&
          SMAP_vodobs(source)%n112,LVT_rc%udef, iret)
!  else
!     call upscaleByAveraging(&
!          SMAP_vodobs(source)%nc*SMAP_vodobs(source)%nr,&
!          LVT_rc%lnc*LVT_rc%lnr,&
!          LVT_rc%udef, &
!          SMAP_vodobs(source)%n112,&
!          vod_data_b, &
!          vod_data, &
!          vodobs_b_ip,&
!          vodobs_ip)
!  endif
!  open(100,file='vodobs.bin',form='unformatted')
!  write(100) vodobs_ip
!  close(100)
!  stop

  deallocate(vod_field)
  deallocate(vod_qa)
  deallocate(ease_row)
  deallocate(ease_col)

!overwrite the input data 
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(vodobs_ip(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
           vodobs_inp(c,r) = & 
                vodobs_ip(c+(r-1)*LVT_rc%lnc)

           SMAP_vodobs(source)%vodtime(c,r) = & 
                time
        endif
     enddo
  enddo

#endif

end subroutine read_SMAPL2vod_data

!BOP
! 
! !ROUTINE: read_SMAPvod
! \label{read_SMAPvod}
!
! !INTERFACE: 
subroutine read_SMAPvod(source, fname, vodobs)
! 
! !USES: 
  use LVT_coreMod,         only : LVT_rc,LVT_domain
  use LVT_logMod
  use SMAP_vodobsMod, only : SMAP_vodobs

#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none

  integer                        :: source 
  character(len=*)               :: fname
  real                           :: vodobs(LVT_rc%lnc,LVT_rc%lnr)
  
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
 
  character*100,   parameter    :: vod_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,   parameter    :: vod_field_name = "vegetation_opacity"

  character*100,    parameter    :: vod_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: vod_field_name_D = "vegetation_opacity"
  character*100,    parameter    :: vod_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: vod_field_name_A = "vegetation_opacity_pm"

  integer(hid_t)                :: file_id, vod_gr_id,vod_field_id
  integer(hid_t)                :: vod_gr_id_D,vod_field_id_D
  integer(hid_t)                :: vod_gr_id_A,vod_field_id_A
  integer(hid_t)                :: dataspace
  integer(hid_t)                :: memspace
  integer                       :: memrank = 2
  integer(hsize_t), allocatable :: dims(:)
  integer(hsize_t), dimension(2) :: dimvod
  integer(hsize_t), dimension(2) :: offset_file
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  real,             allocatable  :: vod_field(:,:)
  real,             allocatable  :: vod_field_D(:,:)
  real,             allocatable  :: vod_field_A(:,:)
  real                           :: vod1d(SMAP_vodobs(source)%nc*SMAP_vodobs(source)%nr)
  real                           :: vod_ip(LVT_rc%lnc*LVT_rc%lnr)
  logical*1                      :: li(SMAP_vodobs(source)%nc*SMAP_vodobs(source)%nr)
  logical*1                      :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: udef
  integer                        :: t, c,r,c1,r1
  integer                        :: iret,status
 
  dimvod      = (/SMAP_vodobs(source)%nc, SMAP_vodobs(source)%nr/)
  count_file = (/SMAP_vodobs(source)%nc, SMAP_vodobs(source)%nr/)
  count_mem  = (/SMAP_vodobs(source)%nc, SMAP_vodobs(source)%nr/)
  
  allocate(vod_field(SMAP_vodobs(source)%nc, SMAP_vodobs(source)%nr))
  allocate(vod_field_D(SMAP_vodobs(source)%nc, SMAP_vodobs(source)%nr))
  allocate(vod_field_A(SMAP_vodobs(source)%nc, SMAP_vodobs(source)%nr))
  allocate(dims(2))

  dims(1) = SMAP_vodobs(source)%nc
  dims(2) = SMAP_vodobs(source)%nr

  if ((SMAP_vodobs(source)%data_designation.eq."SPL3SMP_E") .or. & 
      (SMAP_vodobs(source)%data_designation.eq."SPL3SMP") )then 
     call h5open_f(status)
     call LVT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LVT_verify(status, 'Error opening NASASMAP file ')     

!Read the AM (descending) data     
     call h5gopen_f(file_id,vod_gr_name_D,vod_gr_id_D, status)
     call LVT_verify(status, 'Error opening SM group (AM) in NASASMAP file')

     call h5dopen_f(vod_gr_id_D,vod_field_name_D,vod_field_id_D, status)
     call LVT_verify(status, 'Error opening Veg optical depth field in NASASMAP file')
     
     call h5dget_space_f(vod_field_id_D, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimvod, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vod_field_id_D, H5T_NATIVE_REAL,vod_field_D,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting Veg optical depth (AM) field from NASASMAPfile')

!Read the PM (ascending) data     
     call h5gopen_f(file_id,vod_gr_name_A,vod_gr_id_A, status)
     call LVT_verify(status, 'Error opening SM group (PM) in NASASMAP file')

     call h5dopen_f(vod_gr_id_A,vod_field_name_A,vod_field_id_A, status)
     call LVT_verify(status, 'Error opening Veg optical depth field in NASASMAP file')
     
     call h5dget_space_f(vod_field_id_A, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimvod, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vod_field_id_A, H5T_NATIVE_REAL,vod_field_A,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting Veg optical depth (AM) field from NASASMAPfile')


     call h5dclose_f(vod_field_id_D,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(vod_field_id_A,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(vod_gr_id_D,status)
     call LVT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(vod_gr_id_A,status)
     call LVT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LVT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LVT_verify(status,'Error in H5CLOSE call')

     vod_field = LVT_rc%udef
     !blend the AM and PM overpasses
     do r=1,SMAP_vodobs(source)%nr
        do c=1,SMAP_vodobs(source)%nc
           if(vod_field_D(c,r).ne.LVT_rc%udef) then
              vod_field(c,r) = vod_field_D(c,r)
           endif
        enddo
     enddo
     do r=1,SMAP_vodobs(source)%nr
        do c=1,SMAP_vodobs(source)%nc
           if(vod_field_A(c,r).ne.LVT_rc%udef) then
              if(vod_field(c,r).eq.LVT_rc%udef) then 
                 vod_field(c,r) = vod_field_A(c,r)
              endif
           endif
        enddo
     enddo
  else
     call h5open_f(status)
     call LVT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
     call LVT_verify(status, 'Error opening NASASMAP file ')
     
     call h5gopen_f(file_id,vod_gr_name_D,vod_gr_id_D, status)
     call LVT_verify(status, 'Error opening SM group in NASASMAP file')

     call h5dopen_f(vod_gr_id_D,vod_field_name_D,vod_field_id_D, status)
     call LVT_verify(status, 'Error opening Veg optical depth field in NASASMAP file')
     
     call h5dget_space_f(vod_field_id_D, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimvod, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vod_field_id_D, H5T_NATIVE_REAL,vod_field_D,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting Veg optical depth field from NASASMAPfile')

     call h5gopen_f(file_id,vod_gr_name_A,vod_gr_id_A, status)
     call LVT_verify(status, 'Error opening SM group in NASASMAP file')

     call h5dopen_f(vod_gr_id_A,vod_field_name_A,vod_field_id_A, status)
     call LVT_verify(status, 'Error opening Veg optical depth field in NASASMAP file')
     
     call h5dget_space_f(vod_field_id_A, dataspace, status)
     call LVT_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LVT_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimvod, memspace, status)
     call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')
     
     call h5dread_f(vod_field_id_A, H5T_NATIVE_REAL,vod_field_A,dims,status, &
          memspace, dataspace)
     call LVT_verify(status, 'Error extracting Veg optical depth field from NASASMAPfile')
     
     call h5dclose_f(vod_field_id_D,status)
     call LVT_verify(status,'Error in H5DCLOSE call')

     call h5dclose_f(vod_field_id_A,status)
     call LVT_verify(status,'Error in H5DCLOSE call')
     
     call h5gclose_f(vod_gr_id_D,status)
     call LVT_verify(status,'Error in H5GCLOSE call')

     call h5gclose_f(vod_gr_id_A,status)
     call LVT_verify(status,'Error in H5GCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LVT_verify(status,'Error in H5FCLOSE call')
     
     call h5close_f(status)
     call LVT_verify(status,'Error in H5CLOSE call')
     
     vod_field = LVT_rc%udef
     !blend the AM and PM overpasses
     do r=1,SMAP_vodobs(source)%nr
        do c=1,SMAP_vodobs(source)%nc
           if(vod_field_D(c,r).ne.LVT_rc%udef) then
              vod_field(c,r) = vod_field_D(c,r)
           endif
        enddo
     enddo
     do r=1,SMAP_vodobs(source)%nr
        do c=1,SMAP_vodobs(source)%nc
           if(vod_field_A(c,r).ne.LVT_rc%udef) then
              if(vod_field(c,r).eq.LVT_rc%udef) then 
                 vod_field(c,r) = vod_field_A(c,r)
              endif
           endif
        enddo
     enddo
  endif

  li = .false. 
  t = 1

  do r=1,SMAP_vodobs(source)%nr
     do c=1,SMAP_vodobs(source)%nc        
        vod1d(t) = vod_field(c,r)
        if(vod1d(t).ne.-9999.0) then 
           li(t) = .true.
        endif
        t = t+1
     enddo
  enddo
  
  call neighbor_interp(LVT_rc%gridDesc, li, vod1d, lo, vod_ip, &
       SMAP_vodobs(source)%nc*SMAP_vodobs(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMAP_vodobs(source)%rlat2, SMAP_vodobs(source)%rlon2,&
       SMAP_vodobs(source)%n112,LVT_rc%udef, iret)

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        vodobs(c,r) = vod_ip(c+(r-1)*LVT_rc%lnc)
     enddo
  enddo

  deallocate(vod_field)
  deallocate(dims)

#endif

end subroutine read_SMAPvod

!BOP
! 
! !ROUTINE: SMAP_vod_filename
! \label{SMAP_vod_filename}
!
! !INTERFACE: 
subroutine SMAP_vod_filename(source, name, designation, ndir, yr, mo,da)
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

end subroutine SMAP_vod_filename
