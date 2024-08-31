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
! !ROUTINE: readSMAPTBobs
! \label{readSMAPTBobs}
!
! !INTERFACE: 
subroutine readSMAPTBobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit
  use SMAP_TBobsMod, only : SMAP_TBobs

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
! 
!EOP

  logical           :: alarmcheck, file_exists, readflag
  integer           :: iret
  character*200     :: ls_comm, cmd2, filename
  character*4       :: fyr
  character*2       :: fmo, fda
  real              :: smc(LVT_rc%lnc, LVT_rc%lnr,2)
  integer           :: fsize,k
  real              :: timenow

  smc = LVT_rc%udef

  timenow = float(LVT_rc%dhr(source))*3600 +&
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  if(SMAP_TBobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     SMAP_TBobs(source)%startflag = .false. 

 
     write(unit=fyr, fmt='(i4.4)') LVT_rc%dyr(source)
     write(unit=fmo, fmt='(i2.2)') LVT_rc%dmo(source)
     write(unit=fda, fmt='(i2.2)') LVT_rc%dda(source)
     
     ls_comm = 'ls '//trim(SMAP_TBobs(source)%odir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L1C_TB_*h5 2>&1 2>/dev/null > SMAP_file'

     cmd2 = 'wc -w SMAP_file > SMAP_file.wc'

     call system(ls_comm)
     call system(cmd2)
     
     open(110,file='SMAP_file.wc',form='formatted',action='read') 
     read(110,*) fsize
     close(110)

     open(110,file='SMAP_file',form='formatted',action='read')     
     do k=1,fsize
        read(110,'(a)') filename

        write(LVT_logunit,*) '[INFO] Reading SMAP file ',filename
        call read_SMAPTB(source, filename, smc)
     enddo     
     close(110)

!     open(100,file='test.bin',form='unformatted')
!     write(100) smc(:,:,1)
!     close(100)
!     stop
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source,&
       smc(:,:,1),vlevel=1,units="m3/m3")
 
end subroutine readSMAPTBobs

!BOP
! 
! !ROUTINE: read_SMAPTB
! \label{read_SMAPTB}
!
! !INTERFACE: 
subroutine read_SMAPTB(source, filename, TBobs)
! 
! !USES: 
  use LVT_coreMod,         only : LVT_rc,LVT_domain
  use LVT_logMod
  use SMAP_TBobsMod, only : SMAP_TBobs
  use map_utils
  use easeV2_utils
#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none

  integer                        :: source 
  character(len=*)               :: filename
  real                           :: TBobs(LVT_rc%lnc,LVT_rc%lnr,2)
  
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

  character*100,   parameter    :: sm_gr_name = "Global_Projection"
  character*100,   parameter    :: lat_field_name = "cell_lat"
  character*100,   parameter    :: lon_field_name = "cell_lon"
  character*100,   parameter    :: tbh_a_field_name = "cell_tb_h_aft"
  character*100,   parameter    :: tbv_a_field_name = "cell_tb_v_aft"
  character*100,   parameter    :: tbh_f_field_name = "cell_tb_h_fore"
  character*100,   parameter    :: tbv_f_field_name = "cell_tb_v_fore"
  integer(hid_t)                :: file_id, sm_gr_id
  integer(hid_t)                :: lat_field_id,lon_field_id
  integer(hid_t)                :: tbh_a_field_id,tbv_a_field_id
  integer(hid_t)                :: tbh_f_field_id,tbv_f_field_id
  integer(hsize_t), allocatable :: dims(:)
  integer(hid_t)                :: dataspace
  real, allocatable              :: lat_field(:)
  real, allocatable              :: lon_field(:)
  real, allocatable              :: tbh_a_field(:)
  real, allocatable              :: tbv_a_field(:)
  real, allocatable              :: tbh_f_field(:)
  real, allocatable              :: tbv_f_field(:)
  real                           :: udef
  integer                        :: t, c,r,icol,irow,gindex
  real                           :: col,row
  INTEGER(HSIZE_T) :: dimsr(1), maxdimsr(1)
  integer                        :: ngrid
  real                           :: lat_v, lon_v
  real                           :: tbh_in(SMAP_TBobs(source)%nc*&
       SMAP_TBobs(source)%nr)
  integer                        :: ntbh_in(SMAP_TBobs(source)%nc*&
       SMAP_TBobs(source)%nr)
  real                           :: tbv_in(SMAP_TBobs(source)%nc*&
       SMAP_TBobs(source)%nr)
  integer                        :: ntbv_in(SMAP_TBobs(source)%nc*&
       SMAP_TBobs(source)%nr)
  logical*1                      :: li(SMAP_TBobs(source)%nc*SMAP_TBobs(source)%nr)
  logical*1                      :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: tbh_ip(LVT_rc%lnc*LVT_rc%lnr)
  real                           :: tbv_ip(LVT_rc%lnc*LVT_rc%lnr)
  integer                        :: iret,status
 

  call h5open_f(status)
  call LVT_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(filename),H5F_ACC_RDONLY_F, file_id, status) 
  call LVT_verify(status, 'Error opening SMAP file ')
  
  call h5gopen_f(file_id,sm_gr_name,sm_gr_id, status)
  call LVT_verify(status, 'Error opening SM group in SMAP file')

  call h5dopen_f(sm_gr_id,lat_field_name,lat_field_id, status)
  call LVT_verify(status, 'Error opening lat field in SMAP file')

  call h5dget_space_f(lat_field_id, dataspace, status)
  call LVT_verify(status, 'Error in h5dget_space_f: readSMAPTB')


  CALL h5sget_simple_extent_dims_f(dataspace, dimsr, maxdimsr, status)
!  call LVT_verify(status, 'Error in h5sget_simple_extent_dims_f: readSMAPTB')

  ngrid = dimsr(1)

  allocate(lat_field(ngrid))
  allocate(lon_field(ngrid))
  allocate(tbh_a_field(ngrid))
  allocate(tbv_a_field(ngrid))
  allocate(tbh_f_field(ngrid))
  allocate(tbv_f_field(ngrid))


  call h5dopen_f(sm_gr_id,lon_field_name,lon_field_id, status)
  call LVT_verify(status, 'Error opening lon field in SMAP file')

  call h5dopen_f(sm_gr_id,tbh_a_field_name,tbh_a_field_id, status)
  call LVT_verify(status, 'Error opening TBH field in SMAP file')

  call h5dopen_f(sm_gr_id,tbh_f_field_name,tbh_f_field_id, status)
  call LVT_verify(status, 'Error opening TBH field in SMAP file')

  call h5dopen_f(sm_gr_id,tbv_a_field_name,tbv_a_field_id, status)
  call LVT_verify(status, 'Error opening TBV field in SMAP file')

  call h5dopen_f(sm_gr_id,tbv_f_field_name,tbv_f_field_id, status)
  call LVT_verify(status, 'Error opening TBV field in SMAP file')

  call h5dread_f(lat_field_id, H5T_NATIVE_REAL,lat_field,dimsr,status)
  call LVT_verify(status, 'Error extracting lat field from SMAPfile')

  call h5dread_f(lon_field_id, H5T_NATIVE_REAL,lon_field,dimsr,status)
  call LVT_verify(status, 'Error extracting lon field from SMAPfile')
  
  call h5dread_f(tbh_a_field_id, H5T_NATIVE_REAL,tbh_a_field,dimsr,status)
  call LVT_verify(status, 'Error extracting TBH field from SMAPfile')

  call h5dread_f(tbh_f_field_id, H5T_NATIVE_REAL,tbh_f_field,dimsr,status)
  call LVT_verify(status, 'Error extracting TBH field from SMAPfile')

  call h5dread_f(tbv_a_field_id, H5T_NATIVE_REAL,tbv_a_field,dimsr,status)
  call LVT_verify(status, 'Error extracting TBV field from SMAPfile')

  call h5dread_f(tbv_f_field_id, H5T_NATIVE_REAL,tbv_f_field,dimsr,status)
  call LVT_verify(status, 'Error extracting TBV field from SMAPfile')


  call h5dclose_f(lon_field_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(tbh_a_field_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(tbv_a_field_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(tbh_f_field_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(tbv_f_field_id,status)
  call LVT_verify(status,'Error in H5DCLOSE call')

  call h5gclose_f(sm_gr_id,status)
  call LVT_verify(status,'Error in H5GCLOSE call')

  call h5fclose_f(file_id,status)
  call LVT_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LVT_verify(status,'Error in H5CLOSE call')

  tbh_in = 0
  ntbh_in = 0 

  tbv_in = 0 
  ntbv_in = 0 

  do t=1,ngrid
     call easeV2_convert ('M36', lat_field(t), lon_field(t), col, row)
     
     icol = nint(col)
     irow = nint(row)
     
     if(icol.gt.0.and.irow.gt.0) then 
        gindex = icol + (irow -1)*SMAP_TBobs(source)%nc

        if(tbh_a_field(t).gt.0) then 
           tbh_in(gindex) = tbh_in(gindex) + tbh_a_field(t)
           ntbh_in(gindex) = ntbh_in(gindex) + 1
        endif
        
        if(tbh_f_field(t).gt.0) then 
           tbh_in(gindex) = tbh_in(gindex) + tbh_f_field(t)
           ntbh_in(gindex) = ntbh_in(gindex) + 1
        endif
        
        if(tbv_a_field(t).gt.0) then 
           tbv_in(gindex) = tbv_in(gindex) + tbv_a_field(t)
           ntbv_in(gindex) = ntbv_in(gindex) + 1
        endif
        
        if(tbv_f_field(t).gt.0) then 
           tbv_in(gindex) = tbv_in(gindex) + tbv_f_field(t)
           ntbv_in(gindex) = ntbv_in(gindex) + 1
        endif
     endif
  enddo
  
  li = .false. 
  do r=1,SMAP_TBobs(source)%nr
     do c=1,SMAP_TBobs(source)%nc
        gindex = c + (r -1)*SMAP_TBobs(source)%nc
        if(ntbh_in(gindex).gt.0) then 
           tbh_in(gindex) = tbh_in(gindex)/ntbh_in(gindex)
           li(gindex) = .true.
        else
           tbh_in(gindex) = LVT_rc%udef
        endif
     enddo
  enddo

  udef = -9999.0

  call neighbor_interp(LVT_rc%gridDesc, li, tbh_in, lo, tbh_ip, &
       SMAP_TBobs(source)%nc*SMAP_TBobs(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMAP_TBobs(source)%rlat2, SMAP_TBobs(source)%rlon2,&
       SMAP_TBobs(source)%n112,udef, iret)


  li = .false. 
  do r=1,SMAP_TBobs(source)%nr
     do c=1,SMAP_TBobs(source)%nc
        gindex = c + (r -1)*SMAP_TBobs(source)%nc
        if(ntbv_in(gindex).gt.0) then 
           tbv_in(gindex) = tbv_in(gindex)/ntbv_in(gindex)
           li(gindex) = .true.
        else
           tbv_in(gindex) = LVT_rc%udef
        endif
     enddo
  enddo

  call neighbor_interp(LVT_rc%gridDesc, li, tbv_in, lo, tbv_ip, &
       SMAP_TBobs(source)%nc*SMAP_TBobs(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMAP_TBobs(source)%rlat2, SMAP_TBobs(source)%rlon2,&
       SMAP_TBobs(source)%n112,udef, iret)


  deallocate(lat_field)
  deallocate(lon_field)
  deallocate(tbh_a_field)
  deallocate(tbv_a_field)
  deallocate(tbh_f_field)
  deallocate(tbv_f_field)

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(tbh_ip(c+(r-1)*LVT_rc%lnc).gt.0) then 
           TBobs(c,r,1) = tbh_ip(c+(r-1)*LVT_rc%lnc)
        endif
        if(tbv_ip(c+(r-1)*LVT_rc%lnc).gt.0) then 
           TBobs(c,r,2) = tbv_ip(c+(r-1)*LVT_rc%lnc)
        endif
     enddo
  enddo
#endif

end subroutine read_SMAPTB

subroutine mapLatLonToSMAPindex(lat_in,lon_in,nsize, lats,lons,d_index)
  
  implicit none

  real         :: lat_in
  real         :: lon_in
  integer      :: nsize
  real         :: lats(nsize)
  real         :: lons(nsize)
  integer      :: d_index


  integer      :: k
  real         :: minval, temp

  minval  = 10000.0
  d_index = -1

  do k=1,nsize
     temp = (lats(k)-lat_in)**2 + (lons(k)-lon_in)**2
     
     if(temp.lt.minval) then 
        minval = temp
        d_index = k
     endif
  enddo
  if(temp.gt.10) d_index  = -1
  if(d_index.gt.0) then 
     print*, lat_in, lon_in, lats(d_index), lons(d_index), temp
  endif

end subroutine mapLatLonToSMAPindex
