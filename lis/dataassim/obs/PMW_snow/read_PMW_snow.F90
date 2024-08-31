!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_PMW_snow
! \label{read_PMW_snow}
!
! !REVISION HISTORY:
!  1 Jun 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_PMW_snow(n, OBS_State,OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_pluginIndices, only : LIS_PMWsnowobsId
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use PMW_snow_Mod, only : PMW_snow_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads PMW-based snow observations
!  and packages it into an ESMF State with certain predefined 
!  attributes. The routine reads the snow data at 0z, performs 
!  spatial interpolation to the LIS grid and keeps it in memory. 
!  At given localtime for each grid point, the code then 
!  packages the interpolated observations into the ESMF state. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)              :: snowField,pertfield
  logical                       :: alarmCheck
  logical                       :: data_upd, file_exists
  logical                       :: dataflag(LIS_npes)
  logical                       :: dataflag_local
  integer                       :: c,r, p, t
  character(len=LIS_CONST_PATH_LEN):: obsdir, pmw_filename
  real                          :: lon, lhour, lhour1
  integer                       :: zone
  real                          :: ssdev(LIS_rc%ngrid(n))
  logical*1                     :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real,             pointer     :: obsl(:)
  integer                       :: gid(LIS_rc%ngrid(n))
  integer                       :: assimflag(LIS_rc%ngrid(n))
  integer                       :: status, iret, ierr

  lhour1 = PMW_snow_struc(n)%assim_lhour

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status,'Error in AttributeGet: Data Directory')

!-------------------------------------------------------------------------
!   Read the data at 0z daily. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "PMW snow data read alarm")

  if(alarmCheck.or.PMW_snow_struc(n)%startMode) then 
     PMW_snow_struc(n)%startMode = .false.
     
     PMW_snow_struc(n)%snow = LIS_rc%udef
     
     call PMW_snow_filename(pmw_filename,PMW_snow_struc(n)%data_fn_conv, obsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)       

     inquire(file=pmw_filename,exist=file_exists)
     if(file_exists) then 
        write(LIS_logunit,*)  'Reading PMW snow data from ',trim(pmw_filename)

        if (PMW_snow_struc(n)%data_format .eq. 'HDF4' .and. PMW_snow_struc(n)%data_coordsys .eq. 'EASE') then          
           call read_PMWSnow_HDF4(n, pmw_filename)
  
        elseif (PMW_snow_struc(n)%data_format .eq. 'HDF-EOS' .and. PMW_snow_struc(n)%data_coordsys .eq. 'EASE') then
           call read_PMWSnow_HDFEOS(n, pmw_filename)

        elseif (PMW_snow_struc(n)%data_format .eq. 'HDF5' .and. PMW_snow_struc(n)%data_coordsys .eq. 'LATLON') then
           call read_PMWSnow_HDF5(n, pmw_filename)
        endif 
     endif

        !if (LIS_rc%mo .eq. 11) then
        !open(100,file='snowdata.bin',form='unformatted')
        !write(100) PMW_snow_struc(n)%snow
        !close(100)
        !print*, PMW_snow_struc(n)%snow
        !stop
        !endif

!        call neighbor_interp(LIS_rc%gridDesc(n,:), ibi, li, tsnow_flag,&
!             ibo, lo,ANSAsnow_struc(n)%snwd_flag, &
!             ANSAsnow_struc(n)%mi, LIS_rc%lnc(n)*LIS_rc%lnr(n),&
!             ANSAsnow_struc(n)%rlat,ANSAsnow_struc(n)%rlon,&
!             ANSAsnow_struc(n)%n113,LIS_rc%udef,iret)

!        deallocate(snwd_flag_field)


  endif

!-------------------------------------------------------------------------
!  Update the OBS_State 
!-------------------------------------------------------------------------     

  call ESMF_StateGet(OBS_State,"Observation01",snowfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(snowfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
  obsl = LIS_rc%udef 
  
  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
     
!localtime of this gridcell
           lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(c,r))%lon
           call LIS_localtime(LIS_rc%gmt,lon,lhour,zone)

           !if(lhour.gt.(lhour1-nint(LIS_rc%ts/3600.0)).and.lhour.le.lhour1) then  
           if(lhour.le.lhour1 .and. (lhour1-lhour).lt.LIS_rc%ts/3600.0) then  
              obsl(LIS_domain(n)%gindex(c,r))=&
                   PMW_snow_struc(n)%snow(c+LIS_rc%lnc(n)*(r-1))
           endif

        end if
     end do
  end do

  dataflag_local = .false. 

! LSM-based QC
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_PMWsnowobsId)//char(0), & 
       n, OBS_state)

  call ESMF_StateGet(OBS_State,"Observation01",snowField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(snowField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)

  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           if(obsl(LIS_domain(n)%gindex(c,r)).ne.-9999.0) then 
              dataflag_local = .true. 
           endif
        endif
     end do
  end do


#if (defined SPMD)
  call MPI_ALLGATHER(dataflag_local,1, MPI_LOGICAL, dataflag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, ierr)
#endif
  data_upd = .false.
  
  do p=1,LIS_npes
     data_upd = data_upd.or.dataflag(p)
  enddo

  if(data_upd) then 
     do t=1,LIS_rc%ngrid(n)
        gid(t) = t
        if(obsl(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo
     
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status, 'Error: AttributeSet in Data Update Status')
     

     call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
          rc=status)
     call LIS_verify(status, 'ESMF_StateGet for Observation01 for OBS_Pert_State failed in read_PMW_snow')
     
     if(LIS_rc%ngrid(n).gt.0) then 

!linearly scale the observation err
        ssdev = PMW_snow_struc(n)%ssdev 
        do t=1,LIS_rc%ngrid(n)
           if(obsl(t).ne.-9999.0) then 
              ssdev(t) =  PMW_snow_struc(n)%ssdev !+ 0.05*obsl(t)
!for values adjusted with confidence
!              if(PMW_snow_struc(n)%snwd_flag(t).eq.1) then 
!                 ssdev(t) = 0.1 !assuming multiplicative
!              endif
           endif
        enddo

        call ESMF_AttributeSet(pertField,"Standard Deviation",&
             ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(snowfield,"Grid Number",&
             gid,itemCount=LIS_rc%ngrid(n),rc=status)
        call LIS_verify(status,'Error: AttributeSet in Grid Number')
        
        call ESMF_AttributeSet(snowfield,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
        call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')
        
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status, "Error: AttributeSet Data Update Status")
     return
  end if
  
end subroutine read_PMW_snow


!BOP
! !ROUTINE: read_PMWSnow_HDF4
! \label{read_PMWSnow_HDF4}
! 
! !INTERFACE: 
subroutine read_PMWSnow_HDF4(n,name)
! !USES: 
  use LIS_coreMod,      only : LIS_rc,LIS_domain
  use LIS_logmod,       only : LIS_logunit
  use PMW_snow_Mod, only : PMW_snow_struc
  implicit none

#if (defined USE_HDF4) 
#include "hdf.f90"
#endif
! !ARGUMENTS:   
  integer, intent(in)  :: n 
  character(len=*)     :: name
! 
! !DESCRIPTION: 
!   This routine extracts the SWE or snow depth retrievals and the associated 
!   quality control flags from the HDF4 files. 
! 
!EOP

  real                 :: sb_qc(PMW_snow_struc(n)%mo)
  real                 :: snowobs(PMW_snow_struc(n)%mo)
  
#if (defined USE_HDF4)
  !declare the hdf4 library functions
  integer              :: sfstart, sfselect, sfrdata, sfendacc, sfend
  integer              :: snow_sdsid(2),qc_sdsid(2)
  integer              :: size,igd
  integer              :: file_id,ret
  real,allocatable         :: snow(:)
  real,allocatable      :: qc(:) ! EMK Corrected variable type
  integer,parameter    :: ease_nr=721, ease_nc=721
  integer              :: mi,iret
  integer              :: start(2),edge(2),stride(2)
  logical*1,allocatable      :: li(:)! EMK Corrected variable type
  logical*1              :: lo(PMW_snow_struc(n)%mo) ! EMK Corrected type
  integer              :: i, t, sd_index, sds_id

  sb_qc = LIS_rc%udef
  snowobs = LIS_rc%udef

  ! SDS ids
  snow_sdsid(1)  = PMW_snow_struc(n)%data_nl_sdsid 
  snow_sdsid(2)  = PMW_snow_struc(n)%data_sl_sdsid 
  qc_sdsid(1)   = PMW_snow_struc(n)%flag_nl_sdsid 
  qc_sdsid(2)   = PMW_snow_struc(n)%flag_sl_sdsid 

  !open the hdf4 file for read access

  file_id = sfstart(trim(name),DFACC_READ)
  if (file_id.eq.-1)then
     write(LIS_logunit,*)"Failed to open hdf4 file",trim(name)
     return
  end if

  mi = ease_nr*ease_nc
  allocate(li(mi))

  do igd=PMW_snow_struc(n)%ihemi, PMW_snow_struc(n)%nhemi

     !get snow data SDS
     sd_index = 0
     sds_id = sfselect(snow_sdsid(igd), sd_index)
     if (sds_id .eq. -1) then
        write(LIS_logunit, *) "failed to get snow sds id", snow_sdsid(igd)
        return
     end if

     !retrieve the data
     start(1)=0  !hdf4 lib uses 0-based count
     start(2)=0
     edge(1)=ease_nc
     edge(2)=ease_nr
     stride(1)=1
     stride(2)=1

     allocate(snow(ease_nc*ease_nr))
     ret = sfrdata(sds_id,start,stride,edge,snow)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to get the snow field"
        deallocate(snow)
        deallocate(li)
        return
     end if

     ! terminate access to the snow sds
     ret = sfendacc(sds_id)

     li=.false.
     do t=1,mi
        if(snow(t).ge.0) then
           li(t)=.true.
        endif
     enddo
     call neighbor_interp(LIS_rc%gridDesc(n,:),li,snow,&
          lo,snowobs,mi,PMW_snow_struc(n)%mo, &
          PMW_snow_struc(n)%rlat2(:,igd),PMW_snow_struc(n)%rlon2(:,igd),&
          PMW_snow_struc(n)%n112(:,igd),LIS_rc%udef,iret)
     deallocate(snow)

     if (PMW_snow_struc(n)%use_flag .eq. 1) then
        !get qc sds id
        sd_index = 0
        sds_id = sfselect(qc_sdsid(igd), sd_index)
        if (sds_id .eq. -1) then
           write(LIS_logunit, *) "failed to select sds id", qc_sdsid(igd)
           return
        end if
   
        !retrieve the qc data
        allocate(qc(ease_nc*ease_nr))
        ret = sfrdata(sds_id,start,stride,edge,qc)
        if (ret <0)then
           write(LIS_logunit,*)"Failed to get the snow flag field"
           deallocate(qc)
           deallocate(li)
           return
        end if
   
        ! terminate access to the flag SDS
        ret = sfendacc(sds_id)
   
        li = .true.
        call neighbor_interp(LIS_rc%gridDesc(n,:),li,qc,&
             lo,sb_qc,mi,PMW_snow_struc(n)%mo, &
             PMW_snow_struc(n)%rlat2(:,igd),PMW_snow_struc(n)%rlon2(:,igd),&
             PMW_snow_struc(n)%n112(:,igd),LIS_rc%udef,iret)
   
        deallocate(qc)

     end if 
  end do
  deallocate(li)
  ret=sfend(file_id)
  if (ret <0)then
     write(LIS_logunit,*)"Failed to close file: ",file_id
  end if

  ! qc PMW snow data and save to PMW_snow_struc(n)%snow
  call qc_PMWsnow(n,snowobs,sb_qc, PMW_snow_struc(n)%mo, PMW_snow_struc(n)%snow)

#endif

end subroutine read_PMWSnow_HDF4



!BOP
! !ROUTINE: read_PMWSnow_HDFEOS
! \label{read_PMWSnow_HDFEOS}
! 
! !INTERFACE: 
subroutine read_PMWSnow_HDFEOS(n,name)
! !USES: 
  use LIS_coreMod,      only : LIS_rc,LIS_domain
  use LIS_logmod,       only : LIS_logunit
  use PMW_snow_Mod, only : PMW_snow_struc
  implicit none

#if (defined USE_HDFEOS2)
#include "hdf.f90"
#endif
! !ARGUMENTS:   
  integer, intent(in)  :: n 
  character(len=*)     :: name
! 
! !DESCRIPTION: 
!   This routine extracts the SWE or snow depth retrievals and the associated 
!   quality control flags from the HDF-EOS files. 
! 
!EOP
  real                 :: sb_qc(PMW_snow_struc(n)%mo)
  real                 :: snowobs(PMW_snow_struc(n)%mo)

#if (defined USE_HDFEOS2)
  !declare the hdf-eos library functions 
  integer              :: gdopen,gdattach,gdrdfld
  integer              :: gddetach,gdclose
  character*50         :: grid_name(2),snow_name(2),qc_name(2)
  integer              :: size,igd
  integer              :: file_id,grid_id,ret
  integer*1,allocatable    :: snow(:), qc(:)
  real,allocatable         :: rqc(:)
  integer,parameter    :: ease_nr=721, ease_nc=721
  real,allocatable     :: rsnow(:)
  integer              :: mi,iret
  integer              :: start(2),edge(2),stride(2)
  logical*1,allocatable    :: li(:)
  logical*1            :: lo(PMW_snow_struc(n)%mo)
  integer              :: i, t

  sb_qc = LIS_rc%udef
  snowobs = LIS_rc%udef

  !Grid and field names
  grid_name(1) = PMW_snow_struc(n)%data_nl_grid 
  grid_name(2) = PMW_snow_struc(n)%data_sl_grid 
  snow_name(1)  = PMW_snow_struc(n)%data_nl_sds 
  snow_name(2)  = PMW_snow_struc(n)%data_sl_sds 
  qc_name(1)   = PMW_snow_struc(n)%flag_nl_sds 
  qc_name(2)   = PMW_snow_struc(n)%flag_sl_sds 

  !open the hdf-eos file
  file_id = gdopen(trim(name),DFACC_READ)
  if (file_id.eq.-1)then
     write(LIS_logunit,*)"Failed to open hdf file",trim(name)
     return
  end if
  
  mi = ease_nr*ease_nc
  allocate(li(mi))
  do igd=PMW_snow_struc(n)%ihemi, PMW_snow_struc(n)%nhemi     
    !get the grid id
     grid_id = gdattach(file_id,grid_name(igd))
     if (grid_id.eq.-1)then
        write(LIS_logunit,*)"Failed to attach grid: ",grid_name(igd),trim(name)
        ret = gdclose(file_id)
        deallocate(li)
        return
     end if
     
     start(1)=0  !hdfeos lib uses 0-based count
     start(2)=0
     edge(1)=ease_nc
     edge(2)=ease_nr
     stride(1)=1
     stride(2)=1

     allocate(snow(ease_nc*ease_nr))
     ret = gdrdfld(grid_id,snow_name(igd),start,stride,edge,snow)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to get the snow field"
        ret=gddetach(grid_id)
        ret=gdclose(file_id)
        deallocate(snow)
        deallocate(li)
        return
     end if
        
        !convert short int snow to real 
     allocate(rsnow(ease_nc*ease_nr))
     do i=1,ease_nc*ease_nr
        rsnow(i)=snow(i)
     end do
     li=.false.
     do t=1,mi
        if(rsnow(t).ge.0) then 
           li(t)=.true.
        endif
     enddo
     call neighbor_interp(LIS_rc%gridDesc(n,:),li,rsnow,&
          lo,snowobs,mi,PMW_snow_struc(n)%mo, &
          PMW_snow_struc(n)%rlat2(:,igd),PMW_snow_struc(n)%rlon2(:,igd),&
          PMW_snow_struc(n)%n112(:,igd),LIS_rc%udef,iret)
     
     deallocate(snow)
     deallocate(rsnow)

     if (PMW_snow_struc(n)%use_flag .eq. 1) then
        !get qc
        allocate(qc(ease_nc*ease_nr))
        ret =  gdrdfld(grid_id,qc_name(igd),start,stride,edge,qc)
        !convert to the real number to use the neighbor_interp call
        allocate(rqc(ease_nc*ease_nr))
        do i=1,ease_nc*ease_nr
           rqc(i)=qc(i)
        end do
        deallocate(qc)
        li=.false.
        do t=1,mi
           if(rqc(t).gt.0) then 
              li(t)=.true.
           endif
        enddo
   
        call neighbor_interp(LIS_rc%gridDesc(n,:),li,rqc,&
             lo,sb_qc,mi,PMW_snow_struc(n)%mo, &
             PMW_snow_struc(n)%rlat2(:,igd),PMW_snow_struc(n)%rlon2(:,igd),&
             PMW_snow_struc(n)%n112(:,igd),LIS_rc%udef,iret)
   
        deallocate(qc)
        deallocate(rqc) 
     endif
     ret=gddetach(grid_id)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to detach grid_id: ",grid_id
     end if

  end do
  deallocate(li)
  ret=gdclose(file_id)
  if (ret <0)then
     write(LIS_logunit,*)"Failed to close file: ",file_id
  end if

  ! qc PMW snow data and save to PMW_snow_struc(n)%snow
  call qc_PMWsnow(n,snowobs, sb_qc, PMW_snow_struc(n)%mo, PMW_snow_struc(n)%snow)


#endif
  
end subroutine read_PMWSnow_HDFEOS



!BOP
! !ROUTINE: read_PMWSnow_HDF5
! \label{read_PMWSnow_HDF5}
! 
! !INTERFACE: 
subroutine read_PMWSnow_HDF5(n,name)
! !USES: 
#if (defined USE_HDF5) 
  use hdf5
#endif
  use LIS_coreMod,      only : LIS_rc,LIS_domain
  use LIS_logMod,       only : LIS_logunit, LIS_verify
  use PMW_snow_Mod,     only : PMW_snow_struc
  implicit none

! !ARGUMENTS:   
  integer, intent(in)  :: n 
  character(len=*)     :: name
! 
! !DESCRIPTION: 
!   This routine extracts the SWE or snow depth retrievals and the associated 
!   quality control flags from the LATLON HDF5 files. 
! dim
!EOP
#if (defined USE_HDF5)
   integer(hsize_t), allocatable  :: dims(:)
   integer(hid_t)                 :: file_id, snow_field_id,snow_flag_field_id
   integer, allocatable, target   :: snow_field(:,:),snow_flag_field(:,:)
   real                           :: tsnow(PMW_snow_struc(n)%mi)
   real                           :: tsnow_flag(PMW_snow_struc(n)%mi)
   real                           :: tsnow_qc(PMW_snow_struc(n)%mi)
   integer(hid_t)                 :: dataspace, memspace
   integer                        :: memrank = 2
   integer(hsize_t), dimension(2) :: count_file, count_mem, offset_file
   integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
   integer                        :: status  
   logical*1                      :: li(PMW_snow_struc(n)%mi)
   logical*1                      :: lo(PMW_snow_struc(n)%mo)
   integer                        :: r, c, iret,i

   !dims      = (/PMW_snow_struc(n)%nc, PMW_snow_struc(n)%nr/)
   count_file = (/PMW_snow_struc(n)%nc, PMW_snow_struc(n)%nr/)
   count_mem  = (/PMW_snow_struc(n)%nc, PMW_snow_struc(n)%nr/)
        
   allocate(dims(2))
   dims(1) = PMW_snow_struc(n)%nc
   dims(2) = PMW_snow_struc(n)%nr
   offset_file = (/PMW_snow_struc(n)%offset1, &
        PMW_snow_struc(n)%offset2/)
   
   allocate(snow_field(PMW_snow_struc(n)%nc, &
        PMW_snow_struc(n)%nr))
   allocate(snow_flag_field(PMW_snow_struc(n)%nc, &
        PMW_snow_struc(n)%nr))
   
   !open the fortran interface
   call h5open_f(status)
   call LIS_verify(status, 'Error opening HDF fortran interface')
   
   !open the file
   call h5fopen_f(trim(name),H5F_ACC_RDONLY_F, file_id, status)
   call LIS_verify(status, 'Error opening PMW snow file ')

   !open the snow dataset        
   call h5dopen_f(file_id,PMW_snow_struc(n)%snow_field_name,snow_field_id, status)
   call LIS_verify(status, 'Error opening snow field in PMW snow data file')

   !get dataspace identifier
   call h5dget_space_f(snow_field_id, dataspace, status)
   call LIS_verify(status, 'Error in h5dget_space_f in read_PMWSnow_HDF5')

   !select hyperslab in dataset
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
        start=offset_file, count=count_file, hdferr=status)
   call LIS_verify(status, 'Error setting hyperslab dataspace in read_PMWSnow_HDF5')

   !create memory dataspace
   call h5screate_simple_f(memrank,dims, memspace, status)
   call LIS_verify(status, 'Error in h5create_simple_f in read_PMWSnow_HDF5')

   !select hyperslab in memeory
   call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
        start=offset_mem, count=count_mem, hdferr=status)
   call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_PMWSnow_HDF5')

   !read data from hyperslab in file into hyperslab in memeory
   call h5dread_f(snow_field_id, H5T_NATIVE_INTEGER,snow_field,&
        dims,status, memspace, dataspace)
   call LIS_verify(status, 'Error extracting snow field from PMW snow data file')

   !close dataspace for the dataset
   call h5sclose_f(dataspace, status)
   call LIS_verify(status,'Error in closing dataspace for the snow dataset')

   !close memory space 
   call h5sclose_f(memspace, status)
   call LIS_verify(status,'Error in closing memeroy space')
   
   !close the snow dataset
   call h5dclose_f(snow_field_id,status)
   call LIS_verify(status,'Error in closing the snow dataset')
   
   if(PMW_snow_struc(n)%use_flag .eq. 1) then

      !open snow_flag dataset
      call h5dopen_f(file_id,PMW_snow_struc(n)%snow_flag_field_name,snow_flag_field_id, status)
      call LIS_verify(status, 'Error opening snow flag field in PMW snow data file')
   
      !get dataspace identifier
      call h5dget_space_f(snow_flag_field_id, dataspace, status)
      call LIS_verify(status, 'Error in h5dget_space_f: read_PMWSnow_HDF5')

      !select hyperslab in dataset
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
      call LIS_verify(status, 'Error setting hyperslab dataspace in read_PMWSnow_HDF5')

      !create memory dataspace
      call h5screate_simple_f(memrank,dims, memspace, status)
      call LIS_verify(status, 'Error in h5create_simple_f; read_PMWSnow_HDF5')

      !select hyperslab in memeory
      call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
      call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_PMWSnow_HDF5')

      !read data from hyperslab in file into hyperslab in memeory
      call h5dread_f(snow_flag_field_id, H5T_NATIVE_INTEGER,snow_flag_field,&
         dims,status, memspace, dataspace)
      call LIS_verify(status, 'Error extracting snow flag field from PMW snow data file')

      !close dataspace for the dataset
      call h5sclose_f(dataspace, status)
      call LIS_verify(status,'Error in closing dataspace for the snow flag dataset')

      !close memory space 
      call h5sclose_f(memspace, status)
      call LIS_verify(status,'Error in closing memeroy space') 

      !close snow flag dataset  
      call h5dclose_f(snow_flag_field_id,status)
      call LIS_verify(status,'Error in closing the snow flag dataset')
   
   endif   
   !close file
   call h5fclose_f(file_id, status)
   call LIS_verify(status,'Error in closing the PMW snow data file')

   !close fortran interface
   call h5close_f(status)
   call LIS_verify(status,'Error in closing the PMW snow fortran interface')


   do r=1,PMW_snow_struc(n)%nr
      do c=1,PMW_snow_struc(n)%nc
         tsnow(c+(r-1)*PMW_snow_struc(n)%nc) = real(snow_field(c,r))
         if(PMW_snow_struc(n)%use_flag .eq. 1) then
            tsnow_flag(c+(r-1)*PMW_snow_struc(n)%nc) = real(snow_flag_field(c,r))
         endif
      enddo
   enddo   

   !apply QC before bilinear interpolation  
   call qc_PMWsnow(n, tsnow, tsnow_flag, PMW_snow_struc(n)%mi, tsnow_qc)
       
   li = .false.
   do i=1,PMW_snow_struc(n)%mi
      if (tsnow_qc(i) .ne. LIS_rc%udef) li(i) = .true.
   enddo  
   call bilinear_interp(LIS_rc%gridDesc(n,:),li,tsnow_qc,&
        lo,PMW_snow_struc(n)%snow,&
        PMW_snow_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),&
        LIS_domain(n)%lat, LIS_domain(n)%lon,&
        PMW_snow_struc(n)%w11,PMW_snow_struc(n)%w12, &
        PMW_snow_struc(n)%w21,PMW_snow_struc(n)%w22, &
        PMW_snow_struc(n)%n11,PMW_snow_struc(n)%n12, &
        PMW_snow_struc(n)%n21,PMW_snow_struc(n)%n22, &
        LIS_rc%udef,iret)

   deallocate(dims)
   deallocate(snow_field)
   deallocate(snow_flag_field)

#endif
  
end subroutine read_PMWSnow_HDF5
  
        
!BOP
!
! !ROUTINE: PMW_snow_filename
! \label{PMW_snow_filename}
! 
! !INTERFACE: 
subroutine PMW_snow_filename(name, data_fn_conv, ndir, yr, mo,da)
  
  use LIS_coreMod, only: LIS_rc
  use LIS_logMod, only: LIS_logunit

  implicit none
! !ARGUMENTS: 
  character(len=*)     :: name
  integer           :: yr, mo, da
  character (len=*) :: ndir
  character (len=*) :: data_fn_conv
  character (len=50) :: str1, str2, str0
  integer           :: len1, i1, i2
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped PMW\_snow data filename
!  
!  The arguments are: 
!  \begin{description}
!  \item[name] PMW snow data file name
!  \item[data\_fn\_conv] PMW snow data file name convention
!  \item[ndir] PMW snow data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  character (len=3) :: fdoy
  
  len1 = len_trim(data_fn_conv)
  i1 = index(data_fn_conv, 'YYYYMMDD')
  i2 = index(data_fn_conv, 'YYYYDOY')

  if (i1 .gt. 0) then
    str1 = data_fn_conv(1:i1-1)
    str2 = data_fn_conv(i1+8:len1)
    write(unit=fyr, fmt='(i4.4)') yr
    write(unit=fmo, fmt='(i2.2)') mo
    write(unit=fda, fmt='(i2.2)') da 
    str0 = fyr//fmo//fda
  else if (i2 .gt. 0) then
    str1 = data_fn_conv(1:i2-1)
    str2 = data_fn_conv(i2+7:len1)
    write(unit=fyr, fmt='(i4.4)') yr
    write(unit=fdoy, fmt='(i3.3)') LIS_rc%doy
    str0 = fyr//fdoy
  else
    write(LIS_logunit,*) 'PMW snow data file name convention currently not supported: ', trim(data_fn_conv)
  end if

  name = trim(ndir)//'/'//trim(fyr)//'/'//trim(str1)//trim(str0)//trim(str2)
    
end subroutine PMW_snow_filename


!BOP
! !ROUTINE: qc_PMWSnow
! \label{qc_PMWSnow}
!
! !INTERFACE:
subroutine qc_PMWSnow(n,snowobs, qc_flag, npts, snowobs_qc)
! !USES:
  use LIS_coreMod,      only : LIS_rc
  use PMW_snow_Mod, only : PMW_snow_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in)  :: n
  integer, intent(in)  :: npts
  real, intent(in), dimension(npts)         :: snowobs, qc_flag
  real, intent(inout), dimension(npts)      :: snowobs_qc
  integer              :: i,j

!
! !DESCRIPTION:
!   This routine applies various qc conditions to the PMW snow data
!
!EOP

  do j=1, npts
     snowobs_qc(j) = snowobs(j)  
     if (snowobs_qc(j) .lt. 0) snowobs_qc(j) = LIS_rc%udef
  enddo

! first apply the flags retrieved from PMW data file
  if (PMW_snow_struc(n)%use_flag .eq. 1) then
     if (PMW_snow_struc(n)%flag_n1_invalid_value .gt. 0) then
        do j=1, npts 
           do i=1,PMW_snow_struc(n)%flag_n1_invalid_value                        
              if (abs(qc_flag(j)-PMW_snow_struc(n)%flag_invalid_value(i)) .lt. 1.0e-20) then
                 snowobs_qc(j) = LIS_rc%udef
              endif
           enddo
        enddo 
     endif
  endif

! next mark those PMW data values considered as invalid
  if (PMW_snow_struc(n)%data_n1_invalid_value .gt. 0) then
     do j=1, npts 
        do i=1,PMW_snow_struc(n)%data_n1_invalid_value             
           if (abs(snowobs_qc(j)-PMW_snow_struc(n)%data_invalid_value(i)) .lt. 1.0e-12) then
              snowobs_qc(j) = LIS_rc%udef
           endif
        enddo
     enddo 
  endif

! apply PMW snow data scale factor
  do j=1, npts
     if (snowobs_qc(j) .ne. LIS_rc%udef) then
        snowobs_qc(j) = snowobs_qc(j)*PMW_snow_struc(n)%data_scale
     endif
  enddo

! additional qc with min & max values
  if (PMW_snow_struc(n)%use_minmax .eq. 1) then
     do j=1, npts
        if (snowobs_qc(j).ne.LIS_rc%udef) then
           if (snowobs_qc(j).gt.PMW_snow_struc(n)%data_max .or. &
               snowobs_qc(j).lt.PMW_snow_struc(n)%data_min) then 
               snowobs_qc(j) = LIS_rc%udef
           endif
        endif
     enddo
  endif
        
end subroutine qc_PMWsnow
