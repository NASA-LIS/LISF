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
! !ROUTINE: read_ANSASWEsnow
! \label{read_ANSASWEsnow}
!
! !REVISION HISTORY:
!  1 Jun 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_ANSASWEsnow(n, OBS_State, OBS_Pert_State) 
! !USES: 
#if(defined USE_HDF5) 
  use hdf5
#endif

  use ESMF
  use LIS_mpiMod
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_coreMod,  only : LIS_rc, LIS_domain, LIS_npes, LIS_localPet
  use LIS_timeMgrMod, only : LIS_calendar, LIS_clock
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, LIS_verify
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_pluginIndices, only : LIS_ANSASWEsnowobsId
  use ANSASWEsnow_Mod, only : ANSASWEsnow_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads ANSA snow retrievals
!  and packages it into an ESMF State with certain predefined 
!  attributes. The routine reads the snow data at 0z, performs 
!  spatial interpolation to the LIS grid and keeps it in memory. 
!  At 14:00 localtime for each grid point, the code then 
!  packages the interpolated observations into the ESMF state. 
!  Note that the data files are in HDF5 format and the routine
!  employs the HDF5-subsetting tools to extract and subset data. 
!
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
#if (defined USE_HDF5)
  character*100,   parameter    :: swe_field_name = "ansa_swe_cyl_GB"
  type(ESMF_Field)              :: snowField
  logical                       :: alarmCheck
  logical                       :: data_upd, file_exists
  logical                       :: dataflag(LIS_npes)
  logical                       :: dataflag_local
  integer                       :: c,r, p, t
  character(len=LIS_CONST_PATH_LEN) :: obsdir, ansa_filename
  integer                       :: file_id, swe_field_id
  integer(hsize_t), allocatable :: dims(:)
  integer(hid_t)                :: dataspace
  integer(hid_t)                :: memspace
  integer(hsize_t), dimension(2) :: dimsm 
  integer                       :: memrank = 2
  integer, allocatable, target  :: swe_field(:,:)
  real                          :: tsnow(ANSASWEsnow_struc(n)%nc*ANSASWEsnow_struc(n)%nr)
  logical*1                     :: li(ANSASWEsnow_struc(n)%nc*ANSASWEsnow_struc(n)%nr)
  real                          :: lon, lhour
  integer                       :: zone
  logical*1                     :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real,             pointer     :: obsl(:)
  integer                       :: gid(LIS_rc%ngrid(n))
  integer                       :: assimflag(LIS_rc%ngrid(n))
  integer(hsize_t), dimension(2) :: count_file 
  integer(hsize_t), dimension(2) :: count_mem 
  integer(hsize_t), dimension(2) :: offset_file 
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer                       :: status, iret, ierr


  dimsm      = (/ANSASWEsnow_struc(n)%nc, ANSASWEsnow_struc(n)%nr/)
  count_file = (/ANSASWEsnow_struc(n)%nc, ANSASWEsnow_struc(n)%nr/)
  count_mem  = (/ANSASWEsnow_struc(n)%nc, ANSASWEsnow_struc(n)%nr/)

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status,'Error in AttributeGet: Data Directory')

!-------------------------------------------------------------------------
!   Read the data at 0z daily. 
!-------------------------------------------------------------------------

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "ANSA SWE read alarm")

  if(alarmCheck) then 
     
     ANSASWEsnow_struc(n)%swe = LIS_rc%udef
     
     call ANSAsnow_filename2(ansa_filename,obsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)       

     inquire(file=ansa_filename,exist=file_exists)
     if(file_exists) then 

        write(LIS_logunit,*)  'Reading ANSA SWE data ',trim(ansa_filename)
        
        allocate(dims(2))
        dims(1) = ANSASWEsnow_struc(n)%nc
        dims(2) = ANSASWEsnow_struc(n)%nr

        offset_file = (/ANSASWEsnow_struc(n)%offset1, ANSASWEsnow_struc(n)%offset2/)
        
        allocate(swe_field(ANSASWEsnow_struc(n)%nc, ANSASWEsnow_struc(n)%nr))
        
        call h5open_f(status)
        call LIS_verify(status, 'Error opening HDF fortran interface')
        
        call h5fopen_f(trim(ansa_filename),H5F_ACC_RDONLY_F, file_id, status)
        call LIS_verify(status, 'Error opening ANSA file ')
        
        call h5dopen_f(file_id,swe_field_name,swe_field_id, status)
        call LIS_verify(status, 'Error opening SWE field in ANSA file')
        
        call h5dget_space_f(swe_field_id, dataspace, status)
        call LIS_verify(status, 'Error in h5dget_space_f: readANSASWE')
 
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
             start=offset_file, count=count_file, hdferr=status)
        call LIS_verify(status, 'Error setting hyperslab dataspace in readANSASWE')
        
        call h5screate_simple_f(memrank,dimsm, memspace, status)
        call LIS_verify(status, 'Error in h5create_simple_f; readANSASWE')
        
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
             start=offset_mem, count=count_mem, hdferr=status)
        call LIS_verify(status, 'Error in h5sselect_hyperslab_f: readANSASWE')

        call h5dread_f(swe_field_id, H5T_NATIVE_INTEGER,swe_field,dims,status, &
             memspace, dataspace)
        call LIS_verify(status, 'Error extracting SWE field from ANSA file')
        
        call h5dclose_f(swe_field_id,status)
        call LIS_verify(status,'Error in H5DCLOSE call')
        
        do r=1,ANSASWEsnow_struc(n)%nr
           do c=1,ANSASWEsnow_struc(n)%nc
              tsnow(c+(r-1)*ANSASWEsnow_struc(n)%nc) = swe_field(c,r)
           enddo
        enddo

        li  = .false.
        do c=1,ANSASWEsnow_struc(n)%mi
!           if(tsnow(c).lt.400) then 
!              li(c) = .true. 
!           endif
           if(tsnow(c).ne.494.and.tsnow(c).ne.496.and.&
                tsnow(c).ne.504.and.tsnow(c).ne.506.and.&
                tsnow(c).ne.508.and.tsnow(c).ne.510.and.&
                tsnow(c).ge.0) then 
              li(c) = .true. 
           endif
        enddo

        call bilinear_interp(LIS_rc%gridDesc(n,:),li,tsnow,&
             lo,ANSASWEsnow_struc(n)%swe,&
             ANSASWEsnow_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             ANSASWEsnow_struc(n)%w11,ANSASWEsnow_struc(n)%w12, &
             ANSASWEsnow_struc(n)%w21,ANSASWEsnow_struc(n)%w22, &
             ANSASWEsnow_struc(n)%n11,ANSASWEsnow_struc(n)%n12, &
             ANSASWEsnow_struc(n)%n21,ANSASWEsnow_struc(n)%n22, &
             LIS_rc%udef,iret)
        
        deallocate(swe_field)
           
        call h5fclose_f(file_id,status)
        call LIS_verify(status,'Error in H5FCLOSE call')
        
        call h5close_f(status)
        call LIS_verify(status,'Error in H5CLOSE call')

        deallocate(dims)

!-------------------------------------------------------------------------
!   Impose LSM-based QC on the observations (Flag data for rain, frozen
!   soil, and snow on the ground)
!-------------------------------------------------------------------------     
        
     endif
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

!throw out SWE < 10mm and > 200mm 
           if(ANSASWEsnow_struc(n)%swe(c+LIS_rc%lnc(n)*(r-1)).lt.10.or.&
                ANSASWEsnow_struc(n)%swe(c+LIS_rc%lnc(n)*(r-1)).gt.200) then 
              ANSASWEsnow_struc(n)%swe(c+LIS_rc%lnc(n)*(r-1)) = LIS_rc%udef
           endif
!localtime of this gridcell
           lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(c,r))%lon
           call LIS_localtime(LIS_rc%gmt,lon,lhour,zone)
           
           if(lhour.ge.(14-nint(LIS_rc%ts/3600.0)).and.lhour.le.14) then  
              obsl(LIS_domain(n)%gindex(c,r))=&
                   ANSASWEsnow_struc(n)%swe(c+LIS_rc%lnc(n)*(r-1))
!              if(c.eq.5.and.r.eq.9) print*, c,r,LIS_domain(n)%gindex(c,r),obsl(LIS_domain(n)%gindex(c,r))
           endif

        end if
     end do
  end do

  dataflag_local = .false. 

! LSM-based QC
  call lsmdaqcobsstate(LIS_rc%lsm, LIS_ANSASWEsnowobsId, &
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
     
     call ESMF_AttributeSet(snowfield,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status,'Error: AttributeSet in Grid Number')
     
     call ESMF_AttributeSet(snowfield,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')
           
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status, "Error: AttributeSet Data Update Status")
     return
  end if
  
#endif     
end subroutine read_ANSASWEsnow

!BOP
!
! !ROUTINE: ANSAsnow_filename2
! \label{ANSAsnow_filename2}
! 
! !INTERFACE: 
subroutine ANSAsnow_filename2(name, ndir, yr, mo,da)
  
  implicit none
! !ARGUMENTS: 
  character(len=*)  :: name, ndir
  integer           :: yr, mo, da
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped ANSA filename
!  
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the NESDIS AMSRE soil moisture filename
!  \item[ndir] name of the NESDIS AMSRE soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da 

  name = trim(ndir)//'/'//trim(fyr)//'/ansa_all_'//trim(fyr)//trim(fmo)//trim(fda)//'.h5'
    
end subroutine ANSAsnow_filename2



