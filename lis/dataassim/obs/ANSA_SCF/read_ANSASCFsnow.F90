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
! !ROUTINE: read_ANSASCFsnow
! \label{read_ANSASCFsnow}
!
! !REVISION HISTORY:
!  1 Jun 2009: Sujay Kumar; Initial Specification
!  Aug 2013 Yuqiong Liu; adpated for ANSA SCF assimilation using EnKF or DI

! !INTERFACE: 
subroutine read_ANSASCFsnow(n, OBS_State,OBS_Pert_State) 
! !USES: 
#if(defined USE_HDF5) 
  use hdf5
#endif

  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!  use LIS_pluginIndices, only : LIS_ANSASCFsnowobsId
  use ANSASCFsnow_Mod, only : ANSASCFsnow_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads ANSA snow observations
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
  !character*100,   parameter    :: sca_field_name = "ansa_interpsnow_cyl_GB"

  type(ESMF_Field)              :: snowField,pertfield
  logical                       :: alarmCheck
  logical                       :: data_upd, file_exists
  logical                       :: dataflag(LIS_npes)
  logical                       :: dataflag_local
  integer                       :: c,r, p, t
  character(len=LIS_CONST_PATH_LEN) :: obsdir, ansa_filename
  integer(hid_t)                :: file_id, sca_field_id
  integer(hsize_t), allocatable :: dims(:)
  integer(hid_t)                :: dataspace
  integer(hid_t)                :: memspace
  integer(hsize_t), dimension(2) :: dimsm 
  integer                       :: memrank = 2
  integer, allocatable, target  :: sca_field(:,:)
  real                          :: tsnow(ANSASCFsnow_struc(n)%nc*&
       ANSASCFsnow_struc(n)%nr)
  logical*1                     :: li(ANSASCFsnow_struc(n)%nc*ANSASCFsnow_struc(n)%nr)
  real                          :: lon, lhour, lhour1
  integer                       :: zone
  real                          :: ssdev(LIS_rc%ngrid(n))
  logical*1                     :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real,             pointer     :: obsl(:)
  integer                       :: gid(LIS_rc%ngrid(n))
  integer                       :: assimflag(LIS_rc%ngrid(n))
  integer(hsize_t), dimension(2) :: count_file 
  integer(hsize_t), dimension(2) :: count_mem 
  integer(hsize_t), dimension(2) :: offset_file 
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer                       :: status, iret, ierr
  real                          :: tmp


!  if(LIS_rc%yr.ge.2010) then 
!     ANSASCFsnow_struc(n)%minlat = 0.025
!     ANSASCFsnow_struc(n)%offset1 = &
!          nint((ANSASCFsnow_struc(n)%cornerlon1-&
!          ANSASCFsnow_struc(n)%minlon)/0.05)
!     ANSASCFsnow_struc(n)%offset2 = &
!          nint((ANSASCFsnow_struc(n)%cornerlat1-&
!          ANSASCFsnow_struc(n)%minlat)/0.05)
     
!  endif

  lhour1 = ANSASCFsnow_struc(n)%assim_lhour

  dimsm      = (/ANSASCFsnow_struc(n)%nc, ANSASCFsnow_struc(n)%nr/)
  count_file = (/ANSASCFsnow_struc(n)%nc, ANSASCFsnow_struc(n)%nr/)
  count_mem  = (/ANSASCFsnow_struc(n)%nc, ANSASCFsnow_struc(n)%nr/)

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status,'Error in AttributeGet: Data Directory')

!-------------------------------------------------------------------------
!   Read the data at 0z daily. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "ANSA SCF read alarm")

  if(alarmCheck.or.ANSASCFsnow_struc(n)%startMode) then 
     ANSASCFsnow_struc(n)%startMode = .false.
     
     ANSASCFsnow_struc(n)%sca = LIS_rc%udef
   
     call ANSAscf_filename(ansa_filename,ANSASCFsnow_struc(n)%scf_fn_conv,obsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)       

     inquire(file=ansa_filename,exist=file_exists)
     if(file_exists) then 

        write(LIS_logunit,*)  'Reading ANSA SCF data for EnKF assimilation',trim(ansa_filename)
        allocate(dims(2))

        dims(1) = ANSASCFsnow_struc(n)%nc
        dims(2) = ANSASCFsnow_struc(n)%nr
        offset_file = (/ANSASCFsnow_struc(n)%offset1, &
             ANSASCFsnow_struc(n)%offset2/)
        
        allocate(sca_field(ANSASCFsnow_struc(n)%nc, &
             ANSASCFsnow_struc(n)%nr))
        
        call h5open_f(status)
        call LIS_verify(status, 'Error opening HDF fortran interface')
        
        call h5fopen_f(trim(ansa_filename),H5F_ACC_RDONLY_F, file_id, status)
        call LIS_verify(status, 'Error opening ANSA file ')

        !call h5dopen_f(file_id,sca_field_name,sca_field_id, status)
        call h5dopen_f(file_id,ANSASCFsnow_struc(n)%scaname,sca_field_id, status)
        call LIS_verify(status, 'Error opening SCF field in ANSA file')

        call h5dget_space_f(sca_field_id, dataspace, status)
        call LIS_verify(status, 'Error in h5dget_space_f: read_ANSASCFsnow')
 
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
             start=offset_file, count=count_file, hdferr=status)
        call LIS_verify(status, 'Error setting hyperslab dataspace in read_ANSASCFsnow')
        
        call h5screate_simple_f(memrank,dimsm, memspace, status)
        call LIS_verify(status, 'Error in h5create_simple_f; read_ANSASCFsnow')
        
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
             start=offset_mem, count=count_mem, hdferr=status)
        call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_ANSASCFsnow')

        call h5dread_f(sca_field_id, H5T_NATIVE_INTEGER,sca_field,&
             dims,status, memspace, dataspace)
        call LIS_verify(status, 'Error extracting SCF field from ANSA file')
        
        call h5dclose_f(sca_field_id,status)
        call LIS_verify(status,'Error in H5DCLOSE call')

        do r=1,ANSASCFsnow_struc(n)%nr
           do c=1,ANSASCFsnow_struc(n)%nc
              tsnow(c+(r-1)*ANSASCFsnow_struc(n)%nc) = real(sca_field(c,r))
           enddo
        enddo

        li  = .false.
        do c=1,ANSASCFsnow_struc(n)%mi
            if (tsnow(c).le.100 .and. tsnow(c).ge.0) then
              li(c) = .true. 
           endif
        enddo
        
        call bilinear_interp(LIS_rc%gridDesc(n,:),li,tsnow,&
             lo,ANSASCFsnow_struc(n)%sca,&
             ANSASCFsnow_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             ANSASCFsnow_struc(n)%w11,ANSASCFsnow_struc(n)%w12, &
             ANSASCFsnow_struc(n)%w21,ANSASCFsnow_struc(n)%w22, &
             ANSASCFsnow_struc(n)%n11,ANSASCFsnow_struc(n)%n12, &
             ANSASCFsnow_struc(n)%n21,ANSASCFsnow_struc(n)%n22, &
             LIS_rc%udef,iret)

        deallocate(sca_field)

        call h5fclose_f(file_id,status)
        call LIS_verify(status,'Error in H5FCLOSE call')
        
        call h5close_f(status)
        call LIS_verify(status,'Error in H5CLOSE call')

        deallocate(dims)
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
  dataflag_local = .false. 
  
  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
     
!localtime of this gridcell
           lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(c,r))%lon
           call LIS_localtime(LIS_rc%gmt,lon,lhour,zone)

           if(lhour.le.lhour1 .and. (lhour1-lhour).lt.LIS_rc%ts/3600.0) then
           !if(lhour.eq.ANSASCFsnow_struc(n)%assim_lhour) then            
              obsl(LIS_domain(n)%gindex(c,r))=&
                   ANSASCFsnow_struc(n)%sca(c+LIS_rc%lnc(n)*(r-1))
              if(obsl(LIS_domain(n)%gindex(c,r)).ne.-9999.0) then
                 dataflag_local = .true.
              endif
           endif
        end if
     end do
  end do


! LSM-based QC
!  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
!       //trim(LIS_ANSASCFsnowobsId)//char(0), & 
!       n, OBS_state)

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
     call LIS_verify(status, 'ESMF_StateGet for Observation01 for OBS_Pert_State failed in read_ANSASCFsnow')    

     if(LIS_rc%ngrid(n).gt.0) then 

!linearly scale the observation err
        ssdev = ANSASCFsnow_struc(n)%ssdev 
        if (ANSASCFsnow_struc(n)%obserr .eq. 1) then
           do t=1,LIS_rc%ngrid(n)
              if(obsl(t).ne.-9999.0) then 
                 if (obsl(t).le.50) then
                     ssdev(t) = ssdev(t)*obsl(t)+0.001
                 elseif (obsl(t).gt.50) then
                     ssdev(t) = ssdev(t)*(100-obsl(t))+0.001
                 endif
              endif
           enddo
        endif

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
  
#endif     
end subroutine read_ANSASCFsnow

!BOP
!
! !ROUTINE: ANSAscf_filename
! \label{ANSAscf_filename}
! 
! !INTERFACE: 
subroutine ANSAscf_filename(name, scf_fn_conv, ndir, yr, mo,da)

  use LIS_coreMod, only: LIS_rc
  use LIS_logMod, only: LIS_logunit

  implicit none
! !ARGUMENTS:
  character (len=*) :: name
  integer           :: yr, mo, da
  character (len=*) :: ndir
  character (len=*) :: scf_fn_conv
  integer           :: len1, i1, i2
  character (len=50) :: str1, str2, str0
!
! !DESCRIPTION:
!  This subroutine creates a timestamped ANSA filename
!
!  The arguments are:
!  \begin{description}
!  \item[name] name of the ANSA SCF filename
!  \item[scf\_fn\_conv] ANSA SCF data file name convention
!  \item[ndir] name of the ANSA SCF data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  character (len=3) :: fdoy
 
  len1 = len_trim(scf_fn_conv)
  i1 = index(scf_fn_conv, 'YYYYMMDD')
  i2 = index(scf_fn_conv, 'YYYYDOY')

  if (i1 .gt. 0) then
    str1 = scf_fn_conv(1:i1-1)
    str2 = scf_fn_conv(i1+8:len1)
    write(unit=fyr, fmt='(i4.4)') yr
    write(unit=fmo, fmt='(i2.2)') mo
    write(unit=fda, fmt='(i2.2)') da
    str0 = fyr//fmo//fda
  else if (i2 .gt. 0) then
    str1 = scf_fn_conv(1:i2-1)
    str2 = scf_fn_conv(i2+7:len1)
    write(unit=fyr, fmt='(i4.4)') yr
    write(unit=fdoy, fmt='(i3.3)') LIS_rc%doy
    str0 = fyr//fdoy
  else
    write(LIS_logunit,*) 'ANSA SCF data file name convention currently not supported: ', trim(scf_fn_conv)
  end if

  name = trim(ndir)//'/'//trim(fyr)//'/'//trim(str1)//trim(str0)//trim(str2)

end subroutine ANSAscf_filename

