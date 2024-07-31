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
! !ROUTINE: readANSASNWDobs
! \label{readANSASNWDobs}
!
! !INTERFACE: 
subroutine readANSASNWDobs(source)
! 
! !USES:   
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
!  21 July 2010: Sujay Kumar, Initial Specification
! 
!EOP
!BOP
! 
! 
! 
! !USES:   
#if(defined USE_HDF5) 
  use hdf5
#endif
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use ANSASNWD_obsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: source
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the standard 
! VU soil moisture retrieval product. 
!
!EOP
#if (defined USE_HDF5) 
  character*100,   parameter    :: snwd_field_name = "ansa_swe_depth_cyl_GB"
  character*100,   parameter    :: snwd_field_name1 = "ansa_swe_depth_cyl_NH"
  character*100                 :: ansa_filename
  integer(hsize_t), allocatable :: dims(:)
  integer(hid_t)                :: dataspace
  integer(hid_t)                :: memspace
  integer(hsize_t), dimension(2) :: dimsm 
  integer                       :: memrank = 2
  integer(hid_t)                :: file_id, snwd_field_id
  integer(hsize_t), dimension(2) :: count_file 
  integer(hsize_t), dimension(2) :: count_mem 
  integer(hsize_t), dimension(2) :: offset_file 
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer                       :: c,r
  real                          :: lon, lhour
  integer                       :: zone
  integer, allocatable, target  :: snwd_field(:,:)
  real                          :: tsnow(ANSASNWDobs(source)%nc*ANSASNWDobs(source)%nr)
  integer                       :: iret,status
  logical*1                     :: lo(LVT_rc%lnc*LVT_rc%lnr)
  logical*1                     :: li(ANSASNWDobs(source)%nc*ANSASNWDobs(source)%nr)
  real                          :: snwd(LVT_rc%lnc, LVT_rc%lnr)
  logical                       :: alarmCheck
  logical                       :: file_exists
  real                          :: currtime

  dimsm      = (/ANSASNWDobs(source)%nc, ANSASNWDobs(source)%nr/)
  count_file = (/ANSASNWDobs(source)%nc, ANSASNWDobs(source)%nr/)
  count_mem  = (/ANSASNWDobs(source)%nc, ANSASNWDobs(source)%nr/)

  currtime = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(currtime, 86400.0).eq.0)
  if(ANSASNWDobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     LVT_rc%resetFlag(source) = .false. 

     ANSASNWDobs(source)%startflag = .false. 
     call ANSASNWD_filename(ansa_filename,ANSASNWDobs(source)%odir,&
          LVT_rc%dyr(source),LVT_rc%dmo(source),LVT_rc%dda(source))       
     inquire(file=trim(ansa_filename), exist=file_exists) 

     if(file_exists) then 
        
        write(LVT_logunit,*) '[INFO] Reading ANSA SNWD data', trim(ansa_filename)

        allocate(dims(2))
        dims(1) = ANSASNWDobs(source)%nc
        dims(2) = ANSASNWDobs(source)%nr
        
        offset_file = (/ANSASNWDobs(source)%offset1, ANSASNWDobs(source)%offset2/)
        
        allocate(snwd_field(ANSASNWDobs(source)%nc, ANSASNWDobs(source)%nr))
        
        call h5open_f(status)
        call LVT_verify(status, 'Error opening HDF fortran interface')
        
        call h5fopen_f(trim(ansa_filename),H5F_ACC_RDONLY_F, file_id, status)

        call LVT_verify(status, 'Error opening ANSA file ')

        if(LVT_rc%dyr(source).ge.2010) then 
           call h5dopen_f(file_id,snwd_field_name1,snwd_field_id, status)
           call LVT_verify(status, 'Error opening SNWD field in ANSA file')
        else
           call h5dopen_f(file_id,snwd_field_name,snwd_field_id, status)
           call LVT_verify(status, 'Error opening SNWD field in ANSA file')
        endif

        call h5dget_space_f(snwd_field_id, dataspace, status)
        call LVT_verify(status, 'Error in h5dget_space_f: readANSASNWD')
 
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
             start=offset_file, count=count_file, hdferr=status)
        call LVT_verify(status, 'Error setting hyperslab dataspace in readANSASNWD')
        
        call h5screate_simple_f(memrank,dimsm, memspace, status)
        call LVT_verify(status, 'Error in h5create_simple_f; readANSASNWD')
        
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
             start=offset_mem, count=count_mem, hdferr=status)
        call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASNWD')

        call h5dread_f(snwd_field_id, H5T_NATIVE_INTEGER,snwd_field,dims,status, &
             memspace, dataspace)
        call LVT_verify(status, 'Error extracting SNWD field from ANSA file')
        
        call h5dclose_f(snwd_field_id,status)
        call LVT_verify(status,'Error in H5DCLOSE call')
        
        do r=1,ANSASNWDobs(source)%nr
           do c=1,ANSASNWDobs(source)%nc
              tsnow(c+(r-1)*ANSASNWDobs(source)%nc) = real(snwd_field(c,r))
           enddo
        enddo

        li  = .false.
        do c=1,ANSASNWDobs(source)%mi
!           if(tsnow(c).lt.1500.and.tsnow(c).ge.0) then 
!              li(c) = .true. 
!           endif
           if(tsnow(c).ne.494.and.tsnow(c).ne.496.and.&
                tsnow(c).ne.504.and.tsnow(c).ne.506.and.&
                tsnow(c).ne.508.and.tsnow(c).ne.510.and.&
                tsnow(c).lt.2000.and.tsnow(c).ge.0) then 
              li(c) = .true. 
           endif
        enddo


        li  = .false.
        do c=1,ANSASNWDobs(source)%mi
!           if(tsnow(c).lt.1500.and.tsnow(c).ge.0) then 
!              li(c) = .true. 
!           endif
           if(tsnow(c).ne.494.and.tsnow(c).ne.496.and.&
                tsnow(c).ne.504.and.tsnow(c).ne.506.and.&
                tsnow(c).ne.508.and.tsnow(c).ne.510.and.&
                tsnow(c).lt.2000.and.tsnow(c).ge.0) then 
              li(c) = .true. 
           endif
        enddo

        call bilinear_interp(LVT_rc%gridDesc,li,tsnow,&
             lo,ANSASNWDobs(source)%snwd,&
             ANSASNWDobs(source)%mi,LVT_rc%lnc*LVT_rc%lnr,&
             ANSASNWDobs(source)%rlat,ANSASNWDobs(source)%rlon,&
             ANSASNWDobs(source)%w11,ANSASNWDobs(source)%w12,&
             ANSASNWDobs(source)%w21,ANSASNWDobs(source)%w22,&
             ANSASNWDobs(source)%n11,ANSASNWDobs(source)%n12,&
             ANSASNWDobs(source)%n21,ANSASNWDobs(source)%n22,&
             LVT_rc%udef,iret)

        deallocate(snwd_field)
        call h5fclose_f(file_id,status)
        call LVT_verify(status,'Error in H5FCLOSE call')
        
        call h5close_f(status)
        call LVT_verify(status,'Error in H5CLOSE call')

        deallocate(dims)        
        
        write(LVT_logunit,*) '[INFO] Finished processing ', trim(ansa_filename)
     endif
     
  endif
  
  snwd = LVT_rc%udef
  do r=1,LVT_rc%lnr
     do c=1, LVT_rc%lnc
        if(LVT_domain%gindex(c,r).ne.-1) then 
           if(ANSASNWDobs(source)%snwd(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
              snwd(c,r) = ANSASNWDobs(source)%snwd(c+LVT_rc%lnc*(r-1))/1000.0
           endif
        endif
     enddo
  enddo



 !assume daily analysis
#if 0
  do r=1,LVT_rc%lnr
     do c=1, LVT_rc%lnc
        if(LVT_domain%gindex(c,r).ne.-1) then 
           lon = LVT_domain%grid(LVT_domain%gindex(c,r))%lon
           call LVT_localtime(LVT_rc%gmt, lon, lhour, zone)
           if(lhour.ge.(14-nint(LVT_rc%ts/3600.0)).and.lhour.le.14) then 
              if(ANSASNWDobs(source)%snwd(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
                 snwd(c,r) = ANSASNWDobs(source)%snwd(c+LVT_rc%lnc*(r-1))/1000.0
!                 print*, 'ansa ',c,r,snwd(c,r)
              endif
           endif
        endif
     enddo
  enddo
#endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH, source, snwd,vlevel=1,&
       units="m")
#endif
 
end subroutine readANSASNWDobs

!BOP
! 
! !ROUTINE: ANSAsnwd_filename
! \label{ANSAsnwd_filename}
!
! !INTERFACE: 
subroutine ANSAsnwd_filename(name, ndir, yr, mo,da)
  
  implicit none
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
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
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

!BOP
!
! 
! !ARGUMENTS: 
  character(len=*) :: name
  integer          :: yr, mo, da
  character(len=*) :: ndir
! 
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da 

  name = trim(ndir)//'/'//trim(fyr)//'/ansa_all_'//trim(fyr)//trim(fmo)//trim(fda)//'.h5'
    
end subroutine ANSAsnwd_filename


