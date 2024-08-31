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
! !ROUTINE: readANSASWEobs
! \label{readANSASWEobs}
!
! !INTERFACE: 
subroutine readANSASWEobs(source)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the standard 
! VU soil moisture retrieval product. 
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
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit, LVT_getNextUnitNumber, & 
       LVT_releaseUnitNumber, LVT_verify
  use LVT_timeMgrMod,   only : LVT_get_julss, LVT_localtime
  use ANSASWE_obsMod, only : ANSASWEobs

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: source
!
!EOP
#if (defined USE_HDF5) 
  character*100,   parameter    :: swe_field_name = "ansa_swe_cyl_GB"
  character*100                 :: ansa_filename
  integer(hsize_t), allocatable :: dims(:)
  integer(hid_t)                :: dataspace
  integer(hid_t)                :: memspace
  integer(hsize_t), dimension(2) :: dimsm 
  integer                       :: memrank = 2
  integer(hid_t)                :: file_id, swe_field_id
  integer(hsize_t), dimension(2) :: count_file 
  integer(hsize_t), dimension(2) :: count_mem 
  integer(hsize_t), dimension(2) :: offset_file 
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer                       :: c,r
  real                          :: lon, lhour
  integer                       :: zone
  integer, allocatable, target  :: swe_field(:,:)
  real                          :: tsnow(ANSASWEobs(source)%nc*ANSASWEobs(source)%nr)
  integer                       :: iret,status
  logical*1                     :: lo(LVT_rc%lnc*LVT_rc%lnr)
  logical*1                     :: li(ANSASWEobs(source)%nc*ANSASWEobs(source)%nr)
  real                          :: swe(LVT_rc%lnc, LVT_rc%lnr)
  logical                       :: alarmCheck
  logical                       :: file_exists
  real                          :: currtime

  dimsm      = (/ANSASWEobs(source)%nc, ANSASWEobs(source)%nr/)
  count_file = (/ANSASWEobs(source)%nc, ANSASWEobs(source)%nr/)
  count_mem  = (/ANSASWEobs(source)%nc, ANSASWEobs(source)%nr/)

  currtime = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(currtime, 86400.0).eq.0)
  if(ANSASWEobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     ANSASWEobs(source)%startflag = .false. 
     call ANSASWE_filename(ansa_filename,ANSASWEobs(source)%odir,&
          LVT_rc%dyr(source),LVT_rc%dmo(source),LVT_rc%dda(source))       
     inquire(file=trim(ansa_filename), exist=file_exists) 
     
     if(file_exists) then 
        
        write(LVT_logunit,*) '[INFO] Reading ANSA SWE data', ansa_filename

        allocate(dims(2))
        dims(1) = ANSASWEobs(source)%nc
        dims(2) = ANSASWEobs(source)%nr

        offset_file = (/ANSASWEobs(source)%offset1, ANSASWEobs(source)%offset2/)
        
        allocate(swe_field(ANSASWEobs(source)%nc, ANSASWEobs(source)%nr))
        
        call h5open_f(status)
        call LVT_verify(status, 'Error opening HDF fortran interface')
        
        call h5fopen_f(trim(ansa_filename),H5F_ACC_RDONLY_F, file_id, status)
        call LVT_verify(status, 'Error opening ANSA file ')
        
        call h5dopen_f(file_id,swe_field_name,swe_field_id, status)
        call LVT_verify(status, 'Error opening SWE field in ANSA file')
        
        call h5dget_space_f(swe_field_id, dataspace, status)
        call LVT_verify(status, 'Error in h5dget_space_f: readANSASWE')
 
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
             start=offset_file, count=count_file, hdferr=status)
        call LVT_verify(status, 'Error setting hyperslab dataspace in readANSASWE')
        
        call h5screate_simple_f(memrank,dimsm, memspace, status)
        call LVT_verify(status, 'Error in h5create_simple_f; readANSASWE')
        
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
             start=offset_mem, count=count_mem, hdferr=status)
        call LVT_verify(status, 'Error in h5sselect_hyperslab_f: readANSASWE')

        call h5dread_f(swe_field_id, H5T_NATIVE_INTEGER,swe_field,dims,status, &
             memspace, dataspace)
        call LVT_verify(status, 'Error extracting SWE field from ANSA file')
        
        call h5dclose_f(swe_field_id,status)
        call LVT_verify(status,'Error in H5DCLOSE call')
        
        do r=1,ANSASWEobs(source)%nr
           do c=1,ANSASWEobs(source)%nc
              tsnow(c+(r-1)*ANSASWEobs(source)%nc) = swe_field(c,r)
           enddo
        enddo

        li  = .false.
        do c=1,ANSASWEobs(source)%mi
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
             lo,ANSASWEobs(source)%swe,&
             ANSASWEobs(source)%mi,LVT_rc%lnc*LVT_rc%lnr,&
             ANSASWEobs(source)%rlat,ANSASWEobs(source)%rlon,&
             ANSASWEobs(source)%w11,ANSASWEobs(source)%w12,&
             ANSASWEobs(source)%w21,ANSASWEobs(source)%w22,&
             ANSASWEobs(source)%n11,ANSASWEobs(source)%n12,&
             ANSASWEobs(source)%n21,ANSASWEobs(source)%n22,&
             LVT_rc%udef,iret)

        deallocate(swe_field)
        call h5fclose_f(file_id,status)
        call LVT_verify(status,'Error in H5FCLOSE call')
        
        call h5close_f(status)
        call LVT_verify(status,'Error in H5CLOSE call')

        deallocate(dims)        
       
        write(LVT_logunit,*) '[INFO] Finished processing ', ansa_filename 
     endif
     
  endif
  
  swe = LVT_rc%udef

  do r=1,LVT_rc%lnr
     do c=1, LVT_rc%lnc
        if(LVT_domain%gindex(c,r).ne.-1) then 
           lon = LVT_domain%grid(LVT_domain%gindex(c,r))%lon
           call LVT_localtime(LVT_rc%gmt, lon, lhour, zone)
           if(lhour.ge.(14-nint(LVT_rc%ts/3600.0)).and.lhour.le.14) then 
              if(ANSASWEobs(source)%swe(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
                 swe(c,r) = ANSASWEobs(source)%swe(c+LVT_rc%lnc*(r-1))/1000.0
              endif
           endif
        endif
     enddo
  enddo


  call LVT_logSingleDataStreamVar(LVT_MOC_swe, source, swe,vlevel=1,units="m")

  do r=1,LVT_rc%lnr
     do c=1, LVT_rc%lnc
        if(swe(c,r).ne.-9999.0) then 
           swe(c,r) = swe(c,r)*1000.0
        endif
     enddo
  enddo


  call LVT_logSingleDataStreamVar(LVT_MOC_swe, source, swe,vlevel=1,units="kg/m2")
#endif
 
end subroutine readANSASWEobs

!BOP
! 
! !ROUTINE: ANSAswe_filename
! \label{ANSAswe_filename}
!
! !INTERFACE: 
subroutine ANSAswe_filename(name, ndir, yr, mo,da)
  
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
  character*80      :: name
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da 

  name = trim(ndir)//'/'//trim(fyr)//'/ansa_all_'//trim(fyr)//trim(fmo)//trim(fda)//'.h5'
    
end subroutine ANSAswe_filename


