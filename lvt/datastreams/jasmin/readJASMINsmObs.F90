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
! !ROUTINE: readJASMINsmObs
! \label{readJASMINsmObs}
!
! !INTERFACE: 
subroutine readJASMINsmObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit, LVT_getNextUnitNumber, & 
       LVT_releaseUnitNumber
  use LVT_timeMgrMod,   only : LVT_get_julss
  use JASMINsm_obsMod, only : JASMINsmobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)      :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This subroutine provides the data reader for the standard 
! LPRM soil moisture retrieval product. 
! 
! !NOTES: 
!  The mismatches between the LVT time and LPRM data time is not
!  handled. This is not an issue as long as LVT analysis is not
!  done for sub-daily intervals. 
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
!  25 May 2012: Sujay Kumar, Updated for LPRM version 5. 
! 
!EOP
  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character*100     :: fname
  real              :: smobs(LVT_rc%lnc*LVT_rc%lnr,4)
  real              :: rootsm(LVT_rc%lnc*LVT_rc%lnr)
  real              :: lat,lon

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
       LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  smobs= LVT_rc%udef
  rootsm = LVT_rc%udef

  if(JASMINsmobs(source)%startmode.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     JASMINsmobs(source)%startmode = .false. 

     call create_JASMINsm_filename(JASMINsmobs(source)%odir, &
          LVT_rc%dyr(source),  fname)

     inquire(file=trim(fname),exist=file_exists)

     if(file_exists) then
        write(LVT_logunit,*) '[INFO] Reading JASMIN file: ',trim(fname)
        call read_JASMIN_data(source, fname,smobs)
     else
        write(LVT_logunit,*) '[WARN] Unable to read-in JASMIN file:',trim(fname)
        write(LVT_logunit,*) '[WARN] Note: Check length of filepath, as a long '
        write(LVT_logunit,*) '[WARN] pathname can lead to not reading the file.' 
     endif
     
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(smobs(c+(r-1)*LVT_rc%lnc,1).ne.LVT_rc%udef) then 
              rootsm(c+(r-1)*LVT_rc%lnc) = & 
                   smobs(c+(r-1)*LVT_rc%lnc,1)*JASMINsmobs(source)%rz_wt(1)+&
                   smobs(c+(r-1)*LVT_rc%lnc,2)*JASMINsmobs(source)%rz_wt(2)+&
                   smobs(c+(r-1)*LVT_rc%lnc,3)*JASMINsmobs(source)%rz_wt(3)+&
                   smobs(c+(r-1)*LVT_rc%lnc,4)*JASMINsmobs(source)%rz_wt(4)
                   
           endif
        enddo
     enddo
  endif
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, &
       smobs(:,1),vlevel=1,units="m3/m3")
  call LVT_logSingleDataStreamVar(LVT_MOC_rootmoist, source, &
       rootsm,vlevel=1,units="m3/m3")

end subroutine readJASMINsmObs


!BOP
! 
! !ROUTINE: read_JASMIN_data
! \label(read_JASMIN_data)
!
! !INTERFACE:
subroutine read_JASMIN_data(source, fname,smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use ESMF
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_logMod,   only : LVT_verify
  use map_utils,    only : latlon_to_ij
  use JASMINsm_obsMod, only : JASMINsmobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer              :: source
  character (len=*)    :: fname
  real                 :: smobs_ip(LVT_rc%lnc*LVT_rc%lnr,4)

!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the JASMIN NETCDF file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  land surface temperature is below freezing, if rain is present, if 
!  RFI is present, if residual error is above 0.5 or if optical depth
!  is above 0.8. Finally the routine combines both the C-band and X-band
!  retrievals. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the JASMIN AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LVT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP


  real           :: sm(JASMINsmobs(source)%jasminnc,&
       JASMINsmobs(source)%jasminnr,4)

  real           :: sm_data(JASMINsmobs(source)%jasminnc*&
       JASMINsmobs(source)%jasminnr,4)
  logical*1      :: sm_data_b(JASMINsmobs(source)%jasminnc*JASMINsmobs(source)%jasminnr)
  logical*1      :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)

  type(ESMF_Time) :: currTime
  real, allocatable :: time(:)
  integer        :: c,r,i,j,k,kk
  real           :: rlat,rlon,ri,rj
  integer        :: nid,dimid,tid
  integer        :: ndims
  integer        :: offset, offset_index
  logical        :: find_offset
  integer        :: smid, flagid
  integer        :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  sm = LVT_rc%udef

  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LVT_verify(ios,'Error opening file '//trim(fname))

  ios = nf90_inq_dimid(nid, 'time',dimId)
  call LVT_verify(ios, 'Error nf90_inq_dimid: time')

  ios = nf90_inquire_dimension(nid, dimId, len=ndims)
  call LVT_verify(ios, 'Error nf90_inquire_dimension:')

  ios = nf90_inq_varid(nid,'time',tid)
  call LVT_verify(ios,'Error nf90_inq_varid: time')

  allocate(time(ndims))

  ios = nf90_get_var(nid,tid,time)
  call LVT_verify(ios,'Error nf90_get_var: time')

  ios = nf90_inq_varid(nid, 'sm',smid)
  call LVT_verify(ios, 'Error nf90_inq_varid: sm')

  !find the offset
  
  call ESMF_TimeSet(currTime, yy=LVT_rc%yr,&
       mm = LVT_rc%mo,&
       dd = LVT_rc%da,&
       h = LVT_rc%hr, &
       m = LVT_rc%mn, &
       calendar = LVT_calendar, &
       rc=ios)
  call LVT_verify(ios,'error in timeset: readJASMINsmobs')

  offset = nint((currTime-JASMINsmobs(source)%refTime)/&
       JASMINsmobs(source)%dt)+1

  find_offset = .false. 
  do k=1,ndims
     if(time(k).eq.offset) then 
        find_offset = .true. 
        offset_index = k
     endif
  enddo
  if(find_offset) then 
     ios = nf90_get_var(nid, smid, sm, &
          start=(/1,1,1,offset_index/), &
          count=(/JASMINsmobs(source)%jasminnc,JASMINsmobs(source)%jasminnr,&
          4,1/))
     call LVT_verify(ios, 'Error nf90_get_var: sm')
  endif
     
  ios = nf90_close(ncid=nid)
  call LVT_verify(ios,'Error closing file '//trim(fname))

  do r=1, JASMINsmobs(source)%jasminnr
     do c=1, JASMINsmobs(source)%jasminnc
        sm_data(c+(r-1)*JASMINsmobs(source)%jasminnc,:) = sm(c,r,:)
     enddo
  enddo
  deallocate(time)
  do r=1, JASMINsmobs(source)%jasminnr
     do c=1, JASMINsmobs(source)%jasminnc
        if(sm_data(c+(r-1)*JASMINsmobs(source)%jasminnc,1).gt.0.0.and.&
             sm_data(c+(r-1)*JASMINsmobs(source)%jasminnc,1).lt.1.0) then 
           sm_data_b(c+(r-1)*JASMINsmobs(source)%jasminnc) = .true. 
        else
           sm_data_b(c+(r-1)*JASMINsmobs(source)%jasminnc) = .false.
        endif
     enddo
  enddo
  
  do kk=1,4
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_data(:,kk), smobs_b_ip, smobs_ip(:,kk), &
          JASMINsmobs(source)%jasminnc*JASMINsmobs(source)%jasminnr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          JASMINsmobs(source)%rlat, JASMINsmobs(source)%rlon, &
          JASMINsmobs(source)%w11, JASMINsmobs(source)%w12, &
          JASMINsmobs(source)%w21, JASMINsmobs(source)%w22, &
          JASMINsmobs(source)%n11, JASMINsmobs(source)%n12, &
          JASMINsmobs(source)%n21, JASMINsmobs(source)%n22, &
          LVT_rc%udef, ios)
  enddo

#endif
  
end subroutine read_JASMIN_data

!BOP
! !ROUTINE: create_JASMINsm_filename
! \label{create_JASMINsm_filename}
! 
! !INTERFACE: 
subroutine create_JASMINsm_filename(ndir, yr, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the JASMIN filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the JASMIN soil moisture directory
!  \item[yr]  current year
!  \item[filename] Generated JASMIN filename
! \end{description}
!EOP

  character (len=4) :: fyr
  
  write(unit=fyr, fmt='(i4.4)') yr
  
  filename = trim(ndir)//'/'//&
       'jasmin.vol.smc.'//trim(fyr)//'.nc'

end subroutine create_JASMINsm_filename
