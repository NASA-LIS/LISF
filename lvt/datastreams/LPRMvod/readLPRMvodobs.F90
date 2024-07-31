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
! !ROUTINE: readLPRMvodobs
! \label{readLPRMvodobs}
!
! !INTERFACE: 
subroutine readLPRMvodobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use LPRM_vodobsMod, only : LPRM_vodobs
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

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
! Land Parameter Retrieval Model (LPRM) project
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
  real              :: vod_file(LPRM_vodobs(source)%nc,LPRM_vodobs(source)%nr)
  real              :: vod_inp(LPRM_vodobs(source)%nc*LPRM_vodobs(source)%nr)
  logical*1         :: vod_b_inp(LPRM_vodobs(source)%nc*LPRM_vodobs(source)%nr)
!  real              :: flag(LPRM_vodobs(source)%nc,LPRM_vodobs(source)%nr)
  real              :: vod_out(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: vod_b_out(LVT_rc%lnc*LVT_rc%lnr)
  real              :: lat(LPRM_vodobs(source)%nr)
  real              :: lon(LPRM_vodobs(source)%nc)
  integer           :: fnd 
  integer                :: ftn,ios
  integer                :: vodid,flagid,latid,lonid
  character*100          :: lprm_filename
  integer                :: yr, mo, da, hr, mn, ss
  real                   :: gmt
  integer                :: c,r,c1,r1
  integer                :: doy
  real*8                 :: timenow
  real                   :: vodobs(LVT_rc%lnc,LVT_rc%lnr)

  vodobs = LVT_rc%udef
  vod_out = LVT_rc%udef

  timenow = float(LVT_rc%dhr(source))*3600 +&
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  if(LPRM_vodobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 

     LPRM_vodobs(source)%startflag = .false. 
     
     call create_LPRM_vodfilename(&
          LPRM_vodobs(source)%odir,&
          LPRM_vodobs(source)%data_designation,&
          LVT_rc%dyr(source),&
          LVT_rc%dmo(source),&
          LVT_rc%dda(source),&
          lprm_filename)

     inquire(file=lprm_filename,exist=file_exists)
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] reading ',trim(lprm_filename)
        
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        ios = nf90_open(path=trim(lprm_filename),mode=NF90_NOWRITE,ncid=ftn)
        call LVT_verify(ios,'Error opening file '//trim(lprm_filename))
        
        ios = nf90_inq_varid(ftn, 'lat',latid)
        call LVT_verify(ios, 'Error nf90_inq_varid: lat')

        ios = nf90_inq_varid(ftn, 'lon',lonid)
        call LVT_verify(ios, 'Error nf90_inq_varid: lon')

        ios = nf90_inq_varid(ftn, 'vod',vodid)
        call LVT_verify(ios, 'Error nf90_inq_varid: vod')
        
!        ios = nf90_inq_varid(ftn, 'sensor_flag',flagid)
!        call LVT_verify(ios, 'Error nf90_inq_varid: sensor_flag')
  
        !values
        ios = nf90_get_var(ftn, latid, lat)
        call LVT_verify(ios, 'Error nf90_get_var: lat')

        ios = nf90_get_var(ftn, lonid, lon)
        call LVT_verify(ios, 'Error nf90_get_var: lon')

        ios = nf90_get_var(ftn, vodid, vod_file)
        call LVT_verify(ios, 'Error nf90_get_var: vod')
        
!        ios = nf90_get_var(ftn, flagid, flag)
!        call LVT_verify(ios, 'Error nf90_get_var: flag')
        
        ios = nf90_close(ncid=ftn)
        call LVT_verify(ios,'Error closing file '//trim(lprm_filename))
#endif

        vod_inp = LVT_rc%udef
        vod_b_inp  = .false. 


        do r=1,LPRM_vodobs(source)%nr
           do c=1,LPRM_vodobs(source)%nc
              r1 = nint((lat(r)+89.875)/0.25)+1
              c1 = nint((lon(c)+179.875)/0.25)+1

              if(vod_file(c,r).ne.-999999.0) then 
                 vod_inp(c1+(r1-1)*LPRM_vodobs(source)%nc) = & 
                      vod_file(c,r)
                 vod_b_inp(c1+(r1-1)*LPRM_vodobs(source)%nc) = & 
                      .true. 
              endif
           enddo
        enddo

!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
        call bilinear_interp(LVT_rc%gridDesc(:),&
             vod_b_inp, vod_inp, vod_b_out, vod_out, &
             LPRM_vodobs(source)%nc*LPRM_vodobs(source)%nr, &
             LVT_rc%lnc*LVT_rc%lnr, &
             LPRM_vodobs(source)%rlat, LPRM_vodobs(source)%rlon, &
             LPRM_vodobs(source)%w11, LPRM_vodobs(source)%w12, &
             LPRM_vodobs(source)%w21, LPRM_vodobs(source)%w22, &
             LPRM_vodobs(source)%n11, LPRM_vodobs(source)%n12, &
             LPRM_vodobs(source)%n21, LPRM_vodobs(source)%n22, &
             LVT_rc%udef, ios)

        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              vodobs(c,r) = vod_out(c+(r-1)*LVT_rc%lnc)
           enddo
        enddo
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_VOD, source,&
       vodobs,vlevel=1,units="-")
 
end subroutine readLPRMvodobs

!BOP
! 
! !ROUTINE: create_LPRM_vodfilename
! 
! !INTERFACE: 
subroutine create_LPRM_vodfilename(&
     odir, data_designation,yr, mo, da, filename)
!
! !DESCRIPTION: 
! This routine creates the timestamped filename for the
! LPRM vegetation optical depth dataset
! 
!EOP
  implicit none

  character(len=*) :: odir
  character(len=*) :: data_designation
  integer          :: yr 
  integer          :: mo
  integer          :: da
  character(len=*) :: filename
  
  character*4      :: yyyy
  character*2      :: mm,dd

  write(yyyy,'(i4.4)') yr
  write(mm,'(i2.2)') mo
  write(dd,'(i2.2)') da
  
  filename=trim(odir)//'/'//trim(data_designation)//'/'//trim(yyyy)//&
       '/vodca_v01-0_'//trim(data_designation)//'_'//trim(yyyy)//'-'//&
       trim(mm)//'-'//trim(dd)//'.nc'
       
end subroutine create_LPRM_vodfilename
