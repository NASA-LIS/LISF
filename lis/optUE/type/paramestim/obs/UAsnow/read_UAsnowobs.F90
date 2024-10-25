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
! !ROUTINE: read_UAsnowobs
! \label{read_UAsnowobs}
!
! !REVISION HISTORY:
!  2 May 2020  Sujay Kumar;   Initial Specification
!
! !INTERFACE: 
subroutine read_UAsnowobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_fileIOMod
  use UAsnow_obsMod
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

#if (defined USE_NETCDF3 || defined USE_NETCDF4)  
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_Space] Objective Space
!  \end{description}
!
!EOP
  integer                  :: n
  real,    pointer         :: snod(:)
  type(ESMF_Field)         :: snodField
  real,    pointer         :: swe(:)
  type(ESMF_Field)         :: sweField
  character(len=LIS_CONST_PATH_LEN) :: obsdir
  integer                  :: status
  integer                  :: nid,sweid,snwdid
  integer                  :: c,r,k,nt,iret
  integer                  :: grid_index
  integer                  :: toffset
  type(ESMF_Time)          :: uatime1
  logical                  :: alarmCheck
  logical                  :: file_exists
  logical                  :: data_update
  character(len=LIS_CONST_PATH_LEN) :: fname
  real                     :: dt
  real                     :: timenow
  real, allocatable        :: snwd_out(:)
  real, allocatable        :: swe_out(:)
  logical*1, allocatable   :: lo(:)
  logical*1, allocatable   :: lb(:)
  real, allocatable        :: UA_swe(:,:), UA_snwd(:,:)
  real, allocatable        :: swe_in(:)
  real, allocatable        :: snwd_in(:)

  n = 1

  call ESMF_AttributeGet(Obj_Space,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status, 'Error in ESMF_AttributeGet: Data Directory')
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status, 'Error in ESMF_AttributeGet: Data Update Status')

  allocate(snwd_out(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(swe_out(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  
  snwd_out = LIS_rc%udef
  swe_out = LIS_rc%udef

  timenow = float(LIS_rc%hr*3600+LIS_rc%mn*60+LIS_rc%ss)
  
  alarmCheck = (mod(timenow, 86400.0).eq.43200.0)

! The UA dataset is a daily product, valid approximately
! at 12Z, which is the observation time of many of the
! surface products used to generate the data (personal
! communication with the UA SNOW product generators).
  if(alarmCheck) then 
     
     allocate(UA_swe(UAsnow_obs_struc(n)%nc,&
          UAsnow_obs_struc(n)%nr))
     allocate(UA_snwd(UAsnow_obs_struc(n)%nc,&
          UAsnow_obs_struc(n)%nr))
     UA_swe = LIS_rc%udef
     UA_snwd = LIS_rc%udef

     !UA data is defined for water years, beginning in October
     call ESMF_TimeSet(UAsnow_obs_struc(n)%startTime, &
          yy = LIS_rc%yr, & 
          mm = 10, &
          dd = 1, &
          h = 0, &
          m = 0, & 
          calendar = LIS_calendar, &
          rc=status)
     call LIS_verify(status, 'Error in setting UAsnow start time')

     call create_UAsnowobs_filename(&
          obsdir, &
          LIS_rc%yr, LIS_rc%mo, nt, fname)
     
     inquire(file=trim(fname),exist=file_exists)
     if(file_exists) then 
        write(LIS_logunit,*) '[INFO] Reading UA data ',&
             trim(fname)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        call ESMF_TimeSet(uatime1,yy=LIS_rc%yr,&
             mm=LIS_rc%mo, &
             dd=LIS_rc%da, &
             h=LIS_rc%hr,&
             m=LIS_rc%mn,&
             s=LIS_rc%ss,&
             calendar=LIS_calendar,&
             rc=status)

        toffset=nint((uatime1 - UAsnow_obs_struc(n)%starttime)/&
             UAsnow_obs_struc(n)%timestep)
        if(toffset.lt.0) toffset = toffset + nt+1

        call LIS_verify(nf90_open(path=fname,mode=NF90_NOWRITE,&
             ncid = nid),'Error opening file '//trim(fname))
        call LIS_verify(nf90_inq_varid(nid,'SWE', sweid), &
             'Error nf90_inq_varid: SWE')
        call LIS_verify(nf90_inq_varid(nid,'DEPTH', snwdid),&
             'Error nf90_inq_varid: DEPTH')

        call LIS_verify(nf90_get_var(nid,sweid,UA_swe,&
             start=(/1,1,toffset/),&
             count=(/UAsnow_obs_struc(n)%nc,UAsnow_obs_struc(n)%nr,1/)),&
             'Error in nf90_get_var: SWE')

        call LIS_verify(nf90_get_var(nid,snwdid,UA_snwd,&
             start=(/1,1,toffset/),&
             count=(/UAsnow_obs_struc(n)%nc,UAsnow_obs_struc(n)%nr,1/)),&
             'Error in nf90_get_var: DEPTH')
        call LIS_verify(nf90_close(nid),&
             'Error in nf90_close')
#endif
     else
        write(LIS_logunit,*) '[WARN] UA snow file not found',&
             trim(fname)
     endif
     
     allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
     allocate(lb(UAsnow_obs_struc(n)%nc*UAsnow_obs_struc(n)%nr))
     allocate(swe_in(UAsnow_obs_struc(n)%nc*UAsnow_obs_struc(n)%nr))
     allocate(snwd_in(UAsnow_obs_struc(n)%nc*UAsnow_obs_struc(n)%nr))

     lb = .false. 
     swe_in = LIS_rc%udef
     snwd_in = LIS_rc%udef
     
     do r=1,UAsnow_obs_struc(n)%nr
        do c=1,UAsnow_obs_struc(n)%nc
           if(UA_swe(c,r).ge.0) then 
              swe_in(c+(r-1)*UAsnow_obs_struc(n)%nc) = & 
                   UA_swe(c,r)
              lb(c+(r-1)*UAsnow_obs_struc(n)%nc) = .true. 
           endif
           if(UA_snwd(c,r).ge.0) then 
              snwd_in(c+(r-1)*UAsnow_obs_struc(n)%nc) = & 
                   UA_snwd(c,r)
              lb(c+(r-1)*UAsnow_obs_struc(n)%nc) = .true. 
           endif
        enddo
     enddo
     
     deallocate(UA_swe)
     deallocate(UA_snwd)


     call bilinear_interp(LIS_rc%gridDesc(n,:),&
          lb,swe_in,lo,swe_out,&
          UAsnow_obs_struc(n)%nc*UAsnow_obs_struc(n)%nr,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          UAsnow_obs_struc(n)%w11,&
          UAsnow_obs_struc(n)%w12,&
          UAsnow_obs_struc(n)%w21,&
          UAsnow_obs_struc(n)%w22,&
          UAsnow_obs_struc(n)%n11,&
          UAsnow_obs_struc(n)%n12,&
          UAsnow_obs_struc(n)%n21,&
          UAsnow_obs_struc(n)%n22,&
          LIS_rc%udef,iret)

     call bilinear_interp(LIS_rc%gridDesc(n,:),&
          lb,snwd_in,lo,snwd_out,&
          UAsnow_obs_struc(n)%nc*UAsnow_obs_struc(n)%nr,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          UAsnow_obs_struc(n)%w11,&
          UAsnow_obs_struc(n)%w12,&
          UAsnow_obs_struc(n)%w21,&
          UAsnow_obs_struc(n)%w22,&
          UAsnow_obs_struc(n)%n11,&
          UAsnow_obs_struc(n)%n12,&
          UAsnow_obs_struc(n)%n21,&
          UAsnow_obs_struc(n)%n22,&
          LIS_rc%udef,iret)

     deallocate(lb)
     deallocate(lo)
     deallocate(swe_in)
     deallocate(snwd_in)

  endif
  call ESMF_StateGet(Obj_Space,"UA_SNOD",snodField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: UA_SNOD')

  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status, 'Error in ESMF_FieldGet: snodField')

  call ESMF_StateGet(Obj_Space,"UA_SWE",sweField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: UA_SWE')

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status, 'Error in ESMF_FieldGet: sweField')

  snod = LIS_rc%udef
  swe = LIS_rc%udef

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then
           grid_index =LIS_domain(n)%gindex(c,r)

           if(LIS_rc%hr.eq.12.and.LIS_rc%mn.eq.0.and.LIS_rc%ss.eq.0) then
              snod(grid_index) = &
                   snwd_out(c+(r-1)*LIS_rc%lnc(n))

              swe(grid_index) = &
                   swe_out(c+(r-1)*LIS_rc%lnc(n))

           endif
        endif
     enddo
  enddo

  deallocate(snwd_out)
  deallocate(swe_out)

  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status, 'Error in ESMF_AttributeSet: Data Update Status')

end subroutine read_UAsnowobs

!BOP
!
! !ROUTINE: create_UAsnowobs_filename
! \label(create_UAsnowobs_filename)
!
! !INTERFACE:
subroutine create_UAsnowobs_filename(&
     obsdir, &
     yr, mo, nt, fname)
!
! !USES:
  implicit none
!
! !INPUT PARAMETERS:
  character(len=*), intent(in) :: obsdir
  integer  ,        intent(in) :: yr
  integer  ,        intent(in) :: mo
  integer  ,        intent(out) :: nt
  character(len=*), intent(out) :: fname
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
  character (len=4) :: fyr
  integer           :: tmpyr

  if (mo.ge.10) then
     tmpyr = yr + 1
  else
     tmpyr = yr
  endif
  
  write(fyr, '(i4.4)' ) tmpyr
  
  !leap year
  if (((mod(tmpyr,4).eq.0).and.(mod(tmpyr,100).ne.0)).or. &
       (mod(tmpyr,400).eq.0)) then
     nt = 366
  else
     nt = 365
  endif
  
  fname = trim(obsdir)//'/4km_SWE_Depth_WY'//trim(fyr)//'_v01.nc'

end subroutine create_UAsnowobs_filename


