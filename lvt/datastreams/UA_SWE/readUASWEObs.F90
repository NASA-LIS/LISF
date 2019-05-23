!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readUASWEObs
! \label{readUASWEObs}
!
! !INTERFACE: 
subroutine readUASWEObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,    only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod,     only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_verify, LVT_endrun
  use LVT_timeMgrMod, only : LVT_calendar, LVT_tick
  use UASWE_obsMod,    only : uasweobs
  use map_utils

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer,  intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads and processes the gridded UA SWE
!  data. The data for a given year is read into memory at the 
!  start of a year and indexed into during each day.
! 
! !FILES USED:
! 
!EOP

  integer             :: ftn 
  character*100       :: uafilename
  character*4         :: fyr
  logical             :: file_exists
  real, allocatable   :: swe1(:,:,:),swe2(:,:,:)
  real                :: swe_in(uasweobs(source)%nc*uasweobs(source)%nr)
  logical*1           :: lb(uasweobs(source)%nc*uasweobs(source)%nr)
  logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                :: swe_out(LVT_rc%lnc*LVT_rc%lnr)
  real                :: swe_final(LVT_rc%lnc, LVT_rc%lnr)
  integer             :: k,c,r,nt,jj
  integer             :: nid,varid
  integer             :: iret
  type(ESMF_Time)     :: uatime1
  integer             :: status
  real                :: timenow
  logical             :: alarmCheck

  swe_out = LVT_rc%udef
  swe_final = LVT_rc%udef

if((uasweobs(source)%yr.ne.LVT_rc%dyr(source)).or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     call ESMF_TimeSet(uasweobs(source)%startTime,&
          yy=LVT_rc%dyr(source), &
          mm=1, & 
          dd = 1, &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, 'error in setting scan start time')

     uasweobs(source)%yr = LVT_rc%dyr(source)
     uasweobs(source)%swe = LVT_rc%udef

  do jj=1,2 !once to read the current year and one to read
            !next year -This is done because the UA data is provided
            !in water years. 
     if(jj.eq.1) then 
        call create_UASWE_filename(uasweobs(source)%odir, &
             LVT_rc%dyr(source), uafilename)

             if((mod(LVT_rc%dyr(source),4) .eq. 0 .and. &
                 mod(LVT_rc%dyr(source), 100).ne.0) &!leap year
                 .or.(mod(LVT_rc%dyr(source),400) .eq.0)) then
                 nt = 366
             else
                 nt = 365
             endif

             allocate(swe1(uasweobs(source)%nc,&
                      uasweobs(source)%nr,&    
                      nt))

             inquire(file=trim(uafilename), exist=file_exists)

             if(file_exists) then
               write(LVT_logunit,*) '[INFO] Reading UA data ',&
               trim(uafilename)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
             iret = nf90_open(path=trim(uafilename),mode=NF90_NOWRITE,ncid=nid)
             call LVT_verify(iret, 'Error opening file'//trim(uafilename))

             iret = nf90_inq_varid(nid, 'SWE',varid)
             call LVT_verify(iret, 'Error nf90_inq_varid: SWE')

             iret = nf90_get_var(nid,varid, SWE1)
             call LVT_verify(iret, 'Error nf90_get_var: SWE')

             iret = nf90_close(nid)
             call LVT_verify(iret, 'Error nf90_close')
#endif
             uasweobs(source)%swe(:,:,1:nt-92) = swe1(:,:,93:nt) !Jan.-Sep.
             endif
           
             write(LVT_logunit,*) '[INFO] Finished processing ',trim(uafilename)

     elseif(jj.eq.2) then
        call create_UASWE_filename(uasweobs(source)%odir, &
             LVT_rc%dyr(source)+1, uafilename)

            if((mod(LVT_rc%dyr(source)+1,4) .eq. 0 .and. &
                mod(LVT_rc%dyr(source)+1, 100).ne.0) &!leap year
                .or.(mod(LVT_rc%dyr(source)+1,400) .eq.0)) then
                nt = 366
            else
                nt = 365
            endif

            allocate(swe2(uasweobs(source)%nc,&
                     uasweobs(source)%nr,&    
                     nt))

            inquire(file=trim(uafilename), exist=file_exists)

            if(file_exists) then
              write(LVT_logunit,*) '[INFO] Reading UA data ',&
              trim(uafilename)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
            iret = nf90_open(path=trim(uafilename),mode=NF90_NOWRITE,ncid=nid)
            call LVT_verify(iret, 'Error opening file'//trim(uafilename))

            iret = nf90_inq_varid(nid, 'SWE',varid)
            call LVT_verify(iret, 'Error nf90_inq_varid: SWE')

            iret = nf90_get_var(nid,varid, SWE2)
            call LVT_verify(iret, 'Error nf90_get_var: SWE')

            iret = nf90_close(nid)
            call LVT_verify(iret, 'Error nf90_close')
#endif
            uasweobs(source)%swe(:,:,nt-92+1:nt) = swe2(:,:,1:92) !Oct.-Dec.
            endif
           
            write(LVT_logunit,*) '[INFO] Finished processing ',trim(uafilename)
     endif
  enddo
  deallocate(swe1,swe2)
endif 

  timenow = float(LVT_rc%dhr(source))*3600 + &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  if(alarmCheck) then 
     call ESMF_TimeSet(uatime1, yy=LVT_rc%dyr(source), &
          mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
          h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
          s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'uatime1 set failed')
     
     k = nint((uatime1 - uasweobs(source)%starttime)/&
          uasweobs(source)%timestep)+1
     
     lb = .false. 
     swe_in = LVT_rc%udef
     do r=1,uasweobs(source)%nr
        do c=1,uasweobs(source)%nc
           if(uasweobs(source)%swe(c,r,k).ge.0) then
              swe_in(c+(r-1)*uasweobs(source)%nc) = &
                   uasweobs(source)%swe(c,r,k)
              lb(c+(r-1)*uasweobs(source)%nc) = .true. 
           endif
        enddo
     enddo
     
     call bilinear_interp(LVT_rc%gridDesc,lb,swe_in, &
          lo,swe_out, &
          uasweobs(source)%nc*uasweobs(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr,   &
          uasweobs(source)%rlat, &
          uasweobs(source)%rlon, &
          uasweobs(source)%w11,  &
          uasweobs(source)%w12,  &
          uasweobs(source)%w21,  &
          uasweobs(source)%w22,  &
          uasweobs(source)%n11,  &
          uasweobs(source)%n12,  &
          uasweobs(source)%n21,  &
          uasweobs(source)%n22,  &
          LVT_rc%udef, iret)     
     
     do r=1,LVT_rc%lnr
        do c=1, LVT_rc%lnc
           swe_final(c,r) = swe_out(c+(r-1)*LVT_rc%lnc)
        enddo
     enddo
    
     ! Convert mm to kg/m2s (note that 1 mm = 1 kg/m2)
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(swe_final(c,r).ge.0) then
              swe_final(c,r) = swe_final(c,r)/86400.0 !kg/m2s
           else
              swe_final(c,r) = LVT_rc%udef
           endif
        enddo
     enddo
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,swe_final,vlevel=1,units="kg/m2s")

  ! Now convert from kg/m2s to kg/m2
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(swe_final(c,r).ge.0) then
           swe_final(c,r) = swe_final(c,r)*86400.0 !kg/m2
        else
           swe_final(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,swe_final,vlevel=1,units="kg/m2")

end subroutine readUASWEObs

!BOP
! 
! !ROUTINE: create_UASWE_filename
! \label(create_UASWE_filename)
!
! !INTERFACE:
subroutine create_UASWE_filename(odir, yr,uaname)
! 
! !USES:   
  use LVT_String_Utility
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: uaname
!

  character*4             :: fyr

  write(fyr, '(i4.4)' ) yr

  uaname = trim(odir)//'/4km_SWE_Depth_WY'//trim(fyr)//'_v01.nc'

end subroutine create_UASWE_filename
