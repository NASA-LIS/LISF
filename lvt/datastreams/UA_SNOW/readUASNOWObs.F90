!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
!
! !ROUTINE: readUASNOWObs
! \label{readUASNOWObs}
!
! !INTERFACE:
subroutine readUASNOWObs(source)
!
! !USES:
   use ESMF
   use LVT_histDataMod
   use LVT_coreMod,    only : LVT_rc
   use LVT_timeMgrMod, only : LVT_calendar
   use LVT_logMod,     only : LVT_logunit, LVT_verify
   use UASNOW_obsMod,  only : uasnowobs
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
! This subroutine reads and processes the gridded UA SNOW
! data.  The data for a given year is read into memory at
! the start of a year and indexed into during each day.
!
! !REVISION HISTORY:
! 28 May 2019: Rhae Sung Kim, Initial Specification
! 19 Jun 2019: David Mocko, Set valid time of data to 12Z
! 25 Jun 2019: David Mocko, Only read one file per call
!
!EOP

   integer             :: ftn
   character*100       :: uafilename
   logical             :: file_exists
   real, allocatable   :: swe1(:,:,:),snwd1(:,:,:)
   real                :: swe_in(uasnowobs(source)%nc*uasnowobs(source)%nr)
   real                :: snwd_in(uasnowobs(source)%nc*uasnowobs(source)%nr)
   logical*1           :: lb(uasnowobs(source)%nc*uasnowobs(source)%nr)
   logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
   real                :: swe_out(LVT_rc%lnc*LVT_rc%lnr)
   real                :: snwd_out(LVT_rc%lnc*LVT_rc%lnr)
   real                :: swe_final(LVT_rc%lnc, LVT_rc%lnr)
   real                :: snwd_final(LVT_rc%lnc, LVT_rc%lnr)
   integer             :: c,r,k,nt
   integer             :: nid,varid
   integer             :: iret
   type(ESMF_Time)     :: uatime1
   integer             :: status
   real                :: timenow
   logical             :: alarmCheck

   swe_out = LVT_rc%udef
   swe_final = LVT_rc%udef
   snwd_out = LVT_rc%udef
   snwd_final = LVT_rc%udef

   timenow = float(LVT_rc%dhr(source)*3600 + &
                   LVT_rc%dmn(source)*60 + LVT_rc%dss(source))

! The UA dataset is a daily product, valid approximately
! at 12Z, which is the observation time of many of the
! surface products used to generate the data (personal
! communication with the UA SNOW product generators).
   alarmcheck = (mod(timenow, 86400.0).eq.43200.0)

   if (alarmcheck) then
      LVT_rc%resetFlag(source) = .false.
      call ESMF_TimeSet(uasnowobs(source)%startTime, &
           yy=LVT_rc%dyr(source), &
           mm = 10, &
           dd = 1, &
           h = 0,  &
           m = 0,  &
           calendar = LVT_calendar, &
           rc=status)
      call LVT_verify(status, 'error in setting scan start time')

      uasnowobs(source)%yr = LVT_rc%dyr(source)
      uasnowobs(source)%swe = LVT_rc%udef
      uasnowobs(source)%snwd = LVT_rc%udef

      call create_UASNOW_filename(uasnowobs(source)%odir, &
           LVT_rc%dyr(source), LVT_rc%dmo(source), nt, uafilename)

      inquire(file=trim(uafilename),exist=file_exists)

      if (file_exists) then
         write(LVT_logunit,*) '[INFO] Reading UA data ',&
                                trim(uafilename)

         allocate(swe1(uasnowobs(source)%nc,   &
                       uasnowobs(source)%nr,nt))

         allocate(snwd1(uasnowobs(source)%nc,   &
                        uasnowobs(source)%nr,nt))

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
         iret = nf90_open(path=trim(uafilename),mode=NF90_NOWRITE,ncid=nid)
         call LVT_verify(iret, 'Error opening file'//trim(uafilename))

         iret = nf90_inq_varid(nid, 'SWE',varid)
         call LVT_verify(iret, 'Error nf90_inq_varid: SWE')

         iret = nf90_get_var(nid,varid, SWE1)
         call LVT_verify(iret, 'Error nf90_get_var: SWE')

         iret = nf90_close(nid)
         call LVT_verify(iret, 'Error nf90_close')

         iret = nf90_open(path=trim(uafilename),mode=NF90_NOWRITE,ncid=nid)
         call LVT_verify(iret, 'Error opening file'//trim(uafilename))

         iret = nf90_inq_varid(nid, 'DEPTH',varid)
         call LVT_verify(iret, 'Error nf90_inq_varid: SNWD')

         iret = nf90_get_var(nid,varid, SNWD1)
         call LVT_verify(iret, 'Error nf90_get_var: SNWD')

         iret = nf90_close(nid)
         call LVT_verify(iret, 'Error nf90_close')
#endif

         call ESMF_TimeSet(uatime1, yy=LVT_rc%dyr(source), &
              mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
              h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
              s=LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
         call LVT_verify(status, 'uatime1 set failed')

         k = nint((uatime1 - uasnowobs(source)%starttime) / &
                             uasnowobs(source)%timestep)
         if (k.lt.0) k = k + nt + 1

         uasnowobs(source)%swe(:,:)  = swe1(:,:,k) !Jan. to Sep.
         uasnowobs(source)%snwd(:,:) = snwd1(:,:,k) !Jan. to Sep.
      else
         write(LVT_logunit,*) '[WARN] UA file not found: ',&
                                trim(uafilename)
      endif

      deallocate(swe1)
      deallocate(snwd1)

      lb = .false.
      swe_in = LVT_rc%udef
      snwd_in = LVT_rc%udef

      do r=1,uasnowobs(source)%nr
         do c=1,uasnowobs(source)%nc
            if (uasnowobs(source)%swe(c,r).ge.0) then
               swe_in(c+(r-1)*uasnowobs(source)%nc) = &
                              uasnowobs(source)%swe(c,r)
               lb(c+(r-1)*uasnowobs(source)%nc) = .true.
            endif

            if (uasnowobs(source)%snwd(c,r).ge.0) then
               snwd_in(c+(r-1)*uasnowobs(source)%nc) = &
                               uasnowobs(source)%snwd(c,r)
               lb(c+(r-1)*uasnowobs(source)%nc) = .true.
            endif
         enddo
      enddo

      call bilinear_interp(LVT_rc%gridDesc,lb,swe_in, &
           lo,swe_out, &
           uasnowobs(source)%nc*uasnowobs(source)%nr, &
           LVT_rc%lnc*LVT_rc%lnr,  &
           uasnowobs(source)%rlat, &
           uasnowobs(source)%rlon, &
           uasnowobs(source)%w11,  &
           uasnowobs(source)%w12,  &
           uasnowobs(source)%w21,  &
           uasnowobs(source)%w22,  &
           uasnowobs(source)%n11,  &
           uasnowobs(source)%n12,  &
           uasnowobs(source)%n21,  &
           uasnowobs(source)%n22,  &
           LVT_rc%udef, iret)

      call bilinear_interp(LVT_rc%gridDesc,lb,snwd_in, &
           lo,snwd_out, &
           uasnowobs(source)%nc*uasnowobs(source)%nr, &
           LVT_rc%lnc*LVT_rc%lnr,  &
           uasnowobs(source)%rlat, &
           uasnowobs(source)%rlon, &
           uasnowobs(source)%w11,  &
           uasnowobs(source)%w12,  &
           uasnowobs(source)%w21,  &
           uasnowobs(source)%w22,  &
           uasnowobs(source)%n11,  &
           uasnowobs(source)%n12,  &
           uasnowobs(source)%n21,  &
           uasnowobs(source)%n22,  &
           LVT_rc%udef, iret)

      do r=1,LVT_rc%lnr
         do c=1, LVT_rc%lnc
            swe_final(c,r) = swe_out(c+(r-1)*LVT_rc%lnc)
            snwd_final(c,r) = snwd_out(c+(r-1)*LVT_rc%lnc)
         enddo
      enddo

      do r=1,LVT_rc%lnr
         do c=1,LVT_rc%lnc
            if(swe_final(c,r).ge.0) then
               swe_final(c,r) = swe_final(c,r)          ! (Note that 1 mm = 1 kg/m2)
            else
               swe_final(c,r) = LVT_rc%udef
            endif
            if(snwd_final(c,r).ge.0) then
               snwd_final(c,r) = snwd_final(c,r)/1000.0 ! Convert mm to m
            else
               snwd_final(c,r) = LVT_rc%udef
            endif
         enddo
      enddo
   endif

   call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,swe_final,vlevel=1,units="kg/m2")
   call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH,source,snwd_final,vlevel=1,units="m")

end subroutine readUASNOWObs

!BOP
!
! !ROUTINE: create_UASNOW_filename
! \label(create_UASNOW_filename)
!
! !INTERFACE:
subroutine create_UASNOW_filename(odir,yr,mo,nt,uaname)
!
! !USES:
   use LVT_String_Utility
   implicit none
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)  :: odir
   integer,          intent(in)  :: yr,mo
!
! !OUTPUT PARAMETERS:
   integer, intent(out)          :: nt
   character(len=*), intent(out) :: uaname
!
   character*4                   :: fyr
   integer                       :: tmpyr

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

   uaname = trim(odir)//'/4km_SWE_Depth_WY'//trim(fyr)//'_v01.nc'

end subroutine create_UASNOW_filename
