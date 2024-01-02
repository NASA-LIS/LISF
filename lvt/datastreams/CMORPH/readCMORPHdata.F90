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

!------------------------------------------------------------------------------

subroutine readCMORPHdata(source)

   ! Imports
   use ESMF
   use CMORPH_dataMod
   use LVT_coreMod
   use LVT_histDataMod
   use LVT_logMod
   use LVT_timeMgrMod

   ! Defaults
   implicit none
   
   ! Arguments
   integer,intent(in) :: source

   ! Local variables
   integer,   parameter         :: nc = 1440, nr=480
   character*200                :: filename
   logical                      :: file_exists
   integer                      :: nunit,ufn,iret,ierr
   integer                      :: c,r
   integer                      :: ftn
   real                         :: precip_f(nc*nr)
   integer                      :: nprecip_f(nc*nr)
   real                         :: varfield(LVT_rc%lnc,LVT_rc%lnr)
   integer                      :: yr1, mo1, da1, hr1, mn1, ss1
   integer                      :: yr2, mo2, da2, hr2, mn2, ss2
   type(ESMF_Time)              :: time1
   type(ESMF_Time)              :: time2
   type(ESMF_TimeInterval)      :: lis_ts
   real                         :: var(nc*nr)
   real                         :: tmp(nc,nr)
   integer :: ihr
   integer :: irec
   integer :: c_tmp
   integer :: status,ios

   ! Initialize variables
   precip_f(:)  = 0
   nprecip_f(:) = 0
   
   ! Set the time information.  Note that CMORPH dates indicate the
   ! starting time of the precipitation, so we need look back 3 hours
   ! to get the "present" precipitation.
   yr1 = LVT_rc%dyr(source)
   mo1 = LVT_rc%dmo(source)
   da1 = LVT_rc%dda(source)
   hr1 = LVT_rc%dhr(source)
   mn1 = 0
   ss1 = 0
   call ESMF_TimeSet(time1,yy=yr1, mm=mo1, dd=da1, &
       h=hr1,m=mn1,s=ss1,calendar=LVT_calendar, rc=status)
   call LVT_verify(status)

   call ESMF_TimeIntervalSet(lis_ts, s = LVT_rc%ts, &
        rc=status)
   call LVT_verify(status)  
   time2 = time1 - lis_ts

   call ESMF_TimeGet(time2,yy=yr2, mm=mo2, dd=da2, &
        h=hr2,m=mn2,s=ss2,calendar=LVT_calendar, rc=status)
   call LVT_verify(status)
   
   call create_cmorphdata_filename(source, &
        yr2, mo2, da2, filename)

   inquire(file=trim(filename),exist=file_exists)

   if (.not. file_exists) then
      write(LVT_logunit,*) '[INFO] Cannot find file : ',trim(filename)
   end if

   if(file_exists) then
      ftn = LVT_getNextUnitNumber()
         
      open(unit=ftn, file=trim(filename), access='direct', &
           recl=nc*nr*4, form='unformatted', status='old', &
           convert="little_endian", iostat=ios)
      if (ios .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Stopping routine: File not opened : ',ios
         call LVT_endrun()
      end if
      
      write(LVT_logunit,*) '[INFO] Reading CMORPH file: ',trim(filename)
         
      ! Loop through the file and find the desired time slice.
      irec = 0
      do ihr = 0, 21, 3
         irec = irec + 1
         
         ! Skip if this is the wrong hour.
         if (ihr .ne. hr2) cycle

         read(unit=ftn, rec=irec) var
         
         ! CMORPH grid is from 0 to 360 deg E, but the interpolation code
         ! wants -180 to 180 E.  So we need to reorder the data accordingly.
         do r = 1,nr
            do c = 1, nc
               if ( c .le. (nc/2) ) then
                  c_tmp = (nc/2) + c
               else
                  c_tmp = c - (nc/2)
               end if
               tmp(c_tmp,r) = var(c+(r-1)*nc)
            end do
         end do
         do r = 1,nr
            do c = 1, nc
               var(c+(r-1)*nc) = tmp(c,r)
            end do
         end do
         
         ! The CMORPH GrADS control file lists the units as mm/3hr, but
         ! comparisons with the daily CPC Unified product suggests this is
         ! really a 3-hr accumulation, not a rate.
         ! Copy the valid values, converting from mm to
         ! kg/m2s (1 mm water = 1 kg/m2 of water).
         do r = 1,nr
            do c = 1,nc
               if (var(c+(r-1)*nc) .ge. 0) then
                  precip_f(c+(r-1)*nc) = &
                       var(c+(r-1)*nc) / 10800. ! 3 hr = 10800 sec
                  nprecip_f(c+(r-1)*nc) = 1
               end if
            end do
         end do

         exit ! No need to continue reading the file at this point.
      end do ! ihr

      call LVT_releaseUnitNumber(ftn)
         
      write(LVT_logunit,*)'[INFO] Successfully processed ',trim(filename)
   end if

   ! Now interpolate the data.  First as rate and then as depth.
   call interp_cmorphvar(source, nc, nr, precip_f, nprecip_f, varfield)
   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip, source, varfield, &
        vlevel=1, units="kg/m2s")

   do r = 1,nr
      do c = 1,nc
         if (precip_f(c+(r-1)*nc) .ge. 0) then
            precip_f(c+(r-1)*nc) = precip_f(c+(r-1)*nc) * 10800. ! 3 hr period
         end if
      end do
   end do
   
   call interp_cmorphvar(source, nc, nr, precip_f, nprecip_f, varfield)
   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip, source, varfield, &
        vlevel=1, units="kg/m2")
   
end subroutine readCMORPHdata

!------------------------------------------------------------------------------

subroutine interp_cmorphvar(source, nc, nr, var_input, nvar_input, var_output)

   ! Imports
   use LVT_coreMod,   only : LVT_rc
   use CMORPH_dataMod, only : cmorphdata

   ! Defaults
   implicit none

   ! Arguments
   integer,intent(in)          :: source
   integer,intent(in)          :: nc
   integer,intent(in)          :: nr
   real,intent(inout)          :: var_input(nc*nr)
   integer,intent(in)          :: nvar_input(nc*nr)
   real,intent(out)            :: var_output(LVT_rc%lnc, LVT_rc%lnr)

   ! Local variables
   logical*1          :: lb(nc*nr)
   integer            :: iret
   integer            :: c,r
   logical*1          :: lo(LVT_rc%lnc*LVT_rc%lnr)
   real               :: go(LVT_rc%lnc*LVT_rc%lnr)

   ! Initialize arrays
   var_output(:,:) = LVT_rc%udef
   lb(:) = .false.

   ! Average values and construct logical bit map
   do r=1,nr
     do c=1,nc
        if(nvar_input(c+(r-1)*nc).gt.0) then 
           var_input(c+(r-1)*nc) = &
                var_input(c+(r-1)*nc) / nvar_input(c+(r-1)*nc)
           lb(c+(r-1)*nc) = .true.
        else
           var_input(c+(r-1)*nc) = LVT_rc%udef
        endif
     enddo
  enddo
  
  ! Now interpolate to the LVT grid
  call bilinear_interp(LVT_rc%gridDesc,lb, var_input,&
       lo,go,nc*nr,LVT_rc%lnc*LVT_rc%lnr,&
       cmorphdata(source)%rlat,cmorphdata(source)%rlon,&
       cmorphdata(source)%w11,cmorphdata(source)%w12,&
       cmorphdata(source)%w21,cmorphdata(source)%w22,&
       cmorphdata(source)%n11,cmorphdata(source)%n12,&
       cmorphdata(source)%n21,cmorphdata(source)%n22,LVT_rc%udef,iret)

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        var_output(c,r) = go(c+(r-1)*LVT_rc%lnc)
     enddo
  enddo

end subroutine interp_cmorphvar

!------------------------------------------------------------------------------

subroutine create_cmorphdata_filename(source, yr, mo, da, filename)

   ! Import
   use CMORPH_dataMod

   ! Defaults
   implicit none
   
   ! Arguments
   integer,intent(in) :: source
   integer,intent(in) :: yr
   integer,intent(in) :: mo
   integer,intent(in) :: da
   character(len=*),intent(inout) :: filename

   ! Local variables
   character*4       :: fyr
   character*2       :: fmo
   character*2       :: fda

   write(unit=fyr, fmt='(i4.4)') yr
   write(unit=fmo, fmt='(i2.2)') mo
   write(unit=fda, fmt='(i2.2)') da

!   filename = &
!        trim(cmorphdata(source)%odir) &
!        //'/'//fyr//'/'//fyr//fmo &
!        //'/CMORPH_V1.0_RAW_0.25deg-3HLY_'//fyr//fmo//fda
   filename = &
        trim(cmorphdata(source)%odir) &
        //'/'//fyr//'/'//fyr//fmo &
        //'/CMORPH_V1.0_ADJ_0.25deg-3HLY_'//fyr//fmo//fda

end subroutine create_cmorphdata_filename
