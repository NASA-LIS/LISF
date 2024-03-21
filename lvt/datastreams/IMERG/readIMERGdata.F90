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
!NOTE:  Supports V06x and V07x
subroutine readIMERGdata(source)

   ! Imports
   use ESMF
   use IMERG_dataMod
   use LVT_coreMod
   use LVT_histDataMod
   use LVT_logMod
   use LVT_timeMgrMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif

   ! Defaults
   implicit none
   
   ! Arguments
   integer,intent(in) :: source

   ! Local variables
   integer             :: ftn 
   character*100       :: filename
   logical             :: file_exists
   real                :: prcp_in(imergdata(source)%nc,imergdata(source)%nr)
   real                :: prcp_in1(imergdata(source)%nc*imergdata(source)%nr)
   logical*1           :: lb(imergdata(source)%nc*imergdata(source)%nr)
   logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
   real                :: prcp(LVT_rc%lnc*LVT_rc%lnr)
   real                :: prcp_final(LVT_rc%lnc, LVT_rc%lnr)
   integer             :: t,c,r
   integer             :: iret, ireaderr
   real                :: currTime
   logical             :: alarmCheck
   integer                      :: yr1, mo1, da1, hr1, mn1, ss1
   integer                      :: yr2, mo2, da2, hr2, mn2, ss2
   type(ESMF_Time)              :: time1
   type(ESMF_Time)              :: time2
   type(ESMF_TimeInterval)      :: lis_ts
   integer :: status
   integer :: precipid
   integer :: start(3), count(3)
   integer :: jda2

   ! Initialize variables
   prcp(:) = LVT_rc%udef
   prcp_final(:,:) = LVT_rc%udef
   currTime = float(LVT_rc%dhr(source))*3600+ &
        60*LVT_rc%dmn(source) + LVT_rc%dss(source)
   ! Read every 30 minutes
   alarmCheck = (mod(currtime,1800.0).eq.0)
   ! Read precipitation. IMERG data are available every 30 minutes.
   yr1 = LVT_rc%dyr(source)
   mo1 = LVT_rc%dmo(source)
   da1 = LVT_rc%dda(source)
   hr1 = LVT_rc%dhr(source)
   mn1 = LVT_rc%dmn(source)
   ss1 = 0
   call ESMF_TimeSet(time1,yy=yr1, mm=mo1, dd=da1, &
       h=hr1,m=mn1,s=ss1,calendar=LVT_calendar,rc=status)
   call LVT_verify(status)
   !! If it is 00Z, use previous day's time level
   !if (mod(currtime,86400.0).eq.0) then
   !   call ESMF_TimeIntervalSet(lis_ts, s = 86400, &
   !        rc=status)
   !   call LVT_verify(status)  
   !else
   !   call ESMF_TimeIntervalSet(lis_ts, s = 0, &
   !        rc=status)
   !   call LVT_verify(status)  
   !end if
   ! Use previous IMERG file (IMERG accumulates forward in time)
   call ESMF_TimeIntervalSet(lis_ts, s = 1800, &
        rc=status)
   call LVT_verify(status)
   
   time2 = time1 - lis_ts
   call ESMF_TimeGet(time2,yy=yr2, mm=mo2, dd=da2, &
        h=hr2,m=mn2,s=ss2,calendar=LVT_calendar, &
        dayOfYear=jda2, rc=status)
   call LVT_verify(status)

   if (alarmCheck) then
      !call create_IMERG_filename(imergdata(source)%odir, &
      !                       yr1,mo1,da1,hr1,mn1,filename,imergdata(source)%imergver)
      call create_IMERG_filename(imergdata(source)%odir, &
           yr2,mo2,da2,hr2,mn2,filename,imergdata(source)%imergver, &
           imergdata(source)%imergprd)
      inquire(file=trim(filename),exist=file_exists)

      if(file_exists) then 
         write(LVT_logunit,*) '[INFO] Reading IMERG data ',trim(filename)
         call read_imerghdf(filename, imergdata(source)%nc, &
              imergdata(source)%nr, imergdata(source)%imergver, &
              prcp_in, ireaderr)
         if(ireaderr .eq. 0) then
            ! Use budget-bilinear interpolation if IMERG data are at 
            ! coarser resolution than the analysis grid; otherwise, use
            ! upscale averaging.
            prcp_in1(:) = 0
            lb(:) = .false.
            t = 1
            do r = 1, imergdata(source)%nr
               do c = 1,imergdata(source)%nc
                  prcp_in1(t) = prcp_in(c,r)
                  if (prcp_in1(t) .ge. 0) then
                     lb(t) = .true.
                  end if
                  t = t + 1
               end do ! c
            end do ! r

            if (LVT_isAtAFinerResolution(imergdata(source)%datares)) then
               call conserv_interp(LVT_rc%gridDesc,lb,prcp_in1, &
                     lo,prcp, imergdata(source)%nc*imergdata(source)%nr, &
                     LVT_rc%lnc*LVT_rc%lnr, imergdata(source)%rlat, &
                     imergdata(source)%rlon, &
                     imergdata(source)%w112, imergdata(source)%w122, &
                     imergdata(source)%w212, imergdata(source)%w222, &
                     imergdata(source)%n112, imergdata(source)%n122, &
                     imergdata(source)%n212, imergdata(source)%n222, &
                     LVT_rc%udef, iret)   
            else
               call upscaleByAveraging(&
                     imergdata(source)%nc*imergdata(source)%nr, &
                     LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
                     imergdata(source)%n11, lb, &
                     prcp_in1, lo, prcp)
            endif
            write(LVT_logunit,*) '[INFO] Finished processing ',trim(filename)               
         else
            write(LVT_logunit,*)'[ERR] Read error with IMERG file ', &
                                 trim(filename)
            prcp = LVT_rc%udef
         endif
      else 
         write(LVT_logunit,*)'[ERR] Missing IMERG file ', trim(filename)
         prcp = LVT_rc%udef
      end if ! file_exists

      do r=1,LVT_rc%lnr
         do c=1, LVT_rc%lnc
            prcp_final(c,r) = prcp(c+(r-1)*LVT_rc%lnc)
         end do ! c
      end do ! r
   end if ! alarmCheck
   ! Convert mm/hr to kg/m2s.
   do r=1,LVT_rc%lnr
      do c=1, LVT_rc%lnc
         if (prcp_final(c,r) .ge. 0) then            
            prcp_final(c,r) = prcp_final(c,r)/3600.
         else
            prcp_final(c,r) = LVT_rc%udef
         end if
      end do ! c
   end do ! r
   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final,&
        vlevel=1,units='kg/m2s')
   ! Now convert from kg/m2s to kg/m2
   do r=1,LVT_rc%lnr
      do c=1,LVT_rc%lnc
         if(prcp_final(c,r).ge.0) then
            !prcp_final(c,r) = prcp_final(c,r)*86400.0 !kg/m2
            prcp_final(c,r) = prcp_final(c,r)*1800. ! kg/m2 for 30 minutes
         else
            prcp_final(c,r) = LVT_rc%udef
         endif
      enddo ! c
   enddo ! r
   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final,&
        vlevel=1,units='kg/m2') 
end subroutine readIMERGdata
!------------------------------------------------------------------------------
subroutine read_imerghdf(filename, col, row, version, precipout, ireaderr)
#if(defined USE_HDF5)
      use HDF5
#endif
   use ESMF
   use IMERG_dataMod
   use LVT_coreMod
   use LVT_histDataMod
   use LVT_logMod
   use LVT_timeMgrMod

   implicit none

   character(len=*), intent(in) :: filename
   integer, intent(in)          :: col, row
   character*10, intent(in)     :: version
   integer, intent(out)         :: ireaderr
   real,    intent(out)         :: precipout(col,row)
   !Local variables
   integer :: xsize, ysize
   !character(len=40) :: dsetname='/Grid/precipitationCal'
   character(len=40) :: dsetname
   character(len=40) :: dsetname6='/Grid/precipitationCal' ! V06x
   character(len=40) :: dsetname7='/Grid/precipitation'    ! V07x
   real :: precipin(row,col)
   logical :: bIsError
   integer :: istatus,i
#if(defined USE_HDF5)
     integer(HSIZE_T), dimension(2) :: dims
     integer(HID_T) :: fileid,dsetid
#endif

#if(defined USE_HDF5)
      xsize = col
      ysize = row
      dims(1) = xsize
      dims(2) = ysize
      bIsError=.false.
!open fortran interface
      call h5open_f(istatus)
      if(istatus.ne.0) then
         bIsError=.true.
         write(LVT_logunit,*) 'Error opening HDF5 fortran interface'
         ireaderr = istatus
         return
      endif
!open hdf5 file
      call h5fopen_f(filename,H5F_ACC_RDONLY_F,fileid,istatus)
      if(istatus.ne.0) then
         bIsError=.true.
         write(LVT_logunit,*) 'Error opening IMERG file',trim(filename)
         ireaderr = istatus
         return
      endif
!open dataset
      if (index(trim(version), "V06") .ne. 0) then
         dsetname = dsetname6
      else if (index(trim(version), "V07") .ne. 0) then
         dsetname = dsetname7
      else
         write(LVT_logunit,*)'[ERR] Invalid IMERG version number!'
         write(LVT_logunit,*)'[ERR] Expected V06x or V07x!'
         write(LVT_logunit,*)'[ERR] Received ', trim(version)
         stop
      end if
      call h5dopen_f(fileid,dsetname,dsetid,istatus)
      if(istatus.ne.0) then
         bIsError=.true.
         write(LVT_logunit,*) 'Error opening IMERG dataset',trim(dsetname)
         ireaderr = istatus
         return
      endif
!read dataset 
      call h5dread_f(dsetid,H5T_NATIVE_REAL,precipin,dims,istatus)
      if(istatus.ne.0) then
         bIsError=.true.
         write(LVT_logunit,*) 'Error reading IMERG dataset',trim(dsetname)
         ireaderr = istatus
         return
      endif
!Put the real(1:,1:) on the precipout(0:,0:)
!precipin is (ysize,xsize) starting at (lon=-179.9,lat=-89.9)
      precipout(1:xsize,1:ysize)=transpose(precipin)
!close dataset
      call h5dclose_f(dsetid,istatus)
      if(istatus.ne.0) then
         bIsError=.true.
         write(LVT_logunit,*) 'Error closing IMERG dataset',trim(dsetname)
         ireaderr = istatus
         return
      endif
!close file
      call h5fclose_f(fileid,istatus)
      if(istatus.ne.0) then
         bIsError=.true.
         write(LVT_logunit,*) 'Error closing IMERG file',trim(filename)
         ireaderr = istatus
         return
      endif
!close fortran interface
      call h5close_f(istatus)
      if(istatus.ne.0) then 
         bIsError=.true.
         write(LVT_logunit,*) 'Error closing HDF5 fortran interface'
         ireaderr = istatus
         return
      endif
      ireaderr = istatus
#endif
end subroutine read_imerghdf
!------------------------------------------------------------------------------
subroutine create_IMERG_filename(odir, &
                             yr,mo,da,hr,mn,filename,imVer, imergprd)
   use IMERG_dataMod
   use LVT_logMod

   ! Defaults
   implicit none

   ! Arguments
   character(len=*), intent(in) :: odir, imVer, imergprd
   integer, intent(in) :: yr, mo, da, hr, mn
   character(len=*), intent(out) :: filename

   ! Local variables
   integer :: uyr, umo, uda, uhr, umn, umnadd, umnday, uss
   character*4   :: cyr, cmnday
   character*2   :: cmo, cda, chr, cmn, cmnadd
   character*100 :: fbase, ftimedir, fstem, fext
   
   uyr = yr
   umo = mo
   uda = da
   uhr = 1*(hr/1)  !hour needs to be a multiple of 1 hour
   if (mn .ge. 0 .AND. mn .lt. 30) then
     umn = 0
   else
     umn = 30
   endif
   umnadd = umn + 29
   umnday = uhr*60 + umn
   uss = 0
   write(cyr, '(I4.4)') uyr
   write(cmo, '(I2.2)') umo
   write(cda, '(I2.2)') uda
   write(chr, '(I2.2)') uhr
   write(cmn, '(I2.2)') umn 
   write(cmnadd, '(I2.2)') umnadd
   write(cmnday, '(I4.4)')umnday

   if(imergprd == 'early') then
      fstem = '/3B-HHR-E.MS.MRG.3IMERG.'
      fext = '.RT-H5'
   elseif(imergprd == 'late') then
      fstem = '/3B-HHR-L.MS.MRG.3IMERG.'
      fext = '.RT-H5'
   elseif(imergprd == 'final') then
      fstem = '/3B-HHR.MS.MRG.3IMERG.'
      fext = '.HDF5'
   else
      write(LVT_logunit,*) "[ERR] Invalid IMERG product option was chosen."
      write(LVT_logunit,*) "[ERR] Please choose either 'early', 'late', or 'final'."
      call LVT_endrun()
   endif
   
   filename = trim(odir)//"/"//cyr//cmo//trim(fstem)// &
         cyr//cmo//cda//"-S"//chr//cmn//"00-E"//chr//cmnadd//"59."//cmnday//"."//trim(imVer)//fext
end subroutine create_IMERG_filename
