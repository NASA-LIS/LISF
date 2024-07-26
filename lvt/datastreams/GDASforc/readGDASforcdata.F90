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
! !MODULE: GIMMSAVHRR_NDVIobsMod
! \label(GIMMSAVHRR_NDVIobsMod)
!
! !INTERFACE:
!------------------------------------------------------------------------------

!NOTE:  Currently only 0.05 deg daily GDASforc data are supported.
subroutine readGDASforcdata(source)

   ! Imports
   use ESMF
   use GDASforc_dataMod
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
   character*100       :: name00, name03, name06
   logical             :: F06flag
   logical             :: file_exists1,file_exists2,file_exists3
   logical*1           :: lb(gdasforcdata(source)%ncold*gdasforcdata(source)%nrold)
   logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
   real                :: prcp(LVT_rc%lnc*LVT_rc%lnr)
   real                :: tmp_final(LVT_rc%lnc,LVT_rc%lnr)
   real                :: prcp_final(2,LVT_rc%lnc, LVT_rc%lnr)
   integer             :: t,c,r
   integer             :: iret,fret
   real                :: currTime
   logical             :: alarmCheck
   real*8                       :: timenow, time1_tw, time2_tw
   integer                      :: yr1, mo1, da1, hr1, mn1, ss1
   integer                      :: yr2, mo2, da2, hr2, mn2, ss2
   integer                      :: doy1,doy2,doynow
   real                         :: gmt1, gmt2, gmtnow
   type(ESMF_Time)              :: time1
   type(ESMF_Time)              :: time2
   type(ESMF_TimeInterval)      :: lis_ts
   integer :: status
   integer :: precipid
   integer :: start(3), count(3)
   integer :: jda2
   integer :: ferror2, ferror3
   real    :: gridDesci(50)
   real    :: glbdata_a(3,LVT_rc%lnc,LVT_rc%lnr)
   real    :: glbdata_a_f06(3,LVT_rc%lnc,LVT_rc%lnr)

   ! Initialize variables
   glbdata_a     = LVT_rc%udef
   glbdata_a_f06 = LVT_rc%udef
   prcp(:) = LVT_rc%udef
   tmp_final  = LVT_rc%udef
   prcp_final = LVT_rc%udef

   currTime = float(LVT_rc%dhr(source))*3600+ &
        60*LVT_rc%dmn(source) + LVT_rc%dss(source)
   ! EMK...Only read at 00Z
   alarmCheck = (mod(currtime,10800.0).eq.0)

   yr1 = LVT_rc%dyr(source)
   mo1 = LVT_rc%dmo(source)
   da1 = LVT_rc%dda(source)
   hr1 = LVT_rc%dhr(source)
   mn1 = 0
   ss1 = 0
   call ESMF_TimeSet(time1,yy=yr1, mm=mo1, dd=da1, &
       h=hr1,m=mn1,s=ss1,calendar=LVT_calendar,rc=status)
   call LVT_verify(status)

   ! If it is 00Z, use previous day's time level
   if (mod(currtime,10800.0).eq.0) then
      call ESMF_TimeIntervalSet(lis_ts, s = 10800, &
           rc=status)
      call LVT_verify(status)  
   else
      call ESMF_TimeIntervalSet(lis_ts, s = 0, &
           rc=status)
      call LVT_verify(status)  
   end if
   time2 = time1 - lis_ts

   call ESMF_TimeGet(time2,yy=yr2, mm=mo2, dd=da2, &
        h=hr2,m=mn2,s=ss2,calendar=LVT_calendar, &
        dayOfYear=jda2, rc=status)
   call LVT_verify(status)

   
   call LVT_tick(timenow, doynow, gmtnow, &
        LVT_rc%dyr(source), &
        LVT_rc%dmo(source), &
        LVT_rc%dda(source), &
        LVT_rc%dhr(source), &
        LVT_rc%dmn(source), &
        LVT_rc%dss(source), &
        0)

   call LVT_tick(time1_tw, doy1, gmt1, &
        yr1,mo1,da1,hr1,mn1,ss1,0)

   call LVT_tick(time2_tw, doy2, gmt2, &
        yr2,mo2,da2,hr2,mn2,ss2,0)
   
!= T126:
   if ( time1_tw >= gdasforcdata(source)%griduptime1 .and. &
        time1_tw < gdasforcdata(source)%griduptime2 .and. & 
        gdasforcdata(source)%gridchange1) then 
      
      write(LVT_logunit,*) '[INFO] GDAS -- changing grid to 1991 - 2000 (T126)'
      
      gdasforcdata(source)%ncold = 384
      gdasforcdata(source)%nrold = 190
      gdasforcdata(source)%mi = &
           gdasforcdata(source)%ncold*gdasforcdata(source)%nrold
      !-------------------------------------------------------------------
      ! Reinitialize the weights and neighbors
      !-------------------------------------------------------------------
      gridDesci = 0 
      gridDesci(1) = 4
      gridDesci(2) = 384
      gridDesci(3) = 190
      gridDesci(4) = 89.277
      gridDesci(5) = 0
      gridDesci(6) = 128
      gridDesci(7) = -89.277
      gridDesci(8) = -0.9375
      gridDesci(9) = 0.9375
      gridDesci(10) = 95
      gridDesci(20) = 0


      call gdas_reset_interp_input(source,gridDesci)
      gdasforcdata(source)%gridchange1 = .false.
      gdasforcdata(source)%datares = 0.9375
      != T170:
   elseif ( time1_tw >= gdasforcdata(source)%griduptime2 .and. &
        time1_tw < gdasforcdata(source)%griduptime3 .and. & 
        gdasforcdata(source)%gridchange2) then 
      
      write(LVT_logunit,*) '[INFO] GDAS -- changing grid to 2000 - 2002 (T170)'

      gdasforcdata(source)%ncold = 512
      gdasforcdata(source)%nrold = 256
      gdasforcdata(source)%mi = gdasforcdata(source)%ncold*gdasforcdata(source)%nrold
      !-------------------------------------------------------------------
      ! Reinitialize the weights and neighbors
      !-------------------------------------------------------------------
      gridDesci = 0
      gridDesci(1) = 4
      gridDesci(2) = 512
      gridDesci(3) = 256
      gridDesci(4) = 89.463
      gridDesci(5) = 0
      gridDesci(6) = 128
      gridDesci(7) = -89.463
      gridDesci(8) = -0.703125
      gridDesci(9) = 0.703125
      gridDesci(10) = 128
      gridDesci(20) = 0.0
      
      call gdas_reset_interp_input(source,gridDesci)
      gdasforcdata(source)%gridchange2 = .false.
      gdasforcdata(source)%datares = 0.703125
      
!= T254:
   elseif ( time1_tw >= gdasforcdata(source)%griduptime3 .and. & 
        time1_tw < gdasforcdata(source)%griduptime4 .and. &
        gdasforcdata(source)%gridchange3) then 
      
      write(LVT_logunit,*) '[ERR] GDAS -- changing grid to 2002 - 2005 (T254)'
      
      gdasforcdata(source)%ncold = 768
      gdasforcdata(source)%nrold = 384
      gdasforcdata(source)%mi = gdasforcdata(source)%ncold*gdasforcdata(source)%nrold
      !-------------------------------------------------------------------
      ! Reinitialize the weights and neighbors
      !-------------------------------------------------------------------
      gridDesci = 0
      gridDesci(1) = 4
      gridDesci(2) = 768
      gridDesci(3) = 384
      gridDesci(4) = 89.642
      gridDesci(5) = 0
      gridDesci(6) = 128
      gridDesci(7) = -89.642
      gridDesci(8) = -0.46875
      gridDesci(9) = 0.46875
      gridDesci(10) = 192
      gridDesci(20) = 0.0
      
      call gdas_reset_interp_input(source,gridDesci)
      gdasforcdata(source)%gridchange3 = .false.
      gdasforcdata(source)%datares = 0.46875
!= T382:
   elseif ( time1_tw >= gdasforcdata(source)%griduptime4 .and. &
        time1_tw < gdasforcdata(source)%griduptime5 .and. & 
        gdasforcdata(source)%gridchange4) then 
      
      write(LVT_logunit,*) '[INFO] GDAS -- changing grid to 2005 - 2010 (T382)'
      
      gdasforcdata(source)%ncold = 1152
      gdasforcdata(source)%nrold = 576
      gdasforcdata(source)%mi = gdasforcdata(source)%ncold*gdasforcdata(source)%nrold
      !-------------------------------------------------------------------
      ! Reinitialize the weights and neighbors
     !-------------------------------------------------------------------
      gridDesci = 0
      gridDesci(1) = 4
      gridDesci(2) = 1152
      gridDesci(3) = 576
      gridDesci(4) = 89.761
      gridDesci(5) = 0
      gridDesci(6) = 128
      gridDesci(7) = -89.761
      gridDesci(8) = -0.3125
      gridDesci(9) = 0.3125
      gridDesci(10) = 288
      gridDesci(20) = 0.0
      call gdas_reset_interp_input(source,gridDesci)
      gdasforcdata(source)%gridchange4 = .false.
      gdasforcdata(source)%datares = 0.3125
      != T574:
   elseif ( time1_tw >= gdasforcdata(source)%griduptime5 .and. &
        time1_tw < gdasforcdata(source)%griduptime6 .and. &
        gdasforcdata(source)%gridchange5) then
      
      write(LVT_logunit,*) '[INFO] GDAS -- changing grid to 2010 - 2015 (T574) --'
      
      gdasforcdata(source)%ncold = 1760
      gdasforcdata(source)%nrold = 880
      gdasforcdata(source)%mi = gdasforcdata(source)%ncold*gdasforcdata(source)%nrold
      !-------------------------------------------------------------------
      ! Reinitialize the weights and neighbors
      !-------------------------------------------------------------------
      gridDesci = 0
      gridDesci(1) = 4
      gridDesci(2) = 1760
      gridDesci(3) = 880
      gridDesci(4) = 89.844
      gridDesci(5) = 0
      gridDesci(6) = 128
      gridDesci(7) = -89.844
      gridDesci(8) = -0.204545454545455
      gridDesci(9) = 0.204545454545455
      gridDesci(10) = 440
      gridDesci(20) = 0.0
      
      call gdas_reset_interp_input(source, gridDesci)
      gdasforcdata(source)%gridchange5 = .false.
      gdasforcdata(source)%datares = 0.204545454545455
      != T1534:
   elseif ( time1_tw >= gdasforcdata(source)%griduptime6 .and. &
        gdasforcdata(source)%gridchange6 ) then
      
      write(LVT_logunit,*) '[INFO] GDAS -- changing grid to 2015 (T1534) --'
      
      gdasforcdata(source)%ncold = 3072
      gdasforcdata(source)%nrold = 1536
      gdasforcdata(source)%mi = gdasforcdata(source)%ncold*gdasforcdata(source)%nrold
      !-------------------------------------------------------------------
      ! Reinitialize the weights and neighbors
      !-------------------------------------------------------------------
      gridDesci = 0
      gridDesci(1) = 4
      gridDesci(2) = gdasforcdata(source)%ncold
      gridDesci(3) = gdasforcdata(source)%nrold
      gridDesci(4) = 89.91000
      gridDesci(5) = 0
      gridDesci(6) = 128
      gridDesci(7) = -89.91000
      gridDesci(8) = -0.1171875
      gridDesci(9) = 0.1171875
      gridDesci(10) = 768.0
      gridDesci(20) = 0.0
      
      call gdas_reset_interp_input(source, gridDesci)
      gdasforcdata(source)%gridchange6 = .false.
      gdasforcdata(source)%datares = 0.1171875
   endif

   if (alarmCheck) then      

      call create_GDASforc_filename(gdasforcdata(source)%odir,&
           yr2,mo2,da2,hr2,name00, name03, name06,F06flag)
      
      inquire(file=name00,exist=file_exists1) 
      inquire(file=name03,exist=file_exists2)

      if((.not.file_exists1).or.(.not.file_exists2)) then 

         call create_gdasf9_filename(gdasforcdata(source)%odir,&
              yr2,mo2,da2,hr2,name00, name03,status)

         if(status==1) then 
            write(LVT_logunit,*) '[WARN] backup GDAS f09 files are not ', &
                 'available for this hour'
         else
            inquire(file=name00,exist=file_exists1) 
            inquire(file=name03,exist=file_exists2)
            if ( .not. file_exists1 ) then
               status = 1
            endif
         endif
         if ( .not. file_exists1 ) then
            status = 1
         endif
         if ( .not. file_exists2 ) then
            status = 1
         endif
         if ( file_exists1 .and. file_exists2 ) then 
            status = 0 
         endif
      endif
         
      if ( F06flag ) then
         inquire(file=name06,exist=file_exists3)
         if ( .not. file_exists3 ) then
            status = 1
         endif
      endif
      
      if(status == 0) then 
         write(LVT_logunit,*) '[INFO] Reading GDAS file1 (I) ',trim(name00)
         write(LVT_logunit,*) '[INFO] Reading GDAS file1 (A) ',trim(name03)
         if ( F06flag ) then 
            write(LVT_logunit,*) '[INFO] Reading GDAS file1 (A)',trim(name06)
         endif

!--------------------------------------------------------------------------
! read 3hr forecast for time averaged fields
!--------------------------------------------------------------------------
         call retrieve_gdas_variables(source,name03,glbdata_a, ferror2)

!--------------------------------------------------------------------------
! read 6hr forecast for time averaged fields, if required. 
!--------------------------------------------------------------------------

         if(F06flag) then 
            call retrieve_gdas_variables(source,name06,glbdata_a_f06, ferror3)
         end if

         do r=1,LVT_rc%lnr
            do c=1,LVT_rc%lnc
               if(F06flag) then 
                  tmp_final(c,r) = 2*glbdata_a_F06(1,c,r)-glbdata_a(1,c,r)

                  prcp_final(1,c,r) = 2*glbdata_a_F06(2,c,r)-glbdata_a(2,c,r)
                  if(prcp_final(1,c,r).lt.0) then 
                     prcp_final(1,c,r) =0.0
                  endif
               else
                  tmp_final(c,r) = glbdata_a(1,c,r)
                  prcp_final(1,c,r) = glbdata_a(2,c,r)
               endif
            enddo
         enddo
         
      else
         tmp_final = LVT_rc%udef
         prcp_final = LVT_rc%udef
      endif
   else 
      tmp_final = LVT_rc%udef
      prcp_final = LVT_rc%udef
   end if ! file_exists
   
!   print*, LVT_rc%lnc, LVT_rc%lnr
!   open(100,file='test.bin',form='unformatted')
!   write(100) prcp_final(1,:,:)
!   close(100)
!   stop

   call LVT_logSingleDataStreamVar(LVT_MOC_TAIRFORC,&
        source,tmp_final,&
        vlevel=1,units='K')

   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,&
        source,prcp_final(1,:,:),&
        vlevel=1,units='kg/m2s')
   
   ! Now convert from kg/m2s to kg/m2
   do r=1,LVT_rc%lnr
      do c=1,LVT_rc%lnc
         if(prcp_final(1,c,r).ge.0) then 
            prcp_final(1,c,r) = prcp_final(1,c,r)*10800.0 !kg/m2
         else
            prcp_final(1,c,r) = LVT_rc%udef
         endif
      enddo ! c
   enddo ! r
   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final(1,:,:),&
        vlevel=1,units='kg/m2') 
      
end subroutine readGDASforcdata

!------------------------------------------------------------------------------

!BOP
! 
! !ROUTINE: retrieve_gdas_variables
! \label{retrieve_gdas_variables}
! 
! !INTERFACE: 
subroutine retrieve_gdas_variables(source, fname, glbdata, errorcode)
! !USES: 
  use LVT_coreMod
  use LVT_logMod
  use GDASforc_dataMod

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer               :: source
  character(len=*)      :: fname
  real                  :: glbdata(3,LVT_rc%lnc,LVT_rc%lnr)
  integer               :: errorcode
! 
! !DESCRIPTION: 
!   This subroutine retrieves GDAS forcing variables, and interpolates
!  them to the LVT grid. 
! 
!EOP

  integer               :: ngdas
  real, allocatable :: f(:)
  real, dimension(LVT_rc%lnc, LVT_rc%lnr) :: varfield
  integer :: igrib
  integer :: iv,c,r,t
  real    :: missingValue 
  integer :: iret
  integer :: ftn 
  integer :: pds5(3), pds7(3), pds6(3),pds16(3)
  logical :: var_status(3)
  integer :: pds5_val, pds7_val, pds16_val
  logical*1, allocatable :: lb(:)
  logical :: file_exists
  integer :: kk
  integer :: var_index
  integer :: nvars
  integer :: rc
  logical :: pcp_flag, var_found
  integer :: grid_size

#if(defined USE_GRIBAPI) 
!  pds5 = (/ 011,051,204,205,033,034,001,059,214,084 /) !parameter
!  pds6 = -1
!  pds7 = (/ 002,002,000,000,010,010,000,000,000,000 /) !htlev2
!! index 10 indicates instantaneous, 003 indicates time average
!  pds16 = (/010,010,003,003,010,010,010,003,003,003 /) 

! EMK...Drop last entry (for albedo); not used, and it exceeds array dimension
  pds5 = (/ 011, 059,214 /) !parameter
  pds6 = -1
  pds7 = (/ 002, 000,000 /) !htlev2
! index 10 indicates instantaneous, 003 indicates time average
  pds16 = (/010, 003,003 /) 

  ngdas = (gdasforcdata(source)%ncold*gdasforcdata(source)%nrold)

  varfield = 0 
  errorcode = 1
  var_status = .false. 

  inquire (file=fname, exist=file_exists)
  if (file_exists) then      

     call grib_open_file(ftn,trim(fname),'r',iret)
     if(iret.ne.0) then 
        write(LVT_logunit,*) '[ERR] Could not open file: ',trim(fname)
        errorcode = 0
        return
     endif

     call grib_count_in_file(ftn,nvars,iret)
     call LVT_warning(iret, 'error in grib_count_in_file in read_gdas')
     if(iret.ne.0) then 
        errorcode = 0
        return 
     endif

     allocate(lb(gdasforcdata(source)%ncold*gdasforcdata(source)%nrold))
     allocate(f(gdasforcdata(source)%ncold*gdasforcdata(source)%nrold))
     
     do kk=1,nvars
        call grib_new_from_file(ftn, igrib, iret)
        call LVT_warning(iret, 'error in grib_new_from_file in read_gdas')
        if(iret.ne.0) then 
           write(LVT_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        ! Trap the old "Could not find correct forcing parameter in file"
        ! error from LVT 6.  This error occurred right before a GDAS
        ! grid change.  LVT would try to read ahead, but the new data
        ! would be on the new grid, so LVT would misread it resulting in
        ! the above error message.  The LVT would roll back to the previous
        ! day for GDAS forcing.
        ! Trap this by checking the number of values in one of the
        ! GRIB fields.
        call grib_get_size(igrib,'values',grid_size)
        if ( grid_size /= ngdas ) then
           write(LVT_logunit,*) &
              '[ERR] Number of values does not match expected', trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)
           return
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LVT_verify(rc, 'error in grib_get: indicatorOfParameter in read_gdas')

        call grib_get(igrib,'level',pds7_val,rc)
        call LVT_verify(rc, 'error in grib_get: level in read_gdas')

        call grib_get(igrib,'timeRangeIndicator',pds16_val,rc)
        call LVT_verify(rc, 'error in grib_get: timeRangeIndicator in read_gdas')

        var_found = .false. 
        do iv=1,3
           if((pds5_val.eq.pds5(iv)).and.&
                (pds7_val.eq.pds7(iv)).and.&
                (pds16_val.eq.pds16(iv))) then
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo
        f = -9999.0
        call grib_get(igrib,'values',f,rc)
        call LVT_warning(rc, 'error in grib_get:values in read_gdas')

        if(rc.ne.0) then 
           write(LVT_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'missingValue',missingValue,rc)
        call LVT_verify(rc, 'error in grib_get:missingValue in read_gdas')

        call grib_release(igrib,rc)
        call LVT_verify(rc, 'error in grib_release in read_gdas')
        
        if(var_found) then 
           lb = .false.
           do t=1,ngdas
              if(f(t).ne.missingValue) lb(t) = .true. 
           enddo
           
           call interp_gdas(source,ngdas,f,lb,varfield)

           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 glbdata(var_index,c,r) =&
                      varfield(c,r)
              enddo
           enddo
        endif

     enddo
     call grib_close_file(ftn)

     deallocate(lb)
     deallocate(f)     
         
     do kk=1,2
        if(.not.var_status(kk)) then 
           write(LVT_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           
           return
        endif
     enddo
  else
     write(LVT_logunit,*) &
          '[ERR] Could not find file: ',trim(fname)
     errorcode = 0
  endif
#endif     
end subroutine retrieve_gdas_variables


!BOP
! !ROUTINE: interp_gdas
! \label{interp_gdas}
!
! !INTERFACE:
subroutine interp_gdas(source,ngdas,f,lb,varfield)
! !USES:

  use LVT_coreMod
  use GDASforc_dataMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: source
  integer, intent(in) :: ngdas
  real, intent(out)   :: f(ngdas)
  logical*1           :: lb(ngdas)
  real, intent(out)   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given GDAS field 
!   to the LVT grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[ngdas]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmapgg
! \item[lis\_gds]
!  array description of the LVT grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LVT grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LVT grid
! \item[varfield]
!  output interpolated field
!  \end{description} 
! 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
! \end{description}
!EOP
  integer :: iret
  integer :: count1,i,j,mo

  real, dimension(LVT_rc%lnc*LVT_rc%lnr) :: lis1d

  logical*1 :: lo(LVT_rc%lnc*LVT_rc%lnr)

!=== End variable declarations
!-----------------------------------------------------------------------
! Setting interpolation options (ip=0,bilinear)
! (km=1, one parameter, ibi=1,use undefined bitmap
! (needed for soil moisture and temperature only)
! Use budget bilinear (ip=3) for precip forcing fields
!-----------------------------------------------------------------------
  mo = LVT_rc%lnc*LVT_rc%lnr
!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.
!-----------------------------------------------------------------------  
! Interpolate to LVT grid
  !-----------------------------------------------------------------------
!  print*, gdasforcdata(source)%ncold, gdasforcdata(source)%nrold,ngdas
!  open(100,file='test_inp.bin',form='unformatted')
!  write(100) f
!  close(100)
!  stop

  if (LVT_isAtAFinerResolution(gdasforcdata(source)%datares)) then
     call conserv_interp(LVT_rc%gridDesc,lb,f, &
          lo,lis1d, gdasforcdata(source)%ncold*gdasforcdata(source)%nrold, &
          LVT_rc%lnc*LVT_rc%lnr, gdasforcdata(source)%rlat, &
          gdasforcdata(source)%rlon, &
          gdasforcdata(source)%w112, gdasforcdata(source)%w122, &
          gdasforcdata(source)%w212, gdasforcdata(source)%w222, &
          gdasforcdata(source)%n112, gdasforcdata(source)%n122, &
          gdasforcdata(source)%n212, gdasforcdata(source)%n222, &
          LVT_rc%udef, iret)     
  else
     call upscaleByAveraging(&
          gdasforcdata(source)%ncold*gdasforcdata(source)%nrold, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          gdasforcdata(source)%n11, lb, &
          f, lo, lis1d)
  end if
!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between GDAS & LDAS. For LDAS land 
! points not included in GDAS geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, LVT_rc%lnr
     do i = 1, LVT_rc%lnc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + LVT_rc%lnc
  enddo

end subroutine interp_gdas

!BOP
!
! !ROUTINE: gdas_reset_interp_input
!  \label{gdas_reset_interp_input}
!
! !REVISION HISTORY:
!  01 Feb 2016: James Geiger: Initial specification
!
! !INTERFACE:
subroutine gdas_reset_interp_input(source, gridDesci)
! !USES:

   use LVT_coreMod
   use LVT_logMod
   use GDASforc_dataMod

   implicit none
! !ARGUMENTS: 
   integer, intent(in) :: source
   real, intent(in)    :: gridDesci(50)

! !DESCRIPTION:
! Resets the neighbours and weights arrays used for spatially
! interpolating the GDAS forcing data to the LVT running domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[gridDesci]
!    array of magic numbers describing the GDAS forcing domain
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbors for upscaling by averaging
!   \item[LVT\_isatAfinerResolution](\ref{LVT_isatAfinerResolution}) \newline
!    determines whether LVT' running domain is at a finer resolution
!    than the given resolution
!  \end{description}
!EOP

   integer :: rc

   if(allocated(gdasforcdata(source)%n112)) then 
      deallocate(gdasforcdata(source)%n112, stat=rc)
      deallocate(gdasforcdata(source)%n122, stat=rc)
      deallocate(gdasforcdata(source)%n212, stat=rc)
      deallocate(gdasforcdata(source)%n222, stat=rc)
      deallocate(gdasforcdata(source)%w112, stat=rc)
      deallocate(gdasforcdata(source)%w122, stat=rc)
      deallocate(gdasforcdata(source)%w212, stat=rc)
      deallocate(gdasforcdata(source)%w222, stat=rc)
   endif
   if(allocated(gdasforcdata(source)%n11)) then 
      deallocate(gdasforcdata(source)%n11)
   endif

   if ( LVT_isatAfinerResolution(gridDesci(9)) ) then

            
      allocate(gdasforcdata(source)%n112(LVT_rc%lnc*LVT_rc%lnr,25))
      allocate(gdasforcdata(source)%n122(LVT_rc%lnc*LVT_rc%lnr,25))
      allocate(gdasforcdata(source)%n212(LVT_rc%lnc*LVT_rc%lnr,25))
      allocate(gdasforcdata(source)%n222(LVT_rc%lnc*LVT_rc%lnr,25))
      allocate(gdasforcdata(source)%w112(LVT_rc%lnc*LVT_rc%lnr,25))
      allocate(gdasforcdata(source)%w122(LVT_rc%lnc*LVT_rc%lnr,25))
      allocate(gdasforcdata(source)%w212(LVT_rc%lnc*LVT_rc%lnr,25))
      allocate(gdasforcdata(source)%w222(LVT_rc%lnc*LVT_rc%lnr,25))
      
      call conserv_interp_input(gridDesci,LVT_rc%gridDesc,&
           LVT_rc%lnc*LVT_rc%lnr, &
           gdasforcdata(source)%rlat, gdasforcdata(source)%rlon,&
           gdasforcdata(source)%n112, gdasforcdata(source)%n122, &
           gdasforcdata(source)%n212, gdasforcdata(source)%n222, & 
           gdasforcdata(source)%w112, gdasforcdata(source)%w122, &
           gdasforcdata(source)%w212, gdasforcdata(source)%w222)

   else

      allocate(gdasforcdata(source)%n11(gdasforcdata(source)%mi))
      
      call upscaleByAveraging_input(gridDesci,LVT_rc%gridDesc,&
           gdasforcdata(source)%ncold*gdasforcdata(source)%nrold, &
           LVT_rc%lnc*LVT_rc%lnr, &
           gdasforcdata(source)%n11)
   endif
end subroutine gdas_reset_interp_input


subroutine create_GDASforc_filename(gdasdir,yr,mo,da,hr,name00, name03, name06,&
     F06flag)


  use LVT_timeMgrMod, only : LVT_tick
  ! Defaults
  implicit none

! !ARGUMENTS: 
  character(len=*), intent(in)    :: gdasdir
  integer, intent(in)             :: yr, mo, da, hr
  character(len=*), intent (out)  :: name00
  character(len=*), intent (out)  :: name03
  character(len=*), intent (out)  :: name06
  logical                         :: F06flag
! !DESCRIPTION:
!   This subroutine puts together GDAS file names for 
!   different (0,3,6hr) forecast instances
!
!   LVT uses a mix of instantaneous forcing values and averaged forcing
!   values from the GDAS forcing data.  LVT processes the GDAS forcing at
!   a 3-hourly period.  When it is time to read new GDAS forcing data,
!   LVT reads 3 hours ahead.  So when LVT' clock is at hour hr, LVT
!   processes data for the period ]hr, hr+3].  Here, LVT reads instantaneous
!   forcing valid at hour hr$+3$ and averaged forcing valid for [hr, hr+3]
!
!   GDAS files ``named'' 00.f00, 06.f00, 12.f00, and 18.f00 are the
!   analysis files, and they contain only instantaneous forcing values at
!   hours 00, 06, 12, and 18 respectively.  The files hr.f03, hr.f06,
!   and hr.f09 are, respectively, 3-hour, 6-hour, and 9-hour forecasts
!   from the analysis hour hr, where hr is either 00, 06, 12, or 18.
!
!   The hr.f03 files contain instantaneous forcing values at hour hr+3, and
!   they contain forcing values averaged over the 3-hour period [hr, hr+3].
!
!   The hr.f06 files contain instantaneous forcing values at hour hr+6, and
!   they contain forcing values averaged over the 6-hour period [hr, hr+6].
!
!   The hr.f09 files contain instantaneous forcing values at hour hr+9, and
!   they contain forcing values averaged over the 3-hour period [hr+6, hr+9].
! 
!   Again, LVT processes the GDAS forcing at a 3-hourly period.  When the
!   hour of LVT' clock, hr, is one of the analysis hours (00, 06, 12, 18)
!   then LVT reads ahead to extract instantaneous forcing values from
!   the hr.f03 file, and it reads ahead to extract averaged forcing values
!   from the hr.f03 file.
!
!   When the hour of LVT' clock, hr, is not one of the analysis hours
!   (in this case when hr is 03, 09, 15, or 21), then LVT reads instantaneous
!   forcing values from the from the next analysis hour.
!   E.g., when hr is 03, read from 06.f00; when hr is 09, read from 12.f00.
!   LVT computes averaged forcing values by reading the averaged values
!   from the 3-hour and 6-hour forecasts of the previous analysis hour.
!   E.g., when hr is 03, read from 00.f03 and 00.f06; when hr is 09, read
!   from 06.f03 and 06.f06.
!   For this case, the averaged forcing values are computed by subtracting
!   the 3-hour forecast from the 6-hour forecast both of the previous analysis
!   hour thusly:
!
!   Consider LVT' hour to be 15.  LVT needs averaged forcing valid for the
!   period [15, 18].  The previous analysis hour is then 12.
!   We want the averaged forcing values from 12.f03 and from 12.f06.
!
!   hr.f06 contains averaged forcing, denoted $f_6$, for the 
!   period [hr, hr+6] and hr.f03 contains averaged forcing, denoted $f_3$,
!   for the period [hr, hr+3].  Thus for forcing, $f_6$, $6*f_6$ is the
!   total $f$ at hour hr$+6$ and $3*f_3$ is the total $f$ at hour hr$+3$.
!   $( 6*f_6 - 3*f_3 )$ represents the total $f$ from hr$+3$ to hr$+6$,
!   so $( 6*f_6 - 3*f_3 ) / 3 \equiv 2*f_6 - f_3$ is the averaged $f$ for
!   the period [hr+3, hr+6].
!   
!   So to compute the averaged forcing for period [15, 15+3], LVT must
!   read 12.f03 and 12.f06.  $f_6$ represents averaged $f$ for [12, 18];
!   $f_3$ represents averaged $f$ for [12, 15]; $2*f_6 - f_3$ represents 
!   avergage $f$ for [15, 18], the desired period.
!
!   Note that the variable, F06flag, is used to indicate the above
!   situation when LVT must subtract $f_3$ values from $f_6$ values to
!   obtain the desired averaged forcing values.
!
!   Note that the logic for reading and creating names for bookend1, which
!   is only done at initialization, is the same as described above except
!   imagine that you are at start hour - 3 reading ahead to the start hour.
!
!  The arguments are:
!  \begin{description}
!  \item[option]
!    bookend option (1 or 2)
!  \item[gdasdir]
!    Name of the GDAS directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[name00]
!   name of the GDAS file for reading instantaneous fields
!   \item[name03]
!   name of the 3hr GDAS file for reading time averaged fields
!   \item[name06]
!   name of the 6hr GDAS file for reading time averaged fields
!  \item[F06flag]
!    flag to indicate if 6hr forecast data is required for this interval
!  \end{description}
!
!EOP
  integer :: option
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss
  integer :: uyr0, umo0, uda0, uhr0, umn0, uss0
  logical :: is_analysis_hr
  character(len=2)  :: analysis_hour_inst, analysis_hour_avg, &
                       fcstcode0, fcstcode1, fcstcode2
  character(len=80) :: fbase
  character(len=8)  :: fdir
  character(len=10) :: ftime
  character(len=21) :: fsubs
  real*8      :: dumbtime
  integer     :: doy
  real        :: gmt
!=== End Variable Definition ===============

  option = 2
!=== formats for filename segments
!-----------------------------------------------------------------
!  Make variables for the time used to create the file
!  We don't want these variables being passed out
!-----------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0

  uyr0 = yr
  umo0 = mo
  uda0 = da
  uhr0 = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn0 = 0
  uss0 = 0

  if ( uhr == 0 .or. uhr == 6 .or. uhr == 12 .or. uhr == 18 ) then
     is_analysis_hr = .true.
  else
     !( uhr == 3 .or. uhr == 9 .or. uhr == 15 .or. uhr == 21 )
     is_analysis_hr = .false.
  endif

!-----------------------------------------------------------------
! Cheat sheet
!
! bookend 2: 
!  hour 00 look for 00F3
!  hour 03 look for 00F3 and 00F6
!  hour 06 look for 06F3
!  hour 09 look for 06F3 and 06F6
!  hour 12 look for 12F3
!  hour 15 look for 12F3 and 12F6
!  hour 18 look for 18F3
!  hour 21 look for 18F3 and 18F6
!
! bookend 1: 
!  hour 00 look for 18F3 and 18F6
!  hour 03 look for 00F3
!  hour 06 look for 00F3 and 00F6
!  hour 09 look for 06F3
!  hour 12 look for 06F3 and 06F6
!  hour 15 look for 12F3
!  hour 18 look for 12F3 and 12F6
!  hour 21 look for 18F3
!-----------------------------------------------------------------

  F06flag = .false. 

  if ( option == 2 ) then !bookend 2
     if ( .not. is_analysis_hr ) F06flag = .true. ! need to read two files     
  else
     if ( is_analysis_hr ) F06flag = .true. ! need to read two files     
  endif

  if ( option == 2 ) then !bookend 2 
     if ( uhr == 0 .or. uhr == 6 .or. uhr == 12 .or. uhr == 18 ) then
        ! read averged forcing from analysis_hour file ( 00, 06, 12, 18)
        write(unit=analysis_hour_avg,fmt='(i2.2)') uhr

        ! read instantaneous forcing from analysis_hour file ( 00, 06, 12, 18)
        write(unit=analysis_hour_inst,fmt='(i2.2)') uhr
        fcstcode0 = '03'
     else
        ! read averaged forcing from previous analysis_hour file
        ! (03 -> 00, 09 -> 06, 15 -> 12, 21 -> 18)
        write(unit=analysis_hour_avg,fmt='(i2.2)') uhr - 3

        ! read instantaneous forcing from next analysis_hour file
        ! (03 -> 06, 09 -> 12, 15 -> 18, 21 -> 00)
        write(unit=analysis_hour_inst,fmt='(i2.2)') uhr + 3
        fcstcode0 = '00'
     endif

     ! special case: need to go to the next day
     if ( uhr .eq. 21 ) then
        call LVT_tick(dumbtime,doy,gmt,uyr,umo,uda,uhr,umn,uss,24*60*60)
        analysis_hour_inst = '00'
     endif

  else !bookend 1
     if ( uhr == 0 .or. uhr == 6 .or. uhr == 12 .or. uhr == 18 ) then
        ! read averged forcing from previous analysis_hour file
        ! ( 00 -> -06, 06 -> 00, 12 -> 06, 18 -> 12)
        if(uhr.ne.0) &
             write(unit=analysis_hour_avg,fmt='(i2.2)') uhr - 6

        ! read instantaneous forcing from analysis_hour file ( 00, 06, 12, 18)
        write(unit=analysis_hour_inst,fmt='(i2.2)') uhr
        fcstcode0 = '00'
     else
        ! read averaged forcing from previous analysis_hour file
        ! (03 -> 00, 09 -> 06, 15 -> 12, 21 -> 18)
        write(unit=analysis_hour_avg,fmt='(i2.2)') uhr - 3

        ! read instantaneous forcing from previous analysis_hour file
        ! (03 -> 00, 09 -> 06, 15 -> 12, 21 -> 18)
        write(unit=analysis_hour_inst,fmt='(i2.2)') uhr - 3
        fcstcode0 = '03'
     endif

     ! special case: need to go to the previous day
     if ( uhr0 == 0 ) then
        call LVT_tick(dumbtime,doy,gmt,uyr0,umo0,uda0,uhr0,umn0,uss0,-24*60*60)
        analysis_hour_avg = '18'
     endif

  endif

  fcstcode1 = '03'
  fcstcode2 = '06'
  

  ! construct file names

  !name00
  write(unit=fdir, fmt='(a1, i4.4, i2.2, a1)') '/', uyr, umo, '/'
  write(unit=ftime, fmt='(i4.4, i2.2, i2.2, a2)') uyr, umo, uda, analysis_hour_inst
  write(unit=fsubs, fmt='(a16, a2, a3)') '.gdas1.sfluxgrbf', fcstcode0, '.sg'
  name00 = trim(gdasdir)//fdir//ftime//fsubs
  
  !name03
  write(unit=fdir, fmt='(a1, i4.4, i2.2, a1)') '/', uyr0, umo0, '/'
  write(unit=ftime, fmt='(i4.4, i2.2, i2.2, a2)') uyr0, umo0, uda0, analysis_hour_avg
  write(unit=fsubs, fmt='(a16, a2, a3)') '.gdas1.sfluxgrbf', fcstcode1, '.sg'
  name03 = trim(gdasdir)//fdir//ftime//fsubs

  !name06
  write(unit=fdir, fmt='(a1, i4.4, i2.2, a1)') '/', uyr0, umo0, '/'
  write(unit=ftime, fmt='(i4.4, i2.2, i2.2, a2)') uyr0, umo0, uda0, analysis_hour_avg
  write(unit=fsubs, fmt='(a16, a2, a3)') '.gdas1.sfluxgrbf', fcstcode2, '.sg'
  name06 = trim(gdasdir)//fdir//ftime//fsubs

end subroutine create_GDASforc_filename


subroutine create_gdasf9_filename(gdasdir,yr,mo,da,hr,name00, name03,status)
! !USES: 
  use LVT_timeMgrMod, only : LVT_tick

  implicit none
! !ARGUMENTS: 
  character(len=*), intent(in)    :: gdasdir
  integer, intent(in)             :: yr, mo, da, hr
  character(len=*), intent (out)  :: name00
  character(len=*), intent (out)  :: name03
  integer           :: status

! !DESCRIPTION:
!   This subroutine puts together 9hr forecast GDAS file name
!
!   First read the notes in
!   create\_gdasfilename.F90 (\ref{create_gdasfilename}).
!
!   This routine is used to create the GDAS file name for the backup 9hr 
!   forecast file.  At analysis hours (hr = 00, 06, 12, and 18), LVT
!   reads the hr.f03 GDAS file for both the new instantaneous and the
!   new averaged forcing values (here, name00 == name03). Should this file be
!   missing, the 9-hour forecast of the previous analysis hour contains
!   comparable forcing values.
!
!   An hr.f09 file contains instantaneous forcing values valid at hour hr+9,
!   and it contains 3-hourly averaged forcing values valid for [hr+6, hr+9]
!
!   For example, when LVT is at hour 06, it reads ahead to get instantaneous
!   forcing values at hour 09 and to get averaged forcing data for the
!   period [06, 09].  The GDAS file 06.f03 contains these data.  Should
!   this file be missing, the GDAS file 00.f09 also contains instantaneous
!   forcing values at hour 09 and averged forcing values for the
!   period [06, 09].
!
!   Note that backup 9-hour forecast files exist only for anaylsis hours.
!   Should a name00 or name03 file be missing when LVT is at a non-analysis
!   hour (hr = 03, 09, 15, or 21), then there is no previous hour with a 9-hour
!   forecast that contains comparable forcing values.  Another strategy must
!   be developed to replace those missing files.
!
!   Note that the logic for reading and creating names for bookend1, which
!   is only done at initialization, is the same as described above except
!   imagine that you are at start hour - 3 reading ahead to the start hour.
!
!   Note that the case where LVT is reading name00, name03, and name06
!   data to read ahead is not handled by this routine.
!   Example, consider the case when reading ahead to hour 18.
!   Here LVT is at hour 15, and LVT needs 18f00 for instantaneous values
!   at hour 18, and LVT needs 12f03 and 12f06 for averaged values for the
!   period [15, 18].  12f03 contains averaged values for [12, 15].  If it
!   is missing, LVT could use 06f09 to obtain averaged values for [12, 15].
!   Again, this scenario is not handled by this routine.  LVT will simply
!   roll back to the the previous day to look for data.
! 
!  The arguments are:
!  \begin{description}
!  \item[option]
!    bookend option (1 or 2)
!   \item[name00]
!   name of the timestamped 9hr GDAS file
!   \item[name03]
!   name of the timestamped 9hr GDAS file
!  \item[gdasdir]
!    Name of the GDAS directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!  \item[status]
!   flag indicating if the 9hr forecast is available. 
!  \end{description}
!
!EOP
  integer           :: option
  integer           :: i, c
  integer           :: uyr, umo, uda, uhr, umn, uss
  integer           :: uyr0, umo0, uda0, uhr0, umn0, uss0
  character(len=2)  :: initcode1,  fcstcode1
  character(len=8)  :: fdir
  character(len=10) :: ftime
  character(len=21) :: fsubs
  real*8            :: dumbtime
  integer           :: doy1, doy
  real              :: gmt
!=== End Variable Definition ===============

  option = 2
!=== formats for filename segments
!-----------------------------------------------------------------
!  Make variables for the time used to create the file
!  We don't want these variables being passed out
!-----------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0

  uyr0 = yr
  umo0 = mo
  uda0 = da
  uhr0 = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn0 = 0
  uss0 = 0

!-----------------------------------------------------------------
! Cheat sheet
!
! bookend 2: 
!  hour 00 look for 18F9
!  hour 06 look for 00F9
!  hour 12 look for 06F9
!  hour 18 look for 12F9
!
! bookend 1: 
!  hour 03 look for 18F9
!  hour 09 look for 00F9
!  hour 15 look for 06F9
!  hour 21 look for 12F9
!-----------------------------------------------------------------

  if ( option == 2 ) then !bookend 2
     ! backup f09 files exist only for analysis hours
     if ( uhr == 0 .or. uhr == 6 .or. uhr == 12 .or. uhr == 18 ) then
        status = 0

        ! read forcing from previous analysis_hour file
        write(unit=initcode1,fmt='(i2.2)') uhr-6

        !special case: need to go back to the previous day 18z
        if ( uhr0 == 0 ) then
           initcode1 = '18'
           call LVT_tick(dumbtime, doy, gmt, uyr0, umo0, uda0, &
                uhr0, umn0, uss0, -1*6*60*60)
        endif
     else
        status = 1
     endif
  else !bookend 1
     if (uhr == 3 .or. uhr == 9 .or. uhr == 15 .or. uhr == 21 ) then 
        status = 0

        ! read forcing from previous analysis_hour file
        write(unit=initcode1,fmt='(i2.2)') uhr-9

        !special case: need to go back to the previous day 18z
        if ( uhr0 == 3 ) then
           initcode1 = '18'
           call LVT_tick(dumbtime, doy, gmt, uyr0, umo0, uda0, &
                         uhr0, umn0, uss0, -9*60*60)
        endif
     else
        status = 1
     endif
  endif

  if ( status == 0 ) then
     fcstcode1 = '09'
    
     !name00
     write(UNIT=fdir, fmt='(a1, i4.4, i2.2, a1)') '/', uyr0, umo0, '/'
     write(UNIT=ftime, fmt='(i4.4, i2.2, i2.2, a2)') uyr0, umo0, uda0, initcode1
     write(UNIT=fsubs, fmt='(a16, a2, a3)') '.gdas1.sfluxgrbf', fcstcode1, '.sg'
     name00 = trim(gdasdir)//fdir//ftime//fsubs

     !name03
     name03 = name00
  endif

end subroutine create_gdasf9_filename
