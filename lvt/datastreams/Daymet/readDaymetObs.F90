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
! !ROUTINE:
!
! !INTERFACE:
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
!  29 Sept 2016   Sujay Kumar  Initial Specification
!EOP
!BOP
! 
! !ROUTINE: readDaymetObs
! \label{readDaymetObs}
! 
! !REVISION HISTORY: 
!  23 APR 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readDaymetObs(source)
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use Daymet_obsMod,    only : daymetobs
  use map_utils
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for Daymet data
! LVT expects the data to be organized per calendar year, with 
! each file containing a daily data. Each reported observation is
! assumed to be time averaged. 
!
! The yearly file is read until the data for the current time is 
! found. The data is then upscaled to the LVT analysis grid. 
! 
!EOP
  integer                :: source
  integer                :: yr, mo, da, hr
  integer                :: i,j,t,c,r
  integer                :: stn_col, stn_row
  real                   :: col,row
  character*100          :: daymetfilename
  logical                :: file_exists
  logical                :: readflag
  integer                :: ftn, sweid, ios
  integer                :: status
  type(ESMF_Time)        :: daymettime, daymettime1
  integer                :: stnindex,tind
  real                   :: offset
  real                   :: swe_data
  integer                :: iret
  integer                :: latid, lonid
  real, allocatable      :: lat(:,:),lon(:,:)
  real, allocatable      :: rlat(:),rlon(:) 
  logical*1              :: lb(daymetobs(source)%nc*daymetobs(source)%nr)
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: swe(daymetobs(source)%nr, daymetobs(source)%nc)
  real                   :: swe1d(daymetobs(source)%nc*daymetobs(source)%nr)
  real                   :: swe_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: gridDesci(50)
      
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  swe_ip  = LVT_rc%udef

  call ESMF_TimeSet(daymettime1,yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
       h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'daymettime1 set failed')

  offset = (daymettime1-daymetobs(source)%starttime)/&
       daymetobs(source)%timestep
  
  if((nint(offset)-offset).eq.0) then 
  
     call create_Daymet_swefilename(daymetobs(source)%odir, &
          LVT_rc%dyr(source), &
          daymetfilename)
     
     inquire(file=trim(daymetfilename),exist=file_exists)
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading Daymet data from ',&
             trim(daymetfilename)
        
        if(daymetobs(source)%startMode) then 
           daymetobs(source)%startMode = .false. 
           
           gridDesci = 0 
           gridDesci(1) = 3
           gridDesci(2) = daymetobs(source)%nc
           gridDesci(3) = daymetobs(source)%nr
           griddesci(4) = 6.08138
           gridDesci(5) = -179.996
           gridDesci(6) = 8
           gridDesci(7) = 42.5
           gridDesci(8) = 1
           gridDesci(9) = 1
           gridDesci(10) = 25.0
           gridDesci(11) = -100.0
           gridDesci(20) = 0 
           
           allocate(daymetobs(source)%n11(daymetobs(source)%nc*daymetobs(source)%nr))

           allocate(rlat(daymetobs(source)%nc*daymetobs(source)%nr))
           allocate(rlon(daymetobs(source)%nc*daymetobs(source)%nr))

           allocate(lat(daymetobs(source)%nr, daymetobs(source)%nc))
           allocate(lon(daymetobs(source)%nr, daymetobs(source)%nc))
           
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
           ios=nf90_open(path=trim(daymetfilename),mode=NF90_NOWRITE,ncid=ftn)
           call LVT_verify(ios,'Error in nf90_open: readDaymetObs')
           
           ios=nf90_inq_varid(ftn,'lat',latid)
           call LVT_verify(ios,'Error in nf90_inq_varid: readDaymetObs')
           
           ios=nf90_inq_varid(ftn,'lon',lonid)
           call LVT_verify(ios,'Error in nf90_inq_varid: readDaymetObs')
           
           ios=nf90_get_var(ftn,latid, lat)
           call LVT_verify(ios,'Error in nf90_get_var: readDaymetObs')
           
           ios=nf90_get_var(ftn,lonid, lon)
           call LVT_verify(ios,'Error in nf90_get_var: readDaymetObs')
           
           ios=nf90_close(ftn)
           call LVT_verify(ios,'Error in nf90_close: readDaymetObs')
    

           do r=1,daymetobs(source)%nr
              do c=1,daymetobs(source)%nc
                 rlat(c+(r-1)*daymetobs(source)%nc) = &
                      lat(c,daymetobs(source)%nr-r+1)
                 rlon(c+(r-1)*daymetobs(source)%nc) = &
                      lon(c,daymetobs(source)%nr-r+1)
              enddo
           end do
#endif    
           call upscaleByAveraging_input_with_latlon(&
                gridDesci,LVT_rc%gridDesc(:),&
                daymetobs(source)%nc*daymetobs(source)%nr,&
                LVT_rc%lnc*LVT_rc%lnr,&
                daymetobs(source)%n11,&
                rlat, rlon)

           deallocate(rlat)
           deallocate(rlon)
           deallocate(lat)
           deallocate(lon)

        endif

        ios=nf90_open(path=trim(daymetfilename),mode=NF90_NOWRITE,ncid=ftn)
        call LVT_verify(ios,'Error in nf90_open: readDaymetObs')
        
        ios=nf90_inq_varid(ftn,'swe',sweid)
        call LVT_verify(ios,'Error in nf90_inq_varid: readDaymetObs')

        ios=nf90_get_var(ftn,sweid, swe, start=(/1,1,LVT_rc%doy/),&
             count=(/daymetobs(source)%nr, daymetobs(source)%nc,1/))
        if(ios.ne.0) then 
           swe = LVT_rc%udef
        endif
        
        ios=nf90_close(ftn)
        call LVT_verify(ios,'Error in nf90_close: readDaymetObs')


        lb = .true. 
        do r=1,daymetobs(source)%nr
           do c=1,daymetobs(source)%nc
              swe1d(c+(r-1)*daymetobs(source)%nc) = &
                   swe(c,daymetobs(source)%nr-r+1)
              if(swe(c,daymetobs(source)%nr-r+1).lt.0) &
                   lb(c+(r-1)*daymetobs(source)%nc) = .false. 
           enddo
        end do
        call upscaleByAveraging(&
             daymetobs(source)%nc*daymetobs(source)%nr,&
             LVT_rc%lnc*LVT_rc%lnr,&
             LVT_rc%udef, &
             daymetobs(source)%n11,&
             lb,swe1d,lo,swe_ip)
        
        write(LVT_logunit,*) '[INFO] Successfully processed from ',&
        trim(daymetfilename), 'for day ',LVT_rc%doy
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,&
       swe_ip,vlevel=1,units="kg/m2")

#endif

end subroutine readDaymetObs

!BOP
! 
! !ROUTINE: create_Daymet_swefilename
!  \label(create_Daymet_swefilename)
!
! !INTERFACE: 
subroutine create_Daymet_swefilename(odir, yr,daymetname)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: daymetname
!
! !DESCRIPTION: 
! 
! This routine creates a filename for the Daymet SWE file
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir] CMC snow data base directory
!   \item[yr]   year of data 
!   \item[daymetname]  Name of the Daymet SWE file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  
  write(fyr, '(i4.4)' ) yr

  daymetname = trim(odir)//'/daymet_v3_swe_'&
       //trim(fyr)//'_na.nc4'
  
end subroutine create_Daymet_swefilename
