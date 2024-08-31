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
! !ROUTINE: readOzFluxObs
! \label{readOzFluxObs}
!
! !INTERFACE: 
subroutine readOzFluxObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_logMod,       only : LVT_logunit, LVT_verify
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_calendar, LVT_tick
  use OzFlux_obsMod,      only : ozfluxobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for OzFlux station data. 
! The plugin processes latent, sensible, and ground heat flux data
! from the in-situ measurements. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  1 Apr 2020: Sujay Kumar, Initial Specification
! 
!EOP
!BOP
! !ARGUMENTS: 

!EOP
  type(ESMF_Time)  :: ozfluxtime1
  integer          :: t,c,r
  real             :: gmt
  real                :: time
  type(ESMF_TimeInterval) :: dayInterval
  type(ESMF_Time)         :: initTime
  integer             :: yr, mo, da, hr, mn, ss, doy
  integer             :: status
  integer             :: k
  integer             :: data_index
  real                :: qle(LVT_rc%lnc, LVT_rc%lnr)
  real                :: qh(LVT_rc%lnc, LVT_rc%lnr)
  real                :: qg(LVT_rc%lnc, LVT_rc%lnr)
  
  qle  = LVT_rc%udef 
  qh   = LVT_rc%udef 
  qg   = LVT_rc%udef 


!  time = LVT_rc%dhr(source)*3600+LVT_rc%dmn(source)*60+LVT_rc%dss(source)
!  if((mod(time,86400.0).eq.0.0).or.&
  if(LVT_rc%dyr(source).ne.ozfluxobs(source)%yr) then 

     call ESMF_TimeSet(ozfluxobs(source)%starttime, &
          yy =LVT_rc%dyr(source), &
          mm = LVT_rc%dmo(source), &
          dd = LVT_rc%dda(source), &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status,'error in timeset: readOzfluxobs(Source)')
     ozfluxobs(source)%yr = LVT_rc%dyr(source)

     call process_ozflux_data(source,LVT_rc%dyr(source),&
          LVT_rc%dmo(source),LVT_rc%dda(source))

  endif

  call ESMF_TimeSet(ozfluxtime1,yy=LVT_rc%dyr(source), &
       mm = LVT_rc%dmo(source), &
       dd = LVT_rc%dda(source), &
       h = LVT_rc%dhr(source), &
       m = LVT_rc%dmn(source), &
       calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'error in timeset: readOzfluxobs(Source)')

  data_index = nint((ozfluxtime1 - ozfluxobs(source)%starttime)/&
       ozfluxobs(source)%ts) + 1

  do k=1,ozfluxobs(source)%n_stns
     qle(ozfluxobs(source)%stncol(k),ozfluxobs(source)%stnrow(k)) = &
          ozfluxobs(source)%qle(k,data_index)

     qh(ozfluxobs(source)%stncol(k),ozfluxobs(source)%stnrow(k)) = &
          ozfluxobs(source)%qh(k,data_index)

     qg(ozfluxobs(source)%stncol(k),ozfluxobs(source)%stnrow(k)) = &
          ozfluxobs(source)%qg(k,data_index)


  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_qle,source,qle,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_qh,source,qh,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_qg,source,qg,vlevel=1,units="W/m2",dir="DN")
  
end subroutine readOzfluxobs


!BOP
! 
! !ROUTINE: process_ozflux_data
!  \label{process_ozflux_data}
!
! !INTERFACE:
subroutine process_ozflux_data(source,yr,mo,da)
! 
! !USES:
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod,      only : LVT_logunit
  use OzFlux_obsMod,      only : ozfluxobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads estimates of 
!  sensible,latent and ground heat flux
!  measurements from the Ozflux sites. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[source]  datastream source
!   \item[yr]      year for which the data is being processed
!   \item[mo]      month for which the data is being processed
!   \item[da]      day for which the data is being processed
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer          :: source
  integer          :: yr
  integer          :: mo
  integer          :: da
!
!EOP
  integer          :: i 
  character*100    :: filename
  integer          :: status
  logical          :: file_exists

  do i=1,ozfluxobs(source)%n_stns
     call create_ozflux_filename(ozfluxobs(source)%odir, &
          ozfluxobs(source)%stn_name(i),yr, filename)
     inquire(file=filename,exist=file_exists)
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading OzFlux file   ',trim(filename)
        call read_ozflux_file(source, i,yr,mo,da,filename)
     endif
  enddo

end subroutine process_ozflux_data

!BOP
! 
! !ROUTINE: read_ozflux_file
! \label{read_ozflux_file}
!
! !INTERFACE: 
subroutine read_ozflux_file(source, k, yr,mo,da,filename)
! 
! !USES: 
  use ESMF     
  use LVT_coreMod,      only : LVT_rc
  use LVT_logMod,       only : LVT_logunit, LVT_verify
  use LVT_timeMgrMod,   only : LVT_calendar
  use OzFlux_obsMod,      only : ozfluxobs
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
!
!DESCRIPTION: 
!  This routine reads estimates of 
!  latent, sensible and ground heat fluxes from the OzFlux file, 
!  for a particular station. The data is read from the 
!  native NETCDF files. This routine 
!  also computes the temporal offset of the data relative
!  to the 0z of a particular day. 
!
!  The arguments are: 
!  \begin{description}
!   \item[k]         Index of the OzFlux station
!   \item[yr]        year of OzFlux data 
!   \item[mo]        month of OzFlux data
!   \item[da]        day of OzFlux data 
!   \item[filename]  filename for the ozflux data
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer,  intent(in) :: k 
  integer,  intent(in) :: source
  integer              :: yr
  integer              :: mo
  integer              :: da
  character(len=*)     :: filename 
  integer              :: bend
!EOP
  integer       :: nid, btimeid, latid, lonid, qleid
  integer       :: qhid,qgid, timeid, dimId
  integer       :: ndims
  real          :: lat, lon 
  integer       :: ios
  real, allocatable :: time(:)
  real, allocatable :: qh(:,:,:)
  real, allocatable :: qle(:,:,:)
  real, allocatable :: qg(:,:,:)
  integer           :: kk
  type(ESMF_Time)   :: reftime, datatime!, currtime
  type(ESMF_TimeInterval) :: dt
  integer           :: yr1, mo1, da1, hr1, mn1, ss1
  integer           :: status
  integer           :: data_index

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
  call LVT_verify(ios, 'Error opening file'//trim(filename))

  !variable ids
  ios = nf90_inq_varid(nid,'time',timeid)
  call LVT_verify(ios,'Error in nf90_inq_varid: time')
  
  ios = nf90_inq_varid(nid, 'latitude',latid)
  call LVT_verify(ios, 'Error nf90_inq_varid: latitude')

  ios = nf90_inq_varid(nid, 'longitude',lonid)
  call LVT_verify(ios, 'Error nf90_inq_varid: longitude')

  ios = nf90_inq_varid(nid, 'Fh',qhid)
  call LVT_verify(ios, 'Error nf90_inq_varid: Fh')

  ios = nf90_inq_varid(nid, 'Fe',qleid)
  call LVT_verify(ios, 'Error nf90_inq_varid: qle')

  ios = nf90_inq_varid(nid, 'Fg',qgid)
  call LVT_verify(ios,'Error in nf90_inq_varid: Fg')

!dimensions
  ios = nf90_inq_dimid(nid, 'time',dimId)
  call LVT_verify(ios, 'Error nf90_inq_dimid: time')

  ios = nf90_inquire_dimension(nid, dimId, len=ndims)
  call LVT_verify(ios, 'Error nf90_inquire_dimension:')

!values
  ios = nf90_get_var(nid,latid, lat)
  call LVT_verify(ios, 'Error nf90_get_var: latitude')

  ios = nf90_get_var(nid,lonid, lon)
  call LVT_verify(ios, 'Error nf90_get_var: longitude')

  call ESMF_TimeSet(refTime, yy=yr,mm=mo,dd=da,&
       h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
  call LVT_verify(status, 'error in timeset: readOzFluxObs')

  allocate(time(ndims))
  allocate(qle(1,1,ndims))
  allocate(qh(1,1,ndims))
  allocate(qg(1,1,ndims))

  ios = nf90_get_var(nid,timeid, time)
  call LVT_verify(ios, 'Error nf90_get_var: time')

  ios = nf90_get_var(nid,qhid, qh)
  call LVT_verify(ios, 'Error nf90_get_var: qh')

  ios = nf90_get_var(nid,qleid, qle)
  call LVT_verify(ios, 'Error nf90_get_var: qle')

  ios = nf90_get_var(nid,qgid, qg)
  call LVT_verify(ios, 'Error nf90_get_var: qg')

  ios = nf90_close(nid)
  call LVT_verify(ios, 'Error in nf90_close')
  
  call ESMF_TimeIntervalSet(dt,s=1800,rc=status)
  call LVT_verify(status, 'Error in timeintervalset: readOzFluxobs')

!  reftime = reftime + dt
  datatime = reftime
  
  do kk=1,ndims
!     call ESMF_TimeIntervalSet(dt,s=nint(time(kk)*86400.0),rc=status)
!     call LVT_verify(status, 'Error in timeintervalset: readOzFluxobs')

     
     call ESMF_TimeGet(datatime, yy=yr1,mm=mo1,dd=da1,&
          h=hr1,m=mn1,s=ss1,calendar=LVT_calendar,rc=status)
     call LVT_verify(status, 'error in timeget: readOzFluxObs')

!     call ESMF_TimeSet(currTime, yy=yr, mm=mo, &
!          dd=da, h=0,m=0,s=0,calendar=LVT_calendar,rc=status)
!     call LVT_verify(status, 'error in timeget: readOzFluxObs')

     data_index = nint((datatime - ozfluxobs(source)%starttime)/&
          ozfluxobs(source)%ts) + 1

     ozfluxobs(source)%tindex(k,data_index) = data_index

     ozfluxobs(source)%qh(k,data_index)     = qh(1,1,kk)
     ozfluxobs(source)%qle(k,data_index)    = qle(1,1,kk)
     ozfluxobs(source)%qg(k,data_index)     = qg(1,1,kk)

     datatime = datatime + dt
  enddo

  deallocate(time)
  deallocate(qle)
  deallocate(qh)
  deallocate(qg)
#endif

end subroutine read_ozflux_file


!BOP
! 
! !ROUTINE: create_ozflux_filename
! \label{create_ozflux_filename}
!
! !INTERFACE: 
subroutine create_ozflux_filename(odir, site_id, yr, filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This routine creates a filename for OzFlux in-situ data files. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      OzFlux base directory
!   \item[site_id]   OzFlux site identifier
!   \item[yr]        year of data
!   \item[filename]  Name of the flux file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: site_id
  integer,          intent(in)  :: yr
  character(len=*)              :: filename

!EOP

  character*4       :: fyr

  write(unit=fyr, fmt='(i4.4)') yr

  filename = trim(odir)//'/'//trim(site_id)//&
       '_'//trim(fyr)//'_L3.nc'

end subroutine create_ozflux_filename


