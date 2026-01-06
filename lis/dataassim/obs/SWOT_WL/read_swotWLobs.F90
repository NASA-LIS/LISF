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
#include "LIS_NetCDF_inc.h"
!BOP
! !ROUTINE: read_swotWLobs
! \label{read_swotWLobs}
!
! !REVISION HISTORY:
!  17 Apr 2024: Yeosang Yoon; Initial Specification
!  25 Nov 2025: Yeosang Yoon; Clean up the code
!
! !INTERFACE:
subroutine read_swotWLobs(n, k, OBS_State, OBS_Pert_state)
! !USES:
  use ESMF
  use LIS_mpiMod
  use LIS_historyMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use LIS_dataAssimMod
  use LIS_surfaceModelMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_timeMgrMod
  use swotWLobs_module

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!
!  reads the swot radar altimetry observations
!  The processed data is packaged
!  into an ESMF State for later use within the DA algorithm.
!
!  The arguments are:
!  \begin{description}
!  \item[n]                index of the nest
!  \item[k]                index of the data assimilation instance
!  \item[OBS\_State]       observations state
!  \item[OBS\_Pert\_State] observations perturbation state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: wlfield

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real                :: obs_unsc(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: filename, wlobsdir
  logical             :: data_update
  logical             :: alarmCheck, file_exists, file_exists1
  integer             :: yyyy, mm, dd, hh, mn, ss
  type(ESMF_Time)     :: currTime
  integer             :: fnd
  integer             :: status
  logical             :: data_upd_flag(LIS_npes)
  logical             :: data_upd_flag_local
  logical             :: data_upd
  real                :: wl_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                :: wlobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer             :: t,c,r
  integer             :: time_now, time_start, time_end, dt

  obs_unsc = LIS_rc%udef
  call ESMF_AttributeGet(OBS_State,"Data Directory",wlobsdir,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",data_update,rc=status)
  call LIS_verify(status)

   file_exists = .false.
   file_exists1 = .false.
   alarmCheck = LIS_isAlarmRinging(LIS_rc, "swotWL read alarm")

   ! check current time
   call ESMF_ClockGet(LIS_clock, currTime=currTime, rc=status)
   call ESMF_TimeGet(currTime, yy=yyyy, mm=mm, dd=dd, h=hh, m=mn, s=ss, rc=status)

   if (alarmCheck) then
      call swotWL_filename(filename, wlobsdir, yyyy, mm, dd)

      inquire(file=trim(filename), exist=file_exists)

      if(file_exists) then
         write(LIS_logunit,*)  '[INFO] Reading SWOT WL data ',trim(filename)
         call read_swotWL(filename)
         swot_wl_struc(n)%time_file = &
              calculate_seconds_since_2000(yyyy, mm, dd, hh, mn, ss) !since 2000-01-01
      else
         write(LIS_logunit,*)'[WARN] Cannot find file ',trim(filename)
      end if
   end if

  time_now = calculate_seconds_since_2000(yyyy, mm, dd, hh, mn, ss) !since 2000-01-01
  dt = nint(swot_wl_struc(n)%daInterval)

  if((time_now <= swot_wl_struc(n)%time_file+86400) .and. &
       (mod(time_now, dt).eq.0)) then
     call ESMF_StateGet(OBS_State,"Observation01",wlfield,rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(wlfield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)

     obsl = -1
     wlobs = -9999.0
     fnd = 0

     ! extract SWOT observations from the daily aggregated file for assimilation at the cycle
     time_start = time_now
     time_end = time_now + dt

     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if(swot_wl_struc(n)%lisid(c,r).gt.0) then
              do t=1, swot_wl_struc(n)%ntimes
                 if((swot_wl_struc(n)%time(swot_wl_struc(n)%lisid(c,r),t) >= time_start) &
                    .and. (swot_wl_struc(n)%time(swot_wl_struc(n)%lisid(c,r),t) <= time_end)) then
                    wlobs(c,r) = swot_wl_struc(n)%WLobs(swot_wl_struc(n)%lisid(c,r),t)

                    if(wlobs(c,r).gt.0) then
                       fnd = 1
                    endif
                 endif
              enddo
           endif
        enddo
     enddo

     !TODO: "Normal deviate scaling" is not yet working.
     !Set dascaleoption to 'none'.
     if(LIS_rc%dascaloption(k).eq."Normal deviate scaling".and. &
          fnd.ne.0) then

        call LIS_rescale_with_normal_deviate_scaling(&
             n,k,                         &
             swot_wl_struc(n)%nt,         &
             swot_wl_struc(n)%model_mu,   &
             swot_wl_struc(n)%model_sigma,&
             swot_wl_struc(n)%obs_mu,     &
             swot_wl_struc(n)%obs_sigma,  &
             wlobs)
     endif

     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = wlobs(c,r)
           end if
        end do
     end do

     !  Apply LSM-based QC of observations
     call LIS_checkForValidObs(n, k,obsl,fnd,wl_current)
     call LIS_surfaceModel_DAqcObsState(n,k)
     call LIS_checkForValidObs(n, k,obsl,fnd,wl_current)

     if(fnd.eq.0) then
        data_upd_flag_local = .false.
     else
        data_upd_flag_local = .true.
     endif

#if (defined SPMD)
     call MPI_ALLGATHER(data_upd_flag_local, 1, MPI_LOGICAL, data_upd_flag(:),&
          1, MPI_LOGICAL, LIS_mpi_comm, status)
     data_upd = any(data_upd_flag)
#else
     data_upd = data_upd_flag_local
#endif

     if(data_upd) then
        do t=1,LIS_rc%obs_ngrid(k)
           gid(t) = t
           if(obsl(t).ne.-9999.0) then
              assimflag(t) = 1
              fnd = 1
           else
              assimflag(t) = 0
           endif
        enddo

        call ESMF_AttributeSet(OBS_State,"Data Update Status",.true., rc=status)
        call LIS_verify(status)

        if(LIS_rc%obs_ngrid(k).gt.0) then
           call ESMF_AttributeSet(wlfield,"Grid Number",&
                gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)

           call ESMF_AttributeSet(wlfield,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)

           call ESMF_AttributeSet(wlfield, "Unscaled Obs",&
                obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
           call LIS_verify(status, 'Error in setting Unscaled Obs attribute')
        endif
     else
        call ESMF_AttributeSet(OBS_State,"Data Update Status",.false., rc=status)
        call LIS_verify(status)
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",.false., rc=status)
     call LIS_verify(status)
     return
  endif

contains

   ! Constructs swotWL filename
   subroutine swotWL_filename(filename, dir, yyyy, mm, dd)
      implicit none
      character(len=LIS_CONST_PATH_LEN), intent(inout) :: filename
      character(len=LIS_CONST_PATH_LEN), intent(in)    :: dir
      integer,                           intent(in)    :: yyyy
      integer,                           intent(in)    :: mm
      integer,                           intent(in)    :: dd
      character(4)                                     :: cyyyy
      character(2)                                     :: cmm, cdd

      write(unit=cyyyy, fmt='(i4.4)') yyyy
      write(unit=cmm, fmt='(i2.2)') mm
      write(unit=cdd, fmt='(i2.2)') dd

      filename = trim(dir) // "/SWOT_L2_HR_RiverSP_Node_" &
           // trim(cyyyy) // trim(cmm) // trim(cdd) // ".d01.nc"
   end subroutine swotWL_filename

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   subroutine read_swotWL(filename)
     use netcdf

     implicit none
     character(len=LIS_CONST_PATH_LEN), intent(in) :: filename
     integer                                       :: dim1Id, dim2Id, tId, wlId
     integer                                       :: ftn

     call LIS_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ftn),&
          'Error opening file '//trim(filename))
     call LIS_verify(nf90_inq_dimid(ftn,'times',dim1Id),&
          'Error with nf90_inq_dimid: times')
     call LIS_verify(nf90_inq_dimid(ftn,'lis_grids',dim2Id),&
          'Error with nf90_inq_dimid: lis_grids')

     call LIS_verify(nf90_inquire_dimension(ftn,dim1Id,&
          len=swot_wl_struc(n)%ntimes),'Error with nf90_inquire_dimension: times')
     call LIS_verify(nf90_inquire_dimension(ftn,dim2Id,&
          len=swot_wl_struc(n)%nlisgrids),'Error with nf90_inquire_dimension: lis_grids')

     if(allocated(swot_wl_struc(n)%time)) deallocate(swot_wl_struc(n)%time)
     if(allocated(swot_wl_struc(n)%WLobs)) deallocate(swot_wl_struc(n)%WLobs)

     allocate(swot_wl_struc(n)%time(swot_wl_struc(n)%nlisgrids,swot_wl_struc(n)%ntimes))
     allocate(swot_wl_struc(n)%WLobs(swot_wl_struc(n)%nlisgrids,swot_wl_struc(n)%ntimes))

     call LIS_verify(nf90_inq_varid(ftn,'time',tId),'Error with nf90_iniq_varid: time')
     call LIS_verify(nf90_get_var(ftn,tId,swot_wl_struc(n)%time),'Error with nf90_get_var: time')

     call LIS_verify(nf90_inq_varid(ftn,'wse',wlId),'Error with nf90_inq_varid: wse')
     call LIS_verify(nf90_get_var(ftn,wlId,swot_wl_struc(n)%WLobs),'Error with nf90_get_var: wse')

     call LIS_verify(nf90_close(ftn))
   end subroutine read_swotWL
#endif

   function calculate_seconds_since_2000(year, month, day, hr, mn, ss) result(seconds)
     implicit none
     integer :: year, month, day, hr, mn, ss
     integer :: seconds, days_since_2000
     integer, parameter :: days_per_year = 365, days_per_leap_year = 366
     integer :: y, m

     days_since_2000 = 0
     do y = 2000, year - 1
         if (is_leap_year(y)) then
             days_since_2000 = days_since_2000 + days_per_leap_year
         else
             days_since_2000 = days_since_2000 + days_per_year
         end if
     end do

     do m = 1, month - 1
         days_since_2000 = days_since_2000 + days_in_month(m, year)
     end do

     days_since_2000 = days_since_2000 + (day - 1)

     seconds = days_since_2000 * 86400 + hr * 3600 + mn * 60 + ss
   end function calculate_seconds_since_2000

   function is_leap_year(year) result(leap)
     implicit none
     integer :: year
     logical :: leap
     leap = .false.
     if (mod(year, 4) == 0) then
         if (mod(year, 100) /= 0 .or. mod(year, 400) == 0) then
             leap = .true.
         end if
     end if
   end function is_leap_year

   function days_in_month(month, year) result(days)
     implicit none
     integer :: month, year
     integer :: days
     integer, dimension(12) :: month_days = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

     days = month_days(month)
     if (month == 2 .and. is_leap_year(year)) then
         days = 29
     end if
   end function days_in_month

end subroutine read_swotWLobs

