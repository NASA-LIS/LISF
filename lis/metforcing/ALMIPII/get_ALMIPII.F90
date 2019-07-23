!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: get_ALMIPII
! \label{get_ALMIPII}
!
!
! !REVISION HISTORY:
!  14 Oct 2010: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine get_ALMIPII(n, findex)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf  
#endif
  use ESMF
  use LIS_coreMod, only        : LIS_rc, LIS_domain
  use LIS_timeMgrMod,     only : LIS_calendar
  use LIS_metforcingMod,  only : LIS_forc
  use LIS_logMod,         only : LIS_logunit, LIS_endrun, LIS_verify
  use ALMIPII_forcingMod,   only : ALMIPII_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!
!EOP

  character*4           :: fyr
  character*100         :: ALMIPIIfile
  integer               :: lonid, latid
  integer               :: tair_id, qair_id, swd_id, lwd_id, u_id,v_id
  integer               :: psurf_id, rainf_id, snowf_id
  integer               :: tindex1, tindex2
  real                  :: fvar(ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr)
  real                  :: varfield(LIS_rc%ngrid(n))
  type(ESMF_Time)       :: currTime
  integer               :: status
  
  if(ALMIPII_struc(n)%startFlag.or.(ALMIPII_struc(n)%syr.ne.LIS_rc%yr)) then 

     if((.not.ALMIPII_struc(n)%startFlag).and.&
          (ALMIPII_struc(n)%syr.ne.LIS_rc%yr)) then 
!yr switch during the program
!first close the existing file
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        status = nf90_close(ncid=ALMIPII_struc(n)%ncid)
        call LIS_verify(status, 'Problem closing the ALMIPII file')
#endif     
     endif
     ALMIPII_struc(n)%startFlag = .false.
     ALMIPII_struc(n)%syr = LIS_rc%yr
     write(unit=fyr,fmt='(i4.4)') LIS_rc%yr
     ! time to open a new file. 
     ALMIPIIfile = trim(ALMIPII_struc(n)%dir)//'/'//&
          trim(ALMIPII_struc(n)%filename_prefix)//'_'//trim(fyr)//'.nc'
     
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     write(LIS_logunit,*) 'Opening ALMIPII file ',trim(ALMIPIIfile)
     status = nf90_open(ALMIPIIfile,mode=nf90_nowrite,&
          ncid=ALMIPII_struc(n)%ncid)
     call LIS_verify(status, 'Problem opening the ALMIPII file')
#endif

!   read dimensions
     call LIS_verify(nf90_inq_dimId(ALMIPII_struc(n)%ncid,'lon',lonid),&
          'Problem with nf90_inq_dimId call for lon')
     call LIS_verify(nf90_inquire_dimension(ALMIPII_struc(n)%ncid,lonid,&
          len=ALMIPII_struc(n)%nc),&
          'Problem with nf90_inquire_dimension for nc')
     call LIS_verify(nf90_inq_dimId(ALMIPII_struc(n)%ncid,'lat',latid),&
          'Problem with nf90_inq_dimId call for lat')
     call LIS_verify(nf90_inquire_dimension(ALMIPII_struc(n)%ncid,latid,&
          len=ALMIPII_struc(n)%nr),&
          'Problem with nf90_inquire_dimension for nr')

     call ESMF_TimeSet(ALMIPII_struc(n)%startTime, yy = LIS_rc%yr, &
          mm = 1, &
          dd = 1, &
          h  = 0, &
          m  = 0, & 
          s  = 0, &
          calendar = LIS_calendar, & 
          rc = status)

  endif

! Once the file is opened, we simply index into the two required bookends
! Note that the logic expects the model timestep to be always less
! than or equal to the forcing timestep (which is currently at 30min)

  if(LIS_rc%ss.eq.0) then !only reads to read at some multiple of minutes
     call ESMF_TimeSet(ALMIPII_struc(n)%time1, yy = LIS_rc%yr, &
          mm = LIS_rc%mo, &
          dd = LIS_rc%da, &
          h  = LIS_rc%hr, &
          m  = LIS_rc%mn, & 
          s  = LIS_rc%ss, &
          calendar = LIS_calendar, & 
          rc = status)
     call LIS_verify(status, 'Error in ESMF_TimeSet')
     
     tindex1 = nint((ALMIPII_struc(n)%time1 - ALMIPII_struc(n)%startTime)/&
          ALMIPII_struc(n)%timestep) +1
     tindex2 = tindex1 + 1
     ALMIPII_struc(n)%time2 = ALMIPII_struc(n)%time2 + &
          ALMIPII_struc(n)%timestep 

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

!Tair
     status = nf90_inq_varid(ALMIPII_struc(n)%ncid, "Tair", tair_id)
     call LIS_verify(status, 'Error: nf90_inq_varid for Tair in ALMIPII')

     status = nf90_get_var(ALMIPII_struc(n)%ncid, tair_id,fvar,&
          start = (/1,1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Tair in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)

     ALMIPII_struc(n)%metdata1(1,:) = varfield(:)

     status = nf90_get_var(ALMIPII_struc(n)%ncid, tair_id,fvar,&
          start = (/1,1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Tair in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)
     ALMIPII_struc(n)%metdata2(1,:) = varfield(:)

!Qair
     status = nf90_inq_varid(ALMIPII_struc(n)%ncid, "Qair", qair_id)
     call LIS_verify(status, 'Error: nf90_inq_varid for Qair in ALMIPII')

     status = nf90_get_var(ALMIPII_struc(n)%ncid, qair_id,fvar,&
          start = (/1,1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Qair in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)

     ALMIPII_struc(n)%metdata1(2,:) = varfield(:)

     status = nf90_get_var(ALMIPII_struc(n)%ncid, qair_id,fvar,&
          start = (/1,1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Qair in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)
     ALMIPII_struc(n)%metdata2(2,:) = varfield(:)

!SWdown
     status = nf90_inq_varid(ALMIPII_struc(n)%ncid, "SWdown", swd_id)
     call LIS_verify(status, 'Error: nf90_inq_varid for SWdown in ALMIPII')

     status = nf90_get_var(ALMIPII_struc(n)%ncid, swd_id,fvar,&
          start = (/1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1/))
     call LIS_verify(status, 'Error: nf90_get_var for SWdown in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)

     ALMIPII_struc(n)%metdata1(3,:) = varfield(:)

     status = nf90_get_var(ALMIPII_struc(n)%ncid, swd_id,fvar,&
          start = (/1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1/))
     call LIS_verify(status, 'Error: nf90_get_var for SWdown in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)
     ALMIPII_struc(n)%metdata2(3,:) = varfield(:)

!LWdown
     status = nf90_inq_varid(ALMIPII_struc(n)%ncid, "LWdown", lwd_id)
     call LIS_verify(status, 'Error: nf90_inq_varid for LWdown in ALMIPII')

     status = nf90_get_var(ALMIPII_struc(n)%ncid, lwd_id,fvar,&
          start = (/1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1/))
     call LIS_verify(status, 'Error: nf90_get_var for LWdown in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)

     ALMIPII_struc(n)%metdata1(4,:) = varfield(:)

     status = nf90_get_var(ALMIPII_struc(n)%ncid, lwd_id,fvar,&
          start = (/1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1/))
     call LIS_verify(status, 'Error: nf90_get_var for LWdown in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)
     ALMIPII_struc(n)%metdata2(4,:) = varfield(:)

!Wind_N
     status = nf90_inq_varid(ALMIPII_struc(n)%ncid, "Wind_N", u_id)
     call LIS_verify(status, 'Error: nf90_inq_varid for Wind_N in ALMIPII')

     status = nf90_get_var(ALMIPII_struc(n)%ncid, u_id,fvar,&
          start = (/1,1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Wind_N in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)

     ALMIPII_struc(n)%metdata1(5,:) = varfield(:)

     status = nf90_get_var(ALMIPII_struc(n)%ncid, u_id,fvar,&
          start = (/1,1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Wind_N in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)
     ALMIPII_struc(n)%metdata2(5,:) = varfield(:)

!Wind_E
     status = nf90_inq_varid(ALMIPII_struc(n)%ncid, "Wind_E", v_id)
     call LIS_verify(status, 'Error: nf90_inq_varid for Wind_E in ALMIPII')

     status = nf90_get_var(ALMIPII_struc(n)%ncid, v_id,fvar,&
          start = (/1,1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Wind_E in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)

     ALMIPII_struc(n)%metdata1(6,:) = varfield(:)

     status = nf90_get_var(ALMIPII_struc(n)%ncid, v_id,fvar,&
          start = (/1,1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Wind_E in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)
     ALMIPII_struc(n)%metdata2(6,:) = varfield(:)

!Psurf
     status = nf90_inq_varid(ALMIPII_struc(n)%ncid, "PSurf", psurf_id)
     call LIS_verify(status, 'Error: nf90_inq_varid for PSurf in ALMIPII')

     status = nf90_get_var(ALMIPII_struc(n)%ncid, psurf_id,fvar,&
          start = (/1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1/))
     call LIS_verify(status, 'Error: nf90_get_var for PSurf in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)

     ALMIPII_struc(n)%metdata1(7,:) = varfield(:)

     status = nf90_get_var(ALMIPII_struc(n)%ncid, psurf_id,fvar,&
          start = (/1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1/))
     call LIS_verify(status, 'Error: nf90_get_var for PSurf in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)
     ALMIPII_struc(n)%metdata2(7,:) = varfield(:)

!Rainf
     status = nf90_inq_varid(ALMIPII_struc(n)%ncid, "Rainf", rainf_id)
     call LIS_verify(status, 'Error: nf90_inq_varid for Rainf in ALMIPII')

     status = nf90_get_var(ALMIPII_struc(n)%ncid, rainf_id,fvar,&
          start = (/1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Rainf in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)

     ALMIPII_struc(n)%metdata1(8,:) = varfield(:)

     status = nf90_get_var(ALMIPII_struc(n)%ncid, rainf_id,fvar,&
          start = (/1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Rainf in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)
     ALMIPII_struc(n)%metdata2(8,:) = varfield(:)

!Snowf
     status = nf90_inq_varid(ALMIPII_struc(n)%ncid, "Snowf", snowf_id)
     call LIS_verify(status, 'Error: nf90_inq_varid for Snowf in ALMIPII')

     status = nf90_get_var(ALMIPII_struc(n)%ncid, snowf_id,fvar,&
          start = (/1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Snowf in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)

     ALMIPII_struc(n)%metdata1(9,:) = varfield(:)

     status = nf90_get_var(ALMIPII_struc(n)%ncid, snowf_id,fvar,&
          start = (/1,1,tindex1/),&
          count = (/ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr,1/))
     call LIS_verify(status, 'Error: nf90_get_var for Snowf in ALMIPII')

     call interp_ALMIPII(n, findex, fvar,varfield, .false.)
     ALMIPII_struc(n)%metdata2(9,:) = varfield(:)
#endif 
     
  endif

end subroutine get_ALMIPII







































