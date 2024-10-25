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
! !ROUTINE: read_GlobalLSObsdata
! \label{read_GlobalLSObsdata}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_GlobalLSObsdata(Obj_Space)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod,  only : LIS_rc, LIS_domain, LIS_npes, LIS_localPet
  use LIS_timeMgrMod, only : LIS_calendar
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use LIS_fileIOMod,      only : LIS_readData
  use GlobalLSDataMod, only : globallsobs_struc
  use map_utils

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP
  real,    allocatable    :: lsobs(:,:)
  real,    pointer    :: obsl(:)
  type(ESMF_Field)    :: lsField
  character(len=LIS_CONST_PATH_LEN) :: lsobsdir, name
  logical             :: data_update
  logical             :: file_exists
  logical             :: readflag
  integer             :: status, ierr
  integer             :: fnd
  integer             :: ftn
  integer             :: c,r
  integer             :: istat
  type(ESMF_TimeInterval) :: dt,oneday
  type(ESMF_Time)         :: time1
  integer             :: interval
  integer             :: t
  real                :: col,row
  integer             :: stn_col,stn_row
  logical             :: data_upd, data_upd_flag(LIS_npes)
  integer             :: n,p

  n = 1

  call ESMF_AttributeGet(Obj_Space,"Data Directory",&
       lsobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  
  allocate(lsobs(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  lsobs = 0
  fnd = 1
  do t=1,globallsobs_struc(n)%size
     if(.not.globallsobs_struc(n)%flag(t)) then 
        call ESMF_TimeSet(time1, yy=LIS_rc%yr, mm=LIS_rc%mo, &
             dd=LIS_rc%da,h=LIS_rc%hr,&
             m=LIS_rc%mn,calendar=LIS_calendar,rc=status)
        call ESMF_TimeIntervalSet(oneday,s=86400,rc=status)
        dt = time1 - (globallsobs_struc(n)%time(t)+oneday)
        call ESMF_TimeIntervalGet(dt,s=interval,rc=status)
!        if(abs(interval).le.86400) then 
        if(interval.le.43200.and.interval.ge.0) then !within a 12 hr window 
!for now assume that all obs correspond to the modeling point
#if 0            
           call latlon_to_ij(LIS_domain(n)%lisproj,globallsobs_struc(n)%lat(t),&
                globallsobs_struc(n)%lon(t),col,row)
           stn_col = nint(col)
           stn_row = nint(row)
           globallsobs_struc(n)%flag(t) = .true. 
           if((stn_col.ge.1.and.stn_col.le.LIS_rc%lnc(n)).and.&
                (stn_row.ge.1.and.stn_row.le.LIS_rc%lnr(n))) then 
              lsobs(stn_col,stn_row) = 1.0
           endif
#endif
           globallsobs_struc(n)%flag(t) = .true. 
!           print*, 'found obs ',interval,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr
           lsobs = 1.0
           fnd = 1
        endif
     endif
  enddo

  call ESMF_StateGet(Obj_Space,"Landslide data",lsField,&
       rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(lsField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           obsl(LIS_domain(n)%gindex(c,r)) = lsobs(c,r)
        endif
     enddo
  enddo

  if(fnd.eq.1) then      
     write(LIS_logunit,*) 'Read landslide data'
     call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status)
  else     
     call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
  endif

  deallocate(lsobs)

end subroutine read_GlobalLSObsdata

  
  

  

