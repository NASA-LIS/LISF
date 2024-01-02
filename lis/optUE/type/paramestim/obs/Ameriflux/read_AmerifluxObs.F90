!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_AmerifluxObs
! \label{read_AmerifluxObs}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_AmerifluxObs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_timeMgrMod, only : LIS_calendar, LIS_tick
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_fileIOMod,      only : LIS_readData
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use AmerifluxobsMod, only : AmerifluxObs_struc
  use map_utils,       only : latlon_to_ij

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP
  real,    pointer    :: qle(:)
  real,    pointer    :: qh(:)
  real,    pointer    :: qg(:)
  real,    allocatable    :: sm1(:)
  real,    allocatable    :: sfsm(:)
  real,    allocatable    :: st1(:)
  real,    allocatable    :: sfst(:)
  integer,    allocatable    :: nqle(:)
  integer,    allocatable    :: nqh(:)
  integer,    allocatable    :: nqg(:)
  integer,    allocatable    :: nsm1(:)
  integer,    allocatable    :: nst1(:)

  type(ESMF_Field)    :: qleField
  type(ESMF_Field)    :: qhField
  type(ESMF_Field)    :: qgField
  type(ESMF_Field)    :: sfstField
  type(ESMF_Field)    :: sfsmField
  character(len=LIS_CONST_PATH_LEN) :: obsdir, name
  logical             :: data_update
  logical             :: file_exists
  logical             :: readflag
  integer             :: status
  integer             :: ftn
  integer             :: i
  type(ESMF_Time)     :: amerifluxtime1, amerifluxtime2
  integer             :: t, st, et
  integer             :: yr, mo, da, hr, mn, ss
  real*8              :: lis_prevtime
  integer             :: c,r, gid, stn_row, stn_col
  real                :: col, row
  real                :: gmt
  integer             :: doy
  integer             :: istat
  real                :: gridDesc(6)
  integer             :: n 

  n = 1

  allocate(nqle(LIS_rc%ngrid(n)))
  allocate(nqh(LIS_rc%ngrid(n)))
  allocate(nqg(LIS_rc%ngrid(n)))
  allocate(nsm1(LIS_rc%ngrid(n)))
  allocate(nst1(LIS_rc%ngrid(n)))

  call ESMF_AttributeGet(Obj_Space,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  if(Amerifluxobs_struc(n)%startFlag.or.(Amerifluxobs_struc(n)%yr.ne.LIS_rc%yr)) then 

     call ESMF_TimeSet(AmerifluxObs_struc(n)%starttime, yy=LIS_rc%yr, &
          mm=1, dd=1, h=0, m=0,s = 0, calendar=LIS_calendar, rc=status)
     call LIS_verify(status, 'ameriflux starttime set failed')

     Amerifluxobs_struc(n)%startflag = .false. 
     Amerifluxobs_struc(n)%yr = LIS_rc%yr

     do i = 1,Amerifluxobs_struc(n)%n_stns
        call read_ameriflux_station(n,i)
     end do
  endif

  call ESMF_TimeSet(amerifluxtime1, yy=LIS_rc%yr, &
       mm=LIS_rc%mo, dd=LIS_rc%da, h=LIS_rc%hr, m=LIS_rc%mn, &
       s = LIS_rc%ss, calendar=LIS_calendar, rc=status)
  call LIS_verify(status, 'amerifluxtime1 set failed')

  et = nint((amerifluxtime1 - amerifluxobs_struc(n)%starttime)/amerifluxobs_struc(n)%timestep)+1

  yr = LIS_rc%yr
  mo = LIS_rc%mo
  da = LIS_rc%da
  hr = LIS_rc%hr
  mn = LIS_rc%mn
  ss = LIS_rc%ss

  ! go back one day
  call LIS_tick(lis_prevtime, doy, gmt, yr,mo,da,hr,mn,ss,(-1)*LIS_rc%ts)

  call ESMF_TimeSet(amerifluxtime2, yy=yr, &
       mm = mo, &
       dd = da, &
       h = hr, &
       m = mn, &
       calendar = LIS_calendar, &
       rc=status)
  call LIS_verify(status)

  st = nint((amerifluxtime2 - amerifluxobs_struc(n)%starttime)/amerifluxobs_struc(n)%timestep)+1
  !  print *, 'ending index, starting index: ', et, st

  call ESMF_StateGet(Obj_Space,"Ameriflux Qle",qleField,&
       rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(qleField,localDE=0,farrayPtr=qle,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(Obj_Space,"Ameriflux Qh",qhField,&
       rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(qhField,localDE=0,farrayPtr=qh,rc=status)
  call LIS_verify(status)  

  call ESMF_StateGet(Obj_Space,"Ameriflux Qg",qgField,&
       rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(qgField,localDE=0,farrayPtr=qg,rc=status)
  call LIS_verify(status)  

!  call ESMF_StateGet(Obj_Space,"Ameriflux sfsm",sfsmField,&
!       rc=status)
!  call LIS_verify(status)
  
!  call ESMF_FieldGet(sfsmField,localDE=0,farrayPtr=sfsm,rc=status)
!  call LIS_verify(status)  

!  call ESMF_StateGet(Obj_Space,"Ameriflux sfst",sfstField,&
!       rc=status)
!  call LIS_verify(status)
  
!  call ESMF_FieldGet(sfstField,localDE=0,farrayPtr=sfst,rc=status)
!  call LIS_verify(status)  

  if(st.gt.0.and.et.gt.0) then 
     do t = st+1,et !to avoid double counting from last cycle
        do i=1,amerifluxobs_struc(n)%n_stns
           call latlon_to_ij(LIS_domain(n)%lisproj, amerifluxobs_struc(n)%stnlat(i), &
                amerifluxobs_struc(n)%stnlon(i), col, row)
           stn_col = nint(col)
           stn_row = nint(row)
           gid = LIS_domain(n)%gindex(stn_col,stn_row)

           if((stn_col.ge.1.and.stn_col.le.LIS_rc%lnc(n)).and.&
                (stn_row.ge.1.and.stn_row.le.LIS_rc%lnr(n).and.&
                gid.ne.-1)) then 
              if(amerifluxobs_struc(n)%qle(i,t).ne.LIS_rc%udef) then 
                 qle(gid) = qle(gid) +&
                      amerifluxobs_struc(n)%qle(i,t)
                 nqle(gid) = nqle(gid) + 1          
              end if

              if(amerifluxobs_struc(n)%qh(i,t).ne.LIS_rc%udef) then 
                 qh(gid) = qh(gid) + &
                      amerifluxobs_struc(n)%qh(i,t)
                 nqh(gid) = nqh(gid) + 1          
              end if

              if(amerifluxobs_struc(n)%qg(i, t).ne.LIS_rc%udef) then
                 qg(gid) = qg(gid) + &
                      amerifluxobs_struc(n)%qg(i, t)
                 nqg(gid) = nqg(gid) + 1
              end if

              if(amerifluxobs_struc(n)%sfsm(i, t).ne.LIS_rc%udef) then
                 sm1(gid) = sm1(gid) + &
                      amerifluxobs_struc(n)%sfsm(i, t)
                 nsm1(gid) = nsm1(gid) + 1
              end if

              if(amerifluxobs_struc(n)%sfst(i, t).ne.LIS_rc%udef) then
                 st1(gid) = st1(gid) + &
                      amerifluxobs_struc(n)%sfst(i, t)
                 nst1(gid) = nst1(gid) + 1
              end if
           endif
        end do
     end do
  endif


  do gid=1,LIS_rc%ngrid(n)
     if(nqle(gid).gt.0) then 
        qle(gid) = qle(gid)/nqle(gid)
     else
        qle(gid) = LIS_rc%udef
     end if
     
     if(nqh(gid).gt.0) then 
        qh(gid) = qh(gid)/nqh(gid)
     else
        qh(gid) = LIS_rc%udef
     end if
     
     if(nqg(gid).gt.0) then 
        qg(gid) = qg(gid)/nqg(gid)
     else
        qg(gid) = LIS_rc%udef
     end if
     
!     if(nsm1(gid).gt.0) then 
!        sm1(gid) = sm1(gid)/nsm1(gid)
!     else
!        sm1(gid) = LIS_rc%udef
!     end if
!     if(nst1(gid).gt.0) then 
!        st1(gid) = st1(gid)/nst1(gid)
!     else
!        st1(gid) = LIS_rc%udef
!     end if
  enddo

  
  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status)

!  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
!       .false., rc=status)
!  call LIS_verify(status)
!  return

  deallocate(nqle)
  deallocate(nqh)
  deallocate(nqg)
  deallocate(nsm1)
  deallocate(nst1)
end subroutine read_AmerifluxObs



