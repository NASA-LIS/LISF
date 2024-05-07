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
! !ROUTINE: read_ARSsmobs
! \label{read_ARSsmobs}
!
! !REVISION HISTORY:
!  26 Jan 2018    Soni Yatheendradas;   Initial Specification based on LVT ARSsm datastream
!
! !INTERFACE: 
subroutine read_ARSsmobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_timeMgrMod, only : LIS_calendar, LIS_tick
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_fileIOMod,      only : LIS_readData
  use ARSsm_obsMod, only : ARSsm_obs_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_Space] Objective Space
!  \end{description}
!
!EOP
  real,    pointer    :: smc(:)
!  real,    pointer    :: smstd(:)

  type(ESMF_Field)    :: smcField
!  type(ESMF_Field)    :: smstdField
  character(len=LIS_CONST_PATH_LEN) :: obsdir
  logical             :: data_update
  integer             :: i
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer               :: ios
  integer               :: yr,doy,mo,da,hr,mn,ss
  logical               :: file_exists
  integer               :: status
  type(ESMF_Time)       :: obstime,obstime1
  integer               :: ftn
  integer               :: k,c,r,gid
  integer               :: tind
  real                  :: sfsm, sfsm_std
  integer             :: n 

  n = 1

  call ESMF_AttributeGet(Obj_Space,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  if(ARSsm_obs_struc(n)%yr.ne.LIS_rc%yr) then 

     ARSsm_obs_struc(n)%sm = LIS_rc%udef
     ARSsm_obs_struc(n)%sm_std = LIS_rc%udef

     ARSsm_obs_struc(n)%yr = LIS_rc%yr

     call ESMF_TimeSet(ARSsm_obs_struc(n)%startTime, yy=LIS_rc%yr, &
          mm=1, dd=1, h=0, m=0,s = 0, calendar=LIS_calendar, rc=status) !SY: 2nd argument "mm=" made 1 as against "mm=LIS_rc%mo" done in the LVT datastream
     call LIS_verify(status, 'ARSsm starttime set failed')

     do k = 1,ARSsm_obs_struc(n)%n_stns

        call create_ARSsmobs_filename(ARSsm_obs_struc(n)%odir,&
             ARSsm_obs_struc(n)%stn_name(k), &
             LIS_rc%yr, &
             filename)

        inquire(file=trim(filename), exist=file_exists)
        if(file_exists) then
           write(LIS_logunit,*) '[INFO] Reading ',trim(filename)
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=trim(filename),form='formatted')
           ios = 0
           do while(ios.eq.0)
              read(ftn,*,iostat=ios) yr, mo,da, hr, mn, sfsm, sfsm_std

              if(ios.ne.0) exit
              call ESMF_TimeSet(obstime, yy=yr,mm=mo,dd=da,h=hr,m=mn,s=0,&
                   calendar=LIS_calendar,rc=status)
              call LIS_verify(status, 'ESMF_TimeSet in read_ARSsmObs')

              tind = nint((obstime-ARSsm_obs_struc(n)%startTime)/&
                   ARSsm_obs_struc(n)%timestep(k))+1

              !if(tind.gt.0 .and. sfsm.gt.0 .and. sfsm.lt.0.5) then ! SY see next line 
              if(tind.gt.0 .and. sfsm.ge.0.025 .and. sfsm_std.gt.0) then ! SY 02/27/18: NOTE that this 
                 ! "sfsm.ge.0.025" required condition is mentioned in 
                 ! /discover/nobackup/projects/lis/STN_DATA/ARS_Watersheds/ARS_Data.Readme. 
                 ! Also, I removed the "sfsm.lt.0.5" since I found a few data files having 
                 ! values going above 0.5
                 ! Finally, I consider only values corresponding to sfsm_std > 0. Sujay mentioned to
                 ! consider only non-negative sfsm_std values per email exchange 02/27/18 after I
                 ! noticed [1] hourly files containing negative sfsm_std values, and [2]  hourly 
                 ! Reynolds_Creek*.txt files missing leap year day (Feb 29) time instants and values.  
                 ARSsm_obs_struc(n)%sm(k,tind) = sfsm
                 ! SY: Begin commenting out
                 !if(sfsm_std.gt.0) then
                 !   ARSsm_obs_struc(n)%sm_std(k,tind) = sfsm_std
                 !else
                 !   ARSsm_obs_struc(n)%sm_std(k,tind) = LIS_rc%udef
                 !endif
                 ! SY: End commenting out
                 ARSsm_obs_struc(n)%sm_std(k,tind) = sfsm_std ! SY: this replaces the commenred out 'if' structure immediately above
              endif
           enddo ! do while(ios.eq.0)
           call LIS_releaseUnitNumber(ftn)
           write(LIS_logunit,*) '[INFO] Finished processing ',trim(filename)
        endif !  if(file_exists) then

     end do ! do k = 1,ARSsm_obs_struc(n)%n_stns
  endif !  if(ARSsm_obs_struc(n)%yr.ne.LIS_rc%yr) then 

  call ESMF_TimeSet(obstime1, yy=LIS_rc%yr, &
       mm=LIS_rc%mo, dd=LIS_rc%da, h=LIS_rc%hr, m=LIS_rc%mn, &
       s = LIS_rc%ss, calendar=LIS_calendar, rc=status)
  call LIS_verify(status, 'obstime1 set failed')

!  write(LIS_logunit,*) '[INFO] Here 1 in read_ARSsmobs.F90 '

  call ESMF_StateGet(Obj_Space,"ARS_sm",smcField,&
       rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=smc,rc=status)
  call LIS_verify(status)

!  call ESMF_StateGet(Obj_Space,"ARSsm standard deviation of soil moisture",smstdField,&
!       rc=status)
!  call LIS_verify(status)
!  
!  call ESMF_FieldGet(smstdField,localDE=0,farrayPtr=smstd,rc=status)
!  call LIS_verify(status)  

  smc = LIS_rc%udef ! SY
!  smstd = LIS_rc%udef ! SY
 
!  write(LIS_logunit,*) '[INFO] Here 2 in read_ARSsmobs.F90 '

  do i=1,ARSsm_obs_struc(n)%n_stns

!     write(LIS_logunit,*) '[INFO] stn col and row are:', ARSsm_obs_struc(n)%stn_col(i),',',&
!                          ARSsm_obs_struc(n)%stn_row(i)
     gid = LIS_domain(n)%gindex(ARSsm_obs_struc(n)%stn_col(i),ARSsm_obs_struc(n)%stn_row(i))

     if(ARSsm_obs_struc(n)%stn_col(i).ge.1.and.&
        ARSsm_obs_struc(n)%stn_col(i).le.LIS_rc%lnc(n).and.&
        ARSsm_obs_struc(n)%stn_row(i).ge.1.and.&
        ARSsm_obs_struc(n)%stn_row(i).le.LIS_rc%lnr(n).and.&
        gid.ne.-1) then 

        tind = nint((obstime1 - ARSsm_obs_struc(n)%starttime)/ARSsm_obs_struc(n)%timestep(i))+1
!        write(LIS_logunit,*) '[INFO] tind=',tind,',obstime1=',obstime1,',starttime=',ARSsm_obs_struc(n)%starttime,',timestep = ',ARSsm_obs_struc(n)%timestep(i) 

        smc(gid) = ARSsm_obs_struc(n)%sm(i,tind)

!        smstd(gid) = ARSsm_obs_struc(n)%sm_std(i,tind)

     endif
  end do ! do i=1,ARSsm_obs_struc(n)%n_stns
  
  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status)

!  write(LIS_logunit,*) '[INFO] Finished read_ARSsmobs.F90, smc = ', smc
end subroutine read_ARSsmobs

!BOP
!
! !ROUTINE: create_ARSsmobs_filename
! \label(create_ARSsmobs_filename)
!
! !INTERFACE:
subroutine create_ARSsmobs_filename(odir, stn, yr, filename)
!
! !USES:
  implicit none
!
! !INPUT PARAMETERS:
  character(len=*), intent(in) :: odir
  character(len=*), intent(in) :: stn
  integer  ,        intent(in) :: yr
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
  character(len=*)             :: filename

  character*4 :: fyr

  write(unit=fyr, fmt='(i4.4)') yr

  filename = trim(odir)//'/'//trim(stn)//'_'//trim(fyr)//'.txt'

end subroutine create_ARSsmobs_filename


