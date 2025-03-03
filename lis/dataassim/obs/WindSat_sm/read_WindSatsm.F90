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
! !ROUTINE: read_WindSatsm
! \label{read_WindSatsm}
!
! !REVISION HISTORY:
!  22 Dec 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_WindSatsm(n, OBS_State, OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use map_utils
  use LIS_pluginIndices
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use WindSatsm_Mod, only : WindSatsm_struc


  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the WindSat soil moisture observations
!  and packages it into an ESMF State with certain predefined 
!  attributes. The routine also interpolates the soil moisture
!  data to the LIS resolution
!
!  Current limitations \newline
!  * only neighbor search is supported \newline
!  * the assimilation interval is expected to be 3 hours
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  real,  parameter       :: MAX_SM_VALUE=0.45, MIN_SM_VALUE=0.0001
  integer, parameter  :: windsat_nc=1383, windsat_nr=586
  type(ESMF_Field)    :: smfield
  integer             :: currTime, refTime
  integer             :: yr, mo, da, hour, min, ss
  integer             :: i, p, istat, ierr
  integer             :: col,row,gridid
  real,    pointer    :: obsl(:)
  real                :: sm_current(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                :: obs_unsc(LIS_rc%ngrid(n))
  integer             :: gid(LIS_rc%ngrid(n))
  integer             :: assimflag(LIS_rc%ngrid(n))

  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  logical             :: data_update
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: smname,tmname,tsname,clsname

  logical             :: readflag
  logical             :: data_upd
  logical             :: data_upd_flag(LIS_npes)
  logical             :: data_upd_flag_local
  integer             :: status

  integer             :: t
  integer             :: ftn
  integer             :: index1, iret
  logical*1           :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  logical*1           :: li(windsat_nc*windsat_nr)
  real                :: sm(windsat_nc*windsat_nr)
  real                :: ts(windsat_nc*windsat_nr)
  real*8              :: tm(windsat_nc*windsat_nr)
  real                :: tm1(windsat_nc*windsat_nr)
  integer*2               :: cls(windsat_nc*windsat_nr)
!  integer*1           :: hr(windsat_nc,windsat_nr)
!  integer*1           :: mn(windsat_nc,windsat_nr)
  

  integer             :: binval
  real                :: cdf_obsval
  real                :: smvalue
  real                :: model_delta(LIS_rc%ngrid(n))
  real                :: obs_delta(LIS_rc%ngrid(n))
  logical             :: alarmCheck
  integer             :: c,r
  real                :: dt
  integer             :: fnd

  sm_current    = LIS_rc%udef
  obs_unsc = LIS_rc%udef

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status,'Error in AttributeGet: Data Directory')
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status, 'Error in AttributeGet: Data Update Status')

  data_upd = .false. 

!-------------------------------------------------------------------------
!   Read the data at 0z
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "WindSat sm read alarm")  
  if(alarmCheck.or.WindSatsm_struc(n)%startMode) then 

     WindSatsm_struc(n)%startMode = .false.

     call WindSatsm_filename(smname, tmname,tsname, clsname, smobsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)         

     inquire(file=trim(smname),exist=file_exists)
     if(file_exists) then 
        readflag = .true. 
     else 
        readflag = .false.
     endif
     
     
     if (readflag) then 

        write(LIS_logunit,*)  'Reading WindSat soil moisture data ',&
             trim(smname)
        
        ftn = LIS_getNextUnitNumber()
        open(ftn,file = trim(smname), form='unformatted', status='old',&
             access='direct',recl= WindSatsm_struc(n)%mi*4, iostat = istat)
        if(istat.eq.0) then 
           read(ftn,rec=1) sm
        endif
        call LIS_releaseUnitNumber(ftn)

        ftn = LIS_getNextUnitNumber()
        open(ftn,file = trim(tmname), form='unformatted', status='old',&
             access='direct',recl= WindSatsm_struc(n)%mi*8, iostat = istat)
        if(istat.eq.0) then 
           read(ftn,rec=1) tm
        endif
        call LIS_releaseUnitNumber(ftn)

        ftn = LIS_getNextUnitNumber()
        open(ftn,file = trim(tsname), form='unformatted', status='old',&
             access='direct',recl= WindSatsm_struc(n)%mi*4, iostat = istat)
        if(istat.eq.0) then 
           read(ftn,rec=1) ts
        endif
        call LIS_releaseUnitNumber(ftn)


        ftn = LIS_getNextUnitNumber()
        open(ftn,file = trim(clsname), form='unformatted', status='old',&
             access='direct',recl= WindSatsm_struc(n)%mi*8, iostat = istat)
        if(istat.eq.0) then 
           read(ftn,rec=1) cls
        endif
        call LIS_releaseUnitNumber(ftn)
        
!-------------------------------------------------------------------------
!   Interpolate to the LIS domain using a neighbor-search approach
!-------------------------------------------------------------------------     
        li  = .false.

        do c=1,WindSatsm_struc(n)%mi
           if(sm(c).ne.-999.0) then 
              if(cls(c).eq.80.or.cls(c).eq.90.and.ts(c).gt.0) & 
                   li(c) = .true. 
           endif
        enddo       
        
        call neighbor_interp(LIS_rc%gridDesc(n,:),li,sm,&
             lo,WindSatsm_struc(n)%smobs,&
             WindSatsm_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             WindSatsm_struc(n)%n11,LIS_rc%udef,iret)

!        open(100,file='sm_ip.bin',form='unformatted')
!        write(100) WindSatsm_struc(n)%smobs
!        close(100) 
!        stop
        do c=1,WindSatsm_struc(n)%mi           
           tm1(c) = tm(c)
        enddo

        call neighbor_interp(LIS_rc%gridDesc(n,:),li,tm1,&
             lo,WindSatsm_struc(n)%tmobs,&
             WindSatsm_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             WindSatsm_struc(n)%n11,LIS_rc%udef,iret)
        
     endif
  endif

  call ESMF_StateGet(OBS_State,"Observation01",smfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
     
  obsl = LIS_rc%udef 
  fnd  = 0 

  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if ((LIS_domain(n)%gindex(c,r) .ne. -1).and.&
             (WindSatsm_struc(n)%smobs(c+LIS_rc%lnc(n)*(r-1)).ne.LIS_rc%udef)) then 
           
!Ignore data during winter months.            
           if(LIS_rc%mo.eq.11.or.LIS_rc%mo.eq.12.or.LIS_rc%mo.eq.1.or.LIS_rc%mo.eq.2) then 
              fnd = 0 
              obsl(LIS_domain(n)%gindex(c,r)) = LIS_rc%udef
           else
              yr = 2000
              mo = 1
              da = 1
              hour = 12
              min = 0 
              ss = 0 
              
              call LIS_get_julss(yr, mo, da, hour, min, ss, reftime)
              call LIS_get_julss(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, &
                   LIS_rc%mn, LIS_rc%ss, &
                   currTime)
              dt = WindSatsm_struc(n)%tmobs(c+(r-1)*LIS_rc%lnc(n)) + reftime - & 
                   currTime
              if(WindSatsm_struc(n)%tmobs(c+(r-1)*LIS_rc%lnc(n)).gt.0.0.and.&
                   dt .ge.0.and.dt.le.LIS_rc%ts) then 
                 obsl(LIS_domain(n)%gindex(c,r))=&
                      WindSatsm_struc(n)%smobs(c+LIS_rc%lnc(n)*(r-1))                 
                 fnd = 1
              else
                 fnd = 0 
                 obsl(LIS_domain(n)%gindex(c,r)) = LIS_rc%udef
              endif
           endif
        end if
     end do
  end do

!lsm-based qc

  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_WindSatsmobsId)//char(0),n, OBS_state)

  call ESMF_StateGet(OBS_State,"Observation01",smField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)

  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           sm_current(c,r) = obsl(LIS_domain(n)%gindex(c,r))
        end if
     end do
  end do

!-------------------------------------------------------------------------
!  Store unscaled data for future output
!-------------------------------------------------------------------------     
  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           obs_unsc(LIS_domain(n)%gindex(c,r))=sm_current(c,r)
        end if
     end do
  end do

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     
  if(WindSatsm_struc(n)%scal.ne.0.and.fnd.ne.0) then 
     call LIS_rescale_with_CDF_matching(    &
          n,                                   & 
          WindSatsm_struc(n)%nbins,         & 
          MAX_SM_VALUE,                        & 
          MIN_SM_VALUE,                        & 
          WindSatsm_struc(n)%model_xrange,  &
          WindSatsm_struc(n)%obs_xrange,    &
          WindSatsm_struc(n)%model_cdf,     &
          WindSatsm_struc(n)%obs_cdf,       &
          sm_current)

  endif

  fnd =0 
!  data_upd_flag(LIS_localPet+1) = .false.
  data_upd_flag_local = .false.

  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if(sm_current(c,r).ne.LIS_rc%udef) then
           fnd = 1
        endif
     enddo
  enddo

  call ESMF_StateGet(OBS_State,"Observation01",smField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)
  
  obsl = LIS_rc%udef

  if(fnd.eq.0) then 
     obsl = LIS_rc%udef
     
  else     
     do r =1,LIS_rc%lnr(n)
        do c =1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r) .ne. -1)then
              obsl(LIS_domain(n)%gindex(c,r))=sm_current(c,r)
           end if
        end do
     end do
  endif

  if(fnd.eq.0) then 
     data_upd_flag_local = .false. 
  else
     data_upd_flag_local = .true. 
  endif


#if (defined SPMD)
  call MPI_ALLGATHER(data_upd_flag_local,1, &
       MPI_LOGICAL, data_upd_flag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, ierr)
#endif
  data_upd = .false.
  do p=1,LIS_npes
     data_upd = data_upd.or.data_upd_flag(p)
  enddo
  
  if(data_upd) then 
     do t=1,LIS_rc%ngrid(n)
        gid(t) = t
        if(obsl(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo
  
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true. , rc=status)
     call LIS_verify(status)

     if(LIS_rc%ngrid(n).gt.0) then 
        call ESMF_AttributeSet(smField,"Grid Number",&
             gid,itemCount=LIS_rc%ngrid(n),rc=status)
        call LIS_verify(status)
        
        call ESMF_AttributeSet(smField,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
        call LIS_verify(status)
        
     endif

     call ESMF_AttributeSet(smfield, "Unscaled Obs", &
          obs_unsc,itemCount=LIS_rc%ngrid(n), rc=status)
     call LIS_verify(status, 'Error in setting Unscaled Obs attribute')       

  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif

end subroutine read_WindSatsm

!BOP
!
! !ROUTINE: WindSatsm_filename
! \label{WindSatsm_filename}
! 
! !INTERFACE: 
subroutine WindSatsm_filename(smname, tmname, tsname, clsname, ndir, yr, mo,da)
  
  implicit none
! !ARGUMENTS: 
  character(len=*)     :: smname
  character(len=*)     :: tmname
  character(len=*)     :: tsname
  character(len=*)     :: clsname
  
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the WindSat filenames based the time and date
!  
!  The arguments are: 
!  \begin{description}
!  \item[smname] name of the WindSat soil moisture filename
!  \item[tmname] name of the WindSat time information filename
!  \item[ndir] name of the WindSat soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da 

  smname = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/d'//trim(fyr)//trim(fmo)//trim(fda)//'GEZ25av_d.sm2'

  tmname = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/d'//trim(fyr)//trim(fmo)//trim(fda)//'GEZ25av_d.tm2'
  
  tsname = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/d'//trim(fyr)//trim(fmo)//trim(fda)//'GEZ25av_d.ts2'

  clsname = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/d'//trim(fyr)//trim(fmo)//trim(fda)//'GEZ25av_d.cls2'


end subroutine WindSatsm_filename



