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
! !ROUTINE: read_WindSatCsm
! \label{read_WindSatCsm}
!
! !REVISION HISTORY:
!  22 Dec 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_WindSatCsm(n, OBS_State, OBS_Pert_state) 
! !USES: 
  use ESMF
  use LIS_pluginIndices, only : LIS_WindSatCsmobsId
  use LIS_mpiMod
  use LIS_coreMod,  only : LIS_rc, LIS_domain, LIS_localPet, LIS_npes
  use LIS_timeMgrMod, only : LIS_calendar, LIS_clock, LIS_get_julss
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, & 
       LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use WindSatCsm_Mod, only : WindSatCsm_struc
  use map_utils
  
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
!  Current limitations 
!  * only neighbor search is supported
!  * the assimilation interval is expected to be 3 hours
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  real, parameter     :: MAX_SM_VALUE=0.55
  integer, parameter  :: windsat_nc=1383, windsat_nr=586
  type(ESMF_Field)    :: smfield
  integer             :: currTime

  integer             :: i,p, ierr
  integer             :: col,row,gridid
  real                :: rcol, rrow
  real,    pointer    :: obsl(:)
  real                :: smobs(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                :: obs_unsc(LIS_rc%ngrid(n))
  integer             :: gid(LIS_rc%ngrid(n))
  integer             :: assimflag(LIS_rc%ngrid(n))

  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  logical             :: data_update
  logical             :: file_exists

  logical             :: data_upd
  logical             :: data_upd_flag(LIS_npes)
  logical             :: data_upd_flag_local
  integer             :: status

  integer             :: t
  integer             :: iret
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

  smobs    = LIS_rc%udef
  obs_unsc = LIS_rc%udef

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status,'Error in AttributeGet: Data Directory')
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status, 'Error in AttributeGet: Data Update Status')

!-------------------------------------------------------------------------
!   Read the data at 0z
!-------------------------------------------------------------------------
  data_upd = .false. 


  call ESMF_StateGet(OBS_State,"Observation01",smfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
     
  obsl = LIS_rc%udef 
  fnd  = 0 

!Tolouse location
  call latlon_to_ij(LIS_domain(n)%lisproj,43.25,1.25,rcol,rrow)
  
  c = nint(rcol)
  r = nint(rrow)
  
!Ignore data during winter months.            
  if(LIS_rc%mo.eq.11.or.LIS_rc%mo.eq.12.or.LIS_rc%mo.eq.1.or.LIS_rc%mo.eq.2) then 
     fnd = 0 
     obsl(LIS_domain(n)%gindex(c,r)) = LIS_rc%udef
  else
     call LIS_get_julss(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, 0, 0, &
          currTime)
     do t=1,182
        dt = WindSatCsm_struc(n)%tmobs(t)-currTime
        if(dt.ge.0.and.dt.le.LIS_rc%ts.and.&
             WindSatCsm_struc(n)%smobs(t).gt.0) then 
           obsl(LIS_domain(n)%gindex(c,r)) = &
                WindSatCsm_struc(n)%smobs(t)
           fnd = 1
        endif
     enddo
  endif


!lsm-based qc

  call lsmdaqcobsstate(LIS_rc%lsm, LIS_WindSatCsmobsId , &
       n, OBS_state)

  call ESMF_StateGet(OBS_State,"Observation01",smField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)

  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           smobs(c+LIS_rc%lnc(n)*(r-1)) = obsl(LIS_domain(n)%gindex(c,r))
        end if
     end do
  end do

!-------------------------------------------------------------------------
!  Store unscaled data for future output
!-------------------------------------------------------------------------     
  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           obs_unsc(LIS_domain(n)%gindex(c,r))=smobs(c+LIS_rc%lnc(n)*(r-1))
        end if
     end do
  end do


!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     
  if(WindSatCsm_struc(n)%scal.ne.0) then 
     
     do t=1,LIS_rc%ngrid(n)
        model_delta(t) = WindSatCsm_struc(n)%model_xrange(t,2)-&
             WindSatCsm_struc(n)%model_xrange(t,1)
        obs_delta(t) = WindSatCsm_struc(n)%obs_xrange(t,2)-&
             WindSatCsm_struc(n)%obs_xrange(t,1)
     enddo
     do t=1,LIS_rc%ngrid(n)
        
        col = LIS_domain(n)%grid(t)%col
        row = LIS_domain(n)%grid(t)%row
        
        gridid = col+(row-1)*LIS_rc%lnc(n)
        if(smobs(gridid).ne.-9999.0) then 
           if(obs_delta(t).gt.0) then 
              binval = nint((smobs(gridid)-WindSatCsm_struc(n)%obs_xrange(t,1))/&
                   obs_delta(t))+1
              if(binval.gt.WindSatCsm_struc(n)%nbins) binval = WindSatCsm_struc(n)%nbins
              if(binval.le.0) binval = 1
              cdf_obsval = WindSatCsm_struc(n)%obs_cdf(t,binval)
              if(cdf_obsval.gt.1.0) cdf_obsval = 1.0
              i=1
              do while((WindSatCsm_struc(n)%model_cdf(t,i).lt.cdf_obsval).and.&
                   (i.le.WindSatCsm_struc(n)%nbins))
                 i = i+1
              enddo
              if(i.gt.WindSatCsm_struc(n)%nbins) i = i-1
              smvalue = WindSatCsm_struc(n)%model_xrange(t,i)
              
              if(smvalue.gt.MAX_SM_VALUE) then 
                 write(LIS_logunit,*) 'Problem in scaling SM observations'
                 write(LIS_logunit,*) t, smobs(gridid),smvalue, MAX_SM_VALUE
                 call LIS_endrun()
              endif
              
              if(smvalue.eq.0) then 
                 smvalue = smobs(gridid)
              endif
              smobs(gridid) = smvalue
           else
              print*, 'WARNING: obs_delta for grid point', t, '=0'
              smobs(gridid) = LIS_rc%udef
           endif
        else
           smobs(gridid) = LIS_rc%udef
        endif
     enddo
  endif

  fnd =0 
  data_upd_flag(LIS_localPet+1) = .false.
  do t=1,LIS_rc%lnc(n)*LIS_rc%lnr(n)
     if(smobs(t).ne.LIS_rc%udef) then 
        fnd = 1
     endif        
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
              obsl(LIS_domain(n)%gindex(c,r))=smobs(c+LIS_rc%lnc(n)*(r-1))
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

     call ESMF_AttributeSet(smField,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)
  
     call ESMF_AttributeSet(smField,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(smfield, "Unscaled Obs", &
          obs_unsc, itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status, 'Error in setting Unscaled Obs attribute')       

  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif

end subroutine read_WindSatCsm



