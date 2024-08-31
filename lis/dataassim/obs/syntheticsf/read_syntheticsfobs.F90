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
! !ROUTINE: read_syntheticsfobs
! \label{read_syntheticsfobs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_syntheticsfobs(n, OBS_State, OBS_Pert_state) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_historyMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod
  use LIS_pluginIndices
  use syntheticsfobs_module
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the synthetic streamflow observations 
!  produced from a LIS control run and packages it 
!  into an ESMF State with certain predefined 
!  attributes
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: sffield

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%ngrid(n))
  integer             :: assimflag(LIS_rc%ngrid(n))
  real                :: obs_unsc(LIS_rc%ngrid(n))
  character(len=LIS_CONST_PATH_LEN) :: sfobsdir
  logical             :: data_update
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name
  integer             :: fnd
  integer             :: ftn,p
  logical             :: readflag
  integer             :: status
  logical             :: data_upd_flag(LIS_npes)
  logical             :: data_upd_flag_local
  logical             :: data_upd
  real                :: sf_current(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer             :: t,c,r


  obs_unsc = LIS_rc%udef
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sfobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call synsm_filename(name,sfobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)

  inquire(file=name,exist=file_exists)

  if(file_exists) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif
  
  if (readflag) then 
     write(LIS_logunit,*)  'Reading synthetic data ',trim(name)
     
     call ESMF_StateGet(OBS_State,"Observation01",sffield,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(sffield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
     ftn = LIS_getNextUnitNumber()

     open(ftn,file=name,form='unformatted')
     call LIS_readvar_gridded(ftn, n, obsl)
     call LIS_releaseUnitNumber(ftn)

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     
     
     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
          //trim(LIS_synsmId)//char(0),n, OBS_state)

     call ESMF_StateGet(OBS_State,"Observation01",sffield,&
          rc=status)
     call LIS_verify(status,'ESMF_StateGet failed in read_syntheticsfobs')

     call ESMF_FieldGet(sffield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed in read_syntheticsfobs')

     fnd = 0 
     do r =1,LIS_rc%lnr(n)
        do c =1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r) .ne. -1)then
              sf_current(c,r) = obsl(LIS_domain(n)%gindex(c,r))
              obs_unsc(LIS_domain(n)%gindex(c,r)) = &
                   obsl(LIS_domain(n)%gindex(c,r))
              if(sf_current(c,r).ne.LIS_rc%udef) then 
                 fnd = 1
              endif
           end if
        end do
     end do
     
     if(fnd.eq.0) then 
        data_upd_flag_local = .false.
     else
        data_upd_flag_local = .true.
     endif

#if (defined SPMD)
     call MPI_ALLGATHER(data_upd_flag_local,1, &
          MPI_LOGICAL, data_upd_flag(:),&
          1, MPI_LOGICAL, LIS_mpi_comm, status)
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
              fnd = 1
           else
              assimflag(t) = 0
           endif
        enddo
        
        
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .true., rc=status)
        call LIS_verify(status)
        
        if(LIS_rc%ngrid(n).gt.0) then 
           call ESMF_AttributeSet(sffield,"Grid Number",&
                gid,itemCount=LIS_rc%ngrid(n),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(sffield,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(sffield, "Unscaled Obs",&
                obs_unsc, itemCount=LIS_rc%ngrid(n), rc=status)
           call LIS_verify(status, 'Error in setting Unscaled Obs attribute')      
        endif
     else
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .false., rc=status)
        call LIS_verify(status)
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  end if
  
end subroutine read_syntheticsfobs

subroutine synsm_filename(name, ndir, yr, mo,da,hr,mn)
  
  implicit none
  character(len=*)  :: name
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr,fmn
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  write(unit=fmn, fmt='(i2.2)') mn  
  
  name = trim(ndir)//'/SFOBS_'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
       trim(fmn)//'.bin'
end subroutine synsm_filename



