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
! !ROUTINE: read_simGRACEJPLobs
! \label{read_simGRACEJPLobs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!  11Aug2011: Ben Zaitchik; Modified for simGRACEJPL
!
! !INTERFACE: 
subroutine read_simGRACEJPLobs(n, OBS_State, OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use simGRACEJPLobs_module
  use LIS_fileIOMod
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the simGRACEJPL TWS observations produced using a
!  LIS open loop run and simGRACEJPL TELLUS gridded files, 
!  and packages it into an ESMF State with certain  
!  predefined attributes
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: twsfield

  integer             :: iret

  real,    pointer    :: obsl(:)
  real                :: gracegrid(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                :: gracegrid_glb(LIS_rc%gnc(n),LIS_rc%gnr(n))
  real                :: graceerr(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                :: graceerr_glb(LIS_rc%gnc(n),LIS_rc%gnr(n))
  real                :: gridDesc(6)
  integer             :: gid(LIS_rc%ngrid(n))
  integer             :: assimflag(LIS_rc%ngrid(n))

  character(len=LIS_CONST_PATH_LEN) :: simGRACEJPLobsdir
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name

  integer             :: col,row
  logical             :: data_upd
  logical             :: data_upd_flag_local
  logical             :: data_upd_flag(LIS_npes)
  integer             :: fnd
  logical             :: readflag
  logical             :: alarmCheck, file_status
  integer             :: status,ftn
  real, allocatable       :: ssdev(:)
  integer             :: t,p,r,c
  real             :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/  

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       simGRACEJPLobsdir, rc=status)
  call LIS_verify(status)
  
  call ESMF_AttributeSet(OBS_State,"Data Update Status",&
       .false., rc=status)
  call LIS_verify(status)
  
  call checkStatus_simGRACEJPLfile(simGRACEJPLobsdir,&
       simGRACEJPL_struc(n)%config, &
       file_status)
  
  call ESMF_AttributeSet(OBS_State,"File Status",&
       file_status, rc=status)
  call LIS_verify(status)
  
  alarmCheck = .false. 
  if( LIS_rc%hr.eq.simGRACEJPL_struc(n)%alarmHr.and.&
       LIS_rc%mn.eq.0.and.&
       LIS_rc%ss.eq.0) then 
     call checkSimGRACEJPLdataAlarm(simGRACEJPL_struc(n)%config, &
          alarmCheck)
  endif

  data_upd = .false. 

  if(alarmCheck) then 
     if(LIS_rc%DAincrMode(n).eq.1) then 
        call simGRACEJPL_filename(name,simGRACEJPLobsdir,&
             simGRACEJPL_struc(n)%config, &
             LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
        
        inquire(file=name,exist=file_exists)
        if(file_exists) then 
           readflag = .true. 
        else 
           readflag = .false.
        endif

        if (readflag) then 

           call ESMF_AttributeSet(OBS_State,"File Status",&
                .true., rc=status)
           call LIS_verify(status)

           if(simGRACEJPL_struc(n)%config.eq."GRACE-2") then 
              call ESMF_AttributeSet(OBS_State,&
                   name="Data averaging factor",&
                   value=13.0,rc=status)
              call LIS_verify(status)
           else
              call ESMF_AttributeSet(OBS_State,&
                   name="Data averaging factor",&
                   value=days(LIS_rc%mo),rc=status)
              call LIS_verify(status)
           endif


           write(LIS_logunit,*)  'Reading simGRACEJPL data ',trim(name)

           ftn = LIS_getNextUnitNumber()
           open(ftn, file=trim(name),form='unformatted')
           read(ftn) gracegrid_glb
           read(ftn) graceerr_glb
           call LIS_releaseUnitNumber(ftn)
           write(LIS_logunit,*)  'Done reading simGRACEJPL data ',trim(name)

           gracegrid = gracegrid_glb(&
                LIS_ews_halo_ind(n,LIS_localPet+1):&         
                LIS_ewe_halo_ind(n,LIS_localPet+1), &
                LIS_nss_halo_ind(n,LIS_localPet+1): &
                LIS_nse_halo_ind(n,LIS_localPet+1))

           graceerr = graceerr_glb(&
                LIS_ews_halo_ind(n,LIS_localPet+1):&         
                LIS_ewe_halo_ind(n,LIS_localPet+1), &
                LIS_nss_halo_ind(n,LIS_localPet+1): &
                LIS_nse_halo_ind(n,LIS_localPet+1))

           fnd = 0 
           data_upd_flag_local = .false. 
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(gracegrid(c,r).ne.LIS_rc%udef) then 
                    fnd = 1
                 endif
              enddo
           enddo

           if(simGRACEJPL_struc(n)%useDistErr.eq.1) then 
              if(LIS_rc%ngrid(n).gt.0) then 
                 allocate(ssdev(LIS_rc%ngrid(n)))
                 do r=1,LIS_rc%lnr(n)
                    do c=1,LIS_rc%lnc(n)
                       if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                          ssdev(LIS_domain(n)%gindex(c,r)) = & 
                               graceerr(c,r)
                       endif
                    enddo
                 enddo

                 call ESMF_AttributeSet(OBS_Pert_State,"Standard Deviation",&
                      ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
                 call LIS_verify(status, 'Error with ESMF_AttributeGet in read_simGRACEJPLobs')
                 deallocate(ssdev)
              endif
           endif

           call ESMF_StateGet(OBS_State,"Observation01",twsfield,&
                rc=status)
           call LIS_verify(status, 'ESMF_StateGet failed in read_simGRACEJPLobs')

           call ESMF_FieldGet(twsfield,localDE=0,farrayPtr=obsl,rc=status)
           call LIS_verify(status,'ESMF_FieldGet failed in read_simGRACEJPLobs')

           obsl(:) = -9999.0

           if(fnd.eq.0) then 
              obsl = LIS_rc%udef
           else
              do r=1,LIS_rc%lnr(n)
                 do c=1,LIS_rc%lnc(n)
                    if(LIS_domain(n)%gindex(c,r).ne.-1) then
                       if(gracegrid(c,r).gt.0.0) then
                          obsl(LIS_domain(n)%gindex(c,r)) = gracegrid(c,r)
                       end if
                    endif
                 end do
              end do
           endif

           if(fnd.eq.0) then 
              data_upd_flag_local = .false.
           else
              data_upd_flag_local = .true. 
           endif
#if(defined SPMD)
           call MPI_ALLGATHER(data_upd_flag_local,1,&
                MPI_LOGICAL, data_upd_flag(:),&
                1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
           data_upd = .false.
           do p=1,LIS_npes
              data_upd = data_upd.or.data_upd_flag(p)
           enddo

           write(LIS_logunit,*) 'MSG: read simGRACEJPL data'

           readflag = .false.

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
                   .true., rc=status)
              call LIS_verify(status)

              if(LIS_rc%ngrid(n).gt.0) then 
                 call ESMF_AttributeSet(twsfield,"Grid Number",&
                      gid,itemCount=LIS_rc%ngrid(n),rc=status)
                 call LIS_verify(status)

                 call ESMF_AttributeSet(twsfield,"Assimilation Flag",&
                      assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
                 call LIS_verify(status)
              endif
           endif
        endif

     endif
  end if
  
end subroutine read_simGRACEJPLobs

subroutine simGRACEJPL_filename(name, ndir, config, yr, mo,da)
  
  implicit none
  character(len=*)  :: name
  character (len=*) :: ndir
  character (len=*) :: config
  integer           :: yr, mo,da

  character (len=4) :: fyr
  character (len=2) :: fmo
  character (len=2) :: fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  
  if(config.eq."GRACE-2") then 
     name = trim(ndir)//'/GRACE_obs_'//trim(fyr)//trim(fmo)//trim(fda)//&
          '.bin'
  else
     name = trim(ndir)//'/GRACE_obs_'//trim(fyr)//trim(fmo)//&
          '.bin'
  endif
end subroutine simGRACEJPL_filename


subroutine checkSimGRACEJPLdataAlarm(config, alarmCheck)

  use LIS_coreMod,  only : LIS_rc

  character(len=*) :: config
  logical          :: alarmCheck

  if(config.eq."GRACE-2") then 
     
     if(LIS_rc%mo.eq.1.and. LIS_rc%da.eq.1) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.1.and.LIS_rc%da.eq.14) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.1.and.LIS_rc%da.eq.27) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.2.and.LIS_rc%da.eq.9) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.2.and.LIS_rc%da.eq.22) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.3.and.LIS_rc%da.eq.7) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.3.and.LIS_rc%da.eq.20) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.4.and.LIS_rc%da.eq.2) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.4.and.LIS_rc%da.eq.15) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.4.and.LIS_rc%da.eq.28) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.5.and.LIS_rc%da.eq.11) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.5.and.LIS_rc%da.eq.24) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.6.and.LIS_rc%da.eq.6) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.6.and.LIS_rc%da.eq.19) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.7.and.LIS_rc%da.eq.2) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.7.and.LIS_rc%da.eq.15) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.7.and.LIS_rc%da.eq.28) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.8.and.LIS_rc%da.eq.10) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.8.and.LIS_rc%da.eq.23) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.9.and.LIS_rc%da.eq.5) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.9.and.LIS_rc%da.eq.18) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.10.and.LIS_rc%da.eq.1) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.10.and.LIS_rc%da.eq.14) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.10.and.LIS_rc%da.eq.27) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.11.and.LIS_rc%da.eq.9) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.11.and.LIS_rc%da.eq.22) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.12.and.LIS_rc%da.eq.5) then 
        alarmCheck = .true. 
     elseif(LIS_rc%mo.eq.12.and.LIS_rc%da.eq.18) then 
        alarmCheck = .true. 
     else
        alarmCheck = .false.
     endif
  else
     if(LIS_rc%da.eq.1) then 
        alarmCheck = .true. 
     else
        alarmCheck = .false. 
     endif     
  endif
end subroutine checkSimGRACEJPLdataAlarm

subroutine checkStatus_simGRACEJPLfile(odir,&
     config,file_status)

  use LIS_coreMod, only : LIS_rc

  character(len=*)         :: odir
  character(len=*)         :: config
  logical                  :: file_status

  character(len=LIS_CONST_PATH_LEN) :: name
  integer                  :: da
  character (len=4) :: fyr
  character (len=2) :: fmo
  character (len=2) :: fda


  write(unit=fyr, fmt='(i4.4)') LIS_rc%yr
  write(unit=fmo, fmt='(i2.2)') LIS_rc%mo
  write(unit=fda, fmt='(i2.2)') LIS_rc%da
  
  if(config.eq."GRACE-2") then 
     name = trim(odir)//'/GRACE_obs_'//trim(fyr)//trim(fmo)//trim(fda)//&
          '.bin'
  else
     name = trim(odir)//'/GRACE_obs_'//trim(fyr)//trim(fmo)//&
          '.bin'
  endif
  if(config.eq."GRACE-2") then 
!ideally we should be checking for each file, but in this case we know the files exist. 
     file_status = .true. 
  else
     inquire(file=name,exist=file_status)

  endif

end subroutine checkStatus_simGRACEJPLfile
