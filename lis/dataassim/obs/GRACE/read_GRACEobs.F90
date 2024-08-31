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
! !ROUTINE: read_GRACEobs
! \label{read_GRACEobs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!  11Aug2011: Ben Zaitchik; Modified for GRACE
!
! !INTERFACE: 
subroutine read_GRACEobs(n,k, OBS_State, OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use GRACEobs_module
  use LIS_fileIOMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

!  use smootherDA_runMod, only: smootherDA_increments_mode
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the GRACE TWS observations produced using a
!  LIS open loop run and GRACE TELLUS gridded files, 
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
  real                :: gracegrid(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                :: gracegrid_glb(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k))
  real                :: graceerr(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                :: graceerr_glb(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k))
  real                :: gridDesc(6)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))

  character(len=LIS_CONST_PATH_LEN) :: GRACEobsdir ! BMc, change 200 to LIS_CONST_PATH_LEN
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name ! BMc, change 200 to LIS_CONST_PATH_LEN

  integer             :: col,row
  logical             :: data_upd
  logical             :: data_upd_flag_local
  logical             :: data_upd_flag(LIS_npes)
  integer             :: fnd
  logical             :: readflag
  logical             :: alarmCheck
  integer             :: status,ftn
  real, allocatable       :: ssdev(:)
  integer             :: t,p,r,c
  integer             :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/  !BZ

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       GRACEobsdir, rc=status)
  call LIS_verify(status)
  
  call ESMF_AttributeSet(OBS_State,"Data Update Status",&
       .false., rc=status)
  call LIS_verify(status)
  
  call GRACE_filename(name,GRACEobsdir,LIS_rc%yr,LIS_rc%mo)
  
  inquire(file=name,exist=file_exists)
  if(file_exists) then 
     call ESMF_AttributeSet(OBS_State,"File Status",&
          .true., rc=status)
     call LIS_verify(status)
  else
     call ESMF_AttributeSet(OBS_State,"File Status",&
          .false., rc=status)
     call LIS_verify(status)
  endif

  if(LIS_rc%da.eq.1.and.&
       LIS_rc%hr.eq.GRACE_struc(n)%alarmHr.and.&
       LIS_rc%mn.eq.0.and.&
       LIS_rc%ss.eq.0) then 
     alarmCheck = .true. 
  else
     alarmCheck = .false.
  endif

  data_upd = .false. 

  if(alarmCheck) then 
     if(LIS_rc%DAincrMode(n).eq.1) then 
        call GRACE_filename(name,GRACEobsdir,LIS_rc%yr,LIS_rc%mo)
        
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
           call ESMF_AttributeSet(OBS_State,&
                name="Data averaging factor",&
                value=float(days(LIS_rc%mo)),rc=status)
           call LIS_verify(status)

           write(LIS_logunit,*)  '[INFO] Reading GRACE data ',trim(name)

           ftn = LIS_getNextUnitNumber()
           open(ftn, file=trim(name),form='unformatted')
           read(ftn) gracegrid_glb
           read(ftn) graceerr_glb
           call LIS_releaseUnitNumber(ftn)

           gracegrid = gracegrid_glb(&
                LIS_ews_obs_halo_ind(n,LIS_localPet+1):&         
                LIS_ewe_obs_halo_ind(n,LIS_localPet+1), &
                LIS_nss_obs_halo_ind(n,LIS_localPet+1): &
                LIS_nse_obs_halo_ind(n,LIS_localPet+1))

           graceerr = graceerr_glb(&
                LIS_ews_obs_halo_ind(n,LIS_localPet+1):&         
                LIS_ewe_obs_halo_ind(n,LIS_localPet+1), &
                LIS_nss_obs_halo_ind(n,LIS_localPet+1): &
                LIS_nse_obs_halo_ind(n,LIS_localPet+1))

           fnd = 0 
           data_upd_flag_local = .false. 
           do r=1,LIS_rc%obs_lnr(k)
              do c=1,LIS_rc%obs_lnc(k)
                 if(gracegrid(c,r).ne.LIS_rc%udef) then 
                    fnd = 1
                 endif
              enddo
           enddo

           if(GRACE_struc(n)%useDistErr.eq.1) then 
              if(LIS_rc%obs_ngrid(k).gt.0) then 
                 allocate(ssdev(LIS_rc%obs_ngrid(k)))
                 do r=1,LIS_rc%obs_lnr(k)
                    do c=1,LIS_rc%obs_lnc(k)
                       if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                          ssdev(LIS_obs_domain(n,k)%gindex(c,r)) = & 
                               graceerr(c,r)
                       endif
                    enddo
                 enddo

                 call ESMF_AttributeSet(OBS_Pert_State,"Standard Deviation",&
                      ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                 call LIS_verify(status, &
                      'Error with ESMF_AttributeGet in read_GRACEobs')
                 deallocate(ssdev)
              endif
           endif

           
           call ESMF_StateGet(OBS_State,"Observation01",twsfield,&
                rc=status)
           call LIS_verify(status, 'ESMF_StateGet failed in read_GRACEobs')

           call ESMF_FieldGet(twsfield,localDE=0,farrayPtr=obsl,rc=status)
           call LIS_verify(status,'ESMF_FieldGet failed in read_GRACEobs')

           obsl(:) = -9999.0

           if(fnd.eq.0) then 
              obsl = LIS_rc%udef
           else
              do r=1,LIS_rc%obs_lnr(k)
                 do c=1,LIS_rc%obs_lnc(k)
                    if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                       if(gracegrid(c,r).gt.0.0) then
                          obsl(LIS_obs_domain(n,k)%gindex(c,r)) = gracegrid(c,r)
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

           write(LIS_logunit,*) '[INFO] read GRACE data'

           readflag = .false.

           if(data_upd) then              
              do t=1,LIS_rc%obs_ngrid(k)
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

              if(LIS_rc%obs_ngrid(k).gt.0) then 
                 call ESMF_AttributeSet(twsfield,"Grid Number",&
                      gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                 call LIS_verify(status)

                 call ESMF_AttributeSet(twsfield,"Assimilation Flag",&
                      assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                 call LIS_verify(status)
              endif
           endif
        else
            ! Natt: Let's print it out when LIS cannot find the file
            write(LIS_logunit,*)  '[WARNING] GRACE data not found ',trim(name)
            
        endif
     endif
  end if
  
end subroutine read_GRACEobs

subroutine GRACE_filename(name, ndir, yr, mo)
  
  implicit none
  character(len=*)  :: name
  character (len=*) :: ndir
  integer           :: yr, mo

  character (len=4) :: fyr
  character (len=2) :: fmo
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  
  name = trim(ndir)//'/GRACE_obs_'//trim(fyr)//trim(fmo)//&
       '.bin'
end subroutine GRACE_filename



