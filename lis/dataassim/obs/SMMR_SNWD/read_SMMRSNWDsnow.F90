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
! !ROUTINE: read_SMMRSNWDsnow
! \label{read_SMMRSNWDsnow}
!
! !REVISION HISTORY:
!  16 Oct 2012: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_SMMRSNWDsnow(n,k,OBS_State,OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_pluginIndices, only : LIS_SMMRSNWDsnowobsId
  use SMMRSNWDsnow_Mod, only : SMMRSNWDsnow_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  This routine reads and processes SMMR snow depth observations. 
!  The data is read at 0z every day and is kept in memory. At 
!  each timestep, a subset of data is chosen for use in DA if 
!  the local time of the grid point is 6AM (personal 
!  communication with George Riggs). 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[k]    index of the assimilation instance
!  \item[OBS\_State] observations state
!  \item[OBS\_Pert\_State] observations perturbation state
!  \end{description}
!
!EOP
  type(ESMF_Field)              :: snowField,pertfield
  logical                       :: alarmCheck
  logical                       :: data_upd, file_exists
  logical                       :: dataflag(LIS_npes)
  logical                       :: dataflag_local
  integer                       :: c,r, p, t
  character(len=LIS_CONST_PATH_LEN) :: obsdir, smmr_filename
  real, allocatable             :: snwd_field(:,:)
  real                          :: tsnow(SMMRSNWDsnow_struc(n)%nc*&
       SMMRSNWDsnow_struc(n)%nr)
  logical*1                     :: li(SMMRSNWDsnow_struc(n)%nc*&
       SMMRSNWDsnow_struc(n)%nr)
  integer                       :: ftn
  real                          :: lon, lhour
  real                          :: gmt
  real                          :: dt
  integer                       :: zone
  integer                       :: grid_index
  real                          :: ssdev(LIS_rc%obs_ngrid(k))
  logical*1                     :: lo(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real,             pointer     :: obsl(:)
  integer                       :: gid(LIS_rc%obs_ngrid(k))
  integer                       :: assimflag(LIS_rc%obs_ngrid(k))
  integer                       :: status, iret, ierr
  integer                       :: fnd
  real                          :: snwd_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status,'Error in AttributeGet: Data Directory')

!-------------------------------------------------------------------------
!   Read the data at 0z daily. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SMMR snow depth read alarm")

  if(alarmCheck.or.SMMRSNWDsnow_struc(n)%startMode) then 
     SMMRSNWDsnow_struc(n)%startMode = .false.
     
     SMMRSNWDsnow_struc(n)%snwd = LIS_rc%udef
     SMMRSNWDsnow_struc(n)%snwdtime = -1

     call SMMRsnow_filename3(smmr_filename,obsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)       

     inquire(file=smmr_filename,exist=file_exists)
     if(file_exists) then 

        write(LIS_logunit,*)  '[INFO] Reading SMMR SNWD data ',trim(smmr_filename)
        
        allocate(snwd_field(SMMRSNWDsnow_struc(n)%nc, &
             SMMRSNWDsnow_struc(n)%nr))
        
        ftn = LIS_getNextUnitNumber()
        open(ftn,file=smmr_filename,form='unformatted')
        read(ftn) snwd_field
        call LIS_releaseUnitNumber(ftn)

        do r=1,SMMRSNWDsnow_struc(n)%nr
           do c=1,SMMRSNWDsnow_struc(n)%nc 
              tsnow(c+(r-1)*SMMRSNWDsnow_struc(n)%nc) = snwd_field(c,r)
           enddo
        enddo

        li  = .false.
        do c=1,SMMRSNWDsnow_struc(n)%mi
!-------------------------------------------------------------------------
! assume that 10mm is the threshold of detecting snow for passive microwave
! sensors
!-------------------------------------------------------------------------
           if(tsnow(c).ge.10) then 
              li(c) = .true. 
           endif
        enddo

!-------------------------------------------------------------------------
! use neighbor search approach to interpolate the SMMR data to the 
! observation grid used in assimilation. 
!-------------------------------------------------------------------------
        call neighbor_interp(LIS_rc%obs_gridDesc(k,:),li,tsnow,&
             lo,SMMRSNWDsnow_struc(n)%snwd,&
             SMMRSNWDsnow_struc(n)%mi,LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
             SMMRSNWDsnow_struc(n)%rlat,SMMRSNWDsnow_struc(n)%rlon, &
             SMMRSNWDsnow_struc(n)%n11,LIS_rc%udef,iret)

        deallocate(snwd_field)
!-------------------------------------------------------------------------
! Store the GMT corresponding to 6AM localtime at each grid point
!-------------------------------------------------------------------------
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                 grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
                 lon = LIS_obs_domain(n,k)%lon(grid_index)
                 
                 lhour = 6.0
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 SMMRSNWDsnow_struc(n)%snwdtime(c,r) = gmt

              endif
           enddo
        enddo

     endif
  endif


  call ESMF_StateGet(OBS_State,"Observation01",snowfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(snowfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
  obsl = LIS_rc%udef 

!-------------------------------------------------------------------------
!  Update the OBS_State by subsetting to the local grid time  
!-------------------------------------------------------------------------     

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
           
           dt = (LIS_rc%gmt - SMMRSNWDsnow_struc(n)%snwdtime(c,r))*3600.0
           lon = LIS_obs_domain(n,k)%lon(grid_index)

           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = & 
                   SMMRSNWDsnow_struc(n)%snwd(grid_index)
           endif
           
        endif
     enddo
  enddo

  dataflag_local = .false. 

!-------------------------------------------------------------------------
!  Apply LSM based quality control and screening of observations
!-------------------------------------------------------------------------     

  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_SMMRSNWDsnowobsId)//char(0), & 
       n, k,OBS_state)

  snwd_current = LIS_rc%udef
  call LIS_checkForValidObs(n,k,obsl,fnd,snwd_current)

  if(fnd.eq.0) then 
     dataflag_local = .false. 
  else
     dataflag_local = .true. 
  endif
 
#if (defined SPMD)
  call MPI_ALLGATHER(dataflag_local,1, MPI_LOGICAL, dataflag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, ierr)
#endif
  data_upd = .false.
  
  do p=1,LIS_npes
     data_upd = data_upd.or.dataflag(p)
  enddo

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
     call LIS_verify(status, 'Error: AttributeSet in Data Update Status')
     

     call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
          rc=status)
     call LIS_verify(status, 'ESMF_StateGet for Observation01 for OBS_Pert_State failed in read_SMMRSNWDsnow')
     
     if(LIS_rc%obs_ngrid(k).gt.0) then 

!linearly scale the observation err
        ssdev = SMMRSNWDsnow_struc(n)%ssdev 
        do t=1,LIS_rc%obs_ngrid(k)
           if(obsl(t).ne.-9999.0) then 
              ssdev(t) =  SMMRSNWDsnow_struc(n)%ssdev !+ 0.05*obsl(t)
           endif
        enddo

        call ESMF_AttributeSet(pertField,"Standard Deviation",&
             ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(snowfield,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status,'Error: AttributeSet in Grid Number')
        
        call ESMF_AttributeSet(snowfield,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')
        
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status, "Error: AttributeSet Data Update Status")
     return
  end if
  
end subroutine read_SMMRSNWDsnow

!BOP
!
! !ROUTINE: SMMRsnow_filename3
! \label{SMMRsnow_filename3}
! 
! !INTERFACE: 
subroutine SMMRsnow_filename3(name, ndir, yr, mo,da)
  
  implicit none
! !ARGUMENTS: 
  character(len=*)      :: name
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped SMMR filename
!  
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the SMMR filename
!  \item[ndir] name of the SMMR root directory
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

  name = trim(ndir)//'/'//trim(fyr)//'/SMMR_F08_LL_'//trim(fyr)//trim(fmo)//trim(fda)//'_SD.bin'
    
end subroutine SMMRsnow_filename3



