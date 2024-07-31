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
#include "LIS_NetCDF_inc.h"
!BOP
! !ROUTINE: read_hydrowebWLobs
! \label{read_hydrowebWLobs}
!
! !REVISION HISTORY:
!  17 Jul 2019: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_hydrowebWLobs(n, k, OBS_State, OBS_Pert_state) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_historyMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use LIS_dataAssimMod
  use LIS_surfaceModelMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use hydrowebWLobs_module

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the hydroweb radar altimetry observations 
!  The processed data is packaged 
!  into an ESMF State for later use within the DA algorithm. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]                index of the nest
!  \item[k]                index of the data assimilation instance
!  \item[OBS\_State]       observations state
!  \item[OBS\_Pert\_State] observations perturbation state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: smfield

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real                :: obs_unsc(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: wlobsdir
  logical             :: data_update
  logical             :: file_exists
  integer             :: fnd
  integer             :: dim1Id, dim2Id, tid
  integer             :: wlid
  integer             :: ftn,p
  logical             :: readflag
  integer             :: status
  logical             :: data_upd_flag(LIS_npes)
  logical             :: data_upd_flag_local
  logical             :: data_upd
  real                :: wl_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                :: wlobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer             :: t,c,r
  integer             :: offset
  type(ESMF_Time)     :: startTime, currentTime
  real                :: timenow
  logical             :: alarmCheck


  obs_unsc = LIS_rc%udef
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       wlobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  inquire(file=wlobsdir,exist=file_exists)

  if(file_exists.and.hydroweb_wl_struc(n)%readflag) then 
     
     hydroweb_wl_struc(n)%readflag = .false. 

     write(LIS_logunit,*)  '[INFO] Reading Hydroweb WL data ',trim(wlobsdir)
     
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     call LIS_verify(nf90_open(path=trim(wlobsdir),mode=NF90_NOWRITE,ncid=ftn),&
          'Error opening file '//trim(wlobsdir))
     call LIS_verify(nf90_inq_dimid(ftn,'time',dim1Id),&
          'Error with nf90_inq_dimid: time')
     call LIS_verify(nf90_inq_dimid(ftn,'sites',dim2Id),&
          'Error with nf90_inq_dimid: sites')
     
     call LIS_verify(nf90_inquire_dimension(ftn,dim1Id,&
          len=hydroweb_wl_struc(n)%ntimes),&
          'Error with nf90_inquire_dimension: time')
     call LIS_verify(nf90_inquire_dimension(ftn,dim2Id,&
          len=hydroweb_wl_struc(n)%nsites),&
          'Error with nf90_inquire_dimension: sites')
     
     allocate(hydroweb_wl_struc(n)%WLobs(&
          hydroweb_wl_struc(n)%nsites,&
          hydroweb_wl_struc(n)%ntimes))

     allocate(hydroweb_wl_struc(n)%time(&
          hydroweb_wl_struc(n)%ntimes))

     call LIS_verify(nf90_inq_varid(ftn,'time',tId),&
          'Error with nf90_iniq_varid: time')

     call LIS_verify(nf90_get_var(ftn,tId,hydroweb_wl_struc(n)%time),&
          'Error with nf90_get_var: time')

     call LIS_verify(nf90_inq_varid(ftn,'water_elevation',wlId),&
          'Error with nf90_inq_varid: water_elevation')

     call LIS_verify(nf90_get_var(ftn,wlId,hydroweb_wl_struc(n)%WLobs),&
          'Error with nf90_get_var: water_elevation')
     
     call LIS_verify(nf90_close(ftn))
#endif
  endif

  timenow = float(LIS_rc%hr*3600 + 60*LIS_rc%mn +LIS_rc%ss)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  wlobs = -9999.0

  if(alarmCheck) then 
!-------------------------------------------------------------------------
!  Extract data for the current time
!-------------------------------------------------------------------------     
     call ESMF_StateGet(OBS_State,"Observation01",smfield,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     obsl = -1

     call ESMF_TimeSet(startTime, yy=2016,&
          mm = 7, dd=1, h =0, &
          m = 0, s = 0, calendar = LIS_calendar,&
          rc=status)
     call LIS_verify(status, 'ESMF_TimeSet failed in read_hydrowebWLobs')
     
     call ESMF_TimeSet(currentTime, yy=LIS_rc%yr,&
          mm = LIS_rc%mo, dd=LIS_rc%da, h =LIS_rc%hr, &
          m = LIS_rc%mn, s = 0, calendar = LIS_calendar,&
          rc=status)
     call LIS_verify(status, 'ESMF_TimeSet failed in read_hydrowebWLobs')
     
     offset = nint((currentTime - startTime)/hydroweb_wl_struc(n)%ts) + 1
     
     if(offset.gt.0) then 
        fnd = 0 
        do r =1,LIS_rc%obs_lnr(k)
           do c =1,LIS_rc%obs_lnc(k)
              if(hydroweb_wl_struc(n)%sites_data(c,r).gt.0) then 
                 wlobs(c,r) = hydroweb_wl_struc(n)%WLobs(&
                      nint(hydroweb_wl_struc(n)%sites_data(c,r)),&
                      offset)
                 if(wlobs(c,r).gt.0) then 
                    fnd = 1
                 endif
              endif
           enddo
        enddo
        

        if(LIS_rc%dascaloption(k).eq."Normal deviate scaling".and.fnd.ne.0) then

           call LIS_rescale_with_normal_deviate_scaling(&
                n,k,                               & 
                hydroweb_wl_struc(n)%nt,        & 
                hydroweb_wl_struc(n)%model_mu,  &
                hydroweb_wl_struc(n)%model_sigma,     &
                hydroweb_wl_struc(n)%obs_mu,    &
                hydroweb_wl_struc(n)%obs_sigma,       &
                wlobs)
        endif
        

        do r =1,LIS_rc%obs_lnr(k)
           do c =1,LIS_rc%obs_lnc(k)
              if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
                 obsl(LIS_obs_domain(n,k)%gindex(c,r)) = &
                      wlobs(c,r)
              end if
           end do
        end do
        
        call LIS_checkForValidObs(n, k,obsl,fnd,wl_current)

     
!-------------------------------------------------------------------------
!  Apply LSM-based QC of observations
!-------------------------------------------------------------------------     
     
        call LIS_surfaceModel_DAqcObsState(n,k)
        
        call LIS_checkForValidObs(n, k,obsl,fnd,wl_current)
        
        if(fnd.eq.0) then 
           data_upd_flag_local = .false.
        else
           data_upd_flag_local = .true.
        endif
        
#if (defined SPMD)
        call MPI_ALLGATHER(data_upd_flag_local,1, &
             MPI_LOGICAL, data_upd_flag(:),&
             1, MPI_LOGICAL, LIS_mpi_comm, status)
        data_upd = any(data_upd_flag)
#else
        data_upd = data_upd_flag_local
#endif
        
        if(data_upd) then 
           do t=1,LIS_rc%obs_ngrid(k)
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
           
           if(LIS_rc%obs_ngrid(k).gt.0) then 
              call ESMF_AttributeSet(smfield,"Grid Number",&
                   gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
              call LIS_verify(status)
              
              call ESMF_AttributeSet(smfield,"Assimilation Flag",&
                   assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
              call LIS_verify(status)
              
              call ESMF_AttributeSet(smfield, "Unscaled Obs",&
                   obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
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
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  endif
  
end subroutine read_hydrowebWLobs

