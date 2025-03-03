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
! !ROUTINE: read_syntheticsmobs
! \label{read_syntheticsmobs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_syntheticsmobs(n, k, OBS_State, OBS_Pert_state) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_historyMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod
  use LIS_pluginIndices
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use syntheticsmobs_module
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
!  reads the synthetic soil moisture observations 
!  produced from a LIS control run, applies bias correction, 
!  and quality control. The processed data is packaged 
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
  real,  parameter    :: MAX_SM_VALUE=0.45, MIN_SM_VALUE=0.0001
  type(ESMF_Field)    :: smfield

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real                :: obs_unsc(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  logical             :: data_update
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name
  integer             :: fnd
  integer             :: smid
  integer             :: ftn,p
  logical             :: readflag
  integer             :: status
  logical             :: data_upd_flag(LIS_npes)
  logical             :: data_upd_flag_local
  logical             :: data_upd
  real                :: sm_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                :: smobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer             :: t,c,r


  obs_unsc = LIS_rc%udef
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call synsm_filename(name,smobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)

  inquire(file=name,exist=file_exists)

  if(file_exists) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif
  
  if (readflag) then 
     write(LIS_logunit,*)  '[INFO] Reading syn data ',trim(name)
     
     call ESMF_StateGet(OBS_State,"Observation01",smfield,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     obsl = -1

!     open(ftn,file=name,form='unformatted')
!     call readobsvar_1dgridded(ftn,n,k,obsl)
!     call LIS_releaseUnitNumber(ftn)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     call LIS_verify(nf90_open(path=trim(name),mode=NF90_NOWRITE,ncid=ftn),&
          'Error opening file '//trim(name))
     call LIS_verify(nf90_inq_varid(ftn,'SoilMoist_tavg',smid),&
          'Error nf90_inq_varid: SoilMoist_tavg')
     
     call LIS_verify(nf90_get_var(ftn,smid,smobs),&
          'Error in nf90_get_var')
     call LIS_verify(nf90_close(ftn))
     
     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = &
                   smobs(c,r)
           end if
        end do
     end do

#endif
     
!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     
     obs_unsc = obsl

     call LIS_checkForValidObs(n, k,obsl,fnd,sm_current)
     
     if(LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.eq.1) then  
               
        call LIS_rescale_with_CDF_matching(    &
             n,                                   & 
             k,                                   & 
             synthetic_sm_struc(n)%nbins,         & 
             synthetic_sm_struc(n)%ntimes,        & 
             MAX_SM_VALUE,                        & 
             MIN_SM_VALUE,                        & 
             synthetic_sm_struc(n)%model_xrange,  &
             synthetic_sm_struc(n)%obs_xrange,    &
             synthetic_sm_struc(n)%model_cdf,     &
             synthetic_sm_struc(n)%obs_cdf,       &
             sm_current)
        
        obsl = LIS_rc%udef

        do r =1,LIS_rc%obs_lnr(k)
           do c =1,LIS_rc%obs_lnc(k)
              if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
                 obsl(LIS_obs_domain(n,k)%gindex(c,r))=sm_current(c,r)
              end if
           end do
        end do
     elseif(LIS_rc%dascaloption(k).eq."Linear scaling".and.fnd.eq.1) then  
        call LIS_rescale_with_linear_scaling(    &
             n,                                   & 
             k,                                   & 
             synthetic_sm_struc(n)%nbins,         & 
             synthetic_sm_struc(n)%ntimes,        & 
             synthetic_sm_struc(n)%obs_xrange,    &
             synthetic_sm_struc(n)%obs_cdf,       &
             sm_current)
        
        obsl = LIS_rc%udef

        do r =1,LIS_rc%obs_lnr(k)
           do c =1,LIS_rc%obs_lnc(k)
              if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
                 obsl(LIS_obs_domain(n,k)%gindex(c,r))=sm_current(c,r)
              end if
           end do
        end do
     endif

!-------------------------------------------------------------------------
!  Apply LSM-based QC of observations
!-------------------------------------------------------------------------     

     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
          //trim(LIS_synsmId)//char(0),n, k, OBS_state)

     call LIS_checkForValidObs(n, k,obsl,fnd,sm_current)

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
     write(LIS_logunit,*)  '[INFO] Finished reading syn data ',trim(name)
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  end if
  
end subroutine read_syntheticsmobs

!BOP
!
! !ROUTINE: synsm_filename
! \label{synsm_filename}
!
! !INTERFACE: 
subroutine synsm_filename(name, ndir, yr, mo,da,hr,mn)
  
  implicit none
! !ARGUMENTS: 
  character(len=*) :: name
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This routine creates the filename containing the synthetic soil 
!  moisture data
!
!EOP
  
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr,fmn
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  write(unit=fmn, fmt='(i2.2)') mn  
  
!  name = trim(ndir)//'/SOILM_'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
!       trim(fmn)//'.bin'

  name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/SimObs_'//&
       trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
       trim(fmn)//'.nc'

end subroutine synsm_filename

!BOP
!
! !ROUTINE: readobsvar_1dgridded
! \label{readobsvar_1dgridded}
!
! !INTERFACE: 
subroutine readobsvar_1dgridded(ftn,n,k,var)
! !USES: 
  use LIS_coreMod
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS: 
  integer              :: ftn
  integer              :: n
  integer              :: k
  real                 :: var(LIS_rc%obs_ngrid(k))
!
! !DESCRIPTION: 
!  This routine reads the observation data and subsets to the 
!  local processor's domain, in a 1-d vector formulation. 
!
!EOP

  real,  allocatable   :: gobs(:,:)
  integer              :: nc,c1,r1,c,r,gid

  allocate(gobs(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k)))
  read(ftn) gobs
  
  nc = (LIS_ewe_obs_halo_ind(k,LIS_localPet+1)-&
       LIS_ews_obs_halo_ind(k,LIS_localPet+1))+1
  
  do r=LIS_nss_obs_halo_ind(k,LIS_localPet+1),&
       LIS_nse_obs_halo_ind(k,LIS_localPet+1)
     do c=LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
          LIS_ewe_obs_halo_ind(k,LIS_localPet+1)
        c1 = c-LIS_ews_obs_halo_ind(k,LIS_localPet+1)+1
        r1 = r-LIS_nss_obs_halo_ind(k,LIS_localPet+1)+1
        gid = LIS_obs_domain(n,k)%gindex(c1,r1)
        if(gid.ne.-1) then
           var(gid) = gobs(c,r)
        endif
     enddo
  enddo
  deallocate(gobs)
  
end subroutine readobsvar_1dgridded

