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
! !ROUTINE: read_syntheticsndobs
!  \label{read_syntheticsndobs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_syntheticsndobs(n, k, OBS_State, OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_historyMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
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
!  reads the synthetic SND observations produced from a LIS control run. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: sndField

  integer             :: ftn
  integer             :: c,r
  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))

  integer             :: snodid
  real                :: snodobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  character(len=LIS_CONST_PATH_LEN) :: sndobsdir
  logical             :: data_update
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name

  logical             :: readflag
  integer             :: status

  integer             :: t


  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sndobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call synsnd_filename(name,sndobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)

  inquire(file=name,exist=file_exists)

  if(file_exists) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif

  if (readflag) then 

     call ESMF_StateGet(OBS_State,"Observation01",sndField,&
          rc=status)
     call LIS_verify(status, 'ESMF_StateGet failed in read_syntheticsndobs')
     call ESMF_FieldGet(sndField,localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed in read_syntheticsndobs')
     
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     write(LIS_logunit,*)  '[INFO] Reading syn data ',trim(name)
     
     call LIS_verify(nf90_open(path=trim(name),mode=NF90_NOWRITE,ncid=ftn),&
          'Error opening file '//trim(name))
     call LIS_verify(nf90_inq_varid(ftn,'SnowDepth_tavg',snodid),&
          'Error nf90_inq_varid: SnowDepth_tavg')
     
     call LIS_verify(nf90_get_var(ftn,snodid,snodobs,&
          start=(/LIS_ews_obs_halo_ind(n,LIS_localPet+1),&         
          LIS_nss_obs_halo_ind(n,LIS_localPet+1)/),&
          count = (/LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)/)),&
          'Error in nf90_get_var')
     call LIS_verify(nf90_close(ftn))

     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = &
                   snodobs(c,r)
           end if
        end do
     end do

#endif
     readflag = .false.

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
     call LIS_verify(status,'ESMF_AttributeGet failed in read_syntheticsndobs')

     if(LIS_rc%obs_ngrid(k).gt.0) then 
        call ESMF_AttributeSet(sndField,"Grid Number",&
             gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(sndField,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)
     endif

  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  end if
  
  do t=1,LIS_rc%obs_ngrid(k)
     if(obsl(t).ne.-9999.0) then 
        if(obsl(t).lt.0.0) obsl(t) = 0.0
!        if(obsl(t).gt.200.0 ) obsl(t) = 200.0
     endif
  enddo

end subroutine read_syntheticsndobs

subroutine synsnd_filename(name, ndir, yr, mo,da,hr,mn)
  
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
  
  name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/SimObs_'//&
       trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
       trim(fmn)//'.nc'
end subroutine synsnd_filename

!BOP
!
! !ROUTINE: readobsvar_1dgridded
! \label{readobsvar_1dgridded}
!
! !INTERFACE:
subroutine readobsvar_1dgridded_snd(ftn,n,k,var)
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

end subroutine readobsvar_1dgridded_snd

