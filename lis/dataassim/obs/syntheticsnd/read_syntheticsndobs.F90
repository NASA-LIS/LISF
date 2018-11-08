!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
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

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real, allocatable       :: dummy(:)

  character*100       :: sndobsdir
  logical             :: data_update
  logical             :: file_exists
  character*80        :: name

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
     allocate(dummy(LIS_rc%obs_ngrid(k)))
     write(LIS_logunit,*)  'Reading syn data ',name
     
     call ESMF_StateGet(OBS_State,"Observation01",sndField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(sndField,localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
     open(90, file=trim(name),form='unformatted')
     do t=1,1
        if(t==1) then 
           call readobsvar_1dgridded_snd(90,n,k,obsl)  !Yeosang Yoon
           !call LIS_readvar_gridded(90,n,obsl)
        else 
           call readobsvar_1dgridded_snd(90,n,k,obsl)
           !call LIS_readvar_gridded(90,n,dummy)
        endif
     end do
     close(90)
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
     call LIS_verify(status)

     call ESMF_AttributeSet(sndField,"Grid Number",&
          gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(sndField,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)

     deallocate(dummy)
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
  character*80      :: name
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr,fmn
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  write(unit=fmn, fmt='(i2.2)') mn  
  
  name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
       trim(fmn)//'.d01.gs4r'
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

