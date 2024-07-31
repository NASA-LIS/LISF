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
! !ROUTINE: read_syntheticsweobs
!  \label{read_syntheticsweobs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_syntheticsweobs(n, k, OBS_State, OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_historyMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the synthetic SWE observations produced from a LIS control run. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: sweField

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real, allocatable       :: dummy(:)

  character(len=LIS_CONST_PATH_LEN) :: sweobsdir
  logical             :: data_update
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name

  logical             :: readflag
  integer             :: status

  integer             :: t


  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sweobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call synswe_filename(name,sweobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)

  inquire(file=name,exist=file_exists)

  if(file_exists) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif

  if (readflag) then 
     allocate(dummy(LIS_rc%obs_ngrid(k)))
     write(LIS_logunit,*)  'Reading syn data ',trim(name)
     
     call ESMF_StateGet(OBS_State,"Observation01",sweField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(sweField,localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
     open(90, file=trim(name),form='unformatted')
     do t=1,1
        if(t==1) then 
           call LIS_readvar_gridded(90,n,obsl)
        else 
           call LIS_readvar_gridded(90,n,dummy)
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

     call ESMF_AttributeSet(sweField,"Grid Number",&
          gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(sweField,"Assimilation Flag",&
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

end subroutine read_syntheticsweobs

subroutine synswe_filename(name, ndir, yr, mo,da,hr,mn)
  
  implicit none
  character(len=LIS_CONST_PATH_LEN) :: name
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
end subroutine synswe_filename



