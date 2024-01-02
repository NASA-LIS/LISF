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
! !ROUTINE: read_multisynsmobs
! \label{read_multisynsmobs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_multisynsmobs(n, OBS_State, OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_historyMod, only : LIS_readvar_gridded
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_logunit, LIS_verify
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the synthetic soil moisture observations 
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
  type(ESMF_Field)    :: sm1field, sm2field, sm3field, sm4field

  real,    pointer    :: sm1(:), sm2(:), sm3(:), sm4(:)
  integer             :: gid(LIS_rc%ngrid(n))
  integer             :: assimflag(LIS_rc%ngrid(n))

  character(len=LIS_CONST_PATH_LEN) :: smobsdir, name
  logical             :: data_update
  logical             :: file_exists


  logical             :: readflag
  integer             :: status

  integer             :: t


  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call multismobs_filename(name,smobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)

  inquire(file=name,exist=file_exists)

  if(file_exists) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif
  
  if (readflag) then 
     write(LIS_logunit,*)  'Reading syn data ',trim(name)
     
     call ESMF_StateGet(OBS_State,"Observation01",sm1field,&
          rc=status)
     call LIS_verify(status)

     call ESMF_StateGet(OBS_State,"Observation02",sm2field,&
          rc=status)
     call LIS_verify(status)

     call ESMF_StateGet(OBS_State,"Observation03",sm3field,&
          rc=status)
     call LIS_verify(status)

     call ESMF_StateGet(OBS_State,"Observation04",sm4field,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(sm1field,localDE=0,farrayPtr=sm1,rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(sm2field,localDE=0,farrayPtr=sm2,rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(sm3field,localDE=0,farrayPtr=sm3,rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(sm4field,localDE=0,farrayPtr=sm4,rc=status)
     call LIS_verify(status)
     

     open(90, file=trim(name),form='unformatted')
     call LIS_readvar_gridded(90,n,sm1)
     call LIS_readvar_gridded(90,n,sm2)
     call LIS_readvar_gridded(90,n,sm3)
     call LIS_readvar_gridded(90,n,sm4)
     close(90)

     readflag = .false.
     
     do t=1,LIS_rc%ngrid(n)
        gid(t) = t
        if(sm1(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo

     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(sm1field,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(sm1field,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)

     do t=1,LIS_rc%ngrid(n)
        gid(t) = t
        if(sm2(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo

     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(sm2field,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(sm2field,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)

     do t=1,LIS_rc%ngrid(n)
        gid(t) = t
        if(sm3(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo
     
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(sm3field,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(sm3field,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)

     do t=1,LIS_rc%ngrid(n)
        gid(t) = t
        if(sm4(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo

     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(sm4field,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(sm4field,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  end if
  
end subroutine read_multisynsmobs

subroutine multismobs_filename(name, ndir, yr, mo,da,hr,mn)
  
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
  
  name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
       trim(fmn)//'.d01.gs4r'
end subroutine multismobs_filename



