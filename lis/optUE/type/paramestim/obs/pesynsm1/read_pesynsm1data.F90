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
! !ROUTINE: read_pesynsm1data
! \label{read_pesynsm1data}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_pesynsm1data(Obj_Space)
! !USES: 
  use ESMF
  use LIS_historyMod, only : LIS_readvar_gridded
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_logunit, LIS_verify
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
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
  type(ESMF_Field)    :: smfield
  real,    pointer    :: obsl(:)
  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  logical             :: data_update
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name


  logical             :: readflag
  integer             :: status

  integer             :: t
  integer             :: n 

  n = 1
  call ESMF_AttributeGet(Obj_Space,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call pesynsm1_filename(name,smobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)

  inquire(file=name,exist=file_exists)

  if(file_exists) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif
  
  if (readflag) then 
     write(LIS_logunit,*)  'Reading syn data ',trim(name)
     
     call ESMF_StateGet(Obj_Space,"Synthetic soil moisture",smfield,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     

     open(90, file=trim(name),form='unformatted')
     call LIS_readvar_gridded(90,n,obsl)
     close(90)

     readflag = .false.
     
     call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status)

  else
     call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  end if
  
  do t=1,LIS_rc%ngrid(n)
     if(obsl(t).ne.-9999.0) then 
        if(obsl(t).lt.0.0) obsl(t) = 0.0
        if(obsl(t).gt.100.0 ) obsl(t) = 100.0
     endif
  enddo
end subroutine read_pesynsm1data

subroutine pesynsm1_filename(name, ndir, yr, mo,da,hr,mn)
  
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
end subroutine pesynsm1_filename



