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
! !ROUTINE: read_syntheticSnowTbObs
!  \label{read_syntheticSnowTbObs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!  29Sep2017: Yonghwan Kwon; modified for TB observations
!
! !INTERFACE: 
subroutine read_syntheticSnowTbObs(n, k, OBS_State, OBS_Pert_State) 
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
!  reads the synthetic TB observations produced from a LIS-SVM control run. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: TB18VField, TB18HField, TB36VField, TB36HField

  real,    pointer    :: obsl(:), obs2(:), obs3(:), obs4(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real, allocatable       :: dummy(:)

  character(len=LIS_CONST_PATH_LEN) :: TB18Vobsdir, TB18Hobsdir, TB36Vobsdir, TB36Hobsdir
  logical             :: data_update
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name18V, name18H, name36V, name36H
  character (len=3)   :: fr_channel

  logical             :: readflag
  integer             :: status

  integer             :: t


  call ESMF_AttributeGet(OBS_State,"Data Directory TB_18V",&
       TB18Vobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Directory TB_18H",&
       TB18Hobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Directory TB_36V",&
       TB36Vobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Directory TB_36H",&
       TB36Hobsdir, rc=status)
  call LIS_verify(status)

  write(LIS_logunit,*) 'TB18Vobsdir=', TB18Vobsdir
  write(LIS_logunit,*) 'TB18Hobsdir=', TB18Hobsdir
  write(LIS_logunit,*) 'TB36Vobsdir=', TB36Vobsdir
  write(LIS_logunit,*) 'TB36Hobsdir=', TB36Hobsdir

  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  fr_channel = '18V'
  call synTB_filename(name18V,TB18Vobsdir,fr_channel, &
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)
  fr_channel = '18H'
  call synTB_filename(name18H,TB18Hobsdir,fr_channel, &
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)
  fr_channel = '36V'
  call synTB_filename(name36V,TB36Vobsdir,fr_channel, &
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)
  fr_channel = '36H'
  call synTB_filename(name36H,TB36Hobsdir,fr_channel, &
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)

  write(LIS_logunit,*) 'obs_file_name_18V=', trim(name18V)
  write(LIS_logunit,*) 'obs_file_name_18H=', trim(name18H)
  write(LIS_logunit,*) 'obs_file_name_36V=', trim(name36V)
  write(LIS_logunit,*) 'obs_file_name_36H=', trim(name36H)

  inquire(file=name18V,exist=file_exists)

  if(file_exists) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif

  write(LIS_logunit,*) 'readflag_TB18V=', readflag

  if (readflag) then 
     allocate(dummy(LIS_rc%obs_ngrid(k)))

     !----------------------TB_18V----------------------------------
     write(LIS_logunit,*)  'Reading syn TB_18V data ',trim(name18V)
     
     call ESMF_StateGet(OBS_State,"Observation01",TB18VField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(TB18VField,localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
     open(90, file=trim(name18V),form='unformatted')
     do t=1,1
        if(t==1) then 
           call LIS_readvar_gridded(90,n,obsl)
        else 
           call LIS_readvar_gridded(90,n,dummy)
        endif
     end do
     close(90)
     !----------------------TB_18H----------------------------------
     write(LIS_logunit,*)  'Reading syn TB_18H data ',trim(name18H)

     call ESMF_StateGet(OBS_State,"Observation02",TB18HField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(TB18HField,localDE=0, farrayPtr=obs2,rc=status)
     call LIS_verify(status)

     open(90, file=trim(name18H),form='unformatted')
     do t=1,1
        if(t==1) then
           call LIS_readvar_gridded(90,n,obs2)
        else
           call LIS_readvar_gridded(90,n,dummy)
        endif
     end do
     close(90)
     !----------------------TB_36V----------------------------------
     write(LIS_logunit,*)  'Reading syn TB_36V data ',trim(name36V)

     call ESMF_StateGet(OBS_State,"Observation03",TB36VField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(TB36VField,localDE=0, farrayPtr=obs3,rc=status)
     call LIS_verify(status)

     open(90, file=trim(name36V),form='unformatted')
     do t=1,1
        if(t==1) then
           call LIS_readvar_gridded(90,n,obs3)
        else
           call LIS_readvar_gridded(90,n,dummy)
        endif
     end do
     close(90)
     !----------------------TB_36H----------------------------------
     write(LIS_logunit,*)  'Reading syn TB_36H data ',trim(name36H)

     call ESMF_StateGet(OBS_State,"Observation04",TB36HField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(TB36HField,localDE=0, farrayPtr=obs4,rc=status)
     call LIS_verify(status)

     open(90, file=trim(name36H),form='unformatted')
     do t=1,1
        if(t==1) then
           call LIS_readvar_gridded(90,n,obs4)
        else
           call LIS_readvar_gridded(90,n,dummy)
        endif
     end do
     close(90)
     !-------------------------------------------------------------

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

     !-----------------------TB_18V--------------------------
     call ESMF_AttributeSet(TB18VField,"Grid Number",&
          gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(TB18VField,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)
     !-----------------------TB_18H--------------------------
     call ESMF_AttributeSet(TB18HField,"Grid Number",&
          gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(TB18HField,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)
     !-----------------------TB_36V--------------------------
     call ESMF_AttributeSet(TB36VField,"Grid Number",&
          gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(TB36VField,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)
     !-----------------------TB_36H--------------------------
     call ESMF_AttributeSet(TB36HField,"Grid Number",&
          gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(TB36HField,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
     call LIS_verify(status)
     !-------------------------------------------------------

     deallocate(dummy)
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  end if
  
  do t=1,LIS_rc%obs_ngrid(k)
     if(obsl(t).ne.-9999.0) then 
        if(obsl(t).lt.0.0) obsl(t) = -9999.0
!        if(obsl(t).gt.200.0 ) obsl(t) = 200.0
     endif
     if(obs2(t).ne.-9999.0) then
        if(obs2(t).lt.0.0) obs2(t) = -9999.0
     endif
     if(obs3(t).ne.-9999.0) then
        if(obs3(t).lt.0.0) obs3(t) = -9999.0
     endif
     if(obs4(t).ne.-9999.0) then
        if(obs4(t).lt.0.0) obs4(t) = -9999.0
     endif
  enddo

end subroutine read_syntheticSnowTbObs

subroutine synTB_filename(name, ndir, fr_channel, yr, mo,da,hr,mn)

  use LIS_logMod
  
  implicit none
  character(len=*)  :: name
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr,fmn
  character (len=3) :: fr_channel 
 
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  write(unit=fmn, fmt='(i2.2)') mn  
  
  name = trim(ndir)//'/TB_'//trim(fr_channel)//'/'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
       trim(fmn)//'.d01.gs4r'

  write(LIS_logunit,*) 'synTB_filename_name=', name

end subroutine synTB_filename



