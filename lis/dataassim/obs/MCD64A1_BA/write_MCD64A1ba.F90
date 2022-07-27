!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: write_MCD64A1ba
! \label{write_MCD64A1ba}
! 
! !REVISION HISTORY: 
!  24 Jul 2022  Sujay Kumar; initial specification
! 
! !INTERFACE: 
subroutine write_MCD64A1ba(n, k, OBS_State)
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_historyMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none

! !ARGUMENTS: 

  integer,     intent(in)  :: n 
  integer,     intent(in)  :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION: 
! 
! writes the transformed (interpolated/upscaled/reprojected)  
! MCD64A1 BA observations to a file
! 
!EOP
  type(ESMF_Field)         :: baField
  logical                  :: data_update
  real, pointer            :: baobs(:)
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer                  :: ftn
  integer                  :: status

  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01",baField, &
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(baField, localDE=0, farrayPtr=baobs, rc=status)
     call LIS_verify(status)

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call MCD64A1_baobsname(n,k,obsname)        

        call LIS_create_output_directory('DAOBS')
        open(ftn,file=trim(obsname), form='unformatted')
     endif

     call LIS_writevar_gridded_obs(ftn,n,k,baobs)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  endif  

end subroutine write_MCD64A1ba

!BOP
! !ROUTINE: MCD64A1_baobsname
! \label{MCD64A1_baobsname}
! 
! !INTERFACE: 
subroutine MCD64A1_baobsname(n,k,obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  integer               :: n
  integer               :: k
  character(len=*)      :: obsname
! 
! !DESCRIPTION: 
!
!  writes the assimilated observation to a file.

!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest
!  \item[k] number of observation state
!  \item[obsname] name of the observation
!  \end{description}
!
 
!EOP

  character(len=12) :: cdate1
  character(len=12) :: cdate
  character(len=10) :: cda

  write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
       LIS_rc%yr, LIS_rc%mo, &
       LIS_rc%da, LIS_rc%hr,LIS_rc%mn

  write(unit=cda, fmt='(a2,i2.2)') '.a',k
  write(unit=cdate, fmt='(a2,i2.2)') '.d',n

  obsname = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//&
       '/LISDAOBS_'//cdate1// &
       trim(cda)//trim(cdate)//'.1gs4r'

end subroutine MCD64A1_baobsname
