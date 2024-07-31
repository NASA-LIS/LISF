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
! 
! !ROUTINE: write_syntheticwlobs
! \label{write_syntheticwlobs}
! 
! !REVISION HISTORY: 
! 07 Nov 2019: Sujay Kumar; Initial Specification
! 
! !INTERFACE: 
subroutine write_syntheticwlobs(n, k, OBS_State)
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_historyMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use syntheticwlobs_module
  
  implicit none

! !ARGUMENTS: 

  integer,     intent(in)  :: n 
  integer,     intent(in)  :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION: 
! 
! writes the transformed (interpolated/upscaled/reprojected)  
! synthetic observations to a file
! 
!EOP
  type(ESMF_Field)         :: smField
  logical                  :: data_update
  real, pointer            :: wlobs(:)
  real                     :: wlobs_unsc(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer                  :: ftn
  integer                  :: status

  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Data Update Status failed in write_syntheticwlobs")

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01",smField, &
          rc=status)
     call LIS_verify(status, &
          "ESMF_StateGet failed in write_syntheticwlobs")
     
     call ESMF_FieldGet(smField, localDE=0, farrayPtr=wlobs, rc=status)
     call LIS_verify(status,&
          "ESMF_FieldGet failed in write_syntheticwlobs")

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call synthetic_wlobsname(n,k,obsname)        

        call LIS_create_output_directory('DAOBS')
        open(ftn,file=trim(obsname), form='unformatted')
     endif

     if(LIS_rc%dascaloption(k).ne."none") then 
        call LIS_writevar_gridded_obs(ftn,n,k,wlobs_unsc)
     endif

     call LIS_writevar_gridded_obs(ftn,n,k,wlobs)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  endif  

end subroutine write_syntheticwlobs

!BOP
! !ROUTINE: synthetic_wlobsname
! \label{synthetic_wlobsname}
! 
! !INTERFACE: 
subroutine synthetic_wlobsname(n,k,obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  integer               :: n
  integer               :: k
  character(len=*)      :: obsname
! 
! !DESCRIPTION: 
!  This subroutine creates the name of the synthetic observation
!  data file
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

end subroutine synthetic_wlobsname
