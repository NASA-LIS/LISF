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
! !ROUTINE: write_GRACEobs
! \label{write_GRACEobs}
! 
! !REVISION HISTORY: 
! 25Jan2008: Sujay Kumar; Initial Specification
! 28Sep2011: Ben Zaitchik; modified for GRACE
!
! !INTERFACE: 
subroutine write_GRACEobs(n, k, OBS_State)
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
! GRACE observations to a file
! 
!EOP
  type(ESMF_Field)         :: twsField
  logical                  :: data_update
  real, pointer            :: twsobs(:)
  real                     :: twsobs_unsc(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer                  :: ftn
  integer                  :: status

  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01",twsField, &
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(twsField, localDE=0, farrayPtr=twsobs, rc=status)
     call LIS_verify(status)

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call GRACE_twsobsname(obsname)        

        call LIS_create_output_directory('DAOBS')
        open(ftn,file=trim(obsname), form='unformatted')
     endif

     call LIS_writevar_gridded_obs(ftn,n,k,twsobs)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  endif  

end subroutine write_GRACEobs

!BOP
! !ROUTINE: GRACE_twsobsname
! \label{GRACE_twsobsname}
! 
! !INTERFACE: 
subroutine GRACE_twsobsname(obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  character(len=*)      :: obsname
! 
! !DESCRIPTION: 
! 
!EOP

  character(len=6) :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2)') &
       LIS_rc%yr, LIS_rc%mo

  obsname = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//'/'//cdate1//   &
            '.1gs4r'
  
end subroutine GRACE_twsobsname
