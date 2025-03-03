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
! !ROUTINE: write_AMSRE_SWEobs
! \label{write_AMSRE_SWEobs}
! 
! !REVISION HISTORY: 
! 01 Jul 2010: Sujay Kumar; Initial Specification
! 
! !INTERFACE: 
subroutine write_AMSRE_SWEobs(n, OBS_State)
! !USES: 
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_logMod,     only : LIS_verify, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber
  use LIS_fileIOMod,  only : LIS_create_output_directory
  use LIS_historyMod, only : LIS_writevar_gridded
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none

! !ARGUMENTS: 

  integer,     intent(in)  :: n 
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION: 
! 
! writes the transformed (interpolated/upscaled/reprojected)  
! AMSRE SWE observations to a file
! 
!EOP
  type(ESMF_Field)         :: sweField
  logical                  :: data_update
  real, pointer            :: sweobs(:)
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer                  :: ftn
  integer                  :: status

  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01",sweField, &
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(sweField, localDE=0, farrayPtr=sweobs, rc=status)
     call LIS_verify(status)

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call AMSRE_sweobsname(obsname)        

        call LIS_create_output_directory('DAOBS')
        open(ftn,file=trim(obsname), form='unformatted')
     endif

     call LIS_writevar_gridded(ftn,n,sweobs)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  endif  
end subroutine write_AMSRE_SWEobs

!BOP
! !ROUTINE: AMSRE_SWEobsname
! \label{AMSRE_SWEobsname}
! 
! !INTERFACE: 
subroutine AMSRE_SWEobsname(obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  character(len=*)      :: obsname
! 
! !DESCRIPTION: 
! 
!EOP

  character(len=12) :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
       LIS_rc%yr, LIS_rc%mo, &
       LIS_rc%da, LIS_rc%hr,LIS_rc%mn

  obsname = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//'/'//cdate1//   &
            '.1gs4r'
  
end subroutine AMSRE_SWEobsname
