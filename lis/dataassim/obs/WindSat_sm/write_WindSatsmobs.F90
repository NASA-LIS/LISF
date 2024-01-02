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
! !ROUTINE: write_WindSatsmobs
! \label{write_WindSatsmobs}
! 
! !REVISION HISTORY: 
! 22 Dec 2009: Sujay Kumar; Initial Specification
! 
! !INTERFACE: 
subroutine write_WindSatsmobs(n, OBS_State)
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
! WindSat observations to a file
! 
!EOP
  type(ESMF_Field)         :: smField
  logical                  :: data_update
  real, pointer            :: smobs(:)
  real                     :: smobs_unsc(LIS_rc%ngrid(n))
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer                  :: ftn
  integer                  :: status

  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01",smField, &
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(smField, localDE=0, farrayPtr=smobs, rc=status)
     call LIS_verify(status)

     call ESMF_AttributeGet(smfield,"Unscaled Obs",smobs_unsc,&
          itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status, 'Error in AttributeGet-Unscaled Obs')

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call WindSat_smobsname(obsname)        

        call LIS_create_output_directory('DAOBS')
        open(ftn,file=trim(obsname), form='unformatted')
     endif

!     call LIS_writevar_gridded(ftn,n,smobs_unsc)     
     call LIS_writevar_gridded(ftn,n,smobs)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  endif  

end subroutine write_WindSatsmobs

!BOP
! !ROUTINE: WindSat_smobsname
! \label{WindSat_smobsname}
! 
! !INTERFACE: 
subroutine WindSat_smobsname(obsname)
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
  
end subroutine WindSat_smobsname
