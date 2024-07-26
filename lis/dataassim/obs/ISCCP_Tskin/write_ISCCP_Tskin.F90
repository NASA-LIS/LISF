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
! !ROUTINE: write_ISCCP_Tskin
! \label{write_ISCCP_Tskin}
!
! !REVISION HISTORY:
!  24Jan2008: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine write_ISCCP_Tskin(n, OBS_State)
! !USES: 
  use ESMF
  use LIS_coreMod,    only : LIS_masterproc, LIS_rc
  use LIS_logMod,     only : LIS_verify, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber
  use LIS_fileIOMod,  only : LIS_create_output_directory
  use LIS_historyMod, only : LIS_writevar_gridded
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
! 
! !DESCRIPTION: 
! 
!  writes the transformed (interpolated/upscaled/reprojected) ISCCP Tskin 
!  observations to a file. 
!  
!EOP
  type(ESMF_Field)          :: tskinfield
  logical                   :: data_update
  real, pointer             :: obsl(:)
  character(len=LIS_CONST_PATH_LEN) :: obsname 
  integer                   :: ftn, t
  integer                   :: status

  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01", tskinfield, &
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(tskinfield, localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status)

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call ISCCP_Tskin_obsname(obsname)

        call LIS_create_output_directory('DAOBS') 
        open(ftn,file=trim(obsname), form='unformatted')
     endif

     call LIS_writevar_gridded(ftn, n, obsl)

     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  end if

end subroutine write_ISCCP_Tskin

!BOP
! 
! !ROUTINE: ISCCP_Tskin_obsname
! \label{ISCCP_Tskin_obsname}
! 
! !INTERFACE: 
subroutine ISCCP_Tskin_obsname(obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  character(len=*)     :: obsname
! 
! !DESCRIPTION: 
! 
!EOP

  character(len=12)    :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
       LIS_rc%yr, LIS_rc%mo, &
       LIS_rc%da, LIS_rc%hr,LIS_rc%mn

  obsname = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//'/'//cdate1//   &
            '.1gs4r'
  
end subroutine ISCCP_Tskin_obsname
