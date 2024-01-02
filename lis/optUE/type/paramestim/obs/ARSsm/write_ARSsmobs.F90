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
! !ROUTINE: write_ARSsmobs
! \label{write_ARSsmobs}
!
! !REVISION HISTORY:
!  02 Feb 2018: Soni Yatheendradas; Initial Specification
!
! !INTERFACE: 
subroutine write_ARSsmobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_masterproc
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_fileIOMod,      only : LIS_create_output_directory
  use LIS_historyMod,     only : LIS_writevar_gridded

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  write the ARSsm soil moisture data to disk
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP
  real,    pointer    :: obsl(:)
  type(ESMF_Field)    :: smcField
!  type(ESMF_Field)    :: smstdField
  logical             :: data_update
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer             :: status
  integer             :: ftn
  integer             :: n 

  n = 1

  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(Obj_Space,"ARS_sm",smcField,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(smcField, localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call ARSsm_obsname('smc',obsname)

        call LIS_create_output_directory('PEOBS') 
        open(ftn,file=trim(obsname), form='unformatted')
     endif
     
     call LIS_writevar_gridded(ftn, n, obsl)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

!     call ESMF_StateGet(Obj_Space,"ARSsm standard deviation of soil moisture",smstdField,&
!          rc=status)
!     call LIS_verify(status)
!
!     call ESMF_FieldGet(smstdField, localDE=0,farrayPtr=obsl,rc=status)
!     call LIS_verify(status)
!
!     if(LIS_masterproc) then 
!        ftn = LIS_getNextUnitNumber()
!        call ARSsm_obsname('smstd',obsname)
!
!!        call LIS_create_output_directory('PEOBS') ! SY: already done further above 
!        open(ftn,file=trim(obsname), form='unformatted')
!     endif
!     
!     call LIS_writevar_gridded(ftn, n, obsl)
!     
!     if(LIS_masterproc) then 
!        call LIS_releaseUnitNumber(ftn)
!     endif

  endif

end subroutine write_ARSsmobs

!BOP
! 
! !ROUTINE: ARSsm_obsname
! \label{ARSsm_obsname}
! 
! !INTERFACE: 
subroutine ARSsm_obsname(variabname, obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  character(len=*)      :: obsname
  character(len=*)       :: variabname
! 
! !DESCRIPTION: 
!  This method generates a timestamped filename for the processed
!  ARSsm observations. 
! 
!EOP

  character(len=12)    :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
       LIS_rc%yr, LIS_rc%mo, &
       LIS_rc%da, LIS_rc%hr,LIS_rc%mn

  obsname = trim(LIS_rc%odir)//'/PEOBS/'//trim(variabname)//'_'//cdate1//'.1gs4r'
  
end subroutine ARSsm_obsname



