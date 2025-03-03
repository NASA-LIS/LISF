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
! !ROUTINE: write_SMAPsmobs
! \label{write_SMAPsmobs}
!
! !REVISION HISTORY:
!  21 Sep 2018: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine write_SMAPsmobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_masterproc
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_fileIOMod,      only : LIS_create_output_directory
  use LIS_historyMod,     only : LIS_writevar_gridded
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  write the SMAP soil moisture data to disk
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
     
     call ESMF_StateGet(Obj_Space,"SMAP_sm",smcField,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(smcField, localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call SMAPsm_obsname('smc',obsname)

        call LIS_create_output_directory('PEOBS') 
        open(ftn,file=trim(obsname), form='unformatted')
     endif
     
     call LIS_writevar_gridded(ftn, n, obsl)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  endif

end subroutine write_SMAPsmobs

!BOP
! 
! !ROUTINE: SMAPsm_obsname
! \label{SMAPsm_obsname}
! 
! !INTERFACE: 
subroutine SMAPsm_obsname(variabname, obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  character(len=*)      :: obsname
  character(len=*)       :: variabname
! 
! !DESCRIPTION: 
!  This method generates a timestamped filename for the processed
!  SMAPsm observations. 
! 
!EOP

  character(len=12)    :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
       LIS_rc%yr, LIS_rc%mo, &
       LIS_rc%da, LIS_rc%hr,LIS_rc%mn

  obsname = trim(LIS_rc%odir)//'/PEOBS/'//trim(variabname)//'_'//cdate1//'.1gs4r'
  
end subroutine SMAPsm_obsname



