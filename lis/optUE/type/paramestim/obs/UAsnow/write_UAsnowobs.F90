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
! !ROUTINE: write_UAsnowobs
! \label{write_UAsnowobs}
!
! !REVISION HISTORY:
!  2 May 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine write_UAsnowobs(Obj_Space)
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
!  write the UA snow data used in parameter estimation to disk
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP
  real,    pointer    :: snod(:), swe(:)
  type(ESMF_Field)    :: snodField, sweField
  logical             :: data_update
  character(len=LIS_CONST_PATH_LEN) :: snod_filename, swe_filename 
  integer             :: status
  integer             :: snod_ftn, swe_ftn
  integer             :: n 

  n = 1

  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status, "Error in ESMF_AttributeGet: Data Update Status (UA Snow)")

  if(data_update) then 
     
     call ESMF_StateGet(Obj_Space,"UA_SNOD", snodField, rc=status)
     call LIS_verify(status, "Error in ESMF_StateGet: UA_SNOD")

     call ESMF_FieldGet(snodField, localDE=0, farrayPtr=snod, rc=status)
     call LIS_verify(status, "Error in ESMF_FieldGet: snodField")

     call ESMF_StateGet(Obj_Space, "UA_SWE", sweField, rc=status)
     call LIS_verify(status, "Error in ESMF_StateGet: UA_SWE")

     call ESMF_FieldGet(sweField, localDE=0, farrayPtr=swe, rc=status)
     call LIS_verify(status, "Error in ESMF_FieldGet: sweField")

     if(LIS_masterproc) then 
        snod_ftn = LIS_getNextUnitNumber()
        swe_ftn = LIS_getNextUnitNumber()

        call UAsnow_obsname('UA_SNOD', snod_filename)
        call UAsnow_obsname('UA_SWE', swe_filename)

        call LIS_create_output_directory('PEOBS')

        open(snod_ftn, file=trim(snod_filename), form='unformatted')
        open(swe_ftn, file=trim(swe_filename), form='unformatted')
     endif
     
     call LIS_writevar_gridded(snod_ftn, n, snod)
     call LIS_writevar_gridded(swe_ftn, n, swe)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(snod_ftn)
        call LIS_releaseUnitNumber(swe_ftn)
     endif

  endif

end subroutine write_UAsnowobs

!BOP
! 
! !ROUTINE: UAsnow_obsname
! \label{UAsnow_obsname}
! 
! !INTERFACE: 
subroutine UAsnow_obsname(variabname, obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  character(len=*)       :: obsname
  character(len=*)       :: variabname
! 
! !DESCRIPTION: 
!  This method generates a timestamped filename for the processed
!  UAsnow observations. 
! 
!EOP

  character(len=12)    :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
       LIS_rc%yr, LIS_rc%mo, &
       LIS_rc%da, LIS_rc%hr,LIS_rc%mn

  obsname = trim(LIS_rc%odir)//'/PEOBS/'//trim(variabname)//'_'//cdate1//'.1gs4r'
  
end subroutine UAsnow_obsname



