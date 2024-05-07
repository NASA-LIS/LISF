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
! !ROUTINE: read_USDA_ARSsm_obsdata
! \label{read_USDA_ARSsm_obsdata}
!
! !REVISION HISTORY:
!  11 Jul 2011: Ken Harrison; Initial Specification
!
! !INTERFACE: 
  subroutine read_USDA_ARSsm_obsdata(Obj_Space)
! !USES: 
    use ESMF
    use LIS_mpiMod
    use LIS_coreMod,  only : LIS_rc, LIS_domain, LIS_npes, LIS_localPet
    use LIS_timeMgrMod, only : LIS_calendar
    use LIS_logMod,     only : LIS_logunit, LIS_verify, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_fileIOMod,      only : LIS_readData
    use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
    use USDA_ARSsm_obsMod, only : USDA_ARSsm_obs_struc
    use map_utils

    implicit none
! !ARGUMENTS: 
    type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP
    !    real,    pointer    :: obse(:,:)
    real,    pointer    :: obse(:)
    type(ESMF_Field)    :: smField
    character(len=LIS_CONST_PATH_LEN) :: smobsdir 
    logical             :: data_update
    integer             :: status 
    logical             :: found
    type(ESMF_TimeInterval) :: delta_t
    type(ESMF_Time)         :: lis_time
    integer             :: t
    integer             :: n 
    integer             :: day_index
    integer :: i
    n = 1

    !initialize some values
    found=.false.
    call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
         .false., rc=status)
    call LIS_verify(status)

    if (LIS_rc%hr.eq.0) then
       call ESMF_AttributeGet(Obj_Space,"Data Directory",&
            smobsdir, rc=status)
       call LIS_verify(status)
       call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
            data_update, rc=status)
       call LIS_verify(status)

       allocate(obse(LIS_rc%ngrid(n)))

       call ESMF_TimeSet(lis_time, yy=LIS_rc%yr, mm=LIS_rc%mo, &
            dd=LIS_rc%da,h=LIS_rc%hr,&
            m=LIS_rc%mn,calendar=LIS_calendar,rc=status)

       delta_t = lis_time - USDA_ARSsm_obs_struc(n)%start_time
       call ESMF_TimeIntervalGet(delta_t, d=day_index, &
            calendar=LIS_calendar,rc=status)
       call LIS_verify(status, 'error in timeget: read_CNRSdata.F90')

       !as fortran arrays are zero-based...
       day_index=day_index+1

       call ESMF_StateGet(Obj_Space,"USDA ARS Surface Soil Moisture",smField,&
            rc=status)
       call LIS_verify(status)

       call ESMF_FieldGet(smField,localDE=0,farrayPtr=obse,rc=status)
       call LIS_verify(status)

       obse=USDA_ARSsm_obs_struc(n)%soilmoisture(day_index)
       if(obse(1).ne.LIS_rc%udef) then
          found=.true.
       endif
    endif

    if(found) then      
       write(LIS_logunit,*) 'Read soilmoisture data'
       call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
            .true., rc=status)
       call LIS_verify(status)
    endif

end subroutine read_USDA_ARSsm_obsdata


