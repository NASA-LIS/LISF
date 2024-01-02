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
! !ROUTINE: read_CNRS_em_obsdata
! \label{read_CNRS_em_obsdata}
!
! !REVISION HISTORY:
!  11 Jul 2011: Ken Harrison; Initial Specification
!
! !INTERFACE: 
  subroutine read_CNRS_em_obsdata(Obj_Space)
! !USES: 
    use ESMF
    use LIS_mpiMod
    use LIS_coreMod,  only : LIS_rc, LIS_domain, LIS_npes, LIS_localPet
    use LIS_timeMgrMod, only : LIS_calendar
    use LIS_logMod,     only : LIS_logunit, LIS_verify, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_fileIOMod,      only : LIS_readData
    use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
    use CNRS_em_obsMod, only : CNRS_em_obs_struc
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
    integer, parameter :: numchannels=7
    type(ESMF_Field)    :: emField
    character(len=LIS_CONST_PATH_LEN) :: emobsdir 
    logical             :: data_update
    integer             :: status 
    logical             :: found
    integer             :: c,r
    type(ESMF_TimeInterval) :: delta_t
    type(ESMF_Time)         :: lis_time1
    integer             :: t
    integer             :: n 
!    logical             :: is_ascend_pass_hr
    logical             :: ob_in_curr_hr
    integer             :: day_index
    logical             :: is_overpass_hr
    integer, parameter :: channel_focus=2 ! start with 19H
    integer            :: ipass
    integer :: i
    integer :: nobs_sum, platform_sum
    real    :: emiss_sum, temp, temp_nobs
    real    :: emiss_ave
    n = 1

    !initialize some values
    is_overpass_hr=.false.
    found=.false.
    call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
         .false., rc=status)
    call LIS_verify(status)


    if ((LIS_rc%hr .eq. CNRS_em_obs_struc(n)%overpass_hr_a) &
         .or.(LIS_rc%hr .eq. CNRS_em_obs_struc(n)%overpass_hr_d)) then
       is_overpass_hr=.true.
    endif

    if (is_overpass_hr) then
       call ESMF_AttributeGet(Obj_Space,"Data Directory",&
            emobsdir, rc=status)
       call LIS_verify(status)
       call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
            data_update, rc=status)
       call LIS_verify(status)

       !          allocate(emobs(LIS_rc%lnc(n),LIS_rc%lnr(n),numchannels))
       !          allocate(obse(LIS_rc%ngrid(n),numchannels))
       allocate(obse(LIS_rc%ngrid(n)))

       call ESMF_TimeSet(lis_time1, yy=LIS_rc%yr, mm=LIS_rc%mo, &
            dd=LIS_rc%da,h=LIS_rc%hr,&
            m=LIS_rc%mn,calendar=LIS_calendar,rc=status)

       delta_t = lis_time1 - CNRS_em_obs_struc(n)%start_time
       call ESMF_TimeIntervalGet(delta_t, d=day_index, &
            calendar=LIS_calendar,rc=status)
       call LIS_verify(status, 'error in timeget: read_CNRSdata.F90')

       !as fortran arrays are zero-based...
       day_index=day_index+1

       if (LIS_rc%hr .eq. CNRS_em_obs_struc(n)%overpass_hr_a) then
          ipass=CNRS_em_obs_struc(n)%index_ascend
       elseif (LIS_rc%hr .eq. CNRS_em_obs_struc(n)%overpass_hr_d) then
          ipass=CNRS_em_obs_struc(n)%index_descend
       endif

       call ESMF_StateGet(Obj_Space,"Emissivity",emField,&
            rc=status)
       call LIS_verify(status)

       call ESMF_FieldGet(emField,localDE=0,farrayPtr=obse,rc=status)
       call LIS_verify(status)

       obse=LIS_rc%udef 

       !First pass to count obs and platforms w/ min 3 obs
       nobs_sum=0
       platform_sum=0
       do i=1,3
          temp= CNRS_em_obs_struc(n)%emissivity &
               ( &
               day_index, &
               channel_focus, &
               ipass, &
               i & 
               )
          temp_nobs=CNRS_em_obs_struc(n)%nobs &
               ( &
               day_index, &
               ipass, &
               i &
               )
          if(temp.ne.LIS_rc%udef) then
             nobs_sum=nobs_sum+temp_nobs
             if(temp_nobs.ge.3) then
                platform_sum=platform_sum+1
             endif
          endif
       enddo

       !Second pass to compute average emiss if conditions met
       emiss_sum=0.0
       emiss_ave=LIS_rc%udef
       if(and(nobs_sum.ge.15,platform_sum.ge.2)) then
          do i=1,3
             temp= CNRS_em_obs_struc(n)%emissivity &
                  ( &
                  day_index, &
                  channel_focus, &
                  ipass, &
                  i & 
                  )
             temp_nobs=CNRS_em_obs_struc(n)%nobs &
                  ( &
                  day_index, &
                  ipass, &
                  i &
                  )
             if(temp.ne.LIS_rc%udef) then
                emiss_sum=emiss_sum+real(temp_nobs)*temp
             endif
          enddo
          emiss_ave=emiss_sum/nobs_sum
          found=.true.
       endif

       do t=1,LIS_rc%ngrid(n)
          obse(t)=emiss_ave
       enddo

    !          deallocate(obse)
       if(found) then      
          write(LIS_logunit,*) 'Read emissivity data'
          call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
               .true., rc=status)
          call LIS_verify(status)
       endif
    endif

end subroutine read_CNRS_em_obsdata


