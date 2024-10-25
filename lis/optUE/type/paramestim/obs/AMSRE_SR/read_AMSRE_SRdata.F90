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
! !ROUTINE: read_AMSRE_SR_em_obsdata
! \label{read_AMSRE_SR_em_obsdata}
!
! !REVISION HISTORY:
!  31 Jan 2012: Ken Harrison; Initial Specification
!
! !INTERFACE: 
  subroutine read_AMSRE_SR_em_obsdata(Obj_Space)
    ! !USES: 
    use ESMF
    use LIS_mpiMod
    use LIS_coreMod,  only : LIS_rc, LIS_domain, LIS_npes, LIS_localPet
    use LIS_timeMgrMod, only : LIS_calendar
    use LIS_logMod,     only : LIS_logunit, LIS_verify, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
    use LIS_fileIOMod,      only : LIS_readData
    use AMSRE_SR_em_obsMod, only : SRemobs
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

    real,    pointer    :: obse(:)
    integer, parameter  :: numchannels=6
    integer, parameter  :: numpolarizations=2
    type(ESMF_Field)    :: emField
    character(len=LIS_CONST_PATH_LEN) :: emobsdir 
    logical             :: data_update
    integer             :: status 
    logical             :: found
    type(ESMF_TimeInterval) :: delta_t
    type(ESMF_Time)         :: lis_time1
    integer             :: t
    integer             :: n 
!    logical             :: is_ascend_pass_hr
    logical             :: ob_in_curr_hr
    integer             :: day_index
    logical             :: is_overpass_hr
!    integer, parameter :: freq=1 ! start with 7
!    integer, parameter :: polar=2 ! start with 7H
    integer            :: ipass
    integer                   ::  r,c,gid !table indices
    real                      :: ob
    character*2, parameter :: fnames(6)=(/'07', '11', '19', '24', '37', '89'/)
    character*1, parameter :: pnames(2)=(/'V','H'/)
    integer                   ::  f !freq counter
    integer                   ::  p !polarization counter

    n = 1

    !initialize some values
    is_overpass_hr=.false.
    found=.false.
    call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
         .false., rc=status)
    call LIS_verify(status)


    if ((LIS_rc%hr .eq. SRemobs(n)%overpass_hr_a) &
         .or.(LIS_rc%hr .eq. SRemobs(n)%overpass_hr_d)) then
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

       delta_t = lis_time1 - SRemobs(n)%start_time
       call ESMF_TimeIntervalGet(delta_t, d=day_index, &
            calendar=LIS_calendar,rc=status)
       call LIS_verify(status, 'error in timeget: read_AMSRE_SRdata.F90')

       !as fortran arrays are zero-based...
       day_index=day_index+1

       if (LIS_rc%hr .eq. SRemobs(n)%overpass_hr_a) then
          ipass=SRemobs(n)%index_ascend
       elseif (LIS_rc%hr .eq. SRemobs(n)%overpass_hr_d) then
          ipass=SRemobs(n)%index_descend
       endif

       do f=1,SRemobs(n)%numfreqs
          do p=1,SRemobs(n)%numpolarizations
             
             call ESMF_StateGet(Obj_Space,"Emissivity" // fnames(f) // pnames(p),emField,&
                  rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(emField,localDE=0,farrayPtr=obse,rc=status)
             call LIS_verify(status)
             
             obse=LIS_rc%udef 
             
             do r=1,LIS_rc%lnr(n)
                do c=1,LIS_rc%lnc(n)
                   gid = LIS_domain(n)%gindex(c,r)
                   ob=SRemobs(n)%emissivity(c,r,day_index,f,p,ipass)
                   if (ob.ne.LIS_rc%udef) then
                      obse(gid)=ob
                      found=.true.
                   end if
                enddo
             enddo
          enddo
       enddo
       !          deallocate(obse)
       if(found) then      
          write(LIS_logunit,*) 'Read emissivity data'
          call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
               .true., rc=status)
          call LIS_verify(status)
       endif
    endif

end subroutine read_AMSRE_SR_em_obsdata


