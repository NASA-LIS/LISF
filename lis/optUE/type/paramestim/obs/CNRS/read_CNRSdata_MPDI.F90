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
! !ROUTINE: read_CNRS_mpdi_em_obsdata
! \label{read_CNRS_mpdi_em_obsdata}
!
! !REVISION HISTORY:
!  11 Jul 2011: Ken Harrison; Initial Specification
!
! !INTERFACE: 
  subroutine read_CNRS_mpdi_em_obsdata(Obj_Space)
! !USES: 
    use ESMF
    use LIS_RTMMod, only : LIS_sfcState, LIS_forwardState
    use LIS_mpiMod
    use LIS_coreMod,  only : LIS_rc, LIS_domain, LIS_npes, LIS_localPet
    use LIS_timeMgrMod, only : LIS_calendar
    use LIS_logMod,     only : LIS_logunit, LIS_verify, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_fileIOMod,      only : LIS_readData
    !use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
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
    logical             :: data_update
    integer             :: status 
    logical             :: found
    logical             :: maskout
    integer             :: c,r
    type(ESMF_TimeInterval) :: delta_t
    type(ESMF_Time)         :: lis_time1
    integer             :: t
    integer             :: n 
    logical             :: ob_in_curr_hr
    integer             :: day_index
    logical             :: is_overpass_hr
    integer, parameter :: channel_focus=1 ! start with 19V
    real :: v_em(3) ! 3 platforms
    real :: h_em(3) 
    real :: mpdi(3)
    integer :: nobs(3)
    integer :: ipass
    integer :: i
    integer :: nobs_sum, platform_sum
    real    :: mpdi_sum
    real    :: mpdi_ave
    type(ESMF_Field)      :: varField

    real,    pointer    :: land_temperature(:), soil_temperature(:), &
         canopy_water_content(:), snow_depth(:)

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

    maskout=.false.
    if (is_overpass_hr) then
       !Screen out if snow or frozen ground
       call ESMF_StateGet(LIS_sfcState(n), "Snow Depth", varField, rc=status)
       call LIS_verify(status, 'Error in StateGet: read_CNRSdata_MPDI.F90')
       call ESMF_FieldGet(varField, localDE=0,farrayPtr=snow_depth, rc=status)
       call LIS_verify(status, 'Error in FieldGet: read_CNRSdata_MPDI.F90')
       
       call ESMF_StateGet(LIS_sfcState(n), "Land Temperature", varField, rc=status)
       call LIS_verify(status, 'Error in StateGet: read_CNRSdata_MPDI.F90')
       call ESMF_FieldGet(varField, localDE=0,farrayPtr=land_temperature, rc=status)
       call LIS_verify(status, 'Error in FieldGet: read_CNRSdata_MPDI.F90')
       
       call ESMF_StateGet(LIS_sfcState(n), "Soil Temperature", varField, rc=status)
       call LIS_verify(status, 'Error in StateGet: read_CNRSdata_MPDI.F90')
       call ESMF_FieldGet(varField, localDE=0,farrayPtr=soil_temperature, rc=status)
       call LIS_verify(status, 'Error in FieldGet: read_CNRSdata_MPDI.F90')
       
       do t=1,LIS_rc%ntiles(n)  !for CNRS data, if any one is true then  maskout
          if(snow_depth(t).gt.0.00000001)  then
             maskout=.true.
          endif
          if(land_temperature(t).lt.275.0) then
             maskout=.true.
          endif
          if(soil_temperature(t).lt.275.0) then
             maskout=.true.
          endif
       enddo
    endif

    if (is_overpass_hr.and.(.not.maskout)) then

       call ESMF_TimeSet(lis_time1, yy=LIS_rc%yr, mm=LIS_rc%mo, &
            dd=LIS_rc%da,h=LIS_rc%hr,&
            m=LIS_rc%mn,calendar=LIS_calendar,rc=status)

       delta_t = lis_time1 - CNRS_em_obs_struc(n)%start_time
       call ESMF_TimeIntervalGet(delta_t, d=day_index, &
            calendar=LIS_calendar,rc=status)
       call LIS_verify(status, 'error in timeget: read_CNRSdata_MPDI.F90')

       !as fortran arrays are one-based...
       day_index=day_index+1

       if (day_index.lt.1) then
          found=.false.
       else
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
          
          do i=1,3
          v_em(i)= &
               CNRS_em_obs_struc(n)%emissivity &
               ( &
               day_index, &
               channel_focus, &
               ipass, &
               i & 
               ) 
          h_em(i)= &
               CNRS_em_obs_struc(n)%emissivity &
               ( &
               day_index, &
               channel_focus+1, &   ! note channel_focus + 1
               ipass, &
               i & 
               ) 
          if(and(v_em(i).ne.LIS_rc%udef,h_em(i).ne.LIS_rc%udef)) then
             mpdi(i)= (v_em(i)-h_em(i))/(v_em(i)+h_em(i))
          else
             mpdi(i)=LIS_rc%udef
          endif
          nobs(i)=CNRS_em_obs_struc(n)%nobs &
               ( &
               day_index, &
               ipass, &
               i &
               )
       enddo
       
       !First pass to count obs and platforms w/ min 3 obs
       nobs_sum=0
       platform_sum=0
       do i=1,3 !platform F13, 14, 15
          if (mpdi(i).ne.LIS_rc%udef) then
             nobs_sum=nobs_sum+nobs(i)
             if(nobs(i).ge.3) then
                platform_sum=platform_sum+1
             endif
          endif
       enddo
       
       !Second pass to compute average MPDI if conditions met
       mpdi_sum=0.0
       mpdi_ave=LIS_rc%udef
       if( (nobs_sum.ge.15) .and. (platform_sum.ge.2) ) then
          do i=1,3  !platform F13, 14, 15
             if (mpdi(i).ne.LIS_rc%udef) then
                mpdi_sum=mpdi_sum+real(nobs(i))*mpdi(i)
             endif
          enddo
          mpdi_ave=mpdi_sum/nobs_sum
          found=.true.
       endif
       
       do t=1,LIS_rc%ngrid(n)
          obse(t)=mpdi_ave
       enddo
       

    end if
    if(found) then      
       write(LIS_logunit,*) 'Read emissivity data'
       call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
            .true., rc=status)
       call LIS_verify(status)
    endif
 endif
 
end subroutine read_CNRS_mpdi_em_obsdata

!!!subroutine mygetsfcvar(sfcState, varname, var)
!!!  ! !USES: 
!!!  use ESMF
!!!  use LIS_logMod,  only : LIS_verify
!!!  
!!!  implicit none
!!!  
!!!  type(ESMF_State)      :: sfcState
!!!  type(ESMF_Field)      :: varField
!!!  character(len=*)      :: varname
!!!  real, pointer         :: var(:)
!!!  integer               :: status
!!!
!!!  call ESMF_StateGet(sfcState, trim(varname), varField, rc=status)
!!!  call LIS_verify(status, 'Error in StateGet: read_CNRSdata_MPDI.F90 '//trim(varname))
!!!  call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
!!!  call LIS_verify(status, 'Error in StateGet: read_CNRSdata_MPDI.F90 '//trim(varname))
!!!  
!!!end subroutine mygetsfcvar

