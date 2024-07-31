!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: write_S1_sigmaVVVHSMLAIobs
! \label{write_S1_sigmaVVVHSMLAIobs}
! 
! !REVISION HISTORY: 
! 30 Aug 2019: Hans Lievens; Initial Specification
! 
! !INTERFACE: 
subroutine write_S1_sigmaVVVHSMLAIobs(n, k, OBS_State)
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_historyMod
  use LIS_DAobservationsMod

  implicit none

! !ARGUMENTS: 

  integer,     intent(in)  :: n 
  integer,     intent(in)  :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION: 
! 
! writes the transformed (interpolated/upscaled/reprojected)  
! S1 observations to a file
! 
!EOP
  type(ESMF_Field)         :: s_vvField, s_vhField
  logical                  :: data_update
  real, pointer            :: obs_vv(:), obs_vh(:)
  character*100            :: obsname_vv, obsname_vh
  integer                  :: ftn
  integer                  :: status
  
  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01",s_vvField, &
          rc=status)
     call LIS_verify(status)

     call ESMF_StateGet(OBS_State, "Observation02",s_vhField, &
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(s_vvField, localDE=0, farrayPtr=obs_vv, rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(s_vhField, localDE=0, farrayPtr=obs_vh, rc=status)
     call LIS_verify(status)

!    ! save VV
     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call S1_sigmaVVVHSMLAI_obsname(n,k,obsname_vv,obsname_vh)        

        call LIS_create_output_directory('DAOBS')

        open(ftn,file=trim(obsname_vv), form='unformatted')
     endif
     call LIS_writevar_gridded_obs(ftn,n,k,obs_vv)
     if(LIS_masterproc) then  
       call LIS_releaseUnitNumber(ftn)
     endif
!    ! save VH
     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call S1_sigmaVVVHSMLAI_obsname(n,k,obsname_vv,obsname_vh)        

        call LIS_create_output_directory('DAOBS')

        open(ftn,file=trim(obsname_vh), form='unformatted')
     endif
     call LIS_writevar_gridded_obs(ftn,n,k,obs_vh)
     if(LIS_masterproc) then  
       call LIS_releaseUnitNumber(ftn)
     endif

  endif  

end subroutine write_S1_sigmaVVVHSMLAIobs

!BOP
! !ROUTINE: S1_sigmaVVVHSMLAI_obsname
! \label{S1_sigmaVVVHSMLAI_obsname}
! 
! !INTERFACE: 
subroutine S1_sigmaVVVHSMLAI_obsname(n,k,obsname_vv,obsname_vh)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  integer               :: n
  integer               :: k
  character(len=*)      :: obsname_vv, obsname_vh
! 
! !DESCRIPTION: 
! 
!EOP

  character(len=12) :: cdate1
  character(len=12) :: cdate
  character(len=10) :: cda

  write(unit=cda, fmt='(a2,i2.2)') '.a',k
  write(unit=cdate, fmt='(a2,i2.2)') '.d',n

  write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
       LIS_rc%yr, LIS_rc%mo, &
       LIS_rc%da, LIS_rc%hr,LIS_rc%mn

  obsname_vv = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//&
       '/LISDAOBS_VV_'//cdate1// &
       trim(cda)//trim(cdate)//'.1gs4r'
  obsname_vh = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//&
       '/LISDAOBS_VH_'//cdate1// &
       trim(cda)//trim(cdate)//'.1gs4r'  
end subroutine S1_sigmaVVVHSMLAI_obsname
