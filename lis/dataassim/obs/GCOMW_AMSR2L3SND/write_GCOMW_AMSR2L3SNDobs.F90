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
! !ROUTINE: write_GCOMW_AMSR2L3SNDobs
! \label{write_GCOMW_AMSR2L3SNDobs}
! 
! !REVISION HISTORY: 
! 25Jan2008: Sujay Kumar; Initial Specification
! 
! !INTERFACE: 
subroutine write_GCOMW_AMSR2L3SNDobs(n, k, OBS_State)
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_historyMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none

! !ARGUMENTS: 

  integer,     intent(in)  :: n 
  integer,     intent(in)  :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION: 
! 
! writes the transformed (interpolated/upscaled/reprojected)  
! LPRM AMSRE observations to a file
! 
!EOP
  type(ESMF_Field)         :: snowfield
  logical                  :: data_update
  real, pointer            :: sndobs(:)
  real                     :: imsobs(LIS_rc%obs_ngrid(k))
  real                     :: modisobs(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer                  :: ftn
  integer                  :: status

  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01",snowfield, &
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(snowfield, localDE=0, farrayPtr=sndobs, rc=status)
     call LIS_verify(status)

     if(LIS_rc%obs_ngrid(k).gt.0) then 
        call ESMF_AttributeGet(snowfield,"IMS data",imsobs,&
             itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status, 'Error in AttributeGet-IMS data')
        
        call ESMF_AttributeGet(snowfield,"MODIS data",modisobs,&
             itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status, 'Error in AttributeGet-MODIS data')
     endif

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call GCOMW_AMSR2_sndobsname(n,k,obsname)        

        call LIS_create_output_directory('DAOBS')
        open(ftn,file=trim(obsname), form='unformatted')
     endif

     call LIS_writevar_gridded_obs(ftn,n,k,sndobs)
     call LIS_writevar_gridded_obs(ftn,n,k,imsobs)
     call LIS_writevar_gridded_obs(ftn,n,k,modisobs)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  endif  

end subroutine write_GCOMW_AMSR2L3SNDobs

!BOP
! !ROUTINE: GCOMW_AMSR2_sndobsname
! \label{GCOMW_AMSR2_sndobsname}
! 
! !INTERFACE: 
subroutine GCOMW_AMSR2_sndobsname(n,k,obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  integer               :: n
  integer               :: k
  character(len=*)      :: obsname
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

  obsname = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//&
       '/LISDAOBS_'//cdate1// &
       trim(cda)//trim(cdate)//'.1gs4r'
    
end subroutine GCOMW_AMSR2_sndobsname
