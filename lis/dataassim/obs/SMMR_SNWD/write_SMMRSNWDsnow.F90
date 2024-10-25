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
! !ROUTINE: write_SMMRSNWDsnowobs
! \label{write_SMMRSNWDsnowobs}
! 
! !REVISION HISTORY: 
! 16 Oct 2012: Sujay Kumar; Initial Specification
! 
! !INTERFACE: 
subroutine write_SMMRSNWDsnowobs(n, k, OBS_State)
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
! SMMR observations to a file
! 
!EOP
  type(ESMF_Field)         :: snowField
  logical                  :: data_update
  real, pointer            :: snowobs(:)
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer                  :: ftn
  integer                  :: status
  
!  integer                  :: c,r
!  real                     :: testdata(LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k))

  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01",snowField, &
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(snowField, localDE=0, farrayPtr=snowobs, rc=status)
     call LIS_verify(status)

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call SMMRSNWDsnow_obsname(n,k,obsname)        

        call LIS_create_output_directory('DAOBS')
        open(ftn,file=trim(obsname), form='unformatted')
     endif

#if 0 
     testdata = LIS_rc%udef
     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              testdata(c,r) = snowobs(LIS_obs_domain(n,k)%gindex(c,r))
           endif
        enddo
     enddo
     if(LIS_localPet.eq.241) then 
        open(100,file='test_out.bin',form='unformatted')
        write(100) testdata
        close(100)
        stop
     endif
#endif
     call LIS_writevar_gridded_obs(ftn,n,k,snowobs)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  endif  

end subroutine write_SMMRSNWDsnowobs

!BOP
! !ROUTINE: SMMRSNWDsnow_obsname
! \label{SMMRSNWDsnow_obsname}
! 
! !INTERFACE: 
subroutine SMMRSNWDsnow_obsname(n,k,obsname)
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
  
end subroutine SMMRSNWDsnow_obsname
