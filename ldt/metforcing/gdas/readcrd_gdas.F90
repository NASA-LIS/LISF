!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v2.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_gdas
! \label{readcrd_gdas}
!
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_gdas(findex)

! !USES:
  use ESMF
  use LDT_coreMod,   only : LDT_rc, LDT_config
  use LDT_logMod,    only : LDT_logunit, LDT_verify
  use gdas_forcingMod, only: gdas_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GDAS forcing from 
!  the LDT configuration file. 
!  
!EOP
  implicit none
  integer, intent(in) :: findex
  integer :: n, rc

!  if( LDT_rc%met_ecor(findex) .ne. "none" .and. &
  if( LDT_rc%met_ecor_parms(findex) .ne. "none" .and. &
      LDT_rc%runmode == "LSM parameter processing" ) then

    call ESMF_ConfigFindLabel(LDT_config,"GDAS T126 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%gdasT126elevfile(n),rc=rc)
       call LDT_verify(rc,'GDAS T126 elevation map: not specified')
    end do
    call ESMF_ConfigFindLabel(LDT_config,"GDAS T170 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%gdasT170elevfile(n),rc=rc)
       call LDT_verify(rc,'GDAS T170 elevation map: not specified')
    end do
    call ESMF_ConfigFindLabel(LDT_config,"GDAS T254 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%gdasT254elevfile(n),rc=rc)
       call LDT_verify(rc,'GDAS T254 elevation map: not specified')
    end do
    call ESMF_ConfigFindLabel(LDT_config,"GDAS T382 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%gdasT382elevfile(n),rc=rc)
       call LDT_verify(rc,'GDAS T382 elevation map: not specified')
    end do
    call ESMF_ConfigFindLabel(LDT_config,"GDAS T574 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%gdasT574elevfile(n),rc=rc)
      call LDT_verify(rc,'GDAS T574 elevation map: not specified')
    end do
    call ESMF_ConfigFindLabel(LDT_config,"GDAS T1534 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%gdasT1534elevfile(n),rc=rc)
      call LDT_verify(rc,'GDAS T1534 elevation map: not specified')
    end do
     return

  endif  ! End elev correction check

  call ESMF_ConfigFindLabel(LDT_config,"GDAS forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,gdas_struc(n)%gdasdir,rc=rc)
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) 'GDAS forcing directory :',trim(gdas_struc(n)%GDASDIR)
     gdas_struc(n)%GDASTIME1  = 3000.0
     gdas_struc(n)%GDASTIME2  = 0.0
  enddo

end subroutine readcrd_gdas
