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
! !ROUTINE: readcrd_princeton
! \label{readcrd_princeton}
!
! !REVISION HISTORY:
! 26 Jan 2007; Hiroko Kato, Initial Code adopted from readcrd_princeton.F90
! 25 Jun 2007; Hiroko Kato, upgraded to LDTv5.0
! 25 May 2014; KR Arsenault/Hiroko Kato-Beaudoing, added Princeton elev field
!
! !INTERFACE:    
subroutine readcrd_princeton(findex)

! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod,  only : LDT_logunit
  use princeton_forcingMod, only : princeton_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to PRINCETON forcing from 
!  the LDT configuration file. 
!  
!EOP
  implicit none
  integer, intent(in) :: findex

  integer :: n, rc

  if( LDT_rc%met_ecor_parms(findex) .ne. "none" .and. &
      LDT_rc%runmode == "LSM parameter processing" ) then

    call ESMF_ConfigFindLabel(LDT_config,"PRINCETON elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,princeton_struc(n)%elevfile,rc=rc)
    enddo
    return
  endif

  call ESMF_ConfigFindLabel(LDT_config,"PRINCETON forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,princeton_struc(n)%princetondir,rc=rc)
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) 'PRINCETON forcing directory : ',&
           trim(princeton_struc(n)%princetonDIR)
     princeton_struc(n)%princetontime1 = 3000.0
     princeton_struc(n)%princetontime2 = 0.0
  enddo

end subroutine readcrd_princeton
