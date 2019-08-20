!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_gdasLSWG
! \label{readcrd_gdasLSWG}
!
!
! !REVISION HISTORY:
! 20 Oct 2009; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_gdasLSWG()
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_config
  use LIS_logMod,    only : LIS_logunit, LIS_verify
  use gdasLSWG_forcingMod, only: gdasLSWG_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GDASLSWG forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"GDASLSWG forcing file:",rc=rc)
  call LIS_verify(rc,'GDASLSWG forcing file: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdasLSWG_struc(n)%gdasLSWGfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GDASLSWG domain lower left lat:",rc=rc)
  call LIS_verify(rc,'GDASLSWG domain lower left lat: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdasLSWG_struc(n)%gridDesci(4),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GDASLSWG domain lower left lon:",rc=rc)
  call LIS_verify(rc,'GDASLSWG domain lower left lon: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdasLSWG_struc(n)%gridDesci(5),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GDASLSWG domain upper right lat:",rc=rc)
  call LIS_verify(rc,'GDASLSWG domain upper right lat: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdasLSWG_struc(n)%gridDesci(7),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GDASLSWG domain upper right lon:",rc=rc)
  call LIS_verify(rc,'GDASLSWG domain upper right lon: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdasLSWG_struc(n)%gridDesci(8),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GDASLSWG domain resolution (dx):",rc=rc)
  call LIS_verify(rc,'GDASLSWG domain resolution (dx): not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdasLSWG_struc(n)%gridDesci(9),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GDASLSWG domain resolution (dy):",rc=rc)
  call LIS_verify(rc,'GDASLSWG domain resolution (dy): not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdasLSWG_struc(n)%gridDesci(10),rc=rc)
  enddo

  do n=1,LIS_rc%nnest

     write(LIS_logunit,*) 'Using GDASLSWG forcing'
     write(LIS_logunit,*) 'GDASLSWG forcing file: ',gdasLSWG_struc(n)%GDASLSWGfile
     gdasLSWG_struc(n)%TIME1  = 3000.0
     gdasLSWG_struc(n)%TIME2  = 0.0
     gdasLSWG_struc(n)%startRead = .true. 
  enddo

end subroutine readcrd_gdasLSWG
