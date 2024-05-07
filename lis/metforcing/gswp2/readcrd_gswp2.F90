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
! !ROUTINE: readcrd_gswp2
! \label{readcrd_gswp2}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_gswp2()
! !USES:
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc, LIS_config
  use gswp2_forcingMod, only : gswp2_struc
  use ESMF
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GSWP2 forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc
  call ESMF_ConfigFindLabel(LIS_config,"GSWP2 landmask file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gswp2_struc(n)%mfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GSWP2 2m air temperature map:",rc=rc)
  do n=1,LIS_rc%nnest  
     call ESMF_ConfigGetAttribute(LIS_config,gswp2_struc(n)%tair,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"GSWP2 2m specific humidity map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gswp2_struc(n)%qair,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GSWP2 wind map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gswp2_struc(n)%wind,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GSWP2 surface pressure map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gswp2_struc(n)%psurf,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"GSWP2 convective rainfall rate map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gswp2_struc(n)%rainf_c,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"GSWP2 rainfall rate map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gswp2_struc(n)%rainf,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"GSWP2 snowfall rate map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gswp2_struc(n)%snowf,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"GSWP2 incident shortwave radiation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gswp2_struc(n)%swdown,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"GSWP2 incident longwave radiation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gswp2_struc(n)%lwdown,rc=rc)
  enddo

  write(LIS_logunit,*)'Using GSWP2 forcing'

end subroutine readcrd_gswp2
