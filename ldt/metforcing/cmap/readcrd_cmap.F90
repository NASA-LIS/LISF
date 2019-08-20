!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_cmap
! \label{readcrd_cmap}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_cmap()
! !USES:
  use ESMF 
  use cmap_forcingMod, only : cmap_struc
  use LDT_coreMod,     only : LDT_config,LDT_rc
  use LDT_logMod,      only : LDT_logunit
!
! !DESCRIPTION:
!
!  This routine reads the options specific to CMAP forcing from 
!  the LDT configuration file. 
!  
!EOP
  implicit none
  
  integer :: n, rc

  write(LDT_logunit,*)" Using CMAP forcing"

  call ESMF_ConfigFindLabel(LDT_config, "CMAP forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,cmap_struc(n)%cmapdir,rc=rc)
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*)" CMAP forcing directory : ",trim(cmap_struc(n)%cmapdir)
!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure 
! data is read in during first time step
!------------------------------------------------------------------------
     cmap_struc(n)%cmaptime = 0.0
  enddo

end subroutine readcrd_cmap
