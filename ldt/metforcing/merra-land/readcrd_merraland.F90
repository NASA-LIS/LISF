!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_merraland
! \label{readcrd_merraland}
!
! !REVISION HISTORY:
! 12 Oct 2009; Eric Kemp, Initial code
! 22 Jul 2010: David Mocko, changed to hourly forcing
! 22 Jan 2015: KR Arsenault, added to LDT
!
! !INTERFACE:    
subroutine readcrd_merraland()
! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod
  use merraland_forcingMod, only : merraland_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to MERRA-Land forcing
!  from the LDT configuration file. 
!  
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LDT_config,"MERRA-Land forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,merraland_struc(n)%merralanddir,&
          rc=rc)
     call LDT_verify(rc,&
          'MERRA-Land forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"MERRA-Land use lowest model level forcing:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,merraland_struc(n)%uselml,rc=rc)
     call LDT_verify(rc,&
          'MERRA-Land use lowest model level forcing: not defined')
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) 'Using MERRA-Land forcing'
     write(LDT_logunit,*) 'MERRA-Land forcing directory: ',&
          merraland_struc(n)%MERRALANDDIR

     merraland_struc(n)%nc = 540
     merraland_struc(n)%nr = 361
     merraland_struc(n)%merralandtime1 = 3000.0
     merraland_struc(n)%merralandtime2 = 0.0

  enddo
end subroutine readcrd_merraland
