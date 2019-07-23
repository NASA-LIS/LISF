!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_merraland
! \label{readcrd_merraland}
!
! !REVISION HISTORY:
! 12 Oct 2009; Eric Kemp, Initial code
! 22 Jul 2010: David Mocko, changed to hourly forcing
!
! !INTERFACE:    
subroutine readcrd_merraland()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod
  use merraland_forcingMod, only : merraland_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to MERRA-Land forcing
!  from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"MERRA-Land forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,merraland_struc(n)%merralanddir,&
          rc=rc)
     call LIS_verify(rc,&
          'MERRA-Land forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"MERRA-Land use lowest model level forcing:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,merraland_struc(n)%uselml,rc=rc)
     call LIS_verify(rc,&
          'MERRA-Land use lowest model level forcing: not defined')
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) 'Using MERRA-Land forcing'
     write(LIS_logunit,*) 'MERRA-Land forcing directory: ',&
          merraland_struc(n)%MERRALANDDIR

     merraland_struc(n)%ncold = 540
     merraland_struc(n)%nrold = 361
     merraland_struc(n)%merralandtime1 = 3000.0
     merraland_struc(n)%merralandtime2 = 0.0

  enddo
end subroutine readcrd_merraland
