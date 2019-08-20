!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_pildas
! \label{readcrd_pildas}
!
! !REVISION HISTORY:
! 25 Apr 2013: Sujay Kumar, initial specification
! 14 Jul 2016: Mahdi Navari - Modified for PILDAS
!
! !INTERFACE:    
subroutine readcrd_pildas()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod
  use pildas_forcingMod, only : pildas_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to Merra-Land forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"PILDAS forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,pildas_struc(n)%pildasdir,&
          rc=rc)
     call LIS_verify(rc,&
          'PILDAS forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"PILDAS forcing version:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,pildas_struc(n)%version,&
          rc=rc)
     call LIS_verify(rc,&
          'PILDAS forcing version: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"PILDAS forcing use lowest model level fields:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,pildas_struc(n)%uselml,&
          rc=rc)
     call LIS_verify(rc,&
          'PILDAS forcing use lowest model level fields: not defined')
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) 'Using PILDAS forcing'
     write(LIS_logunit,*) 'PILDAS forcing directory: ',&
          pildas_struc(n)%PILDASDIR

     pildas_struc(n)%nc = 70
     pildas_struc(n)%nr = 50
     pildas_struc(n)%pildastime1 = 3000.0
     pildas_struc(n)%pildastime2 = 0.0

  enddo
end subroutine readcrd_pildas
