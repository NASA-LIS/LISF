!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_coop
! \label{readcrd_coop}
!
! !REVISION HISTORY:
! 13Jul2011; Yuqiong Liu, Initial Specification
!
! !INTERFACE:    
subroutine readcrd_coop()
! !USES:
  use ESMF
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc, LIS_config
  use coop_forcingMod, only : coop_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to COOP station data 
!  forcing from the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"COOP forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,coop_struc(n)%coopdir,&
          rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"COOP metadata file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,coop_struc(n)%metadata,&
          rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"COOP coord file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_Config, coop_struc(n)%coorddata,&
          rc=rc)

     write(LIS_logunit,*)'Using COOP forcing'
     write(LIS_logunit,*) 'COOP forcing directory :',trim(coop_struc(n)%coopdir)
     write(LIS_logunit,*) 'COOP metadata file : ',trim(coop_struc(n)%metadata)
     write(LIS_logunit,*) 'COOP coord file : ',trim(coop_struc(n)%coorddata)

     coop_struc(n)%cooptime2  = 3000.0
     coop_struc(n)%cooptime1  = 0.0
     coop_struc(n)%startRead = .false.
  enddo

end subroutine readcrd_coop
