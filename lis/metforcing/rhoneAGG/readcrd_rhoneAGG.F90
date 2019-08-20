!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_rhoneAGG
! \label{readcrd_rhoneAGG}
!
! !REVISION HISTORY:
! 11 Dec 2003: Sujay Kumar, Initial Code
!
! !INTERFACE:
subroutine readcrd_rhoneAGG()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod, only    : LIS_logunit
  use rhoneAGG_forcingMod, only : rhoneAGG_struc
!
! !DESCRIPTION:
!  Routine to read RHONEAGG specific parameters from the 
!  LIS configuration file.
!
!EOP
  implicit none
  integer  :: n, rc

  call ESMF_ConfigFindLabel(LIS_config, "Rhone AGG forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, rhoneAGG_struc(n)%rhoneAGGdir, rc=rc)
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"Rhone AGG domain x-dimension size:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,rhoneAGG_struc(n)%ncold,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Rhone AGG domain y-dimension size:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,rhoneAGG_struc(n)%nrold,rc=rc)
  enddo
   
  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'Using Rhone AGG forcing'
     write(LIS_logunit,*)'Rhone AGG forcing directory: ',rhoneAGG_struc(n)%RHONEAGGDIR
     rhoneAGG_struc(n)%nmif = 9
     rhoneAGG_struc(n)%rhoneAGGtime1 = 3000.0
     rhoneAGG_struc(n)%rhoneAGGtime2 = 0.0
  enddo
end subroutine readcrd_rhoneAGG

