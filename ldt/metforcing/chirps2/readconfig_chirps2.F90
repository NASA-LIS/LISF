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
! !ROUTINE: readconfig_chirps2
! \label{readconfig_chirps2}
!
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 06 Jul 2015; KR Arsenault, Added support for latest CHIRPS dataset
!
! !INTERFACE:    
subroutine readconfig_chirps2()

! !USES:
  use ESMF 
  use chirps2_forcingMod, only : chirps2_struc
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod,  only : LDT_logunit
!
! !DESCRIPTION:
!
!  This routine reads the options specific to CHIRPS 2.0 forcing 
!   from the LDT configuration file. 
!  
!EOP
  implicit none

  integer :: n,rc

   call ESMF_ConfigFindLabel(LDT_config,"CHIRPS2.0 forcing directory:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,chirps2_struc(n)%directory,rc=rc)
      write(LDT_logunit,*) "CHIRPS 2.0 forcing directory :",&
            trim(chirps2_struc(n)%directory)
   enddo
 
   call ESMF_ConfigFindLabel(LDT_config,"CHIRPS2.0 forcing resolution:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,chirps2_struc(n)%xres,rc=rc)
      chirps2_struc(n)%yres = chirps2_struc(n)%xres
      write(LDT_logunit,*) "CHIRPS 2.0 forcing resolution :",&
            chirps2_struc(n)%xres
   enddo

   chirps2_struc(:)%chirpstime1 = 0.
   chirps2_struc(:)%chirpstime2 = 0.

end subroutine readconfig_chirps2
