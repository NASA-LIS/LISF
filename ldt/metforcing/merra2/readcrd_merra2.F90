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
! !ROUTINE: readcrd_merra2
! \label{readcrd_merra2}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 12 Nov 2015: KR Arsenault, added to LDT
!
! !INTERFACE:    
subroutine readcrd_merra2()
! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod
  use merra2_forcingMod, only : merra2_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to MERRA2 forcing
!  from the LDT configuration file. 
!  
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LDT_config,"MERRA2 forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,merra2_struc(n)%merra2dir,&
          rc=rc)
     call LDT_verify(rc,&
          'MERRA2 forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"MERRA2 use lowest model level forcing:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,merra2_struc(n)%uselml,rc=rc)
     call LDT_verify(rc,&
          'MERRA2 use lowest model level forcing: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"MERRA2 use corrected total precipitation:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,merra2_struc(n)%usecorr,rc=rc)
     call LDT_verify(rc,&
          'MERRA2 use corrected total precipitation: not defined')
  enddo

  ! New! - Added static MERRA-2 terrain height map option:
  call ESMF_ConfigFindLabel(LDT_config,"MERRA2 geopotential terrain height file:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,merra2_struc(n)%merra2hgt_file,rc=rc)
     call LDT_verify(rc,&
          'MERRA2 geopotential terrain height file: not defined')
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) 'Using MERRA2 forcing'
     write(LDT_logunit,*) 'MERRA2 forcing directory: ',&
          merra2_struc(n)%merra2DIR

     merra2_struc(n)%merra2time1 = 3000.0
     merra2_struc(n)%merra2time2 = 0.0

  enddo
end subroutine readcrd_merra2
