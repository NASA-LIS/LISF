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
! !ROUTINE: geowrsi2_writerst
! \label{geowrsi2_writerst}
! 
! !INTERFACE:
subroutine geowrsi2_writerst(n)

! !USES:
  use LIS_coreMod,     only : LIS_rc
  use LIS_logMod,      only : LIS_logunit
  use geowrsi2_lsmMod, only : geowrsi2_CalcSOSlsmRunMode, geowrsi2_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
!
!  Forcing-only (geowrsi2) option for calling the write restart routines. 
!
!  This routine also writes the SOS/SOSa files from the CalcSOS run mode.
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

!- Write out Start-of-season (SOS) files during SOS-calculation mode:
   if( geowrsi2_CalcSOSlsmRunMode .and. geowrsi2_struc(n)%eos_alarm ) then
      write(LIS_logunit,*) "... SOSCALC_MODE: WRITING OUT SOS/SOSa BIL FILES "
      call geowrsi2_writeSOS_Bilfile( n )
      geowrsi2_struc(n)%eos_alarm = .false.
   endif


!- Write out WRSI-model restart files during WRSI-calculation mode:


end subroutine geowrsi2_writerst
