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
! !ROUTINE: hyssib_mapvegc.F90
! \label{hyssib_mapvegc}
!
! !REVISION HISTORY:
! 28 Apr 2002: K Arsenault, Added SSIB LSM to LDAS
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
! 25 Nov 2007: Chuck Alonge, Updates for LIS 5.0
!
! !INTERFACE:
SUBROUTINE HYSSIB_MAPVEGC(VEGT)

   implicit none
! !ARGUMENTS:
   INTEGER, INTENT(INOUT) :: VEGT
!
! !DESCRIPTION: 
!  This subroutine converts the UMD classes to the SIB classes 
!  used by SSIB LSM (v 2.5).
!  (Originally from Dag Lohmann at NCEP)
!
!  The arguments are: 
!  \begin{description}
!  \item[vegt]
!   UMD Vegetation class number to be converted to SIB class
!  \end{description}
!
!EOP
   INTEGER :: SIBVEG

!  Convert UMD Classes to SIB Classes.
   IF (VEGT.EQ.1)  SIBVEG = 4
   IF (VEGT.EQ.2)  SIBVEG = 1
   IF (VEGT.EQ.3)  SIBVEG = 5
   IF (VEGT.EQ.4)  SIBVEG = 2
   IF (VEGT.EQ.5)  SIBVEG = 3
   IF (VEGT.EQ.6)  SIBVEG = 3
   IF (VEGT.EQ.7)  SIBVEG = 6
   IF (VEGT.EQ.8)  SIBVEG = 8
   IF (VEGT.EQ.9)  SIBVEG = 9
   IF (VEGT.EQ.10) SIBVEG = 7
   IF (VEGT.EQ.11) SIBVEG = 12
   IF (VEGT.EQ.12) SIBVEG = 11
   IF (VEGT.EQ.13) SIBVEG = 11
   IF (VEGT.GT.13) THEN
      SIBVEG = 7
   ENDIF

   VEGT = SIBVEG

   return

end
