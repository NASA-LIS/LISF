!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
      module clsmf25_sibalb_coeff

      implicit none

      PRIVATE

      public :: clsmf25_coeffsib

      contains

      FUNCTION clsmf25_COEFFSIB(TABLE, NTABL, LAI ,DX, DY)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NTABL, LAI

      REAL, INTENT(IN) :: DX, DY
      REAL, INTENT(IN), DIMENSION(NTABL,2) :: TABLE
      REAL clsmf25_COEFFSIB

      clsmf25_COEFFSIB = (TABLE(LAI,  1)                                         &
             + (TABLE(LAI  ,2) - TABLE(LAI  ,1)) * DY ) * (1.0-DX)             &
             + (TABLE(LAI+1,1)                                                 &
             + (TABLE(LAI+1,2) - TABLE(LAI+1,1)) * DY ) * DX

      RETURN
      END FUNCTION clsmf25_COEFFSIB

    end module clsmf25_sibalb_coeff
