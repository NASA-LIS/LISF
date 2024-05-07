!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module vic411_module
!BOP
!
! !MODULE: vic411_module.F90
!
! !DESCRIPTION:
!  The code in this file provides a description of the 
!  data structure containing the VIC 4.1.1 1-d variables.
!  The variables specified in the data structure include: 
!
!  \begin{description}
!   \item[min\_Tfactor]
!      Minimum change in temperature due to elevation (C) in each snow
!      elevation band.
!      This value is determined by finding the minimum value of Tfactor from 
!      VIC's soil\_con data structure.
!   \end{description}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 
!EOP
  implicit none

  type vicdec

     real*8 :: min_Tfactor

  end type vicdec

end module vic411_module
