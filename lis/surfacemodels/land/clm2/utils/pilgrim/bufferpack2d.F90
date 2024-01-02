!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !IROUTINE: bufferpack2d --- Pack a ghost region into a buffer
!
! !INTERFACE:  
      subroutine bufferpack2d ( In, xfrom, xto, yfrom, yto,             &
                                x1, x2, y1, y2, buffer )

! !USES:
      use LIS_precisionMod
      implicit none

! !INPUT PARAMETERS:
      integer , intent( in )     :: xfrom,xto                 ! X dims
      integer , intent( in )     :: yfrom,yto                 ! Y dim
      real(r8), intent( in )     :: in(xfrom:xto,yfrom:yto)   ! In array
      integer , intent( in )     :: x1,x2                     ! X limits
      integer , intent( in )     :: y1,y2                     ! Y limits

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent( inout )  :: buffer(x1:x2,y1:y2)       ! Packed buffer

! !DESCRIPTION:
!
!     This routine puts a 2-D ghost region at the end of a buffer, first in
!     X then in Y.
!
! !LOCAL VARIABLES:
      integer  I, J
!
! !REVISION HISTORY:
!   99.09.30   Sawyer     Creation
!   99.10.18   Sawyer     FVCCM3 format (no capitalization)
!   00.07.08   Sawyer     Use precision module
!   01.02.02   Sawyer     Now free format
!
!EOP
!-----------------------------------------------------------------------
!BOC

!
! Fill the buffer, first in X then in Y
!
      do J = y1, y2
        do I = x1, x2
          buffer( I,J ) = in( I, J )
        enddo
      enddo
      return
!EOC
      end subroutine bufferpack2d
!-----------------------------------------------------------------------
