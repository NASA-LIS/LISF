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
! !IROUTINE: bufferunpack4d --- Pack a ghost region into a buffer
!
! !INTERFACE:  
      subroutine bufferunpack4d( buffer, x1, x2, y1, y2,                &
                                 z1, z2, n1, n2, xfrom,                 &
                                 xto, yfrom, yto, zfrom, zto,           &
                                 Nfrom, Nto, Out )

! !USES:
      use LIS_precisionMod
      implicit none

! !INPUT PARAMETERS:
      integer , intent( in )  :: x1,x2                     ! X limits
      integer , intent( in )  :: y1,y2                     ! Y limits
      integer , intent( in )  :: z1,z2                     ! Z limits
      integer , intent( in )  :: n1,n2                     ! N limits
      real(r8), intent( in )  :: buffer(x1:x2,y1:y2,z1:z2,n1:n2)
      integer , intent( in )  :: xfrom,xto                 ! X dims
      integer , intent( in )  :: yfrom,yto                 ! Y dim
      integer , intent( in )  :: zfrom,zto                 ! Z dim
      integer , intent( in )  :: Nfrom,Nto                 ! N dim

! !INPUT/OUTPUT PARAMETERS:
      real(r8),intent(inout)  :: out(xfrom:xto,yfrom:yto,               &
                                     zfrom:zto,nfrom:nto)

! !DESCRIPTION:
!
!     This routine puts a 4-D ghost region at the end of a buffer, 
!     packed in the order of the indices.
!
! !LOCAL VARIABLES:
      integer  I, J, K, L
!
! !REVISION HISTORY:
!   99.09.30   Sawyer     Creation
!   99.10.18   Sawyer     FVCCM3 format (no capitalization)
!   00.07.08   Sawyer     Use precision module
!   01.02.02   Sawyer     Removed SGI directives. OpenMP only; free format
!
!EOP
!-----------------------------------------------------------------------
!BOC

!
! Fill the buffer, first in X then in Y
!
      do L = n1, n2
!$omp  parallel do&
!$omp& default(shared)&
!$omp& private(i,j,k)
        do K = z1, z2
          do J = y1, y2
            do I = x1, x2
              out( I,J,K,L ) = buffer( I,J,K,L )
            enddo
          enddo
        enddo
      enddo
      return
!EOC
      end subroutine bufferunpack4d
!-----------------------------------------------------------------------

