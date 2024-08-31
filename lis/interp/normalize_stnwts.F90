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
! !ROUTINE: normalize_stnwts
! \label{normalize_stnwts}
!
! !REVISION HISTORY: 
!   08Dec04 Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine normalize_stnwts(stndata,nstns,npts,undef,wt)
  implicit none
! !ARGUMENTS: 
  integer    :: nstns
  integer    :: npts
  real       :: stndata(nstns)
  real       :: wt(npts,nstns)
  real       :: undef
!
! !DESCRIPTION: 
!  This subroutine normalizes the interpolation weights to convert 
!  station data into a gridded set. 
!
!  The arguments are: 
!  \begin{description}
!    \item[nstns]
!     number of stations used in interpolation
!    \item[npts]
!     integer maximum number of coordinates
!    \item[stndata]
!     input data (station observations)
!    \item[undef]
!     undefined value used in observations
!    \item[wt] 
!     interpolation weights
!   \end{description}
!
!EOP
  integer :: i,j
  real :: norm(npts)
  norm = 0 
  do i=1,npts
     do j=1,nstns
        if(stndata(j).ne.undef) then
           norm(i) = norm(i)+wt(i,j)
        endif
     enddo
  enddo
  do i=1,npts
     do j=1,nstns
        if(norm(i).eq.0) norm(i) = 1
        wt(i,j) = wt(i,j)/norm(i)
     enddo
  enddo 
end subroutine normalize_stnwts
