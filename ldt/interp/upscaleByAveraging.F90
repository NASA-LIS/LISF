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
! !ROUTINE: upscaleByAveraging
! \label{upscaleByAveraging}
!
! !INTERFACE:    
subroutine upscaleByAveraging(mi, mo, udef, n11, li, gi, lo, go )
!$ use omp_lib
  implicit none
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subprogram upscales scalar data from a finer grid to a coarser
!  grid using linear averaging. 
!  The grids are defined by their grid description arrays. 
!  
!  The grid description arrays are based on the decoding 
!  schemes used by NCEP. 
!   
!  The current code recognizes the following projections: \newline
!             (gridDesc(1)=0) equidistant cylindrical \newline
!            
!  where gridDesc could be defined for either the input grid or the 
!  output grid. 
!
!  The arguments are: 
!  \begin{description}
!    \item[mi] 
!     total number of points in the input grid
!    \item[mo] 
!     total number of points in the output grid
!    \item[udef]
!     undefine value 
!    \item[n11] 
!     array that maps the location of each input grid
!     point in the output grid. 
!    \item[li]
!     input bit mask
!    \item[gi]
!     input data (finer grid)
!    \item[lo]
!     output bit mask
!    \item[go]
!     output data (coarser grid)
!    \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!
! !ARGUMENTS:
  integer       :: mi
  integer       :: mo
  real          :: udef
  integer       :: n11(mi)
  logical*1     :: li(mi)
  real          :: gi(mi)
  logical*1     :: lo(mo)
  real          :: go(mo)
!
!EOP
  integer       :: i 
  integer       :: ngo(mo)

  go = 0.0
  ngo = 0 
!$OMP PARALLEL
!$OMP DO
  do i=1,mi
     if(li(i)) then 
        if(n11(i).ne.0) then            
           go(n11(i)) = go(n11(i)) + gi(i)
           ngo(n11(i)) = ngo(n11(i)) + 1       
        endif
     endif
  enddo
!$OMP END DO
!$OMP DO
  do i=1,mo
     if(ngo(i).gt.0) then 
        go(i) = go(i)/ngo(i)
        lo(i) = .true.
     else
        go(i) = udef
        lo(i) = .false. 
     endif
  enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine upscaleByAveraging

