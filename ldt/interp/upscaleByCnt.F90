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
! !ROUTINE: upscaleByCnt
! \label{upscaleByCnt}
!
! !INTERFACE:    
subroutine upscaleByCnt( mi, mo, nt, udef, n11,li, gi, lo, go )

  implicit none
! 
! !USES:   
!
! !INPUT PARAMETERS: 
  integer,  intent(in) :: mi
  integer,  intent(in) :: mo
  integer,  intent(in) :: nt
  real,     intent(in) :: udef
  integer,  intent(in) :: n11(mi)
  logical*1,intent(in) :: li(mi)
  real,     intent(in) :: gi(mi)

! !OUTPUT PARAMETERS:
  logical*1,intent(out):: lo(mo,nt)
  real,     intent(out):: go(mo,nt)
!
! !DESCRIPTION: 
!  This subprogram upscales scalar data from a finer grid to a coarser
!   grid by counting up the parameter's types or bins. The grids are defined 
!   by their grid description arrays, which are based on the decoding 
!   schemes used by NCEP. 
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
!    \item[nt] 
!     total number of types or bins 
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
! !REVISION HISTORY: 
!  21 Jul 2012; KR Arsenault:  Adapted NCEP code for upscaling
! 
!EOP
!
  integer  :: i, t
  integer  :: itype

  go = 0.0

!- Count types/bins:
   do i = 1, mi
      if(li(i)) then 
         if( n11(i) .ne. 0 .and. gi(i) .ne. udef ) then            
           itype = nint(gi(i))
           if( itype .ne. 0 ) then
              go(n11(i),itype) = go(n11(i),itype) + 1.0
           endif
         endif
      endif
   enddo

   do i = 1, mo
      do t = 1, nt
         if(go(i,t).gt.0) then 
            lo(i,t) = .true.
         else
!           go(i,t) = udef
           go(i,t) = 0.
           lo(i,t) = .false. 
         endif
      enddo
   enddo

end subroutine upscaleByCnt

