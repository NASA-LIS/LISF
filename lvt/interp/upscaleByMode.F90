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
! !ROUTINE: upscaleByMode
! \label{upscaleByMode}
!
! !INTERFACE:    
subroutine upscaleByMode(mi, mo, udef, n11, li, gi, lo, go )

  implicit none
! 
! !USES:   
!
! !DESCRIPTION: 
!  This subprogram upscales scalar data from a finer grid to a coarser
!   grid by finding the dominant type or mode.  The grids are defined 
!   by their grid description arrays, which are based on the decoding 
!   schemes used by NCEP. 
!  Note:  May need some additional testing.
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
!  01 Jun 2012; KR Arsenault:  Adapted NCEP code for mode (dominant type) upscaling
! 
!EOP
!
! !ARGUMENTS:
  integer,  intent(in) :: mi
  integer,  intent(in) :: mo
  real,     intent(in) :: udef
  integer,  intent(in) :: n11(mi)
  logical*1,intent(in) :: li(mi)
  real,     intent(in) :: gi(mi)
  logical*1,intent(out):: lo(mo)
  real,     intent(out):: go(mo)
!
!EOP
  integer  :: i, j, nt, addv
  integer  :: ncts(mo)
  integer  :: max_index(1)
  real     :: min_value
  real, allocatable :: index_values(:,:)
!______________________________________________________

!- Determine min index value to account for 0-indices:
   min_value = minval(gi(:), mask = gi(:).ne.udef )
   if    ( min_value == 0 ) then;  addv = 1
   elseif( min_value == 1 ) then;  addv = 0
   elseif( min_value > 1  ) then;  addv = 0
   endif

!- Determine max number of index values:
   nt = nint(maxval(gi(:)))
   allocate ( index_values(mo,nt+addv) )

   index_values = 0
   ncts = 0;  go = 0.0
!- Count up values:
   do i = 1, mi
      if( li(i) .and. n11(i) .ne. 0 ) then            
          ncts(n11(i)) = ncts(n11(i)) + 1
          index_values(n11(i),nint(gi(i))+addv) = &
               index_values(n11(i),nint(gi(i))+addv) + 1
      endif
   enddo

!- Assign output dominant index value:
   do j = 1, mo
      if(ncts(j).gt.0) then 
         max_index = maxloc(index_values(j,:))
         go(j) = float(max_index(1))-addv 
         lo(j) = .true.
      else
         go(j) = udef
         lo(j) = .false. 
      endif
   enddo
   deallocate( index_values )

end subroutine upscaleByMode

#if 0 
!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: upscaleByDominant
! \label{upscaleByDominant}
!
! !INTERFACE:    
subroutine upscaleByDominant(mi,mo,ncats,&
     udef, n11,li, gi, lo, go)

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
!  grid using dominant lookup. This mode should only be used with 
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
!BOP
!
! !ARGUMENTS:
  integer             :: mi
  integer             :: mo
  integer             :: ncats
  real                :: udef
  integer             :: n11(mi)
  logical*1           :: li(mi)
  real                :: gi(mi)
  logical*1           :: lo(mo)
  real                :: go(mo)

!
!EOP
  integer       :: i 
  integer       :: maxv
  integer       :: ngo(mo,ncats)

  go = udef
  ngo = 0 

  go_t = 0.0

  do i=1,mi
     if(li(i)) then 
        if(n11(i).ne.0) then    
           ngo(n11(i),nint(gi(i))) = ngo(n11(i),nint(gi(i))) + 1
        endif
     endif
  enddo
  
  do i=1,mo
     maxv = -10000
     do k=1,ncats
        if(ngo(i,k).gt.maxv) then 
           maxv = k
        endif
        go(i) = k
     enddo
  enddo
  
  do i=1,mo
     if(go(i).eq.udef) then 
        lo(i) = .false. 
     else
        lo(i) = .true. 
     endif
  enddo

end subroutine upscaleByDominant
#endif
