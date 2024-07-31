!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: create_maskfile
!  \label{create_maskfile}
!
! !REVISION HISTORY:
!  10 Aug  2012: KR Arsenault;  Create mask file and generate
!                                type field
!
! !INTERFACE:
 subroutine create_maskfile(n, nt, gridtransform, vegtype, vegcnt, &
                            localmask )

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_endrun

  implicit none

! !ARGUMENTS: 
  integer, intent(in)     :: n
  integer, intent(in)     :: nt
  character*50, intent(in):: gridtransform
  real,    intent(in)     :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real,    intent(in)     :: vegcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),nt)
  real,    intent(out)    :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine creates landmask data and returns the 
!   mask type arrays.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[nt]
!    number of classes or types
!   \item[vegtype]
!    vegetation type
!   \item[vegcnt]
!    fraction of gridcell
!   \item[localmask]
!    landmask for the region of interest
!   \end{description}
!
!EOP      
  integer  :: c, r, t, i
  real     :: vegsum, totalsum
  real     :: maxv
  integer  :: maxt
!_________________________________________________________________________________

   LDT_rc%nmaskpts = 0.
   localmask = 0.

   write(LDT_logunit,fmt=*) '[INFO] LAT/LON -- Creating mask output'

   if( LDT_rc%inc_water_pts ) then
      write(unit=LDT_logunit,fmt=*) '[INFO] INCLUDING WATER POINTS ...'
   endif

   localmask = 0.0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         maxv = 0
         maxt = 0
         do t = 1, nt
            if(vegcnt(c,r,t).gt.maxv) then 
               maxt = t
               maxv = vegcnt(c,r,t)
            endif
         enddo
         if(maxt.gt.0) localmask(c,r) = maxt
      enddo
   enddo

!- NON-tiled option:
   select case( gridtransform )

     case( "none", "mode", "neighbor" ) 

       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             if( vegtype(c,r) == float(LDT_rc%waterclass) ) then
               localmask(c,r) = 0.0
             else
               localmask(c,r) = 1.0
             endif
           enddo
        enddo

!- Tiled option:
     case( "tile" ) 

       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             vegsum = 0.; totalsum = 0.

          !- Assign land/water mask values 
          !  (need to retain water points for mask-parm checks later):
             do t = 1, nt
                if( t /= LDT_rc%waterclass ) then
                  if( vegcnt(c,r,t) /= LDT_rc%udef) then
                     vegsum = vegsum + vegcnt(c,r,t)
                  endif
                endif
             end do


             totalsum = sum( vegcnt(c,r,1:nt), &
                        mask=vegcnt(c,r,1:nt).ne.LDT_rc%udef )

             ! Check gridcell veg total sum:
             if( totalsum == 0. ) then
                write(LDT_logunit,*) "[WARN] Total vegetation for gridcell, c=",c,", r=",r
                write(LDT_logunit,*) "   has the sum of: ",totalsum
                write(LDT_logunit,*) "  This check performed in routine: create_maskfile "
                ! Set then localmask to 0.0 (water point or undefined)
                localmask(c,r) = 0.0
                cycle
!                call LDT_endrun 
             endif

             if( (vegsum/totalsum) >= LDT_rc%gridcell_water_frac(n) ) then  
               localmask(c,r) = 1.0   ! Designated land points
             else
               localmask(c,r) = 0.0   ! water points
             endif

          end do     ! End nc loop
       end do        ! End nr loop

   end select  ! Type of aggregation condition

!- Generate total number of accounted mask points:
   do r=1,LDT_rc%lnr(n)
      do c=1,LDT_rc%lnc(n)
         if( localmask(c,r) >= 1 ) then
            LDT_rc%nmaskpts(n) = LDT_rc%nmaskpts(n) + 1
         endif
      end do
   end do

end subroutine create_maskfile
