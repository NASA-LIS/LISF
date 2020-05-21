!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: assignirrigtype
!  \label{assignirrigtype}
!
! !REVISION HISTORY:
!  14  Nov 2019: H. Beaudoing; Initial implementation
!
! This routine converts crop water source information from GRIPC to 
! irrigation types of sprinkler, drip, and flood.  The irrigated croplands
! are assigned one of sprinkler, drip, or sprinkler+drip options spcified
! in ldt.config.  For sprinkler+drip, frction of sprinker can be set in
! ldt.config as well.  The sum of three types equals 1.
! Default type is set to flood.
!
!  The GRIPC crop water source legend is:
!    Undefined           = 0
!    Rain-fed croplands  = 1
!    Irrigated croplands = 2
!    Paddy croplands     = 3
!    Not cropped         = 4
! OUTPUT IRRIGTYPE legend is:
!    Sprinkler = 1
!    Drip      = 2
!    Flood     = 3
!
! !INTERFACE:
subroutine assignirrigtype( source, n, typeopt, factor, in_fgrd, &
                           in_num_bins, fgrd, num_bins)

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_endrun
  use LDT_paramDataMod
  use LDT_irrigationMod

  implicit none
  character(len=*),intent(in) :: source
  integer, intent(in) :: n
  character(len=*),intent(in) :: typeopt
  real, intent(in)    :: factor
  integer, intent(in) :: in_num_bins
  real, intent(in)    :: in_fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),in_num_bins)
  integer, intent(in) :: num_bins
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real                :: alltype, FF, XF
  integer             :: c,r

  if ( trim(source) .eq. "GRIPC" ) then
   fgrd = 0.
   where (LDT_LSMparam_struc(n)%landmask%value(:,:,1) > 0.5)
       fgrd(:,:,1) = 0.   ! sprinkler
       fgrd(:,:,2) = 0.   ! drip
       fgrd(:,:,3) = 1.   ! default is flood
   end where
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         if ( LDT_LSMparam_struc(n)%landmask%value(c,r,1) > 0.5 ) then
           if ( in_fgrd(c,r,2) .gt. 0 ) then  ! irrigated
             if ( in_fgrd(c,r,3) .gt. 0 ) then   ! paddy
               alltype = in_fgrd(c,r,2) + in_fgrd(c,r,3)
               FF = in_fgrd(c,r,3) / alltype
               XF = in_fgrd(c,r,2) / alltype
               select case ( typeopt )
                 case ( "sprinkler" ) 
                   fgrd(c,r,1) = XF
                   fgrd(c,r,2) = 0.
                   fgrd(c,r,3) = FF
                 case ( "drip" )
                   fgrd(c,r,1) = 0.
                   fgrd(c,r,2) = XF
                   fgrd(c,r,3) = FF
                 case ( "sprinkler+drip" )
                   fgrd(c,r,1) = XF*factor
                   fgrd(c,r,2) = XF*(1.0-factor)
                   fgrd(c,r,3) = FF
               end select
             else   ! no paddy
               fgrd(c,r,3) = 0.
               select case ( typeopt )
                 case ( "sprinkler" ) 
                   fgrd(c,r,1) = 1.
                   fgrd(c,r,2) = 0.
                 case ( "drip" )
                   fgrd(c,r,1) = 0.
                   fgrd(c,r,2) = 1.
                 case ( "sprinkler+drip" )
                   fgrd(c,r,1) = factor
                   fgrd(c,r,2) = (1.0-factor)
               end select
             endif
           endif   ! irrigated
         endif  ! landmask > 0
      enddo
   enddo

  else
    write(LDT_logunit,*)"[ERR] Assign irrigation types currently not supported"
    write(LDT_logunit,*)"      for ",trim(source)
    write(LDT_logunit,*)"      only GRIPC ... stopping...."
    call LDT_endrun
  endif
  
end subroutine assignirrigtype
