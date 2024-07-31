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
! !ROUTINE: texture
! \label{texture}
! 
! !REVISION HISTORY: 
! 7 Mar 2001 : Matt Rodell, Initial Specification
! 
! !INTERFACE: 
subroutine texture (sand,silt,clay,tex)
! !USES: 
  implicit none
! !ARGUMENTS:   
  real,    intent(in)    :: sand,silt,clay
  integer, intent(inout) :: tex
  
! !DESCRIPTION: 
!    Finds the USDA texture class based on percentages of sand, silt,
!    and clay.
!EOP
  integer :: i
  real :: frac,sa,si,cl
  
  TEX = 0
!    Test for ocean points.
  if (CLAY .lt. 0.00) then
     tex = -9999.0
  else
!      Adjust values at points whose fractions don't sum to 1.00.
     frac = CLAY + & 
          SAND + & 
          SILT
     cl = CLAY / frac
     sa = SAND / frac
     si = SILT / frac
     !      Identify texture class.
     if ((cl .ge. 0.40) .and. (si .lt. 0.40) .and.  &
          (sa .lt. 0.44)) then
        tex = 1  !CLAY
     end if
     
     if ((cl .ge. 0.40) .and. (si .ge. 0.40)) then
        tex = 2  !SILTY CLAY
     end if
     
     if ((cl .ge. 0.36) .and. (sa .ge. 0.44)) then
        tex = 3  !SANDY CLAY
     end if
     
     if ((cl .ge. 0.28) .and. (cl .lt. 0.40) .and.  &
          (sa .ge. 0.20) .and. (sa .lt. 0.44)) then
        tex = 4  !CLAY LOAM
     end if
     
     if ((cl .ge. 0.28) .and. (cl .lt. 0.40) .and.  &
          (sa .lt. 0.20)) then
        tex = 5  !SILTY CLAY LOAM
     end if
     
     if ((cl .ge. 0.20) .and. (cl .lt. 0.36) .and. &
          (si .lt. 0.28) .and. (sa .ge. 0.44)) then
        tex = 6  !SANDY CLAY LOAM
     end if
     
     if ((cl .ge. 0.08) .and. (cl .lt. 0.28) .and.  &
          (si .ge. 0.28) .and. (si .lt. 0.50) &
          .and. (sa .lt. 0.52)) then
        tex = 7  !LOAM
     end if
     
     if ((cl .lt. 0.28) .and. (si .ge. 0.50)) then
        tex = 8  !SILTY LOAM (and SILT)
     end if
     
     if ((sa - cl) .gt. 0.70) then
        if ((sa - (0.3 * cl)) .gt. 0.87) then
           tex = 11  !SAND
        else
           tex = 10  !LOAMY SAND
        end if
     end if
     
     if (tex .eq. 0) then
        tex = 9  !SANDY LOAM
     end if
     
  end if  !CLAY
  return
end subroutine texture
