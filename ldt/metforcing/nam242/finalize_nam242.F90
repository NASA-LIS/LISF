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
! !MODULE: finalize_nam242
! \label{finalize_nam242}
! 
! !REVISION HISTORY: 
!     Sep 2012: NOHRSC/NOAA: Initial specification
! 
! !INTERFACE:
subroutine finalize_nam242(findex)
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use nam242_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for NAM forcing. 
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing
!  \end{description}
!EOP  
  implicit none

  integer :: n 
  integer :: findex
  
  do n=1,LDT_rc%nnest  
     if(LDT_rc%met_gridtransform(findex).eq."bilinear") then 

        deallocate(nam242_struc(n)%n111)
        deallocate(nam242_struc(n)%n121)
        deallocate(nam242_struc(n)%n211)
        deallocate(nam242_struc(n)%n221)
        deallocate(nam242_struc(n)%w111)
        deallocate(nam242_struc(n)%w121)
        deallocate(nam242_struc(n)%w211)
        deallocate(nam242_struc(n)%w221)
     elseif(LDT_rc%met_gridtransform(findex).eq."budget-bilinear") then 

        deallocate(nam242_struc(n)%n111)
        deallocate(nam242_struc(n)%n121)
        deallocate(nam242_struc(n)%n211)
        deallocate(nam242_struc(n)%n221)
        deallocate(nam242_struc(n)%w111)
        deallocate(nam242_struc(n)%w121)
        deallocate(nam242_struc(n)%w211)
        deallocate(nam242_struc(n)%w221)

        deallocate(nam242_struc(n)%n112)
        deallocate(nam242_struc(n)%n122)
        deallocate(nam242_struc(n)%n212)
        deallocate(nam242_struc(n)%n222)
        deallocate(nam242_struc(n)%w112)
        deallocate(nam242_struc(n)%w122)
        deallocate(nam242_struc(n)%w212)
        deallocate(nam242_struc(n)%w222)
     endif
  enddo
  deallocate(nam242_struc)
end subroutine finalize_nam242
