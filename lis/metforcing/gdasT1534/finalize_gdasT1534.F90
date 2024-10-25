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
! !MODULE: finalize_gdasT1534
! \label{finalize_gdasT1534}
! 
! !REVISION HISTORY: 
!  20 June 2014: Sujay Kumar; initial implementation
! 
! !INTERFACE:
subroutine finalize_gdasT1534(findex)
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use gdasT1534_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GDAST1534 forcing. 
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

  do n=1,LIS_rc%nnest  
     if(LIS_rc%met_interp(findex).eq."bilinear") then 

        deallocate(gdasT1534_struc(n)%n111)
        deallocate(gdasT1534_struc(n)%n121)
        deallocate(gdasT1534_struc(n)%n211)
        deallocate(gdasT1534_struc(n)%n221)
        deallocate(gdasT1534_struc(n)%w111)
        deallocate(gdasT1534_struc(n)%w121)
        deallocate(gdasT1534_struc(n)%w211)
        deallocate(gdasT1534_struc(n)%w221)
     elseif(LIS_rc%met_interp(findex).eq."budget-bilinear") then 

        deallocate(gdasT1534_struc(n)%n111)
        deallocate(gdasT1534_struc(n)%n121)
        deallocate(gdasT1534_struc(n)%n211)
        deallocate(gdasT1534_struc(n)%n221)
        deallocate(gdasT1534_struc(n)%w111)
        deallocate(gdasT1534_struc(n)%w121)
        deallocate(gdasT1534_struc(n)%w211)
        deallocate(gdasT1534_struc(n)%w221)

        deallocate(gdasT1534_struc(n)%n112)
        deallocate(gdasT1534_struc(n)%n122)
        deallocate(gdasT1534_struc(n)%n212)
        deallocate(gdasT1534_struc(n)%n222)
        deallocate(gdasT1534_struc(n)%w112)
        deallocate(gdasT1534_struc(n)%w122)
        deallocate(gdasT1534_struc(n)%w212)
        deallocate(gdasT1534_struc(n)%w222)
     endif
  enddo
  deallocate(gdasT1534_struc)
end subroutine finalize_gdasT1534
