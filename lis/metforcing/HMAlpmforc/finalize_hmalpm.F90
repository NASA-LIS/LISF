!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: finalize_hmalpm
! \label{finalize_hmalpm}
!
! !REVISION HISTORY: 
! 23 Dec 2019: Sujay Kumar, initial code
! 
! !INTERFACE:

subroutine finalize_hmalpm(findex)

! !USES:
  use LIS_coreMod,     only : LIS_rc
  use hmalpm_forcingMod, only : hmalpm_struc
!
! !DESCRIPTION:
!  Routine to cleanup ERA5 forcing related memory allocations.   
! 
!EOP
  implicit none

  integer :: findex
  integer :: n

  do n=1,LIS_rc%nnest
    select case( LIS_rc%met_interp(findex) )
     case( "bilinear" )
       deallocate(hmalpm_struc(n)%n111)
       deallocate(hmalpm_struc(n)%n121)
       deallocate(hmalpm_struc(n)%n211)
       deallocate(hmalpm_struc(n)%n221)
       deallocate(hmalpm_struc(n)%w111)
       deallocate(hmalpm_struc(n)%w121)
       deallocate(hmalpm_struc(n)%w211)
       deallocate(hmalpm_struc(n)%w221)

     case( "budget-bilinear" )
       deallocate(hmalpm_struc(n)%n111)
       deallocate(hmalpm_struc(n)%n121)
       deallocate(hmalpm_struc(n)%n211)
       deallocate(hmalpm_struc(n)%n221)
       deallocate(hmalpm_struc(n)%w111)
       deallocate(hmalpm_struc(n)%w121)
       deallocate(hmalpm_struc(n)%w211)
       deallocate(hmalpm_struc(n)%w221)
       deallocate(hmalpm_struc(n)%n112)
       deallocate(hmalpm_struc(n)%n122)
       deallocate(hmalpm_struc(n)%n212)
       deallocate(hmalpm_struc(n)%n222)
       deallocate(hmalpm_struc(n)%w112)
       deallocate(hmalpm_struc(n)%w122)
       deallocate(hmalpm_struc(n)%w212)
       deallocate(hmalpm_struc(n)%w222)

     case( "neighbor" )
       deallocate(hmalpm_struc(n)%n113)
    end select
 enddo
 deallocate(hmalpm_struc)

end subroutine finalize_hmalpm
