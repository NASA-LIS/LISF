!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: apply_conusmask
! \label{apply_conusmask}
!
! !INTERFACE:
subroutine apply_conusmask(n, nc,nr,varfield)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use nldas1_forcingMod,  only : nldas1_struc
  
  implicit none

! !ARGUMENTS:   
  integer, intent(in)   :: n 
  integer, intent(in)   :: nc
  integer, intent(in)   :: nr
  real, intent(inout)   :: varfield(nc,nr)
!
! !DESCRIPTION:
!
!EOP
  
  integer               :: c,r
    
  do r=1,nr
     do c=1,nc
        if(nldas1_struc(n)%conusmask(c,r).eq.0.0) then
           varfield(c,r) = LIS_rc%udef
        endif
     enddo
  enddo

end subroutine apply_conusmask
