!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: finalize_gdas
! \label{finalize_gdas}
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_gdas(findex)
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use gdas_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GDAS forcing. 
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

        deallocate(gdas_struc(n)%n111)
        deallocate(gdas_struc(n)%n121)
        deallocate(gdas_struc(n)%n211)
        deallocate(gdas_struc(n)%n221)
        deallocate(gdas_struc(n)%w111)
        deallocate(gdas_struc(n)%w121)
        deallocate(gdas_struc(n)%w211)
        deallocate(gdas_struc(n)%w221)
     elseif(LDT_rc%met_gridtransform(findex).eq."budget-bilinear") then 

        deallocate(gdas_struc(n)%n111)
        deallocate(gdas_struc(n)%n121)
        deallocate(gdas_struc(n)%n211)
        deallocate(gdas_struc(n)%n221)
        deallocate(gdas_struc(n)%w111)
        deallocate(gdas_struc(n)%w121)
        deallocate(gdas_struc(n)%w211)
        deallocate(gdas_struc(n)%w221)

        deallocate(gdas_struc(n)%n112)
        deallocate(gdas_struc(n)%n122)
        deallocate(gdas_struc(n)%n212)
        deallocate(gdas_struc(n)%n222)
        deallocate(gdas_struc(n)%w112)
        deallocate(gdas_struc(n)%w122)
        deallocate(gdas_struc(n)%w212)
        deallocate(gdas_struc(n)%w222)
     endif
  enddo
  deallocate(gdas_struc)
end subroutine finalize_gdas
