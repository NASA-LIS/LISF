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
! !MODULE: finalize_RFE2gdas
!  \label{finalize_RFE2gdas}
!
! !REVISION HISTORY: 
! 2 June 2010: Soni Yatheendradas; Initial LDT version for FEWSNET
!
! !INTERFACE:
subroutine finalize_RFE2gdas(findex)
! !USES:
  use LDT_coreMod,         only : LDT_rc
  use RFE2gdas_forcingMod, only : RFE2gdas_struc
!
! !DESCRIPTION:
!  Routine to cleanup RFE2gdas forcing related memory allocations.   
! 
!EOP

  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LDT_rc%nnest

    select case( LDT_rc%met_gridtransform(findex) )

     case( "average" )   ! Upscaling 
#if 0
! Soni's older code
       deallocate(RFE2gdas_struc(n)%stc)
       deallocate(RFE2gdas_struc(n)%str)
       deallocate(RFE2gdas_struc(n)%enc)
       deallocate(RFE2gdas_struc(n)%enr)
#endif 
       deallocate(RFE2gdas_struc(n)%n111)

     case( "bilinear" )

       deallocate(RFE2gdas_struc(n)%n111)
       deallocate(RFE2gdas_struc(n)%n121)
       deallocate(RFE2gdas_struc(n)%n211)
       deallocate(RFE2gdas_struc(n)%n221)
       deallocate(RFE2gdas_struc(n)%w111)
       deallocate(RFE2gdas_struc(n)%w121)
       deallocate(RFE2gdas_struc(n)%w211)
       deallocate(RFE2gdas_struc(n)%w221)

     case( "budget-bilinear" )

       deallocate(RFE2gdas_struc(n)%n112)
       deallocate(RFE2gdas_struc(n)%n122)
       deallocate(RFE2gdas_struc(n)%n212)
       deallocate(RFE2gdas_struc(n)%n222)
       deallocate(RFE2gdas_struc(n)%w112)
       deallocate(RFE2gdas_struc(n)%w122)
       deallocate(RFE2gdas_struc(n)%w212)
       deallocate(RFE2gdas_struc(n)%w222)

     case( "neighbor" )

       deallocate(RFE2gdas_struc(n)%n113)

    end select

  enddo

  deallocate(RFE2gdas_struc)

end subroutine finalize_RFE2gdas


