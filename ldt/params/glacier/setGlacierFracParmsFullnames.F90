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
! !ROUTINE: setGlacierFracParmsFullnames
!  \label{setGlacierFracParmsFullnames}i
!
! !REVISION HISTORY:
!  21 Jun 2020: Mahdi Navari; This code is based on the subroutine 
!                setGfracParmsFullnames.F90 initially implemented by K. Arsenault.
!  26 Jan 2021: Mahdi Navari; Modified for glacier fraction
!
! !INTERFACE:
subroutine setGlacierFracParmsFullnames(n,datatype,source)

! !USES:
  use LDT_glacierFractionMod
  use LDT_logMod
  use LDT_paramDataMod

  implicit none
  integer,         intent(in) :: n
  character(len=*),intent(in) :: datatype
  character(len=*),intent(in) :: source

! !ARGUMENTS: 

! !DESCRIPTION:
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[datatype]
!     Irrig data type
!   \item[source]
!     Irrig dataset source
!   \end{description}
!EOP      
!

   select case( datatype )


    case( "glacierfrac" )
      select case( source )
        case( "GLIMS" )
          LDT_glacierfrac_struc(n)%glacierfrac%standard_name =&
              "GLIMS mask glacier gridcell fraction"

      end select

    case default
      write (LDT_logunit, *) "[ERR] Glacier fraction data type not recognized: ",trim(source)
      call LDT_endrun()
   end select

end subroutine setGlacierFracParmsFullnames
