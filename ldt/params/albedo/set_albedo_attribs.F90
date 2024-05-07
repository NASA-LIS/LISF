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
! !ROUTINE: set_albedo_attribs
!  \label{set_albedo_attribs}
!
! !REVISION HISTORY:
!  19 Aug 2014: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine set_albedo_attribs(n,source)

! !USES:
  use LDT_albedoMod

  implicit none

  integer,         intent(in) :: n
  character(len=*),intent(in) :: source

! !ARGUMENTS: 

! !DESCRIPTION:
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[source]
!     Albedo dataset source
!   \end{description}
!EOP      
!

   select case( source )

    case( "NCEP_LIS", "NCEP_Native", "CLSMF2.5" )
      LDT_albedo_struc(n)%albedo%num_bins = 1
      LDT_albedo_struc(n)%albedo%num_times = 12

    case( "NCEP_NativeQtr" )
      LDT_albedo_struc(n)%albedo%num_bins = 1
      LDT_albedo_struc(n)%albedo%num_times = 4

    case( "CONSTANT" )
      LDT_albedo_struc(n)%albedo%num_bins = 1
      LDT_albedo_struc(n)%albedo%num_times = 12

    case default
      write(*,*) "[ERR] Albedo source not recognized: ",trim(source)
      write(*,*) " Please select: NCEP_LIS, NCEP_Native, NCEP_NativeQtr, " 
      write(*,*) "                CLSMF2.5, or CONSTANT"
      write(*,*) " Program stopping ..."
      stop
!      call LDT_endrun
   end select

end subroutine set_albedo_attribs
