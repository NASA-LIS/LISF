!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_metforcScale_pluginMod
!BOP
!
! !MODULE: LDT_metforcScale_pluginModMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining routines to scale meteorological forcing datasets 
!   either temporally and/or spatially (up or down). 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
!   
! !REVISION HISTORY: 
!  11 Dec 2003  Sujay Kumar   Initial Specification
!  11 Dec 2014  KR Arsenault  Added C-function tables for downscaling
! 
!EOP  

  use LDT_pluginIndices

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LDT_timedscale_plugin

 contains

!
!BOP
! !ROUTINE: LDT_timedscale_plugin
!  \label{LDT_timedscale_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing forcing time/space scaling methods.
! The interface mandates that the following routines be implemented
! and registered for each method/approach source. 
!
!  \begin{description}
!  \item[read the forcing downscaling option]      
!      Routines to call methods of metforcing downscaling in time.
!      (to be registered using {\tt applytimedscale} and later called
!       using the generic {\tt applytimedscale} method)
!  \end{description}
! 
!  The user-defined functions are included in the registry using 
!  two indices. 
!  The methods should be defined in the registry as follows:
!  
!    call applytimedscale(ldt%source,ldt%nest)
!
!  \end{verbatim}
!
  subroutine LDT_timedscale_plugin
!
!EOP
    use LDT_SimpleWeight_TdscaleMod, only : LDT_SimpleWeight_Tdscale

! !USES:

!- Apply temporal downscaling (e.g., simple weighting):

    call registerapplytimedscale(trim(LDT_simplewgtId)//char(0),&
         LDT_SimpleWeight_Tdscale)

  end subroutine LDT_timedscale_plugin

end module LDT_metforcScale_pluginMod
