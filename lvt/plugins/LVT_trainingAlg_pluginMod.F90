!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: LVT_trainingAlg_pluginMod
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining routines that initialize various training algorithms
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  14 Jul 2015;   Sujay Kumar  Initial Specification
! 
!EOP
module LVT_trainingAlg_pluginMod

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LVT_trainingAlg_plugin  
contains
!BOP
! !ROUTINE: LVT_trainingAlg_plugin
!  \label{LVT_trainingAlg_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new training algorithm
!  The interface mandates that the following interfaces be implemented
!  and registered for each training algorithm. 
!
!
! !INTERFACE:
  subroutine LVT_trainingAlg_plugin
    use LVT_pluginIndices
    use LinearRegressionMod

    implicit none
!EOP

    call registertraininginit(trim(LVT_LinearRegressionId)//char(0),&
         linearRegression_init) 
    call registertrainingrun(trim(LVT_LinearRegressionId)//char(0),&
         linearRegression_run)

  end subroutine LVT_trainingAlg_plugin
end module LVT_trainingAlg_pluginMod
