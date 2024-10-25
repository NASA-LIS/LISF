!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module RWMCMC_varctl
!BOP
!
! !MODULE: RWMCMC_varctl
!
! !DESCRIPTION:
!  
!  Module for specifying variables used to control the Genetic Algorithm
!  implementation. 
!
! !REVISION HISTORY:
!  04 Feb 2008; Sujay Kumar; Initial Specification
!
!EOP
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  PRIVATE

  type, public ::  mcmcctl
     character(len=LIS_CONST_PATH_LEN) :: decspaceAttribsFile
     character(len=LIS_CONST_PATH_LEN) :: rfile
     character(len=LIS_CONST_PATH_LEN) :: initfile
     integer          :: iterNo
     integer          :: maxiter
     integer          :: nparam
     integer          :: restart
     real             :: minfitness
     real             :: pert_spread
     real,    allocatable :: parmax(:)
     real,    allocatable :: parmin(:)
     character*100, allocatable :: vname(:)
     integer          :: seed   ! initial seed for the random number generator
     logical          :: zerothRun
  end type mcmcctl
  
end module RWMCMC_varctl
