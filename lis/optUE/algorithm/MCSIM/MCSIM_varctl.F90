!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module MCSIM_varctl
!BOP
!
! !MODULE: MCSIM_varctl
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
  implicit none
  PRIVATE

  type, public ::  mcsimctl
     character*100    :: decspaceAttribsFile
     character*100    :: rfile
     character*20     :: restart
     integer          :: iterNo
     integer          :: maxiter
     integer          :: nparam
!     integer          :: restart
     real,    allocatable :: parmax(:)
     real,    allocatable :: parmin(:)
     character*100, allocatable :: vname(:)
     integer          :: seed   ! initial seed for the random number generator
  end type mcsimctl
  
end module MCSIM_varctl
