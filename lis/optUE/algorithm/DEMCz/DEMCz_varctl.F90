!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module DEMCz_varctl
!BOP
!
! !MODULE: DEMCz_varctl
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

  type, public ::  demczctl
     character*100    :: decspaceAttribsFile
     character*100    :: rfile
     character*100    :: optrfile  ! to start from GA run
     integer          :: iterNo
     integer          :: maxiter
     integer          :: nparam
     character*20     :: restart
     real             :: pert_spread
     real             :: minfitness
     real             :: modehopfreq  !eg, 0.1 --> 1 iteration out of every 10 iterations
     real,    allocatable :: parmax(:)
     real,    allocatable :: parmin(:)
     character*100, allocatable :: vname(:)
     integer, allocatable :: useSingleParamSet(:) !use a single set of parameters for 
     integer          :: seed   ! initial seed for the random number generator
     logical          :: overall_check
     logical          :: zerothRun
  end type demczctl
  
end module DEMCz_varctl
