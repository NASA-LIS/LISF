!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module ES_varctl
!BOP
!
! !MODULE: ES_varctl
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

  type, public ::  esctl
     character*100    :: decspaceAttribsFile
     integer          :: nparam
     integer          :: genNo
     real,    allocatable :: parmax(:)
     real,    allocatable :: parmin(:)
     real,    allocatable :: pardel(:)
     integer, allocatable :: npossbl(:)
     integer, allocatable :: parm_count(:)
     integer          :: ntotal
  end type esctl
  
end module ES_varctl
