!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LM_varctl
!BOP
!
! !MODULE: LM_varctl
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

  type, public ::  lmctl
     character(len=LIS_CONST_PATH_LEN) :: decspaceAttribsFile
     character(len=LIS_CONST_PATH_LEN) :: rfile
     integer          :: restart
     integer          :: nparam     !nparam in GA code
     real,    allocatable :: parmax(:)
     real,    allocatable :: parmin(:)
     character*100, allocatable :: vname(:)
     real             :: EPSMCH  !machine precision
     integer          :: Mmax  !initial run will count obs; this will allow for first pass storage
     integer          :: MAXITER
     integer          :: ITER  !lm_struc also has iteration as some grid points may converge before others
     integer          :: MODE
     integer          :: NPRINT  !we will use differently than minpack
     real             :: FTOL    !
     real             :: XTOL    !
     real             :: GTOL    !
     real             :: EPSFCN  
     real             :: FACTOR  
     logical          :: AllTerminated
     logical          :: check
  end type lmctl
  
end module LM_varctl
