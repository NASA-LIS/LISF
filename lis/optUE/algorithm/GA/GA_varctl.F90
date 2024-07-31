!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module GA_varctl
!BOP
!
! !MODULE: GA_varctl
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

  type, public ::  gactl
     character(len=LIS_CONST_PATH_LEN) :: rfile
     character*100,  allocatable     :: vname(:)
     character*20     :: restart
     integer          :: ngens
     integer          :: npopsize
     integer          :: nparam
     integer          :: nchild
     integer          :: genNo
     integer          :: xoverscheme
     real             :: pcross
     real             :: pcreep
     real             :: pmutate
     integer          :: icreep
     integer          :: ielite
     real,    allocatable :: parmax(:)
     real,    allocatable :: parmin(:)
!     real,    allocatable :: pardel(:)
     integer, allocatable :: npossbl(:)
!     real,    allocatable :: lbp(:)     ! lower bound for parameters
!     real,    allocatable :: inp(:)     ! increment for the parameter array
     integer, allocatable :: nbp(:)     ! number of bits per parameter
     integer          :: sum_nbp 
     integer          :: sum_nposs_param
     integer          :: ngenes     ! number of binary bits for each organism
! kwh: "gene" often string of bits used to define parameter
     integer          :: seed   ! initial seed for the random number generator
     real             :: minFitness
     integer, allocatable :: useSingleParamSet(:) !use a single set of parameters for 
                                           !all grid points. 
     integer, allocatable :: useIntegerValues(:)
  end type gactl
  
end module GA_varctl
