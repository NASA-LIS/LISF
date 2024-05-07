!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module SCEUA_varctl
!BOP
!
! !MODULE: SCEUA_varctl
!
! !DESCRIPTION:
!  
!  Module for specifying variables used to control the Shuffled Complex
!  Evolution implementation. 
!
! !REVISION HISTORY:
!  09 Jun 2009; Soni Yatheendradas; Initial Specification
!
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  IMPLICIT NONE
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  type, public ::  sceuactl
     character(len=LIS_CONST_PATH_LEN) :: decspaceAttribsFile ! Decision Space Attributes File
     integer          :: restart ! SCEUA start mode
     character(len=LIS_CONST_PATH_LEN) :: rfile ! SCEUA restart file
     INTEGER          :: MaxNFuncEvals ! Maximum Number of Func. Evals. before Optimization Terminates
     INTEGER          :: NShufflesForMinChange ! Number of Shuffles to Terminate Optimization if Criterion Change less than Minimum
     REAL             :: MinChangeInNShuffles ! Minimum Fractional Criterion Change in Specified Shuffles to Continue Optimization
     INTEGER          :: NComplexes ! Number of SCEUA Optimization Complexes
     INTEGER          :: iseed ! Initial seed should be between 1 to 10 that is then changed to corresponding value from 1st 10 prime numbers. Zero and >10 set to 10th, -ve set to 1st. 
     INTEGER          :: iDefault ! Whether to User-specify the SCEUA Control Parameters
     integer          :: nparam ! Number of optimized parameters
     real*8,  allocatable :: parmax(:)
     real*8,  allocatable :: parmin(:)
     INTEGER          :: NPointsInComplex ! Number of Points in a SCEUA Complex (i.e., Parents)
     INTEGER          :: NPointsInComplexPlus4 ! NPointsInComplexPlus4 is NPointsInComplex plus 4 (i.e., In each Complex, Parents and 4 potential children) 
                                               ! NPointsInComplexPlus4=LIS_rc%nensem(n)/NComplexes if IDefault is 1
                                               ! LIS_rc%nensem(n)=NComplexes*NPointsInComplexPlus4 if IDefault is 0, where NPointsInComplexPlus4=(2*nparam+1)+4
     INTEGER          :: NPointsInSubcomplex ! Number of Points in a Subcomplex 
     INTEGER          :: NStepsBeforeShuffle ! Number of Evolution Steps before Shuffle for a Complex
!     INTEGER          :: NMinComplexes ! This option not implemented in LIS
     INTEGER          :: InitialPointFlag ! Whether Include Initial Point in Population
!     integer          :: npopsize ! This not used, since anyway equal to LIS_rc%nensem(n)
     INTEGER          :: iseed1 ! Corrected initial seed
     real*8,  allocatable :: pardel(:) ! Parameter range (delta)
     integer          :: SumReplacedFuncEvals ! How many REPLACED PARENT function/LSM evaluations yet in each spatial tile?
     integer          :: SumNewMemberFuncEvals ! How many NEW member function/LSM evaluations yet in each spatial tile?
     integer          :: SumAllFuncEvals ! How many OVERALL function/LSM evaluations yet in each spatial tile?
     INTEGER          :: iEvolutionLoop ! Evolution loop number
     INTEGER, POINTER :: isDecSpaceConvergedFlagArray(:) ! Array of sceuastruc%isDecSpaceConvergedFlag over all spatial tiles 
     INTEGER, POINTER :: isObjFuncConvergedFlagArray(:) ! Array of sceuastruc%isObjFuncConvergedFlag over all spatial tiles 
     INTEGER          :: isDecSpaceConvergAllSpatialTiles ! If sceuastruc%isDecSpaceConvergedFlag true for all spatial tiles 
     INTEGER          :: isObjFuncConvergAllSpatialTiles ! If sceuastruc%isObjFuncConvergedFlag true for all spatial tiles 
     REAL*8,  allocatable :: ParentParsets(:,:) ! Parent part of parsets for a given spatial tile 
     REAL*8,  allocatable :: ParentObjFuncs(:) ! Parent part of objfuncs for a given spatial tile
     REAL*8,  allocatable :: complexParsets(:,:) ! Co-ordinates of points in a complex 
     REAL*8,  allocatable :: complexObjFuncs(:) ! Objective function values of points in a complex 
     INTEGER          :: ran1_iff ! iff in ran1
     REAL             :: maxObjFunc ! For use with LIS_feasibleSpace 
  end type sceuactl
  
!EOP

end module SCEUA_varctl
