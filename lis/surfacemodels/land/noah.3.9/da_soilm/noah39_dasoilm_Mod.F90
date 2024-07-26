!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module noah39_dasoilm_Mod
!BOP
!
! !MODULE: noah39_dasoilm_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
! 21Oct2018: Mahdi Navari; Sujay Kumar ; Initial Specification 

! !USES:        
  use ESMF
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noah39_dasoilm_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noah39_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_dec
  
  type(dasm_dec), allocatable :: noah39_dasm_struc(:)

contains
!BOP
! 
! !ROUTINE: noah39_dasoilm_init
! \label{noah39_dasoilm_init}
! 
! !INTERFACE:
  subroutine noah39_dasoilm_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    

    implicit none
    integer                :: k
    integer                :: n 
    character(len=LIS_CONST_PATH_LEN) :: modelcdffile(LIS_rc%nnest)
    integer                :: status
    integer                :: ngrid

    if(.not.allocated(noah39_dasm_struc)) then 
       allocate(noah39_dasm_struc(LIS_rc%nnest))
    endif
    
!TBD: SVK
#if 0 
    if(LIS_rc%dascaloption(k).eq."Linear scaling") then 
       call ESMF_ConfigFindLabel(LIS_config,"Noah.3.9 soil moisture CDF file:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'Noah.3.9 soil moisture CDF file: not defined')
       enddo
       
       do n=1,LIS_rc%nnest
       
!Hardcoded for now.
          noah39_dasm_struc(n)%nbins = 100
          
          call LIS_getCDFattributes(modelcdffile(n),&
               noah39_dasm_struc(n)%ntimes, ngrid)
          
          allocate(noah39_dasm_struc(n)%model_xrange(&
               LIS_rc%ngrid(n), noah39_dasm_struc(n)%ntimes, &
               noah39_dasm_struc(n)%nbins))
          allocate(noah39_dasm_struc(n)%model_cdf(&
               LIS_rc%ngrid(n), noah39_dasm_struc(n)%ntimes, &
               noah39_dasm_struc(n)%nbins))
          
          call LIS_readCDFdata(n,&
               noah39_dasm_struc(n)%nbins, &
               noah39_dasm_struc(n)%ntimes, &
               ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               noah39_dasm_struc(n)%model_xrange,&
               noah39_dasm_struc(n)%model_cdf)
       enddo
    endif
#endif

  end subroutine noah39_dasoilm_init
end module noah39_dasoilm_Mod
