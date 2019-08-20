!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module NoahMP401_dasoilm_Mod
!BOP
!
! !MODULE: NoahMP401_dasoilm_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:

! 15 Dec 2018: Mahdi Navari, Sujay Kumar ; Modified for NoahMP401 !

! !USES:        
  use ESMF
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: NoahMP401_dasoilm_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noahmp401_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_dec
  
  type(dasm_dec), allocatable :: noahmp401_dasm_struc(:)

contains
!BOP
! 
! !ROUTINE: NoahMP401_dasoilm_init
! \label{NoahMP401_dasoilm_init}
! 
! !INTERFACE:
  subroutine NoahMP401_dasoilm_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    

    implicit none
    integer                :: k
    integer                :: n 
    character*100          :: modelcdffile(LIS_rc%nnest)
    integer                :: status
    integer                :: ngrid

    if(.not.allocated(noahmp401_dasm_struc)) then 
       allocate(noahmp401_dasm_struc(LIS_rc%nnest))
    endif
    
!TBD: SVK
#if 0 
    if(LIS_rc%dascaloption(k).eq."Linear scaling") then 
       call ESMF_ConfigFindLabel(LIS_config,"Noah-MP.3.6 soil moisture CDF file:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'Noah-MP.3.6 soil moisture CDF file: not defined')
       enddo
       
       do n=1,LIS_rc%nnest
       
!Hardcoded for now.
          noahmp401_dasm_struc(n)%nbins = 100
          
          call LIS_getCDFattributes(modelcdffile(n),&
               noahmp401_dasm_struc(n)%ntimes, ngrid)
          
          allocate(noahmp401_dasm_struc(n)%model_xrange(&
               LIS_rc%ngrid(n), noahmp401_dasm_struc(n)%ntimes, &
               noahmp401_dasm_struc(n)%nbins))
          allocate(noahmp401_dasm_struc(n)%model_cdf(&
               LIS_rc%ngrid(n), noahmp401_dasm_struc(n)%ntimes, &
               noahmp401_dasm_struc(n)%nbins))
          
          call LIS_readCDFdata(n,&
               noahmp401_dasm_struc(n)%nbins, &
               noahmp401_dasm_struc(n)%ntimes, &
               ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               noahmp401_dasm_struc(n)%model_xrange,&
               noahmp401_dasm_struc(n)%model_cdf)
       enddo
    endif
#endif

  end subroutine NoahMP401_dasoilm_init
end module NoahMP401_dasoilm_Mod
