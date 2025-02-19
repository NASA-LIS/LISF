!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module NoahMP50_dasoilm_Mod
!BOP
!
! !MODULE: NoahMP50_dasoilm_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:

! 15 Dec 2018: Mahdi Navari, Sujay Kumar ; Modified for NoahMP401 !
! May 2023: Cenlin He; modified for refactored NoahMP v5 and later

! !USES:        
  use LIS_coreMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: NoahMP50_dasoilm_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noahmp50_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_dec
  
  type(dasm_dec), allocatable :: noahmp50_dasm_struc(:)

contains
!BOP
! 
! !ROUTINE: NoahMP50_dasoilm_init
! \label{NoahMP50_dasoilm_init}
! 
! !INTERFACE:
  subroutine NoahMP50_dasoilm_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    

    implicit none
    integer                :: k
#if 0
    integer                :: n
    integer                :: status
    integer                :: ngrid
    character(len=LIS_CONST_PATH_LEN) :: modelcdffile(LIS_rc%nnest)
#endif

    if(.not.allocated(noahmp50_dasm_struc)) then 
       allocate(noahmp50_dasm_struc(LIS_rc%nnest))
    endif
    
!TBD: SVK
#if 0 
    if(LIS_rc%dascaloption(k).eq."Linear scaling") then 
       call ESMF_ConfigFindLabel(LIS_config,"Noah-MP.5.0 soil moisture CDF file:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'Noah-MP.5.0 soil moisture CDF file: not defined')
       enddo
       
       do n=1,LIS_rc%nnest
       
!Hardcoded for now.
          noahmp50_dasm_struc(n)%nbins = 100
          
          call LIS_getCDFattributes(modelcdffile(n),&
               noahmp50_dasm_struc(n)%ntimes, ngrid)
          
          allocate(noahmp50_dasm_struc(n)%model_xrange(&
               LIS_rc%ngrid(n), noahmp50_dasm_struc(n)%ntimes, &
               noahmp50_dasm_struc(n)%nbins))
          allocate(noahmp50_dasm_struc(n)%model_cdf(&
               LIS_rc%ngrid(n), noahmp50_dasm_struc(n)%ntimes, &
               noahmp50_dasm_struc(n)%nbins))
          
          call LIS_readCDFdata(n,&
               noahmp50_dasm_struc(n)%nbins, &
               noahmp50_dasm_struc(n)%ntimes, &
               ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               noahmp50_dasm_struc(n)%model_xrange,&
               noahmp50_dasm_struc(n)%model_cdf)
       enddo
    endif
#endif

  end subroutine NoahMP50_dasoilm_init
end module NoahMP50_dasoilm_Mod
