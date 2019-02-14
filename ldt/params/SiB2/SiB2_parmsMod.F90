!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module SiB2_parmsMod
!BOP
!
! !MODULE: SiB2_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read greenness fraction
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  greenness fraction climatology data and allows the users to 
!  specify the frequency of climatology (in months). 
!  The climatological data is temporally interpolated  
!  between months to the current simulation date. 
!
! !REVISION HISTORY:
!
!  04 Nov 2013: K. Arsenault: Added layers for SiB2 model
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: SiB2Parms_init    !allocates memory for required structures
  public :: SiB2Parms_writeHeader
  public :: SiB2Parms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: SiB2_struc

  type, public :: noah_type_dec
!  - SiB2 parameters:
     character*100 :: sib2parmsdir

     ! -  SiB2 LSM-specific:
! -  SiB2 model-specific:
     type(LDT_paramEntry) :: sib2        ! SiB2 model parameters (collective)
     type(LDT_paramEntry) :: itypedef    ! 0 -- Type ??
     type(LDT_paramEntry) :: z2          ! 1   Height of canopy top (m)
     type(LDT_paramEntry) :: z1          ! 2   Height of canopy bottom (m)
     type(LDT_paramEntry) :: vcover      ! 3   Canopy cover fraction (-)
     type(LDT_paramEntry) :: chil        ! 4   Leaf angle distribution factor (-)
     type(LDT_paramEntry) :: sodep       ! 5   Total soil depth (m)
     type(LDT_paramEntry) :: rootd       ! 6   Rooting depth (m)
     type(LDT_paramEntry) :: phc         ! 7   1/2 critical leaf water potential limit (-)
     type(LDT_paramEntry) :: tran1       ! 8   Leaf transmittance (VIS, live)
     type(LDT_paramEntry) :: tran2       ! 9   Leaf transmittance (VIS, dead)
     type(LDT_paramEntry) :: tran3       ! 10  Leaf transmittance (NIR, live)
     type(LDT_paramEntry) :: tran4       ! 11  Leaf transmittance (NIR, dead)
     type(LDT_paramEntry) :: ref1        ! 12  Leaf reflectance (VIS, live)
     type(LDT_paramEntry) :: ref2        ! 13  Leaf reflectance (VIS, dead)
     type(LDT_paramEntry) :: ref3        ! 14  Leaf reflectance (NIR, live)
     type(LDT_paramEntry) :: ref4        ! 15  Leaf reflectance (NIR, dead)
     type(LDT_paramEntry) :: vmax0       ! 16  Max RuBisCO capacity at canopy top (mol m-2 s-1)
     type(LDT_paramEntry) :: effcon      ! 17  Quantum yield of photosynthesis MolCO2 (mol mol-1)
     type(LDT_paramEntry) :: gradm       ! 18  Stomatal conductance slope parameter (-)
     type(LDT_paramEntry) :: binter      ! 19  Minimum stomtal conductance (mol m-2 s-1)
     type(LDT_paramEntry) :: atheta      ! 20  Light and RuBisCO (wc,we) coupling parameter (-)
     type(LDT_paramEntry) :: btheta      ! 21  Light, RuBisCO, CHO sink( wc,we,ws) coupling parameter (-) 
     type(LDT_paramEntry) :: trda        ! 22  Temperature inhibition coef. for leaf respiration (K-1)
     type(LDT_paramEntry) :: trdm        ! 23  1/2 critical temperature inhib. for leaf respiration (K)
     type(LDT_paramEntry) :: trop        ! 24  Q10 temperature coeff. for physiology (298.16K)
     type(LDT_paramEntry) :: respcp      ! 25  Leaf respiration fraction for Vmax (mol m-2 s-1)
     type(LDT_paramEntry) :: slti        ! 26  Slope of low temperature inhibition function (K-1)
     type(LDT_paramEntry) :: shti        ! 27  Slope of high temperature inhibition function (K-1)
     type(LDT_paramEntry) :: hltii       ! 28  1/2 point of low temperature inhibition function (K)
     type(LDT_paramEntry) :: hhti        ! 29  1/2 point of high temperature inhibition function (K)
     type(LDT_paramEntry) :: soref1      ! 30  Soil reflectance (VIS; -)
     type(LDT_paramEntry) :: soref2      ! 31  Soil reflectance (NIR; -)
     type(LDT_paramEntry) :: bee         ! 32  Soil wetness exponent (-)
     type(LDT_paramEntry) :: phsat       ! 33  Soil tension at saturation (m)
     type(LDT_paramEntry) :: satco       ! 34  Hydraulic conductivity at saturation (m s-1)
     type(LDT_paramEntry) :: poros       ! 35  Soil porosity (-)
     type(LDT_paramEntry) :: islope      ! 36  Cosine of mean slope (-)
     type(LDT_paramEntry) :: wopt        ! 37  Optimum soil water perc. for respiration (soiltype dependent) (-)
     type(LDT_paramEntry) :: wsat        ! 38  Respiration rate at soil water saturation coeff. (-)
     type(LDT_paramEntry) :: zm          ! 39  Skewness exponent of respiration vs. soil water (-)

     
  end type noah_type_dec

  type(noah_type_dec), allocatable :: SiB2_struc(:)


contains

!BOP
! 
! !ROUTINE: SiB2Parms_init
! \label{SiB2Parms_init}
! 
! !INTERFACE:
  subroutine SiB2Parms_init
! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify, LDT_endrun
!    use LDT_paramOptCheckMod, only: LDT_noahparmsOptChecks, &
!                       LDT_gridOptChecks

!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the Simple Biosphere 2 (SiB2) model parameter datasets.
!
! Model Reference: 
!  Sellers, P.J., D.A. Randall, G.J. Collatz, J.A. Berry, C.B. Field, 
!   D.A. Dazlich, C. Zhang, G.D. Collelo, and L. Bounoua, 1996. A revised 
!   land surface parameterization (SiB2) for atmospheric GCMs Model 
!   formulation. Journal of Climate, 9, 676.
!
!  The routines invoked are: 
!  \begin{description}
!   \item[sib2Parmssetup](\ref{sib2Parmssetup}) \newline
!    calls the registry to invoke the sib2Parms setup methods. 
!  \end{description}
!
!EOP
    implicit none

     real, allocatable                   :: sib2parms_gridDesc(:,:)
     character*50                    :: sib2parms_proj
     character*50, allocatable           :: sib2parms_gridtransform(:)

    integer  :: n
    integer  :: c,r, m
    integer  :: rc
    type(LDT_fillopts) :: sib2

    ! _____________________________________________________________________

    allocate(SiB2_struc(LDT_rc%nnest))
    do n=1,LDT_rc%nnest
       ! - SiB2 parameters:
       call set_param_attribs(SiB2_struc(n)%sib2,"SIB2")

    enddo
    if( SiB2_struc(1)%sib2%selectOpt == 1 ) then
       write(LDT_logunit,*)" - - - - - - - - - - SiB2 LSM Parameters - - - - - - - - - - - - -"

       allocate(sib2parms_gridDesc(LDT_rc%nnest,20))
       allocate(sib2parms_gridtransform(LDT_rc%nnest))

       sib2parms_gridDesc = 0.

       do n = 1, LDT_rc%nnest
          if( SiB2_struc(n)%sib2%selectOpt == 1 ) then 

             SiB2_struc(n)%sib2%vlevels = SiB2_struc(n)%sib2%num_bins

             ! Fill in derived parameter entries:
             ! ( input_parmattribs -> output_parmattribs ) 
             call populate_param_attribs( "SIB2_Z2", &
                  "Height of canopy top","m", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%z2 )

             allocate(SiB2_struc(n)%z2%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%z2%vlevels))       
             SiB2_struc(n)%z2%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_Z1", &
                  "Height of canopy bottom","m", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%z1 )

             allocate(SiB2_struc(n)%z1%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%z1%vlevels)) 
             SiB2_struc(n)%z1%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_VCOVER", &
                  "Canopy cover fraction","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%vcover )

             allocate(SiB2_struc(n)%vcover%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%vcover%vlevels)) 
             SiB2_struc(n)%vcover%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_CHIL", &
                  "Leaf angle distribution factor","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%chil )

             allocate(SiB2_struc(n)%chil%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%chil%vlevels)) 
             SiB2_struc(n)%chil%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_SODEP", &
                  "Total soil depth","m", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%sodep )

             allocate(SiB2_struc(n)%sodep%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%sodep%vlevels)) 
             SiB2_struc(n)%sodep%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_ROOTD", &
                  "Rooting depth","m", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%rootd )

             allocate(SiB2_struc(n)%rootd%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%rootd%vlevels)) 
             SiB2_struc(n)%rootd%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_PHC", &
                  "1/2 critical leaf water potential limit","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%phc )

             allocate(SiB2_struc(n)%phc%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%phc%vlevels)) 
             SiB2_struc(n)%phc%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_TRAN1", &
                  "Leaf transmittance (VIS, live)","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%tran1 )

             allocate(SiB2_struc(n)%tran1%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%tran1%vlevels)) 
             SiB2_struc(n)%tran1%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_TRAN2", &
                  "Leaf transmittance (VIS, dead)","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%tran2 )

             allocate(SiB2_struc(n)%tran2%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%tran2%vlevels)) 
             SiB2_struc(n)%tran2%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_TRAN3", &
                  "Leaf transmittance (NIR, live)","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%tran3 )

             allocate(SiB2_struc(n)%tran3%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%tran3%vlevels)) 
             SiB2_struc(n)%tran3%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_TRAN4", &
                  "Leaf transmittance (NIR, dead)","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%tran4 )

             allocate(SiB2_struc(n)%tran4%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%tran4%vlevels)) 
             SiB2_struc(n)%tran4%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_REF1", &
                  "Leaf reflectance (VIS, live)","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%ref1 )

             allocate(SiB2_struc(n)%ref1%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%ref1%vlevels)) 
             SiB2_struc(n)%ref1%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_REF2", &
                  "Leaf reflectance (VIS, dead)","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%ref2 )

             allocate(SiB2_struc(n)%ref2%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%ref2%vlevels)) 
             SiB2_struc(n)%ref2%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_REF3", &
                  "Leaf reflectance (NIR, live)","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%ref3 )

             allocate(SiB2_struc(n)%ref3%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%ref3%vlevels)) 
             SiB2_struc(n)%ref3%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_REF4", &
                  "Leaf reflectance (NIR, dead)","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%ref4 )

             allocate(SiB2_struc(n)%ref4%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%ref4%vlevels)) 
             SiB2_struc(n)%ref4%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_VMAX0", &
                  "Max RuBisCO capacity at canopy top","mol m-2 s-1", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%vmax0 )

             allocate(SiB2_struc(n)%vmax0%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%vmax0%vlevels)) 
             SiB2_struc(n)%vmax0%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_EFFCON", &
                  "Quantum yield of photosynthesis MolCO2","mol mol-1", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%effcon )

             allocate(SiB2_struc(n)%effcon%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%effcon%vlevels)) 
             SiB2_struc(n)%effcon%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_GRADM", &
                  "Stomatal conductance slope parameter","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%gradm )

             allocate(SiB2_struc(n)%gradm%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%gradm%vlevels)) 
             SiB2_struc(n)%gradm%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_BINTER", &
                  "Minimum stomtal conductance","mol m-2 s-1", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%binter )

             allocate(SiB2_struc(n)%binter%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%binter%vlevels)) 
             SiB2_struc(n)%binter%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_ATHETA", &
                  "Light and RuBisCO (wc,we) coupling parameter","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%atheta )

             allocate(SiB2_struc(n)%atheta%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%atheta%vlevels)) 
             SiB2_struc(n)%atheta%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_BTHETA", &
                  "Light, RuBisCO, CHO sink( wc,we,ws) coupling parameter","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%btheta )

             allocate(SiB2_struc(n)%btheta%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%btheta%vlevels)) 
             SiB2_struc(n)%btheta%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_TRDA", &
                  "Temperature inhibition coef. for leaf respiration","K-1", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%trda )

             allocate(SiB2_struc(n)%trda%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%trda%vlevels)) 
             SiB2_struc(n)%trda%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_TRDM", &
                  "1/2 critical temperature inhib. for leaf respiration","K", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%trdm )

             allocate(SiB2_struc(n)%trdm%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%trdm%vlevels)) 
             SiB2_struc(n)%trdm%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_TROP", &
                  "Q10 temperature coeff. for physiology (298.16K)","K", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%trop )

             allocate(SiB2_struc(n)%trop%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%trop%vlevels)) 
             SiB2_struc(n)%trop%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_RESPCP", &
                  "Leaf respiration fraction for Vmax","mol m-2 s-1", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%respcp )

             allocate(SiB2_struc(n)%respcp%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%respcp%vlevels)) 
             SiB2_struc(n)%respcp%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_SLTI", &
                  "Slope of low temperature inhibition function","K-1", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%slti )

             allocate(SiB2_struc(n)%slti%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%slti%vlevels)) 
             SiB2_struc(n)%slti%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_SHTI", &
                  "Slope of high temperature inhibition function","K-1", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%shti )

             allocate(SiB2_struc(n)%shti%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%shti%vlevels)) 
             SiB2_struc(n)%shti%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_HLTII", &
                  "1/2 point of low temperature inhibition function","K", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%hltii )

             allocate(SiB2_struc(n)%hltii%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%hltii%vlevels)) 
             SiB2_struc(n)%hltii%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_hhti", &
                  "1/2 point of high temperature inhibition function","K", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%hhti )

             allocate(SiB2_struc(n)%hhti%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%hhti%vlevels)) 
             SiB2_struc(n)%hhti%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_SOREF1", &
                  "Soil reflectance (VIS; -)","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%soref1 )

             allocate(SiB2_struc(n)%soref1%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%soref1%vlevels)) 
             SiB2_struc(n)%soref1%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_SOREF2", &
                  "Soil reflectance (NIR; -)","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%soref2 )

             allocate(SiB2_struc(n)%soref2%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%soref2%vlevels)) 
             SiB2_struc(n)%soref2%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_BEE", &
                  "Soil wetness exponent","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%bee )

             allocate(SiB2_struc(n)%bee%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%bee%vlevels)) 
             SiB2_struc(n)%bee%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_PHSAT", &
                  "Soil tension at saturation","m", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%phsat )

             allocate(SiB2_struc(n)%phsat%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%phsat%vlevels)) 
             SiB2_struc(n)%phsat%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_SATCO", &
                  "Hydraulic conductivity at saturation","m s-1", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%satco )

             allocate(SiB2_struc(n)%satco%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%satco%vlevels)) 
             SiB2_struc(n)%satco%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_POROS", &
                  "Soil porosity","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%poros )

             allocate(SiB2_struc(n)%poros%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%poros%vlevels)) 
             SiB2_struc(n)%poros%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_SLOPE", &
                  "Cosine of mean slope","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%islope )

             allocate(SiB2_struc(n)%islope%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%islope%vlevels)) 
             SiB2_struc(n)%islope%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_WOPT", &
                  "Optimum soil water perc. for respiration","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%wopt )

             allocate(SiB2_struc(n)%wopt%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%wopt%vlevels)) 
             SiB2_struc(n)%wopt%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_WSAT", &
                  "Respiration rate at soil water saturation coeff","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%wsat )

             allocate(SiB2_struc(n)%wsat%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%wsat%vlevels)) 
             SiB2_struc(n)%wsat%value = LDT_rc%udef

             call populate_param_attribs( "SIB2_ZM", &
                  "Skewness exponent of respiration vs. soil water","-", &
                  SiB2_struc(n)%sib2, &
                  SiB2_struc(n)%zm )

             allocate(SiB2_struc(n)%zm%value(&
                  LDT_rc%lnc(n),LDT_rc%lnr(n),&
                  SiB2_struc(n)%zm%vlevels))
             SiB2_struc(n)%zm%value = LDT_rc%udef

          end if
       end do

       !- Read in ldt.config file entries:
       if( SiB2_struc(1)%sib2%selectOpt == 1 ) then 

          !- Read in SiB2 soil parameters from a-priori gridded maps:

          call ESMF_ConfigFindLabel(LDT_config,"SiB2 static parameter directory:",rc=rc)
          do n=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,SiB2_struc(n)%sib2parmsdir,rc=rc)
             call LDT_verify(rc,'SiB2 static parameter directory: not specified')
          enddo

          ! --  Grid info and inputs:
          call ESMF_ConfigFindLabel(LDT_config,"SiB2 parameter spatial transform:",&
               rc=rc)
          do n=1,LDT_rc%nnest
             call ESMF_ConfigGetAttribute(LDT_config,sib2parms_gridtransform(n),&
                  rc=rc)
             call LDT_verify(rc,'SiB2 parameter spatial transform: option not specified in the config file')
          enddo

          sib2%filltype = "none"
          call ESMF_ConfigGetAttribute(LDT_config, sib2%filltype, &
               label="SiB2 parameter fill option:",rc=rc)
          call LDT_verify(rc,"SiB2 parameter fill option: option not specified in the config file")

          if( sib2%filltype == "neighbor" .or. sib2%filltype == "average" ) then
             call ESMF_ConfigGetAttribute(LDT_config, sib2%fillvalue, &
                  label="SiB2 parameter fill value:",rc=rc)
             call LDT_verify(rc,"SiB2 parameter fill value: option not specified in the config file")

             call ESMF_ConfigGetAttribute(LDT_config, sib2%fillradius, &
                  label="SiB2 parameter fill radius:",rc=rc)
             call LDT_verify(rc,"SiB2 parameter fill radius: option not specified in the config file")
          else
             write(LDT_logunit,*) " -- 'NONE' Parameter-Mask Agreement Option Selected for SiB2 parameters"
          end if

          call ESMF_ConfigGetAttribute(LDT_config,sib2parms_proj,&
               label="SiB2 map projection:",rc=rc)
          call LDT_verify(rc,'SiB2 map projection: option not specified in the config file')

          call LDT_readDomainConfigSpecs("SiB2", sib2parms_proj, sib2parms_gridDesc)

       end if   ! End SiB2 parameter array definitions

       !- Loop over all domain nests:
       do n = 1, LDT_rc%nnest
          if( SiB2_struc(n)%sib2%selectOpt == 1 ) then 

             !       call LDT_gridOptChecks( n, "SiB2", LDT_rc%sib2parms_gridtransform, &
             !                       LDT_rc%sib2parms_proj, LDT_rc%sib2parms_gridDesc(n,9) )

             !       call LDT_sib2OptChecks( "SiB2", LDT_rc%sib2parms_proj, LDT_rc%sib2parms_gridtransform )

             ! -------

             write(LDT_logunit,*) "Reading in SiB2 Parms Directory: "//&
                  trim(SiB2_struc(n)%sib2parmsdir)
             call read_sib2_parms(n, SiB2_struc(n)%sib2%num_bins,       &
                  SiB2_struc(n)%z2%value(:,:,:),     &  ! 1
                  SiB2_struc(n)%z1%value(:,:,:),     &  ! 2
                  SiB2_struc(n)%vcover%value(:,:,:), &  ! 3
                  SiB2_struc(n)%chil%value(:,:,:),   &  ! 4
                  SiB2_struc(n)%sodep%value(:,:,:),  &  ! 5
                  SiB2_struc(n)%rootd%value(:,:,:),  &  ! 6
                  SiB2_struc(n)%phc%value(:,:,:),    &  ! 7
                  SiB2_struc(n)%tran1%value(:,:,:),  &  ! 8
                  SiB2_struc(n)%tran2%value(:,:,:),  &  ! 9
                  SiB2_struc(n)%tran3%value(:,:,:),  &  ! 10
                  SiB2_struc(n)%tran4%value(:,:,:),  &  ! 11
                  SiB2_struc(n)%ref1%value(:,:,:),   &  ! 12
                  SiB2_struc(n)%ref2%value(:,:,:),   &  ! 13
                  SiB2_struc(n)%ref3%value(:,:,:),   &  ! 14
                  SiB2_struc(n)%ref4%value(:,:,:),   &  ! 15
                  SiB2_struc(n)%vmax0%value(:,:,:),  &  ! 16
                  SiB2_struc(n)%effcon%value(:,:,:), &  ! 17
                  SiB2_struc(n)%gradm%value(:,:,:),  &  ! 18
                  SiB2_struc(n)%binter%value(:,:,:), &  ! 19
                  SiB2_struc(n)%atheta%value(:,:,:), &  ! 20
                  SiB2_struc(n)%btheta%value(:,:,:), &  ! 21
                  SiB2_struc(n)%trda%value(:,:,:),   &  ! 22
                  SiB2_struc(n)%trdm%value(:,:,:),   &  ! 23
                  SiB2_struc(n)%trop%value(:,:,:),   &  ! 24
                  SiB2_struc(n)%respcp%value(:,:,:), &  ! 25
                  SiB2_struc(n)%slti%value(:,:,:),   &  ! 26
                  SiB2_struc(n)%shti%value(:,:,:),   &  ! 27
                  SiB2_struc(n)%hltii%value(:,:,:),  &  ! 28
                  SiB2_struc(n)%hhti%value(:,:,:),   &  ! 29
                  SiB2_struc(n)%soref1%value(:,:,:), &  ! 30
                  SiB2_struc(n)%soref2%value(:,:,:), &  ! 31
                  SiB2_struc(n)%bee%value(:,:,:),    &  ! 32
                  SiB2_struc(n)%phsat%value(:,:,:),  &  ! 33
                  SiB2_struc(n)%satco%value(:,:,:),  &  ! 34 
                  SiB2_struc(n)%poros%value(:,:,:),  &  ! 35
                  SiB2_struc(n)%islope%value(:,:,:), &  ! 36
                  SiB2_struc(n)%wopt%value(:,:,:),   &  ! 37
                  SiB2_struc(n)%wsat%value(:,:,:),   &  ! 38
                  SiB2_struc(n)%zm%value(:,:,:) )       ! 39

             write(LDT_logunit,*) "Done reading from: "//trim(SiB2_struc(n)%sib2parmsdir)

          end if
       enddo

    endif  ! End SiB2 model check

  end subroutine SiB2Parms_init

  subroutine SiB2Parms_writeHeader(n,ftn,dimID)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
   integer      :: n
   integer      :: ftn
   integer      :: dimID(3)
   integer      :: tdimID(3)

   tdimID(1) = dimID(1)
   tdimID(2) = dimID(2)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

   if( SiB2_struc(1)%sib2%selectOpt == 1 ) then
      
      SiB2_struc(n)%sib2%vlevels = SiB2_struc(n)%sib2%num_bins
      
      call LDT_verify(nf90_def_dim(ftn,'sib2types',&
           SiB2_struc(n)%sib2%vlevels,tdimID(3)))
      
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%z2)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%z1)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%vcover)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%chil)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%sodep)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%rootd)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%phc)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%tran1)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%tran2)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%tran3)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%tran4)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%ref1)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%ref2)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%ref3)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%ref4)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%vmax0)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%effcon)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%gradm)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%binter)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%atheta)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%btheta)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%trda)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%trdm)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%trop)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%respcp)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%slti)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%shti)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%hltii)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%hhti)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%soref1)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%soref2)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%bee)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%phsat)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%satco)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%poros)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%islope)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%wopt)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%wsat)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           SiB2_struc(n)%zm)
      
   endif
#endif
  end subroutine SiB2Parms_writeHeader

  subroutine SiB2Parms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    if( SiB2_struc(1)%sib2%selectOpt == 1 ) then

       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%z2)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%z1)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%vcover)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%chil)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%sodep)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%rootd)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%phc)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%tran1)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%tran2)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%tran3)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%tran4)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%ref1)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%ref2)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%ref3)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%ref4)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%vmax0)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%effcon)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%gradm)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%binter)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%atheta)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%btheta)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%trda)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%trdm)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%trop)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%respcp)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%slti)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%shti)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%hltii)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%hhti)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%soref1)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%soref2)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%bee)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%phsat)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%satco)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%poros)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%islope)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%wopt)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%wsat)
       call LDT_writeNETCDFdata(n,ftn,SiB2_struc(n)%zm)

    endif

  end subroutine SiB2Parms_writeData

!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name)

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name

! ____________________________________________________
    
   
   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = 1
   paramEntry%selectOpt = 1
   paramEntry%source = "SiB2"
   paramEntry%units ="none"
   paramEntry%num_times = 1
   paramEntry%num_bins = 12
   paramEntry%standard_name = trim(short_name)

  end subroutine set_param_attribs

 end module SiB2_parmsMod

