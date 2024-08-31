!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module hyssibveg_module 
!BOP
! !MODULE: hyssibveg_module
!
! !DESCRIPTION:
!  In order to use the HySSiB vegetation table data efficiently, 
!  SiB vegetation parameter data (static and monthly) are read 
!  into this module to reduce model i/o when updating monthly values.
!  The varaiables specified in this data structure include:
!
!  \begin{description}
!   \item[rstpar]
!    par influence on stomatal resist. coefficients (rst)
!   \item[chil]
!    leaf angle distribution factor 
!   \item[topt] 
!    optimum temperature for rst calculation
!   \item[tll]
!    bottom temperature for rst calculation
!   \item[tu]
!    top temperature for rst calculation
!   \item[defac]
!    dew factor for rst calculation
!   \item[ph1]
!    stome slope factor
!   \item[ph2]
!    point at which stomates close
!   \item[rootd]
!    rooting depth for canopy/ground cover
!   \item[bee]
!    Clapp-Hornberger empirical constant
!   \item[phsat]
!    soil moisture potential at saturation
!   \item[satco]
!    saturation hydraulic conductivity
!   \item[poros]
!    soil porosity
!   \item[zdepth]
!    exact independent depth of 3 soil layers
!   \item[slope]
!    avg. topographic slope in %
!   \item[green]
!    green leaf fraction
!   \item[vcover]
!    fraction of vegetation cover [fpar/greenness]
!   \item[zlt]
!    leaf area index
!   \item[z0]
!    roughness height
!   \item[z2]
!    canopy top height
!   \item[z1]
!    canopy base height
!   \item[rdc]
!    ground to canopy air space resistance coeff.
!   \item[rbc]
!    bulk canopy boundary layer resistance coeff.
!  \end{description}
!
! !REVISION HISTORY:
!  01 Dec 2007: Chuck Alonge; Initial Specification
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: hyssibveg_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: hyssibveg
!EOP
  
  type, public ::  hyssibvegdec 
         
    !=== Vegetation Constants used in calculations: =========
    real, dimension(13,2,3)  :: rstpar
    real, dimension(13,2)    :: chil
    real, dimension(13,2)    :: topt
    real, dimension(13,2)    :: tll
    real, dimension(13,2)    :: tu
    real, dimension(13,2)    :: defac
    real, dimension(13,2)    :: ph1
    real, dimension(13,2)    :: ph2
    real, dimension(13,2)    :: rootd
    real, dimension(13)      :: bee
    real, dimension(13)      :: phsat
    real, dimension(13)      :: satco
    real, dimension(13)      :: poros
    real, dimension(13,3)    :: zdepth 
    real, dimension(13)      :: slope

    !=== Monthly Constants used in calculations: =========
    real, dimension(13,12,2) :: green 
    real, dimension(13,12,2) :: vcover 
    real, dimension(13,12,2) :: zlt 
    real, dimension(13,12)   :: z0 
    real, dimension(13,12)   :: dd 
    real, dimension(13,12)   :: z2 
    real, dimension(13,12)   :: z1 
    real, dimension(13,12)   :: rdc
    real, dimension(13,12)   :: rbc
     
  end type hyssibvegdec !hyssibvegdec
  
  type(hyssibvegdec), allocatable  ::  hyssibveg(:)

  save

  contains

!BOP
!
! !ROUTINE: hyssibveg_ini
! \label{hyssibveg_ini}
!
! !INTERFACE:
  subroutine hyssibveg_ini()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use hyssib_lsmMod ! HY-SSiB model control variables
! !DESCRIPTION:
!
!  Reads in runtime HY-SSiB vegetation parameters to populate 
!  constant and monthly varying values Hyssib calculations
!
!EOP
    implicit none 
    integer :: n

    allocate(hyssibveg(LIS_rc%nnest))

    do n=1,LIS_rc%nnest

      open(unit=11,file=hyssib_struc(n)%vfile,status='old', &
                  form='unformatted')
      read(11) hyssibveg(n)%rstpar, hyssibveg(n)%chil, hyssibveg(n)%topt, &
               hyssibveg(n)%tll, hyssibveg(n)%tu, hyssibveg(n)%defac, &
               hyssibveg(n)%ph1, hyssibveg(n)%ph2, hyssibveg(n)%rootd, &
               hyssibveg(n)%bee, hyssibveg(n)%phsat, hyssibveg(n)%satco, &
               hyssibveg(n)%poros, hyssibveg(n)%zdepth, hyssibveg(n)%slope

      read(11) hyssibveg(n)%green, hyssibveg(n)%vcover, hyssibveg(n)%zlt, &
               hyssibveg(n)%z0, hyssibveg(n)%dd, hyssibveg(n)%z2, &
               hyssibveg(n)%z1, hyssibveg(n)%rdc, hyssibveg(n)%rbc

      close(11)

    enddo

  end subroutine hyssibveg_ini

end module hyssibveg_module


