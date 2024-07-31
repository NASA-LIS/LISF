!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module hyssibalb_module 
!BOP
! !MODULE: hyssibalb_module
!
! !DESCRIPTION:
!  In order to use the HySSiB albedo calculations, UMD vegetation types
!  need to be mapped to SSiB vegetation types.  The code in this file provides
!  a description of the data structure used to contain the vegetation specific 
!  constants used in HySSiB's albedo calculations.  The varaiables specified in
!  the data structure include:
!
!  \begin{description}
!   \item[cedfu]
!   \item[cedir] 
!   \item[cedfu1] 
!   \item[cedir1]
!   \item[cedfu2]
!   \item[xmiu]
!   \item[cedir2]
!   \item[clefu]
!   \item[cledir]
!   \item[cether]
!   \item[xmiw]
!  \end{description}
!
! !REVISION HISTORY:
!  01 Dec 2007: Chuck Alonge; Initial Specification
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: hyssibalb_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: hyssibalb 
!EOP
  
  type, public ::  hyssibalbdec 
         
    !=== Constants used in albedo calculations: =========
    real, dimension(13,12,9)       :: cedfu
    real, dimension(13,12,9,3)     :: cedir
    real, dimension(2,13,12,6,3)   :: cedfu1
    real, dimension(2,13,12,6,3,3) :: cedir1
    real, dimension(2,13,12,6,3)   :: cedfu2
    real, dimension(12,3)          :: xmiu
    real, dimension(2,13,12,6,3,3) :: cedir2
    real, dimension(13,12,9)       :: cledfu
    real, dimension(13,12,9,3)     :: cledir
    real, dimension(13,12,2)       :: cether
    real, dimension(12,3)          :: xmiw
     
  end type hyssibalbdec !hyssibalbdec  

  type (hyssibalbdec), allocatable :: hyssibalb(:)
  save

  contains

!BOP
!
! !ROUTINE: hyssibalb_ini
! \label{hyssibalb_ini}
!
! !INTERFACE:
  subroutine hyssibalb_ini()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use hyssib_lsmMod ! HY-SSiB tile variables
! !DESCRIPTION:
!
!  Reads in runtime HY-SSiB albedo parameters to populate constants used
!  in Hyssib albedo calculations
!
!EOP
    implicit none 
    integer :: n

    allocate(hyssibalb(LIS_rc%nnest))
    do n=1,LIS_rc%nnest

      open(unit=12,file=hyssib_struc(n)%afile,status='old', &
                  form='unformatted')
      read(12) hyssibalb(n)%cedfu, hyssibalb(n)%cedir, hyssibalb(n)%cedfu1, &
               hyssibalb(n)%cedir1, hyssibalb(n)%cedfu2, hyssibalb(n)%cedir2, &
               hyssibalb(n)%cledfu, hyssibalb(n)%cledir, hyssibalb(n)%xmiu, &
               hyssibalb(n)%cether, hyssibalb(n)%xmiw
      close(12)

    enddo

  end subroutine hyssibalb_ini

end module hyssibalb_module

