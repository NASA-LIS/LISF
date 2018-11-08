!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module LISrunoffdataMod
!BOP
! 
! !MODULE: LISrunoffdataMod
! 
! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!  8 Jan 2016: Sujay Kumar, initial implementation
! 17 Mar 2016: Augusto Getirana, Save in memory input file name and surface runoff and baseflow variables - this will reduce the number of times input files are read
! 
! !USES: 
  use ESMF
  
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LISrunoffdata_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: LISrunoffdata_struc
  
  type, public :: LISrunoffdatadec
     
     real                :: outInterval 
     character*50        :: odir
 
     !ag - 17Mar2016
     character*100       :: previous_filename
     real, allocatable   :: qs(:,:),qsb(:,:),evap(:,:)
     
  end type LISrunoffdatadec

  type(LISrunoffdatadec), allocatable :: LISrunoffdata_struc(:)

contains
 
!BOP
!
! !ROUTINE: LISrunoffdata_init
! \label{LISrunoffdata_init}
! 
  subroutine LISrunoffdata_init
    !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod

    integer              :: n 
    integer              :: status
    character*10         :: time

    allocate(LISrunoffdata_struc(LIS_rc%nnest))
       
    call ESMF_ConfigFindLabel(LIS_config,&
         "LIS runoff output directory:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,LISrunoffdata_struc(n)%odir,rc=status)
       call LIS_verify(status,&
            "LIS runoff output directory: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "LIS runoff output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "LIS runoff output interval: not defined")

       call LIS_parseTimeString(time,LISrunoffdata_struc(n)%outInterval)
    enddo
    
    !ag - 17Mar2016
    do n=1, LIS_rc%nnest
      LISrunoffdata_struc(n)%previous_filename='none'
      allocate(LISrunoffdata_struc(n)%qs(LIS_rc%gnc(n),LIS_rc%gnr(n)))
      allocate(LISrunoffdata_struc(n)%qsb(LIS_rc%gnc(n),LIS_rc%gnr(n)))
      allocate(LISrunoffdata_struc(n)%evap(LIS_rc%gnc(n),LIS_rc%gnr(n)))
      LISrunoffdata_struc(n)%qs=LIS_rc%udef
      LISrunoffdata_struc(n)%qsb=LIS_rc%udef
      LISrunoffdata_struc(n)%evap=LIS_rc%udef
    enddo
    
  end subroutine LISrunoffdata_init
end module LISrunoffdataMod
