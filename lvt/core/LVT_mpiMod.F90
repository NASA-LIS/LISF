!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!
!
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: LVT_mpiMod
! \label(LVT_mpiMod)
!
! !INTERFACE:
module LVT_mpiMod
! 
! !USES:
#if (defined SPMD)
#if (defined USE_INCLUDE_MPI)
  include 'mpif.h' 
#else
  use mpi
#endif
#endif
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! Data and parameters used for MPI. Some shorthand variables 
! with shorter names than the standard MPI parameters. 
! Also some variables used for heap management. Adopted from CLM
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 02 Jan 2002 : Sujay Kumar: Initial Version
! 
!EOP

end module LVT_mpiMod
