!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
! NASA GSFC Land Information System (LIS)
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_mpiMod
!BOP
!
! !MODULE: LDT_mpiMod
!
! !REVISION HISTORY: 
! 02 Jan 2002 : Sujay Kumar: Initial Version
! 
! !USES:
#if (defined SPMD)
#if (defined USE_INCLUDE_MPI)
  include 'mpif.h' 
#else
  use mpi
#endif
#endif
! !DESCRIPTION: 
!
! Data and parameters used for MPI. Some shorthand variables 
! with shorter names than the standard MPI parameters. 
! Also some variables used for heap management. Adopted from CLM
!
!EOP

end module LDT_mpiMod
