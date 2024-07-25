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
!BOP
! !ROUTINE: read_orig_nldas20_elevdiff
! \label{read_orig_nldas20_elevdiff}
!
! !REVISION HISTORY:
! 11 Jul 2024: David Mocko, Initial Specification
!                           (derived from read_orig_nldas2_elevdiff.F90)
!
! !INTERFACE:
subroutine read_orig_nldas20_elevdiff(n)
! !USES:
  use LIS_coreMod
  use nldas20_forcingMod, only : nldas20_struc
  use LIS_logmod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
!
! !DESCRIPTION:
!     Open and read in original NLDAS elevation difference file.
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \end{description}
!
!EOP
  integer :: i, err
  logical :: file_exists
  integer :: nid,elevId
  integer :: c,r
  real    :: elevdiff(nldas20_struc(n)%ncold,nldas20_struc(n)%nrold)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=LIS_rc%paramfile(n),exist=file_exists)

  if (file_exists) then
     write(LIS_logunit,*)                                          &
          "[INFO] Reading NLDAS-2 original elevation difference file ..."

     call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),           &
          mode=NF90_NOWRITE,ncid=nid),                             &
          "nf90_open failed in read_orig_nldas20_elevdiff")
     call LIS_verify(nf90_inq_varid(nid,"ELEVDIFF_NLDAS2",elevId), &
          "nf90_inq_varid failed in read_orig_nldas20_elevdiff")
     call LIS_verify(nf90_get_var(nid,elevId,elevdiff),            &
          "nf90_get_var failed in read_orig_nldas20_elevdiff")
     call LIS_verify(nf90_close(nid))

     do r = 1,nldas20_struc(n)%nrold
        do c = 1,nldas20_struc(n)%ncold
           nldas20_struc(n)%orig_ediff(c+(r-1)*nldas20_struc(n)%ncold) = elevdiff(c,r)
        enddo
     enddo
  endif

  write(LIS_logunit,*) "[INFO] Finished reading original NLDAS-2"
  write(LIS_logunit,*) "[INFO] elevation difference data."
#endif

end subroutine read_orig_nldas20_elevdiff

