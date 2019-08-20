!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_orig_nldas1_elevdiff
! \label{read_orig_nldas1_elevdiff}
!
! !REVISION HISTORY: 
!  7 Nov 2012: Sujay Kumar - initial specification based on the NLDAS1
!                            implementation
!
! !INTERFACE:
subroutine read_orig_nldas1_elevdiff(n)
! !USES:
  use LIS_coreMod
  use nldas1_forcingMod, only : nldas1_struc
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
   real    :: elevdiff(nldas1_struc(n)%ncold,nldas1_struc(n)%nrold)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

   inquire(file=LIS_rc%paramfile(n), exist=file_exists)
   if(file_exists) then 

      write(LIS_logunit,*) "Reading NLDAS1 original elevation difference file ..."

      call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
           mode=NF90_NOWRITE,ncid=nid),&
           'nf90_open failed in read_orig_nldas1_elevdiff')
      call LIS_verify(nf90_inq_varid(nid,"ELEVDIFF_NLDAS1",elevId),&
           'nf90_inq_varid failed in read_orig_nldas1_elevdiff')
      call LIS_verify(nf90_get_var(nid,elevId,elevdiff),&
           'nf90_get_var failed in read_orig_nldas1_elevdiff')      
      call LIS_verify(nf90_close(nid))

      do r=1,nldas1_struc(n)%nrold
         do c=1,nldas1_struc(n)%ncold            
            nldas1_struc(n)%orig_ediff(c+(r-1)*nldas1_struc(n)%ncold) = & 
                 elevdiff(c,r)
         enddo
      enddo
   endif

   write(LIS_logunit,*) & 
        "Finish reading original NLDAS1 elevation difference data"
#endif
 end subroutine read_orig_nldas1_elevdiff
