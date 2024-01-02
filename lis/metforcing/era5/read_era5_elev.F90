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
! !ROUTINE: read_era5_elev
!  \label{read_era5_elev}
!
! !REVISION HISTORY:
!
!  23 dec 2019; Sujay Kumar; Initial Specificaton
!
! !INTERFACE:
subroutine read_era5_elev(n,findex)
! !USES:
  use LIS_coreMod
  use LIS_metforcingMod
  use LIS_logMod
  use era5_forcingMod
  use LIS_fileIOMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
! !DESCRIPTION:
!
!  Opens, reads, and interpolates ERA5 model elevation to the LIS
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
! 
!  The routines invoked are: 
!   \begin{description}
!    \item[ij\_to\_latlon](\ref{ij_to_latlon}) \newline
!     computes the lat lon values in LIS grid projection
!   \end{description}
!EOP
   integer :: i, err
   logical :: file_exists
   integer :: nid,elevId
   integer :: c,r
   real    :: elev(LIS_rc%gnc(n),LIS_rc%gnr(n))
   real    :: elev_subset(LIS_rc%lnc(n),LIS_rc%lnr(n))
  
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

   inquire(file=LIS_rc%paramfile(n), exist=file_exists)
   if(file_exists) then 

      write(LIS_logunit,*) " Reading ERA5 elevation data ... "

      call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
           mode=NF90_NOWRITE,ncid=nid),&
           'nf90_open failed in read_era5_elev')
      call LIS_verify(nf90_inq_varid(nid,"ELEV_ERA5",elevId),&
           'nf90_inq_varid failed in read_era5_elev')
      call LIS_verify(nf90_get_var(nid,elevId,elev),&
           'nf90_get_var failed in read_era5_elev')      
      call LIS_verify(nf90_close(nid))

      elev_subset(:,:) = elev(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1))

      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            if(LIS_domain(n)%gindex(c,r).ne.-1) then 
               LIS_forc(n,findex)%modelelev(LIS_domain(n)%gindex(c,r)) =&
                    elev_subset(c,r)
            endif
         enddo
      enddo
   endif

   write(LIS_logunit,*) & 
        "Finish reading original ERA5 elevation data"
#endif

end subroutine read_era5_elev
 
