!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: get_cdf_params
! \label{get_cdf_params}
!
! !REVISION HISTORY:
! 01 May 2023: Yeosang Yoon, initial code
!
! !INTERFACE:
subroutine get_cdf_params (n, fname, month, param_a, param_b, mean, std)

  !USES:
  use LIS_coreMod
  use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  !ARGUMENTS:
  integer, intent(in)           :: n
  character(len=*), intent(in)  :: fname
  integer,          intent(in)  :: month            ! month of year (1-12)

  !EOP
  logical                        :: file_exists
  integer                        :: ftn
  integer                        :: param_a_id, param_b_id, mean_id, &
       std_id, ngrid_id, nlead_id
  integer                        :: ngrid, nlead
  integer                        :: r, c, nc, c1, r1, gid, l
  real, allocatable              :: param_hires(:,:)
  real, allocatable              :: param_hires_2d(:,:)
  real                           :: param_a(LIS_rc%ngrid(n),8) !8:lead time
  real                           :: param_b(LIS_rc%ngrid(n),8)
  real                           :: mean(LIS_rc%ngrid(n),8)
  real                           :: std(LIS_rc%ngrid(n),8)

  ! check file
  inquire(file=fname, exist=file_exists)
  if (.not. file_exists) then
     write(LIS_logunit,*) '[ERR] ', trim(fname) // ' does not exist'
     call LIS_endrun()
  endif

  write(LIS_logunit,*) &
       '[INFO] Getting MOGREPS-G bias correction parameters ', trim(fname)
  call LIS_verify(nf90_open(path=trim(fname), mode=NF90_NOWRITE, &
       ncid=ftn), 'nf90_open failed in get_cdf_params')

  call LIS_verify(nf90_inq_dimid(ftn, "ngrid", ngrid_id), &
       'nf90_inq_dimid failed for ngrid in get_cdf_params')
  call LIS_verify(nf90_inquire_dimension(ftn, ngrid_id, len=ngrid),&
       'nf90_inquire_dimension failed for ngrid in get_cdf_params')

  if (LIS_rc%gnc(n)*LIS_rc%gnr(n) .ne. ngrid) then
     write(LIS_logunit,*) '[ERR] The input dimensions of the '//trim(fname)
     write(LIS_logunit,*) '(',ngrid,')'
     write(LIS_logunit,*) &
          'does not match the dimensions in the LIS parameter file'
     write(LIS_logunit,*) '(',LIS_rc%gnc(n)*LIS_rc%gnr(n),')'
     call LIS_endrun()
  endif

  call LIS_verify(nf90_inq_dimid(ftn, "lead_time", nlead_id), &
       'nf90_inq_dimid failed for nlead in get_cdf_params')
  call LIS_verify(nf90_inquire_dimension(ftn, nlead_id, len=nlead),&
       'nf90_inquire_dimension failed for nlead in get_cdf_params')

  allocate(param_hires(LIS_rc%gnc(n)*LIS_rc%gnr(n),nlead))
  allocate(param_hires_2d(LIS_rc%gnc(n),LIS_rc%gnr(n)))

  param_hires = -9999.0
  param_hires_2d = -9999.0

  ! read param_a
  call LIS_verify(nf90_inq_varid(ftn, 'cdf_param_a', param_a_id), &
       'nf90_inq_varid failed for cdf_param_a in get_cdf_params')
  call LIS_verify(nf90_get_var(ftn, param_a_id, param_hires, &
       start=(/1,1,month/), count=(/ngrid,nlead,1/)), &
       'nf90_get_var failed for cdf_param_a in get_cdf_params')

  do l = 1, nlead
     ! 1D -> 2D
     do r = 1, LIS_rc%gnr(n)
        do c = 1, LIS_rc%gnc(n)
           param_hires_2d(c,r) = param_hires(c+(r-1)*LIS_rc%gnc(n),l)
        enddo
     enddo

     !subsets the data for each processor's domain
     nc = (LIS_ewe_halo_ind(n,LIS_localPet+1)-LIS_ews_halo_ind(n,LIS_localPet+1))+1
     do r = LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c = LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           c1 = c - LIS_ews_halo_ind(n,LIS_localPet+1)+1
           r1 = r - LIS_nss_halo_ind(n,LIS_localPet+1)+1
           gid = LIS_domain(n)%gindex(c1,r1)
           if (gid .ne. -1) then
              param_a(gid,l) = param_hires_2d(c,r)
           endif
        enddo
     enddo
  enddo

  ! read param_b
  call LIS_verify(nf90_inq_varid(ftn, 'cdf_param_b', param_b_id), &
       'nf90_inq_varid failed for cdf_param_b in get_cdf_params')
  call LIS_verify(nf90_get_var(ftn, param_b_id, param_hires, &
       start=(/1,1,month/), count=(/ngrid,nlead,1/)),&
       'nf90_get_var failed for cdf_param_b in get_cdf_params')

  do l = 1, nlead
     ! 1D -> 2D
     do r = 1, LIS_rc%gnr(n)
        do c = 1,LIS_rc%gnc(n)
           param_hires_2d(c,r) = param_hires(c+(r-1)*LIS_rc%gnc(n),l)
        enddo
     enddo

     !subsets the data for each processor's domain
     nc = (LIS_ewe_halo_ind(n,LIS_localPet+1)-LIS_ews_halo_ind(n,LIS_localPet+1))+1
     do r = LIS_nss_halo_ind(n,LIS_localPet+1), LIS_nse_halo_ind(n,LIS_localPet+1)
        do c = LIS_ews_halo_ind(n,LIS_localPet+1), LIS_ewe_halo_ind(n,LIS_localPet+1)
           c1 = c - LIS_ews_halo_ind(n,LIS_localPet+1) + 1
           r1 = r - LIS_nss_halo_ind(n,LIS_localPet+1) + 1
           gid = LIS_domain(n)%gindex(c1,r1)
           if (gid.ne.-1) then
              param_b(gid,l) = param_hires_2d(c,r)
           endif
        enddo
     enddo
  enddo

  ! read mean
  call LIS_verify(nf90_inq_varid(ftn, 'mean', mean_id), &
       'nf90_inq_varid failed for mean in get_cdf_params')
  call LIS_verify(nf90_get_var(ftn, mean_id, param_hires, &
       start=(/1,1,month/), count=(/ngrid,nlead,1/)),&
       'nf90_get_var failed for mean in get_cdf_params')

  do l = 1, nlead
     ! 1D -> 2D
     do r = 1, LIS_rc%gnr(n)
        do c = 1, LIS_rc%gnc(n)
           param_hires_2d(c,r) = param_hires(c+(r-1)*LIS_rc%gnc(n),l)
        enddo
     enddo

     !subsets the data for each processor's domain
     nc = (LIS_ewe_halo_ind(n,LIS_localPet+1)-LIS_ews_halo_ind(n,LIS_localPet+1))+1
     do r = LIS_nss_halo_ind(n,LIS_localPet+1), LIS_nse_halo_ind(n,LIS_localPet+1)
        do c = LIS_ews_halo_ind(n,LIS_localPet+1), LIS_ewe_halo_ind(n,LIS_localPet+1)
           c1 = c - LIS_ews_halo_ind(n,LIS_localPet+1)+1
           r1 = r - LIS_nss_halo_ind(n,LIS_localPet+1)+1
           gid = LIS_domain(n)%gindex(c1,r1)
           if (gid.ne.-1) then
              mean(gid,l) = param_hires_2d(c,r)
           endif
        enddo
     enddo
  enddo

  ! read std
  call LIS_verify(nf90_inq_varid(ftn, 'std', std_id), &
       'nf90_inq_varid failed for std in get_cdf_params')
  call LIS_verify(nf90_get_var(ftn, std_id, param_hires, &
       start=(/1,1,month/), count=(/ngrid,nlead,1/)),&
       'nf90_get_var failed for std in get_cdf_params')

  do l = 1, nlead
     ! 1D -> 2D
     do r = 1, LIS_rc%gnr(n)
        do c = 1, LIS_rc%gnc(n)
           param_hires_2d(c,r) =param_hires(c+(r-1)*LIS_rc%gnc(n),l)
        enddo
     enddo

     !subsets the data for each processor's domain
     nc = (LIS_ewe_halo_ind(n,LIS_localPet+1)-LIS_ews_halo_ind(n,LIS_localPet+1))+1
     do r = LIS_nss_halo_ind(n,LIS_localPet+1), LIS_nse_halo_ind(n,LIS_localPet+1)
        do c = LIS_ews_halo_ind(n,LIS_localPet+1), LIS_ewe_halo_ind(n,LIS_localPet+1)
           c1 = c - LIS_ews_halo_ind(n,LIS_localPet+1) + 1
           r1 = r - LIS_nss_halo_ind(n,LIS_localPet+1) + 1
           gid = LIS_domain(n)%gindex(c1,r1)
           if (gid.ne.-1) then
              std(gid,l) = param_hires_2d(c,r)
           endif
        enddo
     enddo
  enddo

  deallocate(param_hires)
  deallocate(param_hires_2d)

  call LIS_verify(nf90_close(ftn),'failed to close in get_cdf_params')

  write(LIS_logunit,*) &
       '[INFO] Done reading MOGREPS-G bias correction parameters data '

end subroutine get_cdf_params

