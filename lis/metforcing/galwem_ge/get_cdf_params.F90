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
!
! !ROUTINE: get_cdf_params
! \label{get_cdf_params}
!
! !REVISION HISTORY:
! 01 May 2023: Yeosang Yoon, initial code
!
! !INTERFACE:
subroutine get_cdf_params (n, rootdir, month, model_cdf, ref_cdf, landmask)

  !USES:
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_coreMod
  use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  !ARGUMENTS:
  integer, intent(in)           :: n
  character(len=*), intent(in)  :: rootdir
  integer,          intent(in)  :: month            ! month of year (1-12)

  !EOP
  logical                        :: file_exists
  integer                        :: ftn
  integer                        :: ngrid_id, nlead_id, nperc_id
  integer                        :: m_cdf_id, r_cdf_id, mask_id
  integer                        :: ngrid, nlead, nperc
  integer                        :: r, c, c1, r1, gid, l
  real, allocatable              :: param_hires1(:,:,:)
  real, allocatable              :: param_hires2(:,:,:)
  real                           :: model_cdf(LIS_rc%ngrid(n),101,8) !npercentile=101, nlead=8
  real                           :: ref_cdf(LIS_rc%ngrid(n),101,8)
  real                           :: landmask(LIS_rc%ngrid(n))
  real, allocatable              :: param_tmp1(:)
  real, allocatable              :: param_tmp2(:,:)

  character(len=LIS_CONST_PATH_LEN) :: fname
  character(len=LIS_CONST_PATH_LEN) :: filename
  character(len=2)               :: fmn

  ! get file name
  fname = 'mogrepsg_cdf_params_'
  write (UNIT=fmn, FMT='(i2.2)') month
  filename = trim(rootdir) // '/' // trim(fname) // fmn // '.nc'

  ! check file
  inquire(file=filename, exist=file_exists)
  if (.not. file_exists) then
     write(LIS_logunit,*) '[ERR] ', trim(filename) // ' does not exist'
     call LIS_endrun()
  endif

  write(LIS_logunit,*) &
       '[INFO] Getting GALWEM-GE bias correction parameters ', &
       trim(filename)
  call LIS_verify(nf90_open(path=trim(filename), mode=NF90_NOWRITE, &
       ncid=ftn), 'nf90_open failed in get_cdf_params')

  call LIS_verify(nf90_inq_dimid(ftn, "ngrid", ngrid_id), &
       'nf90_inq_dimid failed for ngrid in get_cdf_params')
  call LIS_verify(nf90_inquire_dimension(ftn, ngrid_id, len=ngrid),&
       'nf90_inquire_dimension failed for ngrid in get_cdf_params')

  if (LIS_rc%gnc(n)*LIS_rc%gnr(n) .ne. ngrid) then
     write(LIS_logunit,*) &
          '[ERR] The input dimensions of the ' // trim(fname)
     write(LIS_logunit,*) '(',ngrid,')'
     write(LIS_logunit,*) &
          'does not match the dimensions in the LIS parameter file'
     write(LIS_logunit,*) '(', LIS_rc%gnc(n)*LIS_rc%gnr(n), ')'
     call LIS_endrun()
  endif

  call LIS_verify(nf90_inq_dimid(ftn, "lead_time", nlead_id), &
       'nf90_inq_dimid failed for nlead in get_cdf_params')
  call LIS_verify(nf90_inquire_dimension(ftn, nlead_id, len=nlead),&
       'nf90_inquire_dimension failed for nlead in get_cdf_params')

  call LIS_verify(nf90_inq_dimid(ftn, "npercentile", nperc_id), &
       'nf90_inq_dimid failed for npercentile in get_cdf_params')
  call LIS_verify(nf90_inquire_dimension(ftn, nperc_id, len=nperc),&
       'nf90_inquire_dimension failed for npercentile in get_cdf_params')

  allocate(param_hires1(LIS_rc%gnc(n)*LIS_rc%gnr(n),nperc,nlead))
  allocate(param_hires2(LIS_rc%gnc(n),LIS_rc%gnr(n),nperc))
  allocate(param_tmp1(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
  allocate(param_tmp2(LIS_rc%gnc(n),LIS_rc%gnr(n)))

  param_hires1 = -9999.0
  param_hires2 = -9999.0
  param_tmp1 = -9999.0
  param_tmp2 = -9999.0

  ! read landmask
  call LIS_verify(nf90_inq_varid(ftn, 'landmask', mask_id), &
       'nf90_inq_varid failed for landmask in get_cdf_params')
  call LIS_verify(nf90_get_var(ftn, mask_id, param_tmp1), &
       'nf90_get_var failed for landmask in get_cdf_params')

  do r = 1, LIS_rc%gnr(n)
     do c = 1, LIS_rc%gnc(n)
        param_tmp2(c,r) = param_tmp1(c+(r-1)*LIS_rc%gnc(n))
     enddo
  enddo

  !subsets the data for each processor's domain
  do r = LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
     do c = LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
        c1 = c - LIS_ews_halo_ind(n,LIS_localPet+1)+1
        r1 = r - LIS_nss_halo_ind(n,LIS_localPet+1)+1
        gid = LIS_domain(n)%gindex(c1,r1)
        if (gid .ne. -1) then
           landmask(gid) = param_tmp2(c,r)
        endif
     enddo
  enddo

  ! read model cdf
  call LIS_verify(nf90_inq_varid(ftn, 'model_cdf', m_cdf_id), &
       'nf90_inq_varid failed for model_cdf in get_cdf_params')
  call LIS_verify(nf90_get_var(ftn, m_cdf_id, param_hires1,&
       start=(/1,1,1/), count=(/ngrid,nperc,nlead/)), &
       'nf90_get_var failed for model_cdf in get_cdf_params')

  do l = 1, nlead
     do r = 1, LIS_rc%gnr(n)
        do c = 1, LIS_rc%gnc(n)
           param_hires2(c,r,:) = param_hires1(c+(r-1)*LIS_rc%gnc(n),:,l)
        enddo
     enddo

     !subsets the data for each processor's domain
     do r = LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c = LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           c1 = c - LIS_ews_halo_ind(n,LIS_localPet+1)+1
           r1 = r - LIS_nss_halo_ind(n,LIS_localPet+1)+1
           gid = LIS_domain(n)%gindex(c1,r1)
           if (gid .ne. -1) then
              model_cdf(gid,:,l) = param_hires2(c,r,:)
           endif
        enddo
     enddo
  enddo

  ! read ref cdf
  call LIS_verify(nf90_inq_varid(ftn, 'ref_cdf', r_cdf_id), &
       'nf90_inq_varid failed for ref_cdf in get_cdf_params')
  call LIS_verify(nf90_get_var(ftn, r_cdf_id, param_hires1), &
       'nf90_get_var failed for ref_cdf in get_cdf_params')

  do l = 1, nlead
     do r = 1, LIS_rc%gnr(n)
        do c = 1, LIS_rc%gnc(n)
           param_hires2(c,r,:) = param_hires1(c+(r-1)*LIS_rc%gnc(n),:,l)
        enddo
     enddo

     !subsets the data for each processor's domain
     do r = LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c = LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           c1 = c - LIS_ews_halo_ind(n,LIS_localPet+1)+1
           r1 = r - LIS_nss_halo_ind(n,LIS_localPet+1)+1
           gid = LIS_domain(n)%gindex(c1,r1)
           if (gid .ne. -1) then
              ref_cdf(gid,:,l) = param_hires2(c,r,:)
           endif
        enddo
     enddo
  enddo

  call LIS_verify(nf90_close(ftn), 'failed to close in get_cdf_params')

  deallocate(param_hires1)
  deallocate(param_hires2)
  deallocate(param_tmp1)
  deallocate(param_tmp2)

  write(LIS_logunit,*) &
       '[INFO] Done reading GALWEM-GE bias correction parameters data '

end subroutine get_cdf_params

