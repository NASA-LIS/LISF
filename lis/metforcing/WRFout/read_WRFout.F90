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
! !ROUTINE: read_WRFout
!  \label{read_WRFout}
!
! !REVISION HISTORY:
!  14 Mar 2013: Sujay Kumar; Initial Specification in LIS
! 
! !INTERFACE:
subroutine read_WRFout(n, findex, order, fname, ferror)
! !USES:
  use LIS_coreMod
  use LIS_metforcingMod, only : LIS_forc
  use LIS_logMod
  use WRFout_forcingMod, only : WRFout_struc
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: findex
  integer, intent(in)          :: order
  character(len=*), intent(in) :: fname
  integer, intent(out)         :: ferror
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  WRF output files, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing source
!  \item[order]
!    flag indicating which data to be read (order=1, read for the previous 
!    1hr bookend, order=2, read for the next 1hr bookend)
!  \item[fname]
!    name of the WRF output file
!  \item[ferror]
!    flag to indicate success of the call (=1 indicates success)
!  \end{description}
!EOP
  integer                 :: c,r,gindex
  integer                 :: ftn
  logical                 :: file_exists, rainc_exists
  integer                 :: t2Id, q2Id, swdownId, glwId
  integer                 :: u10Id, v10Id, psfcId, rainncId
  integer                 :: raincID
  real                    :: gvar(LIS_rc%gnc(n),LIS_rc%gnr(n),1)
  real                    :: t2(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: q2(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: swdown(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: glw(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: u10(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: v10(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: psfc(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: rainnc(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: rainc(LIS_rc%lnc(n),LIS_rc%lnr(n))

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  inquire (file=trim(fname), exist=file_exists)
  if (file_exists) then      
     ferror = 1
     call LIS_verify(nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn), 'nf90_open failed in read_WRFout')
     call LIS_verify(nf90_inq_varid(ftn,'T2',t2Id),&
          'nf90_inq_varid failed for T2 in read_WRFout')
     call LIS_verify(nf90_inq_varid(ftn,'Q2',q2Id),&
          'nf90_inq_varid failed for Q2 in read_WRFout')
     call LIS_verify(nf90_inq_varid(ftn,'SWDOWN',swdownId),&
          'nf90_inq_varid failed for SWDOWN in read_WRFout')
     call LIS_verify(nf90_inq_varid(ftn,'GLW',glwId),&
          'nf90_inq_varid failed for GLW in read_WRFout')
     call LIS_verify(nf90_inq_varid(ftn,'U10',u10Id),&
          'nf90_inq_varid failed for U10 in read_WRFout')
     call LIS_verify(nf90_inq_varid(ftn,'V10',v10Id),&
          'nf90_inq_varid failed for V10 in read_WRFout')
     call LIS_verify(nf90_inq_varid(ftn,'PSFC',psfcId),&
          'nf90_inq_varid failed for PSFC in read_WRFout')
     call LIS_verify(nf90_inq_varid(ftn,'RAINNC',rainncId),&
          'nf90_inq_varid failed for RAINNC in read_WRFout')
     if ( nf90_inq_varid(ftn,'RAINC',raincId) == 0 ) then
        rainc_exists = .false.
        write(LIS_logunit,*) '[INFO] read_WRFout -- optional RAINC not present'
     else
        rainc_exists = .true.
     endif

     call LIS_verify(nf90_get_var(ftn,t2id,gvar),&
          'nf90_get_var failed for t2 in read_WRFout')
     t2(:,:) = gvar(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1),1)

     call LIS_verify(nf90_get_var(ftn,q2id,gvar),&
          'nf90_get_var failed for q2 in read_WRFout')
     q2(:,:) = gvar(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1),1)

     call LIS_verify(nf90_get_var(ftn,swdownid,gvar),&
          'nf90_get_var failed for swdown in read_WRFout')
     swdown(:,:) = gvar(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1),1)

     call LIS_verify(nf90_get_var(ftn,glwid,gvar),&
          'nf90_get_var failed for glw in read_WRFout')
     glw(:,:) = gvar(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1),1)

     call LIS_verify(nf90_get_var(ftn,u10id,gvar),&
          'nf90_get_var failed for u10 in read_WRFout')
     u10(:,:) = gvar(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1),1)

     call LIS_verify(nf90_get_var(ftn,v10id,gvar),&
          'nf90_get_var failed for v10 in read_WRFout')
     v10(:,:) = gvar(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1),1)

     call LIS_verify(nf90_get_var(ftn,psfcid,gvar),&
          'nf90_get_var failed for psfc in read_WRFout')
     psfc(:,:) = gvar(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1),1)

     call LIS_verify(nf90_get_var(ftn,rainncid,gvar),&
          'nf90_get_var failed for rainnc in read_WRFout')
     rainnc(:,:) = gvar(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1),1)

     if ( rainc_exists ) then
        call LIS_verify(nf90_get_var(ftn,raincid,gvar),&
             'nf90_get_var failed for rainc in read_WRFout')
        rainc(:,:) = gvar(LIS_ews_halo_ind(n,LIS_localPet+1):&         
             LIS_ewe_halo_ind(n,LIS_localPet+1), &
             LIS_nss_halo_ind(n,LIS_localPet+1): &
             LIS_nse_halo_ind(n,LIS_localPet+1),1)
     endif

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if ( LIS_domain(n)%gindex(c,r) /= -1 ) then 
              gindex = LIS_domain(n)%gindex(c,r)
              if ( order == 1 ) then
                 WRFout_struc(n)%metdata1(1,gindex) = t2(c,r)
                 WRFout_struc(n)%metdata1(2,gindex) = q2(c,r)
                 WRFout_struc(n)%metdata1(3,gindex) = swdown(c,r)
                 WRFout_struc(n)%metdata1(4,gindex) = glw(c,r)
                 WRFout_struc(n)%metdata1(5,gindex) = u10(c,r)
                 WRFout_struc(n)%metdata1(6,gindex) = v10(c,r)
                 WRFout_struc(n)%metdata1(7,gindex) = psfc(c,r)
                 if ( rainc_exists ) then
                    WRFout_struc(n)%metdata1(8,gindex) = rainnc(c,r)+ &
                                                            rainc(c,r)
                 else
                    WRFout_struc(n)%metdata1(8,gindex) = rainnc(c,r)
                 endif
              else
                 WRFout_struc(n)%metdata2(1,gindex) = t2(c,r)
                 WRFout_struc(n)%metdata2(2,gindex) = q2(c,r)
                 WRFout_struc(n)%metdata2(3,gindex) = swdown(c,r)
                 WRFout_struc(n)%metdata2(4,gindex) = glw(c,r)
                 WRFout_struc(n)%metdata2(5,gindex) = u10(c,r)
                 WRFout_struc(n)%metdata2(6,gindex) = v10(c,r)
                 WRFout_struc(n)%metdata2(7,gindex) = psfc(c,r)
                 if ( rainc_exists ) then
                    WRFout_struc(n)%metdata2(8,gindex) = rainnc(c,r)+ &
                                                            rainc(c,r)
                 else
                    WRFout_struc(n)%metdata2(8,gindex) = rainnc(c,r)
                 endif
              endif
           endif
        enddo
     enddo

  else
     write(LIS_logunit,*) '[ERR] Forcing file '//trim(fname)//' not found'
     ferror = 0 
  endif
#else
  write(LIS_logunit,*) '[ERR] read_WRFout requires NetCDF'
  write(LIS_logunit,*) '[ERR] please recompile LIS'
  call LIS_endrun
#endif
end subroutine read_WRFout

