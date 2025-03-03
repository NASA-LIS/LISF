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
! !ROUTINE: read_COAMPSout
!  \label{read_COAMPSout}
!
! !REVISION HISTORY:
!  14 Mar 2013: Sujay Kumar; Initial Specification in LIS
!  10 Jun 2021: Mahdi Navari; several bugs fixed
! 
! !INTERFACE:
subroutine read_COAMPSout(n, findex, order, fname, ferror)
! !USES:
  use LIS_coreMod
  use LIS_metforcingMod, only : LIS_forc
  use LIS_logMod
  use COAMPSout_forcingMod, only : COAMPSout_struc
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
!  COAMPS output files, transforms into 9 LIS forcing 
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
!    name of the COAMPS output file
!  \item[ferror]
!    flag to indicate success of the call (=1 indicates success)
!  \end{description}
!EOP
  integer                 :: c,r,gindex
  integer                 :: ftn, ios
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

  real                    :: t2scal,q2scal,u10scal,v10scal,swdscal,glwscal
  real                    :: psfcscal,rainncscal
  real                    :: t2offset,q2offset,u10offset,v10offset
  real                    :: psfcoffset,rainncoffset,swdoffset,glwoffset

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  inquire (file=trim(fname), exist=file_exists)
  if (file_exists) then      
     ferror = 1
     call LIS_verify(nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn), 'nf90_open failed in read_COAMPSout')
     call LIS_verify(nf90_inq_varid(ftn,'air_temp_2m',t2Id),&
          'nf90_inq_varid failed for air_temp_2m in read_COAMPSout')
     call LIS_verify(nf90_get_att(ftn,t2Id,"scale_factor",t2scal),&
          'nf90_get_att failed for nf90_get_att')
     call LIS_verify(nf90_get_att(ftn,t2Id,"add_offset",t2offset),&
          'nf90_get_att failed for nf90_get_att')

     call LIS_verify(nf90_inq_varid(ftn,'spec_hum_2m',q2Id),&
          'nf90_inq_varid failed for spec_hum_2m in read_COAMPSout')
     call LIS_verify(nf90_get_att(ftn,q2Id,"scale_factor",q2scal),&
          'nf90_get_att failed for nf90_get_att')
     call LIS_verify(nf90_get_att(ftn,q2Id,"add_offset",q2offset),&
          'nf90_get_att failed for nf90_get_att')

     !call LIS_verify(nf90_inq_varid(ftn,'sw_rad_down',swdownId),&
     !     'nf90_inq_varid failed for sw_rad_down in read_COAMPSout')
     ios = nf90_inq_varid(ftn,'sw_rad_down',swdownId)
     if ( ios /= 0) then
         ios = nf90_inq_varid(ftn,'sw_flux_dn',swdownId)
     endif
     call LIS_verify(ios,'nf90_inq_varid failed for sw_rad_down or sw_flux_dn in read_COAMPSout') 
     call LIS_verify(nf90_get_att(ftn,swdownId,"scale_factor",swdscal),&
          'nf90_get_att failed for nf90_get_att')
     call LIS_verify(nf90_get_att(ftn,swdownId,"add_offset",swdoffset),&
          'nf90_get_att failed for nf90_get_att')

     !call LIS_verify(nf90_inq_varid(ftn,'lw_rad_down',glwId),&
     !     'nf90_inq_varid failed for lw_rad_down in read_COAMPSout')
     ios = nf90_inq_varid(ftn,'lw_rad_down',glwId)
     if ( ios /= 0) then
         ios = nf90_inq_varid(ftn,'lw_flux_dn',glwId)
     endif
     call LIS_verify(ios,'nf90_inq_varid failed for lw_rad_down or lw_flux_dn in read_COAMPSout')
     call LIS_verify(nf90_get_att(ftn,glwId,"scale_factor",glwscal),&
          'nf90_get_att failed for nf90_get_att')
     call LIS_verify(nf90_get_att(ftn,glwId,"add_offset",glwoffset),&
          'nf90_get_att failed for nf90_get_att')

     call LIS_verify(nf90_inq_varid(ftn,'wind_10m_x',u10Id),&
          'nf90_inq_varid failed for wind_10m_x in read_COAMPSout')
     call LIS_verify(nf90_get_att(ftn,u10Id,"scale_factor",u10scal),&
          'nf90_get_att failed for nf90_get_att')
     call LIS_verify(nf90_get_att(ftn,u10Id,"add_offset",u10offset),&
          'nf90_get_att failed for nf90_get_att')

     call LIS_verify(nf90_inq_varid(ftn,'wind_10m_y',v10Id),&
          'nf90_inq_varid failed for wind_10m_y in read_COAMPSout')
     call LIS_verify(nf90_get_att(ftn,v10Id,"scale_factor",v10scal),&
          'nf90_get_att failed for nf90_get_att')
     call LIS_verify(nf90_get_att(ftn,v10Id,"add_offset",v10offset),&
          'nf90_get_att failed for nf90_get_att')

     call LIS_verify(nf90_inq_varid(ftn,'surf_atm_pres',psfcId),&
          'nf90_inq_varid failed for surf_atm_pres in read_COAMPSout')
     call LIS_verify(nf90_get_att(ftn,psfcId,"scale_factor",psfcscal),&
          'nf90_get_att failed for nf90_get_att')
     call LIS_verify(nf90_get_att(ftn,psfcId,"add_offset",psfcoffset),&
          'nf90_get_att failed for nf90_get_att')

     call LIS_verify(nf90_inq_varid(ftn,'ttl_prcp',rainncId),&
          'nf90_inq_varid failed for ttl_prcp in read_COAMPSout')
     call LIS_verify(nf90_get_att(ftn,rainncId,"scale_factor",rainncscal),&
          'nf90_get_att failed for nf90_get_att')
     call LIS_verify(nf90_get_att(ftn,rainncId,"add_offset",rainncoffset),&
          'nf90_get_att failed for nf90_get_att')

     rainc_exists = .false.
     rainc = 0.0 ! initialize to 0 
     call LIS_verify(nf90_get_var(ftn,t2id,gvar),&
          'nf90_get_var failed for t2 in read_COAMPSout')
     do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           if(gvar(c,r,1).gt.-32768) then               
              t2(c-LIS_ews_halo_ind(n,LIS_localPet+1)+1,&
                   r-LIS_nss_halo_ind(n,LIS_localPet+1)+1) = &
                   gvar(c,r,1)*t2scal+t2offset
           endif
        enddo
     enddo
     
     call LIS_verify(nf90_get_var(ftn,q2id,gvar),&
          'nf90_get_var failed for q2 in read_COAMPSout')

     do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           if(gvar(c,r,1).gt.-32768) then               
              q2(c-LIS_ews_halo_ind(n,LIS_localPet+1)+1,&
                   r-LIS_nss_halo_ind(n,LIS_localPet+1)+1) = &
                   gvar(c,r,1)*q2scal+q2offset
           endif
        enddo
     enddo

     call LIS_verify(nf90_get_var(ftn,swdownid,gvar),&
          'nf90_get_var failed for swdown in read_COAMPSout')

     do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           if(gvar(c,r,1).gt.-32768) then               
              swdown(c-LIS_ews_halo_ind(n,LIS_localPet+1)+1,&
                   r-LIS_nss_halo_ind(n,LIS_localPet+1)+1) = &
                   gvar(c,r,1)*swdscal+swdoffset
           endif
        enddo
     enddo

     call LIS_verify(nf90_get_var(ftn,glwid,gvar),&
          'nf90_get_var failed for glw in read_COAMPSout')

     do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           if(gvar(c,r,1).gt.-32768) then               
              glw(c-LIS_ews_halo_ind(n,LIS_localPet+1)+1,&
                   r-LIS_nss_halo_ind(n,LIS_localPet+1)+1) = &
                   gvar(c,r,1)*glwscal+glwoffset
           endif
        enddo
     enddo

     call LIS_verify(nf90_get_var(ftn,u10id,gvar),&
          'nf90_get_var failed for u10 in read_COAMPSout')
     do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           if(gvar(c,r,1).gt.-32768) then               
              u10(c-LIS_ews_halo_ind(n,LIS_localPet+1)+1,&
                   r-LIS_nss_halo_ind(n,LIS_localPet+1)+1) = &
                   gvar(c,r,1)*u10scal+u10offset
           endif
        enddo
     enddo

     call LIS_verify(nf90_get_var(ftn,v10id,gvar),&
          'nf90_get_var failed for v10 in read_COAMPSout')
     do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           if(gvar(c,r,1).gt.-32768) then               
              v10(c-LIS_ews_halo_ind(n,LIS_localPet+1)+1,&
                   r-LIS_nss_halo_ind(n,LIS_localPet+1)+1) = &
                   gvar(c,r,1)*v10scal+v10offset
           endif
        enddo
     enddo

     call LIS_verify(nf90_get_var(ftn,psfcid,gvar),&
          'nf90_get_var failed for psfc in read_COAMPSout')
     do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           if(gvar(c,r,1).gt.-32768) then               
              psfc(c-LIS_ews_halo_ind(n,LIS_localPet+1)+1,&
                   r-LIS_nss_halo_ind(n,LIS_localPet+1)+1) = &
                   gvar(c,r,1)*psfcscal+psfcoffset
           endif
        enddo
     enddo

     call LIS_verify(nf90_get_var(ftn,rainncid,gvar),&
          'nf90_get_var failed for rainnc in read_COAMPSout')
     do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
        do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
           if(gvar(c,r,1).gt.-32768) then               
              rainnc(c-LIS_ews_halo_ind(n,LIS_localPet+1)+1,&
                   r-LIS_nss_halo_ind(n,LIS_localPet+1)+1) = &
                   gvar(c,r,1)*rainncscal+rainncoffset
           endif
        enddo
     enddo

     call LIS_verify(nf90_close(ftn), &
          'failed to close file in read_COAMPSout')

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if ( LIS_domain(n)%gindex(c,r) /= -1 ) then 
              gindex = LIS_domain(n)%gindex(c,r)
              if ( order == 1 ) then
                 COAMPSout_struc(n)%metdata1(1,gindex) = t2(c,r)
                 COAMPSout_struc(n)%metdata1(2,gindex) = q2(c,r)
                 COAMPSout_struc(n)%metdata1(3,gindex) = swdown(c,r)
                 COAMPSout_struc(n)%metdata1(4,gindex) = glw(c,r)
                 COAMPSout_struc(n)%metdata1(5,gindex) = u10(c,r)
                 COAMPSout_struc(n)%metdata1(6,gindex) = v10(c,r)
                 COAMPSout_struc(n)%metdata1(7,gindex) = psfc(c,r)
                 if ( rainc_exists ) then
                    COAMPSout_struc(n)%metdata1(8,gindex) = rainnc(c,r)+ &
                                                            rainc(c,r)
                 else
                    COAMPSout_struc(n)%metdata1(8,gindex) = rainnc(c,r)
                 endif
              else
                 COAMPSout_struc(n)%metdata2(1,gindex) = t2(c,r)
                 COAMPSout_struc(n)%metdata2(2,gindex) = q2(c,r)
                 COAMPSout_struc(n)%metdata2(3,gindex) = swdown(c,r)
                 COAMPSout_struc(n)%metdata2(4,gindex) = glw(c,r)
                 COAMPSout_struc(n)%metdata2(5,gindex) = u10(c,r)
                 COAMPSout_struc(n)%metdata2(6,gindex) = v10(c,r)
                 COAMPSout_struc(n)%metdata2(7,gindex) = psfc(c,r)
                 if ( rainc_exists ) then
                    COAMPSout_struc(n)%metdata2(8,gindex) = rainnc(c,r)+ &
                                                            rainc(c,r)
                 else
                    COAMPSout_struc(n)%metdata2(8,gindex) = rainnc(c,r)
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
  write(LIS_logunit,*) '[ERR] read_COAMPSout requires NetCDF'
  write(LIS_logunit,*) '[ERR] please recompile LIS'
  call LIS_endrun
#endif
end subroutine read_COAMPSout
