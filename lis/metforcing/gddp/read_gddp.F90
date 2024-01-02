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
! !ROUTINE: read_gddp
!  \label{read_gddp}
!
! !REVISION HISTORY:
!  03 Feb 2022: Sujay Kumar, Initial specification
! 
! !INTERFACE:

subroutine read_gddp(n, findex, order, year, doy, &
     names,ref_dailyclimofile, &
     ref_hourlyclimofile, ferror)
  
! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_constantsMod
  use LIS_metforcingMod
  use gddp_forcingMod
  
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: findex  ! Forcing index
  integer, intent(in)          :: order
  integer, intent(in)          :: year 
  integer, intent(in)          :: doy
  character(len=*), intent(in) :: names(gddp_struc(n)%nmodels,7)
  character(len=*), intent(in) :: ref_dailyclimofile
  character(len=*), intent(in) :: ref_hourlyclimofile(7)
  integer, intent(out)         :: ferror
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  GDDP data, transforms into LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  A (MERRA2-based) climatology temporal disaggregation of the daily
!  data to hourly is performed.
!  The code expects the availability of the climatology of the reference
!  data generated at hourly and daily timescales.
!
!    The temporal disaggreation is performed as:
!
!      GDDP_new = GDDP_orig * REF_H/REF_D
!
!   where GDDP_new is the temporally disaggreated data
!         GDDP_orig is the original GDDP data
!         REF_H is the hourly climatology data
!         REF_D is the daily climatology data
!  
!EOP

  integer                   :: iv, ftn
  integer                   :: ngddp
  integer                   :: k,t,m,c,r,iret,rc
  integer                   :: var_index
  real                      :: missingValue 
  integer                   :: nvars
  integer                   :: tasid,hussid,rsdsid,rldsid
  integer                   :: sfcWindid,prid,hursid,psid
  integer                   :: tdimId, tdims
  integer                   :: igrib
  logical                   :: pcp_flag, var_found
  logical                   :: var_status(11)
  logical                   :: file_exists
  logical*1, allocatable    :: lb(:)
  real, allocatable         :: f(:)
  logical                   :: found
  integer                   :: inpnc,inpnr,rad
  integer                   :: c1,c2,r1,r2,i,j
  integer                   :: doy_upd
  real                      :: tas_c, es, e
  logical                   :: read_flag
  logical                   :: leap_year
  real                      :: tas_inp(gddp_struc(n)%nc,gddp_struc(n)%nr)
  real                      :: huss_inp(gddp_struc(n)%nc,gddp_struc(n)%nr)
  real                      :: hurs_inp(gddp_struc(n)%nc,gddp_struc(n)%nr)
  real                      :: rsds_inp(gddp_struc(n)%nc,gddp_struc(n)%nr)
  real                      :: rlds_inp(gddp_struc(n)%nc,gddp_struc(n)%nr)
  real                      :: sfcWind_inp(gddp_struc(n)%nc,gddp_struc(n)%nr)
  real                      :: pr_inp(gddp_struc(n)%nc,gddp_struc(n)%nr)
  real                      :: psurf_inp(gddp_struc(n)%nc,gddp_struc(n)%nr)

  real                      :: tair_climo1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                      :: qair_climo1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                      :: swdown_climo1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                      :: lwdown_climo1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                      :: wind_climo1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                      :: psurf_climo1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                      :: prcp_climo1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  
  real                      :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))

  real                      :: tmpval
  ferror = 1
  iv = 0

  inpnc = gddp_struc(n)%nc
  inpnr = gddp_struc(n)%nr
  
  ngddp = (gddp_struc(n)%nc*gddp_struc(n)%nr)
  
!  allocate(gddp_forcing(gddp_struc(n)%nc*gddp_struc(n)%nr,11))

  if((mod(year,4) .eq. 0 .and. mod(year, 100).ne.0) &!leap year
         .or.(mod(year,400) .eq.0)) then 
     leap_year = .true. 
  else 
     leap_year = .false. 
  endif
  
  varfield = 0 
  ferror = 1
  var_status = .false.

  read_flag = .false. 
  if(order.eq.1) then
     if(gddp_struc(n)%day_check1.ne.doy) then 
        read_flag = .true.
        gddp_struc(n)%day_check1 = doy
     endif
  else
     if(gddp_struc(n)%day_check2.ne.doy) then 
        read_flag = .true.
        gddp_struc(n)%day_check2 = doy
     endif
  endif

#if (defined USE_NETCDF3 || defined USE_NETCDF4)  

  inquire (file=ref_hourlyclimofile(1), exist=file_exists)

  if (file_exists) then

     call LIS_verify(nf90_open(path=ref_hourlyclimofile(1),&
          mode=nf90_nowrite,ncid=ftn),&
          'Error in opening file '//trim(ref_hourlyclimofile(1)))
     call LIS_verify(nf90_inq_varid(ftn,'TAIR',tasid),&
          'Error in nf90_inq_varid: TAIR')
     call LIS_verify(nf90_get_var(ftn,tasid,tair_climo1,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'Error in nf90_get_var:TAIR')
     call LIS_verify(nf90_close(ftn))

     call LIS_verify(nf90_open(path=ref_hourlyclimofile(2),&
          mode=nf90_nowrite,ncid=ftn),&
          'Error in opening file '//trim(ref_hourlyclimofile(2)))     
     call LIS_verify(nf90_inq_varid(ftn,'QAIR',hussid),&
          'Error in nf90_inq_varid: QAIR')
     call LIS_verify(nf90_get_var(ftn,hussid,qair_climo1,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'Error in nf90_get_var:QAIR')
     call LIS_verify(nf90_close(ftn))

     call LIS_verify(nf90_open(path=ref_hourlyclimofile(3),&
          mode=nf90_nowrite,ncid=ftn),&
          'Error in opening file '//trim(ref_hourlyclimofile(3)))
     call LIS_verify(nf90_inq_varid(ftn,'SWDOWN',rsdsid),&
          'Error in nf90_inq_varid: SWDOWN')
     call LIS_verify(nf90_get_var(ftn,rsdsid,swdown_climo1,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'Error in nf90_get_var:SWDOWN')
     call LIS_verify(nf90_close(ftn))

     call LIS_verify(nf90_open(path=ref_hourlyclimofile(4),&
          mode=nf90_nowrite,ncid=ftn),&
          'Error in opening file '//trim(ref_hourlyclimofile(4)))
     call LIS_verify(nf90_inq_varid(ftn,'LWDOWN',rldsid),&
          'Error in nf90_inq_varid: LWDOWN')
     call LIS_verify(nf90_get_var(ftn,rldsid,lwdown_climo1,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&          
          'Error in nf90_get_var:LWDOWN')
     call LIS_verify(nf90_close(ftn))

     call LIS_verify(nf90_open(path=ref_hourlyclimofile(5),&
          mode=nf90_nowrite,ncid=ftn),&
          'Error in opening file '//trim(ref_hourlyclimofile(5)))     
     call LIS_verify(nf90_inq_varid(ftn,'WIND',sfcwindid),&
          'Error in nf90_inq_varid: WIND')
     call LIS_verify(nf90_get_var(ftn,sfcwindid,wind_climo1,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'Error in nf90_get_var:WIND')
     call LIS_verify(nf90_close(ftn))

     call LIS_verify(nf90_open(path=ref_hourlyclimofile(6),&
          mode=nf90_nowrite,ncid=ftn),&
          'Error in opening file '//trim(ref_hourlyclimofile(6)))     
     call LIS_verify(nf90_inq_varid(ftn,'PSURF',psid),&
          'Error in nf90_inq_varid: PSURF')
     call LIS_verify(nf90_get_var(ftn,psid,psurf_climo1,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&          
          'Error in nf90_get_var:PSURF')
     call LIS_verify(nf90_close(ftn))

     call LIS_verify(nf90_open(path=ref_hourlyclimofile(7),&
          mode=nf90_nowrite,ncid=ftn),&
          'Error in opening file '//trim(ref_hourlyclimofile(7)))          
     call LIS_verify(nf90_inq_varid(ftn,'PRCP',prid),&
          'Error in nf90_inq_varid: PRCP')
     call LIS_verify(nf90_get_var(ftn,prid,prcp_climo1,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'Error in nf90_get_var:PRCP')     
     call LIS_verify(nf90_close(ftn))
     
     
  else
     write(LIS_logunit,*) '[ERR] climatology file '//trim(ref_hourlyclimofile(1))
     write(LIS_logunit,*) '[ERR] does not exist '
     call LIS_endrun()
  endif

  if(read_flag) then 
     inquire (file=ref_dailyclimofile, exist=file_exists)
     if (file_exists) then
        
        call LIS_verify(nf90_open(path=ref_dailyclimofile,mode=nf90_nowrite,ncid=ftn),&
             'Error in opening file '//trim(ref_dailyclimofile))
        call LIS_verify(nf90_inq_varid(ftn,'TAIR',tasid),&
             'Error in nf90_inq_varid: TAIR')
        call LIS_verify(nf90_get_var(ftn,tasid,gddp_struc(n)%tair_climo,&
             start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
             LIS_nss_halo_ind(n,LIS_localPet+1)/),&
             count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
             'Error in nf90_get_var:TAIR')
        
        call LIS_verify(nf90_inq_varid(ftn,'QAIR',hussid),&
             'Error in nf90_inq_varid: QAIR')
        call LIS_verify(nf90_get_var(ftn,hussid,gddp_struc(n)%qair_climo,&
             start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
             LIS_nss_halo_ind(n,LIS_localPet+1)/),&
             count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
             'Error in nf90_get_var:QAIR')
        
        call LIS_verify(nf90_inq_varid(ftn,'SWDOWN',rsdsid),&
             'Error in nf90_inq_varid: SWDOWN')
        call LIS_verify(nf90_get_var(ftn,rsdsid,gddp_struc(n)%swdown_climo,&
             start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
             LIS_nss_halo_ind(n,LIS_localPet+1)/),&
             count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
             'Error in nf90_get_var:SWDOWN')
        
        call LIS_verify(nf90_inq_varid(ftn,'LWDOWN',rldsid),&
             'Error in nf90_inq_varid: LWDOWN')
        call LIS_verify(nf90_get_var(ftn,rldsid,gddp_struc(n)%lwdown_climo,&
             start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
             LIS_nss_halo_ind(n,LIS_localPet+1)/),&
             count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&          
             'Error in nf90_get_var:LWDOWN')
        
        call LIS_verify(nf90_inq_varid(ftn,'WIND',sfcwindid),&
             'Error in nf90_inq_varid: WIND')
        call LIS_verify(nf90_get_var(ftn,sfcwindid,gddp_struc(n)%wind_climo,&
             start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
             LIS_nss_halo_ind(n,LIS_localPet+1)/),&
             count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
             'Error in nf90_get_var:WIND')
        
        call LIS_verify(nf90_inq_varid(ftn,'PSURF',psid),&
             'Error in nf90_inq_varid: PSURF')
        call LIS_verify(nf90_get_var(ftn,psid,gddp_struc(n)%psurf_climo,&
             start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
             LIS_nss_halo_ind(n,LIS_localPet+1)/),&
             count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&          
             'Error in nf90_get_var:PSURF')
        
        call LIS_verify(nf90_inq_varid(ftn,'PRCP',prid),&
             'Error in nf90_inq_varid: PRCP')
        call LIS_verify(nf90_get_var(ftn,prid,gddp_struc(n)%prcp_climo,&
             start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
             LIS_nss_halo_ind(n,LIS_localPet+1)/),&
             count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
             'Error in nf90_get_var:PRCP')
        
        call LIS_verify(nf90_close(ftn))
        
        
     else
        write(LIS_logunit,*) '[ERR] climatology file '//trim(ref_dailyclimofile)
        write(LIS_logunit,*) '[ERR] does not exist '
        call LIS_endrun()
     endif
     
     do m=1,gddp_struc(n)%nmodels
        inquire (file=names(m,1), exist=file_exists)
        if (file_exists) then
           write(LIS_logunit,*)'[INFO] reading.. ',trim(names(m,1))
           call LIS_verify(nf90_open(path=names(m,1),mode=nf90_nowrite,ncid=ftn),&
                'Error in opening file '//trim(names(m,1)))
           call LIS_verify(nf90_inq_dimid(ftn,'time',tdimID),&
                'Error in opening nf90_inq_dimid')
           call LIS_verify(nf90_inquire_dimension(ftn,tdimId,len=tdims),&
                'Error in nf90_inquire_dimension')
           
           call LIS_verify(nf90_inq_varid(ftn,'tas',tasid),&
                'Error in nf90_inq_varid: tas')
           doy_upd = doy
           if(leap_year.and.tdims.ne.366.and.tdims.eq.365) then
              if(doy.gt.60) then !feb 29th onwards
                 doy_upd = doy-1   
              endif
           endif
           if(tdims.eq.360) then
              if(leap_year) then !distribute 6 days
                 if(doy.ge.366) then
                    doy_upd = doy-6
                 elseif(doy.ge.305) then
                    doy_upd = doy-5
                 elseif(doy.ge.244) then
                    doy_upd = doy-4
                 elseif(doy.ge.183) then
                    doy_upd = doy-3
                 elseif(doy.ge.122) then
                    doy_upd = doy-2
                 elseif(doy.ge.61) then
                    doy_upd = doy-1
                 endif
              else !distribute 5 days
                 if(doy.ge.365) then
                    doy_upd = doy-5
                 elseif(doy.ge.292) then
                    doy_upd = doy-4
                 elseif(doy.ge.219) then
                    doy_upd = doy-3
                 elseif(doy.ge.146) then
                    doy_upd = doy-2
                 elseif(doy.ge.73) then
                    doy_upd = doy-1
                 endif                 
                 
              endif
           elseif(tdims.eq.364) then
              if(leap_year) then !distribute 2 days
                 if(doy.ge.366) then
                    doy_upd = doy-2
                 elseif(doy.ge.183) then
                    doy_upd = doy-1
                 endif
              else !distribute 1 days
                 if(doy.ge.365) then
                    doy_upd = doy-1
                 endif                                  
              endif              
           endif
           
           call LIS_verify(nf90_get_var(ftn,tasid,tas_inp,&
                start=(/1,1,doy_upd/),count=(/inpnc,inpnr,1/)),&
                'Error in nf90_get_var:tas')
           call LIS_verify(nf90_close(ftn))
           
           
           call LIS_verify(nf90_open(path=names(m,2),mode=nf90_nowrite,ncid=ftn),&
                'Error in opening file '//trim(names(m,2)))
           call LIS_verify(nf90_inq_varid(ftn,'huss',hussid),&
                'Error in nf90_inq_varid: huss')
           call LIS_verify(nf90_get_var(ftn,hussid,huss_inp,&
                start=(/1,1,doy_upd/),count=(/inpnc,inpnr,1/)),&
                'Error in nf90_get_var:huss')
           call LIS_verify(nf90_close(ftn))
           
           call LIS_verify(nf90_open(path=names(m,3),mode=nf90_nowrite,ncid=ftn),&
                'Error in opening file '//trim(names(m,3)))
           call LIS_verify(nf90_inq_varid(ftn,'rsds',rsdsid),&
                'Error in nf90_inq_varid: rsds')
           call LIS_verify(nf90_get_var(ftn,rsdsid,rsds_inp,&
                start=(/1,1,doy_upd/),count=(/inpnc,inpnr,1/)),&
                'Error in nf90_get_var:rsds')
           call LIS_verify(nf90_close(ftn))
                      
           call LIS_verify(nf90_open(path=names(m,4),mode=nf90_nowrite,ncid=ftn),&
                'Error in opening file '//trim(names(m,4)))
           call LIS_verify(nf90_inq_varid(ftn,'rlds',rldsid),&
                'Error in nf90_inq_varid: rlds')
           call LIS_verify(nf90_get_var(ftn,rldsid,rlds_inp,&
                start=(/1,1,doy_upd/),count=(/inpnc,inpnr,1/)),&
                'Error in nf90_get_var:rlds')
           call LIS_verify(nf90_close(ftn))
           
           
           call LIS_verify(nf90_open(path=names(m,5),mode=nf90_nowrite,ncid=ftn),&
                'Error in opening file '//trim(names(m,5)))
           call LIS_verify(nf90_inq_varid(ftn,'sfcWind',sfcWindid),&
                'Error in nf90_inq_varid: sfcWind')
           call LIS_verify(nf90_get_var(ftn,sfcWindid,sfcWind_inp,&
                start=(/1,1,doy_upd/),count=(/inpnc,inpnr,1/)),&
                'Error in nf90_get_var:sfcWind')
           call LIS_verify(nf90_close(ftn))
           
           
           call LIS_verify(nf90_open(path=names(m,6),mode=nf90_nowrite,ncid=ftn),&
                'Error in opening file '//trim(names(m,6)))
           call LIS_verify(nf90_inq_varid(ftn,'pr',prid),&
                'Error in nf90_inq_varid: pr')
           call LIS_verify(nf90_get_var(ftn,prid,pr_inp,&
                start=(/1,1,doy_upd/),count=(/inpnc,inpnr,1/)),&
                'Error in nf90_get_var:pr')             
           call LIS_verify(nf90_close(ftn))
           
           call LIS_verify(nf90_open(path=names(m,7),mode=nf90_nowrite,ncid=ftn),&
                'Error in opening file '//trim(names(m,7)))
           call LIS_verify(nf90_inq_varid(ftn,'hurs',hursid),&
                'Error in nf90_inq_varid: hurs')
           call LIS_verify(nf90_get_var(ftn,hursid,hurs_inp,&
                start=(/1,1,doy_upd/),count=(/inpnc,inpnr,1/)),&
                'Error in nf90_get_var:hurs')             
           call LIS_verify(nf90_close(ftn))

           !checking for unphysical specific humidity values
           do r=1,inpnr
              do c=1,inpnc
                 if(huss_inp(c,r).ne.1e+20) then 
                    if(huss_inp(c,r).le.0.) then
                       rad = 1
                       found = .false.
                       do while(.not.found)
                          c1 = max(1,c-rad)
                          c2 = min(inpnc,c+rad)
                          r1 = max(1,r-rad)
                          r2 = min(inpnr,r+rad)
                          
                          do j=r1,r2
                             do i=c1,c2
                                if(huss_inp(i,j).ne.0) then
                                   huss_inp(c,r) = huss_inp(i,j)
                                   found = .true.
                                   exit
                                endif
                             enddo
                             if(found) exit
                          enddo
                          if(.not.found) rad = rad+1
                       enddo
                    elseif(huss_inp(c,r).lt.0) then
                       write(LIS_logunit,*)'[ERR] Found unphysical Qair values that'
                       write(LIS_logunit,*)'[ERR] could not be fixed with neighbor search'
                       call LIS_endrun()
                    endif
                 endif
              enddo
           enddo

           !check for unphysical lwdown values
           do r=1,inpnr
              do c=1,inpnc
                 if(rlds_inp(c,r).ne.1e+20) then 
                    if(rlds_inp(c,r).lt.0) then 
                       rad = 1
                       found = .false.
                       do while(.not.found)
                          c1 = max(1,c-rad)
                          c2 = min(inpnc,c+rad)
                          r1 = max(1,r-rad)
                          r2 = min(inpnr,r+rad)
                          
                          do j=r1,r2
                             do i=c1,c2
                                if(rlds_inp(i,j).ne.1e+20.and.&
                                     rlds_inp(i,j).ge.0) then
                                   rlds_inp(c,r) = rlds_inp(i,j)
                                   found = .true.
                                   exit
                                endif
                             enddo
                             if(found) exit
                          enddo
                          if(.not.found) rad = rad+1
                       enddo
                    elseif(rlds_inp(c,r).lt.0) then 
                       write(LIS_logunit,*)'[ERR] Found unphysical LWdown values that'
                       write(LIS_logunit,*)'[ERR] could not be fixed with neighbor search'
                       call LIS_endrun()
                    endif
                 endif
              enddo
           enddo
           !check for unphysical SWdown values
           do r=1,inpnr
              do c=1,inpnc
                 if(rsds_inp(c,r).ne.1e+20) then 
                    if(rsds_inp(c,r).gt.LIS_CONST_SOLAR.or.&
                         rsds_inp(c,r).lt.0) then 
                       rad = 1
                       found = .false.
                       do while(.not.found)
                          c1 = max(1,c-rad)
                          c2 = min(inpnc,c+rad)
                          r1 = max(1,r-rad)
                          r2 = min(inpnr,r+rad)
                          
                          do j=r1,r2
                             do i=c1,c2
                                if(rsds_inp(i,j).le.LIS_CONST_SOLAR.and.&
                                     rsds_inp(i,j).ne.1e+20.and.&
                                     rsds_inp(i,j).ge.0) then
                                   rsds_inp(c,r) = rsds_inp(i,j)
                                   found = .true.
                                   exit
                                endif
                             enddo
                             if(found) exit
                          enddo
                          if(.not.found) rad = rad+1
                       enddo
                    elseif(rsds_inp(c,r).gt.LIS_CONST_SOLAR) then
                       write(LIS_logunit,*)'[ERR] Found unphysical SWdown values that'
                       write(LIS_logunit,*)'[ERR] could not be fixed with neighbor search'
                       call LIS_endrun()
                    endif
                 endif
              enddo
           enddo

           do r=1,inpnr
              do c=1,inpnc
                 if(tas_inp(c,r).ne.1e+20) then
                    tas_c = tas_inp(c,r)-273.15
                    !August-Roche-Magnus formula
                    es = 6.112*exp((17.67*tas_c)/(tas_c+243.04))
                    e = hurs_inp(c,r)*es/100.0
                    psurf_inp(c,r) = 100*(0.622+0.378*huss_inp(c,r))*e/&
                         (huss_inp(c,r)) !in Pascals                       
                 else
                    psurf_inp(c,r) = 1e+20
                 endif
              enddo
           enddo

           do r=1,inpnr
              do c=1,inpnc
                 if(psurf_inp(c,r).gt.200000.or.psurf_inp(c,r).lt.1000) then 
                    rad = 1
                    found = .false.
                    do while(.not.found)
                       c1 = max(1,c-rad)
                       c2 = min(inpnc,c+rad)
                       r1 = max(1,r-rad)
                       r2 = min(inpnr,r+rad)
                       
                       do j=r1,r2
                          do i=c1,c2
                             if(psurf_inp(i,j).lt.200000.and.&
                                  psurf_inp(i,j).gt.1000.and.&
                                  psurf_inp(i,j).ne.1e+20) then
                                psurf_inp(c,r) = psurf_inp(i,j)
                                tas_inp(c,r) = tas_inp(i,j)
                                huss_inp(c,r) = huss_inp(i,j)
                                hurs_inp(c,r) = hurs_inp(i,j)
                                found = .true.
                                exit
                             endif
                          enddo
                          if(found) exit
                       enddo
                       if(.not.found) rad = rad+1
                    enddo
                 endif
              enddo
           enddo           
        
           allocate(lb(gddp_struc(n)%nc*gddp_struc(n)%nr))
           
           !tair
           pcp_flag = .false.
           
           call interp_gddp(n, 1,findex,  pcp_flag, inpnc,inpnr,&
                tas_inp, & 
                lb, LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield )
           
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    gddp_struc(n)%metdata(m,1,&
                         LIS_domain(n)%gindex(c,r)) &
                         = varfield(c,r)
                 endif
              end do
           enddo
           
           !qair        
           pcp_flag = .false. 
           call interp_gddp(n, 2, findex,  pcp_flag, inpnc,inpnr,&
                huss_inp, & 
                lb, LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield )
           
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    gddp_struc(n)%metdata(m,2,&
                         LIS_domain(n)%gindex(c,r)) &
                         = varfield(c,r)
                 endif
              end do
           enddo
           
           !swdown        
           pcp_flag = .false. 
           call interp_gddp(n, 3, findex,  pcp_flag, inpnc,inpnr,&
                rsds_inp, & 
                lb, LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield )
           
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then
                    if(varfield(c,r).lt.0) varfield(c,r) = 0.0
                    gddp_struc(n)%metdata(m,3,&
                         LIS_domain(n)%gindex(c,r)) &
                         = varfield(c,r)                          
                 endif
              end do
           enddo
           !lwdown        
           pcp_flag = .false. 
           call interp_gddp(n, 4, findex,  pcp_flag, inpnc,inpnr,&
                rlds_inp, & 
                lb, LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield )
           
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    gddp_struc(n)%metdata(m,4,&
                         LIS_domain(n)%gindex(c,r)) &
                         = varfield(c,r)
                 endif
              end do
           enddo
           !wind
           pcp_flag = .false. 
           call interp_gddp(n, 5, findex,  pcp_flag, inpnc,inpnr,&
                sfcwind_inp, & 
                lb, LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield )
           
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    gddp_struc(n)%metdata(m,5,&
                         LIS_domain(n)%gindex(c,r)) &
                         = varfield(c,r)
                 endif
              end do
           enddo
           !psurf        
           pcp_flag = .false. 
           call interp_gddp(n, 6, findex,  pcp_flag, inpnc,inpnr,&
                psurf_inp, & 
                lb, LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield )
           
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    gddp_struc(n)%metdata(m,6,&
                         LIS_domain(n)%gindex(c,r)) &
                         = varfield(c,r)
                 endif
              end do
           enddo
           !precip        
           pcp_flag = .false. 
           call interp_gddp(n, 7, findex,  pcp_flag, inpnc,inpnr,&
                pr_inp, & 
                lb, LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield )
           
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then
                    gddp_struc(n)%metdata(m,7,&
                         LIS_domain(n)%gindex(c,r)) &
                         = varfield(c,r)
                 endif
              end do
           enddo
           
           deallocate(lb)
        else
           write(LIS_logunit,*) &
                '[ERR] Could not find file: ',trim(names(m,1))
           ferror = 0
           call LIS_endrun()
                      
        end if
     enddo
  endif

  !disaggregate to hourly
  do m=1,gddp_struc(n)%nmodels 
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 gddp_struc(n)%metdata1(m,1,&
                      LIS_domain(n)%gindex(c,r)) &
                      = gddp_struc(n)%metdata(m,1,&
                      LIS_domain(n)%gindex(c,r))*&
                      tair_climo1(c,r)/&
                      gddp_struc(n)%tair_climo(c,r)

                 gddp_struc(n)%metdata1(m,2,&
                      LIS_domain(n)%gindex(c,r)) &
                      = gddp_struc(n)%metdata(m,2,&
                      LIS_domain(n)%gindex(c,r))*&
                      qair_climo1(c,r)/&
                      gddp_struc(n)%qair_climo(c,r)

                 if(gddp_struc(n)%swdown_climo(c,r).ne.0) then 
                    gddp_struc(n)%metdata1(m,3,&
                         LIS_domain(n)%gindex(c,r)) &
                         = gddp_struc(n)%metdata(m,3,&
                         LIS_domain(n)%gindex(c,r))*&
                         swdown_climo1(c,r)/&
                         gddp_struc(n)%swdown_climo(c,r)                 
                 else
                    
                    gddp_struc(n)%metdata1(m,3,&
                         LIS_domain(n)%gindex(c,r)) = 0
                    
                 endif
                 gddp_struc(n)%metdata1(m,4,&
                      LIS_domain(n)%gindex(c,r)) &
                      = gddp_struc(n)%metdata(m,4,&
                      LIS_domain(n)%gindex(c,r))*&
                      lwdown_climo1(c,r)/&
                      gddp_struc(n)%lwdown_climo(c,r)
                 
                 gddp_struc(n)%metdata1(m,5,&
                      LIS_domain(n)%gindex(c,r)) &
                      = gddp_struc(n)%metdata(m,5,&
                      LIS_domain(n)%gindex(c,r))*&
                      wind_climo1(c,r)/&
                      gddp_struc(n)%wind_climo(c,r)                 

                 gddp_struc(n)%metdata1(m,6,&
                      LIS_domain(n)%gindex(c,r)) &
                      =gddp_struc(n)%metdata(m,6,&
                      LIS_domain(n)%gindex(c,r))*&
                      psurf_climo1(c,r)/&
                      gddp_struc(n)%psurf_climo(c,r)

                 gddp_struc(n)%metdata1(m,7,&
                      LIS_domain(n)%gindex(c,r)) &
                      = gddp_struc(n)%metdata(m,7,&
                      LIS_domain(n)%gindex(c,r)) *&
                      prcp_climo1(c,r)/&
                      gddp_struc(n)%prcp_climo(c,r)
                 
                 
              elseif(order.eq.2) then 
                 gddp_struc(n)%metdata2(m,1,&
                      LIS_domain(n)%gindex(c,r))&
                      = gddp_struc(n)%metdata(m,1,&
                      LIS_domain(n)%gindex(c,r))*&
                      tair_climo1(c,r)/&
                      gddp_struc(n)%tair_climo(c,r)

                 gddp_struc(n)%metdata2(m,2,&
                      LIS_domain(n)%gindex(c,r))&
                      = gddp_struc(n)%metdata(m,2,&
                      LIS_domain(n)%gindex(c,r))*&
                      qair_climo1(c,r)/&
                      gddp_struc(n)%qair_climo(c,r)
                                  
                 if(gddp_struc(n)%swdown_climo(c,r).ne.0) then 
                    tmpval = &                          
                         gddp_struc(n)%metdata(m,3,&
                         LIS_domain(n)%gindex(c,r))*&
                         swdown_climo1(c,r)/&
                         gddp_struc(n)%swdown_climo(c,r)

                    gddp_struc(n)%metdata2(m,3,&
                         LIS_domain(n)%gindex(c,r)) = tmpval
                 else
                    gddp_struc(n)%metdata2(m,3,&
                         LIS_domain(n)%gindex(c,r)) = 0
                 endif
                 gddp_struc(n)%metdata2(m,4,&
                      LIS_domain(n)%gindex(c,r)) &
                      = gddp_struc(n)%metdata(m,4,&
                      LIS_domain(n)%gindex(c,r))*&
                      lwdown_climo1(c,r)/&
                      gddp_struc(n)%lwdown_climo(c,r)

                 gddp_struc(n)%metdata2(m,5,&
                      LIS_domain(n)%gindex(c,r)) &
                      = gddp_struc(n)%metdata(m,5,&
                      LIS_domain(n)%gindex(c,r))*&
                      wind_climo1(c,r)/&
                      gddp_struc(n)%wind_climo(c,r)

                 gddp_struc(n)%metdata2(m,6,&
                      LIS_domain(n)%gindex(c,r)) &
                      =gddp_struc(n)%metdata(m,6,&
                      LIS_domain(n)%gindex(c,r))*&
                      psurf_climo1(c,r)/&
                      gddp_struc(n)%psurf_climo(c,r)

                 if(gddp_struc(n)%prcp_climo(c,r).ne.0) then 
                    gddp_struc(n)%metdata2(m,7,&
                         LIS_domain(n)%gindex(c,r)) &
                         = gddp_struc(n)%metdata(m,7,&
                         LIS_domain(n)%gindex(c,r)) *&
                         prcp_climo1(c,r)/&
                         gddp_struc(n)%prcp_climo(c,r)
                 else
                    gddp_struc(n)%metdata2(m,7,&
                         LIS_domain(n)%gindex(c,r)) = 0.0
                 endif
                 
              endif
           endif
        end do
     enddo
  enddo
    
#endif

end subroutine read_gddp




!BOP
! !ROUTINE: interp_gddp
! \label{interp_gddp}
!
! !INTERFACE:
subroutine interp_gddp(n,vid, findex, pcp_flag, &
     input_nc, input_nr,input_data,input_bitmap,&
     lis_gds,nc,nr, &
     output_2d)
! !USES:
  use LIS_coreMod
  use LIS_spatialDownscalingMod
  use LIS_logMod
  use gddp_forcingMod, only :gddp_struc
 
  implicit none

! !ARGUMENTS:   
  integer, intent(in)   :: n 
  integer, intent(in)   :: vid
  integer, intent(in)   :: findex
  logical, intent(in)   :: pcp_flag
  integer, intent(in)   :: input_nc
  integer, intent(in)   :: input_nr
  real                  :: input_data(input_nc,input_nr)
  logical*1             :: input_bitmap(input_nc*input_nr)
  real, intent(in)      :: lis_gds(50)
  integer, intent(in)   :: nc
  integer, intent(in)   :: nr
  real, intent(inout)   :: output_2d(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given GDDP field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[kpds]
!  grib decoding array
! \item[input\_size]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[input\_bitmap]
!  input bitmap
! \item[lis\_gds]
!  array description of the LIS grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LIS grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LIS grid
! \item[output\_2d]
!  output interpolated field
!  \end{description} 
! 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using nearest neighbor interpolation
! \end{description}
!EOP
  integer   :: iret
  integer   :: maxrad = 25
  integer   :: c,r,mo
  integer   :: rad, c1,r1,c2,r2
  logical   :: found
  integer   :: i,j
  real      :: input_data_1d(input_nc*input_nr)
  real      :: output_data(nc*nr)
  logical*1 :: output_bitmap(nc*nr)

!=== End variable declarations

  mo = nc*nr
  input_bitmap = .false.
  output_2d = LIS_rc%udef
  
  do r=1,input_nr
     do c=1,input_nc/2
        if(input_data(c,r).eq.1e+20) input_data(c,r) = LIS_rc%udef
        c1 = c+input_nc/2
        input_data_1d(c1+(r-1)*input_nc) = input_data(c,r)
        if(input_data(c,r).ne.LIS_rc%udef) then
           input_bitmap(c1+(r-1)*input_nc) = .true.
        endif
     enddo
  enddo

  do r=1,input_nr
     do c=input_nc/2,input_nc
        if(input_data(c,r).eq.1e+20) input_data(c,r) = LIS_rc%udef
        c1 = c-input_nc/2+1
        input_data_1d(c1+(r-1)*input_nc) = input_data(c,r)
        if(input_data(c,r).ne.LIS_rc%udef) then
           input_bitmap(c1+(r-1)*input_nc) = .true.
        endif
     enddo
  enddo
  
!-----------------------------------------------------------------------
! Initialize output bitmap. 
!-----------------------------------------------------------------------
  output_bitmap = .true.

!-----------------------------------------------------------------------  
! Interpolate to LIS grid
!-----------------------------------------------------------------------
  select case( LIS_rc%met_interp(findex) )

    case( "bilinear" )
     call bilinear_interp(lis_gds,input_bitmap,input_data_1d,&
          output_bitmap,&
          output_data,gddp_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gddp_struc(n)%w111, gddp_struc(n)%w121,&
          gddp_struc(n)%w211,gddp_struc(n)%w221,&
          gddp_struc(n)%n111,gddp_struc(n)%n121,&
          gddp_struc(n)%n211,gddp_struc(n)%n221,LIS_rc%udef,iret)

    case( "budget-bilinear" )
     if (pcp_flag) then     
        call conserv_interp(lis_gds,input_bitmap,input_data_1d,&
             output_bitmap,&
             output_data,gddp_struc(n)%mi,mo,& 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gddp_struc(n)%w112,gddp_struc(n)%w122,&
             gddp_struc(n)%w212,gddp_struc(n)%w222,&
             gddp_struc(n)%n112,gddp_struc(n)%n122,&
             gddp_struc(n)%n212,gddp_struc(n)%n222,LIS_rc%udef,iret)
     else 
        call bilinear_interp(lis_gds,input_bitmap,input_data_1d,&
             output_bitmap,&
             output_data,gddp_struc(n)%mi,mo,&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gddp_struc(n)%w111,gddp_struc(n)%w121,&
             gddp_struc(n)%w211,gddp_struc(n)%w221,&
             gddp_struc(n)%n111,gddp_struc(n)%n121,&
             gddp_struc(n)%n211,gddp_struc(n)%n221,LIS_rc%udef,iret)
     endif
     
  case( "neighbor" )
     call neighbor_interp(lis_gds,input_bitmap,input_data_1d,&
          output_bitmap,&
          output_data,gddp_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gddp_struc(n)%n113,LIS_rc%udef,iret)

  end select

!-----------------------------------------------------------------------    
! convert the interpolated data to 2d. 
!-----------------------------------------------------------------------    

  do j = 1, nr
     do i = 1, nc
        output_2d(i,j) = output_data(i+(j-1)*nc)
     enddo
  enddo

  do r=1,nr
     do c=1,nc
        found = .true. 
        if(LIS_domain(n)%gindex(c,r).ne.-1.and.output_2d(c,r).eq.-9999.0) then
           rad = 1
           found = .false.

           do while(.not.found)
              c1 = max(1,c-rad)
              c2 = min(nc,c+rad)
              r1 = max(1,r-rad)
              r2 = min(nr,r+rad)
              
              do j=r1,r2
                 do i=c1,c2
                    if(output_2d(i,j).ne.-9999.0) then
                       output_2d(c,r) = output_2d(i,j)
                       found = .true.
                       
                       exit
                    endif
                 enddo
                 if(found) exit
              enddo

              if(.not.found) rad = rad +1
              if(rad.gt.maxrad) exit
           enddo
           
        endif
        if(.not.found) then
           write(LIS_logunit,*) '[ERR]: unable to fill in GDDP:', LIS_localPet, vid,  c,r,output_2d(c,r),found
           write(LIS_logunit,*) '[ERR]: ',c+LIS_ews_halo_ind(n,LIS_localPet+1),&
                r+LIS_nss_halo_ind(n,LIS_localPet+1)
           call LIS_endrun()
        endif
     enddo
  enddo
  
   
end subroutine interp_gddp
