!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readVIIRS_GVFObs
! \label{readVIIRS_GVFObs}
! 
! !REVISION HISTORY: 
!  8 Oct 2021: Yonghwan Kwon, Initial Specification
! 
! !INTERFACE: 
subroutine readVIIRS_GVFObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_DAobsDataMod
  use VIIRSGVFobsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for 
! VIIRS Green Vegetation Fraction (GVF) data
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r
  character*100     :: fname
  integer           :: vtype
  real              :: gvfobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  VIIRSgvfobs(n)%gvfobs = LDT_rc%udef
  gvfobs= LDT_rc%udef

  call create_VIIRSgvf_filename(VIIRSgvfobs(n)%odir, &
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da,&
       fname,vtype)

  inquire(file=trim(fname),exist=file_exists)
  if(file_exists) then

     write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
     call read_VIIRS_GVF_data(n, fname, gvfobs, vtype)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(gvfobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then
              VIIRSgvfobs(n)%gvfobs(c,r) = gvfobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%gvf_obs,&
       VIIRSgvfobs(n)%gvfobs,vlevel=1)

end subroutine readVIIRS_GVFObs

!BOP
! 
! !ROUTINE: read_VIIRS_GVF_data
! \label{read_VIIRS_GVF_data}
!
! !INTERFACE:
subroutine read_VIIRS_GVF_data(n, fname, gvfobs_ip, vtype)
! 
! !USES:

  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use VIIRSGVFobsMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n
  integer                       :: vtype
  character (len=*)             :: fname
  real                          :: gvfobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine reads the VIIRS GVF file 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the VIIRS GVF file
!  \item[gvtobs\_ip]   GVF data processed to the LIS domain
! \end{description}
!
!
!EOP

! !USES:   
  integer,  parameter     :: nc=10000, nr=5000
  integer                 :: gvf_raw(VIIRSgvfobs(n)%nc,VIIRSgvfobs(n)%nr)
  real                    :: gvf_in(VIIRSgvfobs(n)%nc*VIIRSgvfobs(n)%nr)
  logical*1               :: gvf_data_b(VIIRSgvfobs(n)%nc*VIIRSgvfobs(n)%nr)
  logical*1               :: gvfobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer                 :: gvfid
  integer                 :: ios, nid
  integer                 :: c,r,c1,r1

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LDT_verify(ios,'Error opening file '//trim(fname))

  if (vtype == 1) then
     ios = nf90_inq_varid(nid, '4km_gvf',gvfid)
     call LDT_verify(ios, 'Error nf90_inq_varid: 4km_gvf')
  else
     ios = nf90_inq_varid(nid, 'gvf_4km',gvfid)
     call LDT_verify(ios, 'Error nf90_inq_varid: gvf_4km')
  endif

  !values

  gvf_data_b = .false.

  ios = nf90_get_var(nid, gvfid, gvf_raw)
  call LDT_verify(ios, 'Error nf90_get_var: gvf')

  ios = nf90_close(ncid=nid)
  call LDT_verify(ios,'Error closing file '//trim(fname))

#endif

  do r=1,VIIRSgvfobs(n)%nr
     do c=1,VIIRSgvfobs(n)%nc
        c1 = c
        r1 = nr - r + 1 

        if (gvf_raw(c1,r1)>=0.and.&
           gvf_raw(c1,r1)<=100) then

           gvf_in(c+(r-1)*VIIRSgvfobs(n)%nc) = gvf_raw(c1,r1)*0.01
           gvf_data_b(c+(r-1)*VIIRSgvfobs(n)%nc) = .true.
        else
           gvf_in(c+(r-1)*VIIRSgvfobs(n)%nc) = LDT_rc%udef
           gvf_data_b(c+(r-1)*VIIRSgvfobs(n)%nc) = .false.
        endif 
     enddo
  enddo 

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  if(LDT_isLDTatAfinerResolution(n,0.036)) then
     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          gvf_data_b, gvf_in, gvfobs_b_ip, gvfobs_ip, &
          VIIRSgvfobs(n)%nc*VIIRSgvfobs(n)%nr, &
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          VIIRSgvfobs(n)%w11,VIIRSgvfobs(n)%w12,&
          VIIRSgvfobs(n)%w21,VIIRSgvfobs(n)%w22,&
          VIIRSgvfobs(n)%n11,VIIRSgvfobs(n)%n12,&
          VIIRSgvfobs(n)%n21,VIIRSgvfobs(n)%n22,LDT_rc%udef,ios)
  else
     call upscaleByAveraging(VIIRSgvfobs(n)%nc*VIIRSgvfobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_rc%udef, VIIRSgvfobs(n)%n11,&
          gvf_data_b,gvf_in, gvfobs_b_ip, gvfobs_ip)
  endif

end subroutine read_VIIRS_GVF_data


!BOP
! !ROUTINE: create_VIIRSgvf_filename
! \label{create_VIIRSgvf_filename}
! 
! !INTERFACE:
subroutine create_VIIRSgvf_filename(ndir,yr,mo,da,filename,vtype)
! !USES:

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
  integer           :: vtype
! 
! !DESCRIPTION: 
!  This subroutine creates the VIIRS GVF filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the VIIRS GVF data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated VIIRS GVF filename
! \end{description}
!EOP

  character*4             :: yyyy
  character*2             :: mm,dd

  write(unit=yyyy, fmt='(i4.4)') yr
  write(unit=mm, fmt='(i2.2)') mo
  write(unit=dd, fmt='(i2.2)') da
  
!v1r0: 2015, 2016, 2017, to 2018-09-26
!v2r1: from 2018-09-27, to 2019-09-17, 
!v2r3: from 2019-09-18 to 2021
!/2015/GVF-WKL-GLB_v1r0_npp_e20151231.nc
!/2018/GVF-WKL-GLB_v2r1_npp_e20181231.nc
!/2020/GVF-WKL-GLB_v2r3_npp_e20201230.nc

  if (yr<=2017) then
     filename = trim(ndir)//'/'//trim(yyyy)//'/GVF-WKL-GLB_v1r0_npp_e'//&
               trim(yyyy)//trim(mm)//trim(dd)//'.nc'
     vtype = 1;
  elseif (yr==2018.and.mo<=9.and.da<=26) then
     filename = trim(ndir)//'/'//trim(yyyy)//'/GVF-WKL-GLB_v1r0_npp_e'//&
               trim(yyyy)//trim(mm)//trim(dd)//'.nc'
     vtype = 1;
  elseif (yr==2018.and.mo>=9.and.da>=27) then
     filename = trim(ndir)//'/'//trim(yyyy)//'/GVF-WKL-GLB_v2r1_npp_e'//&
               trim(yyyy)//trim(mm)//trim(dd)//'.nc'
     vtype = 2;
  elseif (yr==2019.and.mo<=9.and.da<=17) then
     filename = trim(ndir)//'/'//trim(yyyy)//'/GVF-WKL-GLB_v2r1_npp_e'//&
               trim(yyyy)//trim(mm)//trim(dd)//'.nc'
     vtype = 2;
  elseif (yr==2019.and.mo>=9.and.da>=18) then
     filename = trim(ndir)//'/'//trim(yyyy)//'/GVF-WKL-GLB_v2r3_npp_e'//&
               trim(yyyy)//trim(mm)//trim(dd)//'.nc'
     vtype = 2;
  else !yr >= 2020
     filename = trim(ndir)//'/'//trim(yyyy)//'/GVF-WKL-GLB_v2r3_npp_e'//&
               trim(yyyy)//trim(mm)//trim(dd)//'.nc'
     vtype = 2;
  endif

end subroutine create_VIIRSgvf_filename





