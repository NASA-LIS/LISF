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
! !ROUTINE: read_VIIRSgvf
! \label{read_VIIRSgvf}
!
! !REVISION HISTORY:
!  15 Oct 2021    Yonghwan Kwon; initial specification
!
! !INTERFACE: 
subroutine read_VIIRSgvf(n, k, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use map_utils
  use LIS_pluginIndices
  use VIIRSgvf_Mod, only : VIIRSgvf_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the VIIRS green vegetation fraction (GVF) observations 
!  and rescales the data to a reference LAI data.
!  This routine essentially converts the GVF 
!  datasets into the LAI space.  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  real, parameter        ::  minssdev = 0.05
  real, parameter        ::  maxssdev = 0.11
  real,  parameter       :: MAX_LAI_VALUE=10.0, MIN_LAI_VALUE=0.0001
  integer                :: status
  integer                :: grid_index
  character*100          :: gvfobsdir
  character*100          :: fname
  integer                :: vtype
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: smfield, pertField
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  real                   :: obs_unsc(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: gvfobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: lai_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer                :: fnd
  real, allocatable      :: ssdev(:)
  !integer                :: cyr, cmo, cda, chr, cmn, css
  !integer                :: nyr, nmo, nda, nhr, nmn, nss

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       gvfobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false.
  obs_unsc = LIS_rc%udef

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "VIIRS GVF read alarm")

  gvfobs = LIS_rc%udef

  !cyr = LIS_rc%yr
  !cmo = LIS_rc%mo
  !cda = LIS_rc%da
  !chr = LIS_rc%hr
  !cmn = LIS_rc%mn
  !css = LIS_rc%ss

  !call LIS_tick(time1,doy,gmt,cyr,cmo,cda,chr,cmn,css,0.0)
  !nyr = LIS_rc%yr
  !nmo = LIS_rc%mo
  !nda = LIS_rc%da
  !nhr = LIS_rc%hr
  !nmn = LIS_rc%mn
  !nss = LIS_rc%ss

  !call LIS_tick(time2,doy,gmt,nyr,nmo,nda,nhr,nmn,nss,3600.0)
  !nyr = LIS_rc%yr
  !nmo = LIS_rc%mo
  !nda = LIS_rc%da
  !nhr = LIS_rc%hr
  !nmn = LIS_rc%mn
  !nss = LIS_rc%ss

  !call LIS_tick(time3,doy,gmt,nyr,nmo,nda,nhr,nmn,nss,LIS_rc%ts)

  if(alarmCheck.or.VIIRSgvf_struc(n)%startMode) then
     VIIRSgvf_struc(n)%startMode = .false.

     call create_VIIRSgvf_filename(gvfobsdir, &
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname, vtype)

     inquire(file=fname,exist=file_exists)

     if(file_exists) then
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
        call read_VIIRS_GVF_data(n,k, fname, gvfobs, vtype)
     else
        write(LIS_logunit,*) '[WARN] Missing VIIRS GVF file: ',trim(fname)
     endif

     VIIRSgvf_struc(n)%gvfobs  = LIS_rc%udef

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           grid_index = LIS_obs_domain(n,k)%gindex(c,r)
           if(grid_index.ne.-1) then
              if(gvfobs(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then
                 VIIRSgvf_struc(n)%gvfobs(c,r) = &
                         gvfobs(c+(r-1)*LIS_rc%obs_lnc(k))
              endif
           endif
        enddo
     enddo
  endif ! alarm

  if(alarmCheck) then
     call ESMF_StateGet(OBS_State, "Observation01", smfield, &
                        rc=status)
     call LIS_verify(status, 'Error: StateGet Observation01')

     call ESMF_FieldGet(smfield, localDE=0, farrayPtr=obsl, rc=status)
     call LIS_verify(status, 'Error: FieldGet')

     fnd = 0
     lai_current = LIS_rc%udef

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
              lai_current(c,r) = VIIRSgvf_struc(n)%gvfobs(c,r)

              if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                 obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = &
                          lai_current(c,r)
              endif
              if(lai_current(c,r).ne.LIS_rc%udef) then
                 fnd = 1
              endif
           endif
        enddo
     enddo

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------

     ! Read monthly CDF (only for the current month)

     if (VIIRSgvf_struc(n)%ntimes.gt.1.and.VIIRSgvf_struc(n)%cdf_read_opt.eq.1) then
        if (.not. VIIRSgvf_struc(n)%cdf_read_mon .or. LIS_rc%da .eq. 1 .and. LIS_rc%hr .eq. 0 .and. &
             LIS_rc%mn .eq. 0 .and. LIS_rc%ss .eq. 0) then
           call LIS_readMeanSigmaData(n, k, &
                VIIRSgvf_struc(n)%ntimes, &
                LIS_rc%obs_ngrid(k), &
                VIIRSgvf_struc(n)%modelcdffile, &
                "LAI", &
                VIIRSgvf_struc(n)%model_mu, &
                VIIRSgvf_struc(n)%model_sigma, &
                LIS_rc%mo)

           call LIS_readMeanSigmaData(n, k, &
                VIIRSgvf_struc(n)%ntimes, &
                LIS_rc%obs_ngrid(k), &
                VIIRSgvf_struc(n)%obscdffile, &
                "GVF", &
                VIIRSgvf_struc(n)%obs_mu, &
                VIIRSgvf_struc(n)%obs_sigma, &
                LIS_rc%mo)

           call LIS_readCDFdata(n, k, &
                VIIRSgvf_struc(n)%nbins, &
                VIIRSgvf_struc(n)%ntimes, &
                LIS_rc%obs_ngrid(k), &
                VIIRSgvf_struc(n)%modelcdffile, &
                "LAI", &
                VIIRSgvf_struc(n)%model_xrange, &
                VIIRSgvf_struc(n)%model_cdf, &
                LIS_rc%mo)

           call LIS_readCDFdata(n, k, &
                VIIRSgvf_struc(n)%nbins, &
                VIIRSgvf_struc(n)%ntimes, &
                LIS_rc%obs_ngrid(k), &
                VIIRSgvf_struc(n)%obscdffile, &
                "GVF", &
                VIIRSgvf_struc(n)%obs_xrange, &
                VIIRSgvf_struc(n)%obs_cdf, &
                LIS_rc%mo)

           VIIRSgvf_struc(n)%cdf_read_mon = .true.
        endif
     endif

     if (LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then
        if (VIIRSgvf_struc(n)%ntimes.gt.1.and.VIIRSgvf_struc(n)%cdf_read_opt.eq.1) then
           call LIS_rescale_with_CDF_matching( &
                n, k, &
                VIIRSgvf_struc(n)%nbins, &
                1, &
                MAX_LAI_VALUE, &
                MIN_LAI_VALUE, &
                VIIRSgvf_struc(n)%model_xrange, &
                VIIRSgvf_struc(n)%obs_xrange, &
                VIIRSgvf_struc(n)%model_cdf, &
                VIIRSgvf_struc(n)%obs_cdf, &
                lai_current)
        else
           call LIS_rescale_with_CDF_matching( &
                n, k, &
                VIIRSgvf_struc(n)%nbins, &
                VIIRSgvf_struc(n)%ntimes, &
                MAX_LAI_VALUE, &
                MIN_LAI_VALUE, &
                VIIRSgvf_struc(n)%model_xrange, &
                VIIRSgvf_struc(n)%obs_xrange, &
                VIIRSgvf_struc(n)%model_cdf, &
                VIIRSgvf_struc(n)%obs_cdf, &
                lai_current)
        endif
     endif

     obsl = LIS_rc%udef
     do r=1, LIS_rc%obs_lnr(k)
        do c=1, LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=lai_current(c,r)
           endif
        enddo
     enddo
     !-------------------------------------------------------------------------
     !  Apply LSM based QC and screening of observations
     !-------------------------------------------------------------------------
     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
          //trim(LIS_VIIRSgvfobsId)//char(0),n,k,OBS_state)

     call LIS_checkForValidObs(n,k,obsl,fnd,lai_current)

     if(fnd.eq.0) then
        data_upd_flag_local = .false.
     else
        data_upd_flag_local = .true.
     endif

#if (defined SPMD)
     call MPI_ALLGATHER(data_upd_flag_local,1, &
          MPI_LOGICAL, data_upd_flag(:),&
          1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
     data_upd = .false.
     do p=1,LIS_npes
        data_upd = data_upd.or.data_upd_flag(p)
     enddo

     if(data_upd) then

        do t=1,LIS_rc%obs_ngrid(k)
           gid(t) = t
           if(obsl(t).ne.-9999.0) then
              assimflag(t) = 1
           else
              assimflag(t) = 0
           endif
        enddo

        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .true. , rc=status)
        call LIS_verify(status)

        if(LIS_rc%obs_ngrid(k).gt.0) then
           call ESMF_AttributeSet(smField,"Grid Number",&
                gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)

           call ESMF_AttributeSet(smField,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)

           call ESMF_AttributeSet(smfield, "Unscaled Obs",&
                obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
           call LIS_verify(status, 'Error in setting Unscaled Obs attribute')

        endif

        if(LIS_rc%dascaloption(k).eq."CDF matching") then
           if(VIIRSgvf_struc(n)%useSsdevScal.eq.1) then
              call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                   rc=status)
              call LIS_verify(status, 'Error: StateGet Observation01')

              allocate(ssdev(LIS_rc%obs_ngrid(k)))
              ssdev = VIIRSgvf_struc(n)%ssdev_inp

              if(VIIRSgvf_struc(n)%ntimes.eq.1) then
                 jj = 1
              else
                 jj = LIS_rc%mo
              endif
              do t=1,LIS_rc%obs_ngrid(k)
                 if(VIIRSgvf_struc(n)%obs_sigma(t,jj).gt.0) then
                    ssdev(t) = ssdev(t)*VIIRSgvf_struc(n)%model_sigma(t,jj)/&
                         VIIRSgvf_struc(n)%obs_sigma(t,jj)
                    if(ssdev(t).lt.minssdev) then
                       ssdev(t) = minssdev
                    endif
                 endif
              enddo

              if(LIS_rc%obs_ngrid(k).gt.0) then
                 call ESMF_AttributeSet(pertField,"Standard Deviation",&
                      ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                 call LIS_verify(status)
              endif
              deallocate(ssdev)
           endif
        endif
     else
         call ESMF_AttributeSet(OBS_State, "Data Update Status", &
                                .false., rc=status)
         call LIS_verify(status)
     endif
  else   !alarmCheck == F
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
  endif

end subroutine read_VIIRSgvf

!BOP
! 
! !ROUTINE: read_VIIRS_GVF_data
! \label{read_VIIRS_GVF_data}
!
! !INTERFACE:
subroutine read_VIIRS_GVF_data(n,k,fname,gvfobs_ip,vtype)
! 
! !USES:   

  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use VIIRSgvf_Mod, only : VIIRSgvf_struc

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n
  integer                  :: k
  integer                  :: vtype
  character (len=*)        :: fname
  real                     :: gvfobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
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
  integer                 :: gvf_raw(VIIRSgvf_struc(n)%nc,VIIRSgvf_struc(n)%nr)
  real                    :: gvf_in(VIIRSgvf_struc(n)%nc*VIIRSgvf_struc(n)%nr)
  logical*1               :: gvf_data_b(VIIRSgvf_struc(n)%nc*VIIRSgvf_struc(n)%nr)
  logical*1               :: gvfobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                 :: gvfid
  integer                 :: ios, nid
  integer                 :: c,r,c1,r1

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LIS_verify(ios,'Error opening file '//trim(fname))

  if (vtype == 1) then
     ios = nf90_inq_varid(nid, '4km_gvf',gvfid)
     call LIS_verify(ios, 'Error nf90_inq_varid: 4km_gvf')
  else
     ios = nf90_inq_varid(nid, 'gvf_4km',gvfid)
     call LIS_verify(ios, 'Error nf90_inq_varid: gvf_4km')
  endif

  !values

  gvf_data_b = .false.

  ios = nf90_get_var(nid, gvfid, gvf_raw)
  call LIS_verify(ios, 'Error nf90_get_var: gvf')

  ios = nf90_close(ncid=nid)
  call LIS_verify(ios,'Error closing file '//trim(fname))

#endif

  do r=1,VIIRSgvf_struc(n)%nr
     do c=1,VIIRSgvf_struc(n)%nc
        c1 = c
        r1 = nr - r + 1

        if (gvf_raw(c1,r1)>=0.and.&
           gvf_raw(c1,r1)<=100) then

           gvf_in(c+(r-1)*VIIRSgvf_struc(n)%nc) = gvf_raw(c1,r1)*0.01
           gvf_data_b(c+(r-1)*VIIRSgvf_struc(n)%nc) = .true.
        else
           gvf_in(c+(r-1)*VIIRSgvf_struc(n)%nc) = LIS_rc%udef
           gvf_data_b(c+(r-1)*VIIRSgvf_struc(n)%nc) = .false.
        endif
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!--------------------------------------------------------------------------
  if(LIS_rc%obs_gridDesc(k,10).le.0.036) then
     call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
          gvf_data_b, gvf_in, gvfobs_b_ip, gvfobs_ip, &
          VIIRSgvf_struc(n)%nc*VIIRSgvf_struc(n)%nr, &
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          VIIRSgvf_struc(n)%rlat,VIIRSgvf_struc(n)%rlon,&
          VIIRSgvf_struc(n)%w11,VIIRSgvf_struc(n)%w12,&
          VIIRSgvf_struc(n)%w21,VIIRSgvf_struc(n)%w22,&
          VIIRSgvf_struc(n)%n11,VIIRSgvf_struc(n)%n12,&
          VIIRSgvf_struc(n)%n21,VIIRSgvf_struc(n)%n22,LIS_rc%udef,ios)
  else
     call upscaleByAveraging(VIIRSgvf_struc(n)%nc*VIIRSgvf_struc(n)%nr,&
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          LIS_rc%udef, VIIRSgvf_struc(n)%n11,&
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
  use LIS_logMod

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
