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
! !ROUTINE: read_CDFSgvf
! \label{read_CDFSgvf}
!
! !REVISION HISTORY:
!  04 Apr 2022    Yonghwan Kwon; initial specification
!
! !INTERFACE: 
subroutine read_CDFSgvf(n, k, OBS_State, OBS_Pert_State)
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
  use CDFSgvf_Mod, only : CDFSgvf_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the CDFS green vegetation fraction (GVF) observations 
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

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "CDFS GVF read alarm")

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

  if(alarmCheck.or.CDFSgvf_struc(n)%startMode) then
     CDFSgvf_struc(n)%startMode = .false.

     call create_CDFSgvf_filename(gvfobsdir, &
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname)

     inquire(file=fname,exist=file_exists)

     if(file_exists) then
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
        call read_CDFS_GVF_data(n,k, fname, gvfobs)
     else
        write(LIS_logunit,*) '[WARN] Missing CDFS GVF file: ',trim(fname)
     endif

     CDFSgvf_struc(n)%gvfobs  = LIS_rc%udef

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           grid_index = LIS_obs_domain(n,k)%gindex(c,r)
           if(grid_index.ne.-1) then
              if(gvfobs(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then
                 CDFSgvf_struc(n)%gvfobs(c,r) = &
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
              lai_current(c,r) = CDFSgvf_struc(n)%gvfobs(c,r)

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

     if (CDFSgvf_struc(n)%ntimes.gt.1.and.CDFSgvf_struc(n)%cdf_read_opt.eq.1) then
        if (.not. CDFSgvf_struc(n)%cdf_read_mon .or. LIS_rc%da .eq. 1 .and. LIS_rc%hr .eq. 0 .and. &
             LIS_rc%mn .eq. 0 .and. LIS_rc%ss .eq. 0) then
           call LIS_readMeanSigmaData(n, k, &
                CDFSgvf_struc(n)%ntimes, &
                LIS_rc%obs_ngrid(k), &
                CDFSgvf_struc(n)%modelcdffile, &
                "LAI", &
                CDFSgvf_struc(n)%model_mu, &
                CDFSgvf_struc(n)%model_sigma, &
                LIS_rc%mo)

           call LIS_readMeanSigmaData(n, k, &
                CDFSgvf_struc(n)%ntimes, &
                LIS_rc%obs_ngrid(k), &
                CDFSgvf_struc(n)%obscdffile, &
                "GVF", &
                CDFSgvf_struc(n)%obs_mu, &
                CDFSgvf_struc(n)%obs_sigma, &
                LIS_rc%mo)

           call LIS_readCDFdata(n, k, &
                CDFSgvf_struc(n)%nbins, &
                CDFSgvf_struc(n)%ntimes, &
                LIS_rc%obs_ngrid(k), &
                CDFSgvf_struc(n)%modelcdffile, &
                "LAI", &
                CDFSgvf_struc(n)%model_xrange, &
                CDFSgvf_struc(n)%model_cdf, &
                LIS_rc%mo)

           call LIS_readCDFdata(n, k, &
                CDFSgvf_struc(n)%nbins, &
                CDFSgvf_struc(n)%ntimes, &
                LIS_rc%obs_ngrid(k), &
                CDFSgvf_struc(n)%obscdffile, &
                "GVF", &
                CDFSgvf_struc(n)%obs_xrange, &
                CDFSgvf_struc(n)%obs_cdf, &
                LIS_rc%mo)

           CDFSgvf_struc(n)%cdf_read_mon = .true.
        endif
     endif

     if (LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then
        if (CDFSgvf_struc(n)%ntimes.gt.1.and.CDFSgvf_struc(n)%cdf_read_opt.eq.1) then
           call LIS_rescale_with_CDF_matching( &
                n, k, &
                CDFSgvf_struc(n)%nbins, &
                1, &
                MAX_LAI_VALUE, &
                MIN_LAI_VALUE, &
                CDFSgvf_struc(n)%model_xrange, &
                CDFSgvf_struc(n)%obs_xrange, &
                CDFSgvf_struc(n)%model_cdf, &
                CDFSgvf_struc(n)%obs_cdf, &
                lai_current)
        else
           call LIS_rescale_with_CDF_matching( &
                n, k, &
                CDFSgvf_struc(n)%nbins, &
                CDFSgvf_struc(n)%ntimes, &
                MAX_LAI_VALUE, &
                MIN_LAI_VALUE, &
                CDFSgvf_struc(n)%model_xrange, &
                CDFSgvf_struc(n)%obs_xrange, &
                CDFSgvf_struc(n)%model_cdf, &
                CDFSgvf_struc(n)%obs_cdf, &
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
          //trim(LIS_CDFSgvfobsId)//char(0),n,k,OBS_state)

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
           if(CDFSgvf_struc(n)%useSsdevScal.eq.1) then
              call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                   rc=status)
              call LIS_verify(status, 'Error: StateGet Observation01')

              allocate(ssdev(LIS_rc%obs_ngrid(k)))
              ssdev = CDFSgvf_struc(n)%ssdev_inp

              if(CDFSgvf_struc(n)%ntimes.eq.1) then
                 jj = 1
              else
                 jj = LIS_rc%mo
              endif
              do t=1,LIS_rc%obs_ngrid(k)
                 if(CDFSgvf_struc(n)%obs_sigma(t,jj).gt.0) then
                    ssdev(t) = ssdev(t)*CDFSgvf_struc(n)%model_sigma(t,jj)/&
                         CDFSgvf_struc(n)%obs_sigma(t,jj)
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

end subroutine read_CDFSgvf

!BOP
! 
! !ROUTINE: read_CDFS_GVF_data
! \label{read_CDFS_GVF_data}
!
! !INTERFACE:
subroutine read_CDFS_GVF_data(n,k,fname,gvfobs_ip)
! 
! !USES:   

  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use CDFSgvf_Mod, only : CDFSgvf_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n
  integer                  :: k
  character (len=*)        :: fname
  real                     :: gvfobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
!
! !DESCRIPTION:
!  This subroutine reads the CDFS GVF file 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the CDFS GVF file
!  \item[gvtobs\_ip]   GVF data processed to the LIS domain
! \end{description}
!
!
!EOP

! !USES:   
  integer,  parameter     :: nc=7200, nr=3600
  real*4                  :: gvf_raw(CDFSgvf_struc(n)%nc,CDFSgvf_struc(n)%nr)
  real                    :: gvf_in(CDFSgvf_struc(n)%nc*CDFSgvf_struc(n)%nr)
  logical*1               :: gvf_data_b(CDFSgvf_struc(n)%nc*CDFSgvf_struc(n)%nr)
  logical*1               :: gvfobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                 :: gvfid
  integer                 :: ios, nid
  integer                 :: c,r
  integer                 :: ftn1

  ftn1 = LIS_getNextUnitNumber()
  open(unit=ftn1,file=fname,form='unformatted',access='direct',convert='little_endian',recl=4*nc*nr,status='old')
  read(ftn1, rec=1) gvf_raw
  close(1)
  call LIS_releaseUnitNumber(ftn1)

  gvf_data_b = .false.

  do r=1,CDFSgvf_struc(n)%nr
     do c=1,CDFSgvf_struc(n)%nc
        if (gvf_raw(c,r)>=0.and.&
           gvf_raw(c,r)<=100) then

           gvf_in(c+(r-1)*CDFSgvf_struc(n)%nc) = gvf_raw(c,r)
           gvf_data_b(c+(r-1)*CDFSgvf_struc(n)%nc) = .true.
        else
           gvf_in(c+(r-1)*CDFSgvf_struc(n)%nc) = LIS_rc%udef
           gvf_data_b(c+(r-1)*CDFSgvf_struc(n)%nc) = .false.
        endif
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!--------------------------------------------------------------------------
  if(LIS_rc%obs_gridDesc(k,10).le.0.05) then
     call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
          gvf_data_b, gvf_in, gvfobs_b_ip, gvfobs_ip, &
          CDFSgvf_struc(n)%nc*CDFSgvf_struc(n)%nr, &
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          CDFSgvf_struc(n)%rlat,CDFSgvf_struc(n)%rlon,&
          CDFSgvf_struc(n)%w11,CDFSgvf_struc(n)%w12,&
          CDFSgvf_struc(n)%w21,CDFSgvf_struc(n)%w22,&
          CDFSgvf_struc(n)%n11,CDFSgvf_struc(n)%n12,&
          CDFSgvf_struc(n)%n21,CDFSgvf_struc(n)%n22,LIS_rc%udef,ios)
  else
     call upscaleByAveraging(CDFSgvf_struc(n)%nc*CDFSgvf_struc(n)%nr,&
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          LIS_rc%udef, CDFSgvf_struc(n)%n11,&
          gvf_data_b,gvf_in, gvfobs_b_ip, gvfobs_ip)

  endif

end subroutine read_CDFS_GVF_data

!BOP
! !ROUTINE: create_CDFSgvf_filename
! \label{create_CDFSgvf_filename}
! 
! !INTERFACE:
subroutine create_CDFSgvf_filename(ndir,yr,mo,da,filename)
! !USES:
  use LIS_logMod

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the CDFS GVF filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the CDFS GVF data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated CDFS GVF filename
! \end{description}
!EOP

  character*4             :: yyyy
  character*2             :: mm,dd

  write(unit=yyyy, fmt='(i4.4)') yr
  write(unit=mm, fmt='(i2.2)') mo
  write(unit=dd, fmt='(i2.2)') da

  filename = trim(ndir)//'/'//trim(yyyy)//'/green.'//&
             trim(yyyy)//trim(mm)//trim(dd)//'.1gd4r'

end subroutine create_CDFSgvf_filename
