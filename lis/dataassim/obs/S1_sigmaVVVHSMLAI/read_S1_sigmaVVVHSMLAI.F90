!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_S1_sigmaVVVHSMLAI
! \label{read_S1_sigmaVVVHSMLAI}
!
! !REVISION HISTORY:
!  29 Aug 2019: Hans Lievens; Initial Specification for snow depth
!  9 Mar 2021:  Isis Brangers, Michel Bechtold; adaptation to backscatter
!  29 Jun 2022:  Louise Busschaert; getting domain dims from files
! !INTERFACE: 
subroutine read_S1_sigmaVVVHSMLAI(n,k,OBS_State,OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use map_utils
  use LIS_DAobservationsMod
  use LIS_pluginIndices, only : LIS_S1_sigmaVVVHSMLAI_obsId
  use S1_sigmaVVVHSMLAI_Mod, only : S1_sigma_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  This routine reads Sentinel-1 backscatter observations in netcdf4 format
!  The data is read at 0z every day and is kept in memory. At 
!  each timestep, a subset of data is chosen for use in DA if 
!  the local time of the grid point is 6AM. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[k]    index of the assimilation instance
!  \item[OBS\_State] observations state
!  \item[OBS\_Pert\_State] observations perturbation state
!  \end{description}
!
!EOP
  type(ESMF_Field)              :: s_vvField,s_vhField,pertfield
  logical                       :: alarmCheck
  logical                       :: data_upd, file_exists
  logical                       :: dataflag(LIS_npes)
  logical                       :: dataflag_local
  integer                       :: c,r, p, t
  character*100                 :: obsdir
  character*80                  :: S1_filename
  integer                       :: ftn
  real                          :: lon, lhour
  real                          :: gmt
  real                          :: dt
  integer                       :: zone
  integer                       :: grid_index
  real                          :: ssdev(LIS_rc%obs_ngrid(k))
  real,             pointer     :: s_vv(:), s_vh(:)
  integer                       :: gid(LIS_rc%obs_ngrid(k))
  integer                       :: assimflag_vv(LIS_rc%obs_ngrid(k))
  integer                       :: assimflag_vh(LIS_rc%obs_ngrid(k))
  integer                       :: status, iret, ierr
  integer                       :: fnd
  real                          :: svv_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                          :: svh_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer                       :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/  !BZ



  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status,'Error in AttributeGet: Data Directory')

!-------------------------------------------------------------------------
!   Read the data at 0z daily. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "S1 backscatter read alarm")

  if(alarmCheck.or.S1_sigma_struc(n)%startMode) then 
     S1_sigma_struc(n)%startMode = .false.
     
     S1_sigma_struc(n)%s_vv = LIS_rc%udef
     S1_sigma_struc(n)%s_vh = LIS_rc%udef
     S1_sigma_struc(n)%sigmatime = -1

     call S1_sigmaVVVHSMLAI_filename(S1_filename,obsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)       

     inquire(file=S1_filename,exist=file_exists)
     if(file_exists) then 

        write(LIS_logunit,*)  '[INFO] Reading S1 sigma data ',S1_filename
        call read_S1_sigmaVVVHSMLAI_data(n,k, S1_filename, S1_sigma_struc(n)%s_vv, &
             S1_sigma_struc(n)%s_vh)


!-------------------------------------------------------------------------
! Store the GMT corresponding to 6AM localtime at each grid point
!-------------------------------------------------------------------------
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                 grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
                 lon = LIS_obs_domain(n,k)%lon(grid_index)
                 
                 lhour = 6.0
                 call LIS_localtime2gmt(gmt,lon,lhour,zone)
                 S1_sigma_struc(n)%sigmatime(c,r) = gmt

              endif
           enddo
        enddo

     endif
  endif

!-------------------------------------------------------------------------
! Update the OBS_State
!-------------------------------------------------------------------------

  call ESMF_StateGet(OBS_State,"Observation01",s_vvField,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')

  call ESMF_StateGet(OBS_State,"Observation02",s_vhField,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation02')
  
  call ESMF_FieldGet(s_vvField,localDE=0,farrayPtr=s_vv,rc=status)
  call LIS_verify(status, 'Error: FieldGet')

  call ESMF_FieldGet(s_vhField,localDE=0,farrayPtr=s_vh,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
  s_vv = LIS_rc%udef 
  s_vh = LIS_rc%udef 

!-------------------------------------------------------------------------
!  Update the OBS_State by subsetting to the local grid time  
!-------------------------------------------------------------------------     

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
           
           dt = (LIS_rc%gmt - S1_sigma_struc(n)%sigmatime(c,r))*3600.0
           lon = LIS_obs_domain(n,k)%lon(grid_index)

           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
              s_vv(LIS_obs_domain(n,k)%gindex(c,r)) = & 
                   S1_sigma_struc(n)%s_vv(c,r)
              s_vh(LIS_obs_domain(n,k)%gindex(c,r)) = & 
                   S1_sigma_struc(n)%s_vh(c,r)
           endif
           
        endif
     enddo
  enddo

  dataflag_local = .false. 

!-------------------------------------------------------------------------
!  Apply LSM based quality control and screening of observations
!-------------------------------------------------------------------------     

  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_S1_sigmaVVVHSMLAI_obsId)//char(0), & 
       n, k,OBS_state)

  svv_current = LIS_rc%udef
  svh_current = LIS_rc%udef
  call LIS_checkForValidObs(n,k,s_vv,fnd,svv_current) 
  call LIS_checkForValidObs(n,k,s_vh,fnd,svh_current) 

  if(fnd.eq.0) then 
     dataflag_local = .false. 
  else
     dataflag_local = .true. 
  endif
 
#if (defined SPMD)
  call MPI_ALLGATHER(dataflag_local,1, MPI_LOGICAL, dataflag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, ierr)
#endif
  data_upd = .false.
  
  do p=1,LIS_npes
     data_upd = data_upd.or.dataflag(p)
  enddo

  if(data_upd) then 
     do t=1,LIS_rc%obs_ngrid(k)
        gid(t) = t
        if(s_vv(t).ne.-9999.0) then 
           assimflag_vv(t) = 1
        else
           assimflag_vv(t) = 0
        endif
     enddo
     do t=1,LIS_rc%obs_ngrid(k)
        gid(t) = t
        if(s_vh(t).ne.-9999.0) then 
           assimflag_vh(t) = 1
        else
           assimflag_vh(t) = 0
        endif
     enddo
     
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status, 'Error: AttributeSet in Data Update Status')
     

     call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
          rc=status)
     call LIS_verify(status, 'ESMF_StateGet for Observation01 for OBS_Pert_State failed in read_S1_sigma')

     call ESMF_StateGet(OBS_Pert_State,"Observation02",pertfield,&
          rc=status)
     call LIS_verify(status, 'ESMF_StateGet for Observation02 for OBS_Pert_State failed in read_S1_sigma')
     
     if(LIS_rc%obs_ngrid(k).gt.0) then 

        !linearly scale the observation err
        ssdev = S1_sigma_struc(n)%ssdev_inp
        do t=1,LIS_rc%obs_ngrid(k)
           if(s_vv(t).ne.-9999.0) then 
              ssdev(t) =  S1_sigma_struc(n)%ssdev_inp !+ 0.05*obsl(t)
           endif
        enddo

        call ESMF_AttributeSet(pertField,"Standard Deviation",&
             ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(s_vvField,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status,'Error: AttributeSet in Grid Number')
        
        call ESMF_AttributeSet(s_vvField,"Assimilation Flag",&
             assimflag_vv,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')

        call ESMF_AttributeSet(s_vhField,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status,'Error: AttributeSet in Grid Number')
        
        call ESMF_AttributeSet(s_vhField,"Assimilation Flag",&
             assimflag_vh,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')
        
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status, "Error: AttributeSet Data Update Status")
     return
  end if
 

end subroutine read_S1_sigmaVVVHSMLAI



!BOP
!
! !ROUTINE: read_S1_sigmaVVVHSMLAI_data
! \label{read_S1_sigmaVVVHSMLAI_data}
!
! !INTERFACE:
subroutine read_S1_sigmaVVVHSMLAI_data(n, k, fname, svv_ip, svh_ip)
!
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod
  use map_utils,    only : latlon_to_ij
  use S1_sigmaVVVHSMLAI_Mod, only : S1_sigma_struc

  implicit none
!
! !INPUT PARAMETERS:
!
  integer                       :: n
  integer                       :: k
  character (len=*)             :: fname

! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine reads the S1 NETCDF files
!
!  The arguments are:
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the S1 sigma file
!  \item[svvobs\_ip]   sigma depth data processed to the LIS domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
  real                        :: s_vv(S1_sigma_struc(n)%nr,S1_sigma_struc(n)%nc)
  real                        :: s_vh(S1_sigma_struc(n)%nr,S1_sigma_struc(n)%nc)
  real                        :: lat_nc(S1_sigma_struc(n)%nr)
  real                        :: lat_nc_fold(S1_sigma_struc(n)%nr)
  real                        :: lon_nc(S1_sigma_struc(n)%nc)
  real                        :: svv_ip(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                        :: svh_ip(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer                     :: ns_ip(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  logical                     :: file_exists
  integer                     :: c,r,i,j
  integer                     :: stn_col,stn_row
  real                        :: col,row
  integer                     :: nid
  integer                     :: svvId,svhId,latId,lonId
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

!  s_vv = LIS_rc%udef
!  lat_nc = LIS_rc%udef
!  lon_nc = LIS_rc%udef
!  svv_ip = LIS_rc%udef

  inquire(file=fname, exist=file_exists)
  if(file_exists) then
     write(LIS_logunit,*) 'Reading ',trim(fname)
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(fname))
 
     ! variables
     ios = nf90_inq_varid(nid, 'g0vv',svvid)
     if(ios.lt.0) then
       ios = nf90_inq_varid(nid, 's0vv',svvid)
     endif
     call LIS_verify(ios, 'Error nf90_inq_varid: backscatter data')

     ios = nf90_inq_varid(nid, 'g0vh',svhid)
     if(ios.lt.0) then
       ios = nf90_inq_varid(nid, 's0vh',svhid)
     endif
     call LIS_verify(ios, 'Error nf90_inq_varid: backscatter data')

     ios = nf90_inq_varid(nid, 'lat',latid)
     call LIS_verify(ios, 'Error nf90_inq_varid: latitude data')

     ios = nf90_inq_varid(nid, 'lon',lonid)
     call LIS_verify(ios, 'Error nf90_inq_varid: longitude data')

     !values
     ios = nf90_get_var(nid, svvid, s_vv)
     call LIS_verify(ios, 'Error nf90_get_var: s0vv')

     ios = nf90_get_var(nid, svhid, s_vh)
     call LIS_verify(ios, 'Error nf90_get_var: s0vh')

     ios = nf90_get_var(nid, latid, lat_nc_fold)
     call LIS_verify(ios, 'Error nf90_get_var: lat')
     !new g0 dataset of Hans with flipped lat variable
     do i=1,size(lat_nc_fold)
       lat_nc(i)=lat_nc_fold(size(lat_nc_fold)+1-i)
     end do 
     ios = nf90_get_var(nid, lonid, lon_nc)
     call LIS_verify(ios, 'Error nf90_get_var: lon')

     ! close file
     ios = nf90_close(ncid=nid)
     call LIS_verify(ios,'Error closing file '//trim(fname))



!    ! Initialize 
     svv_ip = 0
     svh_ip = 0
     ns_ip = 0

     ! Interpolate the data by averaging 
     do i=1,S1_sigma_struc(n)%nr
        do j=1,S1_sigma_struc(n)%nc
              
           call latlon_to_ij(LIS_domain(n)%lisproj,&
                !lat_nc(i),lon_nc(j),col,row)
                lat_nc(S1_sigma_struc(n)%nr-(i-1)),lon_nc(j),col,row)
           stn_col = nint(col)
           stn_row = nint(row)

           if(s_vv(i,j).ge.-999.and.&
                stn_col.gt.0.and.stn_col.le.LIS_rc%obs_lnc(k).and.&
                stn_row.gt.0.and.stn_row.le.LIS_rc%obs_lnr(k)) then
              ! add in linear scale for later averaging
                svv_ip(stn_col,stn_row) = svv_ip(stn_col,stn_row) + 10.**(s_vv(i,j)/10.)
              svh_ip(stn_col,stn_row) = svh_ip(stn_col,stn_row) + 10.**(s_vh(i,j)/10.)
              ns_ip(stn_col,stn_row) = ns_ip(stn_col,stn_row) + 1
           endif
        enddo
     enddo

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if(ns_ip(c,r).ne.0.and.svv_ip(c,r).gt.0.0.and.svh_ip(c,r).gt.0.0) then
              ! average in linear scale
              svv_ip(c,r) = svv_ip(c,r)/ns_ip(c,r)
              svh_ip(c,r) = svh_ip(c,r)/ns_ip(c,r)
           else
              svv_ip(c,r) = LIS_rc%udef
              svh_ip(c,r) = LIS_rc%udef
           endif
        enddo
     enddo


     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if(svv_ip(c,r).ne.LIS_rc%udef) then
              ! back to log scale, dB
              svv_ip(c,r) = 10.*LOG10(svv_ip(c,r));
            endif
           if(svh_ip(c,r).ne.LIS_rc%udef) then
              ! back to log scale, dB
              svh_ip(c,r) = 10.*LOG10(svh_ip(c,r));
           endif
        enddo
     enddo

  endif

#endif

end subroutine read_S1_sigmaVVVHSMLAI_data



!BOP
!
! !ROUTINE: S1_sigmaVVVHSMLAI_filename
! \label{S1_SWND_filename}
! 
! !INTERFACE: 
subroutine S1_sigmaVVVHSMLAI_filename(filename, ndir, yr, mo, da)
  
  implicit none
! !ARGUMENTS: 
  character*80      :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped S1 filename
!  
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the S1 filename
!  \item[ndir] name of the S1 root directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da 
  filename = trim(ndir)//'/S1_g0_'//trim(fyr)//trim(fmo)//trim(fda)//'.nc'
    
end subroutine S1_sigmaVVVHSMLAI_filename



