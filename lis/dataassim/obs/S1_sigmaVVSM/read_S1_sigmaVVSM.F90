!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_S1_sigmaVVSM
! \label{read_S1_sigmaVVSM}
!
! !REVISION HISTORY:
!  29 Aug 2019: Hans Lievens; Initial Specification for snow depth
!  9 Mar 2021:  Isis Brangers, Michel Bechtold; adaptation to backscatter
! 29 Jun 2022:  Louise Busschaert; getting domain dims from files
!
! !INTERFACE: 
subroutine read_S1_sigmaVVSM(n,k,OBS_State,OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use map_utils
  use LIS_DAobservationsMod
  use LIS_pluginIndices, only : LIS_S1_sigmaVVSM_obsId
  use S1_sigmaVVSM_Mod, only : S1_sigma_struc

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
  type(ESMF_Field)              :: sigmaField,pertfield
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
  real,             pointer     :: obsl(:)
  integer                       :: gid(LIS_rc%obs_ngrid(k))
  integer                       :: assimflag(LIS_rc%obs_ngrid(k))
  integer                       :: status, iret, ierr
  integer                       :: fnd
  real                          :: sigma_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
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
     
     S1_sigma_struc(n)%sigma = LIS_rc%udef
     S1_sigma_struc(n)%sigmatime = -1

     call S1_sigmaVVSM_filename(S1_filename,obsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)       

     inquire(file=S1_filename,exist=file_exists)
     if(file_exists) then 

        write(LIS_logunit,*)  '[INFO] Reading S1 sigma data ',S1_filename
        call read_S1_sigmaVVSM_data(n,k, S1_filename, S1_sigma_struc(n)%sigma)


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

  call ESMF_StateGet(OBS_State,"Observation01",sigmafield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(sigmafield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
  obsl = LIS_rc%udef 

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
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = & 
                   S1_sigma_struc(n)%sigma(c,r)
           endif
           
        endif
     enddo
  enddo

  dataflag_local = .false. 

!-------------------------------------------------------------------------
!  Apply LSM based quality control and screening of observations
!-------------------------------------------------------------------------     

  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_S1_sigmaVVSM_obsId)//char(0), & 
       n, k,OBS_state)

  sigma_current = LIS_rc%udef
  call LIS_checkForValidObs(n,k,obsl,fnd,sigma_current)

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
        if(obsl(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo
     
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status, 'Error: AttributeSet in Data Update Status')
     

     call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
          rc=status)
     call LIS_verify(status, 'ESMF_StateGet for Observation01 for OBS_Pert_State failed in read_S1_sigma')
     
     if(LIS_rc%obs_ngrid(k).gt.0) then 

!linearly scale the observation err
        ssdev = S1_sigma_struc(n)%ssdev 
        do t=1,LIS_rc%obs_ngrid(k)
           if(obsl(t).ne.-9999.0) then 
              ssdev(t) =  S1_sigma_struc(n)%ssdev !+ 0.05*obsl(t)
           endif
        enddo

        call ESMF_AttributeSet(pertField,"Standard Deviation",&
             ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(sigmafield,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status,'Error: AttributeSet in Grid Number')
        
        call ESMF_AttributeSet(sigmafield,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')
        
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status, "Error: AttributeSet Data Update Status")
     return
  end if
 

end subroutine read_S1_sigmaVVSM



!BOP
!
! !ROUTINE: read_S1_sigmaVVSM_data
! \label{read_S1_sigmaVVSM_data}
!
! !INTERFACE:
subroutine read_S1_sigmaVVSM_data(n, k, fname, sigma_ip)
!
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod
  use map_utils,    only : latlon_to_ij
  use S1_sigmaVVSM_Mod, only : S1_sigma_struc

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
!  \item[s0vvobs\_ip]   sigma depth data processed to the LIS domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
  real                        :: sigma(S1_sigma_struc(n)%nr,S1_sigma_struc(n)%nc)
  real                        :: lat_nc(S1_sigma_struc(n)%nr)
  real                        :: lat_nc_fold(S1_sigma_struc(n)%nr)
  real                        :: lon_nc(S1_sigma_struc(n)%nc)
  real                        :: sigma_ip(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                        :: sigma_fill(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer                     :: nsigma_ip(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  logical                     :: file_exists
  integer                     :: c,r,i,j
  integer                     :: stn_col,stn_row
  real                        :: col,row
  integer                     :: nid
  integer                     :: sigmaId,latId,lonId
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

!  sigma = LIS_rc%udef
!  lat_nc = LIS_rc%udef
!  lon_nc = LIS_rc%udef
!  sigma_ip = LIS_rc%udef

  inquire(file=fname, exist=file_exists)
  if(file_exists) then
     write(LIS_logunit,*) 'Reading ',trim(fname)
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(fname))
 
     ! variables
     ios = nf90_inq_varid(nid, 'g0vv',sigmaid)
     if(ios.lt.0) then
       ios = nf90_inq_varid(nid, 's0vv',sigmaid)
     endif
     call LIS_verify(ios, 'Error nf90_inq_varid: backscatter data')

     ios = nf90_inq_varid(nid, 'lat',latid)
     call LIS_verify(ios, 'Error nf90_inq_varid: latitude data')

     ios = nf90_inq_varid(nid, 'lon',lonid)
     call LIS_verify(ios, 'Error nf90_inq_varid: longitude data')

     !values
     ios = nf90_get_var(nid, sigmaid, sigma)
     call LIS_verify(ios, 'Error nf90_get_var: s0vv')

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
     sigma_ip = 0
     nsigma_ip = 0

     ! Interpolate the data by averaging 
     do i=1,S1_sigma_struc(n)%nr
        do j=1,S1_sigma_struc(n)%nc
              
           call latlon_to_ij(LIS_domain(n)%lisproj,&
                !lat_nc(i),lon_nc(j),col,row)
                lat_nc(S1_sigma_struc(n)%nr-(i-1)),lon_nc(j),col,row)
           stn_col = nint(col)
           stn_row = nint(row)

           if(sigma(i,j).ge.-999.and.&
                stn_col.gt.0.and.stn_col.le.LIS_rc%obs_lnc(k).and.&
                stn_row.gt.0.and.stn_row.le.LIS_rc%obs_lnr(k)) then
              ! add in linear scale for later averaging
              sigma_ip(stn_col,stn_row) = sigma_ip(stn_col,stn_row) + 10.**(sigma(i,j)/10.)
              nsigma_ip(stn_col,stn_row) = nsigma_ip(stn_col,stn_row) + 1
           endif
        enddo
     enddo

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if(nsigma_ip(c,r).ne.0) then
              ! average in linear scale
              sigma_ip(c,r) = sigma_ip(c,r)/nsigma_ip(c,r)
           else
              sigma_ip(c,r) = LIS_rc%udef
           endif
        enddo
     enddo


     ! Fix stripes in obs caused by regridding (missing obs for LIS grid cell)
     sigma_fill = 0
     nsigma_ip = 0
     do r=2,LIS_rc%obs_lnr(k)-1
        do c=2,LIS_rc%obs_lnc(k)-1
           if (sigma_ip(c,r).ne.LIS_rc%udef) then
              sigma_fill(c-1,r-1)= sigma_fill(c-1,r-1) + sigma_ip(c,r)
              sigma_fill(c-1,r)= sigma_fill(c-1,r) + sigma_ip(c,r)
              sigma_fill(c-1,r+1)= sigma_fill(c-1,r+1) + sigma_ip(c,r)
              sigma_fill(c,r-1)= sigma_fill(c,r-1) + sigma_ip(c,r)
              sigma_fill(c,r+1)= sigma_fill(c,r+1) + sigma_ip(c,r)
              sigma_fill(c+1,r-1)= sigma_fill(c+1,r-1) + sigma_ip(c,r)
              sigma_fill(c+1,r)= sigma_fill(c+1,r) + sigma_ip(c,r)
              sigma_fill(c+1,r+1)= sigma_fill(c+1,r+1) + sigma_ip(c,r)
              
              nsigma_ip(c-1,r-1)= nsigma_ip(c-1,r-1) + 1
              nsigma_ip(c-1,r)= nsigma_ip(c-1,r) + 1
              nsigma_ip(c-1,r+1)= nsigma_ip(c-1,r+1) + 1
              nsigma_ip(c,r-1)= nsigma_ip(c,r-1) + 1
              nsigma_ip(c,r+1)= nsigma_ip(c,r+1) + 1
              nsigma_ip(c+1,r-1)= nsigma_ip(c+1,r-1) + 1
              nsigma_ip(c+1,r)= nsigma_ip(c+1,r) + 1
              nsigma_ip(c+1,r+1)= nsigma_ip(c+1,r+1) + 1
           endif
        enddo
     enddo

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if (sigma_ip(c,r).eq.LIS_rc%udef.and.nsigma_ip(c,r).gt.0) then
              ! average in linear scale
              sigma_ip(c,r) = sigma_fill(c,r)/nsigma_ip(c,r) 
           endif
           if (sigma_ip(c,r).ne.LIS_rc%udef) then
             ! back to log scale, dB
             sigma_ip(c,r) = 10.*LOG10(sigma_ip(c,r))
           endif
        enddo
     enddo

  endif

#endif

end subroutine read_S1_sigmaVVSM_data



!BOP
!
! !ROUTINE: S1_sigmaVVSM_filename
! \label{S1_SWND_filename}
! 
! !INTERFACE: 
subroutine S1_sigmaVVSM_filename(filename, ndir, yr, mo, da)
  
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
    
end subroutine S1_sigmaVVSM_filename



