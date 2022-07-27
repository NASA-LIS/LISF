!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_MCD64A1ba
! \label{read_MCD64A1ba}
!
! !REVISION HISTORY:
!  24 Jul 2022  Sujay Kumar; initial specification
!
! !INTERFACE: 
subroutine read_MCD64A1ba(n, k, OBS_State, OBS_Pert_State)
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
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use MCD64A1ba_Mod, only : MCD64A1ba_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the MCD64A1 burn area observations from NETCDF files.

! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[k] number of observation state
!  \item[OBS\_State] observations state
!  \item[OBS\_Pert\_State] observation perturbations state
!  \end{description}
!
!EOP
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: baobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname1,fname2, climofile1, climofile2
  integer                :: cyr, cmo, cda, chr,cmn,css,cdoy
  real                   :: wt1, wt2,ts
  integer                :: count
  real                   :: cgmt
  real*8                 :: time
  logical                :: alarmCheck, file_exists,dataCheck
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: bafield, pertField
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: baobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                :: fnd
  real                   :: timenow

  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       baobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "MCD64A1 burned area read alarm")

  if(alarmCheck.or.MCD64A1ba_struc(n)%startMode) then 
     MCD64A1ba_struc(n)%startMode = .false.

     call LIS_tick(time,cdoy,cgmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
          LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,0.0)
     
     call create_MCD64A1ba_filename(baobsdir, &
          MCD64A1ba_struc(n)%version, &
          LIS_rc%yr, cdoy, fname1)
         
     inquire(file=fname1,exist=file_exists)

     if(file_exists) then 
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname1)
        call read_MCD64A1_BA_data(n,k, fname1,baobs)
        fnd = 1
        dataCheck = .true.
     else
        fnd = 0
        dataCheck = .false.
        write(LIS_logunit,*) '[WARN] Missing BA file: ',trim(fname1)
     endif     
  else
     dataCheck = .false.
     fnd = 0 
     baobs = LIS_rc%udef
  endif


  if(dataCheck) then 
        
     call ESMF_StateGet(OBS_State,"Observation01",bafield,&
          rc=status)
     call LIS_verify(status, 'Error: StateGet Observation01')
     
     call ESMF_FieldGet(bafield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status, 'Error: FieldGet')
     
     obsl = LIS_rc%udef 
     do r=1, LIS_rc%obs_lnr(k)
        do c=1, LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                   baobs(c+(r-1)*LIS_rc%obs_lnc(k))
           endif
        enddo
     enddo
     
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
           call ESMF_AttributeSet(bafield,"Grid Number",&
                gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(bafield,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
        endif
     endif
  end if
  call ESMF_AttributeSet(OBS_State,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status)     

!  else
!     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
!          .false., rc=status)
!     call LIS_verify(status)     
!  endif
end subroutine read_MCD64A1ba

!BOP
! 
! !ROUTINE: read_MCD64A1_BA_data
! \label{read_MCD64A1_BA_data}
!
! !INTERFACE:
subroutine read_MCD64A1_BA_data(n, k, fname, baobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use MCD64A1ba_Mod, only : MCD64A1ba_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: baobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real*8                        :: cornerlat(2), cornerlon(2)
  character*3                   :: fdoy

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the MCD64A1 BA file and applies the data
!  quality flags to filter the data. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[k]            number of observation state
!  \item[k]            number of observation state
!  \item[fname]        name of the MCD64A1 BA file
!  \item[climofile]    Generated MCD152AH BA climatology file
!  \item[baobs\_ip]   MCD64A1 BA data processed to the LIS domain
!  \end{description}
!
!
!EOP

!--------------Wanshu -----------------------
  integer,  parameter     :: nc=86400, nr=43200
  integer                 :: lat_off, lon_off
  integer                 :: ba(MCD64A1ba_struc(n)%nc,MCD64A1ba_struc(n)%nr)
  integer                 :: flag(MCD64A1ba_struc(n)%nc,MCD64A1ba_struc(n)%nr)
  real                    :: ba_flagged(MCD64A1ba_struc(n)%nc,MCD64A1ba_struc(n)%nr)
  real                    :: ba_in(MCD64A1ba_struc(n)%nc*MCD64A1ba_struc(n)%nr)
  logical*1               :: ba_data_b(MCD64A1ba_struc(n)%nc*MCD64A1ba_struc(n)%nr)
  logical*1               :: baobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                 :: c,r,t
  integer                 :: nid
  integer                 :: baid, flagid
  integer                 :: ios
  
  integer, dimension(nf90_max_var_dims) :: dimIDs
  integer                                :: numLons, numLats
  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LIS_verify(ios,'Error opening file '//trim(fname))
  
  ios = nf90_inq_varid(nid, "Burn Date",baid)
  call LIS_verify(ios, 'Error nf90_inq_varid: Burn Date')
  
  !values
  
  cornerlat(1)=MCD64A1ba_struc(n)%gridDesci(4)
  cornerlon(1)=MCD64A1ba_struc(n)%gridDesci(5)
  cornerlat(2)=MCD64A1ba_struc(n)%gridDesci(7)
  cornerlon(2)=MCD64A1ba_struc(n)%gridDesci(8)
  
  ba_data_b = .false.
  
  lat_off = nint((cornerlat(1)+89.9979167)/0.00416667)+1
  lon_off = nint((cornerlon(1)+179.9979167)/0.00416667)+1


  ios = nf90_get_var(nid, baid, ba, &
       start=(/lon_off,lat_off/), &
       count=(/MCD64A1ba_struc(n)%nc,MCD64A1ba_struc(n)%nr/)) 
  
  call LIS_verify(ios, 'Error nf90_get_var: Burn Date')
  
  ios = nf90_close(ncid=nid)
  call LIS_verify(ios,'Error closing file '//trim(fname))
  
  do r=1, MCD64A1ba_struc(n)%nr
     do c=1, MCD64A1ba_struc(n)%nc
        
        if(ba(c,r).gt.0.and.ba(c,r).le.400) then
           ba_flagged(c,r) = ba(c,r)
        else
           ba_flagged(c,r) = LIS_rc%udef
        endif
        
     end do
  end do


  do r=1, MCD64A1ba_struc(n)%nr
     do c=1, MCD64A1ba_struc(n)%nc
        ba_in(c+(r-1)*MCD64A1ba_struc(n)%nc) = ba_flagged(c,r)
        if(ba_flagged(c,r).ne.LIS_rc%udef) then
           ba_data_b(c+(r-1)*MCD64A1ba_struc(n)%nc) = .true.
        else
           ba_data_b(c+(r-1)*MCD64A1ba_struc(n)%nc) = .false.
        endif
     enddo
  enddo

  if(LIS_rc%obs_gridDesc(k,10).le.0.00416667) then 
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
          ba_data_b, ba_in, baobs_b_ip, baobs_ip, &
          MCD64A1ba_struc(n)%nc*MCD64A1ba_struc(n)%nr, &
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          MCD64A1ba_struc(n)%rlat,MCD64A1ba_struc(n)%rlon,&
          MCD64A1ba_struc(n)%w11,MCD64A1ba_struc(n)%w12,&
          MCD64A1ba_struc(n)%w21,MCD64A1ba_struc(n)%w22,&
          MCD64A1ba_struc(n)%n11,MCD64A1ba_struc(n)%n12,&
          MCD64A1ba_struc(n)%n21,MCD64A1ba_struc(n)%n22,LIS_rc%udef,ios)
  else
     call upscaleByAveraging(MCD64A1ba_struc(n)%nc*MCD64A1ba_struc(n)%nr,&
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          LIS_rc%udef, MCD64A1ba_struc(n)%n11,&
          ba_data_b,ba_in, baobs_b_ip, baobs_ip)

  endif

#endif

end subroutine read_MCD64A1_BA_data


!BOP
! !ROUTINE: create_MCD64A1ba_filename
! \label{create_MCD64A1ba_filename}
! 
! !INTERFACE: 
subroutine create_MCD64A1ba_filename(ndir, version, yr, doy, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: version
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the MCD64A1 BA filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the MCD64A1 BA data directory
!  \item[version] version of the MCD64A1 BA data
!  \item[yr]  current year
!  \item[doy]  current day of the year
!  \item[filename] Generated MCD64A1 BA filename
!  \item[climofile] Generated MCD152AH BA climatology file
!  \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  if(version.eq."006") then
     filename = trim(ndir)//'/'//trim(fyr)//'/MCD64A1.006_'//&
          trim(fyr)//trim(fdoy)//'.nc4'
  endif


end subroutine create_MCD64A1ba_filename





