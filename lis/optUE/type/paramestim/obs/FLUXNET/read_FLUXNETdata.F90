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
! !ROUTINE: read_FLUXNETdata
! \label{read_FLUXNETdata}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_FLUXNETdata(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_fileIOMod,      only : LIS_readData
  use LIS_timeMgrMod,     only : LIS_calendar, LIS_tick
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use map_utils
  use FLUXNETdata_module,     only : FLUXNETdata_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)  
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP
  integer                  :: n 
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer                  :: ftn, qleId,qhId
  logical                  :: file_exists
  integer                  :: c,r,t
  integer                  :: status
  integer                  :: cmo
  integer                  :: iret
  real,  pointer           :: qle_obs(:),qh_obs(:)
  type(ESMF_Field)         :: qleField,qhField
  real, allocatable            :: qle_ip(:),qh_ip(:)
  logical*1,allocatable        :: li(:)
  logical*1,allocatable        :: lo(:)
  real,     allocatable        :: qle(:,:,:),qh(:,:,:)
  real,     allocatable        :: qle1d(:),qh1d(:)
  integer                  :: nyr,nmo,nda,nhr,nmn,nss
  type(ESMF_Time)          :: currTime
  type(ESMF_TimeInterval)  :: ts

#if (defined USE_NETCDF3 || defined USE_NETCDF4)    
  n = 1
  allocate(qle_ip(LIS_rc%ngrid(n)))
  qle_ip = LIS_rc%udef
  allocate(qh_ip(LIS_rc%ngrid(n)))
  qh_ip = LIS_rc%udef

! read data when year changes. 
!
  if(FLUXNETdata_struc(n)%yr.ne.LIS_rc%yr) then 

     FLUXNETdata_struc(n)%yr = LIS_rc%yr

     call create_fluxnet_lh_filename(FLUXNETdata_struc(n)%odir,&
          LIS_rc%yr, filename)
     inquire(file=trim(filename),exist=file_exists)

     if(file_exists) then 
        write(LIS_logunit,*) 'Reading FLUXNET LH file ',trim(filename)
        
        call LIS_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,&
             ncid=ftn),'nf90_open failed in readFLUXETdata')
        call LIS_verify(nf90_inq_varid(ftn,'EnsembleLEcor_May09',qleId),&
             'nf90_inq_varid failed in read_FLUXNETdata')
        
        allocate(qle(FLUXNETdata_struc(n)%nc,FLUXNETdata_struc(n)%nr,12))
        allocate(qle1d(FLUXNETdata_struc(n)%nc*FLUXNETdata_struc(n)%nr))
        allocate(li(FLUXNETdata_struc(n)%nc*FLUXNETdata_struc(n)%nr))
        allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

        call LIS_verify(nf90_get_var(ftn,qleId,qle),&
             'nf90_get_var for qle failed in read_FLUXNETdata')
        call LIS_verify(nf90_close(ftn))

        FLUXNETdata_struc(n)%qle = LIS_rc%udef

        do t=1,12
           li = .false. 
           do r=1,FLUXNETdata_struc(n)%nr
              do c=1,FLUXNETdata_struc(n)%nc
                 if(isNaN(qle(c,FLUXNETdata_struc(n)%nr-r+1,t))) then 
                    qle1d(c+(r-1)*FLUXNETdata_struc(n)%nc)= &
                         LIS_rc%udef
                 elseif(qle(c,FLUXNETdata_struc(n)%nr-r+1,t).lt.0) then 
                    qle1d(c+(r-1)*FLUXNETdata_struc(n)%nc)= &
                         LIS_rc%udef
                 else
                    qle1d(c+(r-1)*FLUXNETdata_struc(n)%nc)= &
                         qle(c,FLUXNETdata_struc(n)%nr-r+1,t)*0.01*1E6/86400.0
                    li(c+(r-1)*FLUXNETdata_struc(n)%nc) = .true. 
                 endif
              enddo
           enddo

           call neighbor_interp(LIS_rc%gridDesc(n,:),li,&
                qle1d,lo,FLUXNETdata_struc(n)%qle(:,t),&
                FLUXNETdata_struc(n)%nc*FLUXNETdata_struc(n)%nr,&
                LIS_rc%lnc(n)*LIS_rc%lnr(n),&
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                FLUXNEtdata_struc(n)%n11,LIS_rc%udef, iret)
        enddo
        
        deallocate(li)
        deallocate(lo)
        deallocate(qle)
        deallocate(qle1d)
     endif
     call create_fluxnet_sh_filename(FLUXNETdata_struc(n)%odir,&
          LIS_rc%yr, filename)
     inquire(file=trim(filename),exist=file_exists)

     if(file_exists) then 
        write(LIS_logunit,*) 'Reading FLUXNET SH file ',trim(filename)
        
        call LIS_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,&
             ncid=ftn),'nf90_open failed in readFLUXETdata')
        call LIS_verify(nf90_inq_varid(ftn,'EnsembleHcor_Aug09',qhId),&
             'nf90_inq_varid failed in read_FLUXNETdata')
        
        allocate(qh(FLUXNETdata_struc(n)%nc,FLUXNETdata_struc(n)%nr,12))
        allocate(qh1d(FLUXNETdata_struc(n)%nc*FLUXNETdata_struc(n)%nr))
        allocate(li(FLUXNETdata_struc(n)%nc*FLUXNETdata_struc(n)%nr))
        allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

        call LIS_verify(nf90_get_var(ftn,qhId,qh),&
             'nf90_get_var for qh failed in read_FLUXNETdata')
        call LIS_verify(nf90_close(ftn))

        FLUXNETdata_struc(n)%qh = LIS_rc%udef

        do t=1,12
           li = .false. 
           do r=1,FLUXNETdata_struc(n)%nr
              do c=1,FLUXNETdata_struc(n)%nc
                 if(isNaN(qh(c,FLUXNETdata_struc(n)%nr-r+1,t))) then 
                    qh1d(c+(r-1)*FLUXNETdata_struc(n)%nc)= &
                         LIS_rc%udef
                 elseif(qh(c,FLUXNETdata_struc(n)%nr-r+1,t).lt.0) then 
                    qh1d(c+(r-1)*FLUXNETdata_struc(n)%nc)= &
                         LIS_rc%udef
                 else
                    qh1d(c+(r-1)*FLUXNETdata_struc(n)%nc)= &
                         qh(c,FLUXNETdata_struc(n)%nr-r+1,t)*0.01*1E6/86400.0
                    li(c+(r-1)*FLUXNETdata_struc(n)%nc) = .true. 
                 endif
              enddo
           enddo
 
           call neighbor_interp(LIS_rc%gridDesc(n,:),li,&
                qh1d,lo,FLUXNETdata_struc(n)%qh(:,t),&
                FLUXNETdata_struc(n)%nc*FLUXNETdata_struc(n)%nr,&
                LIS_rc%lnc(n)*LIS_rc%lnr(n),&
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                FLUXNEtdata_struc(n)%n11,LIS_rc%udef, iret)

        enddo
        
        deallocate(li)
        deallocate(lo)
        deallocate(qh)
        deallocate(qh1d)
     endif
  endif
!log data at the change of month

  call ESMF_TimeSet(currTime,  yy=LIS_rc%yr, &
       mm = LIS_rc%mo, &
       dd = LIS_rc%da, &
       h = LIS_rc%hr, &
       m = LIS_rc%mn, &
       s = LIS_rc%ss,&
       calendar = LIS_calendar, &
       rc=status)
  call LIS_verify(status,'error in ESMF_TimeSet:read_FLUXNETdata')
  
  call ESMF_TimeIntervalSet(ts,s=nint(LIS_rc%ts),rc=status)
  call LIS_verify(status,'error in ESMF_TimeSet:read_FLUXNETdata')
  
  currTime = currTime + ts
  call ESMF_TimeGet(currTime,  yy=nyr, &
       mm = nmo, &
       dd = nda, &
       h = nhr, &
       m = nmn, &
       s = nss,&
       calendar = LIS_calendar, &
       rc=status)
  call LIS_verify(status,'error in ESMF_TimeSet:read_FLUXNETdata')
  if(nmo.ne.FLUXNETdata_struc(n)%mo) then 

     write(LIS_logunit,*) 'Reading FLUXNET data for month ',&
          nmo

     FLUXNETdata_struc(n)%mo = nmo
     cmo = FLUXNETdata_struc(n)%mo
     
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           t = LIS_domain(n)%gindex(c,r)
           if(t.ne.-1) then 
              qle_ip(t) = FLUXNETdata_struc(n)%qle(&
                   c+(r-1)*LIS_rc%lnc(n),cmo)
           endif
        enddo
     enddo

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           t = LIS_domain(n)%gindex(c,r)
           if(t.ne.-1) then 
              qh_ip(t) = FLUXNETdata_struc(n)%qh(&
                   c+(r-1)*LIS_rc%lnc(n),cmo)
           endif
        enddo
     enddo
  else
     qle_ip = LIS_rc%udef
     qh_ip = LIS_rc%udef
  endif

  call ESMF_StateGet(Obj_Space,"FLUXNET Latent Heat Flux",qleField,&
       rc=status)
  call LIS_verify(status, "StateGet failed for FLUXNET Latent Heat Flux in read_FLUXNETdata")
  call ESMF_FieldGet(qleField,localDE=0,farrayPtr=qle_obs,rc=status)
  call LIS_verify(status, "FieldGet failed for FLUXNET Latent Heat Flux in read_FLUXNETdata") 
  qle_obs = qle_ip

  call ESMF_StateGet(Obj_Space,"FLUXNET Sensible Heat Flux",qhField,&
       rc=status)
  call LIS_verify(status, "StateGet failed for FLUXNET Sensible Heat Flux in read_FLUXNETdata")
  call ESMF_FieldGet(qhField,localDE=0,farrayPtr=qh_obs,rc=status)
  call LIS_verify(status, "FieldGet failed for FLUXNET Sensible Heat Flux in read_FLUXNETdata") 
  qh_obs = qh_ip

! Set this to true even if the data logged is undefined. This is to ensure that 
! the code will enter the update objective function call at each timestep 
! so that it can perform computations of the model climatology (if optimization
! is done against a climatology dataset such as the gridded fluxnet data, 
! which is a monthly average). Setting this flag to true will enable
! the continuous update and calculation of monthly averages of the LIS
! model simulated fluxes. 

  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status)
#endif

end subroutine read_FLUXNETdata


!BOP
! 
! !ROUTINE: create_fluxnet_lh_filename
! \label{create_fluxnet_lh_filename}
! 
! !INTERFACE: 
subroutine create_fluxnet_lh_filename(odir,yr,filename)
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for FLUXNET latent heat 
! data files based on the given date
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the FLUXNET\_SH file
!  \end{description}
!EOP
  character(len=*)         :: odir
  integer                  :: yr
  character(len=*)         :: filename
  
  character*4            :: fyr

  write(unit=fyr,fmt='(i4.4)') yr
  filename = trim(odir)//'/EnsembleLEcor_May09_'//trim(fyr)//'.nc'

end subroutine create_fluxnet_lh_filename


!BOP
! 
! !ROUTINE: create_fluxnet_sh_filename
! \label{create_fluxnet_sh_filename}
! 
! !INTERFACE: 
subroutine create_fluxnet_sh_filename(odir,yr,filename)
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for FLUXNET latent heat 
! data files based on the given date
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the FLUXNET\_SH file
!  \end{description}
! 
!EOP
  character(len=*)         :: odir
  integer                  :: yr
  character(len=*)         :: filename
  
  character*4            :: fyr

  write(unit=fyr,fmt='(i4.4)') yr
  filename = trim(odir)//'/EnsembleHcor_Aug09_'//trim(fyr)//'.nc'

end subroutine create_fluxnet_sh_filename
