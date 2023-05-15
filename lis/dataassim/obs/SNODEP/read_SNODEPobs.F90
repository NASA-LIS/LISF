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
! !ROUTINE: read_SNODEPobs
! \label{read_SNODEPobs}
!
! !REVISION HISTORY:
!  21 Aug 2008: Sujay Kumar; Initial Specification
!  25 Jan 2012: Sujay Kumar; Updated to use grib api
!  13 Jul 2016: Sujay Kumar; Updated the processing to the Observation space
! 
! !INTERFACE: 
subroutine read_SNODEPobs(n, k, OBS_State,OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_historyMod, only : LIS_readvar_gridded
  use LIS_timeMgrMod, only : LIS_clock
  use LIS_logMod,     only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use SNODEPobs_Mod,  only : SNODEP_obs_obj
  use LIS_mpiMod
#if (defined USE_GRIBAPI)
  use grib_api
#endif
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the AFWA/557 SNODEP observations and processes it
!  for for later use within the DA algorithm. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[k] index of the data assimilation instance
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: snowfield

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))

  character(len=LIS_CONST_PATH_LEN) :: obsdir
  logical             :: data_update
  logical             :: file_exists1,file_exists2
  character(len=LIS_CONST_PATH_LEN) :: name_nh, name_sh

  real                :: value
  logical             :: readflag
  integer             :: rc,status
  integer             :: t,c,r
  integer             :: ftn1,ftn2,iret
  integer             :: npts
  logical*1           :: lb1(SNODEP_obs_obj(n)%pmax*SNODEP_obs_obj(n)%pmax)
  logical*1           :: lb2(SNODEP_obs_obj(n)%pmax*SNODEP_obs_obj(n)%pmax)
  real                :: f1(SNODEP_obs_obj(n)%pmax*SNODEP_obs_obj(n)%pmax)
  real                :: f2(SNODEP_obs_obj(n)%pmax*SNODEP_obs_obj(n)%pmax)
  integer             :: hemi
  real                :: gi(2,SNODEP_obs_obj(n)%pmax,SNODEP_obs_obj(n)%pmax)
  real                :: varfield(LIS_rc%obs_lnc(n), LIS_rc%obs_lnr(n))
  logical             :: alarmCheck
  real                :: missingValue
  integer             :: igrib
  integer             :: pds5,pds7,pds5_val,pds7_val
  integer             :: yr,mo,da,hr,mn,ss
  type(ESMF_Time)     :: currTime
  type(ESMF_TimeInterval) :: deltaT

#if(defined USE_GRIBAPI)
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  file_exists1 = .false.
  file_exists2 = .false.
  name_nh = ''
  name_sh = ''

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SNODEP read alarm")    

  if(alarmCheck) then 
     hemi = 1
     call ESMF_TimeIntervalSet(deltaT,s=nint(LIS_rc%ts))
     call ESMF_ClockGet(LIS_clock,currTime=currTime,rc=status)       
     currTime = currTime + deltaT

     call ESMF_TimeGet(currTime,yy=yr,mm=mo,dd=da,h=hr,m=mn,s=ss,rc=status)
     call SNODEP_filename(name_nh, SNODEP_obs_obj(n)%mesh, hemi, obsdir,&
          yr, mo, da, hr, mn, SNODEP_obs_obj(n)%conv)
     hemi = 2
     call SNODEP_filename(name_sh, SNODEP_obs_obj(n)%mesh, hemi, obsdir,&
          yr, mo, da, hr, mn, SNODEP_obs_obj(n)%conv)
     
     inquire(file=name_nh,exist=file_exists1)
     inquire(file=name_sh,exist=file_exists2)
  endif

  if(file_exists1.and.file_exists2) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif
  
  if (readflag) then 
     
     call ESMF_StateGet(OBS_State,"Observation01",snowfield,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(snowfield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)

     if(file_exists1) then 
        npts = SNODEP_obs_obj(n)%pmax*SNODEP_obs_obj(n)%pmax
        if(SNODEP_obs_obj(n)%mesh.eq.8) then 
           pds5 = 174 
           pds7 = 0
        else
           pds5 = 66 
           pds7 = 0
        endif
        write(LIS_logunit,*) '[INFO] SNODEP file ',trim(name_nh)

        call grib_open_file(ftn1,trim(name_nh),'r',iret)
        call LIS_verify(iret,'error in grib_open_file in read_SNODEPobs')

        call grib_new_from_file(ftn1,igrib,iret)
        call LIS_verify(iret,'error in grib_new_from_file in read_SNODEPobs')

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LIS_verify(rc,'error in grib_get:indicatorOfParameter in read_SNODEPobs')

        call grib_get(igrib,'level',pds7_val,rc)
        call LIS_verify(rc,'error in grib_get:level in read_SNODEPobs')

        f1 = -9999.0
        if((pds5.eq.pds5_val).and.(pds7.eq.pds7_val)) then

           call grib_get(igrib,'values',f1,rc)
           call LIS_verify(rc,'error in grib_get:values in read_SNODEPobs')

           call grib_get(igrib,'missingValue',missingValue,rc)
           call LIS_verify(rc, 'error in grib_get:missingValue in read_SNODEPobs')
           
           call grib_release(igrib,rc)
           call LIS_verify(rc, 'error in grib_release in read_SNODEPobs')
        else
           write(LIS_logunit,*) '[ERR] Unable to retrieve values from ',trim(name_nh)
           call LIS_endrun()

        endif
        call grib_close_file(ftn1)

        do r=1,SNODEP_obs_obj(n)%pmax
           do c=1,SNODEP_obs_obj(n)%pmax
              value = f1(c+(r-1)*SNODEP_obs_obj(n)%pmax)
              gi(1,c,r) = value
              if(value.gt.408.9) gi(1,c,r) = LIS_rc%udef
           enddo
        enddo
     else
        gi(1,:,:) = LIS_rc%udef
     endif
     
     if(file_exists2) then 
        npts = SNODEP_obs_obj(n)%pmax*SNODEP_obs_obj(n)%pmax
        if(SNODEP_obs_obj(n)%mesh.eq.8) then 
           pds5 = 174 
           pds7 = 0
        else
           pds5 = 66 
           pds7 = 0
        endif
        write(LIS_logunit,*) '[INFO] SNODEP file ',trim(name_sh)

        call grib_open_file(ftn2,trim(name_sh),'r',iret)
        call LIS_verify(iret,'error in grib_open_file in read_SNODEPobs')

        call grib_new_from_file(ftn2,igrib,iret)
        call LIS_verify(iret,'error in grib_new_from_file in read_SNODEPobs')

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LIS_verify(rc,'error in grib_get:indicatorOfParameter in read_SNODEPobs')

        call grib_get(igrib,'level',pds7_val,rc)
        call LIS_verify(rc,'error in grib_get:level in read_SNODEPobs')

        f2 = -9999.0
        if((pds5.eq.pds5_val).and.(pds7.eq.pds7_val)) then

           call grib_get(igrib,'values',f2,rc)
           call LIS_verify(rc,'error in grib_get:values in read_SNODEPobs')

           call grib_get(igrib,'missingValue',missingValue,rc)
           call LIS_verify(rc, 'error in grib_get:missingValue in read_SNODEPobs')
           
           call grib_release(igrib,rc)
           call LIS_verify(rc, 'error in grib_release in read_SNODEPobs')
        else
           write(LIS_logunit,*) '[ERR] Unable to retrieve values from ',trim(name_sh)
           call LIS_endrun()

        endif
        call grib_close_file(ftn2)

        do r=1,SNODEP_obs_obj(n)%pmax
           do c=1,SNODEP_obs_obj(n)%pmax
              value = f2(c+(r-1)*SNODEP_obs_obj(n)%pmax)
              gi(2,c,r) = value
              if(value.gt.408.9) gi(2,c,r) = LIS_rc%udef
           enddo
        enddo
     else
        gi(2,:,:)=LIS_rc%udef
     endif
     
     if(file_exists1.or.file_exists2) then 
        call interp_SNODEPfield(n, k, 1,gi,LIS_rc%udef,varfield)
     else
        varfield = LIS_rc%udef
     endif

     do r =1,LIS_rc%obs_lnr(n)
        do c =1,LIS_rc%obs_lnc(n)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=varfield(c,r)
           end if
        end do
     end do
     
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true. , rc=status)
     call LIS_verify(status)

     do t=1,LIS_rc%obs_ngrid(k)
        gid(t) = t
        if(obsl(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo
     
     if(LIS_rc%obs_ngrid(k).gt.0) then 
        call ESMF_AttributeSet(snowfield,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)
        
        call ESMF_AttributeSet(snowfield,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)
     endif

  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  endif

#endif
  
end subroutine read_SNODEPobs

subroutine  SNODEP_filename(name, mesh, hemi, ndir, yr, mo, da, hr, mn, conv)

  implicit none

  character*100     :: name
  integer           :: hemi
  integer           :: mesh
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr,fmn
  character (len=*) :: conv
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  write(unit=fmn, fmt='(i2.2)') mn  

  if(mesh.eq.8) then 
     if(hemi.eq.1) then 
        if ( conv == "LIS" ) then
           name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//&
                  '/SNODEP/depth_nh.'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)
        else
           name = trim(ndir)//'/depth_nh.'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)
        endif
     elseif(hemi.eq.2) then 
        if ( conv == "LIS" ) then
           name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//&
                  '/SNODEP/depth_sh.'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)
        else
           name = trim(ndir)//'/depth_sh.'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)
        endif
     endif
  elseif(mesh.eq.16) then 
     if(hemi.eq.1) then 
        if ( conv == "LIS" ) then
           name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//&
                  '/SNODEP/SNODEP_16_NH_'//trim(fyr)//trim(fmo)//&
                  trim(fda)//trim(fhr)//'.GR1'
        else
           name = trim(ndir)//'/SNODEP_16_NH_'//trim(fyr)//trim(fmo)//&
                  trim(fda)//trim(fhr)//'.GR1'
        endif
     elseif(hemi.eq.2) then 
        if ( conv == "LIS" ) then
           name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//&
                  '/SNODEP/SNODEP_16_SH_'//trim(fyr)//trim(fmo)//&
                  trim(fda)//trim(fhr)//'.GR1'
        else
           name = trim(ndir)//'/SNODEP_16_SH_'//trim(fyr)//trim(fmo)//&
                  trim(fda)//trim(fhr)//'.GR1'
        endif
     endif
  endif
end subroutine SNODEP_filename

!BOP
! !ROUTINE: interp_SNODEPfield
!  \label{interp_SNODEPfield}
! 
! !INTERFACE:
subroutine interp_SNODEPfield(n,k,ip,gi,udef,varfield)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use SNODEPobs_Mod,     only : SNODEP_obs_obj

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n 
  integer, intent(in)    :: k
  integer, intent(in)    :: ip
  real, intent(in)       :: gi(2,SNODEP_obs_obj(n)%pmax,SNODEP_obs_obj(n)%pmax)
  real, intent(in)       :: udef
  real, intent(inout)    :: varfield(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
!
! !DESCRIPTION:
!   This subroutine interpolates a given AFWA field 
!   to the LIS domain (from polar stereographic grid to the LIS grid)
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[ip]
!    interpolation option
!  \item[gi]
!    input AGRMET field (both hemispheres)
!  \item[udef]
!    undefined value in the input field
!  \item[varfield]
!    interpolated field in the LIS grid
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using neighbor interpolation
! \end{description}
!EOP
  real                   :: gi_temp(SNODEP_obs_obj(n)%mi)
  logical*1              :: li(SNODEP_obs_obj(n)%mi)
  integer                :: ihemi,i,j,iret
!  real                   :: gridDesco(200)
  logical*1, allocatable :: lo_nh(:),lo_sh(:)
  real, allocatable      :: go_nh(:),go_sh(:)
  real, allocatable      :: comb_data(:)

!TBD - needs updating and testing. 
  if(ip.eq.1) then 
     if(SNODEP_obs_obj(n)%gridspan.eq.1) then 
        allocate(lo_nh(SNODEP_obs_obj(n)%mo1))
        allocate(go_nh(SNODEP_obs_obj(n)%mo1))
     elseif(SNODEP_obs_obj(n)%gridspan.eq.2) then 
        allocate(lo_sh(SNODEP_obs_obj(n)%mo2))
        allocate(go_sh(SNODEP_obs_obj(n)%mo2))
     else
        allocate(lo_nh(SNODEP_obs_obj(n)%mo1))
        allocate(lo_sh(SNODEP_obs_obj(n)%mo2))
        allocate(go_nh(SNODEP_obs_obj(n)%mo1))
        allocate(go_sh(SNODEP_obs_obj(n)%mo2))
     endif
     allocate(comb_data(SNODEP_obs_obj(n)%mo1+SNODEP_obs_obj(n)%mo2))
     
     do ihemi = SNODEP_obs_obj(n)%shemi,SNODEP_obs_obj(n)%nhemi

        li = .false.
        gi_temp = LIS_rc%udef

        do i=1,SNODEP_obs_obj(n)%pmax
           do j=1,SNODEP_obs_obj(n)%pmax
              if(gi(ihemi,i,j).ne.udef) then 
                 li(i+(j-1)*SNODEP_obs_obj(n)%pmax) = .true.
                 gi_temp(i+(j-1)*SNODEP_obs_obj(n)%pmax) = gi(ihemi,i,j)
!                 if(gi(ihemi,i,j).gt.30) then 
!                    print*, ihemi,i,j,gi(ihemi,i,j)
!                 endif
              endif
           enddo
        enddo
!        gridDesco(1) = 0 
!        gridDesco(2) = SNODEP_obs_obj(n)%hemi_nc(ihemi)
!        gridDesco(3) = SNODEP_obs_obj(n)%hemi_nr(ihemi)
!        gridDesco(5) = LIS_rc%obs_gridDesc(k,5)
!        gridDesco(8) = LIS_rc%obs_gridDesc(k,8)
!        gridDesco(6) = LIS_rc%obs_gridDesc(k,6)
!        gridDesco(9) = LIS_rc%obs_gridDesc(k,9)
!        gridDesco(10) = LIS_rc%obs_gridDesc(k,10)
!        gridDesco(20) = 0 
!        if(SNODEP_obs_obj(n)%gridspan.eq.1.or.SNODEP_obs_obj(n)%gridspan.eq.2) then 
!           gridDesco(4) = LIS_rc%obs_gridDesc(k,4)
!           gridDesco(7) = LIS_rc%obs_gridDesc(k,7)
!        else
!           if(ihemi.eq.1) then 
!              gridDesco(4) = LIS_rc%obs_gridDesc(k,9)/2
!              gridDesco(7) = LIS_rc%obs_gridDesc(k,7)
!           else
!              gridDesco(4) = LIS_rc%obs_gridDesc(k,4)
!              gridDesco(7) = -LIS_rc%obs_gridDesc(k,9)/2
!           endif
!     endif
        if(ihemi.eq.1) then 
           call neighbor_interp(SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                li,gi_temp,&
                lo_nh,go_nh,SNODEP_obs_obj(n)%mi,SNODEP_obs_obj(n)%mo1,&
                SNODEP_obs_obj(n)%rlat1_nh,SNODEP_obs_obj(n)%rlon1_nh,&
                SNODEP_obs_obj(n)%n111_nh, & 
                LIS_rc%udef,iret)
#if 0 
           call bilinear_interp(SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                li,gi_temp,&
                lo_nh,go_nh,SNODEP_obs_obj(n)%mi,SNODEP_obs_obj(n)%mo1,&
                SNODEP_obs_obj(n)%rlat1_nh,SNODEP_obs_obj(n)%rlon1_nh,&
                SNODEP_obs_obj(n)%w111_nh,SNODEP_obs_obj(n)%w121_nh,&
                SNODEP_obs_obj(n)%w211_nh,SNODEP_obs_obj(n)%w221_nh,&
                SNODEP_obs_obj(n)%n111_nh,SNODEP_obs_obj(n)%n121_nh,&
                SNODEP_obs_obj(n)%n211_nh,SNODEP_obs_obj(n)%n221_nh,&
                LIS_rc%udef,iret)
#endif
        elseif(ihemi.eq.2) then 
           call neighbor_interp(SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                li,gi_temp,&
                lo_sh,go_sh,SNODEP_obs_obj(n)%mi,SNODEP_obs_obj(n)%mo2,&
                SNODEP_obs_obj(n)%rlat1_sh, SNODEP_obs_obj(n)%rlon1_sh,&
                SNODEP_obs_obj(n)%n111_sh,&
                LIS_rc%udef,iret)
#if 0 
           call bilinear_interp(SNODEP_obs_obj(n)%gridDesco(ihemi,:),&
                li,gi_temp,&
                lo_sh,go_sh,SNODEP_obs_obj(n)%mi,SNODEP_obs_obj(n)%mo2,&
                SNODEP_obs_obj(n)%rlat1_sh, SNODEP_obs_obj(n)%rlon1_sh,&
                SNODEP_obs_obj(n)%w111_sh,SNODEP_obs_obj(n)%w121_sh,&
                SNODEP_obs_obj(n)%w211_sh,SNODEP_obs_obj(n)%w221_sh,&
                SNODEP_obs_obj(n)%n111_sh,SNODEP_obs_obj(n)%n121_sh,&
                SNODEP_obs_obj(n)%n211_sh,SNODEP_obs_obj(n)%n221_sh,&
                LIS_rc%udef,iret)
#endif
        endif
     end do
     if(SNODEP_obs_obj(n)%gridspan.eq.1) then
        comb_data(1:SNODEP_obs_obj(n)%mo2) = LIS_rc%udef
        comb_data(SNODEP_obs_obj(n)%mo2+1:SNODEP_obs_obj(n)%mo1+SNODEP_obs_obj(n)%mo2) = go_nh(:)
     elseif(SNODEP_obs_obj(n)%gridspan.eq.2) then 
        comb_data(1:SNODEP_obs_obj(n)%mo2) = go_sh(:)
        comb_data(SNODEP_obs_obj(n)%mo2+1:SNODEP_obs_obj(n)%mo1+SNODEP_obs_obj(n)%mo2) = LIS_rc%udef
     else
        comb_data(1:SNODEP_obs_obj(n)%mo2) = go_sh(:)
        comb_data(SNODEP_obs_obj(n)%mo2+1:SNODEP_obs_obj(n)%mo1+SNODEP_obs_obj(n)%mo2) = go_nh(:)
     endif
     varfield = RESHAPE(comb_data(1:SNODEP_obs_obj(n)%mo1+SNODEP_obs_obj(n)%mo2),(/LIS_rc%obs_lnc(n),LIS_rc%obs_lnr(n)/))

     if(SNODEP_obs_obj(n)%gridspan.eq.1) then 
        deallocate(lo_nh)
        deallocate(go_nh)
     elseif(SNODEP_obs_obj(n)%gridspan.eq.2) then 
        deallocate(lo_sh)
        deallocate(go_sh)
     else
        deallocate(lo_nh)
        deallocate(lo_sh)
        deallocate(go_nh)
        deallocate(go_sh)
     endif
     deallocate(comb_data)     
  endif

end subroutine interp_SNODEPfield
