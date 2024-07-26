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
! !ROUTINE: read_ARMdata
! \label{read_ARMdata}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_ARMdata(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_fileIOMod,      only : LIS_readData
  use LIS_timeMgrMod,     only : LIS_calendar, LIS_tick
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use map_utils
  use ARMdata_module,     only : ARMdata_struc

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
  
  integer                   :: n 
  integer                   :: yr, mo, da, hr, mn, ss,doy
  real*8                    :: lis_prevtime
  integer                   :: status
  character(len=LIS_CONST_PATH_LEN) :: filename
  real                      :: time, gmt
  integer                   :: i,t,gid,c,r,st,et
  type(ESMF_TimeInterval)   :: dayInterval
  type(ESMF_Time)           :: initTime
  type(ESMF_Field)          :: qleField
  type(ESMF_Field)          :: qhField
  type(ESMF_Field)          :: qgField
  type(ESMF_Field)          :: sfsmField
  type(ESMF_Field)          :: sfstField
  real                      :: col, row
  integer                   :: stn_col, stn_row
  type(ESMF_Time)           :: armtime1, armtime2
  real, allocatable             :: qle(:)
  real, allocatable             :: qh(:)
  real, allocatable             :: qg(:)
  real, allocatable             :: sfsm(:)
  real, allocatable             :: sfst(:)
  integer, allocatable          :: nqle(:)
  integer, allocatable          :: nqh(:)
  integer, allocatable          :: nqg(:)
  integer, allocatable          :: nsfsm(:)
  integer, allocatable          :: nsfst(:)
  real,    pointer          :: qle_obs(:)
  real,    pointer          :: qh_obs(:)
  real,    pointer          :: qg_obs(:)
  real,    allocatable          :: sfsm_obs(:)
  real,    allocatable          :: sfst_obs(:)

  n = 1

  allocate(qle(LIS_rc%ngrid(n)))
  allocate(qh(LIS_rc%ngrid(n)))
  allocate(qg(LIS_rc%ngrid(n)))
  allocate(sfsm(LIS_rc%ngrid(n)))
  allocate(sfst(LIS_rc%ngrid(n)))

  allocate(nqle(LIS_rc%ngrid(n)))
  allocate(nqh(LIS_rc%ngrid(n)))
  allocate(nqg(LIS_rc%ngrid(n)))
  allocate(nsfsm(LIS_rc%ngrid(n)))
  allocate(nsfst(LIS_rc%ngrid(n)))

  qle = 0 
  qh = 0 
  qg = 0 
  sfsm = 0 
  sfst = 0 

  nqle = 0 
  nqh = 0 
  nqg = 0 
  nsfsm = 0 
  nsfst = 0 

  time = LIS_rc%hr*3600+LIS_rc%mn*60+LIS_rc%ss

  if((mod(time,86400.0).eq.0.0).or.(LIS_rc%da.ne.ARMdata_struc(n)%da)) then 

     call ESMF_TimeSet(ARMdata_struc(n)%starttime, yy=LIS_rc%yr, &
          mm = LIS_rc%mo, &
          dd = LIS_rc%da, &
          h = 0, &
          m = 0, &
          calendar = LIS_calendar, &
          rc=status)
     call LIS_verify(status,'error in timeset: read_ARMdata')
     ARMdata_struc(n)%da = LIS_rc%da

     call process_ecor_flux_data(n,2,LIS_rc%yr,LIS_rc%mo,LIS_rc%da)     
     call process_baebbr_flux_data(n,2,LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
     call process_ebbr_data(n,2,LIS_rc%yr, LIS_rc%mo, LIS_rc%da)

     if(.not.ARMdata_struc(n)%startflag) then 
        if(ARMdata_struc(n)%baebbr_select.eq.1) then 
           ARMdata_struc(n)%baebbr_qle_p = ARMdata_struc(n)%baebbr_qle_c
           ARMdata_struc(n)%baebbr_qh_p = ARMdata_struc(n)%baebbr_qh_c
           ARMdata_struc(n)%baebbr_qg_p = ARMdata_struc(n)%baebbr_qg_c
        endif
        if(ARMdata_struc(n)%ecor_select.eq.1) then 
           ARMdata_struc(n)%ecor_qle_p = ARMdata_struc(n)%ecor_qle_c
           ARMdata_struc(n)%ecor_qh_p = ARMdata_struc(n)%ecor_qh_c
        endif
        if(ARMdata_struc(n)%ebbr_select.eq.1) then 
           ARMdata_struc(n)%ebbr_sfsm_p = ARMdata_struc(n)%ebbr_sfsm_c
           ARMdata_struc(n)%ebbr_sfst_p = ARMdata_struc(n)%ebbr_sfst_c
        endif
     else !read first bookend
        ARMdata_struc(n)%startflag = .false. 
        
        yr = LIS_rc%yr
        mo = LIS_rc%mo
        da = LIS_rc%da
        hr = LIS_rc%hr
        mn = LIS_rc%mn
        ss = LIS_rc%ss
        
        !set back by one day. 
        call LIS_tick(lis_prevtime, doy, gmt, yr,mo,da,hr,mn,ss,-1.0*86400)

        call process_ecor_flux_data(n,1,yr,mo,da)
        call process_baebbr_flux_data(n,1,yr,mo,da)
        call process_ebbr_data(n,1,yr,mo,da)

     end if
  endif

  call ESMF_TimeSet(armtime1,yy=LIS_rc%yr, &
       mm = LIS_rc%mo, &
       dd = LIS_rc%da, &
       h = LIS_rc%hr, &
       m = LIS_rc%mn, &
       calendar = LIS_calendar, &
       rc=status)
  call LIS_verify(status, 'error in timeset: readARMobs')
  yr = LIS_rc%yr
  mo = LIS_rc%mo
  da = LIS_rc%da
  hr = LIS_rc%hr
  mn = LIS_rc%mn
  ss = LIS_rc%ss
  
  call LIS_tick(lis_prevtime, doy, gmt, yr,mo,da,hr,mn,ss,(-1)*LIS_rc%ts)
  call ESMF_TimeSet(armtime2, yy=yr, &
       mm = mo, &
       dd = da, &
       h = hr, &
       m = mn, &
       calendar = LIS_calendar, &
       rc=status)

  call LIS_verify(status, 'error in timeset: readARMobs')

  et = nint((armtime1 - ARMdata_struc(n)%starttime)/&
       ARMdata_struc(n)%baebbr_ts) + 1
  st = nint((armtime2 - ARMdata_struc(n)%starttime)/&
       ARMdata_struc(n)%baebbr_ts) + 1

  if(mod(time,1800.0).eq.0.0) then 
     if(st.gt.0) then 
        call tavg_ecor_flux_data(n,st+1,et,2,qle,nqle,qh,nqh)
     else
        !now add previous day
        call ESMF_TimeIntervalSet(dayInterval,s=86400,rc=status)
        call LIS_verify(status,'error in timeintervalset: readARMobs')
        initTime = ARMdata_struc(n)%startTime - dayInterval
        st = nint((armtime2 - initTime)/ARMdata_struc(n)%baebbr_ts) + 1
        et = 48
        call tavg_ecor_flux_data(n,st,et,1,qle,nqle,qh,nqh)
     endif
  else
     qle = LIS_rc%udef
     nqle = 0 
     qh = LIS_rc%udef
     nqh = 0
  endif

!process the flux data only at 1/2 hour intervals
  if(mod(time,1800.0).eq.0.0) then 
     if(st.gt.0) then 
        call tavg_baebbr_flux_data(n,st+1,et,2,qle,nqle,qh,nqh,qg,nqg)
     else
        !now add previous day
        call ESMF_TimeIntervalSet(dayInterval,s=86400,rc=status)
        call LIS_verify(status,'error in timeintervalset: readARMobs')
        initTime = ARMdata_struc(n)%startTime - dayInterval
        st = nint((armtime2 - initTime)/ARMdata_struc(n)%baebbr_ts) + 1
        et = 48
        call tavg_baebbr_flux_data(n,st,et,1,qle,nqle,qh,nqh,qg,nqg)
     endif
  else
     qle = LIS_rc%udef
     nqle = 0 
     qh = LIS_rc%udef
     nqh = 0
     qg = LIS_rc%udef
     nqg = 0 
  endif

  et = nint((armtime1 - ARMdata_struc(n)%starttime)/ARMdata_struc(n)%ebbr_ts) +1
  st = nint((armtime2 - ARMdata_struc(n)%starttime)/ARMdata_struc(n)%ebbr_ts) +1

  if(mod(time,1800.0).eq.0.0) then 
     if(st.gt.0) then 
        call tavg_ebbr_data(n,st+1,et,2,sfsm,nsfsm,sfst,nsfst)
     else
        !now add previous day
        call ESMF_TimeIntervalSet(dayInterval,s=86400,rc=status)
        call LIS_verify(status,'error in timeintervalset: readARMobs')
        initTime = ARMdata_struc(n)%startTime - dayInterval
        st = nint((armtime2 - initTime)/ARMdata_struc(n)%baebbr_ts) + 1
        et = 48
        call tavg_ebbr_data(n,st,et,1,sfsm,nsfsm,sfst,nsfst)
     endif
  else
     sfsm = LIS_rc%udef
     nsfsm = 0
     sfst = LIS_rc%udef
     nsfst = 0
  endif
 
  do t=1,LIS_rc%ngrid(n)
     if(nqle(t).gt.0) then 
        qle(t) = qle(t)/nqle(t)
     else
        qle(t) = LIS_rc%udef
     endif
  enddo
  
  do t=1,LIS_rc%ngrid(n)
     if(nqh(t).gt.0) then 
        qh(t) = qh(t)/nqh(t)
     else
        qh(t) = LIS_rc%udef
     endif
  enddo

  do t=1,LIS_rc%ngrid(n)
     if(nqg(t).gt.0) then 
        qg(t) = qg(t)/nqg(t)
     else
        qg(t) = LIS_rc%udef
     endif
  enddo
  do t=1,LIS_rc%ngrid(n)
     if(nsfsm(t).gt.0) then 
        sfsm(t) = sfsm(t)/nsfsm(t)
     else
        sfsm(t) = LIS_rc%udef
     endif
  enddo

  do t=1,LIS_rc%ngrid(n)
     if(nsfst(t).gt.0) then 
        sfst(t) = sfst(t)/nsfst(t)
     else
        sfst(t) = LIS_rc%udef
     endif
  enddo

  call ESMF_StateGet(Obj_Space,"ARM Latent Heat Flux",qleField,&
       rc=status)
  call LIS_verify(status, "StateGet failed for ARM Latent Heat Flux in read_ARMdata")
  call ESMF_FieldGet(qleField,localDE=0,farrayPtr=qle_obs,rc=status)
  call LIS_verify(status, "FieldGet failed for ARM Latent Heat Flux in read_ARMdata") 
  qle_obs = qle

  call ESMF_StateGet(Obj_Space,"ARM Sensible Heat Flux",qhField,&
       rc=status)
  call LIS_verify(status, "StateGet failed for ARM Sensible Heat Flux in read_ARMdata")
  call ESMF_FieldGet(qhField,localDE=0,farrayPtr=qh_obs,rc=status)
  call LIS_verify(status, "FieldGet failed for ARM Sensible Heat Flux in read_ARMdata") 
  qh_obs = qh

  call ESMF_StateGet(Obj_Space,"ARM Ground Heat Flux",qgField,&
       rc=status)
  call LIS_verify(status, "StateGet failed for ARM Ground Heat Flux in read_ARMdata")
  call ESMF_FieldGet(qgField,localDE=0,farrayPtr=qg_obs,rc=status)
  call LIS_verify(status, "FieldGet failed for ARM Ground Heat Flux in read_ARMdata") 
  qg_obs = qg

!  call ESMF_StateGet(Obj_Space,"Surface Soil Moisture",sfsmField,&
!       rc=status)
!  call LIS_verify(status, "StateGet failed for SFSM in read_ARMdata")
!  call ESMF_FieldGet(sfsmField,localDE=0,farrayPtr=sfsm_obs,rc=status)
!  call LIS_verify(status, "FieldGet failed for SFSM in read_ARMdata") 
!  sfsm_obs = sfsm

!  call ESMF_StateGet(Obj_Space,"Surface Soil Temperature",sfstField,&
!       rc=status)
!  call LIS_verify(status, "StateGet failed for SFST in read_ARMdata")
!  call ESMF_FieldGet(sfstField,localDE=0,farrayPtr=sfst_obs,rc=status)
!  call LIS_verify(status, "FieldGet failed for SFST in read_ARMdata") 
!  sfst_obs = sfst

  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status)

  deallocate(qle)
  deallocate(qh)
  deallocate(qg)
  deallocate(sfsm)
  deallocate(sfst)

  deallocate(nqle)
  deallocate(nqh)
  deallocate(nqg)
  deallocate(nsfsm)
  deallocate(nsfst)

end subroutine read_ARMdata

!BOP
!
! !ROUTINE: process_ebbr_data
!  \label{process_ebbr_data}
!
! !INTERFACE:
subroutine process_ebbr_data(n,tindex,yr,mo,da)
! !USES:
  use LIS_coreMod,     only : LIS_rc
  use LIS_logMod,      only : LIS_logunit
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use ARMdata_module,     only : ARMdata_struc

  implicit none
! !ARGUMENTS: 
  integer          :: n 
  integer          :: tindex
  integer          :: yr
  integer          :: mo
  integer          :: da
!
! !DESCRIPTION: 
!  This subroutine reads the bulk aerodynamic estimates of 
!  sensible,latent and ground heat flux
!  measurements from the Energy Balance Bowen Ratio (BAEBBR) 
!  system. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[tindex]  value indicating the bookend to be read
!                  into (2-current, 1-previous)     
!   \item[yr]      year for which the data is being processed
!   \item[mo]      month for which the data is being processed
!   \item[da]      day for which the data is being processed
!  \end{description}
!EOP
  integer          :: i 
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer          :: status

  if(ARMdata_struc(n)%ebbr_select.eq.1) then 
  ! Read the baebbr data
     do i=1,ARMdata_struc(n)%n_stns
        call create_arm_ebbr_filename(ARMdata_struc(n)%odir, &
             ARMdata_struc(n)%site_id, ARMdata_struc(n)%stn_name(i), &
             yr, mo, da, filename,status)
        if(status.eq.0) then 
           write(LIS_logunit,*) 'Reading ebbr file ',trim(filename)
           call read_ebbr_file(n,i,yr, mo, da, &
                filename,tindex)
        else
           if(tindex.eq.2) then 
              ARMdata_struc(n)%ebbr_sfsm_c(i,:) = LIS_rc%udef
              ARMdata_struc(n)%ebbr_sfst_c(i,:) = LIS_rc%udef
           else
              ARMdata_struc(n)%ebbr_sfsm_p(i,:) = LIS_rc%udef
              ARMdata_struc(n)%ebbr_sfst_p(i,:) = LIS_rc%udef
           endif
        endif
     enddo
  endif

end subroutine process_ebbr_data

!BOP
!
! !ROUTINE: process_baebbr_flux_data
!  \label{process_baebbr_flux_data}
!
! !INTERFACE:
subroutine process_baebbr_flux_data(n,tindex,yr,mo,da)
! !USES:
  use LIS_coreMod,     only : LIS_rc
  use LIS_logMod,      only : LIS_logunit
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use ARMdata_module,     only : ARMdata_struc

  implicit none
! !ARGUMENTS: 
  integer          :: n
  integer          :: tindex
  integer          :: yr
  integer          :: mo
  integer          :: da
!
! !DESCRIPTION: 
!  This subroutine reads the bulk aerodynamic estimates of 
!  sensible,latent and ground heat flux
!  measurements from the Energy Balance Bowen Ratio (BAEBBR) 
!  system. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[tindex]  value indicating the bookend to be read
!                  into (2-current, 1-previous)     
!   \item[yr]      year for which the data is being processed
!   \item[mo]      month for which the data is being processed
!   \item[da]      day for which the data is being processed
!  \end{description}
!EOP
  integer          :: i 
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer          :: status

  if(ARMdata_struc(n)%baebbr_select.eq.1) then 
  ! Read the baebbr data
     do i=1,ARMdata_struc(n)%n_stns
        call create_arm_baebbr_flux_filename(ARMdata_struc(n)%odir, &
             ARMdata_struc(n)%site_id, ARMdata_struc(n)%stn_name(i), &
             yr, mo, da, filename,status)
        if(status.eq.0) then 
           write(LIS_logunit,*) 'Reading baebbr file   ',trim(filename)
           call read_baebbr_flux_file(n,i,yr, mo, da, &
                filename,tindex)
        else
           if(tindex.eq.2) then 
              ARMdata_struc(n)%baebbr_qle_c(i,:) = LIS_rc%udef
              ARMdata_struc(n)%baebbr_qh_c(i,:) = LIS_rc%udef
              ARMdata_struc(n)%baebbr_qg_c(i,:) = LIS_rc%udef
           else
              ARMdata_struc(n)%baebbr_qle_p(i,:) = LIS_rc%udef
              ARMdata_struc(n)%baebbr_qh_p(i,:) = LIS_rc%udef
              ARMdata_struc(n)%baebbr_qg_p(i,:) = LIS_rc%udef
           endif
        endif
     enddo
  endif

end subroutine process_baebbr_flux_data

!BOP
!
! !ROUTINE: process_ecor_flux_data
! \label{process_ecor_flux_data}
!
! !INTERFACE:
subroutine process_ecor_flux_data(n,tindex,yr,mo,da)
  ! !USES:
  use LIS_coreMod,     only : LIS_rc
  use LIS_logMod,      only : LIS_logunit
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use ARMdata_module,     only : ARMdata_struc

  implicit none

  integer          :: n 
  integer          :: tindex
  integer          :: yr
  integer          :: mo
  integer          :: da

!
! !DESCRIPTION: 
!  This subroutine reads the sensible and latent heat flux
!  measurements from the Eddy Correlation (ECOR) flux 
!  measurement system.
!
!  The arguments are: 
!  \begin{description}
!   \item[tindex]  value indicating the bookend to be read
!                  into (2-current, 1-previous)     
!   \item[yr]      year for which the data is being processed
!   \item[mo]      month for which the data is being processed
!   \item[da]      day for which the data is being processed
!  \end{description}
!EOP 
  integer          :: i 
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer          :: status

  if(ARMdata_struc(n)%ecor_select.eq.1) then 
  ! read the ecor data
     do i=1,ARMdata_struc(n)%n_stns
        call create_arm_ecor_flux_filename(ARMdata_struc(n)%odir, &
             ARMdata_struc(n)%site_id, ARMdata_struc(n)%stn_name(i), &
             yr, mo, da, filename,status)
        if(status.eq.0) then 
           write(LIS_logunit,*) 'Reading ecor file  ',trim(filename)
           !exclude E14 ECOR data. 
           if(trim(ARMdata_struc(n)%stn_name(i)).ne."E14") then 
              call read_ecor_flux_file(n,i, yr, mo, da, &
                   filename,tindex)
           endif
        else
           if(tindex.eq.2) then
              ARMdata_struc(n)%ecor_qh_c(i,:) = LIS_rc%udef
              ARMdata_struc(n)%ecor_qle_c(i,:) = LIS_rc%udef
           else
              ARMdata_struc(n)%ecor_qh_p(i,:) = LIS_rc%udef
              ARMdata_struc(n)%ecor_qle_p(i,:) = LIS_rc%udef
           endif
        endif
     enddo
  endif
end subroutine process_ecor_flux_data

!BOP
! !ROUTINE: tavg_ebbr_data
! \label{tavg_ebbr_data}
! 
! !INTERFACE: 
subroutine tavg_ebbr_data(n,st,et,tindex,sfsm,nsfsm,sfst,nsfst)
! !USES:
  use LIS_coreMod,     only : LIS_rc, LIS_domain
  use LIS_logMod,      only : LIS_logunit
  use ARMdata_module,     only : ARMdata_struc
  use map_utils
  
  implicit none

! !ARGUMENTS: 
  integer          :: n 
  integer          :: st,et
  integer          :: tindex
  real             :: sfsm(LIS_rc%ngrid(n))
  integer          :: nsfsm(LIS_rc%ngrid(n))
  real             :: sfst(LIS_rc%ngrid(n))
  integer          :: nsfst(LIS_rc%ngrid(n))

!
! !DESCRIPTION: 
!  This subroutine temporally aggregates the BAEBBR flux measurement
!  data up to the LIS output timestep 
!
!  The arguments are: 
!  \begin{description}
!   \item[st]  starting time index                
!   \item[et]  ending time index   
!   \item[tindex]  value indicating the bookend to be read
!                  into (2-current, 1-previous)        
!   \item[qle]  array containting aggregated latent heat flux  
!   \item[nqle] array containting temporal counts of latent heat flux  
!   \item[qh]   array containting aggregated sensible heat flux    
!   \item[nqh]  array containting temporal counts of sensible heat flux  

!  \end{description}
!EOP

  integer          :: c,t,gid,stn_col,stn_row
  real             :: col,row
  
  if(ARMdata_struc(n)%ebbr_select.eq.1) then 
     do t=st,et
        do c=1,ARMdata_struc(n)%n_stns
           call latlon_to_ij(LIS_domain(n)%lisproj, ARMdata_struc(n)%stnlat(c), &
                ARMdata_struc(n)%stnlon(c), col, row)
           stn_col = nint(col)
           stn_row = nint(row)
           
           gid = -1
           if(stn_col.ge.1.and.stn_col.le.LIS_rc%lnc(n).and.&
                stn_row.ge.1.and.stn_row.le.LIS_rc%lnr(n)) then 
              gid = LIS_domain(n)%gindex(stn_col,stn_row)
           endif
           
           if(tindex.eq.2.and.gid.ne.-1) then  
              if(ARMdata_struc(n)%ebbr_sfsm_c(c,t).ne.LIS_rc%udef) then 
                 sfsm(gid) = sfsm(gid) + & 
                      ARMdata_struc(n)%ebbr_sfsm_c(c,t) 
                 nsfsm(gid) = nsfsm(gid) + 1
              endif
              if(ARMdata_struc(n)%ebbr_sfst_c(c,t).ne.LIS_rc%udef) then 
                 sfst(gid) = sfst(gid) + & 
                      ARMdata_struc(n)%ebbr_sfst_c(c,t) 
                 nsfst(gid) = nsfst(gid) + 1
              endif
           elseif(tindex.eq.1.and.gid.ne.-1) then 
              if(ARMdata_struc(n)%ebbr_sfsm_p(c,t).ne.LIS_rc%udef) then 
                 sfsm(gid) = sfsm(gid) + & 
                      ARMdata_struc(n)%ebbr_sfsm_p(c,t) 
                 nsfsm(gid) = nsfsm(gid) + 1
              endif
              if(ARMdata_struc(n)%ebbr_sfst_p(c,t).ne.LIS_rc%udef) then 
                 sfst(gid) = sfst(gid) + & 
                      ARMdata_struc(n)%ebbr_sfst_p(c,t) 
                 nsfst(gid) = nsfst(gid) + 1
              endif
           endif
        enddo
     enddo
  endif
end subroutine tavg_ebbr_data
   
!BOP
! !ROUTINE: tavg_baebbr_flux_data
! \label{tavg_baebbr_flux_data}
! 
! !INTERFACE: 
subroutine tavg_baebbr_flux_data(n,st,et,tindex,qle,nqle,qh,nqh,qg,nqg)
! !USES:
  use LIS_coreMod,     only : LIS_rc, LIS_domain
  use ARMdata_module,     only : ARMdata_struc
  use LIS_logMod,      only : LIS_logunit
  use map_utils

  implicit none

! !ARGUMENTS: 
  integer          :: n 
  integer          :: st,et
  integer          :: tindex
  real             :: qle(LIS_rc%ngrid(n))
  integer          :: nqle(LIS_rc%ngrid(n))
  real             :: qh(LIS_rc%ngrid(n))
  integer          :: nqh(LIS_rc%ngrid(n))
  real             :: qg(LIS_rc%ngrid(n))
  integer          :: nqg(LIS_rc%ngrid(n))
!
! !DESCRIPTION: 
!  This subroutine temporally aggregates the BAEBBR flux measurement
!  data up to the LIS output timestep 
!
!  The arguments are: 
!  \begin{description}
!   \item[st]  starting time index                
!   \item[et]  ending time index   
!   \item[tindex]  value indicating the bookend to be read
!                  into (2-current, 1-previous)        
!   \item[qle]  array containting aggregated latent heat flux  
!   \item[nqle] array containting temporal counts of latent heat flux  
!   \item[qh]   array containting aggregated sensible heat flux    
!   \item[nqh]  array containting temporal counts of sensible heat flux  
!   \item[qg]   array containting aggregated ground heat flux  
!   \item[nqg]  array containting temporal counts of ground heat flux          
!  \end{description}
!EOP

  integer          :: c,t,stn_col,stn_row,gid
  real             :: col,row
  
  if(ARMdata_struc(n)%baebbr_select.eq.1) then 
     do t=st,et
        do c=1,ARMdata_struc(n)%n_stns
           call latlon_to_ij(LIS_domain(n)%lisproj, ARMdata_struc(n)%stnlat(c), &
                ARMdata_struc(n)%stnlon(c), col, row)
           stn_col = nint(col)
           stn_row = nint(row)
           
           gid = -1
           if(stn_col.ge.1.and.stn_col.le.LIS_rc%lnc(n).and.&
                stn_row.ge.1.and.stn_row.le.LIS_rc%lnr(n)) then 
              gid = LIS_domain(n)%gindex(stn_col,stn_row)
           endif
           
           if(tindex.eq.2.and.gid.ne.-1) then      
              if(ARMdata_struc(n)%baebbr_qle_c(c,t).ne.LIS_rc%udef) then 
                 
                 qle(gid) = qle(gid) + & 
                      (-1)*ARMdata_struc(n)%baebbr_qle_c(c,t) 
                 nqle(gid) = nqle(gid) + 1
              endif
              
              if(ARMdata_struc(n)%baebbr_qh_c(c,t).ne.LIS_rc%udef) then 
                 qh(gid) = qh(gid) + & 
                      (-1)*ARMdata_struc(n)%baebbr_qh_c(c,t) 
                 nqh(gid) = nqh(gid) + 1
              endif
              
              if(ARMdata_struc(n)%baebbr_qg_c(c,t).ne.LIS_rc%udef) then 
                 qg(gid) = qg(gid) + & 
                      ARMdata_struc(n)%baebbr_qg_c(c,t) 
                 nqg(gid) = nqg(gid) + 1
              endif
           elseif(tindex.eq.1.and.gid.ne.-1) then 
              if(ARMdata_struc(n)%baebbr_qle_p(c,t).ne.LIS_rc%udef) then 
                 
                 qle(gid) = qle(gid) + & 
                      (-1)*ARMdata_struc(n)%baebbr_qle_p(c,t) 
                 nqle(gid) = nqle(gid) + 1
              endif
              
              if(ARMdata_struc(n)%baebbr_qh_p(c,t).ne.LIS_rc%udef) then 
                 qh(gid) = qh(gid) + & 
                      (-1)*ARMdata_struc(n)%baebbr_qh_p(c,t) 
                 nqh(gid) = nqh(gid) + 1
              endif
              
              if(ARMdata_struc(n)%baebbr_qg_p(c,t).ne.LIS_rc%udef) then 
                 qg(gid) = qg(gid) + & 
                      ARMdata_struc(n)%baebbr_qg_p(c,t) 
                 nqg(gid) = nqg(gid) + 1
              endif
           endif
        enddo
     enddo
  endif
end subroutine tavg_baebbr_flux_data
   
!BOP
! !ROUTINE: tavg_ecor_flux_data
! \label{tavg_ecor_flux_data}
! 
! !INTERFACE: 
subroutine tavg_ecor_flux_data(n,st,et,tindex,qle,nqle,qh,nqh)
! !USES:
  use LIS_coreMod,     only : LIS_rc, LIS_domain
  use ARMdata_module,     only : ARMdata_struc
  use LIS_logMod,      only : LIS_logunit
  use map_utils
  
  implicit none

! !ARGUMENTS: 
  integer          :: n 
  integer          :: st,et
  integer          :: tindex
  real             :: qle(LIS_rc%ngrid(n))
  integer          :: nqle(LIS_rc%ngrid(n))
  real             :: qh(LIS_rc%ngrid(n))
  integer          :: nqh(LIS_rc%ngrid(n))

!
! !DESCRIPTION: 
!  This subroutine temporally aggregates the ECOR flux measurement
!  data up to the LIS output timestep 
!
!  The arguments are: 
!  \begin{description}
!   \item[st]  starting time index                
!   \item[et]  ending time index   
!   \item[tindex]  value indicating the bookend to be read
!                  into (2-current, 1-previous)        
!   \item[qle]  array containting aggregated latent heat flux  
!   \item[nqle] array containting temporal counts of latent heat flux  
!   \item[qh]   array containting aggregated sensible heat flux    
!   \item[nqh]  array containting temporal counts of sensible heat flux  
!  \end{description}
!EOP

  integer          :: c,t,stn_col,stn_row,gid
  real             :: col,row

  if(ARMdata_struc(n)%ecor_select.eq.1) then 
     do t=st,et
        do c=1,ARMdata_struc(n)%n_stns
           call latlon_to_ij(LIS_domain(n)%lisproj, ARMdata_struc(n)%stnlat(c), &
                ARMdata_struc(n)%stnlon(c), col, row)
           stn_col = nint(col)
           stn_row = nint(row)
           
           gid = -1
           if(stn_col.ge.1.and.stn_col.le.LIS_rc%lnc(n).and.&
                stn_row.ge.1.and.stn_row.le.LIS_rc%lnr(n)) then 
              gid = LIS_domain(n)%gindex(stn_col,stn_row)
           endif

           if(tindex.eq.2.and.gid.ne.-1) then 
              if(ARMdata_struc(n)%ecor_qle_c(c,t).ne.LIS_rc%udef) then 
                 qle(gid) = qle(gid) + & 
                      ARMdata_struc(n)%ecor_qle_c(c,t) 
                 nqle(gid) = nqle(gid) + 1
              endif
              
              if(ARMdata_struc(n)%ecor_qh_c(c,t).ne.LIS_rc%udef) then 
                 qh(gid) = qh(gid) + & 
                      ARMdata_struc(n)%ecor_qh_c(c,t) 
                 nqh(gid) = nqh(gid) + 1
              endif
           elseif(tindex.eq.1.and.gid.ne.-1) then 
              if(ARMdata_struc(n)%ecor_qle_p(c,t).ne.LIS_rc%udef) then 
                 qle(gid) = qle(gid) + & 
                      ARMdata_struc(n)%ecor_qle_p(c,t) 
                 nqle(gid) = nqle(gid) + 1
              endif
              
              if(ARMdata_struc(n)%ecor_qh_p(c,t).ne.LIS_rc%udef) then 
                 qh(gid) = qh(gid) + & 
                      ARMdata_struc(n)%ecor_qh_p(c,t) 
                 nqh(gid) = nqh(gid) + 1
              endif
              
           endif
        enddo
     enddo
  endif
end subroutine tavg_ecor_flux_data

!BOP
! 
! !ROUTINE: read_ebbr_file
! \label{read_ebbr_file}
! 
! !INTERFACE: 
subroutine read_ebbr_file(n,k, yr, mo, da, filename,bend)
! !USES: 
  use ESMF     
  use LIS_coreMod,      only : LIS_rc
  use LIS_logMod,       only : LIS_logunit, LIS_verify
  use LIS_timeMgrMod,   only : LIS_calendar
  use ARMdata_module,     only : ARMdata_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer,  intent(in) :: n
  integer,  intent(in) :: k 
  integer              :: yr
  integer              :: mo
  integer              :: da
  character(len=*)     :: filename 
  integer              :: bend

!DESCRIPTION: 
! 
!  The arguments are: 
!  \begin{description}
!   \item[k]         Index of the ARM station
!   \item[yr]        year of ARM data 
!   \item[mo]        month of ARM data
!   \item[da]        day of ARM data
!   \item[filename]  Name of the BAEBBR file
!   \item[bend]      bookend index (2-current, 1-previous)
!  \end{description}
!EOP
  integer       :: nid, btimeid, latid, lonid
  integer       :: timeid, dimId
  integer       :: sm1id,sm2id,sm3id,sm4id,sm5id
  integer       :: qc_sm1id,qc_sm2id,qc_sm3id,qc_sm4id,qc_sm5id
  integer       :: st1id,st2id,st3id,st4id,st5id
  integer       :: qc_st1id,qc_st2id,qc_st3id,qc_st4id,qc_st5id
  integer       :: ndims
  real          :: lat, lon 
  integer       :: base_time
  integer       :: ios
  real, allocatable :: time(:)
  real, allocatable :: ebbr_sm1(:),ebbr_st1(:)
  real, allocatable :: ebbr_sm2(:),ebbr_st2(:)
  real, allocatable :: ebbr_sm3(:),ebbr_st3(:)
  real, allocatable :: ebbr_sm4(:),ebbr_st4(:)
  real, allocatable :: ebbr_sm5(:),ebbr_st5(:)
  
  integer, allocatable :: ebbr_qc_sm1(:),ebbr_qc_st1(:)
  integer, allocatable :: ebbr_qc_sm2(:),ebbr_qc_st2(:)
  integer, allocatable :: ebbr_qc_sm3(:),ebbr_qc_st3(:)
  integer, allocatable :: ebbr_qc_sm4(:),ebbr_qc_st4(:)
  integer, allocatable :: ebbr_qc_sm5(:),ebbr_qc_st5(:)
  logical           :: qc_sm,qc_st
  integer           :: kk
  type(ESMF_Time)   :: reftime, datatime, currtime
  type(ESMF_TimeInterval) :: dt
  real              :: sfsm,sfst
  integer           :: yr1, mo1, da1, hr1, mn1, ss1
  integer           :: status
  integer           :: data_index

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
  call LIS_verify(ios, 'Error opening file'//trim(filename))

!variable ids
  ios = nf90_inq_varid(nid, 'base_time',btimeid)
  call LIS_verify(ios, 'Error nf90_inq_varid: base_time')

  ios = nf90_inq_varid(nid, 'lat',latid)
  call LIS_verify(ios, 'Error nf90_inq_varid: lat')

  ios = nf90_inq_varid(nid, 'lon',lonid)
  call LIS_verify(ios, 'Error nf90_inq_varid: lon')

  ios = nf90_inq_varid(nid, 'time_offset',timeid)
  call LIS_verify(ios, 'Error nf90_inq_varid: time')

  ios = nf90_inq_varid(nid, 'soil_moisture_1',sm1id)
  call LIS_verify(ios, 'Error nf90_inq_varid: sm1')

  ios = nf90_inq_varid(nid, 'qc_soil_moisture_1',qc_sm1id)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_sm1')

  ios = nf90_inq_varid(nid, 'soil_moisture_2',sm2id)
  call LIS_verify(ios, 'Error nf90_inq_varid: sm2')

  ios = nf90_inq_varid(nid, 'qc_soil_moisture_2',qc_sm2id)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_sm2')

  ios = nf90_inq_varid(nid, 'soil_moisture_3',sm3id)
  call LIS_verify(ios, 'Error nf90_inq_varid: sm3')

  ios = nf90_inq_varid(nid, 'qc_soil_moisture_3',qc_sm3id)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_sm3')

  ios = nf90_inq_varid(nid, 'soil_moisture_4',sm4id)
  call LIS_verify(ios, 'Error nf90_inq_varid: sm4')

  ios = nf90_inq_varid(nid, 'qc_soil_moisture_4',qc_sm4id)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_sm4')

  ios = nf90_inq_varid(nid, 'soil_moisture_5',sm5id)
  call LIS_verify(ios, 'Error nf90_inq_varid: sm5')

  ios = nf90_inq_varid(nid, 'qc_soil_moisture_5',qc_sm5id)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_sm5')

  ios = nf90_inq_varid(nid, 'soil_temp_1',st1id)
  call LIS_verify(ios, 'Error nf90_inq_varid: st1')

  ios = nf90_inq_varid(nid, 'qc_soil_temp_1',qc_st1id)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_st1')

  ios = nf90_inq_varid(nid, 'soil_temp_2',st2id)
  call LIS_verify(ios, 'Error nf90_inq_varid: st2')

  ios = nf90_inq_varid(nid, 'qc_soil_temp_2',qc_st2id)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_st2')

  ios = nf90_inq_varid(nid, 'soil_temp_3',st3id)
  call LIS_verify(ios, 'Error nf90_inq_varid: st3')

  ios = nf90_inq_varid(nid, 'qc_soil_temp_3',qc_st3id)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_st3')

  ios = nf90_inq_varid(nid, 'soil_temp_4',st4id)
  call LIS_verify(ios, 'Error nf90_inq_varid: st4')

  ios = nf90_inq_varid(nid, 'qc_soil_temp_4',qc_st4id)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_st4')

  ios = nf90_inq_varid(nid, 'soil_temp_5',st5id)
  call LIS_verify(ios, 'Error nf90_inq_varid: st5')

  ios = nf90_inq_varid(nid, 'qc_soil_temp_5',qc_st5id)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_st5')
!dimensions
  ios = nf90_inq_dimid(nid, 'time',dimId)
  call LIS_verify(ios, 'Error nf90_inq_dimid: time')

  ios = nf90_inquire_dimension(nid, dimId, len=ndims)
  call LIS_verify(ios, 'Error nf90_inquire_dimension:')

!values
  ios = nf90_get_var(nid,latid, lat)
  call LIS_verify(ios, 'Error nf90_get_var: lat')

  ios = nf90_get_var(nid,lonid, lon)
  call LIS_verify(ios, 'Error nf90_get_var: lon')

  ios = nf90_get_var(nid,btimeid, base_time)
  call LIS_verify(ios, 'Error nf90_get_var: base_time')

  call ESMF_TimeSet(refTime, yy=1970,mm=1,dd=1,&
       h=0,m=0,s=0,calendar=LIS_calendar,rc=status)
  call LIS_verify(status, 'error in timeset: readARMObs')

!  ARMdata_struc(n)%stnlat(k) = lat
!  ARMdata_struc(n)%stnlon(k) = lon

  allocate(time(ndims))
  allocate(ebbr_sm1(ndims))
  allocate(ebbr_sm2(ndims))
  allocate(ebbr_sm3(ndims))
  allocate(ebbr_sm4(ndims))
  allocate(ebbr_sm5(ndims))
  allocate(ebbr_st1(ndims))
  allocate(ebbr_st2(ndims))
  allocate(ebbr_st3(ndims))
  allocate(ebbr_st4(ndims))
  allocate(ebbr_st5(ndims))

  allocate(ebbr_qc_sm1(ndims))
  allocate(ebbr_qc_sm2(ndims))
  allocate(ebbr_qc_sm3(ndims))
  allocate(ebbr_qc_sm4(ndims))
  allocate(ebbr_qc_sm5(ndims))
  allocate(ebbr_qc_st1(ndims))
  allocate(ebbr_qc_st2(ndims))
  allocate(ebbr_qc_st3(ndims))
  allocate(ebbr_qc_st4(ndims))
  allocate(ebbr_qc_st5(ndims))

  ebbr_sm1 = 0
  ebbr_sm2 = 0 
  ebbr_sm3 = 0 
  ebbr_sm4 = 0 
  ebbr_sm5 = 0 

  ebbr_qc_sm1 = -1
  ebbr_qc_sm2 = -1
  ebbr_qc_sm3 = -1
  ebbr_qc_sm4 = -1
  ebbr_qc_sm5 = -1

  ebbr_st1 = 0
  ebbr_st2 = 0 
  ebbr_st3 = 0 
  ebbr_st4 = 0 
  ebbr_st5 = 0 

  ebbr_qc_st1 = -1
  ebbr_qc_st2 = -1
  ebbr_qc_st3 = -1
  ebbr_qc_st4 = -1
  ebbr_qc_st5 = -1


  ios = nf90_get_var(nid,timeid, time)
  call LIS_verify(ios, 'Error nf90_get_var: time')

  ios = nf90_get_var(nid,sm1id, ebbr_sm1)
  call LIS_verify(ios, 'Error nf90_get_var: sm1')

  ios = nf90_get_var(nid,qc_sm1id, ebbr_qc_sm1)
  call LIS_verify(ios, 'Error nf90_get_var: qc_sm1')

  ios = nf90_get_var(nid,sm2id, ebbr_sm2)
  call LIS_verify(ios, 'Error nf90_get_var: sm2')

  ios = nf90_get_var(nid,qc_sm2id, ebbr_qc_sm2)
  call LIS_verify(ios, 'Error nf90_get_var: qc_sm2')

  ios = nf90_get_var(nid,sm3id, ebbr_sm3)
  call LIS_verify(ios, 'Error nf90_get_var: sm3')

  ios = nf90_get_var(nid,qc_sm3id, ebbr_qc_sm3)
  call LIS_verify(ios, 'Error nf90_get_var: qc_sm3')

  ios = nf90_get_var(nid,sm4id, ebbr_sm4)
  call LIS_verify(ios, 'Error nf90_get_var: sm4')

  ios = nf90_get_var(nid,qc_sm4id, ebbr_qc_sm4)
  call LIS_verify(ios, 'Error nf90_get_var: qc_sm4')

  ios = nf90_get_var(nid,sm5id, ebbr_sm5)
  call LIS_verify(ios, 'Error nf90_get_var: sm5')

  ios = nf90_get_var(nid,qc_sm5id, ebbr_qc_sm5)
  call LIS_verify(ios, 'Error nf90_get_var: qc_sm5')

  ios = nf90_get_var(nid,qc_st1id, ebbr_qc_st1)
  call LIS_verify(ios, 'Error nf90_get_var: qc_st1')

  ios = nf90_get_var(nid,st2id, ebbr_st2)
  call LIS_verify(ios, 'Error nf90_get_var: st2')

  ios = nf90_get_var(nid,qc_st2id, ebbr_qc_st2)
  call LIS_verify(ios, 'Error nf90_get_var: qc_st2')

  ios = nf90_get_var(nid,st3id, ebbr_st3)
  call LIS_verify(ios, 'Error nf90_get_var: st3')

  ios = nf90_get_var(nid,qc_st3id, ebbr_qc_st3)
  call LIS_verify(ios, 'Error nf90_get_var: qc_st3')

  ios = nf90_get_var(nid,st4id, ebbr_st4)
  call LIS_verify(ios, 'Error nf90_get_var: st4')

  ios = nf90_get_var(nid,qc_st4id, ebbr_qc_st4)
  call LIS_verify(ios, 'Error nf90_get_var: qc_st4')

  ios = nf90_get_var(nid,st5id, ebbr_st5)
  call LIS_verify(ios, 'Error nf90_get_var: st5')

  ios = nf90_get_var(nid,qc_st5id, ebbr_qc_st5)
  call LIS_verify(ios, 'Error nf90_get_var: qc_st5')

  ios = nf90_close(nid)
  call LIS_verify(ios, 'Error in nf90_close')
  
  call ESMF_TimeIntervalSet(dt,s=base_time,rc=status)
  call LIS_verify(status, 'Error in timeintervalset: readARMobs')

  reftime = reftime + dt

  do kk=1,ndims
     call ESMF_TimeIntervalSet(dt,s=nint(time(kk)),rc=status)
     call LIS_verify(status, 'Error in timeintervalset: readARMobs')

     datatime = reftime + dt
     
     call ESMF_TimeGet(datatime, yy=yr1,mm=mo1,dd=da1,&
          h=hr1,m=mn1,s=ss1,calendar=LIS_calendar,rc=status)
     call LIS_verify(status, 'error in timeget: readARMObs')

! data index from 0z with 30mn timestep. All this extensive processing is
! required because some of the ARM files report data at times different 
! from 0z. 
     
     call ESMF_TimeSet(currTime, yy=yr, mm=mo, &
          dd=da, h=0,m=0,s=0,calendar=LIS_calendar,rc=status)
     call LIS_verify(status, 'error in timeget: readARMObs')

     data_index = (datatime - currtime)/ARMdata_struc(n)%ebbr_ts + 1
     ARMdata_struc(n)%ebbr_tindex(k,data_index) = data_index
     if(bend.eq.2) then 
        
        sfsm = LIS_rc%udef
        sfst = LIS_rc%udef
        qc_sm = .true. 
        qc_sm = qc_sm .and. &
             (ebbr_qc_sm1(kk).eq.0).and.&
             (ebbr_qc_sm2(kk).eq.0).and.&
             (ebbr_qc_sm3(kk).eq.0).and.&
             (ebbr_qc_sm4(kk).eq.0).and.&
             (ebbr_qc_sm5(kk).eq.0)
        
        if(qc_sm) then
           sfsm = 0 
           if(ARMdata_struc(n)%stnbd(k).ne.-9999.0) then 
              sfsm = sfsm + &
                   (((ebbr_sm1(kk) + &
                   ebbr_sm2(kk) + &
                   ebbr_sm3(kk) + &
                   ebbr_sm4(kk) + &
                   ebbr_sm5(kk))/5.0)*ARMdata_struc(n)%stnbd(k))/100.0
           else
              sfsm = LIS_rc%udef
           endif
        endif

        qc_st = .true. 
        qc_st = qc_st .and. &
             (ebbr_qc_st1(kk).eq.0).and.&
             (ebbr_qc_st2(kk).eq.0).and.&
             (ebbr_qc_st3(kk).eq.0).and.&
             (ebbr_qc_st4(kk).eq.0).and.&
             (ebbr_qc_st5(kk).eq.0).and.&
             (ebbr_st1(kk).gt.-100).and.&
             (ebbr_st2(kk).gt.-100).and.&
             (ebbr_st3(kk).gt.-100).and.&
             (ebbr_st4(kk).gt.-100).and.&
             (ebbr_st5(kk).gt.-100).and.&
             (ebbr_st1(kk).lt.100).and.&
             (ebbr_st2(kk).lt.100).and.&
             (ebbr_st3(kk).lt.100).and.&
             (ebbr_st4(kk).lt.100).and.&
             (ebbr_st5(kk).lt.100)
        if(qc_st) then
           sfst = 0 
           sfst = sfst + &
                ((ebbr_st1(kk) + & 
                ebbr_st2(kk) + & 
                ebbr_st3(kk) + & 
                ebbr_st4(kk) + & 
                ebbr_st5(kk))/5)+273.15
           if(sfst.lt.0) then 
              print*, ebbr_st1(kk),ebbr_qc_st1(kk),ebbr_st2(kk),&
                   ebbr_st3(kk),ebbr_st4(kk),&
                   ebbr_st5(kk)
              stop
           endif
        else
           sfst = LIS_rc%udef
        endif

        ARMdata_struc(n)%ebbr_sfsm_c(k,data_index)     = sfsm
        ARMdata_struc(n)%ebbr_sfst_c(k,data_index)     = sfst

     elseif(bend.eq.1) then 
        
        sfsm = LIS_rc%udef
        sfst = LIS_rc%udef
        qc_sm = .true. 
        qc_sm = qc_sm .and. &
             (ebbr_qc_sm1(kk).eq.0).and.&
             (ebbr_qc_sm2(kk).eq.0).and.&
             (ebbr_qc_sm3(kk).eq.0).and.&
             (ebbr_qc_sm4(kk).eq.0).and.&
             (ebbr_qc_sm5(kk).eq.0)
        
        if(qc_sm) then
           sfsm = 0 
           if(ARMdata_struc(n)%stnbd(k).ne.-9999.0) then 
              sfsm = sfsm + &
                   (((ebbr_sm1(kk) + &
                   ebbr_sm2(kk) + &
                   ebbr_sm3(kk) + &
                   ebbr_sm4(kk) + &
                   ebbr_sm5(kk))/5.0)*ARMdata_struc(n)%stnbd(k))/100.0
           else
              sfsm = LIS_rc%udef
           endif
        endif
        
        qc_st = .true. 
        qc_st = qc_st .and. &
             (ebbr_qc_st1(kk).eq.0).and.&
             (ebbr_qc_st2(kk).eq.0).and.&
             (ebbr_qc_st3(kk).eq.0).and.&
             (ebbr_qc_st4(kk).eq.0).and.&
             (ebbr_qc_st5(kk).eq.0).and.&
             (ebbr_st1(kk).gt.-100).and.&
             (ebbr_st2(kk).gt.-100).and.&
             (ebbr_st3(kk).gt.-100).and.&
             (ebbr_st4(kk).gt.-100).and.&
             (ebbr_st5(kk).gt.-100).and.&
             (ebbr_st1(kk).lt.100).and.&
             (ebbr_st2(kk).lt.100).and.&
             (ebbr_st3(kk).lt.100).and.&
             (ebbr_st4(kk).lt.100).and.&
             (ebbr_st5(kk).lt.100)
        if(qc_st) then
           sfst = 0 
           sfst = sfst + &
                ((ebbr_st1(kk) + & 
                ebbr_st2(kk) + & 
                ebbr_st3(kk) + & 
                ebbr_st4(kk) + & 
                ebbr_st5(kk))/5)+273.15
           if(sfst.lt.0) then 
              print*, ebbr_st1(kk),ebbr_qc_st1(kk),ebbr_st2(kk),&
                   ebbr_st3(kk),ebbr_st4(kk),&
                   ebbr_st5(kk)
              stop
           endif
        else
           sfst = LIS_rc%udef
        endif
        ARMdata_struc(n)%ebbr_sfsm_p(k,data_index)     = sfsm
        ARMdata_struc(n)%ebbr_sfst_p(k,data_index)     = sfst
     endif
  enddo

  deallocate(time)
  deallocate(ebbr_sm1)
  deallocate(ebbr_sm2)
  deallocate(ebbr_sm3)
  deallocate(ebbr_sm4)
  deallocate(ebbr_sm5)
  deallocate(ebbr_st1)
  deallocate(ebbr_st2)
  deallocate(ebbr_st3)
  deallocate(ebbr_st4)
  deallocate(ebbr_st5)

  deallocate(ebbr_qc_sm1)
  deallocate(ebbr_qc_sm2)
  deallocate(ebbr_qc_sm3)
  deallocate(ebbr_qc_sm4)
  deallocate(ebbr_qc_sm5)
  deallocate(ebbr_qc_st1)
  deallocate(ebbr_qc_st2)
  deallocate(ebbr_qc_st3)
  deallocate(ebbr_qc_st4)
  deallocate(ebbr_qc_st5)

#endif

end subroutine read_ebbr_file


!BOP
! 
! !ROUTINE: read_ecor_flux_file
! \label{read_ecor_flux_file}
! 
! !INTERFACE: 
subroutine read_ecor_flux_file(n, k, yr, mo, da, filename, bend)
! !USES: 
  use ESMF     
  use LIS_coreMod,      only : LIS_rc
  use LIS_logMod,       only : LIS_logunit, LIS_verify
  use LIS_timeMgrMod,   only : LIS_calendar
  use ARMdata_module,   only : ARMdata_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

!DESCRIPTION: 
!  This routine reads the eddy correlation estimates of 
!  latent and sensible fluxes from the ECOR file, 
!  for a particular station. The data is read from the 
!  native NETCDF files, and is processed by applying the 
!  quality control flags. Only data points with 'good'
!  classification (qc value=0) is chosen. This routine 
!  also computes the temporal offset of the data relative
!  to the 0z of a particular day. 
!
!  The arguments are: 
!  \begin{description}
!   \item[k]         Index of the ARM station
!   \item[yr]        year of ARM data 
!   \item[mo]        month of ARM data
!   \item[da]        day of ARM data
!   \item[filename]  Name of the ECOR file
!   \item[bend]      bookend index (2-current, 1-previous)
!  \end{description}
!EOP


  integer,  intent(in) :: k 
  integer,  intent(in) :: n 
  character(len=*)  :: filename 
  integer           :: yr, mo, da
  integer          :: bend

  integer       :: nid, btimeid, latid, lonid, qleid
  integer       :: qhid,timeid, dimId
  integer       :: qc_qhid, qc_qleid
  integer       :: ndims
  real          :: lat, lon 
  integer       :: base_time
  integer       :: ios
  real, allocatable :: time(:)
  real, allocatable :: ecor_qh(:)
  real, allocatable :: ecor_qle(:)
  integer, allocatable :: ecor_qc_qh(:)
  integer, allocatable :: ecor_qc_qle(:)
  integer           :: kk
  type(ESMF_Time)   :: reftime, datatime, currtime
  type(ESMF_TimeInterval) :: dt,ts
  integer           :: yr1, mo1, da1, hr1, mn1, ss1
  integer           :: status
  integer           :: data_index

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
  call LIS_verify(ios, 'Error opening file'//trim(filename))

!variable ids
  ios = nf90_inq_varid(nid, 'base_time',btimeid)
  call LIS_verify(ios, 'Error nf90_inq_varid: base_time')

  ios = nf90_inq_varid(nid, 'lat',latid)
  call LIS_verify(ios, 'Error nf90_inq_varid: lat')

  ios = nf90_inq_varid(nid, 'lon',lonid)
  call LIS_verify(ios, 'Error nf90_inq_varid: lon')

  ios = nf90_inq_varid(nid, 'time_offset',timeid)
  call LIS_verify(ios, 'Error nf90_inq_varid: time')

  ios = nf90_inq_varid(nid, 'h',qhid)
  call LIS_verify(ios, 'Error nf90_inq_varid: qh')

  ios = nf90_inq_varid(nid, 'qc_h',qc_qhid)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_qh')

  ios = nf90_inq_varid(nid, 'lv_e',qleid)
  call LIS_verify(ios, 'Error nf90_inq_varid: qle')

  ios = nf90_inq_varid(nid, 'qc_lv_e',qc_qleid)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_qle')

!dimensions
  ios = nf90_inq_dimid(nid, 'time',dimId)
  call LIS_verify(ios, 'Error nf90_inq_dimid: time')

  ios = nf90_inquire_dimension(nid, dimId, len=ndims)
  call LIS_verify(ios, 'Error nf90_inquire_dimension:')

!values
  ios = nf90_get_var(nid,latid, lat)
  call LIS_verify(ios, 'Error nf90_get_var: lat')

  ios = nf90_get_var(nid,lonid, lon)
  call LIS_verify(ios, 'Error nf90_get_var: lon')

  ios = nf90_get_var(nid,btimeid, base_time)
  call LIS_verify(ios, 'Error nf90_get_var: base_time')

  call ESMF_TimeSet(refTime, yy=1970,mm=1,dd=1,&
       h=0,m=0,s=0,calendar=LIS_calendar,rc=status)
  call LIS_verify(status, 'error in timeset: readARMObs')

!  ARMdata_struc(n)%stnlat(k) = lat
!  ARMdata_struc(n)%stnlon(k) = lon

  allocate(time(ndims))
  allocate(ecor_qle(ndims))
  allocate(ecor_qh(ndims))

  allocate(ecor_qc_qle(ndims))
  allocate(ecor_qc_qh(ndims))

  ios = nf90_get_var(nid,timeid, time)
  call LIS_verify(ios, 'Error nf90_get_var: time')

  ios = nf90_get_var(nid,qhid, ecor_qh)
  call LIS_verify(ios, 'Error nf90_get_var: qh')

  ios = nf90_get_var(nid,qc_qhid, ecor_qc_qh)
  call LIS_verify(ios, 'Error nf90_get_var: qc_qh')

  ios = nf90_get_var(nid,qleid, ecor_qle)
  call LIS_verify(ios, 'Error nf90_get_var: qle')

  ios = nf90_get_var(nid,qc_qleid, ecor_qc_qle)
  call LIS_verify(ios, 'Error nf90_get_var: qc_qle')

  ios = nf90_close(nid)
  call LIS_verify(ios, 'Error in nf90_close')
  
  call ESMF_TimeIntervalSet(dt,s=base_time,rc=status)
  call LIS_verify(status, 'Error in timeintervalset: readARMobs')

  reftime = reftime + dt

  do kk=1,ndims
     call ESMF_TimeIntervalSet(dt,s=nint(time(kk)),rc=status)
     call LIS_verify(status, 'Error in timeintervalset: readARMobs')

     datatime = reftime + dt
     
     call ESMF_TimeGet(datatime, yy=yr1,mm=mo1,dd=da1,&
          h=hr1,m=mn1,s=ss1,calendar=LIS_calendar,rc=status)
     call LIS_verify(status, 'error in timeget: readARMObs')

! data index from 0z with 30mn timestep. All this extensive processing is
! required because some of the ARM files report data at times different 
! from 0z. 
     
     call ESMF_TimeSet(currTime, yy=yr, mm=mo, &
          dd=da, h=0,m=0,s=0,calendar=LIS_calendar,rc=status)
     call LIS_verify(status, 'error in timeget: readARMObs')

     data_index = (datatime - currtime)/ARMdata_struc(n)%ecor_ts + 1
     ARMdata_struc(n)%ecor_tindex(k,data_index) = data_index
     if(bend.eq.2) then 
        if(ecor_qc_qh(kk).eq.0) then 
           ARMdata_struc(n)%ecor_qh_c(k,data_index)     = ecor_qh(kk)
        else
           ARMdata_struc(n)%ecor_qh_c(k,data_index)     = LIS_rc%udef
        endif
        if(ecor_qc_qle(kk).eq.0) then 
           ARMdata_struc(n)%ecor_qle_c(k,data_index)    = ecor_qle(kk)
        else
           ARMdata_struc(n)%ecor_qle_c(k,data_index)    = LIS_rc%udef
        endif
     else
        if(ecor_qc_qh(kk).eq.0) then 
           ARMdata_struc(n)%ecor_qh_p(k,data_index)     = ecor_qh(kk)
        else
           ARMdata_struc(n)%ecor_qh_p(k,data_index)     = LIS_rc%udef
        endif
        if(ecor_qc_qle(kk).eq.0) then 
           ARMdata_struc(n)%ecor_qle_p(k,data_index)    = ecor_qle(kk)
        else
           ARMdata_struc(n)%ecor_qle_p(k,data_index)    = LIS_rc%udef
        endif
     endif     
  enddo

  deallocate(time)
  deallocate(ecor_qle)
  deallocate(ecor_qh)
  deallocate(ecor_qc_qle)
  deallocate(ecor_qc_qh)

#endif

end subroutine read_ecor_flux_file

!BOP
! 
! !ROUTINE: read_baebbr_flux_file
! \label{read_baebbr_flux_file}
! 
! !INTERFACE: 
subroutine read_baebbr_flux_file(n, k, yr, mo, da, filename,bend)
! !USES: 
  use ESMF     
  use LIS_coreMod,      only : LIS_rc
  use LIS_logMod,       only : LIS_logunit, LIS_verify
  use LIS_timeMgrMod,   only : LIS_calendar
  use ARMdata_module,   only : ARMdata_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

!DESCRIPTION: 
!  This routine reads the bulk aerodynamic estimates of 
!  latent, sensible and ground heat fluxes from the BAEBBR file, 
!  for a particular station. The data is read from the 
!  native NETCDF files, and is processed by applying the 
!  quality control flags. Only data points with 'good'
!  classification (qc value=0) is chosen. This routine 
!  also computes the temporal offset of the data relative
!  to the 0z of a particular day. 
!
!  The arguments are: 
!  \begin{description}
!   \item[k]         Index of the ARM station
!   \item[yr]        year of ARM data 
!   \item[mo]        month of ARM data
!   \item[da]        day of ARM data
!   \item[filename]  Name of the BAEBBR file
!   \item[bend]      bookend index (2-current, 1-previous)
!  \end{description}
!EOP


  integer,  intent(in) :: k 
  integer,  intent(in) :: n
  integer              :: yr
  integer              :: mo
  integer              :: da
  character(len=*)     :: filename 
  integer              :: bend

  integer       :: nid, btimeid, latid, lonid, qleid
  integer       :: qhid,qgid, timeid, dimId
  integer       :: qc_qhid, qc_qleid, qc_qgid
  integer       :: ndims
  real          :: lat, lon 
  integer       :: base_time
  integer       :: ios
  real, allocatable :: time(:)
  real, allocatable :: baebbr_qh(:)
  real, allocatable :: baebbr_qle(:)
  real, allocatable :: baebbr_qg(:)
  integer, allocatable :: baebbr_qc_qh(:)
  integer, allocatable :: baebbr_qc_qle(:)
  integer, allocatable :: baebbr_qc_qg(:)
  integer           :: kk
  type(ESMF_Time)   :: reftime, datatime, currtime
  type(ESMF_TimeInterval) :: dt,ts
  integer           :: yr1, mo1, da1, hr1, mn1, ss1
  integer           :: status
  integer           :: data_index

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
  call LIS_verify(ios, 'Error opening file'//trim(filename))

!variable ids
  ios = nf90_inq_varid(nid, 'base_time',btimeid)
  call LIS_verify(ios, 'Error nf90_inq_varid: base_time')

  ios = nf90_inq_varid(nid, 'lat',latid)
  call LIS_verify(ios, 'Error nf90_inq_varid: lat')

  ios = nf90_inq_varid(nid, 'lon',lonid)
  call LIS_verify(ios, 'Error nf90_inq_varid: lon')

  ios = nf90_inq_varid(nid, 'time_offset',timeid)
  call LIS_verify(ios, 'Error nf90_inq_varid: time')

  ios = nf90_inq_varid(nid, 'be_sensible_heat_flux',qhid)
  call LIS_verify(ios, 'Error nf90_inq_varid: qh')

  ios = nf90_inq_varid(nid, 'qc_be_sensible_heat_flux',qc_qhid)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_qh')

  ios = nf90_inq_varid(nid, 'be_latent_heat_flux',qleid)
  call LIS_verify(ios, 'Error nf90_inq_varid: qle')

  ios = nf90_inq_varid(nid, 'qc_be_latent_heat_flux',qc_qleid)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_qle')

  ios = nf90_inq_varid(nid, 'surface_soil_heat_flux_avg',qgid)
  call LIS_verify(ios, 'Error nf90_inq_varid: qg')

  ios = nf90_inq_varid(nid, 'qc_surface_soil_heat_flux_avg',qc_qgid)
  call LIS_verify(ios, 'Error nf90_inq_varid: qc_qg')

!dimensions
  ios = nf90_inq_dimid(nid, 'time',dimId)
  call LIS_verify(ios, 'Error nf90_inq_dimid: time')

  ios = nf90_inquire_dimension(nid, dimId, len=ndims)
  call LIS_verify(ios, 'Error nf90_inquire_dimension:')

!values
  ios = nf90_get_var(nid,latid, lat)
  call LIS_verify(ios, 'Error nf90_get_var: lat')

  ios = nf90_get_var(nid,lonid, lon)
  call LIS_verify(ios, 'Error nf90_get_var: lon')

  ios = nf90_get_var(nid,btimeid, base_time)
  call LIS_verify(ios, 'Error nf90_get_var: base_time')

  call ESMF_TimeSet(refTime, yy=1970,mm=1,dd=1,&
       h=0,m=0,s=0,calendar=LIS_calendar,rc=status)
  call LIS_verify(status, 'error in timeset: readARMObs')

!  ARMdata_struc(n)%stnlat(k) = lat
!  ARMdata_struc(n)%stnlon(k) = lon

  allocate(time(ndims))
  allocate(baebbr_qle(ndims))
  allocate(baebbr_qh(ndims))
  allocate(baebbr_qg(ndims))

  allocate(baebbr_qc_qle(ndims))
  allocate(baebbr_qc_qh(ndims))
  allocate(baebbr_qc_qg(ndims))

  ios = nf90_get_var(nid,timeid, time)
  call LIS_verify(ios, 'Error nf90_get_var: time')

  ios = nf90_get_var(nid,qhid, baebbr_qh)
  call LIS_verify(ios, 'Error nf90_get_var: qh')

  ios = nf90_get_var(nid,qc_qhid, baebbr_qc_qh)
  call LIS_verify(ios, 'Error nf90_get_var: qc_qh')

  ios = nf90_get_var(nid,qleid, baebbr_qle)
  call LIS_verify(ios, 'Error nf90_get_var: qle')

  ios = nf90_get_var(nid,qc_qleid, baebbr_qc_qle)
  call LIS_verify(ios, 'Error nf90_get_var: qc_qle')

  ios = nf90_get_var(nid,qgid, baebbr_qg)
  call LIS_verify(ios, 'Error nf90_get_var: qg')

  ios = nf90_get_var(nid,qc_qgid, baebbr_qc_qg)
  call LIS_verify(ios, 'Error nf90_get_var: qc_qg')

  ios = nf90_close(nid)
  call LIS_verify(ios, 'Error in nf90_close')
  
  call ESMF_TimeIntervalSet(dt,s=base_time,rc=status)
  call LIS_verify(status, 'Error in timeintervalset: readARMobs')

  reftime = reftime + dt

  do kk=1,ndims
     call ESMF_TimeIntervalSet(dt,s=nint(time(kk)),rc=status)
     call LIS_verify(status, 'Error in timeintervalset: readARMobs')

     datatime = reftime + dt
     
     call ESMF_TimeGet(datatime, yy=yr1,mm=mo1,dd=da1,&
          h=hr1,m=mn1,s=ss1,calendar=LIS_calendar,rc=status)
     call LIS_verify(status, 'error in timeget: readARMObs')

! data index from 0z with 30mn timestep. All this extensive processing is
! required because some of the ARM files report data at times different 
! from 0z. 
     
     call ESMF_TimeSet(currTime, yy=yr, mm=mo, &
          dd=da, h=0,m=0,s=0,calendar=LIS_calendar,rc=status)
     call LIS_verify(status, 'error in timeget: readARMObs')

     data_index = (datatime - currtime)/ARMdata_struc(n)%baebbr_ts + 1
     ARMdata_struc(n)%baebbr_tindex(k,data_index) = data_index
     if(bend.eq.2) then 
        if(baebbr_qc_qh(kk).eq.0) then        
           ARMdata_struc(n)%baebbr_qh_c(k,data_index)     = baebbr_qh(kk)
        else
           ARMdata_struc(n)%baebbr_qh_c(k,data_index)     = LIS_rc%udef
        endif
        if(baebbr_qc_qle(kk).eq.0) then 
           ARMdata_struc(n)%baebbr_qle_c(k,data_index)    = baebbr_qle(kk)
        else
           ARMdata_struc(n)%baebbr_qle_c(k,data_index)    = LIS_rc%udef
        endif
        
        if(baebbr_qc_qg(kk).eq.0) then
           ARMdata_struc(n)%baebbr_qg_c(k,data_index)     = baebbr_qg(kk)
        else
           ARMdata_struc(n)%baebbr_qg_c(k,data_index)     = LIS_rc%udef
        endif
     elseif(bend.eq.1) then 
        if(baebbr_qc_qh(kk).eq.0) then        
           ARMdata_struc(n)%baebbr_qh_p(k,data_index)     = baebbr_qh(kk)
        else
           ARMdata_struc(n)%baebbr_qh_p(k,data_index)     = LIS_rc%udef
        endif
        if(baebbr_qc_qle(kk).eq.0) then 
           ARMdata_struc(n)%baebbr_qle_p(k,data_index)    = baebbr_qle(kk)
        else
           ARMdata_struc(n)%baebbr_qle_p(k,data_index)    = LIS_rc%udef
        endif
        
        if(baebbr_qc_qg(kk).eq.0) then 
           ARMdata_struc(n)%baebbr_qg_p(k,data_index)     = baebbr_qg(kk)
        else
           ARMdata_struc(n)%baebbr_qg_p(k,data_index)     = LIS_rc%udef
        endif
     endif
  enddo

  deallocate(time)
  deallocate(baebbr_qle)
  deallocate(baebbr_qh)
  deallocate(baebbr_qg)

  deallocate(baebbr_qc_qle)
  deallocate(baebbr_qc_qh)
  deallocate(baebbr_qc_qg)

#endif

end subroutine read_baebbr_flux_file


!BOP
! 
! !ROUTINE: create_arm_ecor_flux_filename
! \label{create_arm_ecor_flux_filename}
! 
! !INTERFACE: 
subroutine create_arm_ecor_flux_filename(odir, site_id, stnid, &
     yr, mo, da, filename, rc)
  
  use LIS_coreMod,  only : LIS_localPet
  implicit none

! !ARGUMRENTS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: site_id
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
  character(len=*), intent(out) :: filename
  integer,          intent(out) :: rc 

! !DESCRIPTION: 
! 
! This routine creates a filename for ARM in-situ data files 
! from the eddy correlation (ECOR) flux measurement system. The 
! system produces 30mn estimates of vertical fluxes of 
! sensible and latent heat fluxes. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      ARM base directory
!   \item[site\_id]  ARM site identifier
!   \item[stnid]     Station ID 
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the ECOR file
!   \item[rc]        return code (0-success, 1-failed to generate)
!  \end{description}
!EOP

  character*4       :: fyr
  character*2       :: fmo
  character*2       :: fda
  character*4       :: fproc
  
  integer           :: fsize
  character*100     :: ls_comm, cmd2

  rc = 1 !fail to find the file

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fproc,fmt='(i4.4)') LIS_localPet

  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//&
       trim(site_id)//'30ecor'//&
       trim(adjustl(stnid))//'.'//'b1.'//&
       trim(fyr)//trim(fmo)//trim(fda)//&
       '*cdf 2>&1 2>/dev/null > ecor_file.'//trim(fproc)
  
  cmd2 = 'wc -w ecor_file.'//trim(fproc)//' > ecor_file_wc.'//trim(fproc)

  call system(ls_comm)
  call system(cmd2)

  open(110,file='ecor_file_wc.'//trim(fproc),form='formatted',action='read') 
  read(110,*) fsize
  close(110)

  if(fsize.gt.0) then 
     open(110,file='ecor_file.'//trim(fproc),form='formatted',action='read')
     read(110,'(a)') filename
     close(110)
     rc=0
  endif
  
end subroutine create_arm_ecor_flux_filename


!BOP
! 
! !ROUTINE: create_arm_baebbr_flux_filename
! \label{create_arm_baebbr_flux_filename}
! 
! !INTERFACE: 
subroutine create_arm_baebbr_flux_filename(odir, site_id, stnid, &
     yr, mo, da, filename, rc)
  
  use LIS_coreMod, only : LIS_localPet

  implicit none

! !ARGUMRENTS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: site_id
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
  character(len=*), intent(out) :: filename
  integer,          intent(out) :: rc 

! !DESCRIPTION: 
! 
! This routine creates a filename for ARM in-situ data files 
! from the Energy Balance Bowen Ratio (BAEBBR) system. The 
! system produces 30mn estimates of vertical fluxes of 
! sensible and latent heat fluxes. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      ARM base directory
!   \item[site\_id]  ARM site identifier
!   \item[stnid]     Station ID 
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the BAEBBR file
!   \item[rc]        return code (0-success, 1-failed to generate)
!  \end{description}
!EOP

  character*4       :: fyr
  character*2       :: fmo
  character*2       :: fda
  character*4       :: fproc
  
  integer           :: fsize
  character*100     :: ls_comm, cmd2

  rc = 1 !fail to find the file

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fproc,fmt='(i4.4)') LIS_localPet

  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//&
       trim(site_id)//'30baebbr'//&
       trim(adjustl(stnid))//'.'//'s1.'//&
       trim(fyr)//trim(fmo)//trim(fda)//&
       '*cdf 2>&1 2>/dev/null > baebbr_file.'//trim(fproc)
  
  cmd2 = 'wc -w baebbr_file.'//trim(fproc)//' > baebbr_file_wc.'//trim(fproc)

  call system(ls_comm)
  call system(cmd2)

  open(110,file='baebbr_file_wc.'//trim(fproc),form='formatted',action='read') 
  read(110,*) fsize
  close(110)

  if(fsize.gt.0) then 
     open(110,file='baebbr_file.'//trim(fproc),form='formatted',action='read')
     read(110,'(a)') filename
     close(110)
     rc=0
  endif
  
end subroutine create_arm_baebbr_flux_filename


!BOP
! 
! !ROUTINE: create_arm_ebbr_filename
! \label{create_arm_ebbr_filename}
! 
! !INTERFACE: 
subroutine create_arm_ebbr_filename(odir, site_id, stnid, &
     yr, mo, da, filename, rc)
! !USES: 
  use LIS_coreMod, only : LIS_localPet
  
  implicit none

! !ARGUMRENTS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: site_id
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
  character(len=*), intent(out) :: filename
  integer,          intent(out) :: rc 

! !DESCRIPTION: 
! 
! This routine creates a filename for ARM in-situ data files 
! from the Energy Balance Bowen Ratio (EBBR) system.
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      ARM base directory
!   \item[site\_id]  ARM site identifier
!   \item[stnid]     Station ID 
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the EBBR file
!   \item[rc]        return code (0-success, 1-failed to generate)
!  \end{description}
!EOP

  character*4       :: fyr
  character*2       :: fmo
  character*2       :: fda
  character*4       :: fproc

  integer           :: fsize
  character*100     :: ls_comm, cmd2

  rc = 1 !fail to find the file

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fproc,fmt='(i4.4)') LIS_localPet

  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//&
       trim(site_id)//'30ebbr'//&
       trim(adjustl(stnid))//'.'//'b1.'//&
       trim(fyr)//trim(fmo)//trim(fda)//&
       '*cdf 2>&1 2>/dev/null > ebbr_file.'//trim(fproc)

  cmd2 = 'wc -w ebbr_file.'//trim(fproc)//' > ebbr_file_wc.'//trim(fproc)  
  
  call system(ls_comm)
  call system(cmd2)

  open(110,file='ebbr_file_wc.'//trim(fproc),form='formatted',action='read') 
  read(110,*) fsize
  close(110)

  if(fsize.gt.0) then 
     open(110,file='ebbr_file.'//trim(fproc),form='formatted',action='read')
     read(110,'(a)') filename
     close(110)
     rc=0
  endif

end subroutine create_arm_ebbr_filename


