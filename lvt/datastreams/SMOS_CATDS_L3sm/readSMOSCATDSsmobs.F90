!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readSMOSCATDSsmobs
! \label{readSMOSCATDSsmobs}
!
! !INTERFACE: 
subroutine readSMOSCATDSsmobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit
  use SMOSCATDS_smobsMod, only : SMOSCATDS_smobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in)       :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the 
! Level 3 SMOS CATDS product. SMOS overpass time 
! ascending - 6AM and descending - 6pm 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  28 Jan 2017: Sujay Kumar, Initial Specification
! 
!EOP

  logical           :: alarmcheck, file_exists, readflag
  integer           :: iret
  character*200     :: name
  real              :: smc(LVT_rc%lnc, LVT_rc%lnr)
  real              :: smc_A(LVT_rc%lnc, LVT_rc%lnr)
  real              :: smc_D(LVT_rc%lnc, LVT_rc%lnr)
  integer           :: fnd,c,r
  real              :: timenow

  smc   = LVT_rc%udef
  smc_A = LVT_rc%udef
  smc_D = LVT_rc%udef

  timenow = float(LVT_rc%dhr(source))*3600 +&
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  if(SMOSCATDS_smobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     SMOSCATDS_smobs(source)%startflag = .false. 

     call SMOSCATDS_sm_filename(source,name,&
          SMOSCATDS_smobs(source)%odir, 'A',& 
        LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source))
                 
     inquire(file=name, exist=file_exists) 

     if(file_exists) then 
        readflag = .true. 
     else
        readflag = .false. 
     endif
     
     if(readflag) then 
        write(LVT_logunit,*) '[INFO] Reading SMOSCATDS file ',name
        call read_SMOSCATDSsm(source, name, smc_A)
        
     endif

     call SMOSCATDS_sm_filename(source,name,&
          SMOSCATDS_smobs(source)%odir, 'D',& 
        LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source))
                 
     inquire(file=name, exist=file_exists) 

     if(file_exists) then 
        readflag = .true. 
     else
        readflag = .false. 
     endif
     
     if(readflag) then 
        write(LVT_logunit,*) '[INFO] Reading SMOSCATDS file ',name
        call read_SMOSCATDSsm(source, name, smc_D)
        
     endif
     
     smc = smc_A
     
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(smc_D(c,r).ne.-9999.0) then 
              smc(c,r) = smc_D(c,r)
           endif
        enddo
     enddo

  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source,&
       smc,vlevel=1,units="m3/m3")
 
end subroutine readSMOSCATDSsmobs

!BOP
! 
! !ROUTINE: read_SMOSCATDSsm
! \label{read_SMOSCATDSsm}
!
! !INTERFACE: 
subroutine read_SMOSCATDSsm(source, filename, smobs)
! 
! !USES: 
  use LVT_coreMod,         only : LVT_rc,LVT_domain
  use LVT_logMod
  use SMOSCATDS_smobsMod, only : SMOSCATDS_smobs

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  integer                        :: source
  real                           :: smobs(LVT_rc%lnc,LVT_rc%lnr)
  character(len=*)               :: filename
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  integer                        :: iret,smid,rfiid,ios,nid,c,r,c1,r1,t
  real                           :: udef,sm_offset, sm_scalef, rfi_scalef, rfi_offset
  real                           :: rfi_prob_val

  logical*1                  :: li(SMOSCATDS_smobs(source)%nc*SMOSCATDS_smobs(source)%nr)
  integer*2                  :: sm(SMOSCATDS_smobs(source)%nc, SMOSCATDS_smobs(source)%nr)
  integer*1                  :: rfi_prob(SMOSCATDS_smobs(source)%nc,SMOSCATDS_smobs(source)%nr)
  real                       :: sm1d(SMOSCATDS_smobs(source)%nc*SMOSCATDS_smobs(source)%nr)

  logical*1                  :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                       :: sm_ip(LVT_rc%lnc*LVT_rc%lnr)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
  call LVT_verify(ios, 'Error opening file'//trim(filename))

! Soil moisture
  ios = nf90_inq_varid(nid, 'Soil_Moisture',smid)
  call LVT_verify(ios, 'Error nf90_inq_varid: Soil_Moisture')

  ios = nf90_get_var(nid,smid, sm)
  call LVT_verify(ios, 'Error nf90_get_var: Soil_Moisture')

  ios = nf90_get_att(nid,smid,'add_offset',sm_offset)
  call LVT_verify(ios, 'Error nf90_get_att: add_offset Soil_Moisture')

  ios = nf90_get_att(nid,smid,'scale_factor',sm_scalef)
  call LVT_verify(ios, 'Error nf90_get_att: scale_factor Soil_Moisture')

!RFI probability
  ios = nf90_inq_varid(nid, 'Rfi_Prob',rfiid)
  call LVT_verify(ios, 'Error nf90_inq_varid: Rfi_prob')

  ios = nf90_get_var(nid,rfiid, rfi_prob)
  call LVT_verify(ios, 'Error nf90_get_var: rfi_prob')

  ios = nf90_get_att(nid,rfiid,'add_offset',rfi_offset)
  call LVT_verify(ios, 'Error nf90_get_att: add_offset rfi_prob')

  ios = nf90_get_att(nid,rfiid,'scale_factor',rfi_scalef)
  call LVT_verify(ios, 'Error nf90_get_att: scale_factor rfi_prob')

  ios = nf90_close(nid)
  call LVT_verify(ios, 'Error nf90_close')

  li = .false. 

  do r=1,SMOSCATDS_smobs(source)%nr
     do c=1,SMOSCATDS_smobs(source)%nc       
        t = (c+(SMOSCATDS_smobs(source)%nr-r+1-1)*&
             SMOSCATDS_smobs(source)%nc)        
        rfi_prob_val = rfi_prob(c,r)*rfi_scalef+rfi_offset
        if(sm(c,r)*sm_scalef+sm_offset.gt.0) then
           if(rfi_prob(c,r).ne.-128) then 
              if(rfi_prob_val.lt.0.05) then
                 sm1d(t) = &
                      sm(c,r)*sm_scalef+sm_offset
                 li(t) = .true.
              endif
           endif
        else
           sm1d(t) = -9999.0
        endif
     enddo
  enddo

  udef = -9999.0

  call neighbor_interp(LVT_rc%gridDesc, li, sm1d, lo, sm_ip, &
       SMOSCATDS_smobs(source)%nc*SMOSCATDS_smobs(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMOSCATDS_smobs(source)%rlat2, SMOSCATDS_smobs(source)%rlon2,&
       SMOSCATDS_smobs(source)%n112,udef, iret)

  smobs = -9999.0

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        smobs(c,r) = sm_ip(c+(r-1)*LVT_rc%lnc)
     enddo
  enddo

#endif

end subroutine read_SMOSCATDSsm

!BOP
! 
! !ROUTINE: SMOSCATDS_sm_filename
! \label{SMOSCATDS_sm_filename}
!
! !INTERFACE: 

subroutine SMOSCATDS_sm_filename(source, name, ndir, pass, yr, mo,da)
! 
! !USES:   
  use LVT_coreMod,only : LVT_rc
  use LVT_logMod, only : LVT_logunit

  implicit none
!
! !ARGUMENTS: 
  integer            :: source
  character(len=*)   :: name
  character(len=*)   :: pass
  integer            :: yr, mo, da, hr,mn
  character (len=*)  :: ndir
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine creates the SMOS CATDS L3 filename based on the time and date 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  logical           :: oper_flag

  write(unit=fyr, fmt='(i4.4)') LVT_rc%dyr(source)
  write(unit=fmo, fmt='(i2.2)') LVT_rc%dmo(source)
  write(unit=fda, fmt='(i2.2)') LVT_rc%dda(source)
  
  oper_flag = .true. 
  if(yr.gt.2015) then 
     oper_flag = .true. 

  elseif(yr.eq.2015) then 
     if(mo.gt.5) then 
        oper_flag = .true. 
     elseif(mo.eq.5) then 
        if(da.ge.6) then 
           oper_flag = .true. 
        else
           oper_flag = .false. 
        endif
     else
        oper_flag = .false.
     endif
        
  else
     oper_flag = .false.
  endif

  if(oper_flag) then 
     name = trim(ndir)//'/'//trim(fyr)//'/SM_OPER_MIR_CLF31' &
          //trim(pass)//'_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          'T000000_'//trim(fyr)//trim(fmo)//trim(fda)//&
          'T235959_300_001_7.DBL.nc'

  else
     name = trim(ndir)//'/'//trim(fyr)//'/SM_RE04_MIR_CLF31' &
          //trim(pass)//'_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          'T000000_'//trim(fyr)//trim(fmo)//trim(fda)//&
          'T235959_300_001_7.DBL.nc'
  endif
end subroutine SMOSCATDS_sm_filename
