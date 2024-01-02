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
! !ROUTINE: readSMOSNESDISsmobs
! \label{readSMOSNESDISsmobs}
!
! !INTERFACE: 
subroutine readSMOSNESDISsmobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit
  use SMOSNESDIS_smobsMod, only : SMOSNESDIS_smobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in)       :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the standard 
! NASA soil moisture retrieval product. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
!EOP

  logical           :: alarmcheck, file_exists, readflag
  integer           :: iret
  character*200     :: name
  real              :: smc(LVT_rc%lnc, LVT_rc%lnr)
  integer           :: fnd 
  real              :: timenow

  smc = LVT_rc%udef

  timenow = float(LVT_rc%dhr(source))*3600 +&
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  if(SMOSNESDIS_smobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     SMOSNESDIS_smobs(source)%startflag = .false. 
     call SMOSNESDIS_sm_filename(source,name,&
          SMOSNESDIS_smobs(source)%odir, & 
          LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source))
                 
     inquire(file=name, exist=file_exists) 

     if(file_exists) then 
        readflag = .true. 
     else
        readflag = .false. 
     endif
     
     if(readflag) then 
        write(LVT_logunit,*) '[INFO] Reading SMOSNESDIS file ',name
        call read_SMOSNESDISsm(source, name, smc)
        
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source,&
       smc,vlevel=1,units="m3/m3")
 
end subroutine readSMOSNESDISsmobs

!BOP
! 
! !ROUTINE: read_SMOSNESDISsm
! \label{read_SMOSNESDISsm}
!
! !INTERFACE: 
subroutine read_SMOSNESDISsm(source, name, smobs)
! 
! !USES: 
  use LVT_coreMod,         only : LVT_rc,LVT_domain
  use LVT_logMod
  use SMOSNESDIS_smobsMod, only : SMOSNESDIS_smobs

  implicit none

  integer                        :: source
  character(len=*)               :: name
  real                           :: smobs(LVT_rc%lnc,LVT_rc%lnr)
  
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

  integer*2           :: sm_raw(SMOSNESDIS_smobs(source)%nc,&
       SMOSNESDIS_smobs(source)%nr)
  integer*2           :: sm_raw_time(SMOSNESDIS_smobs(source)%nc,&
       SMOSNESDIS_smobs(source)%nr)
  real                :: sm_data(SMOSNESDIS_smobs(source)%nc*&
       SMOSNESDIS_smobs(source)%nr)
  integer*1           :: sm_data_hr(SMOSNESDIS_smobs(source)%nc,&
       SMOSNESDIS_smobs(source)%nr)
  integer*1           :: sm_data_mn(SMOSNESDIS_smobs(source)%nc,&
       SMOSNESDIS_smobs(source)%nr)
  integer*1           :: sm_data_hr1(SMOSNESDIS_smobs(source)%nc*&
       SMOSNESDIS_smobs(source)%nr)
  integer*1           :: sm_data_mn1(SMOSNESDIS_smobs(source)%nc*&
       SMOSNESDIS_smobs(source)%nr)
  logical*1           :: sm_data_b(SMOSNESDIS_smobs(source)%nc*&
       SMOSNESDIS_smobs(source)%nr)
  logical*1           :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)
  real                :: sm_ip(LVT_rc%lnc*LVT_rc%lnr)
  integer             :: hr_val, mn_val
  integer             :: julss
  integer             :: c,r,c1,r1,i,j,kk,ios
  integer             :: ftn,iret,igrib,nvars
  integer             :: param_num
  logical             :: var_found
  real                :: err, ql, udef
 
  ftn = LVT_getNextUnitNumber()
  open(ftn,file=trim(name),form='unformatted',access='direct',&
       recl=SMOSNESDIS_smobs(source)%nc*SMOSNESDIS_smobs(source)%nr*2,&
       convert='little_endian')
  read(ftn,rec=1) sm_raw
  read(ftn,rec=2) sm_raw_time
  call LVT_releaseUnitNumber(ftn)
  
  do r=1,SMOSNESDIS_smobs(source)%nr
     do c=1,SMOSNESDIS_smobs(source)%nc
        if((sm_raw(c,r)*0.0001.lt.0.01).or.sm_raw_time(c,r).lt.0) then 
           
           sm_data_hr(c,r) = -1
           sm_data_mn(c,r) = -1
        else
           sm_data_hr(c,r) = sm_raw_time(c,r)/100
           sm_data_mn(c,r) = sm_raw_time(c,r) - &
                sm_data_hr(c,r)*100
        endif
     enddo
  enddo

  do r=1,SMOSNESDIS_smobs(source)%nr
     do c=1,SMOSNESDIS_smobs(source)%nc
        if((sm_raw(c,r)*0.0001.lt.0.01).or.sm_raw_time(c,r).lt.0) then 
           sm_data(c+((SMOSNESDIS_smobs(source)%nr-r+1)-1)*&
                SMOSNESDIS_smobs(source)%nc) = LVT_rc%udef
           sm_data_hr1(c+((SMOSNESDIS_smobs(source)%nr-r+1)-1)*&
                SMOSNESDIS_smobs(source)%nc) = sm_data_hr(c,r)
           sm_data_mn1(c+((SMOSNESDIS_smobs(source)%nr-r+1)-1)*&
                SMOSNESDIS_smobs(source)%nc) = sm_data_mn(c,r)
           sm_data_b(c+((SMOSNESDIS_smobs(source)%nr-r+1)-1)*&
                SMOSNESDIS_smobs(source)%nc) = .false.
        else
           sm_data(c+((SMOSNESDIS_smobs(source)%nr-r+1)-1)*&
                SMOSNESDIS_smobs(source)%nc) = &
                sm_raw(c,r)*0.0001
           sm_data_hr1(c+((SMOSNESDIS_smobs(source)%nr-r+1)-1)*&
                SMOSNESDIS_smobs(source)%nc) =  sm_data_hr(c,r)
           sm_data_mn1(c+((SMOSNESDIS_smobs(source)%nr-r+1)-1)* &
                SMOSNESDIS_smobs(source)%nc) =  sm_data_mn(c,r)
           sm_data_b(c+((SMOSNESDIS_smobs(source)%nr-r+1)-1)*&
                SMOSNESDIS_smobs(source)%nc) = .true. 

        endif
     enddo
  enddo

  udef = -9999.0

  call neighbor_interp(LVT_rc%gridDesc, &
       sm_data_b, sm_data, smobs_b_ip, sm_ip, &
       SMOSNESDIS_smobs(source)%nc*SMOSNESDIS_smobs(source)%nr,&
       LVT_rc%lnc*LVT_rc%lnr,&
       SMOSNESDIS_smobs(source)%rlat2, SMOSNESDIS_smobs(source)%rlon2,&
       SMOSNESDIS_smobs(source)%n112,udef, iret)

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        c1 = c
        r1 = r
        smobs(c1,r1) = sm_ip(c+(r-1)*LVT_rc%lnc)
     enddo
  enddo

end subroutine read_SMOSNESDISsm

!BOP
! 
! !ROUTINE: SMOSNESDIS_sm_filename
! \label{SMOSNESDIS_sm_filename}
!
! !INTERFACE: 
subroutine SMOSNESDIS_sm_filename(source, filename, ndir, yr, mo,da)
! 
! !USES:   
  use LVT_coreMod,only : LVT_rc
  use LVT_logMod, only : LVT_logunit

  implicit none
!
! !ARGUMENTS: 
  integer            :: source
  character(len=*)   :: filename
  integer            :: yr, mo, da, hr,mn
  character (len=*)  :: ndir
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine creates the NASA AMSRE filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the NASA AMSRE soil moisture filename
!  \item[ndir] name of the NASA AMSRE soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  logical           :: oper_flag
  
 oper_flag = .false.

  if(yr.gt.2011) then
     if(yr.eq.2011) then 
        if(mo.eq.12) then 
           oper_flag = .true. 
        elseif(mo.eq.11) then 
           if(da.eq.30) then 
              oper_flag = .true. 
           endif
        endif
     else
        oper_flag = .true. 
     endif
  endif

  
  write(unit=fyr, fmt='(i4.4)') LVT_rc%dyr(source)
  write(unit=fmo, fmt='(i2.2)') LVT_rc%dmo(source)
  write(unit=fda, fmt='(i2.2)') LVT_rc%dda(source)
  
  if(oper_flag) then 
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/SMOS_OPER_MIR_SMUDP2_SoilMoisture_'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.dat'
  else
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/SMOS_REPR_MIR_SMUDP2_SoilMoisture_'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.dat'
  endif

end subroutine SMOSNESDIS_sm_filename
