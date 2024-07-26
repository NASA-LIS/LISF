!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readSMOSNESDISsmObs
! \label{readSMOSNESDISsmObs}
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readSMOSNESDISsmObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod,      only : LDT_rc
  use LDT_timeMgrMod,   only : LDT_get_julss
  use LDT_logMod,       only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use SMOSNESDISsm_obsMod, only : SMOSNESDISsmobs
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the ESACCI
! soil moisture retrieval product. 
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character(len=LDT_CONST_PATH_LEN)     :: fname
  real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  SMOSNESDISsmobs(n)%smobs = LDT_rc%udef
  smobs= LDT_rc%udef

  call create_SMOSNESDISsm_filename(SMOSNESDISsmobs(n)%odir, &
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)
  
  inquire(file=trim(fname),exist=file_exists)
  if(file_exists) then
     
     write(LDT_logunit,*) '[INFO] Reading ..',trim(fname)
     call read_SMOSNESDIS_data(n, fname, smobs)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(smobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
              SMOSNESDISsmobs(n)%smobs(c,r) = smobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%soilmoist_obs,&
       SMOSNESDISsmobs(n)%smobs,vlevel=1)

end subroutine readSMOSNESDISsmObs


!BOP
! 
! !ROUTINE: read_SMOSNESDIS_data
! \label(read_SMOSNESDIS_data)
!
! !INTERFACE:
subroutine read_SMOSNESDIS_data(n, fname, smobs_ip)
! 
! !USES:   
  use LDT_coreMod
  use LDT_logMod
  use map_utils
  use SMOSNESDISsm_obsMod, only : SMOSNESDISsmobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname
  real                          :: smobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SMOS NESDIS file
!  \item[smobs\_ip]    soil moisture data processed to the LDT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  integer*2         :: sm(SMOSNESDISsmobs(n)%nc,SMOSNESDISsmobs(n)%nr)
  real              :: sm_combined(SMOSNESDISsmobs(n)%nc,SMOSNESDISsmobs(n)%nr)
  real              :: sm_data(SMOSNESDISsmobs(n)%nc*SMOSNESDISsmobs(n)%nr)
  logical*1         :: sm_data_b(SMOSNESDISsmobs(n)%nc*SMOSNESDISsmobs(n)%nr)
  logical*1         :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  integer           :: c,r,i,j
  integer           :: ftn
  integer           :: smid
  integer           :: ios

  ftn = LDT_getNextUnitNumber()
  open(ftn,file=trim(fname), form='unformatted',access='direct',&
       recl=SMOSNESDISsmobs(n)%nc*SMOSNESDISsmobs(n)%nr*2,&
       convert='little_endian')
  read(ftn,rec=1) sm
  call LDT_releaseUnitNumber(ftn)

  do r=1, SMOSNESDISsmobs(n)%nr
     do c=1, SMOSNESDISsmobs(n)%nc
        if(sm(c,r).lt.0.001) then 
           sm_combined(c,SMOSNESDISsmobs(n)%nr-r+1) = LDT_rc%udef
        else
           sm_combined(c,SMOSNESDISsmobs(n)%nr-r+1) = sm(c,r)*0.0001
        endif
     enddo
  enddo
 
  do r=1, SMOSNESDISsmobs(n)%nr
     do c=1, SMOSNESDISsmobs(n)%nc
        sm_data(c+(r-1)*SMOSNESDISsmobs(n)%nc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LDT_rc%udef) then 
           sm_data_b(c+(r-1)*SMOSNESDISsmobs(n)%nc) = .true. 
        else
           sm_data_b(c+(r-1)*SMOSNESDISsmobs(n)%nc) = .false.
        endif
        if(sm_combined(c,r).gt.0.5) then 
           sm_combined(c,r) = LDT_rc%udef
           sm_data_b(c+(r-1)*SMOSNESDISsmobs(n)%nc) = .false.
        endif
     enddo
  enddo
  
!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LDT_rc%gridDesc(n,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       SMOSNESDISsmobs(n)%nc*SMOSNESDISsmobs(n)%nr, &
       LDT_rc%lnc(n)*LDT_rc%lnr(n), &
       LDT_domain(n)%lat, LDT_domain(n)%lon,&
       SMOSNESDISsmobs(n)%n11,&
       LDT_rc%udef, ios)

end subroutine read_SMOSNESDIS_data

!BOP
! !ROUTINE: create_SMOSNESDISsm_filename
! \label{create_SMOSNESDISsm_filename}
! 
! !INTERFACE: 
subroutine create_SMOSNESDISsm_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the SMOS NESDIS filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the SMOS NESDIS soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated SMOS NESDIS filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  logical           :: oper_flag
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

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

  
  if(oper_flag) then 
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/SMOS_OPER_MIR_SMUDP2_SoilMoisture_'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.dat'
  else
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/SMOS_REPR_MIR_SMUDP2_SoilMoisture_'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.dat'
  endif

end subroutine create_SMOSNESDISsm_filename
