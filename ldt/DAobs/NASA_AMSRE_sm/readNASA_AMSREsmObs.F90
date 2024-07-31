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
! !ROUTINE: readNASA_AMSREsmObs
! \label{readNASA_AMSREsmObs}
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readNASA_AMSREsmObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod,      only : LDT_rc
  use LDT_logMod,       only : LDT_logunit
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use NASA_AMSREsm_obsMod, only : NASA_AMSREsmobs

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the standard 
! NASA soil moisture retrieval product. 
!
!EOP
  
  logical           :: alarmcheck, file_exists, readflag
  integer           :: iret
  character(len=LDT_CONST_PATH_LEN)     :: name
  real              :: smc(LDT_rc%lnc(n), LDT_rc%lnr(n))
  integer           :: fnd 
  integer           :: c,r
  real              :: timenow

!  alarmCheck = ESMF_AlarmIsRinging(NASA_AMSREsmobs%readAlarm, rc=iret)
!  if(alarmCheck) then 
!     call ESMF_AlarmRingerOff(NASA_AMSREsmobs%readAlarm, rc=iret)

  timenow = float(LDT_rc%hr)*3600 + 60*LDT_rc%mn + LDT_rc%ss
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  if(NASA_AMSREsmobs%startflag.or.alarmCheck) then 
     
     NASA_AMSREsmobs%startflag = .false. 
     call NASA_AMSREsm_filename(name,NASA_AMSREsmobs%odir, & 
        LDT_rc%yr, LDT_rc%mo, LDT_rc%da)
                 
     inquire(file=name, exist=file_exists) 

     if(file_exists) then 
        readflag = .true. 
     else
        readflag = .false. 
     endif
     
     if(readflag) then 
        write(LDT_logunit,*) 'Reading NASA AMSRE file ',name
        call read_AMSREsm(n,name)
        call maskObs_basedonQC(LDT_rc%lnc(n)*LDT_rc%lnr(n), &
             NASA_AMSREsmobs%smobs, &
             NASA_AMSREsmobs%smqc)
     endif
  endif

  fnd = 0 
  smc = LDT_rc%udef
  call maskObs_basedonTime(n, LDT_rc%lnc(n)*LDT_rc%lnr(n), &
       NASA_AMSREsmobs%smobs, &
       NASA_AMSREsmobs%smtime, smc,fnd)

!  open(100,file='smc.bin',form='unformatted')
!  write(100) smc
!  close(100) 
!  stop
   call LDT_logSingleDAobs(n,LDT_DAobsData(1)%soilmoist_obs, smc,vlevel=1)

end subroutine readNASA_AMSREsmObs


!BOP
! !ROUTINE: read_AMSREsm
! \label{read_AMSREsm}
! 
! !INTERFACE: 
subroutine read_AMSREsm(n, name)
! !USES: 
  use LDT_coreMod,         only : LDT_rc,LDT_domain
  use LDT_logmod,          only : LDT_logunit
  use NASA_AMSREsm_obsMod, only : NASA_AMSREsmobs
  implicit none

#if (defined USE_HDFEOS2)
#include "hdf.f90"
#endif
! !ARGUMENTS:   
  integer, intent(in)  :: n
  character(len=*)     :: name
! 
! !DESCRIPTION: 
! 
!EOP
  real               :: sb_rqc(NASA_AMSREsmobs%mo)

#if (defined USE_HDFEOS2)
  !declare the hdf-eos library functions 
  integer              :: gdopen,gdattach,gddefboxreg,gdrdfld
  integer              :: gdgetpix,gdextreg,gddetach,gdclose
  character*50         :: grid_name(2),sm_name(2),tm_name(2),qc_name(2)
  integer              :: ltime
  integer              :: ntype,rank,dims(2),size,igd
  integer              :: file_id,grid_id,region_id,ret
  integer*2,allocatable  :: smc(:)
  integer*2,allocatable  :: qc(:)
  real,allocatable       :: rqc(:)
  integer,parameter    :: ease_nr=586
  integer,parameter    :: ease_nc=1383
  real,allocatable     :: rsmc(:)
  real*8,allocatable   :: tm(:)
  real,allocatable     :: r4tm(:)
  real*8               :: upleftpt(2),lowrightpt(2)
  real*8               :: cornerlon(2),cornerlat(2)
  integer              :: mi,iret
  integer              :: start(2),edge(2),stride(2)
  logical*1,allocatable  :: li(:)
  logical*1            :: lo(NASA_AMSREsmobs%mo)
  real                 :: udef
  integer              :: ir,ic,ij,i
  integer              :: im,jm
  integer              :: t

  !Grid and field names
  grid_name(1) ="Ascending_Land_Grid"
  grid_name(2) ="Descending_Land_Grid"
  sm_name(1)   ="A_Soil_Moisture"
  sm_name(2)   ="D_Soil_Moisture"
  tm_name(1)   ="A_Time"
  tm_name(2)   ="D_Time"
  qc_name(1)   ="A_Inversion_QC_Flag"
  qc_name(2)   ="D_Inversion_QC_Flag"

  !open the hdf file

  file_id = gdopen(trim(name),DFACC_READ)
  if (file_id.eq.-1)then
     write(LDT_logunit,*)"Failed to open hdf file",name
     stop
     return
  end if
  
  mi = ease_nr*ease_nc
  allocate(li(mi))
  do igd=1,2
     
    !get the grid id
     grid_id = gdattach(file_id,grid_name(igd))
     if (grid_id.eq.-1)then
        write(LDT_logunit,*)"Failed to attach grid: ",grid_name(igd),name
        ret = gdclose(file_id)
        deallocate(li)
        return
     end if
     
     !retrieve the entire global grid
     start(1)=0  !hdfeos lib uses 0-based count
     start(2)=0
     edge(1)=ease_nc
     edge(2)=ease_nr
     stride(1)=1
     stride(2)=1
     allocate(tm(ease_nc*ease_nr))
     ret = gdrdfld(grid_id,tm_name(igd),start,stride,edge,tm)
     if (ret <0)then
        write(LDT_logunit,*)"Failed to get the time field"
        ret=gddetach(grid_id)
        ret=gdclose(file_id)
        deallocate(tm)
        deallocate(li)
        return
     end if
     
     !project to the subset lat/lon domain
     li=.false.
     udef=-9999.0
     do t=1,mi
        if(tm(t).gt.0.and.tm(t).ne.9999.0) then 
           li(t)=.true.
        endif
     enddo

     !convert double precision time to single precision
     allocate(r4tm(ease_nc*ease_nr))
     do i=1,mi
        r4tm(i)=tm(i)
     end do
     deallocate(tm)

     NASA_AMSREsmobs%smtime(:,igd) = 0

     call neighbor_interp(LDT_rc%gridDesc,li,r4tm,&
          lo,NASA_AMSREsmobs%smtime(:,igd),mi,NASA_AMSREsmobs%mo, &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          NASA_AMSREsmobs%n112,udef,iret)
     deallocate(r4tm)
!     call mask_tm(NASA_AMSREsmobs%smtime,NASA_AMSREsmobs%mo,ltime)
     
     !if found matching time
!     if (ltime==1)then
        !get smc
     allocate(smc(ease_nc*ease_nr))
     ret = gdrdfld(grid_id,sm_name(igd),start,stride,edge,smc)
     if (ret <0)then
        write(LDT_logunit,*)"Failed to get the smc field"
        ret=gddetach(grid_id)
        ret=gdclose(file_id)
        deallocate(smc)
        deallocate(li)
        return
     end if
        
        !convert short int smc to real 
     allocate(rsmc(ease_nc*ease_nr))
     do i=1,ease_nc*ease_nr
        rsmc(i)=smc(i)
     end do
     deallocate(smc)
     li=.false.
     do t=1,mi
        if(rsmc(t).gt.0.and.rsmc(t).ne.9999.0) then 
           li(t)=.true.
        endif
     enddo

     call neighbor_interp(LDT_rc%gridDesc,li,rsmc,&
          lo,NASA_AMSREsmobs%smobs(:,igd),mi,NASA_AMSREsmobs%mo, &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          NASA_AMSREsmobs%n112,udef,iret)

     deallocate(rsmc)

        !get qc
     allocate(qc(ease_nc*ease_nr))
     ret = gdrdfld(grid_id,qc_name(igd),start,stride,edge,qc)
     !convert to the real number to use the neighbor_interp call
     if (ret <0)then
        write(LDT_logunit,*)"Failed to get the qc field"
        ret=gddetach(grid_id)
        ret=gdclose(file_id)
        deallocate(qc)
        deallocate(li)
        return
     end if

     allocate(rqc(ease_nc*ease_nr))
     do i=1,ease_nc*ease_nr
        rqc(i)=qc(i)
     end do

     deallocate(qc)
     li=.false.
     do t=1,mi
        if(rqc(t).gt.0.and.rqc(t).ne.9999.0) li(t)=.true.
     enddo
     call neighbor_interp(LDT_rc%gridDesc,li,rqc,&
          lo,sb_rqc,mi,NASA_AMSREsmobs%mo, &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          NASA_AMSREsmobs%n112,udef,iret)
     deallocate(rqc) 
        !convert real qc back to integer for bit reading
     do i=1,NASA_AMSREsmobs%mo
        NASA_AMSREsmobs%smqc(i,igd)=nint(sb_rqc(i))
     end do
        
     ret=gddetach(grid_id)
     if (ret <0)then
        write(LDT_logunit,*)"Failed to detach grid_id: ",grid_id
     end if

  end do
  deallocate(li)
  ret=gdclose(file_id)
  if (ret <0)then
     write(LDT_logunit,*)"Failed to close file: ",file_id
  end if
#endif
  
end subroutine read_AMSREsm

!BOP
! !ROUTINE: maskObs_basedonQC
! \label{maskObs_basedonQC}
! 
! !INTERFACE: 
subroutine maskObs_basedonQC(npts, smc, qc)

  implicit none
!
! !ARGUMENTS: 
  integer        :: npts
  real           :: smc(npts,2)
  integer*2      :: qc(npts,2)
!
! !DESCRIPTION: 
!   This subroutine masks the observations based on the QC flags. 
!   data points with ice, dense vegetation and RFI are masked out. 
!
!    http://nsidc.org/data/docs/daac/ae_land3_l3_soil_moisture.gd.html
!EOP

  integer        :: k, i, iq

  do k=1,2
     do i = 1,npts
        if (smc(i,k)>0)then
           iq=qc(i,k)
              
           if (btest(iq,0).or.btest(iq,1).or.btest(iq,2)&
                .or.btest(iq,3).or.btest(iq,4).or.btest(iq,5)&
                .or.btest(iq,6))then
              smc(i,k)=-9999.0
           else
              smc(i,k)=smc(i,k)/1000.0
           end if
        else
           smc(i,k)=-9999.0
        end if
     end do
  enddo


end subroutine maskObs_basedonQC


!BOP
! 
! !ROUTINE: maskObs_basedonTime
! \label{maskObs_basedonTime}
! 
! !INTERFACE: 
subroutine maskObs_basedonTime(n, npts, smc, smtime, rsmc,fnd)
! !USES: 
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit
  use LDT_timeMgrMod, only : LDT_get_julhr

  implicit none

!
! !ARGUMENTS: 
  integer        :: n 
  integer        :: npts
  real           :: smc(npts,2)
  real           :: smtime(npts,2)
  real           :: rsmc(LDT_rc%lnc(n), LDT_rc%lnr(n))
  integer        :: fnd 
  real           :: rsmc1(npts)
  integer        :: c,r
! 
! !DESCRIPTION: 
!  This subroutine finds the observations in the AMSRE data that matches the 
!  current LDT timestep 
!EOP
  integer        :: k,i
  integer        :: ltime
  integer        :: lis_julhr,julhr1993
  real dt

  ltime=0
  !convert lis.time to sec since 01/01/1993 which is the begining time of the amsr-e measurements

  do k=1,2
     call LDT_get_julhr(LDT_rc%yr,LDT_rc%mo,LDT_rc%da,LDT_rc%hr,0,0,lis_julhr) 
     call LDT_get_julhr(1993,1,1,0,0,0,julhr1993)
     do i = 1,npts
        dt = smtime(i,k)+julhr1993*3600-(lis_julhr*3600+(LDT_rc%mn)*60+LDT_rc%ss)
        if ( smtime(i,k)>0.0 .and. dt >= 0 .and. dt <= LDT_rc%ts) then
           fnd = 1
           rsmc1(i) = smc(i,k)
        else
           rsmc1(i) = LDT_rc%udef
        end if
     end do
  enddo

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        rsmc(c,r) = rsmc1(c+(r-1)*LDT_rc%lnc(n)) 
     enddo
  enddo
    

end subroutine maskObs_basedonTime


!BOP
! !ROUTINE: NASA_AMSREsm_filename
! \label{NASA_AMSREsm_filename}
! 
! !INTERFACE: 
subroutine NASA_AMSREsm_filename(name, ndir, yr, mo,da)
! !USES:   
  use LDT_coreMod,only : LDT_rc
  use LDT_logMod, only : LDT_logunit

  implicit none
! !ARGUMENTS: 
  character(len=*)      :: name
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
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
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') LDT_rc%yr
  write(unit=fmo, fmt='(i2.2)') LDT_rc%mo
  write(unit=fda, fmt='(i2.2)') LDT_rc%da
  
  name = trim(ndir)//'/'//trim(fyr)//'/'//'AMSR_E_L3_DailyLand_V06_'&
         //trim(fyr)//trim(fmo)//trim(fda)//'.hdf'

end subroutine NASA_AMSREsm_filename
