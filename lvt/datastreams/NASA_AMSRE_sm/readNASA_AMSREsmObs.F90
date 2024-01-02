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
! !ROUTINE: readNASA_AMSREsmObs
! \label{readNASA_AMSREsmObs}
!
! !INTERFACE: 
subroutine readNASA_AMSREsmObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit
  use NASA_AMSREsm_obsMod, only : NASA_AMSREsmobs

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

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  if(NASA_AMSREsmobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     NASA_AMSREsmobs(source)%startflag = .false. 
     call NASA_AMSREsm_filename(source,name,NASA_AMSREsmobs(source)%odir, & 
        LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source))
                 
     inquire(file=name, exist=file_exists) 

     if(file_exists) then 
        readflag = .true. 
     else
        readflag = .false. 
     endif
     
     if(readflag) then 
        write(LVT_logunit,*) '[INFO] Reading NASA AMSRE file ',name
        call read_AMSREsm(source, name)
        call maskObs_basedonQC(LVT_rc%lnc*LVT_rc%lnr, &
             NASA_AMSREsmobs(source)%smobs, &
             NASA_AMSREsmobs(source)%smqc)
     endif
  endif

  fnd = 0 
  smc = LVT_rc%udef
  call maskObs_basedonTime(source,LVT_rc%lnc*LVT_rc%lnr, &
       NASA_AMSREsmobs(source)%smobs, &
       NASA_AMSREsmobs(source)%smtime, smc,fnd)

  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source,smc,vlevel=1,units="m3/m3")
 
end subroutine readNASA_AMSREsmObs

!BOP
! 
! !ROUTINE: read_AMSREsm
! \label{read_AMSREsm}
!
! !INTERFACE: 
subroutine read_AMSREsm(source, name)
! 
! !USES: 
  use LVT_coreMod,         only : LVT_rc,LVT_domain
  use LVT_logmod,          only : LVT_logunit
  use NASA_AMSREsm_obsMod, only : NASA_AMSREsmobs
  implicit none
#if (defined USE_HDFEOS2)
#include "hdf.f90"
#endif
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
!BOP
! !ARGUMENTS:   
  integer              :: source 
  character(len=*)     :: name
!EOP
  real                 :: sb_rqc(NASA_AMSREsmobs(source)%mo)

#if (defined USE_HDFEOS2)
  !declare the hdf-eos library functions 
  integer              :: gdopen,gdattach,gddefboxreg,gdrdfld
  integer              :: gdgetpix,gdextreg,gddetach,gdclose
  character*50         :: grid_name(2),sm_name(2),tm_name(2),qc_name(2)
  integer              :: ltime
  integer              :: ntype,rank,dims(2),size,igd
  integer              :: file_id,grid_id,region_id,ret
  integer*2,allocatable    :: smc(:)
  integer*2,allocatable    :: qc(:)
  real,allocatable         :: rqc(:)
  integer,parameter    :: ease_nr=586
  integer,parameter    :: ease_nc=1383
  real,allocatable     :: rsmc(:)
  real*8,allocatable   :: tm(:)
  real,allocatable     :: r4tm(:)
  real*8               :: upleftpt(2),lowrightpt(2)
  real*8               :: cornerlon(2),cornerlat(2)
  integer              :: mi,iret
  integer              :: start(2),edge(2),stride(2)
  logical*1,allocatable    :: li(:)
  logical*1            :: lo(NASA_AMSREsmobs(source)%mo)
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
     write(LVT_logunit,*) "[ERR] Failed to open hdf file",name
     stop
     return
  end if
  
  mi = ease_nr*ease_nc
  allocate(li(mi))
  do igd=1,2
     
    !get the grid id
     grid_id = gdattach(file_id,grid_name(igd))
     if (grid_id.eq.-1)then
        write(LVT_logunit,*) "[ERR] Failed to attach grid: ",grid_name(igd),name
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
        write(LVT_logunit,*) "[ERR] Failed to get the time field"
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
        if(tm(t).gt.0) li(t)=.true.
     enddo
     !convert double precision time to single precision
     allocate(r4tm(ease_nc*ease_nr))
     do i=1,mi
        r4tm(i)=tm(i)
     end do
     deallocate(tm)

     NASA_AMSREsmobs(source)%smtime(:,igd) = 0

     call neighbor_interp(LVT_rc%gridDesc,li,r4tm,&
          lo,NASA_AMSREsmobs(source)%smtime(:,igd),mi,NASA_AMSREsmobs(source)%mo, &
          NASA_AMSREsmobs(source)%rlat2,NASA_AMSREsmobs(source)%rlon2,&
          NASA_AMSREsmobs(source)%n112,udef,iret)
     deallocate(r4tm)
!     call mask_tm(NASA_AMSREsmobs(source)%smtime,NASA_AMSREsmobs(source)%mo,ltime)
     
     !if found matching time
!     if (ltime==1)then
        !get smc
     allocate(smc(ease_nc*ease_nr))
     ret = gdrdfld(grid_id,sm_name(igd),start,stride,edge,smc)
     if (ret <0)then
        write(LVT_logunit,*) "[ERR] Failed to get the smc field"
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
        if(rsmc(t).gt.0.and.rsmc(t).ne.9999.0) li(t)=.true.
     enddo
     call neighbor_interp(LVT_rc%gridDesc,li,rsmc,&
          lo,NASA_AMSREsmobs(source)%smobs(:,igd),mi,NASA_AMSREsmobs(source)%mo, &
          NASA_AMSREsmobs(source)%rlat2,NASA_AMSREsmobs(source)%rlon2,&
          NASA_AMSREsmobs(source)%n112,udef,iret)
     
     deallocate(rsmc)

        !get qc
     allocate(qc(ease_nc*ease_nr))
     ret = gdextreg(grid_id,region_id,qc_name(igd),qc)
     !convert to the real number to use the neighbor_interp call
     allocate(rqc(ease_nc*ease_nr))
     do i=1,ease_nc*ease_nr
        rqc(i)=qc(i)
     end do
     deallocate(qc)
     li=.false.
     do t=1,mi
        if(rqc(t).gt.0) li(t)=.true.
     enddo
     call neighbor_interp(LVT_rc%gridDesc,li,rqc,&
          lo,sb_rqc,mi,NASA_AMSREsmobs(source)%mo, &
          NASA_AMSREsmobs(source)%rlat2,NASA_AMSREsmobs(source)%rlon2,&
          NASA_AMSREsmobs(source)%n112,udef,iret)
     deallocate(rqc) 
        !convert real qc back to integer for bit reading
     do i=1,NASA_AMSREsmobs(source)%mo
        NASA_AMSREsmobs(source)%smqc(i,igd)=nint(sb_rqc(i))
     end do
        
     ret=gddetach(grid_id)
     if (ret <0)then
        write(LVT_logunit,*) "[ERR] Failed to detach grid_id: ",grid_id
     end if

  end do
  deallocate(li)
  ret=gdclose(file_id)
  if (ret <0)then
     write(LVT_logunit,*) "[ERR] Failed to close file: ",file_id
  end if
#endif
  
end subroutine read_AMSREsm

!BOP
! 
! !ROUTINE: maskObs_basedonQC
! \label{maskObs_basedonQC}
!
! !INTERFACE: 
subroutine maskObs_basedonQC(npts, smc, qc)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine masks the observations based on the QC flags. 
!   data points with ice, dense vegetation and RFI are masked out. 
!
!    http://nsidc.org/data/docs/daac/ae_land3_l3_soil_moisture.gd.html
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer        :: npts
  real           :: smc(npts,2)
  integer*2      :: qc(npts,2)
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
subroutine maskObs_basedonTime(source, npts, smc, smtime, rsmc,fnd)
! 
! !USES: 
  use LVT_coreMod,  only : LVT_rc
  use LVT_logMod,   only : LVT_logunit
  use LVT_timeMgrMod, only : LVT_get_julhr

  implicit none
!
! !AGRUMENTS: 
  integer,        intent(in) :: source
  integer        :: npts
  real           :: smc(npts,2)
  real           :: smtime(npts,2)
  real           :: rsmc(LVT_rc%lnc, LVT_rc%lnr)
  integer        :: fnd 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine finds the observations in the AMSRE data that matches the 
!  current LVT timestep 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  real           :: rsmc1(npts)
  integer        :: c,r
!EOP
  integer        :: k,i
  integer        :: ltime
  integer        :: lis_julhr,julhr1993
  real dt
  ltime=0
  !convert lis.time to sec since 01/01/1993 which is the begining time of the amsr-e measurements
  do k=1,2
     
     call LVT_get_julhr(LVT_rc%dyr(source),LVT_rc%dmo(source),LVT_rc%dda(source),LVT_rc%dhr(source),0,0,lis_julhr) 
     call LVT_get_julhr(1993,1,1,0,0,0,julhr1993)
     do i = 1,npts
        dt = smtime(i,k)+julhr1993*3600-(lis_julhr*3600+(LVT_rc%dmn(source))*60+LVT_rc%dss(source))
        if ( smtime(i,k)>0.0 .and. dt >= 0 .and. dt <= LVT_rc%ts) then
           fnd = 1
           rsmc1(i) = smc(i,k)
        else
           rsmc1(i) = LVT_rc%udef
        end if
     end do
  enddo

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        rsmc(c,r) = rsmc1(c+(r-1)*LVT_rc%lnc) 
     enddo
  enddo
    

end subroutine maskObs_basedonTime


!BOP
! 
! !ROUTINE: NASA_AMSREsm_filename
! \label{NASA_AMSREsm_filename}
!
! !INTERFACE: 
subroutine NASA_AMSREsm_filename(source, name, ndir, yr, mo,da)
! 
! !USES:   
  use LVT_coreMod,only : LVT_rc
  use LVT_logMod, only : LVT_logunit

  implicit none
!
! !ARGUMENTS: 
  integer            :: source
  character*200      :: name
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
  
  write(unit=fyr, fmt='(i4.4)') LVT_rc%dyr(source)
  write(unit=fmo, fmt='(i2.2)') LVT_rc%dmo(source)
  write(unit=fda, fmt='(i2.2)') LVT_rc%dda(source)
  
  name = trim(ndir)//'/'//trim(fyr)//'/'//'AMSR_E_L3_DailyLand_V06_'&
         //trim(fyr)//trim(fmo)//trim(fda)//'.hdf'

end subroutine NASA_AMSREsm_filename
