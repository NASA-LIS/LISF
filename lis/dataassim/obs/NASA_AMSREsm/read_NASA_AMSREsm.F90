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
! !ROUTINE: read_NASA_AMSREsm
! \label{read_NASA_AMSREsm}
!
! !REVISION HISTORY:
!  06Nov2007: Bailing Li; Initial Specification
!  17Jun2010: Sujay Kumar; Updated for use with NASA AMSRE Version 6, 
!                          added LSM-based QC, generic handling of obs scaling
!
! !INTERFACE: 
subroutine read_NASA_AMSREsm(n, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod,    only : LIS_rc, LIS_domain, &
       LIS_masterproc, LIS_npes, LIS_masterproc, LIS_localPet
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod,     only : LIS_logunit, LIS_verify
  use LIS_timeMgrMod, only : LIS_calendar, LIS_clock
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, & 
       LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
  use NASA_AMSREsm_Mod, only : NASA_AMSREsm_struc
  use LIS_pluginIndices 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the amsr-e soil moisture observations 
!  from HDF-EOS files and packages it 
!  into an ESMF State with certain predefined 
!  attributes
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  real, parameter     :: MAX_SM_VALUE=0.55
  type(ESMF_Field)    :: smField

  integer             :: iret

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%ngrid(n))
  integer             :: assimflag(LIS_rc%ngrid(n))

  real                :: smobs(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                :: obs_unsc(LIS_rc%ngrid(n))
  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  logical             :: data_update
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name

  logical             :: alarmCheck

  logical             :: readflag
  integer             :: status
  logical             :: data_upd_flag(LIS_npes)
  logical             :: data_upd_flag_local
  logical             :: data_upd
  integer             :: t,c,r,p

  integer             :: col, row, i, gridid
  integer             :: binval
  real                :: cdf_obsval
  real                :: smvalue
  real                :: model_delta(LIS_rc%ngrid(n))
  real                :: obs_delta(LIS_rc%ngrid(n))

  integer             :: fnd
  integer             :: ierr


  smobs    = LIS_rc%udef
  obs_unsc = LIS_rc%udef 

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)
  data_upd = .false. 
!-------------------------------------------------------------------------
!   Read the data at 0z
!-------------------------------------------------------------------------

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "NASA AMSR-E read alarm")
  if(alarmCheck) then 
     call NASA_AMSREsm_filename(name,smobsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)   
     
     inquire(file=name,exist=file_exists)
     if(file_exists) then 
        readflag = .true. 
     else 
        readflag = .false.
     endif

     if (readflag) then 
        write(LIS_logunit,*) 'Reading NASA AMSRE file ',trim(name)

        call read_AMSREsm(n,name)
        
        call maskObs_basedonQC(LIS_rc%lnc(n)*LIS_rc%lnr(n), &
             NASA_AMSREsm_struc(n)%smobs, &
             NASA_AMSREsm_struc(n)%smqc)
     endif
  endif

  call ESMF_StateGet(OBS_State,"Observation01",smField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)
  
  obsl = LIS_rc%udef
  fnd = 0 
  call maskObs_basedonTime(LIS_rc%lnc(n)*LIS_rc%lnr(n), &
       NASA_AMSREsm_struc(n)%smobs, &
       NASA_AMSREsm_struc(n)%smtime, smobs,fnd)

  if(fnd.eq.0) then 
     obsl = LIS_rc%udef
     
  else     

  !map tmp_obsl to obsl.
     do r =1,LIS_rc%lnr(n)
        do c =1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r) .ne. -1)then
              obsl(LIS_domain(n)%gindex(c,r))=smobs(c+LIS_rc%lnc(n)*(r-1))
           end if
        end do
     end do

  endif

!lsm-based qc

  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_NASA_AMSREsmobsId)//char(0), &
       n, OBS_state)

  call ESMF_StateGet(OBS_State,"Observation01",smField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)

  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           smobs(c+LIS_rc%lnc(n)*(r-1)) = obsl(LIS_domain(n)%gindex(c,r))
        end if
     end do
  end do
  

!-------------------------------------------------------------------------
!  Store unscaled data for future output
!-------------------------------------------------------------------------     
  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           obs_unsc(LIS_domain(n)%gindex(c,r))=smobs(c+LIS_rc%lnc(n)*(r-1))
        end if
     end do
  end do
        

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     
  if(NASA_AMSREsm_struc(n)%scal.ne.0) then 
     
     do t=1,LIS_rc%ngrid(n)
        model_delta(t) = NASA_AMSREsm_struc(n)%model_xrange(t,2)-&
             NASA_AMSREsm_struc(n)%model_xrange(t,1)
        obs_delta(t) = NASA_AMSREsm_struc(n)%obs_xrange(t,2)-&
             NASA_AMSREsm_struc(n)%obs_xrange(t,1)
     enddo
     do t=1,LIS_rc%ngrid(n)
        
        col = LIS_domain(n)%grid(t)%col
        row = LIS_domain(n)%grid(t)%row
        
        gridid = col+(row-1)*LIS_rc%lnc(n)
        if(smobs(gridid).ne.-9999.0) then 
           if(obs_delta(t).gt.0) then 
              binval = nint((smobs(gridid)-NASA_AMSREsm_struc(n)%obs_xrange(t,1))/&
                   obs_delta(t))+1
              if(binval.gt.NASA_AMSREsm_struc(n)%nbins) binval = NASA_AMSREsm_struc(n)%nbins
              if(binval.le.0) binval = 1
              cdf_obsval = NASA_AMSREsm_struc(n)%obs_cdf(t,binval)
              if(cdf_obsval.gt.1.0) cdf_obsval = 1.0
              i=1
              do while((NASA_AMSREsm_struc(n)%model_cdf(t,i).lt.cdf_obsval).and.&
                   (i.le.NASA_AMSREsm_struc(n)%nbins))
                 i = i+1
              enddo
              if(i.gt.NASA_AMSREsm_struc(n)%nbins) i = i-1
              smvalue = NASA_AMSREsm_struc(n)%model_xrange(t,i)
              
              if(smvalue.gt.MAX_SM_VALUE) then 
                 write(LIS_logunit,*) 'Problem in scaling SM observations'
                 write(LIS_logunit,*) t, smobs(gridid),smvalue, MAX_SM_VALUE
                 call LIS_endrun()
              endif
              
              if(smvalue.eq.0) then 
                 smvalue = smobs(gridid)
              endif
              smobs(gridid) = smvalue
           else
              print*, 'WARNING: obs_delta for grid point', t, '=0'
              smobs(gridid) = LIS_rc%udef
           endif
        else
           smobs(gridid) = LIS_rc%udef
        endif
     enddo
  endif


  fnd =0 
  data_upd_flag_local = .false.
  do t=1,LIS_rc%lnc(n)*LIS_rc%lnr(n)
     if(smobs(t).ne.LIS_rc%udef) then 
        fnd = 1
     endif        
  enddo

  call ESMF_StateGet(OBS_State,"Observation01",smField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)
  
  obsl = LIS_rc%udef

  if(fnd.eq.0) then 
     obsl = LIS_rc%udef
     
  else     
     do r =1,LIS_rc%lnr(n)
        do c =1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r) .ne. -1)then
              obsl(LIS_domain(n)%gindex(c,r))=smobs(c+LIS_rc%lnc(n)*(r-1))
           end if
        end do
     end do
  endif

  if(fnd.eq.0) then 
     data_upd_flag_local = .false. 
  else
     data_upd_flag_local = .true. 
  endif


#if (defined SPMD)
  call MPI_ALLGATHER(data_upd_flag_local,1, &
       MPI_LOGICAL, data_upd_flag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, ierr)
#endif
  data_upd = .false.
  do p=1,LIS_npes
     data_upd = data_upd.or.data_upd_flag(p)
  enddo
  
  if(data_upd) then 
     do t=1,LIS_rc%ngrid(n)
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

     call ESMF_AttributeSet(smField,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)
  
     call ESMF_AttributeSet(smField,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(smfield, "Unscaled Obs",&
          obs_unsc, itemCount=LIS_rc%ngrid(n), rc=status)
     call LIS_verify(status, 'Error in setting Unscaled Obs attribute')       

  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif
end subroutine read_NASA_AMSREsm

!BOP
! !ROUTINE: read_AMSREsm
! \label{read_AMSREsm}
! 
! !INTERFACE: 
subroutine read_AMSREsm(n,name)
! !USES: 
  use LIS_coreMod,      only : LIS_rc,LIS_domain
  use LIS_logmod,       only : LIS_logunit
  use NASA_AMSREsm_Mod, only : NASA_AMSREsm_struc
  implicit none

#if (defined USE_HDF4) 
#include "hdf.f90"
#endif
! !ARGUMENTS:   
  integer, intent(in)  :: n 
  character(len=*)     :: name
! 
! !DESCRIPTION: 
! 
!EOP
  real                 :: sb_rqc(NASA_AMSREsm_struc(n)%mo)

#if (defined USE_HDF4)
  !declare the hdf-eos library functions 
  integer              :: gdopen,gdattach,gdrdfld
  integer              :: gddetach,gdclose
  character*50         :: grid_name(2),sm_name(2),tm_name(2),qc_name(2)
  integer              :: igd
  integer              :: file_id,grid_id,ret
  integer*2,allocatable    :: smc(:)
  integer*2,allocatable    :: qc(:)
  real,allocatable         :: rqc(:)
  integer,parameter    :: ease_nr=586
  integer,parameter    :: ease_nc=1383
  real,allocatable     :: rsmc(:)
  real*8,allocatable   :: tm(:)
  real,allocatable     :: r4tm(:)
  integer              :: mi,iret
  integer              :: start(2),edge(2),stride(2)
  logical*1,allocatable    :: li(:)
  logical*1            :: lo(NASA_AMSREsm_struc(n)%mo)
  real                 :: udef
  integer              :: i
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
     write(LIS_logunit,*)"Failed to open hdf file",trim(name)
     return
  end if
  
  mi = ease_nr*ease_nc
  allocate(li(mi))
  do igd=1,2
     
    !get the grid id
     grid_id = gdattach(file_id,grid_name(igd))
     if (grid_id.eq.-1)then
        write(LIS_logunit,*)"Failed to attach grid: ",grid_name(igd),trim(name)
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
        write(LIS_logunit,*)"Failed to get the time field"
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
        if(tm(t).gt.0.and.tm(t).ne.9999.0) li(t)=.true.
     enddo
     !convert double precision time to single precision
     allocate(r4tm(ease_nc*ease_nr))
     do i=1,mi
        r4tm(i)=tm(i)
     end do
     deallocate(tm)

     NASA_AMSREsm_struc(n)%smtime(:,igd) = 0
     call neighbor_interp(LIS_rc%gridDesc(n,:),li,r4tm,&
          lo,NASA_AMSREsm_struc(n)%smtime(:,igd),mi,NASA_AMSREsm_struc(n)%mo, &
          LIS_domain(n)%lat,LIS_domain(n)%lon,&
          NASA_AMSREsm_struc(n)%n112,udef,iret)
     deallocate(r4tm)
!     call mask_tm(NASA_AMSREsm_struc(n)%smtime,NASA_AMSREsm_struc(n)%mo,ltime)
     
     !if found matching time
!     if (ltime==1)then
        !get smc
     allocate(smc(ease_nc*ease_nr))
     ret = gdrdfld(grid_id,sm_name(igd),start,stride,edge,smc)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to get the smc field"
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
     call neighbor_interp(LIS_rc%gridDesc(n,:),li,rsmc,&
          lo,NASA_AMSREsm_struc(n)%smobs(:,igd),mi,NASA_AMSREsm_struc(n)%mo, &
          LIS_domain(n)%lat,LIS_domain(n)%lon,&
          NASA_AMSREsm_struc(n)%n112,udef,iret)
     
     deallocate(rsmc)

        !get qc
     allocate(qc(ease_nc*ease_nr))
!     ret = gdextreg(grid_id,region_id,qc_name(igd),qc)
     ret = gdrdfld(grid_id, qc_name(igd), start, stride, edge, qc)
     !convert to the real number to use the neighbor_interp call
     allocate(rqc(ease_nc*ease_nr))
     do i=1,ease_nc*ease_nr
        rqc(i)=qc(i)
     end do
     deallocate(qc)
     li=.false.
     do t=1,mi
        if(rqc(t).gt.0.and.rqc(t).ne.9999.0) li(t)=.true.
     enddo
     call neighbor_interp(LIS_rc%gridDesc(n,:),li,rqc,&
          lo,sb_rqc,mi,NASA_AMSREsm_struc(n)%mo, &
          LIS_domain(n)%lat,LIS_domain(n)%lon,&
          NASA_AMSREsm_struc(n)%n112,udef,iret)
     deallocate(rqc) 
        !convert real qc back to integer for bit reading
     do i=1,NASA_AMSREsm_struc(n)%mo
        NASA_AMSREsm_struc(n)%smqc(i,igd)=nint(sb_rqc(i))
     end do
        
        !mask the smc based on qc and convert smc to 0~0.5
!        call mask_obs(NASA_AMSREsm_struc(n)%smtime,NASA_AMSREsm_struc(n)%smobs,NASA_AMSREsm_struc(n)%smqc,NASA_AMSREsm_struc(n)%mo)
        
!        write(LIS_logunit,*)"NASA_AMSREsm_struc(n)%smobs after masking",NASA_AMSREsm_struc(n)%smobs       
!     deallocate(li)
     ret=gddetach(grid_id)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to detach grid_id: ",grid_id
     end if

  end do
  deallocate(li)
  ret=gdclose(file_id)
  if (ret <0)then
     write(LIS_logunit,*)"Failed to close file: ",file_id
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
!    http://nsidc.org/data/docs/daac/ae\_land3\_l3\_soil\_moisture.gd.html
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
subroutine maskObs_basedonTime(npts, smc, smtime, rsmc,fnd)
! !USES: 
  use LIS_coreMod,  only : LIS_rc
  use LIS_logMod,   only : LIS_logunit
  use LIS_timeMgrMod, only : LIS_get_julhr

  implicit none

!
! !ARGUMENTS: 
  integer        :: npts
  real           :: smc(npts,2)
  real           :: smtime(npts,2)
  real           :: rsmc(npts)
  integer        :: fnd 
! 
! !DESCRIPTION: 
!  This subroutine finds the observations in the AMSRE data that matches the 
!  current LIS timestep 
!EOP
  integer        :: k,i
  integer        :: ltime
  integer        :: lis_julhr,julhr1993
  real dt
  ltime=0
  !convert lis.time to sec since 01/01/1993 which is the begining time of the amsr-e measurements
  do k=1,2
     
     call LIS_get_julhr(LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,0,0,lis_julhr) 
     call LIS_get_julhr(1993,1,1,0,0,0,julhr1993)
     do i = 1,npts
        dt = smtime(i,k)+julhr1993*3600-(lis_julhr*3600+(LIS_rc%mn)*60+LIS_rc%ss)
        if ( smtime(i,k)>0.0 .and. dt >= 0 .and. dt <= LIS_rc%ts) then          
!do not overwrite ascending and descending passes 
           if(rsmc(i).eq.LIS_rc%udef) then 
              fnd = 1
              rsmc(i) = smc(i,k)
           endif
        else
           rsmc(i) = LIS_rc%udef
        end if
     end do
  enddo

end subroutine maskObs_basedonTime


!BOP
! !ROUTINE: NASA_AMSREsm_filename
! \label{NASA_AMSREsm_filename}
! 
! !INTERFACE: 
subroutine NASA_AMSREsm_filename(name, ndir, yr, mo,da)
! !USES:   
  use LIS_coreMod,only : LIS_rc
  use LIS_logMod, only : LIS_logunit

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: name
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
  
  write(unit=fyr, fmt='(i4.4)') LIS_rc%yr
  write(unit=fmo, fmt='(i2.2)') LIS_rc%mo
  write(unit=fda, fmt='(i2.2)') LIS_rc%da
  
  name = trim(ndir)//'/'//trim(fyr)//'/'//'AMSR_E_L3_DailyLand_V06_'&
         //trim(fyr)//trim(fmo)//trim(fda)//'.hdf'

end subroutine NASA_AMSREsm_filename




