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
! !ROUTINE: read_AMSRE_SWE
! \label{read_AMSRE_SWE}
!
! !REVISION HISTORY:
!  01 Jul 2010: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_AMSRE_SWE(n, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_coreMod,    only : LIS_rc, LIS_domain, &
       LIS_masterproc, LIS_npes, LIS_masterproc, LIS_localPet
  use LIS_logMod,     only : LIS_logunit, LIS_verify
  use LIS_timeMgrMod, only : LIS_calendar, LIS_clock
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, & 
       LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use AMSRE_SWE_Mod, only : AMSRE_SWE_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the Level 3 AMSR-E snow water equivalent observations 
!  from the HDF-EOS files and packages into an ESMF state object. 
!  The routine reads the snow data at 0z, performs spatial 
!  interpolation to the LIS grid and keeps it in memory. 
!  At 10.30 AM localtime for each grid point, the code then 
!  packages the interpolated observations into an ESMF state object. 
!  This routine handles the AMSR-E SWE retrievals available in the 
!  HDF-EOS format. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: sweField

  integer             :: iret

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%ngrid(n))
  integer             :: assimflag(LIS_rc%ngrid(n))

  character(len=LIS_CONST_PATH_LEN) :: sweobsdir, name
  logical             :: data_update
  logical             :: file_exists

  logical             :: alarmCheck

  logical             :: readflag
  integer             :: status
  logical             :: dataflag(LIS_npes)
  logical             :: dataflag_local
  logical             :: data_upd
  integer             :: t,c,r,p
  real                :: lon, lhour
  integer             :: zone
  integer             :: ierr


  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sweobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 

!-------------------------------------------------------------------------
!   Read the data at 0z
!-------------------------------------------------------------------------

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "AMSR-E SWE read alarm")  

  if(alarmCheck) then 

     AMSRE_SWE_struc(n)%sweobs = LIS_rc%udef
     AMSRE_SWE_struc(n)%sweqc = LIS_rc%udef
     
     call AMSRE_SWE_filename(name,sweobsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)   
     
     inquire(file=name,exist=file_exists)
     if(file_exists) then 
        readflag = .true. 
     else 
        readflag = .false.
     endif

     if (readflag) then 
        write(LIS_logunit,*) 'Reading NASA AMSRE file ',trim(name)
        call read_AMSREswe(n,name)
        
        call maskSWEobs_basedonQC(LIS_rc%lnc(n)*LIS_rc%lnr(n), &
             AMSRE_SWE_struc(n)%sweobs, &
             AMSRE_SWE_struc(n)%sweqc)
        
     endif
  endif

  call ESMF_StateGet(OBS_State,"Observation01",sweField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)
  
  obsl = LIS_rc%udef 
  dataflag_local = .false. 
  
  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           
!localtime of this gridcell
           lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(c,r))%lon
           call LIS_localtime(LIS_rc%gmt,lon,lhour,zone)
           
           if(lhour.eq.10.5) then            
              obsl(LIS_domain(n)%gindex(c,r))=&
                   AMSRE_SWE_struc(n)%sweobs(c+LIS_rc%lnc(n)*(r-1))
              if(obsl(LIS_domain(n)%gindex(c,r)).ne.-9999.0) then 
                 dataflag_local = .true.
              endif
           endif
           
        end if
     end do
  end do

#if (defined SPMD)
  call MPI_ALLGATHER(dataflag_local,1, MPI_LOGICAL, dataflag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, ierr)
#endif
  data_upd = .false.
  
  do p=1,LIS_npes
     data_upd = data_upd.or.dataflag(p)
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
          .true., rc=status)
     call LIS_verify(status, 'Error: AttributeSet in Data Update Status')
     
     call ESMF_AttributeSet(sweField,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status,'Error: AttributeSet in Grid Number')
     
     call ESMF_AttributeSet(sweField,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')
           
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status, "Error: AttributeSet Data Update Status")
     return
  end if
  
end subroutine read_AMSRE_SWE

!BOP
! !ROUTINE: read_AMSREswe
! \label{read_AMSREswe}
! 
! !INTERFACE: 
subroutine read_AMSREswe(n,name)
! !USES: 
  use LIS_coreMod,      only : LIS_rc,LIS_domain
  use LIS_logmod,       only : LIS_logunit
  use AMSRE_SWE_Mod, only : AMSRE_SWE_struc
  implicit none

#if (defined USE_HDF4) 
#include "hdf.f90"
#endif
! !ARGUMENTS:   
  integer, intent(in)  :: n 
  character(len=*)     :: name
! 
! !DESCRIPTION: 
!   This routine extracts the SWE retrievals and the associated 
!   quality control flags from the Level 3 AMSR-E SWE HDF-EOS files. 
! 
!EOP
  real                 :: sb_rqc(AMSRE_SWE_struc(n)%mo)
  real                 :: sweobs(AMSRE_SWE_struc(n)%mo)
#if (defined USE_HDF4)
  !declare the hdf-eos library functions 
  integer              :: gdopen,gdattach,gdrdfld
  integer              :: gddetach,gdclose
  character*50         :: grid_name(2),swe_name(2),qc_name(2)
  integer              :: size,igd
  integer              :: file_id,grid_id,ret
  integer*1,allocatable    :: swe(:)
  integer*1,allocatable    :: qc(:)
  real,allocatable         :: rqc(:)
  integer,parameter    :: ease_nr=721
  integer,parameter    :: ease_nc=721
  real,allocatable     :: rswe(:)
  integer              :: mi,iret
  integer              :: start(2),edge(2),stride(2)
  logical*1,allocatable    :: li(:)
  logical*1            :: lo(AMSRE_SWE_struc(n)%mo)
  integer              :: i
  integer              :: t

  sb_rqc = LIS_rc%udef
  sweobs = LIS_rc%udef

  !Grid and field names
  grid_name(1) ="Northern Hemisphere"
  grid_name(2) ="Southern Hemisphere"
  swe_name(1)   ="SWE_NorthernDaily"
  swe_name(2)   ="SWE_SouthernDaily"
  qc_name(1)   ="Flags_NorthernDaily"
  qc_name(2)   ="Flags_SouthernDaily"

  !open the hdf file

  file_id = gdopen(trim(name),DFACC_READ)
  if (file_id.eq.-1)then
     write(LIS_logunit,*)"Failed to open hdf file",trim(name)
     return
  end if
  
  mi = ease_nr*ease_nc
  allocate(li(mi))

  do igd=AMSRE_SWE_struc(n)%ihemi, AMSRE_SWE_struc(n)%nhemi     
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

     allocate(swe(ease_nc*ease_nr))
     ret = gdrdfld(grid_id,swe_name(igd),start,stride,edge,swe)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to get the swe field"
        ret=gddetach(grid_id)
        ret=gdclose(file_id)
        deallocate(swe)
        deallocate(li)
        return
     end if
        
        !convert short int swe to real 
     allocate(rswe(ease_nc*ease_nr))
     do i=1,ease_nc*ease_nr
        rswe(i)=swe(i)
     end do
     deallocate(swe)
     li=.false.
     do t=1,mi
        if(rswe(t).ge.0) then 
           li(t)=.true.
        endif
     enddo
     call neighbor_interp(LIS_rc%gridDesc(n,:),li,rswe,&
          lo,sweobs,mi,AMSRE_SWE_struc(n)%mo, &
          AMSRE_SWE_struc(n)%rlat2(:,igd),AMSRE_SWE_struc(n)%rlon2(:,igd),&
          AMSRE_SWE_struc(n)%n112(:,igd),LIS_rc%udef,iret)
     
     do t=1,LIS_rc%lnc(n)*LIS_rc%lnr(n)
        if(sweobs(t).ne.LIS_rc%udef) then
           AMSRE_SWE_struc(n)%sweobs(t) = sweobs(t)*2.0
        endif
     enddo

     !get qc
     allocate(qc(ease_nc*ease_nr))
     ret =  gdrdfld(grid_id,qc_name(igd),start,stride,edge,qc)
     !convert to the real number to use the neighbor_interp call
     allocate(rqc(ease_nc*ease_nr))
     do i=1,ease_nc*ease_nr
        rqc(i)=qc(i)
     end do
     deallocate(qc)
     li=.false.
     do t=1,mi
        if(rqc(t).gt.0) then 
           li(t)=.true.
        endif
     enddo

     call neighbor_interp(LIS_rc%gridDesc(n,:),li,rqc,&
          lo,sb_rqc,mi,AMSRE_SWE_struc(n)%mo, &
          AMSRE_SWE_struc(n)%rlat2(:,igd),AMSRE_SWE_struc(n)%rlon2(:,igd),&
          AMSRE_SWE_struc(n)%n112(:,igd),LIS_rc%udef,iret)

     deallocate(rswe)
     deallocate(rqc) 

        !convert real qc back to integer for bit reading
     do i=1,AMSRE_SWE_struc(n)%mo        
        if(sb_rqc(i).ne.LIS_rc%udef) then 
           AMSRE_SWE_struc(n)%sweqc(i)=sb_rqc(i)
        endif
     end do

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
  
end subroutine read_AMSREswe

!BOP
! !ROUTINE: maskSWEobs_basedonQC
! \label{maskSWEobs_basedonQC}
! 
! !INTERFACE: 
subroutine maskSWEobs_basedonQC(npts, swe, qc)

  implicit none
!
! !ARGUMENTS: 
  integer        :: npts
  real           :: swe(npts)
  real           :: qc(npts)
!
! !DESCRIPTION: 
!   This subroutine masks the observations based on the QC flags. 
!   non-validated, off-earth, ice-sheet, water, and missing data
!   points are excluded in this routine. 
! 
!EOP

  integer        :: i
  
  do i = 1,npts
     if (swe(i)>=0)then
        if(qc(i).eq.241.0.or.qc(i).eq.248.0.or.&
             qc(i).eq.252.or.qc(i).eq.253.0.or.&
             qc(i).eq.254.or.qc(i).eq.255.0) then 
           swe(i)=-9999.0
        else
           swe(i)=swe(i)
        end if
     else
        swe(i)=-9999.0
     end if
  end do

end subroutine maskSWEobs_basedonQC



!BOP
! !ROUTINE: AMSRE_SWE_filename
! \label{AMSRE_SWE_filename}
! 
! !INTERFACE: 
subroutine AMSRE_SWE_filename(name, ndir, yr, mo,da)
! !USES:   
  use LIS_coreMod,only : LIS_rc
  use LIS_logMod, only : LIS_logunit

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: name, ndir
  integer           :: yr, mo, da, hr,mn
! 
! !DESCRIPTION: 
!  This subroutine creates the Level 3 AMSRE SWE filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the AMSRE SWE filename
!  \item[ndir] name of the AMSRE SWE root directory
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
  
  name = trim(ndir)//'/'//trim(fyr)//'/'//'AMSR_E_L3_DailySnow_V09_'&
         //trim(fyr)//trim(fmo)//trim(fda)//'.hdf'

end subroutine AMSRE_SWE_filename




