!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: readNASMDObs
! \label{readNASMDObs}
!
! !INTERFACE: 
subroutine readNASMDObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_calendar, LVT_tick
  use LVT_logMod,       only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_endrun, LVT_verify
  use NASMD_obsMod,      only : nasmdobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer,       intent(in)   :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for NASMD station data, 
! to process the soil moisture and soil temperature measurements.  
! At the start of the program, an entire year's worth of 
! data is read and stored
! based on the time location and the station it corresponds to. 
! At future times, the read routine simply indexes into the right
! location. Depending upon the frequency of computing output statistics,
! the routine also computes time average (between different LIS 
! output intervals). The routine also interpolates the soil moisture
! and soil temperature data to the model's vertical resolution. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  16 Feb 2008: Sujay Kumar, Initial Specification
! 
!EOP
  integer             :: ftn
  integer             :: c,r,stn_col,stn_row
  real                :: col,row 
  character*10        :: id
  integer             :: i,t,k,kk,st,et,ios
  logical             :: readflag,file_exists
  character*100       :: nasmdfilename
  real                :: sf_wtsum, rz_wtsum
  integer             :: yr,mo,da,hr,mn,ss,doy
  real                :: gmt
  integer             :: status
  integer             :: nasmd_yr,nasmd_mo,nasmd_da,nasmd_hr,nasmd_mn,nasmd_ss
  real*8              :: lis_prevtime
  type(ESMF_Time)     :: startTime,nasmdtime,nasmdtime1,nasmdtime2,initTime
  type(ESMF_TimeInterval) :: dayInterval
  real,  allocatable      :: sm(:)
  real,  allocatable      :: sf_wt(:),rz_wt(:)
  real                :: sfsm, rzsm
  integer             :: obs_nlayers
  real,  allocatable      :: obs_d(:)
  real                :: smc(LVT_rc%lnc, LVT_rc%lnr)
  real                :: dummy(LVT_rc%lnc, LVT_rc%lnr)
  real                :: rootsm(LVT_rc%lnc, LVT_rc%lnr)

  smc     = LVT_rc%udef
  rootsm  = LVT_rc%udef

  yr = LVT_rc%dyr(source)
  mo = LVT_rc%dmo(source)
  da = LVT_rc%dda(source)
  hr = LVT_rc%dhr(source)
  mn = LVT_rc%dmn(source)
  ss = LVT_rc%dss(source) 

  if((nasmdobs(source)%yr.ne.LVT_rc%dyr(source)).or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     call ESMF_TimeSet(nasmdobs(source)%startTime,  yy=LVT_rc%dyr(source), &
          mm = 1, &
          dd = 1, &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, 'error in setting nasmd start time')

     nasmdobs(source)%yr = LVT_rc%dyr(source)
     nasmdobs(source)%sfsm = LVT_rc%udef
     nasmdobs(source)%rzsm = LVT_rc%udef

     do i=1,nasmdobs(source)%nstns
        call create_NASMD_filename(nasmdobs(source)%odir, LVT_rc%dyr(source), &
             nasmdobs(source)%stnid(i), nasmdfilename)
        
        inquire(file=trim(nasmdfilename), exist=file_exists) 
        
        obs_nlayers = nasmdobs(source)%stnnlayers(i)
        if(file_exists.and.obs_nlayers.gt.0) then 
           if(.not.allocated(obs_d)) then 
              allocate(obs_d(obs_nlayers))
           endif
           obs_d = real(nasmdobs(source)%stndepths(i,:))/100.0
           write(LVT_logunit,*) '[INFO] Reading NASMD file ',trim(nasmdfilename)
           
           ftn=LVT_getNextUnitNumber()
           open(ftn,file=trim(nasmdfilename),form='formatted')
           readflag = .true.         
           allocate(sm(obs_nlayers))
           allocate(sf_wt(obs_nlayers))
           allocate(rz_wt(obs_nlayers))
           sm = -9999.0
           do while(readflag) 
              read(ftn,fmt=*,iostat=ios) &
                   id, yr, mo, da, doy, &
                   (sm(k),k=1,obs_nlayers)

              if(ios.eq.0) then 
                 
                 call ESMF_TimeSet(nasmdtime, yy=yr,&
                      mm = mo, dd=da, h=0, m =0, s=0,&
                      calendar=LVT_calendar, &
                      rc=status)
                 call LVT_verify(status)
                 
                 t = nint((nasmdtime-nasmdobs(source)%starttime)/&
                      nasmdobs(source)%timestep) + 1
                 
                 call compute_vinterp_weights(&
                      obs_nlayers,LVT_rc%lis_sf_d,&
                      LVT_rc%lis_rz_d,&
                      obs_d,sf_wt,rz_wt)
                 sfsm = 0 
                 rzsm = 0
                 sf_wtsum = 0 
                 rz_wtsum = 0
                 do kk=1,obs_nlayers
                    if(sf_wt(kk).ne.0.and.sm(kk).lt.0.01) then 
                       sfsm = LVT_rc%udef
                       exit
                    else
                       sfsm = sfsm + sf_wt(kk)*sm(kk)
                       sf_wtsum = sf_wtsum + sf_wt(kk)
                    endif
                    if(rz_wt(kk).ne.0.and.sm(kk).lt.0.01) then 
                       rzsm = LVT_rc%udef
                       exit
                    else
                       rzsm = rzsm + rz_wt(kk)*sm(kk)
                       rz_wtsum = rz_wtsum + rz_wt(kk)
                    endif
                 enddo
                 
                 if(sfsm.ne.LVT_rc%udef.and.&
                      abs(sf_wtsum-1.0).lt.0.001) then 
                    nasmdobs(source)%sfsm(i,t) = sfsm
                 else
                    nasmdobs(source)%sfsm(i,t) = LVT_rc%udef
                 endif
                 if(rzsm.ne.LVT_rc%udef.and.&
                      abs(rz_wtsum-1.0).lt.0.001) then 
                    nasmdobs(source)%rzsm(i,t) = rzsm
                 else
                    nasmdobs(source)%rzsm(i,t) = LVT_rc%udef
                 endif                 
              endif
              if(ios.ne.0) then 
                 readflag = .false. 
              endif
           enddo
           deallocate(sm)
           deallocate(sf_wt)
           deallocate(rz_wt)
           call LVT_releaseUnitNumber(ftn)

        endif
     enddo
  endif
  
  call ESMF_TimeSet(nasmdtime1, yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
       h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'nasmdtime1 set failed')

  t = nint((nasmdtime1 - nasmdobs(source)%starttime)/&
       nasmdobs(source)%timestep)+1

  do i=1,nasmdobs(source)%nstns
     call latlon_to_ij(LVT_domain%lvtproj, nasmdobs(source)%stnlat(i), &
          nasmdobs(source)%stnlon(i),&
          col, row)
     stn_col = nint(col)
     stn_row = nint(row)
     
     if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
          stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
        
        if(nasmdobs(source)%sfsm(i,t).ne.LVT_rc%udef) then 
           if(nasmdobs(source)%sfsm(i,t).lt.0) then 
              print*, 'Error in sm ',stn_col, stn_row, i, t, &
                   nasmdobs(source)%sfsm(i,t), LVT_rc%udef
              stop
           endif
           smc(stn_col, stn_row) = &
                nasmdobs(source)%sfsm(i,t)
        endif
        if(nasmdobs(source)%rzsm(i,t).ne.LVT_rc%udef) then 
           rootsm(stn_col, stn_row) =&
                nasmdobs(source)%rzsm(i,t)
           if(nasmdobs(source)%rzsm(i,t).lt.0) then 
              print*, 'Error in rzsm ',stn_col, stn_row, i, t, &
                   nasmdobs(source)%rzsm(i,t)
              stop
           endif
        endif
     endif
  enddo
     
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(smc(c,r).gt.0.5) then !ignore
           smc(c,r) = LVT_rc%udef
        endif
        if(smc(c,r).lt.0.001) then !ignore
           smc(c,r) = LVT_rc%udef
        endif
        if(rootsm(c,r).gt.0.5) then !ignore
           rootsm(c,r) = LVT_rc%udef
        endif
        if(rootsm(c,r).lt.0.001) then !ignore
           rootsm(c,r) = LVT_rc%udef
        endif
     enddo
  enddo
  dummy = LVT_rc%udef

  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source,smc,&
       vlevel=1,units="m3/m3")

  do c=2,LVT_rc%nsmlayers
     call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source,&
          dummy,vlevel=c,units="m3/m3")
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_rootmoist, source,rootsm,&
       vlevel=1,units="m3/m3")
 
end subroutine readNASMDObs

!BOP
! 
! !ROUTINE: create_NASMD_filename
! \label{create_NASMD_filename}
!
! !INTERFACE: 
subroutine create_NASMD_filename(odir, yr, stnid, nasmdname)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
  character(len=*), intent(in)  :: stnid
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: nasmdname
!
! !DESCRIPTION: 
! 
! This routine creates a filename for the NASMD station
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      NASMD base directory
!   \item[yr]        year of data
!   \item[stnid]     Station ID 
!   \item[nasmdname]  Name of the NASMD file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4       :: fyr
  
  write(unit=fyr, fmt='(i4.4)') yr
  
  nasmdname = trim(odir)//'/'//trim(fyr)//'/'//trim(stnid)//&
       '_'//trim(fyr)//'.dat'
  
  
end subroutine create_NASMD_filename
