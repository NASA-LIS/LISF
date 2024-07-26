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
! !ROUTINE: readGHCNObs
! \label{readGHCNObs}
!
! !INTERFACE: 
subroutine readGHCNObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use LVT_logMod
  use GHCN_obsMod
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for GHCN station data. 
! The code expects yearly GHCN files. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  16 Feb 2008: Sujay Kumar, Initial Specification
! 
!EOP

  integer                 :: status
  integer                 :: gh_yr, gh_mo, gh_da, gh_hr, gh_mn, gh_ss
  integer                 :: yr, mo, da, hr, mn, ss
  integer                 :: doy
  real                    :: gmt
  real*8                  :: lis_prevtime
  integer                 :: st,et
  integer                 :: i,t,c,r,kk
  integer                 :: stn_col, stn_row
  real                    :: col, row
  real                    :: offset  
  character*100           :: ghcnname
  type(ESMF_TimeInterval) :: dayInterval
  type(ESMF_Time)         :: startTime, initTime
  type(ESMF_Time)         :: ghcntime1, ghcntime2
  real                    :: snowdepth(LVT_rc%lnc,LVT_rc%lnr)
  integer                 :: nsnowdepth(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: prcp(LVT_rc%lnc,LVT_rc%lnr)
  integer                 :: nprcp(LVT_rc%lnc,LVT_rc%lnr)
  
  snowdepth  = 0.0
  nsnowdepth = 0.0

  prcp = 0
  nprcp = 0 
  
  if((ghcnobs(source)%yr.ne.LVT_rc%dyr(source)).or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     ghcnobs(source)%yr = LVT_rc%dyr(source)
     ghcnobs(source)%snod = LVT_rc%udef

     ghcnobs(source)%prcp = LVT_rc%udef
     
     call ESMF_TimeSet(ghcnobs(source)%startTime,  yy=LVT_rc%dyr(source), &
          mm = 1, &
          dd = 1, &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, 'error in setting ghcn time')
     
     call  create_ghcnsnwd_filename(ghcnobs(source)%odir, LVT_rc%dyr(source), &
          ghcnname)
     
     call read_ghcndata(source, ghcnname)  

  endif

  call ESMF_TimeSet(ghcntime1, yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), h=LVT_rc%dhr(source),&
       m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'ghcntime1 set failed')

  t = nint((ghcntime1-ghcnobs(source)%startTime)/ghcnobs(source)%timestep) + 1

  do i=1,ghcnobs(source)%nstns
     call latlon_to_ij(LVT_domain%lvtproj, &
          ghcnobs(source)%stnlat(i), ghcnobs(source)%stnlon(i),&
          col,row)
     stn_col = nint(col)
     stn_row = nint(row)
     
     if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
          stn_row.ge.1.and.stn_row.le.LVT_rc%lnr.and.&
          ghcnobs(source)%snod(i,t).gt.0) then 
        snowdepth(stn_col,stn_row) = &
             snowdepth(stn_col,stn_row) + & 
             ghcnobs(source)%snod(i,t)
        nsnowdepth(stn_col,stn_row) = & 
             nsnowdepth(stn_col,stn_row) + 1
     endif
     
     if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
          stn_row.ge.1.and.stn_row.le.LVT_rc%lnr.and.&
          ghcnobs(source)%prcp(i,t).gt.0) then 
        prcp(stn_col,stn_row) = &
             ghcnobs(source)%prcp(i,t)
     endif
  enddo

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nsnowdepth(c,r).gt.0) then 
           snowdepth(c,r) = snowdepth(c,r)/nsnowdepth(c,r)
        else
           snowdepth(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_snowdepth,source,snowdepth,vlevel=1,units="m")
  call LVT_logSingleDataStreamVar(LVT_MOC_rainf,source,prcp,vlevel=1,units='kg/m2')
  call LVT_logSingleDataStreamVar(LVT_MOC_rainfforc,source,prcp,vlevel=1,units='kg/m2')
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp,vlevel=1,units='kg/m2')

  do r=1, LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(prcp(c,r).ne.LVT_rc%udef) then 
           prcp(c,r) = prcp(c,r)/86400.0
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_rainf,source,prcp,vlevel=1,units='kg/m2s')
  call LVT_logSingleDataStreamVar(LVT_MOC_rainfforc,source,prcp,vlevel=1,units='kg/m2s')
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp,vlevel=1,units='kg/m2s')
end subroutine readGhcnobs

!BOP
! 
! !ROUTINE: read_ghcndata
! \label{read_ghcndata}
!
! !INTERFACE: 
subroutine read_ghcndata(source, filename)
! 
! !USES: 
  use ESMF
  use LVT_coreMod,    only : LVT_rc, LVT_domain
  use LVT_timeMgrMod, only : LVT_calendar
  use LVT_logMod,     only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_verify, LVT_endrun
  use GHCN_obsMod,    only : ghcnobs

  implicit none

  character(len=*)        :: filename
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads the GHCN data files and categorizes the data
!  based on station index and the temporal location of the data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 

  integer                 :: i   
  integer                 :: source 
  integer          :: ftn
  integer          :: ios,iloc
  real             :: lat,lon,col,row,value
  integer          :: stn_col, stn_row
  character*50     :: stnname
  character*200    :: cline
  type(ESMF_Time)  :: ghcntime
  integer          :: status
  integer          :: t
  integer          :: year, month, day, hour, minute
  real             :: maxtemp, mintemp, prcp, snowf, snod
  logical          :: file_exists

  inquire(file=trim(filename),exist=file_exists)

  if(file_exists) then 
     write(LVT_logunit,*) '[INFO] Reading GHCN file ',trim(filename)
     ftn = LVT_getNextUnitNumber()
     open(ftn,file=trim(filename),form='formatted')

!skip the first line
     read(ftn, *)

     ios = 0 
     do while(ios.eq.0)
        read(ftn,fmt='(a)',iostat=ios) cline
        if(ios.ne.0) exit
        iloc = index(cline,",")
        read(cline(1:iloc-1),*) stnname
        cline = cline(iloc+1:len(cline))

        iloc = index(cline,",")
        read(cline(1:iloc-1),*) year
        cline = cline(iloc+1:len(cline))

        iloc = index(cline,",")
        read(cline(1:iloc-1),*) month
        cline = cline(iloc+1:len(cline))

        iloc = index(cline,",")
        read(cline(1:iloc-1),*) day
        cline = cline(iloc+1:len(cline))

        iloc = index(cline,":")
        read(cline(1:iloc-1),*) hour
        cline = cline(iloc+1:len(cline))

        iloc = index(cline,",")
        read(cline(1:iloc-1),*) minute
        cline = cline(iloc+1:len(cline))

        iloc = index(cline,",")
        read(cline(1:iloc-1),*) mintemp
        cline = cline(iloc+1:len(cline))

        iloc = index(cline,",")
        read(cline(1:iloc-1),*) maxtemp
        cline = cline(iloc+1:len(cline))
        
        iloc = index(cline,",")
        read(cline(1:iloc-1),*) prcp
        cline = cline(iloc+1:len(cline))

        iloc = index(cline,",")
        read(cline(1:iloc-1),*) snowf
        cline = cline(iloc+1:len(cline))

        read(cline,*) snod
        
        if(month.ne.99.and.day.ne.99.and.hour.ne.99.and.minute.ne.99.and.&
             snod.gt.0) then 
           
           call getGHCNstnIndex(source, stnname, i)
           
           if(i.gt.0) then 
              call ESMF_TimeSet(ghcntime, yy=year, mm=month,&
                   dd=day, calendar=LVT_calendar, rc=status)
              call LVT_verify(status, 'ESMF_TimeSet failed in read_ghcndata')
              
              t = nint((ghcntime - ghcnobs(source)%startTime)/&
                   ghcnobs(source)%timestep) + 1

              ghcnobs(source)%snod(i,t) = snod/1000.0
              ghcnobs(source)%prcp(i,t) = prcp

           endif
        endif
     enddo

     call LVT_releaseUnitNumber(ftn)
  endif

end subroutine read_ghcndata

subroutine getGHCNstnIndex(source, stnname,stnIndex)

  use GHCN_obsMod,    only : ghcnobs
  
  implicit none
  
  integer              :: source 
  character(len=*)     :: stnname
  integer              :: stnIndex
  integer              :: i

  stnIndex = -1

  do i=1,ghcnobs(source)%nstns  
     if(stnname.eq.ghcnobs(source)%stnid(i)) then 
        stnIndex = i
        exit;
     endif
  enddo
end subroutine getGHCNstnIndex
!BOP
! 
! !ROUTINE: create_ghcnsnwd_filename
! \label{create_ghcnsnwd_filename}
!
! !INTERFACE: 
subroutine create_ghcnsnwd_filename(odir, yr, ghcnname)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: ghcnname
!
! !DESCRIPTION: 
! 
! This routine creates a filename for the GHCN station
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir] location of the GHCN files
!   \item[yr]   year of the data
!   \item[stnid] station id
!   \item[ghcnname] the name of the GHCN file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4             :: fyr
  
  write(fyr, '(i4.4)' ) yr

  ghcnname = trim(odir)//'/'//trim(fyr)//'/ghcn-proc_'//&
       trim(fyr)//'.txt'
  
end subroutine create_ghcnsnwd_filename


