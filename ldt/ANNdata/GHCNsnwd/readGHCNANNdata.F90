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
! !ROUTINE: readGHCNANNdata
! \label{readGHCNANNdata}
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readGHCNANNdata(n,iomode,sindex,eindex)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_ANNMod
  use GHCN_ANNdataMod
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: iomode
  integer, intent(in) :: sindex
  integer, intent(in) :: eindex
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the synthetic
! soil moisture retrieval product. 
!
!EOP

  integer                 :: i,t,c,r,kk
  character(len=LDT_CONST_PATH_LEN) :: ghcnname
  type(ESMF_Time)         :: ghcntime1, ghcntime2
  real                    :: snowdepth(LDT_rc%lnc(n),LDT_rc%lnr(n))
  integer                 :: stn_col, stn_row
  real                    :: col, row
  integer                 :: status

  snowdepth  = -9999.0

  if(GHCNobs(n)%yr.ne.LDT_rc%yr) then 
     GHCNobs(n)%yr = LDT_rc%yr
     GHCNobs(n)%snod = LDT_rc%udef

     call ESMF_TimeSet(GHCNobs(n)%startTime,  yy=LDT_rc%yr, &
          mm = 1, &
          dd = 1, &
          h = 0, &
          m = 0, &
          calendar = LDT_calendar, &
          rc=status)
     call LDT_verify(status, 'error in setting ghcn time')
     
     call  create_ghcnsnwd_filename(GHCNobs(n)%odir, LDT_rc%yr, &
          ghcnname)
     
     call read_ghcndata(n, ghcnname)  

  endif

  call ESMF_TimeSet(ghcntime1, yy=LDT_rc%yr, &
       mm=LDT_rc%mo, dd=LDT_rc%da, h=LDT_rc%hr,&
       m=LDT_rc%mn, &
       s = LDT_rc%ss, calendar=LDT_calendar, rc=status)
  call LDT_verify(status, 'ghcntime1 set failed')

  t = nint((ghcntime1-GHCNobs(n)%startTime)/GHCNobs(n)%timestep) + 1

  do i=1,GHCNobs(n)%nstns
     call latlon_to_ij(LDT_domain(n)%ldtproj, &
          GHCNobs(n)%stnlat(i), GHCNobs(n)%stnlon(i),&
          col,row)
     stn_col = nint(col)
     stn_row = nint(row)
     
     if(stn_col.ge.1.and.stn_col.le.LDT_rc%lnc(n).and.&
          stn_row.ge.1.and.stn_row.le.LDT_rc%lnr(n).and.&
          GHCNobs(n)%snod(i,t).gt.0) then 
        snowdepth(stn_col,stn_row) = &
             GHCNobs(n)%snod(i,t)
     endif
     
  enddo
  if(LDT_rc%mo.eq.6.or.LDT_rc%mo.eq.7.or.LDT_rc%mo.eq.8) then 
     snowdepth = 0.0
  endif

  call LDT_logSingleANNdata(n,&
       snowdepth,  &
       pindex=sindex, &
       iomode = iomode, &
       name = "SnowDepth",&
       units="K")
  
end subroutine readGHCNANNdata

!BOP
! 
! !ROUTINE: read_ghcndata
! \label{read_ghcndata}
!
! !INTERFACE: 
subroutine read_ghcndata(n, filename)
! 
! !USES: 
  use ESMF
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use GHCN_ANNdataMod

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
  integer                 :: n 
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
     write(LDT_logunit,*) '[INFO] Reading GHCN file ',trim(filename)
     ftn = LDT_getNextUnitNumber()
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
           
           call getGHCNstnIndex(n, stnname, i)
           
           if(i.gt.0) then 
              call ESMF_TimeSet(ghcntime, yy=year, mm=month,&
                   dd=day, calendar=LDT_calendar, rc=status)
              call LDT_verify(status, 'ESMF_TimeSet failed in read_ghcndata')
              
              t = nint((ghcntime - ghcnobs(n)%startTime)/&
                   ghcnobs(n)%timestep) + 1

              ghcnobs(n)%snod(i,t) = snod/1000.0

           endif
        endif
     enddo

     call LDT_releaseUnitNumber(ftn)
  endif

end subroutine read_ghcndata

subroutine getGHCNstnIndex(n, stnname,stnIndex)

  use GHCN_ANNdataMod,    only : ghcnobs
  
  implicit none
  
  integer              :: n 
  character(len=*)     :: stnname
  integer              :: stnIndex
  integer              :: i

  stnIndex = -1

  do i=1,ghcnobs(n)%nstns  
     if(stnname.eq.ghcnobs(n)%stnid(i)) then 
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


