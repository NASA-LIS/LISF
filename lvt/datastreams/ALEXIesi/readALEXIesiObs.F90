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
! !ROUTINE: readALEXIesiObs
! \label{readALEXIesiObs}
!
! !INTERFACE: 
subroutine readALEXIesiObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_logMod
  use LVT_histDataMod
  use ALEXIesi_obsMod

  implicit none
  integer,   intent(in)   :: source
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  
! This subroutine provides the data reader for the 
! ESI estimates from the ALEXI data. The daily ESI
! files are read and are spatially interpolated/upscaled
! to match the analysis grid/domain. 
! 
! Note that currently, ESI data over CONUS only is 
! supported. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  08 Feb 2020: Sujay Kumar, Initial Specification
! 
!EOP
  real                   :: currTime
  logical                :: alarmCheck 
  logical                :: file_exists
  character*100          :: esi_filename
  integer                :: c,r,ios,ftn
  logical*1              :: li(ALEXIesiobs(source)%nc*ALEXIesiobs(source)%nr)
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: esit(ALEXIesiobs(source)%nc,ALEXIesiobs(source)%nr)
  real                   :: esi_in(ALEXIesiobs(source)%nc*ALEXIesiobs(source)%nr)
  real                   :: esi_out(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: varfield(LVT_rc%lnc, LVT_rc%lnr)

  varfield = LVT_rc%udef
  esi_out = LVT_rc%udef

  currTime = float(LVT_rc%dhr(source))*3600+ &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmCheck = (mod(currtime,86400.0).eq.0)

  if(ALEXIesiobs(source)%startFlag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     ALEXIesiobs(source)%startFlag = .false. 
     
     call create_ALEXIesi_filename(ALEXIesiobs(source)%odir, &
          LVT_rc%dyr(source), &
          LVT_rc%ddoy(source), esi_filename)

     inquire(file=esi_filename, exist=file_exists) 
     if(file_exists) then 
        ftn = LVT_getNextUnitNumber()
        write(LVT_logunit,*) '[INFO] Reading ',trim(esi_filename)
        open(ftn, file=esi_filename, form='unformatted',access='direct',&
             convert='little_endian',&
             recl=ALEXIesiobs(source)%nc*ALEXIesiobs(source)%nr*4)
        read(ftn,rec=1) esit
        call LVT_releaseUnitNumber(ftn)
        
        li = .false. 
        do r=1,ALEXIesiobs(source)%nr
           do c=1,ALEXIesiobs(source)%nc
              esi_in(c+(r-1)*ALEXIesiobs(source)%nc) =&
                   esit(c,ALEXIesiobs(source)%nr-r+1)
              if(esit(c,ALEXIesiobs(source)%nr-r+1).ne.LVT_rc%udef.and.&
                   .not.isNaN(esit(c,ALEXIesiobs(source)%nr-r+1))) then 
                 li(c+(r-1)*ALEXIesiobs(source)%nc) = .true. 
              else
                 esi_in(c+(r-1)*ALEXIesiobs(source)%nc) = LVT_rc%udef
              endif
           enddo
        enddo
        if(LVT_isAtAfinerResolution(ALEXIesiobs(source)%datares)) then
           call neighbor_interp(LVT_rc%gridDesc,li,esi_in,&
                lo, esi_out, ALEXIesiobs(source)%nc*ALEXIesiobs(source)%nr, &
                LVT_rc%lnc*LVT_rc%lnr,&
                ALEXIesiobs(source)%rlat, ALEXIesiobs(source)%rlon, &
                ALEXIesiobs(source)%n11,LVT_rc%udef, ios)
        else
           call upscaleByAveraging(ALEXIesiobs(source)%nc*ALEXIesiobs(source)%nr,&
                LVT_rc%lnc*LVT_rc%lnr, &
                LVT_rc%udef, ALEXIesiobs(source)%n11,&
                li,esi_in, lo, esi_out)
           
        endif
     else
        esi_out = LVT_rc%udef
     endif
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(esi_out(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then
              varfield(c,r) = esi_out(c+(r-1)*LVT_rc%lnc)
           else
              varfield(c,r) = LVT_rc%udef
           endif
        enddo
     enddo
  endif
  call LVT_logSingleDataStreamVar(LVT_MOC_ESI,source,varfield,&
       vlevel=1,units="-")

  
end subroutine readALEXIesiObs

!BOP
! 
! !ROUTINE: create_ALEXIesi_filename
! \label{create_ALEXIesi_filename}
!
! !INTERFACE: 
subroutine create_ALEXIesi_filename(odir,yr,doy,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for ALEXI ESI data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      ALEXI ESI  base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the ALEXI ESI data file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: doy
  character(len=*)             :: filename
!EOP

  character*4             :: fyr
  character*3             :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy,fmt='(i3.3)') doy

  filename = trim(odir)//'/'//trim(fyr)//'/'//&
       'FPPM_FLT_'//trim(fyr)//trim(fdoy)//'.dat'

end subroutine create_ALEXIesi_filename

