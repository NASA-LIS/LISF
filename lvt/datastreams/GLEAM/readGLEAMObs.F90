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
! !ROUTINE: readGLEAMObs
! \label{readGLEAMObs}
!
! !INTERFACE: 
subroutine readGLEAMObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_logMod
  use LVT_histDataMod
  use GLEAM_obsMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
  integer,   intent(in)   :: source
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The GLEAM output is available at daily intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
!
!   Currently the implementation only supports the processing of
!   GLEAM evaporation data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  14 Feb 2017: Sujay Kumar, Initial Specification
! 
!EOP
  real                   :: currTime
  logical                :: alarmCheck 
  logical                :: file_exists
  character*100          :: lh_filename
  integer                :: c,r,ios,ftn, lhid
  logical*1              :: li(GLEAMobs(source)%nc*GLEAMobs(source)%nr)
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: lht(GLEAMobs(source)%nc*GLEAMobs(source)%nr)
  real                   :: lh(GLEAMobs(source)%nr, GLEAMobs(source)%nc)
  real                   :: qle(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: varfield(LVT_rc%lnc, LVT_rc%lnr)

  currTime = float(LVT_rc%dhr(source))*3600+ &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmCheck = (mod(currtime,86400.0).eq.0)

  if(GLEAMobs(source)%startFlag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     LVT_rc%resetFlag(source) = .false. 

     GLEAMobs(source)%startFlag = .false. 
     
     call create_GLEAM_lh_filename(GLEAMobs(source)%odir, &
          GLEAMobs(source)%version, &
          LVT_rc%dyr(source), &
          lh_filename)

     inquire(file=lh_filename, exist=file_exists) 

     if(file_exists) then 
        
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

        write(LVT_logunit,*) '[INFO] reading ',trim(lh_filename)

        ios = nf90_open(path=trim(lh_filename),mode=NF90_NOWRITE,ncid=ftn)
        call LVT_verify(ios, 'Error opening file, '//trim(lh_filename))

        ios = nf90_inq_varid(ftn,'E',lhid)
        call LVT_verify(ios,'readGLEAMObs: Error in nf90_inq_varid: E')
        
        ios = nf90_get_var(ftn,lhid,lh, start=(/1,1,LVT_rc%ddoy(source)/), &
             count=(/GLEAMobs(source)%nr, GLEAMobs(source)%nc,1/))
        call LVT_verify(ios,'readGLEAMObs: Error in nf90_get_var: E')
        
        ios = nf90_close(ftn)
        call LVT_verify(ios,'readGLEAMObs: Error in nf90_close')

        do r=1,GLEAMobs(source)%nr
           do c=1,GLEAMobs(source)%nc
              lht(c+(r-1)*GLEAMobs(source)%nc) = lh(GLEAMobs(source)%nr-r+1,c)
           enddo
        enddo
#endif
        li = .false. 
        do r=1,GLEAMobs(source)%nr
           do c=1,GLEAMobs(source)%nc
              if(lht(c+(r-1)*GLEAMobs(source)%nc).ne.-999.0) then 
                 li(c+(r-1)*GLEAMobs(source)%nc) = .true. 
              else
                 lht(c+(r-1)*GLEAMobs(source)%nc) = LVT_rc%udef
              endif
           enddo
        enddo
        call neighbor_interp(LVT_rc%gridDesc,li,lht,&
             lo, qle, GLEAMobs(source)%nc*GLEAMobs(source)%nr, &
             LVT_rc%lnc*LVT_rc%lnr,&
             GLEAMobs(source)%rlat, GLEAMobs(source)%rlon, &
             GLEAMobs(source)%n11,LVT_rc%udef, ios)
     else
        qle = LVT_rc%udef
     endif
     
  else
     write(LVT_logunit,*)'[WARN] GLEAM TotalEvap file missing: ',trim(lh_filename)
     write(LVT_logunit,*)'[WARN]  OR incorrectly entered.'
     qle = LVT_rc%udef
  endif

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(qle(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then
           varfield(c,r) = qle(c+(r-1)*LVT_rc%lnc)*2454000.0/86400.0 ! mm/day to W/m2
        else
           varfield(c,r) = LVT_rc%udef
        endif
     enddo
  enddo
  
!  print*, 'here'
!  open(100,file='test.bin',form='unformatted')
!  write(100) varfield
!  close(100)
!  stop

  call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,varfield,vlevel=1,units="W/m2")

  
end subroutine readGLEAMObs

!BOP
! 
! !ROUTINE: create_GLEAM_lh_filename
! \label{create_GLEAM_lh_filename}
!
! !INTERFACE: 
subroutine create_GLEAM_lh_filename(odir,version,yr,filename)
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
! This routine creates a timestamped filename for GLEAM_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the GLEAM_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  character(len=*)             :: version
  integer                      :: yr
  integer                      :: res
  integer                      :: doy
  character(len=*)             :: filename
!EOP

  character*4             :: fyr
  
  write(unit=fyr, fmt='(i4.4)') yr

  ! Yeosang Yoon
  filename = trim(odir)//'/v'//trim(version)//"/"//trim(fyr)//'/'//&
       'E_'//trim(fyr)//'_GLEAM_v'//trim(version)//'.nc'
!  filename = trim(odir)//'/'//trim(fyr)//'/'//&
!       'E_'//trim(fyr)//'_GLEAM_v'//trim(version)//'.nc'
  
end subroutine create_GLEAM_lh_filename

