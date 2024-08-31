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
! !ROUTINE: get_PALSmetdata
! \label{get_PALSmetdata}
!
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
! 10 Sept 2018: Jules Kouatchou
!    - Fixed a bug in the routine read_PALSvar where
!      (tindex2-tindex1)+1 values were read in instead of
!      only  tindex2-tindex1 values.
!    - (Lines 109-112) Use the current date/time (instead of January 1st of the current year)
!      for the first calculation of time1. With the new setting, an experiment
!      can begin at any date/time after the reference Date/Time of the PALS data.
! 
! !INTERFACE:
subroutine get_PALSmetdata(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_metforcingMod,  only : LIS_forc
  use LIS_logMod
  use PALSmetdata_forcingMod,  only : PALSmetdata_struc
  use LIS_constantsMod,        only : LIS_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and assigns 1/2 hourly, PALS met data as forcing. 
!  At the beginning of a simulation or at the change of years, 
!  the code reads the data (for that year and keeps it in memory). 
!  past data (nearest 1./2 hourly interval), and the nearest future data.
!  
!  NOTES: No spatial interpolation is applied. The LIS domain is expected
!  to be setup over the station of interest. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[PALSmetdatafile](\ref{PALSmetdatafile}) \newline
!    Generates the PALS forcing file name
!  \item[read\_PALSvar](\ref{read_PALSvar}) \newline
!    Reads a variable from PALS file, applies the QC and 
!    assigns it to the LIS data structures.
!  \end{description}
!
!EOP
  character(len=LIS_CONST_PATH_LEN) :: name
  logical                       :: file_exists
  integer                       :: ftn
  integer                       :: t
  integer                       :: tindex1, tindex2, tindex3, tindex_max
  integer                       :: timeId
  type(ESMF_Time)               :: time1, time2, currTime
  integer                       :: tairId,qairId,swdownId, lwdownId
  integer                       :: windId,psurfId,rainfId 
  integer                       :: status
  real,       allocatable           :: varfield(:)
  integer                       :: yr1,mo1,da1,hr1,mn1,ss1,ts1,doy1
  real                          :: gmt1
  
! At the start of the year, read the current year and the next year's data
! and keep in memory. 
! 
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  if(PALSmetdata_struc(n)%startFlag.or.&
       (PALSmetdata_struc(n)%yr.ne.LIS_rc%yr)) then 

     PALSmetdata_struc(n)%startFlag = .false. 
     PALSmetdata_struc(n)%yr = LIS_rc%yr

     call PALSmetdatafile(name,PALSmetdata_struc(n)%PALSmetdatadir,&     
          PALSmetdata_struc(n)%stn_name)

     inquire(file=name,exist=file_exists) 

     if(file_exists) then 
        write(LIS_logunit,*) 'Reading '//trim(name)
        call LIS_verify(nf90_open(name, mode=nf90_nowrite,&
             ncid=ftn), 'nf90_open failed in get_PALSmetdata')     
        
        call LIS_verify(nf90_inq_dimid(ftn,'time',timeId),&
             'nf90_inq_dimid failed for time')
        call LIS_verify(nf90_inquire_dimension(ftn,timeId,len=tindex_max),&
             'nf90_inquire_dimension failed for time')

        call ESMF_TimeSet(time1, yy=LIS_rc%yr, &
             mm=LIS_rc%mo,dd=LIS_rc%da,h=LIS_rc%hr,m=LIS_rc%mn,s=LIS_rc%ss,calendar=LIS_calendar, &
             !mm=1,dd=1,h=0,m=0,s=0,calendar=LIS_calendar, &
             rc=status)
        
        call ESMF_TimeSet(time2, yy=LIS_rc%yr, &
             mm=12,dd=31,h=24,m=0,s=0,calendar=LIS_calendar, &
             rc=status)
        call LIS_verify(status, 'error in ESMF_TimeSet in get_PALSmetdata')
        
        tindex1 = nint((time1-PALSmetdata_struc(n)%reftime)/&
             PALSmetdata_struc(n)%timestep) + 1 
        tindex2 = nint((time2-PALSmetdata_struc(n)%reftime)/&
             PALSmetdata_struc(n)%timestep) + 1 
        
        PALSmetdata_struc(n)%tindex1 = tindex1
        PALSmetdata_struc(n)%tindex2 = tindex2 

        if(tindex1.le.0) then 
           write(LIS_logunit,*) 'LIS start time is greater than the '
           write(LIS_logunit,*) 'start time of PALS data ..'
           call LIS_endrun()
        endif

        if(tindex1.gt.tindex_max.or.tindex2.gt.tindex_max) then 
           write(LIS_logunit,*) 'LIS time is exceeds the  '
           write(LIS_logunit,*) 'PALS time frame '
           call LIS_endrun()
        endif

!Tair        
        allocate(varfield(tindex2-tindex1))
        
        call read_PALSvar(ftn,"Tair", tindex1, tindex2,varfield,&
             .true.)
        PALSmetdata_struc(n)%tair1(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)
!Qair
        call read_PALSvar(ftn,"Qair", tindex1, tindex2,varfield,&
             .true.)
        PALSmetdata_struc(n)%qair1(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)
!SWdown
        call read_PALSvar(ftn,"SWdown", tindex1, tindex2,varfield,&
             .false.)
        PALSmetdata_struc(n)%swdown1(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)
!Lwdown
        call read_PALSvar(ftn,"LWdown", tindex1, tindex2,varfield,&
             .false.)
        PALSmetdata_struc(n)%lwdown1(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)
!Wind
        call read_PALSvar(ftn,"Wind", tindex1, tindex2,varfield,&
             .true.)
        PALSmetdata_struc(n)%wind1(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)

!Psurf
         call read_PALSvar(ftn,"PSurf", tindex1, tindex2,varfield,&
              .false.)
        PALSmetdata_struc(n)%psurf1(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)

!Rainf
        call read_PALSvar(ftn,"Rainf", tindex1, tindex2,varfield,&
             .false.)
        PALSmetdata_struc(n)%rainf1(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)
        deallocate(varfield)

        call ESMF_TimeSet(time1, yy=LIS_rc%yr+1, &
             mm=1,dd=1,h=0,m=0,s=0,calendar=LIS_calendar, &
             rc=status)
        
        call ESMF_TimeSet(time2, yy=LIS_rc%yr+1, &
             mm=12,dd=31,h=24,m=0,s=0,calendar=LIS_calendar, &
             rc=status)
        
        tindex1 = nint((time1-PALSmetdata_struc(n)%reftime)/&
             PALSmetdata_struc(n)%timestep) + 1 
        tindex2 = nint((time2-PALSmetdata_struc(n)%reftime)/&
             PALSmetdata_struc(n)%timestep) + 1 
        
        if(tindex1.le.0) then 
           write(LIS_logunit,*) 'The LIS start time is greater than the '
           write(LIS_logunit,*) 'start time of PALS data ..'
           call LIS_endrun()
        endif

        if(tindex1.gt.tindex_max.or.tindex2.gt.tindex_max) then 

           PALSmetdata_struc(n)%tair2 = LIS_rc%udef

        else

           allocate(varfield(tindex2-tindex1))
           call read_PALSvar(ftn,"Tair", tindex1, tindex2,varfield,&
                .true.)
           PALSmetdata_struc(n)%tair2(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)

           !Qair
           call read_PALSvar(ftn,"Qair", tindex1, tindex2,varfield,&
                .true.)
           PALSmetdata_struc(n)%qair2(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)   

           !SWdown
           call read_PALSvar(ftn,"SWdown", tindex1, tindex2,varfield,&
                .false.)
           PALSmetdata_struc(n)%swdown2(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)     

           !Lwdown           
           call read_PALSvar(ftn,"LWdown", tindex1, tindex2,varfield,&
                .false.)
           PALSmetdata_struc(n)%lwdown2(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)     
           !Wind
           call read_PALSvar(ftn,"Wind", tindex1, tindex2,varfield,&
                .true.)
           PALSmetdata_struc(n)%wind2(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)    
           !Psurf
           call read_PALSvar(ftn,"PSurf", tindex1, tindex2,varfield,&
                .false.)
           PALSmetdata_struc(n)%psurf2(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)   

           !Rainf
           call read_PALSvar(ftn,"Rainf", tindex1, tindex2,varfield,&
                .false.)
           PALSmetdata_struc(n)%rainf2(1:tindex2-tindex1) = varfield(1:tindex2-tindex1)
           deallocate(varfield)

        endif
        
        call LIS_verify(nf90_close(ftn))
     else
        write(LIS_logunit,*) 'Forcing file '//trim(name)//' does not exist'
        call LIS_endrun()
     endif
  endif
! At all intermediate times except the start of an year, index into the 
! data

  call ESMF_TimeSet(currTime, yy=LIS_rc%yr, mm=LIS_rc%mo, &
       dd=LIS_rc%da, h=LIS_rc%hr, m=LIS_rc%mn, s=LIS_rc%ss,&
       calendar = LIS_calendar, rc=status)
  tindex1 = (nint((currTime - PALSmetdata_struc(n)%reftime)/&
       PALSmetdata_struc(n)%timestep) + 1 ) - PALSmetdata_struc(n)%tindex1 + 1

  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=LIS_rc%ss
  ts1=0
  call LIS_tick(PALSmetdata_struc(n)%fcsttime1,doy1,gmt1,&
       yr1,mo1,da1,hr1,mn1,ss1,real(ts1))

  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=LIS_rc%ss
  ts1=1800
  call LIS_tick(PALSmetdata_struc(n)%fcsttime2,doy1,gmt1,&
       yr1,mo1,da1,hr1,mn1,ss1,real(ts1))

  PALSmetdata_struc(n)%metdata1(1,:) = PALSmetdata_struc(n)%tair1(tindex1)
  PALSmetdata_struc(n)%metdata1(2,:) = PALSmetdata_struc(n)%qair1(tindex1)
  PALSmetdata_struc(n)%metdata1(3,:) = PALSmetdata_struc(n)%swdown1(tindex1)
  PALSmetdata_struc(n)%metdata1(4,:) = PALSmetdata_struc(n)%lwdown1(tindex1)
  PALSmetdata_struc(n)%metdata1(5,:) = PALSmetdata_struc(n)%wind1(tindex1)
  PALSmetdata_struc(n)%metdata1(6,:) = 0.0
  PALSmetdata_struc(n)%metdata1(7,:) = PALSmetdata_struc(n)%psurf1(tindex1)
  PALSmetdata_struc(n)%metdata1(8,:) = PALSmetdata_struc(n)%rainf1(tindex1)
  
  if(tindex1.le.PALSmetdata_struc(n)%tindex2) then 
     PALSmetdata_struc(n)%metdata2(1,:) = PALSmetdata_struc(n)%tair1(tindex1)
     PALSmetdata_struc(n)%metdata2(2,:) = PALSmetdata_struc(n)%qair1(tindex1)
     PALSmetdata_struc(n)%metdata2(3,:) = PALSmetdata_struc(n)%swdown1(tindex1)
     PALSmetdata_struc(n)%metdata2(4,:) = PALSmetdata_struc(n)%lwdown1(tindex1)
     PALSmetdata_struc(n)%metdata2(5,:) = PALSmetdata_struc(n)%wind1(tindex1)
     PALSmetdata_struc(n)%metdata2(6,:) = 0.0
     PALSmetdata_struc(n)%metdata2(7,:) = PALSmetdata_struc(n)%psurf1(tindex1)
     PALSmetdata_struc(n)%metdata2(8,:) = PALSmetdata_struc(n)%rainf1(tindex1)
  else !2nd bookend should be from year 2, first data. 
     tindex3 = tindex1 -PALSmetdata_struc(n)%tindex2
     PALSmetdata_struc(n)%metdata2(1,:) = PALSmetdata_struc(n)%tair2(tindex3)
     PALSmetdata_struc(n)%metdata2(2,:) = PALSmetdata_struc(n)%qair2(tindex3)
     PALSmetdata_struc(n)%metdata2(3,:) = PALSmetdata_struc(n)%swdown2(tindex3)
     PALSmetdata_struc(n)%metdata2(4,:) = PALSmetdata_struc(n)%lwdown2(tindex3)
     PALSmetdata_struc(n)%metdata2(5,:) = PALSmetdata_struc(n)%wind2(tindex3)
     PALSmetdata_struc(n)%metdata2(6,:) = 0.0
     PALSmetdata_struc(n)%metdata2(7,:) = PALSmetdata_struc(n)%psurf2(tindex3)
     PALSmetdata_struc(n)%metdata2(8,:) = PALSmetdata_struc(n)%rainf2(tindex3)
  endif

#endif  
end subroutine get_PALSmetdata

!BOP
! !ROUTINE: PALSmetdatafile
! \label{PALSmetdatafile}
!
! !INTERFACE:
subroutine PALSmetdatafile(name,PALSmetdatadir,stn_name)
  
  implicit none
! !ARGUMENTS: 
  character(len=*), intent(out) :: name
  character(len=*), intent(in)  :: PALSmetdatadir
  character(len=*), intent(in)  :: stn_name

! !DESCRIPTION:
!  This subroutine puts together PALS met forcing file name 
! 
!  The arguments are:
!  \begin{description}
!  \item[name]
!   name of the PALS met forcing file
!  \item[PALSmetdatadir]
!    Name of the PALS data directory
!  \item[stn\_name]
!   name of the station
!  \end{description}
!
!EOP
  name = trim(PALSmetdatadir)//'/'//trim(stn_name)//'Fluxnet.1.4_met.nc'

end subroutine PALSmetdatafile

!BOP
! !ROUTINE: read_PALSvar
! \label{read_PALSvar}
!
! !INTERFACE: 
subroutine read_PALSvar(ftn, varname, tindex1, tindex2, varfield, vardim)
! !USES: 
  use LIS_coreMod
  use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! 
!
! !ARGUMENTS: 
  integer              :: ftn
  character(len=*)     :: varname
  integer              :: tindex1
  integer              :: tindex2
  real                 :: varfield(tindex2-tindex1)
  logical              :: vardim ! true means 3d else 2d
!
! !DESCRIPTION: 
!  This routine reads a variable from the station data file and 
!  applies the quality control filters. 
! 
!  The arguments are:
!  \begin{description}
!  \item[ftn]
!   netcdf file handle
!  \item[varname]
!   name of the variable to be retrieved
!  \item[tindex1]
!   starting time index of the data
!  \item[tindex2]
!   ending time index of the data
!  \item[varfield]
!   retrieved variable field
!  \item[vardim]
!   logical variable indicating if the variable to be 
!   read is 2d or 3d (.true. = 3d, .false. = 2d)
!  \end{description}
!EOP
  integer              :: varqcfield(tindex2-tindex1)
  integer              :: t
  integer              :: varId, varqcId

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 

  call LIS_verify(nf90_inq_varid(ftn,varname,varId),&
       'nf90_inq_varid failed for '//trim(varname)//' in get_PALSmetdata')
  call LIS_verify(nf90_inq_varid(ftn,trim(varname)//'_qc',varqcId),&
       'nf90_inq_varid failed for '//trim(varname)//'_qc in get_PALSmetdata')
  if(vardim) then 
     call LIS_verify(nf90_get_var(ftn,varId, &
          varfield,start=(/1,1,1,tindex1/), &
          count=(/1,1,1,tindex2-tindex1/)),&
          'nf90_get_var failed for '//trim(varname))
  else
     call LIS_verify(nf90_get_var(ftn,varId, &
          varfield,start=(/1,1,tindex1/), &
          count=(/1,1,tindex2-tindex1/)),&
          'nf90_get_var failed for '//trim(varname))
  endif
  call LIS_verify(nf90_get_var(ftn,varqcId, &
       varqcfield,start=(/1,1,tindex1/), &
       count=(/1,1,tindex2-tindex1/)),&
       'nf90_get_var failed for '//trim(varname)//'_qc')

!  do t=1,tindex2-tindex1
!     if(varqcfield(t).ne.1) then 
!        varfield(t) = LIS_rc%udef
!     endif
!  enddo
#endif

end subroutine read_PALSvar
