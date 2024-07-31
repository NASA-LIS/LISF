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
! !ROUTINE: readASCATTUWsmObs
! \label{readASCATTUWsmObs}
! 
! !REVISION HISTORY: 
!  8 May 2013: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readASCATTUWsmObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod,      only : LDT_rc, LDT_domain
  use LDT_timeMgrMod,   only : LDT_get_julss
  use LDT_logMod,       only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use ASCATTUWsm_obsMod, only : ASCATTUWsmobs
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the ASCATTUW
! soil moisture retrieval product. The data has many layers 
! and the reader provides the options to select the 
! desired layer(s).
! 

!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------

  timenow = float(LDT_rc%hr)*3600 + 60*LDT_rc%mn + LDT_rc%ss
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  ASCATTUWsmobs(n)%smobs = LDT_rc%udef
  smobs= LDT_rc%udef
        
  if(ASCATTUWsmobs(n)%startmode.or.alarmCheck) then 
     
     ASCATTUWsmobs(n)%startmode = .false. 

     call read_ASCATTUW_data(n, ASCATTUWsmobs(n)%odir,&
          LDT_rc%yr, LDT_rc%mo, LDT_rc%da, smobs)

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(smobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
              ASCATTUWsmobs(n)%smobs(c,r) = smobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  else
     ASCATTUWsmobs(n)%smobs = LDT_rc%udef
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%soilmoist_obs,&
       ASCATTUWsmobs(n)%smobs,vlevel=1)

end subroutine readASCATTUWsmObs


!BOP
! 
! !ROUTINE: read_ASCATTUW_data
! \label(read_ASCATTUW_data)
!
! !INTERFACE:
subroutine read_ASCATTUW_data(n, odir, yr,mo,da,sm_data)
! 
! !USES:   
#if (defined USE_GRIBAPI)
  use grib_api
#endif
  use LDT_coreMod,  only : LDT_rc, LDT_domain
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use map_utils,    only : latlon_to_ij
  use ASCATTUWsm_obsMod, only : ASCATTUWsmobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: odir
  integer                       :: yr,mo,da
  real                          :: sm_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))

! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the ASCATTUW grib2 file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  the estimated error is above a predefined threshold (the recommeded
!  value is 5%). 

!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the ASCATTUW AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LDT domain
! \end{description}
!
!
!EOP
  real, parameter        :: err_threshold = 0.03
  character*200          :: ls_comm, cmd2
  character(len=LDT_CONST_PATH_LEN)          :: fname
  integer                :: ftn1, ftn2
  integer                :: fsize,n_data
  integer                :: i,c,r,t,k
  integer                :: iret
  real                   :: col,row
  character (len=4)      :: fyr
  character (len=2)      :: fmo,fda
  real,    allocatable       :: lat(:)
  real,    allocatable       :: lon(:)
  real,    allocatable       :: sds(:)
  real,    allocatable       :: err(:)
  logical*1              :: sm_data_b(ASCATTUWsmobs(n)%nc*ASCATTUWsmobs(n)%nr)
  logical*1              :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                   :: sm_ascat(ASCATTUWsmobs(n)%nc*ASCATTUWsmobs(n)%nr)

#if (defined USE_GRIBAPI)
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  
  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'.'//trim(fmo)//'/'//&
       'SDS_'//trim(fyr)//trim(fmo)//trim(fda)//&
       '*bin 2>&1 2>/dev/null > ascat_file'
  cmd2 = 'wc -w ascat_file > ascat_file.wc'

  call system(ls_comm)
  call system(cmd2)
  
  ftn1 = LDT_getNextUnitNumber()
  open(ftn1,file='ascat_file.wc',form='formatted',action='read')
  read(ftn1,*) fsize
  call LDT_releaseUnitNumber(ftn1)
  
  sm_ascat = LDT_rc%udef

  ftn1 = LDT_getNextUnitNumber()     
  open(ftn1,file='ascat_file',form='formatted',action='read')  

  do k=1,fsize
     read(ftn1,'(a)') fname    
     
     ftn2 = LDT_getNextUnitNumber()
     write(LDT_logunit,*) 'Reading ',trim(fname)
     open(ftn2,file=trim(fname),form='unformatted')
     read(ftn2) n_data
     
     if(n_data.gt.0) then 
        allocate(lat(n_data))
        allocate(lon(n_data))
        allocate(sds(n_data))
        allocate(err(n_data))
        
        read(ftn2) lon
        read(ftn2) lat
        read(ftn2) sds
        read(ftn2) err
        
        do i=1,n_data
           
           call latlon_to_ij(ASCATTUWsmobs(n)%ascattuwproj, &
                lat(i),lon(i),col,row)
           
           c = nint(col)
           r = nint(row)
           if(err(i).le.err_threshold) then 
              sm_ascat(c+(r-1)*ASCATTUWsmobs(n)%nc) = sds(i)
           endif
        enddo

        deallocate(lat)
        deallocate(lon)
        deallocate(sds)
        deallocate(err)
     endif
     call LDT_releaseUnitNumber(ftn2)
  enddo

  call LDT_releaseUnitNumber(ftn1)
  
  sm_data_b = .false. 
  do t=1,ASCATTUWsmobs(n)%nc*ASCATTUWsmobs(n)%nr
     if(sm_ascat(t).ne.LDT_rc%udef) then 
        sm_data_b(t) = .true. 
     endif
  enddo
!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp( LDT_rc%gridDesc(n,:),&
       sm_data_b, sm_ascat, smobs_b_ip, sm_data, &
       ASCATTUWsmobs(n)%nc*ASCATTUWsmobs(n)%nr, &
       LDT_rc%lnc(n)*LDT_rc%lnr(n), &
       LDT_domain(n)%lat, LDT_domain(n)%lon,&
       ASCATTUWsmobs(n)%n11, LDT_rc%udef, iret)

!  open(100,file='test_ip.bin',form='unformatted')
!  write(100) sm_data
!  close(100)
!  stop
#endif
end subroutine read_ASCATTUW_data

