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
! !ROUTINE: readAquariusL2smObs
! \label{readAquariusL2smObs}
! 
! !REVISION HISTORY: 
!  21 July 2014: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readAquariusL2smObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod,      only : LDT_rc
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use AquariusL2sm_obsMod, only : AquariusL2smobs
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
!  
!
!EOP

  real                    :: timenow
  logical                 :: alarmCheck
  logical                 :: file_exists
  integer                 :: c,r,i,j
  character(len=LDT_CONST_PATH_LEN)           :: fname
  character(len=LDT_CONST_PATH_LEN)           :: aquarius_filename
  character*7             :: yyyyddd
  character*4             :: fyr
  character*2             :: fmo,fda
  character*200           :: list_files
  integer                 :: sind
  integer                 :: yr,doy,mo,da,hr,mn,ss
  integer                 :: ftn
  type(ESMF_Time)         :: currTime, prevTime
  type(ESMF_Time)         :: obsTime
  type(ESMF_TimeInterval) :: dayInt
  integer                 :: status
  integer                 :: ierr
  real                    :: smobs(LDT_rc%lnc(n),LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  call ESMF_TimeintervalSet(dayInt, d=1,rc=status)
  call LDT_verify(status, 'ESMF_TimeIntervalSet failed in readAquariusL2smObs')

  timenow = float(LDT_rc%hr)*3600 + 60*LDT_rc%mn + LDT_rc%ss
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  AquariusL2smobs(n)%smobs = LDT_rc%udef
  smobs= LDT_rc%udef
        
  if(AquariusL2smobs(n)%startmode.or.alarmCheck) then 
     
     AquariusL2smobs(n)%startmode = .false. 

! dump the list of files for the current date to a file (note that
! we assume a flat organization of the files under the Aquarius observation
! directory. 

     write(yyyyddd,'(i4.4,i3.3)') &
          LDT_rc%yr, &
          LDT_rc%doy

     write(fyr,'(i4.4)') & 
          LDT_rc%yr
     write(fmo,'(i2.2)') & 
          LDT_rc%mo
     write(fda,'(i2.2)') & 
          LDT_rc%da

     list_files = 'ls '//trim(AquariusL2smobs(n)%odir)//'/'//&
          trim(fyr)//'.'//trim(fmo)//'.'//trim(fda)//'/*'//&
          trim(yyyyddd)&
          //'*_V3.0 > Aquarius_filelist.dat'
     call system(trim(list_files))

     ftn = LDT_getNextUnitNumber()
     open(ftn,file="./Aquarius_filelist.dat",status='old',iostat=ierr)
     do while(ierr.eq.0) 
        read(ftn,'(a)',iostat=ierr) aquarius_filename
        if(ierr.ne.0) then 
           exit
        else
! check first if the observation time is within the LDT time window (daily)
           call ESMF_TimeSet(currTime, yy=LDT_rc%yr,&
                mm = LDT_rc%mo, dd=LDT_rc%da, h = LDT_rc%hr, &
                m = LDT_rc%mn, s = LDT_rc%ss, calendar = LDT_calendar,&
                rc=status)
           call LDT_verify(status, 'ESMF_TimeSet failed in readAquariusL2smObs')
           prevTime = currTime 
           currTime = currTime + dayInt
           
           sind = index(aquarius_filename, "Q")
           sind = sind+1
           read(aquarius_filename(sind:sind+3),'(i4)') yr
           read(aquarius_filename(sind+4:sind+6),'(i3)') doy
           read(aquarius_filename(sind+7:sind+8),'(i2)') hr
           read(aquarius_filename(sind+9:sind+10),'(i2)') mn
           read(aquarius_filename(sind+11:sind+12),'(i2)') ss
           
!           call LDT_doy2moda(yr, doy, mo, da)
           call LDT_doy2date(yr, doy, mo, da)

           call ESMF_TimeSet(obsTime, yy=yr,&
                mm = mo, dd=da, h = hr, &
                m = mn, s = ss, calendar = LDT_calendar,&
                rc=status)
           call LDT_verify(status, 'ESMF_TimeSet failed in readAquariusL2smObs')
 
           if((obsTime.gt.prevTime).and.(obsTime.le.currTime)) then 
              call read_AquariusL2_data(n, aquarius_filename, smobs)
           endif
        endif
     enddo

     call LDT_releaseUnitNumber(ftn)

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(smobs(c,r).ne.-9999.0) then 
              AquariusL2smobs(n)%smobs(c,r) = smobs(c,r)
           endif
        enddo
     enddo
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%soilmoist_obs,&
       AquariusL2smobs(n)%smobs,vlevel=1)

end subroutine readAquariusL2smObs


!BOP
! 
! !ROUTINE: read_AquariusL2_data
! \label(read_AquariusL2_data)
!
! !INTERFACE:
subroutine read_AquariusL2_data(n, fname, smobs)
! 
! !USES:  
#if(defined USE_HDF5) 
  use hdf5
#endif 

  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use map_utils
  use AquariusL2sm_obsMod, only : AquariusL2smobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname
  real                          :: smobs(LDT_rc%lnc(n),LDT_rc%lnr(n))


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the AquariusL2 file
!  \item[smobs]    soil moisture data processed to the LDT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
#if(defined USE_HDF5) 
!need to be read from "Number of Blocks"
!  integer, parameter               :: aq_nc=4083,aq_nr = 3
  integer                          :: aq_nc,aq_nr
  integer(hid_t)                   :: file_id
  integer(hid_t)                   :: attr_id
  integer                          :: t,i,c,r,c1,c2,r1,r2
  real                             :: col1,row1,col2,row2
  integer(hid_t)                   :: rad_sm_id
  integer(hid_t)                   :: lat_id
  integer(hid_t)                   :: lon_id
  integer(hid_t)                   :: flg_id
  integer                          :: status
  real                             :: dx
  integer(hsize_t), allocatable    :: dims(:)
  real, allocatable                :: rad_sm(:,:)
  real, allocatable                :: lat(:,:)
  real, allocatable                :: lon(:,:)
  integer, allocatable             :: flag(:,:)
  logical                          :: file_exists
  integer(hsize_t), dimension(1) :: adims 
  integer, dimension(1)::  attr_data 

  inquire(file=fname,exist=file_exists)
  
  if(file_exists) then 
     aq_nr = 3

     write(LDT_logunit,*) 'Reading '//trim(fname)

     call h5open_f(status)
     call LDT_verify(status, 'Error opening HDF fortran interface')
     
     call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status)
     call LDT_verify(status, 'Error opening Aquarius file ')

     call h5dopen_f(file_id,"Aquarius Data/rad_sm",rad_sm_id,status)
     call LDT_verify(status,'h5dopen_f failed for rad_sm')

     call h5dopen_f(file_id,"Navigation/beam_clat",lat_id,status)
     call LDT_verify(status,'h5dopen_f failed for beam_clat')

     call h5dopen_f(file_id,"Navigation/beam_clon",lon_id,status)
     call LDT_verify(status,'h5dopen_f failed for beam_clon')

     call h5dopen_f(file_id,"Aquarius Flags/radiometer_flags",flg_id,status)
     call LDT_verify(status,'h5dopen_f failed for radiometer_flags')

     call h5aopen_f(file_id, "Number of Blocks", attr_id,status)
     call LDT_verify(status,'h5aopen_f failed for Number of Blocks')
     
     adims(1) = 1
     call h5aread_f(attr_id, H5T_NATIVE_INTEGER,attr_data, &
          adims,status)
     call LDT_verify(status,'h5aread_f failed for Number of Blocks')

     aq_nc = attr_data(1)

     allocate(rad_sm(aq_nc,aq_nr))
     allocate(lat(aq_nc,aq_nr))
     allocate(lon(aq_nc,aq_nr))
     allocate(flag(aq_nc,aq_nr))
     allocate(dims(2))
     dims(1) = aq_nc
     dims(2) = aq_nr

     call h5dread_f(rad_sm_id, H5T_NATIVE_REAL,rad_sm,&
          dims,status)
     call LDT_verify(status, 'h5dread_f failed for rad_sm')
     
     call h5dread_f(lat_id, H5T_NATIVE_REAL,lat,&
          dims,status)
     call LDT_verify(status, 'h5dread_f failed for lat')

     call h5dread_f(lon_id, H5T_NATIVE_REAL,lon,&
          dims,status)
     call LDT_verify(status, 'h5dread_f failed for lon')

     call h5dread_f(flg_id, H5T_NATIVE_INTEGER,flag,&
          dims,status)
     call LDT_verify(status, 'h5dread_f failed for flag')
     
!     dx = 0.25 !assuming a conversative 50km swath 
     dx = 0.10
     do i=1,aq_nr
        do t=1,aq_nc
           if(rad_sm(t,i).gt.0) then 
              if(.not.(btest(flag(t,i),0)).and.& !soil moisture retrieval performed
                   .not.(btest(flag(t,i),1)).and.& !Tb within acceptable range
                   .not.(btest(flag(t,i),3)).and.& !RFI acceptable
                   .not.(btest(flag(t,i),5)).and.& !ground not frozen
                   .not.(btest(flag(t,i),6)).and.& !SWE < 10 kg/m2
                   .not.(btest(flag(t,i),7)).and.& !ice fraction < 0.1
                   .not.(btest(flag(t,i),9)).and.& !veg water content < 5kg/m2
                   .not.(btest(flag(t,i),12))) then !not water
                 
                 
                 call latlon_to_ij(LDT_domain(n)%ldtproj,lat(t,i)-dx,lon(t,i)-dx,&
                      col1,row1)
                 call latlon_to_ij(LDT_domain(n)%ldtproj,lat(t,i)+dx,lon(t,i)+dx,&
                      col2,row2)
                 c1 = nint(col1)
                 r1 = nint(row1)
                 c2 = nint(col2)
                 r2 = nint(row2)
                 do r=r1,r2
                    do c=c1,c2
                       if((c.ge.1.and.c.le.LDT_rc%lnc(n)).and.&
                            (r.ge.1.and.r.le.LDT_rc%lnr(n))) then 
                          smobs(c,r) = rad_sm(t,i)
                       endif
                    enddo
                 enddo
              endif
           endif
        enddo
     enddo

     call h5dclose_f(rad_sm_id,status)
     call LDT_verify(status,'Error in H5DCLOSE call')
     
     call h5fclose_f(file_id,status)
     call LDT_verify(status, 'Error in h5fclose_f')
     
     call h5close_f(status)
     call LDT_verify(status, 'Error in h5close_f')
     deallocate(dims)
     
     deallocate(rad_sm)
     deallocate(lat)
     deallocate(lon)
     deallocate(flag)
  endif
#endif

end subroutine read_AquariusL2_data

