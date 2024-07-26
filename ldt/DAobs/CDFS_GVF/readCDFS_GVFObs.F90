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
! !ROUTINE: readCDFS_GVFObs
! \label{readCDFS_GVFObs}
! 
! !REVISION HISTORY: 
!  4 Mar 2022: Yonghwan Kwon, Initial Specification
! 
! !INTERFACE: 
subroutine readCDFS_GVFObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_DAobsDataMod
  use CDFSGVFobsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for 
! CDFS Green Vegetation Fraction (GVF) data
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r
  character*100     :: fname
  real              :: gvfobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  CDFSgvfobs(n)%gvfobs = LDT_rc%udef
  gvfobs= LDT_rc%udef

  call create_CDFSgvf_filename(CDFSgvfobs(n)%odir, &
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da,&
       fname)

  inquire(file=trim(fname),exist=file_exists)
  if(file_exists) then

     write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
     call read_CDFS_GVF_data(n, fname, gvfobs)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(gvfobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then
              CDFSgvfobs(n)%gvfobs(c,r) = gvfobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%gvf_obs,&
       CDFSgvfobs(n)%gvfobs,vlevel=1)

end subroutine readCDFS_GVFObs

!BOP
! 
! !ROUTINE: read_CDFS_GVF_data
! \label{read_CDFS_GVF_data}
!
! !INTERFACE:
subroutine read_CDFS_GVF_data(n, fname, gvfobs_ip)
! 
! !USES:

  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use CDFSGVFobsMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n
  character (len=*)             :: fname
  real                          :: gvfobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine reads the CDFS GVF file 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the CDFS GVF file
!  \item[gvtobs\_ip]   GVF data processed to the LIS domain
! \end{description}
!
!
!EOP

! !USES:   
  integer,  parameter     :: nc=7200, nr=3600
  real*4                  :: gvf_raw(CDFSgvfobs(n)%nc,CDFSgvfobs(n)%nr)
  real                    :: gvf_in(CDFSgvfobs(n)%nc*CDFSgvfobs(n)%nr)
  logical*1               :: gvf_data_b(CDFSgvfobs(n)%nc*CDFSgvfobs(n)%nr)
  logical*1               :: gvfobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer                 :: gvfid
  integer                 :: ios, nid
  integer                 :: c,r
  integer                 :: ftn1

  ftn1 = LDT_getNextUnitNumber()
  open(unit=ftn1,file=fname,form='unformatted',access='direct',convert='little_endian',recl=4*nc*nr,status='old')
  read(ftn1, rec=1) gvf_raw
  close(1)
  call LDT_releaseUnitNumber(ftn1)

  do r=1,CDFSgvfobs(n)%nr
     do c=1,CDFSgvfobs(n)%nc
        if (gvf_raw(c,r)>=0.and.&
           gvf_raw(c,r)<=100) then
           gvf_in(c+(r-1)*CDFSgvfobs(n)%nc) = gvf_raw(c,r)
           gvf_data_b(c+(r-1)*CDFSgvfobs(n)%nc) = .true.
        else
           gvf_in(c+(r-1)*CDFSgvfobs(n)%nc) = LDT_rc%udef
           gvf_data_b(c+(r-1)*CDFSgvfobs(n)%nc) = .false.
        endif 
     enddo
  enddo 

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  if(LDT_isLDTatAfinerResolution(n,0.05)) then
     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          gvf_data_b, gvf_in, gvfobs_b_ip, gvfobs_ip, &
          CDFSgvfobs(n)%nc*CDFSgvfobs(n)%nr, &
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          CDFSgvfobs(n)%w11,CDFSgvfobs(n)%w12,&
          CDFSgvfobs(n)%w21,CDFSgvfobs(n)%w22,&
          CDFSgvfobs(n)%n11,CDFSgvfobs(n)%n12,&
          CDFSgvfobs(n)%n21,CDFSgvfobs(n)%n22,LDT_rc%udef,ios)
  else
     call upscaleByAveraging(CDFSgvfobs(n)%nc*CDFSgvfobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_rc%udef, CDFSgvfobs(n)%n11,&
          gvf_data_b,gvf_in, gvfobs_b_ip, gvfobs_ip)
  endif

end subroutine read_CDFS_GVF_data


!BOP
! !ROUTINE: create_CDFSgvf_filename
! \label{create_CDFSgvf_filename}
! 
! !INTERFACE:
subroutine create_CDFSgvf_filename(ndir,yr,mo,da,filename)
! !USES:

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the CDFS GVF filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the CDFS GVF data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated CDFS GVF filename
! \end{description}
!EOP

  character*4             :: yyyy
  character*2             :: mm,dd

  write(unit=yyyy, fmt='(i4.4)') yr
  write(unit=mm, fmt='(i2.2)') mo
  write(unit=dd, fmt='(i2.2)') da
  
  filename = trim(ndir)//'/'//trim(yyyy)//'/green.'//&
             trim(yyyy)//trim(mm)//trim(dd)//'.1gd4r'

end subroutine create_CDFSgvf_filename





