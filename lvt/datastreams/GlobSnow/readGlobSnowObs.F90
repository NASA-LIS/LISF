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
! !ROUTINE: readGlobSnowObs
! \label{readGlobSnowObs}
!
! !INTERFACE: 
subroutine readGlobSnowObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod,      only : LVT_logunit, LVT_verify
  use LVT_histDataMod
  use GlobSnow_obsMod, only : GlobSnowObs

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)  :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!   This subroutine provides the data reader for the GlobSnow daily 
!   data in NETCDF format (L3A daily SWE). The routine reads the
!   data from the NETCDF file and spatially interpolates it to the
!   LIS model grid and resolution using the neighbor interpolation. 
!
!  NOTES: 
!   The GlobSnow output is available at daily intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  10 Dec 2010: Sujay Kumar, Initial Specification
! 
!EOP

  character*100          :: filename
  logical                :: file_exists
  integer                :: nid, ios
  integer                :: sweid, rowId, colId
  integer                :: nc,nr
  real,  allocatable     :: swe(:,:)
  real,  allocatable     :: swe1d(:)
  logical*1, allocatable :: li(:)
  integer                :: c,r
  real                   :: swe1d_ip(LVT_rc%lnc*LVT_rc%lnr)
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: varfield(LVT_rc%lnc, LVT_rc%lnr)

  swe1d_ip = LVT_rc%udef

  call create_globsnow_filename(GlobSnowObs(source)%odir, &
       LVT_rc%dyr(source), LVT_rc%dmo(source), &
       LVT_rc%dda(source), filename)
  inquire(file=trim(filename),exist=file_exists) 
  
  if(file_exists) then 
     write(LVT_logunit,*) '[INFO] Reading GlobSnow file ',trim(filename)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
     ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios, 'Error opening file '//trim(filename))
     
     ios = nf90_inq_varid(nid,'SWE',sweid)
     call LVT_verify(ios, 'Error nf90_inq_varid: SWE')
! dimensions
     ios = nf90_inq_dimid(nid,'rows',rowId)
     call LVT_verify(ios, 'Error nf90_inq_dimid: rows')

     ios = nf90_inquire_dimension(nid,rowId, len=nr)
     call LVT_verify(ios, 'Error nf90_inquire_dimension: row')

     ios = nf90_inq_dimid(nid,'cols',colId)
     call LVT_verify(ios, 'Error nf90_inq_dimid: cols')

     ios = nf90_inquire_dimension(nid,colId, len=nc)
     call LVT_verify(ios, 'Error nf90_inquire_dimension: col')

     allocate(swe(nc,nr))
     allocate(swe1d(nc*nr))
     allocate(li(nc*nr))
     li = .true. 
     
!values
     ios = nf90_get_var(nid,sweid, swe)
     call LVT_verify(ios, 'Error nf90_get_var: swe')

     ios = nf90_close(nid)
     call LVT_verify(ios, 'Error in nf90_close')

     do r=1,nr
        do c=1,nc
           if(isNaN(swe(c,r))) then 
              swe1d(c+(r-1)*nc) = LVT_rc%udef
              li(c+(r-1)*nc) = .false. 
           elseif(swe(c,r).lt.0) then 
              swe1d(c+(r-1)*nc) = LVT_rc%udef
              li(c+(r-1)*nc) = .false. 
           else
              swe1d(c+(r-1)*nc) = swe(c,r)/1000.0
              li(c+(r-1)*nc) = .true. 
           endif
        enddo
     enddo

     call neighbor_interp(LVT_rc%gridDesc,li,swe1d,&
          lo, swe1d_ip, nc*nr, LVT_rc%lnc*LVT_rc%lnr,&
          GlobSnowObs(source)%rlat, GlobSnowObs(source)%rlon, &
          GlobSnowObs(source)%n11,LVT_rc%udef, ios)

     deallocate(swe1d)
     deallocate(swe)
     deallocate(li)
#endif
  endif
  
  varfield = LVT_rc%udef

  do r=1, LVT_rc%lnr
     do c=1, LVT_rc%lnc
        if(isNaN(swe1d_ip(c+(r-1)*LVT_rc%lnc))) then 
           varfield(c,r) = LVT_rc%udef
        endif
        if(swe1d_ip(c+(r-1)*LVT_rc%lnc).lt.0.or.&
             swe1d_ip(c+(r-1)*LVT_rc%lnc).gt.1000) then 
           varfield(c,r) = LVT_rc%udef
        else
           varfield(c,r) = swe1d_ip(c+(r-1)*LVT_rc%lnc)
           varfield(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,varfield,vlevel=1,units="m")

  do r=1, LVT_rc%lnr
     do c=1, LVT_rc%lnc
        if(varfield(c,r).ne.-9999.0) then 
           varfield(c,r) = varfield(c,r)*1000.0
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,varfield,vlevel=1,units="kg/m2")
  
end subroutine readGlobSnowObs

!BOP
! 
! !ROUTINE: create_globsnow_filename
! \label{create_globsnow_filename}
!
! !INTERFACE:
subroutine create_globsnow_filename(odir,yr,mo,da,filename)
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
! This routine creates a timestamped filename for GlobSnow data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the GlobSnow file
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
  integer                      :: mo
  integer                      :: da
  character(len=*)             :: filename
!
!EOP

  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  
  filename = trim(odir)//'/'//trim(fyr)//'/L3A_daily_SWE/GlobSnow_SWE_L3A_'//&
       trim(fyr)//trim(fmo)//trim(fda)//'_v1.0.nc'
  
end subroutine create_globsnow_filename

