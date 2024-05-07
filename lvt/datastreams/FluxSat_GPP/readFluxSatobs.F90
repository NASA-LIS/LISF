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
! !ROUTINE: readFluxSatObs
! \label{readFluxSatObs}
!
! !INTERFACE: 
subroutine readFluxSatObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_timeMgrMod,  only : LVT_calendar
  use LVT_logMod,      only : LVT_logunit, LVT_verify
  use LVT_histDataMod
  use FluxSat_obsMod, only : FluxSatObs

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The FluxSat output is available at daily intervals.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 Feb 2021   Sujay Kumar  Initial Specification
! 
!EOP

  character*100           :: filename
  logical                 :: file_exists
  integer                 :: nid, ios
  integer                 :: gppid, rowId, colId
  integer                 :: neeid
  integer                 :: nc,nr
  real,  allocatable      :: gpp(:,:)
  real,  allocatable      :: gpp1d(:)
  real,  allocatable      :: nee(:,:)
  real,  allocatable      :: nee1d(:)
  logical*1, allocatable  :: li(:)
  integer                 :: c,r,t,kk
  type(ESMF_Time)         :: currTime
  type(ESMF_TimeInterval) :: ts
  integer                 :: cyr, cmo, cda, chr, cmn, css
  integer                 :: status
  logical*1               :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                    :: varfield1d(LVT_rc%lnc*LVT_rc%lnr)
  real                    :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  logical                 :: alarmcheck
  real                    :: timenow
 

  timenow = float(LVT_rc%dhr(source))*3600 +&
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  varfield = -9999.0

  if(FluxSatobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then
    
     LVT_rc%resetFlag(source) = .false. 

     call create_fluxsat_gpp_filename(FluxSatobs(Source)%odir, &
          LVT_rc%dyr(source),LVT_rc%dmo(source),filename)

     inquire(file=trim(filename),exist=file_exists) 


     if(file_exists) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)  
 
        write(LVT_logunit,*) '[INFO] Reading FluxSat GPP file ',trim(filename)
                
        ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
        call LVT_verify(ios, 'Error opening file '//trim(filename))
        
        ios = nf90_inq_varid(nid,'GPP',gppid)
        call LVT_verify(ios, 'Error nf90_inq_varid: GPP')
        ! dimensions
        ios = nf90_inq_dimid(nid,'lat',rowId)
        call LVT_verify(ios, 'Error nf90_inq_dimid: lat')
        
        ios = nf90_inquire_dimension(nid,rowId, len=nr)
        call LVT_verify(ios, 'Error nf90_inquire_dimension: lat')
        
        ios = nf90_inq_dimid(nid,'lon',colId)
        call LVT_verify(ios, 'Error nf90_inq_dimid: lon')
        
        ios = nf90_inquire_dimension(nid,colId, len=nc)
        call LVT_verify(ios, 'Error nf90_inquire_dimension: lon')

        allocate(gpp(nc,nr))
        allocate(gpp1d(nc*nr))
        allocate(li(nc*nr))
        li = .true. 
        
        !values
        ios = nf90_get_var(nid,gppid, gpp,&
             start=(/1,1,LVT_rc%dda(source)/),&
             count=(/nc,nr,1/))
        call LVT_verify(ios, 'Error nf90_get_var: GPP')
        
        ios = nf90_close(nid)
        call LVT_verify(ios, 'Error in nf90_close')

#endif      

        gpp1d = -9999.0
        do r=1,nr
           do c=1,nc
              if(isNaN(gpp(c,nr-r+1))) then 
                 gpp1d(c+(r-1)*nc) = LVT_rc%udef
                 li(c+(r-1)*nc) = .false. 
              elseif(gpp(c,nr-r+1).eq.-9999.0) then 
                 gpp1d(c+(r-1)*nc) = LVT_rc%udef
                 li(c+(r-1)*nc) = .false. 
              else
                 gpp1d(c+(r-1)*nc) = gpp(c,nr-r+1)
                 li(c+(r-1)*nc) = .true. 
              endif
           enddo
        enddo

        call neighbor_interp(LVT_rc%gridDesc,li,gpp1d,&
             lo, varfield1d, nc*nr, &
             LVT_rc%lnc*LVT_rc%lnr,&
             FluxSatobs(Source)%rlat, FluxSatobs(Source)%rlon, &
             FluxSatobs(Source)%n11,LVT_rc%udef, ios)
        
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(varfield1d(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then
!convert to g/m2s
                 varfield(c,r) = varfield1d(c+(r-1)*LVT_rc%lnc)/86400
              else
                 varfield(c,r) = LVT_rc%udef
              endif
           enddo
        enddo

        deallocate(gpp)
        deallocate(gpp1d)
        deallocate(li)
     endif
  endif
  
  call LVT_logSingleDataStreamVar(LVT_MOC_GPP,source,varfield,&
       vlevel=1,units="g/m2s")
  
end subroutine readFluxSatObs

!BOP
! 
! !ROUTINE: create_fluxsat_gpp_filename
! \label{create_fluxsat_gpp_filename}
!
! !INTERFACE: 
subroutine create_fluxsat_gpp_filename(odir,yr,mo,filename)
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
! This routine creates a timestamped filename for FluxSat data files 
! based on the given date (year, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      Fluxsat base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the FluxSat file
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
  character(len=*)             :: filename
!EOP

  character*4             :: fyr
  character*4             :: fmo
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  
  filename = trim(odir)//'/'//trim(fyr)//'/GPP_FluxSat_daily_v2.0_'//&
       trim(fyr)//trim(fmo)//'.nc'
  
end subroutine create_fluxsat_gpp_filename

