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
! !ROUTINE: readGOME2_SIFObs
! \label{readGOME2_SIFObs}
!
! !INTERFACE: 
subroutine readGOME2_SIFObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use GOME2_SIFobsMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
          
  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  6 Jan 2016: Sujay Kumar, Initial Specification
! 
!EOP

  integer                :: c,r, tindex
  integer                :: flag
  integer                :: ftn
  character*100          :: fname
  logical*1              :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  logical*1              :: input_bitmap(&
       GOME2sifobs(source)%nc*GOME2sifobs(source)%nr)
  real                   :: ndvi_out(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: ndvi_inp(&
       GOME2sifobs(source)%nc*GOME2sifobs(source)%nr)
  real                   :: ndvi_temp(&
       GOME2sifobs(source)%nc,GOME2sifobs(source)%nr)
  real                   :: ndvi_ip(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: sif_out(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: sif_inp(&
       GOME2sifobs(source)%nc*GOME2sifobs(source)%nr)
  real                   :: sif_temp(&
       GOME2sifobs(source)%nc,GOME2sifobs(source)%nr)
  real                   :: sif_ip(LVT_rc%lnc,LVT_rc%lnr)
  logical                :: file_exists
  integer                :: nid,ndviid,sifid,iret
  integer                :: nc, nr
  
  ndvi_ip = LVT_rc%udef
  sif_ip  = LVT_rc%udef

  nc = GOME2sifobs(source)%nc
  nr = GOME2sifobs(source)%nr

  if(GOME2sifobs(source)%mo.ne.LVT_rc%d_nmo(source).or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     if(GOME2sifobs(source)%startFlag) then 
        GOME2sifobs(source)%startFlag = .false. 
     endif
     GOME2sifobs(source)%yr = LVT_rc%d_nyr(source)
     GOME2sifobs(source)%mo = LVT_rc%d_nmo(source)

     call create_GOME2sif_filename(GOME2sifobs(Source)%odir, &
          GOME2sifobs(source)%version,&
          LVT_rc%dyr(source),LVT_rc%dmo(source),&
          fname)
     
     inquire(file=trim(fname),exist=file_exists) 

     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading GOME-2 SIF file ',trim(fname)
            
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

        call LVT_verify(nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid),&
             'Error opening file'//trim(fname))

        call LVT_verify(nf90_inq_varid(nid, 'NDVI',ndviid),&
             'error nf90_inq_varid')
        call LVT_verify(nf90_inq_varid(nid, 'SIF',sifid),&
             'error nf90_inq_varid')

        call LVT_verify(nf90_get_var(nid,ndviid, ndvi_temp),&
             'Error nf90_get_var: NDVI')
        call LVT_verify(nf90_get_var(nid,sifid, sif_temp),&
             'Error nf90_get_var: SIF')

        call LVT_verify(nf90_close(nid),&
             'Error in nf90_close')

#endif
     else
        ndvi_temp  = LVT_rc%udef
        sif_temp   = LVT_rc%udef
     endif
        
     input_bitmap = .false. 
     
     do r=1,nr
        do c=1,nc
           ndvi_inp(c+(r-1)*nc) = ndvi_temp(c,nr-r+1)
           if(.not.isNaN(ndvi_inp(c+(r-1)*nc)).and.&
                ndvi_inp(c+(r-1)*nc).gt.0) then 
              input_bitmap(c+(r-1)*nc) = .true. 
           else
              ndvi_inp(c+(r-1)*nc) = LVT_rc%udef
           endif

        enddo
     enddo

     if(LVT_isAtAfinerResolution(gome2sifobs(source)%datares)) then
        call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
             ndvi_inp, output_bitmap, ndvi_out, &
             nc*nr, &
             LVT_rc%lnc*LVT_rc%lnr, &
             gome2sifobs(source)%rlat, & 
             gome2sifobs(source)%rlon, &
             gome2sifobs(source)%n11, &
             LVT_rc%udef, iret)
        
     else
        call upscaleByAveraging(&
             nc*nr, &
             LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
             gome2sifobs(source)%n11, input_bitmap, &
             ndvi_inp, output_bitmap, ndvi_out)
        
     endif
     
     do r=1, LVT_rc%lnr
        do c=1, LVT_rc%lnc
           if(ndvi_out(c+(r-1)*LVT_rc%lnc).gt.0) then 
              ndvi_ip(c,r) = ndvi_out(c+(r-1)*LVT_rc%lnc)
           else
              ndvi_ip(c,r) = LVT_rc%udef
           endif
        enddo
     enddo

     input_bitmap = .false. 
     
     do r=1,nr
        do c=1,nc
           sif_inp(c+(r-1)*nc) = sif_temp(c,nr-r+1)
           if(.not.isNaN(sif_inp(c+(r-1)*nc)).and.&
                sif_inp(c+(r-1)*nc).gt.0) then 
              input_bitmap(c+(r-1)*nc) = .true. 
           else
              sif_inp(c+(r-1)*nc) = LVT_rc%udef
           endif

        enddo
     enddo

     if(LVT_isAtAfinerResolution(gome2sifobs(source)%datares)) then
        call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
             sif_inp, output_bitmap, sif_out, &
             nc*nr, &
             LVT_rc%lnc*LVT_rc%lnr, &
             gome2sifobs(source)%rlat, & 
             gome2sifobs(source)%rlon, &
             gome2sifobs(source)%n11, &
             LVT_rc%udef, iret)
        
     else
        call upscaleByAveraging(&
             nc*nr, &
             LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
             gome2sifobs(source)%n11, input_bitmap, &
             sif_inp, output_bitmap, sif_out)
        
     endif
     
     do r=1, LVT_rc%lnr
        do c=1, LVT_rc%lnc
           if(sif_out(c+(r-1)*LVT_rc%lnc).gt.0) then 
              sif_ip(c,r) = sif_out(c+(r-1)*LVT_rc%lnc)
           else
              sif_ip(c,r) = LVT_rc%udef
           endif
        enddo
     enddo
     
  else
     ndvi_ip  = LVT_rc%udef
     sif_ip  = LVT_rc%udef
  endif
  
  call LVT_logSingleDataStreamVar(LVT_MOC_NDVI,source,ndvi_ip,&
       vlevel=1,units="-")

  call LVT_logSingleDataStreamVar(LVT_MOC_SIF,source,sif_ip,&
       vlevel=1,units="mW/m^2/nm/sr")
  
end subroutine readGOME2_SIFObs

!BOP
! 
! !ROUTINE: create_GOME2sif_filename
! \label{create_GOME2sif_filename}
!
! !INTERFACE: 
subroutine create_GOME2sif_filename(odir,version,yr,mo,filename)
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
! This routine creates a timestamped filename for GOME2 SIF data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GOME2 SIF base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the GOME2 SIF file
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
  integer                      :: mo
  character(len=*)             :: filename

!
!EOP
  
  character*4             :: fyr
  character*2             :: fmo

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  
!  filename = trim(odir)//'/'//trim(fyr)//&
!       '/ret_f_nr5_nsvd12_v26_waves734_nolog.grid_SIF__v20_1x1_'//&
!       trim(fyr)//trim(fmo)//'01_31.nc'
  filename = trim(odir)//'/'//trim(fyr)//&
       '/ret_f_nr5_nsvd12_v26_waves734_nolog.grid_SIF_'//trim(version)//'_'//&
       trim(fyr)//trim(fmo)//'01_31.nc'
  
end subroutine create_GOME2sif_filename


