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
! !ROUTINE: readOCO2_SIFObs
! \label{readOCO2_SIFObs}
!
! !INTERFACE: 
subroutine readOCO2_SIFObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use OCO2_SIFobsMod
  use map_utils
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
!  This subroutine reads the OCO2 SIF data (from 757 and 771 nm channels). 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  20 Apr 2018: Sujay Kumar, Initial Specification
! 
!EOP

  real                   :: currTime
  logical                :: alarmCheck 
  integer                :: c,r, k
  real                   :: col,row
  integer                :: col_min,col_max,row_min,row_max
  real                   :: lat_min,lat_max,lon_min,lon_max
  integer                :: flag
  integer                :: ftn
  character*100          :: fname
  real, allocatable      :: sif_temp(:),lat(:), lon(:),flat(:,:),flon(:,:)
  real                   :: sif_ip(LVT_rc%lnc,LVT_rc%lnr)
  integer                :: nsif_ip(LVT_rc%lnc,LVT_rc%lnr)
  logical                :: file_exists
  integer                :: dim1id,sifid,latid, lonid,flatid,flonid
  integer                :: ndata
  
  nsif_ip = 0.0
  sif_ip  = 0.0

  currTime = float(LVT_rc%dhr(source))*3600+ &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmCheck = (mod(currtime,86400.0).eq.0)

  if(OCO2sifobs(source)%startFlag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 

     if(OCO2sifobs(source)%startFlag) then 
        OCO2sifobs(source)%startFlag = .false. 
     endif

     call create_OCO2sif_filename(OCO2sifobs(Source)%odir, &
          LVT_rc%dyr(source),&
          LVT_rc%dmo(source),&
          LVT_rc%dda(source),&
          fname)
     
     inquire(file=trim(fname),exist=file_exists) 

     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading OCO-2 SIF file ',trim(fname)
            
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

        call LVT_verify(nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=ftn),&
             'Error opening file'//trim(fname))

        call LVT_verify(nf90_inq_dimid(ftn, 'sounding_dim',dim1Id),&
             'Error nf90_inq_dimid: sounding_dim')

        call LVT_verify(nf90_inquire_dimension(ftn, dim1Id, len=ndata), &
             'Error nf90_inquire_dimension: ')

        allocate(sif_temp(ndata))
        allocate(lat(ndata))
        allocate(lon(ndata))
        allocate(flat(4,ndata))
        allocate(flon(4,ndata))

        call LVT_verify(nf90_inq_varid(ftn, 'latitude',latid),&
             'error nf90_inq_varid: latitude')

        call LVT_verify(nf90_get_var(ftn,latid, lat),&
             'Error nf90_get_var: latitude')

        call LVT_verify(nf90_inq_varid(ftn, 'longitude',lonid),&
             'error nf90_inq_varid: longitude')

        call LVT_verify(nf90_get_var(ftn,lonid, lon),&
             'Error nf90_get_var: longitude')


        call LVT_verify(nf90_inq_varid(ftn, &
             'footprint_vertex_latitude',flatid),&
             'error nf90_inq_varid: footprint_vertex_latitude')

        call LVT_verify(nf90_get_var(ftn,flatid, flat),&
             'Error nf90_get_var: footprint_vertex_latitude')

        call LVT_verify(nf90_inq_varid(ftn,&
             'footprint_vertex_longitude',flonid),&
             'error nf90_inq_varid: footprint_vertex_longitude')

        call LVT_verify(nf90_get_var(ftn,flonid, flon),&
             'Error nf90_get_var: footprint_vertex_longitude')

        
        if(OCO2SIFobs(source)%channel.eq."757nm") then 
           call LVT_verify(nf90_inq_varid(ftn, 'SIF_757nm',sifid),&
                'error nf90_inq_varid')
           
           call LVT_verify(nf90_get_var(ftn,sifid, sif_temp),&
                'Error nf90_get_var: SIF_757nm')
        elseif(OCO2SIFobs(source)%channel.eq."771nm") then 
           call LVT_verify(nf90_inq_varid(ftn, 'SIF_771nm',sifid),&
                'error nf90_inq_varid')
           
           call LVT_verify(nf90_get_var(ftn,sifid, sif_temp),&
                'Error nf90_get_var: SIF_771nm')

        endif
        call LVT_verify(nf90_close(ftn),&
             'Error in nf90_close')
        
#endif
!Assuming a 2.25kmx1.3km grid
        do k=1,ndata

           lat_min = minval(flat(:,k))
           lat_max = maxval(flat(:,k))
           lon_min = minval(flon(:,k))
           lon_max = maxval(flon(:,k))

           call latlon_to_ij(LVT_domain%lvtproj,&
                lat_min,lon_min,col,row)

           col_min = nint(col)
           row_min = nint(row)

           call latlon_to_ij(LVT_domain%lvtproj,&
                lat_max,lon_max,col,row)
           col_max = nint(col)
           row_max = nint(row)
        
           if((col_min.ge.1.and.col_min.le.LVT_rc%lnc).and.&
                (col_max.ge.1.and.col_max.le.LVT_rc%lnc).and.&
                (row_min.ge.1.and.row_min.le.LVT_rc%lnr).and.& 
                (row_max.ge.1.and.row_max.le.LVT_rc%lnr)) then 

              if(.not.isNaN(sif_temp(k))) then 

                 do r=row_min, row_max
                    do c=col_min,col_max
                       sif_ip(c,r) = sif_ip(c,r) + sif_temp(k)
                       nsif_ip(c,r) = nsif_ip(c,r) + 1
                    enddo
                 enddo
              endif
           endif
        enddo

        deallocate(sif_temp)
        deallocate(lat)
        deallocate(lon)

     endif

  end if

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nsif_ip(c,r).gt.0) then 
           sif_ip(c,r) = sif_ip(c,r)/nsif_ip(c,r)
        else
           sif_ip(c,r) = LVT_rc%udef
        endif
     enddo
  enddo
  
  call LVT_logSingleDataStreamVar(LVT_MOC_SIF,source,sif_ip,&
       vlevel=1,units="mW/m^2/nm/sr")
  
end subroutine readOCO2_SIFObs

!BOP
! 
! !ROUTINE: create_OCO2sif_filename
! \label{create_OCO2sif_filename}
!
! !INTERFACE: 
subroutine create_OCO2sif_filename(odir,yr,mo,da,filename)
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
! This routine creates a timestamped filename for OCO2 SIF data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      OCO2 SIF base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the OCO2 SIF file
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
! The standard filenames are : 
!  oco2_LtSIF_[AcquisitionDate]_[ShortBuildId]_[ProductionDateTime][Source].nc4
!
! For the sake of simplicity the archived files have been renamed by 
! stripping off the '[ProductionDateTime]' part. 
! 
!EOP
  
  character*4             :: fyr
  character*1             :: fyr2(4)
  character*2             :: fmo
  character*2             :: fda

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  
  filename = trim(odir)//'/'//trim(fyr)//&
       '/oco2_LtSIF_'//&
       trim(fyr(3:3))//trim(fyr(4:4))//trim(fmo)//trim(fda)//'_B8100r.nc4'
  
end subroutine create_OCO2sif_filename


