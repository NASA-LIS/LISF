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
! !ROUTINE: readGOES_LSTObs
! \label{readGOES_LSTObs}
!
! !INTERFACE: 
subroutine readGOES_LSTObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use map_utils
  use GOES_LSTobsMod, only : GOESLSTObs
          
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
!   The GOES LST data is available at daily intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
!
!  Main qc flag 
!   0 - LST produced, good quality, not necessary to examine more detailed QA
!   1 - LST produced, other quality, recommend examination of more detailed QA
!   2 - LST not produced due to cloud effects
!   3 - LST not produced primarily due to reasons other than cloud
!   
!  LST error flag key: 
!   0 = average LST error <=1k
!   1 = average LST error <=2k
!   2 = average LST error <=3k
!   3 = average LST error > 3k
!
!  Emissivity error flag key: 
!   0 = average LST error <=0.01
!   1 = average LST error <=0.02
!   2 = average LST error <=0.04
!   3 = average LST error > 0.04
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015: Sujay Kumar, Initial Specification
! 
!EOP

  character*100          :: GEfile, GWfile
  integer                :: ftn 
  logical                :: file_exists
  integer                :: ios, iret,ivar1, ivar2, ivar3
  real                   :: time
  integer*2, allocatable      :: lst(:,:)
  integer*2, allocatable      :: lat(:,:)
  integer*2, allocatable      :: lon(:,:)
  real, allocatable      :: lst_in(:,:) 
  real, allocatable      :: lat_in(:,:)
  real, allocatable      :: lon_in(:,:)
  real                   :: lat1,lon1
  integer                :: nc,nr
  real                   :: col,row
  integer                :: c,r,c1,r1,t,kk
  integer                :: lstid, latid, lonid,latdimid, londimid
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  integer                :: nvarfield(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: lst_ip_ge(LVT_rc%lnc,LVT_rc%lnr)
  integer                :: nlst_ip_ge(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: lst_ip_gw(LVT_rc%lnc,LVT_rc%lnr)
  integer                :: nlst_ip_gw(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: scale_lat, scale_lon, scale_lst
  real                   :: offset_lat, offset_lon, offset_lst
  real                   :: fill_lat, fill_lon, fill_lst
  integer                :: yr, mo, da, doy,hr
  
  GOESlstobs(Source)%lst = LVT_rc%udef
  nvarfield = 0 
  varfield = 0 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  
  time = LVT_rc%dhr(source)*3600+LVT_rc%dmn(source)*60+LVT_rc%dss(source)
  yr = LVT_rc%dyr(source)
  mo = LVT_rc%dmo(source)
  da = LVT_rc%dda(source)
  hr = LVT_rc%dhr(source)
  doy = LVT_rc%ddoy(source)
  
  if(mod(time,10800.0).eq.0.0.or.LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     call create_GWDISK_filename(GOESlstobs(Source)%odir, &
          yr, mo, da, doy, hr, GWfile)
     
     inquire(file=trim(GWfile),exist=file_exists) 
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading GOES LST file ',trim(GWfile)
        
        call LVT_verify(nf90_open(path=trim(GWfile),mode=NF90_NOWRITE, &
             ncid=ftn), &
             'error opening file in readGOES_LSTobs')
        
        iret = nf90_inq_dimid(ftn,"nlon_grid",londimid)
        iret = nf90_inq_dimid(ftn,"nlat_grid",latdimid)
        iret = nf90_inquire_dimension(ftn,londimid,len=nc)
        iret = nf90_inquire_dimension(ftn,latdimid,len=nr)
        
        allocate(lst(nc,nr))
        allocate(lat(nc,nr))
        allocate(lon(nc,nr))
        
        allocate(lst_in(nc,nr))
        allocate(lat_in(nc,nr))
        allocate(lon_in(nc,nr))
        
        ivar1 = nf90_inq_varid(ftn,'LST',lstid)
        ivar2 = nf90_inq_varid(ftn,'Pixel_Longitude_scaled',lonid)
        ivar3 = nf90_inq_varid(ftn,'Pixel_Latitude_scaled',latid)
        
        if(ivar1.eq.0.and.ivar2.eq.0.and.ivar3.eq.0) then 
           call LVT_verify(nf90_get_var(ftn,lstid,lst),&
                'Error in nf90_get_var: lst')
           call LVT_verify(nf90_get_var(ftn,latid,lat),&
                'Error in nf90_get_var: lat')
           call LVT_verify(nf90_get_var(ftn,lonid,lon),&
                'Error in nf90_get_var: lon')              
           
           call LVT_verify(nf90_get_att(ftn,lstid,"Factor",&
                scale_lst), "Error in nf90_get_att Factor")
           call LVT_verify(nf90_get_att(ftn,lstid,"Offset",&
                offset_lst), "Error in nf90_get_att Offset")
           call LVT_verify(nf90_get_att(ftn,lstid,"FillValue",&
                fill_lst), "Error in nf90_get_att Offset")
           
           call LVT_verify(nf90_get_att(ftn,latid,"Factor",&
                scale_lat), "Error in nf90_get_att Factor")
           call LVT_verify(nf90_get_att(ftn,latid,"Offset",&
                offset_lat), "Error in nf90_get_att Offset")
           call LVT_verify(nf90_get_att(ftn,latid,"FillValue",&
                fill_lat), "Error in nf90_get_att Offset")
           
           call LVT_verify(nf90_get_att(ftn,lonid,"Factor",&
                scale_lon), "Error in nf90_get_att Factor")
           call LVT_verify(nf90_get_att(ftn,lonid,"Offset",&
                offset_lon), "Error in nf90_get_att Offset")
           call LVT_verify(nf90_get_att(ftn,lonid,"FillValue",&
                fill_lon), "Error in nf90_get_att Offset")
           
           do r=1,nr
              do c=1,nc
                 lat_in(c,r) = (lat(c,r)+65536) * scale_lat + offset_lat
                 lon_in(c,r) = (lon(c,r)+65536)* scale_lon + offset_lon
                 lst_in(c,r) = lst(c,r) * scale_lst + offset_lst
              enddo
           enddo
           
        else
           lst_in = LVT_rc%udef
           lat_in = LVT_rc%udef
           lon_in = LVT_rc%udef
        endif
        deallocate(lst)
        deallocate(lat)
        deallocate(lon)

        call LVT_verify(nf90_close(ftn),&
             'Error in nf90_close')
        
        lst_ip_gw = 0
        nlst_ip_gw = 0 
        
        do r=1,nr
           do c=1,nc
              if(lst_in(c,r).gt.0) then 
                 call latlon_to_ij(LVT_domain%lvtproj, &
                      lat_in(c,r), lon_in(c,r),&
                      col, row)
                 c1 = nint(col)
                 r1 = nint(row)
                 if(c1.ge.1.and.c1.le.LVT_rc%lnc.and.&
                      r1.ge.1.and.r1.le.LVT_rc%lnr) then 
                    lst_ip_gw(c1,r1) = lst_ip_gw(c1,r1) + lst_in(c,r)
                    nlst_ip_gw(c1,r1) = nlst_ip_gw(c1,r1) + 1
                 endif
              endif
           enddo
        enddo
        
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(nlst_ip_gw(c,r).gt.0) then 
                 lst_ip_gw(c,r) = lst_ip_gw(c,r)/nlst_ip_gw(c,r)
              else
                 lst_ip_gw(c,r) = LVT_rc%udef
              endif
           enddo
        enddo

        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(lst_ip_gw(c,r).ne.LVT_rc%udef) then 
                 varfield(c,r) = varfield(c,r) + lst_ip_gw(c,r)
                 nvarfield(c,r) = nvarfield(c,r) + 1
              endif
           enddo
        enddo

        deallocate(lat_in)
        deallocate(lon_in)
        deallocate(lst_in)
     endif

     call create_GEDISK_filename(GOESlstobs(Source)%odir, &
          yr, mo, da, doy, hr, GEfile)

     inquire(file=trim(GEfile),exist=file_exists) 
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading GOES LST file ',trim(GEfile)
        
        call LVT_verify(nf90_open(path=trim(GEfile),mode=NF90_NOWRITE, &
             ncid=ftn), &
             'error opening file in readGOES_LSTobs')
        
        iret = nf90_inq_dimid(ftn,"nlon_grid",londimid)
        iret = nf90_inq_dimid(ftn,"nlat_grid",latdimid)
        iret = nf90_inquire_dimension(ftn,londimid,len=nc)
        iret = nf90_inquire_dimension(ftn,latdimid,len=nr)
        
        allocate(lst(nc,nr))
        allocate(lat(nc,nr))
        allocate(lon(nc,nr))
        
        allocate(lst_in(nc,nr))
        allocate(lat_in(nc,nr))
        allocate(lon_in(nc,nr))
        
        ivar1 = nf90_inq_varid(ftn,'LST',lstid)
        ivar2 = nf90_inq_varid(ftn,'Pixel_Longitude_scaled',lonid)
        ivar3 = nf90_inq_varid(ftn,'Pixel_Latitude_scaled',latid)
        
        if(ivar1.eq.0.and.ivar2.eq.0.and.ivar3.eq.0) then 
           call LVT_verify(nf90_get_var(ftn,lstid,lst),&
                'Error in nf90_get_var: lst')
           call LVT_verify(nf90_get_var(ftn,latid,lat),&
                'Error in nf90_get_var: lat')
           call LVT_verify(nf90_get_var(ftn,lonid,lon),&
                'Error in nf90_get_var: lon')              
           
           call LVT_verify(nf90_get_att(ftn,lstid,"Factor",&
                scale_lst), "Error in nf90_get_att Factor")
           call LVT_verify(nf90_get_att(ftn,lstid,"Offset",&
                offset_lst), "Error in nf90_get_att Offset")
           call LVT_verify(nf90_get_att(ftn,lstid,"FillValue",&
                fill_lst), "Error in nf90_get_att Offset")
           
           call LVT_verify(nf90_get_att(ftn,latid,"Factor",&
                scale_lat), "Error in nf90_get_att Factor")
           call LVT_verify(nf90_get_att(ftn,latid,"Offset",&
                offset_lat), "Error in nf90_get_att Offset")
           call LVT_verify(nf90_get_att(ftn,latid,"FillValue",&
                fill_lat), "Error in nf90_get_att Offset")
           
           call LVT_verify(nf90_get_att(ftn,lonid,"Factor",&
                scale_lon), "Error in nf90_get_att Factor")
           call LVT_verify(nf90_get_att(ftn,lonid,"Offset",&
                offset_lon), "Error in nf90_get_att Offset")
           call LVT_verify(nf90_get_att(ftn,lonid,"FillValue",&
                fill_lon), "Error in nf90_get_att Offset")
           
           do r=1,nr
              do c=1,nc
                 lat_in(c,r) = (lat(c,r)+65536) * scale_lat + offset_lat
                 lon_in(c,r) = (lon(c,r)+65536)* scale_lon + offset_lon
                 lst_in(c,r) = lst(c,r) * scale_lst + offset_lst
              enddo
           enddo
           
        else
           lst_in = LVT_rc%udef
           lat_in = LVT_rc%udef
           lon_in = LVT_rc%udef
        endif
        deallocate(lst)
        deallocate(lat)
        deallocate(lon)
        call LVT_verify(nf90_close(ftn),&
             'Error in nf90_close')
        
        lst_ip_ge = 0
        nlst_ip_ge = 0 
        
        do r=1,nr
           do c=1,nc
              if(lst_in(c,r).gt.0) then 
                 call latlon_to_ij(LVT_domain%lvtproj, &
                      lat_in(c,r), lon_in(c,r),&
                      col, row)
                 c1 = nint(col)
                 r1 = nint(row)
                 if(c1.ge.1.and.c1.le.LVT_rc%lnc.and.&
                      r1.ge.1.and.r1.le.LVT_rc%lnr) then 
                    lst_ip_ge(c1,r1) = lst_ip_ge(c1,r1) + lst_in(c,r)
                    nlst_ip_ge(c1,r1) = nlst_ip_ge(c1,r1) + 1
                 endif
              endif
           enddo
        enddo
        
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(nlst_ip_ge(c,r).gt.0) then 
                 lst_ip_ge(c,r) = lst_ip_ge(c,r)/nlst_ip_ge(c,r)
              else
                 lst_ip_ge(c,r) = LVT_rc%udef
              endif
           enddo
        enddo

        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(lst_ip_ge(c,r).ne.LVT_rc%udef) then 
                 varfield(c,r) = varfield(c,r) + lst_ip_ge(c,r)
                 nvarfield(c,r) = nvarfield(c,r) + 1
              endif
           enddo
        enddo

        deallocate(lat_in)
        deallocate(lon_in)
        deallocate(lst_in)
     end if
     
  endif
#endif
  
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nvarfield(c,r).gt.0) then 
           varfield(c,r) = varfield(c,r)/nvarfield(c,r)
        else
           varfield(c,r) = LVT_rc%udef
        endif
     enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT,source,varfield,vlevel=1,units="K")
  
end subroutine readGOES_LSTObs

!BOP
! 
! !ROUTINE: create_GWDISK_filename
! \label{create_GWDISK_filename}
!
! !INTERFACE: 
subroutine create_GWDISK_filename(odir,yr,mo,da, doy,hr,filename)
! 
! !USES:   
  use LVT_logMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for GOES LST data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the GOESLST_LH file
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
  integer                      :: doy
  integer                      :: hr
  character(len=*)             :: filename
!
!EOP

  character*200           :: list_file
  character*4             :: fyr
  character*2             :: fmo
  character*3             :: fdoy
  character*2             :: fhr
  integer                 :: ftn
  integer                 :: ierr

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fdoy, fmt='(i3.3)') doy
  write(unit=fhr, fmt='(i2.2)') hr
  
  list_file = 'ls '//trim(odir)//'/'//trim(fyr)//'/*GWDISK_*'//&
       trim(fyr)//trim(fdoy)//'_'//trim(fhr)//'*.nc > GOES_file'

  call system(trim(list_file))
  ftn=LVT_getNextUnitNumber()
  open(ftn,file='GOES_file',status='old',iostat=ierr)
  do while(ierr.eq.0) 
     read(ftn,'(a)',iostat=ierr) filename
     if(ierr.ne.0) then 
        exit
     endif
  enddo
  call LVT_releaseUnitNumber(ftn)

end subroutine create_GWDISK_filename

!BOP
! 
! !ROUTINE: create_GEDISK_filename
! \label{create_GEDISK_filename}
!
! !INTERFACE: 
subroutine create_GEDISK_filename(odir,yr,mo,da,doy,hr,filename)
! 
! !USES:   
  use LVT_logMod
  use LVT_timeMgrMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for GOES LST data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the GOESLST_LH file
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
  integer                      :: doy
  integer                      :: hr
  character(len=*)             :: filename
!
!EOP

  character*200           :: list_file
  character*4             :: fyr
  character*2             :: fmo
  character*3             :: fdoy
  character*2             :: fhr
  integer                 :: ftn
  integer                 :: ierr

  integer                 :: ts, doy1, mn, ss
  real*8                  :: time
  real                    :: gmt

  ts = -900
  mn = 0 
  ss = 0 

  call LVT_tick(time, doy1, gmt, yr,mo,da,hr,mn,ss,ts)

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fdoy, fmt='(i3.3)') doy
  write(unit=fhr, fmt='(i2.2)') hr
  
  list_file = 'ls '//trim(odir)//'/'//trim(fyr)//'/*GEDISK_*'//&
       trim(fyr)//trim(fdoy)//'_'//trim(fhr)//'*.nc > GOES_file'

  call system(trim(list_file))
  ftn=LVT_getNextUnitNumber()
  open(ftn,file='GOES_file',status='old',iostat=ierr)
  do while(ierr.eq.0) 
     read(ftn,'(a)',iostat=ierr) filename
     if(ierr.ne.0) then 
        exit
     endif
  enddo
  call LVT_releaseUnitNumber(ftn)

end subroutine create_GEDISK_filename


