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
! !ROUTINE: readUSDMObs
! \label{readUSDMObs}
!
! !INTERFACE: 
subroutine readUSDMObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod
  use LVT_histDataMod
  use USDM_obsMod, only : USDMObs
  use map_utils

#if (defined USE_GDAL)
  use, INTRINSIC:: iso_c_binding
  use fortranc
  use gdal
#endif
  implicit none
!
! !INPUT PARAMETERS: 
  integer,     intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The USDM output is available at weekly intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  8 Mar 2017   Sujay Kumar, Initial Specification
! 
!EOP

  character*100          :: filename
  integer                :: ftn 
  logical                :: file_exists
  integer                :: nid, ios
  logical*1, allocatable :: li(:)
  real                   :: lat1,lon1
  integer                :: c,r,t,kk
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: gridDesc(6)
  integer                :: c1,r1
  real                   :: lat,lon, dres
  real                   :: udef
  real                   :: cornerlat1, cornerlat2
  real                   :: cornerlon1, cornerlon2
  real                   :: timenow
  logical                :: alarmCheck
#if (defined USE_GDAL)
  TYPE(gdaldriverh)      :: driver
  TYPE(gdaldataseth)     :: ds
  TYPE(gdalrasterbandh)  :: band
  INTEGER(kind=c_int)    :: xsize, ysize
  REAL(kind=c_double)    :: x1, y1, x2, y2, gt(6)
  INTEGER(kind=c_int)    :: i1, j1, k1, i, j, k, ierr
  REAL,ALLOCATABLE       :: usdm_inp(:),zval(:,:)

#endif

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
       LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  if(USDMobs(source)%startflag.or.alarmCheck) then 

#if (defined USE_GDAL) 
     call GDALAllRegister()

     call create_USDM_filename(USDMobs(source)%odir,&
          LVT_rc%dyr(source),&
          LVT_rc%dmo(source),&
          LVT_rc%dda(source),&
          filename)

     inquire(file=trim(filename),exist=file_exists) 
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading USDM file ',trim(filename)
        
        driver =  gdalgetdriverbyname('GTiff'//CHAR(0))
        ds = gdalopen(TRIM(filename)//CHAR(0), GA_ReadOnly)
        if (.not.gdalassociated(ds)) then
           write(LVT_logunit,*) '[ERR] opening dataset on file ',TRIM(filename), ' failed'
           call LVT_endrun()
        ENDIF

        ierr = gdalgetgeotransform(ds, gt)
        xsize = gdalgetrasterxsize(ds)
        ysize = gdalgetrasterysize(ds)

        CALL gdalapplygeotransform(gt, 0.5_c_double, 0.5_c_double, x2, y1)
        
        band = gdalgetrasterband(ds, 1)
        if (.NOT.gdalassociated(band)) THEN
           write(LVT_logunit,*) '[ERR] failed getting raster band from GeoTIFF dataset on file ',TRIM(filename)
           call LVT_endrun()           
        endif

        allocate(zval(xsize, ysize))
        ierr = gdalrasterio_f(band, GF_Read, 0, 0, zval)
        if (ierr /= 0) THEN
           write(LVT_logunit,*) '[ERR] reading data from GeoTIFF dataset on file ',TRIM(filename)
           call LVT_endrun()
        endif
        call gdalclose(ds)

        dres = abs(gt(6))

        x1 = (xsize-1)*dres + x2
        y2 = y1 - (ysize-1)*dres

!        print*, gt
!        print*, x1, x2
!        print*, y1, y2
!        print*, xsize, ysize

        
        if(USDMobs(source)%startFlag) then 
           USDMobs(source)%startflag = .false.
           
           
           cornerlat1 = max(y2,nint((LVT_rc%gridDesc(4)-y2)/dres)*dres+y2-2*dres)
           cornerlon1 = max(x2,nint((LVt_rc%gridDesc(5)-x2)/dres)*dres+x2-2*dres)
           cornerlat2 = min(y1,nint((LVT_rc%gridDesc(7)-y2)/dres)*dres+y2+2*dres)
           cornerlon2 = min(x1,nint((LVT_rc%gridDesc(8)-x2)/dres)*dres+x2+2*dres)
                      
           USDMobs(source)%nr = nint((cornerlat2 - cornerlat1)/dres)+1
           USDMobs(source)%nc = nint((cornerlon2 - cornerlon1)/dres)+1

           allocate(USDMObs(source)%n11(USDMobs(source)%nc*USDMobs(source)%nr))
           
           !filling the items needed by the interpolation library
           USDMobs(source)%gridDesc(1) = 0  !input is EASE grid
           USDMobs(source)%gridDesc(2) = USDMobs(source)%nc
           USDMobs(source)%gridDesc(3) = USDMobs(source)%nr
           USDMobs(source)%gridDesc(4) = cornerlat1
           USDMobs(source)%gridDesc(5) = cornerlon1
           USDMobs(source)%gridDesc(7) = cornerlat2
           USDMobs(source)%gridDesc(8) = cornerlon2
           USDMobs(source)%gridDesc(6) = 128
           USDMobs(source)%gridDesc(9) = dres
           USDMobs(source)%gridDesc(10) = dres
           USDMobs(source)%gridDesc(20) = 64
           
           call upscaleByAveraging_input(USDMobs(source)%gridDesc,&
                LVT_rc%gridDesc,USDMobs(source)%nc*USDMobs(source)%nr,&
                LVT_rc%lnc*LVT_rc%lnr,USDMobs(source)%n11)
           

           call map_set(PROJ_LATLON,&
                USDMobs(source)%gridDesc(4),&
                USDMobs(source)%gridDesc(5), &
                0.0, &
                USDMobs(source)%gridDesc(9),&
                USDMobs(source)%gridDesc(10), &
                0.0, &
                USDMobs(source)%nc, &
                USDMobs(source)%nr, &
                USDMobs(source)%map_proj)           
        endif
        
        allocate(usdm_inp(USDMobs(source)%nc*USDMobs(source)%nr))
        allocate(li(USDMobs(source)%nc*USDMobs(source)%nr))

        do r=1,USDMobs(source)%nr
           do c=1,USDMobs(source)%nc
              call ij_to_latlon(USDMobs(source)%map_proj, & 
                   float(c), float(r), & 
                   lat, lon)

              r1 = (ysize - (nint((lat-y2)/dres)+1) + 1)
              c1 = nint((lon-x2)/dres) + 1

!              usdm_inp(c+((USDMobs(source)%nr-r+1)-1)*&
!                   USDMobs(source)%nc) = & 
!                   zval(c1,r1)

              usdm_inp(c+(r-1)*&
                   USDMobs(source)%nc) = & 
                   zval(c1,r1)
           enddo
        enddo
              
        deallocate(zval)

!        print*, ''
!        print*, USDMobs(source)%nc, USDMobs(source)%nr
!        print*, cornerlat1, cornerlon1

!        open(100,file='test.bin',form='unformatted')
!        write(100) usdm_inp
!       close(100)

        li = .false. 
        udef = 9.0

        do r=1,USDMobs(source)%nr
           do c=1,USDMobs(source)%nc
              if(usdm_inp(c+(r-1)*USDMobs(source)%nc).ne.udef) then 
                 li(c+(r-1)*USDMobs(source)%nc) = .true. 
              else
                 usdm_inp(c+(r-1)*USDMobs(source)%nc) = -9999.0
                 li(c+(r-1)*USDMobs(source)%nc) = .false. 
              endif
           enddo
        enddo
        call upscaleByMode(USDMobs(source)%nc*USDMobs(source)%nr,&
             LVT_rc%lnc*LVT_rc%lnr,LVT_rc%udef,&
             USDMobs(source)%n11,li,usdm_inp,lo,USDMobs(source)%drcategory)

        do r=1, LVT_rc%lnr
           do c=1, LVT_rc%lnc
              if(lo(c+(r-1)*LVT_rc%lnc)) then
! Adding 1 to make the drought categories to go from 1 to 5 
! instead of 0 to 4 (D0 to D4)
! 1 represents D0 and 5 represents D4
 
                 varfield(c,r) = USDMobs(source)%drcategory(&
                      c+(r-1)*LVT_rc%lnc)+ 1                         
              else
                 varfield(c,r) = LVT_rc%udef
              endif
           enddo
        enddo
        
!        open(100,file='test_out.bin',form='unformatted')
!        write(100) varfield
!        close(100)
!        stop
     else
        varfield  = LVT_rc%udef
     endif
  else
     varfield  = LVT_rc%udef
#endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_DR_CATEGORY,source,varfield,&
       vlevel=1,units="-")

end subroutine readUSDMObs

!BOP
! 
! !ROUTINE: create_USDM_filename
! \label{create_USDM_filename}
!
! !INTERFACE: 
subroutine create_USDM_filename(odir,yr,mo,da,filename)
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
! This routine creates a timestamped filename for USDM data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      USDM base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the USDM_LH file
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
  character*2             :: fda
  character*2             :: fmo

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  
  filename = trim(odir)//'/USDM_'//trim(fyr)//trim(fmo)//&
       trim(fda)//'.tif'

end subroutine create_USDM_filename


