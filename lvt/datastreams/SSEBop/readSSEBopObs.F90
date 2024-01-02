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
! !ROUTINE: readSSEBopObs
! \label{readSSEBopObs}
!
! !INTERFACE: 
subroutine readSSEBopObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod
  use LVT_histDataMod
  use SSEBop_obsMod, only : SSEBopObs
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
!   The SSEBop output is available at monthly intervals. So 
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
#if (defined USE_GDAL)
  TYPE(gdaldriverh)      :: driver
  TYPE(gdaldataseth)     :: ds
  TYPE(gdalrasterbandh)  :: band
  INTEGER(kind=c_int)    :: xsize, ysize
  REAL(kind=c_double)    :: x1, y1, x2, y2, gt(6)
  INTEGER(kind=c_int)    :: i1, j1, k1, i, j, k, ierr
  REAL,ALLOCATABLE       :: sseb_inp(:),zval(:,:)

#endif

  if((Ssebopobs(Source)%mo.ne.LVT_rc%d_nmo(source)).or.&
       LVT_rc%resetFlag(source)) then
     LVT_rc%resetFlag(source) = .false. 
     if(ssebopobs(source)%startFlag) then 
        Ssebopobs(Source)%yr = LVT_rc%dyr(source)
        Ssebopobs(Source)%mo = LVT_rc%dmo(source)
     endif

     Ssebopobs(Source)%et_var = LVT_rc%udef

#if (defined USE_GDAL) 
     call GDALAllRegister()

     call create_ssebop_filename(Ssebopobs(Source)%odir,&
          Ssebopobs(Source)%use_anomaly, &
          Ssebopobs(Source)%yr,&
          Ssebopobs(Source)%mo,filename)

     inquire(file=trim(filename),exist=file_exists) 
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading SSEBop file ',trim(filename)
        
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

        
        if(ssebopobs(source)%startFlag) then 
           ssebopobs(source)%startflag = .false.
           
           
           cornerlat1 = max(y2,nint((LVT_rc%gridDesc(4)-y2)/dres)*dres+y2-2*dres)
           cornerlon1 = max(x2,nint((LVt_rc%gridDesc(5)-x2)/dres)*dres+x2-2*dres)
           cornerlat2 = min(y1,nint((LVT_rc%gridDesc(7)-y2)/dres)*dres+y2+2*dres)
           cornerlon2 = min(x1,nint((LVT_rc%gridDesc(8)-x2)/dres)*dres+x2+2*dres)
                      
           ssebopobs(source)%nr = nint((cornerlat2 - cornerlat1)/dres)+1
           ssebopobs(source)%nc = nint((cornerlon2 - cornerlon1)/dres)+1

           allocate(SSEBopObs(source)%n11(ssebopobs(source)%nc*ssebopobs(source)%nr))
           
           !filling the items needed by the interpolation library
           ssebopobs(source)%gridDesc(1) = 0  !input is EASE grid
           ssebopobs(source)%gridDesc(2) = ssebopobs(source)%nc
           ssebopobs(source)%gridDesc(3) = ssebopobs(source)%nr
           ssebopobs(source)%gridDesc(4) = cornerlat1
           ssebopobs(source)%gridDesc(5) = cornerlon1
           ssebopobs(source)%gridDesc(7) = cornerlat2
           ssebopobs(source)%gridDesc(8) = cornerlon2
           ssebopobs(source)%gridDesc(6) = 128
           ssebopobs(source)%gridDesc(9) = dres
           ssebopobs(source)%gridDesc(10) = dres
           ssebopobs(source)%gridDesc(20) = 64
           
           call upscaleByAveraging_input(ssebopobs(source)%gridDesc,&
                LVT_rc%gridDesc,ssebopobs(source)%nc*ssebopobs(source)%nr,&
                LVT_rc%lnc*LVT_rc%lnr,ssebopobs(source)%n11)
           

           call map_set(PROJ_LATLON,&
                ssebopobs(source)%gridDesc(4),&
                ssebopobs(source)%gridDesc(5), &
                0.0, &
                ssebopobs(source)%gridDesc(9),&
                ssebopobs(source)%gridDesc(10), &
                0.0, &
                ssebopobs(source)%nc, &
                ssebopobs(source)%nr, &
                ssebopobs(source)%map_proj)           
        endif
        
        allocate(sseb_inp(ssebopobs(source)%nc*ssebopobs(source)%nr))
        allocate(li(ssebopobs(source)%nc*ssebopobs(source)%nr))

        do r=1,ssebopobs(source)%nr
           do c=1,ssebopobs(source)%nc
              call ij_to_latlon(ssebopobs(source)%map_proj, & 
                   float(c), float(r), & 
                   lat, lon)

              r1 = (ysize - (nint((lat-y2)/dres)+1) + 1)
              c1 = nint((lon-x2)/dres) + 1

!              sseb_inp(c+((ssebopobs(source)%nr-r+1)-1)*&
!                   ssebopobs(source)%nc) = & 
!                   zval(c1,r1)

              sseb_inp(c+(r-1)*&
                   ssebopobs(source)%nc) = & 
                   zval(c1,r1)
           enddo
        enddo
              
        deallocate(zval)

!        print*, ''
!        print*, ssebopobs(source)%nc, ssebopobs(source)%nr
!        print*, cornerlat1, cornerlon1

!        open(100,file='test.bin',form='unformatted')
!        write(100) sseb_inp
!       close(100)

        li = .false. 
        if(ssebopobs(Source)%use_anomaly.eq.1) then 
           udef = 0.0
        else
           udef = 255.0
        endif

        do r=1,ssebopobs(source)%nr
           do c=1,ssebopobs(source)%nc
              if(sseb_inp(c+(r-1)*ssebopobs(source)%nc).ne.udef) then 
                 li(c+(r-1)*ssebopobs(source)%nc) = .true. 
              else
                 sseb_inp(c+(r-1)*ssebopobs(source)%nc) = -9999.0
                 li(c+(r-1)*ssebopobs(source)%nc) = .false. 
              endif
           enddo
        enddo
        
        call upscaleByAveraging(ssebopobs(source)%nc*ssebopobs(source)%nr,&
             LVT_rc%lnc*LVT_rc%lnr,LVT_rc%udef,&
             Ssebopobs(Source)%n11,li,sseb_inp,lo,Ssebopobs(Source)%et_var)

!        open(100,file='test_inp.bin',form='unformatted')
!        write(100) Ssebopobs(Source)%et_var
!        close(100)
!        stop
        
     endif

     do r=1, LVT_rc%lnr
        do c=1, LVT_rc%lnc
           if(lo(c+(r-1)*LVT_rc%lnc)) then 
              varfield(c,r) = Ssebopobs(Source)%et_var(c+(r-1)*LVT_rc%lnc)          
           else
              varfield(c,r) = LVT_rc%udef
           endif
        enddo
     enddo
     
     Ssebopobs(Source)%yr = LVT_rc%d_nyr(source)
     Ssebopobs(Source)%mo = LVT_rc%d_nmo(source)
  else
     varfield  = LVT_rc%udef
#endif
  endif

  if(ssebopobs(Source)%use_anomaly.eq.1) then 
     call LVT_logSingleDataStreamVar(LVT_MOC_ETa,source,varfield,&
          vlevel=1,units="%")
  else
     call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,varfield,&
          vlevel=1,units="W/m2")
  endif

end subroutine readSSEBopObs

!BOP
! 
! !ROUTINE: create_ssebop_filename
! \label{create_ssebop_filename}
!
! !INTERFACE: 
subroutine create_ssebop_filename(odir,use_anomaly, yr,mo,filename)
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
! This routine creates a timestamped filename for SSEBop_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the SSEBop_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: use_anomaly
  integer                      :: yr
  integer                      :: mo
  character(len=*)             :: filename
!
!EOP

  integer                 :: yr1, yr2
  character*4             :: fyr
  character*2             :: fyr1
  character*2             :: fmo

  write(unit=fyr, fmt='(i4.4)') yr
  read(unit=fyr,fmt='(i2.2,i2.2)') yr2, yr1
  write(unit=fyr1,fmt='(i2.2)') yr1

  write(unit=fmo, fmt='(i2.2)') mo
  
  if(use_anomaly.eq.0) then 
  
     filename = trim(odir)//'/'//trim(fyr)//'/m'//trim(fyr1)//&
          trim(fmo)//'modisSSEBopET.tif'
  else
     filename = trim(odir)//'/'//trim(fyr)//'/ma'//trim(fyr1)//&
          trim(fmo)//'.modisSSEBopET.tif'
  endif
end subroutine create_ssebop_filename


