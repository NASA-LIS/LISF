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
! !ROUTINE: readASOSWEObs
! \label{readASOSWEObs}
!
! !INTERFACE:
subroutine readASOSWEObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use ASOSWE_obsMod
  use UTM_utils
  use map_utils
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for ASOSWE SWE data.
! LVT expects the data to be provided in a timestamped manner. The raw 
! ASOSWE binary files are read, and the data is interpolated to the 
! LIS model output. The data is interpolated using the bilinear 
! averaging technique. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  6 May 2010: Sujay Kumar, Initial Specification
! 
!EOP
  integer                :: source
  real                   :: timenow
  logical                :: alarmCheck
  integer                :: yr, mo, da, hr
  integer                :: i,j,t,c,r
  integer                :: stn_col, stn_row
  real                   :: col,row
  character*100          :: asoswe_filename
  logical                :: file_exists
  integer                :: ftn, ios
  integer                :: status
  integer                :: stnindex,tind
  real                   :: offset
  integer                :: dim1id, dim2id, xid, yid, sweid, nx,ny
  integer                :: iret
  integer                :: zone
  real                   :: latdeg, londeg
  real,  allocatable     :: var(:,:)
  real,  allocatable     :: xval(:),yval(:)
  real                   :: var_ip(LVT_rc%lnc,LVT_rc%lnr)
  integer                :: nvar_ip(LVT_rc%lnc,LVT_rc%lnr)
  
!  real                   :: minlat,minlon, maxlat, maxlon

  var_ip  = LVT_rc%udef
  
  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
       LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  
  if(alarmCheck) then 
     call create_ASOSWE_swefilename(&
          asosweobs(source)%odir, LVT_rc%dyr(source), &
          LVT_rc%dmo(source), LVT_rc%dda(source), &
          asoswe_filename)
     
     inquire(file=trim(asoswe_filename),exist=file_exists)
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading ASO SWE SWE data ',&
             trim(asoswe_filename)
        
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

        ios = nf90_open(path=trim(asoswe_filename),mode=NF90_NOWRITE,ncid=ftn)
        call LVT_verify(ios, 'Error opening file, '//trim(asoswe_filename))

        ios = nf90_inq_dimid(ftn, 'x',dim1Id)
        call LVT_verify(ios, 'Error nf90_inq_dimid: x')
        
        ios = nf90_inquire_dimension(ftn, dim1Id, len=nx)
        call LVT_verify(ios, 'Error nf90_inquire_dimension:x')        

        ios = nf90_inq_dimid(ftn, 'y',dim2Id)
        call LVT_verify(ios, 'Error nf90_inq_dimid: y')
        
        ios = nf90_inquire_dimension(ftn, dim2Id, len=ny)
        call LVT_verify(ios, 'Error nf90_inquire_dimension:y')   

        allocate(var(nx,ny))

        allocate(xval(nx))
        allocate(yval(ny))

        ios = nf90_inq_varid(ftn, 'x',xid)
        call LVT_verify(ios, 'Error nf90_inq_varid: x')   

        ios = nf90_get_var(ftn, xid, xval)
        call LVT_verify(ios, 'Error nf90_get_var: x')

        ios = nf90_inq_varid(ftn, 'y',yid)
        call LVT_verify(ios, 'Error nf90_inq_varid: y')   

        ios = nf90_get_var(ftn, yid, yval)
        call LVT_verify(ios, 'Error nf90_get_var: y')
  
        ios = nf90_inq_varid(ftn, 'Band1',sweid)
        call LVT_verify(ios, 'Error nf90_inq_varid: Band1')     
        
        ios = nf90_get_var(ftn, sweid, var)
        call LVT_verify(ios, 'Error nf90_get_var: Band1')

        zone = 11
        var_ip = 0 
        nvar_ip = 0 

!        minlat = 10000
!        maxlat = -10000
!        minlon = 10000
!        maxlon = -10000

        do i=1,nx
           do j=1,ny
              call UTM2Geo(zone,yval(j), xval(i),latdeg,londeg)

!              if(latdeg.lt.minlat) minlat = latdeg
!              if(londeg.lt.minlon) minlon = londeg
!              if(latdeg.gt.maxlat) maxlat = latdeg
!              if(londeg.gt.maxlon) maxlon = londeg

              call latlon_to_ij(LVT_domain%lvtproj, &
                   latdeg, londeg, col, row)

              stn_col = nint(col)
              stn_row = nint(row)

              if(.not.isNaN(var(i,j)).and.var(i,j).ge.0) then 
                 var_ip(stn_col,stn_row) = var_ip(stn_col,stn_row) +& 
                      var(i,j)
                 nvar_ip(stn_col,stn_row) = nvar_ip(stn_col,stn_row) +1  
              endif
           enddo
        enddo
        
        ios = nf90_close(ftn)
        call LVT_verify(ios,'readASOSWEobs: Error in nf90_close')

        deallocate(var)
        deallocate(xval)
        deallocate(yval)

        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(nvar_ip(c,r).ne.0) then 
                 var_ip(c,r) = var_ip(c,r)/nvar_ip(c,r)
              else
                 var_ip(c,r) = LVT_rc%udef
              endif
           enddo
        enddo

#endif
        

        write(LVT_logunit,*) '[INFO] Successfully processed SWE from ',trim(asoswe_filename)
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,var_ip,vlevel=1,units="m")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(var_ip(c,r).ne.-9999.0) then 
           var_ip(c,r) = var_ip(c,r)*1000.0
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,var_ip,vlevel=1,units="kg/m2")


end subroutine readASOSWEObs

!BOP
! 
! !ROUTINE: create_ASOSWE_swefilename
! \label{create_ASOSWE_swefilename}
!
! !INTERFACE:
subroutine create_ASOSWE_swefilename(odir, yr,mo,da, asoswename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: asoswename
!
! !DESCRIPTION: 
! 
! This routine creates a timestamped ASOSWE filename. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir] ASOSWE base directory
!   \item[yr]   year of data 
!   \item[mo]   month of data
!   \item[da]   day of data
!   \item[asoswename]  Name of the ASOSWE file  
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  
  write(fyr, '(i4.4)' ) yr
  write(fmo, '(i2.2)' ) mo
  write(fda, '(i2.2)' ) da

  asoswename = trim(odir)//'/'//&
       '/ASO_swe_' &
       //trim(fyr)//trim(fmo)//trim(fda)//'.nc'
  
end subroutine create_ASOSWE_swefilename

