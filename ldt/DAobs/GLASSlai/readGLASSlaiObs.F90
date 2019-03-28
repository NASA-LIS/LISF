!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readGLASSlaiObs
! \label{readGLASSlaiObs}
! 
! !REVISION HISTORY: 
!  26 Mar 2019: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readGLASSlaiObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_DAobsDataMod
  use GLASSlai_obsMod
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the GLASS
! LAI data product.
!
!EOP 

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character*100     :: fname
  real              :: laiobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  GLASSlaiobs(n)%laiobs = LDT_rc%udef
  laiobs= LDT_rc%udef

  call create_GLASSlai_filename(GLASSlaiobs(n)%odir, &
       GLASSlaiobs(n)%source,&
       LDT_rc%yr, LDT_rc%doy, fname)
  
  inquire(file=trim(fname),exist=file_exists)
  if(file_exists) then
     
     write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
     call read_GLASSlai_data(n, fname, laiobs)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(laiobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
              GLASSlaiobs(n)%laiobs(c,r) = laiobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%lai_obs,&
       GLASSlaiobs(n)%laiobs,vlevel=1)

end subroutine readGLASSlaiObs

!BOP
! 
! !ROUTINE: read_GLASSlai_data
! \label{read_GLASSlai_data}
!
! !INTERFACE:
subroutine read_GLASSlai_data(n, fname, laiobs_ip)
! 
! !USES:   

  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use GLASSlai_obsMod

  implicit none
#if (defined USE_HDF4) 
#include "hdf.f90"
#endif
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname
  real                          :: laiobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real*8                        :: cornerlat(2), cornerlon(2)


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the GLASS LAI file and applies the data
!  quality flags to filter the data. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTGLASS AMSR-E file
!  \item[laiobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDF4)
  integer,  parameter     :: nc=7200, nr=3600
  integer                 :: gdopen,gdattach,gdrdfld
  integer                 :: gddetach,gdclose
  integer                 :: file_id,grid_id,region_id,iret,c,r,c1,r1
  character*50            :: grid_name,lai_name
  integer                 :: start(2), edge(2), stride(2)
  integer*2, allocatable  :: lai_raw_avhrr(:)
  integer*1, allocatable  :: lai_raw_modis(:)
  integer                 :: lat_off, lon_off
  real                    :: lai_in(GLASSlaiobs(n)%nc*GLASSlaiobs(n)%nr)
  logical*1               :: lai_data_b(GLASSlaiobs(n)%nc*GLASSlaiobs(n)%nr)
  logical*1               :: laiobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  real                    :: test(nc,nr)
  if(GLASSlaiobs(n)%source.eq."AVHRR") then 
     grid_name ="GLASS01B02"
  elseif(GLASSlaiobs(n)%source.eq."MODIS") then 
     grid_name ="GLASS01B01"
  endif

  file_id = gdopen(trim(fname),DFACC_READ)
  if (file_id.eq.-1)then
     write(LDT_logunit,*) "[ERR] Failed to open hdf file",fname
  end if
  
  lai_name = "LAI"

  grid_id = gdattach(file_id,grid_name)

  start(1)=0  !hdfeos lib uses 0-based count
  start(2)=0
  edge(1)=nc
  edge(2)=nr
  stride(1)=1
  stride(2)=1
  
  cornerlat(1)=GLASSlaiobs(n)%gridDesci(4)
  cornerlon(1)=GLASSlaiobs(n)%gridDesci(5)
  cornerlat(2)=GLASSlaiobs(n)%gridDesci(7)
  cornerlon(2)=GLASSlaiobs(n)%gridDesci(8)

  if(GLASSlaiobs(n)%source.eq."AVHRR") then 
     allocate(lai_raw_avhrr(nc*nr))
     iret = gdrdfld(grid_id,lai_name,start,stride,edge,lai_raw_avhrr)
     
     lai_data_b = .false. 
     
     lat_off = nint((cornerlat(1)+89.975)/0.05)+1
     lon_off = nint((cornerlon(1)+179.975)/0.05)+1
     
     do r=1,GLASSlaiobs(n)%nr
        do c=1,GLASSlaiobs(n)%nc
           c1 = c + lon_off
           r1 = nr - (r + lat_off) + 1
           
           if(lai_raw_avhrr(c1+(r1-1)*nc).gt.0.and.&
                lai_raw_avhrr(c1+(r1-1)*nc).ne.2550) then 
              lai_in(c+(r-1)*GLASSlaiobs(n)%nc) =&
                   lai_raw_avhrr(c1+(r1-1)*nc)*0.01
              lai_data_b(c+(r-1)*GLASSlaiobs(n)%nc) =  .true. 
           else
              lai_in(c+(r-1)*GLASSlaiobs(n)%nc) = -9999.0
              lai_data_b(c+(r-1)*GLASSlaiobs(n)%nc) = .false. 
           endif
        enddo
     enddo
     deallocate(lai_raw_avhrr)
  elseif(GLASSlaiobs(n)%source.eq."MODIS") then 
     allocate(lai_raw_modis(nc*nr))
     iret = gdrdfld(grid_id,lai_name,start,stride,edge,lai_raw_modis)
     lai_data_b = .false. 
     
     lat_off = nint((cornerlat(1)+89.975)/0.05)+1
     lon_off = nint((cornerlon(1)+179.975)/0.05)+1
     
     do r=1,GLASSlaiobs(n)%nr
        do c=1,GLASSlaiobs(n)%nc
           c1 = c + lon_off
           r1 = nr - (r + lat_off) + 1
           
           if(lai_raw_modis(c1+(r1-1)*nc).gt.0.and.&
                lai_raw_modis(c1+(r1-1)*nc).ne.2550) then 
              lai_in(c+(r-1)*GLASSlaiobs(n)%nc) =&
                   lai_raw_modis(c1+(r1-1)*nc)*0.1
              lai_data_b(c+(r-1)*GLASSlaiobs(n)%nc) =  .true. 
           else
              lai_in(c+(r-1)*GLASSlaiobs(n)%nc) = -9999.0
              lai_data_b(c+(r-1)*GLASSlaiobs(n)%nc) = .false. 
           endif
        enddo
     enddo
     deallocate(lai_raw_modis)

  endif

!  open(100,file='test_inp.bin',form='unformatted')
!  write(100) lai_in
!  close(100)


  iret=gddetach(grid_id)
  iret=gdclose(file_id)

  if(LDT_isLDTatAfinerResolution(n,0.05)) then

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
          GLASSlaiobs(n)%nc*GLASSlaiobs(n)%nr, &
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          GLASSlaiobs(n)%w11,GLASSlaiobs(n)%w12,&
          GLASSlaiobs(n)%w21,GLASSlaiobs(n)%w22,&
          GLASSlaiobs(n)%n11,GLASSlaiobs(n)%n12,&
          GLASSlaiobs(n)%n21,GLASSlaiobs(n)%n22,LDT_rc%udef,iret)
  else
     call upscaleByAveraging(GLASSlaiobs(n)%nc*GLASSlaiobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_rc%udef, GLASSlaiobs(n)%n11,&
          lai_data_b,lai_in, laiobs_b_ip, laiobs_ip)
  endif
  
#endif

end subroutine read_GLASSlai_data


!BOP
! !ROUTINE: create_GLASSlai_filename
! \label{create_GLASSlai_filename}
! 
! !INTERFACE: 
subroutine create_GLASSlai_filename(ndir, source, yr, doy, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: source
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the GLASS LAI filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the GLASS LAI data directory
!  \item[yr]  current year
!  \item[mo]  current doy
!  \item[filename] Generated GLASS LAI filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  if(source.eq."AVHRR") then 
     filename = trim(ndir)//'/'//trim(fyr)//'/GLASS01B02.V04.A'//&
          trim(fyr)//trim(fdoy)//'.hdf'
  elseif(source.eq."MODIS") then 
     filename = trim(ndir)//'/'//trim(fyr)//'/GLASS01B01.V04.A'//&
          trim(fyr)//trim(fdoy)//'.hdf'

  endif
end subroutine create_GLASSlai_filename


