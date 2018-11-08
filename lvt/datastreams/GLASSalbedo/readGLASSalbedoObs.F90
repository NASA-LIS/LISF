!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readGLASSalbedoObs
! \label{readGLASSalbedoObs}
!
! !INTERFACE: 
subroutine readGLASSalbedoObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use GLASSalbedoobsMod
          
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
!  7 Mar 2015: Sujay Kumar, Initial Specification
! 14 Sep 2018   Abheera Hazra - Adapted GLASSlai to GLASSalbedo 
! 
!EOP

  integer                :: c,r, tindex
  integer                :: flag
  integer                :: ftn
  character*200          :: fname
  logical                :: file_exists
  integer                :: iret  
  integer                :: nc, nr
  real                   :: varfield(6,LVT_rc%lnc,LVT_rc%lnr)

  varfield = LVT_rc%udef

  nc = GLASSalbedoobs(source)%nc
  nr = GLASSalbedoobs(source)%nr

  LVT_rc%resetFlag(source) = .false. 
  GLASSalbedoobs(source)%startFlag = .false. 
  
  call create_GLASSalbedo_filename(glassalbedoobs(Source)%odir, &
       LVT_rc%dyr(source),LVT_rc%ddoy(source),&
       fname)
  
  inquire(file=trim(fname),exist=file_exists) 
  
  ! Check if both files exist:
  if( file_exists ) then 
     write(LVT_logunit,*) '[INFO] Reading GLASS ALBEDO file ',trim(fname)
     if(GLASSalbedoobs(source)%source.eq."AVHRR") then
        call read_GLASS_ALBEDO_data(source, fname, varfield)
     elseif(GLASSalbedoobs(source)%source.eq."MODIS") then
        write(LVT_logunit,*)'[WARN] GLASS ALBEDO for MODIS is not yet supported by LVT.'
        stop 'End of code' 
    endif 
  else
     write(LVT_logunit,*)'[WARN] GLASS ALBEDO file missing: ',trim(fname)
     varfield  = LVT_rc%udef
  endif
  
  if(LVT_MOC_VISDIRALBEDO(source).ge.1) then
     call LVT_logSingleDataStreamVar(LVT_MOC_VISDIRALBEDO,source,varfield(1,:,:),vlevel=1,units="-")
  endif

  if(LVT_MOC_VISDIFALBEDO(source).ge.1) then 
     call LVT_logSingleDataStreamVar(LVT_MOC_VISDIFALBEDO,source,varfield(2,:,:),vlevel=1,units="-")
  endif
 
  if(LVT_MOC_NIRDIRALBEDO(source).ge.1) then
     call LVT_logSingleDataStreamVar(LVT_MOC_NIRDIRALBEDO,source,varfield(3,:,:),vlevel=1,units="-")
  endif

  if(LVT_MOC_NIRDIFALBEDO(source).ge.1) then 
     call LVT_logSingleDataStreamVar(LVT_MOC_NIRDIFALBEDO,source,varfield(4,:,:),vlevel=1,units="-")
  endif

  if(LVT_MOC_SWDIRALBEDO(source).ge.1) then
     call LVT_logSingleDataStreamVar(LVT_MOC_SWDIRALBEDO,source,varfield(5,:,:),vlevel=1,units="-")
  endif

  if(LVT_MOC_SWDIFALBEDO(source).ge.1) then 
     call LVT_logSingleDataStreamVar(LVT_MOC_SWDIFALBEDO,source,varfield(6,:,:),vlevel=1,units="-")
  endif



end subroutine readGLASSalbedoObs

!BOP
! 
! !ROUTINE: read_GLASS_ALBEDO_data
! \label{read_GLASS_ALBEDO_data}
!
! !INTERFACE:
subroutine read_GLASS_ALBEDO_data(source, fname, albedoobs_ip)
! 
! !USES:   

  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use GLASSalbedoobsMod

  implicit none
#if (defined USE_HDF4) 
#include "hdf.f90"
#endif
!
! !INPUT PARAMETERS: 
! 
  integer                       :: source
  character (len=*)             :: fname
  real                 :: albedoobs_ip(6,LVT_rc%lnc*LVT_rc%lnr)

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the GLASS ALBEDO file and applies the data
!  quality flags to filter the data. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTGLASS AMSR-E file
!  \item[albedoobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDF4)
  integer                 :: k
  real*8                  :: cornerlat(2), cornerlon(2)
  integer,  parameter     :: nc=7200, nr=3600
  integer                 :: gdopen,gdattach,gdrdfld
  integer                 :: gddetach,gdclose
  integer                 :: file_id,grid_id,region_id,c,r,c1,r1,i
  integer		  :: iret,iret1,iret2,iret3,iret4,iret5,iret6
  character*50            :: grid_name,albedo_name
  character*50            :: albedo_name_ws_vis,albedo_name_bs_vis,albedo_name_ws_nir,albedo_name_bs_nir,albedo_name_ws_sw,albedo_name_bs_sw
  integer                 :: start(2), edge(2), stride(2)
  integer*2, allocatable  :: albedo_raw_avhrr_bs_vis(:),albedo_raw_avhrr_ws_vis(:),albedo_raw_avhrr_bs_nir(:),albedo_raw_avhrr_ws_nir(:),albedo_raw_avhrr_bs_sw(:),albedo_raw_avhrr_ws_sw(:)
  integer*2, allocatable  :: albedo_raw_avhrr(:),albedo_raw_modis(:)
  integer                 :: lat_off, lon_off
  real                    :: albedo_in(GLASSalbedoobs(source)%nc*GLASSalbedoobs(source)%nr)
  real                    :: albedo_in1(6,GLASSalbedoobs(source)%nc*GLASSalbedoobs(source)%nr)

  logical*1               :: albedo_data_b(GLASSalbedoobs(source)%nc*GLASSalbedoobs(source)%nr)
  logical*1               :: albedo_data_b1(6,GLASSalbedoobs(source)%nc*GLASSalbedoobs(source)%nr)
  logical*1               :: albedoobs_b_ip(6,LVT_rc%lnc*LVT_rc%lnr)



  if(GLASSalbedoobs(source)%source.eq."AVHRR") then 
     grid_name ="GLASS02B05"
  elseif(GLASSalbedoobs(source)%source.eq."MODIS") then 
     grid_name ="GLASS02A06"
   endif

  file_id = gdopen(trim(fname),DFACC_READ)
  if (file_id.eq.-1)then
     write(LVT_logunit,*) "[ERR] Failed to open hdf file",fname
     call LVT_endrun()
  endif

     albedo_name_bs_vis = "ABD_BSA_VIS"

     albedo_name_ws_vis = "ABD_WSA_VIS"

     albedo_name_bs_nir = "ABD_BSA_NIR"

     albedo_name_ws_nir = "ABD_WSA_NIR"
  
     albedo_name_bs_sw = "ABD_BSA_shortwave"

     albedo_name_ws_sw = "ABD_WSA_shortwave"
  

  grid_id = gdattach(file_id,grid_name)

  start(1)=0  !hdfeos lib uses 0-based count
  start(2)=0
  edge(1)=nc
  edge(2)=nr
  stride(1)=1
  stride(2)=1
  
  cornerlat(1)=GLASSalbedoobs(source)%gridDesc(4)
  cornerlon(1)=GLASSalbedoobs(source)%gridDesc(5)
  cornerlat(2)=GLASSalbedoobs(source)%gridDesc(7)
  cornerlon(2)=GLASSalbedoobs(source)%gridDesc(8)
  
  if(GLASSalbedoobs(source)%source.eq."AVHRR") then 
     allocate(albedo_raw_avhrr_bs_vis(nc*nr))
     allocate(albedo_raw_avhrr_ws_vis(nc*nr))
     allocate(albedo_raw_avhrr_bs_nir(nc*nr))
     allocate(albedo_raw_avhrr_ws_nir(nc*nr))
     allocate(albedo_raw_avhrr_bs_sw(nc*nr))
     allocate(albedo_raw_avhrr_ws_sw(nc*nr))

     iret1 = gdrdfld(grid_id,albedo_name_bs_vis,start,stride,edge,albedo_raw_avhrr_bs_vis)
     iret2 = gdrdfld(grid_id,albedo_name_ws_vis,start,stride,edge,albedo_raw_avhrr_ws_vis)
     iret3 = gdrdfld(grid_id,albedo_name_bs_nir,start,stride,edge,albedo_raw_avhrr_bs_nir)
     iret4 = gdrdfld(grid_id,albedo_name_ws_nir,start,stride,edge,albedo_raw_avhrr_ws_nir)
     iret5 = gdrdfld(grid_id,albedo_name_bs_sw,start,stride,edge,albedo_raw_avhrr_bs_sw)
     iret6 = gdrdfld(grid_id,albedo_name_ws_sw,start,stride,edge,albedo_raw_avhrr_ws_sw)         
 
     albedo_data_b1 = .false. 
     
     lat_off = nint((cornerlat(1)+89.975)/0.05)+1
     lon_off = nint((cornerlon(1)+179.975)/0.05)+1
     
     do r=1,GLASSalbedoobs(source)%nr
        do c=1,GLASSalbedoobs(source)%nc
           c1 = c + lon_off
           r1 = nr - (r + lat_off) + 1

             if(albedo_raw_avhrr_bs_vis(c1+(r1-1)*nc).gt.0.and.&
                albedo_raw_avhrr_bs_vis(c1+(r1-1)*nc).ne.2550) then 
              albedo_in1(1,c+(r-1)*GLASSalbedoobs(source)%nc) =&
                   albedo_raw_avhrr_bs_vis(c1+(r1-1)*nc)*0.01
              albedo_data_b1(1,c+(r-1)*GLASSalbedoobs(source)%nc) =  .true. 
             else
              albedo_in1(1,c+(r-1)*GLASSalbedoobs(source)%nc) = -9999.0
              albedo_data_b1(1,c+(r-1)*GLASSalbedoobs(source)%nc) = .false. 
             endif
!
             if(albedo_raw_avhrr_ws_vis(c1+(r1-1)*nc).gt.0.and.&
                albedo_raw_avhrr_ws_vis(c1+(r1-1)*nc).ne.2550) then 
              albedo_in1(2,c+(r-1)*GLASSalbedoobs(source)%nc) =&
                   albedo_raw_avhrr_ws_vis(c1+(r1-1)*nc)*0.01
              albedo_data_b1(2,c+(r-1)*GLASSalbedoobs(source)%nc) =  .true. 
             else
              albedo_in1(2,c+(r-1)*GLASSalbedoobs(source)%nc) = -9999.0
              albedo_data_b1(2,c+(r-1)*GLASSalbedoobs(source)%nc) = .false. 
             endif
!
             if(albedo_raw_avhrr_bs_nir(c1+(r1-1)*nc).gt.0.and.&
                albedo_raw_avhrr_bs_nir(c1+(r1-1)*nc).ne.2550) then 
              albedo_in1(3,c+(r-1)*GLASSalbedoobs(source)%nc) =&
                   albedo_raw_avhrr_bs_nir(c1+(r1-1)*nc)*0.01
              albedo_data_b1(3,c+(r-1)*GLASSalbedoobs(source)%nc) =  .true. 
             else
              albedo_in1(3,c+(r-1)*GLASSalbedoobs(source)%nc) = -9999.0
              albedo_data_b1(3,c+(r-1)*GLASSalbedoobs(source)%nc) = .false. 
             endif
!
             if(albedo_raw_avhrr_ws_nir(c1+(r1-1)*nc).gt.0.and.&
                albedo_raw_avhrr_ws_nir(c1+(r1-1)*nc).ne.2550) then 
              albedo_in1(4,c+(r-1)*GLASSalbedoobs(source)%nc) =&
                   albedo_raw_avhrr_ws_nir(c1+(r1-1)*nc)*0.01
              albedo_data_b1(4,c+(r-1)*GLASSalbedoobs(source)%nc) =  .true. 
             else
              albedo_in1(4,c+(r-1)*GLASSalbedoobs(source)%nc) = -9999.0
              albedo_data_b1(4,c+(r-1)*GLASSalbedoobs(source)%nc) = .false. 
             endif
!
             if(albedo_raw_avhrr_bs_sw(c1+(r1-1)*nc).gt.0.and.&
                albedo_raw_avhrr_bs_sw(c1+(r1-1)*nc).ne.2550) then 
              albedo_in1(5,c+(r-1)*GLASSalbedoobs(source)%nc) =&
                   albedo_raw_avhrr_bs_sw(c1+(r1-1)*nc)*0.01
              albedo_data_b1(5,c+(r-1)*GLASSalbedoobs(source)%nc) =  .true. 
             else
              albedo_in1(5,c+(r-1)*GLASSalbedoobs(source)%nc) = -9999.0
              albedo_data_b1(5,c+(r-1)*GLASSalbedoobs(source)%nc) = .false. 
             endif
!
             if(albedo_raw_avhrr_ws_sw(c1+(r1-1)*nc).gt.0.and.&
                albedo_raw_avhrr_ws_sw(c1+(r1-1)*nc).ne.2550) then 
              albedo_in1(6,c+(r-1)*GLASSalbedoobs(source)%nc) =&
                   albedo_raw_avhrr_ws_sw(c1+(r1-1)*nc)*0.01
              albedo_data_b1(6,c+(r-1)*GLASSalbedoobs(source)%nc) =  .true. 
             else
              albedo_in1(6,c+(r-1)*GLASSalbedoobs(source)%nc) = -9999.0
              albedo_data_b1(6,c+(r-1)*GLASSalbedoobs(source)%nc) = .false. 
             endif
!
        enddo
     enddo
     deallocate(albedo_raw_avhrr_bs_vis)
     deallocate(albedo_raw_avhrr_ws_vis)
     deallocate(albedo_raw_avhrr_bs_nir)
     deallocate(albedo_raw_avhrr_ws_nir)
     deallocate(albedo_raw_avhrr_bs_sw)
     deallocate(albedo_raw_avhrr_ws_sw)

  endif
  
  iret=gddetach(grid_id)
  iret=gdclose(file_id)


  if(GLASSalbedoobs(source)%source.eq."AVHRR") then 
    if(LVT_isAtAfinerResolution(GLASSalbedoobs(source)%datares)) then
!--------------------------------------------------------------------------
! Interpolate to the LVT analysis domain
!-------------------------------------------------------------------------- 
      do i=1,6
        call bilinear_interp(LVT_rc%gridDesc,albedo_data_b1(i,:), &
          albedo_in1(i,:),albedoobs_b_ip(i,:),albedoobs_ip(i,:), &
          GLASSalbedoobs(source)%nc*GLASSalbedoobs(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          GLASSalbedoobs(source)%rlat,GLASSalbedoobs(source)%rlon,&
          GLASSalbedoobs(source)%w11,GLASSalbedoobs(source)%w12,&
          GLASSalbedoobs(source)%w21,GLASSalbedoobs(source)%w22,&
          GLASSalbedoobs(source)%n11,GLASSalbedoobs(source)%n12,&
          GLASSalbedoobs(source)%n21,GLASSalbedoobs(source)%n22,LVT_rc%udef,iret)
      enddo
    else
      do i=1,6
        call upscaleByAveraging(GLASSalbedoobs(source)%nc*GLASSalbedoobs(source)%nr,&
          LVT_rc%lnc*LVT_rc%lnr, &
          LVT_rc%udef, GLASSalbedoobs(source)%n11,&
          albedo_data_b1(i,:),albedo_in1(i,:), albedoobs_b_ip(i,:), albedoobs_ip(i,:))
      enddo
    endif
  endif



#endif

end subroutine read_GLASS_ALBEDO_data


!BOP
! !ROUTINE: create_GLASSalbedo_filename
! \label{create_GLASSalbedo_filename}
! 
! !INTERFACE: 
subroutine create_GLASSalbedo_filename(ndir,  yr, doy, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the GLASS ALBEDO filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the GLASS ALBEDO data directory
!  \item[yr]  current year
!  \item[mo]  current doy
!  \item[filename] Generated GLASS ALBEDO filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  filename = trim(ndir)//'/'//trim(fyr)//'/GLASS02B05.V04.A'//&
       trim(fyr)//trim(fdoy)//'.hdf'

end subroutine create_GLASSalbedo_filename


