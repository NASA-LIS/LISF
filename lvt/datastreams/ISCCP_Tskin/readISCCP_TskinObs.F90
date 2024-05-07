!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: readISCCP_Tskinobs
! \label(readISCCP_Tskinobs)
!
! !INTERFACE:
subroutine readISCCP_Tskinobs(source)
! !USES:   
  use LVT_coreMod
  use LVT_historyMod
  use LVT_logMod
  use LVT_histDataMod
  use ISCCP_TskinobsMod, only : ISCCP_Tskin_obs

  implicit none
!
! !INPUT PARAMETERS: 
  integer,    intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  real, parameter     :: tskin_min = 200.0
  real, parameter     :: tskin_max = 400.0
  integer, parameter  :: N_gswp2_compressed = 15238
  integer, parameter  :: N_gswp2            = 64800

  ! land_i_gswp2 and land_j_gswp2 as stored in 
  ! ISCCP_Tskin_GSWP2_grid_V1 files (by Sarith) follow the GSWP2 convention
  ! for grid orientation, that is counting from north-to-south
  ! and from west-to-east

  real, parameter     :: minlon_gswp2 = -180.5
  real, parameter     :: maxlat_gswp2 =   90.5
  
  real, parameter     :: dx_gswp2 = 1.
  real, parameter     :: dy_gswp2 = 1.
  
  integer             :: land_i_gswp2, land_j_gswp2
  real                :: tsclr
!  real, dimension(N_gswp2_compressed) :: tmp_obs, tmp_lat, tmp_lon
  character*100       :: name
  real, dimension(N_gswp2) :: tmp_obs
  integer             :: c,r,j,i,istat
  real                :: tskin(LVT_rc%lnc, LVT_rc%lnr)
  integer             :: ftn
  integer             :: index1, iret
  logical*1           :: li(N_gswp2)
  logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                :: go(LVT_rc%lnc*LVT_rc%lnr)
  logical             :: readflag, file_exists


  call ISCCP_Tskin_filename(name, ISCCP_Tskin_obs(source)%odir, &
       LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), LVT_rc%dhr(source))
  
  inquire(file=name,exist=file_exists)
  
  if(file_exists.and.LVT_rc%dmn(source).eq.0) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif
  
  tskin = LVT_rc%udef 

  if (readflag) then 
     write(LVT_logunit,*)  '[INFO] Reading ISCCP Tskin data ',name
     
     li = .false. 
     tmp_obs = -9999.0

     ftn = LVT_getNextUnitNumber()
     open(ftn,file=trim(name), form='unformatted', status='old', &
          iostat = istat)

     if(istat.eq.0) then 
        j = 0 
        do i=1,N_gswp2_compressed
           read(ftn) land_i_gswp2, land_j_gswp2, tsclr
           
           if((tsclr > tskin_min ) .and. &
                (tsclr < tskin_max)) then 
              j = j + 1
              index1 = land_i_gswp2 + (180 - land_j_gswp2)*360
              tmp_obs(index1) = tsclr
              li(index1) = .true. 
           endif
        enddo

        call neighbor_interp(LVT_rc%gridDesc,  li, tmp_obs, &
             lo, go, &
             isccp_tskin_obs(source)%mi, LVT_rc%lnc*LVT_rc%lnr, &
             isccp_tskin_obs(source)%rlat, isccp_tskin_obs(source)%rlon, &
             isccp_tskin_obs(source)%n11, LVT_rc%udef, iret)

!        open(100,file='isccp.bin',form='unformatted')
!        write(100) tmp_obs
!        close(100)
!        stop

        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              tskin(c,r) = go(c+LVT_rc%lnc*(r-1))
           enddo
        enddo

     endif
     
     call LVT_releaseUnitNumber(ftn)
  endif  

  call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT,source,tskin,vlevel=1,units="K")

end subroutine readISCCP_Tskinobs


!BOP
! 
! !ROUTINE: ISCCP_Tskin_filename
! \label(ISCCP_Tskin_filename)
!
! !INTERFACE:
subroutine ISCCP_Tskin_filename(name, ndir, yr, mo,da,hr)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  
  implicit none
  character*80      :: name
  integer           :: yr, mo, da, hr
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  
  name = trim(ndir)//'/Y'//trim(fyr)//'/M'//trim(fmo)//'/isccpdx_tskin.'//trim(fyr)//trim(fmo)//trim(fda)//'_'//trim(fhr)//'z.bin'

end subroutine ISCCP_Tskin_filename

