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
! !MODULE: SNODEP_obsMod
!  \label(SNODEP_obsMod)
!
! !INTERFACE:
module SNODEP_obsMod
! 
! !USES: 
  use ESMF

  implicit none
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
!  23 Apr 2010   Sujay Kumar  Initial Specification
! 
!EOP
!
! 
!
! 
!

  PRIVATE 

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SNODEP_obsinit !Initializes structures for reading SNODEP data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SNODEPobs !Object to hold SNODEP observation attributes
!EOP
  type, public :: snodepobsdec
     real, allocatable      :: rlat1_nh(:)
     real, allocatable      :: rlat1_sh(:)
     real, allocatable      :: rlon1_nh(:)
     real, allocatable      :: rlon1_sh(:)
     integer, allocatable   :: n111_nh(:)
     integer, allocatable   :: n111_sh(:)
     integer, allocatable   :: n121_nh(:)
     integer, allocatable   :: n121_sh(:)
     integer, allocatable   :: n211_nh(:)
     integer, allocatable   :: n211_sh(:)
     integer, allocatable   :: n221_nh(:)
     integer, allocatable   :: n221_sh(:)
     real, allocatable      :: w111_nh(:)
     real, allocatable      :: w121_nh(:)
     real, allocatable      :: w111_sh(:)
     real, allocatable      :: w121_sh(:)
     real, allocatable      :: w211_nh(:)
     real, allocatable      :: w221_nh(:)
     real, allocatable      :: w211_sh(:)
     real, allocatable      :: w221_sh(:)     
     integer                :: gridspan
     integer                :: shemi
     integer                :: nhemi
     integer                :: mo1
     integer                :: mo2
     integer                :: hemi_nc(2)
     integer                :: hemi_nr(2)
     integer                :: mi
     integer                :: pmax
     character*100               :: odir
     integer                     :: mesh
     logical                     :: startFlag
  end type snodepobsdec

  type(snodepobsdec), allocatable :: snodepobs(:)

contains
  
!BOP
! 
! !ROUTINE: SNODEP_obsInit
! \label{SNODEP_obsInit}
!
! !INTERFACE: 
  subroutine SNODEP_obsinit(source)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading SNODEP data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    
    integer            :: ts
    integer            :: status, rc
    real               :: gridDesci(50)
    real               :: gridDesco(50)
    integer            :: ihemi
    real, parameter :: xmeshl1 = 47.625
    real, parameter :: xmeshl2 = 23.812
    real, parameter :: xpnmcaf1 = 257
    real, parameter :: ypnmcaf1 = 257
    real, parameter :: xpnmcaf2 = 513
    real, parameter :: ypnmcaf2 = 513
    real :: xmesh, orient,xi1,xj1
    real :: alat1,alon1
    
    if(.not.allocated(snodepobs)) then 
       allocate(snodepobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, snodepobs(source)%odir, &
         label='SNODEP observation directory:',rc=status)
    call LVT_verify(status, 'SNODEP observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, snodepobs(source)%mesh, &
         label='SNODEP mesh resolution:',rc=status)
    call LVT_verify(status, 'SNODEP mesh resolution: not defined')

    if(snodepobs(source)%mesh.eq.8) then 
       snodepobs(source)%pmax = 512
       snodepobs(source)%mi = 512*512
    elseif(snodepobs(source)%mesh.eq.16) then 
       snodepobs(source)%pmax = 1024
       snodepobs(source)%mi = 1024*1024
    endif

    snodepobs(source)%startFlag = .true. 
    ts = 86400
    call LVT_update_timestep(LVT_rc, 86400)

    if(LVT_rc%gridDesc(1).eq.0) then 
       if(LVT_rc%gridDesc(4).ge.0.and.LVT_rc%gridDesc(7).ge.0) then 
          snodepobs(source)%gridspan = 1
          snodepobs(source)%shemi = 1
          snodepobs(source)%nhemi = 1
       elseif(LVT_rc%gridDesc(4).le.0.and.LVT_rc%gridDesc(7).le.0) then 
          snodepobs(source)%gridspan = 2
          snodepobs(source)%shemi = 2
          snodepobs(source)%nhemi = 2
       else
          snodepobs(source)%gridspan = 3
          snodepobs(source)%shemi = 1
          snodepobs(source)%nhemi = 2
       endif
    endif

    if(LVT_rc%gridDesc(1) .eq.0) then !latlon domain
       if(LVT_rc%gridDesc(4).ge.0.and.LVT_rc%gridDesc(7).ge.0) then 
          snodepobs(source)%gridspan = 1 ! NH only
          snodepobs(source)%shemi = 1
          snodepobs(source)%nhemi = 1
       elseif(LVT_rc%gridDesc(4).le.0.and.LVT_rc%gridDesc(7).le.0) then 
          snodepobs(source)%gridspan = 2 ! SH only 
          snodepobs(source)%shemi = 2
          snodepobs(source)%nhemi = 2
       else
          snodepobs(source)%gridspan = 3 ! NH and SH
          snodepobs(source)%shemi = 1
          snodepobs(source)%nhemi = 2
       endif
       snodepobs(source)%hemi_nc = nint((LVT_rc%gridDesc(8)-&
            LVT_rc%gridDesc(5))&
            /LVT_rc%gridDesc(9))+1
       if(snodepobs(source)%gridspan.eq.1) then 
          snodepobs(source)%hemi_nr(1) = nint((LVT_rc%gridDesc(7)-&
               LVT_rc%gridDesc(4))/LVT_rc%gridDesc(10))+1
          snodepobs(source)%hemi_nr(2) = 0 
          snodepobs(source)%mo1 = snodepobs(source)%hemi_nc(1)*&
               snodepobs(source)%hemi_nr(1)
          snodepobs(source)%mo2 = 0 
       elseif(snodepobs(source)%gridspan.eq.2) then 
          snodepobs(source)%hemi_nr(1) = 0 
          snodepobs(source)%hemi_nr(2) = nint((LVT_rc%gridDesc(7)-&
               LVT_rc%gridDesc(4))/LVT_rc%gridDesc(10)+1)
          snodepobs(source)%mo1 = 0 
          snodepobs(source)%mo2 = snodepobs(source)%hemi_nc(2)*&
               snodepobs(source)%hemi_nr(2)
       else
          snodepobs(source)%hemi_nr(1) = nint((LVT_rc%gridDesc(7)-&
               LVT_rc%gridDesc(10)/2)/LVT_rc%gridDesc(10)+1)
          snodepobs(source)%hemi_nr(2) = nint((-LVT_rc%gridDesc(10)/2-&
               LVT_rc%gridDesc(4))/LVT_rc%gridDesc(10)+1) 
          snodepobs(source)%mo1 = snodepobs(source)%hemi_nc(1)*&
               snodepobs(source)%hemi_nr(1)
          snodepobs(source)%mo2 = snodepobs(source)%hemi_nc(2)*&
               snodepobs(source)%hemi_nr(2)
       endif
       
       gridDesco = 0 
       gridDesci = 0 
       do ihemi = snodepobs(source)%shemi,snodepobs(source)%nhemi
          gridDesco(1) = 0 
          gridDesco(2) = snodepobs(source)%hemi_nc(ihemi)
          gridDesco(3) = snodepobs(source)%hemi_nr(ihemi)
          gridDesco(5) = LVT_rc%gridDesc(5)
          gridDesco(8) = LVT_rc%gridDesc(8)
          gridDesco(6) = LVT_rc%gridDesc(6)
          gridDesco(9) = LVT_rc%gridDesc(9)
          gridDesco(10) = LVT_rc%gridDesc(10)
          gridDesco(20) = 255
          if(snodepobs(source)%gridspan.eq.1.or.&
               snodepobs(source)%gridspan.eq.2) then 
             gridDesco(4) = LVT_rc%gridDesc(4)
             gridDesco(7) = LVT_rc%gridDesc(7)
          else
             if(ihemi.eq.1) then 
                gridDesco(4) = LVT_rc%gridDesc(9)/2
                gridDesco(7) = LVT_rc%gridDesc(7)
             else
                gridDesco(4) = LVT_rc%gridDesc(4)
                gridDesco(7) = -LVT_rc%gridDesc(9)/2
             endif
          endif
          
          if(snodepobs(source)%mesh.eq.8) then 
             if(ihemi.eq.1) then 
                xmesh = xmeshl1
                orient = 100.0
             else
                xmesh = -xmeshl1
                orient = 280.0
             endif
             xj1 = float(1)-ypnmcaf1
             xi1 = float(1)-xpnmcaf1
             
             call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
          elseif(snodepobs(source)%mesh.eq.16) then 
             if(ihemi.eq.1) then 
                xmesh = xmeshl2
                orient = 100.0
             else
                xmesh = -xmeshl2
                orient = 280.0
             endif
             xj1 = float(1)-ypnmcaf2
             xi1 = float(1)-xpnmcaf2
             
             call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
          endif
          
          gridDesci = 0 
          gridDesci(1) = 5
          gridDesci(2) = snodepobs(source)%pmax
          gridDesci(3) = snodepobs(source)%pmax
          gridDesci(4) = alat1
          gridDesci(5) = alon1
          gridDesci(6) = 8
          gridDesci(7) = orient
          gridDesci(8) = xmesh
          gridDesci(9) = xmesh
          gridDesci(10) = 0.0       
          if(ihemi .eq.2) then 
             gridDesci(20) = 128
             gridDesci(11) = 128
          endif
          gridDesci(13) = 1  !global grid
          gridDesci(20) = 0 
          
          if(snodepobs(source)%gridspan.eq.1) then 
             allocate(snodepobs(source)%rlat1_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%rlon1_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%n111_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%n121_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%n211_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%n221_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%w111_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%w121_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%w211_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%w221_nh(snodepobs(source)%mo1))
             call bilinear_interp_input(gridDesci,gridDesco,&
                  snodepobs(source)%mo1,snodepobs(source)%rlat1_nh,&
                  snodepobs(source)%rlon1_nh,snodepobs(source)%n111_nh,&
                  snodepobs(source)%n121_nh,snodepobs(source)%n211_nh,&
                  snodepobs(source)%n221_nh,snodepobs(source)%w111_nh,&
                  snodepobs(source)%w121_nh,snodepobs(source)%w211_nh,&
                  snodepobs(source)%w221_nh)
          elseif(snodepobs(source)%gridspan.eq.2) then 
             allocate(snodepobs(source)%rlat1_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%rlon1_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%n111_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%n121_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%n211_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%n221_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%w111_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%w121_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%w211_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%w221_sh(snodepobs(source)%mo2))
             call bilinear_interp_input(gridDesci,gridDesco,&
                  snodepobs(source)%mo2,snodepobs(source)%rlat1_sh,&
                  snodepobs(source)%rlon1_sh,&
                  snodepobs(source)%n111_sh,snodepobs(source)%n121_sh,&
                  snodepobs(source)%n211_sh,snodepobs(source)%n221_sh,&
                  snodepobs(source)%w111_sh,snodepobs(source)%w121_sh,&
                  snodepobs(source)%w211_sh,snodepobs(source)%w221_sh)
          else
             if(ihemi.eq.1) then 
                allocate(snodepobs(source)%rlat1_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%rlon1_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%n111_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%n121_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%n211_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%n221_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%w111_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%w121_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%w211_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%w221_nh(snodepobs(source)%mo1))
                call bilinear_interp_input(gridDesci,gridDesco,&
                     snodepobs(source)%mo1,snodepobs(source)%rlat1_nh,&
                     snodepobs(source)%rlon1_nh,snodepobs(source)%n111_nh,&
                     snodepobs(source)%n121_nh,snodepobs(source)%n211_nh,&
                     snodepobs(source)%n221_nh,snodepobs(source)%w111_nh,&
                     snodepobs(source)%w121_nh,snodepobs(source)%w211_nh,&
                     snodepobs(source)%w221_nh)
             elseif(ihemi.eq.2) then 
                allocate(snodepobs(source)%rlat1_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%rlon1_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%n111_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%n121_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%n211_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%n221_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%w111_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%w121_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%w211_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%w221_sh(snodepobs(source)%mo2))
                call bilinear_interp_input(gridDesci,gridDesco,&
                     snodepobs(source)%mo2,snodepobs(source)%rlat1_sh,&
                     snodepobs(source)%rlon1_sh,snodepobs(source)%n111_sh,&
                     snodepobs(source)%n121_sh,snodepobs(source)%n211_sh,&
                     snodepobs(source)%n221_sh,snodepobs(source)%w111_sh,&
                     snodepobs(source)%w121_sh,snodepobs(source)%w211_sh,&
                     snodepobs(source)%w221_sh)
             endif
          endif
       enddo
    elseif(LVT_rc%gridDesc(1).ne.0) then 
       if(LVT_rc%gridDesc(4).ge.0.and.LVT_rc%gridDesc(10).ge.0) then 
          snodepobs(source)%gridspan = 1 ! NH only
          snodepobs(source)%shemi = 1
          snodepobs(source)%nhemi = 1
       elseif(LVT_rc%gridDesc(4).le.0.and.LVT_rc%gridDesc(10).ge.0) then 
          snodepobs(source)%gridspan = 2 ! SH only 
          snodepobs(source)%shemi = 2
          snodepobs(source)%nhemi = 2
       elseif(LVT_rc%gridDesc(10).eq.-100.and.LVT_rc%gridDesc(20).eq.-100) then 
!-----------------------------------------------------------------------------
! Global grid in polar stereographic projection. No interpolation 
! will be done.
!-----------------------------------------------------------------------------
          snodepobs(source)%gridspan = 3
          snodepobs(source)%shemi = 1
          snodepobs(source)%nhemi = 2
       else
          write(*,*) 'Currently spanning across hemispheres is'
          write(*,*) 'not supported for SNODEP data'
          call LVT_endrun()
       endif
       snodepobs(source)%hemi_nc = LVT_rc%gridDesc(2)
       if(snodepobs(source)%gridspan.eq.1) then 
          snodepobs(source)%hemi_nr(1) = LVT_rc%gridDesc(3)
          snodepobs(source)%hemi_nr(2) = 0 
          snodepobs(source)%mo1 = snodepobs(source)%hemi_nc(1)*&
               snodepobs(source)%hemi_nr(1)
          snodepobs(source)%mo2 = 0 
       elseif(snodepobs(source)%gridspan.eq.2) then 
          snodepobs(source)%hemi_nr(1) = 0 
          snodepobs(source)%hemi_nr(2) = LVT_rc%gridDesc(3)
          snodepobs(source)%mo1 = 0 
          snodepobs(source)%mo2 = snodepobs(source)%hemi_nc(2)*&
               snodepobs(source)%hemi_nr(2)
       else
          snodepobs(source)%hemi_nr(1) = LVT_rc%gridDesc(3)
          snodepobs(source)%hemi_nr(2) = LVT_rc%gridDesc(3)
          snodepobs(source)%mo1 = snodepobs(source)%hemi_nc(1)*&
               snodepobs(source)%hemi_nr(1)
          snodepobs(source)%mo2 = snodepobs(source)%hemi_nc(2)*&
               snodepobs(source)%hemi_nr(2)             
       endif
       gridDesco = 0 
       gridDesci = 0 
       do ihemi = snodepobs(source)%shemi,snodepobs(source)%nhemi
          gridDesco(1) = LVT_rc%gridDesc(1)
          gridDesco(2) = snodepobs(source)%hemi_nc(ihemi)
          gridDesco(3) = snodepobs(source)%hemi_nr(ihemi)
          gridDesco(5) = LVT_rc%gridDesc(5)
          gridDesco(8) = LVT_rc%gridDesc(8)
          gridDesco(6) = LVT_rc%gridDesc(6)
          gridDesco(9) = LVT_rc%gridDesc(9)
          gridDesco(10) = LVT_rc%gridDesc(10)
          gridDesco(11) = LVT_rc%gridDesc(11)
          gridDesco(20) = 255
          if(snodepobs(source)%gridspan.eq.1.or.&
               snodepobs(source)%gridspan.eq.2) then 
             gridDesco(4) = LVT_rc%gridDesc(4)
             gridDesco(7) = LVT_rc%gridDesc(7)
          endif
          if(snodepobs(source)%mesh.eq.8) then 
             if(ihemi.eq.1) then 
                xmesh = xmeshl1
                orient = 100.0
             else
                xmesh = -xmeshl1
                orient = 280.0
             endif
             xj1 = float(1)-ypnmcaf1
             xi1 = float(1)-xpnmcaf1
             
             call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
          elseif(snodepobs(source)%mesh.eq.16) then 
             if(ihemi.eq.1) then 
                xmesh = xmeshl2
                orient = 100.0
             else
                xmesh = -xmeshl2
                orient = 280.0
             endif
             xj1 = float(1)-ypnmcaf2
             xi1 = float(1)-xpnmcaf2
             
             call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
          endif
          
          gridDesci = 0 
          gridDesci(1) = 5
          gridDesci(2) = snodepobs(source)%pmax
          gridDesci(3) = snodepobs(source)%pmax
          gridDesci(4) = alat1
          gridDesci(5) = alon1
          gridDesci(6) = 8
          gridDesci(7) = orient
          gridDesci(8) = xmesh
          gridDesci(9) = xmesh
          gridDesci(10) = 0.0
          if(ihemi .eq.2) then 
             gridDesci(20) = 128
             gridDesci(11) = 128
          endif
          
          gridDesci(20) = 0 
          gridDesci(13) = 1  !global grid
          
          if(snodepobs(source)%gridspan.eq.1) then 
             allocate(snodepobs(source)%rlat1_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%rlon1_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%n111_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%n121_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%n211_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%n221_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%w111_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%w121_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%w211_nh(snodepobs(source)%mo1))
             allocate(snodepobs(source)%w221_nh(snodepobs(source)%mo1))
             call bilinear_interp_input(gridDesci,gridDesco,&
                  snodepobs(source)%mo1,snodepobs(source)%rlat1_nh,&
                  snodepobs(source)%rlon1_nh,snodepobs(source)%n111_nh,&
                  snodepobs(source)%n121_nh,snodepobs(source)%n211_nh,&
                  snodepobs(source)%n221_nh,snodepobs(source)%w111_nh,&
                  snodepobs(source)%w121_nh,snodepobs(source)%w211_nh,&
                  snodepobs(source)%w221_nh)
          elseif(snodepobs(source)%gridspan.eq.2) then 
             
             allocate(snodepobs(source)%rlat1_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%rlon1_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%n111_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%n121_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%n211_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%n221_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%w111_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%w121_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%w211_sh(snodepobs(source)%mo2))
             allocate(snodepobs(source)%w221_sh(snodepobs(source)%mo2))
             call bilinear_interp_input(gridDesci,gridDesco,&
                  snodepobs(source)%mo2,snodepobs(source)%rlat1_sh,&
                  snodepobs(source)%rlon1_sh,&
                  snodepobs(source)%n111_sh,snodepobs(source)%n121_sh,&
                  snodepobs(source)%n211_sh,snodepobs(source)%n221_sh,&
                  snodepobs(source)%w111_sh,snodepobs(source)%w121_sh,&
                  snodepobs(source)%w211_sh,snodepobs(source)%w221_sh)
          elseif(snodepobs(source)%gridspan.eq.3) then 
             !-------------------------------------------------------------------------
             !   No interpolation is being done. So no weights will be computed. 
             !-------------------------------------------------------------------------
             if(ihemi.eq.1) then 
                allocate(snodepobs(source)%rlat1_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%rlon1_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%n111_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%n121_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%n211_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%n221_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%w111_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%w121_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%w211_nh(snodepobs(source)%mo1))
                allocate(snodepobs(source)%w221_nh(snodepobs(source)%mo1))
             elseif(ihemi.eq.2) then 
                allocate(snodepobs(source)%rlat1_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%rlon1_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%n111_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%n121_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%n211_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%n221_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%w111_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%w121_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%w211_sh(snodepobs(source)%mo2))
                allocate(snodepobs(source)%w221_sh(snodepobs(source)%mo2))
             endif
          endif
       enddo
    endif

  end subroutine SNODEP_obsinit


end module SNODEP_obsMod
