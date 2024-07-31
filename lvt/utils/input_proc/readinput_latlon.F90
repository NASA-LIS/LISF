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
! !ROUTINE: readinput_latlon
!  \label{readinput_latlon}
!  
! !REVISION HISTORY:
! 
!  11 Apr 2000: Brian Cosgrove; Added Elevation correction and Forcing
!               Mask read statements
!  14 Nov 2003: Sujay Kumar; Modified card file that includes regional 
!               modeling options
! !INTERFACE:
subroutine readinput_latlon

! !USES:
  use ESMF 
  use preprocMod
  use map_utils

! !DESCRIPTION: 
!
!EOP
  implicit none

  integer              :: i, j,count
  integer              :: k,n
  real, allocatable    :: run_dd(:,:)  ! LIS Target grid/domain
  integer              :: lnc,lnr
  integer              :: nc, nr
  integer              :: ierr
  integer              :: rc
  integer              :: ips, ipe, jps, jpe
  integer              :: Px, Py, P
  integer              :: mytask_x, mytask_y
  real(ESMF_KIND_R8)   :: stlat, stlon, dx, dy
  integer              :: status
! __________________________________________________________________________

  allocate(run_dd(LVT_rc%nnest,20))

  call ESMF_ConfigFindLabel(LVT_config,"run domain lower left lat:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,1),rc=rc)
     call LVT_verify(rc, 'run domain lower left lat: not defined')
  enddo

  call ESMF_ConfigFindLabel(LVT_config,"run domain lower left lon:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,2),rc=rc)
     call LVT_verify(rc, 'run domain lower left lon: not defined')
  enddo
  call ESMF_ConfigFindLabel(LVT_config,"run domain upper right lat:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,3),rc=rc)
     call LVT_verify(rc, 'run domain upper right lat: not defined')
  enddo

  call ESMF_ConfigFindLabel(LVT_config,"run domain upper right lon:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,4),rc=rc)
     call LVT_verify(rc, 'run domain upper right lon: not defined')
  enddo

  call ESMF_ConfigFindLabel(LVT_config,"run domain resolution (dx):",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,5),rc=rc)
     call LVT_verify(rc, 'run domain resolution (dx): not defined')
  enddo

  call ESMF_ConfigFindLabel(LVT_config,"run domain resolution (dy):",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,6),rc=rc)
     call LVT_verify(rc, 'run domain resolution (dy): not defined')
  enddo

  do n=1,LVT_rc%nnest
     stlat = run_dd(n,1)
     stlon = run_dd(n,2)
     dx    = run_dd(n,5)
     dy    = run_dd(n,6)
     nc    = (nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,5))) + 1
     nr    = (nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,6))) + 1
     
     LVT_rc%gridDesc(n,4) = run_dd(n,1)
     LVT_rc%gridDesc(n,5) = run_dd(n,2)
     LVT_rc%gridDesc(n,7) = run_dd(n,3)
     LVT_rc%gridDesc(n,8) = run_dd(n,4)


     LVT_rc%gnc(n) = nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,5))+1
     LVT_rc%gnr(n) = nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,6))+1

     LVT_rc%lnc(n) = LVT_rc%gnc(n)
     LVT_rc%lnr(n) = LVT_rc%gnr(n)

     LVT_rc%gridDesc(n,1) = 0
     LVT_rc%gridDesc(n,9) = run_dd(n,5)
     
     if(LVT_rc%gridDesc(n,1).eq.0) then 
        LVT_rc%gridDesc(n,10) = run_dd(n,6)
        LVT_rc%gridDesc(n,6)  = 128
        LVT_rc%gridDesc(n,11) = 64
        LVT_rc%gridDesc(n,20) = 64
     endif
     if(LVT_rc%gridDesc(n,7).lt.LVT_rc%gridDesc(n,4)) then
        write(unit=*,fmt=*) 'lat2 must be greater than lat1'
        write(*,*) LVT_rc%gridDesc(n,7),LVT_rc%gridDesc(n,4)
        write(unit=*,fmt=*) 'Stopping run...'
        stop
     endif
     if(LVT_rc%gridDesc(n,8).lt.LVT_rc%gridDesc(n,5)) then
        write(unit=*,fmt=*) 'lon2 must be greater than lon1'
        write(*,*) LVT_rc%gridDesc(n,8),LVT_rc%gridDesc(n,5)
        write(unit=*,fmt=*) 'Stopping run...'
        stop
     endif
     
     LVT_rc%gridDesc(n,2) = nint((LVT_rc%gridDesc(n,8)-LVT_rc%gridDesc(n,5))&
          /LVT_rc%gridDesc(n,10))+ 1
     LVT_rc%gridDesc(n,3) = nint((LVT_rc%gridDesc(n,7)-LVT_rc%gridDesc(n,4))&
          /LVT_rc%gridDesc(n,9)) + 1  

     write(*,*)'-------------------- LVT/LIS Domain ----------------------'
     do k=1,10
        write(unit=*,fmt=*) '(',k,',',LVT_rc%gridDesc(n,k),')'
     enddo
     
     call map_set( PROJ_LATLON, LVT_rc%gridDesc(n,4),LVT_rc%gridDesc(n,5), &
                   0.0, LVT_rc%gridDesc(n,9),LVT_rc%gridDesc(n,10), 0.0,   &
                   LVT_rc%lnc(n),LVT_rc%lnr(n),LVT_domain(n)%lvtproj )
     
     write(*,*)'----------------------------------------------------------'

  enddo  ! end nest loop
  
  deallocate(run_dd)

end subroutine readinput_latlon
      
