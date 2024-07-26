!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: readinput_easev2
!  \label{readinput_easev2}
!  
! !REVISION HISTORY:
! 
!  14 Sept 2016: Sujay Kumar; Initial specification
!               
! !INTERFACE:
subroutine readinput_easev2(nest)

! !USES:
  use ESMF 
  use LDT_coreMod
  use LDT_domainMod
  use LDT_logMod
  use LDT_fileIOMod
  use map_utils
  use easeV2_utils

  implicit none

  integer, intent(in) :: nest

! !DESCRIPTION: 
!
!  This routine reads the options specified in the LIS configuration 
!  file to determine the region of interest in the LIS simulation, in 
!  an EASE grid configuration. The routine reads the extents of the 
!  running domain. The land surface parameters used for the simulation
!  are read from a file that spans the area of interest. 
!  Based on the number of processors and the processor layout specified, 
!  this routine also performs the domain decomposition. 
!
! !NOTES: The code is not currently tested for parallel domain decomposition
!
!  The routines invoked are: 
!  \begin{description}
!  \item[listask\_for\_point](\ref{LDT_mpDecomp}) \newline
!   routine to perform domain decomposition
!  \end{description}
!EOP

  integer              :: i, j, k, n, count
  real, allocatable    :: run_dd(:,:)  ! LIS Target grid/domain
  integer              :: lnc,lnr
  integer              :: nc, nr
  integer              :: ierr, rc
  integer              :: ips, ipe, jps, jpe
  integer              :: Px, Py, P
  integer              :: mytask_x, mytask_y
  character*40         :: run_dd_grid_desgn
  real                 :: stlat, stlon, enlat,enlon, dx, dy
  integer              :: status
  real                 :: c1, r1, c2, r2
  integer              :: c,r
  real                 :: lat1,lon1
! __________________________________________________________________________

  allocate(run_dd(LDT_rc%nnest,20))

  LDT_rc%lis_map_resfactor(nest) = 1.

  call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lat:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,1),rc=rc)
     call LDT_verify(rc, 'Run domain lower left lat: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lon:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,2),rc=rc)
     call LDT_verify(rc, 'Run domain lower left lon: not defined')
  enddo
  call ESMF_ConfigFindLabel(LDT_config,"Run domain upper right lat:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,3),rc=rc)
     call LDT_verify(rc, 'Run domain upper right lat: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain upper right lon:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,4),rc=rc)
     call LDT_verify(rc, 'Run domain upper right lon: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain grid designation:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd_grid_desgn,rc=rc)
     call LDT_verify(rc, 'Run domain grid designation: not defined')
  enddo

  do n=1,nest
     stlat = run_dd(n,1)
     stlon = run_dd(n,2)
     enlat = run_dd(n,3)
     enlon = run_dd(n,4)

     if(run_dd_grid_desgn.eq."M09") then 
!        nc    = 3856
!        nr    = 1624
        call easeV2_convert ("M09", stlat, stlon, c1,r1)
        call easeV2_convert ("M09", enlat, enlon, c2,r2)
        
!        c1 = real(nint(c1))
!        r1 = real(nint(r1))
!        c2 = real(nint(c2))
!        r2 = real(nint(r2))
!        call easeV2_inverse("M36", c1,r1,stlat, stlon)
!        call easeV2_inverse("M36", c2,r2,enlat, enlon)

        nc = nint(c2)-nint(c1)+1
        nr = r1-r2+1

        LDT_rc%gridDesc(n,9) = 5
        LDT_rc%gridDesc(n,10) = 0.09
     elseif(run_dd_grid_desgn.eq."M36") then 

        call easeV2_convert ("M36", stlat, stlon, c1,r1)
        call easeV2_convert ("M36", enlat, enlon, c2,r2)
        
!        c1 = real(nint(c1))
!        r1 = real(nint(r1))
!        c2 = real(nint(c2))
!        r2 = real(nint(r2))
!        call easeV2_inverse("M36", c1,r1,stlat, stlon)
!        call easeV2_inverse("M36", c2,r2,enlat, enlon)

        nc = nint(c2)-nint(c1)+1
        nr = r1-r2+1

!        call easeV2_inverse("M36", 1.0,1.0,stlat, stlon)
!        call easeV2_inverse("M36", 964.0,406.0,enlat, enlon)
!        print*, stlat, stlon, enlat, enlon
!        do c=1,964
!           call easeV2_inverse("M36", real(c),1.0,stlat, stlon)
!           print*, c,stlat, stlon
!        enddo

!        stop
!        nc    = 964
!        nr    = 406
        LDT_rc%gridDesc(n,9) = 4
        LDT_rc%gridDesc(n,10) = 0.36
     else
        write(LDT_logunit,*) '[ERR] target grid designation '//&
             trim(run_dd_grid_desgn)//' not supported'
        call LDT_endrun()
     endif

     call LDT_quilt_domain(n,nc,nr)

     LDT_rc%gridDesc(n,4) = stlat
     LDT_rc%gridDesc(n,5) = stlon
     LDT_rc%gridDesc(n,7) = enlat
     LDT_rc%gridDesc(n,8) = enlon
     
     write(unit=LDT_logunit,fmt=*) 'local domain ',&
          LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,7),&
          LDT_rc%gridDesc(n,5),LDT_rc%gridDesc(n,8)
     
     LDT_rc%gnc(n) = nc
     LDT_rc%gnr(n) = nr
     LDT_rc%lnc(n) = nc
     LDT_rc%lnr(n) = nr
     LDT_rc%gridDesc(n,1) = 9     
     LDT_rc%gridDesc(n,20) = 64

     
     LDT_rc%gridDesc(n,2) = nc
     LDT_rc%gridDesc(n,3) = nr

     write(LDT_logunit,*)'-------------------- LDT/LIS Domain ----------------------'
     do k=1,13
        write(unit=LDT_logunit,fmt=*) '(',k,',',LDT_rc%gridDesc(n,k),')'
     enddo
     
   ! local domain of each processor
     call map_set( PROJ_EASEV2, LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,5), &
                   LDT_rc%gridDesc(n,9), LDT_rc%gridDesc(n,9),&
                   LDT_rc%gridDesc(n,9), 0.0,   &
                   LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_domain(n)%ldtproj )


   ! global domain
     call map_set( PROJ_EASEV2,LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,5), &
                   LDT_rc%gridDesc(n,9), LDT_rc%gridDesc(n,9),&
                   LDT_rc%gridDesc(n,9), 0.0,   &
                   LDT_rc%gnc(n),LDT_rc%gnr(n),LDT_domain(n)%ldtglbproj )
!     
     write(LDT_logunit,*)'----------------------------------------------------------'

  enddo  ! end nest loop
  
  deallocate(run_dd)

end subroutine readinput_easev2
      
