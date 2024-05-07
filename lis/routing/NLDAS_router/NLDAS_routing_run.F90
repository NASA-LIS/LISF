!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine NLDAS_routing_run(n)

  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_masterproc, LIS_localPet
  use LIS_timeMgrMod,   only : LIS_isAlarmRinging
  use LIS_routingMod,   only : LIS_runoff_state
  use LIS_logMod,       only : LIS_logunit, LIS_verify
  use LIS_historyMod,   only : LIS_grid2patch, LIS_tile2grid
  use LIS_histDataMod
  use NLDAS_routingMod, only : NLDAS_routing_struc
  
  implicit none
  
  integer, intent(in)   :: n  
  type(ESMF_Field)      :: sf_runoff_field
  type(ESMF_Field)      :: sf_runoff_count_field
  type(ESMF_Field)      :: baseflow_field
  real,   allocatable       :: surface_runoff(:,:)
  real,   allocatable       :: baseflow(:,:)
  real,   pointer       :: surface_runoff_t(:)
  real,   pointer       :: baseflow_t(:)
  real,   allocatable   :: streamflow_lvec(:)
  logical               :: alarmCheck
  integer               :: status
  integer               :: c,r,t

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "NLDAS router model alarm")
  if(alarmCheck) then 
     allocate(streamflow_lvec(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     if(LIS_masterproc) then
!run the routing model at 1 hour output interval
     
        call ESMF_StateGet(LIS_runoff_state(n),"Surface Runoff",sf_runoff_field,&
             rc=status)
        call LIS_verify(status, "ESMF_StateGet failed for Surface Runoff")
        
        call ESMF_StateGet(LIS_runoff_state(n),"Subsurface Runoff",&
             baseflow_field, rc=status)
        call LIS_verify(status, "ESMF_StateGet failed for Subsurface Runoff")
        
        call ESMF_FieldGet(sf_runoff_field,localDE=0,&
             farrayPtr=surface_runoff_t,&
             rc=status)
        call LIS_verify(status, "ESMF_FieldGet failed for Surface Runoff")
        
        call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow_t,&
             rc=status)
        call LIS_verify(status, "ESMF_FieldGet failed for Subsurface Runoff")
        
        allocate(surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n)))
        allocate(baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n)))
        
        surface_runoff = 0.0
        baseflow = 0.0 

        call LIS_tile2grid(n,surface_runoff,surface_runoff_t,1)
        call LIS_tile2grid(n,baseflow, baseflow_t,1)

        do r=1, LIS_rc%gnr(n)
           do c=1,LIS_rc%gnc(n) 
              surface_runoff(c,r) = (surface_runoff(c,r))*3600000 !mm/hr
              baseflow(c,r) = (baseflow(c,r))*3600000 
           enddo
        enddo
        !     stop
        call rout_it(LIS_rc%gnc(n),LIS_rc%gnr(n),&
             NLDAS_routing_struc(n)%luh,&
             NLDAS_routing_struc(n)%ltr,&
             surface_runoff,baseflow, & 
             NLDAS_routing_struc(n)%streamflow,&
             NLDAS_routing_struc(n)%runoff_intern,&
             NLDAS_routing_struc(n)%runoff_trans,&
             NLDAS_routing_struc(n)%uh_intern,&
             NLDAS_routing_struc(n)%uh_trans,&
             NLDAS_routing_struc(n)%order,&
             NLDAS_routing_struc(n)%order_n,&
             NLDAS_routing_struc(n)%area)

! the following is done to distribute the global arrays to the individual 
! processors, so that they can use the generic LIS interfaces for writing
! output -- which works based on the assumption of multiprocessors. 
! When the NLDAS routine itself is parallelized, this can (should) be taken
! out. 
        deallocate(surface_runoff)
        deallocate(baseflow)
     endif

     call LIS_grid2patch(n,LIS_rc%lsm_index, NLDAS_routing_struc(n)%streamflow,&
             streamflow_lvec)
        
     do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
        call LIS_diagnoseRoutingOutputVar(n, t,LIS_MOC_STREAMFLOW,&
             value=streamflow_lvec(t),vlevel=1,unit="m3 s-1",&
             direction="-")
     enddo
     
     deallocate(streamflow_lvec)
  endif

end subroutine NLDAS_routing_run

SUBROUTINE ROUT_IT(NX,NY,LUH,LTR,SURFACE_RUNOFF,BASEFLOW,& 
     STREAMFLOW,RUNOFF_INTERN,RUNOFF_TRANS,UH_INTERN,UH_TRANS,&
     ORDER,ORDER_N,AREA)

  IMPLICIT NONE

!     NX    -- grid points in west-east direction
!     NY    -- grid points in south-north direction
!     DT    -- time step in seconds 
!     LUH   -- length of the internal unit-hydrograph in DT
!     LTR   -- length of transport unit-hydrograph in DT
! surface runoff and baseflow are in units of mm/h
! mm/h*km2 = 1000*1000/1000*3600 = 3.6
! streamflow in units of m3 s-1

  INTEGER NX
  INTEGER NY
  INTEGER LUH
  INTEGER LTR
  INTEGER NOB
  INTEGER ORDER_N
  
  INTEGER ORDER(4,NX*NY)
  
  REAL SURFACE_RUNOFF(NX,NY)
  REAL BASEFLOW(NX,NY)
  REAL STREAMFLOW(NX,NY)
  REAL UH_INTERN(LUH,NX,NY)
  REAL UH_TRANS(LTR,NX,NY)
  REAL RUNOFF_INTERN(LUH,NX,NY)
  REAL RUNOFF_TRANS(LTR,NX,NY)
  REAL RUNOFF_IN(NX,NY)
  REAL AREA(NX,NY)
  REAL DT
  
  INTEGER N, I, J, IX, IY, IXX, IYY
  
  DO J = 1, NY
     DO I = 1, NX
        if(surface_runoff(i,j).ne.-9999.0) then 
           BASEFLOW(I,J) = BASEFLOW(I,J) * AREA(I,J)/3.6
           SURFACE_RUNOFF(I,J) = SURFACE_RUNOFF(I,J) * AREA(I,J)/3.6
        endif
     END DO
  END DO

!      WRITE(*,*) 'ORDER_N = ', ORDER_N
  DO N = 1,ORDER_N       
     IX  = ORDER(1,N)
     IY  = ORDER(2,N)
     RUNOFF_IN(IX,IY) = 0.0
  END DO
  STREAMFLOW = -9999.0

  DO N = 1,ORDER_N         
     IX  = ORDER(1,N)
     IY  = ORDER(2,N)
     IXX = ORDER(3,N)
     IYY = ORDER(4,N)
     
     IF (BASEFLOW(IX,IY) .LT. 0.0) THEN
        !            WRITE(*,*) IX, IY,BASEFLOW(IX,IY) 
        BASEFLOW(IX,IY) = 0.0
     END IF
     
     IF (SURFACE_RUNOFF(IX,IY) .LT. 0.0) THEN
!        WRITE(*,*) 'srunoff ',IX, IY, SURFACE_RUNOFF(IX,IY)
        SURFACE_RUNOFF(IX,IY) = 0.0
     END IF

     DO I = 1,LUH
        RUNOFF_INTERN(I,IX,IY) = RUNOFF_INTERN(I,IX,IY) & 
             + (SURFACE_RUNOFF(IX,IY) + BASEFLOW(IX,IY)) & 
             * UH_INTERN(I,IX,IY)
     END DO
     RUNOFF_IN(IXX,IYY) = RUNOFF_IN(IXX,IYY) +  & 
          RUNOFF_INTERN(1,IX,IY)

     DO I = 1,LTR
        RUNOFF_TRANS(I,IX,IY) = RUNOFF_TRANS(I,IX,IY) +  & 
             UH_TRANS(I,IX,IY) * & 
             RUNOFF_IN(IX,IY)
     END DO
     RUNOFF_IN(IXX,IYY) = RUNOFF_IN(IXX,IYY) +  & 
          RUNOFF_TRANS(1,IX,IY)
  END DO
  
  DO N = 1,ORDER_N      
     IX  = ORDER(1,N)
     IY  = ORDER(2,N)
     STREAMFLOW(IX,IY) = RUNOFF_TRANS(1,IX,IY) + RUNOFF_INTERN(1,IX,IY)
     
     DO I = 2,LUH
        RUNOFF_INTERN(I-1,IX,IY) = RUNOFF_INTERN(I,IX,IY)
     END DO
     RUNOFF_INTERN(LUH,IX,IY) = 0.0

     DO I = 2,LTR
        RUNOFF_TRANS(I-1,IX,IY) = RUNOFF_TRANS(I,IX,IY)
     END DO
     RUNOFF_TRANS(LTR,IX,IY) = 0.0
         
  END DO

END SUBROUTINE ROUT_IT

