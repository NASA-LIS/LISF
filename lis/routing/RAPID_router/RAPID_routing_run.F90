!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! 
! !DESCRIPTION: 
! 
!  Reference: 
!
! !REVISION HISTORY: 
! 18 Mar 2021: Yeosang Yoon: Initial implementation in LIS 
! 25 Oct 2022: Yeosang Yoon: Support to run with LSM ensemble mean runoff variables
! 27 Apr 2023: Eric Kemp: Updated length of output file.
! !USES: 
subroutine RAPID_routing_run(n)

  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_routingMod
  use LIS_logMod
  use LIS_historyMod
  use LIS_histDataMod
  use LIS_constantsMod
  use LIS_fileIOMod
  use RAPID_routingMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_mpiMod
  
  implicit none
  
  integer, intent(in)   :: n  
  type(ESMF_Field)      :: sf_runoff_field
  type(ESMF_Field)      :: baseflow_field
  real,   allocatable   :: surface_runoff(:,:)
  real,   pointer       :: surface_runoff_t(:)
  real,   pointer       :: baseflow_t(:)
  real,   allocatable   :: baseflow(:,:)
  integer               :: status
  logical               :: alarmCheck

  !character*200         :: qout_filename
  character(len=LIS_CONST_PATH_LEN) :: qout_filename ! EMK
  ! for mpi
  real,   allocatable   :: runoff1_t(:)
  real,   allocatable   :: runoff2_t(:)
  real,   allocatable   :: gvar1(:,:)
  real,   allocatable   :: gvar2(:,:)
  real,   allocatable   :: gtmp1(:)
  real,   allocatable   :: gtmp2(:)
  integer               :: count1 ,c,r,ntiles,t,gid,stid,tid,l,i,m

  real,   allocatable   :: meanv1(:) 
  real,   allocatable   :: meanv2(:)
! _______________________________________________

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "RAPID router model alarm")
  if(alarmCheck) then 
     call ESMF_StateGet(LIS_runoff_state(n),"Surface Runoff",sf_runoff_field,&
          rc=status)
     call LIS_verify(status, "ESMF_StateGet failed for Surface Runoff")
     
     call ESMF_StateGet(LIS_runoff_state(n),"Subsurface Runoff",&
          baseflow_field, rc=status)
     call LIS_verify(status, "ESMF_StateGet failed for Subsurface Runoff")
              
     call ESMF_FieldGet(sf_runoff_field,localDE=0,farrayPtr=surface_runoff_t,&
          rc=status)
     call LIS_verify(status, "ESMF_FieldGet failed for Surface Runoff")
              
     call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow_t,&
          rc=status)
     call LIS_verify(status, "ESMF_FieldGet failed for Subsurface Runoff")
                         
     allocate(surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n)))
     allocate(baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n)))

     if(LIS_npes==1) then    !single core
        call LIS_tile2grid(n,surface_runoff,surface_runoff_t,1)
        call LIS_tile2grid(n,baseflow,baseflow_t,1)   
     else                    !mpi
        allocate(runoff1_t(LIS_rc%ntiles(n)))
        allocate(runoff2_t(LIS_rc%ntiles(n)))
       
        runoff1_t=surface_runoff_t
        runoff2_t=baseflow_t

        allocate(gtmp1(LIS_rc%glbngrid(n)))            ! temporary memory for mpi_gather
        allocate(gtmp2(LIS_rc%glbngrid(n)))
        allocate(gvar1(LIS_rc%gnc(n),LIS_rc%gnr(n)))   ! temporary memory for reorganizing variable structure
        allocate(gvar2(LIS_rc%gnc(n),LIS_rc%gnr(n)))
        gtmp1 = 0.0
        gtmp2 = 0.0
        gvar1 = 0.0
       
        ! only support LSM ensmble mean runoff variables (2022/10/25)
        if(RAPID_routing_struc(n)%useens==1) then     
           allocate(meanv1(LIS_rc%ngrid(n)))
           allocate(meanv2(LIS_rc%ngrid(n)))
           meanv1 = 0
           meanv2 = 0
           do i=1,LIS_rc%ntiles(n), LIS_rc%nensem(n)
              c=LIS_domain(n)%tile(i)%index
              do m=1, LIS_rc%nensem(n)
                 t = i+m-1
                 if (runoff1_t(t) == -9999.0 ) then ! surface runoff
                    meanv1(c) = -9999.0
                 else  ! make ensemble mean
                    meanv1(c) = meanv1(c) + runoff1_t(t)*LIS_domain(n)%tile(t)%fgrd*&
                       LIS_domain(n)%tile(t)%pens
                 endif

                 if (runoff2_t(t) == -9999.0 ) then ! baseflow
                    meanv2(c) = -9999.0
                 else
                    meanv2(c) = meanv2(c) + runoff2_t(t)*LIS_domain(n)%tile(t)%fgrd*&
                       LIS_domain(n)%tile(t)%pens
                 endif
              enddo
           enddo

#if (defined SPMD)
           call MPI_GATHERV(meanv1,LIS_deltas(n,LIS_localPet),MPI_REAL,&   ! surface runoff  
                            gtmp1,LIS_deltas(n,:),LIS_offsets(n,:),MPI_REAL,0,LIS_mpi_comm,status)
           call MPI_GATHERV(meanv2,LIS_deltas(n,LIS_localPet),MPI_REAL,&   ! baseflow
                            gtmp2,LIS_deltas(n,:),LIS_offsets(n,:),MPI_REAL,0,LIS_mpi_comm,status)
#endif
           deallocate(meanv1)
           deallocate(meanv2)
        else if(RAPID_routing_struc(n)%useens==0) then ! OL loop
#if (defined SPMD)
           call MPI_GATHERV(runoff1_t,LIS_deltas(n,LIS_localPet),MPI_REAL,&   ! surface runoff
                            gtmp1,LIS_deltas(n,:),LIS_offsets(n,:),MPI_REAL,0,LIS_mpi_comm,status)
           call MPI_GATHERV(runoff2_t,LIS_deltas(n,LIS_localPet),MPI_REAL,&   ! baseflow
                            gtmp2,LIS_deltas(n,:),LIS_offsets(n,:),MPI_REAL,0,LIS_mpi_comm,status)
#endif
        endif

        ! reorganizing variable structure (only for root)
        if(LIS_masterproc) then
           gvar1 = LIS_rc%udef
           gvar2 = LIS_rc%udef
           count1=1
           do l=1,LIS_npes
              do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                 do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                    gid=c+(r-1)*LIS_rc%gnc(n)
                    ntiles=LIS_domain(n)%ntiles_pergrid(gid)
                    stid=LIS_domain(n)%str_tind(gid)

                    if(ntiles.ne.0) then
                       if(r.ge.LIS_nss_ind(n,l).and.&
                            r.le.LIS_nse_ind(n,l).and.&
                            c.ge.LIS_ews_ind(n,l).and.&
                            c.le.LIS_ewe_ind(n,l)) then !points not in halo
                          gvar1(c,r) = gtmp1(count1) ! surface runoff
                          gvar2(c,r) = gtmp2(count1) ! baseflow
                       endif
                       count1 = count1 + 1
                    endif
                 enddo
              enddo          
           enddo
        endif
       
       ! The message is sent from the root process to all processes in the group
#if (defined SPMD)
       call MPI_BCAST(gvar1,LIS_rc%glbngrid(n),MPI_REAL,0, &
            LIS_mpi_comm, status)
       call MPI_BCAST(gvar2,LIS_rc%glbngrid(n),MPI_REAL,0, &
            LIS_mpi_comm, status)
#endif
        
        surface_runoff = gvar1
        baseflow = gvar2
        
        deallocate(runoff1_t)
        deallocate(runoff2_t)
        deallocate(gtmp1)
        deallocate(gtmp2)
        deallocate(gvar1)
        deallocate(gvar2)
     endif

     ! output file name
     call LIS_create_output_directory('ROUTING')
     call LIS_create_output_filename(n,qout_filename,model_name='ROUTING', &
             writeint=RAPID_routing_struc(n)%outInterval)

     ! run RAPID
#ifdef PETSc
     call RAPID_model_main (n,RAPID_routing_struc(n)%bQinit,RAPID_routing_struc(n)%bQfinal,RAPID_routing_struc(n)%bV,             &
                            RAPID_routing_struc(n)%bhum,RAPID_routing_struc(n)%bfor,RAPID_routing_struc(n)%bdam,                  &
                            RAPID_routing_struc(n)%binfluence,RAPID_routing_struc(n)%buq,                                         &
                            RAPID_routing_struc(n)%run_opt,RAPID_routing_struc(n)%routing_opt,RAPID_routing_struc(n)%phi_opt,     &
                            RAPID_routing_struc(n)%connectfile,RAPID_routing_struc(n)%max_reach,RAPID_routing_struc(n)%n_riv_tot, &
                            RAPID_routing_struc(n)%weightfile,RAPID_routing_struc(n)%n_wei_table,                                 &
                            RAPID_routing_struc(n)%basinIDfile,RAPID_routing_struc(n)%n_riv_bas,                                  &
                            RAPID_routing_struc(n)%kfile,RAPID_routing_struc(n)%xfile,                                            &
                            RAPID_routing_struc(n)%nmlfile,qout_filename,                                                         &
                            LIS_rc%gnc(n),LIS_rc%gnr(n),surface_runoff,baseflow,RAPID_routing_struc(n)%initCheck,                 &
                            RAPID_routing_struc(n)%dt,RAPID_routing_struc(n)%routingInterval)
#endif
     deallocate(surface_runoff)
     deallocate(baseflow)
     
  endif
end subroutine RAPID_routing_run
