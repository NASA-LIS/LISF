!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
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
! 
! !USES: 
subroutine RAPID_routing_run(n)

  use ESMF
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

  character*200         :: qout_filename

  ! for mpi
  real,   allocatable   :: runoff1_t(:)
  real,   allocatable   :: runoff2_t(:)
  real,   allocatable   :: gvar1(:)
  real,   allocatable   :: gvar2(:)
  real,   allocatable   :: gtmp1(:)
  real,   allocatable   :: gtmp2(:)
  integer               :: count1 ,c,r,ntiles,t,gid,stid,tid,l

  !temp
  character*200         :: lsm_filename
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

        ! temporary memory for mpi_gather
        allocate(gtmp1(LIS_rc%glbntiles(n)))
        allocate(gtmp2(LIS_rc%glbntiles(n)))
        
        ! temporary memory for reorganizing variable structure
        allocate(gvar1(LIS_rc%glbntiles(n)))
        allocate(gvar2(LIS_rc%glbntiles(n)))

#if (defined SPMD)
        call MPI_GATHERV(runoff1_t,LIS_deltas(n,LIS_localPet),MPI_REAL,&       ! surface runoff  
                         gtmp1,LIS_deltas(n,:),LIS_offsets(n,:),MPI_REAL,0,LIS_mpi_comm,status)
        call MPI_GATHERV(runoff2_t,LIS_deltas(n,LIS_localPet),MPI_REAL,&       ! baseflow
                         gtmp2,LIS_deltas(n,:),LIS_offsets(n,:),MPI_REAL,0,LIS_mpi_comm,status)
#endif

        ! reorganizing variable structure (only for root)
        if(LIS_masterproc) then
           count1=1
           do l=1,LIS_npes
              do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                 do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                    gid=c+(r-1)*LIS_rc%gnc(n)
                    ntiles=LIS_domain(n)%ntiles_pergrid(gid)
                    stid=LIS_domain(n)%str_tind(gid)
            
                    do t=1,ntiles
                       tid=stid + t-1
                       gvar1(tid)=gtmp1(count1)  ! surface runoff
                       gvar2(tid)=gtmp2(count1)  ! baseflow
                       count1=count1 + 1         
                    enddo
                 enddo
              enddo          
           enddo
        endif
       
       ! The message is sent from the root process to all processes in the group
#if (defined SPMD)
       call MPI_BCAST(gvar1,LIS_rc%glbntiles(n),MPI_REAL,0, &
            LIS_mpi_comm, status)
       call MPI_BCAST(gvar2,LIS_rc%glbntiles(n),MPI_REAL,0, &
            LIS_mpi_comm, status)
#endif

        surface_runoff=reshape(gvar1,(/LIS_rc%gnc(n),LIS_rc%gnr(n)/))
        baseflow=reshape(gvar2,(/LIS_rc%gnc(n),LIS_rc%gnr(n)/))

        deallocate(runoff1_t)
        deallocate(runoff2_t)
        deallocate(gtmp1)
        deallocate(gtmp2)
        deallocate(gvar1)
        deallocate(gvar2)
     endif

     !TODO
     ! output file name
     call LIS_create_output_directory('ROUTING')
     call LIS_create_output_filename(n,qout_filename,model_name='ROUTING', &
             writeint=RAPID_routing_struc(n)%outInterval)

     ! output file name (LSM)
     call LIS_create_output_filename(n,lsm_filename,model_name='SURFACEMODEL', &
             writeint=RAPID_routing_struc(n)%outInterval)

     ! run RAPID
     call RAPID_model_main (RAPID_routing_struc(n)%bQinit,RAPID_routing_struc(n)%bQfinal,RAPID_routing_struc(n)%bV,               &
                            RAPID_routing_struc(n)%bhum,RAPID_routing_struc(n)%bfor,RAPID_routing_struc(n)%bdam,                  &
                            RAPID_routing_struc(n)%binfluence,RAPID_routing_struc(n)%buq,                                         &
                            RAPID_routing_struc(n)%run_opt,RAPID_routing_struc(n)%routing_opt,RAPID_routing_struc(n)%phi_opt,     &
                            RAPID_routing_struc(n)%connectfile,RAPID_routing_struc(n)%max_reach,RAPID_routing_struc(n)%n_riv_tot, &
                            RAPID_routing_struc(n)%weightfile,RAPID_routing_struc(n)%n_wei_table,                                 &
                            RAPID_routing_struc(n)%basinIDfile,RAPID_routing_struc(n)%n_riv_bas,                                  &
                            RAPID_routing_struc(n)%kfile,RAPID_routing_struc(n)%xfile,                                            &
                            RAPID_routing_struc(n)%nmlfile,qout_filename,                                                         &
                            LIS_rc%gnc(n),LIS_rc%gnr(n),surface_runoff,baseflow,RAPID_routing_struc(n)%initCheck,                 &
                            RAPID_routing_struc(n)%rst_Qout,RAPID_routing_struc(n)%startMode,                                     &
                            RAPID_routing_struc(n)%dt,RAPID_routing_struc(n)%routingInterval,lsm_filename)

     deallocate(surface_runoff)
     deallocate(baseflow)
     
  endif
end subroutine RAPID_routing_run
