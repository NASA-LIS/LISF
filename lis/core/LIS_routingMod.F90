!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_routingMod
!BOP
! 
! !MODULE: LIS_routingMod
! 
! !DESCRIPTION: 
!  The code in this file provides the top level calls to manage the 
!  operation of different runoff-routing algorithms and models. 
! 
! !REVISION HISTORY: 
!  6 May 2011: Sujay Kumar, Initial implementation
!  1 Oct 2022: Yeosang Yoon; excluded RAPID from sws and DA modules
! 
! !USES: 
  use LIS_coreMod
  use LIS_PRIV_tileMod
  use LIS_logMod
  use LIS_histDataMod
  use LIS_DAobs_pluginMod
  use LIS_DAobservationsMod
  use LIS_mpiMod
  use map_utils
  use ESMF

  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_routing_init         ! initialize the routing model
  public :: LIS_routing_readrestart  ! read the routing model restart file
  public :: LIS_routing_run           ! execute the routing model 
  public :: LIS_routing_writeoutput   ! write the output file
  public :: LIS_routing_writerestart   ! write the restart file

  public :: LIS_routing_perturb_states
  public :: LIS_routing_DAGetObsPred
  public :: LIS_routing_DAGetStateVar
  public :: LIS_routing_DASetStateVar
  public :: LIS_routing_DAScaleStateVar
  public :: LIS_routing_DADescaleStateVar
  public :: LIS_routing_DAUpdateState
  public :: LIS_routing_DAQCState
  public :: LIS_routing_DAgetStateSpaceSize
  public :: LIS_routing_DAextractStateVector
  public :: LIS_routing_DAsetFreshIncrementsStatus
  public :: LIS_routing_DAgetFreshIncrementsStatus
  public :: LIS_routing_DAsetAnlysisUpdates
  public :: LIS_routing_DAmapTileSpaceToObsSpace
  public :: LIS_routing_DAgetStateVarNames
  public :: LIS_routing_getlatlons
  public :: LIS_routing_DAobsTransform
  public :: LIS_routing_DAmapObsToRouting
  public :: LIS_routing_DAqcObsState
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LIS_runoff_state
  public :: LIS_routing_state         !state vector in DA
  public :: LIS_routing_pert_state
  public :: LIS_routing_incr_state
!EOP
  
  type(ESMF_State), allocatable :: LIS_runoff_state(:)
  type(ESMF_State), allocatable :: LIS_routing_state(:,:)
  type(ESMF_State), allocatable :: LIS_routing_pert_state(:,:)
  type(ESMF_State), allocatable :: LIS_routing_incr_state(:,:)


  contains
!BOP
! !ROUTINE: LIS_routing_init
! \label{LIS_routing_init}
! 
! !INTERFACE: 
  subroutine LIS_routing_init
! !USES: 
! 
! !DESCRIPTION: 
!EOP
    integer           :: rc
    integer           :: i
    integer           :: n
    integer           :: status
    character*1       :: nestid(2)
    character*100     :: temp

    integer              :: c,r,j,k,kk,m,l
    integer              :: cg,rg
    integer              :: ftn
    integer              :: ierr
    integer              :: count
    integer, allocatable :: ntiles_pergrid(:)
    integer, allocatable :: gtmp(:,:)
    integer, allocatable :: gtmp1(:)

    call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%routingmodel, &
         label="Routing model:",default="none", rc=rc)

    if(LIS_rc%routingmodel.eq."RAPID router") then

       allocate(LIS_runoff_state(LIS_rc%nnest))

       do n=1, LIS_rc%nnest

          write(unit=temp,fmt='(i2.2)') n
          read(unit=temp,fmt='(2a1)') nestid

          LIS_runoff_state(n) = ESMF_StateCreate(name="LIS Runoff State"//&
               nestid(1)//nestid(2), rc=status)
          call LIS_verify(status, "ESMF_StateCreate failed in LIS_routing_init")

       enddo

       if(LIS_rc%lsm.eq."none") then
          call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%runoffdatasource, &
               label="External runoff data source:", rc=rc)
          call LIS_verify(rc,"External runoff data source: not defined")
       endif
      
       call routinginit(trim(LIS_rc%routingmodel)//char(0))
    else if((LIS_rc%routingmodel.eq."HYMAP router") .or. &
            (LIS_rc%routingmodel.eq."HYMAP2 router") .or. &
            (LIS_rc%routingmodel.eq."NLDAS router")) then

       allocate(LIS_routing(LIS_rc%nnest))

       do n=1,LIS_rc%nnest
          allocate(LIS_routing(n)%dommask(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(LIS_routing(n)%nextx(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(LIS_routing(n)%gindex(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(LIS_routing(n)%ntiles_pergrid(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
          allocate(ntiles_pergrid(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       enddo

       allocate(LIS_rc%nroutinggrid(LIS_rc%nnest))
       allocate(LIS_rc%glbnroutinggrid(LIS_rc%nnest))
       allocate(LIS_rc%glbnroutinggrid_red(LIS_rc%nnest))
       
       allocate(LIS_routing_gdeltas(LIS_rc%nnest,0:LIS_npes-1))
       allocate(LIS_routing_goffsets(LIS_rc%nnest,0:LIS_npes-1))

       allocate(LIS_rc%Routing_DAinst_valid(LIS_rc%ndas))
       
       do i=1,LIS_rc%ndas
          call LIS_isDAinstanceValid(LIS_rc%daset(i),&
               "Routing",LIS_rc%Routing_DAinst_valid(i))
       enddo

       allocate(LIS_runoff_state(LIS_rc%nnest))

       do n=1, LIS_rc%nnest

          write(unit=temp,fmt='(i2.2)') n
          read(unit=temp,fmt='(2a1)') nestid

          LIS_runoff_state(n) = ESMF_StateCreate(name="LIS Runoff State"//&
               nestid(1)//nestid(2), rc=status)
          call LIS_verify(status, "ESMF_StateCreate failed in LIS_routing_init")
               
       enddo

       if(LIS_rc%lsm.eq."none") then 
          call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%runoffdatasource, &
               label="External runoff data source:", rc=rc)
          call LIS_verify(rc,"External runoff data source: not defined")
       endif

       call routinginit(trim(LIS_rc%routingmodel)//char(0))

       do n=1,LIS_rc%nnest
          allocate(LIS_routing(n)%grid(LIS_rc%nroutinggrid(n)))

          allocate(LIS_routing(n)%tile(LIS_rc%nroutinggrid(n)*&
               LIS_rc%nensem(n)))

          call LIS_routingHistDataInit(n,LIS_rc%nroutinggrid(n)*&
               LIS_rc%nensem(n))
          
          k = 1
          kk = 1
          LIS_routing(n)%gindex = -1
          ntiles_pergrid = 0 

          do c=1,LIS_rc%lnc(n)
             do r=1,LIS_rc%lnr(n)
                cg = c+LIS_ews_halo_ind(n,LIS_localPet+1)-1
                rg = r+LIS_nss_halo_ind(n,LIS_localPet+1)-1
                if(LIS_routing(n)%dommask(c,r).gt.0.and.&
                     LIS_routing(n)%nextx(cg,rg).ne.LIS_rc%udef) then      

                   LIS_routing(n)%grid(kk)%lat = &
                        LIS_domain(n)%lat(c+(r-1)*LIS_rc%lnc(n))
                   LIS_routing(n)%grid(kk)%lon = &
                        LIS_domain(n)%lon(c+(r-1)*LIS_rc%lnc(n))

                   LIS_routing(n)%grid(kk)%col = c
                   LIS_routing(n)%grid(kk)%row = r

                   do m=1,LIS_rc%nensem(n)
                      LIS_routing(n)%tile(k)%ensem = m
                      LIS_routing(n)%tile(k)%row = r
                      LIS_routing(n)%tile(k)%col = c
                      LIS_routing(n)%tile(k)%index = kk
                      ntiles_pergrid(r+(c-1)*LIS_rc%lnr(n)) = &
                           ntiles_pergrid(r+(c-1)*LIS_rc%lnr(n)) + 1
                      k = k + 1
                   enddo
                   LIS_routing(n)%gindex(c,r) = kk

                   kk = kk + 1
                endif
             enddo
          enddo
          
          if(LIS_masterproc) then 
             allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
             allocate(gtmp1(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
          else
             allocate(gtmp1(1))
          endif
#if (defined SPMD)
          call MPI_GATHERV(ntiles_pergrid,LIS_deltas(n,LIS_localPet),&
               MPI_INTEGER,gtmp1,&
               LIS_deltas(n,:),LIS_offsets(n,:),&
               MPI_INTEGER,0,LIS_mpi_comm,ierr)
#else
          gtmp1 = ntiles_pergrid
#endif
          if(LIS_masterproc) then 
             count=1
             do l=1,LIS_npes
                do c=LIS_ews_ind(n,l), LIS_ewe_ind(n,l)
                   do r=LIS_nss_ind(n,l), LIS_nse_ind(n,l)
                      gtmp(c,r) = gtmp1(count)
                      count = count+1
                   enddo
                enddo
             enddo

             count=1
             do c=1,LIS_rc%gnc(n)
                do r=1,LIS_rc%gnr(n)
                   LIS_routing(n)%ntiles_pergrid(count) = gtmp(c,r)
                   count = count+1
                enddo
             enddo

             deallocate(gtmp)
             deallocate(gtmp1)
          else
             deallocate(gtmp1)
          endif

#if (defined SPMD)
          call MPI_Bcast(LIS_routing(n)%ntiles_pergrid,&
               LIS_rc%gnc(n)*LIS_rc%gnr(n),&
               MPI_INTEGER,0,LIS_mpi_comm,ierr)
#endif   
          do k=1,LIS_rc%nroutinggrid(n)*LIS_rc%nensem(n)
             LIS_routing(n)%tile(k)%pens = 1.0/float(LIS_rc%nensem(n))
             LIS_routing(n)%tile(k)%fgrd = 1.0  !no subgrid routing tiling for now
          enddo

       end do
    else if (LIS_rc%routingmodel.eq."none") then

       allocate(LIS_rc%Routing_DAinst_valid(LIS_rc%ndas))

       do i=1,LIS_rc%ndas
          call LIS_isDAinstanceValid(LIS_rc%daset(i),&
               "Routing",LIS_rc%Routing_DAinst_valid(i))
       enddo

    end if

  end subroutine LIS_routing_init


!BOP
! !ROUTINE: LIS_routing_readrestart
! \label{LIS_routing_readrestart}
! 
! !INTERFACE: 
  subroutine LIS_routing_readrestart
! !USES: 

! 
! !DESCRIPTION: 
!EOP
    if(LIS_rc%routingmodel.ne."none") then 
       call routingreadrestart(trim(LIS_rc%routingmodel)//char(0))
    endif
  
  end subroutine LIS_routing_readrestart

!BOP
! !ROUTINE: LIS_routing_run
! \label{LIS_routing_run}
! 
! !INTERFACE: 
  subroutine LIS_routing_run(n)
! !USES: 

! 
! !DESCRIPTION: 
!EOP
    integer, intent(in) :: n

    if(LIS_rc%routingmodel.ne."none") then
       if(LIS_rc%lsm.ne."none") then 
          call lsmroutinggetrunoff(trim(LIS_rc%lsm)//"+"//&
               trim(LIS_rc%routingmodel)//char(0),n)
       endif
       if(LIS_rc%glaciermodel.ne."none") then
          call glacierroutinggetrunoff(trim(LIS_rc%glaciermodel)//"+"//&
               trim(LIS_rc%routingmodel)//char(0),n)
       endif

       call routingrun(trim(LIS_rc%routingmodel)//char(0),n)
       if(LIS_rc%lsm.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then 
          call lsmroutinggetsws(trim(LIS_rc%lsm)//"+"//&
               trim(LIS_rc%routingmodel)//char(0),n)
       endif
    endif

  end subroutine LIS_routing_run

!BOP
! !ROUTINE: LIS_routing_writeoutput
! \label{LIS_routing_writeoutput}
! 
! !INTERFACE: 
  subroutine LIS_routing_writeoutput(n)
! !USES: 

! 
! !DESCRIPTION: 
!EOP

    integer, intent(in) :: n
  
    if(LIS_rc%routingmodel.ne."none") then 
       call routingoutput(trim(LIS_rc%routingmodel)//char(0),n)
    endif

  end subroutine LIS_routing_writeoutput


!BOP
! !ROUTINE: LIS_routing_writerestart
! \label{LIS_routing_writerestart}
! 
! !INTERFACE: 
  subroutine LIS_routing_writerestart(n)
! !USES: 

! 
! !DESCRIPTION: 
!EOP

    integer, intent(in) :: n
  
    if(LIS_rc%wopt_rst.ne.0) then 
       if(LIS_rc%routingmodel.ne."none") then 
          call routingwriterestart(trim(LIS_rc%routingmodel)//char(0),n)
       endif
    endif

  end subroutine LIS_routing_writerestart

!BOP
! 
! !ROUTINE: LIS_routing_perturb_states
! \label{LIS_routing_perturb_states}
! 
! !INTERFACE: 
  subroutine LIS_routing_perturb_states(n)
! !USES: 

! !ARGUMENTS: 
    integer,    intent(IN)  :: n 
!
! !DESCRIPTION:
! This interface provides the entry point to the Routing routines that 
! computes the perturbations on the Routing prognostic state variables
! 
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[perturbmethod] (\ref{perturbmethod}) \newline
!    invokes the abstract method to invoke the perturbation
!    algorithm to perturb Routing prognostic variables
!  \item[routingdagetstatevar] (\ref{routingdagetstatevar}) \newline
!    obtains the list of prognostic variables
!  \item[applyRoutingPert] (\ref{applyRoutingPert}) \newline
!    applies the specified perturbations to the Routing state
!  \item[routingdaqcstate] (\ref{routingdaqcstate}) \newline
!   performs the QC of the perturbed Routing state
!  \item[routingdasetstatevar] (\ref{routingdasetstatevar}) \newline
!   assigns the prognostic variables back to the model states
! \end{description}
!EOP

    real                    :: curr_time
    integer                 :: k 


    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then 
       do k=1, LIS_rc%nperts
          if(LIS_rc%Routing_DAinst_valid(k)) then
             if(LIS_rc%perturb_state(k).ne."none") then 
                curr_time = float(LIS_rc%hr)*3600+&
                     60*float(LIS_rc%mn)+float(LIS_rc%ss)
                if(mod(curr_time,real(LIS_rc%pertstateInterval(k))).eq.0) then
!------------------------------------------------------------------------
!   Returns the perturbed state based on the chosen algorithm
!------------------------------------------------------------------------
                   call perturbmethod(trim(LIS_rc%perturb_state(k))//char(0),&
                        4, n,k,&
                        LIS_Routing_State(n,k),LIS_Routing_Pert_State(n,k))
!------------------------------------------------------------------------
!   Propagate step or applying the perturbations to prognostic states
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!   apply the routing perturbations to the the Routing state
!------------------------------------------------------------------------
                   call routingdagetstatevar(trim(LIS_rc%routingModel)//"+"//&
                        trim(LIS_rc%daset(k))//char(0), n, LIS_Routing_State(n,k))
                   call applyRoutingPert(n, k, &
                        LIS_Routing_State(n,k),LIS_Routing_Pert_State(n,k))
!------------------------------------------------------------------------
!   Diagnose the perturbed state (updates the model prognostic states)
!------------------------------------------------------------------------
                   call routingdaqcstate(trim(LIS_rc%routingModel)//"+"//&
                        trim(LIS_rc%daset(k))//char(0), n, LIS_Routing_State(n,k))
                   call routingdasetstatevar(trim(LIS_rc%routingModel)//"+"//&
                        trim(LIS_rc%daset(k))//char(0), n, LIS_Routing_State(n,k)) 
                endif
             endif
          endif
       enddo
    endif
  end subroutine LIS_routing_perturb_states

!BOP
! 
! !ROUTINE: applyRoutingPert
! \label{applyRoutingPert}
! 
! !INTERFACE: 
  subroutine applyRoutingPert(n, k, LIS_Routing_State, LIS_Routing_Pert_State)
! !USES: 

! !ARGUMENTS:     
    integer, intent(IN) :: n 
    integer, intent(IN) :: k
    type(ESMF_State)    :: LIS_Routing_State
    type(ESMF_State)    :: LIS_Routing_Pert_State
!
! !DESCRIPTION:
! 
!  This routine applies the specified perturbations to the 
!  Routing prognostic state variables. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]        index of the nest
!   \item[Routing\_State]   ESMF State with prognostic variables
!   \item[Routing\_Pert\_State]    ESMF State with prognostic variable 
!                       perturbations
!  \end{description}
!EOP

    integer                   :: i,j,t,f,m,t_unpert
    integer                   :: stateSize
    integer                   :: lsm_state_count
    integer                   :: status
    character*100,    allocatable :: lsm_state_objs(:)
    type(ESMF_Field), allocatable :: lsm_field(:)
    real, pointer             :: lsm_temp(:)
    real                      :: delta(LIS_rc%nstvars(k))
    integer, allocatable          :: typ(:)
    real, allocatable         :: stvar(:,:)
    real, allocatable         :: stvar_up(:,:)
    real, allocatable         :: pert(:,:)
    integer                   :: nxx

    call LIS_routing_DAgetStateSpaceSize(n,k,stateSize)

    allocate(stvar(LIS_rc%nstvars(k),&
         stateSize))
    allocate(stvar_up(LIS_rc%nstvars(k),&
         stateSize))
    allocate(pert(LIS_rc%nstvars(k),&
         stateSize))

    call ESMF_StateGet(LIS_Routing_State,itemCount=lsm_state_count,rc=status)
    call LIS_verify(status, &
         "ESMF_StateGet failed in applyRoutingPert")
    
    allocate(lsm_state_objs(lsm_state_count))
    allocate(lsm_field(lsm_state_count))
    
    call ESMF_StateGet(LIS_Routing_State,itemNameList=lsm_state_objs,rc=status)
    call LIS_verify(status,&
         "ESMF_StateGet failed in applyRoutingPert")        
    
    do i=1,lsm_state_count
       call ESMF_StateGet(LIS_Routing_State,lsm_state_objs(i),lsm_field(i),&
            rc=status)
       call LIS_verify(status,&
            "ESMF_StateGet failed in applyRoutingPert")
       call ESMF_FieldGet(lsm_field(i), localDE=0,farrayPtr=lsm_temp,rc=status)
       call LIS_verify(status,&
            "ESMF_FieldGet failed in applyRoutingPert")
       stvar(i,:) = lsm_temp(:)
       stvar_up(i,:) = lsm_temp(:)
    enddo
    deallocate(lsm_field)

    allocate(lsm_field(lsm_state_count))
    allocate(typ(lsm_state_count))
    
    do i=1,lsm_state_count
       call ESMF_StateGet(LIS_Routing_Pert_State,&
            trim(lsm_state_objs(i)),&
            lsm_field(i),rc=status)
       call LIS_verify(status,&
            "ESMF_StateGet failed in applyRoutingPert")

       call ESMF_AttributeGet(lsm_field(i),"Perturbation Type",typ(i),&
            rc=status)
       call LIS_verify(status,&
            "ESMF_AttributeGet failed in applyRoutingPert")

       call ESMF_FieldGet(lsm_field(i), localDE=0,farrayPtr=lsm_temp,rc=status)
       call LIS_verify(status,&
            "ESMF_FieldGet failed in applyRoutingPert")
       pert(i,:) = lsm_temp(:)
    enddo

    if(LIS_rc%perturb_state(k).eq."uniform") then 

       nxx = LIS_rc%nensem(n)/LIS_rc%nmetforc
       do i=1,stateSize/LIS_rc%nensem(n)
          do f=1,LIS_rc%nmetforc
             t_unpert = (i-1)*LIS_rc%nensem(n) + f*nxx
             do j=1,LIS_rc%nstvars(k)
                do m=1,nxx
                   t=(i-1)*LIS_rc%nensem(n)+(f-1)*nxx + m
!                if(typ(j).eq.0) then 
                   if(LIS_rc%pert_bias_corr.eq.1) then 
                      if(m.ne.nxx) then 
                         if(j.eq.m) then 
                            stvar(j,t) = stvar(j,t) + pert(j,t)
                         endif
                      endif
                   else
                      stvar(j,t) = stvar(j,t) + pert(j,t)
                   endif
                enddo
             enddo
          enddo
       enddo
    else

! If perturbation bias correction scheme is turned on, then 
! apply the perturbations to (nensem-1) ensemble members. Keep 
! the last ensemble member for each tile unperturbed. 
! This will be used later to apply the Ryu et al. JHM (2009) 
! perturbation bias correction
!
       nxx = LIS_rc%nensem(n)/LIS_rc%nmetforc
       do i=1,stateSize/LIS_rc%nensem(n)
          do f=1,LIS_rc%nmetforc
             t_unpert = (i-1)*LIS_rc%nensem(n) + f*nxx
             do j=1,LIS_rc%nstvars(k)
                do m=1,nxx
                   t=(i-1)*LIS_rc%nensem(n)+(f-1)*nxx + m
                   if(typ(j).eq.0) then 
                      if(LIS_rc%pert_bias_corr.eq.1) then 
                         if(m.ne.nxx) then 
                            stvar(j,t) = stvar(j,t) + pert(j,t)
                         endif
                      else
                         stvar(j,t) = stvar(j,t) + pert(j,t)
                      endif
                   elseif(typ(j).eq.1) then 
                      if(LIS_rc%pert_bias_corr.eq.1) then 
                         if(m.ne.nxx) then 
                            stvar(j,t) = stvar(j,t) * pert(j,t)
                         endif
                      else
                         stvar(j,t) = stvar(j,t) * pert(j,t)
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
       
       if(LIS_rc%pert_bias_corr.eq.1) then 
          do i=1,stateSize/LIS_rc%nensem(n)
             do f=1,LIS_rc%nmetforc
                t_unpert = (i-1)*LIS_rc%nensem(n) + f*nxx
                do j=1,LIS_rc%nstvars(k)
                   delta(j) = 0.0
                   do m=1,nxx-1
                      t = (i-1)*LIS_rc%nensem(n)+(f-1)*nxx+m
                      delta(j) = delta(j) + (stvar(j,t)-stvar(j,t_unpert))
                   enddo
                enddo
                
                do j=1,LIS_rc%nstvars(k)
                   delta(j) = delta(j)/(nxx-1)
                   do m=1,nxx-1
                      t = (i-1)*LIS_rc%nensem(n)+(f-1)*nxx+m
                      stvar(j,t) = stvar(j,t) - delta(j)
                   enddo
                enddo
             end do
          enddo
       endif
    endif

    deallocate(lsm_state_objs)
    deallocate(lsm_field)
    deallocate(typ)
    
!set the perturbed variables 
   
   allocate(lsm_state_objs(lsm_state_count))
   allocate(lsm_field(lsm_state_count))

   call ESMF_StateGet(LIS_Routing_State,itemNameList=lsm_state_objs,rc=status)
   call LIS_verify(status,&
        "ESMF_StateGet failed in applyRoutingPert")       
 
   do i=1,lsm_state_count
      call ESMF_StateGet(LIS_Routing_State,lsm_state_objs(i),lsm_field(i),&
           rc=status)
      call LIS_verify(status,&
           "ESMF_StateGet failed in applyRoutingpert")
      call ESMF_FieldGet(lsm_field(i), localDE=0,farrayPtr=lsm_temp,rc=status)
      call LIS_verify(status,&
           "ESMF_FieldGet failed in applyRoutingPert")
      lsm_temp(:) = stvar(i,:) 
   enddo

   do i=1,lsm_state_count
      call ESMF_StateGet(LIS_Routing_Pert_State,&
           trim(lsm_state_objs(i)),&
           lsm_field(i),rc=status)
      call LIS_verify(status,&
           "ESMF_StateGet failed in applyRoutingPert")
      
      call ESMF_FieldGet(lsm_field(i), localDE=0,farrayPtr=lsm_temp,rc=status)
      call LIS_verify(status,&
           "ESMF_FieldGet failed in applyRoutingPert")
      lsm_temp(:) = stvar(i,:) - stvar_up(i,:)
    enddo

   deallocate(lsm_state_objs)
   deallocate(lsm_field)

   deallocate(stvar)
   deallocate(stvar_up)
   deallocate(pert)

 end subroutine applyRoutingPert

  subroutine LIS_routing_DAGetObsPred(n,k,Obs_Pred)
    
    integer                :: n
    integer                :: k
    real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))

    integer                :: m

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%Routing_DAinst_valid(k)) then
          call routingdagetobspred(trim(LIS_rc%routingmodel)//"+"//&
               trim(LIS_rc%daset(k))//char(0),n, k, Obs_pred)
       endif
    endif
  end subroutine LIS_routing_DAGetObsPred

!BOP
!
!ROUTINE: LIS_routing_DAGetStateVar
! \label{LIS_routing_DAGetStateVar}
!
! !INTERFACE:
  subroutine LIS_routing_DAGetStateVar(n,k)

! !ARGUMENTS:
    integer                :: n
    integer                :: k

!
! !DESCRIPTION:
! 
!  This interface invokes the land model to return the state vector
!  used in data assimilation. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then 
       if(LIS_rc%routing_DAinst_valid(k)) then
          call routingdagetstatevar(trim(LIS_rc%routingmodel)//"+"//&
               trim(LIS_rc%daset(k))//char(0), n, LIS_Routing_State(n,k))
       endif
    endif
  end subroutine LIS_routing_DAGetStateVar
  
!BOP
!
!ROUTINE: LIS_routing_DASetStateVar
! \label{LIS_routing_DASetStateVar}
!
! !INTERFACE:
  subroutine LIS_routing_DASetStateVar(n,k)

! !ARGUMENTS:
    integer                :: n
    integer                :: k

!
! !DESCRIPTION:
! 
!  This interface invokes the land model to set the state vector
!  used in data assimilation. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          call routingdasetstatevar(trim(LIS_rc%routingmodel)//"+"//&
               trim(LIS_rc%daset(k))//char(0), n, LIS_Routing_State(n,k))
       endif
    endif
  end subroutine LIS_routing_DASetStateVar

!BOP
!
!ROUTINE: LIS_routing_DAScaleStateVar
! \label{LIS_routing_DAScaleStateVar}
!
! !INTERFACE:
  subroutine LIS_routing_DAScaleStateVar(n,k)

! !ARGUMENTS:
    integer                :: n
    integer                :: k
!
! !DESCRIPTION:
! 
!  This interface invokes the land model to allow the scaling
!  of the state vector variables
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          call routingdascalestatevar(trim(LIS_rc%routingmodel)//"+"//&
               trim(LIS_rc%daset(k))//char(0), n, LIS_Routing_State(n,k))
       endif
    endif
  end subroutine LIS_routing_DAScaleStateVar

!BOP
!
!ROUTINE: LIS_routing_DAScaleStateVar
! \label{LIS_routing_DAScaleStateVar}
!
! !INTERFACE:
  subroutine LIS_routing_DADescaleStateVar(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
!
! !DESCRIPTION:
! 
!  This interface invokes the land model to allow the descaling
!  of the state vector variables
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          call routingdadescalestatevar(trim(LIS_rc%routingmodel)//"+"//&
               trim(LIS_rc%daset(k))//char(0), n, LIS_Routing_State(n,k), &
               LIS_Routing_Incr_State(n,k))
       endif
    endif
  end subroutine LIS_routing_DADescaleStateVar

!BOP
!
!ROUTINE: LIS_routing_DAUpdateState
! \label{LIS_routing_DAUpdateState}
!
! !INTERFACE:
  subroutine LIS_routing_DAUpdateState(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
!
! !DESCRIPTION:
! 
!  This interface invokes the land model to enable the update of 
!  of the state vector variables based on the increments from the
!  DA algorithm. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          call routingdaupdatestate(trim(LIS_rc%routingmodel)//"+"//&
               trim(LIS_rc%daset(k))//char(0), n, LIS_Routing_State(n,k), &
               LIS_Routing_Incr_State(n,k))
       endif
    endif
  end subroutine LIS_routing_DAUpdateState

!BOP
!
!ROUTINE: LIS_routing_DAQCState
! \label{LIS_routing_DAQCState}
!
! !INTERFACE:
  subroutine LIS_routing_DAQCState(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k    
!
! !DESCRIPTION:
! 
!  This interface invokes the land model to enable the QC of
!  state vector used in DA
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          call routingdaqcstate(trim(LIS_rc%routingmodel)//"+"//&
               trim(LIS_rc%daset(k))//char(0),n,LIS_Routing_State(n,k))
       endif
    endif
  end subroutine LIS_routing_DAQCState


!BOP
!
!ROUTINE: LIS_routing_DAextractStateVector
! \label{LIS_routing_DAextractStateVector}
!
! !INTERFACE:
  subroutine LIS_routing_DAextractStateVector(n,k,state_size,stvar)

! !ARGUMENTS:
    integer                :: n
    integer                :: k
    integer                :: state_size
    real                   :: stvar(LIS_rc%nstvars(k),state_size)
!
! !DESCRIPTION:
! 
!  This interface extracts the state vector variables from the 
!  ESMF state object. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]     index of the nest
!   \item[k]     index of the data assimilation instance   
!   \item[stvar] state vector from the Routing state
!  \end{description}
!EOP
    integer                :: status
    integer                :: v,t
    character*100,    allocatable     :: routing_state_objs(:)
    type(ESMF_Field)                  :: routing_field(LIS_rc%nstvars(k))
    real,         pointer             :: stdata(:)
    real,         pointer             :: stincrdata(:)

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          allocate(routing_state_objs(LIS_rc%nstvars(k)))
          
          call ESMF_StateGet(LIS_Routing_State(n,k),itemNameList=routing_state_objs,&
               rc=status)
          call LIS_verify(status, &
               "ESMF_StateGet failed in enkf_increments")
          
          do v=1,LIS_rc%nstvars(k)
             call ESMF_StateGet(LIS_Routing_State(n,k),trim(routing_state_objs(v)),&
                  routing_field(v),rc=status)
             call LIS_verify(status, &
                  "ESMF_StateGet failed in enkf_increments")
             
             call ESMF_FieldGet(routing_field(v),localDE=0, farrayPtr=stdata,rc=status)
             call LIS_verify(status,&
                  "ESMF_FieldGet failed in enkf_increments")
             
             
             do t=1,state_size
                
                stvar(v,t) = stdata(t)     
             enddo
          enddo
          deallocate(routing_state_objs)
       endif
    endif
  end subroutine LIS_routing_DAextractStateVector
   
!BOP
!
!ROUTINE: LIS_routing_DAgetFreshIncrementsStatus
! \label{LIS_routing_DAgetFreshIncrementsStatus}
!
! !INTERFACE: 
  subroutine  LIS_routing_DAgetFreshIncrementsStatus(n,k,setStatus)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
    logical                :: setStatus
!
! !DESCRIPTION:
! 
!  This interface gets the 'fresh increments status' attribute of
!  the Routing increments object. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]     index of the nest
!   \item[k]     index of the data assimilation instance   
!   \item[setStatus] fresh increments status attribute
!  \end{description}
!EOP
    integer                :: status

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          call ESMF_AttributeGet(LIS_Routing_Incr_State(n,k),&
               "Fresh Increments Status", setstatus, rc=status)
          call LIS_verify(status,&
               'ESMF_AttributeSet: Fresh Increments Status failed in enkf_increments')   
       endif
    endif

  end subroutine LIS_routing_DAgetFreshIncrementsStatus

!BOP
!
!ROUTINE: LIS_routing_DAsetFreshIncrementsStatus
! \label{LIS_routing_DAsetFreshIncrementsStatus}
!
! !INTERFACE:  
  subroutine  LIS_routing_DAsetFreshIncrementsStatus(n,k,setStatus)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
    logical                :: setStatus
!
! !DESCRIPTION:
! 
!  This interface sets the 'fresh increments status' attribute in
!  the Routing increments object. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]     index of the nest
!   \item[k]     index of the data assimilation instance   
!   \item[setStatus] fresh increments status attribute
!  \end{description}
!EOP
    integer                :: status

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          call ESMF_AttributeSet(LIS_Routing_Incr_State(n,k),&
               "Fresh Increments Status", setstatus, rc=status)
          call LIS_verify(status,&
               'ESMF_AttributeSet: Fresh Increments Status failed in enkf_increments')   
       endif
    endif
  end subroutine LIS_routing_DAsetFreshIncrementsStatus

!BOP
!
!ROUTINE: LIS_routing_DAsetAnlysisUpdates
! \label{LIS_routing_DAsetAnlysisUpdates}
!
! !INTERFACE: 
  subroutine LIS_routing_DAsetAnlysisUpdates(n,k,state_size,stvar,stincr)

! !ARGUMENTS:
    integer                :: n
    integer                :: k
    integer                :: state_size
    real                   :: stvar(LIS_rc%nstvars(k),state_size)
    real                   :: stincr(LIS_rc%nstvars(k),state_size)

!
! !DESCRIPTION:
! 
!  This interface extracts the variables from the state vector
!  and state increments vector objects. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]      index of the nest
!   \item[k]      index of the data assimilation instance   
!   \item[stvar]  state vector variables
!   \item[stincr] state increments vector variables
!  \end{description}
!EOP

    integer                :: status
    integer                :: v,t
    character*100,    allocatable     :: routing_state_objs(:)
    type(ESMF_Field)                  :: routing_field(LIS_rc%nstvars(k))
    type(ESMF_Field)                  :: routing_incr_field(LIS_rc%nstvars(k))
    real,         pointer             :: stdata(:)
    real,         pointer             :: stincrdata(:)

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          allocate(routing_state_objs(LIS_rc%nstvars(k)))
          
          call ESMF_StateGet(LIS_Routing_State(n,k),itemNameList=routing_state_objs,&
               rc=status)
          call LIS_verify(status, &
               "ESMF_StateGet failed in enkf_increments")
          
          do v=1,LIS_rc%nstvars(k)
             call ESMF_StateGet(LIS_Routing_State(n,k),trim(routing_state_objs(v)),&
                  routing_field(v),rc=status)
             call LIS_verify(status, &
                  "ESMF_StateGet failed in enkf_increments")
             
             call ESMF_StateGet(LIS_Routing_Incr_State(n,k),trim(routing_state_objs(v)),&
                  routing_incr_field(v),rc=status)
             call LIS_verify(status, &
                  "ESMF_StateGet failed in enkf_increments")
             
             call ESMF_FieldGet(routing_field(v),localDE=0, farrayPtr=stdata,rc=status)
             call LIS_verify(status,&
                  "ESMF_FieldGet failed in enkf_increments")
             
             call ESMF_FieldGet(routing_incr_field(v),localDE=0,farrayPtr=stincrdata,&
                  rc=status)
             call LIS_verify(status, &
                  'ESMF_FieldGet failed in enkf_increments')
             
             do t=1,state_size
                stdata(t) =  stvar(v,t)
                stincrdata(t) = stincr(v,t)
             enddo
             
          enddo
          
          deallocate(routing_state_objs)
       endif
    endif
  end subroutine LIS_routing_DAsetAnlysisUpdates

!BOP
!
!ROUTINE: LIS_routing_DAmapTileSpaceToObsSpace
! \label{LIS_routing_DAmapTileSpaceToObsSpace}
!
! !INTERFACE:
  subroutine LIS_routing_DAmapTileSpaceToObsSpace(n,k,tileid,st_id,en_id)

! !DESCRIPTION: 
! This routine derives the observation space location that maps to the
! input tile space of a given patch space
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [k]
!     index of the DA instance
!   \item [tileid]
!     location in the tile space 
!   \item [st\_id]
!     starting index of the observation space location
!   \item [en\_id]
!     ending index of the observation space location
!  \end{description} 
!EOP

    integer                         :: n
    integer                         :: k
    integer                         :: tileid
    integer                         :: st_id
    integer                         :: en_id

    real                            :: lat, lon, col,row
    integer                         :: gid,c,r

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          
          c = LIS_routing(n)%tile(tileid)%col
          r = LIS_routing(n)%tile(tileid)%row
          
          gid = -1
          if(c.ge.1.and.c.le.LIS_rc%obs_lnc(k).and.&
               r.ge.1.and.r.le.LIS_rc%obs_lnr(k)) then 
             gid = LIS_obs_domain(n,k)%gindex(c,r)
             
          endif
          call ij_to_latlon(LIS_obs_domain(n,k)%lisproj,real(c),real(r),&
               lat,lon)
          st_id = gid
          en_id = gid

#if 0 
          lat = LIS_domain(n)%grid(LIS_domain(n)%gindex( & 
               LIS_surface(n,LIS_rc%lsm_index)%tile(tileid)%col,&
               LIS_surface(n,LIS_rc%lsm_index)%tile(tileid)%row))%lat
          lon = LIS_domain(n)%grid(LIS_domain(n)%gindex( & 
               LIS_surface(n,LIS_rc%lsm_index)%tile(tileid)%col,&
               LIS_surface(n,LIS_rc%lsm_index)%tile(tileid)%row))%lon
          
          call latlon_to_ij(LIS_obs_domain(n,k)%lisproj,lat,lon,&
               col,row)
          c = nint(col)
          r = nint(row)
          
          gid = -1
          if(c.ge.1.and.c.le.LIS_rc%obs_lnc(k).and.&
               r.ge.1.and.r.le.LIS_rc%obs_lnr(k)) then 
             gid = LIS_obs_domain(n,k)%gindex(c,r)
             
          endif
          call ij_to_latlon(LIS_obs_domain(n,k)%lisproj,real(c),real(r),&
               lat,lon)
          st_id = gid
          en_id = gid
#endif
       endif
    endif

  end subroutine LIS_routing_DAmapTileSpaceToObsSpace

!BOP
!
!ROUTINE: LIS_routing_DAgetStateVarNames
! \label{LIS_routing_DAgetStateVarNames}
!
! !INTERFACE:
  subroutine LIS_routing_DAgetStateVarNames(n,k,stateNames)
! !ARGUMENTS:
    integer            :: n
    integer            :: k
    character(len=*)   :: stateNames(LIS_rc%nstVars(k))
! 
! !DESCRIPTION:
! 
!  This interface extracts the variable names from the state vector object
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]      index of the nest
!   \item[k]      index of the data assimilation instance   
!   \item[stateNames  state vector names
!  \end{description}
!EOP
    integer            :: status

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          call ESMF_StateGet(LIS_Routing_State(n,k),itemNameList=stateNames,&
               rc=status)
       endif
    endif
  end subroutine LIS_routing_DAgetStateVarNames


!BOP
!
!ROUTINE: LIS_routing_getlatlons
! \label{LIS_routing_getlatlons}
!
! !INTERFACE:
  subroutine LIS_routing_getlatlons(n,k,state_size,lats,lons)
! !ARGUMENTS:
    integer             :: n    
    integer             :: k
    integer             :: state_size
    real                :: lats(state_size)
    real                :: lons(state_size)
! 
! !DESCRIPTION:
! 
!  This routine extracts the lat/lon values corresponding to the
!  state vector space used by the Routing model. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]      index of the nest
!   \item[k]      index of the data assimilation instance   
!   \item[stateNames  state vector names
!  \end{description}
!EOP

    integer             :: i,c,r,gid

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%routing_DAinst_valid(k)) then 
          do i=1,state_size
             c = LIS_routing(n)%tile(i)%col
             r = LIS_routing(n)%tile(i)%row
             gid = LIS_routing(n)%gindex(c,r)

             lats(i) = LIS_routing(n)%grid(gid)%lat
             lons(i) = LIS_routing(n)%grid(gid)%lon
             
          enddo
       endif
    endif

  end subroutine LIS_routing_getlatlons
        
  subroutine LIS_routing_DAobsTransform(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%Routing_DAinst_valid(k)) then
          call routingdaobstransform(trim(LIS_rc%routingmodel)//"+"//&
               trim(LIS_rc%daset(k))//char(0), n, LIS_OBS_State(n,k))
       endif
    endif
  end subroutine LIS_routing_DAobsTransform

  subroutine LIS_routing_DAmapObsToRouting(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%Routing_DAinst_valid(k)) then
          call routingdamapobstorouting(trim(LIS_rc%routingmodel)//"+"//&
               trim(LIS_rc%daset(k))//char(0), n, k, LIS_OBS_State(n,k),&
               LIS_Routing_Incr_State(n,k))

       endif
    endif
  end subroutine LIS_routing_DAmapObsToRouting

  subroutine LIS_routing_DAqcObsState(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k

    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%Routing_DAinst_valid(k)) then

          call routingdaqcobsstate(trim(LIS_rc%routingmodel)//"+"&
               //trim(LIS_rc%daset(k))//char(0),n, k, LIS_OBS_State(n,k))

       endif
    endif
  end subroutine LIS_routing_DAqcObsState

!BOP
!
!ROUTINE: LIS_routing_DAgetStateSpaceSize
! \label{LIS_routing_DAgetStateSpaceSize}
!
! !INTERFACE:
  subroutine LIS_routing_DAgetStateSpaceSize(n,k,size)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
    integer                :: size
!
! !DESCRIPTION:
! 
!  This routine returns the spatial size of the state vector
!  used in DA
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!  \end{description}
!EOP
    if(LIS_rc%routingmodel.ne."none" .and. LIS_rc%routingmodel.ne."RAPID router") then
       if(LIS_rc%Routing_DAinst_valid(k)) then
          call routingdagetstatespacesize(trim(LIS_rc%routingmodel)//"+"//&
               trim(LIS_rc%daset(k))//char(0), n, size)
       endif
    endif
  end subroutine LIS_routing_DAgetStateSpaceSize

end module LIS_routingMod
