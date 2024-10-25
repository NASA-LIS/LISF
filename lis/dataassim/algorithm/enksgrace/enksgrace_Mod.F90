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
module enksgrace_Mod
!BOP
!
! !MODULE: enksgrace_Mod
!
! !DESCRIPTION:
!   This module contains interfaces and subroutines that control
!   the incorporation of a data set using the ensemble kalman filter
!   (EnKF) method, into a land surface model. 
!
!  The EnKF algorithm is based on the work of Rolf Reichle at the NASA
!  Global Modeling and Assimilation Office (GMAO) at the NASA GSFC.
!  
! !REVISION HISTORY:
!   27 Feb 2005: Sujay Kumar; Initial Specification
! 
! !USES: 
  use ESMF
  use enks_types
  use LIS_coreMod
  use LIS_DAobservationsMod
  use LIS_logMod
  use LIS_lsmMod
  use LIS_fileIOMod
  use LIS_timeMgrMod
  use LIS_historyMod
  use enks_general
  use enks_matrix_functions
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: enksgrace_init  ! Initialization for EnKF 
  public :: enksgrace_setup
  public :: enksgrace_increments ! compute analysis increments
  public :: enksgrace_update ! apply the analysis increments
  public :: enksgrace_diagnostics ! write EnKF related diagnostics
  public :: enksgrace_final ! Finalization for EnKF
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: enksgrace_struc ! data structure containing EnKF diagnostics
!EOP

  type, public ::  enksgrace_dec
     logical     :: fileOpen
     real, allocatable :: innov(:)
     real, allocatable :: forecast_var(:) !HPHt
     real, allocatable :: norm_innov(:)
     real, allocatable :: k_gain(:,:)

  end type enksgrace_dec
!EOP  

  type(enksgrace_dec), allocatable :: enksgrace_struc(:,:)

contains

!BOP
! 
! !ROUTINE: enksgrace_init
! \label{enksgrace_init}
!  
! !INTERFACE: 
  subroutine enksgrace_init()
! !USES: 


!
! !DESCRIPTION: 
!  This method performs the required initializations for the 
!  GMAO EnKF method. The method reads the runtime settings from 
!  the LIS configuration file. 
! 
!EOP
    allocate(enksgrace_struc(LIS_rc%nnest, LIS_rc%ndas))
   
  end subroutine enksgrace_init
!BOP
! 
! !ROUTINE: enksgrace_setup
! \label{enksgrace_setup}
!  
! !INTERFACE: 
  subroutine enksgrace_setup(k)
! !USES: 

!
! !DESCRIPTION: 
!  This method performs the required initializations for the 
!  GMAO EnKF method. The method reads the runtime settings from 
!  the LIS configuration file. 
! 
!EOP
    integer                      :: n
    integer                      :: k
    integer                      :: status
    integer                      :: Nobjs,N_obs_size

    do n=1,LIS_rc%nnest
       if(LIS_rc%nensem(n).le.1) then
          write(LIS_logunit,*) '[ERR] Please set number of ensembles '
          write(LIS_logunit,*) '[ERR] to greater than 1 for ENKS',LIS_rc%nensem(n)
          call LIS_endrun
       endif
       enksgrace_struc(n,k)%fileOpen = .false.
    enddo
    if(mod(1800.0,LIS_rc%ts).ne.0) then 
       write(LIS_logunit,*) '[ERR] Please make sure that 30 minutes is a '
       write(LIS_logunit,*) '[ERR] multiple of the LIS clock timestep in GRACE-DA mode'
       call LIS_endrun()
    endif

    do n=1,LIS_rc%nnest
       call ESMF_StateGet(LIS_OBS_State(n,k),itemCount=Nobjs,rc=status)
       call LIS_verify(status, 'ESMF_StateGet failed in enkf_Mod')
       
       call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Number Of Observations",&
            value=N_obs_size,rc=status)
       call LIS_verify(status, 'ESMF_AttributeGet failed in enkf_Mod')
       
       if(LIS_rc%winnov(k).eq.1) then
          allocate(enksgrace_struc(n,k)%norm_innov(Nobjs*N_obs_size))
          allocate(enksgrace_struc(n,k)%innov(Nobjs*N_obs_size))
          allocate(enksgrace_struc(n,k)%forecast_var(Nobjs*N_obs_size))
          allocate(enksgrace_struc(n,k)%k_gain(LIS_rc%npatch(n,LIS_rc%lsm_index), &
               LIS_rc%nstvars(k)))
       endif
    enddo
  end subroutine enksgrace_setup

!BOP
! 
! !ROUTINE: enksgrace_increments
! \label{enksgrace_increments}
!  
! !INTERFACE: 
  subroutine enksgrace_increments(n,k)
! !USES: 

! !ARGUMENTS: 
    integer, intent(IN)    :: n 
    integer, intent(IN)    :: k
! 
! !DESCRIPTION: 
!  This routine computes the analysis increments for Enks. The state variable
!  increments are computed and stored int the Increment State. 
! 
!  The arguments and variables are: 
!  \begin{description}
!   \item[n]         index of the nest 
!   \item[data\_status]  flag indicating if the observations are new    
!   \item[N\_state]       Number of state prognostic variables
!   \item[N\_ens]         Number of ensembles
!   \item[N\_obs\_size]   total number of observations   
!   \item[N\_selected\_obs] number of selected observations   
!   \item[st\_id,en\_id]  indices of the selected observations  
!   \item[stvar]          state prognostic variables
!   \item[Obs\_pert]      observation perturbations     
!   \item[Obs\_cov]       observations error covariances   
!   \item[Observations]    observations object used in the GMAO routines
!  \end{description}
! 
!
!  The methods invoked are: 
!  \begin{description}
!   \item[generateObservations]\ref{generateObservations_enksgrace} \newline
!    obtain the observations
!   \item[lsmdagetobspred](\ref{lsmdagetobspred}) \newline
!    obtain model's estimate of the observations
!   \item[getObsPert](\ref{getObsPert_enksgrace}) \newline
!    obtain the observation perturbations
!   \item[generateObsparam](\ref{generateObsparam_enksgrace}) \newline
!    generate the 'obsparam' (metadata for observations) \newline
!   \item[lsmdagetstatevar]\ref{lsmdagetstatevar} \newline
!    obtain the specified LSM prognostic variables. 
!   \item[assemble\_obs\_cov](\ref{assemble_obs_cov_enksgrace})
!    assembles the observation error covariance
!   \item{getSelectedObsNumber}(\ref{getselectedobsnumber}) \newline
!    obtain the number of selected observations for each 
!    modeling point.
!   \item[enks\_analysis](\ref{enks_analysis}) \newline
!    apply the Enks filter to obtain the prognostic variable
!    increments. 
!   \item[row\_variance](\ref{row_variance_enksgrace}) \newline
!    computes the row variance HPH' 
!   \item[lsmdascalestatvar](\ref{lsmdascalestatevar}) \newline
!    scales the state variables for computational stability \newline
!   \item[lsmdadescalestatevar](\ref{lsmdadescalestatevar}) \newline
!    descales the state variables to the original state \newline
!  \end{description}
! 
!EOP
    logical                       :: data_status
    integer                       :: status
    integer                       :: Nobjs
    integer                       :: Nobs
    integer                       :: N_obs_size
    integer                       :: N_selected_obs
    integer                       :: N_ens
    integer                       :: N_state
    type(obs_type), allocatable       :: Observations(:)
    type(obs_type), allocatable       :: obs_da(:)
    type(obs_param_type), allocatable :: obs_param(:)
    real,         allocatable         :: Obs_pred(:,:)
    real,         allocatable         :: obspred_da(:,:)
    real,         allocatable         :: Obs_pert(:,:)
    real,         allocatable         :: Obs_cov(:,:)
    real,         allocatable         :: obspert_da(:,:)
    integer                       :: i,v,tileid
    integer                       :: st_id, en_id, sid,eid
    character*100,    allocatable     :: lsm_state_objs(:)
    type(ESMF_Field)              :: lsm_field(LIS_rc%nstvars(k))
    type(ESMF_Field)              :: lsm_incr_field(LIS_rc%nstvars(k))
    real,         allocatable         :: stvar(:,:)
    real,         pointer         :: stdata(:)
    real,         pointer         :: stincrdata(:)
    real,         allocatable         :: state_incr(:,:)
    real                          :: innov,std_innov(1)
    integer                       :: kk,m
    logical                       :: assim
    integer                       :: gid, t,p
    real                          :: avg_factor
!    integer                       :: days(12), mindex !BZ    
!    data days /31,28,31,30,31,30,31,31,30,31,30,31/  !BZ
!----------------------------------------------------------------------------
!  Check if the observation state is updated or not. If it is updated,
!  increments are calculated (but calc is decoupled from application in ENKS). 
!----------------------------------------------------------------------------

    call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Data Update Status",&
         value=data_status,rc=status)
    call LIS_verify(status)

    call ESMF_AttributeSet(LIS_LSM_Incr_State(n,k),&
         "Fresh Increments Status", .false., rc=status)
    call LIS_verify(status)    

    if(data_status) then  
       write(LIS_logunit,*) '[INFO] Assimilating Observations using GMAO EnKS'

       call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Data averaging factor",&
            value=avg_factor,rc=status)
       call LIS_verify(status, 'Data averaging factor: attribute not found')

       N_state = LIS_rc%nstvars(k)
       N_ens = LIS_rc%nensem(n)
!----------------------------------------------------------------------------
! It is assumed that the obs_state in this subroutine is a superset of
! the required observations for each point in the processor's modeling
! domain. 
!----------------------------------------------------------------------------
       call ESMF_StateGet(LIS_OBS_State(n,k),itemCount=Nobjs,rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Number Of Observations",&
            value=N_obs_size,rc=status)
       call LIS_verify(status)

       Nobs = Nobjs*N_obs_size
       allocate(Observations(Nobs))
       call generateObservations(n, k, Nobjs, Nobs, LIS_OBS_State(n,k), &
            LIS_OBS_Pert_State(n,k),Observations)

!----------------------------------------------------------------------------
!  Retrieve Obs_pred : model's estimate of the observations
!----------------------------------------------------------------------------

       allocate(Obs_pred(Nobs,N_ens))      
!       write(LIS_logunit,*)'calling to getObsPred'
!       write(LIS_logunit,*)'Nobs=',Nobs,'N_ens=',N_ens,'n=',n
!       write(LIS_logunit,*)'LIS_rc%lsm = ',LIS_rc%lsm,'LIS_rc%daset(k) = ',LIS_rc%daset(k)
       call lsmdagetobspred(&
            trim(LIS_rc%lsm)//"+"//trim(LIS_rc%daset(k))//char(0),&
            n, k, Obs_pred)
!----------------------------------------------------------------------------
!  Retrieve Obs_pert : observation perturbations
!---------------------------------------------------------------------------- 
       allocate(Obs_pert(Nobs,N_ens))
       call getObsPert(LIS_OBS_Pert_State(n,k),N_obs_size,&
            Nobs, N_ens, Obs_pert)

!     Implemented Ryu et al(2009) Ens_Pert Bias Removal Approach: 
!       if(LIS_rc%pert_bias_corr == 1 ) then
!     !-- Remove perturbation to observation for N-ensemble member ... 
!          do p=1,Nobs
!             Obs_pert(p,N_ens) = 0.
!          end do
!       endif
!----------------------------------------------------------------------------
!  Assemble observation covariances. 
!----------------------------------------------------------------------------
       allocate(obs_param(LIS_rc%nobtypes(k)))
       call generateObsparam(Nobjs, LIS_OBS_Pert_State(n,k),obs_param)
       
!----------------------------------------------------------------------------
! retrieve the state variables
!----------------------------------------------------------------------------
       allocate(stvar(N_state,LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(state_incr(N_state,LIS_rc%npatch(n,LIS_rc%lsm_index)))

       call lsmdagetstateVar(&
            trim(LIS_rc%lsm)//"+"//trim(LIS_rc%daset(k))//char(0),&
            n, LIS_LSM_State(n,k))

       call lsmdascalestateVar(&
            trim(LIS_rc%lsm)//"+"//trim(LIS_rc%daset(k))//char(0),&
            n, LIS_LSM_State(n,k))
       
       allocate(lsm_state_objs(LIS_rc%nstvars(k)))

       call ESMF_StateGet(LIS_LSM_State(n,k),itemNameList=lsm_state_objs,&
            rc=status)
       call LIS_verify(status)
       do v=1,LIS_rc%nstvars(k)
          call ESMF_StateGet(LIS_LSM_State(n,k),trim(lsm_state_objs(v)),&
               lsm_field(v),rc=status)
          call LIS_verify(status)

          call ESMF_StateGet(LIS_LSM_Incr_State(n,k),trim(lsm_state_objs(v)),&
               lsm_incr_field(v),rc=status)
          call LIS_verify(status)

          call ESMF_FieldGet(lsm_field(v),localDE=0, farrayPtr=stdata,rc=status)
          call LIS_verify(status)

          call ESMF_FieldGet(lsm_incr_field(v),localDE=0, farrayPtr=stincrdata,&
               rc=status)
          call LIS_verify(status)

          do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
             stvar(v,t) = stdata(t)     
             state_incr(v,t) = stdata(t)
          enddo
       enddo
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
          tileid = (i-1)*LIS_rc%nensem(n)+1

          call LIS_lsm_DAmapTileSpaceToObsSpace(n, k, &
               tileid, st_id, en_id)

!          call getSelectedObsNumber(trim(LIS_rc%daset(k))//char(0),n,&
!               LIS_surface(n,LIS_rc%lsm_index)%tile(tileid)%index,st_id,en_id)
          if(st_id.lt.0.or.en_id.lt.0) then 
             assim = .false. 
          else          
             N_selected_obs = (en_id-st_id+1)*Nobjs
             assim = .true.
             do kk=st_id,en_id
                assim = assim .and.Observations(kk)%assim
             enddo

             allocate(obs_da(N_selected_obs))
             allocate(obspred_da(N_selected_obs,N_ens))
             allocate(obspert_da(N_selected_obs,N_ens))
             allocate(obs_cov(N_selected_obs, N_selected_obs))
             
             kk = 1
             do while(kk.le.N_selected_obs)
                sid = st_id + (kk-1)*N_obs_size
                eid = en_id + (kk-1)*N_obs_size
                obs_da(kk:kk+(en_id-st_id))       = Observations(sid:eid)
                obspred_da(kk:kk+(en_id-st_id),:) = Obs_pred(sid:eid,:)
                obspert_da(kk:kk+(en_id-st_id),:) = Obs_pert(sid:eid,:)
                kk = kk+(en_id-st_id+1)
             enddo

             call assemble_obs_cov(LIS_rc%nobtypes(k), N_selected_obs, &
                  obs_param, obs_da,Obs_cov)             
          endif
          if(assim) then  
             
!for serial run
!             if(i.eq.113) then 
!                print*, 'obs', obs_da
!                print*, 'obspred',obspred_da
!                print*, 'obspert ',obspert_da
!                print*, 'obscov',obs_cov
!                print*, 'state', state_incr(:, ((i-1)*N_ens+1):((i-1)*N_ens+N_ens))
!             endif
!for 4 process run
!             if(i.eq.66.and. LIS_localPet.eq.0) then 
!             if(LIS_localPet.eq.0) then 
!                if(LIS_domain(n)%grid(i)%lat.eq.42.125.and.&
!                     LIS_domain(n)%grid(i)%lon.eq.-88.375) then 
!                   print*, 'i ',i, LIS_domain(n)%grid(i)%lat, LIS_domain(n)%grid(i)%lon
!                   print*, 'obs', obs_da
!                   print*, 'obspred',obspred_da
!                   print*, 'obspert ',obspert_da
!                   print*, 'obscov',obs_cov
!                   print*, 'state', state_incr(:, ((i-1)*N_ens+1):((i-1)*N_ens+N_ens))
!                endif
!             endif
             call enks_analysis(gid,N_state,N_selected_obs, N_ens, &
                  obs_da, & 
                  obspred_da, &
                  obspert_da, &
                  obs_cov, &
                  state_incr(:, ((i-1)*N_ens+1):((i-1)*N_ens+N_ens)))  
          else
             state_incr(:,(i-1)*N_ens+1:(i-1)*N_ens+N_ens) = 0.0            
          endif
          if(.not.(st_id.lt.0.or.en_id.lt.0)) then 
             deallocate(obs_da)
             deallocate(obspred_da)
             deallocate(obspert_da)
             deallocate(obs_cov)
          endif
       enddo
       call ESMF_AttributeSet(LIS_LSM_Incr_State(n,k),&
            "Fresh Increments Status", .true., rc=status)
       call LIS_verify(status)
       
       if(LIS_rc%winnov(k).eq.1) then 
          do i=1,Nobs
             if(Observations(i)%assim) then
                innov = Observations(i)%value - &
                     sum(Obs_pred(i,:))/real(LIS_rc%nensem(n))
             ! compute diag(HPHt), put it into std_innov
                call row_variance(1,LIS_rc%nensem(n),Obs_pred(i,:),std_innov(1))
             !  add diag (R)
                enksgrace_struc(n,k)%forecast_var(i) = std_innov(1)
                std_innov = std_innov+(Observations(i)%std)**2          
                std_innov = sqrt(std_innov)
                enksgrace_struc(n,k)%norm_innov(i) = innov/std_innov(1)
                enksgrace_struc(n,k)%innov(i) = innov
             else
                enksgrace_struc(n,k)%norm_innov(i) = LIS_rc%udef
                enksgrace_struc(n,k)%innov(i) = LIS_rc%udef
                enksgrace_struc(n,k)%forecast_var(i) = LIS_rc%udef
             endif
          enddo
       endif
!----------------------------------------------------------------------------
! Updating LIS_LSM_State and LSM_Incr_State
!----------------------------------------------------------------------------
       do v=1,LIS_rc%nstvars(k)
          call ESMF_FieldGet(lsm_field(v),localDE=0,farrayPtr=stdata,rc=status)
          call LIS_verify(status)
          
          call ESMF_FieldGet(lsm_incr_field(v),localDE=0,farrayPtr=stincrdata,&
               rc=status)
          call LIS_verify(status)

          do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

             stdata(t) =  stvar(v,t)
             !BZ divide by number of days in month
!             if((mod(LIS_rc%yr,4) .eq. 0 .and. &
!                  mod(LIS_rc%yr, 100).ne.0) &!leap year
!                  .or.(mod(LIS_rc%yr,400) .eq.0)) then
!                days(2) = 29
!             else
!                days(2) = 28
!             endif
!             mindex = LIS_rc%mo

!             stincrdata(t) = state_incr(v,t)/days(mindex)  ! BZ simple smoother
             stincrdata(t) = state_incr(v,t)/avg_factor

          enddo
       enddo
       call lsmdadescalestateVar(trim(LIS_rc%lsm)//"+"&
            //trim(LIS_rc%daset(k))//char(0),&
            n, LIS_LSM_State(n,k), &
            LIS_LSM_Incr_State(n,k))
!----------------------------------------------------------------------------
! Cleanup
!----------------------------------------------------------------------------
       deallocate(lsm_state_objs)
       deallocate(obs_param)
       deallocate(stvar)
       deallocate(State_incr)
       deallocate(Observations)
       deallocate(Obs_pred)
       deallocate(Obs_pert)
    end if
    
  end subroutine enksgrace_increments

!BOP
! 
! !ROUTINE: enksgrace_update
! \label{enksgrace_update}
!
! !INTERFACE: 
  subroutine enksgrace_update(n,k)
! !USES: 

    implicit none
    
    integer,       intent(in)    :: n
    integer,       intent(in)    :: k
! 
! !DESCRIPTION: 
!  This routine updates the model prognostics using the analysis
!  increments computed earlier. 
!
!EOP
    character*100,    allocatable     :: lsm_state_objs(:)
    integer                       :: status,m
!BZ    type(ESMF_Logical)            :: fresh_incr
! BZ add following
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    character(len=LIS_CONST_PATH_LEN) :: GRACEobsdir
    logical             :: file_exists
    character(len=LIS_CONST_PATH_LEN) :: name
    character (len=4) :: fyr
    character (len=2) :: fmo    

    allocate(lsm_state_objs(LIS_rc%nstvars(k)))
    
    call ESMF_StateGet(LIS_LSM_State(n,k),itemNameList=lsm_state_objs,&
         rc=status)
    call LIS_verify(status)
    

    call ESMF_AttributeGet(LIS_OBS_State(n,k),name="File Status",&
         value=file_exists,rc=status)

! BZ apply smoother logic here [using simple daily applications for now]
    if(file_exists.and.(LIS_rc%DAincrMode(n).eq.1).and.(LIS_rc%hr.eq.12)&
         .and.(LIS_rc%mn.eq.30)) then

      if(LIS_rc%incroption(k).eq.0) then !include analysis increments
         call lsmdaupdatestate(trim(LIS_rc%lsm)//"+"//&
              trim(LIS_rc%daset(k))//char(0),&
              n, LIS_LSM_State(n,k), &
              LIS_LSM_Incr_State(n,k))
!----------------------------------------------------------------------
!  Update the LSM's state variables
!----------------------------------------------------------------------       
    !BZ     call qclsmstate(LIS_rc%lsm,LIS_rc%daset(k),n,LIS_LSM_State(n,k))
         call lsmdasetstatevar(trim(LIS_rc%lsm)//"+"//&
              trim(LIS_rc%daset(k))//char(0),&
              n, LIS_LSM_State(n,k))
      endif

!----------------------------------------------------------------------
!  Cleanup
!---------------------------------------------------------------------
       call ESMF_AttributeSet(LIS_OBS_State(n,k),"Data Assimilate Status",&
         .true., rc=status)
       call LIS_verify(status)

       write(LIS_logunit,*) '[INFO] Finished assimilating Observations using GMAO EnKS'
    else
       call ESMF_AttributeSet(LIS_OBS_State(n,k),"Data Assimilate Status",&
            .false., rc=status)
       call LIS_verify(status)
    endif

end subroutine enksgrace_update


!BOP
! 
! !ROUTINE: enksgrace_diagnostics
!  \label{enksgrace_diagnostics}
! 
! !INTERFACE:
  subroutine enksgrace_diagnostics(n,k)
! !USES:

! !ARGUMENTS:
    integer, intent(IN)    :: n 
    integer, intent(IN)    :: k
!  
! !DESCRIPTION:  
!  This subroutine generates the Enks related diagnostics and outputs
!  it to a file. This includes a text output of selected ensemble
!  members, their mean, standard deviation, observations, and normalized
!  innovations. The frequency of diagnostic outputs can be specified
!  in the LIS configuration file.
!
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!    call to create output directory for DA statistics
!   \item[getLSMvarnames](\ref{getLSMvarnames_enksgrace}) \newline
!    retrieve the names of the LSM prognostic variables
!   \item[pruneVarname](\ref{pruneVarname_enksgrace}) \newline
!    trims the variable name, eliminating white spaces
!   \item[LIS\_create\_stats\_filename](\ref{LIS_create_stats_filename}) \newline
!    creates the filename for statistics 
!   \item[LIS\_create\_innov\_filename](\ref{LIS_create_innov_filename}) \newline
!    creates the name of the innovations file
!   \item[lsmdagetstatevar](\ref{lsmdagetstatevar}) \newline
!    retrieve the lsm state variables
!   \item[getLSMdata](\ref{getLSMData_enksgrace}) \newline
!    unpack the LSM state and retrive the data
!  \end{description}
!EOP

    integer                :: v
    integer                :: status
    logical                :: assim_status
    logical                :: alarmCheck
    character*3            :: fda

    write(fda,'(i3.3)') k
    alarmCheck = LIS_isAlarmRinging(LIS_rc,"LIS DA output "//trim(fda))
       
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    call ESMF_AttributeGet(LIS_OBS_State(n,k),name="Data Assimilate Status",&
         value=assim_status,rc=status)
    
    if(assim_status) then 
!--------------------------------------------------------------------------
! write innovations file with the following entries for all observation 
! types in each data assimilation instance. 
! 
! * Normalized innovations
! * Raw innovations
! * Ensemble spread
! * Analysis increments
! * Forecast variance
!--------------------------------------------------------------------------
       if(.not.enksgrace_struc(n,k)%fileopen.and.LIS_masterproc) then
          
          call LIS_create_output_directory("EnKF")
          
       endif
       
       call writeInnovationOutput(n,k)

    endif

    if(alarmCheck) then 
       if(.not.enksgrace_struc(n,k)%fileopen.and.LIS_masterproc) then
          
          call LIS_create_output_directory("EnKF")
          
       endif

       call writeEnsembleSpread(n,k)

    endif


#endif

  end subroutine enksgrace_diagnostics

!BOP
! 
! !ROUTINE: writeInnovationOutput
! \label{writeInnovationOutput}
!
! !INTERFACE: 
  subroutine writeInnovationOutput(n,k)
! !ARGUMENTS:
    integer,  intent(in)    :: n 
    integer,  intent(in)    :: k 
!
! !DESCRIPTION: 
! 
!  This routine writes the innovation values (observation minus the model 
!  forecast) to an external file. 
! 
!  The arguments and variables are: 
!  \begin{description}
!   \item[n]    index of the nest 
!   \item[k]    index of the observation datastream
!  \end{description}
! 
!
!EOP
    integer                :: ftn
    character(len=LIS_CONST_PATH_LEN) :: innovfile, gainfile, incrfile
    integer                :: shuffle, deflate, deflate_level
    integer                :: dimID(3), ares_Id, ninnov_Id, innov_id
    character*100          :: varname, vardimname, standard_name
    character*2            :: finst
    integer                :: status

    if(LIS_rc%winnov(k).eq.1) then 

       shuffle = 1
       deflate = 1
       deflate_level =9

       if(LIS_masterproc) then 
          call LIS_create_innov_filename(n,k, innovfile,&
               'EnKF')
          
#if (defined USE_NETCDF4)
          status = nf90_create(path=innovfile,cmode=nf90_hdf5,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(innovfile)//&
               ' failed in enksgrace_Mod')
#endif
#if (defined USE_NETCDF3)
          status = nf90_create(path=innovfile,cmode=nf90_clobber,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(innovfile)//&
               ' failed in enksgrace_Mod')
#endif
          
          call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%obs_gnc(k),&
               dimID(1)),'nf90_def_dim for east_west failed in enksgrace_mod')
          call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%obs_gnr(k),&
               dimID(2)),'nf90_def_dim for north_south failed in enksgrace_mod')

          call LIS_verify(nf90_put_att(ftn,&
               NF90_GLOBAL,"missing_value", LIS_rc%udef),&
               'nf90_put_att failed for missing_value in enksgrace_mod')
          
!--------------------------------------------------------------------------
!  Normalized innovations -meta data
!--------------------------------------------------------------------------
          write(unit=finst, fmt='(i2.2)') k
          varname = "ninnov_"//trim(finst)
          vardimname = "ninnov_"//trim(finst)//"_levels"
          standard_name = "Normalized_innovations_for_DA_instance_"//&
               trim(finst)
          
          call LIS_verify(nf90_def_dim(ftn,&
               vardimname,LIS_rc%nobtypes(k),dimId(3)),&
               'nf90_def_dim failed for ninnov_'//trim(finst))
          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=ninnov_Id),&
               'nf90_def_var failed for '//trim(varname)//' in enksgrace_mod')
          
#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               ninnov_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate failed for ninnov in enksgrace_mod')      
#endif
          call LIS_verify(nf90_put_att(ftn,ninnov_Id,&
               "standard_name",standard_name),&
               'nf90_put_att failed for '//trim(standard_name)//' in enksgrace_mod')
          
!--------------------------------------------------------------------------
!  Innovations -meta data
!--------------------------------------------------------------------------
          varname = "innov_"//trim(finst)
          vardimname = "innov_"//trim(finst)//"_levels"
          standard_name = "Innovations_for_DA_instance_"//&
               trim(finst)
          
          call LIS_verify(nf90_def_dim(ftn,&
               vardimname,LIS_rc%nobtypes(k),dimId(3)),&
               'nf90_def_dim failed for innov_'//trim(finst))

          call LIS_verify(nf90_def_var(ftn,varname,&
               nf90_float,&
               dimids = dimID(1:2), varID=innov_Id),&
               'nf90_def_var failed for innov')
          
#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               innov_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for innov failed in enksgrace_mod')             
#endif
          call LIS_verify(nf90_put_att(ftn,innov_Id,&
               "standard_name",standard_name),&
               'nf90_put_att failed for innov in enksgrace_mod')
          
       endif
       
       call LIS_writevar_innov(ftn,n, k, ninnov_id, &
            enksgrace_struc(n,k)%norm_innov)
       call LIS_writevar_innov(ftn,n, k, innov_id, &
            enksgrace_struc(n,k)%innov)

       if(LIS_masterproc) then 
          call LIS_verify(nf90_close(ftn),&
               'nf90_close failed in enksgrace_mod')
       endif

    endif
  end subroutine writeInnovationOutput


!BOP
! 
! !ROUTINE: writeEnsembleSpread
! \label{writeEnsembleSpread}
!
! !INTERFACE: 
  subroutine writeEnsembleSpread(n,k)
!
! !DESCRIPTION: 
!  This routine writes the ensemble spread (standard deviation) 
!  of the model state vector. 
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[n]    index of the nest 
!   \item[k]    index of the observation datastream
!  \end{description}

!EOP
    integer,  intent(in)    :: n 
    integer,  intent(in)    :: k 

    integer                :: ftn 
    integer                :: v
    character(len=LIS_CONST_PATH_LEN) :: spreadfile
    integer                :: shuffle, deflate, deflate_level
    integer                :: dimID(3)
    integer                :: ensspread_id(LIS_rc%nstvars(k))
    character*100          :: varname, vardimname, standard_name
    character*2            :: finst
    integer                :: status
    real, allocatable      :: stvar(:,:)
    character*100,    allocatable     :: lsm_state_objs(:)

    if(LIS_rc%wensems(k).eq.1) then 

       shuffle = 1
       deflate = 1
       deflate_level =9
       
       if(LIS_masterproc) then 
          call LIS_create_daspread_filename(n,k,spreadfile,&
               'EnKF')
          
#if (defined USE_NETCDF4)
          status = nf90_create(path=spreadfile,cmode=nf90_hdf5,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(spreadfile)//&
               ' failed in enksgrace_Mod')
#endif
#if (defined USE_NETCDF3)
          status = nf90_create(path=spreadfile,cmode=nf90_clobber,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(spreadfile)//&
               ' failed in enksgrace_Mod')
#endif
          
          if(LIS_rc%wopt.eq."1d gridspace") then 
             call LIS_verify(nf90_def_dim(ftn,'ngrid',&
                  LIS_rc%glbngrid_red(n),&
                  dimID(1)),'nf90_def_dim for ngrid failed in enksgrace_mod')
          elseif(LIS_rc%wopt.eq."2d gridspace") then 
             call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
                  dimID(1)),'nf90_def_dim for east_west failed in enksgrace_mod')
             call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
                  dimID(2)),'nf90_def_dim for north_south failed in enksgrace_mod')
          endif

          call LIS_verify(nf90_put_att(ftn,&
               NF90_GLOBAL,"missing_value", LIS_rc%udef),&
               'nf90_put_att for missing_value failed in enksgrace_mod')
          
!--------------------------------------------------------------------------
!  Ensemble spread -meta data
!--------------------------------------------------------------------------
          allocate(lsm_state_objs(LIS_rc%nstvars(k)))          
          call ESMF_StateGet(LIS_LSM_State(n,k),itemNameList=lsm_state_objs,&
               rc=status)

          do v = 1, LIS_rc%nstvars(k)
             write(unit=finst, fmt='(i2.2)') k
             varname = "ensspread_"//trim(lsm_state_objs(v))//"_"//trim(finst)
             vardimname = "ensspread_"//trim(lsm_state_objs(v))//&
                  "_"//trim(finst)//"_levels"
             standard_name = "Ensemble_spread_for_DA_instance_"//&
                  trim(lsm_state_objs(v))//"_"//&
                  trim(finst)

             if(LIS_rc%wopt.eq."1d gridspace") then           
                call LIS_verify(nf90_def_var(ftn,varname,&
                     nf90_float,&
                     dimids = dimID(1), varID=ensspread_Id(v)),&
                     'nf90_def_var for ensspread failed in enksgrace_mod')
                
             elseif(LIS_rc%wopt.eq."2d gridspace") then 
                call LIS_verify(nf90_def_var(ftn,varname,&
                     nf90_float,&
                     dimids = dimID(1:2), varID=ensspread_Id(v)),&
                     'nf90_def_var for ensspread failed in enksgrace_mod') 
             endif
          
#if(defined USE_NETCDF4)
             call LIS_verify(nf90_def_var_deflate(ftn,&
                  ensspread_Id(v),&
                  shuffle, deflate, deflate_level),&
                  'nf90_def_var_deflate for ensspread failed in enksgrace_mod')             
#endif
             call LIS_verify(nf90_put_att(ftn,ensspread_Id(v),&
                  "standard_name",standard_name),&
                  'nf90_put_att for ensspread failed in enksgrace_mod')
             call LIS_verify(nf90_enddef(ftn),&
                  'nf90_enddef failed in enksgrace_mod')
          end do
          deallocate(lsm_state_objs)          
       endif
       
       allocate(stvar(LIS_rc%nstvars(k),&
            LIS_rc%npatch(n,LIS_rc%lsm_index)))
       
       
       call lsmdagetstatevar(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k))
       call getLSMdata(LIS_LSM_State(n,k), LIS_rc%nstvars(k), &
            LIS_rc%npatch(n,LIS_rc%lsm_index), stvar)
       
       do v=1,LIS_rc%nstvars(k)
          call LIS_writevar_spread(ftn,n,k,ensspread_id(v), &
               stvar(v,:),v)
       enddo
       
       deallocate(stvar)

       if(LIS_masterproc) then 
          call LIS_verify(nf90_close(ftn),&
               'nf90_close failed in enksgrace_mod')
       endif
    endif

  end subroutine writeEnsembleSpread


!BOP
! 
! !ROUTINE: pruneVarname
! \label{pruneVarname_enksgrace}
!
! !INTERFACE:
  subroutine pruneVarname(varname)
! !ARGUMENTS: 
    character(len=*), intent(INOUT) :: varname
! 
! !DESCRIPTION:
! 
!  This routine generates a filename based on the names of the state
!  prognostic variables
! 
!  \begin{description}
!  \item[varname]  name of the variable filename
!  \end{description}
!EOP
    character*100                   :: temp
    character*1                     :: ftemp(100)
    character*1                     :: fdir(9)
    integer                         :: i,c
    
    
    write(unit=temp,fmt='(A100)') varname
    read(unit=temp,fmt='(100a1)') (ftemp(i),i=1,100)

    write(unit=temp,fmt='(A9)') 'enksgrace/'
    read(unit=temp,fmt='(9a1)') (fdir(i),i=1,9)
    
    do i=100,1,-1
       if(ftemp(i).ne.(' ')) then 
          c = i
          exit
       endif
    enddo

    do i=1,c
       if(ftemp(i).eq.(' ')) ftemp(i) = '_'
    enddo
    write(unit=temp,fmt='(100a1)')(fdir(i),i=1,9),(ftemp(i),i=1,c)
    read(unit=temp,fmt='(a100)') varname

  end subroutine pruneVarname


!BOP
! 
! !ROUTINE: assemble_obs_cov
! \label{assemble_obs_cov_enksgrace}
!
! !INTERFACE:
  subroutine assemble_obs_cov(nob, N_obs, obs_param, &
       Observations, Obs_cov)

    implicit none
! !ARGUMENTS:    
    integer, intent(in)                              :: nob, N_obs
    type(obs_param_type), dimension(nob), intent(in) :: obs_param
    type(obs_type), dimension(N_obs), intent(in)     :: Observations  
    real, intent(out), dimension(N_obs,N_obs)        :: Obs_cov
! 
! !DESCRIPTION:
!
! assemble measurements error covariance (reichle, 27 Jul 2005)
! 
! \begin{description}
!  \item[nob] number of observation types
!  \item[N\_obs] number of observations
!  \item[obs\_param] observation settings
!  \item[Observations] Observations data type
!  \item[Obs\_cov]     observation error covariance
! \end{description}
!EOP
    
    integer :: i, j, i_species, j_species   !! inum, jnum
    
    real :: fac, xcorr_tmp, ycorr_tmp
    
    ! -------------------------------------------------------------
    
    ! assemble measurement error covariance 
    
    ! initialize

    Obs_cov = 0.

    ! diagonal elements

    do i=1,N_obs

       Obs_cov(i,i) = Observations(i)%std**2
       
    end do
    
    ! off-diagonal elements
    
    do i=1,N_obs
       do j=(i+1),N_obs
          
          i_species = Observations(i)%species
          j_species = Observations(j)%species
          
          ! have non-zero correlation only between observations of same type
          
          if (i_species == j_species) then
             
             xcorr_tmp = obs_param(i_species)%xcorr
             ycorr_tmp = obs_param(i_species)%ycorr
             
             ! check for zero correlation distance 
             
             if (xcorr_tmp>0. .and. ycorr_tmp>0.) then
                
                ! compute correlation between observation locations
                
                !!inum = Observations(i)%catnum
                !!jnum = Observations(j)%catnum 
                
                ! compute Gaussian correlation
                
                !!fac =  & 
                !!  ((tile_coord(inum)%com_lon-tile_coord(jnum)%com_lon)**2 &
                !!  /xcorr_tmp**2 )                                         &
                !!  +                                                       &
                !!  ((tile_coord(inum)%com_lat-tile_coord(jnum)%com_lat)**2 & 
                !!  /ycorr_tmp**2 )               
                
                fac =  & 
                     ((Observations(i)%lon-Observations(j)%lon)**2 &
                     /xcorr_tmp**2 )                                         &
                     +                                                       &
                     ((Observations(i)%lat-Observations(j)%lat)**2 & 
                     /ycorr_tmp**2 )               

                fac = exp(-.5*fac)
                
                Obs_cov(i,j) = Observations(i)%std * Observations(j)%std * fac
                
                Obs_cov(j,i) = Obs_cov(i,j)
                
             end if
          end if
          
       end do
    end do
        
  end subroutine assemble_obs_cov


!BOP
! 
! !ROUTINE: enksgrace_final
! \label{enksgrace_final}
! 
! !INTERFACE: 
  subroutine enksgrace_final
! 
! !DESCRIPTION: 
!  This method performs the finalization for all Enks
!  related structures and subroutines. 
!
!EOP
    deallocate(enksgrace_struc)
  end subroutine enksgrace_final

!BOP
! 
! !ROUTINE: getLSMvarnames
! \label{getLSMvarnames_enksgrace}
! 
! !INTERFACE: 
  subroutine getLSMvarnames(LIS_LSM_State, dim1, varname)
! !USES:

! !ARGUMENTS:         
    type(ESMF_State)      :: LIS_LSM_State
    integer               :: dim1
    character(len=*)      :: varname(dim1)
!
! !DESCRIPTION:
! 
!  This routine retrieves the names of the state prognostic variables
!  from the LSM state. 
!EOP    
    integer               :: status

    call ESMF_StateGet(LIS_LSM_State,itemNameList=varname,rc=status)
    call LIS_verify(status)        

  end subroutine getLSMvarnames

!BOP
! 
! !ROUTINE: getObsPert
! \label{getObsPert_enksgrace}
!
! !INTERFACE: 
  subroutine getObsPert(OBS_Pert_State, gsize,dim1, dim2, pert)
! !USES: 

! !ARGUMENTS: 
    type(ESMF_State)      :: OBS_Pert_State
    integer               :: gsize
    integer               :: dim1, dim2
    real                  :: pert(dim1,dim2)
! 
! !DESCRIPTION: 
! 
!  This routine retrieves the perturbations to be applied to the
!  observations from the Observation perturbation state
! 
! The arguments are: 
! \begin{description}
!  \item[OBS\_Pert\_State]  Observation perturbation state 
!  \item[dim1,dim2]         dimensions of the perturbations array
!  \item[pert]              perturbation array
! \end{description}
!EOP
    real, pointer         :: obs_temp(:,:)
    integer               :: i 
    integer               :: obs_state_count
    integer               :: status
    character*100,allocatable     :: obs_state_objs(:)
    type(ESMF_Field), allocatable :: obs_field(:)

    
    call ESMF_StateGet(OBS_Pert_State,itemCount=obs_state_count,rc=status)
    call LIS_verify(status)
    
    allocate(obs_state_objs(obs_state_count))
    allocate(obs_field(obs_state_count))
    
    call ESMF_StateGet(OBS_Pert_State,itemNameList=obs_state_objs,rc=status)
    call LIS_verify(status)        
    
    do i=1,obs_state_count
       call ESMF_StateGet(OBS_Pert_State,obs_state_objs(i),&
            obs_field(i),rc=status)
       call LIS_verify(status)
       call ESMF_FieldGet(obs_field(i),localDE=0,farrayPtr=obs_temp,rc=status)
       call LIS_verify(status)
!       pert(:,:) = obs_temp(:,:)
       pert((i-1)*gsize + 1:i*gsize,:) = obs_temp(:,:)
    enddo
    deallocate(obs_state_objs)
    deallocate(obs_field)

  end subroutine getObsPert

!BOP
! 
! !ROUTINE: generateObservations
! \label{generateObservations_enksgrace}
! 
! !INTERFACE: 
  subroutine generateObservations(n, k, Nobjs, N_obs_size, LIS_OBS_State, &
       LIS_OBS_Pert_State, Observations )
! !USES: 

! !ARGUMENTS: 
    integer,     intent(IN)        :: n 
    integer,     intent(IN)        :: k
    integer,     intent(IN)        :: Nobjs
    integer,     intent(IN)        :: N_obs_size
    type(ESMF_State)               :: LIS_OBS_State
    type(ESMF_State)               :: LIS_OBS_Pert_State
    type(obs_type)                 :: Observations(N_obs_size)
! 
! !DESCRIPTION: 
!   This subroutine unpacks the Observation state and packages them into the
!   'obs\_type' datastructure. 
!  
!   \begin{description}
!   \item[n]  index of the nest
!   \item[Nobjs] Number of observation types
!   \item[OBS\_State] ESMF State containing observations
!   \item[Observations] Observations datastructure being returned. 
!   \end{description}
!EOP

    integer                        :: typ
    integer                        :: gval, gsize
    integer                        :: count1
    character*100                  :: temp
    integer                        :: i, t,g
    character*1                    :: vid(2)
    integer                        :: status
    type(ESMF_Field)               :: valfield
    type(ESMF_Field)               :: pertField
    real, pointer                  :: value1(:)
    integer, allocatable               :: gid1(:)
    real,    allocatable               :: obsstd(:)
    integer                        :: pert_type
    integer, allocatable               :: obsassimflag(:)
    integer                        :: counts, col, row
    integer                        :: pertyp(N_obs_size)
    real                           :: value(N_obs_size)
    integer                        :: gid(N_obs_size)
    real                           :: std(N_obs_size)
    integer                        :: species(N_obs_size)
    integer                        :: assimflag(N_obs_size)

    count1 = 1
    typ    = 1
    do i=1,Nobjs

       write(unit=temp,fmt='(i2.2)') i
       read(unit=temp,fmt='(2a1)') vid
       
       call ESMF_StateGet(LIS_OBS_State,"Observation"//vid(1)//vid(2),&
            valfield,rc=status)
       call LIS_verify(status)

       call ESMF_StateGet(LIS_OBS_Pert_State,&
            "Observation"//vid(1)//vid(2),&
            pertField,rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField, "Perturbation Type",&
            pert_type, rc=status)
       call LIS_verify(status, &
            "Perturbation Type attribute not found in enksgrace_Mod")


       call ESMF_AttributeGet(LIS_OBS_State,name="Number Of Observations",&
            value = counts,rc=status)
       call LIS_verify(status)

       call ESMF_FieldGet(valfield,localDE=0,farrayPtr=value1,rc=status)
       call LIS_verify(status)
       
       allocate(gid1(counts))
       allocate(obsassimflag(counts))
       allocate(obsstd(counts))

       if(counts.gt.0) then 
          call ESMF_AttributeGet(pertField,"Standard Deviation",obsstd,&
               rc=status)
!               itemCount=counts, rc=status)
          call LIS_verify(status)
          

          call ESMF_AttributeGet(valfield,"Grid Number",gid1,&
               rc=status)
!               itemCount=counts,rc=status)
          call LIS_verify(status)
       
          call ESMF_AttributeGet(valfield,"Assimilation Flag",&
               obsassimflag,rc=status)
          !itemCount=counts,rc=status)
          
          call LIS_verify(status)
       endif

       value(count1: count1+(counts-1))  = value1(:)          
       gid(count1:count1+(counts-1))     = gid1(:)       
       pertyp(count1:count1+(counts-1))     = pert_type
       if(pert_type.eq.0) then 
          std(count1:count1+(counts-1))     = obsstd(:)
       else
          do g = count1,count1+(counts-1)
             if(value1(g).ne.-9999.0) then 
                std(g) = obsstd(g)*value1(g)
                if(std(g).eq.0) std(g) = 0.1
             else
                std(g) = -9999.0
             endif
          enddo
       endif

       species(count1:count1+(counts-1)) = typ
       assimflag(count1:count1+(counts-1)) = obsassimflag(:)
       count1 = count1 + counts
       typ    = typ+1
       
       deallocate(gid1)
       deallocate(obsassimflag)
       deallocate(obsstd)
    enddo

!----------------------------------------------------------------------------
! Map obs_state to Observations object
!----------------------------------------------------------------------------
    gsize = N_obs_size/Nobjs
    do t=1,N_obs_size
       gval = t - nint(real((t-1)/gsize))*gsize
       Observations(t)%species = species(t)
       Observations(t)%catnum = gid(t)

       col = LIS_obs_domain(n,k)%col(t)
       row = LIS_obs_domain(n,k)%row(t)
       Observations(t)%lon = LIS_obs_domain(n,k)%lon(&
            col+(row-1)*LIS_rc%obs_lnc(k))
       Observations(t)%lat = LIS_obs_domain(n,k)%lat(&
            col+(row-1)*LIS_rc%obs_lnc(k))

       Observations(t)%value = value(t)
       Observations(t)%std = std(t)
       Observations(t)%pert_type = pertyp(t)
       if(assimflag(t).eq.1) then 
          Observations(t)%assim = .true.
       else
          Observations(t)%assim = .false.
       endif
    enddo
    
  end subroutine generateObservations

!BOP
! 
! !ROUTINE: generateObsparam
! \label{generateObsparam_enksgrace}
!
! !INTERFACE: 
  subroutine generateObsparam(Nobjs, OBS_Pert_State, obs_param)
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer                       :: Nobjs
    type(ESMF_State)              :: OBS_Pert_State
    type(obs_param_type)          :: obs_param(Nobjs)
!
! !DESCRIPTION: 
!   This routine obtains the metadata information associated with 
!   the observations and puts them into the obs\_param datastructure
!  
!   The arguments are: 
!   \begin{description}
!    \item[Nobjs]            Number of observation types \newline
!    \item[OBS\_Pert\_State] ESMF state containing observation 
!                            perturbations
!    \item[obs\_param]       obs\_param datastructure being 
!                            returned
!   \end{description}
!
!EOP    
    integer                       :: i 
    character*100                 :: temp
    character*1                   :: vid(2)
    type(ESMF_Field)              :: pertField
    integer                       :: status
    real                          :: std_normal_max(Nobjs)
    real                          :: xcorr(Nobjs)
    real                          :: ycorr(Nobjs)

    do i=1,Nobjs
       write(unit=temp,fmt='(i2.2)') i
       read(unit=temp,fmt='(2a1)') vid
       
       call ESMF_StateGet(OBS_Pert_State,"Observation"//vid(1)//vid(2),&
            pertField,rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField,"Std Normal Max",std_normal_max(i),&
            rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField,"X Correlation Scale",xcorr(i),&
            rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField,"Y Correlation Scale",ycorr(i),&
            rc=status)
       call LIS_verify(status)          
       
       obs_param(i)%species = i
       obs_param(i)%assim = .TRUE.
       obs_param(i)%std_normal_max = std_normal_max(i)
       obs_param(i)%xcorr = xcorr(i)
       obs_param(i)%ycorr = ycorr(i)
    enddo
  end subroutine generateObsparam
!BOP
! 
! !ROUTINE: getLSMData
! \label{getLSMData_enksgrace}
!
! !INTERFACE:
  subroutine getLSMData(LIS_LSM_State, dim1, dim2, lsmdata)
! !USES: 

! !ARGUMENTS:         
    type(ESMF_State)      :: LIS_LSM_State
    integer               :: dim1, dim2
    real                  :: lsmdata(dim1, dim2)
! 
! !DESCRIPTION:
! 
! This routine retrives the LSM state prognostic variables from the 
! LSM state. 
!
! The arguments are:
! 
! \begin{description} 
!  \item[LSM\_State]  LSM State containing prognostic variables
!  \item[dim1,dim2]   dimensions of the prognostic variable array
!  \item[lsmdata]     prognostic variable array
! \end{description}
!EOP
    real, pointer         :: lsm_temp(:)
    integer               :: i 
    integer               :: lsm_state_count
    integer               :: status
    character*100,allocatable     :: lsm_state_objs(:)
    type(ESMF_Field), allocatable :: lsm_field(:)

    
    call ESMF_StateGet(LIS_LSM_State,itemCount=lsm_state_count,rc=status)
    call LIS_verify(status)
    
    allocate(lsm_state_objs(lsm_state_count))
    allocate(lsm_field(lsm_state_count))
    
    call ESMF_StateGet(LIS_LSM_State,itemNameList=lsm_state_objs,rc=status)
    call LIS_verify(status)        
    
    do i=1,lsm_state_count
       call ESMF_StateGet(LIS_LSM_State,lsm_state_objs(i),lsm_field(i),&
            rc=status)
       call LIS_verify(status)
       call ESMF_FieldGet(lsm_field(i),localDE=0,farrayPtr=lsm_temp,rc=status)
       call LIS_verify(status)
       lsmdata(i,:) = lsm_temp(:)
    enddo
    deallocate(lsm_state_objs)
    deallocate(lsm_field)

  end subroutine getLSMData
end module enksgrace_Mod
