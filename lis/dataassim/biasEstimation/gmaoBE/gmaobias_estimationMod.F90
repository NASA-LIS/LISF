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
! !MODULE: gmaobias_estimationMod
!
! !DESCRIPTION:
!   This module contains interfaces and subroutines that control
!   the computation of dynamic bias estimation routines using the 
!   approach of Bosilovich et al. (2007)
!
!  The implementation of the algorithm is based on the work of Rolf Reichle at
!  the NASA Global Modeling and Assimilation Office (GMAO) at the NASA GSFC.
!  
! !REVISION HISTORY:
!   26 Nov 2007, Sujay Kumar  : Initial Specification in LIS
! 
!EOP

module gmaobias_estimationMod
  
! ---------------------------------------------------------------------------
!
! Nparam indicates how many bias parameters are estimated per model field:
!
! Nparam = 0 - no bias correction
!
! Nparam = 1 - constant bias corr (w/o diurnal cycle)
!
! Nparam = 3 - diurnal sine/cosine bias corr
! Nparam = 5 - semi-diurnal sine/cosine bias corr
!
! Nparam = 2 - "time-of-day" bias corr with 2 separate bias estimates per day
! Nparam = 4 - "time-of-day" bias corr with 4 separate bias estimates per day
! Nparam = 8 - "time-of-day" bias corr with 8 separate bias estimates per day
!
! The bias estimate is updated from analysis increments whenever observations 
! are available.  The bias time scale relative to the temporal spacing of 
! the observations is described by "tconst_bias".
!
! tconst_bias = dimensionless bias time constant
!
! DEFINITION:  P_bias = tconst_bias * P_forecast   (P is forecast error cov)
!
!   tconst_bias is a.k.a. "gamma" in Dee's 2003 ECWMF proceedings paper
!   
! ASSUMPTION:              tconst << 1
!
! CRUDE APPROXIMATION:     tconst = dt_obs/tcorr
!
!   where dt_obs is the interval between the updates from obs
!   and tcorr is the bias time scale
!
! IN MORE DETAIL:
! 
!   starting from Dee & Todling, MWR, 2000, equation (9) we get
! 
!   b_k = (1-lambda)*b_(k-1) - lambda*innov
!
!   => AR(1) time correlation exp(-dt/tcorr)=(1-lambda)
!
!   use lambda = tconst*var_f/(var_f+var_o) to get
!
!   tconst = (var_f+var_o)/var_f*(1-exp(-dt/tcorr))
!
!   note that (1-exp(-dt/tcorr)) = dt/tcorr  for dt/tcorr<<1
!
! !NOTES: TBD: need to switch the dimensions over to the patch space. A new 
!  implementation of LIS_readvar_reduced_tilespace for the patch space
!  is required -- svk  

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: gmaobiasestimation_init
  PUBLIC :: gmaobiasestimation_setup
  PUBLIC :: gmaobiasestimation_calc
  PUBLIC :: gmaobiasestimation_update
  PUBLIC :: gmaobiasestimation_write_restart
  PUBLIC :: gmaobiasestimation_finalize
  
  integer,   parameter          :: N_bias_param_max = 8

  type, public :: bias_dec_type
     integer             :: Nparam
     real                :: tconst
     real                :: trelax
     real, allocatable       :: param(:)
  end type bias_dec_type

  type, public :: bias_struc_type
     type(bias_dec_type), allocatable :: bias_type(:,:)     
  end type bias_struc_type


  type(bias_struc_type) , allocatable :: bias_struc(:,:)
  

contains

  subroutine gmaobiasestimation_init()

    use LIS_coreMod,      only : LIS_rc

    implicit none

    allocate(bias_struc(LIS_rc%nnest, LIS_rc%ndas))

  end subroutine gmaobiasestimation_init

  subroutine gmaobiasestimation_setup(k)

    use LIS_coreMod,      only : LIS_rc
    use LIS_logMod,       only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun
    use LIS_historyMod,   only : LIS_readvar_reduced_tilespace

    implicit none

    integer                 :: k 
    integer                 :: N_size
    integer                 :: N_state
    integer                 :: N_ens
    integer                 :: v,t,n
    integer                 :: ftn 
    
    integer, allocatable        :: Nparam(:)
    real,    allocatable        :: tconst(:)
    real,    allocatable        :: trelax(:)
    character*40, allocatable   :: vname(:)
    logical                 :: file_exists
    integer                 :: m
    real,    allocatable        :: tempvar(:)

    do n=1,LIS_rc%nnest

       N_size  = LIS_rc%ntiles(n)/LIS_rc%nensem(n)
       N_state = LIS_rc%nstvars(k)
       N_ens   = LIS_rc%nensem(n)
          
       allocate(bias_struc(n,k)%bias_type(N_state, N_size))
       
       allocate(Nparam(N_state))
       allocate(tconst(N_state))
       allocate(trelax(N_state))
       allocate(vname(N_state))
       
       do v=1,N_state
          do t=1,N_size
             allocate(bias_struc(n,k)%bias_type(v,t)%param(N_bias_param_max))
             bias_struc(n,k)%bias_type(v,t)%param(:) = 0.0
          enddo
       enddo
       ftn = LIS_getNextUnitNumber()
       
       write(LIS_logunit,*) 'Opening attributes for bias estimation '
       open(ftn,file=LIS_rc%biasOptionsFile(k), form='formatted', status='old')
       read(ftn,*)
       do v=1,N_state
          read(ftn,fmt='(a40)') vname(v)
          write(LIS_logunit,*) vname(v)
          read(ftn,*) Nparam(v), tconst(v), trelax(v)
          write(LIS_logunit,*) Nparam(v), tconst(v), trelax(v)
          
          bias_struc(n,k)%bias_type(v,:)%Nparam = Nparam(v)
          bias_struc(n,k)%bias_type(v,:)%tconst = tconst(v)
          bias_struc(n,k)%bias_type(v,:)%trelax = trelax(v)
       enddo
       
       call LIS_releaseUnitNumber(ftn)
       
       deallocate(Nparam)
       deallocate(tconst)
       deallocate(trelax)
       deallocate(vname)

       if(trim(LIS_rc%biasrst(n)).ne."none") then ! use the bias restart file
          inquire(file=trim(LIS_rc%biasrstfile(k)), exist=file_exists) 

          if(.not.file_exists) then 
             write(LIS_logunit,*) 'Bias restart file ',trim(LIS_rc%biasrstfile(k)), & 
                  'does not exist '
             write(LIS_logunit,*) 'Program stopping ... '
             call LIS_endrun()          
          endif
          
          ftn=LIS_getNextUnitNumber()
          write(LIS_logunit,*) 'Reading bias restart file ',trim(LIS_rc%biasrstfile(k))
          open(ftn,file=trim(LIS_rc%biasrstfile(k)),form='unformatted')
          
          allocate(tempvar(N_size))
          do v=1,LIS_rc%nstvars(k)
             call LIS_readvar_reduced_tilespace(ftn, n, tempvar)
             do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                bias_struc(n,k)%bias_type(v,t)%Nparam = tempvar(t)
             enddo

             call LIS_readvar_reduced_tilespace(ftn, n, tempvar)
             do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                bias_struc(n,k)%bias_type(v,t)%tconst = tempvar(t)
             enddo
             
             call LIS_readvar_reduced_tilespace(ftn, n, tempvar)
             do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                bias_struc(n,k)%bias_type(v,t)%trelax = tempvar(t)
             enddo
             
             do m=1, N_bias_param_max       
                call LIS_readvar_reduced_tilespace(ftn, n, tempvar)
                do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                   bias_struc(n,k)%bias_type(v,t)%param(m) = tempvar(t)
                enddo
                
             enddo
          enddo
          deallocate(tempvar)
       endif
    enddo
  end subroutine gmaobiasestimation_setup
  
!BOP
! !ROUTINE: gmaobiasestimation_calc
!  \label{gmaobiasestimation_calc}
! 
! !INTERFACE: 
  subroutine gmaobiasestimation_calc(n,k)
! !USES:     
    use ESMF
    use LIS_coreMod, only : LIS_rc
    use LIS_lsmMod,  only : LIS_LSM_State
    use LIS_logMod,  only : LIS_verify, LIS_logunit

    implicit none

! !ARGUMENTS: 
    integer,    intent(in) :: n 
    integer,    intent(in) :: k 
! 
! !DESCRIPTION: 
! 
! 
!EOP
    integer                :: N_size, N_ens, N_tile, N_state
    integer                :: v
    real,          pointer :: stdata(:)
    real                   :: stvars(LIS_rc%nstvars(k),LIS_rc%ntiles(n))
!    real,          allocatable :: stvars_vec(:,:)
    character*100, allocatable :: lsm_state_objs(:)
    type(ESMF_Field)       :: lsm_field(LIS_rc%nstvars(k))
    integer                :: status
    
    allocate(lsm_state_objs(LIS_rc%nstvars(k)))

    call ESMF_StateGet(LIS_LSM_State(n,k),itemNameList=lsm_state_objs,&
         rc=status)
    call LIS_verify(status, &
         'ESMF_StateGet failed in gmaobiasestimation_calc')

    do v=1,LIS_rc%nstvars(k)
       call ESMF_StateGet(LIS_LSM_State(n,k),trim(lsm_state_objs(v)), &
            lsm_field(v),rc=status)
       call LIS_verify(status,&
            'ESMF_StateGet failed in gmaobiasestimation_calc')
       
       call ESMF_FieldGet(lsm_field(v),localDE=0,farrayPtr=stdata, rc=status)
       call LIS_verify(status,&
            'ESMF_FieldGet failed in gmaobiasestimation_calc')
       
       stvars(v,:)     = stdata(:)
    enddo

    N_tile = LIS_rc%ntiles(n)
    N_ens  = LIS_rc%nensem(n)
    N_size = N_tile/N_ens
    N_state = LIS_rc%nstvars(k)
    
!    allocate(stvars_vec(N_state, N_size))
    ! --------------------------------------------------------------
    
    ! apply bias correction to cat_progn, relax bias parameters
    
    call bias_corr( n, k, N_tile, N_state, N_size, N_ens, Stvars)
    
    ! ------------------------------------------------------------------
    !
    ! recompute diagnostics 
    !
    ! NOTE: subroutine recompute_diagnostics() is based on update_type,
    !       not on cat_bias(:)%Nparam.  This means that the relevant 
    !       diagnostics are recomputed even if no bias correction was
    !       applied.  To avoid this, for example for soil moisture 
    !       assimilation, it is best to compile the code without the bias
    !       module whenever bias is not estimated.
    
!    call get_ens_avg(N_tile, N_size, N_state, N_ens, &
!         Stvars, Stvars_vec)
!    call get_cat_progn_ens_avg( N_catd, N_ens, cat_progn, cat_progn_vec )
    
!    call qc and set from recompute? 
!    call recompute_diagnostic( N_catd, update_type, &
!         cat_param, cat_progn_vec, cat_diagn)
    
    ! ------------------------------------------------------------------

!----------------------------------------------------------------------------
! Updating LSM_State
!----------------------------------------------------------------------------
    do v=1,LIS_rc%nstvars(k)
       call ESMF_FieldGet(lsm_field(v),localDE=0,farrayPtr=stdata,rc=status)
       call LIS_verify(status,&
            'ESMF_FieldGet failed in gmaobiasestimation_calc')
       stdata(:) = stvars(v,:)
    enddo

    call lsmdaqcstate(LIS_rc%lsm,LIS_rc%daset(k),n,LIS_LSM_State(n,k))
    call lsmdasetstatevar(LIS_rc%lsm,LIS_rc%daset(k),n, LIS_LSM_State(n,k))

!    deallocate(stvars_vec)
    deallocate(lsm_state_objs)

  end subroutine gmaobiasestimation_calc


!BOP
! 
! !ROUTINE: gmaobiasestimation_update
! \label{gmaobiasestimation_update}
! 
! !INTERFACE: 
  subroutine gmaobiasestimation_update(n,k)
! !USES: 
    use ESMF
    use LIS_coreMod,     only : LIS_rc
    use LIS_lsmMod,      only : LIS_LSM_Incr_State
    use LIS_logMod,      only : LIS_verify
    
    implicit none
! !ARGUMENTS: 
    integer,        intent(in)  :: n 
    integer,        intent(in)  :: k
! 
! !DESCRIPTION: 
! 
!EOP
    integer                     :: N_size, N_ens, N_tile, N_state
    integer                     :: v
    character*100,    allocatable   :: lsm_state_objs(:)
    real,             pointer   :: stincrdata(:)
    real,           allocatable :: State_incr_vec(:,:)
    real                        :: state_incr(LIS_rc%nstvars(k), LIS_rc%ntiles(n))
    logical                     :: fresh_incr
    type(ESMF_Field)            :: lsm_incr_field(LIS_rc%nstvars(k))
    integer                     :: status

    call ESMF_AttributeGet(LIS_LSM_Incr_State(n,k), &
         "Fresh Increments Status", value = fresh_incr, rc=status)
    call LIS_verify(status,&
         'ESMF_AttributeGet: Fresh Increments Status failed in gmaobiasestimation_update')

    allocate(lsm_state_objs(LIS_rc%nstvars(k)))

    call ESMF_StateGet(LIS_LSM_Incr_State(n,k), itemNameList=lsm_state_objs,&
         rc=status)
    call LIS_verify(status,&
         'ESMF_StateGet failed in gmaobiasestimation_update')

    do v=1,LIS_rc%nstvars(k)
       call ESMF_StateGet(LIS_LSM_Incr_State(n,k),trim(lsm_state_objs(v)), &
            lsm_incr_field(v), rc=status)
       call LIS_verify(status,&
            'ESMF_StateGet failed in gmaobiasestimation_update')

       call ESMF_FieldGet(lsm_incr_field(v),localDE=0,farrayPtr=stincrdata, rc=status)
       call LIS_verify(status,&
            'ESMF_FieldGet failed in gmaobiasestimation_update')

       state_incr(v,:) = stincrdata(:)
       
    enddo
    
    N_tile = LIS_rc%ntiles(n)
    N_ens  = LIS_rc%nensem(n)
    N_size = N_tile/N_ens
    N_state = LIS_rc%nstvars(k)

    allocate(State_incr_vec(N_state, N_size))

    if(fresh_incr) then 
       
       ! update bias parameters from ensemble average increments
       
       call get_ens_avg(N_tile, N_size, N_state, N_ens, &
            State_incr, State_incr_vec)

!       call get_cat_progn_ens_avg(N_catd, N_ens, cat_progn_incr, &
!            cat_progn_incr_vec)
       
       call bias_update(n, k, N_state, N_size, State_incr_vec)
!       call bias_update(date_time, model_dtstep_real, N_catd, &
!            cat_progn_incr_vec, cat_bias )
       
    end if
!----------------------------------------------------------------------------
! Updating LSM_Incr_State
!----------------------------------------------------------------------------
    do v=1,LIS_rc%nstvars(k)

       call ESMF_FieldGet(lsm_incr_field(v),localDE=0,farrayPtr=stincrdata,rc=status)
       call LIS_verify(status,&
            'ESMF_FieldGet failed in gmaobiasestimation_update')

       stincrdata(:) = state_incr(v,:)       
    enddo


    deallocate(lsm_state_objs)
    deallocate(State_incr_vec)

  end subroutine gmaobiasestimation_update

!BOP
! !ROUTINE: gmaobiasestimation_restart
!  \label{gmaobiasestimation_restart}
! 
! !INTERFACE: 
  subroutine gmaobiasestimation_write_restart(n,k)
! !USES:    
    use LIS_historyMod, only    : LIS_writevar_reduced_tilespace
    use LIS_coreMod,    only    : LIS_rc, LIS_masterproc
    use LIS_logMod,     only    : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber
    use LIS_fileIOMod,  only    : LIS_create_output_directory
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
! 
! !ARGUMENTS: 
!     
    integer,   intent(in)       :: n
    integer,   intent(in)       :: k
! 
! !DESCRIPTION: 
!   This subroutine generates the bias restart file based on the
!   specified output frequency. 
!   
!EOP
    integer       :: ftn
    integer       :: v, t, m
    real          :: tempvar(LIS_rc%ntiles(n)/LIS_rc%nensem(n))
    real          :: curr_time
    character(len=LIS_CONST_PATH_LEN) :: rstname

    curr_time = float(LIS_rc%hr)*3600+60*float(LIS_rc%mn)+float(LIS_rc%ss)

    if(mod(curr_time,real(LIS_rc%biasrstInterval(n))).eq.0) then

       if(LIS_masterproc) then 
          ftn = LIS_getNextUnitNumber()
          call gmaobias_restart_filename(rstname)
          
          call LIS_create_output_directory('GMAOBE')
          write(LIS_logunit,*) 'Writing bias restart ',trim(rstname)
          open(ftn,file=trim(rstname),form='unformatted')
       endif

       do v=1,LIS_rc%nstvars(k)
          
          do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
             tempvar(t) = bias_struc(n,k)%bias_type(v,t)%Nparam
!             print*, t, tempvar(t)
          enddo
          call LIS_writevar_reduced_tilespace(ftn,n,tempvar)

          do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
             tempvar(t) = bias_struc(n,k)%bias_type(v,t)%tconst
          enddo
          call LIS_writevar_reduced_tilespace(ftn,n,tempvar)
          do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
             tempvar(t) = bias_struc(n,k)%bias_type(v,t)%trelax
          enddo
          call LIS_writevar_reduced_tilespace(ftn,n,tempvar)
          
          do m=1, N_bias_param_max       
             do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                tempvar(t) = bias_struc(n,k)%bias_type(v,t)%param(m)
             enddo
             call LIS_writevar_reduced_tilespace(ftn,n,tempvar)
          enddo
       enddo
       
       if(LIS_masterproc) then 
          call LIS_releaseUnitNumber(ftn)
       endif
    endif

  end subroutine gmaobiasestimation_write_restart

!BOP
! !ROUTINE: gmaobias_restart_filename
! \label{gmaobias_restart_filename}
!
! !INTERFACE: 
  subroutine gmaobias_restart_filename(filename)
! !USES: 
    use LIS_coreMod,   only : LIS_rc
! !ARGUMENTS:     
    character(len=*)       :: filename
!
! !DESCRIPTION: 
! This subroutine generates a time stamped name of the bias restart file
!
!EOP
    character(len=12)      :: cdate1
    
    write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
         LIS_rc%yr, LIS_rc%mo, &
         LIS_rc%da, LIS_rc%hr,LIS_rc%mn

    filename = trim(LIS_rc%odir)//&
         '/GMAOBE/Bias_rst_'//cdate1//'.1gs4r'
    
  end subroutine gmaobias_restart_filename
  
!BOP
! 
! !ROUTINE: bias_update
! \label{bias_update}
! 
! !REVISION HISTORY: 
! 19 Oct 2005, Rolf Reichle : Initial Specification
! 26 Nov 2007, Sujay Kumar  : Included in LIS
! 19 Aug 2008, Rolf Reichle : Added call to get_bias_const_param()
!
! !INTERFACE: 
  subroutine bias_update(n, k, N_state, N_size, progn_incr)
! !USES: 
    use LIS_coreMod,         only      : LIS_rc
    
    implicit none

! !ARGUMENTS: 
    integer,           intent(in)    :: n 
    integer,           intent(in)    :: k 
    integer,           intent(in)    :: N_state
    integer,           intent(in)    :: N_size
    real,              intent(in)    :: progn_incr(N_state, N_size)
! 
! !DESCRIPTION: 
! 
!EOP
    integer                          :: v,t 
    real                             :: const_param(N_bias_param_max)
    
    call bias_get_const_param(const_param)

    do v=1,N_state
       do t=1,N_size
          if( bias_struc(n,k)%bias_type(v,t)%Nparam.gt.0 ) then 
             call bias_update_helper(real(LIS_rc%ts), &
                  bias_struc(n,k)%bias_type(v,t)%Nparam,&
                  const_param, &
                  bias_struc(n,k)%bias_type(v,t)%tconst, &
                  progn_incr(v,t), &
                  bias_struc(n,k)%bias_type(v,t)%param(1:bias_struc(n,k)%bias_type(v,t)%Nparam))
          endif
       enddo
    enddo
  end subroutine bias_update

!BOP
! 
! !ROUTINE: bias_update_helper
! \label{bias_update_helper}
! 
! !REVISION HISTORY: 
!  17 Oct 2005 Rolf Reichle:  Initial Specification
!  26 Nov 2006 Sujay Kumar :  Included in LIS
!
! !INTERFACE: 
  subroutine bias_update_helper( dtstep, Nparam, const_param, tconst, &
       increment, bias_param )
    
    implicit none
! !ARGUMENTS:     
    real, intent(in) :: dtstep, tconst, increment
    integer, intent(in) :: Nparam    
    real, dimension(Nparam), intent(in) :: const_param
    real, dimension(Nparam), intent(inout) :: bias_param
!
! !DESCRIPTION: 
! update bias parameters for a single tile and a single field of cat\_progn
!
!EOP    
    
    ! local variables

    real, dimension(Nparam) :: const_param_tmp

    real :: tmp_real
    
    integer :: i, ind_start, ind_end
    
    ! -------------------------------------------------------------
    
    if (Nparam > 0) then
       
       call bias_options_helper( Nparam, const_param, &
            ind_start, ind_end, const_param_tmp )
       
       ! update bias parameters
       
       tmp_real = - tconst * increment / dtstep
       
       do i=ind_start,ind_end         
          bias_param(i) = bias_param(i) + tmp_real*const_param_tmp(i)
       end do
    end if
    
  end subroutine bias_update_helper

!BOP
! 
! !ROUTINE: get_ens_avg
! \label{get_ens_avg}
! 
! !INTERFACE: 
subroutine get_ens_avg(N_tile, N_size, N_state, N_ens, &
     State_incr, State_incr_vec)
! 
! !DESCRIPTION: 
! 
!EOP
  implicit none

  integer,        intent(in)    :: N_tile
  integer,        intent(in)    :: N_size
  integer,        intent(in)    :: N_state
  integer,        intent(in)    :: N_ens
  real,           intent(in)    :: State_incr(N_state, N_tile)
  real,           intent(inout) :: State_incr_vec(N_state, N_size)

  integer                       :: i, m, tid
  
  do i=1,N_size
     State_incr_vec(:,i) = 0 
     do m=1,N_ens           
        tid = (i-1)*N_ens+m
        State_incr_vec(:,i) = State_incr_vec(:,i) + State_incr(:,tid)
     enddo
     State_incr_vec(:,i) = State_incr_vec(:,i)/real(N_ens)
  enddo
  
end subroutine get_ens_avg
  
!BOP
! 
! !ROUTINE: bias_corr
! \label{bias_corr}
! 
! !REVISION HISTORY: 
! 19 Oct 2005, Rolf Reichle : Initial Specification
! 26 Nov 2007, Sujay Kumar  : Included in LIS
! 19 Aug 2008, Rolf Reichle : Added call to get_bias_const_param()
! 
! !INTERFACE: 
  subroutine bias_corr(n, k, N_tile, N_state, N_size, N_ens, progn_field)
! !USES: 
    use LIS_coreMod, only : LIS_rc

! 
! !DESCRIPTION: 
! apply bias correction to prognostic variables, relaxes bias 
! parameters
! 
!EOP
    implicit none
    integer,        intent(in)      :: n 
    integer,        intent(in)      :: k
    integer,        intent(in)      :: N_tile
    integer,        intent(in)      :: N_state
    integer,        intent(in)      :: N_size
    integer,        intent(in)      :: N_ens
    real,           intent(INOUT)   :: progn_field(N_state, N_tile)


    integer                          :: v,t 
    real                             :: const_param(N_bias_param_max)
!    real                             :: bias_param(N_bias_param_max)

    call bias_get_const_param(const_param)
    
    do v=1,N_state
       do t=1,N_size
          call bias_corr_helper(real(LIS_rc%ts),&
               N_ens, &
               bias_struc(n,k)%bias_type(v,t)%Nparam,&
               const_param, &
               bias_struc(n,k)%bias_type(v,t)%trelax, &
               progn_field(v,(t-1)*N_ens+1:(t-1)*N_ens+N_ens),&
               bias_struc(n,k)%bias_type(v,t)%param(1:bias_struc(n,k)%bias_type(v,t)%Nparam))
       enddo
    enddo

  end subroutine bias_corr

!BOP
! 
! !ROUTINE: bias_corr_helper
! \label{bias_corr_helper}
!
! !REVISION HISTORY: 
!  17 Oct 2005 Rolf Reichle:  Initial Specification
!  26 Nov 2006 Sujay Kumar :  Included in LIS
! 
! !INTERFACE: 
  subroutine bias_corr_helper( dtstep, N_ens, Nparam,   &
       const_param, trelax, progn_field, bias_param )
    
    implicit none

!
! !DESCRIPTION:     
! diagnose/apply bias flux and relax bias parameters for a single 
! tile and a single field of progn
!
!EOP

    real,    intent(in) :: dtstep, trelax

    integer, intent(in) :: N_ens, Nparam
    
    real, dimension(Nparam), intent(in) :: const_param
    
    real, dimension(N_ens), intent(inout) :: progn_field
    
    real, dimension(Nparam), intent(inout) :: bias_param
    
    ! local variables
    
    real :: tmpflux, relax_fac
    
    integer :: i, ind_start, ind_end
    
    real, dimension(Nparam) :: const_param_tmp
    
    ! -----------------------------------------------------------------    

    if (Nparam > 0) then

       call bias_options_helper( Nparam, const_param, &
            ind_start, ind_end, const_param_tmp )

       ! diagnose bias flux from bias parameters
       
       tmpflux = 0.
       
       do i=ind_start,ind_end
          
          tmpflux = tmpflux + bias_param(i)*const_param_tmp(i)
          
       end do
       
       ! apply bias flux
       
       progn_field(:) = progn_field(:) - tmpflux*dtstep
       
       ! relax bias parameters
       
       relax_fac = (1. - dtstep/trelax)
       
       do i=ind_start,ind_end
          
          bias_param(i) = bias_param(i) * relax_fac
          
       end do
       
    end if
    
  end subroutine bias_corr_helper

!BOP
! 
! !ROUTINE: bias_get_const_param
! \label{bias_get_const_param}
! 
! !REVISION HISTORY: 
!  18 Aug 2008 Rolf Reichle:  Initial Specification
!  19 Aug 2008 Rolf Reichle:  Included in LIS
!
! !INTERFACE: 
  subroutine bias_get_const_param(const_param)
! !USES: 
    use LIS_coreMod,         only      : LIS_rc
    use LIS_constantsMod,    only      : LIS_CONST_PI
    
    implicit none

! !ARGUMENTS: 
    real,           intent(out)   :: const_param(N_bias_param_max) 
! 
! !DESCRIPTION: 
! 
!EOP
    real,       parameter            :: omega1 = 2*LIS_CONST_PI/86400.0
    real,       parameter            :: omega2 = 2*LIS_CONST_PI/43200.0
    real                             :: secs_in_day
    real                             :: om1_t
    real                             :: om2_t

    ! -----------------------------------------------------------------    

    secs_in_day = real(LIS_rc%hr)*3600.0+real(LIS_rc%mn)*60.0+&
         real(LIS_rc%ss)
    
    om1_t    = secs_in_day*omega1
    om2_t    = secs_in_day*omega2
    
    const_param(1) = 1.0
    
    const_param(2) = cos(om1_t)
    const_param(3) = sin(om1_t)
    
    const_param(4) = cos(om2_t)
    const_param(5) = sin(om2_t)

  end subroutine bias_get_const_param



!BOP
! 
! !ROUTINE: bias_time_of_day_index
! \label{bias_time_of_day_index}
! 
! !REVISION HISTORY: 
!  18 Aug 2008 Rolf Reichle:  Initial Specification
!  19 Aug 2008 Rolf Reichle:  Included in LIS
!
! !INTERFACE: 
  integer function bias_time_of_day_index( Nparam )
! !USES: 
    use LIS_coreMod,         only      : LIS_rc
    
    implicit none

! !ARGUMENTS: 
    integer, intent(in)  :: Nparam             
! 
! !DESCRIPTION: 
! 
!EOP
    
    integer :: tmpint, dtstep_tmp
    
    ! -----------------------------------------------
    
    dtstep_tmp = 86400/Nparam
    
    ! seconds-in-day 
    
    tmpint = LIS_rc%hr*3600 + LIS_rc%mn*60 + LIS_rc%ss 
    
    ! add half of dtstep_assim (to center interval around nominal time)
    
    tmpint = tmpint + dtstep_tmp/2
    
    ! modulus calculation (yields tmpint ranging from 1 to N+1)
    
    tmpint = tmpint/dtstep_tmp + 1
    
    ! fix "N+1" with another modulus calculation
    
    bias_time_of_day_index = mod( tmpint-1, 86400/dtstep_tmp) + 1
    
  end function bias_time_of_day_index
  
  ! --------------------------------------------------------


!BOP
! 
! !ROUTINE: bias_options_helper
! \label{bias_options_helper}
! 
! !REVISION HISTORY: 
!  18 Aug 2008 Rolf Reichle:  Initial Specification
!  19 Aug 2008 Rolf Reichle:  Included in LIS
!
! !INTERFACE: 
  subroutine bias_options_helper( Nparam, const_param_in, &
       ind_start, ind_end, const_param_out )
! !USES: 
    
    implicit none

! !ARGUMENTS: 
    integer,                    intent(in)  :: Nparam     
    real,    dimension(Nparam), intent(in)  :: const_param_in
    
    integer,                    intent(out) :: ind_start, ind_end    
    real,    dimension(Nparam), intent(out) :: const_param_out
! 
! !DESCRIPTION: 
! 
!EOP
    
    select case (Nparam)
       
    case (1,3,5)          ! constant or sine/cosine (semi-)diurnal bias corr
       
       ind_start = 1
       ind_end   = Nparam
       
       const_param_out = const_param_in
       
    case (2,4,8)          ! separate "time-of-day" bias correction
       
       ind_start = bias_time_of_day_index( Nparam )
       ind_end   = ind_start
       
       const_param_out = 1.
       
    case default
       
       write (*,*) 'bias_options_helper(): unknown Nparam = ', Nparam
       
    end select
    
  end subroutine bias_options_helper
  
  ! --------------------------------------------------------

  subroutine gmaobiasestimation_finalize
    
    deallocate(bias_struc)
  end subroutine gmaobiasestimation_finalize

end module gmaobias_estimationMod
