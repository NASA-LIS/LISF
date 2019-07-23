!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: LVT_coreMod
! \label(LVT_coreMod)
!
! !INTERFACE:
module LVT_coreMod
! 
! !USES:       
  use ESMF
  use LVT_timeMgrMod
  use LVT_PRIV_rcMod
  use LVT_PRIV_gridMod
  use LVT_PRIV_tileMod
  use LVT_logMod
  use LVT_mpiMod
  use map_utils

  PRIVATE
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  The code in this file contains the basic datastructures and 
!  control routines for the operation of LIS verification toolkit
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_configinit
  public :: LVT_ticktime
  public :: LVT_endofrun
  public :: LVT_isAtAfinerResolution
  public :: LVT_core_init
  public :: LVT_557post_alarm_is_on ! EMK
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LVT_rc
  public :: LVT_LIS_rc
  public :: LVT_vm
  public :: LVT_config
  public :: LVT_localPet
  public :: LVT_npes
  public :: LVT_masterproc
  public :: LVT_domain 
  public :: LVT_LIS_domain 
  public :: LVT_surface

  public :: LVT_ews_ind   ! East-West starting index of the LVT domain
  public :: LVT_ewe_ind   ! East-West ending index of the LVT domain
  public :: LVT_nss_ind   ! North-South starting index of the LVT domain
  public :: LVT_nse_ind   ! North-South ending index of the LVT domain    

  public :: LVT_lis_ews_ind   ! East-West starting index of the LIS domain
  public :: LVT_lis_ewe_ind   ! East-West ending index of the LIS domain
  public :: LVT_lis_nss_ind   ! North-South starting index of the LIS domain
  public :: LVT_lis_nse_ind   ! North-South ending index of the LIS domain     

  public :: LVT_ews_halo_ind ! East-West starting index of the LVT domain with halo   
  public :: LVT_ewe_halo_ind ! East-West ending index of the LVT domain with halo   
  public :: LVT_nss_halo_ind ! North-South starting index of the LVT domain with halo    
  public :: LVT_nse_halo_ind ! North-South ending index of the LVT domain with halo       

  public :: LVT_lis_ews_halo_ind ! East-West starting index of the LIS domain with halo     
  public :: LVT_lis_ewe_halo_ind ! East-West ending index of the LIS domain with halo     
  public :: LVT_lis_nss_halo_ind ! North-South starting index of the LIS domain with halo        
  public :: LVT_lis_nse_halo_ind ! North-South ending index of the LIS domain with halo          

  public :: LVT_deltas     
  public :: LVT_offsets   
  public :: LVT_lis_deltas     
  public :: LVT_lis_offsets   
  public :: LVT_tdeltas    
  public :: LVT_toffsets   
  public :: LVT_lis_tdeltas    
  public :: LVT_lis_toffsets   
  public :: LVT_gdeltas    
  public :: LVT_goffsets   
  public :: LVT_lis_gdeltas    
  public :: LVT_lis_goffsets   
  public :: LVT_odeltas    
  public :: LVT_ooffsets   
  public :: LVT_lis_patch_deltas
  public :: LVT_lis_patch_offsets
  public :: LVT_ntiless       
  public :: LVT_ngrids     
  public :: LVT_lis_ntiless       
  public :: LVT_lis_ngrids     
  public :: LVT_lis_npatches
  public :: LVT_vecTile   
  public :: LVT_vecPatch
  public :: LVT_vecGrid   
  public :: LVT_ensOnGrid 
!EOP

  type, public :: lvt_domain_type 
     type(proj_info)            :: lvtproj
     type(proj_info)            :: lvtparamproj
     integer, allocatable       :: gindex(:,:)
     type(griddec), allocatable :: grid(:)
     real                       :: minLat, maxLat, minLon, maxLon
  end type lvt_domain_type

  type, public :: lis_domain_type
     type(proj_info)            :: proj
     integer, allocatable       :: gindex(:,:)
     type(griddec), allocatable :: grid(:)
     type(tiledec), allocatable :: tile(:)
     integer, allocatable       :: ntiles_pergrid(:)
     integer, allocatable       :: str_tind(:) !starting tile id for each grid
  end type lis_domain_type

  type, public :: lvt_domain_sf_type 
     type(tiledec), allocatable :: lis_tile(:)
     integer,       allocatable :: npatch_pergrid(:)
     integer,       allocatable :: str_patch_ind(:)
  end type lvt_domain_sf_type

  type(lvtrcdec),                  save :: LVT_rc
  type(lisrcdec),                  save :: LVT_LIS_rc(3)
  type(lvt_domain_type),           save :: LVT_domain
  type(lis_domain_type),           save :: LVT_LIS_domain(3)
  type(ESMF_VM),                   save :: LVT_vm
  type(ESMF_Config),               save :: LVT_config
  integer                               :: LVT_localPet
  integer                               :: LVT_npes
  logical                               :: LVT_masterproc
  type(lvt_domain_sf_type), allocatable :: LVT_surface(:,:)

  integer, allocatable                  :: LVT_ews_ind(:),LVT_ewe_ind(:)
  integer, allocatable                  :: LVT_nss_ind(:),LVT_nse_ind(:)
  integer, allocatable                  :: LVT_ews_halo_ind(:),LVT_ewe_halo_ind(:)
  integer, allocatable                  :: LVT_nss_halo_ind(:),LVT_nse_halo_ind(:)
  integer, allocatable                  :: LVT_deltas(:),LVT_offsets(:)
  integer, allocatable                  :: LVT_tdeltas(:),LVT_toffsets(:)
  integer, allocatable                  :: LVT_gdeltas(:),LVT_goffsets(:)
  integer, allocatable                  :: LVT_odeltas(:),LVT_ooffsets(:)
  integer, allocatable                  :: LVT_ntiless(:),LVT_ngrids(:)
 
  integer, allocatable                  :: LVT_lis_deltas(:,:),LVT_lis_offsets(:,:)
  integer, allocatable                  :: LVT_lis_tdeltas(:,:),LVT_lis_toffsets(:,:)
  integer, allocatable                  :: LVT_lis_gdeltas(:,:),LVT_lis_goffsets(:,:)
  integer, allocatable                  :: LVT_lis_ews_ind(:,:),LVT_lis_ewe_ind(:,:)
  integer, allocatable                  :: LVT_lis_nss_ind(:,:),LVT_lis_nse_ind(:,:)
  integer, allocatable                  :: LVT_lis_ews_halo_ind(:,:),LVT_lis_ewe_halo_ind(:,:)
  integer, allocatable                  :: LVT_lis_nss_halo_ind(:,:),LVT_lis_nse_halo_ind(:,:)
  integer, allocatable                  :: LVT_lis_ntiless(:,:),LVT_lis_ngrids(:,:)
  integer, allocatable                  :: LVT_lis_npatches(:,:,:)
  integer, allocatable                  :: LVT_lis_patch_deltas(:,:,:), LVT_lis_patch_offsets(:,:,:)
  type(ESMF_Grid), save                 :: LVT_vecTile(3)
  type(ESMF_Grid), allocatable          :: LVT_vecPatch(:,:)
  type(ESMF_Grid), save                 :: LVT_vecGrid(3)
  type(ESMF_Grid), save                 :: LVT_ensOnGrid(3)
  
contains
!BOP
! 
! !ROUTINE: LVT_configinit
! \label{LVT_configinit}
!
! !INTERFACE: 
  subroutine LVT_configinit(configfile)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  Performs the initialization of LVT runtime configuration for an 
!  offline (uncoupled simulation).   
! 
!  The arguments are: 
!  \begin{description}
!    \item[configfile]
!       name of the lvt configuration file
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!    \item[spmd\_init\_offline] (\ref{spmd_init_offline}) \newline
!      performs MPI initialization
!    \item[LVT\_log\_init] (\ref{LVT_log_init}) \newline
!      initialize the LVT logging modules
!    \item[LVT\_readConfig](\ref{LVT_readConfig}) \newline
!      reads the LVT configuration file
!    \item[spmd\_setup] (\ref{spmd_setup}) \newline
!      allocates memory for the variables to track domain decomposition
!    \item[LVT\_timemgr\_init](\ref{LVT_timemgr_init}) \newline
!      initializes the LVT time manager
!   \end{description}
! 
!EOP
    implicit none
    character(len=*) :: configfile
    integer :: k
    integer :: status

    call spmd_init_offline(LVT_vm)
    call LVT_log_init(LVT_getNextUnitNumber())
    call LVT_readConfig(configfile)
    call spmd_setup()
    call LVT_timemgr_init(LVT_rc)

    allocate(LVT_surface(2,LVT_rc%max_model_types))
    do k=1,2
       allocate(LVT_LIS_rc(k)%glbnpatch(LVT_rc%max_model_types))
       allocate(LVT_LIS_rc(k)%glbnpatch_red(LVT_rc%max_model_types))
    enddo

    allocate(LVT_lis_npatches(2,LVT_rc%max_model_types,0:LVT_npes-1))
    allocate(LVT_lis_patch_offsets(2,LVT_rc%max_model_types,0:LVT_npes-1))
    allocate(LVT_lis_patch_deltas(2,LVT_rc%max_model_types,0:LVT_npes-1))

    allocate(LVT_vecPatch(2,LVT_rc%max_model_types))

    LVT_rc%endtime = 0
    LVT_rc%gridDesc = 0 

  end subroutine LVT_configinit

!BOP
! 
! !ROUTINE: LVT_ticktime
! \label{LVT_ticktime}
!
! !INTERFACE: 
  subroutine LVT_ticktime()
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine calls the time manager to increment the runtime 
!  clock by sthe model timestep. LVT keeps track of the current
!  time (LVT_rc%yr, LVT_rc%mo, ...), the time at the next timestep
!  (LVT_rc%nyr, LVT_rc%nmo, ...), time based on the specified
!  temporal lag (LVT_rc%lyr, LVT_rc%lmo, ...) and the time at
!  the next timestep for in the lagged space (LVT_rc%l_nyr, 
!  LVT_rc%l_nmo, ...). Finally vairiables with the 'd' prefix
!  (LVT_rc%dyr, LVT_rc%dmo, are used to keep track of the 
!  time based on the source of the data. Datastream source 1
!  is given the current time and datastream 2 is given the 
!  lagged time. 
!  
!
!  The routine also determines if the temporal averaging interval
!  is met. This is computed differently for time averaged and 
!  instantanteous variables. For instantaneous variables, the 
!  'next' time is used to determine if the temporal averaging 
!  interval is met. For the time averaged variables, the 'current'
!  time values are used since the temporal averaging is assumed to be
!  back-averaged. 
! 
!  A special timestepping for 'dekadal' intervals are also included. 
!  It is, however, possible to use the regular timestepping even 
!  if the LIS outputs are only available at dekadal intervals. LVT
!  would simply skip over temporal locations where the LIS output
!  doesn't exist. 
!
!  The Calling sequence is : 
!  \begin{description}
!   \item[LVT\_advance\_timestep] (\ref{LVT_advance_timestep}) \newline
!    advances the clock 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer                 :: yr,mo,da,hr,mn,ss
    integer                 :: status
    type(ESMF_Time)         :: currTime, currTime1, lagTime, advTime
    type(ESMF_Time)         :: start_dekad1, start_dekad2, start_dekad3
    type(ESMF_TimeInterval) :: ts, lagTS
    integer*4, parameter    :: gDEKAD_TIME_SCHEMA = 10
    integer*4               :: dekOy
    integer                 :: days(12)
    integer                 :: tfact
    data days /31,28,31,30,31,30,31,31,30,31,30,31/

    call LVT_advance_timestep(LVT_rc)

    LVT_rc%computeFlag = .false. 
    if(LVT_rc%tsconv.eq."regular") then 
       LVT_rc%chkTS = .true. 
       
! The current time is advanced by one timestep (LIS output frequency)
! to see if it is time to compute the metrics and write time series 
! outputs. 

       call ESMF_TimeSet(currTime,  yy=LVT_rc%yr, &
            mm = LVT_rc%mo, &
            dd = LVT_rc%da, &
            h = LVT_rc%hr, &
            m = LVT_rc%mn, &
            s = LVT_rc%ss,&
            calendar = LVT_calendar, &
            rc=status)
       call LVT_verify(status, 'error in LIS_ticktime')
       
       currTime1 = currTime
       call ESMF_TimeIntervalSet(ts,s=LVT_rc%ts,rc=status)
       call LVT_verify(status, 'error in LIS_ticktime')
       
       currTime = currTime + ts
       call ESMF_TimeGet(currTime,  yy=LVT_rc%nyr, &
            mm = LVT_rc%nmo, &
            dd = LVT_rc%nda, &
            h = LVT_rc%nhr, &
            m = LVT_rc%nmn, &
            s = LVT_rc%nss,&
            calendar = LVT_calendar, &
            rc=status)
       call LVT_verify(status, 'error in LIS_ticktime')       
       
       
       if(mod(LVT_rc%tlag,2592000).ne.0) then 
          call ESMF_TimeIntervalSet(lagTS,s=LVT_rc%tlag,rc=status)
          call LVT_verify(status, 'error in LIS_ticktime')
          lagTime = currTime1 + lagTS
       else
          tfact = LVT_rc%tlag/2592000
          LVT_rc%lmo = LVT_rc%mo + tfact
          LVT_rc%lyr = LVT_rc%yr
!-------------------------------------------------------------------------
! Assuming here that the month interval will not span more than 12 months
!-------------------------------------------------------------------------
          if(LVT_rc%lmo.gt.12) then 
             LVT_rc%lyr = LVT_rc%lyr + 1
             LVT_rc%lmo = LVT_rc%lmo -12
          endif
          if(LVT_rc%lmo.lt.1) then 
             LVT_rc%lyr = LVT_rc%lyr - 1
             LVT_rc%lmo = 12 - LVT_rc%lmo 
          endif
          if(LVT_rc%lmo.eq.2) then 
             if((mod(LVT_rc%lyr,4) .eq. 0 .and. &
                  mod(LVT_rc%lyr, 100).ne.0) &!leap year
                  .or.(mod(LVT_rc%lyr,400) .eq.0)) then 
                if(LVT_rc%da.gt.29) then 
                   LVT_rc%lda = 29
                else
                   LVT_rc%lda = LVT_rc%da
                endif
             else 
                if(LVT_rc%da.gt.28) then 
                   LVT_rc%lda = 28
                else
                   LVT_rc%lda = LVT_rc%da
                endif
             endif
          else
             LVT_rc%lda = LVT_rc%da
          endif

          LVT_rc%lhr = LVT_rc%hr
          LVT_rc%lmn = LVT_rc%mn 
          LVT_rc%lss = LVT_rc%ss 
       endif

       call ESMF_TimeGet(lagTime,  yy=LVT_rc%lyr, &
            mm = LVT_rc%lmo, &
            dd = LVT_rc%lda, &
            h = LVT_rc%lhr, &
            m = LVT_rc%lmn, &
            s = LVT_rc%lss,&
            calendar = LVT_calendar, &
            rc=status)
       call LVT_date2time(LVT_rc%ltime,LVT_rc%ldoy,LVT_rc%lgmt,&
            LVT_rc%lyr,LVT_rc%lmo,LVT_rc%lda,LVT_rc%lhr,LVT_rc%lmn,LVT_rc%lss)


       lagTime = lagTime + ts
       call ESMF_TimeGet(currTime,  yy=LVT_rc%l_nyr, &
            mm = LVT_rc%l_nmo, &
            dd = LVT_rc%l_nda, &
            h = LVT_rc%l_nhr, &
            m = LVT_rc%l_nmn, &
            s = LVT_rc%l_nss,&
            calendar = LVT_calendar, &
            rc=status)
       call LVT_verify(status, 'error in LIS_ticktime')  
       
       LVT_rc%dyr(1)   = LVT_rc%yr
       LVT_rc%dmo(1)   = LVT_rc%mo
       LVT_rc%dda(1)   = LVT_rc%da
       LVT_rc%dhr(1)   = LVT_rc%hr
       LVT_rc%dmn(1)   = LVT_rc%mn
       LVT_rc%dss(1)   = LVT_rc%ss
       LVT_rc%ddoy(1)  = LVT_rc%doy
       LVT_rc%dgmt(1)  = LVT_rc%gmt
       LVT_rc%dtime(1) = LVT_rc%time

       LVT_rc%d_nyr(1)   = LVT_rc%yr
       LVT_rc%d_nmo(1)   = LVT_rc%mo
       LVT_rc%d_nda(1)   = LVT_rc%da
       LVT_rc%d_nhr(1)   = LVT_rc%hr
       LVT_rc%d_nmn(1)   = LVT_rc%mn
       LVT_rc%d_nss(1)   = LVT_rc%ss

       LVT_rc%dyr(2)   = LVT_rc%lyr
       LVT_rc%dmo(2)   = LVT_rc%lmo
       LVT_rc%dda(2)   = LVT_rc%lda
       LVT_rc%dhr(2)   = LVT_rc%lhr
       LVT_rc%dmn(2)   = LVT_rc%lmn
       LVT_rc%dss(2)   = LVT_rc%lss
       LVT_rc%ddoy(2)  = LVT_rc%ldoy
       LVT_rc%dgmt(2)  = LVT_rc%lgmt
       LVT_rc%dtime(2) = LVT_rc%ltime

       LVT_rc%d_nyr(2)   = LVT_rc%l_nyr
       LVT_rc%d_nmo(2)   = LVT_rc%l_nmo
       LVT_rc%d_nda(2)   = LVT_rc%l_nda
       LVT_rc%d_nhr(2)   = LVT_rc%l_nhr
       LVT_rc%d_nmn(2)   = LVT_rc%l_nmn
       LVT_rc%d_nss(2)   = LVT_rc%l_nss

       LVT_rc%dyr(3)   = LVT_rc%lyr
       LVT_rc%dmo(3)   = LVT_rc%lmo
       LVT_rc%dda(3)   = LVT_rc%lda
       LVT_rc%dhr(3)   = LVT_rc%lhr
       LVT_rc%dmn(3)   = LVT_rc%lmn
       LVT_rc%dss(3)   = LVT_rc%lss
       LVT_rc%ddoy(3)  = LVT_rc%ldoy
       LVT_rc%dgmt(3)  = LVT_rc%lgmt
       LVT_rc%dtime(3) = LVT_rc%ltime

       LVT_rc%d_nyr(3)   = LVT_rc%l_nyr
       LVT_rc%d_nmo(3)   = LVT_rc%l_nmo
       LVT_rc%d_nda(3)   = LVT_rc%l_nda
       LVT_rc%d_nhr(3)   = LVT_rc%l_nhr
       LVT_rc%d_nmn(3)   = LVT_rc%l_nmn
       LVT_rc%d_nss(3)   = LVT_rc%l_nss

       if(LVT_rc%timeAvgOpt.eq.0) then !instantaneous variable 
          if(mod(LVT_rc%tavgInterval,31536000).eq.0)then 
             if(LVT_rc%nmo.eq.LVT_rc%use_shift_mo) then 
                if(LVT_rc%nyr.ne.LVT_rc%prev_yr_tavg) then 
                   LVT_rc%prev_yr_tavg = LVT_rc%nyr
                   LVT_rc%computeFlag = .true. 
                endif
             endif
          elseif(mod(LVT_rc%tavgInterval,15552000).eq.0)then !6 monthly 
             if(LVT_rc%nmo.ne.LVT_rc%prev_mo_tavg) then 
                LVT_rc%prev_mo_tavg = LVT_rc%nmo
                LVT_rc%monthCount = LVT_rc%monthCount + 1
                if(LVT_rc%monthCount.eq.6) then 
                   LVT_rc%computeFlag = .true. 
                   LVT_rc%monthCount = 0 
                endif
             endif
          elseif(mod(LVT_rc%tavgInterval,7776000).eq.0)then !3 monthly 
             if(LVT_rc%nmo.ne.LVT_rc%prev_mo_tavg) then 
                LVT_rc%prev_mo_tavg = LVT_rc%nmo
                LVT_rc%monthCount = LVT_rc%monthCount + 1
                if(LVT_rc%monthCount.eq.3) then 
                   LVT_rc%computeFlag = .true. 
                   LVT_rc%monthCount = 0 
                endif
             endif
          elseif(mod(LVT_rc%tavgInterval,2592000).eq.0)then              
             if(LVT_rc%nmo.ne.LVT_rc%prev_mo_tavg) then 
                LVT_rc%prev_mo_tavg = LVT_rc%nmo
                LVT_rc%computeFlag = .true. 
             endif
             
          elseif(mod(LVT_rc%tavgInterval,604800).eq.0) then 
             if(mod(real(LVT_rc%nhr)*3600+60*real(LVT_rc%nmn)+&
                  float(LVT_rc%nss),&
                  real(LVT_rc%tavgInterval)).eq.0) then        
                LVT_rc%dayCount = LVT_rc%dayCount + 1
                if(LVT_rc%dayCount.eq.7) then 
                   LVT_rc%computeFlag = .true. 
                   LVT_rc%dayCount = 0 
                endif
             endif
          elseif(LVT_rc%tavgInterval.le.86400) then 
             if(mod(real(LVT_rc%nhr)*3600+60*real(LVT_rc%nmn)+&
                  float(LVT_rc%nss),&
                  real(LVT_rc%tavgInterval)).eq.0) then        
                LVT_rc%computeFlag = .true.
             else
                LVT_rc%computeFlag = .false. 
             endif
          else
             write(LVT_logunit,*) '[ERR] The support for the specified averaging'
             write(LVT_logunit,*) '[ERR] is not supported currently. '             
             write(LVT_logunit,*) '[ERR] Please contact the LVT development team.'
             call LVT_endrun()
          endif
          if(LVT_rc%endtime.eq.1) then 
             LVT_rc%computeFlag = .true. 
          endif

       elseif(LVT_rc%timeAvgOpt.eq.1.or.LVT_rc%timeAvgOpt.eq.3) then !time averaged variable

          if(mod(LVT_rc%tavgInterval,31536000).eq.0)then 
             if(LVT_rc%mo.eq.LVT_rc%use_shift_mo) then                 
                if(LVT_rc%yr.ne.LVT_rc%prev_yr_tavg) then 
                   LVT_rc%prev_yr_tavg = LVT_rc%yr
                   LVT_rc%computeFlag = .true. 
                endif
             endif
          elseif(mod(LVT_rc%tavgInterval,15552000).eq.0)then !6 monthly 
             if(LVT_rc%mo.ne.LVT_rc%prev_mo_tavg) then 
                LVT_rc%prev_mo_tavg = LVT_rc%mo
                LVT_rc%monthCount = LVT_rc%monthCount + 1
                if(LVT_rc%monthCount.eq.6) then 
                   LVT_rc%computeFlag = .true. 
                   LVT_rc%monthCount = 0 
                endif
             endif
          elseif(mod(LVT_rc%tavgInterval,7776000).eq.0)then !3 monthly 
             if(LVT_rc%mo.ne.LVT_rc%prev_mo_tavg) then 
                LVT_rc%prev_mo_tavg = LVT_rc%mo
                LVT_rc%monthCount = LVT_rc%monthCount + 1
                if(LVT_rc%monthCount.eq.3) then 
                   LVT_rc%computeFlag = .true. 
                   LVT_rc%monthCount = 0 
                endif
             endif
          elseif(mod(LVT_rc%tavgInterval,2592000).eq.0)then              
             if(LVT_rc%mo.ne.LVT_rc%prev_mo_tavg) then 
                LVT_rc%prev_mo_tavg = LVT_rc%mo
                LVT_rc%computeFlag = .true. 
             endif
             
          elseif(mod(LVT_rc%tavgInterval,604800).eq.0) then 
             if(mod(real(LVT_rc%hr)*3600+60*real(LVT_rc%mn)+float(LVT_rc%ss),&
                  real(LVT_rc%tavgInterval)).eq.0) then        
                LVT_rc%dayCount = LVT_rc%dayCount + 1
                if(LVT_rc%dayCount.eq.7) then 
                   LVT_rc%computeFlag = .true. 
                   LVT_rc%dayCount = 0 
                endif
             endif
          elseif(LVT_rc%tavgInterval.le.86400) then 
             if(mod(real(LVT_rc%hr)*3600+60*real(LVT_rc%mn)+float(LVT_rc%ss),&
                  real(LVT_rc%tavgInterval)).eq.0) then        
                LVT_rc%computeFlag = .true.
             else
                LVT_rc%computeFlag = .false. 
             endif
          else
             write(LVT_logunit,*) '[ERR] The support for the specified averaging'
             write(LVT_logunit,*) '[ERR] is not supported currently. '             
             write(LVT_logunit,*) '[ERR] Please contact the LVT development team.'
             call LVT_endrun()
          endif
          if(LVT_rc%endtime.eq.1) then 
             LVT_rc%computeFlag = .true. 
          endif
       endif
    elseif(LVT_rc%tsconv.eq."dekad") then !dekad
       !DOUBLE CHECK THIS to be consistent with the way time is advanced at
       ! the beginning of the loop
       !          LVT_rc%chkTS = checkTS(secsFrmBegYr=secsFrmBegYr(LVT_rc%doy, &
       !               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss), &
       !               yr=LVT_rc%yr, timeSchema_arg=gDEKAD_TIME_SCHEMA, &
       !               dekOy_arg=dekOy)
       !          
       !          if(LVT_rc%chkTS) then 
       !             LVT_rc%computeFlag = .true.
       !          endif
       LVT_rc%chkTS = .true. 
       call ESMF_TimeSet(currTime,  yy=LVT_rc%yr, &
            mm = LVT_rc%mo, &
            dd = LVT_rc%da, &
            h = LVT_rc%hr, &
            m = LVT_rc%mn, &
            s = LVT_rc%ss,&
            calendar = LVT_calendar, &
            rc=status)
       call LVT_verify(status, 'error in LVT_ticktime')
       
       call ESMF_TimeIntervalSet(ts,s=LVT_rc%ts,rc=status)
       call LVT_verify(status, 'error in LVT_ticktime')
       
       if((mod(LVT_rc%yr,4) .eq. 0 .and. mod(LVT_rc%yr, 100).ne.0) &!leap year
            .or.(mod(LVT_rc%yr,400) .eq.0)) then 
          days(2) = 29
       else 
          days(2) = 28
       endif
       
       call ESMF_TimeGet(currTime,yy=yr,mm=mo,dd=da,h=hr,m=mn,s=ss,rc=status)

       currTime = currTime + ts
       call ESMF_TimeGet(currTime,yy=yr,mm=mo,dd=da,h=hr,m=mn,s=ss,rc=status)

       call ESMF_TimeSet(start_dekad1, yy=yr, mm=mo, & 
            dd=1, h=0,  m=0, s=0, rc=status)       
       call LVT_verify(status, 'error in LVT_ticktime')
       call ESMF_TimeSet(start_dekad2, yy=yr, mm=mo, &
            dd=11, h=0, m=0, s=0, rc=status)       
       call LVT_verify(status, 'error in LVT_ticktime')
       call ESMF_TimeSet(start_dekad3, yy=yr, mm=mo, &
            dd=21, h=0, m=0, s=0, rc=status)
       call LVT_verify(status, 'error in LVT_ticktime')
       

       if(currTime < start_dekad1) then 
          call ESMF_TimeSet(advTime, yy=yr, mm=mo, dd=1, &
               h=hr, m=mn, s=ss, rc=status)
          call LVT_verify(status, 'error in LVT_ticktime')
       elseif(currTime .ge. start_dekad1 .and. currTime .lt. start_dekad2) then 
          call ESMF_TimeSet(advTime, yy=yr, mm=mo, dd=10, &
               h=hr, m=mn, s=ss, rc=status)
          call LVT_verify(status, 'error in LVT_ticktime')
       elseif(currTime .ge. start_dekad2 .and. currTime .lt. start_dekad3) then 
          call ESMF_TimeSet(advTime, yy=yr, mm=mo, dd=20, &
               h=hr, m=mn, s=ss, rc=status)
          call LVT_verify(status, 'error in LVT_ticktime')
       elseif(currTime.ge.start_dekad3) then 
          call ESMF_TimeSet(advTime, yy=yr, mm=mo, &
               dd=days(mo), &
               h=hr, m=mn, s=ss, rc=status)
          call LVT_verify(status, 'error in LVT_ticktime')
       endif
       call ESMF_TimeGet(advTime,  yy=LVT_rc%nyr, &
            mm = LVT_rc%nmo, &
            dd = LVT_rc%nda, &
            h = LVT_rc%nhr, &
            m = LVT_rc%nmn, &
            s = LVT_rc%nss,&
            calendar = LVT_calendar, &
            rc=status)
       call LVT_verify(status, 'error in LVT_ticktime')       
    
       if(LVT_rc%endtime.eq.1) then 
          LVT_rc%computeFlag = .true. 
       endif
       
       LVT_rc%computeFlag = .true.
       if(LVT_rc%tavgInterval.ne.864000) then 
          write(LVT_logunit,*) '[ERR] Only dekad interval for time averaging '
          write(LVT_logunit,*) '[ERR] is supported when reading LIS outputs '
          write(LVT_logunit,*) '[ERR] in dekad intervals'
          call LVT_endrun()
       endif
    end if

    ! EMK...For 557 post mode, alarm should be set when ready to write.
    if (LVT_rc%runmode .eq. "557 post") then
       LVT_rc%computeFlag = LVT_557post_alarm_is_on()
    end if
  end subroutine LVT_ticktime

!BOP
! 
! !ROUTINE: secsFrmBegYr
!  \label{secsFrmBegYr}
! 
! !INTERFACE: 
  function secsFrmBegYr(doy, hr, mm, ss)
! !ARGUMENTS:     
    integer*4 :: secsFrmBegYr    
    integer, intent(in) :: doy
    integer, intent(in) :: hr
    integer, intent(in) :: mm
    integer, intent(in) :: ss
!
! !DESCRIPTION: 
!  This routine returns the total number of seconds from the beginning
!  of the year. 
! 
!  The arguments are: 
!  \begin{description}
!    \item[doy]
!       julian day of the year
!    \item[hr]
!       hour value
!    \item[mm]
!       minute value
!    \item[ss]
!       second value
!    \item[secsFrmBegYr]
!       the returned seconds from the beginning of the year value
!  \end{description}
!
!EOP 
    
    secsFrmBegYr = (86400*doy) + (3600*hr) + (60*mm) + ss
    
  end function secsFrmBegYr

!BOP
! 
! !ROUTINE: checkTS
! \label{checkTS}
!
! !INTERFACE: 
  function checkTS(secsFrmBegYr, yr, &
       timeSchema_arg, dekOy_arg)
! !ARGUMENTS:     
   logical                            :: checkTS
   integer, intent(in)                :: yr
   integer*4, intent(in), optional    :: timeSchema_arg
   integer*4, intent(inout), optional :: dekOy_arg
!
! !DESCRIPTION: 
! 
!EOP

   integer*4, parameter :: gDEKAD_TIMESTEPS_PER_YEAR = 36
   integer*4, parameter :: gDEKAD_TIME_SCHEMA = 10

   integer*4, intent(in) :: secsFrmBegYr



   integer*4 :: timeSchema 
   integer*4 :: ts
   integer*4, save :: dekOy = (0-1)
!!! Dekad time schema variables begin !!!
! Remove the guess work re. ts start hour:  We should by default always be using offset 0600 per corresp. Jul 2011
#ifdef NO_OFFSET_BY_0600
#undef NO_OFFSET_BY_0600
#endif
   ! Dekad time step seconds from beginning of the year non- leap year
   integer*4, parameter, dimension(gDEKAD_TIMESTEPS_PER_YEAR) :: dekTSsfbyNLY = (/&
#ifndef NO_OFFSET_BY_0600
885600  ,&
1749600 ,&
2700000 ,&
3564000 ,&
4428000 ,&
5119200 ,&
5983200 ,&
6847200 ,&
7797600 ,&
8661600 ,&
9525600 ,&
10389600,&
11253600,&
12117600,&
13068000,&
13932000,&
14796000,&
15660000,&
16524000,&
17388000,&
18338400,&
19202400,&
20066400,&
21016800,&
21880800,&
22744800,&
23608800,&
24472800,&
25336800,&
26287200,&
27151200,&
28015200,&
28879200,&
29743200,&
30607200,&
31557600 &
#else
864000  ,&
1728000 ,&
2678400 ,&
3542400 ,&
4406400 ,&
5097600 ,&
5961600 ,&
6825600 ,&
7776000 ,&
8640000 ,&
9504000 ,&
10368000,&
11232000,&
12096000,&
13046400,&
13910400,&
14774400,&
15638400,&
16502400,&
17366400,&
18316800,&
19180800,&
20044800,&
20995200,&
21859200,&
22723200,&
23587200,&
24451200,&
25315200,&
26265600,&
27129600,&
27993600,&
28857600,&
29721600,&
30585600,&
31536000 &
#endif
/)
   ! Dekad time step seconds from beginning of the year leap year
   integer*4, parameter, dimension(gDEKAD_TIMESTEPS_PER_YEAR) :: dekTSsfbyLY = (/&
#ifndef NO_OFFSET_BY_0600
885600  ,&
1749600 ,&
2700000 ,&
3564000 ,&
4428000 ,&
5205600 ,&
6069600 ,&
6933600 ,&
7884000 ,&
8748000 ,&
9612000 ,&
10476000,&
11340000,&
12204000,&
13154400,&
14018400,&
14882400,&
15746400,&
16610400,&
17474400,&
18424800,&
19288800,&
20152800,&
21103200,&
21967200,&
22831200,&
23695200,&
24559200,&
25423200,&
26373600,&
27237600,&
28101600,&
28965600,&
29829600,&
30693600,&
31644000 &
#else
864000  ,&
1728000 ,&
2678400 ,&
3542400 ,&
4406400 ,&
5184000 ,&
6048000 ,&
6912000 ,&
7862400 ,&
8726400 ,&
9590400 ,&
10454400,&
11318400,&
12182400,&
13132800,&
13996800,&
14860800,&
15724800,&
16588800,&
17452800,&
18403200,&
19267200,&
20131200,&
21081600,&
21945600,&
22809600,&
23673600,&
24537600,&
25401600,&
26352000,&
27216000,&
28080000,&
28944000,&
29808000,&
30672000,&
31622400 &
#endif
/)
!!! Dekad time schema variables end !!!


   checkTS = .false.
   timeSchema = gDEKAD_TIME_SCHEMA ! assume dekad-time-schema* unless specified otherwise
   if(present(timeSchema_arg)) timeSchema = timeSchema_arg

   select case(timeSchema)
      case (gDEKAD_TIME_SCHEMA)

! From Tamuka Magadzire e-mail dated February 16, 2011:
!> For each month:
!> Dekad 1: the first 10 days of the month Dekad 2: the 11th to the 20th 
!> of the month Dekad 3: whatever days are left in the month after the 
!> 20th, whether its 8 days, 9, 10 or 11 days.

! Accordingly, see 'timestepTablesDekads2011.xls' (Excel 97-2003 workbook) 
! which was created to give the above table/s.

         if(  (mod(yr,4) /= 0) .or. ( (mod(yr,100) == 0) .and. (mod(yr,400) /= 0) )  ) then
            ! non- leap year

            do ts=1, gDEKAD_TIMESTEPS_PER_YEAR, 1
              if( secsFrmBegYr == dekTSsfbyNLY(ts) ) then
                checkTS = .true.
                dekOy = ts
              endif
            end do

         else
            ! is leap year

            do ts=1, gDEKAD_TIMESTEPS_PER_YEAR, 1
              if( secsFrmBegYr == dekTSsfbyLY(ts) ) then
                checkTS = .true.
                dekOy = ts
              endif
            end do

         endif
! Seconds from beginning of year is used instead of other temporal representations because:
!   1) this was given the most thought and therefore closest to completion at time of this writing
!   2) this representation represents the minimal requirements for any code that may use a time loop*
!   3) this exacting approach should require less logic in LIS for discrete time steps
!  *including both LIS and standalone model codes

      case default
         ! unhandled time schema
   end select

   if(present(dekOy_arg)) dekOy_arg = dekOy

end function checkTS 


!BOP
! 
! !ROUTINE: LVT_endofrun
! \label{LVT_endofrun}
!
! !INTERFACE:            
  function LVT_endofrun() result(finish)
! !ARGUMENTS:
    logical :: finish
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This function checks to see if the runtime clock has reached the
!  specified stop time of the simulation. 
! 
!   The arguments are:
!   \begin{description}
!   \item [finish]
!     boolean value indicating if the end of simulation is reached. 
!   \end{description}
!    
!  The calling sequence is: 
!  \begin{description}
!   \item[LVT\_is\_last\_step] (\ref{LVT_is_last_step}) \newline
!    check if the clock has reached the stop time
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer :: ierr

    if(LVT_masterproc) then
       finish = LVT_is_last_step(LVT_rc)
    endif
#if (defined SPMD)      
    call MPI_BCAST(finish, 1, MPI_LOGICAL, 0, & 
         MPI_COMM_WORLD, ierr)
#endif
  end function LVT_endofrun  

!BOP
! 
! !ROUTINE: spmd_init_offline
! \label{spmd_init_offline} 
!
! !INTERFACE:
  subroutine spmd_init_offline(lvtvm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Initializes MPI and retrieves the number of CPUs and the 
!  processors IDs. The MPI initialization is done in an
!  offline mode. In a coupled mode, the initialized
!  environment is passed to LVT from a parent component. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none
    type(ESMF_VM)   :: lvtvm
    integer         :: ier        ! return error status      '

    call ESMF_Initialize(vm=lvtvm,&
         defaultCalKind=ESMF_CALKIND_GREGORIAN,&
         logkindflag=ESMF_LOGKIND_NONE,rc=ier)
    call ESMF_VMGet(lvtvm,localPet=LVT_localPet,petCount=LVT_npes,&
         rc=ier)

    if (LVT_localPet==0) then 
       LVT_masterproc = .true.
    else
       LVT_masterproc = .false.
    end if
  end subroutine spmd_init_offline

!BOP
! 
! !ROUTINE: spmd_setup
! \label{spmd_setup} 
!
! !INTERFACE: 
  subroutine spmd_setup()
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  Allocates memory for the variables describing 
!  domain decomposition. 
!  
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    allocate(LVT_ews_ind(LVT_npes))
    allocate(LVT_ewe_ind(LVT_npes))
    allocate(LVT_nss_ind(LVT_npes))
    allocate(LVT_nse_ind(LVT_npes))

    allocate(LVT_lis_ews_ind(2,LVT_npes))
    allocate(LVT_lis_ewe_ind(2,LVT_npes))
    allocate(LVT_lis_nss_ind(2,LVT_npes))
    allocate(LVT_lis_nse_ind(2,LVT_npes))

    allocate(LVT_ews_halo_ind(LVT_npes))
    allocate(LVT_ewe_halo_ind(LVT_npes))
    allocate(LVT_nss_halo_ind(LVT_npes))
    allocate(LVT_nse_halo_ind(LVT_npes))

    allocate(LVT_lis_ews_halo_ind(2,LVT_npes))
    allocate(LVT_lis_ewe_halo_ind(2,LVT_npes))
    allocate(LVT_lis_nss_halo_ind(2,LVT_npes))
    allocate(LVT_lis_nse_halo_ind(2,LVT_npes))

    allocate(LVT_deltas(0:LVT_npes-1))
    allocate(LVT_offsets(0:LVT_npes-1))
    allocate(LVT_lis_deltas(2,0:LVT_npes-1))
    allocate(LVT_lis_offsets(2,0:LVT_npes-1))
    allocate(LVT_tdeltas(0:LVT_npes-1))
    allocate(LVT_toffsets(0:LVT_npes-1))
    allocate(LVT_lis_tdeltas(2,0:LVT_npes-1))
    allocate(LVT_lis_toffsets(2,0:LVT_npes-1))
    allocate(LVT_gdeltas(0:LVT_npes-1))
    allocate(LVT_goffsets(0:LVT_npes-1))
    allocate(LVT_lis_gdeltas(2,0:LVT_npes-1))
    allocate(LVT_lis_goffsets(2,0:LVT_npes-1))
    allocate(LVT_odeltas(0:LVT_npes-1))
    allocate(LVT_ooffsets(0:LVT_npes-1))
    allocate(LVT_ntiless(0:LVT_npes-1))
    allocate(LVT_ngrids(0:LVT_npes-1))
    allocate(LVT_lis_ntiless(2,0:LVT_npes-1))
    allocate(LVT_lis_ngrids(2,0:LVT_npes-1))
    
    LVT_ews_ind = 0
    LVT_ewe_ind = 0
    LVT_nss_ind = 0 
    LVT_nse_ind = 0 

    LVT_ews_halo_ind = 0
    LVT_ewe_halo_ind = 0
    LVT_nss_halo_ind = 0 
    LVT_nse_halo_ind = 0 

    LVT_lis_ews_halo_ind = 0
    LVT_lis_ewe_halo_ind = 0
    LVT_lis_nss_halo_ind = 0 
    LVT_lis_nse_halo_ind = 0 

    LVT_deltas = 0
    LVT_offsets = 0
    LVT_lis_deltas = 0
    LVT_lis_offsets = 0
    LVT_tdeltas = 0
    LVT_toffsets = 0
    LVT_lis_tdeltas = 0
    LVT_lis_toffsets = 0
    LVT_gdeltas = 0
    LVT_goffsets = 0
    LVT_lis_gdeltas = 0
    LVT_lis_goffsets = 0
    LVT_odeltas = 0
    LVT_ooffsets = 0
    LVT_ntiless = 0
    LVT_ngrids = 0
    LVT_lis_ntiless = 0
    LVT_lis_ngrids = 0
  end subroutine spmd_setup

!BOP
! 
! !ROUTINE: LVT_isAtAfinerResolution
! 
! !INTERFACE:
  function LVT_isAtAfinerResolution(datares) result(finish)
! !ARGUMENTS:
    logical :: finish
!
! !DESCRIPTION:
! 
!  This function checks if the given resolution is finer than 
!  the LIS/target resolution. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[datares]
!     input data resolution (in degrees)
!  \item[finish]
!     logical flag (true if the input resolution is finer than the 
!     LVT analysis resolution; else false)
!  \end{description}
! 

!EOP
    real    :: datares

    
    if(LVT_rc%gridDesc(1).eq.0) then 
       if(LVT_rc%gridDesc(10).le.datares) then 
          finish = .true.
       else
          finish = .false. 
       endif
    else
       write(LVT_logunit,*) '[ERR] LVT_isAtAfinerResolution check not supported '
       write(LVT_logunit,*) '[ERR] for this map projection'
       call LVT_endrun()
    endif
  end function LVT_isAtAfinerResolution

!BOP
! !ROUTINE: LVT_core_init
! \label{LVT_core_init}
!
! !INTERFACE: 
  subroutine LVT_core_init
! !DESCRIPTION: 
!
! Completes the initialization of the LVT' clock and alarms.
!EOP

    call LVT_update_clock(float(LVT_rc%ts))

    if(LVT_rc%tavgInterval.lt.LVT_rc%ts) then 
       write(LVT_logunit,*) '[ERR] Time averaging interval is less than the LIS output interval'
       call LVT_endrun()
    endif
  end subroutine LVT_core_init

  ! EMK...Return logical indicating if alarm should ring.
  ! Used by "557 post" runmode.
  logical function LVT_557post_alarm_is_on() result(alarmCheck)
     use LVT_timeMgrMod,      only : LVT_get_julhr
     implicit none

     logical, save :: firstTime = .true.
     integer, save :: starttime = 0
     integer :: curtime
     integer :: difftime

     alarmCheck = .false.
     if (firstTime) then
        call LVT_get_julss(LVT_rc%syr, LVT_rc%smo, LVT_rc%sda, &
             LVT_rc%shr, LVT_rc%smn, LVT_rc%sss, &
             starttime)
        firstTime = .false.
     end if
     call LVT_get_julss(LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
          LVT_rc%hr, LVT_rc%mn, LVT_rc%ss, &
          curtime)
     difftime = curtime - starttime
     if (difftime .gt. 0) then
        if (mod(difftime,LVT_rc%statswriteint).eq.0) then
           alarmCheck = .true.
        end if
     end if

     return
  end function LVT_557post_alarm_is_on

end module LVT_coreMod
