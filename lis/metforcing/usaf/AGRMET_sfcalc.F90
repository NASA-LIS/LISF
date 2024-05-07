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
! !ROUTINE: AGRMET_sfcalc
! \label{AGRMET_sfcalc}
! 
! !INTERFACE: 
subroutine AGRMET_sfcalc(n) 
! !USES: 
  use LIS_coreMod,         only : LIS_rc, LIS_masterproc, LIS_localPet, &
       LIS_domain, LIS_gdeltas, LIS_goffsets, LIS_npes, &
       LIS_ews_halo_ind, LIS_ewe_halo_ind, &
       LIS_nss_halo_ind, LIS_nse_halo_ind, &
       LIS_ews_ind, LIS_ewe_ind, LIS_nss_ind, LIS_nse_ind
  use LIS_LMLCMod,         only : LIS_LMLC
  use LIS_topoMod,         only : LIS_topo
  use LIS_timeMgrMod,      only : LIS_get_julhr, LIS_julhr_date
  use LIS_logMod,          only : LIS_logunit
  use LIS_mpiMod 
  use LIS_historyMod,      only : LIS_gather_2d_local_to_global
!  use LIS_fileIOMod,       only : LIS_putget
  use AGRMET_forcingMod,   only : agrmet_struc
  use USAF_bratsethMod, only: USAF_ObsData, USAF_createObsData, &
       USAF_setbratsethscreenstats, USAF_interpbacktotypeobsdata, &
       USAF_analyzescreen, USAF_multobsdata, USAF_analyzescreen
  use USAF_OBAMod, only: OBA, newOBA,makefilename,writetofile,&
       destroyOBA

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! This routine generates the surface fields of temperature, pressure,
! winds and relative humidity. It reads the GFS data, interpolates 
! first guess height, temperature, moisture, and wind data to the 
! AGRMET grid. These fields are vertically interpolated to each point's
! elevation. Finally the routine blends the surface observations of 
! temperature, relative humidity and surface wind speed with the 
! first guess fields using a barnes analysis. 
!
! !REVISION HISTORY:
! 30 AUG 2010 Modified to include coastal obs that were previously ignored
!             due to resolution of land mask (not resolving coastal cities and
!             obs locations as land, so flagging them as water); also,
!             modified to fix bug in conditional for range of Barnes spread
!             vector indices, and to adapt spread to land mask and resolution
!             flexibility (i.e., spreading of values over similar distances
!             regardless of number of grid points, e.g.)...Michael Shaw/WXE
! 10 SEP 2010 Modified to use variable size arrays for surface observation 
!             data.................................Chris Franks/16WS/WXE/SEMS
! 14 Jun 2017 Added GFS/GALWEM ground height, 2-m T, RH........Eric Kemp/GSFC
! 16 Dec 2021 Replaced julhr with YYYYMMDD in log.........Eric Kemp/NASA/SSAI
! 
! The arguments and variables are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[yr1,mo1,da1,hr1,mn1,ss1,yr,mo,da,hr,mn,ss]
!    time/date variables
!  \item[kprs]
!    number of isobaric levels for first guess data
!  \item[prslvls]
!    isobaric levels for first guess data
!  \item[alert\_number]
!    number of alerts that occur in the program
!  \item[timetoReadGFS]
!    flag to check for GFS reading frequency
!  \item[radrlh]
!    radius of influence- for each AGRMET grid point-
!    of wind speed observations in barnes analysis
!  \item[radspd]
!    radius of influence- for each AGRMET grid point-
!    of wind speed observations in barnes analysis
!  \item[radtmp]
!    radius of influence- for each AGRMET grid point-
!    of temperature observations in barnes analysis
!  \item[lokrlh]
!    standard radii of influence of relative humidity
!    observations in barnes analysis
!  \item[lokspd]
!    standard radii of influence of wind speed 
!    observations in barnes analysis
!  \item[loktmp]
!    standard radii of influence of temperature 
!    observations in barnes analysis
!  \item[ri]
!    i coordinate of surface observations on the AGRMET grid
!  \item[rj]
!    j coordinate of surface observations on the AGRMET grid
!  \item[obsrlh]
!    array of rlh observations
!  \item[obsspd]
!    array of wind speed observations
!  \item[obstmp]
!    array of temperature observations
!  \item[obscnt]
!    number of retrieved surface observations
!  \item[jul,julhr,cur\_jul,prev\_jul]
!   julian hour variables
!  \item[i,j]
!   loop indices
!  \item[order]
!   flag to indicate which data is to be read
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[find\_agrfld\_starttime](\ref{find_agrfld_starttime}) \newline
!    computes the start time for the FLDBLD processing.
!  \item[AGRMET\_fldbld] (\ref{AGRMET_fldbld}) \newline
!    generates the surface fields 
!  \item[LIS\_get\_julhr] (\ref{LIS_get_julhr}) \newline
!    converts the date to a julian hour
!  \item[AGRMET\_fgfill](\ref{AGRMET_fgfill}) \newline
!    determine a first guess field for the current time
!  \item[AGRMET\_getsfc](\ref{AGRMET_getsfc}) \newline
!    retrieve surface observations of temperature, relative
!    humidity, and wind speed from CDMS 
!  \item[AGRMET\_sfcalc\_barnes](\ref{AGRMET_sfcalc_barnes}) \newline
!    barnes optimal interpolation 
!  \end{description}
!EOP

  integer                 :: kprs 
  integer                 :: alert_number = 0
  logical                 :: timeToReadGFS
  integer,   allocatable  :: radrlh      ( : , : )
  integer,   allocatable  :: radspd      ( : , : )
  integer,   allocatable  :: radtmp      ( : , : )
  integer                 :: lokrlh      ( 2:6 )
  integer                 :: lokspd      ( 2:6 )
  integer                 :: loktmp      ( 2:6 )
  real,      allocatable  :: ri          ( : )
  real,      allocatable  :: rj          ( : )
  real,      allocatable  :: obsrlh      ( : )
  real,      allocatable  :: obsspd      ( : )
  real,      allocatable  :: obstmp      ( : )  
  integer                 :: obscnt
  integer                 :: julhr
  integer                 :: istart
  integer                 :: julend
  integer                 :: i, j
  integer                 :: order
  integer                 :: step
  integer                 :: step_dum
  integer                 :: jultmp
  real                    :: half_degree_ratio
  logical                 :: bounds_check
  ! EMK
  type(USAF_ObsData) :: t2mObs, rh2mObs, spd10mObs
!<rm -- jim merge>
#if 0
  real, allocatable :: t2m_back(:,:), rh2m_back(:,:), spd10m_back(:,:)
  real, allocatable :: tmp_back_1d_local(:), tmp_back_1d_glb(:)
  real, allocatable :: topo_glb(:,:)
#endif
!</rm -- jim merge>
  type(OBA) :: t2mOBA, rh2mOBA, spd10mOBA
  character(len=50) :: t2mPathOBA,rh2mPathOBA,spd10mPathOBA
  character(len=120) :: obaFilename
  character(len=10) :: yyyymmddhh
  integer :: ierr
  integer :: r,c, L
  character(len=32) :: type
  integer :: gdeltas, gid, ntiles
  integer :: count1
  logical :: found_inq
  integer :: rc
  integer, external :: LIS_create_subdirs
  logical :: use_wigos_sfcobs

  data lokspd     / 15, 25, 30, 40, 50 /
  data lokrlh     / 10, 15, 25, 35, 40 /
  data loktmp     / 10, 15, 20, 25, 30 /

  t2mPathOBA = "./t2m_OBA" ! EMK
  rh2mPathOBA = "./rh2m_OBA" ! EMK
  spd10mPathOBA = "./spd10m_OBA" ! EMK
  type = 'ALL' ! EMK

  ! See if subdirectories exist
  if (agrmet_struc(n)%oba_switch .eq. 1 .or. &
       agrmet_struc(n)%oba_switch .eq. 2) then
     if (LIS_masterproc) then
        inquire(file=trim(t2mPathOBA), exist=found_inq)
        if (.not. found_inq) then
           rc = lis_create_subdirs(len_trim(t2mPathOBA), trim(t2mPathOBA))
           if (rc .ne. 0) then
              write(LIS_logunit, *) &
                   '[WARN] Cannot create directory ', trim(t2mPathOBA)
           end if
        end if

        inquire(file=trim(rh2mPathOBA), exist=found_inq)
        if (.not. found_inq) then
           rc = lis_create_subdirs(len_trim(rh2mPathOBA), trim(rh2mPathOBA))
           if (rc .ne. 0) then
              write(LIS_logunit, *) &
                   '[WARN] Cannot create directory ', trim(rh2mPathOBA)
           end if
        end if

        inquire(file=trim(spd10mPathOBA), exist=found_inq)
        if (.not. found_inq) then
           rc = lis_create_subdirs(len_trim(spd10mPathOBA), &
                trim(spd10mPathOBA))
           if (rc .ne. 0) then
              write(LIS_logunit, *) &
                   '[WARN] Cannot create directory ', trim(spd10mPathOBA)
           end if
        end if

     end if
   end if
   
  half_degree_ratio=.5/LIS_rc%gridDesc(n,9)

  call LIS_get_julhr(LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,&
       0,0,jultmp)
  if(jultmp .gt.agrmet_struc(n)%lastSfcalcHour) then 
    timeToReadGFS = .true. 
 else
    timeToReadGFS = .false.
  endif

  if(timeToReadGFS) then 

!------------------------------------------------------------------    
! Find the time to start the processing from
!------------------------------------------------------------------    

     call find_agrfld_starttime(LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,istart)

     julend = istart+6
     
     agrmet_struc(n)%lastSfcalcHour = julend

!     ------------------------------------------------------------------
!    check if the current instance is a restart or continuation of a run
!    if the run is a cold start, read current and previous 6 hour data
!    else copy the current to previous and read the new current data. 
!     ------------------------------------------------------------------

     if(agrmet_struc(n)%findtime1.eq.1.and.agrmet_struc(n)%findtime2.eq.1) then
        order = 2
        call AGRMET_fldbld(n,order,istart)
        
        order = 1
        call AGRMET_fldbld(n,order,julend)
     else
        agrmet_struc(n)%agr_bgrd_src_p = agrmet_struc(n)%agr_bgrd_src_c
        agrmet_struc(n)%agr_tmp_p = agrmet_struc(n)%agr_tmp_c
        agrmet_struc(n)%agr_hgt_p = agrmet_struc(n)%agr_hgt_c
        agrmet_struc(n)%agr_rh_p = agrmet_struc(n)%agr_rh_c
        agrmet_struc(n)%agr_tmp_sfc_p = agrmet_struc(n)%agr_tmp_sfc_c
        agrmet_struc(n)%agr_hgt_sfc_p = agrmet_struc(n)%agr_hgt_sfc_c
        agrmet_struc(n)%agr_rh_sfc_p = agrmet_struc(n)%agr_rh_sfc_c
        agrmet_struc(n)%agr_wspd_p = agrmet_struc(n)%agr_wspd_c

        order = 1
        call AGRMET_fldbld(n,order,julend)
     endif
 
!EMK      
!     allocate ( radrlh (LIS_rc%lnc(n), LIS_rc%lnr(n)))
!     allocate ( radspd (LIS_rc%lnc(n), LIS_rc%lnr(n)))
!     allocate ( radtmp (LIS_rc%lnc(n), LIS_rc%lnr(n)))
  
!     ------------------------------------------------------------------
!         loop thru points in hemi. if grid pt is land, set spread
!         radii based on mean obs count (irad), otherwise set spread
!         radii to zero and reset altitude to zero
!     ------------------------------------------------------------------      
!      do j = 1,LIS_rc%lnr(n)
!         do i = 1,LIS_rc%lnc(n)
!            bounds_check = .false.
!            if(LIS_LMLC(n)%landmask(i,j).gt.0) then 
!               bounds_check = .true. 
!            endif
!            if(i.lt.LIS_rc%lnc(n)) then 
!               if(LIS_LMLC(n)%landmask(i+1,j).gt.0) then 
!                  bounds_check = .true.
!               endif
!            endif
!            if(j.lt.LIS_rc%lnr(n)) then 
!               if(LIS_LMLC(n)%landmask(i,j+1).gt.0) then 
!                  bounds_check = .true.
!               endif
!            endif
!            if(i.lt.LIS_rc%lnc(n).and.j.lt.LIS_rc%lnr(n)) then 
!               if(LIS_LMLC(n)%landmask(i+1,j+1).gt.0) then 
!                  bounds_check = .true.
!               endif
!            endif
!            if(i.gt.1) then 
!               if(LIS_LMLC(n)%landmask(i-1,j).gt.0) then 
!                  bounds_check = .true.
!               endif
!            endif
!            if(j.gt.1) then 
!               if(LIS_LMLC(n)%landmask(i,j-1).gt.0) then 
!                  bounds_check = .true.
!               endif
!            endif
!            if(i.gt.1.and.j.gt.1) then 
!               if(LIS_LMLC(n)%landmask(i-1,j-1).gt.0) then 
!                  bounds_check = .true.
!               endif
!            endif
! !           if( LIS_LMLC(n)%landmask(i,j)    .gt. 0 &
! !          .or. i .lt. LIS_rc%lnc(n) .and. LIS_LMLC(n)%landmask(i+1,j)                              .gt. 0 &
! !          .or. j .lt. LIS_rc%lnr(n) .and. LIS_LMLC(n)%landmask(i,j+1)                              .gt. 0 &
! !          .or. i .lt. LIS_rc%lnc(n) .and. j .lt. LIS_rc%lnr(n) .and. LIS_LMLC(n)%landmask(i+1,j+1) .gt. 0 &
! !          .or. i .gt. 1 .and. LIS_LMLC(n)%landmask(i-1,j)                                          .gt. 0 &
! !          .or. j .gt. 1 .and. LIS_LMLC(n)%landmask(i,j-1)                                          .gt. 0 &
! !          .or. i .gt. 1 .and. j .gt. 1 .and. LIS_LMLC(n)%landmask(i-1,j-1)                         .gt. 0 )then
!            if(bounds_check) then 
!               if( (agrmet_struc(n)%irad(i,j)/half_degree_ratio .lt. 2) .or. & !  .and. (agrmet_struc(n)%irad(i,j) .gt. 0) .or. &
!                   (agrmet_struc(n)%irad(i,j)/half_degree_ratio .ge. 6) )then
!                  radspd(i,j) = 1
!                  radrlh(i,j) = 1
!                  radtmp(i,j) = 1
!               else
!                  radspd(i,j) = half_degree_ratio*lokspd(nint(agrmet_struc(n)%irad(i,j)/half_degree_ratio))
!                  radrlh(i,j) = half_degree_ratio*lokrlh(nint(agrmet_struc(n)%irad(i,j)/half_degree_ratio))
!                  radtmp(i,j) = half_degree_ratio*loktmp(nint(agrmet_struc(n)%irad(i,j)/half_degree_ratio))
!               endif
!            else
! !              agrmet_struc(n)%alt(i,j)    = 0.0
!               radspd(i,j) = 0
!               radrlh(i,j) = 0
!               radtmp(i,j) = 0
!            endif
!         enddo
!      enddo

!     ------------------------------------------------------------------
!             allocate arrays for sfc observations
!     ------------------------------------------------------------------

     allocate ( ri (agrmet_struc(n)%max_sfcobs))
     allocate ( rj (agrmet_struc(n)%max_sfcobs))
     allocate ( obsrlh (agrmet_struc(n)%max_sfcobs))
     allocate ( obsspd (agrmet_struc(n)%max_sfcobs))
     allocate ( obstmp (agrmet_struc(n)%max_sfcobs))

!     ------------------------------------------------------------------
!             loop for the surface obs files
!     ------------------------------------------------------------------

     step = 1        
     do julhr= julend-5,julend

        call AGRMET_julhr_date10(julhr, yyyymmddhh) ! EMK
        write(LIS_logunit,*)' '
        write(LIS_logunit,*)'---------------------------- '
        !write(LIS_logunit,*)'- PROCESSING-SFC JULHR ', julhr
        write(LIS_logunit,*)'[INFO] PROCESSING-SFC YYYYMMDDHH ', yyyymmddhh
        write(LIS_logunit,*)'---------------------------- '

!     ------------------------------------------------------------------
!             fill and/or process the first guess fields
!     ------------------------------------------------------------------
           
!<rm -- jim merge>
#if 0
       allocate(tmp_back_1d_local(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(tmp_back_1d_glb(LIS_rc%gnc(n)*LIS_rc%gnr(n)))

       ! Create global topography array
       call MPI_Barrier(LIS_MPI_COMM, ierr)
       tmp_back_1d_local(:) = -9999.0
       do r = 1, LIS_rc%lnr(n)
          do c = 1, LIS_rc%lnc(n)
             if (LIS_domain(n)%gindex(c,r) .eq. -1) cycle
             tmp_back_1d_local(LIS_domain(n)%gindex(c,r)) = &
                  LIS_topo(n)%elevation(c,r,1)
          end do ! c
       end do ! r
       tmp_back_1d_glb(:) = -9999.0
       gdeltas = LIS_gdeltas(n,LIS_localPet)
       call MPI_ALLGATHERV(tmp_back_1d_local,gdeltas,MPI_REAL, &
            tmp_back_1d_glb,LIS_gdeltas(n,:),LIS_goffsets(n,:), MPI_REAL, &
            LIS_MPI_COMM,ierr)
       count1=1
       allocate(topo_glb(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       topo_glb(:,:) = -9999.0
       do L=1,LIS_npes
          do r=LIS_nss_halo_ind(n,L),LIS_nse_halo_ind(n,L)
             do c=LIS_ews_halo_ind(n,L),LIS_ewe_halo_ind(n,L)
                gid = c+(r-1)*LIS_rc%gnc(n)
                ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                if (ntiles.ne.0) then
                   if(r.ge.LIS_nss_ind(n,l).and.&
                        r.le.LIS_nse_ind(n,l).and.&
                        c.ge.LIS_ews_ind(n,l).and.&
                        c.le.LIS_ewe_ind(n,l)) then !points not in halo
                      topo_glb(c,r) = tmp_back_1d_glb(count1)
                   end if
                   count1 = count1 + 1
                end if
             end do ! c
          end do ! r
       end do ! L       
#endif
!</rm -- jim merge>

        kprs = 13
        step_dum = step

        call AGRMET_fgfill( n, LIS_topo(n)%elevation(:,:,1), &
             agrmet_struc(n)%sfctmp, &
             agrmet_struc(n)%sfcrlh, &
             agrmet_struc(n)%sfcspd, &
             agrmet_struc(n)%sfcprs, &
!             LIS_LMLC(n)%landmask, & ! EMK...Removed land mask
             julend, &
             agrmet_struc(n)%lastmp, &
             agrmet_struc(n)%lasrlh, &
             agrmet_struc(n)%lasspd, &
             agrmet_struc(n)%lasprs, &
             step, order,agrmet_struc(n)%minwnd, &
             LIS_rc%lnc(n), LIS_rc%lnr(n),&
             kprs, agrmet_struc(n)%agr_tmp_c, &
             agrmet_struc(n)%agr_hgt_c,&
             agrmet_struc(n)%agr_rh_c,&
             agrmet_struc(n)%agr_tmp_sfc_c,&
             agrmet_struc(n)%agr_hgt_sfc_c,&
             agrmet_struc(n)%agr_rh_sfc_c,&
             agrmet_struc(n)%agr_wspd_c,&
             agrmet_struc(n)%agr_pres_c,&
             agrmet_struc(n)%agr_tmp_p, &
             agrmet_struc(n)%agr_hgt_p,&
             agrmet_struc(n)%agr_rh_p,&
             agrmet_struc(n)%agr_tmp_sfc_p,&
             agrmet_struc(n)%agr_hgt_sfc_p,&
             agrmet_struc(n)%agr_rh_sfc_p,&
             agrmet_struc(n)%agr_wspd_p,&
             agrmet_struc(n)%agr_pres_p,&
             trim(agrmet_struc(n)%agrmetdir))
!<rm -- jim merge>
#if 0
        deallocate(topo_glb)
#endif
!</rm -- jim merge>

        do  i = 1, 6
           call LIS_gather_2d_local_to_global(n, agrmet_struc(n)%sfctmp(i,:,:), &
                                                 agrmet_struc(n)%sfctmp_glb(i,:,:))
           call LIS_gather_2d_local_to_global(n, agrmet_struc(n)%sfcrlh(i,:,:), &
                                                 agrmet_struc(n)%sfcrlh_glb(i,:,:))
           call LIS_gather_2d_local_to_global(n, agrmet_struc(n)%sfcspd(i,:,:), &
                                                 agrmet_struc(n)%sfcspd_glb(i,:,:))
        enddo


!     ------------------------------------------------------------------
!             retrieve any observed temps, rh's, and wind speeds for
!             this hemisphere and time from the surface database
!     ------------------------------------------------------------------

        ! EMK Specify the source of the background field before
        ! getting obs.  This is required so that the appropriate observation
        ! error variances are assigned.
        call USAF_setBratsethScreenStats(agrmet_struc(n)%agr_bgrd_src_p,n)

        ! EMK Added structures for surface obs
        !call USAF_createObsData(t2mObs,n,maxobs=20000)
        !call USAF_createObsData(rh2mObs,n,maxobs=20000)
        !call USAF_createObsData(spd10mObs,n,maxobs=20000)
        ! EMK Patch...Dynamically set the max obs sizes
        call USAF_createObsData(t2mObs, n, maxobs=agrmet_struc(n)%max_sfcobs)
        call USAF_createObsData(rh2mObs, n, maxobs=agrmet_struc(n)%max_sfcobs)
        call USAF_createObsData(spd10mObs, n, &
             maxobs=agrmet_struc(n)%max_sfcobs)

        if (agrmet_struc(n)%sfcobsfmt == 1) then
           use_wigos_sfcobs = .false.
        else
           use_wigos_sfcobs = .true.
        end if
        call AGRMET_getsfc( n, julhr, t2mObs, rh2mObs, spd10mObs, &
             ri, rj, obstmp, obsrlh, obsspd, &
             obscnt, agrmet_struc(n)%max_sfcobs, agrmet_struc(n)%minwnd, &
             alert_number, LIS_rc%lnc(n), LIS_rc%lnr(n),&
             agrmet_struc(n)%agrmetdir,agrmet_struc(n)%cdmsdir,&
             agrmet_struc(n)%use_timestamp, use_wigos_sfcobs)

!        call MPI_Barrier(LIS_mpi_comm, ierr)
!        stop
!     ------------------------------------------------------------------
!             for all parameters except surface pressure... use barnes
!             technique to fold observed values of this parameter into
!             the first guess surface field
!     ------------------------------------------------------------------
#if 0 
           tempvar = -9999.0
           do k=1,obscnt
              tempvar(nint(ri(k)),nint(rj(k))) = obsspd(k)
           enddo
           open(100,file='tempvar.bin',form='unformatted')
           write(100) tempvar
           close(100)
           stop
#endif

! EMK...Comment out
!         if( obscnt .gt. 0 )then

!            write(LIS_logunit,*)' '
!            write(LIS_logunit,*)'- PERFORM BARNES ANALYSIS FOR TEMPERATURE.'
           
!            call AGRMET_sfcalc_barnes( obscnt, obstmp, ri, rj, &
!                 agrmet_struc(n)%sfctmp(step_dum,:,:),radtmp, &
!                 'tp', LIS_LMLC(n)%landmask(:,:),&
!                 agrmet_struc(n)%minwnd,&
!                 agrmet_struc(n)%max_sfcobs, LIS_rc%lnc(n), LIS_rc%lnr(n))
           
!            write(LIS_logunit,*)' '
!            write(LIS_logunit,*)'- PERFORM BARNES ANALYSIS FOR RH.'
           
!            call AGRMET_sfcalc_barnes( obscnt, obsrlh, ri, rj, &
!                 agrmet_struc(n)%sfcrlh(step_dum,:,:),radrlh, &
!                 'rh', LIS_LMLC(n)%landmask, &
!                 agrmet_struc(n)%minwnd, &
!                 agrmet_struc(n)%max_sfcobs, LIS_rc%lnc(n), LIS_rc%lnr(n))

!            write(LIS_logunit,*)' '
!            write(LIS_logunit,*)'- PERFORM BARNES ANALYSIS FOR WIND SPEED.'
!            write(LIS_logunit,*)' '
              
!            call AGRMET_sfcalc_barnes( obscnt, obsspd, ri, rj, &
!                 agrmet_struc(n)%sfcspd(step_dum,:,:), radspd, &
!                 'ws', LIS_LMLC(n)%landmask(:,:), &
!                 agrmet_struc(n)%minwnd,&
!                 agrmet_struc(n)%max_sfcobs, LIS_rc%lnc(n), LIS_rc%lnr(n))
!        endif

!        ! EMK
!        call MPI_Barrier(LIS_MPI_COMM, ierr)
!        tmp_back_1d_local(:) = 0
!        do r = 1, LIS_rc%lnr(n)
!           do c = 1, LIS_rc%lnc(n)
!              if (LIS_domain(n)%gindex(c,r) .eq. -1) cycle
!              tmp_back_1d_local(LIS_domain(n)%gindex(c,r)) = &
!                   agrmet_struc(n)%sfctmp(step_dum,c,r)
!           end do ! c
!        end do ! r
!        tmp_back_1d_glb(:) = 0
!        gdeltas = LIS_gdeltas(n,LIS_localPet)
!        call MPI_ALLGATHERV(tmp_back_1d_local,gdeltas,MPI_REAL, &
!             tmp_back_1d_glb,LIS_gdeltas(n,:),LIS_goffsets(n,:), MPI_REAL, &
!             LIS_MPI_COMM,ierr)
!        count1=1
!        allocate(t2m_back(LIS_rc%gnc(n),LIS_rc%gnr(n)))
!        t2m_back(:,:) = 0
!        do L=1,LIS_npes
!           do r=LIS_nss_halo_ind(n,L),LIS_nse_halo_ind(n,L)
!              do c=LIS_ews_halo_ind(n,L),LIS_ewe_halo_ind(n,L)
!                 gid = c+(r-1)*LIS_rc%gnc(n)
!                 ntiles = LIS_domain(n)%ntiles_pergrid(gid)
!                 if (ntiles.ne.0) then
!                    if(r.ge.LIS_nss_ind(n,l).and.&
!                         r.le.LIS_nse_ind(n,l).and.&
!                         c.ge.LIS_ews_ind(n,l).and.&
!                         c.le.LIS_ewe_ind(n,l)) then !points not in halo
!                       t2m_back(c,r) = tmp_back_1d_glb(count1)
!                    end if
!                 end if
!              end do ! c
!           end do ! r
!        end do ! L       
       call USAF_interpBackToTypeObsData(t2mObs,n,LIS_rc%gnc(n),LIS_rc%gnr(n), &
             agrmet_struc(n)%sfctmp_glb(step_dum,:,:), type)

!        call MPI_Barrier(LIS_MPI_COMM, ierr)
!        tmp_back_1d_local(:) = 0
!        do r = 1, LIS_rc%lnr(n)
!           do c = 1, LIS_rc%lnc(n)
!              if (LIS_domain(n)%gindex(c,r) .eq. -1) cycle
!              tmp_back_1d_local(LIS_domain(n)%gindex(c,r)) = &
!                   agrmet_struc(n)%sfcrlh(step_dum,c,r)
!           end do ! c
!        end do ! r
!        tmp_back_1d_glb(:) = 0
!        gdeltas = LIS_gdeltas(n,LIS_localPet)
!        call MPI_ALLGATHERV(tmp_back_1d_local,gdeltas,MPI_REAL, &
!             tmp_back_1d_glb,LIS_gdeltas(n,:),LIS_goffsets(n,:), MPI_REAL, &
!             LIS_MPI_COMM,ierr)
!        count1=1
!        allocate(rh2m_back(LIS_rc%gnc(n),LIS_rc%gnr(n)))
!        rh2m_back(:,:) = 0
!        do L=1,LIS_npes
!           do r=LIS_nss_halo_ind(n,L),LIS_nse_halo_ind(n,L)
!              do c=LIS_ews_halo_ind(n,L),LIS_ewe_halo_ind(n,L)
!                 gid = c+(r-1)*LIS_rc%gnc(n)
!                 ntiles = LIS_domain(n)%ntiles_pergrid(gid)
!                 if (ntiles.ne.0) then
!                    if(r.ge.LIS_nss_ind(n,l).and.&
!                         r.le.LIS_nse_ind(n,l).and.&
!                         c.ge.LIS_ews_ind(n,l).and.&
!                         c.le.LIS_ewe_ind(n,l)) then !points not in halo
!                       rh2m_back(c,r) = tmp_back_1d_glb(count1)
!                    end if
!                 end if
!              end do ! c
!           end do ! r
!        end do ! L       
       call USAF_interpBackToTypeObsData(rh2mObs,n,LIS_rc%gnc(n),LIS_rc%gnr(n),&
            agrmet_struc(n)%sfcrlh_glb(step_dum,:,:),type)

!       print*,'EMK: LIS_localPet,maxval(rh2m_back) = ', &
!            LIS_localPet,maxval(rh2m_back)
!       print*,'EMK: LIS_localPet,minval(rh2m_back) = ', &
!            LIS_localPet,minval(rh2m_back)

!        call MPI_Barrier(LIS_MPI_COMM, ierr)
!        tmp_back_1d_local(:) = 0
!        do r = 1, LIS_rc%lnr(n)
!           do c = 1, LIS_rc%lnc(n)
!              if (LIS_domain(n)%gindex(c,r) .eq. -1) cycle
!              tmp_back_1d_local(LIS_domain(n)%gindex(c,r)) = &
!                   agrmet_struc(n)%sfcspd(step_dum,c,r)
!           end do ! c
!        end do ! r
!        tmp_back_1d_glb(:) = 0
!        gdeltas = LIS_gdeltas(n,LIS_localPet)
!        call MPI_ALLGATHERV(tmp_back_1d_local,gdeltas,MPI_REAL, &
!             tmp_back_1d_glb,LIS_gdeltas(n,:),LIS_goffsets(n,:), MPI_REAL, &
!             LIS_MPI_COMM,ierr)
!        count1=1
!        allocate(spd10m_back(LIS_rc%gnc(n),LIS_rc%gnr(n)))
!        spd10m_back(:,:) = 0
!        do L=1,LIS_npes
!           do r=LIS_nss_halo_ind(n,L),LIS_nse_halo_ind(n,L)
!              do c=LIS_ews_halo_ind(n,L),LIS_ewe_halo_ind(n,L)
!                 gid = c+(r-1)*LIS_rc%gnc(n)
!                 ntiles = LIS_domain(n)%ntiles_pergrid(gid)
!                 if (ntiles.ne.0) then
!                    if(r.ge.LIS_nss_ind(n,l).and.&
!                         r.le.LIS_nse_ind(n,l).and.&
!                         c.ge.LIS_ews_ind(n,l).and.&
!                         c.le.LIS_ewe_ind(n,l)) then !points not in halo
!                       spd10m_back(c,r) = tmp_back_1d_glb(count1)
!                    end if
!                 end if
!              end do ! c
!           end do ! r
!        end do ! L       
       call USAF_interpBackToTypeObsData(spd10mObs,n,LIS_rc%gnc(n), &
            LIS_rc%gnr(n),agrmet_struc(n)%sfcspd_glb(step_dum,:,:), type)

#if (defined SPMD)
       call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif

!<rm -- jim merge>
#if 0
       deallocate(tmp_back_1d_local)
       deallocate(tmp_back_1d_glb)
#endif
!</rm -- jim merge>

       ! EMK Run Bratseth scheme
       call AGRMET_julhr_date10(julhr, yyyymmddhh)
       write(LIS_logunit,*) &
            '[INFO] RUNNING BRATSETH TEMPERATURE ANALYSIS FOR ',yyyymmddhh
       call USAF_analyzeScreen(t2mObs,n,&
            agrmet_struc(n)%sfctmp_glb(step_dum,:,:),&
            agrmet_struc(n)%bratseth_t2m_back_sigma_b_sqr, &
            agrmet_struc(n)%bratseth_t2m_max_dist, &
            agrmet_struc(n)%bratseth_t2m_back_err_scale_length, &
            agrmet_struc(n)%sfctmp(step_dum,:,:),t2mOBA)

       ! Output OBA data
       if (agrmet_struc(n)%oba_switch .eq. 1 .or. &
            agrmet_struc(n)%oba_switch .eq. 2) then
          if (LIS_masterproc) then
             call AGRMET_julhr_date10(julhr, yyyymmddhh)
             call makeFilename(t2mPathOBA,yyyymmddhh,1,obaFilename)
             call writeToFile(t2mOBA,obaFilename)
          end if          
          call destroyOBA(t2mOBA)
       end if

       call AGRMET_julhr_date10(julhr, yyyymmddhh)
       write(LIS_logunit,*) &
            '[INFO] RUNNING BRATSETH RELATIVE HUMIDITY ANALYSIS FOR ',&
            yyyymmddhh
       ! Note:  For analysis, we use RH in percent.
       call USAF_multObsData(rh2mObs,100.)
       agrmet_struc(n)%sfcrlh_glb(step_dum,:,:) = &
            agrmet_struc(n)%sfcrlh_glb(step_dum,:,:)*100.
       call USAF_analyzeScreen(rh2mObs,n,&
            agrmet_struc(n)%sfcrlh_glb(step_dum,:,:),&
            agrmet_struc(n)%bratseth_rh2m_back_sigma_b_sqr, &
            agrmet_struc(n)%bratseth_rh2m_max_dist, &
            agrmet_struc(n)%bratseth_rh2m_back_err_scale_length, &
            agrmet_struc(n)%sfcrlh(step_dum,:,:),rh2mOBA)
       ! Make sure RH is between 0 and 1 (convert from percent).
       do r = 1, LIS_rc%lnr(n)
          do c = 1, LIS_rc%lnc(n)
             agrmet_struc(n)%sfcrlh(step_dum,c,r) = &
                  max(0.,min(1.,0.01*agrmet_struc(n)%sfcrlh(step_dum,c,r)))
          end do ! c
       end do ! r

       ! Output OBA data
       if (agrmet_struc(n)%oba_switch .eq. 1 .or. &
            agrmet_struc(n)%oba_switch .eq. 2) then
          if (LIS_masterproc) then
             call AGRMET_julhr_date10(julhr, yyyymmddhh)
             call makeFilename(rh2mPathOBA,yyyymmddhh,1,obaFilename)
             call writeToFile(rh2mOBA,obaFilename)
          end if
          call destroyOBA(rh2mOBA)       
       end if

       call AGRMET_julhr_date10(julhr, yyyymmddhh)
       write(LIS_logunit,*) &
            '[INFO] RUNNING BRATSETH WIND SPEED ANALYSIS FOR ',yyyymmddhh
       call USAF_analyzeScreen(spd10mObs,n,&
            agrmet_struc(n)%sfcspd_glb(step_dum,:,:),&
            agrmet_struc(n)%bratseth_spd10m_back_sigma_b_sqr, &
            agrmet_struc(n)%bratseth_spd10m_max_dist, &
            agrmet_struc(n)%bratseth_spd10m_back_err_scale_length, &
            agrmet_struc(n)%sfcspd(step_dum,:,:),spd10mOBA)

       ! Output OBA data
       if (agrmet_struc(n)%oba_switch .eq. 1 .or. &
            agrmet_struc(n)%oba_switch .eq. 2) then
          if (LIS_masterproc) then
             call AGRMET_julhr_date10(julhr, yyyymmddhh)
             call makeFilename(spd10mPathOBA,yyyymmddhh,1,obaFilename)
             call writeToFile(spd10mOBA,obaFilename)
          end if
          call destroyOBA(spd10mOBA)
       end if

!<rm -- jim merge>
#if 0
       ! Make sure pressure is copied to local array
       agrmet_struc(n)%sfcprs(step_dum,:,:) = &
            agrmet_struc(n)%sfcprs_glb(step_dum, &
            LIS_ews_halo_ind(n,LIS_localPet+1):LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1):LIS_nse_halo_ind(n,LIS_localPet+1))
#endif
!</rm -- jim merge>

       ! EMK Clean up
!       deallocate(t2m_back)
!       deallocate(rh2m_back)
!       deallocate(spd10m_back)
!       call USAF_destroyObsData(t2mObs)
!       call USAF_destroyObsData(rh2mObs)
!       call USAF_destroyObsData(spd10mObs)

#if 0 
!     ------------------------------------------------------------------
!    Save the analysis
!     ------------------------------------------------------------------     

        if(LIS_masterproc) then 
           write(LIS_logunit,*) 'Writing surface analysis to disk'
           call LIS_julhr_date(julhr,yr2,mo2,da2,hr2)
           call agrmet_sfctmp_filename(name_nh,agrmet_struc(n)%agrmetdir,&
                agrmet_struc(n)%sfcalcdir,agrmet_struc(n)%use_timestamp,&
                hemi,yr2,mo2,da2,hr2)
           call LIS_putget( agrmet_struc(n)%sfctmp(hemi,step_dum,:,:), &
                'w', name_nh, &
                routine_name, &
                agrmet_struc(n)%imax, agrmet_struc(n)%jmax )        
           call agrmet_sfcrlh_filename(name_nh,agrmet_struc(n)%agrmetdir,&
                agrmet_struc(n)%sfcalcdir,agrmet_struc(n)%use_timestamp,&
                hemi,yr2,mo2,da2,hr2)
           call LIS_putget( agrmet_struc(n)%sfcrlh(hemi,step_dum,:,:), &
                'w', name_nh, &
                routine_name ,&
                agrmet_struc(n)%imax, agrmet_struc(n)%jmax )        
           call agrmet_sfcspd_filename(name_nh,agrmet_struc(n)%agrmetdir,&
                agrmet_struc(n)%sfcalcdir,agrmet_struc(n)%use_timestamp,&
                hemi,yr2,mo2,da2,hr2)
           call LIS_putget( agrmet_struc(n)%sfcspd(hemi,step_dum,:,:), &
                'w', name_nh, &
                routine_name ,&
                agrmet_struc(n)%imax, agrmet_struc(n)%jmax )        
           call agrmet_sfcprs_filename(name_nh,agrmet_struc(n)%agrmetdir,&
                agrmet_struc(n)%sfcalcdir,agrmet_struc(n)%use_timestamp,&
                hemi,yr2,mo2,da2,hr2)
           call LIS_putget( agrmet_struc(n)%sfcprs(hemi,step_dum,:,:), &
                'w', name_nh, &
                routine_name ,&
                agrmet_struc(n)%imax, agrmet_struc(n)%jmax )        
        endif
#endif
     end do

!     ------------------------------------------------------------------
!     Deallocate observation arrays
!     ------------------------------------------------------------------

     deallocate ( ri )
     deallocate ( rj )
     deallocate ( obsrlh )
     deallocate ( obsspd )
     deallocate ( obstmp )
!     deallocate ( radrlh ) 
!     deallocate ( radspd ) 
!     deallocate ( radtmp )      

  endif
  
end subroutine AGRMET_sfcalc

!BOP
! 
! !ROUTINE: find_agrfld_starttime
! \label{find_agrfld_starttime}
! 
! !INTERFACE:
subroutine find_agrfld_starttime(yr,mo,da,hr,julbeg)        
! !USES: 
  use LIS_timeMgrMod, only : LIS_tick, LIS_get_julhr

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: yr
  integer, intent(in) :: mo
  integer, intent(in) :: da
  integer, intent(in) :: hr
  integer, intent(inout) :: julbeg
! 
! !DESCRIPTION: 
!  This routine finds the julian hour to start the AGRMET
!  precip processing from, based on the current input time. 
! 
!  The arguments are:
!  \begin{description}
!  \item[yr]
!    the current year
!  \item[mo]
!    the current month
!  \item[da]
!    the current day
!  \item[hr]
!    the current hour
!  \item[julbeg]
!    output starting julian hour
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_LIS\_tick] (\ref{LIS_tick}) \newline
!    computes previous 6 hour time. 
!  \item[LIS\_get\_julhr] (\ref{LIS_get_julhr}) \newline
!    converts the date to a julian hour
!  \end{description}
!EOP
  
  integer          :: yr1,mo1,da1,hr1,mn1,ss1,ts1
  real*8           :: time
  integer          :: doy
  real             :: gmt

  yr1 = yr
  mo1 = mo
  da1 = da
  mn1 = 0 
  ss1 = 0 

  if(hr.ge.1.and.hr.le.6) then 
     hr1 = 0
  elseif(hr.ge.7.and.hr.le.12) then 
     hr1 = 6
  elseif(hr.ge.13.and.hr.le.18) then 
     hr1 = 12
  elseif(hr.ge.19.and.hr.le.23) then 
     hr1 = 18
  elseif(hr.eq.0) then 
     hr1 = 0
     ts1 = -6*60*60
     call LIS_tick(time,doy,gmt,yr1,mo1,da1,hr1,mn1,ss1,real(ts1))
  endif
  call LIS_get_julhr(yr1,mo1,da1,hr1,mn1,ss1,julbeg)

end subroutine find_agrfld_starttime

