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

! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
!BOP
!
! !ROUTINE: readagrmetpcpforcing
! \label{readagrmetpcpforcing}
!
! !REVISION HISTORY:
! 29Jul2005; Sujay Kumar, Initial Code
! 27Nov2007; Yudong Tian, Now it is called by getagrmet.F90 at every time step 
!    in real-time mode, instead of every hour before. 
! 27Dec2007: Marv Freimund, use trim on filenames
! 1 Aug2008: Switched processing from PS to a generic cartesian grid with 
!            built-in implicit parallelism
! 14Aug2008: Chris Franks, added call to readblacklist
! 21Apr2009: Sujay Kumar, Removed the intermediate I/O routines
! 29Jun2009: Chris Franks, increase spreading radii for 1/4 degree
! 16Nov2010: Chris Franks, remove call to AGRMET_readblacklist
! 19Jun2012: Ryan Ruhge, Add support for CMORPH
! 20Oct2017: Eric Kemp, switched to Bratseth scheme
! 04Jan2018: Eric Kemp, disable 6-hr cycle unless 12-hr is not run.
! 05Sep2018: Eric Kemp, code clean-up and addition of IMERG.
! 07Jan2020: Eric Kemp, added creation of precip_OBA directory.
!
! !INTERFACE:    
subroutine readagrmetpcpforcing(n,findex, order)

   ! !USES:
   use AGRMET_forcingMod, only : agrmet_struc
#ifdef ESMF_TRACE
   use ESMF, only: ESMF_TraceRegionEnter, ESMF_TraceRegionExit
#endif
   use LIS_coreMod, only         : LIS_rc, LIS_domain, LIS_masterproc
   use LIS_logMod, only          : LIS_logunit
   use LIS_mpiMod, only          : LIS_mpi_comm
   use LIS_pluginIndices, only   : LIS_agrmetrunId
   use LIS_timeMgrMod, only      : LIS_julhr_date, LIS_get_julhr
   use USAF_bratsethMod, only: USAF_obsData, USAF_createObsData, &
        USAF_getbacknwp, &
        USAF_split6hrgaugeobsdata,USAF_destroyobsdata,USAF_getcmorphobsdata, &
        USAF_interpbacktotypeobsdata, USAF_getssmiobsdata, &
        USAF_getgeoprecipobsdata, USAF_setbratsethprecipstats, &
        USAF_analyzeprecip, USAF_split12hrgaugeobsdata, USAF_analyzeprecip, &
        USAF_countGoodObs,USAF_waterQC,USAF_filterObsData,USAF_dupQC, &
        USAF_snowQC,USAF_snowDepthQC,USAF_backQC,USAF_superstatQC
   !use USAF_ImergHHMod, only: newImergHHPrecip, destroyImergHHPrecip, &
   !     update30minImergHHPrecip,copyToObsDataImergHHPrecip, &
   !     calc3hrImergHHPrecip, count3hrObsImergHHPrecip, &
   !     create_Imerg_HH_filename, fetch3hrImergHH
   use USAF_ImergHHMod, only: newImergHHPrecip, destroyImergHHPrecip, &
        update30minImergHHPrecip,copyToObsDataImergHHPrecip, &
        count3hrObsImergHHPrecip, &
        create_Imerg_HH_filename, fetch3hrImergHH
   use USAF_OBAMod, only: OBA,destroyOBA,makeFilename,writeToFile

   implicit none

   ! !ARGUMENTS: 
   integer, intent(in) :: n
   integer, intent(in) :: findex
   integer, intent(in) :: order
!
! !DESCRIPTION:
!  This routine calls the AGRMET precipitation analysis for the 
!  current time. The analysis includes the use of observed
!  precipitation (rain gauge reports) from AFWA's global surface
!  observation database, estimates using methods based on 
!  remotely-sensed and climatological data. The analysis also
!  blends real and estimated precip amounts using a barnes technique. 
!  A hierarchy of these estimates is established to ensure
!  the most reliable data is used. 
!
!  The arguments and variables are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!  \item[n]
!    index of the nest
!  \item[c,r,i,t]
!    looping and indexing variables
!  \item[use\_twelve]
!   flag to use the 12-hourly precip amounts or the 6 hourly 
!   amounts 
!  \item[julbeg]
!    starting julian hour
!  \item[julend]
!    ending julian hour
!  \item[j3hr]
!    3 hourly julian time
!  \item[p6]
!    6 hourly rain-gauge precip amounts (mm) on the AGRMET grid
!  \item[estpcp]
!    array of precip estimates (last two 3-hourly periods)
!  \item[estp6]
!    estimated 6hr precip array (mm/6hr)
!  \item[source]
!    array that contains source-of-the-estimate information 
!    with specified source flags. 
!  \item[cdfs2est]
!    CDFS2 cloud based precip estimate (mm/3hr)
!  \item[prcpwe]
!    present/past weather estimate on the grid
!  \item[p12]
!    12 hourly precip amounts (mm) on the AGRMET grid
!  \item[mrgp6]
!    6 hourly merged precip amounts (mm)
!  \item[relp]
!    parsed 3 hour real precip array (mm/3hr)
!  \item[relp6]
!    parsed 6hr real precip array (mm/6hr)
!  \item[cdfsii6]
!    6hr CDFS2 based precip estimate (mm/6hr)
!  \item[srcwts]
!    value weight placed on each source
!  \item[addrad]
!    added radius value
!  \item[maxrad]
!    maximum radius value
!  \item[minrad]
!    minimum radius value
!  \item[varrad]
!    radius change amount due to the real precip variance
!  \item[tmp]
!    temporary array to hold weighted values
!  \item[wgt]
!    weighted array point values
!  \item[curr\_time]
!    current time 
!  \item[pcp\_process]
!    flag to determine if the precip analysis should be performed
!    in the current hour
!  \item[varfield]
!    interpolated variable
!  \item[yr1,mo1,da1,hr1,mn1,ss1]
!    time/date specific variables 
!  \item[alert\_number]
!    alert message number
!  \item[ip]
!    interpolation option
!  \item[c,r,k]
!   looping and indexing variables
!  \item[pathpcp]
!   Path to agrmet precip files
!  \item[cmorphpixel]
!   Indicates if pixels are fed with cmorph data
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[find\_agrpcp\_starttime](\ref{find_agrpcp_starttime}) \newline
!    computes the start time for the precip processing.
!  \item[julhr\_date] (\ref{LIS_julhr_date}) \newline
!    converts julian hour to a date format
!  \item[AGRMET\_getpcpobs](\ref{AGRMET_getpcpobs}) \newline
!    read AFWA precip surface observations
!  \item[AGRMET\_makest](\ref{AGRMET_makest}) \newline
!    create an consolidated precip amount based on 
!    several possible precipitation estimates
!  \item[AGRMET\_phsrel12](\ref{AGRMET_phsrel12}) \newline
!    retrieve phase 12 hourly rain gauge amounts
!  \item[AGRMET\_phsrel6](\ref{AGRMET_phsrel6}) \newline
!    retrieve phase 6 hourly raing gauge amounts
!  \item[AGRMET\_valid](\ref{AGRMET_valid}) \newline
!    validate precip amounts
!  \item[AGRMET\_pcp\_barnes](\ref{AGRMET_pcp_barnes}) \newline
!    perform barnes analysis on the precip data
!  \item[AGRMET\_make03](\ref{AGRMET_make03}) \newline
!    make 3 hourly merged precipitation
!  \item[interp\_agrmetvar](\ref{interp_agrmetvar}) \newline
!    spatial interpolation of an AGRMET variable to LIS grid
!  \item[AGRMET\_fillgaps](\ref{AGRMET_fillgaps}) \newline
!    fills the gaps in the interpolated field due to mismatches in 
!    LIS and AGRMET masks
!  \end{description}
!EOP
   logical        :: use_twelve 
   integer        :: julhr
   integer        :: julbeg
   integer        :: julend_lis
   integer        :: j3hr
   real           :: p6(LIS_rc%lnc(n), LIS_rc%lnr(n))   
   real           :: estpcp(LIS_rc%lnc(n), LIS_rc%lnr(n), 4)
   real           :: estp6(LIS_rc%lnc(n),LIS_rc%lnr(n))
   integer        :: source(LIS_rc%lnc(n), LIS_rc%lnr(n), 4)
   logical        :: cmorphpixel(LIS_rc%lnc(n), LIS_rc%lnr(n), 4)
   real           :: cdfs2est(LIS_rc%lnc(n), LIS_rc%lnr(n), 4)
   real           :: prcpwe(LIS_rc%lnc(n), LIS_rc%lnr(n), 4)   
   real           :: p12(LIS_rc%lnc(n), LIS_rc%lnr(n))  
   real           :: mrgp6(LIS_rc%lnc(n), LIS_rc%lnr(n))  
   real           :: relp(LIS_rc%lnc(n), LIS_rc%lnr(n), 4)   
   real           :: relp6(LIS_rc%lnc(n), LIS_rc%lnr(n))  
   real           :: cdfsii6(LIS_rc%lnc(n), LIS_rc%lnr(n))  
   real           :: srcwts(8)
   integer        :: addrad(8)
   integer        :: maxrad(6)
   integer        :: minrad(6)
   integer        :: varrad(4)
   integer        :: k
   real           :: curr_time
   real           :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))
   integer        :: c,r
   integer        :: yr1,mo1,da1,hr1
   integer        :: alert_number
   character*100       :: pathpcp
   character*10        :: date10_03
   character*4         :: fyr
   character*2         :: fmo,fda
   character*100       :: ofil
   character*100       :: ifil
   logical             :: exists
   integer             :: ftn
   real                :: gridDesc(6)
   integer             :: ierr
   integer             :: icount

   ! EMK NEW
   integer, external :: LIS_create_subdirs
   integer :: rc
   real, allocatable :: back4(:,:,:)
   integer :: nobs_good
   type(USAF_obsData) :: precip_6hr_gage_tmp, precip_6hr_gage
   type(USAF_obsData) :: precip_12hr_gage_tmp, precip_12hr_gage
   type(USAF_obsData) :: precip_3hrly_gage_tmp(4), precip_3hrly_gage(4)
   type(USAF_obsData) :: precip_3hrly_ssmi_tmp(4), precip_3hrly_ssmi(4)
   type(USAF_obsData) :: precip_3hrly_geop_tmp(4), precip_3hrly_geop(4)
   type(USAF_obsData) :: precip_3hrly_cmorph_tmp(4), precip_3hrly_cmorph(4)
   type(USAF_obsData) :: precip_3hrly_imerg_tmp(4), precip_3hrly_imerg(4)
   type(USAF_obsData) :: precip_3hrly_all
   integer :: julbeg_adj, ii, nobs_good_extra

   !type(USAF_obsData) :: precip06gage, precip012gage
   !type(USAF_obsData) :: precip06gage_tmp, precip012gage_tmp
   !type(USAF_obsData) :: precip3gage, precip6gage, precip9gage, precip12gage
   !type(USAF_obsData) :: precip3ssmi, precip6ssmi, precip9ssmi, precip12ssmi
   !type(USAF_obsData) :: precip3ssmi_tmp, precip6ssmi_tmp, precip9ssmi_tmp, &
   !     precip12ssmi_tmp
   !type(USAF_obsData) :: precip3geo, precip6geo, precip9geo, precip12geo
   !type(USAF_obsData) :: precip3geo_tmp, precip6geo_tmp, precip9geo_tmp, &
   !     precip12geo_tmp
   !type(USAF_obsData) :: precip3cmorph, precip6cmorph, precip9cmorph, &
   !     precip12cmorph
   !type(USAF_obsData) :: precip3cmorph_tmp, precip6cmorph_tmp, &
   !     precip9cmorph_tmp, &
   !     precip12cmorph_tmp
   !type(USAF_obsData) :: precip3Imerg, precip6Imerg, precip9Imerg, &
   !     precip12Imerg
   !type(USAF_obsData) :: precip3Imerg_tmp, precip6Imerg_tmp, precip9Imerg_tmp, &
   !     precip12Imerg_tmp
   integer :: hourindex
   integer :: k1,k2,k3,k4
   character(len=32) :: type
   character(len=10) :: yyyymmddhh
   character(len=50) :: pathOBA
   logical :: found_inq
   character(len=120) :: obaFilename
   type(OBA) :: precipOBA ! EMK
   character(len=6) :: pcp_src(4)
   character(len=255) :: imerg_datadir
   character(len=20) :: imerg_product
   character(len=20) :: imerg_version
   integer*2 :: imerg_plp_thresh
   real :: imerg_sigmaOSqr
   real :: imerg_oErrScaleLength
   character(len=32) :: imerg_net, imerg_platform
   real :: sigmaBSqr
   character(len=32) :: new_name
   integer :: use_expanded_station_ids
   data alert_number / 0 / 
   data srcwts /100.0,50.0,4.0,4.0,1.0,1.0,60.0,1.0/
   data addrad /0, -1, -2, -4, -5, -5, 0, -2/
   ! data addrad /0, -1, -2, -4, -5, -5, -4, -6/
   ! data maxrad /15, 12, 12, 6, 3, 3/
   data maxrad /15, 1, 12, 1, 1, 1/
   data minrad /4, 1, 2, 1, 1, 1/
   data varrad /0, -2, -4, -6/
   
   pathOBA = "./precip_OBA" ! EMK 

   ! See if subdirectory exists
   if (agrmet_struc(n)%oba_switch .eq. 1 .or. &
        agrmet_struc(n)%oba_switch .eq. 2) then
      if (LIS_masterproc) then
         inquire(file=trim(pathOBA), exist=found_inq)
         if (.not. found_inq) then
            rc = lis_create_subdirs(len_trim(pathOBA), trim(pathOBA))
            if (rc .ne. 0) then
               write(LIS_logunit, *) &
                    '[WARN] Cannot create directory ', trim(pathOBA)
            end if
         end if
      end if
   end if

   ! YDT 9/26/07 save merged precip in analysis dir, instead of forcing dir
   !  to avoid overwriting the forcing archive. 
   !  pathpcp = trim(agrmet_struc(n)%agrmetdir)
   pathpcp = trim(agrmet_struc(n)%analysisdir)
   curr_time = float(LIS_rc%hr)*60+float(LIS_rc%mn)+float(LIS_rc%ss/60)
   
   ! return immediatly if it is not at 
   !  3hr + ts, i.e. 0:15, 3:15, 6:15, 9:15, 12:15, 15:15, 18:15 and 21:15
   ! which are the time to read in precp data for the upcoming 3 hrs, and it 
   ! is not a fresh start (pcp_ready = false for fresh start). 

   if ( mod(curr_time, 180.0) .ne. LIS_rc%ts/60.0 .and. &
        agrmet_struc(n)%pcp_ready ) return 

   TRACE_ENTER("agrmet_readpcpforc")

   p6 = LIS_rc%udef
   p12 = LIS_rc%udef
   prcpwe = LIS_rc%udef

   k1 = 1
   k2 = 2
   k3 = 3
   k4 = 4
   pcp_src(:) = 'NULL'
   
   !determines which of the two time segments we are in: (0, 12], or (12, 24].
   ! need to get 3hr ahead because precip at 3Z is used for run in (0Z, 3Z]. 
   ! Note this subroutine is only called hourly  or at 3hr + ts
   ! when LIS is on [12, 15, 18, 21, 24), it is time to prepare seg II. 
   if (curr_time .GE. 0 .and. curr_time .LT. 12*60 ) then 
      call LIS_get_julhr(LIS_rc%yr,LIS_rc%mo,LIS_rc%da,&
           0, LIS_rc%mn,LIS_rc%ss,julbeg)
   else  !   prepare seg II.  
      call LIS_get_julhr(LIS_rc%yr,LIS_rc%mo,LIS_rc%da,&
           12, LIS_rc%mn,LIS_rc%ss,julbeg)
   end if
   call LIS_get_julhr(LIS_rc%eyr,LIS_rc%emo,LIS_rc%eda,&
        LIS_rc%ehr, LIS_rc%emn,LIS_rc%ess,julend_lis)

   ! when start from fresh (in LIS_init) or first step into a new segment, 
   !  precip is not ready. Need to do pcp analysis 
   if (curr_time .EQ. LIS_rc%ts/60 .or. curr_time .EQ. 12*60.0+LIS_rc%ts/60 ) &
        agrmet_struc(n)%pcp_ready = .false. 
  
   call AGRMET_julhr_date10(julbeg, date10_03)
   write(LIS_logunit,*)'[INFO] Entering pcp proc. pcp_ready=', &
        agrmet_struc(n)%pcp_ready, date10_03

   !====================================================================== 
   if ( .not. agrmet_struc(n)%pcp_ready ) then  ! do pcp analysis for the 
                                                ! whole seg 

      ! EMK...Need entire background field to interpolate to all observations
      allocate(back4(LIS_rc%gnc(n),LIS_rc%gnr(n),4))
      back4 = LIS_rc%udef ! EMK
      sigmaBSqr = agrmet_struc(n)%bratseth_precip_back_sigma_b_sqr
      
      !EMK...We need to decide if 6-hr or 12-hr cycle should be run.  NOT BOTH,
      !since the 12-hr cycle values overwrite the 6-hr ones.
      use_twelve = .false.

      ! EMK...The runmode check is probably for legacy testing.  But I will
      ! leave it alone for now.
      if (( (julend_lis - julbeg) .ge. 12) .or. &
           (LIS_rc%runmode .ne. LIS_agrmetrunId) ) then
         use_twelve = .true.
      end if

      ! EMK...6-hr cycle
      if (.not. use_twelve) then

         write(LIS_logunit,*)'[INFO] PCP Processing 6-hr cycle starts .... ', &
              date10_03
         !first comes the 6-hour cycle ---------------------------
         use_twelve = .false.

         ! Get the NWP background field
         call USAF_getBackNWP(n,back4,pcp_src,use_twelve,julbeg, findex)

         ! Handle rain gage data.
         if(agrmet_struc(n)%pcpobswch.eq.1) then

            ! First fetch 6hr accumulations,  Only 6hr values are fetched
            ! since use_twelve is false.
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            write(LIS_logunit,*) &
                 '[INFO] Fetching 6-hr gage data'
            call USAF_createObsData(precip_6hr_gage_tmp,n,maxobs=15000)
            !call AGRMET_getpcpobs(n, julbeg, LIS_rc%mo, prcpwe, &
            !     use_twelve, p6, p12, alert_number, precip_6hr_gage_tmp,  &
            !     precip_12hr_gage_tmp, &
            !     pcp_src)
            use_expanded_station_ids = agrmet_struc(n)%pcpobsfmt - 1
            call USAF_getpcpobs(n, julbeg, LIS_rc%mo, use_twelve, &
                 pcp_src, use_expanded_station_ids, alert_number, &
                 precip_6hr_gage_tmp, precip_12hr_gage_tmp)

            ! Reject data over water
            write(LIS_logunit,*) &
                 '[INFO] Running waterQC on 6 hr gauge observations'
            call USAF_waterQC(precip_6hr_gage_tmp,n)
            nobs_good = USAF_countGoodObs(precip_6hr_gage_tmp)
            nobs_good_extra = nint(nobs_good*1.10)
            call USAF_createObsData(precip_6hr_gage,n,maxobs=nobs_good_extra)
            call USAF_filterObsData(precip_6hr_gage,precip_6hr_gage_tmp)
            call USAF_destroyObsData(precip_6hr_gage_tmp)

            ! Split the 6hr accumulations into 3hr totals.  This also
            ! interpolates background field to the gages.
            call USAF_split6hrGaugeObsData(precip_6hr_gage,n, &
                 LIS_rc%gnc(n), LIS_rc%gnr(n), back4, agrmet_struc(n)%pcap, &
                 precip_3hrly_gage_tmp(1), precip_3hrly_gage_tmp(2)) ! EMK
            call USAF_destroyObsData(precip_6hr_gage)

            ! Now work with 3hrly values.  Run additional QC.
            do ii = 1,2
               ! Handle duplicate 3-hr reports.
               write(LIS_logunit,*) &
                    '[INFO] Running dupQC on 3hrly gage obs, set ',ii
               call USAF_dupQC(precip_3hrly_gage_tmp(ii))

               hourindex = calc_hourindex(ii)

               ! Reject apparent snowfall.
               write(LIS_logunit,*) &
                    '[INFO] Running snowQC on 3hrly gage obs, set ',ii
               call USAF_snowQC(precip_3hrly_gage_tmp(ii),n,hourindex, &
                    threshold=275.)

               ! Compare with background field
               if (agrmet_struc(n)%skip_backqc .ne. 1) then
                  write(LIS_logunit,*) &
                       '[INFO] Running backQC on 3hrly gage obs, set ',ii
                  call USAF_backQC(precip_3hrly_gage_tmp(ii),sigmaBSqr)
               end if

               ! Create "superobservations" from close gauges
               if (agrmet_struc(n)%skip_superstatqc .ne. 1) then
                  new_name = "SUPERGAGE"
                  write(LIS_logunit,*) &
                       '[INFO] Running superstatQC on 3hrly gage obs, set ',ii
                  call USAF_superstatQC(precip_3hrly_gage_tmp(ii),n,new_name)
                  type = new_name
                  call USAF_interpBackToTypeObsData(precip_3hrly_gage_tmp(ii),&
                       n, &
                       LIS_rc%gnc(n),LIS_rc%gnr(n),back4(:,:,ii),type)
               end if

               ! Finally, filter out the bad obs
               nobs_good = USAF_countGoodObs(precip_3hrly_gage_tmp(ii))
               call USAF_createObsData(precip_3hrly_gage(ii),n,&
                    maxobs=nobs_good)
               call USAF_filterObsData(precip_3hrly_gage(ii), &
                    precip_3hrly_gage_tmp(ii))
               call USAF_destroyObsData(precip_3hrly_gage_tmp(ii))
            end do           
         else
            ! Create dummy objects with no obs.
            do ii = 1,2
               call USAF_createObsData(precip_3hrly_gage(ii),n,maxobs=1)
            end do
            p6 = LIS_rc%udef
            p12 = LIS_rc%udef
            prcpwe = LIS_rc%udef
         endif         
         call AGRMET_julhr_date10(julbeg, date10_03)
#if (defined SPMD)
         call MPI_Barrier(LIS_mpi_comm,ierr)       
#endif

         ! Handle SSMI
         if (agrmet_struc(n)%raswch .eq. 1) then
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            write(LIS_logunit,*) &
                 '[INFO] Fetching 6 hours of SSMI data'
            do ii = 1, 2
               nobs_good_extra = nint(1024*1024*2*1.10)
               call USAF_createObsData(precip_3hrly_ssmi_tmp(ii),n,&
                    maxobs=nobs_good_extra)
            end do
            ! Only 0-3hr and 3-6hr values fetched since use_twelve is false
            call USAF_getSSMIObsData(n,julbeg,use_twelve, &
                 precip_3hrly_ssmi_tmp(1),precip_3hrly_ssmi_tmp(2), &
                 precip_3hrly_ssmi_tmp(3),precip_3hrly_ssmi_tmp(4),pcp_src)

            do ii = 1,2
               ! Reject obs over water
               write(LIS_logunit,*) &
                    '[INFO] Running waterQC on 3hrly SSMI obs, set ',&
                    ii
               call USAF_waterQC(precip_3hrly_ssmi_tmp(ii),n, &
                    silent_rejects=.true.)

               ! Interpolate background field to remaining obs
               type = "SSMI"
               call USAF_interpBackToTypeObsData(precip_3hrly_ssmi_tmp(ii),n, &
                    LIS_rc%gnc(n), &
                    LIS_rc%gnr(n),back4(:,:,ii),type)

               hourindex = calc_hourindex(ii)

               ! Reject apparent snowfall.
               write(LIS_logunit,*) &
                    '[INFO] Running snowQC on 3hrly SSMI obs, set ',ii
               call USAF_snowQC(precip_3hrly_ssmi_tmp(ii),n,hourindex, &
                    threshold=278., &
                    silent_rejects=.true.)

               ! Reject satellite estimates over snow
               write(LIS_logunit,*) &
                    '[INFO] Running snowDepthQC on 3hrly SSMI obs, set ',ii
               call USAF_snowDepthQC(precip_3hrly_ssmi_tmp(ii),n, &
                    silent_rejects=.true.)

               ! Compare with background field
               !if (agrmet_struc(n)%skip_backqc .ne. 1) then
               !   write(LIS_logunit,*) &
               !        '[INFO] Running backQC on 3hrly SSMI obs, set ',ii
               !   call USAF_backQC(precip_3hrly_ssmi_tmp(ii),sigmaBSqr, &
               !        silent_rejects=.true.)
               !end if

               ! Create "superobservations" from close SSMI retrievals
               !if (agrmet_struc(n)%skip_superstatqc .ne. 1) then
               !   new_name = "SUPERSSMI"
               !   write(LIS_logunit,*) &
               !     '[INFO] Running superstatQC on 3hrly SSMI obs, set ',ii
               !   call USAF_superstatQC(precip_3hrly_SSMI_tmp(ii),n,new_name, &
               !     silent_rejects=.true.)
               !   type = new_name
               !   call USAF_interpBackToTypeObsData(precip_3hrly_ssmi_tmp(ii),&
               !        n, &
               !        LIS_rc%gnc(n),LIS_rc%gnr(n),back4(:,:,ii),type)
               !end if

               ! Filter out the bad obs
               nobs_good = USAF_countGoodObs(precip_3hrly_ssmi_tmp(ii))
               call USAF_createObsData(precip_3hrly_ssmi(ii),n, &
                    maxobs=nobs_good)
               call USAF_filterObsData(precip_3hrly_ssmi(ii), &
                    precip_3hrly_ssmi_tmp(ii))
               call USAF_destroyObsData(precip_3hrly_ssmi_tmp(ii))

            end do

         else
            ! Create dummy objects with no obs.
            do ii = 1,2
               call USAF_createObsData(precip_3hrly_ssmi(ii),n,maxobs=1)
            end do
         end if

         ! Handle GEOPRECIP
         if (agrmet_struc(n)%geoswch .eq. 1 .or. &
              agrmet_struc(n)%geoswch .eq. 2) then
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            write(LIS_logunit,*) &
                 '[INFO] Fetching 6 hours of GEOPRECIP data'
            do ii = 1,2
               nobs_good_extra = nint(1024*1024*2*1.10)
               call USAF_createObsData(precip_3hrly_geop_tmp(ii),n, &
                    maxobs=nobs_good_extra)
            end do
            call USAF_getGeoPrecipObsData(n,julbeg,use_twelve, &
                 precip_3hrly_geop_tmp(1), &
                 precip_3hrly_geop_tmp(2), &
                 precip_3hrly_geop_tmp(3), &
                 precip_3hrly_geop_tmp(4),pcp_src)

            do ii = 1,2
               ! Reject obs over water
               write(LIS_logunit,*) &
                    '[INFO] Running waterQC on 3 hrly GEOPRECIP obs, set ',ii
               call USAF_waterQC(precip_3hrly_geop_tmp(ii),n,&
                    silent_rejects=.true.)

               ! Interpolate background field
               type = "GEOPRECIP"
               call USAF_interpBackToTypeObsData(precip_3hrly_geop_tmp(ii),n, &
                    LIS_rc%gnc(n), &
                    LIS_rc%gnr(n),back4(:,:,ii),type)

               hourindex = calc_hourindex(ii)

               ! Reject apparent snowfall.
               write(LIS_logunit,*) &
                    '[INFO] Running snowQC on 3hrly GEOPRECIP obs, set ',ii
               call USAF_snowQC(precip_3hrly_geop_tmp(ii),n,hourindex, &
                    threshold=real(agrmet_struc(n)%geomaxthresh), &
                    silent_rejects=.true.)

               ! Reject satellite estimates over snow
               write(LIS_logunit,*) &
                    '[INFO] Running snowDepthQC on 3hrly GEOPRECIP obs, set ',&
                    ii
               call USAF_snowDepthQC(precip_3hrly_geop_tmp(ii),n, &
                    silent_rejects=.true.)

               ! Compare with background field
               !if (agrmet_struc(n)%skip_backqc .ne. 1) then
               !   write(LIS_logunit,*) &
               !        '[INFO] Running backQC on 3hrly GEOPRECIP obs, set ',ii
               !   call USAF_backQC(precip_3hrly_geop_tmp(ii),sigmaBSqr, &
               !        silent_rejects=.true.)
               !end if

               ! Create "superobservations" from close obs
               ! if (agrmet_struc(n)%skip_superstatqc .ne. 1) then
               !    new_name = "SUPERGEO"
               !    !write(LIS_logunit,*) &
               !    !     '[INFO] Running superstatQC on 3hrly GEOPRECIP obs, set ',&
               !    !     ii
               !    !call USAF_superstatQC(precip_3hrly_geop_tmp(ii),n,new_name,&
               !    !     silent_rejects=.true.)
               !    type = new_name
               !    call USAF_interpBackToTypeObsData(precip_3hrly_geop_tmp(ii),&
               !         n, &
               !         LIS_rc%gnc(n),LIS_rc%gnr(n),back4(:,:,ii),type)
               ! end if

               ! Filter the bad obs
               nobs_good = USAF_countGoodObs(precip_3hrly_geop_tmp(ii))
               call USAF_createObsData(precip_3hrly_geop(ii),n, &
                    maxobs=nobs_good)
               call USAF_filterObsData(precip_3hrly_geop(ii), &
                    precip_3hrly_geop_tmp(ii))
               call USAF_destroyObsData(precip_3hrly_geop_tmp(ii))

            end do
         else
            ! Create dummy objects with no obs.
            do ii = 1,2
               call USAF_createObsData(precip_3hrly_geop(ii),n,maxobs=1)
            end do
         end if

         ! Handle CMORPH
         if (agrmet_struc(n)%cmorswch .eq. 1) then
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            write(LIS_logunit,*) &
                 '[INFO] Fetching 6 hours of CMORPH data'
            do ii = 1,2
               nobs_good_extra = nint(4948*1649*1.10)
               call USAF_createObsData(precip_3hrly_cmorph_tmp(ii),n, &
                    maxobs=nobs_good_extra)
            end do
            ! Only 0-3hr and 3-6hr values fetched since use_twelve is false
            call USAF_getCMORPHObsData(n,julbeg,use_twelve, &
                 precip_3hrly_cmorph_tmp(1), &
                 precip_3hrly_cmorph_tmp(2), &
                 precip_3hrly_cmorph_tmp(3), &
                 precip_3hrly_cmorph_tmp(4),pcp_src)

            do ii = 1,2
               ! Reject obs over water
               write(LIS_logunit,*) &
                    '[INFO] Running waterQC on 3hrly CMORPH obs, set ',ii
               call USAF_waterQC(precip_3hrly_cmorph_tmp(ii),n, &
                    silent_rejects=.true.)

               ! Interpolate the background field
               type = "CMORPH"
               call USAF_interpBackToTypeObsData(precip_3hrly_cmorph_tmp(ii),&
                    n,&
                    LIS_rc%gnc(n), &
                    LIS_rc%gnr(n),back4(:,:,ii),type)

               hourindex = calc_hourindex(ii)

               ! Reject apparent snowfall.
               write(LIS_logunit,*) &
                    '[INFO] Running snowQC on 3hrly CMORPH obs, set ',ii
               call USAF_snowQC(precip_3hrly_cmorph_tmp(ii),n,hourindex, &
                    threshold=real(agrmet_struc(n)%cmormaxthresh), &
                    silent_rejects=.true.)

               ! Reject satellite estimates over snow
               write(LIS_logunit,*) &
                    '[INFO] Running snowDepthQC on 3hrly CMORPH obs, set ',&
                    ii
               call USAF_snowDepthQC(precip_3hrly_cmorph_tmp(ii),n, &
                    silent_rejects=.true.)

               ! ! Compare with background field
               ! if (agrmet_struc(n)%skip_backqc .ne. 1) then
               !    write(LIS_logunit,*) &
               !         '[INFO] Running backQC on 3hrly CMORPH obs, set ',ii
               !    call USAF_backQC(precip_3hrly_cmorph_tmp(ii),sigmaBSqr, &
               !         silent_rejects=.true.)
               ! end if

               ! Create "superobservations" from close CMORPH retrievals
               ! if (agrmet_struc(n)%skip_superstatqc .ne. 1) then
               !    new_name = "SUPERCMRPH"
               !    !write(LIS_logunit,*) &
               !    !     '[INFO] Running superstatQC on 3hrly CMORPH obs, set ',ii
               !    !call USAF_superstatQC(precip_3hrly_cmorph_tmp(ii),n,new_name, &
               !         silent_rejects=.true.)
               !    type = new_name
               !    call USAF_interpBackToTypeObsData(precip_3hrly_cmorph_tmp(ii),&
               !         n, &
               !         LIS_rc%gnc(n),LIS_rc%gnr(n),back4(:,:,ii),type)
               ! end if

               ! Filter out the bad obs
               nobs_good = USAF_countGoodObs(precip_3hrly_cmorph_tmp(ii))
               call USAF_createObsData(precip_3hrly_cmorph(ii),n, &
                    maxobs=nobs_good)
               call USAF_filterObsData(precip_3hrly_cmorph(ii), &
                    precip_3hrly_cmorph_tmp(ii))
               call USAF_destroyObsData(precip_3hrly_cmorph_tmp(ii))

            end do
         else
            ! Create dummy objects with no obs.
            do ii = 1,2
               call USAF_createObsData(precip_3hrly_cmorph(ii),n,maxobs=1)
            end do
         end if

         ! Handle IMERG
         if (agrmet_struc(n)%imerg_swch .eq. 1) then
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            imerg_datadir = trim(agrmet_struc(n)%imerg_dir)
            imerg_product = trim(agrmet_struc(n)%imerg_product)
            imerg_version = trim(agrmet_struc(n)%imerg_version)
            imerg_plp_thresh = agrmet_struc(n)%imerg_plp_thresh
            imerg_sigmaOSqr = &
                 agrmet_struc(n)%bratseth_precip_imerg_sigma_o_sqr
            imerg_oErrScaleLength = &
                 agrmet_struc(n)%bratseth_precip_imerg_err_scale_length
            imerg_net = "IMERG"
            imerg_platform = "IMERG"            
            do ii = 1,2
               write(LIS_logunit,*) &
                    '[INFO] Fetching 6 hours of IMERG data, set ',ii
               julbeg_adj = (ii-1)*3
               call fetch3hrImergHH(julbeg+julbeg_adj, &
                    imerg_datadir,imerg_product, &
                    imerg_version,imerg_plp_thresh,n,imerg_sigmaOSqr,&
                    imerg_oErrScaleLength,imerg_net,imerg_platform,&
                    precip_3hrly_imerg_tmp(ii))

               ! Reject observations over water
               write(LIS_logunit,*) &
                    '[INFO] Running waterQC on 3hrly IMERG obs, set ',ii
               call USAF_waterQC(precip_3hrly_imerg_tmp(ii),n, &
                    silent_rejects=.true.)

               ! Interpolate background field
               type = "IMERG"
               call USAF_interpBackToTypeObsData(precip_3hrly_imerg_tmp(ii),&
                    n, &
                    LIS_rc%gnc(n), &
                    LIS_rc%gnr(n),back4(:,:,ii),type)

               hourindex = calc_hourindex(ii)

               ! Reject apparent snowfall.
               write(LIS_logunit,*) &
                    '[INFO] Running snowQC on 3hrly IMERG obs, set ',ii
               call USAF_snowQC(precip_3hrly_imerg_tmp(ii),n,hourindex, &
                    threshold=real(agrmet_struc(n)%imerg_t_thresh), &
                    silent_rejects=.true.)

               ! Reject satellite estimates over snow
               write(LIS_logunit,*) &
                    '[INFO] Running snowDepthQC on 3hrly IMERG obs, set ',&
                    ii
               call USAF_snowDepthQC(precip_3hrly_imerg_tmp(ii),n, &
                    silent_rejects=.true.)

               ! ! Compare with background field
               ! if (agrmet_struc(n)%skip_backqc .ne. 1) then
               !    write(LIS_logunit,*) &
               !         '[INFO] Running backQC on 3hrly IMERG obs, set ',ii
               !    call USAF_backQC(precip_3hrly_imerg_tmp(ii),sigmaBSqr, &
               !         silent_rejects=.true.)
               ! end if

               ! Create "superobservations" from close IMERG retrievals
               ! if (agrmet_struc(n)%skip_superstatqc .ne. 1) then
               !    new_name = "SUPERIMERG"
               !    !write(LIS_logunit,*) &
               !    !     '[INFO] Running superstatQC on 3hrly IMERG obs, set ',ii
               !    !call USAF_superstatQC(precip_3hrly_imerg_tmp(ii),n,new_name, &
               !    !     silent_rejects=.true.)
               !    type = new_name
               !    call USAF_interpBackToTypeObsData(precip_3hrly_imerg_tmp(ii),&
               !         n, &
               !         LIS_rc%gnc(n),LIS_rc%gnr(n),back4(:,:,ii),type)
               ! end if

               ! Filter out bad obs
               nobs_good = USAF_countGoodObs(precip_3hrly_imerg_tmp(ii))
               call USAF_createObsData(precip_3hrly_imerg(ii),n, &
                    maxobs=nobs_good)
               call USAF_filterObsData(precip_3hrly_imerg(ii), &
                    precip_3hrly_imerg_tmp(ii))
               call USAF_destroyObsData(precip_3hrly_imerg_tmp(ii))

            end do
         else
            ! Create dummy objects with no obs.
            do ii = 1,2
               call USAF_createObsData(precip_3hrly_imerg(ii),n,maxobs=1)
            end do
         end if

#if (defined SPMD)
         call MPI_Barrier(LIS_mpi_comm,ierr)       
#endif

         estp6 = LIS_rc%udef
         relp6   = LIS_rc%udef
         cdfsii6 = LIS_rc%udef
         mrgp6   = LIS_rc%udef
         k = 1
         icount = 0
         THREE_H : Do  j3hr = julbeg + 3, julbeg + 6, 3 

            agrmet_struc(n)%mrgp(:,:,k) = LIS_rc%udef
                        
            hourindex = calc_hourindex(k)

            call AGRMET_julhr_date10(j3hr, yyyymmddhh)       
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            write(LIS_logunit,*) &
                 '[INFO] RUNNING BRATSETH PRECIPITATION ANALYSIS FOR ', &
                 yyyymmddhh

            ! Set Bratseth error statistics based on source of background 
            ! field.
            call USAF_setBratsethPrecipStats(pcp_src(k),n)

            icount = icount + 1
            
            ! Consolidate all good obs into one structure.  This also
            ! destroys the source structures.
            call create_precipAll(n, &
                 precip_3hrly_gage(k), &
                 precip_3hrly_ssmi(k), &
                 precip_3hrly_geop(k), &
                 precip_3hrly_cmorph(k), &
                 precip_3hrly_imerg(k), &
                 precip_3hrly_all)

            call USAF_analyzePrecip(precip_3hrly_all, n, &
                 back4(:,:,k),hourindex, &
                 agrmet_struc(n)%mrgp(:,:,k), precipOBA)
            
            !-----------------------------------------------------------
            !   write 3-hrly precip amounts to file.  accumulate 3-hrly
            !   merged precip amounts into 6-hrly merged precip amounts.
            !-----------------------------------------------------------
            call AGRMET_make03(n, agrmet_struc(n)%mrgp(:,:,k), mrgp6)
            
            ! EMK:  Write OBA data to file
            if (agrmet_struc(n)%oba_switch .eq. 1 .or. &
                 agrmet_struc(n)%oba_switch .eq. 2) then
               if (LIS_masterproc) then
                  call AGRMET_julhr_date10(j3hr, yyyymmddhh)
                  call makeFilename(pathOBA,yyyymmddhh,6,obaFilename)
                  call writeToFile(precipOBA,obaFilename)
               end if
               call destroyOBA(precipOBA)
            end if

#if (defined SPMD)
            call MPI_Barrier(LIS_mpi_comm,ierr)       
#endif
            k = k+1
         enddo THREE_H

#if (defined SPMD)
         call MPI_Barrier(LIS_mpi_comm,ierr)       
#endif
       
      ! EMK 12-hr cycle
      else

         !then comes the 12-hour cycle ---------------------------
         write(LIS_logunit,*) &
              '[INFO] PCP Processing 12-hr cycle starts .... ', date10_03
         use_twelve = .true.

         ! Get NWP field.  Call twice to get 4 time slices instead of two.
         call USAF_getBackNWP(n,back4,pcp_src,.false., julbeg, findex)
         call USAF_getBackNWP(n,back4,pcp_src,.true., julbeg+6, findex)

         ! Handle rain gages
         !first get obs and est at (6, 12Z] or (18, 24Z]
         ! julbeg = 0Z or 12Z all the time 
         if(agrmet_struc(n)%pcpobswch.eq.1) then
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            write(LIS_logunit,*) &
                 '[INFO] Fetching 12 hour gage data'
            call USAF_createObsData(precip_6hr_gage_tmp,n,maxobs=15000) 
            call USAF_createObsData(precip_12hr_gage_tmp,n,maxobs=15000)
            ! EMK...Call this twice to ensure we get obs for first two time
            ! levels.
            !call AGRMET_getpcpobs(n, julbeg, LIS_rc%mo, prcpwe, &
            !     .false., p6, p12, alert_number,precip_6hr_gage_tmp, &
            !     precip_12hr_gage_tmp, &
            !     pcp_src)          
            !call AGRMET_getpcpobs(n, julbeg+6, LIS_rc%mo, prcpwe, &
            !     .true., p6, p12, alert_number,precip_6hr_gage_tmp, &
            !     precip_12hr_gage_tmp, &
            !     pcp_src)
            use_expanded_station_ids = agrmet_struc(n)%pcpobsfmt - 1
            call USAF_getpcpobs(n, julbeg, LIS_rc%mo, .false., pcp_src, &
                 use_expanded_station_ids, alert_number, &
                 precip_6hr_gage_tmp, precip_12hr_gage_tmp)
            call USAF_getpcpobs(n, julbeg+6, LIS_rc%mo, .true., pcp_src, &
                 use_expanded_station_ids, alert_number, &
                 precip_6hr_gage_tmp, precip_12hr_gage_tmp)

            ! Not needed at this point since we have the 12hr accum
            call USAF_destroyObsData(precip_6hr_gage_tmp)

            ! Reject and filter out gage reports over water
            write(LIS_logunit,*) &
                 '[INFO] Running waterQC on 12 hr gauge observations'
            call USAF_waterQC(precip_12hr_gage_tmp,n)
            nobs_good = USAF_countGoodObs(precip_12hr_gage_tmp)
            nobs_good_extra = nint(nobs_good*1.10)
            call USAF_createObsData(precip_12hr_gage,n,maxobs=nobs_good_extra)
            call USAF_filterObsData(precip_12hr_gage,precip_12hr_gage_tmp)
            call USAF_destroyObsData(precip_12hr_gage_tmp)

            ! Split the 12hr accumulations into 3hr totals.  This also
            ! interpolates background field to the gages.
            call USAF_split12hrGaugeObsData(precip_12hr_gage,n, &
                 LIS_rc%gnc(n), LIS_rc%gnr(n), back4, agrmet_struc(n)%pcap, &
                 precip_3hrly_gage_tmp(1), &
                 precip_3hrly_gage_tmp(2), &
                 precip_3hrly_gage_tmp(3), &
                 precip_3hrly_gage_tmp(4)) ! EMK
            call USAF_destroyObsData(precip_12hr_gage)
         
            ! Now work on 3hrly data.  Apply additional QC tests.
            do ii = 1,4

               ! Handle duplicate 3-hr reports.
               write(LIS_logunit,*) &
                    '[INFO] Running dupQC on 3hrly gage obs, set ',ii
               call USAF_dupQC(precip_3hrly_gage_tmp(ii))

               hourindex = calc_hourindex(ii)

               ! Reject apparent snowfall.
               write(LIS_logunit,*) &
                    '[INFO] Running snowQC on 3hrly gage obs, set ',ii
               call USAF_snowQC(precip_3hrly_gage_tmp(ii),n,hourindex, &
                    threshold=275.)

               ! Compare with background field
               if (agrmet_struc(n)%skip_backqc .ne. 1) then
                  write(LIS_logunit,*) &
                       '[INFO] Running backQC on 3hrly gage obs, set ',ii
                  call USAF_backQC(precip_3hrly_gage_tmp(ii),sigmaBSqr)
               end if

               ! Create "superobservations" from close gage reports
               if (agrmet_struc(n)%skip_superstatqc .ne. 1) then
                  new_name = "SUPERGAGE"
                  write(LIS_logunit,*) &
                       '[INFO] Running superstatQC on 3hrly gage obs, set ',ii
                  call USAF_superstatQC(precip_3hrly_gage_tmp(ii),n,new_name)
                  type = new_name
                  call USAF_interpBackToTypeObsData(precip_3hrly_gage_tmp(ii),&
                       n, &
                       LIS_rc%gnc(n),LIS_rc%gnr(n),back4(:,:,ii),type)
               end if

               ! Finally, filter out the bad obs
               nobs_good = USAF_countGoodObs(precip_3hrly_gage_tmp(ii))
               call USAF_createObsData(precip_3hrly_gage(ii),n, &
                    maxobs=nobs_good)
               call USAF_filterObsData(precip_3hrly_gage(ii), &
                    precip_3hrly_gage_tmp(ii))
               call USAF_destroyObsData(precip_3hrly_gage_tmp(ii))
            end do           
         else
            ! Create dummy objects with no obs.
            do ii = 1,4
               call USAF_createObsData(precip_3hrly_gage(ii),n,maxobs=1)
            end do
            p6 = LIS_rc%udef
            p12 = LIS_rc%udef
            prcpwe = LIS_rc%udef
         endif       
         call AGRMET_julhr_date10(julbeg + 6, date10_03)

#if (defined SPMD)
         call MPI_Barrier(LIS_mpi_comm,ierr)
#endif

         ! Handle SSMI
         if (agrmet_struc(n)%raswch .eq. 1) then
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            write(LIS_logunit,*) &
                 '[INFO] Fetching 12 hours of SSMI data'
            do ii = 1,4
               nobs_good_extra = nint(2*1024*1024*1.10)            
               call USAF_createObsData(precip_3hrly_ssmi_tmp(ii),n, &
                    maxobs=nobs_good_extra)
            end do
            use_twelve=.false.
            call USAF_getSSMIObsData(n,julbeg,use_twelve,&
                 precip_3hrly_ssmi_tmp(1), &
                 precip_3hrly_ssmi_tmp(2), &
                 precip_3hrly_ssmi_tmp(3), &
                 precip_3hrly_ssmi_tmp(4),pcp_src)
            use_twelve=.true.
            call USAF_getSSMIObsData(n,julbeg+6,use_twelve,&
                 precip_3hrly_ssmi_tmp(1), &
                 precip_3hrly_ssmi_tmp(2), &
                 precip_3hrly_ssmi_tmp(3), &
                 precip_3hrly_ssmi_tmp(4),pcp_src)

            do ii = 1,4
               ! Reject obs over water
               write(LIS_logunit,*) &
                    '[INFO] Running waterQC on 3 hrly SSMI obs, set ',ii
               call USAF_waterQC(precip_3hrly_ssmi_tmp(ii),n, &
                    silent_rejects=.true.)

               ! Interpolate background field
               type = "SSMI"
               call USAF_interpBackToTypeObsData(precip_3hrly_ssmi_tmp(ii),n, &
                    LIS_rc%gnc(n),&
                    LIS_rc%gnr(n),back4(:,:,ii),type)

               hourindex = calc_hourindex(ii)

               ! Reject apparent snowfall.
               write(LIS_logunit,*) &
                    '[INFO] Running snowQC on 3hrly SSMI obs, set ',ii
               call USAF_snowQC(precip_3hrly_ssmi_tmp(ii),n,hourindex, &
                    threshold=278., &
                    silent_rejects=.true.)

               ! Reject satellite estimates over snow
               write(LIS_logunit,*) &
                    '[INFO] Running snowDepthQC on 3hrly SSMI obs, set ',&
                    ii
               call USAF_snowDepthQC(precip_3hrly_ssmi_tmp(ii),n, &
                    silent_rejects=.true.)

               ! ! Compare with background field
               ! if (agrmet_struc(n)%skip_backqc .ne. 1) then
               !    write(LIS_logunit,*) &
               !         '[INFO] Running backQC on 3hrly SSMI obs, set ',ii
               !    call USAF_backQC(precip_3hrly_ssmi_tmp(ii),sigmaBSqr, &
               !         silent_rejects=.true.)
               ! end if

               ! Create "superobservations" from close SSMI reports
               ! if (agrmet_struc(n)%skip_superstatqc .ne. 1) then
               !    new_name = "SUPERSSMI"
               !    !write(LIS_logunit,*) &
               !    !     '[INFO] Running superstatQC on 3hrly SSMI obs, set ',ii
               !    !call USAF_superstatQC(precip_3hrly_ssmi_tmp(ii),n,new_name, &
               !    !    silent_rejects=.true.)
               !    type = new_name
               !    call USAF_interpBackToTypeObsData(precip_3hrly_ssmi_tmp(ii),&
               !         n, &
               !         LIS_rc%gnc(n),LIS_rc%gnr(n),back4(:,:,ii),type)
               ! end if

               ! Filter out the bad obs
               nobs_good = USAF_countGoodObs(precip_3hrly_ssmi_tmp(ii))
               call USAF_createObsData(precip_3hrly_ssmi(ii),n, &
                    maxobs=nobs_good)
               call USAF_filterObsData(precip_3hrly_ssmi(ii), &
                    precip_3hrly_ssmi_tmp(ii))
               call USAF_destroyObsData(precip_3hrly_ssmi_tmp(ii))

            end do
         else
            ! Create dummy objects with no obs.
            do ii = 1,4
               call USAF_createObsData(precip_3hrly_ssmi(ii),n,maxobs=1)
            end do
         end if
        
         ! Handle GEOPRECIP
         if (agrmet_struc(n)%geoswch .eq. 1 .or. &
              agrmet_struc(n)%geoswch .eq. 2) then
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            write(LIS_logunit,*) &
                 '[INFO] Fetching 12 hours of GEOPRECIP data'
            do ii = 1,4
               nobs_good_extra = nint(2*1024*1024*1.10)            
               call USAF_createObsData(precip_3hrly_geop_tmp(ii),n,&
                    maxobs=nobs_good_extra)
            end do
            
            use_twelve=.false.
            call USAF_getGeoPrecipObsData(n,julbeg,use_twelve, &
                 precip_3hrly_geop_tmp(1),&
                 precip_3hrly_geop_tmp(2),&
                 precip_3hrly_geop_tmp(3),&
                 precip_3hrly_geop_tmp(4),pcp_src)
            use_twelve=.true.
            call USAF_getGeoPrecipObsData(n,julbeg+6,use_twelve, &
                 precip_3hrly_geop_tmp(1),&
                 precip_3hrly_geop_tmp(2),&
                 precip_3hrly_geop_tmp(3),&
                 precip_3hrly_geop_tmp(4),pcp_src)

            do ii = 1,4
               ! Reject obs over water
               write(LIS_logunit,*) &
                    '[INFO] Running waterQC on 3hrly GEOPRECIP obs, set ', &
                    ii
               call USAF_waterQC(precip_3hrly_geop_tmp(ii),n, &
                    silent_rejects=.true.)

               ! Interpolate background field
               type = "GEOPRECIP"
               call USAF_interpBackToTypeObsData(precip_3hrly_geop_tmp(ii),n, &
                    LIS_rc%gnc(n),&
                    LIS_rc%gnr(n),back4(:,:,ii),type)

               hourindex = calc_hourindex(ii)

               ! Reject apparent snowfall.
               write(LIS_logunit,*) &
                    '[INFO] Running snowQC on 3hrly GEOPRECIP obs, set ',ii
               call USAF_snowQC(precip_3hrly_geop_tmp(ii),n,hourindex, &
                    threshold=real(agrmet_struc(n)%geomaxthresh), &
                    silent_rejects=.true.)

               ! Reject satellite estimates over snow
               write(LIS_logunit,*) &
                    '[INFO] Running snowDepthQC on 3hrly GEOPRECIP obs, set ',&
                    ii
               call USAF_snowDepthQC(precip_3hrly_geop_tmp(ii),n, &
                    silent_rejects=.true.)

               ! ! Compare with background field
               ! if (agrmet_struc(n)%skip_backqc .ne. 1) then
               !    write(LIS_logunit,*) &
               !         '[INFO] Running backQC on 3hrly GEOPRECIP obs, set ',ii
               !    call USAF_backQC(precip_3hrly_geop_tmp(ii),sigmaBSqr, &
               !         silent_rejects=.true.)
               ! end if

               ! Create "superobservations" from close GEOPRECIP reports
               ! if (agrmet_struc(n)%skip_superstatqc .ne. 1) then
               !    new_name = "SUPERGEO"
               !    !write(LIS_logunit,*) &
               !    !     '[INFO] Running superstatQC on 3hrly GEOPRECIP obs, set ',&
               !    !     ii
               !    !call USAF_superstatQC(precip_3hrly_geop_tmp(ii),n,new_name, &
               !    !     silent_rejects=.true.)
               !    type = new_name
               !    call USAF_interpBackToTypeObsData(precip_3hrly_geop_tmp(ii),&
               !         n, &
               !         LIS_rc%gnc(n),LIS_rc%gnr(n),back4(:,:,ii),type)
               ! end if

               ! Filter out bad obs
               nobs_good = USAF_countGoodObs(precip_3hrly_geop_tmp(ii))
               call USAF_createObsData(precip_3hrly_geop(ii),n, &
                    maxobs=nobs_good)
               call USAF_filterObsData(precip_3hrly_geop(ii), &
                    precip_3hrly_geop_tmp(ii))
               call USAF_destroyObsData(precip_3hrly_geop_tmp(ii))

            end do
         else
            ! Create dummy objects with no obs.
            do ii = 1,4
               call USAF_createObsData(precip_3hrly_geop(ii),n,maxobs=1)
            end do
         end if

         ! Handle CMORPH
         if (agrmet_struc(n)%cmorswch .eq. 1) then
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            write(LIS_logunit,*) &
                 '[INFO] Fetching 12 hours of CMORPH data'
            do ii = 1,4
               nobs_good_extra = nint(4948*1649*1.10)            
               call USAF_createObsData(precip_3hrly_cmorph_tmp(ii),n, &
                    maxobs=nobs_good_extra)
            end do
            use_twelve=.false.
            call USAF_getCMORPHObsData(n,julbeg,use_twelve,&
                 precip_3hrly_cmorph_tmp(1), &
                 precip_3hrly_cmorph_tmp(2), &
                 precip_3hrly_cmorph_tmp(3), &
                 precip_3hrly_cmorph_tmp(4),pcp_src)
            use_twelve=.true.
            call USAF_getCMORPHObsData(n,julbeg+6,use_twelve,&
                 precip_3hrly_cmorph_tmp(1), &
                 precip_3hrly_cmorph_tmp(2), &
                 precip_3hrly_cmorph_tmp(3), &
                 precip_3hrly_cmorph_tmp(4),pcp_src)

            do ii = 1,4
               ! Reject obs over water
               write(LIS_logunit,*) &
                    '[INFO] Running waterQC on 3hrly CMORPH obs, set ',ii
               call USAF_waterQC(precip_3hrly_cmorph_tmp(ii),n, &
                    silent_rejects=.true.)

               ! Interpolate background field
               type = "CMORPH"
               call USAF_interpBackToTypeObsData(precip_3hrly_cmorph_tmp(ii), &
                    n, &
                    LIS_rc%gnc(n),&
                    LIS_rc%gnr(n),back4(:,:,ii),type)

               hourindex = calc_hourindex(ii)

               ! Reject apparent snowfall.
               write(LIS_logunit,*) &
                    '[INFO] Running snowQC on 3hrly CMORPH obs, set ',ii
               call USAF_snowQC(precip_3hrly_cmorph_tmp(ii),n,hourindex, &
                    threshold=real(agrmet_struc(n)%cmormaxthresh), &
                    silent_rejects=.true.)

               ! Reject satellite estimates over snow
               write(LIS_logunit,*) &
                    '[INFO] Running snowDepthQC on 3hrly CMORPH obs, set ',&
                    ii
               call USAF_snowDepthQC(precip_3hrly_cmorph_tmp(ii),n, &
                    silent_rejects=.true.)

               ! ! Compare with background field
               ! if (agrmet_struc(n)%skip_backqc .ne. 1) then
               !    write(LIS_logunit,*) &
               !         '[INFO] Running backQC on 3hrly CMORPH obs, set ',ii
               !    call USAF_backQC(precip_3hrly_cmorph_tmp(ii),sigmaBSqr, &
               !         silent_rejects=.true.)
               ! end if

               ! Create "superobservations" from close CMORPH reports
               ! if (agrmet_struc(n)%skip_superstatqc .ne. 1) then
               !    new_name = "SUPERCMRPH"
               !    !write(LIS_logunit,*) &
               !    !     '[INFO] Running superstatQC on 3hrly CMORPH obs, set ',&
               !    !     ii
               !    !call USAF_superstatQC(precip_3hrly_cmorph_tmp(ii),n,new_name, &
               !    !     silent_rejects=.true.)
               !    type = new_name
               !    call USAF_interpBackToTypeObsData(precip_3hrly_cmorph_tmp(ii),&
               !         n, &
               !         LIS_rc%gnc(n),LIS_rc%gnr(n),back4(:,:,ii),type)
               ! end if

               ! Filter out bad obs
               nobs_good = USAF_countGoodObs(precip_3hrly_cmorph_tmp(ii))
               call USAF_createObsData(precip_3hrly_cmorph(ii),n, &
                    maxobs=nobs_good)
               call USAF_filterObsData(precip_3hrly_cmorph(ii), &
                    precip_3hrly_cmorph_tmp(ii))
               call USAF_destroyObsData(precip_3hrly_cmorph_tmp(ii))

            end do
         else
            ! Create dummy objects with no obs.
            do ii = 1,4
               call USAF_createObsData(precip_3hrly_cmorph(ii),n,maxobs=1)
            end do
         end if

         ! Handle IMERG
         if (agrmet_struc(n)%imerg_swch .eq. 1) then
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            imerg_datadir = trim(agrmet_struc(n)%imerg_dir)
            imerg_product = trim(agrmet_struc(n)%imerg_product)
            imerg_version = trim(agrmet_struc(n)%imerg_version)
            imerg_plp_thresh = agrmet_struc(n)%imerg_plp_thresh
            imerg_sigmaOSqr = &
                 agrmet_struc(n)%bratseth_precip_imerg_sigma_o_sqr
            imerg_oErrScaleLength = &
                 agrmet_struc(n)%bratseth_precip_imerg_err_scale_length
            imerg_net = "IMERG"
            imerg_platform = "IMERG"                  
            do ii = 1,4
               write(LIS_logunit,*) &
                    '[INFO] Fetching 12 hours of IMERG data, set ',ii
               julbeg_adj = (ii-1)*3
               call fetch3hrImergHH(julbeg+julbeg_adj, &
                    imerg_datadir,imerg_product, &
                    imerg_version,imerg_plp_thresh,n,imerg_sigmaOSqr,&
                    imerg_oErrScaleLength,imerg_net,imerg_platform, &
                    precip_3hrly_imerg_tmp(ii))

               ! Reject obs over water
               write(LIS_logunit,*) &
                    '[INFO] Running waterQC on 3hrly IMERG obs, set ', &
                    ii
               call USAF_waterQC(precip_3hrly_imerg_tmp(ii),n,&
                    silent_rejects=.true.)

               ! Interpolate background field
               type = "IMERG"
               call USAF_interpBackToTypeObsData(precip_3hrly_imerg_tmp(ii),&
                    n, &
                    LIS_rc%gnc(n), &
                    LIS_rc%gnr(n),back4(:,:,ii),type)

               hourindex = calc_hourindex(ii)

               ! Reject apparent snowfall.
               write(LIS_logunit,*) &
                    '[INFO] Running snowQC on 3hrly IMERG obs, set ',ii
               call USAF_snowQC(precip_3hrly_imerg_tmp(ii),n,hourindex, &
                    threshold=real(agrmet_struc(n)%imerg_t_thresh), &
                    silent_rejects=.true.)

               ! Reject satellite estimates over snow
               write(LIS_logunit,*) &
                    '[INFO] Running snowDepthQC on 3hrly IMERG obs, set ',&
                    ii
               call USAF_snowDepthQC(precip_3hrly_imerg_tmp(ii),n, &
                    silent_rejects=.true.)

               ! ! Compare with background field
               ! if (agrmet_struc(n)%skip_backqc .ne. 1) then
               !    write(LIS_logunit,*) &
               !         '[INFO] Running backQC on 3hrly IMERG obs, set ',ii
               !    call USAF_backQC(precip_3hrly_imerg_tmp(ii),sigmaBSqr, &
               !         silent_rejects=.true.)
               ! end if

               ! Create "superobservations" from close IMERG reports
               ! if (agrmet_struc(n)%skip_superstatqc .ne. 1) then
               !    new_name = "SUPERIMERG"
               !    !write(LIS_logunit,*) &
               !    !     '[INFO] Running superstatQC on 3hrly IMERG obs, set ',&
               !    !     ii
               !    !call USAF_superstatQC(precip_3hrly_imerg_tmp(ii),n,new_name,&
               !    !     silent_rejects=.true.)
               !    type = new_name
               !    call USAF_interpBackToTypeObsData(precip_3hrly_imerg_tmp(ii),&
               !         n, &
               !         LIS_rc%gnc(n),LIS_rc%gnr(n),back4(:,:,ii),type)
               ! end if

               ! Filter out the bad obs
               nobs_good = USAF_countGoodObs(precip_3hrly_imerg_tmp(ii))
               call USAF_createObsData(precip_3hrly_imerg(ii),n,&
                    maxobs=nobs_good)
               call USAF_filterObsData(precip_3hrly_imerg(ii), &
                    precip_3hrly_imerg_tmp(ii))
               call USAF_destroyObsData(precip_3hrly_imerg_tmp(ii))

            end do
         else
            ! Create dummy objects with no obs.
            do ii = 1,4
               call USAF_createObsData(precip_3hrly_imerg(ii),n,maxobs=1)
            end do
         end if
            
#if (defined SPMD)
         call MPI_Barrier(LIS_mpi_comm,ierr)
#endif

         estp6 = LIS_rc%udef
         relp6   = LIS_rc%udef
         cdfsii6 = LIS_rc%udef
         mrgp6   = LIS_rc%udef
         
         k = 1
         
         icount = 0
         !---------** (3, 6, 9, 12) or (15, 18, 21, 24) 
         SIX_H : Do  j3hr = julbeg + 3, julbeg + 12, 3  
            
            agrmet_struc(n)%mrgp(:,:,k) = LIS_rc%udef
            
            hourindex = calc_hourindex(k)

            call AGRMET_julhr_date10(j3hr, yyyymmddhh)       
            write(LIS_logunit,*) &
                 '------------------------------------------------------------'
            write(LIS_logunit,*) &
                 '[INFO] RUNNING BRATSETH PRECIPITATION ANALYSIS FOR ', &
                 yyyymmddhh

            ! Set Bratseth error statistics based on source of background 
            ! field.
            call USAF_setBratsethPrecipStats(pcp_src(k),n)

            ! Consolidate all good obs into one structure.  This also
            ! destroys the source structures.
            call create_precipAll(n, &
                 precip_3hrly_gage(k), &
                 precip_3hrly_ssmi(k), &
                 precip_3hrly_geop(k), &
                 precip_3hrly_cmorph(k), &
                 precip_3hrly_imerg(k), &
                 precip_3hrly_all)

            call USAF_analyzePrecip(precip_3hrly_all, n, &
                 back4(:,:,k),hourindex, &
                 agrmet_struc(n)%mrgp(:,:,k), precipOBA)
            
            !-----------------------------------------------------------
            !   write 3-hrly precip amounts to file.  accumulate 3-hrly
            !   merged precip amounts into 6-hrly merged precip amounts.
            !-----------------------------------------------------------

            call AGRMET_make03(n, agrmet_struc(n)%mrgp(:,:,k), mrgp6)
            
            ! EMK:  Write OBA data to file
            if (agrmet_struc(n)%oba_switch .eq. 1 .or. &
                 agrmet_struc(n)%oba_switch .eq. 2) then
               if (LIS_masterproc) then
                  call AGRMET_julhr_date10(j3hr, yyyymmddhh)
                  call makeFilename(pathOBA,yyyymmddhh,12,obaFilename)
                  call writeToFile(precipOBA,obaFilename)
               end if
               call destroyOBA(precipOBA)
            end if


#if (defined SPMD)
            call MPI_Barrier(LIS_mpi_comm,ierr)       
#endif

            k = k+1
         enddo SIX_H
      end if

      agrmet_struc(n)%pcp_ready = .true. 
   end if

 
#if (defined SPMD)
   call MPI_Barrier(LIS_mpi_comm,ierr)
#endif

   write(LIS_logunit,*)'[INFO] Finished with Bratseth precip analysis'

   ! Clean up
   if (allocated(back4)) deallocate(back4)

   !====================================================================== 

   !******** now pcp is ready for model run *******************************
  
   if ( mod(curr_time, 180.0).eq.LIS_rc%ts/60.0 ) then 
      agrmet_struc(n)%pcp_start = .false.  !once booted up, no need
     
      ! reading the nearest reading time
      call find_agrpcp_readtime(LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,julhr) 
     
      call LIS_julhr_date(julhr,yr1,mo1,da1,hr1)
      call AGRMET_julhr_date10( julhr, date10_03 )
      if(LIS_rc%hr.ge.0.and.LIS_rc%hr.lt.3) then 
         k = 1
      elseif(LIS_rc%hr.ge.3.and.LIS_rc%hr.lt.6) then 
         k = 2
      elseif(LIS_rc%hr.ge.6.and.LIS_rc%hr.lt.9) then 
         k = 3
      elseif(LIS_rc%hr.ge.9.and.LIS_rc%hr.lt.12) then 
         k = 4
      elseif(LIS_rc%hr.ge.12.and.LIS_rc%hr.lt.15) then 
         k = 1
      elseif(LIS_rc%hr.ge.15.and.LIS_rc%hr.lt.18) then 
         k = 2
      elseif(LIS_rc%hr.ge.18.and.LIS_rc%hr.lt.21) then 
         k = 3
      else
         k = 4
      endif

      !simply index into the rigth data. 
      varfield = agrmet_struc(n)%mrgp(:,:,k)
      
      do c =1, LIS_rc%lnc(n)
         do r = 1,LIS_rc%lnr(n)
            if (LIS_domain(n)%gindex(c,r).ne. -1) then 
               agrmet_struc(n)%metdata2(8,LIS_domain(n)%gindex(c,r)) = &
                    varfield(c,r)

            endif
         end do
      end do
   end if
   

#if (defined SPMD)
   call MPI_Barrier(LIS_mpi_comm,ierr)
#endif
   TRACE_EXIT("agrmet_readpcpforc")
       
contains

   integer function calc_hourindex(k)
      implicit none
      integer,intent(in) :: k
      if ( (k .eq. 1) .or. (k .eq. 3) ) then
         calc_hourindex = 2
      else if ( (k .eq. 2) .or. (k .eq. 4) ) then
         calc_hourindex=5
      endif
   end function calc_hourindex

   ! Internal subroutine to collect all good obs into same data structure.
   subroutine create_precipAll(n, &
        precip_3hrly_gage, &
        precip_3hrly_ssmi, &
        precip_3hrly_geop, &
        precip_3hrly_cmorph, &
        precip_3hrly_imerg, &
        precip_3hrly_all)

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: n
      type(USAF_ObsData),intent(inout) :: precip_3hrly_gage
      type(USAF_ObsData),intent(inout) :: precip_3hrly_ssmi
      type(USAF_ObsData),intent(inout) :: precip_3hrly_geop
      type(USAF_ObsData),intent(inout) :: precip_3hrly_cmorph
      type(USAF_ObsData),intent(inout) :: precip_3hrly_imerg
      type(USAF_ObsData),intent(out)   :: precip_3hrly_all

      ! Local variables
      integer :: good_obs_gage
      integer :: good_obs_ssmi
      integer :: good_obs_geop
      integer :: good_obs_cmorph
      integer :: good_obs_imerg
      integer :: good_obs

      ! Get counts
      good_obs_gage   = USAF_countGoodObs(precip_3hrly_gage)
      good_obs_ssmi   = USAF_countGoodObs(precip_3hrly_ssmi)
      good_obs_geop   = USAF_countGoodObs(precip_3hrly_geop)
      good_obs_cmorph = USAF_countGoodObs(precip_3hrly_cmorph)
      good_obs_imerg  = USAF_countGoodObs(precip_3hrly_imerg)

      good_obs = good_obs_gage + good_obs_ssmi + good_obs_geop + &
           good_obs_cmorph + good_obs_imerg

      write(LIS_logunit,*)'[INFO] Will work with following observation counts:'
      write(LIS_logunit,*)'  Gauges:       ',good_obs_gage
      write(LIS_logunit,*)'  SSMI:         ',good_obs_ssmi
      write(LIS_logunit,*)'  GEOPRECIP:    ',good_obs_geop
      write(LIS_logunit,*)'  CMORPH:       ',good_obs_cmorph
      write(LIS_logunit,*)'  IMERG:        ',good_obs_imerg
      write(LIS_logunit,*)'  Total count:  ',good_obs

      ! Collect the good obs together
      call USAF_createObsData(precip_3hrly_all,n,maxobs=good_obs)
      call USAF_filterObsData(precip_3hrly_all,precip_3hrly_gage)
      call USAF_filterObsData(precip_3hrly_all,precip_3hrly_ssmi)
      call USAF_filterObsData(precip_3hrly_all,precip_3hrly_geop)
      call USAF_filterObsData(precip_3hrly_all,precip_3hrly_cmorph)
      call USAF_filterObsData(precip_3hrly_all,precip_3hrly_imerg)

      ! Clean up
      call USAF_destroyObsData(precip_3hrly_gage)
      call USAF_destroyObsData(precip_3hrly_ssmi)
      call USAF_destroyObsData(precip_3hrly_geop)
      call USAF_destroyObsData(precip_3hrly_cmorph)
      call USAF_destroyObsData(precip_3hrly_imerg)

   end subroutine create_precipAll
end subroutine readagrmetpcpforcing

!BOP
! 
! !ROUTINE: find_agrpcp_starttime
! \label{find_agrpcp_starttime}
! 
! !INTERFACE:
subroutine find_agrpcp_starttime(yr,mo,da,hr,julbeg)        
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
!  \item[LIS\_tick] (\ref{LIS_tick}) \newline
!    computes previous 6 hour time. 
!  \item[get\_julhr] (\ref{LIS_get_julhr}) \newline
!    converts the date to a julian hour
!  \end{description}
!EOP
  
   integer          :: yr1,mo1,da1,hr1,mn1,ss1
   real*8           :: time1
   integer          :: doy
   real             :: gmt,ts1

   yr1 = yr
   mo1 = mo
   da1 = da
   hr1 = 0     
   mn1 = 0 
   ss1 = 0 
   
   if(hr.ge.1.and.hr.le.6) then 
      hr1 =0
   elseif(hr.ge.7.and.hr.le.12) then 
      hr1 = 6
   elseif(hr.ge.13.and.hr.le.18) then
      hr1 = 12
   elseif(hr.ge.19.and.hr.le.23) then 
      hr1 = 18 
   elseif(hr.eq.0) then 
      hr1 = 0 
      ts1 = -6*60*60
      call LIS_tick(time1,doy,gmt,yr1,mo1,da1,hr1,mn1,ss1,ts1)
   endif
   call LIS_get_julhr(yr1,mo1,da1,hr1,mn1,ss1,julbeg)
end subroutine find_agrpcp_starttime

!BOP
! 
! !ROUTINE: find_agrpcp_readtime
! \label{find_agrpcp_readtime}
! 
! !INTERFACE:
subroutine find_agrpcp_readtime(yr,mo,da,hr,julhr)        
! !USES: 
   use LIS_timeMgrMod, only : LIS_tick, LIS_get_julhr

   implicit none
! !ARGUMENTS: 
   integer, intent(in) :: yr
   integer, intent(in) :: mo
   integer, intent(in) :: da
   integer, intent(in) :: hr
   integer, intent(inout) :: julhr
! 
! !DESCRIPTION: 
!  This routine finds the julian hour to read the AGRMET
!  precip from, based on the current input time. 
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
!    output AGRMET precip reading time
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick] (\ref{LIS_tick}) \newline
!    computes previous 6 hour time. 
!  \item[get\_julhr] (\ref{LIS_get_julhr}) \newline
!    converts the date to a julian hour
!  \end{description}
!EOP
   
   integer          :: yr1,mo1,da1,hr1,mn1,ss1
   real*8           :: time1
   integer          :: doy
   real             :: gmt,ts1
   
   yr1 = yr
   mo1 = mo
   da1 = da
   mn1 = 0 
   ss1 = 0 
   
   if(hr.ge.0.and.hr.lt.3) then 
      hr1 = 3
   elseif(hr.ge.3.and.hr.lt.6) then 
      hr1 = 6
   elseif(hr.ge.6.and.hr.lt.9) then 
      hr1 = 9
   elseif(hr.ge.9.and.hr.lt.12) then 
      hr1 = 12
   elseif(hr.ge.12.and.hr.lt.15) then 
      hr1 = 15
   elseif(hr.ge.15.and.hr.lt.18) then 
      hr1 = 18
   elseif(hr.ge.18.and.hr.lt.21) then 
      hr1 = 21
   elseif(hr.ge.21.and.hr.le.23) then 
      hr1 = 0
! YDT 10/4/07
!     ts1 = 3*60*60
      ts1 = 24*60*60
      call LIS_tick(time1,doy,gmt,yr1,mo1,da1,hr1,mn1,ss1,ts1)
   endif
   call LIS_get_julhr(yr1,mo1,da1,hr1,mn1,ss1,julhr)

end subroutine find_agrpcp_readtime

