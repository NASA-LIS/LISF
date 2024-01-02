!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
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
! !ROUTINE: readagrmetforcing
! \label{readagrmetforcing}
!
! !REVISION HISTORY:
! 29Jul2005; Sujay Kumar, Initial Code
! 13Dec2007: Marv Freimund, Do NOT abort for missing 10 files
! 31Dec2007: Marv Freimund, simplify file creation code; use trim on filenames
! 1 Aug2008: Switched processing from PS to a generic cartesian grid with 
!            built-in implicit parallelism
! 10 MAR 2010 Added some variables to pass in calls to interp_agrmet, loadcloud,
!             tr_coeffs.................................Michael Shaw/WXE
!
! !INTERFACE:
subroutine readagrmetforcing(n,findex, order)
! !USES:
  use LIS_coreMod, only         : LIS_rc, LIS_domain, LIS_localPet
  use LIS_timeMgrMod,only       : LIS_get_julhr,LIS_tick,LIS_time2date
  use LIS_LMLCMod,   only       : LIS_LMLC
  use LIS_albedoMod, only       : LIS_alb
  use LIS_vegDataMod, only      : LIS_gfrac
  use LIS_snowMod, only         : LIS_snow_struc
  use LIS_logMod, only          : LIS_logunit
  use AGRMET_forcingMod, only : agrmet_struc

  use LIS_mpiMod
#ifdef ESMF_TRACE
   use ESMF
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: order
!
! !DESCRIPTION:
!  This routine generates the meteorological forcing variables for the
!  current time. The surface observations along with GFS data is
!  processed to generate estimates of temperature, pressure, wind and
!  relative humidity. The WWMCA cloud data is read in to generate the
!  downward shortwave and longwave radiation fields. This routine also
!  performs the spatial interpolation to the LIS grid
!  and filling any gaps introduced due to mismatches in the LIS
!  and AGRMET masks.
!
!  The indices of glbdata1 and glbdata2 correspond to the following
!  variables:
!  \begin{verbatim}
!   1 - 2m air temp
!   2 - 2m relative humidity
!   3 - shortwave
!   4 - longwave
!   5 - u wind
!   6 - v wind
!   7 - surface pressure
!   8 - precip
!  \end{verbatim}
!
!  The arguments and variables are:
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous
!    hourly instance, order=2, read the next hourly instance)
!  \item[n]
!    index of the nest
!  \item[filename]
!    generated filename
!  \item[hemi]
!    index of hemisphere loops
!  \item[c,r,i,t]
!    looping and indexing variables
!  \item[istat]
!    io error stat variable
!  \item[file\_exists]
!    check for archived file
!  \item[varfield]
!    interpolated variable
!  \item[tair]
!    interpolated surface temperature
!  \item[psurf]
!    interpolated surface pressure
!  \item[rlh]
!    interpolated relative humidity
!  \item[vapor\_pres]
!    vapor pressure value
!  \item[ip]
!    choice of interpoation algorithm
!  \item[julhr]
!    julian hour value
!  \item[amounts]
!    retrieved CDFS2 cloud amounts
!  \item[types]
!    retrieved CDFS2 cloud types
!  \item[tops]
!    retrieved CDFS2 cloud tops
!  \item[times]
!    retrieved CDFS2 pixel times
!  \item[cldamt\_nh, cldamt\_sh]
!    retrieved cloud amounts for each hemisphere
!  \item[cldtyp\_nh, cldtyp\_sh]
!    retrieved cloud types for each hemisphere
!  \item[fog\_nh,fog\_sh]
!    fog present flags for each hemisphere
!  \item[thres]
!    cloud top thresholds
!  \item[q2sat]
!    saturation mixing rati at the 1st model level
!  \item[esat]
!    saturation vapor pressure
!  \item[udef]
!    undefined variable used for spatial interpolation
!  \item[doy1,yr1,mo1,da1,hr1,mn1,ss1,try,ts1,backtime1,gmt1]
!    time/date specific variables
!  \item[fog]
!    fog present flag
!  \item[bare]
!    flag to eliminate greenness at a point
!  \item[albedo]
!    snow free albedo
!  \item[snup]
!    snow depth thresholds
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[agrmet\_sfcalc](\ref{AGRMET_sfcalc}) \newline
!    calls the routines to produce an analysis of temperature, pressure
!    winds and relative humidity
!  \item[find\_agrsfc\_dataindex](\ref{find_agrsfc_dataindex}) \newline
!    computes the data index for the current hour
!  \item[interp\_agrmetvar](\ref{interp_agrmetvar}) \newline
!    spatial interpolation of an AGRMET variable to LIS grid
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the input times of different files
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!  \item[agrmet\_cdfs\_type\_filename](\ref{agrmet_cdfs_type_filename}) \newline
!    generates the filename to read CDFS2 cloud types
!  \item[agrmet\_cdfs\_pcts\_filename](\ref{agrmet_cdfs_pcts_filename}) \newline
!    generates the filename to read CDFS2 cloud amounts
!  \item[agrmet\_cdfs\_hgts\_filename](\ref{agrmet_cdfs_hgts_filename}) \newline
!    generates the filename to read CDFS2 cloud tops
!  \item[agrmet\_cdfs\_pixltime\_filename](\ref{agrmet_cdfs_pixltime_filename}) \newline
!    generates the filename to read CDFS2 pixel times
!  \item[get\_julhr] (\ref{LIS_get_julhr}) \newline
!    converts the date to a julian hour
!  \item[AGRMET\_loadcloud](\ref{AGRMET_loadcloud}) \newline
!   converts CDFS2 data for use in radiation routines
!  \item[AGRMET\_svp](\ref{AGRMET_svp}) \newline
!   computes saturation vapor pressure and saturation specific humidity
!  \item[AGRMET\_longwv](\ref{AGRMET_longwv}) \newline
!   computes the downward longwave radiation
!  \item[AGRMET\_calc\_albedo](\ref{AGRMET_calc_albedo}) \newline
!   computes the snow free albedo for current day
!  \item[AGRMET\_tr\_coeffs](\ref{AGRMET_tr_coeffs}) \newline
!   computes the transmissivity and reflectivity coefficients
!  \item[AGRMET\_solar](\ref{AGRMET_solar}) \newline
!   computes the downward shortwave radiation value
!  \end{description}
!EOP
  character*200    :: filename
  integer          :: hemi
  integer          :: c,r,t,gid,sftype
  integer          :: istat
  logical          :: file_exists
  real,allocatable :: longwv(:,:)
  real,allocatable :: coszen_ps(:,:,:)
  real,allocatable :: swdown(:,:)
  real,allocatable :: varfield(:,:)
  real,allocatable :: coszen(:,:)
  real,allocatable :: tair(:,:)
  real             :: vapor_pres
  real,allocatable :: psurf(:,:)
  real,allocatable :: rlh(:,:)
  integer          :: ip
  integer          :: julhr
  real,allocatable             :: cldamt(:,:,:)
  integer,allocatable          :: cldamt_nh( :,:,:)
  integer,allocatable          :: cldtyp_nh( :,:,:)
  integer,allocatable          :: cldamt_sh( :,:,:)
  integer,allocatable          :: cldtyp_sh( :,:,:)
  logical,allocatable          :: fog_nh( :,:)
  logical,allocatable          :: fog_sh( :,:)
  integer          :: thres(5) !needs to be read in.
  real             :: q2sat
  real             :: e,esat
  real             :: udef
  integer          :: doy1,yr1,mo1,da1,hr1,mn1,ss1,try
  real*8           :: backtime1
  real             :: gmt1,ts1
  logical          :: fog, bare
  real,allocatable             :: albedo(:)
  real,allocatable             :: albedo_tmp(:)
  integer,allocatable          :: ncount(:)
  real             :: snup(24)
  integer          :: dindex

  real,allocatable            :: r1_ps(:,:,:)
  real,allocatable            :: r2_ps(:,:,:)
  real,allocatable            :: r3_ps(:,:,:)
  real,allocatable            :: t1_ps(:,:,:)
  real,allocatable            :: t2_ps(:,:,:)
  real,allocatable            :: t3_ps(:,:,:)

  real,allocatable            :: t1(:,:)
  real,allocatable            :: t2(:,:)
  real,allocatable            :: t3(:,:)
  real,allocatable            :: r1(:,:)
  real,allocatable            :: r2(:,:)
  real,allocatable            :: r3(:,:)

  real,allocatable            :: cod(:,:,:)

  data snup /.080, .040, .040, .040, .040, .040, .040, .040, &
             .040, .080, .080, .080, .080, .080, .080, .080, &
             .040, .080, .025, .040, .040, .040, .025, .025/


  data THRES /12600, 12300, 12000, 11700, 11400/

  TRACE_ENTER("agrmet_readforc")
  allocate(longwv(LIS_rc%lnc(n), LIS_rc%lnr(n)))
  allocate(varfield(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  allocate(tair(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  allocate(psurf(LIS_rc%lnc(n), LIS_rc%lnr(n)))
  allocate(rlh(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  allocate(coszen_ps(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
  allocate(coszen(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  allocate(cldamt(3, LIS_rc%lnc(n), LIS_rc%lnr(n)))

  allocate(t1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  allocate(t2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  allocate(t3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  allocate(r1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  allocate(r2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  allocate(r3(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  if ( agrmet_struc(n)%compute_radiation == 'cloud types' ) then
     allocate(cldamt_nh( 3,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
     allocate(cldtyp_nh( 3,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
     allocate(cldamt_sh( 3,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
     allocate(cldtyp_sh( 3,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
     allocate(fog_nh( agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
     allocate(fog_sh( agrmet_struc(n)%imax,agrmet_struc(n)%jmax))

     allocate(r1_ps(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
     allocate(r2_ps(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
     allocate(r3_ps(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
     allocate(t1_ps(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
     allocate(t2_ps(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
     allocate(t3_ps(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
  else
     allocate(cod(3, LIS_rc%lnc(n), LIS_rc%lnr(n)))
  endif

  allocate(albedo(LIS_rc%ngrid(n)))
  allocate(albedo_tmp(LIS_rc%ntiles(n)))
  allocate(ncount(LIS_rc%ngrid(n)))


  call AGRMET_sfcalc(n)

  call find_agrsfc_dataindex(LIS_rc%hr,dindex)


!  open(100,file='tmp.bin',form='unformatted')
!  write(100) agrmet_struc(n)%sfctmp(dindex,:,:)
!  close(100) 
!  print*, 'Finished with SFCALC '
!  stop
  if(LIS_rc%run_model) then

     tair  = agrmet_struc(n)%sfctmp(dindex,:,:)
     psurf = agrmet_struc(n)%sfcprs(dindex,:,:)
     rlh = agrmet_struc(n)%sfcrlh(dindex,:,:)

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              agrmet_struc(n)%metdata2(1,LIS_domain(n)%gindex(c,r)) = &
                   agrmet_struc(n)%sfctmp(dindex,c,r)
              agrmet_struc(n)%metdata2(7,LIS_domain(n)%gindex(c,r)) = &
                   agrmet_struc(n)%sfcprs(dindex,c,r)
              agrmet_struc(n)%metdata2(5,LIS_domain(n)%gindex(c,r)) = &
                   agrmet_struc(n)%sfcspd(dindex,c,r)
              agrmet_struc(n)%metdata2(6,LIS_domain(n)%gindex(c,r)) = 0.0

           endif
        enddo
     end do

!     print*,'EMK: maxval(agrmet_struc(n)%sfcrlh(dindex,:,:) = ', &
!          maxval(agrmet_struc(n)%sfcrlh(dindex,:,:))
!     print*,'EMK: minval(agrmet_struc(n)%sfcrlh(dindex,:,:) = ', &
!          minval(agrmet_struc(n)%sfcrlh(dindex,:,:))

     do c =1, LIS_rc%lnc(n)
        do r = 1,LIS_rc%lnr(n)
           if (LIS_domain(n)%gindex(c,r).ne. -1) then
              vapor_pres = 611*exp((17.3*(tair(c,r)-273.16))/&
                   (tair(c,r)-35.86))*agrmet_struc(n)%sfcrlh(dindex,c,r)
              varfield(c,r) =  0.622*vapor_pres/&
                   (psurf(c,r)-0.378*vapor_pres)
           endif
        enddo
     enddo

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              agrmet_struc(n)%metdata2(2,LIS_domain(n)%gindex(c,r)) = varfield(c,r)
           endif
        enddo
     enddo

     if ( agrmet_struc(n)%compute_radiation == 'cloud types' ) then
        call compute_type_based_clouds(n, cldamt_nh, cldamt_sh, cldamt, &
                                       cldtyp_nh, cldtyp_sh, fog_nh, fog_sh)
     else
        call compute_cod_based_clouds(n, cod, cldamt)
     endif

     yr1 = LIS_rc%yr
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = LIS_rc%hr
     mn1 = LIS_rc%mn
     ss1 = LIS_rc%ss

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           ! EMK...Calculate for all valid grid indices
!           if(LIS_LMLC(n)%landmask(c,r).eq.1) then 
           if ( LIS_domain(n)%gindex(c,r) /= -1 ) then
              call AGRMET_svp(q2sat,esat,&
                   psurf(c,r), tair(c,r))
              e = esat*rlh(c,r)
              call AGRMET_longwv(tair(c,r),e,cldamt(:,c,r),longwv(c,r))
           else
              longwv(c,r) = LIS_rc%udef
           endif
        enddo
     enddo

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              agrmet_struc(n)%metdata2(4,LIS_domain(n)%gindex(c,r)) = longwv(c,r)
           endif
        enddo
     enddo

     fog = .false.
     if(agrmet_struc(n)%findtime1.eq.1.and.agrmet_struc(n)%findtime2.eq.1) then
#if 0 
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(c,r).ne.-1) then
                 t = LIS_domain(n)%gindex(c,r)
! The snup value is looked up from a table based on the vegtype.
! This is a little tricky when subgrid tiling is used since each tile
! in a gridcell has a different vegtype. For now proceeding with a
! simple 1st tile lookup.
!
                 if  (LIS_domain(n)%tile(t)%vegt .eq. LIS_rc%bareclass  .or. &
                      LIS_domain(n)%tile(t)%vegt .eq. LIS_rc%urbanclass .or. &
                      LIS_domain(n)%tile(t)%vegt .eq. LIS_rc%snowclass)  then
                    bare = .true.
                 else
                    bare = .false.
                 endif

                 call AGRMET_calc_albedo(LIS_alb(n)%albsf(t),albedo(t), &
                      LIS_gfrac(n)%greenness(t), LIS_alb(n)%mxsnalb(t), &
                      LIS_snow_struc(n)%sneqv(t),&
                      snup(LIS_domain(n)%tile(t)%vegt),&
                      agrmet_struc(n)%salp, bare)
              endif
           enddo
        enddo
        LIS_alb(n)%albedo(:) = albedo(:)
#endif        
        do t=1,LIS_rc%ntiles(n)
           sftype = LIS_domain(n)%tile(t)%sftype
           if(sftype.eq.1) then 
              if  (LIS_domain(n)%tile(t)%vegt .eq. LIS_rc%bareclass  .or. &
                   LIS_domain(n)%tile(t)%vegt .eq. LIS_rc%urbanclass .or. &
                   LIS_domain(n)%tile(t)%vegt .eq. LIS_rc%snowclass)  then
                 bare = .true.
              else
                 bare = .false.
              endif

              call AGRMET_calc_albedo(LIS_alb(n)%albsf(t),albedo_tmp(t), &
                   LIS_gfrac(n)%greenness(t), LIS_alb(n)%mxsnalb(t), &
                   LIS_snow_struc(n)%sneqv(t),&
                   snup(LIS_domain(n)%tile(t)%vegt),&
                   agrmet_struc(n)%salp, bare)
           else
!ocean - set a fixed albedo
!for small zenith angle 0.03-0.10
! for large zenith angle 0.10 -1.0
              albedo_tmp(t) = 0.10
           endif
        enddo        
        LIS_alb(n)%albedo(:) = albedo_tmp(:)
     else
        albedo_tmp(:)   = LIS_alb(n)%albedo(:)
     endif

     albedo = 0 
     ncount = 0 
     do t=1,LIS_rc%ntiles(n)
        c= LIS_domain(n)%tile(t)%col
        r= LIS_domain(n)%tile(t)%row
        gid = LIS_domain(n)%gindex(c,r)
        if(gid.ne.-1) then 
           albedo(gid) = albedo(gid) + albedo_tmp(t)
           ncount(gid) = ncount(gid) + 1
        endif
     enddo

     do t=1,LIS_rc%ngrid(n)
        if(ncount(t).gt.0) then 
           albedo(t) = albedo(t)/ncount(t)
        endif
     enddo


     if ( agrmet_struc(n)%compute_radiation == 'cloud types' ) then
        r1_ps = 0.0
        r2_ps = 0.0
        r3_ps = 0.0

        t1_ps = 0.0
        t2_ps = 0.0
        t3_ps = 0.0

        coszen_ps = -1.0

        do hemi=1,2
           do c=1,agrmet_struc(n)%imax
              do r=1,agrmet_struc(n)%jmax
                 if(agrmet_struc(n)%land(c,r,hemi).gt.0) then
                    if(hemi.eq.1) then
                       call AGRMET_tr_coeffs(c,r,hemi,&
                            fog_nh(c,r), cldtyp_nh(:,c,r),&
                            cldamt_nh(:,c,r), &
                            r1_ps(hemi,c,r),r2_ps(hemi,c,r),r3_ps(hemi,c,r),&
                            t1_ps(hemi,c,r),t2_ps(hemi,c,r),t3_ps(hemi,c,r),&
                            coszen_ps(hemi,c,r),yr1,mo1,da1,hr1,n)
                    else
                       call AGRMET_tr_coeffs(c,r,hemi,&
                            fog_sh(c,r), cldtyp_sh(:,c,r),&
                            cldamt_sh(:,c,r), &
                            r1_ps(hemi,c,r),r2_ps(hemi,c,r),r3_ps(hemi,c,r),&
                            t1_ps(hemi,c,r),t2_ps(hemi,c,r),t3_ps(hemi,c,r),&
                            coszen_ps(hemi,c,r),yr1,mo1,da1,hr1,n)
                    endif
                 end if
              enddo
           enddo
        enddo

        ip = 1
        call interp_agrmetvar(n,ip,r1_ps,0.0,r1,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
        call interp_agrmetvar(n,ip,r2_ps,0.0,r2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
        call interp_agrmetvar(n,ip,r3_ps,0.0,r3,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)

        call interp_agrmetvar(n,ip,t1_ps,0.0,t1,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
        call interp_agrmetvar(n,ip,t2_ps,0.0,t2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
        call interp_agrmetvar(n,ip,t3_ps,0.0,t3,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)

        call interp_agrmetvar(n,ip,coszen_ps,-1.0,coszen,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
     else
        call compute_tr_from_cod(n, t1, r1, cod(1,:,:))
        call compute_tr_from_cod(n, t2, r2, cod(2,:,:))
        call compute_tr_from_cod(n, t3, r3, cod(3,:,:))

        coszen_ps = -1.0
        do hemi=1,2
           do c=1,agrmet_struc(n)%imax
              do r=1,agrmet_struc(n)%jmax
                 if(agrmet_struc(n)%land(c,r,hemi).gt.0) then
                    call AGRMET_coszen(yr1, mo1, da1, hr1, hemi, c, r, &
                                       int(agrmet_struc(n)%imax/64),   &
                                       coszen_ps(hemi,c,r))
                 endif
              enddo
           enddo
        enddo

        ip = 1
        call interp_agrmetvar(n,ip,coszen_ps,-1.0,coszen,&
                              agrmet_struc(n)%imax,agrmet_struc(n)%jmax)

        ! See AGRMET_tr_coeffs
        where ( coszen < 0.01 )
           t1 = 0.0; t2 = 0.0; t3 = 0.0
           r1 = 0.0; r2 = 0.0; r3 = 0.0
        endwhere

        ! When any ti equals 0, then swdown will be 0.  See AGRMET_solar.
        ! However, when t1 = 0 = t2 or when t2 = 0 = t3, then backscatter
        ! becomes 0, which leads to a divide-by-zero error.  Thus wherever
        ! any level (low, middle, high) of transmittance is 0, set all
        ! levels of transmittance and reflectance to 0.  This still yields
        ! swdown = 0, but avoids any divide-by-zero errors.
        where ( t1 == 0.0 .or. t2 == 0.0 .or. t3 == 0.0 )
           t1 = 0.0; t2 = 0.0; t3 = 0.0
           r1 = 0.0; r2 = 0.0; r3 = 0.0
        endwhere
     endif

     allocate(swdown(LIS_rc%lnc(n),LIS_rc%lnr(n)))

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if ( LIS_domain(n)%gindex(c,r) /= -1 ) then
              t = LIS_domain(n)%gindex(c,r)
              if(r1(c,r).eq.LIS_rc%udef) r1(c,r) = 0.0
              if(r2(c,r).eq.LIS_rc%udef) r2(c,r) = 0.0
              if(r3(c,r).eq.LIS_rc%udef) r3(c,r) = 0.0
              if(t1(c,r).eq.LIS_rc%udef) t1(c,r) = 0.0
              if(t2(c,r).eq.LIS_rc%udef) t2(c,r) = 0.0
              if(t3(c,r).eq.LIS_rc%udef) t3(c,r) = 0.0

              call AGRMET_solar(coszen(c,r),&
                   albedo(t),r1(c,r),r2(c,r),r3(c,r),&
                   t1(c,r),t2(c,r),t3(c,r),&
                   swdown(c,r),yr1,mo1,da1,hr1)
           else
              swdown(c,r) = LIS_rc%udef
           endif
        enddo
     enddo

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              agrmet_struc(n)%metdata2(3,LIS_domain(n)%gindex(c,r)) = swdown(c,r)
           endif
        enddo
     enddo

     deallocate(swdown)
     deallocate(longwv)
     deallocate(varfield)
     deallocate(tair)
     deallocate(psurf)
     deallocate(rlh)

     deallocate(coszen_ps)
     deallocate(coszen)

     deallocate(cldamt)

     deallocate(t1)
     deallocate(t2)
     deallocate(t3)
     deallocate(r1)
     deallocate(r2)
     deallocate(r3)

     if ( agrmet_struc(n)%compute_radiation == 'cloud types' ) then
        deallocate(cldamt_nh)
        deallocate(cldtyp_nh)
        deallocate(cldamt_sh)
        deallocate(cldtyp_sh)
        deallocate(fog_nh)
        deallocate(fog_sh)

        deallocate(r1_ps)
        deallocate(r2_ps)
        deallocate(r3_ps)
        deallocate(t1_ps)
        deallocate(t2_ps)
        deallocate(t3_ps)
     else
        deallocate(cod)
     endif

     deallocate(albedo)
     deallocate(albedo_tmp)
     deallocate(ncount)
  end if
  TRACE_EXIT("agrmet_readforc")
end subroutine readagrmetforcing
!BOP
!
! !ROUTINE: agrmet_cdfs_pcts_filename
! \label{agrmet_cdfs_pcts_filename}
!
! !INTERFACE:
subroutine agrmet_cdfs_pcts_filename(fname,rootdir,dir,&
                                     use_timestamp,hemi,yr,mo,da,hr)

  implicit none
! !ARGUMENTS:
  character(*)        :: fname
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: hemi
  integer, intent(in) :: yr,mo,da,hr
!
! !DESCRIPTION:
!  This routines generates the name of the CDFS2 cloud amounts file
!  by appending the hemisphere and timestamps to the root directory.
!
!  The arguments are:
!  \begin{description}
!   \item[fname]
!    created filename
!   \item[rootdir]
!    path to the root directory containing the data
!   \item[dir]
!    name of the subdirectory containing the data
!   \item[use\_timestamp]
!    flag to indicate whether the directories
!    should be timestamped or not
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[yr]
!    4 digit year
!   \item[mo]
!    integer value of month (1-12)
!   \item[da]
!    day of the month
!   \item[hr]
!    hour of the day
!  \end{description}
!EOP

  character( 3), parameter :: fhemi(2) = (/'NH_','SH_'/)
  character(10) :: ftime1, ftime2

  write (UNIT=ftime2, FMT='(i4, i2.2, i2.2, i2.2)') yr,mo,da,hr

  if ( use_timestamp == 1 ) then
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/',yr,mo,da,'/'
     fname = trim(rootdir)//ftime1//trim(dir)//'/WWMCA_LYRD_CLOUD_PCTS_'// &
             fhemi(hemi) // ftime2
  else
     fname = trim(rootdir)//trim(dir)//'/WWMCA_LYRD_CLOUD_PCTS_'// &
             fhemi(hemi) // ftime2
  endif

end subroutine agrmet_cdfs_pcts_filename

!BOP
!
! !ROUTINE: agrmet_cdfs_type_filename
! \label{agrmet_cdfs_type_filename}
!
! !INTERFACE:
subroutine agrmet_cdfs_type_filename(fname,rootdir,dir,&
                                     use_timestamp,hemi,yr,mo,da,hr)

  implicit none
! !ARGUMENTS:
  character(*)        :: fname
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: hemi
  integer, intent(in) :: yr,mo,da,hr
!
! !DESCRIPTION:
!  This routines generates the name of the CDFS2 cloud types file
!  the hemisphere and timestamps to the root directory.
!
!  The arguments are:
!  \begin{description}
!   \item[fname]
!    created filename
!   \item[rootdir]
!    path to the root directory containing the data
!   \item[dir]
!    name of the subdirectory containing the data
!   \item[use\_timestamp]
!    flag to indicate whether the directories
!    should be timestamped or not
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[yr]
!    4 digit year
!   \item[mo]
!    integer value of month (1-12)
!   \item[da]
!    day of the month
!   \item[hr]
!    hour of the day
!  \end{description}
!EOP

  character( 3), parameter :: fhemi(2) = (/'NH_','SH_'/)
  character(10) :: ftime1, ftime2

  write (UNIT=ftime2, fmt='(i4, i2.2, i2.2, i2.2)') yr,mo,da,hr

  if ( use_timestamp == 1 ) then
     write (UNIT=ftime1, fmt='(a1, i4, i2.2, i2.2, a1)') '/',yr,mo,da,'/'
     fname = trim(rootdir)//ftime1//trim(dir)//'/WWMCA_LYRD_CLOUD_TYPE_'//&
             fhemi(hemi) // ftime2
  else
     fname = trim(rootdir)//trim(dir)//'/WWMCA_LYRD_CLOUD_TYPE_'// &
             fhemi(hemi) // ftime2
  endif

end subroutine agrmet_cdfs_type_filename

!BOP
!
! !ROUTINE: agrmet_cdfs_hgts_filename
! \label{agrmet_cdfs_hgts_filename}
!
! !INTERFACE:
subroutine agrmet_cdfs_hgts_filename(fname,rootdir,dir,&
                                     use_timestamp,hemi,yr,mo,da,hr)

  implicit none
! !ARGUMENTS:
  character(*)        :: fname
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: hemi
  integer, intent(in) :: yr,mo,da,hr
!
! !DESCRIPTION:
!  This routines generates the name of the CDFS2 cloud tops file
!  with hemisphere and timestamps to the root directory.
!
!  The arguments are:
!  \begin{description}
!   \item[fname]
!    created filename
!   \item[rootdir]
!    path to the root directory containing the data
!   \item[dir]
!    name of the subdirectory containing the data
!   \item[use\_timestamp]
!    flag to indicate whether the directories
!    should be timestamped or not
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[yr]
!    4 digit year
!   \item[mo]
!    integer value of month (1-12)
!   \item[da]
!    day of the month
!   \item[hr]
!    hour of the day
!  \end{description}
!EOP

  character( 3), parameter :: fhemi(2) = (/'NH_','SH_'/)
  character(len=10) :: ftime1, ftime2

  write (UNIT=ftime2, fmt='(i4, i2.2, i2.2, i2.2)') yr,mo,da,hr

  if (use_timestamp == 1) then
     write (UNIT=ftime1, fmt='(a1, i4, i2.2, i2.2, a1)') '/',yr,mo,da,'/'
     fname = trim(rootdir)//ftime1//trim(dir)//'/WWMCA_LYRD_TOPS_HGHTS_'// &
           fhemi(hemi) // ftime2
  else
     fname = trim(rootdir)//trim(dir)//'/WWMCA_LYRD_TOPS_HGHTS_'// &
             fhemi(hemi) // ftime2
  endif

end subroutine agrmet_cdfs_hgts_filename

!BOP
!
! !ROUTINE: agrmet_cdfs_pixltime_filename
! \label{agrmet_cdfs_pixltime_filename}
!
! !INTERFACE:
subroutine agrmet_cdfs_pixltime_filename(fname,rootdir,dir,&
                                         use_timestamp,hemi,yr,mo,da,hr)

  implicit none
! !ARGUMENTS:
  character(*)        :: fname
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: hemi
  integer, intent(in) :: yr,mo,da,hr
!
! !DESCRIPTION:
!  This routines generates the name of the CDFS2 pixel times file.
!  the hemisphere and timestamps to the root directory.
!
!  The arguments are:
!  \begin{description}
!   \item[fname]
!    created filename
!   \item[rootdir]
!    path to the root directory containing the data
!   \item[dir]
!    name of the subdirectory containing the data
!   \item[use\_timestamp]
!    flag to indicate whether the directories
!    should be timestamped or not
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[yr]
!    4 digit year
!   \item[mo]
!    integer value of month (1-12)
!   \item[da]
!    day of the month
!   \item[hr]
!    hour of the day
!  \end{description}
!EOP

  character( 3), parameter :: fhemi(2) = (/'NH_','SH_'/)
  character(10) :: ftime1, ftime2

  write (UNIT=ftime2, FMT='(i4, i2.2, i2.2, i2.2)') yr,mo,da,hr

  if (use_timestamp .eq. 1 ) then
  write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/',yr,mo,da,'/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/WWMCA_PIXL_TIME_MEANS_' // &
             fhemi(hemi) // ftime2
  else
     fname = trim(rootdir) //  '/'  // trim(dir) // '/WWMCA_PIXL_TIME_MEANS_' // &
             fhemi(hemi) // ftime2
  endif

end subroutine agrmet_cdfs_pixltime_filename

!BOP
!
! !ROUTINE: find_agrsfc_dataindex
! \label{find_agrsfc_dataindex}
!
! !INTERFACE:
subroutine find_agrsfc_dataindex(hr,dindex)

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: hr
  integer, intent(inout) :: dindex
!
! !DESCRIPTION:
!  This routine finds index to be used in the SFCALC
!  data for the current instance. The SFCALC processing
!  is done every 6 hours, which produces products
!  for the upcoming hourly intervals. For every hour within
!  a 6 hour window, to retrieve the data, the program
!  will simply use the index returned by this routine,
!  instead of redoing the processing.
!
!  The arguments are:
!  \begin{description}
!  \item[hr]
!    the current hour
!  \item[dindex]
!    output data index
!  \end{description}
!
!EOP
  if(hr.ge.1.and.hr.le.6) then
     dindex = hr
  elseif(hr.ge.7.and.hr.le.12) then
     dindex = hr-6
  elseif(hr.ge.13.and.hr.le.18) then
     dindex = hr-12
  elseif(hr.ge.19.and.hr.le.23) then
     dindex = hr-18
  elseif(hr.eq.0) then
     dindex = 6
  endif

end subroutine find_agrsfc_dataindex

subroutine getforcfile(fname,var,hemi,yr,mo,da,hr)

  implicit none

  character(*)       :: fname
  integer,intent(in) :: var
  integer,intent(in) :: hemi, yr,mo,da,hr

  character( 4), parameter :: fhemi(2) = (/'_NH_','_SH_'/)
  character( 3), parameter :: fname2(4) = (/'tmp','rlh','prs','wnd'/)
  character(10) :: ftime2

  write (UNIT=ftime2, FMT='(i4, i2.2, i2.2, i2.2)') yr,mo,da,hr

  fname = 'lis_' // fname2(var) // fhemi(hemi) // ftime2

end subroutine getforcfile
