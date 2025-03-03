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
! !ROUTINE: get_gfs
!  \label{get_gfs}
!
! !REVISION HISTORY:
!  16 Mar 2008: Sujay Kumar; Initial specification
!  29 Apr 2010: Sujay Kumar: Fixed the problems in mixing the instantaneous
!               and time averaged GFS forecast files. 
!
! !INTERFACE:
subroutine get_gfs(n,findex)
! !USES:
  use LDT_coreMod,       only : LDT_rc
  use LDT_metforcingMod, only : LDT_forc
  use LDT_timeMgrMod,    only : LDT_tick, LDT_get_nstep
  use LDT_logMod,        only : LDT_logunit, LDT_endrun
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use gfs_forcingMod,    only : gfs_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, GFS forcing. 
!  The idea is to open either the 00, 03 and 06 forecast file associated with
!  the most recent GFS assimilation (available every 6 hours).
!  Precipitation rates and radiation fluxes will be taken from the F03 and 
!  F06 files, since averages are provided. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 3 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The strategy for missing data is to first look for 9hr forecasts and if not
!  go backwards up to 10 days to get forcing at the same time of day.
!  Note that when run in forecast mode (for future days), the F00 files
!  are not available and the code defaults to using F03 and F06 files.  
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    determines the GFS data times
!  \item[create\_gfsfilename](\ref{create_gfsfilename}) \newline
!    Puts together appropriate file name for 3 hour intervals
!  \item[create\_gfs\_f0backup\_filename](\ref{create_gfs_f0backup_filename}) \newline
!    Puts together the backup filenames when F00 forecasts are not available. 
!  \item[read_gfs](\ref{read_gfs}) \newline
!      Interpolates GFS data to LDT grid
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[read\_gfs\_elev](\ref{read_gfs_elev}) \newline
!    reads the native elevation of the GFS data to be used
!    for topographic adjustments to the forcing 
!  \end{description}
!EOP

  integer :: t,f
  integer :: ferror
  integer :: try
  integer, parameter :: ndays = 10  ! # of days to look back for forcing data
  integer :: order     ! 1 indicates lesser interpolation boundary time
                       ! 2 indicates greater interpolation boundary time
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  integer :: movetime  ! 1=move time2 into time1
  integer :: nforce    ! Number of forcing variables for model init. option
  real*8  :: timenow, time1, time2
  real*8  :: dumbtime1, dumbtime2
  real    :: gmt1, gmt2,ts1,ts2
  character(len=LDT_CONST_PATH_LEN) :: name00, name03, name06
  logical :: file_exists, file_exists1, file_exists2
  integer :: option
  real :: gridDesci(20)
  integer :: nstep
  logical :: F06flag
  integer :: status

  gfs_struc(n)%findtime1=0
  gfs_struc(n)%findtime2=0

  movetime=0

!-----------------------------------------------------------------
! Determine the correct number of forcing variables
!-----------------------------------------------------------------
  nstep = LDT_get_nstep(LDT_rc,n)
  nforce = gfs_struc(n)%nmif

  if ( nstep == 0 .or. LDT_rc%rstflag(n)== 1) then
    gfs_struc(n)%findtime1=1
    gfs_struc(n)%findtime2=1
    movetime=0        ! movetime is not properly set at time-step = 1
    LDT_rc%rstflag(n) = 0
  endif 

!-----------------------------------------------------------------
! Determine required GFS data times 
! (previous assimilation, current & future assimilation hours)
! The adjustment of the hour and the direction will be done
! in the subroutines that generate the names
!-----------------------------------------------------------------
  yr1 = LDT_rc%yr  !current time
  mo1 = LDT_rc%mo
  da1 = LDT_rc%da
  hr1 = LDT_rc%hr
  mn1 = LDT_rc%mn
  ss1 = 0
  ts1 = 0
  call LDT_tick( timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

  yr1 = LDT_rc%yr  !previous assimilation/forecast hour
  mo1 = LDT_rc%mo
  da1 = LDT_rc%da
  hr1 = 3*(int(real(LDT_rc%hr)/3.0))
  mn1 = 0
  ss1 = 0
  ts1 = 0
  call LDT_tick( time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

  yr2 = LDT_rc%yr  !next assimilation/forecast hour
  mo2 = LDT_rc%mo
  da2 = LDT_rc%da
  hr2 = 3*(int(real(LDT_rc%hr)/3.0))
  mn2 = 0
  ss2 = 0
  ts2 = 3*60*60
  call LDT_tick( time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )
!-----------------------------------------------------------------
! Use these if need to roll back time.
!-----------------------------------------------------------------
  dumbtime1 = time1
  dumbtime2 = time2
  if ( nstep == 0 .or. nstep == 1 ) then
     gfs_struc(n)%gfstime1 = time1
     gfs_struc(n)%gfstime2 = time2
  endif
!-----------------------------------------------------------------
! Check to see if current time (timenow) has crossed past gfstime2,
! requiring that both gfstime2 be assigned to gfstime1 and a new
! gfstime2 be set 3 or 6 hours ahead of the current gfstime2.
!-----------------------------------------------------------------
  if ( timenow > gfs_struc(n)%gfstime2 ) then
     movetime  = 1
     gfs_struc(n)%findtime2 = 1
  end if
  
  if ( time1 >= gfs_struc(n)%griduptime1 .and. &
       time1 < gfs_struc(n)%griduptime2 .and. & 
       gfs_struc(n)%gridchange1) then 

     write(LDT_logunit,*) 'MSG: get_gfs -- changing grid to 1991 -- 2000'

     gfs_struc(n)%nc = 384
     gfs_struc(n)%nr = 190
     gfs_struc(n)%mi = gfs_struc(n)%nc*gfs_struc(n)%nr
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
     gridDesci = 0 
     gridDesci(1) = 4
     gridDesci(2) = 384
     gridDesci(3) = 190
     gridDesci(4) = 89.277
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.277
     gridDesci(8) = -0.938
     gridDesci(9) = 0.938
     gridDesci(10) = 95
     gridDesci(20) = 0

     if(trim(LDT_rc%met_ecor(findex)).ne."none") then 
        call read_gfs_elev(n,findex, 1)
     endif

     if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 
        call bilinear_interp_input(n, gridDesci, &
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221)
     elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
        call bilinear_interp_input(n, gridDesci, &
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221)
        call conserv_interp_input(n, gridDesci, &
             gfs_struc(n)%n112,gfs_struc(n)%n122,&
             gfs_struc(n)%n212,gfs_struc(n)%n222,&
             gfs_struc(n)%w112,gfs_struc(n)%w122,&
             gfs_struc(n)%w212,gfs_struc(n)%w222)
     endif
     gfs_struc(n)%gridchange1 = .false.

  elseif ( time1 >= gfs_struc(n)%griduptime2 .and. &
           time1 < gfs_struc(n)%griduptime3 .and. & 
           gfs_struc(n)%gridchange2) then 

     write(LDT_logunit,*) 'MSG: get_gfs -- changing grid to 2000 -- 2002'

     gfs_struc(n)%nc = 512
     gfs_struc(n)%nr = 256
     gfs_struc(n)%mi = gfs_struc(n)%nc*gfs_struc(n)%nr
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 512
     gridDesci(3) = 256
     gridDesci(4) = 89.463
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.463
     gridDesci(8) = -0.703
     gridDesci(9) = 0.703
     gridDesci(10) = 128
     gridDesci(20) = 0.0

     if(trim(LDT_rc%met_ecor(findex)).ne."none") then 
        call read_gfs_elev(n,findex, 2)
     endif

     if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 
        call bilinear_interp_input(n, gridDesci, &
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221)

     elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
        call bilinear_interp_input(n, gridDesci,&
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221)
        call conserv_interp_input(n, gridDesci, &
             gfs_struc(n)%n112,gfs_struc(n)%n122,&
             gfs_struc(n)%n212,gfs_struc(n)%n222,&
             gfs_struc(n)%w112,gfs_struc(n)%w122,&
             gfs_struc(n)%w212,gfs_struc(n)%w222)
     endif
     gfs_struc(n)%gridchange2 = .false.

  elseif ( time1 >= gfs_struc(n)%griduptime3 .and. & 
           time1 < gfs_struc(n)%griduptime4 .and. &
           gfs_struc(n)%gridchange3) then 

     write(LDT_logunit,*) 'MSG: get_gfs -- changing grid to 2002 -- 2005'

     gfs_struc(n)%nc = 768
     gfs_struc(n)%nr = 384
     gfs_struc(n)%mi = gfs_struc(n)%nc*gfs_struc(n)%nr
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 768
     gridDesci(3) = 384
     gridDesci(4) = 89.642
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.642
     gridDesci(8) = -0.469
     gridDesci(9) = 0.469
     gridDesci(10) = 192
     gridDesci(20) = 0.0

     if(trim(LDT_rc%met_ecor(findex)).ne."none") then 
        call read_gfs_elev(n,findex, 3)
     endif

     if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 
        call bilinear_interp_input(n, gridDesci, &
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221)
     elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
        call bilinear_interp_input(n, gridDesci, &
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221)
        call conserv_interp_input(n, gridDesci, &
             gfs_struc(n)%n112,gfs_struc(n)%n122,&
             gfs_struc(n)%n212,gfs_struc(n)%n222,&
             gfs_struc(n)%w112,gfs_struc(n)%w122,&
             gfs_struc(n)%w212,gfs_struc(n)%w222)
     endif
     gfs_struc(n)%gridchange3 = .false.

  elseif ( time1 >= gfs_struc(n)%griduptime4 .and. &
           time1 < gfs_struc(n)%griduptime5 .and. &
           gfs_struc(n)%gridchange4) then

     write(LDT_logunit,*) 'MSG: get_gfs -- changing grid to 2005 -- 2010'

     gfs_struc(n)%nc = 1152
     gfs_struc(n)%nr = 576
     gfs_struc(n)%mi = gfs_struc(n)%nc*gfs_struc(n)%nr
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 1152
     gridDesci(3) = 576
     gridDesci(4) = 89.761
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.761
     gridDesci(8) = -0.313
     gridDesci(9) = 0.313
     gridDesci(10) = 288
     gridDesci(20) = 0.0

     if(trim(LDT_rc%met_ecor(findex)).ne."none") then 
        call read_gfs_elev(n,findex, 4)
     endif

     if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 
        call bilinear_interp_input(n, gridDesci, &
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221)
     elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
        call bilinear_interp_input(n, gridDesci, &
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221)
        call conserv_interp_input(n, gridDesci, &
             gfs_struc(n)%n112,gfs_struc(n)%n122,&
             gfs_struc(n)%n212,gfs_struc(n)%n222,&
             gfs_struc(n)%w112,gfs_struc(n)%w122,&
             gfs_struc(n)%w212,gfs_struc(n)%w222)
     endif
     gfs_struc(n)%gridchange4 = .false.

  elseif ( time1 >= gfs_struc(n)%griduptime5 .and. &
           gfs_struc(n)%gridchange5 ) then

     write(LDT_logunit,*) 'MSG: get_gfs -- changing grid to 2010 --'

     gfs_struc(n)%nc = 1760
     gfs_struc(n)%nr = 880
     gfs_struc(n)%mi = gfs_struc(n)%nc*gfs_struc(n)%nr
!-------------------------------------------------------------------
! Reinitialize the weights and neighbors
!-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 1760
     gridDesci(3) = 880
     gridDesci(4) = 89.844
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.844
     gridDesci(8) = -0.205
     gridDesci(9) = 0.205
     gridDesci(10) = 440
     gridDesci(20) = 0.0
  
     if(trim(LDT_rc%met_ecor(findex)).ne."none") then
        call read_gfs_elev(n,findex, 5)
     endif
  
     if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then
        call bilinear_interp_input(n, gridDesci, &
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221)
     elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then
        call bilinear_interp_input(n, gridDesci, &
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221) 
        call conserv_interp_input(n, gridDesci, &
             gfs_struc(n)%n112,gfs_struc(n)%n122,&
             gfs_struc(n)%n212,gfs_struc(n)%n222,&
             gfs_struc(n)%w112,gfs_struc(n)%w122,&
             gfs_struc(n)%w212,gfs_struc(n)%w222)
     endif
     gfs_struc(n)%gridchange5 = .false.
  endif
!-----------------------------------------------------------------
! Establish bookend 1: The code sets up the filenames required to 
! read both instantaneous and time average fields.  
!-----------------------------------------------------------------
  if ( gfs_struc(n)%findtime1 == 1 ) then  !get new time1 from the past
     write(LDT_logunit,*) 'Getting new time1 data'
     ferror = 0
     order = 1
     try = 0
     ts1 = -24*60*60
     status = 0 
     do
        if ( ferror /= 0 ) exit
        try = try+1

        call create_gfsfilename( order, name00, name03, name06, F06flag, &
             gfs_struc(n)%gfsdir, yr1, mo1, da1, hr1 )

        inquire(file=name00,exist=file_exists1) 
        inquire(file=name03,exist=file_exists2)
!-----------------------------------------------------------------
! look for backup files if F0 forecast is not available. 
!-----------------------------------------------------------------
        if(.not.file_exists1) then !f0 forecast does not exist
           call create_gfs_f0backup_filename(order, name00, &
                gfs_struc(n)%gfsdir, yr1,mo1,da1,hr1,status)
           inquire(file=name00,exist=file_exists1) 
           if(.not.file_exists1) then 
              status = 1
           endif
        endif
               
        if(status.eq.0) then 
           write(LDT_logunit,*) 'Reading GFS file1 (I) ',trim(name00)
           write(LDT_logunit,*) 'Reading GFS file1 (A) ',trim(name03)
           if(F06flag) then 
              write(LDT_logunit,*) 'Reading GFS file1 (A)',trim(name06)
           endif
           call read_gfs( order, n, name00, name03, name06, F06flag, ferror, try)
           if ( ferror == 1 ) then  
!-----------------------------------------------------------------
! successfully retrieved forcing data
!-----------------------------------------------------------------
              gfs_struc(n)%gfstime1 = time1
           else  
!-----------------------------------------------------------------
! ferror still=0, so roll back one day & start again
!-----------------------------------------------------------------
              call LDT_tick( dumbtime1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
           end if
        endif
        if ( try > ndays ) then  
!-----------------------------------------------------------------
! data gap exceeds 10 days so stop
!-----------------------------------------------------------------
           write(LDT_logunit,*) 'ERROR: GFS data gap exceeds 10 days on file 1'
           call LDT_endrun
        endif
     end do
  end if
!-----------------------------------------------------------------
! Establish gfstime2
!-----------------------------------------------------------------
  if ( movetime == 1 ) then  
     gfs_struc(n)%gfstime1 = gfs_struc(n)%gfstime2
     gfs_struc(n)%findtime2 = 1  
     do f = 1, LDT_rc%met_nf(findex)
        do t = 1, LDT_rc%ngrid(n)
           LDT_forc(n,findex)%metdata1(f,t) = LDT_forc(n,findex)%metdata2(f,t)
        end do
     end do
  end if
  
  if ( gfs_struc(n)%findtime2 == 1 ) then  
     write(LDT_logunit,*) 'Getting new time2 data'
  
     ferror = 0
     order = 2
     try = 0
     ts2 = -24*60*60
!-----------------------------------------------------------------
! determine the required forecast files. 
!-----------------------------------------------------------------
     do
        if ( ferror /= 0 ) exit
        try = try+1
        call create_gfsfilename( order, name00, name03,name06, &
             F06flag, gfs_struc(n)%gfsdir, &
             yr1, mo1, da1, hr1 )
        inquire(file=name00,exist=file_exists1) 
        inquire(file=name03,exist=file_exists2)
!-----------------------------------------------------------------
! look for backup files if F0 forecast is not available. 
!-----------------------------------------------------------------
        if(.not.file_exists1) then !f0 forecast does not exist
           call create_gfs_f0backup_filename(order, name00, &
                gfs_struc(n)%gfsdir, yr1,mo1,da1,hr1,status)
           inquire(file=name00,exist=file_exists1) 
           if(.not.file_exists1) then 
              status = 1
           endif
        endif
        
        if(file_exists1.and.file_exists2) then 
           status = 0 
           write(LDT_logunit,*) 'Reading GFS file2 (I) ',trim(name00)
           write(LDT_logunit,*) 'Reading GFS file2 (A) ',trim(name03)
        else
           status = 1
        endif

        if(F06flag) then 
           write(LDT_logunit,*) 'Reading GFS file2 (A)',trim(name06)
        endif       
        call read_gfs( order, n, name00, name03, name06, F06flag, ferror, try)    
        if ( ferror == 1 ) then  
!-----------------------------------------------------------------
! successfully retrieved forcing data
!-----------------------------------------------------------------
           gfs_struc(n)%gfstime2 = time2
        else  
!-----------------------------------------------------------------
! ferror still=0, so roll back one day & start again
!-----------------------------------------------------------------
           call LDT_tick( dumbtime2, doy2, gmt2, yr1, mo1, da1, hr1, mn1, ss1, ts2 )
        end if
        if ( try > ndays ) then  
!-----------------------------------------------------------------
! data gap exceeds 10 days so stop
!-----------------------------------------------------------------
           write(LDT_logunit,*) 'ERROR: GFS data gap exceeds 10 days on file 2'
           call LDT_endrun
        end if
     end do
  end if
!-----------------------------------------------------------------
! Check for negative values in radiation 
!-----------------------------------------------------------------
  do f = 1,LDT_rc%met_nf(findex)  
     if ( (f == 3) .or. (f == 4) ) then
        do t=1,LDT_rc%ngrid(n)
           if ( (LDT_forc(n,findex)%metdata2(f,t) /= -9999.9) .and.  &
                (LDT_forc(n,findex)%metdata2(f,t) < 0)) then
              LDT_forc(n,findex)%metdata2(f,t) = (-1) * &
                   LDT_forc(n,findex)%metdata2(f,t)
           endif
           if ( (LDT_forc(n,findex)%metdata1(f,t) /= -9999.9) .and.  &
                (LDT_forc(n,findex)%metdata1(f,t) < 0)) then
              LDT_forc(n,findex)%metdata1(f,t) = (-1) * &
                   LDT_forc(n,findex)%metdata1(f,t)
           endif
        enddo
     endif
  enddo
end subroutine get_gfs

 
