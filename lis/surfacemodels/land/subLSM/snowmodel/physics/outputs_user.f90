! outputs_user.f90

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine OUTPUTS_USER(nx,ny,iter,Tair_grid,rh_grid,&
     &  uwind_grid,vwind_grid,windspd_grid,winddir_grid,&
     &  Qsi_grid,Qli_grid,prec_grid,Tsfc,Qle,Qh,Qe,Qc,Qm,Qf,&
     &  e_balance,snow_depth,xro_snow,swe_depth,ro_nsnow,&
     &  runoff,rain,sprec,sum_prec,sum_runoff,w_balance,&
     &  snow_d,topo_land,wbal_qsubl,sum_sprec,wbal_salt,&
     &  wbal_susp,ro_snow_grid,sum_Qcs,canopy_int,Qcs,&
     &  iyear,imonth,iday,xhour,undef,deltax,xmn,ymn,&
     &  wbal_subgrid,canopy_unload,sum_qsubl,sum_trans,&
     &  sum_unload,sum_glacmelt,glacier_melt,swemelt,&
     &  sfc_pressure,sum_swemelt,albedo,nrecs_max,&
     &  icorr_factor_loop,swesublim,vegtype,iter_start,&
     &  seaice_run,print_inc,cloud_frac_grid,&
     &  output_path_wo_assim,output_path_wi_assim,print_var,&
     &  print_outvars,Qsubl_depth,Qsalt,Qsusp)

! This subroutine is available to provide user-defined outputs.
!   These might be special-case situations, like just writing out
!   data at the end of every day, writing out a few grid cells,
!   saving each data arrays to individual files, etc.

      use snowmodel_inc
      implicit none

      integer i,j,nx,ny,iter,max_poly,num_poly,iyear,imonth,iday,&
     &  icorr_factor_loop,iter_start,icorr_loop_new,&
     &  individual_files,k,irec

      real Tair_grid(nx,ny),rh_grid(nx,ny),&
     &  uwind_grid(nx,ny),vwind_grid(nx,ny),&
     &  windspd_grid(nx,ny),winddir_grid(nx,ny),&
     &  Qsi_grid(nx,ny),Qli_grid(nx,ny),&
     &  prec_grid(nx,ny),Tsfc(nx,ny),&
     &  Qle(nx,ny),Qh(nx,ny),Qe(nx,ny),&
     &  Qc(nx,ny),Qm(nx,ny),Qf(nx,ny),&
     &  e_balance(nx,ny),snow_depth(nx,ny),&
     &  xro_snow(nx,ny),swe_depth(nx,ny),&
     &  ro_nsnow(nx,ny),runoff(nx,ny),&
     &  rain(nx,ny),sprec(nx,ny),&
     &  sum_prec(nx,ny),sum_runoff(nx,ny),&
     &  w_balance(nx,ny),snow_d(nx,ny),&
     &  topo_land(nx,ny),wbal_qsubl(nx,ny),&
     &  sum_sprec(nx,ny),wbal_salt(nx,ny),&
     &  wbal_susp(nx,ny),ro_snow_grid(nx,ny),&
     &  sum_Qcs(nx,ny),canopy_int(nx,ny),&
     &  Qcs(nx,ny),wbal_subgrid(nx,ny),&
     &  canopy_unload(nx,ny),sum_qsubl(nx,ny),&
     &  sum_trans(nx,ny),glacier_melt(nx,ny),&
     &  sum_unload(nx,ny),sum_glacmelt(nx,ny),&
     &  swemelt(nx,ny),sfc_pressure(nx,ny),&
     &  sum_swemelt(nx,ny),swesublim(nx,ny),&
     &  vegtype(nx,ny),albedo(nx,ny),&
     &  cloud_frac_grid(nx,ny),Qsubl_depth(nx,ny),&
     &  Qsalt(nx,ny),Qsusp(nx,ny)

      real undef,xhour,deltax,pi,rad2deg,seaice_run,print_inc
      double precision xmn,ymn
      double precision nrecs_max,nrecs

      real uwnd(nx,ny)
      real vwnd(nx,ny)

! Define the output variable data block.  Note that the print_var
!   "yes/no" array was generated in readparam_code.f.
!      real vars(nx_max,ny_max,n_print_vars)
      real vars(nx,ny,n_print_vars)  ! KRA

      character*80 output_path_wo_assim,output_path_wi_assim
      character*1 print_var(n_print_vars)
      character*4 print_outvars(n_print_vars)

! This was defined in preprocess_code.f, in subroutine mk_ctl_files.
!     data print_outvars /'tair','relh','wspd','qsin','qlin',
!    &                    'qlem','albd','wdir','prec','rpre',
!    &                    'spre','smlt','ssub','roff','glmt',
!    &                    'snod','sden','swed','sspr','ssmt',
!    &                    'cldf','var1','var2','var3','var4',
!    &                    'var5','var6','var7','var8','var9'/
      
! These now come in from the .par file.
!     character path1*(*) 
!     character path2*(*) 

      integer i_trailing_blanks,trailing_blanks,i_len_wo,i_len_wi

! Calculate how long the paths are.
      i_len_wo = 80 - trailing_blanks(output_path_wo_assim)
      i_len_wi = 80 - trailing_blanks(output_path_wi_assim)
!     print *, i_len_wo,i_len_wi
!     print *, output_path_wo_assim(1:i_len_wo)
!     print *, output_path_wi_assim(1:i_len_wi)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! BEGIN USER EDIT SECTION.

! Define which variables you want to save, whether they will be
!   output at every time step or some time-step increment, which
!   directories the data will be put in, etc.

! Define the output file locations (paths).
! These now come in from the .par file.
!     parameter (path1 =
!    &  'outputs/wo_assim/') 
!     parameter (path2 =
!    &  'outputs/wi_assim/') 

! Write a seperate file for each variable (individual_files = 1).
!   No other option has been implemented here.
      individual_files = 1

! Define the number of time steps you are going to sum or average
!   over.  If you want to output data at every model time step, set
!   print_inc = 1.0.  For run with an hourly time step and data
!   writes once a day, print_inc = 24.0.  For a run with 3-hourly
!   time steps and data writes once a day, print_inc = 8.0.
! This now comes in from the .par file.
!     print_inc = 24.0
!     print_inc = 1.0
!     print_inc = 8.0

! Define the variables you want to save.  The following are the
!   variables this subroutine is currently set up to output.
!   Listed are the output variable name and the corresponding model
!   variable name.
! Note that these "yes/no" flags are now defined in snowmodel.par.

! VALUES AVERAGED OVER THE PERIOD.
!    1   tair(i,j) = Tair_grid(i,j) - 273.15
!    2   relh(i,j) = rh_grid(i,j)
!    3   wspd(i,j) = windspd_grid(i,j)
!    4   qsin(i,j) = Qsi_grid(i,j)
!    5   qlin(i,j) = Qli_grid(i,j)
!    6   qlem(i,j) = Qle(i,j)
!    7   albd(i,j) = albedo(i,j)
!    8   wdir(i,j) = from uwind_grid(i,j) and vwind_grid(i,j)

! VALUES SUMMED OVER THE PERIOD.
!    9   prec(i,j) = prec_grid(i,j)
!   10   rpre(i,j) = rain(i,j)
!   11   spre(i,j) = sprec(i,j)
!   12   smlt(i,j) = swemelt(i,j)
!   13   ssub(i,j) = swesublim(i,j)
!   14   roff(i,j) = runoff(i,j)
!   15   glmt(i,j) = glacier_melt(i,j)

! VALUES SAVED AT THE END OF THE PERIOD.
!   16   snod(i,j) = snow_depth(i,j)
!   17   sden(i,j) = xro_snow(i,j)
!   18   swed(i,j) = swe_depth(i,j)
!   19   sspr(i,j) = sum_sprec(i,j)
!   20   ssmt(i,j) = sum_swemelt(i,j)

! NEW VARIABLES ADDED TO THE OUTPUT DATA BLOCK.
!   (you have to modify the code below to
!    do the processing you want for these)
!   21   cldf(i,j) = cloud_fraction_grid(i,j)

! Define which variables you want to save by placing a yes = 'y' or
!   no = 'n' in front of the variable number.

! VALUES AVERAGED OVER THE PERIOD.
! VALUES AVERAGED OVER THE PERIOD.
! 1 = tair(i,j) = Tair_grid(i,j) - 273.15
!     print_var(1)  = 'y'

! 2 = relh(i,j) = rh_grid(i,j)
!     print_var(2)  = 'n'

! 3 = wspd(i,j) = windspd_grid(i,j)
!     print_var(3)  = 'n'

! 4 = qsin(i,j) = Qsi_grid(i,j)
!     print_var(4)  = 'n'

! 5 = qlin(i,j) = Qli_grid(i,j)
!     print_var(5)  = 'n'

! 6 = qlem(i,j) = Qle(i,j)
!     print_var(6)  = 'n'

! 7 = albd(i,j) = albedo(i,j)
!     print_var(7)  = 'n'

! 8 = wdir(i,j) = from uwind_grid(i,j) and vwind_grid(i,j)
!     print_var(8)  = 'n'

! VALUES SUMMED OVER THE PERIOD.
! VALUES SUMMED OVER THE PERIOD.
!  9 = prec(i,j) = prec_grid(i,j)
!     print_var(9)  = 'y'

! 10 = rpre(i,j) = rain(i,j)
!     print_var(10) = 'y'

! 11 = spre(i,j) = sprec(i,j)
!     print_var(11) = 'y'

! 12 = smlt(i,j) = swemelt(i,j)
!     print_var(12) = 'n'

! 13 = ssub(i,j) = swesublim(i,j)
!     print_var(13) = 'n'

! 14 = roff(i,j) = runoff(i,j)
!     print_var(14) = 'n'

! 15 = glmt(i,j) = glacier_melt(i,j)
!     print_var(15) = 'n'

! VALUES SAVED AT THE END OF THE PERIOD.
! VALUES SAVED AT THE END OF THE PERIOD.
! 16 = snod(i,j) = snow_depth(i,j)
!     print_var(16) = 'y'

! 17 = sden(i,j) = xro_snow(i,j)
!     print_var(17) = 'y'

! 18 = swed(i,j) = swe_depth(i,j)
!     print_var(18) = 'y'

! 19 = sspr(i,j) = sum_sprec(i,j)
!     print_var(19) = 'y'

! 20 = ssmt(i,j) = sum_swemelt(i,j)
!     print_var(20) = 'y'

! NEW VARIABLES ADDED TO THE OUTPUT DATA BLOCK.
! NEW VARIABLES ADDED TO THE OUTPUT DATA BLOCK.
!   (you have to modify the code below to
!    do the processing you want for these)
! 21 = cldf(i,j) = cloud_fraction_grid(i,j)
!     print_var(21) = 'n'

! Extra variables.
!     print_var(22) = 'n'
!     print_var(23) = 'n'
!     print_var(24) = 'n'
!     print_var(25) = 'n'
!     print_var(26) = 'n'
!     print_var(27) = 'n'
!     print_var(28) = 'n'
!     print_var(29) = 'n'
!     print_var(30) = 'n'

! Note that this data output implementation is currently configured
!   to mask out the ocean points (vegtype.eq.24.0) if this is a
!   land run (seaice_run = 0.0); mask out all land points (vegtype.
!   ne.24.0) if this is an ocean/sea ice run (seaice_run = 1.0, 3.0,
!   or 4.0); and to not mask out anything if this is a combined land
!   and sea ice run (seaice_run = 2.0).

! END USER EDIT SECTION.

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Define the constants used in the wind-direction averaging.
      pi = 2.0 * acos(0.0)
      rad2deg = 180.0 / pi

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      uwnd = 0.0
      vwnd = 0.0

! Use individual output files for each variable.
      if (individual_files.eq.1) then

        if (iter.eq.iter_start) then
          nrecs = nx * ny
          if (nrecs.gt.nrecs_max) then
            print *,'Your simulation domain has too many grid cells'
            print *,'to print the .gdat files like the write statements'
            print *,'are currently configured.  You will have to change'
            print *,'them to look like:'
            print *,'    do j=1,ny'
            print *,'      write (51,rec=j) (var(i,j),i=1,nx)'
            print *,'    enddo'
            stop
          endif
        endif

! Open individual output files for each variable.
        if (iter.eq.iter_start) then
          if (icorr_factor_loop.eq.1) then
            do k=1,n_print_vars
              if (print_var(k).eq.'y') then
                open (220+k,&
     &file=output_path_wo_assim(1:i_len_wo)//print_outvars(k)//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*ny,&
     &            status='replace')
              endif
            enddo
          endif

          if (icorr_factor_loop.eq.2) then
            do k=1,n_print_vars
              if (print_var(k).eq.'y') then
                open (320+k,&
     &file=output_path_wi_assim(1:i_len_wi)//print_outvars(k)//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*ny,&
     &            status='replace')
              endif
            enddo
          endif
        endif

        if (iter.eq.iter_start) then
! Initialize the averaging and summing arrays.
          do j=1,ny
            do i=1,nx
              do k=1,n_print_vars
                vars(i,j,k) = 0.0
              enddo
              uwnd(i,j) = 0.0
              vwnd(i,j) = 0.0
            enddo
          enddo
        endif

! Perform the avaraging, summing, etc.
        do j=1,ny
          do i=1,nx
! Values averaged over the period.
            vars(i,j,1) = vars(i,j,1) + (Tair_grid(i,j) - 273.15) / &
     &        print_inc
            vars(i,j,2) = vars(i,j,2) + rh_grid(i,j) / print_inc
            vars(i,j,3) = vars(i,j,3) + windspd_grid(i,j) / print_inc
            vars(i,j,4) = vars(i,j,4) + Qsi_grid(i,j) / print_inc
            vars(i,j,5) = vars(i,j,5) + Qli_grid(i,j) / print_inc
            vars(i,j,6) = vars(i,j,6) + Qle(i,j) / print_inc
            vars(i,j,7) = vars(i,j,7) + albedo(i,j) / print_inc

            uwnd(i,j) = uwnd(i,j) + uwind_grid(i,j) / print_inc
            vwnd(i,j) = vwnd(i,j) + vwind_grid(i,j) / print_inc

! Some compilers do not allow both u and v to be 0.0 in
!   the atan2 computation.
            if (abs(uwnd(i,j)).lt.1e-10) uwnd(i,j) = 1e-10

            vars(i,j,8) = rad2deg * atan2(uwnd(i,j),vwnd(i,j))
            if (vars(i,j,8).ge.180.0) then
              vars(i,j,8) = vars(i,j,8) - 180.0
            else
              vars(i,j,8) = vars(i,j,8) + 180.0
            endif

! Values summed over the period.
            vars(i,j,9) = vars(i,j,9) + prec_grid(i,j)
            vars(i,j,10) = vars(i,j,10) + rain(i,j)
            vars(i,j,11) = vars(i,j,11) + sprec(i,j)
            vars(i,j,12) = vars(i,j,12) + swemelt(i,j)
            vars(i,j,13) = vars(i,j,13) + swesublim(i,j)
            vars(i,j,14) = vars(i,j,14) + runoff(i,j)
            vars(i,j,15) = vars(i,j,15) + glacier_melt(i,j)

! Values saved at the end of the day.
            vars(i,j,16) = snow_depth(i,j)
            vars(i,j,17) = xro_snow(i,j)
            vars(i,j,18) = swe_depth(i,j)
            vars(i,j,19) = sum_sprec(i,j)
            vars(i,j,20) = sum_swemelt(i,j)

! Values averaged over the period.
          vars(i,j,21) = vars(i,j,21) + cloud_frac_grid(i,j) / print_inc

! New variables.
            vars(i,j,22) = vars(i,j,22) + Qsubl_depth(i,j)
            vars(i,j,23) = sum_trans(i,j)
            vars(i,j,24) = vars(i,j,24) + Qsalt(i,j)
            vars(i,j,25) = vars(i,j,25) + Qsusp(i,j)

          enddo
        enddo

! Check to see whether this is the data-write time step.
        if (mod(iter,nint(print_inc)).eq.0) then

! Mask out the ocean points (vegtype.eq.24.0) if this is
!   a land run (seaice_run = 0.0).  Mask out all land points
!   (vegtype.ne.24.0) if this is an ocean/sea ice run
!   (seaice_run = 1.0, 3.0, or 4.0).  Do not mask out anything
!   if this this is a combined land and sea ice run (seaice_run
!   = 2.0).
          if (seaice_run.eq.0.0) then
            do j=1,ny
              do i=1,nx
                if (vegtype(i,j).eq.24.0) then
                  do k=1,n_print_vars
                    vars(i,j,k) = undef
                  enddo
                endif
              enddo
            enddo
          elseif (seaice_run.eq.1.0 .or. seaice_run.eq.3.0 .or. &
     &      seaice_run.eq.4.0) then
            do j=1,ny
              do i=1,nx
                if (vegtype(i,j).ne.24.0) then
                  do k=1,n_print_vars
                    vars(i,j,k) = undef
                  enddo
                endif
              enddo
            enddo
          endif

! Write out the data.
          irec = iter / nint(print_inc)
          if (icorr_factor_loop.eq.1) then
            do k=1,n_print_vars
              if (print_var(k).eq.'y') then
                write (220+k,rec=irec) ((vars(i,j,k),i=1,nx),j=1,ny)
              endif
            enddo
          elseif (icorr_factor_loop.eq.2) then
            do k=1,n_print_vars
              if (print_var(k).eq.'y') then
                write (320+k,rec=irec) ((vars(i,j,k),i=1,nx),j=1,ny)
              endif
            enddo
          endif

! Reinitialize the averaging and summing arrays.
          do j=1,ny
            do i=1,nx
              do k=1,n_print_vars
                vars(i,j,k) = 0.0
              enddo
              uwnd(i,j) = 0.0
              vwnd(i,j) = 0.0
            enddo
          enddo
        endif

! Use more than one variable in an output file.
      else

        print *,'Use more than one variable in an output file:'
        print *,'  THIS HAS NOT BEEN IMPLEMENTED YET'

      endif

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! THIS IS AN EXAMPLE OF SAVING DATA IN ASCII/TEXT FORMAT.

! I have completely removed this example; I now do this with
!   improved codes as part of post-processing steps.  It is
!   just too slow to do it as part of the model simulation.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THE CODE BELOW WAS USED TO SAVE AVERAGES OVER POLYGONS.

! I have completely removed this example; if I were to do this
!   again I would do it as a post-processing step.  If you really
!   want to see what I did here, you can look in one of the pre-
!   2018 code distributions.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

