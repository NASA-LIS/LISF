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
! !ROUTINE: read_era5
! \label{read_era5}
! 
! !REVISION HISTORY:
! 23 dec 2019: Sujay Kumar, initial code 
! 15 apr 2025: Hiroko Beaudoing, adopted ERA5 routines for the public CDS
!                               data format
!
! !INTERFACE:
subroutine read_era5cds(n, kk, order, year, month, day, hour, read_flag, findex,&
     instfile, avgfile, lmlfile, prevavgfile, ferror)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_forc
  use era5cds_forcingMod, only : era5cds_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: kk
  integer, intent(in)          :: order ! lower(1) or upper(2) time interval bdry
  integer, intent(in)          :: year
  integer, intent(in)          :: month
  integer, intent(in)          :: day
  integer, intent(in)          :: hour
  logical, intent(in)          :: read_flag
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: instfile, avgfile, lmlfile, prevavgfile
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  ERA5 data, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain. \newline
!
!  ERA5 model output variables used to force LIS are provided in 3 files:
!  1) inst, Instantaneous values, available every hour.
!  2) accum, Time integrated values (accumulations) over the past hour.
!  3) lmlinst, Lowest Model Level Instantaneous values, available every hour.
!
!  ERA5 FORCING VARIABLES:
!  1. T2M       inst,  Atmospheric temperature (K) (2m)
!  2. D2M       inst,  Dewpoint temperature (K) (2m)
!  3. SSRD      accum, surface solar radiation downwards [J m**-2]
!  4. STRD      accum, surface thermal radiation downwards [J m**-2]
!  5. U10M      inst, zonal wind [m/s] (10m)
!  6. V10M      inst, meridional wind [m/s] (10m)
!  7. SP        inst, surface pressure [Pa]
!  8. LSP       accum, large scale precipitation [m]
!  9. CP        accum, convective precipitation [m]
!  ERA5 LOWEST MODEL LEVEL VARIABLES:
!  1. T         inst,  air temperature (K)
!  2. Q         inst,  specific_humidity (kg kg**-1)
!  5. U         inst, zonal wind [m/s] 
!  6. V         inst, meridional wind [m/s]
!
!  also inst file includes:
!  Z: Geopotential (m**2 s**-2)
!
!  NOTE 1: be aware that ECMWF outputs large-scale and convective precipitation
!  separately.  For total precipitation, need to sum the two fields,
!  LSP+CP=TP. 
!  NOTE 2: SW & LW flux accumulations from time1 are interpolated via new zterp.
!  NOTE 3: accum files starts at 7z on the 1st day of month and ends at 6z 
!          on the 1st day of next month.

!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, initial read, 
!    order=2, assign the next 1 hourly instance )
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the 1 hour ERA5 analysis file
!  \item[tscount]
!    time step count
!  \item[ferror]
!    return error code (1 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[era5grid\_2\_lisgrid](\ref{era5grid_2_lisgrid}) \newline
!    reverse south-north direction and shift 180 deg east-west
!  \item[interp\_era5cds\_var](\ref{interp_era5cds_var}) \newline
!    spatially interpolate the forcing data using interpolation of choice
!  \item[assign\_processed\_era5cdsf](\ref{assign_processed_era5cdsf}) \newline
!    assigns interpolated forcing data to the module data structure 
!  \end{description}
!EOP
  
  integer   :: ftn
  integer   :: tmpId, qId, uwindId, vwindId, lwdId, psId, rainfId, crainfId
  integer   :: swdId, timeId
  integer   :: c,r,t,k,l,i,j,ll
  integer   :: tindex,atindex
  integer   :: mo,rec_size, prev_rec_size
  integer   :: start_time_index
  logical   :: file_exists
  real      :: missingValue

  real, allocatable      :: tair(:,:,:)
  real, allocatable      :: qair(:,:,:)
  real, allocatable      :: swd(:,:,:)
  real, allocatable      :: lwd(:,:,:)
  real, allocatable      :: uwind(:,:,:)
  real, allocatable      :: vwind(:,:,:)
  real, allocatable      :: ps(:,:,:)
  real, allocatable      :: rainf(:,:,:)
  real, allocatable      :: crainf(:,:,:)
  real, allocatable      :: data4d(:,:,:,:)
  real, allocatable      :: var1(:,:,:)

  real, allocatable, dimension(:,:)  :: datalis

  integer, parameter :: n_last_steps = 7
  integer :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/
  real, parameter :: epsln = 0.622        !Constants' values taken from
  real, parameter :: A = 2.53E8, B=5.42E3 !Roger&Yau, A Short Course
                                          !in Cloud Physics, pp.12-17
                                          ! A in kPa, B in K.
  real            :: p_kPa
! __________________________________________________________________________

  ferror = 0 ! 1 success

  if(read_flag) then
#if (defined USE_NETCDF4)

     if((mod(year,4) .eq. 0 .and. mod(year, 100).ne.0) &!leap year
          .or.(mod(year,400) .eq.0)) then
        days(2) = 29
     else
        days(2) = 28
     endif

     mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
     rec_size = days(month)*24 

     allocate(datalis(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold))
     !=== If using the lowest model level forcing (inst and mlinst):
     !=== "t","q","u","v", and "sp " ====
     if(era5cds_struc(n)%uselml.eq.1) then
       inquire(file=lmlfile,exist=file_exists)
       if(file_exists) then
          era5cds_struc(n)%ps     = LIS_rc%udef
          era5cds_struc(n)%tair   = LIS_rc%udef
          era5cds_struc(n)%qair   = LIS_rc%udef
          era5cds_struc(n)%uwind  = LIS_rc%udef
          era5cds_struc(n)%vwind  = LIS_rc%udef
          write(LIS_logunit,*)'[INFO] Reading ERA5 file (bookend,', order,' -',trim(instfile), ')'

          call LIS_verify(nf90_open(path=trim(instfile), mode=NF90_NOWRITE, &
               ncid=ftn), 'nf90_open failed in read_era5cds')

          allocate(ps(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,rec_size))

          call LIS_verify(nf90_inq_varid(ftn,'sp',psId), &
               'nf90_inq_varid failed for psurf in read_era5cds')
          call LIS_verify(nf90_get_var(ftn,psId, ps),&
               'nf90_get_var failed for ps in read_era5cds')
          call LIS_verify(nf90_get_att(ftn,psId, "_FillValue", missingValue), &
               'nf90_get_att failed in read_era5cds for sp FillValue')
          call LIS_verify(nf90_close(ftn), &
               'failed to close file in read_era5cds inst in lml')

          do l=1,rec_size
            call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                   era5cds_struc(n)%nrold,ps(:,:,l),datalis)
            call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                   era5cds_struc(n)%ps(:,l))
          enddo
          deallocate(ps)

          write(LIS_logunit,*)'[INFO] Reading ERA5 file (bookend,', order,' -',trim(lmlfile), ')'
          call LIS_verify(nf90_open(path=trim(lmlfile), mode=NF90_NOWRITE, &
               ncid=ftn), 'nf90_open failed in read_era5cds')

          ! model level data are in 4-dimentional
          allocate(data4d(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,1,rec_size))
          call LIS_verify(nf90_inq_varid(ftn,'t',tmpId), &
               'nf90_inq_varid failed for t in read_era5cds')
          call LIS_verify(nf90_get_var(ftn,tmpId, data4d),&
               'nf90_get_var failed for t in read_era5cds')
          call LIS_verify(nf90_get_att(ftn,tmpId, "_FillValue", missingValue), &
               'nf90_get_att failed in read_era5cds for t FillValue')

          do l=1,rec_size
            call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                   era5cds_struc(n)%nrold,data4d(:,:,1,l),datalis)
            call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                   era5cds_struc(n)%tair(:,l))
          enddo

          call LIS_verify(nf90_inq_varid(ftn,'q',qId), &
               'nf90_inq_varid failed for q in read_era5cds')
          call LIS_verify(nf90_get_var(ftn,qId, data4d),&
               'nf90_get_var failed for q in read_era5cds')
          call LIS_verify(nf90_get_att(ftn,qId, "_FillValue", missingValue), &
               'nf90_get_att failed in read_era5cds for q FillValue')

          do l=1,rec_size
            call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                   era5cds_struc(n)%nrold,data4d(:,:,1,l),datalis)
            call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                   era5cds_struc(n)%qair(:,l))
          enddo

          call LIS_verify(nf90_inq_varid(ftn,'u',uwindId), &
               'nf90_inq_varid failed for u in read_era5cds')
          call LIS_verify(nf90_get_var(ftn,uwindId, data4d),&
               'nf90_get_var failed for u in read_era5cds')
          call LIS_verify(nf90_get_att(ftn,uwindId, "_FillValue", missingValue), &
               'nf90_get_att failed in read_era5cds for u FillValue')

          do l=1,rec_size
            call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                   era5cds_struc(n)%nrold,data4d(:,:,1,l),datalis)
            call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                   era5cds_struc(n)%uwind(:,l))
          enddo

          call LIS_verify(nf90_inq_varid(ftn,'v',vwindId), &
               'nf90_inq_varid failed for v in read_era5cds')
          call LIS_verify(nf90_get_var(ftn,vwindId, data4d),&
               'nf90_get_var failed for v in read_era5cds')
          call LIS_verify(nf90_get_att(ftn,vwindId, "_FillValue", missingValue), &
               'nf90_get_att failed in read_era5cds for v FillValue')

          do l=1,rec_size
            call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                   era5cds_struc(n)%nrold,data4d(:,:,1,l),datalis)
            call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                   era5cds_struc(n)%vwind(:,l))
          enddo
          deallocate(data4d)

          call LIS_verify(nf90_close(ftn), &
               'failed to close file in read_era5cds lml')
          ferror = 1
       else
         write(LIS_logunit,*) '[ERR] ',trim(lmlfile)//' does not exist'
         call LIS_endrun()
       endif

     else  !=== If using 2m and 10m forcing (inst):
           !=== "t2m","d2m","u10","v10","sp " ====

       inquire(file=instfile,exist=file_exists)
       if(file_exists) then
          era5cds_struc(n)%ps     = LIS_rc%udef
          era5cds_struc(n)%tair   = LIS_rc%udef
          era5cds_struc(n)%qair   = LIS_rc%udef
          era5cds_struc(n)%uwind  = LIS_rc%udef
          era5cds_struc(n)%vwind  = LIS_rc%udef

          write(LIS_logunit,*)'[INFO] Reading ERA5 2m10m (bookend,', order,' -',trim(instfile), ')'

          call LIS_verify(nf90_open(path=trim(instfile), mode=NF90_NOWRITE, &
               ncid=ftn), 'nf90_open failed in read_era5cds')

          allocate(ps(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,rec_size))

          call LIS_verify(nf90_inq_varid(ftn,'sp',psId), &
               'nf90_inq_varid failed for psurf in read_era5cds')
          call LIS_verify(nf90_get_var(ftn,psId, ps),&
               'nf90_get_var failed for ps in read_era5cds')
          call LIS_verify(nf90_get_att(ftn,psId, "_FillValue", missingValue), &
               'nf90_get_att failed in read_era5cds for sp FillValue')

          do l=1,rec_size
            call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                   era5cds_struc(n)%nrold,ps(:,:,l),datalis)
            call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                   era5cds_struc(n)%ps(:,l))
          enddo
          deallocate(ps)

          allocate(qair(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,rec_size))
          call LIS_verify(nf90_inq_varid(ftn,'d2m',qId), &
               'nf90_inq_varid failed for d2m in read_era5cds')
          call LIS_verify(nf90_get_var(ftn,qId, qair),&
               'nf90_get_var failed for d2m in read_era5cds')
          call LIS_verify(nf90_get_att(ftn,qId, "_FillValue", missingValue), &
               'nf90_get_att failed in read_era5cds for d2m FillValue')
          do l=1,rec_size
            call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                   era5cds_struc(n)%nrold,qair(:,:,l),datalis)
            call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                   era5cds_struc(n)%qair(:,l))
          enddo
          deallocate(qair)
          !-----------------------------------------------------------------
          ! Calculate specific humidity from Dew point Temp. & Sfc Pressure.
          ! Approximate: q~epsilon*e/p, p~p_sfc, e=e_s(Td), e_s=A*exp(-B/T)
          !-----------------------------------------------------------------
          do l=1,rec_size
            do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 k=c+(r-1)*LIS_rc%lnc(n)
                 if(era5cds_struc(n)%ps(k,l).ne.LIS_rc%udef.and.&
                     era5cds_struc(n)%qair(k,l).ne.LIS_rc%udef) then
                    p_kPa = era5cds_struc(n)%ps(k,l) / 1000.0
                    era5cds_struc(n)%qair(k,l) = epsln*(A * EXP(-B/era5cds_struc(n)%qair(k,l))) / p_kPa
                 endif
              enddo
            enddo
          enddo

          allocate(tair(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,rec_size))
          call LIS_verify(nf90_inq_varid(ftn,'t2m',tmpId), &
               'nf90_inq_varid failed for t2m in read_era5cds')
          call LIS_verify(nf90_get_var(ftn,tmpId, tair),&
               'nf90_get_var failed for t2m in read_era5cds')
          call LIS_verify(nf90_get_att(ftn,tmpId, "_FillValue", missingValue), &
               'nf90_get_att failed in read_era5cds for t2m FillValue')

          do l=1,rec_size
            call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                   era5cds_struc(n)%nrold,tair(:,:,l),datalis)
            call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                   era5cds_struc(n)%tair(:,l))
          enddo
          deallocate(tair)

          allocate(uwind(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,rec_size))
          call LIS_verify(nf90_inq_varid(ftn,'u10',uwindId), &
               'nf90_inq_varid failed for u10 in read_era5cds')
          call LIS_verify(nf90_get_var(ftn,uwindId, uwind),&
               'nf90_get_var failed for u10 in read_era5cds')
          call LIS_verify(nf90_get_att(ftn,uwindId, "_FillValue", missingValue), &
               'nf90_get_att failed in read_era5cds for u10 FillValue')

          do l=1,rec_size
            call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                   era5cds_struc(n)%nrold,uwind(:,:,l),datalis)
            call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                   era5cds_struc(n)%uwind(:,l))
          enddo
          deallocate(uwind)

          allocate(vwind(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,rec_size))
          call LIS_verify(nf90_inq_varid(ftn,'v10',vwindId), &
               'nf90_inq_varid failed for v10 in read_era5cds')
          call LIS_verify(nf90_get_var(ftn,vwindId, vwind),&
               'nf90_get_var failed for v10 in read_era5cds')
          call LIS_verify(nf90_get_att(ftn,vwindId, "_FillValue", missingValue), &
               'nf90_get_att failed in read_era5cds for v10 FillValue')

          do l=1,rec_size
            call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                   era5cds_struc(n)%nrold,vwind(:,:,l),datalis)
            call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                   era5cds_struc(n)%vwind(:,l))
          enddo
          deallocate(vwind)

          call LIS_verify(nf90_close(ftn), &
               'failed to close file in read_era5cds inst')
          ferror = 1
        else
          write(LIS_logunit,*) '[ERR] ',trim(instfile)//' does not exist'
          call LIS_endrun()
        endif

     endif  !=== era5cds_struc(n)%uselml

     !=== Read precipitation and radiation fields in accum file
     !=== "ssrd","strd","lsp","cp" ====
     inquire(file=avgfile,exist=file_exists)
     if(file_exists) then
        era5cds_struc(n)%rainf  = LIS_rc%udef
        era5cds_struc(n)%crainf  = LIS_rc%udef
        era5cds_struc(n)%swd    = LIS_rc%udef
        era5cds_struc(n)%lwd    = LIS_rc%udef

        write(LIS_logunit,*)'[INFO] Reading ERA5 file (bookend,', order,' -',trim(avgfile), ')'
        call LIS_verify(nf90_open(path=trim(avgfile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_era5cds')

        allocate(rainf(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'lsp',rainfId), &
             'nf90_inq_varid failed for lsp in read_era5cds')
        call LIS_verify(nf90_get_var(ftn,rainfId, rainf),&
             'nf90_get_var failed for lsp in read_era5cds')
        call LIS_verify(nf90_get_att(ftn,rainfId, "_FillValue", missingValue), &
             'nf90_get_att failed in read_era5cds for lsp FillValue')

        ! large scale rainfall: lsp [m] -> mm/s
        do l=1,rec_size
          do j=1,era5cds_struc(n)%nrold
            do i=1,era5cds_struc(n)%ncold
               if(rainf(i,j,l).ne.missingValue) then
                  rainf(i,j,l) = rainf(i,j,l)*1000.0/(1.0*60*60) 
               endif
            enddo
          enddo
        enddo

        do l=1,rec_size
          call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                 era5cds_struc(n)%nrold,rainf(:,:,l),datalis)
          call interp_era5cds_var(n,findex,month,datalis,missingValue,.true.,&
                 era5cds_struc(n)%rainf(:,l))
        enddo
        deallocate(rainf)

        allocate(crainf(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'cp',crainfId), &
             'nf90_inq_varid failed for cp in read_era5cds')
        call LIS_verify(nf90_get_var(ftn,crainfId, crainf),&
             'nf90_get_var failed for cp in read_era5cds')
        call LIS_verify(nf90_get_att(ftn,crainfId, "_FillValue", missingValue), &
             'nf90_get_att failed in read_era5cds for cp FillValue')

        ! convective rainfall: cp [m] -> mm/s
        do l=1,rec_size
          do j=1,era5cds_struc(n)%nrold
            do i=1,era5cds_struc(n)%ncold
               if(crainf(i,j,l).ne.missingValue) then
                  crainf(i,j,l) = crainf(i,j,l)*1000.0/(1.0*60*60)
               endif
            enddo
          enddo
        enddo

        do l=1,rec_size
          call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                 era5cds_struc(n)%nrold,crainf(:,:,l),datalis)
          call interp_era5cds_var(n,findex,month,datalis,missingValue,.true.,&
                 era5cds_struc(n)%crainf(:,l))
        enddo
        deallocate(crainf)

        !-----------------------------------------------------------------
        ! rain = lsp + cp [mm/s]
        !-----------------------------------------------------------------
        do l=1,rec_size
          do r=1,LIS_rc%lnr(n)
            do c=1,LIS_rc%lnc(n)
               k=c+(r-1)*LIS_rc%lnc(n)
               if(era5cds_struc(n)%rainf(k,l).ne.LIS_rc%udef.and.&
                   era5cds_struc(n)%crainf(k,l).ne.LIS_rc%udef) then
                  era5cds_struc(n)%rainf(k,l) = era5cds_struc(n)%rainf(k,l) + &
                                                era5cds_struc(n)%crainf(k,l)
               endif
            enddo
          enddo
        enddo
        write(LIS_logunit,*) 'HKB rainf and crainf done...'

        allocate(swd(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'ssrd',swdId), &
             'nf90_inq_varid failed for ssrd in read_era5cds')
        call LIS_verify(nf90_get_var(ftn,swdId, swd),&
             'nf90_get_var failed for swd in read_era5cds')
        call LIS_verify(nf90_get_att(ftn,swdId, "_FillValue", missingValue), &
             'nf90_get_att failed in read_era5cds for swd FillValue')

        ! swd - convert accumulated field to rate
        do l=1,rec_size
          do j=1,era5cds_struc(n)%nrold
            do i=1,era5cds_struc(n)%ncold
               if(swd(i,j,l).ne.missingValue) then
                  swd(i,j,l) = swd(i,j,l)/(1.0*60*60)
               endif
            enddo
          enddo
        enddo
        do l=1,rec_size
          call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                 era5cds_struc(n)%nrold,swd(:,:,l),datalis)
          call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                 era5cds_struc(n)%swd(:,l))
        enddo
        deallocate(swd)

        allocate(lwd(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'strd',lwdId), &
             'nf90_inq_varid failed for strd in read_era5cds')
        call LIS_verify(nf90_get_var(ftn,lwdId, lwd),&
             'nf90_get_var failed for lwd in read_era5cds')
        call LIS_verify(nf90_get_att(ftn,lwdId, "_FillValue", missingValue), &
             'nf90_get_att failed in read_era5cds for lwd FillValue')

        ! lwd - convert accumulated field to rate
        do l=1,rec_size
          do j=1,era5cds_struc(n)%nrold
            do i=1,era5cds_struc(n)%ncold
               if(lwd(i,j,l).ne.missingValue) then
                  lwd(i,j,l) = lwd(i,j,l)/(1.0*60*60)
               endif
            enddo
          enddo
        enddo
        do l=1,rec_size
          call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                 era5cds_struc(n)%nrold,lwd(:,:,l),datalis)
          call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                 era5cds_struc(n)%lwd(:,l))
        enddo
        deallocate(lwd)

        call LIS_verify(nf90_close(ftn), &
             'failed to close avg file in read_era5cds')
        ferror = 1
     else
        write(LIS_logunit,*) '[ERR] ',trim(avgfile)//' does not exist'
        call LIS_endrun()

     endif
     !=== Read previous precipitation and radiation fields in accum file
     !=== populate prev* arrays
      inquire(file=prevavgfile,exist=file_exists)
      if(file_exists) then
         ! backward fill: read in last 7 time steps of prevavgfile
         write(LIS_logunit,*)'[INFO] Reading prev ERA5 file (bookend,', order,' -',trim(prevavgfile), ')'
         call LIS_verify(nf90_open(path=trim(prevavgfile), mode=NF90_NOWRITE, &
                ncid=ftn), 'nf90_open failed in read_era5cds')
        
         call LIS_verify(nf90_inq_dimid(ftn,'valid_time',timeId), &
              'nf90_inq_dimid failed for timeId in read_era5cds')
         call LIS_verify(nf90_inquire_dimension(ftn,timeId,len=prev_rec_size), &
              'nf90_inquire_dimension failed for timeId in read_era5cds')
         start_time_index = prev_rec_size - n_last_steps + 1

         allocate(var1(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,n_last_steps))
         call LIS_verify(nf90_inq_varid(ftn,'lsp',rainfId), &
              'nf90_inq_varid failed for lsp in read_era5cds')
         call LIS_verify(nf90_get_var(ftn,rainfId, var1, &
              start=(/1,1,start_time_index/), &
              count=(/era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,n_last_steps/)),&
              'nf90_get_var failed for lsp in read_era5cds')
         call LIS_verify(nf90_get_att(ftn,rainfId, "_FillValue", missingValue), &
              'nf90_get_att failed in read_era5cds for lsp FillValue')

         ! lsp [m] -> mm/s
         do l=1,n_last_steps
           do j=1,era5cds_struc(n)%nrold
             do i=1,era5cds_struc(n)%ncold
                if(var1(i,j,l).ne.missingValue) then
                   var1(i,j,l) = var1(i,j,l)*1000.0/(1.0*60*60) 
                endif
             enddo
           enddo
         enddo

         do l=1,n_last_steps
           call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                  era5cds_struc(n)%nrold,var1(:,:,l),datalis)
           call interp_era5cds_var(n,findex,month,datalis,missingValue,.true.,&
                  era5cds_struc(n)%prev_rainf(:,l))
         enddo

         call LIS_verify(nf90_inq_varid(ftn,'cp',crainfId), &
              'nf90_inq_varid failed for cp in read_era5cds')
         call LIS_verify(nf90_get_var(ftn,crainfId, var1, &
              start=(/1,1,start_time_index/), &
              count=(/era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,n_last_steps/)),&
              'nf90_get_var failed for cp in read_era5cds')

         ! cp [m] -> mm/s
         do l=1,n_last_steps
           do j=1,era5cds_struc(n)%nrold
             do i=1,era5cds_struc(n)%ncold
                if(var1(i,j,l).ne.missingValue) then
                   var1(i,j,l) = var1(i,j,l)*1000.0/(1.0*60*60)
                endif
             enddo
           enddo
         enddo

         do l=1,n_last_steps
           call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                  era5cds_struc(n)%nrold,var1(:,:,l),datalis)
           call interp_era5cds_var(n,findex,month,datalis,missingValue,.true.,&
                  era5cds_struc(n)%prev_crainf(:,l))
         enddo

         !-----------------------------------------------------------------
         ! rain = lsp + cp [mm/s]
         !-----------------------------------------------------------------
         do l=1,n_last_steps
           do r=1,LIS_rc%lnr(n)
             do c=1,LIS_rc%lnc(n)
                k=c+(r-1)*LIS_rc%lnc(n)
                if(era5cds_struc(n)%prev_rainf(k,l).ne.LIS_rc%udef.and.&
                   era5cds_struc(n)%prev_crainf(k,l).ne.LIS_rc%udef) then
                   era5cds_struc(n)%prev_rainf(k,l) =  &
                                      era5cds_struc(n)%prev_rainf(k,l) + &
                                      era5cds_struc(n)%prev_crainf(k,l)
                endif
             enddo
           enddo
         enddo
         write(LIS_logunit,*) 'HKB prev rainf and crainf done...'

         call LIS_verify(nf90_inq_varid(ftn,'ssrd',swdId), &
               'nf90_inq_varid failed for ssrd in read_era5cds')
         call LIS_verify(nf90_get_var(ftn,swdId, var1, &
              start=(/1,1,start_time_index/), &
              count=(/era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,n_last_steps/)),&
               'nf90_get_var failed for swd in read_era5cds')
         call LIS_verify(nf90_get_att(ftn,swdId, "_FillValue", missingValue), &
             'nf90_get_att failed in read_era5cds for swd FillValue')

         ! swd - convert accumulated field to rate
         do l=1,n_last_steps
           do j=1,era5cds_struc(n)%nrold
             do i=1,era5cds_struc(n)%ncold
                if(var1(i,j,l).ne.missingValue) then
                   var1(i,j,l) = var1(i,j,l)/(1.0*60*60)
                endif
             enddo
           enddo
         enddo
         do l=1,n_last_steps
           call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                  era5cds_struc(n)%nrold,var1(:,:,l),datalis)
           call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                  era5cds_struc(n)%prev_swd(:,l))
         enddo

         call LIS_verify(nf90_inq_varid(ftn,'strd',lwdId), &
              'nf90_inq_varid failed for strd in read_era5cds')
         call LIS_verify(nf90_get_var(ftn,lwdId, var1, &
              start=(/1,1,start_time_index/), &
              count=(/era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,n_last_steps/)),&
              'nf90_get_var failed for lwd in read_era5cds')
         call LIS_verify(nf90_get_att(ftn,lwdId, "_FillValue", missingValue), &
              'nf90_get_att failed in read_era5cds for lwd FillValue')

         ! lwd - convert accumulated field to rate
         do l=1,n_last_steps
           do j=1,era5cds_struc(n)%nrold
             do i=1,era5cds_struc(n)%ncold
                if(var1(i,j,l).ne.missingValue) then
                   var1(i,j,l) = var1(i,j,l)/(1.0*60*60)
                endif
             enddo
           enddo
         enddo
         do l=1,n_last_steps
           call era5grid_2_lisgrid(era5cds_struc(n)%ncold, &
                  era5cds_struc(n)%nrold,var1(:,:,l),datalis)
           call interp_era5cds_var(n,findex,month,datalis,missingValue,.false.,&
                  era5cds_struc(n)%prev_lwd(:,l))
         enddo
         deallocate(var1)

         call LIS_verify(nf90_close(ftn), &
              'failed to close avg file in read_era5cds')
      else  ! prevavgfile exists
         ! at start of ERA5 in Jan 1940, no data for 0-6z 
         write(LIS_logunit,*) '[INFO] ',trim(prevavgfile)//' not available!'
         do l = 1, n_last_steps
            ll = rec_size - n_last_steps + l
            era5cds_struc(n)%prev_rainf(:,l)  = LIS_rc%udef
            era5cds_struc(n)%prev_crainf(:,l) = LIS_rc%udef
            era5cds_struc(n)%prev_swd(:,l)    = LIS_rc%udef
            era5cds_struc(n)%prev_lwd(:,l)    = LIS_rc%udef
         enddo
      endif

     deallocate(datalis)
#endif
  endif !if(read_flag)

  tindex = (day - 1)*24 + hour + 1
  atindex = tindex - n_last_steps

  call assign_processed_era5cdsf(n,kk,order,1,era5cds_struc(n)%tair(:,tindex))
  call assign_processed_era5cdsf(n,kk,order,2,era5cds_struc(n)%qair(:,tindex))
  call assign_processed_era5cdsf(n,kk,order,5,era5cds_struc(n)%uwind(:,tindex))
  call assign_processed_era5cdsf(n,kk,order,6,era5cds_struc(n)%vwind(:,tindex))
  call assign_processed_era5cdsf(n,kk,order,7,era5cds_struc(n)%ps(:,tindex))
  if ( atindex .le. 0 ) then
    atindex = atindex + n_last_steps
    call assign_processed_era5cdsf(n,kk,order,3,era5cds_struc(n)%prev_swd(:,atindex))
    call assign_processed_era5cdsf(n,kk,order,4,era5cds_struc(n)%prev_lwd(:,atindex))
    call assign_processed_era5cdsf(n,kk,order,8,era5cds_struc(n)%prev_rainf(:,atindex))
    call assign_processed_era5cdsf(n,kk,order,9,era5cds_struc(n)%prev_crainf(:,atindex))
  else
    call assign_processed_era5cdsf(n,kk,order,3,era5cds_struc(n)%swd(:,atindex))
    call assign_processed_era5cdsf(n,kk,order,4,era5cds_struc(n)%lwd(:,atindex))
    call assign_processed_era5cdsf(n,kk,order,8,era5cds_struc(n)%rainf(:,atindex))
    call assign_processed_era5cdsf(n,kk,order,9,era5cds_struc(n)%crainf(:,atindex))
  endif
  if ( .not. read_flag .and. ferror == 0 ) ferror = 1
  write(LIS_logunit,*) 'tindex=',tindex,' atindex=',atindex, read_flag, order

end subroutine read_era5cds

!BOP
! 
! !ROUTINE: interp_era5cds_var
! \label{interp_era5cds_var}
! 
! !INTERFACE: 
subroutine interp_era5cds_var(n,findex, month, input_var, missingValue, &
     pcp_flag, output_var)

! !USES: 
  use LIS_coreMod
  use LIS_logMod
  use LIS_spatialDownscalingMod
  use era5cds_forcingMod, only : era5cds_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: month
  real,    intent(in)    :: input_var(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold)
  real,    intent(in)    :: missingValue
  logical, intent(in)    :: pcp_flag
  real,    intent(out)   :: output_var(LIS_rc%lnc(n)*LIS_rc%lnr(n))

!
! !DESCRIPTION: 
!  This subroutine spatially interpolates a ERA5CDS field
!  to the LIS running domain
! 
!EOP

  integer   :: t,c,r,k,iret
  integer   :: doy
  integer   :: ftn
  real      :: f(era5cds_struc(n)%ncold*era5cds_struc(n)%nrold)
  logical*1 :: lb(era5cds_struc(n)%mi)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer   :: input_size
  integer   :: input_nc, input_nr
  integer   :: count1,nrec
! _____________________________________________________________

  input_size = era5cds_struc(n)%mi
  output_var = LIS_rc%udef
  lo = .true.
  input_nc = era5cds_struc(n)%ncold
  input_nr = era5cds_struc(n)%nrold
  
!-----------------------------------------------------------------------
! Apply corrections
!-----------------------------------------------------------------------

  lb = .false.
  do r=1,era5cds_struc(n)%nrold
     do c=1,era5cds_struc(n)%ncold
        k= c+(r-1)*era5cds_struc(n)%ncold
        f(k) = input_var(c,r)
        if(f(k) .ne. missingValue) then
           lb(k) = .true.
        endif
     enddo
  enddo
!-----------------------------------------------------------------------    
! Apply downscaling
!-----------------------------------------------------------------------    
     
  if(pcp_flag.and.LIS_rc%pcp_downscale(findex).ne.0) then 
     call LIS_generatePcpClimoRatioField(n,findex,"ERA5CDS",&
          month, & 
          input_size, &
          f, &
          lb)     
  endif
          
  if(pcp_flag.and.&
       trim(era5cds_struc(n)%met_interp).eq."budget-bilinear") then 
     
     call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var, &
          era5cds_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),& 
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          era5cds_struc(n)%w112,era5cds_struc(n)%w122,&
          era5cds_struc(n)%w212,era5cds_struc(n)%w222,&
          era5cds_struc(n)%n112,era5cds_struc(n)%n122,&
          era5cds_struc(n)%n212,era5cds_struc(n)%n222,&
          LIS_rc%udef, iret)
     
  elseif(trim(era5cds_struc(n)%met_interp).eq."bilinear".or.&
       trim(era5cds_struc(n)%met_interp).eq."budget-bilinear") then 

     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var, &
          era5cds_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n), & 
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          era5cds_struc(n)%w111,era5cds_struc(n)%w121,&
          era5cds_struc(n)%w211,era5cds_struc(n)%w221,&
          era5cds_struc(n)%n111,era5cds_struc(n)%n121,&
          era5cds_struc(n)%n211,era5cds_struc(n)%n221,&
          LIS_rc%udef, iret)
     
  elseif(trim(era5cds_struc(n)%met_interp).eq."neighbor") then 
     call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var,era5cds_struc(n)%mi,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          era5cds_struc(n)%n113,LIS_rc%udef,iret)

  elseif(trim(era5cds_struc(n)%met_interp).eq."average") then 
     call upscaleByAveraging(era5cds_struc(n)%mi, &
          LIS_rc%lnc(n)*LIS_rc%lnr(n), LIS_rc%udef, &
          era5cds_struc(n)%n111, lb, f, lo, output_var)

  elseif(trim(era5cds_struc(n)%met_interp).eq."none") then 
     ! at 0.25 degree, no interpolation needed. input_var(input_nc*input_rc)
     ! If input grid matches the output grid size:
       if( era5cds_struc(n)%subset_nc == input_nc .and. &
           era5cds_struc(n)%subset_nr == input_nr ) then
           output_var = f
     ! Otherwise, subset for same domain resolutions:
       else
         count1 = 0
         do r = 1, era5cds_struc(n)%subset_nr
           do c = 1, era5cds_struc(n)%subset_nc
              count1 = count1 + 1
              nrec = (era5cds_struc(n)%lat_line(c,r)-1)*input_nc + era5cds_struc(n)%lon_line(c,r)
              output_var(count1) = f(nrec)
           enddo ! columns
         enddo ! rows
       endif

  else
     write(LIS_logunit,*) '[ERR] Spatial interpolation option '//&
          trim(LIS_rc%met_interp(findex))//&
          ' not supported for ERA5CDS'
     call LIS_endrun()
  endif
  
  if( pcp_flag.and.LIS_rc%pcp_downscale(findex).ne.0 ) then 
     
     call LIS_pcpClimoDownscaling(n, findex, month,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n), output_var(:), lo)
     
  endif

end subroutine interp_era5cds_var
!BOP
!
! !ROUTINE: era5grid_2_lisgrid
!
! !DESCRIPTION:
! Changes grid_data from ECMWF data convention to GLDAS convention
!
! ERA5: North-to-South around Greenwich Meridian
! Global grid. Data are written as flat binary from "upper left to
! lower right" starting at 0.25-degree grid point center coordinates:
! 0.125E,89.875N and going to 0.125W,89.875S.
!
! LIS: South-to-North around Date Line
! Full global grid.  Starts at the southernmost latitude and date line,
! going east and then north.
!
! !REVISION HISTORY:
!  10 Apr 2002: Urszula Jambor;  Code adapted from
!               ecmwfgrid_2_grid2catgrid, by R. Reichle
!  21 Apr 2025: Hiroko Beaudoing; modified berggrid_2_gldasgrid.F90
!                                 to work for odd ny.
!
! !INTERFACE:
subroutine era5grid_2_lisgrid( nx, ny, grid_data, out_data )
!EOP
  implicit none

  integer, intent(in) :: nx, ny
  real, intent(in), dimension(nx,ny) :: grid_data
  real, intent(out), dimension(nx*ny):: out_data
  real, dimension(nx,ny):: lis_data

  integer :: i, j, m, n, c
  real :: tmp, tmp_data1(nx)

  ! ------------------------------------------------------------------
  ! some checks

  if ((nx /= 1440) .or. (ny /= 720)) then
     write (*,*) 'era5grid_2_gldasgrid(): This routine has only been'
     write (*,*) 'checked for nx=1440 and ny=720. STOPPING'
     stop
  end if
  if ((mod(nx,2) /= 0)) then
     write (*,*) 'era5grid_2_gldasgrid(): This routine can only work'
     write (*,*) 'for even nx. STOPPING.'
     stop
  end if

  !-------------------------------------------------------------------

  do j=1,ny

     ! swap latitude bands (North-to-South becomes South-to-North)
     n = ny-j+1
     tmp_data1      = grid_data(:,j)

     do i=1,nx/2

        ! shift longitudes (wrapping around Greenwhich Meridian becomes
        !  wrapping around Date Line)
        m = i + nx/2
        tmp          = tmp_data1(i)
        tmp_data1(i) = tmp_data1(m)
        tmp_data1(m) = tmp

     end do

     lis_data(:,n) = tmp_data1

  end do

  ! return output in 1D
  c=0
  do j=1,ny
     do i=1,nx
        c=c+1
        out_data(c) = lis_data(i,j)
     end do
  end do

end subroutine era5grid_2_lisgrid

!BOP
!
! !ROUTINE: assign_processed_era5cdsf
! \label{assign_processed_era5cdsf}
!
! !INTERFACE:
subroutine assign_processed_era5cdsf(n,kk,order,var_index,era5forc)
! !USES:
  use LIS_coreMod
  use era5cds_forcingMod, only : era5cds_struc
!
! !DESCRIPTION:
!  This routine assigns the interpolated ERA5 forcing data
!  to the module data structures to be used later for
!  time interpolation
!
!EOP
  implicit none

  integer :: n
  integer :: kk
  integer :: order
  integer :: var_index
  real    :: era5forc(LIS_rc%lnc(n)*LIS_rc%lnr(n))


  integer :: c,r

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then
           if(order.eq.1) then
              era5cds_struc(n)%metdata1(kk,var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   era5forc(c+(r-1)*LIS_rc%lnc(n))
           elseif(order.eq.2) then
              era5cds_struc(n)%metdata2(kk,var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   era5forc(c+(r-1)*LIS_rc%lnc(n))
           endif
        endif
     enddo
  enddo
end subroutine assign_processed_era5cdsf
