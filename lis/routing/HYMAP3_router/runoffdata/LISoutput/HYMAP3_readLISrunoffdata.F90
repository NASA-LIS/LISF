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
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  7 Jan 2016: Sujay Kumar, Initial implementation
! 17 Mar 2016: Augusto Getirana, changes in input file name generation
! and surface runoff and baseflow variables - this will reduce the
! number of times input files are read

subroutine HYMAP3_readLISrunoffdata(n, surface_runoff, baseflow)

! !USES:
  use HYMAP3_LISrunoffdataMod
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_coreMod
  use LIS_fileIOMod
  use LIS_logMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  integer,  intent(in) :: n
  real, intent(out)    :: surface_runoff(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real, intent(out)    :: baseflow(LIS_rc%lnc(n),LIS_rc%lnr(n))

  !caveats
  ! 1) assumes the LIS outputs are the same output interval as that of
  ! the HYMAP model timestep. No temporal aggregation is done.
  !
  ! 3) Assumes that the units of Qs and Qsb are in kg/m2s
  !
  ! 4) Assumes that the LIS outputs are in the same grid/map projection/
  ! resolution.
  !
  ! 5) LIS outputs are in NETCDF format.
  !
  !Added total evapotranspiration (Evap)
  integer                       :: c,r
  !ag - 17Mar2016
  real                  :: qs2d(HYMAP3_LISrunoffdata_struc(n)%nc, &
       HYMAP3_LISrunoffdata_struc(n)%nr)
  real                  :: qsb2d(HYMAP3_LISrunoffdata_struc(n)%nc, &
       HYMAP3_LISrunoffdata_struc(n)%nr)
  logical*1             :: lb(&
       HYMAP3_LISrunoffdata_struc(n)%nc*HYMAP3_LISrunoffdata_struc(n)%nr)
  real                  :: var_input( &
       HYMAP3_LISrunoffdata_struc(n)%nc*HYMAP3_LISrunoffdata_struc(n)%nr)
  logical*1             :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                  :: var_out(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer               :: ios, nid,qsid,qsbid !,evapid
  character(LIS_CONST_PATH_LEN) :: filename
  logical               :: file_exists
  logical               :: check_Flag

  external :: neighbor_interp
  external :: upscaleByAveraging

  !create LIS filename
  call LIS_create_output_filename(n, &
       filename, "netcdf", check_flag, &
       model_name='SURFACEMODEL', &
       odir=HYMAP3_LISrunoffdata_struc(n)%odir,&
       writeint=HYMAP3_LISrunoffdata_struc(n)%outInterval)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  if(trim(HYMAP3_LISrunoffdata_struc(n)%previous_filename) /= &
       trim(filename))then

     HYMAP3_LISrunoffdata_struc(n)%previous_filename=filename

     inquire(file=filename, exist=file_exists)
     if(file_exists) then
        write(LIS_logunit,*) '[INFO] Reading '//trim(filename)
        !ag - 17Mar2016
        ios = nf90_open(path=filename,&
             mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios, &
             'Error opening file in HYMAP3_readLISrunoffdata')

        ios = nf90_inq_varid(nid,'Qs_tavg',qsid)
        call LIS_verify(ios, &
             'failed to read Qs_tavg field in HYMAP3_readLISrunoffdata')

        ios = nf90_inq_varid(nid,'Qsb_tavg',qsbid)
        call LIS_verify(ios, &
             'failed to read Qsb_tavg field in HYMAP3_readLISrunoffdata')

        if(HYMAP3_LISrunoffdata_struc(n)%domainCheck) then
           if(LIS_rc%wopt.eq."2d gridspace") then
              ios = nf90_get_var(nid,qsid, &
                   HYMAP3_LISrunoffdata_struc(n)%qs, &
                   start=(/LIS_ews_halo_ind(n,LIS_localPet+1), &
                   LIS_nss_halo_ind(n,LIS_localPet+1)/),&
                   count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
              call LIS_verify(ios, &
                   'failed to read Qs_tavg field in HYMAP3_readLISrunoffdata')

              ios = nf90_get_var(nid,qsbid, &
                   HYMAP3_LISrunoffdata_struc(n)%qsb,&
                   start=(/LIS_ews_halo_ind(n,LIS_localPet+1), &
                   LIS_nss_halo_ind(n,LIS_localPet+1)/),&
                   count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
              call LIS_verify(ios, &
                   'failed to read Qsb_tavg field in HYMAP3_readLISrunoffdata')

           else
              write(LIS_logunit,*) &
                   "[ERR] Stand-alone HYMAP3 is only supported for '2d gridspace' outputs currently"
              call LIS_endrun()
           endif

           call LIS_verify(nf90_close(nid))

        else
           if(LIS_rc%wopt.eq."2d gridspace") then

              ios = nf90_get_var(nid,qsid,qs2d)
              call LIS_verify(ios, &
                   'failed to read Qs_tavg field in HYMAP3_readLISrunoffdata')

              ios = nf90_get_var(nid,qsbid,qsb2d)
              call LIS_verify(ios, &
                   'failed to read Qsb_tavg field in HYMAP3_readLISrunoffdata')

              if(LIS_isAtAfinerResolution(n, &
                   HYMAP3_LISrunoffdata_struc(n)%datares)) then

                 lb = .true.
                 do r=1,HYMAP3_LISrunoffdata_struc(n)%nr
                    do c=1,HYMAP3_LISrunoffdata_struc(n)%nc
                       var_input(c+(r-1)* &
                            HYMAP3_LISrunoffdata_struc(n)%nc) = qs2d(c,r)
                       if(qs2d(c,r).lt.0) then
                          lb(c+(r-1)*HYMAP3_LISrunoffdata_struc(n)%nc) = &
                               .false.
                       endif
                    enddo
                 enddo

                 call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
                      lo,var_out, &
                      HYMAP3_LISrunoffdata_struc(n)%nc* &
                      HYMAP3_LISrunoffdata_struc(n)%nr,&
                      LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
                      LIS_domain(n)%lat, LIS_domain(n)%lon,  &
                      HYMAP3_LISrunoffdata_struc(n)%n11,     &
                      LIS_rc%udef,ios)

                 do r=1,LIS_rc%lnr(n)
                    do c=1,LIS_rc%lnc(n)
                       HYMAP3_LISrunoffdata_struc(n)%qs(c,r) = &
                            var_out(c+(r-1)*LIS_rc%lnc(n))
                    enddo
                 enddo

                 lb = .true.
                 do r=1,HYMAP3_LISrunoffdata_struc(n)%nr
                    do c=1,HYMAP3_LISrunoffdata_struc(n)%nc
                       var_input(c+(r-1)* &
                            HYMAP3_LISrunoffdata_struc(n)%nc) = qsb2d(c,r)
                       if(qsb2d(c,r).lt.0) then
                          lb(c+(r-1)* &
                               HYMAP3_LISrunoffdata_struc(n)%nc) = .false.
                       endif
                    enddo
                 enddo

                 call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
                      lo,var_out, &
                      HYMAP3_LISrunoffdata_struc(n)%nc* &
                      HYMAP3_LISrunoffdata_struc(n)%nr,&
                      LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
                      LIS_domain(n)%lat, LIS_domain(n)%lon,  &
                      HYMAP3_LISrunoffdata_struc(n)%n11,     &
                      LIS_rc%udef,ios)

                 do r=1,LIS_rc%lnr(n)
                    do c=1,LIS_rc%lnc(n)
                       HYMAP3_LISrunoffdata_struc(n)%qsb(c,r) = &
                            var_out(c+(r-1)*LIS_rc%lnc(n))
                    enddo
                 enddo

              else

                 lb = .true.
                 do r=1,HYMAP3_LISrunoffdata_struc(n)%nr
                    do c=1,HYMAP3_LISrunoffdata_struc(n)%nc
                       var_input(c+(r-1)* &
                            HYMAP3_LISrunoffdata_struc(n)%nc) = qs2d(c,r)
                       if(qs2d(c,r).lt.0) then
                          lb(c+(r-1)* &
                               HYMAP3_LISrunoffdata_struc(n)%nc) = .false.
                       endif
                    enddo
                 enddo

                 call upscaleByAveraging(&
                      HYMAP3_LISrunoffdata_struc(n)%nc* &
                      HYMAP3_LISrunoffdata_struc(n)%nr, &
                      LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                      LIS_rc%udef, &
                      HYMAP3_LISrunoffdata_struc(n)%n11, lb, &
                      var_input, lo, var_out)

                 do r=1,LIS_rc%lnr(n)
                    do c=1,LIS_rc%lnc(n)
                       HYMAP3_LISrunoffdata_struc(n)%qs(c,r) = &
                            var_out(c+(r-1)*LIS_rc%lnc(n))
                    enddo
                 enddo

                 lb = .true.
                 do r=1,HYMAP3_LISrunoffdata_struc(n)%nr
                    do c=1,HYMAP3_LISrunoffdata_struc(n)%nc
                       var_input(c+(r-1)* &
                            HYMAP3_LISrunoffdata_struc(n)%nc) = qsb2d(c,r)
                       if(qsb2d(c,r).lt.0) then
                          lb(c+(r-1)*HYMAP3_LISrunoffdata_struc(n)%nc) = &
                               .false.
                       endif
                    enddo
                 enddo

                 call upscaleByAveraging(&
                      HYMAP3_LISrunoffdata_struc(n)%nc* &
                      HYMAP3_LISrunoffdata_struc(n)%nr, &
                      LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                      LIS_rc%udef, &
                      HYMAP3_LISrunoffdata_struc(n)%n11, lb, &
                      var_input, lo, var_out)

                 do r=1,LIS_rc%lnr(n)
                    do c=1,LIS_rc%lnc(n)
                       HYMAP3_LISrunoffdata_struc(n)%qsb(c,r) = &
                            var_out(c+(r-1)*LIS_rc%lnc(n))
                    enddo
                 enddo

              endif
           else
              write(LIS_logunit,*) &
                   "[ERR] Stand-alone HYMAP3 is only supported for '2d gridspace' outputs currently"
              call LIS_endrun()
           endif

           call LIS_verify(nf90_close(nid))
        endif
     else
        write(LIS_logunit,*) '[ERR] Failed to find '//trim(filename)
        call LIS_endrun()
     endif
  endif

#endif

  !ag - 17Mar2016
  where(HYMAP3_LISrunoffdata_struc(n)%qs/=LIS_rc%udef)
     surface_runoff = HYMAP3_LISrunoffdata_struc(n)%qs
     baseflow = HYMAP3_LISrunoffdata_struc(n)%qsb
  else where
     surface_runoff = 0.0
     baseflow = 0.0
  end where

end subroutine HYMAP3_readLISrunoffdata
