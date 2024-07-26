!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! ROUTINE: USAFsfcval
!
! REVISION HISTORY: 
! 22 Jun 2017  Initial version based on AGRMET_sfcval..Eric Kemp/SSAI/NASA
! 07 Sat 2018  Renamed to USAF_sfcval..................Eric Kemp/SSAI/NASA
!
! DESCRIPTION:
!  
! Vertically interpolate input fldbld fields (isobaric and near-surface)
! to LIS model surface.  The fldbld fields come from GRIB (GFS or GALWEM).
! Using the near-surface fields is a technique borrowed from the WRF REAL 
! preprocessor, and it greatly reduces bias in high terrain.
!-------------------------------------------------------------------------

subroutine USAF_sfcval( alt, &
     minwnd, sfcprs, sfcrlh, &
     sfcspd, sfctmp, imax, jmax, &
     kprs, prslvls, kprs1,tmp, hgt, rlh, &
     tmp_sfc, hgt_sfc, rlh_sfc, wndspd, pres)

   ! Imports

   ! Defaults
   implicit none

   ! Arguments
   integer,  intent(in)       :: kprs
   integer,  intent(in)       :: kprs1
   integer,  intent(in)       :: imax
   integer,  intent(in)       :: jmax
   real,     intent(in)       :: alt  ( imax, jmax )
   real,     intent(in)       :: hgt   ( imax, jmax, kprs1 )
   real,     intent(in)       :: rlh    ( imax, jmax, kprs1 )
   real,     intent(in)       :: tmp   ( imax, jmax, kprs1 )
   real,     intent(in)       :: hgt_sfc   ( imax, jmax )
   real,     intent(in)       :: rlh_sfc    ( imax, jmax )
   real,     intent(in)       :: tmp_sfc   ( imax, jmax )   
   real,     intent(in)       :: wndspd  ( imax,jmax)
   real,     intent(in)       :: pres  ( imax,jmax)
   integer,  intent(in)       :: prslvls   ( 30 )
   real,     intent(in)       :: minwnd
   real,     intent(out)      :: sfcprs    ( imax,jmax )
   real,     intent(out)      :: sfcrlh    ( imax,jmax )
   real,     intent(out)      :: sfcspd    ( imax,jmax )
   real,     intent(out)      :: sfctmp    ( imax,jmax )

   ! Local variables
   integer                       :: i
   integer                       :: j 
   integer                       :: kk   
   real                          :: aalt
   real                          :: abas
   real                          :: athick
   real                          :: bashgt
   real                          :: basp  
   real                          :: basrh
   real                          :: bastmp
   real                          :: calprs
   real                          :: h1
   real                          :: linint
   real                          :: p 
   real                          :: t1
   real                          :: t2
   real                          :: tbas
   real                          :: thick 
   real                          :: toprh
   real                          :: toptmp
   real                          :: ttop
   real                          :: za
   real                          :: zb
   real       :: alt_t  ( imax, jmax )
   logical :: found

   !------------------------------------------------------------------
   ! Define internal functions h1 and calprs for calculation of
   ! final surface pressure using the hypsometric equation.  Also 
   ! define internal function linint for vertically interpolating
   ! from isobaric levels to the model surface. 
   !------------------------------------------------------------------

   h1(t1, t2) = 287.04 * ((t1 + t2) * 0.5) / 9.80665
   
   calprs(t1, t2, za, zb, p) = max (min  &
        (p * exp(-((za - zb) /  &
        (h1 (t1, t2)))) * 100.0, &
        110000.0), 40000.0)
   
   linint(ttop, tbas, athick, aalt, abas) = tbas + ((ttop - tbas) *   &
        (aalt - abas) / athick)
   
   !------------------------------------------------------------------
   ! Loop through each model grid point.
   !------------------------------------------------------------------

   alt_t = alt
   JGRID : do j = 1, jmax
      IGRID : do i = 1, imax

         ! Sanity check to handle missing LIS terrain
         if(alt_t(i,j).eq.-9999.0) alt_t(i,j) = 0.0

         if (alt_t(i,j) < hgt_sfc(i,j)) then
            ! Case 1: LIS elevation is below GRIB ground surface.
            ! Use the GRIB surface values for RH and sfcprs as-is, and
            ! use standard -6.5 K/km lapse rate to extrapolate temperature.

            sfctmp(i,j) = tmp_sfc(i,j) + &
                 (0.0065 * (hgt_sfc(i,j) - alt_t(i,j)))
            sfctmp(i,j) = max( min( sfctmp(i,j), 350.0 ), 200.0 )
            sfcprs(i,j) = pres(i,j)
            sfcrlh(i,j) = rlh_sfc(i,j)
            sfcrlh(i,j) = max( min( sfcrlh(i,j), 1.0 ), 0.1 )

         else if (alt_t(i,j) == hgt_sfc(i,j)) then
            ! Case 2: LIS elevation is exactly on GRIB ground surface.
            ! Just use the GRIB (near) surface values.

            sfctmp(i,j) = tmp_sfc(i,j)
            sfctmp(i,j) = max( min( sfctmp(i,j), 350.0 ), 200.0 )
            sfcprs(i,j) = pres(i,j)
            sfcrlh(i,j) = rlh_sfc(i,j)
            sfcrlh(i,j) = max( min( sfcrlh(i,j), 1.0 ), 0.1 )
            
         else
            ! Case 3: LIS elevation is above GRIB ground surface.
            ! Interpolate between GRIB surface and first isobaric level above
            ! the LIS elevation, and use hypsometric equation for surface
            ! pressure.

            found = .false.

            SEARCH : do kk = 1, kprs

               if ( alt_t(i,j) .le.  hgt(i,j,kk) ) then
                  
                  basp   = pres(i,j) * 0.01 ! Pa to hPa
                  bastmp = tmp_sfc(i,j)
                  toptmp = tmp(i,j,kk)
                  bashgt = hgt_sfc(i,j)
                  thick  = hgt(i,j,kk) - hgt_sfc(i,j)

                  sfctmp(i,j) = linint( toptmp, bastmp, thick, &
                       alt_t(i,j), bashgt )
                  sfctmp(i,j) = max( min( sfctmp(i,j), 350.0 ), 200.0 )
                    
                  sfcprs(i,j) = calprs( bastmp, toptmp, alt_t(i,j),  &
                       bashgt, basp )

                  basrh       = rlh_sfc(i,j)
                  toprh       = rlh(i,j,kk)
                  sfcrlh(i,j) = linint( toprh, basrh, thick,  &
                       alt_t(i,j), bashgt )
                  sfcrlh(i,j) = max( min( sfcrlh(i,j), 1.0 ), 0.1 )
                  
                  found = .true.
                  exit ! Get out of kk loop
               end if
            end do SEARCH ! kk
               
            ! Contingency:  If LIS terrain is above all isobaric data, just
            ! use the highest level.
            if (.not. found) then
               sfctmp(i,j) = tmp(i,j,kprs)
               sfctmp(i,j) = max( min( sfctmp(i,j), 350.0 ), 200.0 )
               sfcprs(i,j) = float(prslvls(kprs)) * 100. ! hPa to Pa
               sfcrlh(i,j) = rlh(i,j,kprs)
               sfcrlh(i,j) = max( min( sfcrlh(i,j), 1.0 ), 0.1 )
            end if
         end if
           
         !------------------------------------------------------------------
         ! Set the output surface wind speed to the input resultant
         ! wind speed.  Limit the values between a max and min.
         !------------------------------------------------------------------
           
         sfcspd(i,j) = wndspd(i,j)
         sfcspd(i,j) = max( min( sfcspd(i,j), 75.0 ), minwnd )           

      enddo IGRID
   enddo JGRID

end subroutine USAF_sfcval
