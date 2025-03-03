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
! !ROUTINE: agrmet_cdfs_cod_filename
! \label{agrmet_cdfs_cod_filename}
!
! !INTERFACE:
subroutine agrmet_cdfs_cod_filename(fname,rootdir,dir,&
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
!  This routines generates the name of the CDFS2 cloud optical depth file
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

  character(len=5), parameter :: fhemi(2) = (/'N-HEM','S-HEM'/)
  character(len=8)            :: ftime1
  character(len=10)           :: ftime2
  character(len=2)            :: fhr
  character(len=41)           :: fname1
  character(len=13)           :: fname2

! PS.AFWA_SC.U_DI.A_GP.WWMCA-M_GR.P24KM_AR.N-HEM_PA.WWMCA_CY.00_FH.000_20130601000000

  fname1 = 'PS.AFWA_SC.U_DI.A_GP.WWMCA-M_GR.P24KM_AR.'
  fname2 = '_PA.WWMCA_CY.'

  write (UNIT=fhr, FMT='(i2.2)') hr
  write (UNIT=ftime2, fmt='(i4, i2.2, i2.2, i2.2)') yr,mo,da,hr

  if ( use_timestamp == 1 ) then
     write (UNIT=ftime1, fmt='(i4, i2.2, i2.2)') yr,mo,da
     fname = trim(rootdir) // '/' // ftime1 // '/' // trim(dir) // '/' // fname1 // &
             fhemi(hemi) // fname2 // fhr // '_FH.000_' // ftime2 // '0000'
  else
     fname = trim(rootdir) // '/' // trim(dir) // '/' // fname1 // &
             fhemi(hemi) // fname2 // fhr // '_FH.000_' // ftime2 // '0000'
  endif

end subroutine agrmet_cdfs_cod_filename
!BOP
!
! !ROUTINE: read_cod
! \label{read_cod}
!
! !REVISION HISTORY:
! 31 Mar 2016; James Geiger, Initial Code
!
! !INTERFACE:
subroutine read_cod(n, codfile, hr, cod, cloud_base, cloud_top, cloud_pcts, &
                    cloud_tot_pcts, cloud_times)
   ! !USES:
   use LIS_logMod,        only : LIS_logunit, LIS_verify, LIS_warning
   use AGRMET_forcingMod, only : agrmet_struc

#if (defined USE_GRIBAPI)
   use grib_api
#endif

   implicit none
   ! !ARGUMENTS:
   integer, intent(in)                       :: n
   character(len=*), intent(in)              :: codfile
   integer, intent(in)                       :: hr
   real, dimension(4,1024,1024), intent(out) :: cod
   real, dimension(4,1024,1024), intent(out) :: cloud_base
   real, dimension(4,1024,1024), intent(out) :: cloud_top
   real, dimension(4,1024,1024), intent(out) :: cloud_pcts
   real, dimension(1024,1024),   intent(out) :: cloud_tot_pcts
   real, dimension(1024,1024),   intent(out) :: cloud_times
!
! !DESCRIPTION:
!  This routine reads the CDFS II cloud optical depth.
!
!EOP
   integer, parameter     :: NV=18
   integer                :: ftn
   logical                :: file_exists
   integer                :: iv, iv_total
   integer                :: k
   integer                :: nvars
   integer                :: igrib
   integer                :: var_index
   logical                :: var_found
   logical                :: var_status(NV)
   real                   :: missingValue
   real,      allocatable :: f(:)
   integer                :: indparam(NV), level(NV), indlevel(NV)
   integer                :: datatime
   integer                :: indparam_val, level_val, indlevel_val, datatime_val
   integer                :: rc

#if ( defined USE_GRIBAPI )
   !--------------------------------------------------------------------------
   ! Set the GRIB parameter specifiers
   !--------------------------------------------------------------------------

   indparam = (/ 185, 185, 185, 185, 227, 227, 227, 227, 228, 228, 228, 228, 163, 163, 163, 163, 71, 183 /)
   indlevel = (/ 109, 109, 109, 109, 109 ,109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 200, 200 /)
   level    = (/ 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 0, 0 /)
   !datatime = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
   datatime = 100*hr

   iv_total = NV
   inquire (file=codfile, exist=file_exists)
   if (file_exists) then

      call grib_open_file(ftn,trim(codfile),'r',rc)
      if ( rc /= 0 ) then 
         write(LIS_logunit,*) 'Could not open file: ',trim(codfile)
         return
      endif

      call grib_count_in_file(ftn,nvars,rc)
      call LIS_verify(rc, 'error in grib_count_in_file in read_cod')

      allocate(f(1024*1024))

      var_status = .false.
      do k = 1, nvars
         call grib_new_from_file(ftn, igrib, rc)
         call LIS_warning(rc, 'error in grib_new_from_file in read_cod')
         if ( rc /= 0 ) then 
            write(LIS_logunit,*) &
               'Could not retrieve entries in file: ',trim(codfile)
            deallocate(f)
            return
         endif

         call grib_get(igrib,'indicatorOfParameter',indparam_val,rc)
         call LIS_verify(rc, &
                 'error in grib_get: indicatorOfParameter in read_cod')

         call grib_get(igrib,'level',level_val,rc)
         call LIS_verify(rc, 'error in grib_get: level in read_cod')

         call grib_get(igrib,'indicatorOfTypeOfLevel',indlevel_val,rc)
         call LIS_verify(rc, &
                 'error in grib_get: indicatorOfTypeOfLevel in read_cod')

         call grib_get(igrib,'dataTime',datatime_val,rc)
         call LIS_verify(rc, 'error in grib_get: dataTime in read_cod')

         var_found = .false.
         do iv = 1, iv_total
            if ( (indparam_val == indparam(iv)) .and. &
                 (level_val == level(iv))       .and. &
                 (indlevel_val == indlevel(iv)) .and. &
                 (datatime_val == datatime) ) then
               var_found = .true.
               var_index = iv
               var_status(iv) = .true. 
               exit
            endif
         enddo

         if ( var_found ) then
            call grib_get(igrib,'values',f,rc)
            call LIS_warning(rc, 'error in grib_get:values in read_cod')

            if ( rc /= 0 ) then
               write(LIS_logunit,*) &
                  'Could not retrieve entries in file: ',trim(codfile)
               write(LIS_logunit,*) 'for variables ',k
               deallocate(f)
               return
            endif

            call grib_get(igrib,'missingValue',missingValue,rc)
            call LIS_verify(rc, 'error in grib_get:missingValue in read_cod')

            ! layer cloud optical depth was multiplied by 10
            ! when written to disk.
            !
            ! layer cloud particle size was multiplied by 10
            ! when written to disk.
            if ( indparam_val == 185 ) then
               where ( f /= missingValue )
                  f = f / 10.0
               endwhere
               cod(level_val,:,:) = reshape(f, (/1024,1024/))
            elseif ( indparam_val == 227 ) then
               cloud_base(level_val,:,:) = reshape(f, (/1024,1024/))
            elseif ( indparam_val == 228 ) then
               cloud_top(level_val,:,:) = reshape(f, (/1024,1024/))
            elseif ( indparam_val == 163 ) then
               cloud_pcts(level_val,:,:) = reshape(f, (/1024,1024/))
            elseif ( indparam_val == 71 ) then
               cloud_tot_pcts(:,:) = reshape(f, (/1024,1024/))
            elseif ( indparam_val == 183 ) then
               where ( f /= missingValue )
                  f = f / 60.0 ! convert from minutes to hours
               endwhere
               cloud_times(:,:) = reshape(f, (/1024,1024/))
            endif
         endif

         call grib_release(igrib,rc)
         call LIS_verify(rc, 'error in grib_release in read_cod')
      enddo

      call grib_close_file(ftn)

      deallocate(f)

      do k = 1, iv_total
         if ( .not. var_status(k) ) then
            write(LIS_logunit,*) &
               'Could not retrieve entries in file: ',trim(codfile)
            return
         endif
      enddo
   else
      write(LIS_logunit,*) 'Could not find file: ',trim(codfile)
   endif
#endif
end subroutine read_cod

