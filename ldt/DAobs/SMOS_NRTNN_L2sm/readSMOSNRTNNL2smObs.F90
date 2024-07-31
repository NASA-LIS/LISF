!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: readSMOSNRTNNL2smObs
! \label{readSMOSNRTNNL2smObs}
!
! !REVISION HISTORY:
!  06 Jan. 2021: Yonghwan Kwon, Initial Specification
!  21 Feb. 2021: Mahdi Navari, Before the bugfix, code did
!                not work with optimization level -1 and -2
!  21 Feb. 2021: Mahdi Navari, code modified to write the DGG 
!                lookup table into a netCDF file
!  3  Apr. 2023: Mahdi Navari, Bugfix: issue in grid finding 
!                period (Issue #1301)
! 
! !INTERFACE:
subroutine readSMOSNRTNNL2smObs(n)
! !USES:

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif
   use ESMF
   use LDT_coreMod
   use LDT_logMod
   use LDT_constantsMod, only : LDT_CONST_PATH_LEN
   use LDT_timeMgrMod
   use LDT_DAobsDataMod
   use SMOSNRTNNL2sm_obsMod
   use map_utils

   implicit none
! !ARGUMENTS:
   integer, intent(in) :: n
!
! !DESCRIPTION:
!
! This subroutine provides the data reader for the SMOS NRT NN L2
! soil moisture retrieval product. 
!
!EOP

   real*8             :: timenow
   logical            :: alarmCheck
   logical            :: file_exists
   integer            :: c, r, i, j
   character(len=LDT_CONST_PATH_LEN)      :: fname
   character(len=LDT_CONST_PATH_LEN) :: nc_filename
   integer            :: mn_ind
   integer            :: yr, mo, da, hr, mn, ss
   integer            :: doy
   integer            :: ftn
   integer            :: ierr
   integer            :: ios
   integer            :: nid
   real               :: gmt
   character*8        :: yyyymmdd
   character*4        :: yyyy
   character*2        :: mm, dd, hh
   character*100      :: list_files
   character(len=LDT_CONST_PATH_LEN)      :: smos_filename(10)
   real               :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   integer            :: lat_varid, lon_varid, sm_varid, dim_ids(2)

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations.
!-----------------------------------------------------------------------

   SMOSNRTNNsmobs(n)%smobs = LDT_rc%udef
   smobs = LDT_rc%udef

   if (LDT_rc%ts .gt. 3600) then
      write (LDT_logunit, *) '[ERR] Please set the LDT timestep to 1hr or less'
      write (LDT_logunit, *) '[ERR] This is required for SMOS NRT NN L2 data processing'
      call LDT_endrun()
   endif

   write (yyyymmdd, '(i4.4,2i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da
   write (yyyy, '(i4.4)') LDT_rc%yr
   write (mm, '(i2.2)') LDT_rc%mo
   write (dd, '(i2.2)') LDT_rc%da
   write (hh, '(i2.2)') LDT_rc%hr


   if (SMOSNRTNNsmobs(n)%start_day.ne.LDT_rc%da) then
      SMOSNRTNNsmobs(n)%count_day = SMOSNRTNNsmobs(n)%count_day + 1
      SMOSNRTNNsmobs(n)%start_day = LDT_rc%da
   endif


   list_files = 'ls '//trim(SMOSNRTNNsmobs(n)%odir)// &
                '/'//trim(yyyy)//'.'//trim(mm)//'.'//dd// &
                '/W_XX-ESA,SMOS,NRTNN_C_LEMM_*_' &
                //trim(yyyymmdd)//trim(hh)//'*_'//trim(yyyy) &
                //"*.nc > SMOS_filelist.dat"

   call system(trim(list_files))

   i = 1
   ftn = LDT_getNextUnitNumber()
   open (ftn, file="./SMOS_filelist.dat", &
         status='old', iostat=ierr)

   do while (ierr .eq. 0)
      read (ftn, '(a)', iostat=ierr) fname
      if (ierr .ne. 0) then
         exit
      endif
      ss = 0

      smos_filename(i) = fname

      write (LDT_logunit, *) '[INFO] reading ', trim(smos_filename(i))

      call read_SMOSNRTL2sm_data(n, smos_filename(i), SMOSNRTNNsmobs(n)%smobs)

      i = i + 1
   enddo
   call LDT_releaseUnitNumber(ftn)

   call LDT_logSingleDAobs(n, LDT_DAobsData(n)%soilmoist_obs, &
                           SMOSNRTNNsmobs(n)%smobs, vlevel=1)

end subroutine readSMOSNRTNNL2smObs

!BOP
! 
! !ROUTINE: read_SMOSNRTL2sm_data
! \label{read_SMOSNRTL2sm_data}
!
! !INTERFACE:
subroutine read_SMOSNRTL2sm_data(n, fname, smobs_inp)
! 
! !USES:

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use SMOSNRTNNL2sm_obsMod
  use write_lookup_table

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer              :: n
  character (len=200)  :: fname
  real                 :: smobs_inp(LDT_rc%lnc(n),LDT_rc%lnr(n))
  !real*8               :: time

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!
!
!EOP
   integer              :: ios
   integer              :: nid, smid, dimid, lonid, latid, smunctid, dggid
   integer              :: dgg_length
   real, allocatable    :: sm(:), smunct(:)
   real, allocatable    :: lon_dgg(:), lat_dgg(:)
   integer, allocatable :: DGG_id_number(:)
   integer              :: max_dgg_id_number
   real                 :: smobs_sum(LDT_rc%lnc(n),LDT_rc%lnr(n))
   integer              :: count_smobs(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real                 :: lat2d(LDT_rc%lnc(n),LDT_rc%lnr(n)), lon2d(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real                 :: dx, dy
   real                 :: max_lon_dgg, min_lon_dgg, max_lat_dgg, min_lat_dgg
   character (len = 20) :: dimname 
   integer              :: c,r,i,j
   integer, allocatable :: i_dgg(:)
   real, allocatable    :: sm_dgg(:), smunct_dgg(:)
   integer              :: total_length
   real, allocatable    :: dgg_indices_1d(:)


#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
   call LDT_verify(ios,'Error opening file '//trim(fname))

   ios = nf90_inq_varid(nid, 'longitude',lonid)
   call LDT_verify(ios, 'Error nf90_inq_varid: longitude')

   ios = nf90_inq_varid(nid, 'latitude',latid)
   call LDT_verify(ios, 'Error nf90_inq_varid: latitude')

   ios = nf90_inq_varid(nid, 'soil_moisture',smid)
   call LDT_verify(ios, 'Error nf90_inq_varid: soil_moisture')

   ios = nf90_inq_varid(nid, 'soil_moisture_uncertainty',smunctid)
   call LDT_verify(ios, 'Error nf90_inq_varid: soil_moisture_uncertainty')

   ios = nf90_inq_varid(nid, 'DGG_id_number',dggid)
   call LDT_verify(ios, 'Error nf90_inq_varid: DGG_id_number')

   !dimension
   ios = nf90_inq_dimid(nid, 'DGG_id_number', dimid)
   call LDT_verify(ios, 'Error nf90_inq_varid: soil_moisture dimension id')

   ios = nf90_inquire_dimension(nid, dimid, dimname, dgg_length)
   call LDT_verify(ios, 'Error nf90_inq_varid: soil_moisture dimension')

   !values
   if (dgg_length > 0) then 
      allocate(sm(dgg_length))
      allocate(smunct(dgg_length))
      allocate(lon_dgg(dgg_length))
      allocate(lat_dgg(dgg_length)) 
      allocate(DGG_id_number(dgg_length))
      allocate(i_dgg(dgg_length))

      ios = nf90_get_var(nid, smid, sm)
      call LDT_verify(ios, 'Error nf90_get_var: soil_moisture')

      ios = nf90_get_var(nid, smunctid, smunct)
      call LDT_verify(ios, 'Error nf90_get_var: soil_moisture_uncertainty')

      ios = nf90_get_var(nid, lonid, lon_dgg)
      call LDT_verify(ios, 'Error nf90_get_var: longitude') 

      ios = nf90_get_var(nid, latid, lat_dgg)
      call LDT_verify(ios, 'Error nf90_get_var: latitude')       

      ios = nf90_get_var(nid, dggid, DGG_id_number)
      call LDT_verify(ios, 'Error nf90_get_var: DGG_id_number')

      ! find max/min lon_dgg/lat_dgg
      max_lon_dgg = maxval(lon_dgg)
      min_lon_dgg = minval(lon_dgg)
      max_lat_dgg = maxval(lat_dgg)
      min_lat_dgg = minval(lat_dgg)

      ! find min/max DGG_id_number in a file
      max_dgg_id_number = maxval(DGG_id_number)

      if (max_dgg_id_number > 0) then
         allocate(sm_dgg(max_dgg_id_number))
         allocate(smunct_dgg(max_dgg_id_number))

         sm_dgg = LDT_rc%udef
         smunct_dgg = LDT_rc%udef

         ! put sm to corresponding DGG_id_number
         sm_dgg(DGG_id_number(:)) = sm(:)
         smunct_dgg(DGG_id_number(:)) = smunct(:)

         ! Put the SMOS soil moisture data in LIS input grid
         smobs_sum = 0 !initialize
         count_smobs = 0 !initialize
         i_dgg = (/(j, j = 1,dgg_length) /)

         do r=1,LDT_rc%lnr(n)
            do c=1,LDT_rc%lnc(n)
               lon2d(c,r) = LDT_domain(n)%lon(c+(r-1)*LDT_rc%lnc(n))
               lat2d(c,r) = LDT_domain(n)%lat(c+(r-1)*LDT_rc%lnc(n))

               !find dx and dy
               if (r==1.and.c==1) then
                  lon2d(c+1,r) = LDT_domain(n)%lon((c+1)+(r-1)*LDT_rc%lnc(n))
                  lat2d(c,r+1) = LDT_domain(n)%lat(c+((r+1)-1)*LDT_rc%lnc(n))
                  dx = lon2d(c+1,r) - lon2d(c,r)
                  dy = lat2d(c,r+1) - lat2d(c,r)
               endif
 
               if (lon2d(c,r)+dx/2 >= min_lon_dgg.and.&
                   lon2d(c,r)-dx/2 <= max_lon_dgg.and.&
                   lat2d(c,r)+dy/2 >= min_lat_dgg.and.&
                   lat2d(c,r)-dy/2 <= max_lat_dgg) then


                  if (SMOSNRTNNsmobs(n)%count_day <= 150) then
                     ! assume that during 150 days after the simulation start date
                     ! all land grids have assigned dgg_id_number
                     ! in order to avoid looping over ocean grids for every files
                     ! originally count_day was set to 30 days, but it seems that this 
                     ! period was not enough to assign all SMOS grids to the LIS grids.  
                     if (.not. SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_assign) then 
                        call find_SMOS_Dgg_id_number(n,c,r, lon_dgg, lat_dgg, &
                                                     dgg_length, i_dgg, DGG_id_number, &
                                                     min_lon_dgg, max_lon_dgg, &
                                                     min_lat_dgg, max_lat_dgg, &
                                                     dx, dy, lon2d(c,r), lat2d(c,r))
                     endif
                  endif

                  if (SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_assign) then
                     if (size(SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices) > 0) then
                        do i=1,size(SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices)
                        ! soil moisture uncertainty < 0.07 m3m-3
                    
                           if (SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices(i) <= max_dgg_id_number) then
                              if (smunct_dgg(SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices(i)) < 0.07) then
                                 if (sm_dgg(SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices(i)) >= 0) then
                                    smobs_sum(c,r) = smobs_sum(c,r) + sm_dgg(SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices(i))
                                    count_smobs(c,r) = count_smobs(c,r) + 1
                                 endif
                              endif
                           endif
                        enddo
                     endif
                  endif
               endif

            enddo
         enddo

         do r=1,LDT_rc%lnr(n)
            do c=1,LDT_rc%lnc(n)
               lon2d(c,r) = LDT_domain(n)%lon(c+(r-1)*LDT_rc%lnc(n))
               lat2d(c,r) = LDT_domain(n)%lat(c+(r-1)*LDT_rc%lnc(n))

               if (lon2d(c,r)+dx/2 >= min_lon_dgg.and.&
                   lon2d(c,r)-dx/2 <= max_lon_dgg.and.&
                   lat2d(c,r)+dy/2 >= min_lat_dgg.and.&
                   lat2d(c,r)-dy/2 <= max_lat_dgg) then

                  if (count_smobs(c,r) > 0) then
                     smobs_inp(c,r) = smobs_sum(c,r)/count_smobs(c,r)
  
                  else ! count_smobs(c,r) == 0
                     if (r > 1) then
                        if (count_smobs(c,r-1) > 0) then
                           smobs_sum(c,r) = smobs_sum(c,r) + smobs_sum(c,r-1)/count_smobs(c,r-1)
                           count_smobs(c,r) = count_smobs(c,r) + 1
                        endif
                     endif
                     if (c > 1) then
                        if (count_smobs(c-1,r) > 0) then
                           smobs_sum(c,r) = smobs_sum(c,r) + smobs_sum(c-1,r)/count_smobs(c-1,r)
                           count_smobs(c,r) = count_smobs(c,r) + 1
                        endif
                     endif
                     if (c < LDT_rc%lnc(n)) then
                        if (count_smobs(c+1,r) > 0) then
                           smobs_sum(c,r) = smobs_sum(c,r) + smobs_sum(c+1,r)/count_smobs(c+1,r)
                           count_smobs(c,r) = count_smobs(c,r) + 1
                        endif
                     endif
                     if (r < LDT_rc%lnr(n)) then
                        if (count_smobs(c,r+1) > 0) then
                           smobs_sum(c,r) = smobs_sum(c,r) + smobs_sum(c,r+1)/count_smobs(c,r+1)
                           count_smobs(c,r) = count_smobs(c,r) + 1
                        endif
                     endif

                     ! get average smobs only if the number of neighbor grids having soil moisutre values
                     ! is 3 or greater
                     if (count_smobs(c,r) >= 3) then
                        smobs_inp(c,r) = smobs_sum(c,r)/count_smobs(c,r)
    
                     endif
                  endif
               endif
            enddo
         enddo

         deallocate(sm_dgg)
         deallocate(smunct_dgg)

      endif !max_dgg_id_number > 0

      deallocate(sm)
      deallocate(lon_dgg)
      deallocate(lat_dgg)
      deallocate(smunct)
      deallocate(DGG_id_number)
      deallocate(i_dgg)
   endif

   ios = nf90_close(ncid=nid)
   call LDT_verify(ios,'Error closing file '//trim(fname))


  ! assume that during 150 days after the simulation start date
  ! all land grids have assigned dgg_id_number

   if (SMOSNRTNNsmobs(n)%count_day .eq. 151) then

     total_length = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             total_length = total_length + 1
             if (SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_assign) then
               total_length = total_length + size(SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices)
             end if
          enddo
       enddo
       ! Allocate memory and store
       allocate(dgg_indices_1d(total_length))
       dgg_indices_1d = 0 ! Initialize memory
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             if (SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_assign) then
                dgg_indices_1d(i) = size(SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices)
                do j = 1, size(SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices)
                   i = i + 1
                   dgg_indices_1d(i) = SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices(j)
                end do
             end if
          enddo
       enddo
    
    SMOSNRTNNsmobs(n)%dgg_lookup_1d = dgg_indices_1d
    deallocate(dgg_indices_1d)

      call LDT_SMOS_DGG_lookup(n)
   endif

#endif

end subroutine read_SMOSNRTL2sm_data


subroutine find_SMOS_Dgg_id_number(n,c,r, lon_dgg, lat_dgg, & 
                                   dgg_length, i_dgg, DGG_id_number, &
                                   min_lon_dgg, max_lon_dgg, &
                                   min_lat_dgg, max_lat_dgg, &
                                   dx, dy, lon2d, lat2d)
! 
! !USES:
  use SMOSNRTNNL2sm_obsMod

  implicit none

! !INPUT PARAMETERS: 
! 
  integer              :: n,c,r
  integer              :: dgg_length
  real                 :: lon_dgg(dgg_length), lat_dgg(dgg_length)
  real                 :: max_lon_dgg, min_lon_dgg, max_lat_dgg, min_lat_dgg
  integer              :: DGG_id_number(dgg_length)
!
! !OUTPUT PARAMETERS:
!
!EOP
  integer              :: i_dgg(dgg_length)
  real                 :: lat2d, lon2d
  real                 :: dx, dy
  integer, allocatable :: indices(:)

  indices = pack([i_dgg], lon_dgg >= lon2d-dx/2.and.&
                          lon_dgg <= lon2d+dx/2.and.&
                          lat_dgg >= lat2d-dy/2.and.&
                          lat_dgg <= lat2d+dy/2)

  if (size(indices) > 0) then
     allocate(SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices(size(indices)))

     SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_indices = DGG_id_number(indices)
     SMOSNRTNNsmobs(n)%SMOS_lookup(c,r)%dgg_assign = .true.

  endif

  deallocate(indices) 

end subroutine find_SMOS_Dgg_id_number


