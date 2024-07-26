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
! !ROUTINE: compute_cod_based_clouds
! \label{compute_cod_based_clouds}
!
! !REVISION HISTORY:
! 28 Apr 2016: James Geiger; Initial specification
!
! !INTERFACE:
subroutine compute_cod_based_clouds(n, cod, cloud_pcts)
!  !USES:
   use LIS_coreMod!,       only : LIS_rc
   use LIS_logMod,        only : LIS_logunit
   use LIS_timeMgrMod,    only : LIS_get_julhr,LIS_tick,LIS_time2date
   use AGRMET_forcingMod, only : agrmet_struc
!<debug -- jim testing>
use LIS_historyMod
use LIS_mpiMod
!</debug -- jim testing>

   implicit none
!  !ARGUMENTS:
   integer, intent(in) :: n
   real, intent(out)   :: cod(3, LIS_rc%lnc(n), LIS_rc%lnr(n))
   real, intent(out)   :: cloud_pcts(3, LIS_rc%lnc(n), LIS_rc%lnr(n))

   integer             :: hemi
   integer             :: yr1,mo1,da1,hr1,mn1,ss1,julhr,doy1
   character*200       :: filename
   logical             :: file_exists
   integer             :: ip
   real                :: udef
   integer             :: try
   real*8              :: backtime1
   real                :: ts1,gmt1
!<debug -- jim testing>
integer :: istat
character*200       :: cfilename
real, allocatable, dimension(:,:) :: temp_cod
!</debug -- jim testing>

! !DESCRIPTION:
! This routine processes the four-layer CDFS II cloud optical depth into
! High, Middle, and Low cloud amounts for computing longwave radiation.
!
!  The arguments and variables are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[cod]
!    processed high, middle, low cloud optical depth
!  \item[cloud\_pcts]
!    processed high, middle, low cloud amounts
!  \item[hemi]
!    index of hemisphere loops
!  \item[yr1,mo1,da1,hr1,mn1,ss1,ts1,julhr,backtime1,gmt1,doy1]
!    time/date specific variables
!  \item[try]
!    counter for rolling back
!  \item[filename]
!    generated filename
!  \item[file\_exists]
!    check for archived file
!  \item[ip]
!    choice of interpoation algorithm
!  \item[udef]
!    undefined variable used for spatial interpolation
!  \item[cod\_4\_nh, cod\_4\_sh]
!    hemispheric cloud optical depth (each contains four layers)
!  \item[cod\_lmh\_nh, cod\_lmh\_sh]
!    hemispheric cloud optical depth for low, middle, and high layers
!  \item[cod\_1, cod\_2, cod\_3]
!    cloud optical depth for low, middle, and high layers
!    (each contains both hemispheres)
!  \item[cod]
!    spatially interpolated cloud optical depth for low, middle, and high layers
!  \item[cloud\_base\_4\_nh, cloud\_base\_4\_sh]
!    hemispheric cloud base (each contains four layers)
!  \item[cloud\_top\_4\_nh, cloud\_top\_4\_sh]
!    hemispheric cloud top (each contains four layers)
!  \item[cloud\_pcts\_4\_nh, cloud\_pcts\_4\_sh]
!    hemispheric cloud amounts (each contains four layers)
!  \item[cloud\_pcts\_lmh\_nh, cloud\_pcts\_lmh\_sh]
!    hemispheric cloud amounts for low, middle, and high layers
!  \item[cloud\_pcts\_1, cloud\_pcts\_2, cloud\_pcts\_3]
!    cloud amounts for low, middle, and high layers
!    (each contains both hemispheres)
!  \item[cloud\_pcts]
!    spatially interpolated cloud amounts for low, middle, and high layers
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[agrmet\_cdfs\_cod\_filename](\ref{agrmet_cdfs_cod_filename}) \newline
!    generates the filename to read CDFS II cloud optical depth
!  \item[read\_cod](\ref{read_cod}) \newline
!    reads the CDFS II cloud optical depth
!  \item[interp\_agrmetvar](\ref{interp_agrmetvar}) \newline
!    spatial interpolation of an AGRMET variable to LIS grid
!  \end{description}
!EOP

   real,allocatable      :: cod_4_nh(:,:,:)
   real,allocatable      :: cod_4_sh(:,:,:)
   real,allocatable      :: cod_lmh_nh(:,:,:)
   real,allocatable      :: cod_lmh_sh(:,:,:)
   real,allocatable      :: cod_1(:,:,:)
   real,allocatable      :: cod_2(:,:,:)
   real,allocatable      :: cod_3(:,:,:)

   real,allocatable      :: cloud_base_4_nh(:,:,:)
   real,allocatable      :: cloud_base_4_sh(:,:,:)

   real,allocatable      :: cloud_top_4_nh(:,:,:)
   real,allocatable      :: cloud_top_4_sh(:,:,:)

   real,allocatable      :: cloud_pcts_4_nh(:,:,:)
   real,allocatable      :: cloud_pcts_4_sh(:,:,:)
   real,allocatable      :: cloud_pcts_lmh_nh(:,:,:)
   real,allocatable      :: cloud_pcts_lmh_sh(:,:,:)
   real,allocatable      :: cloud_pcts_1(:,:,:)
   real,allocatable      :: cloud_pcts_2(:,:,:)
   real,allocatable      :: cloud_pcts_3(:,:,:)

   real,allocatable      :: cloud_tot_pcts_nh(:,:)
   real,allocatable      :: cloud_tot_pcts_sh(:,:)

   real,allocatable      :: cloud_times(:,:)

   allocate(cod_4_nh(4,1024,1024))
   allocate(cod_4_sh(4,1024,1024))
   allocate(cod_lmh_nh(3,1024,1024))
   allocate(cod_lmh_sh(3,1024,1024))
   allocate(cod_1(2,1024,1024))
   allocate(cod_2(2,1024,1024))
   allocate(cod_3(2,1024,1024))

   allocate(cloud_base_4_nh(4,1024,1024))
   allocate(cloud_base_4_sh(4,1024,1024))

   allocate(cloud_top_4_nh(4,1024,1024))
   allocate(cloud_top_4_sh(4,1024,1024))

   allocate(cloud_pcts_4_nh(4,1024,1024))
   allocate(cloud_pcts_4_sh(4,1024,1024))
   allocate(cloud_pcts_lmh_nh(3,1024,1024))
   allocate(cloud_pcts_lmh_sh(3,1024,1024))
   allocate(cloud_pcts_1(2,1024,1024))
   allocate(cloud_pcts_2(2,1024,1024))
   allocate(cloud_pcts_3(2,1024,1024))

   allocate(cloud_tot_pcts_nh(1024,1024))
   allocate(cloud_tot_pcts_sh(1024,1024))

   allocate(cloud_times(1024,1024))

!<debug -- jim testing>
   allocate(temp_cod(LIS_rc%lnc(n), LIS_rc%lnr(n)))
!</debug -- jim testing>

   ! Low    = 1
   ! Middle = 2
   ! High   = 3

   cod_lmh_nh = 0.0
   cod_lmh_sh = 0.0
   cloud_pcts_lmh_nh = 0.0
   cloud_pcts_lmh_sh = 0.0

   do hemi = 1,2

      yr1 = LIS_rc%yr
      mo1 = LIS_rc%mo
      da1 = LIS_rc%da
      hr1 = LIS_rc%hr
      mn1 = LIS_rc%mn
      ss1 = LIS_rc%ss

      try = 0

      do ! read-or-rollback loop
         call LIS_get_julhr(yr1, mo1, da1, hr1, 0, 0, julhr)

         call agrmet_cdfs_cod_filename(filename,                      &
                                       agrmet_struc(n)%agrmetdir,     &
                                       agrmet_struc(n)%clouddir,      &
                                       agrmet_struc(n)%use_timestamp, &
                                       hemi, yr1, mo1, da1, hr1)

         inquire (file=trim(filename), exist=file_exists)
         if ( file_exists ) then
            write(LIS_logunit,*)'[INFO] OPENING CDFS II DATA ', trim(filename)
            if ( hemi == 1 ) then
               call read_cod(n, filename, hr1, cod_4_nh, cloud_base_4_nh, &
                  cloud_top_4_nh, cloud_pcts_4_nh,                        &
                  cloud_tot_pcts_nh, cloud_times)
               call process_cloud_layers(cod_4_nh, cloud_base_4_nh,       &
                  cloud_top_4_nh, cloud_pcts_4_nh, &
                  cloud_tot_pcts_nh,               &
                  julhr, cloud_times,              &
                  cod_lmh_nh, cloud_pcts_lmh_nh)
            else
               call read_cod(n, filename, hr1, cod_4_sh, cloud_base_4_sh, &
                  cloud_top_4_sh, cloud_pcts_4_sh,                        &
                  cloud_tot_pcts_sh, cloud_times)
               call process_cloud_layers(cod_4_sh, cloud_base_4_sh,       &
                  cloud_top_4_sh, cloud_pcts_4_sh, &
                  cloud_tot_pcts_sh,               &
                  julhr, cloud_times,              &
                  cod_lmh_sh, cloud_pcts_lmh_sh)
            endif
            exit ! read file; exit read-or-rollback loop
         else !rolling back 1 hour at a time
            write(LIS_logunit,*) &
               '[WARN] File missing, shifting to previous hour ',trim(filename)
            ts1 = -60*60
            call LIS_tick(backtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
            call LIS_time2date(backtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1)
            try = try + 1
            if ( try == 10 ) then
               write(LIS_logunit,*) '[WARN] Missing file ',trim(filename)
               write(LIS_logunit,*) '[WARN] Leaving cod set to zero...'
               exit ! cannot file file; give up; exit read-or-rollback loop
            endif
         endif
      enddo
   enddo

   ip = 3
   udef = -1.0

   cod_1(1,:,:) = cod_lmh_nh(1,:,:)
   cod_1(2,:,:) = cod_lmh_sh(1,:,:)
   cod_2(1,:,:) = cod_lmh_nh(2,:,:)
   cod_2(2,:,:) = cod_lmh_sh(2,:,:)
   cod_3(1,:,:) = cod_lmh_nh(3,:,:)
   cod_3(2,:,:) = cod_lmh_sh(3,:,:)

   call interp_agrmetvar(n,ip,cod_1,udef,cod(1,:,:),1024,1024)
   call interp_agrmetvar(n,ip,cod_2,udef,cod(2,:,:),1024,1024)
   call interp_agrmetvar(n,ip,cod_3,udef,cod(3,:,:),1024,1024)

   cloud_pcts_1(1,:,:) = cloud_pcts_lmh_nh(1,:,:)
   cloud_pcts_1(2,:,:) = cloud_pcts_lmh_sh(1,:,:)
   cloud_pcts_2(1,:,:) = cloud_pcts_lmh_nh(2,:,:)
   cloud_pcts_2(2,:,:) = cloud_pcts_lmh_sh(2,:,:)
   cloud_pcts_3(1,:,:) = cloud_pcts_lmh_nh(3,:,:)
   cloud_pcts_3(2,:,:) = cloud_pcts_lmh_sh(3,:,:)

   call interp_agrmetvar(n,ip,cloud_pcts_1,udef,cloud_pcts(1,:,:),1024,1024)
   call interp_agrmetvar(n,ip,cloud_pcts_2,udef,cloud_pcts(2,:,:),1024,1024)
   call interp_agrmetvar(n,ip,cloud_pcts_3,udef,cloud_pcts(3,:,:),1024,1024)
!<debug -- jim testing>
if ( LIS_localPet == 0 ) then
write(UNIT=cfilename,FMT='(a,i4,i2.2,i2.2,i2.2,a)') 'cod_levels.',yr1,mo1,da1,hr1,'.bin'
open(unit=666,file=trim(cfilename),access='direct',recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
endif
call LIS_writevar_bin(666, n, cod(1,:,:), 1)
call LIS_writevar_bin(666, n, cod(2,:,:), 2)
call LIS_writevar_bin(666, n, cod(3,:,:), 3)

call MPI_Barrier(LIS_mpi_comm, istat)
if ( LIS_localPet == 0 ) then
close(666)
endif


if ( LIS_localPet == 0 ) then
write(UNIT=cfilename,FMT='(a,i4,i2.2,i2.2,i2.2,a)') 'cod_layers.',yr1,mo1,da1,hr1,'.bin'
open(unit=666,file=trim(cfilename),access='direct',recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
endif

cod_1(1,:,:) = cod_4_nh(1,:,:)
cod_1(2,:,:) = cod_4_sh(1,:,:)
call interp_agrmetvar(n,ip,cod_1,udef,temp_cod(:,:),1024,1024)
call LIS_writevar_bin(666, n, temp_cod, 1)

cod_1(1,:,:) = cod_4_nh(2,:,:)
cod_1(2,:,:) = cod_4_sh(2,:,:)
call interp_agrmetvar(n,ip,cod_1,udef,temp_cod(:,:),1024,1024)
call LIS_writevar_bin(666, n, temp_cod, 2)

cod_1(1,:,:) = cod_4_nh(3,:,:)
cod_1(2,:,:) = cod_4_sh(3,:,:)
call interp_agrmetvar(n,ip,cod_1,udef,temp_cod(:,:),1024,1024)
call LIS_writevar_bin(666, n, temp_cod, 3)

cod_1(1,:,:) = cod_4_nh(4,:,:)
cod_1(2,:,:) = cod_4_sh(4,:,:)
call interp_agrmetvar(n,ip,cod_1,udef,temp_cod(:,:),1024,1024)
call LIS_writevar_bin(666, n, temp_cod, 4)

call MPI_Barrier(LIS_mpi_comm, istat)
if ( LIS_localPet == 0 ) then
close(666)
!stop 666
endif
!</debug -- jim testing>

   deallocate(cod_4_nh)
   deallocate(cod_4_sh)
   deallocate(cod_lmh_nh)
   deallocate(cod_lmh_sh)
   deallocate(cod_1)
   deallocate(cod_2)
   deallocate(cod_3)

   deallocate(cloud_base_4_nh)
   deallocate(cloud_base_4_sh)

   deallocate(cloud_top_4_nh)
   deallocate(cloud_top_4_sh)

   deallocate(cloud_pcts_4_nh)
   deallocate(cloud_pcts_4_sh)
   deallocate(cloud_pcts_lmh_nh)
   deallocate(cloud_pcts_lmh_sh)
   deallocate(cloud_pcts_1)
   deallocate(cloud_pcts_2)
   deallocate(cloud_pcts_3)

   deallocate(cloud_tot_pcts_nh)
   deallocate(cloud_tot_pcts_sh)

   deallocate(cloud_times)
!<debug -- jim testing>
   deallocate(temp_cod)
!</debug -- jim testing>

end subroutine compute_cod_based_clouds


!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: compute_tr_from_cod
! \label{compute_tr_from_cod}
!
! !REVISION HISTORY:
! 12 Jul 2016: James Geiger; Initial specification
!
! !INTERFACE:
subroutine compute_tr_from_cod(n, t, r, cod)
!  !USES:
   use LIS_coreMod, only : LIS_rc

   implicit none
!  !ARGUMENTS:
   integer, intent(in) :: n
   real, intent(out)   :: t(LIS_rc%lnc(n), LIS_rc%lnr(n))
   real, intent(out)   :: r(LIS_rc%lnc(n), LIS_rc%lnr(n))
   real, intent(in)    :: cod(LIS_rc%lnc(n), LIS_rc%lnr(n))

! !DESCRIPTION:
! This routine computes transmittance and reflectance from cloud optical depth.
!
!  The arguments and variables are:
!  \begin{description}
!    \item[n]
!       nest index
!    \item[t]
!       transmittance
!    \item[r]
!       reflectance
!    \item[cod]
!       cloud optical depth
!  \end{description}
!EOP

   integer :: i,j

   t = 0.0
   r = 0.0

   do j = 1, LIS_rc%lnr(n)
      do i = 1, LIS_rc%lnc(n)
         if ( cod(i,j) /= LIS_rc%udef ) then
            t(i,j) = exp(-1.0*cod(i,j))

            ! quality check
            if ( t(i,j) < 0.01 ) then
               t(i,j) = 0.0
            endif
         endif
      enddo
   enddo

   r = 1.0 - t

end subroutine compute_tr_from_cod
