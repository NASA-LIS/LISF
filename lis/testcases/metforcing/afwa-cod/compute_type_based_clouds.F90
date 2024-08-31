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
! !ROUTINE: compute_type_based_clouds
! \label{compute_type_based_clouds}
!
! !REVISION HISTORY:
! 20 Apr 2016; James Geiger, Initial specification
!
! !INTERFACE:
subroutine compute_type_based_clouds(n, cldamt_nh, cldamt_sh, cldamt, &
                                     cldtyp_nh, cldtyp_sh, fog_nh, fog_sh)
!  !USES:
   use LIS_coreMod,       only : LIS_rc, LIS_domain, LIS_localPet
   use LIS_logMod,        only : LIS_logunit
   use LIS_timeMgrMod,    only : LIS_get_julhr,LIS_tick,LIS_time2date
   use AGRMET_forcingMod, only : agrmet_struc
!<debug -- jim testing>
use LIS_mpiMod
use LIS_historyMod
!</debug -- jim testing>

   implicit none
!  !ARGUMENTS:
   integer,intent(in)  :: n
   integer,intent(out) :: cldamt_nh(3,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
   integer,intent(out) :: cldamt_sh(3,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
   real,   intent(out) :: cldamt(3, LIS_rc%lnc(n), LIS_rc%lnr(n))
   integer,intent(out) :: cldtyp_nh(3,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
   integer,intent(out) :: cldtyp_sh(3,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
   logical,intent(out) :: fog_nh(agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
   logical,intent(out) :: fog_sh(agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
!EOP

   integer               :: hemi
   byte,allocatable      :: amounts   (:,:,:)
   byte,allocatable      :: types     (:,:,:)
   integer*2,allocatable :: tops      (:,:,:)
   integer*4,allocatable :: times     (:,:)
   integer               :: doy1,yr1,mo1,da1,hr1,mn1,ss1,try
   character*200         :: filename
   logical               :: file_exists
   integer               :: istat
   real*8                :: backtime1
   real                  :: gmt1,ts1
   integer               :: julhr
   integer               :: thres(5) !needs to be read in.
   integer               :: ip
   real                  :: udef

   real,allocatable      :: cldamt1(:,:,:)
   real,allocatable      :: cldamt2(:,:,:)
   real,allocatable      :: cldamt3(:,:,:)

   data THRES /12600, 12300, 12000, 11700, 11400/

   allocate(amounts(4,1024,1024))
   allocate(types(4,1024,1024))
   allocate(tops(4,1024,1024))
   allocate(times(1024,1024))

   allocate(cldamt1(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
   allocate(cldamt2(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))
   allocate(cldamt3(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax))

   do hemi = 1,2
      types   = 25
      amounts = 0
      tops    = 0
      times   = 0

      yr1 = LIS_rc%yr
      mo1 = LIS_rc%mo
      da1 = LIS_rc%da
      hr1 = LIS_rc%hr
      mn1 = LIS_rc%mn
      ss1 = LIS_rc%ss

      call agrmet_cdfs_type_filename(filename,                      &
                                     agrmet_struc(n)%agrmetdir,     &
                                     agrmet_struc(n)%clouddir,      &
                                     agrmet_struc(n)%use_timestamp, &
                                     hemi, yr1, mo1, da1, hr1)

      inquire (file=trim(filename), exist=file_exists)
      if(file_exists) then
         write(LIS_logunit,*)'[INFO] OPENING CDFS II DATA ', trim(filename)
         open(8, file=trim(filename), access='direct', recl=1024*1024*4, &
              status='old', iostat=istat)
         read(8, rec=1, iostat=istat) types
         close(8)
      else !rolling back 1 hour at a time
         try  = 1
         do while(try.lt.10)
            write(LIS_logunit,*) &
               '[WARN] File missing, shifting to previous hour'
            ts1 = -60*60
            call LIS_tick(backtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
            call LIS_time2date(backtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1)
            call agrmet_cdfs_type_filename(filename,                      &
                                           agrmet_struc(n)%agrmetdir,     &
                                           agrmet_struc(n)%clouddir,      &
                                           agrmet_struc(n)%use_timestamp, &
                                           hemi, yr1, mo1, da1, hr1)
            inquire (file=trim(filename), exist=file_exists)
            if(file_exists) then
               write(LIS_logunit,*)'[INFO] OPENING CDFS II DATA ',trim(filename)
               open(8, file=trim(filename), access='direct', recl=1024*1024*4, &
                    status='old', iostat=istat)
               read(8, rec=1, iostat=istat) types
               close(8)
               exit;
            else
               try = try+1
            endif
         enddo
         if(try.ge.10) then
            write(LIS_logunit,*) '[WARN] Missing file ',trim(filename)
            write(LIS_logunit,*) '[WARN] Leaving types set to zero...'
         endif
      endif

      yr1 = LIS_rc%yr
      mo1 = LIS_rc%mo
      da1 = LIS_rc%da
      hr1 = LIS_rc%hr
      mn1 = LIS_rc%mn
      ss1 = LIS_rc%ss

      call agrmet_cdfs_pcts_filename(filename,                      &
                                     agrmet_struc(n)%agrmetdir,     &
                                     agrmet_struc(n)%clouddir,      &
                                     agrmet_struc(n)%use_timestamp, &
                                     hemi, yr1, mo1, da1, hr1)

      inquire (file=trim(filename), exist=file_exists)
      if(file_exists) then
         write(LIS_logunit,*)'[INFO] OPENING CDFS II DATA ', trim(filename)
         open(8, file=trim(filename), access='direct', recl=1024*1024*4, &
              status='old', iostat=istat)
         read(8, rec=1, iostat=istat) amounts
         close(8)
      else !rolling back 1 hour at a time
         try  = 1
         do while(try.lt.10)
            write(LIS_logunit,*) &
               '[WARN] File missing, shifting to previous hour'
            ts1 = -60*60
            call LIS_tick(backtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
            call LIS_time2date(backtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1)
            call agrmet_cdfs_pcts_filename(filename,                      &
                                           agrmet_struc(n)%agrmetdir,     &
                                           agrmet_struc(n)%clouddir,      &
                                           agrmet_struc(n)%use_timestamp, &
                                           hemi, yr1, mo1, da1, hr1)
            inquire (file=trim(filename), exist=file_exists)
            if(file_exists) then
               write(LIS_logunit,*)'[INFO] OPENING CDFS II DATA ',trim(filename)
               open(8, file=trim(filename), access='direct', recl=1024*1024*4, &
                    status='old', iostat=istat)
               read(8, rec=1, iostat=istat) amounts
               close(8)
               exit;
            else
               try = try+1
            endif
         enddo
         if(try.ge.10) then
            write(LIS_logunit,*) '[WARN] Missing file ',trim(filename)
            write(LIS_logunit,*) '[WARN] Leaving amounts set to zero...'
         endif
      endif

      yr1 = LIS_rc%yr
      mo1 = LIS_rc%mo
      da1 = LIS_rc%da
      hr1 = LIS_rc%hr
      mn1 = LIS_rc%mn
      ss1 = LIS_rc%ss

      call agrmet_cdfs_hgts_filename(filename,                      &
                                     agrmet_struc(n)%agrmetdir,     &
                                     agrmet_struc(n)%clouddir,      &
                                     agrmet_struc(n)%use_timestamp, &
                                     hemi, yr1, mo1, da1, hr1)

      inquire (file=trim(filename), exist=file_exists)
      if(file_exists) then
         write(LIS_logunit,*)'[INFO] OPENING CDFS II DATA ', trim(filename)
         open(8, file=trim(filename), access='direct', recl=1024*1024*8, &
              status='old', iostat=istat)
         read(8, rec=1, iostat=istat) tops
         close(8)
      else !rolling back 1 hour at a time
         try  = 1
         do while(try.lt.10)
            write(LIS_logunit,*) &
               '[WARN] File missing, shifting to previous hour'
            ts1 = -60*60
            call LIS_tick(backtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
            call LIS_time2date(backtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1)
            call agrmet_cdfs_hgts_filename(filename,                      &
                                           agrmet_struc(n)%agrmetdir,     &
                                           agrmet_struc(n)%clouddir,      &
                                           agrmet_struc(n)%use_timestamp, &
                                           hemi, yr1, mo1, da1, hr1)
            inquire (file=trim(filename), exist=file_exists)
            if(file_exists) then
               write(LIS_logunit,*)'[INFO] OPENING CDFS II DATA ',trim(filename)
               open(8, file=trim(filename), access='direct', recl=1024*1024*8, &
                    status='old', iostat=istat)
               read(8, rec=1, iostat=istat) tops
               close(8)
               exit;
            else
               try = try+1
            endif
         enddo
         if(try.ge.10) then
            write(LIS_logunit,*) '[WARN] Missing file ',trim(filename)
            write(LIS_logunit,*) '[WARN] Leaving tops set to zero...'
         endif
      endif

      yr1 = LIS_rc%yr
      mo1 = LIS_rc%mo
      da1 = LIS_rc%da
      hr1 = LIS_rc%hr
      mn1 = LIS_rc%mn
      ss1 = LIS_rc%ss

      call agrmet_cdfs_pixltime_filename(filename,                      &
                                         agrmet_struc(n)%agrmetdir,     &
                                         agrmet_struc(n)%clouddir,      &
                                         agrmet_struc(n)%use_timestamp, &
                                         hemi, yr1, mo1, da1, hr1)
      inquire (file=trim(filename), exist=file_exists)
      if(file_exists) then
         write(LIS_logunit,*) '[INFO] READING CDFS II DATA ',trim(filename)
         open(8, file=trim(filename), access='direct', recl=1024*1024*4, &
              status='old', iostat=istat)
         read(8, rec=1, iostat=istat) times
         close(8)
      else !rolling back 1 hour at a time
         try  = 1
         do while(try.lt.10)
            write(LIS_logunit,*) &
               '[WARN] File missing, shifting to previous hour'
            ts1 = -60*60
            call LIS_tick(backtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
            call LIS_time2date(backtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1)
            call agrmet_cdfs_pixltime_filename(filename,                      &
                                               agrmet_struc(n)%agrmetdir,     &
                                               agrmet_struc(n)%clouddir,      &
                                               agrmet_struc(n)%use_timestamp, &
                                               hemi, yr1, mo1, da1, hr1)
            inquire (file=trim(filename), exist=file_exists)
            if(file_exists) then
               write(LIS_logunit,*)'[INFO] OPENING CDFS II DATA ',trim(filename)
               open(8, file=trim(filename), access='direct', recl=1024*1024*4, &
                    status='old', iostat=istat)
               read(8, rec=1, iostat=istat) times
               close(8)
               exit;
            else
               try = try+1
            endif
         enddo
         if(try.ge.10) then
            write(LIS_logunit,*) '[WARN] Missing file ',trim(filename)
            write(LIS_logunit,*) '[WARN] Leaving times set to zero...'
         endif
      endif

      yr1 = LIS_rc%yr
      mo1 = LIS_rc%mo
      da1 = LIS_rc%da
      hr1 = LIS_rc%hr
      mn1 = LIS_rc%mn
      ss1 = LIS_rc%ss

      times = times/60

      call LIS_get_julhr(yr1, mo1, da1, hr1, &
         0, 0,julhr)
!<debug -- jim testing>
write(LIS_logunit,*) 'GREP: julhr ', julhr
write(LIS_logunit,*) 'GREP: times ', minval(times), maxval(times)
call flush(LIS_logunit)
call MPI_Barrier(LIS_mpi_comm, istat)
!</debug -- jim testing>
      if(hemi.eq.1) then
         call AGRMET_loadcloud(hemi,agrmet_struc(n)%land(:,:,hemi),thres, &
                               times,amounts,tops,types,                  &
                               cldtyp_nh,cldamt_nh,fog_nh,                &
                               julhr, agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
      else
         call AGRMET_loadcloud(hemi,agrmet_struc(n)%land(:,:,hemi), &
                               thres,times,amounts,tops,types,      &
                               cldtyp_sh,cldamt_sh,fog_sh,          &
                               julhr, agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
      endif
   end do

   cldamt1(1,:,:) = real(cldamt_nh(1,:,:))
   cldamt1(2,:,:) = real(cldamt_sh(1,:,:))
   cldamt2(1,:,:) = real(cldamt_nh(2,:,:))
   cldamt2(2,:,:) = real(cldamt_sh(2,:,:))
   cldamt3(1,:,:) = real(cldamt_nh(3,:,:))
   cldamt3(2,:,:) = real(cldamt_sh(3,:,:))

   ip =1
   udef = -1.0

   call interp_agrmetvar(n,ip,cldamt1,udef,cldamt(1,:,:), &
                         agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
   call interp_agrmetvar(n,ip,cldamt2,udef,cldamt(2,:,:), &
                         agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
   call interp_agrmetvar(n,ip,cldamt3,udef,cldamt(3,:,:), &
                         agrmet_struc(n)%imax,agrmet_struc(n)%jmax)

!<debug -- jim testing>
if ( LIS_localPet == 0 ) then
open(unit=667,file='cldamt.bin',form='unformatted',access='direct',recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
endif
call LIS_writevar_bin(667, n, cldamt(1,:,:),1)
call LIS_writevar_bin(667, n, cldamt(2,:,:),2)
call LIS_writevar_bin(667, n, cldamt(3,:,:),3)
call MPI_Barrier(LIS_mpi_comm, istat)
if ( LIS_localPet == 0 ) then
flush(667)
close(667)
endif
!</debug -- jim testing>
   deallocate(amounts)
   deallocate(types)
   deallocate(tops)
   deallocate(times)

   deallocate(cldamt1)
   deallocate(cldamt2)
   deallocate(cldamt3)

end subroutine compute_type_based_clouds
