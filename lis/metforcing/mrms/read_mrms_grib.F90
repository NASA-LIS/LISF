!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_mrms
! \label{read_mrms}
!
! !REVISION HISTORY:
!  25 May 2006: Kristi Arsenault;  Data and code implementation
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!  12 Feb 2015: Jonathan Case; revised routine for MRMS
!  05 Sep 2017; Jessica Erlingis; modified routine for operational MRMS
!  29 Mar 2018: Jessica Erlingis; add options for monthly-variying static mask
!  22 Feb 2019; Jessica Erlingis; added options to config file for mask
!
! !INTERFACE:
subroutine read_mrms_grib( n, fname, findex, order, yr, mo, da, ferror_mrms_grib )

! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod,         only : LIS_logunit, LIS_verify
  use LIS_metforcingMod,  only : LIS_forc
  use mrms_grib_forcingMod,    only : mrms_grib_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=120)   :: fname
  character(len=200)  :: maskname          
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_mrms_grib
  integer             :: yr, mo, da
  real                :: maxprec, maxinterp
  character(2)        :: cmon
! !DESCRIPTION:
!  For the given time, reads parameters from MRMS operational datasets
!  and interpolates to a designated user-domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the hourly MRMS file
!  \item[ferror\_mrms]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_mrms_grib](\ref{interp_mrms_grib})\\
!    spatially interpolates the MRMS data
!  \end{description}
!EOP

  integer            :: i, j, t, iret  ! Loop indicies and error flags
  integer            :: ftn, ndata, index1, ftn2
  real, allocatable  :: precip_regrid(:,:)
!  real               :: precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n))     
!  real               :: mrms_grib_in(mrms_grib_struc(n)%ncol*mrms_grib_struc(n)%nrow)
!  real               :: mrms_mask(mrms_grib_struc(n)%ncol*mrms_grib_struc(n)%nrow)
!  logical*1          :: lb(mrms_grib_struc(n)%ncol*mrms_grib_struc(n)%nrow)
  real, allocatable  :: mrms_grib_in(:)
  real, allocatable  :: mrms_mask(:)
  logical*1, allocatable :: lb(:)

  logical            :: file_exists            ! Check MRMS file status
  logical            :: mask_file_exists 
  integer            :: igrib,igrib2
  real               :: missingValue =-1 !MRMS missing value for precip fields
  real               :: noCovValue = -3  !MRMS value for no coverage field
  real               :: upperbound = 350 !upper bound for precip [mm]
  real               :: maskthresh = 50. !default upper  bound for static RQI mask
! See http://www.nssl.noaa.gov/projects/mrms/operational/tables.php for missing and 
! no coverage values for various MRMS fields 

!=== End Variable Definition =======================

  if(order.eq.1) then 
     mrms_grib_struc(n)%metdata1 = LIS_rc%udef
  elseif(order.eq.2) then 
     mrms_grib_struc(n)%metdata2 = LIS_rc%udef
  endif

!-- Set necessary parameters for call to interp_mrms_grib   
  precip_regrid = LIS_rc%udef

!-- Check initially if file exists:
  inquire (file=fname, exist=file_exists ) ! Check if file exists
  if (.not. file_exists)  then 
     if (LIS_masterproc) write(LIS_logunit,*) "** Missing MRMS precipitation file: ", fname
     ferror_mrms_grib = 1
     return
  endif

#if(defined USE_GRIBAPI)

  allocate(mrms_grib_in(mrms_grib_struc(n)%ncol*mrms_grib_struc(n)%nrow))
  
  ndata = mrms_grib_struc(n)%ncol * mrms_grib_struc(n)%nrow
  mrms_grib_in = missingValue !LIS_rc%udef

  call grib_open_file(ftn,trim(fname),'r',iret)
  call LIS_verify(iret,'error grib_open_file in read_mrms_grib')

  call grib_new_from_file(ftn,igrib,iret)
  call LIS_verify(iret,'error in grib_new_from_file in read_mrms_grib')

  call grib_get(igrib,'values',mrms_grib_in,iret)
  call LIS_verify(iret, 'error in grib_get: values in read_mrms_grib')

  !write(LIS_logunit,*) 'Max MRMS value: ', maxval(mrms_grib_in)
  !write(LIS_logunit,*) 'Min MRMS value: ', minval(mrms_grib_in)

  ! Use this block to apply a monthly mask of average RQI

  if (mrms_grib_struc(n)%mrms_mask_opt.eq.1) then

    maskthresh = mrms_grib_struc(n)%mrms_mask_thresh

    allocate(mrms_mask(mrms_grib_struc(n)%ncol*mrms_grib_struc(n)%nrow))

    mrms_mask = 0.0

    if ( mo < 10 ) write ( cmon, '(a1,i1)' ) "0", mo
    if ( mo >= 10) write ( cmon, '(i2)' ) mo

    maskname = trim(mrms_grib_struc(n)%mrms_mask_dir)//'AvgRQI_'//cmon//'_conus_sm.grib2'
    write(LIS_logunit,*) 'Using mask ',maskname
    !write(LIS_logunit,*) 'With cutoff theshold ',maskthresh

    call grib_open_file(ftn2,trim(maskname),'r',iret)
    call LIS_verify(iret,'error grib_open_file for RQI mask')
  
    call grib_new_from_file(ftn2,igrib2,iret)
    call LIS_verify(iret,'error in grib_new_from_file for RQI mask')
 
    call grib_get(igrib2,'values',mrms_mask,iret)
    call LIS_verify(iret,'error in grib_get: values for RQI mask')

    !write(LIS_logunit,*) 'Maximum mask value: ', maxval(mrms_mask)
    !write(LIS_logunit,*) 'Minimum mask value: ', minval(mrms_mask)

  endif

  allocate(lb(mrms_grib_struc(n)%ncol*mrms_grib_struc(n)%nrow))

  lb = .false.
  do t = 1, ndata
    lb(t)=.false.
    if (mrms_grib_struc(n)%mrms_mask_opt.eq.1) then
      if (mrms_mask(t).lt.maskthresh) then
         mrms_grib_in(t)=LIS_rc%udef
         lb(t) = .false.
      endif
    endif
    !write(LIS_logunit,*) 'Current Mask Value: ',mrms_mask(t)
    if ((mrms_grib_in(t) .gt. upperbound)) then
      mrms_grib_in(t) =  LIS_rc%udef !missingValue !0.0
      lb(t) = .false.
    !elseif ((mrms_grib_struc(n)%mrms_mask_opt.eq.1).and.(mrms_mask(t).lt.maskthresh)) then
    !    !write(LIS_logunit,*) 'Does not meet mask value '
    !  mrms_grib_in(t) = 0.0 !LIS_rc%udef !missingValue !0.0
    !  lb(t) = .false.
    elseif (( mrms_grib_in(t) .ne. missingValue ) .and. ( mrms_grib_in(t) .ne. noCovValue) &
          .and. (mrms_grib_in(t) .ge. 0) .and. (mrms_grib_in(t).ne.LIS_rc%udef) .and. &
          (abs(mrms_grib_in(t)) .le. upperbound)) then
        !write(LIS_logunit,*) 'Good value: ',mrms_grib_in(t)
        lb(t) = .true. 
    endif
  enddo

  if (mrms_grib_struc(n)%mrms_mask_opt.eq.1) then    
    deallocate(mrms_mask)
    call grib_release(igrib2,iret)
    call LIS_verify(iret,'error in grib_release in mask')
    call grib_close_file(ftn2)
  endif

  allocate(precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  precip_regrid(:,:) = LIS_rc%udef
     
  write(LIS_logunit,*) "interpolating MRMS data -- "
  call interp_mrms_grib( n, findex, ndata, mrms_grib_in, lb, LIS_rc%gridDesc(n,:), &
       LIS_rc%lnc(n), LIS_rc%lnr(n), precip_regrid )
        
  do j = 1, LIS_rc%lnr(n)
    do i = 1, LIS_rc%lnc(n)
           
       if (( precip_regrid(i,j) .ne. LIS_rc%udef ) .and. &
         (precip_regrid(i,j) .le. upperbound)) then       !Sanity check to prevent bad values

         index1 = LIS_domain(n)%gindex(i,j)
         if(index1 .ne. -1) then
           if(order .eq. 1) then 
             mrms_grib_struc(n)%metdata1(1,index1) = precip_regrid(i,j) 
           elseif(order .eq. 2) then 
             mrms_grib_struc(n)%metdata2(1,index1) = precip_regrid(i,j) 
           endif
         endif
       endif
    enddo
  enddo

! Deallocate local variables

  deallocate(mrms_grib_in)
  deallocate(precip_regrid)
  deallocate(lb)

  call grib_release(igrib,iret)
  call LIS_verify(iret,'error in grib_release in read_mrms_grib')

  ferror_mrms_grib = 0 
  write(LIS_logunit,*) "- Obtained MRMS precipitation data from: ", trim(fname)
  write(LIS_logunit,*) " ------------------------------------- "

  call grib_close_file(ftn)

#endif

end subroutine read_mrms_grib
