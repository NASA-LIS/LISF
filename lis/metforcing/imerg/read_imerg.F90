!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_imerg
! \label{read_imerg}
!
! !REVISION HISTORY:
!  17 Jul 2001: Jon Gottschalck; Initial code
!  09 Mar 2015: Jon Case;  Added IMERG precipitation reader
!
! !INTERFACE:
subroutine read_imerg (n, kk, name_imerg, findex, order, ferror_imerg )
! !USES:
  use LIS_coreMod,only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber, &
                         LIS_releaseUnitNumber
  use LIS_metforcingMod,only : LIS_forc
  use imerg_forcingMod, only : imerg_struc

  implicit none
! !ARGUMENTS:   
  integer, intent(in) :: n 
  integer, intent(in) :: kk
  character(len=*)    :: name_imerg
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_imerg
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  IMERG data and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[name\_imerg]
!    name of the IMERG file
!  \item[ferror\_imerg]
!    flag to indicate success of the call (=0 indicates success)
!  \item[iflg]
!    flag indicating which 1/2 hour to read
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_imerg](\ref{interp_imerg}) \newline
!    spatially interpolates the IMERG data
!  \end{description}
!EOP

  integer :: index1

!==== Local Variables=======================

  integer :: r,ios
  integer :: i,j,xd,yd
  parameter(xd=3600,yd=1800)    ! Dimension of original IMERG 0.1-deg data
  character*1  precip(xd,yd),timestamp(xd,yd)
  character*1  staid(xd,yd)
  character*1  testout1(xd,yd)                             ! Reconfigured original precip array
  real :: realprecip(xd,yd)
  real :: testout(xd,yd)
  real, allocatable :: precip_regrid(:,:)                      ! Interpolated precip array
  character(len=99) :: fname, zname                  ! Filename variables
  logical           :: file_exists
  integer           :: ftn
  integer           :: ireaderr
!=== End Variable Definition =======================

  allocate (precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  precip_regrid = -1.0
  realprecip = -1.0
  fname = name_imerg
  if(order.eq.1) then 
     imerg_struc(n)%metdata1 = -1.0
  elseif(order.eq.2) then 
     imerg_struc(n)%metdata2 = -1.0
  endif

 ferror_imerg = 1
 inquire(file=fname, EXIST=file_exists)
 if (file_exists) then
   if(LIS_masterproc) write(LIS_logunit,*) &
        "[INFO] Reading HDF5 IMERG precipitation data from ", fname
   call read_imerghdf(fname, xd, yd, realprecip, ireaderr)
   if (ireaderr .ne. 0) then
     if(LIS_masterproc) write(LIS_logunit,*) &
        "[WARN] Error reading IMERG file ",fname
     ferror_imerg = 0
   endif
 else
   if(LIS_masterproc) write(LIS_logunit,*) &
      "[WARN] Missing IMERG precipitation data:: ",fname
   ferror_imerg = 0
 endif

 if(ferror_imerg.eq.1) then 
    call interp_imerg(n, xd, yd, realprecip, LIS_rc%gridDesc(n,:), &
         LIS_rc%lnc(n),LIS_rc%lnr(n),precip_regrid)
    do j = 1,LIS_rc%lnr(n)
       do i = 1,LIS_rc%lnc(n)
          if (precip_regrid(i,j) .ge. 0.0) then
             index1 = LIS_domain(n)%gindex(i,j)
             if(index1 .ne. -1) then
                if(order.eq.1) then 
                   imerg_struc(n)%metdata1(kk,1,index1) = &
                        precip_regrid(i,j)   !here is mm/h
                elseif(order.eq.2) then 
                   imerg_struc(n)%metdata2(kk,1,index1) = &
                        precip_regrid(i,j)   !here is mm/h
                endif
             endif
          endif
       enddo
    enddo
    if (LIS_masterproc) &
      write(LIS_logunit,*) "Obtained IMERG precipitation data from ", trim(fname)
 endif

 deallocate (precip_regrid)
    
end subroutine read_imerg


! J.Case (3/9/2015) -- below will be the HDF5 reader subroutine.
subroutine read_imerghdf(filename, xsize, ysize, precipout, istatus)
! !USES:
#if (defined USE_HDF5)
  use hdf5
#endif
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber

  implicit none

! ARGUMENTS
  character(len=99)    :: filename
  integer, intent(in)  :: xsize, ysize

  character(len=40) :: dsetname='/Grid/precipitationCal'
  real :: precipin(ysize,xsize)
  real :: precipout(xsize,ysize)
  logical :: bIsError
  integer :: istatus,i
#if (defined USE_HDF5)
  integer(HSIZE_T), dimension(2) :: dims
  integer(HID_T) :: fileid,dsetid

  dims(1) = xsize
  dims(2) = ysize

  bIsError=.false.
!open fortran interface
  call h5open_f(istatus)
  if(istatus.ne.0) then
    bIsError=.true.
    if (LIS_masterproc) write (LIS_logunit,*) 'Error opening HDF5 fortran interface'
    return
  endif
!open hdf5 file
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,fileid,istatus)
  if(istatus.ne.0) then
    bIsError=.true.
    if (LIS_masterproc) write (LIS_logunit,*) 'Error opening IMERG file',trim(filename)
    return
  endif
!open dataset
  call h5dopen_f(fileid,dsetname,dsetid,istatus)
  if(istatus.ne.0) then
    bIsError=.true.
    if (LIS_masterproc) write (LIS_logunit,*) 'Error opening IMERG dataset',trim(dsetname)
    return
  endif
!read dataset 
  call h5dread_f(dsetid,H5T_NATIVE_REAL,precipin,dims,istatus)
  if(istatus.ne.0) then
    bIsError=.true.
    if (LIS_masterproc) write (LIS_logunit,*) 'Error reading IMERG dataset',trim(dsetname)
    return
  endif
!Put the real(1:,1:) on the precipout(0:,0:)
!precipin is (ysize,xsize) starting at (lon=-179.9,lat=-89.9)
  precipout(1:xsize,1:ysize)=transpose(precipin)

!close dataset
  call h5dclose_f(dsetid,istatus)
  if(istatus.ne.0) then
    bIsError=.true.
    if (LIS_masterproc) write (LIS_logunit,*) 'Error closing IMERG dataset',trim(dsetname)
    return
  endif
!close file
  call h5fclose_f(fileid,istatus)
  if(istatus.ne.0) then
    bIsError=.true.
    if (LIS_masterproc) write (LIS_logunit,*) 'Error closing IMERG file',trim(filename)
    return
  endif
!close fortran interface
  call h5close_f(istatus)
  if(istatus.ne.0) then
    bIsError=.true.
    if (LIS_masterproc) write (LIS_logunit,*) 'Error closing HDF5 fortran interface'
    return
  endif

#endif

end subroutine read_imerghdf
