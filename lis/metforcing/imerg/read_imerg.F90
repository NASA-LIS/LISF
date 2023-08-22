!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
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
!  07 Aug 2023: Jessica Erlingas; Added support for IMERG V07
!  16 Aug 2023: Eric Kemp; Improved graceful HDF5 library cleanup if
!               error is detected.
!
! !INTERFACE:
subroutine read_imerg (n, kk, name_imerg, findex, order, ferror_imerg )
! !USES:
  use LIS_coreMod,only : LIS_rc, LIS_domain
  use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber, &
                         LIS_releaseUnitNumber
  use LIS_metforcingMod,only : LIS_forc
  use imerg_forcingMod, only : imerg_struc
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

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
  character(len=LIS_CONST_PATH_LEN) :: fname ! Filename variables
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
   write(LIS_logunit,*) &
        "[INFO] Reading HDF5 IMERG precipitation data from ", trim(fname)
   call read_imerghdf(n, fname, xd, yd, realprecip, ireaderr)
   if (ireaderr .ne. 0) then
     write(LIS_logunit,*) &
        "[WARN] Error reading IMERG file ",trim(fname)
     ferror_imerg = 0
   endif
 else
   write(LIS_logunit,*) &
      "[WARN] Missing IMERG precipitation data:: ",trim(fname)
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
    write(LIS_logunit,*) "Obtained IMERG precipitation data from ", trim(fname)
 endif

 deallocate (precip_regrid)
    
end subroutine read_imerg


! J.Case (3/9/2015) -- below will be the HDF5 reader subroutine.
subroutine read_imerghdf(n, filename, xsize, ysize, precipout, istatus)
! !USES:
#if (defined USE_HDF5)
  use hdf5
#endif
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber
  use imerg_forcingMod, only : imerg_struc

  implicit none

! ARGUMENTS
  character(len=*)    :: filename
  integer, intent(in)  :: xsize, ysize
  character(len=40) :: dsetname
  real :: precipin(ysize,xsize)
  real :: precipout(xsize,ysize)
  integer :: istatus,i
  character(len=4) :: vlname
  integer         :: vnum, n
#if (defined USE_HDF5)
  integer(HSIZE_T), dimension(2) :: dims
  integer(HID_T) :: fileid, dsetid
  integer :: ierr

  dims(1) = xsize
  dims(2) = ysize

  ! Variable names changed in IMERG V07
  vlname = trim(imerg_struc(n)%imergver)
  read( vlname(2:3), '(I2)') vnum
  if (vnum .ge. 7) then
     dsetname='/Grid/precipitation'
  else
     dsetname='/Grid/precipitationCal'
  endif

  istatus = 0 ! istatus is returned from this subroutine
  ierr = 0    ! ierr is strictly a local variable

  ! Open Fortran interface
  call h5open_f(istatus)
  if (istatus .ne. 0) then
     ! No need to proceed, something is wrong with HDF5.
     write (LIS_logunit,*) '[WARN] Error opening HDF5 Fortran interface'
     return
  end if

  ! Open HDF5 file.
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, fileid, istatus)
  if (istatus .ne. 0) then
     ! We couldn't open the file.  Close the Fortran interface and return.
     ! Don't overwrite the istatus value.
     write (LIS_logunit,*) '[WARN] Error opening IMERG file', &
          trim(filename)
     call h5close_f(ierr)
     return
  end if

  ! Open dataset
  call h5dopen_f(fileid, dsetname, dsetid, istatus)
  if (istatus .ne. 0) then
     ! We can't open the dataset.  Close the file and Fortran interface.
     ! Don't overwrite the istatus value.
     write (LIS_logunit,*) '[WARN] Error opening IMERG dataset', &
             trim(dsetname)
     call h5fclose_f(fileid, ierr)
     call h5close_f(ierr)
     return
  end if

  ! Read dataset
  call h5dread_f(dsetid, H5T_NATIVE_REAL, precipin, dims, istatus)
  if (istatus .ne. 0) then
     ! We can't read the dataset.  Close the dataset, the file, and
     ! the Fortran interface.  Don't overwrite the istatus value.
     write (LIS_logunit,*) '[WARN] Error reading IMERG dataset', &
             trim(dsetname)
     call h5dclose_f(dsetid, ierr)
     call h5fclose_f(fileid, ierr)
     call h5close_f(ierr)
     return
  end if

  ! Put the real(1:,1:) on the precipout(0:,0:)
  ! precipin is (ysize,xsize) starting at (lon=-179.9,lat=-89.9)
  precipout(1:xsize,1:ysize)=transpose(precipin)

  ! Close the dataset, file, and Fortran interface.  We already have
  ! the IMERG data in the array, so we don't need to overwrite the
  ! istatus value.  But we'll still log problems.
  call h5dclose_f(dsetid, ierr)
  if (ierr .ne. 0) then
     write (LIS_logunit,*) '[WARN] Error closing IMERG dataset', &
          trim(dsetname)
  end if
  call h5fclose_f(fileid, ierr)
  if (ierr .ne. 0) then
     write (LIS_logunit,*) '[WARN] Error closing IMERG file', &
          trim(filename)
  end if
  call h5close_f(ierr)
  if (ierr .ne. 0) then
     write (LIS_logunit,*) '[WARN] Error closing HDF5 Fortran interface'
  end if

#endif

end subroutine read_imerghdf
