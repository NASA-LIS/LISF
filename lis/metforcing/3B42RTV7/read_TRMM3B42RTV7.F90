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
! !ROUTINE: read_TRMM3B42RTV7
! \label{read_TRMM3B42RTV7}
!
! !REVISION HISTORY: 
!  17 Jul 2001: Jon Gottschalck; Initial code
!  04 Feb 2002: Jon Gottschalck; Added necessary code to use global precip
!               observations with LIS_domain 3 (2x2.5)
!  06 Jan 2015: KR Arsenault; Added support for latest V7 RT dataset
!
! !INTERFACE:
subroutine read_TRMM3B42RTV7 (n, kk, filename_TRMM3B42RT, findex, &
                              order, ferror_TRMM3B42RT)
! !USES:
 use LIS_coreMod, only       : LIS_rc, LIS_domain
 use LIS_logMod, only        : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber
 use LIS_metforcingMod, only : LIS_forc
 use TRMM3B42RTV7_forcingMod, only : TRMM3B42RTV7_struc
 use LIS_constantsMod,        only : LIS_CONST_PATH_LEN
 
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: kk     ! Forecast ensemble member
  character(len=*)  :: filename_TRMM3B42RT
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_TRMM3B42RT
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  TRMM 3B42RT V7 data and interpolates to the LIS domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[name\_TRMM3B42RT V7]
!    name of the 3 hour TRMM 3B42RT V7 forecast file
!  \item[ferror\_TRMM3B42RTV7]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_TRMM3B42RTV7](\ref{interp_TRMM3B42RTV7}) \newline
!    spatially interpolates the TRMM 3B42RTV7 data
!  \end{description}
!EOP

  integer :: index1

!==== Local Variables=======================
       
  integer :: ios
  integer :: i, j, xd, yd
  real    :: precip(TRMM3B42RTV7_struc(n)%nc, TRMM3B42RTV7_struc(n)%nr)
  real    :: precip2(TRMM3B42RTV7_struc(n)%nc, TRMM3B42RTV7_struc(n)%nr)
  real    :: tmp(TRMM3B42RTV7_struc(n)%nc, TRMM3B42RTV7_struc(n)%nr)   
  real, allocatable  :: precip_regrid(:,:)   ! Interpolated precipitation array
 
  character(len=LIS_CONST_PATH_LEN) :: filename             ! Filename variables
  character(len=LIS_CONST_PATH_LEN) :: dirfile 
  integer            :: ftn
  logical            :: file_exists
 
!=== End Variable Definition =======================

  filename = filename_TRMM3B42RT
  xd = TRMM3B42RTV7_struc(n)%nc
  yd = TRMM3B42RTV7_struc(n)%nr

!------------------------------------------------------------------------
! Fill necessary arrays to assure not using old TRMM 3B42RT V7 data
!------------------------------------------------------------------------
! J.Case (4/22/2013) -- Make consistent with Stg4/NMQ routines
  if(order.eq.1) then 
     TRMM3B42RTV7_struc(n)%metdata1 = LIS_rc%udef 
  elseif(order.eq.2) then 
     TRMM3B42RTV7_struc(n)%metdata2 = LIS_rc%udef
  endif
  allocate (precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  precip_regrid = -1.0     ! J.Case

!------------------------------------------------------------------------
! Find TRMM 3B42RT V7 precip data, read it in and assign to forcing precip array.
! Must reverse grid in latitude dimension to be consistent with LIS grid
!------------------------------------------------------------------------

! Check for and read in a gzipped file:
  dirfile = trim(filename)//".bin.gz"
  inquire(file=dirfile, EXIST=file_exists) 
  if (file_exists) then 
     write(LIS_logunit, *) trim(filename)//".bin.gz"
     call read_3B42RTV7_gzip(dirfile, precip, xd, yd) 
  else 
   ! Check for and read in native binary file:
     dirfile = trim(filename)//".bin"
     inquire(file=dirfile, EXIST=file_exists) 
     if (file_exists) then 
        write(LIS_logunit, *) trim(filename)//".bin"
        call read_3B42RTV7_bin(dirfile, precip, xd, yd) 
     else 
      ! Check for and read in processed 1gd4r binary file:
        dirfile = trim(filename)//".1gd4r" 
        inquire(file=dirfile, EXIST=file_exists) 
        if (file_exists) then 
          write(LIS_logunit, *) trim(filename)//".1gd4r"
          call read_3B42RTV7_1gd4r(dirfile, precip, xd, yd) 
        else 
          dirfile = trim(filename)//" [.bin.gz|.bin|.1gd4r] file missing"
          precip = -1.0
        end if
     end if
  end if

 ! Because raw RT data goes from 0.125 to 359.875, need to swap 
 ! Western/Eastern Hemisphere to make it go: -179.875 to 179.875

  precip2 = precip

  Do j=1, yd
    Do i=1, xd/2  ! 720
       tmp(i, j) = precip(i+(xd/2), j)
       precip(i+(xd/2), j) = precip(i, j)
       precip(i, j) = tmp(i, j)
    End Do 
  End Do

!  print *,"after hemi-switch: ",dirfile(41:50),", x=1133, y=84:",&
!          precip2(1133,397),precip(413,397)

!------------------------------------------------------------------------
! Interpolating to desired LIS_domain and resolution
! Global precip datasets not used currently to force NLDAS
!------------------------------------------------------------------------

   call interp_TRMM3B42RTV7(n, xd, yd, precip, LIS_rc%gridDesc(n,:), &
        LIS_rc%lnc(n),LIS_rc%lnr(n),precip_regrid, findex)
    
   do j = 1,LIS_rc%lnr(n)
      do i = 1,LIS_rc%lnc(n)
         if (precip_regrid(i,j) .ne. -1.0) then
            index1 = LIS_domain(n)%gindex(i,j)
            if(index1 .ne. -1) then
               if(order.eq.1) then 
                  TRMM3B42RTV7_struc(n)%metdata1(kk,1,index1) = precip_regrid(i,j)   !here is mm/h
               elseif(order.eq.2) then 
                  TRMM3B42RTV7_struc(n)%metdata2(kk,1,index1) = precip_regrid(i,j)   !here is mm/h
               endif
            endif
         endif
      enddo
   enddo

   ferror_TRMM3B42RT = 1

   deallocate (precip_regrid)

 end subroutine read_TRMM3B42RTV7


!========== Read .bin file ===========

 subroutine read_3B42RTV7_bin(dirfile, precip, xd, yd)

  use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber, &
                         LIS_releaseUnitNumber
  use LIS_coreMod, only : LIS_rc

   implicit none
   character(len=*), intent(in) :: dirfile
   integer,       intent(in) :: xd, yd 
   real,       intent(inout) :: precip(xd,yd)

   integer   :: i, j, ftn
   real      :: output(xd, yd)
   integer*2 :: rr(xd, yd)
! _____________________________________________

   ftn = LIS_getNextUnitNumber()
   open(unit=ftn,file=dirfile, status='old', &
        access='direct',recl=xd*2, form='unformatted') 
   Do j=1, yd 
      read (ftn,rec=(j+1)) rr(:, j)  !skip 2880-byte header
   End Do 

   close(ftn)
   call LIS_releaseUnitNumber(ftn)

!- Convert integer to real, flip N-S, and set undef values
   Do j=1, yd
      Do i = 1, xd 
         if( rr(i, j).GE.0 ) then
           output(i, yd-j+1) = rr(i, j)*0.01
         else
           output(i, yd-j+1) = LIS_rc%udef   ! -9999.0
         end if
      End Do
   End Do

#if 0
! Subset to 50 N/S (where data are only valid):
!   Do j=41, 440
   Do j=41, (xd/2)
      precip(:, j-40) = output(:, j)     
   End Do 
#endif

   precip = output

 end subroutine read_3B42RTV7_bin 

!============ Read .bin.gz file =================

 subroutine read_3B42RTV7_gzip(zipfile, output, xd, yd)

  use LIS_coreMod, only : LIS_rc

  character(len=*), intent(in) :: zipfile
  integer, intent(in) :: xd, yd
  real, intent(inout) :: output(xd, yd)

  integer, parameter :: nc=1440, nr=480

  integer*2   :: input(xd,yd), itmp(nc, nr+1), rtmp(xd)  ! tmp includes header
  character*1 :: array(nc*(nr+1)*2), ct    ! buffer space
  integer     :: readzipf, dlen, rdlen, i, j, l
  equivalence (itmp, array)

! _________________________________________________________

  dlen=xd*(yd+1)*2

  rdlen = readzipf(trim(zipfile)//char(0), array, dlen)
! swap the bytes for big-endian representation
  Do l=1, xd*(yd+1)*2-1, 2
     ct = array(l)
     array(l) = array(l+1)
     array(l+1) = ct
  End Do
  if(rdlen .ne. dlen) then
     write(*, *) "rdlen=", rdlen, " File reading error ..."
     write(*, *) "Fill array with undef"
     itmp = -999
  end if

  Do j=1, yd 
     input(:, j) = itmp(:, j+1)
  End Do

!- Convert integer to real, flip N-S, and set undef values
   Do j=1, yd
!   Do j=41, 440
      Do i = 1, xd
         if( input(i, j).GE.0 ) then
           output(i, yd-j+1) = input(i, j)*0.01
!           output(i, yd-j+1-40) = input(i, j)*0.01
         else
           output(i, yd-j+1) = LIS_rc%udef
!           output(i, yd-j+1-40) = -9999.0
         end if
      End Do
   End Do

end subroutine read_3B42RTV7_gzip

!========== Read .1gd4r file ===========

 subroutine read_3B42RTV7_1gd4r(dirfile, precip, xd, yd)

  use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber, &
                         LIS_releaseUnitNumber

  character(len=*), intent(in) :: dirfile
  integer, intent(in)       :: xd,yd
  real,    intent(inout)    :: precip(xd,yd)

  integer :: i, j
  integer :: ftn

  ftn = LIS_getNextUnitNumber()
  open(unit=ftn,file=dirfile, status='old', &
       access='direct',recl=xd*yd*4, form='unformatted')
  read (ftn,rec=1) precip

  close(ftn)
  call LIS_releaseUnitNumber(ftn)

 end subroutine read_3B42RTV7_1gd4r

! ========================================
