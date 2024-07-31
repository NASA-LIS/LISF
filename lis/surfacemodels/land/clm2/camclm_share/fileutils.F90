!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <LIS_misc.h>
module fileutils

!public methods

    public :: get_filename  
    public :: set_filename
    public :: opnfil
    public :: getfil
    public :: putfil
    public :: getavu
    public :: relavu
    private:: shell_cmd

    logical, public :: lsmiou(99)  !I/O file unit numbers (1 to 99)

!=======================================================================
CONTAINS
!=======================================================================

  character(len=256) function get_filename (fulpath)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! returns filename given full pathname
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

! ------------------------ arguments --------------------------------
    character(len=*), intent(in)  :: fulpath !full pathname
! -------------------------------------------------------------------

! ------------------------ local variables --------------------------
    integer i               !loop index
    integer klen            !length of fulpath character string
! -------------------------------------------------------------------

    klen = len_trim(fulpath)
    do i = klen, 1, -1
       if (fulpath(i:i) == '/') go to 10
    end do
    i = 0
10  get_filename = fulpath(i+1:klen)
    
    return
  end function get_filename
  
!=======================================================================

  character(len=256) function set_filename (rem_dir, loc_fn)
  
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set remote full path filename
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

! ------------------------ arguments ------------------------------
    character(len=*), intent(in)  :: rem_dir !remote directory
    character(len=*), intent(in)  :: loc_fn  !local full path filename
! -----------------------------------------------------------------
    
! ------------------------ local variables ------------------------
    integer :: i   !integer
! -----------------------------------------------------------------

    set_filename = ' '
    do i = len_trim(loc_fn), 1, -1
       if (loc_fn(i:i)=='/') go to 10
    end do
    i = 0
10  set_filename = trim(rem_dir) // loc_fn(i+1:len_trim(loc_fn))
    
  end function set_filename

!=======================================================================

   subroutine getfil (fulpath, locfn, iflag)
     use LIS_logMod, only : LIS_logunit
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Obtain local copy of file
! First check current working directory
! Next check full pathname[fulpath] on disk
! Finally check full pathname[fulpath] on mass store
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
 
! ------------------------ arguments -----------------------------------
   character(len=*), intent(in)  :: fulpath !MSS or permanent disk full pathname
   character(len=*), intent(out) :: locfn   !output local file name
   integer, optional, intent(in) :: iflag   !0=>abort if file not found 1=>do not abort
! --------------------------------------------------------------------
 
! ------------------------ local variables ---------------------------
   integer i               !loop index
   integer klen            !length of fulpath character string
   integer ierr            !error status
   logical lexist          !true if local file exists
   character(len=256) text !mswrite command
! --------------------------------------------------------------------
 
 
! get local file name from full name: start at end. look for first "/"
 
   klen = len_trim(fulpath)
   do i = klen, 1, -1
      if (fulpath(i:i).eq.'/') go to 100
   end do
   i = 0
  100 locfn = fulpath(i+1:klen)
   if (len_trim(locfn) == 0) then
      write(LIS_logunit,*)'(GETFIL): local filename has zero length'
      call endrun
   else
      write(LIS_logunit,*)'(GETFIL): attempting to find local file ',          &
     &     trim(locfn)
   endif
 
! first check if file is in current working directory.
 
   inquire (file=locfn,exist=lexist)
   if (lexist) then
      write (LIS_logunit,*) '(GETFIL): using ',trim(locfn),                    &
     &     ' in current working directory'
      RETURN
   endif
 
! second check for full pathname on disk
 
   inquire(file=fulpath,exist=lexist)
   if (lexist) then
      locfn = trim(fulpath)
      write(LIS_logunit,*)'(GETFIL): using ',trim(fulpath)
      return
   endif
 
! finally check on mass store
 
   text='msread '//trim(locfn)//' '//trim(fulpath)
   call shell_cmd(text, ierr)
   if (ierr==0) then
      write(LIS_logunit,*)'(GETFIL): File ',trim(locfn),' read from MSS'
   else  ! all tries to get file have been unsuccessful
      write(LIS_logunit,*)'(GETFIL): failed cmd=',trim(text)
      if (present(iflag) .and. iflag==0) then
         call endrun
      else
         RETURN
      endif
   end if
 
   return
   end subroutine getfil
 
!=======================================================================
 
   subroutine putfil(locfn, mssfpn, pass, irt, lremov)
     use LIS_logMod, only : LIS_logunit
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Dispose to Mass Store only if nonzero retention period.
! 
! Method: 
! Put mswrite command in background for asynchronous behavior.
! The string put into 'cmd' below needs to be changed to 
! the appropriate archival command for the users system 
! if a shell command 'mswrite' does not exist.
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
 
!------------------------------Arguments--------------------------------
   character(len=*), intent(in) :: locfn   ! Local filename
   character(len=*), intent(in) :: mssfpn  ! Mass Store full pathname
   character(len=*), intent(in) :: pass    ! write password
   integer, intent(in) :: irt              ! Mass Store retention time
   logical, intent(in) :: lremov           ! true=>remove local file
!-----------------------------------------------------------------------
 
!---------------------------Local workspace-----------------------------
   character(len=256) cmd     ! Command string
   character(len=256) cmdtem  ! Temporary for command string
   character(len=  4) crt     ! Retention time as characters
   character(len= 16) wpass   ! Write password
   integer ier                ! error number
!-----------------------------------------------------------------------
 
 
   if (irt/=0) then
      wpass = ' '
      if (pass(1:1) /= ' ') wpass = ' -w ' // trim(pass)
      write (crt,'(i4)') irt
      write (cmd,'(100a)') 'mswrite ',' -t ',crt,trim(wpass),' ',&
           trim(locfn),' ',trim(mssfpn)
      if (lremov) then
         cmdtem = '('//trim(cmd)//'; /bin/rm '//trim(locfn)//' )&'
      else
         cmdtem = '('//trim(cmd)//' )&'
      end if
      write(LIS_logunit,*)'(PUTFIL): Issuing shell cmd:',trim(cmdtem)
      call shell_cmd(cmdtem, ier)
      if (ier /= 0) then
         write(LIS_logunit,*)'(PUTFIL): Error from shell cmd'
         call endrun
      end if
   endif

   return
   end subroutine putfil
 
!=======================================================================
 
   subroutine opnfil (locfn, iun, form)
     use LIS_logMod, only : LIS_logunit
!----------------------------------------------------------------------- 
! 
! Purpose: 
! open file locfn in unformatted or formatted form on unit iun
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
 
! ------------------------ input variables ---------------------------
   character(len=*), intent(in):: locfn  !file name
   integer, intent(in):: iun             !fortran unit number
   character(len=1), intent(in):: form   !file format: u = unformatted. f = formatted
! --------------------------------------------------------------------
 
! ------------------------ local variables ---------------------------
   integer ioe             !error return from fortran open
   character(len=11) ft    !format type: formatted. unformatted
! --------------------------------------------------------------------
 
   if (len_trim(locfn) == 0) then
      write(LIS_logunit,*)'(OPNFIL): local filename has zero length'
      call endrun
   endif
   if (form=='u' .or. form=='U') then
      ft = 'unformatted'
   else
      ft = 'formatted  '
   end if
   open (unit=iun,file=locfn,status='unknown',form=ft,iostat=ioe)
   if (ioe /= 0) then
      write(LIS_logunit,*)'(OPNFIL): failed to open file ',trim(locfn),        &
     &     ' on unit ',iun,' ierr=',ioe
      call endrun
   else
      write(LIS_logunit,*)'(OPNFIL): Successfully opened file ',trim(locfn),   &
     &     ' on unit= ',iun
   end if
 
   return
   end subroutine opnfil
 
!=======================================================================

  integer function getavu()
     use LIS_logMod, only : LIS_logunit
!----------------------------------------------------------------------- 
! 
! Purpose: 
! get next available Fortran unit number
!
! Method: 
! Get next available Fortran unit number itst. Set lsmiou(itst), in 
! lsmio common block, true. If coupled to CAM, use CAM function navu
! to get available unit number, in which case lsmiou is not needed.
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------

#if (defined COUP_CAM)
    use units     !CAM units module
#endif

! ------------------------ local variables ------------------------
    integer itst  !Fortran unit number
! -----------------------------------------------------------------

#if (defined COUP_CAM)
    getavu = getunit()
    RETURN
#else
    do itst = 1, 99
       if (.not.lsmiou(itst)) then
          getavu = itst
          lsmiou(itst) = .true.
          RETURN
       end if
    end do
    write (LIS_logunit,*) 'GETAVU error: ran out of Fortran unit numbers'
    call endrun
#endif
  end function getavu

!=======================================================================

  subroutine relavu (iunit)
     use LIS_logMod, only : LIS_logunit
!----------------------------------------------------------------------- 
! 
! Purpose: 
! close and release Fortran unit no longer in use
!
! Method: 
! Close and release Fortran unit number iunit. Set lsmiou(iunit) to 
! false. If coupled to cam, use cam function relunit to close/release 
! unit number.
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------

#if (defined COUP_CAM)
    use units     !CAM units module
#endif

! ------------------------ arguments ------------------------------
    integer, intent(in) :: iunit    !Fortran unit number
! -----------------------------------------------------------------

#if (defined COUP_CAM)
    close(iunit)
    call freeunit(iunit)
#else
    if (.not.lsmiou(iunit)) then
       write (LIS_logunit,*) 'RELAVU eror: unit ',iunit,' is not flagged as in use'
       call endrun
    end if
    if (iunit<1 .or. iunit>99) then
       write (LIS_logunit,*) 'RELAVU error: attempt to return out of range unit'
       call endrun
    end if
    close(iunit)
    lsmiou(iunit) = .false.
#endif
    return
  end subroutine relavu

!=======================================================================

   subroutine shell_cmd(text, ier)
 
! ------------------------ arguments -----------------------------------
   character(len=*), intent(in) :: text
   integer         , intent(out):: ier
! ----------------------------------------------------------------------
 
! ------------------------ local variables -----------------------------
#if ( defined CRAY )
   integer, external :: ishell ! System routine, execute shell command
#elif (!defined AIX)
   integer, external :: system ! System routine, execute shell command
#endif
! ----------------------------------------------------------------------
 
#if ( defined CRAY )
   ier = ishell(trim(text))
#elif ( defined AIX )
   call system(trim(text), ier)
   ier = 0               ! Set ier to zero
#elif (!defined CRAY) && (!defined AIX)
   ier = system(trim(text))
#endif
 
   return
   end subroutine shell_cmd

!=======================================================================

end module fileutils

