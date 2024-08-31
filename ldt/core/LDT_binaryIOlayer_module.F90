!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_binaryIOlayer_module

!BOP
!
! !MODULE: LDT_binaryIOlayer_module
!
! !DESCRIPTION:
!  The routines in this module contain methods to read and handle
!  different binary format type files.
!
!  Origins of each subroutine vary with organization and source.
!  
!
! !REVISION HISTORY:
  implicit none

  PRIVATE

  public :: LDT_read_bilfile
  public :: LDT_read_bsqfile
!  public :: LDT_read_xmrgfile

!EOP
  contains

!BOP
!
! !ROUTINE: LDT_read_bilfile
! \label{LDT_read_bilfile}
!
! !DESCRIPTION:
! !MODULE:  LDT_read_bilfile
!
! !DESCRIPTION:
!  This routine contains methods to read and handle
!   band-interleaved (BIL) binary data array.
!
! !REVISION HISTORY:
!  28 Oct 2013: KR Arsenault : Specification in LDT.
!
! !INTERFACE:
 SUBROUTINE LDT_read_bilfile( array_len, list )
! EOP

  implicit none

!- Inputs::
   integer, intent (in)                         :: array_len
   real, dimension (array_len), intent(inout)   :: list

 end subroutine LDT_read_bilfile

! =========================================================

!BOP
!
! !ROUTINE: LDT_read_bsqfile
! \label{LDT_read_bsqfile}
!
! !DESCRIPTION:
! !MODULE:  LDT_read_bsqfile
!
! !DESCRIPTION:
!     This routine contains methods to read and unpack
!     band-sequential (BSQ) binary data array.
!
!     This subroutine extracts a specified window (rectangular subset)
!     from a band-sequential (BSQ) binary data array and places it into
!     an INTEGER*2 array.  In a band-sequential array, all data for the
!     first band come first, followed by all data for the second band,
!     and so on until the last band is complete.  Within each band,
!     Within each band, the data are organized by rows
!
!     The 'bands' of the array may represent spectral bands, geographic
!     layers, successive time intervals, or the like.
!
!     For multi-byte data values, this subroutine assumes that the input
!     array is in "big endian" format; i.e., the most-significant byte
!     of each value precedes the least-significant byte.  If the host
!     computer being used to read the data uses "little endian" format,
!     the compile-time parameter BIGEND should be changed -- see the
!     PARAMETER statements below.
!
!  ORIGINAL USAGE:
!     CALL READBSQ (lun, inband, inrow, incol, inbyte,
!    +   iband, irow, icol, nband, nrow, ncol, outarray)
!
!  Arguments:
!     INTEGER*4     lun         Logical unit number to be used for
!                               accessing the input data.
!     INTEGER*4     inband      Number of bands in input data array; for
!                               two-dimensional array, inband = 1.
!     INTEGER*4     inrow       Number of rows in input data array.
!     INTEGER*4     incol       Number of columns per row in input data
!                               array.
!     INTEGER*4     inbyte      Number of bytes per input data value.
!                               Allowed values are 1 and 2.
!     INTEGER*4     iband       Initial band for which data are desired.
!     INTEGER*4     irow        Initial row for which data are desired:
!                               rows are numbered from top to bottom
!                               (north to south), beginning with row 1.
!     INTEGER*4     icol        Initial column for which data are
!                               desired:  columns are numbered from left
!                               to right (west to east), beginning with
!                               column 1.
!     INTEGER*4     nband       Number of consecutive bands to be
!                               extracted from input data array
!     INTEGER*4     nrow        Number of consecutive rows to be
!                               extracted from each band of data array.
!     INTEGER*4     ncol        Number of consecutive columns to be
!                               extracted from each row of data array.
!     INTEGER*2     outarray(ncol,nrow,nband)
!                               Array to receive the data from input
!                               file.  If nband = 1, the array may be
!                               declared as outarray(ncol,nrow).
!
!  Compile-time parameters:
!     BIGEND          Byte offset in memory from start of an INTEGER*2
!                     value to the most significant byte of the value.
!                     Set to 1 if host computer stores the least
!                     significant byte of integer data values at a
!                     lower memory address than the most significant
!                     byte ("little-endian"); set to 0 otherwise.
!     MAXINBYTE       Maximum number of bytes per row in input data
!                     array, MAXINBYTE >= in_col.
!
!  Addtional variables and arrays used:
!     band            Current input band.
!     bytebuf         Buffer for moving bytes into integers.
!     col             Current input column.
!     inbuf           Buffer to receive one record of input data.
!                     declared as type CHARACTER.
!     intbuf          Buffer for moving bytes into integers, declared
!                     INTEGER*2, equivalenced to  bytebuf  array.
!     row             Current input row.
!
! !REVISION HISTORY:
!     Jan 1996: R.A. White, PSU/ESSC: Initial version.
!  28 Oct 2014: KR Arsenault: Specification in LDT.
!
! !INTERFACE:
 SUBROUTINE LDT_read_bsqfile( &
            lun, inband, inrow, incol, inbyte, &
            iband, irow, icol, nband, nrow, ncol, outarray)
!
! EOP
!
! USES:
  use LDT_logMod, only : LDT_logunit, LDT_endrun

  implicit none

  integer*4, intent(in) :: lun
  integer*4, intent(in) :: inband
  integer*4, intent(in) :: inrow
  integer*4, intent(in) :: incol
  integer*4, intent(in) :: inbyte  ! Number of bytes per input data value;
                                   !  Allowed values are 1 and 2.
  integer*4, intent(in) :: iband
  integer*4, intent(in) :: irow
  integer*4, intent(in) :: icol
  integer*4, intent(in) :: nband
  integer*4, intent(in) :: nrow
  integer*4, intent(in) :: ncol
!  integer*2, intent(inout):: outarray(ncol, 1, nband)
  integer*2, intent(inout):: outarray(ncol, nrow, nband)
!(VK) This is original definition:    INTEGER*2 outarray(ncol, nrow, nband)

!- Local parameters:

! -- Use BIGEND = 0 for big-endian hosts, such as Sun SPARC
!   integer, parameter :: BIGEND = 0

! -- Use BIGEND = 1 for little-endian hosts, such as VAX and Intel
!     80x86-compatible PCs
   integer, parameter :: BIGEND = 1
   integer, parameter :: MAXINBYTE = 16384

! -- In the following format statement, the repeat count must match the
!     value of MAXINBYTE
 1000 FORMAT(16384A1)

  integer     :: irec, i
  integer*4   :: band, row, col
  integer*2   :: intbuf
  character*1 :: bytebuf(2)
  character*1 :: inbuf(2*MAXINBYTE)

  equivalence(intbuf, bytebuf)

! ---------------------------------------------------------

  intbuf = 0

  do band = iband, iband+nband-1
     do row = irow, irow+nrow-1

        irec= row+(band-1)*inrow
        read( lun, REC=row+(band-1)*inrow, FMT=1000) &
            ( inbuf(i), i=1,incol*inbyte )

        do col = icol, icol+ncol-1
!
!  Unpack 1-byte data to 2-byte integers
!
           IF(inbyte .GT. 1) GOTO 20

           bytebuf(2-BIGEND) = inbuf(col)
           GOTO 40
!
!  Extract two-byte data from input buffer           
!
 20        CONTINUE
           bytebuf(BIGEND + 1) = inbuf(2*col-1)
           bytebuf(2 - BIGEND) = inbuf(2*col)

 40        CONTINUE
           outarray(col-icol+1, row-irow+1, band-iband+1) = intbuf

        end do
     end do
  end do
  return

 end subroutine LDT_read_bsqfile

! =========================================================

! Other subroutines ...


end module LDT_binaryIOlayer_module
