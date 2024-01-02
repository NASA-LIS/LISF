!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module fbil_module

! Simple, generic ESRI bil I/O library for reading and writing files
!  in the binary interleaved data file format.
!
! Written during FLDAS-WRSI development for support of that project 
! for FEWSNET/USGS, Jim Verdin, Christa Peters-Lidard, et al. 
! by Brad Wind of ESSIC/UMD.
!
! * FEWSNET_bil I/O library Disclaimer:
!
! This file reading infrastructure handles ESRI bil 'binary interleaved' files.
! Specifically, those of the format utilized by FEWSNET, and more specifically by 
! the USGS Santa Barbara application "GeoWRSI."  The following code began its 
! life in support of the needs of translating the WRSI model work to Fortran from
! Visual Basic, however every attempt was made to permit this I/O framework to be 
! flexible if not robust enough to accommodate reading of files of the ESRI bil 
! format generally, with minimal modification, for any purpose.
!
! Accordingly, successful support for files of the formats .hdr, .blw, and .bil 
! has been tested for the WRSI project.  Untested functionality exists herein for 
! files of the 'Float' format.  Similarly, untested code exists for 1) 'nbits=gBIL_BIT' 
! bit reading and for 2) 'fillSpatialGaps=.true.' as well as for resampling without 
! correcting for spatial overflow, all features of which were provided for 
! in the library source code upon which this code was based but none of which came 
! complete with both data file examples and use cases in the context of the WRSI 
! development effort sufficient to declare these features well-tested by any 
! stretch.  Such scenarios may well exist even in the context of ongoing WRSI work, 
! however these processing paths remain largely untested, and are likely so in 
! the original library upon which this code is based as well.  As in that library, 
! so-called .stx statistics files are not natively supported at time of this 
! writing nor are any other related extraneous data header or file types other 
! than those explicitly mentioned, once again lacking both for examples and project 
! requirements to now.  However, addition of these features is encouraged in the 
! spirit of an efficient, consistent ESRI .bil I/O framework.
!
! Endian-resolution of the reader code was successfully tested.  The in-code 
! compiler-independent endian resolution permits the code to detect the endian 
! of the data file and the endian of the system and properly interpret read data.  
! Accordingly, endian tests were completed both on little endian Linux platforms 
! using gfortran and ifort; and on the big endian Power PC MacOSX platform 
! using Absoft Of course, while 'm' and 'msbfirst' and 'motorola' are the 
! keywords used for case-insensitive comparison of the standard ESRI header 
! field 'byteorder' this keyword list represents only those keywords used for 
! endian considerations in practice by FEWSNET in the WRSI project files provided
! as well as those keywords resulting from a brief Google search and may by 
! no means be complete.  The list is intended to be added-to as needs arise, as 
! is the case with any other project application -specific encoding or 
! interpretation of header fields that may prove necessary.
!
! Because this infrastructure was intended first and foremost to provide reading 
! support, the 'O' of the I/O herein is as of this writing light compared with 
! the 'I' or input functionality.  Enhancement of standardized write capabilities 
! is strongly encouraged.
!
! Within the code below, 150 lines of code remain completely untested.  
!  The lines are those pertaining to the two features enumerated above.  
!  These features represent significatly different file processing paths, which 
!  have been marked-off with the precompiler flag 'BRAVE_'.  
!  Pound-define BRAVE_ (#define BRAVE_) to use them but please test them first.  
!  While reasonable validation was performed for the rest of the code, especially
!  nbits=gBIL_UNSIGNEDBYTE and nbits=gBIL_REAL4, none was done for the BRAVE_ sections.
!
! +Brad Wind(BW) March, 2011
! ___________________________________________________________________
!
! FEWSNET BIL I/O library I/O precompiler flag directives
!#include "LIS_misc.h"
#define HEADER_FILE_READ_METHOD_1_
#define ASSUME_BIG_ENDIAN_FEWSNET_BIL_READ_

! Define the first of the following two if the code is built w/ ifort -convert big_endian
! which instructs byte order interpretation as big_endian regardless of system byteorder.
!#define ASSUME_BIG_ENDIAN_FEWSNET_BIL_READ_ 1
#ifndef ASSUME_BIG_ENDIAN_FEWSNET_BIL_READ_
!#define ASSUME_LITTLE_ENDIAN_FEWSNET_BIL_READ_ 1
#endif

! On ifort and Absoft/xlf this was the supported method for reading 
! an ascii file character-by-character
!#define HEADER_FILE_READ_METHOD_1_ 1 
! On gfortran it was necessary to use this alternative method instead
!#define HEADER_FILE_READ_METHOD_2_ 0

! Represent diagnostics for file I/O for the model code:
#define DIAGNOSTICS_ 1

! Reverse data in the y dimension. Some .bil data goes from North to South 
! whereas LIS and GrADS go from S. to N.
#define YREV_ 1
! As of this writing (10/2011) this YREV_ option works for a single band

! "Brave" code sections never tested (See FEWSNET_bil I/O library Disclaimer)
!#define BRAVE_
! ________________________________________________________________________

! !USES:
  use LIS_logMod,       only : LIS_logunit, LIS_endrun

 implicit none

 type charN
   ! This struct is used so that we can easily change 
   !  the # characters allocated per str
   character*300, pointer  :: str
 end type charN

! FORTRAN_INTRINSIC_SUCCESS_VALUES_
  integer, parameter :: gALLOCATE_SUCCESS = 0
  integer, parameter :: gFILEOPEN_SUCCESS = 0
  integer, parameter :: gFGETC_SUCCESS = 0
  integer, parameter :: gFILEREAD_SUCCESS = 0

! This struct is used so that we can easily access coordinate space information
 type coord_spc
  ! Double *ulxmap = Math.Round(minLon, 5) 'AnalysisMinLL is actually the upper-left corner, 
  !  the origin of the program coord system
   real*8    :: minLon 
   real*8    :: maxLon
   real*8    :: pixLon   ! *xdim=pixLon
   real*8    :: minLat
   real*8    :: maxLat   ! *ulymap = Math.Round(maxLat, 5)
   real*8    :: pixLat   ! Double *ydim=pixLat
  ! *These assignments and comment noted from GetHeaderFromRegionSet function in GeoWRSI

  ! The following derived from the above.  Useful to compute once; store for later use
  !  nCols = CInt4(((.maxLon - .minLon) / .PixLon) + 1) ' lon associated with x +BW
  !  nRows = CInt4(((.maxLat - .minLat) / .PixLat) + 1) ' lat associated with y +BW
  ! These are used in memory-buffering the data file. For subsquent model-array sizing, 
  ! see gNcolsModelArr, gNrowsModelArr or the like elsewhere which may be defined 
  ! differently.
   integer*4 :: maxX ! Integer ! this, immediately set as nCols upon config file read
   integer*4 :: maxY ! Integer ! this, immediately set as nRows upon config file read

 ! These last two vars set & used exclusively by the 'gCoords' and for two purposes
 ! to size arrays: ((gCoords%maxX-gCoords%minX) +1), ((gCoords%maxY-gCoords%minY) +1)
 ! to loop: xCols=gCoords%minX,gCoords%maxX,1;yRows=gCoords%minY,gCoords%maxY,1
   integer*4 :: minX ! Integer
   integer*4 :: minY ! Integer
 end type coord_spc


! BIL_PARAMETERS_
  logical,   parameter   :: gReadAllFilesAsThoughWrittenRowMajor = .true.
  character*3, parameter :: gBIL_HEADERFILENAMEEXT = 'hdr'
  character*3, parameter :: gBIL_BLWFILENAMEEXT = 'blw'
  integer*4,   parameter :: gBIL_UNSIGNEDBYTE = 8
  integer*4,   parameter :: gBIL_BIT = 1
  integer*4,   parameter :: gBIL_INT2 = 16
  integer*4,   parameter :: gBIL_INT4 = 32
  integer*4,   parameter :: gBIL_REAL4 = 33

!! Following two variables represent modalities encountered when reading tokens in BIL 
!! header files & must merely be different from one another in ways shown
  integer*4,   parameter :: gBIL_identifierLabel = 1
  integer*4,   parameter :: gBIL_idVal = 2 
  integer*4,   parameter :: gBIL_MAX_HEADERFILE_LINELENGTH = 300

! The following two variables from String_Utility.f90
  CHARACTER( * ), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz' 
  CHARACTER( * ), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 

! The following four variables are used variously to establish correct reading of hdr
  integer*4, parameter :: gInitValNcols = 0
  real*8,    parameter :: gInitValyllcorner = (0-1)
  real*4,    parameter :: gInitValXdim = (0-1)
  integer*4, parameter :: gInitValnbands = (0-1)

 type FEWSNET_bil__header
   ! This type was originally GeoWRSI's GIS I/O library, Public Structure BilHeader.
   ! Here, also supports alternative, similar 'Float' file format header info.

 ! Calculated during header file read:
   integer*4                            :: numEntriesRead ! Added by BW

 ! Read, or computed during header file read, directly from header file:
   character(len=gBIL_MAX_HEADERFILE_LINELENGTH) &
                                        :: byteorder ! String
   character(len=gBIL_MAX_HEADERFILE_LINELENGTH) &
                                        :: layout    ! String
   integer*4                            :: nbits     ! Integer
   integer*4                            :: skipbytes ! Integer
   real*4                               :: xdim      ! Single
   real*4                               :: ydim      ! Single
   integer*4                            :: ncols     ! Integer
   integer*4                            :: nrows     ! Integer
   integer*4                            :: nbands    ! Integer
   real*8                               :: ulxmap    ! Double
   real*8                               :: ulymap    ! Double
   integer*4                            :: missingvalue  ! Integer
   integer*8                            :: bandrowbytes  ! Integer
   integer*8                            :: totalrowbytes ! Integer
   integer*8                            :: bandgapbytes  ! Integer 
   real*8                               :: yllcorner ! Double ! I added this +BW

 ! Useful for read process of, and hard-wired for, a particular header file (sub)type:
   integer*4                            :: minimumNumRequiredEntries ! Added by BW

 ! For handling potential mis-match between read- and target- coordinate systems 
 !  and used for setting derived variables below:
   logical                              :: refCoordsIsNull
   type(coord_spc), pointer             :: refCoords
   logical :: abortFileReadIfHeaderDataRegionDoesntSubsumeAnalysisRegion

 ! Derived variables, evaluated secondarily post-read from read-in vars
   integer*4                            :: minX
   integer*4                            :: maxX
   integer*4                            :: minY
   integer*4                            :: maxY
 ! The same, corrected for spatial overflow
   integer*4                            :: minX_corr_spatl_overflow
   integer*4                            :: maxX_corr_spatl_overflow
   integer*4                            :: minY_corr_spatl_overflow
   integer*4                            :: maxY_corr_spatl_overflow

 ! These, the original read-in variables, adjusted per derived variables
 ! This group adjusted per derived variables not corrected for spatial overflow
   real*8                               :: ulxmap_adj
   real*8                               :: ulymap_adj
   integer*4                            :: ncols_adj
   integer*4                            :: nrows_adj
 ! These adjusted using derived variables corrected for spatial overflow
   real*8                               :: ulxmap_corr_spatl_overflow
   real*8                               :: ulymap_corr_spatl_overflow
   integer*4                            :: ncols_corr_spatl_overflow
   integer*4                            :: nrows_corr_spatl_overflow

 end type FEWSNET_bil__header

  integer*4, parameter :: gDEFAULT_SIGDIGSRTOFDECPT = 5
  integer*4, parameter :: gDEFAULT_SIGDIGSTOTAL = 5

  real*8,    parameter :: gCOORD_SPECIAL_VALUE_1 = (0-9999)
  integer*2, parameter :: gSPATIAL_GAP_FILL_VALUE_1 = 255
  integer*4, parameter :: gDATA_SPECIAL_VALUE_1 = (0-9999)


!================begin application-layer input file reading================

! The following the user-supplied spatial bounding coordinates, and pixel width & height
! I.e. the Model Spatial Domain
  type(coord_spc), pointer   :: gCoords

!! Following four variables represent 2d file-read modalities in FEWSNET_bil &
!!  must merely be different from one another in ways shown:
  integer*4, parameter :: gUSE_CURRENT_TIMESTEP_PPT                 = 1
  integer*4, parameter :: gUSE_LONG_TERM_AVERAGE_CLIMATOLOGICAL_PPT = 3
  integer*4, parameter :: gUSE_CURRENT_TIMESTEP_PET                 = 2
  integer*4, parameter :: gUSE_LONG_TERM_AVERAGE_CLIMATOLOGICAL_PET = 3

!================end application-layer input file reading================

CONTAINS


 subroutine initialize_FEWSNET_bil__header(headerInfo, &
   initialize_float_hdrFileFormat)

   type(FEWSNET_bil__header), intent(inout) :: headerInfo
   logical, intent(in),            optional :: initialize_float_hdrFileFormat

   logical                                  :: init_flt_hdr
! ___________________________________________________________________

   init_flt_hdr = .false.
   if ( present(initialize_float_hdrFileFormat) ) then
      init_flt_hdr = initialize_float_hdrFileFormat
   endif

   headerInfo%numEntriesRead = 0
   write(headerInfo%byteorder, '(a)') ''
   write(headerInfo%layout,    '(a)') ''
   if( init_flt_hdr .eqv. .true.) write(headerInfo%layout, '(a)') 'BIL'
   headerInfo%nbits = (0-1)
   if( init_flt_hdr .eqv. .true.) headerInfo%nbits = gBIL_REAL4
   headerInfo%skipbytes = (0-1)
   headerInfo%xdim = gInitValXdim
   headerInfo%ydim = (0-1)

 ! Original GeoWRSI/GIS codes based decisions on ncols remaining '0' +BW
   headerInfo%ncols = gInitValNcols
   headerInfo%nrows = (0-1)
   headerInfo%nbands = gInitValnbands 
   headerInfo%ulxmap = (0-1)
   headerInfo%ulymap = (0-1)
   headerInfo%missingvalue = (0-1)
   headerInfo%bandrowbytes = (0-1)
   headerInfo%totalrowbytes = (0-1)
   headerInfo%bandgapbytes = (0-1)
   headerInfo%yllcorner = gInitValyllcorner
   headerInfo%minimumNumRequiredEntries = (0-1)

   if( init_flt_hdr .eqv. .true. ) headerInfo%minimumNumRequiredEntries = 5
   ! Unique are:
   !  headerInfo%refCoordsIsNull
   !  headerInfo%refCoords
   !  headerInfo%abortFileReadIfHeaderDataRegionDoesntSubsumeAnalysisRegion
   ! in that they are never initialized; rather 
   ! they are presumed to have been set elsewhere outside of the read process
   headerInfo%minX = (0-1)
   headerInfo%maxX = (0-1)
   headerInfo%minY = (0-1)
   headerInfo%maxY = (0-1)
   headerInfo%minX_corr_spatl_overflow = (0-1)
   headerInfo%maxX_corr_spatl_overflow = (0-1)
   headerInfo%minY_corr_spatl_overflow = (0-1)
   headerInfo%maxY_corr_spatl_overflow = (0-1)
   headerInfo%ulxmap_adj = (0-1)
   headerInfo%ulymap_adj = (0-1)
   headerInfo%ncols_adj = (0-1)
   headerInfo%nrows_adj = (0-1)
   headerInfo%ulxmap_corr_spatl_overflow = (0-1)
   headerInfo%ulymap_corr_spatl_overflow = (0-1)
   headerInfo%ncols_corr_spatl_overflow = (0-1)
   headerInfo%nrows_corr_spatl_overflow = (0-1)


end subroutine initialize_FEWSNET_bil__header


function fileReadFailedBasedOnValXremainingForVarY(&
     int4VarY, int4ValX,   &
     real4VarY, real4ValX, &
     real8VarY, real8ValX)

   logical :: fileReadFailedBasedOnValXremainingForVarY

   integer*4, intent(in), optional :: int4VarY
   integer*4, intent(in), optional :: int4ValX
   real*4, intent(in),    optional :: real4VarY
   real*4, intent(in),    optional :: real4ValX
   real*8, intent(in),    optional :: real8VarY
   real*8, intent(in),    optional :: real8ValX

   fileReadFailedBasedOnValXremainingForVarY = .false.

   if(present(int4VarY) .and. present(int4ValX)) then
      if(int4VarY == int4ValX)   fileReadFailedBasedOnValXremainingForVarY = .true.
   endif
   if(present(real4VarY) .and. present(real4ValX)) then
      if(real4VarY == real4ValX) fileReadFailedBasedOnValXremainingForVarY = .true.
   endif
   if(present(real8VarY) .and. present(real8ValX)) then
      if(real8VarY == real8ValX) fileReadFailedBasedOnValXremainingForVarY = .true.
   endif

end function fileReadFailedBasedonValXremainingForVarY


function checkFileExistence(filename)

   integer*4                 :: checkfileExistence

   character(*), intent(in)  :: filename
   logical                   :: exist
   character(len=11)         :: readability


   checkFileExistence = 0

   inquire( file = filename,    &
            exist= exist,       &
            read = readability)
   if (.not. exist) then
      checkFileExistence = 1
   endif
  if (readability == 'NO') then
     checkFileExistence = 2
  endif

end function checkFileExistence


function getFreeFileHandle()

   integer       :: getFreeFileHandle

   logical       :: isOpen

   do getFreeFileHandle=100, 1000, 1
      inquire (unit=getFreeFileHandle, opened=isOpen)
      if (.not.(isOpen)) return
   end do

   getFreeFileHandle = (0-1)

end function getFreeFileHandle


FUNCTION StrLowCase ( Input_String ) RESULT ( Output_String ) 
  ! This function from String_Utility.f90
  ! -- Argument and result 
  CHARACTER( * ), INTENT( IN )     :: Input_String 
  CHARACTER( LEN( Input_String ) ) :: Output_String 
  ! -- Local variables 
  INTEGER :: i, n 

  ! -- Copy input string 
  Output_String = Input_String 
  ! -- Loop over string elements 
  DO i = 1, LEN( Output_String ) 
    ! -- Find location of letter in upper case constant string 
    n = INDEX( UPPER_CASE, Output_String( i:i ) ) 
    ! -- If current substring is an upper case letter, make it lower case 
    IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n ) 
  END DO 

END FUNCTION StrLowCase


! The following must be done prior to calling this subroutine
!  headerInfo%refCoords => reference_coordinate_system
!  headerInfo%refCoordsIsNull = .true. or .false. depending on whether its set
!  - set headerInfo%abortFileReadIfHeaderDataRegionDoesntSubsumeAnalysisRegion

subroutine readFEWSNET_bil_or_float_headerFile(mapFilename, headerInfo)

   type(charN), intent(in)                      :: mapFilename
   type(FEWSNET_bil__header), intent(inout)     :: headerInfo


   call readFEWSNET_bil_headerFile(mapFilename, headerInfo)

   ! In original GIS lib code where FEWSNET_bil and float header file reading
   ! codes were not united, the following check was utilized to establish a 
   ! failed bil header file read and subsequent go-ahead for 'float' header read  
   !if(fileReadFailedBasedOnValXremainingForVarY(&
   !   real4VarY=headerInfo%xdim, real4ValX=gInitValXdim) .eqv. .true.) then

   ! Now that the reading of these two header file types was unified, the check is:
   if(.not. (fileReadFailedBasedOnValXremainingForVarY(&
      real8VarY=headerInfo%yllcorner, real8ValX=gInitValyllcorner) .eqv. .true.)) then
 
     ! If the file was not in the expected format, try another format:
      call readFEWSNET_bil_headerFile(mapFilename, headerInfo, &
                  read_float_hdrFileFormat=.true.)

      if(fileReadFailedBasedOnValXremainingForVarY(&
         real4VarY=headerInfo%xdim, real4ValX=gInitValXdim) .eqv. .true.) then
        ! If you can't read even the other format, then exit the subroutine.
        ! Error: Could not read the Bil header file ...
         call initialize_FEWSNET_bil__header(headerInfo)
      endif

   endif
!    call print_header_info(headerInfo)    ! KRA

end subroutine readFEWSNET_bil_or_float_headerFile


subroutine readFEWSNET_bil_headerFile(mapFilename, headerInfo, read_float_hdrFileFormat)

! This subroutine assumes:
! * mapFilename is just the root filename with no filename extension
!   .bil, .hdr, .blw will be concatenated for input read & output write
!   Therefore, in LIS configuration file specification and others, specify root filenames only
! * There are no spaces in the filename or path lest special handling is added for such
! * The following was done prior to calling this subroutine:
!     headerInfo%refCoords => reference_coordinate_system
!     headerInfo%refCoordsIsNull = .true. or .false. depending on whether its set
!   Set headerInfo%abortFileReadIfHeaderDataRegionDoesntSubsumeAnalysisRegion

   type(charN), intent(in)                  :: mapFilename
   type(FEWSNET_bil__header), intent(inout) :: headerInfo
   logical, intent(in),            optional :: read_float_hdrFileFormat

   logical                                  :: initialize_fltHdrFileFormat
   type(charN)                              :: hdrFilename
   integer                                  :: fileHandle
   integer                                  :: pos
   integer                                  :: filestatus
   type(charN)                              :: blwFilename
   real*4                                   :: tmpReal4
!____________________________________________________________________________

   filestatus = -1

   initialize_fltHdrFileFormat = .false.
   if (present(read_float_hdrFileFormat) ) then
      initialize_fltHdrFileFormat = read_float_hdrFileFormat
   endif

   call initialize_FEWSNET_bil__header(headerInfo, &
      initialize_float_hdrFileFormat=initialize_fltHdrFileFormat)

   if (mapFilename%str == "") return

   allocate(hdrFilename%str)
   hdrFilename%str = trim(mapFilename%str) // '.' // gBIL_HEADERFILENAMEEXT

   allocate(blwFilename%str)
   blwFilename%str = trim(mapFilename%str) // '.' // gBIL_BLWFILENAMEEXT

!- Read BIL header:
!   print *, "(bilheader_reader): "  ! KRA

   if (checkFileExistence(trim(hdrFilename%str)) /= 0) return ! File does not exist

   fileHandle = getFreeFileHandle()
   if(fileHandle == (0-1) ) return   ! Couldn't obtain a free file handle

#ifdef HEADER_FILE_READ_METHOD_1_
   open (unit=fileHandle, file=trim(hdrFilename%str), access='stream', &
         iostat=filestatus)
#else
#ifdef HEADER_FILE_READ_METHOD_2_
   open (unit=fileHandle, file=trim(hdrFilename%str), iostat=filestatus)
#endif
#endif

   if(filestatus /= gFILEOPEN_SUCCESS) then ! File couldn't be opened
      return
   endif

   pos=1
   call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
              file_handle=fileHandle, pos=pos, status=filestatus, &
              headerInfo=headerInfo)

   if ( present(read_float_hdrFileFormat) ) then
      if ( read_float_hdrFileFormat .eqv. .true. ) then

         if(headerInfo%byteorder /= '') then
            if(trim(StrLowCase(headerInfo%byteorder)) == 'lsbfirst') then
              write(headerInfo%byteorder, '(a)') 'I'
            else
              write(headerInfo%byteorder, '(a)') 'M'
            endif
         endif
         if( (headerInfo%nrows /= (0-1)) .and. &
             (headerInfo%ydim /= (0-1))  .and. &
             (headerInfo%yllcorner /= (0-1)) ) then
            headerInfo%ulymap = headerInfo%yllcorner + &
                                (headerInfo%nrows - 1) * headerInfo%ydim
            headerInfo%numEntriesRead = headerInfo%numEntriesRead + 1
         endif
      endif
   endif
   close(fileHandle)

!-- Read blw file, if present: --
   fileHandle = getFreeFileHandle()
   if ( (checkFileExistence(trim(blwFilename%str)) == 0) .and. &
        (fileHandle /= (0-1)) ) then

      ! The following is adapted from the ArcView help file
      ! When the blw file is present, ArcView does a 6-parameter affine 
      !  transformation from image coordinates to world coordinates as follows:
      !  x1 = Ax + By + C
      !  y1 = Dx + Ey + F
      ! where
      !  x1 = world x-coordinate of pixel
      !  y1 =  world y-coordinate of pixel
      !   x = column number of pixel in image (starting at 0 in UL corner)
      !   y = row number of pixel in image (starting at 0 in UL corner)
      !   A = size of pixel in world units in x direction
      !B, D = rotation parameters
      !C, F = x,y world coordinates of center of upper-left pixel (for translation)
      !   E = -1 * size of pixel in world units in y direction
      ! The parameters are stored in the world file with the following order:
      ! 0.1 - A
      ! 0.0 - D
      ! 0.0 - B
      ! -0.1 - E
      !-20.0 - C
      ! 40.0 - F

#ifdef HEADER_FILE_READ_METHOD_1_
      open (unit=fileHandle, file=trim(blwFilename%str), access='stream', &
            iostat=filestatus)
#else
#ifdef HEADER_FILE_READ_METHOD_2_
      open (unit=fileHandle, file=trim(blwFilename%str), iostat=filestatus)
#endif
#endif
      if(filestatus /= gFILEOPEN_SUCCESS) then ! File couldn't be opened
         return
      endif

      pos=1
      call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
         file_handle=fileHandle, pos=pos, status=filestatus,  &
         reading_token_of_this_class=gBIL_idVal,          &
         real4Val=headerInfo%xdim,                        &
         no_round=.true.)  ! Explicit rounding wasn't done for any of these here +BW

      !! The next two lines were and are effectively skipped
         call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
            file_handle=fileHandle, pos=pos, status=filestatus, &
            reading_token_of_this_class=gBIL_idVal,         &
            real4Val=tmpReal4, no_round=.true.)
         call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
            file_handle=fileHandle, pos=pos, status=filestatus, &
            reading_token_of_this_class=gBIL_idVal,         &
            real4Val=tmpReal4, no_round=.true.)

         call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
            file_handle=fileHandle, pos=pos, status=filestatus, &
            reading_token_of_this_class=gBIL_idVal,         &
            real4Val=headerInfo%ydim, no_round=.true.)
         headerInfo%ydim = (0-1) * headerInfo%ydim

         call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
            file_handle=fileHandle, pos=pos, status=filestatus, &
            reading_token_of_this_class=gBIL_idVal,         &
            real4Val=tmpReal4, no_round=.true.)
         headerInfo%ulxmap = tmpReal4 ! ulxmap and ulymap were read as 
                                      ! aline = LineInput(filehandle1)
                                      ! ahdr.ulxmap = CSng(aline)
         call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
            file_handle=fileHandle, pos=pos, status=filestatus, &
            reading_token_of_this_class=gBIL_idVal,         &
            real4Val=tmpReal4, no_round=.true.)
         headerInfo%ulymap = tmpReal4
      close(fileHandle)

   endif   ! End BLW file read

! --- END BLW ---

   ! if read_float_hdrFileFormat .eqv. .false. then
   !
   ! Regular first-order bil read.  A sanity check here 
   ! should be based on some other criterion than # entries read
   ! lots of good bil files have less than 10 entries..
   !
   ! ..However, with the alternative 'Float' format, the following was done..
   if ( present(read_float_hdrFileFormat) ) then
      if (read_float_hdrFileFormat .eqv. .true.) then
         if (headerInfo%numEntriesRead < &
             headerInfo%minimumNumRequiredEntries) then
            ! Not enough entries could be read from the header file:
            call initialize_FEWSNET_bil__header(headerInfo, &
               initialize_float_hdrFileFormat=initialize_fltHdrFileFormat)
         endif
      endif
   endif

 ! The header file read has gone ok to this point ...
 !  Do some calculations of derived header file data members immediately,
 !  so that they are available to any and all further processing; 
 !  quit via reinitializing header info upon failure to perform adjustment - 
 !  this is effectively the same as returning success or failure b/c subsequent 
 !  checks examine whether header info remains in initialized state to proceed.

   if(headerInfo%refCoordsIsNull .neqv. .true.) then

      if( adjust_bil_headerFile_variables_post_read(headerInfo) .neqv. .true.) then  ! FALSE ...
         call initialize_FEWSNET_bil__header(headerInfo, &
               initialize_float_hdrFileFormat=initialize_fltHdrFileFormat)
        if( headerInfo%ncols == 0 ) then
           print *, " param header file has ncols reset to 0 ..."
           stop 
        else
           print *, " "
           print *, " "
        end if
      end if

   endif

end subroutine readFEWSNET_bil_headerFile


recursive subroutine FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
      file_handle, status, pos, headerInfo, numEntriesRead,  &
      reading_token_of_this_class, charStrVal, int4Val,      & 
      real4Val, real8Val, int8Val, no_round   )

! This subroutine represents recursive read of actual characters of the header file. 
! Spacing patterns therein, we are told by the I/O library how they are to be
! handled as though arbitrary.

   integer*4, intent(in)                     :: file_handle
   integer, intent(inout)                    :: status
   integer, intent(inout)                    :: pos ! req'd for stream read only
  ! Next 8 optional because initial and recursive calls require only subsets 
   type(FEWSNET_bil__header), intent(inout), optional :: headerInfo
   integer*4, intent(inout),                 optional :: numEntriesRead
   integer*4, intent(in),           optional :: reading_token_of_this_class
   character(gBIL_MAX_HEADERFILE_LINELENGTH), &
                     intent(inout), optional :: charStrVal
   integer*4, intent(inout),        optional :: int4Val
   real*4, intent(inout),           optional :: real4Val
   real*8, intent(inout),           optional :: real8Val
   integer*8, intent(inout),        optional :: int8Val
   logical, intent(in),             optional :: no_round

   character :: c
   integer*4 :: d
   integer*4, parameter :: nDelimeters = 4
   integer*4, parameter, dimension(nDelimeters) :: delimeters = (/ 32, 10, 13, 9 /)
  ! FORTRAN ASCII
  ! 32 SPACE - Token to "Split" in Tamuka Magadzire's VB GIS lib translated here
  ! 10 LF - These two necessary in fortran at end of line 's
  ! 13 CR - These two necessary in fortran at end of line 's
  ! 09 TAB - No guidance provided re ths. Incl for completeness. Cld be others too 
   logical   :: cIsDelim
   character(len=gBIL_MAX_HEADERFILE_LINELENGTH) :: stringToken
   integer*4 :: sig_digits_right_of_dec_pt = gDEFAULT_SIGDIGSRTOFDECPT
   real*8    :: tmpDbl
   logical   :: no_round_val
   logical   :: read_gbil_identifierlabel
! ________________________________________________________________________

   write(stringToken, '(a)') ''
!- Loop over all characters in header file to extract header information:
   do

#ifdef HEADER_FILE_READ_METHOD_1_
      read(file_handle, pos=pos, iostat=status) c
#else
#ifdef HEADER_FILE_READ_METHOD_2_
      call fgetc(file_handle, c, status)
#endif
#endif
      if(status /= gFGETC_SUCCESS) return  ! File not accessible
      pos = pos + 1

      cIsDelim = .false.
      do d = 1, nDelimeters
         if (ichar(c) == delimeters(d)) then
            cIsDelim = .true. ; exit
         endif
      end do
!      print*, 'c=', c, ';cIsDelim=', cIsDelim

   !- No delimiter present:
      if(cIsDelim .neqv. .true.) then
         stringToken = trim(stringToken) // c ! accumulate the token

   !- Delimiter reached; a complete "token" has now been assembled:
      else

         if( len_trim(stringToken) == 0 ) cycle

         read_gbil_identifierlabel = .true.
         if ( present(reading_token_of_this_class) ) then
            if (reading_token_of_this_class /= gBIL_identifierLabel) then
               read_gbil_identifierlabel = .false.
            endif
         endif
      !- Header information present to be assigned to headerInfo pointer array:
         !if( ( (.not.(present(reading_token_of_this_class))) .or.    &
         !      (reading_token_of_this_class == gBIL_identifierLabel) &
         !  ) ) then
         if ( read_gbil_identifierlabel .eqv. .true. ) then

          ! Here we handle the FEWSNET_bil header file-specific identifier 
          ! labels and route the values to the appropriate struct data members:
          ! E.g.,
          !  byteorder I
          !  layout BIL
          !  nbits 8
          !  xdim 0.1
          !  ydim 0.1
          !  ncols 751
          !  nrows 801
          !  nbands 1
          !  ulxmap -20
          !  ulymap 40
          !
          ! We also handle the alternative, similar FEWSNET float header  
          ! file-specific entries:
          ! E.g.,
          !  byteorder lsbfirst
          !  cellsize 0.1
          !  nodata_value -32768
          !  ncols 421
          !  nrows 441
          !  nbands 1
          !  xllcorner 10
          !  yllcorner -37

            select case (trim(StrLowCase(stringToken)))
             ! 'byteorder' also in the alternative FEWSNET 'float' file format
             !  but will need to be post-processed if read-in from that file type

               case ('byteorder')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     charStrVal=headerInfo%byteorder)

             ! 'cellsize' from the alternative FEWSNET 'float' file format
               case ('cellsize')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     real4Val=headerInfo%xdim,                        &
                     no_round=.true.)
                  headerInfo%ydim=headerInfo%xdim

               case ('layout')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     charStrVal=headerInfo%layout)

               case ('nbits')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     int4Val=headerInfo%nbits)

               case ('xdim')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     real4Val=headerInfo%xdim)

               case ('ydim')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     real4Val=headerInfo%ydim)

               case ('ncols')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     int4Val=headerInfo%ncols)

               case ('nrows')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     int4Val=headerInfo%nrows)

               case ('nbands')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     int4Val=headerInfo%nbands)

              ! 'nodata_value', 'xllcorner', and 'yllcorner' all from 
              ! the alternative FEWSNET 'float' file format
               case ('nodata_value')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     int4Val=headerInfo%missingvalue)

               case ('xllcorner')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     real8Val=headerInfo%ulxmap,                      &
                     no_round=.true.)

              ! 'yllcorner' from the alternative 'float' file format is used 
              ! together with 'nrows' and 'ydim' in post-process to arrive at 'ulymap'
               case ('yllcorner')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     reading_token_of_this_class=gBIL_idVal,          &
                     real8Val=headerInfo%yllcorner,                   &
                     no_round=.true.)

               case ('ulxmap')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     real8Val=headerInfo%ulxmap)

               case ('ulymap')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     real8Val=headerInfo%ulymap)

              ! Note next 3 in the VB were: ahdr.bandrowbytes = CDbl(temp_arr(1))
              ! ahdr.totalrowbytes and ahdr.bandgapbytes both = CDbl(temp_arr(1))
              ! The translation here therefore assumes a potentially large value 
              ! for each accordingly, hence their designation as integer*8 values.
               case ('bandrowbytes')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     int8Val=headerInfo%bandrowbytes)

             ! In any case, these next two don't even appear to be used ? (yet?)
               case ('totalrowbytes')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     int8Val=headerInfo%totalrowbytes)

               case ('bandgapbytes')
                  call FEWSNET_bil_parse_arbitrarily_delineated_ascii(&
                     file_handle=file_handle, pos=pos, status=status, &
                     numEntriesRead=headerInfo%numEntriesRead,        &
                     reading_token_of_this_class=gBIL_idVal,          &
                     int8Val=headerInfo%bandgapbytes)

               case default
                ! unsupported entry; skip and proceed with reading

            end select

            if(status /= gFGETC_SUCCESS) return
            write(stringToken, '(a)') ''

      !- Alternative Header information - file handling exceptions:
         else
            ! I.e.
            !    if( (present(reading_token_of_this_class )) .and. &
            !        (reading_token_of_this_class == gBIL_idVal) ) then
            !
            ! Here we handle the FEWSNET_bil(& alt. float) header file entry-specific 
            ! value reading, in some places conversion, and memory assignment.

            if(present(charStrVal)) then
               write(charStrVal, '(a)') trim(stringToken)
            endif

            if(present(int4Val)) then
               ! int4Val's all were read in and interpreted as follows:
               ! int4Val's all      ahdr.ncols = CInt(temp_arr(1))
               read(stringToken, *) int4Val
            endif

            if ( present(real4Val) ) then
               no_round_val = .false.
               if ( present(no_round) ) then
                  if ( no_round .eqv. .true.) then
                     no_round_val = .true.
                  endif
               endif

               if ( no_round_val .eqv. .true. ) then
                     read(stringToken, *) real4Val
               else
                  ! real4Val's all were read in and interpreted as follows:
                  ! ahdr.xdim = Math.Round(CDbl(temp_arr(1)), decimalplaces)
                  read(stringToken, *) tmpDbl
                  real4Val = real(round_real8(tmpDbl, &
                  sig_digits_right_of_dec_pt_arg=sig_digits_right_of_dec_pt),4)
               endif
            endif

            if ( present(real8Val) ) then
               no_round_val = .false.
               if ( present(no_round) ) then
                  if ( no_round .eqv. .true.) then
                     no_round_val = .true.
                  endif
               endif

               if ( no_round_val .eqv. .true. ) then
                  read(stringToken, *) real8Val
               else
                  ! real8Val's all were read in and interpreted as follows:
                  ! ahdr.ulxmap = Math.Round(CDbl(temp_arr(1)), decimalplaces)
                  read(stringToken, *) tmpDbl
                  real8Val = round_real8(tmpDbl, &
                     sig_digits_right_of_dec_pt_arg=sig_digits_right_of_dec_pt)
               endif
            endif

            if(present(int8Val)) then
               ! int8Val's all were read in and interpreted as follows:
               ! ahdr.bandrowbytes = CDbl(temp_arr(1))
               read(stringToken, *) int8Val
            endif

            if ( present(numEntriesRead) ) then
               numEntriesRead = numEntriesRead + 1
            endif
            return
         endif

      endif
   end do

end subroutine FEWSNET_bil_parse_arbitrarily_delineated_ascii


! Don't call this function if headerInfo%refCoordsIsNull is .true.
function adjust_bil_headerFile_variables_post_read(headerInfo)

   logical :: adjust_bil_headerFile_variables_post_read ! success is .true.

   type(FEWSNET_bil__header), intent(inout) :: headerInfo

   real*8   :: refMinLon
   real*8   :: refMaxLon
   real*8   :: refMinLat
   real*8   :: refMaxLat
   real*8   :: refMaxLon_LR
   real*8   :: refMinLat_LR
   real*8   :: hdrMinLon
   real*8   :: hdrMaxLat
   real*8   :: hdrMaxLon_LR
   real*8   :: hdrMinLat_LR
   integer  :: decprec_checklat, decprec_checklon
   integer*4, parameter :: dPrec = gDEFAULT_SIGDIGSRTOFDECPT
   integer*4, parameter :: nRnd = (0-1)
   logical, parameter   :: decrement = .true.
! _______________________________________________________________

   adjust_bil_headerFile_variables_post_read = .false.

   refMinLon    = headerInfo%refCoords%minLon
   refMaxLat    = headerInfo%refCoords%maxLat
   refMaxLon_LR = calcEndPt(startPt=refMinLon, &
                      nSteps=(headerInfo%refCoords%maxX - 1), &
                      stepSz=headerInfo%refCoords%pixLon, &
                      precision=dPrec)
   refMinLat_LR = calcEndPt(refMaxLat, (headerInfo%refCoords%maxY - 1), &
                      headerInfo%refCoords%pixLat, dPrec, decrement)

   hdrMinLon    = headerInfo%ulxmap
   hdrMaxLat    = headerInfo%ulymap
   hdrMaxLon_LR = calcEndPt(hdrMinLon, (headerInfo%ncols - 1), real(headerInfo%xdim, 8), dPrec)
   hdrMinLat_LR = calcEndPt(hdrMaxLat,           &
                      (headerInfo%nrows - 1),    &
                       real(headerInfo%ydim, 8), &
                       dPrec, decrement)
   
!   print *, "ref:",refMinLon, refMaxLat, refMaxLon_LR, refMinLat_LR
!   print *, "hdr:",hdrMinLon, hdrMaxLat, hdrMaxLon_LR, hdrMinLat_LR  ! KRA - HERE!!

   if(headerInfo%abortFileReadIfHeaderDataRegionDoesntSubsumeAnalysisRegion &
      .eqv. .true.) then

   !- The check below is included to deal with differences in grid array data types:
   !-  LIS gridDesc array elements == real*4
   !-  BIL header grid array elements == real*8:
      decprec_checklat=0; decprec_checklon=0
      if( abs(hdrMinLon-refMinLon) < 0.01 .and. abs(hdrMaxLat-refMaxLat) < 0.01 ) then
        decprec_checklat = 1
      endif
      if( abs(hdrMinLon-refMinLon) < 0.01 .and. abs(hdrMaxLat-refMaxLat) < 0.01 ) then
        decprec_checklon = 1
      endif

    ! Constraint: Successful continuation of reading is contingent on hdr's space 
    !  containing at least the data points of the reference coordinate space ...

      if (.not. ( &
         (hdrMinLon    <= refMinLon)    .and. &
         (hdrMaxLon_LR >= refMaxLon_LR) .and. &
         (hdrMinLat_LR <= refMinLat_LR) .and. &
         (hdrMaxLat    >= refMaxLat)          &
               ) ) then
!          return  ! Function value remains .false., i.e. no success
         if( decprec_checklat .ne. 1 .or. decprec_checklon .ne. 1 ) then
           write(LIS_logunit,*) "ERR MSG: LIS run domain is just outside a WRSI-parameter file domain. "
           if(hdrMinLon > refMinLon)  write(LIS_logunit,*) &
               " WRSI parm min lon:",hdrMinLon, "; LIS run min lon: ",refMinLon
           if( hdrMaxLon_LR < refMaxLon_LR )  write(LIS_logunit,*)  &
               " WRSI parm max lon:",hdrMaxLon_LR,"; LIS run max lon: ",refMaxLon_LR
           if( hdrMinLat_LR > refMinLat_LR ) write(LIS_logunit,*) &
               " WRSI parm min lat:",hdrMinLat_LR,"; LIS run min lat: ",refMinLat_LR
           if( hdrMaxLat < refMaxLat ) write(LIS_logunit,*) &
               " WRSI parm max lat:",hdrMaxLat, "; LIS run max lat: ",refMaxLat
           write(LIS_logunit,*) "      -- Please check your LIS domain or parameter file  "
           write(LIS_logunit,*) "     (sometimes shifting run domain slightly can help address this issue.)"
           write(LIS_logunit,*) " --  "
         endif
      endif
   endif

! gCOORD_SPECIAL_VALUE_1 ==> set as undefined, -9999:
   if (refMinLon == gCOORD_SPECIAL_VALUE_1) refMinLon = hdrMinLon
   refMaxLon = headerInfo%refCoords%maxLon
   if (refMaxLon == gCOORD_SPECIAL_VALUE_1) refMaxLon = &
      calcEndPt(hdrMinLon, (headerInfo%ncols - 1), real(headerInfo%xdim, 8), nRnd)
   refMinLat = headerInfo%refCoords%minLat
   if (refMinLat == gCOORD_SPECIAL_VALUE_1) refMinLat = &
      calcEndPt(hdrMaxLat, (headerInfo%nrows - 1), real(headerInfo%ydim, 8), nRnd, decrement)
   if (refMaxLat == gCOORD_SPECIAL_VALUE_1) refMaxLat = hdrMaxLat

 ! Calculate for the headerInfo, minX, maxX, minY, maxY:
   call calcXyBounds( &
                     minLon=refMinLon,              &
                     ulxmap=headerInfo%ulxmap,      &
                     xdim=real(headerInfo%xdim, 8), &
                     maxLon=refMaxLon,              &
                     ulymap=headerInfo%ulymap,      &
                     minLat=refMinLat,              &
                     ydim=real(headerInfo%ydim, 8), &
                     maxLat=refMaxLat,              &
                     minX=headerInfo%minX,          &
                     maxX=headerInfo%maxX,          &
                     minY=headerInfo%minY,          &
                     maxY=headerInfo%maxY           ) 

   headerInfo%minX_corr_spatl_overflow = restrictRangeInt4( & 
      headerInfo%minX, lowend=0, highend=headerInfo%ncols, outInt4vL=0, outInt4vH=0)
   headerInfo%maxX_corr_spatl_overflow = &
      restrictRangeInt4(headerInfo%maxX, 0, headerInfo%ncols, &
      (headerInfo%ncols - 1), (headerInfo%ncols - 1))
   headerInfo%minY_corr_spatl_overflow = & 
      restrictRangeInt4(headerInfo%minY, 0, headerInfo%nrows, 0, 0)
   headerInfo%maxY_corr_spatl_overflow = &
      restrictRangeInt4(headerInfo%maxY, 0, headerInfo%nrows, &
      (headerInfo%nrows - 1), (headerInfo%nrows - 1))

   headerInfo%ulxmap_adj = calcEndPt( &
      headerInfo%ulxmap, headerInfo%minX, real(headerInfo%xdim, 8), nRnd)
   headerInfo%ulymap_adj = calcEndPt( &
      headerInfo%ulymap, headerInfo%minY, real(headerInfo%ydim, 8), nRnd, decrement)
   headerInfo%ulxmap_adj = round_real8(headerInfo%ulxmap_adj, &
      sig_digits_total=gDEFAULT_SIGDIGSTOTAL)
   headerInfo%ulymap_adj = round_real8(headerInfo%ulymap_adj, &
      sig_digits_total=gDEFAULT_SIGDIGSTOTAL)
   headerInfo%ncols_adj  = headerInfo%maxX - headerInfo%minX + 1
   headerInfo%nrows_adj  = headerInfo%maxY - headerInfo%minY + 1

   headerInfo%ulxmap_corr_spatl_overflow = calcEndPt( &
              headerInfo%ulxmap, headerInfo%minX_corr_spatl_overflow, &
              real(headerInfo%xdim, 8), nRnd)
   headerInfo%ulymap_corr_spatl_overflow = calcEndPt( &
              headerInfo%ulymap, headerInfo%minY_corr_spatl_overflow, &
              real(headerInfo%ydim, 8), nRnd, decrement)
   headerInfo%ulxmap_corr_spatl_overflow = round_real8( &
              headerInfo%ulxmap_corr_spatl_overflow, sig_digits_total=gDEFAULT_SIGDIGSTOTAL)
   headerInfo%ulymap_corr_spatl_overflow = round_real8( &
              headerInfo%ulymap_corr_spatl_overflow, sig_digits_total=gDEFAULT_SIGDIGSTOTAL)
   headerInfo%ncols_corr_spatl_overflow  = headerInfo%maxX_corr_spatl_overflow - &
              headerInfo%minX_corr_spatl_overflow + 1
   headerInfo%nrows_corr_spatl_overflow  = headerInfo%maxY_corr_spatl_overflow - &
              headerInfo%minY_corr_spatl_overflow + 1

 ! Adjust array index variables for fortran:
   headerInfo%minX = headerInfo%minX +1
   headerInfo%maxX = headerInfo%maxX +1
   headerInfo%minY = headerInfo%minY +1
   headerInfo%maxY = headerInfo%maxY +1
   headerInfo%minX_corr_spatl_overflow = headerInfo%minX_corr_spatl_overflow +1
   headerInfo%maxX_corr_spatl_overflow = headerInfo%maxX_corr_spatl_overflow +1
   headerInfo%minY_corr_spatl_overflow = headerInfo%minY_corr_spatl_overflow +1
   headerInfo%maxY_corr_spatl_overflow = headerInfo%maxY_corr_spatl_overflow +1

   adjust_bil_headerFile_variables_post_read = .true.

end function adjust_bil_headerFile_variables_post_read


function calcEndPt(startPt, nSteps, stepSz, precision, decrement)

   real*8                :: calcEndPt

   real*8,    intent(in) :: startPt
   integer*4, intent(in) :: nSteps
   real*8,    intent(in) :: stepSz
   integer*4             :: precision
   logical, optional     :: decrement

   integer*4 :: countUp1DownNeg1

   countUp1DownNeg1 = 1
   if ( present(decrement) ) then
      if ( decrement .eqv. .true. ) then
         countUp1DownNeg1 = (0 - 1)
      endif
   endif

   calcEndPt = startPt + countUp1DownNeg1*(nSteps * stepSz)
   if ( precision > 0 ) then
      calcEndPt = round_real8(calcEndPt, &
                              sig_digits_right_of_dec_pt_arg=precision)
   endif

end function calcEndPt


subroutine calcXyBounds(addOffset_arg, ifMidPointRoundUp, &
      minLon, ulxmap, xdim, maxLon, ulymap, minLat, ydim, maxLat, &
      minX, maxX, minY, maxY)

   integer*4, intent(in), optional :: addOffset_arg
   logical, intent(in),   optional :: ifMidPointRoundUp
   real*8, intent(in) :: minLon
   real*8, intent(in) :: ulxmap
   real*8, intent(in) :: xdim
   real*8, intent(in) :: maxLon
   real*8, intent(in) :: ulymap
   real*8, intent(in) :: minLat
   real*8, intent(in) :: ydim
   real*8, intent(in) :: maxLat
   integer*4, intent(inout), optional :: minX
   integer*4, intent(inout), optional :: maxX
   integer*4, intent(inout), optional :: minY
   integer*4, intent(inout), optional :: maxY

   integer*4 :: addOffset
   logical   :: ifMPalwaysRoundUp
! ____________________________________________

   addOffset = 0
   if(present(addOffset_arg)) then
      addOffset = addOffset_arg
   endif
   ! Round up option should not be used if min's are not known in advance to 
   ! evaluate to non-midpoint values.  This option is intended for max calc.
   ifMPalwaysRoundUp = .false.
   if(present(ifMidPointRoundUp)) then
      ifMPalwaysRoundUp = ifMidPointRoundUp
   endif

   if(present(minX)) minX = CInt4((((minLon - ulxmap) / xdim) + addOffset), &
      ifMidPointRoundUp=ifMPalwaysRoundUp)
   if(present(maxX)) maxX = CInt4((((minLon - ulxmap) / xdim + (maxLon - minLon) / xdim) &
      + addOffset), ifMidPointRoundUp=ifMPalwaysRoundUp)
   if(present(maxY)) maxY = CInt4((((ulymap - minLat) / ydim) + addOffset), &
      ifMidPointRoundUp=ifMPalwaysRoundUp)
   if(present(minY)) minY = CInt4((((ulymap - minLat) / ydim - (maxLat - minLat) / ydim) &
      + addOffset), ifMidPointRoundUp=ifMPalwaysRoundUp)

end subroutine calcXyBounds


subroutine xy2ll(x, y, ulxmap, ulymap, xdim, ydim, longitude, latitude)

   integer*4, intent(in) :: x
   integer*4, intent(in) :: y
   real*8, intent(in)    :: ulxmap
   real*8, intent(in)    :: ulymap
   real*8, intent(in)    :: xdim
   real*8, intent(in)    :: ydim
   real*8, intent(inout) :: longitude
   real*8, intent(inout) :: latitude

   ! Adjust array index variables for fortran
   real*8, parameter :: arrIndexAdj = (0 +1)

   latitude = ulymap - (y - arrIndexAdj) * ydim
   longitude = ulxmap + (x - arrIndexAdj) * xdim

end subroutine xy2ll


subroutine ll2xy(longitude, latitude, ulxmap, ulymap, xdim, ydim, x, y)

   real*8,    intent(in)    :: longitude
   real*8,    intent(in)    :: latitude
   real*8,    intent(in)    :: ulxmap
   real*8,    intent(in)    :: ulymap
   real*8,    intent(in)    :: xdim
   real*8,    intent(in)    :: ydim
   integer*4, intent(inout) :: x
   integer*4, intent(inout) :: y

   y = CInt4((ulymap - latitude) / ydim) 
  ! In this case, rows start from 0 to (nrows-1). 
  ! For them to start from 1, add "+ 1" to this statement.
   x = CInt4((longitude - ulxmap) / xdim)
  ! Adjust array index variables for fortran
   y = y +1
   x = x +1

end subroutine ll2xy


! This function was called ReadSubWindowGeoCoords
! The window is cut from (minlon,minlat) to (maxlon,maxlat) of the analysis region
!  to which the output/destination ('Dst') array is presumed to be sized in advance.
!  Don't call this function if headerInfo%refCoordsIsNull is .true.
function memCpFEWSNET_bilDataSrcToDst( &
   real4ptr1dSrc, int2ptr2dSrc, &
   headerInfo, int2ptr2dDst, real4ptr2dDst, &
   check_if_can_copy_but_dont_copy, &
   correctSpatialOverflow_arg, &
   resample_arg, fillSpatialGaps_arg, &
   spatialGapFillValue_arg, &
   missingValueFrom_arg, &
   missingValueTo_arg, &
   missingValueOut_arg &
   )
 
   logical :: memCpFEWSNET_bilDataSrcToDst

   real*4, intent(in), dimension(:),     optional :: real4ptr1dSrc ! inparray(VB), int4ptr1dL(caller)
   integer*2, intent(in), dimension(:,:),optional :: int2ptr2dSrc  ! inparray(VB), int2ptr2dL(caller)
   type(FEWSNET_bil__header), intent(in)          :: headerInfo
   integer*2, pointer, dimension(:,:),   optional :: int2ptr2dDst  ! outarray(VB), int2ptr2d(caller), gPPT3d(ts, :, :)
   real*4, pointer, dimension(:,:),      optional :: real4ptr2dDst ! added initially for e1v parameters
   logical, intent(in),   optional :: check_if_can_copy_but_dont_copy
   logical, intent(in),   optional :: correctSpatialOverflow_arg
   logical, intent(in),   optional :: resample_arg
   logical, intent(in),   optional :: fillSpatialGaps_arg
   integer*2, intent(in), optional :: spatialGapFillValue_arg
   integer*4, intent(in), optional :: missingValueFrom_arg
   integer*4, intent(in), optional :: missingValueTo_arg
   integer*4, intent(in), optional :: missingValueOut_arg

   ! Control variables
   logical   :: correctSpatialOverflow
   logical   :: fillSpatialGaps
   logical   :: resample
   ! Source memory location variables
   integer*2 :: spatialGapFillValue
   integer*4 :: missingValueFrom
   integer*4 :: missingValueTo
   integer*4 :: missingValueOut
   real*8    :: ulymap_unadjusted
   real*4    :: ydim_unadjusted
   integer*4 :: ncols_unadjusted
   integer*4 :: nrows_unadjusted
   real*8    :: ulxmap_unadjusted
   real*4    :: xdim_unadjusted
   integer*4 :: minXsrc
   integer*4 :: maxXsrc
   integer*4 :: minYsrc
   integer*4 :: maxYsrc
   real*4    :: xdimSrc
   real*4    :: ydimSrc
   real*8    :: ulxmapSrc
   real*8    :: ulymapSrc
   integer*4 :: ncolsSrc
   integer*4 :: nrowsSrc
   ! Destination memory location variables
   integer*4 :: ncolsDst
   integer*4 :: nrowsDst
   real*8    :: pixLonDst
   real*8    :: pixLatDst
   real*8    :: minLonDst
   real*8    :: maxLatDst
   integer*4 :: minXdst
   integer*4 :: maxXdst
   integer*4 :: minYdst
   integer*4 :: maxYdst
   ! other function local variables
   integer*4 :: xDst
   integer*4 :: yDst
   integer   :: status
   integer*4 :: rsmplFactorHalved
   real*8    :: rsmplLat
   integer*4 :: oldY
   real*8    :: pixULy ! type 'Variant' (in the VB)
   real*8    :: pixLRy ! type 'Variant' (in the VB)
   real*8    :: rsmplLon
   integer*4 :: oldX
   real*8    :: pixULx ! type 'Variant' (in the VB)
   real*8    :: pixLRx ! type 'Variant' (in the VB)
   real*8    :: avgPixValue
   integer*4 :: numPoints
   integer*4 :: yRsmpl
   real*8    :: oldLat
   integer*4 :: xRsmpl
   real*8    :: oldLon
   real*4    :: real4v
   integer*4 :: xSrc
   integer*4 :: ySrc
   integer*4 :: xDstSz
   integer*4 :: yDstSz
   real*8    :: lon
   real*8    :: lat
! ______________________________

   memCpFEWSNET_bilDataSrcToDst = .false. ! default fail, until that is we succeed=.true.

   correctSpatialOverflow = .false. ! default
   if(present(correctSpatialOverflow_arg)) &
      correctSpatialOverflow = correctSpatialOverflow_arg
   resample = .false. ! default
   if(present(resample_arg)) &
      resample = resample_arg
   fillSpatialGaps = .false. ! default
   if(present(fillSpatialGaps_arg)) &
      fillSpatialGaps = fillSpatialGaps_arg
   spatialGapFillValue = gSPATIAL_GAP_FILL_VALUE_1 ! default
   if(present(spatialGapFillValue_arg)) &
      spatialGapFillValue = spatialGapFillValue_arg
   missingValueFrom = gDATA_SPECIAL_VALUE_1 ! default
   if(present(missingValueFrom_arg)) &
      missingValueFrom = missingValueFrom_arg
   missingValueTo = gDATA_SPECIAL_VALUE_1 ! default
   if(present(missingValueTo_arg)) &
      missingValueTo = missingValueTo_arg
   missingValueOut = gDATA_SPECIAL_VALUE_1 ! default
   if(present(missingValueOut_arg)) &
      missingValueOut = missingValueOut_arg

   ulymap_unadjusted = headerInfo%ulymap
   ydim_unadjusted   = headerInfo%ydim
   ncols_unadjusted  = headerInfo%ncols
   nrows_unadjusted  = headerInfo%nrows
   ulxmap_unadjusted = headerInfo%ulxmap
   xdim_unadjusted   = headerInfo%xdim

   minXsrc = headerInfo%minX
   maxXsrc = headerInfo%maxX
   minYsrc = headerInfo%minY
   maxYsrc = headerInfo%maxY
   xdimSrc = headerInfo%xdim
   ydimSrc = headerInfo%ydim
   ulxmapSrc = headerInfo%ulxmap
   ulymapSrc = headerInfo%ulymap
   ncolsSrc = headerInfo%ncols
   nrowsSrc = headerInfo%nrows
   !ulxmapSrc = headerInfo%ulxmap_adj
   !ulymapSrc = headerInfo%ulymap_adj
   !ncolsSrc = headerInfo%ncols_adj
   !nrowsSrc = headerInfo%nrows_adj
   if( (correctSpatialOverflow .eqv. .true.) .and. &
       (fillSpatialGaps .neqv. .true.) ) then
      minXsrc = headerInfo%minX_corr_spatl_overflow
      maxXsrc = headerInfo%maxX_corr_spatl_overflow
      minYsrc = headerInfo%minY_corr_spatl_overflow
      maxYsrc = headerInfo%maxY_corr_spatl_overflow
      !ulxmapSrc = headerInfo%ulxmap_corr_spatl_overflow
      !ulymapSrc = headerInfo%ulymap_corr_spatl_overflow
      !ncolsSrc = headerInfo%ncols_corr_spatl_overflow
      !nrowsSrc = headerInfo%ncols_corr_spatl_overflow
   endif
   pixLonDst = headerInfo%refCoords%pixLon
   pixLatDst = headerInfo%refCoords%pixLat
   minLonDst = headerInfo%refCoords%minLon
   maxLatDst = headerInfo%refCoords%maxLat

 ! Calculate ncolsDst, nrowsDst
   call calcXyBounds( &
                     addOffset_arg=1,                    &
                     minLon=minLonDst,                   &
                     ulxmap=minLonDst,                   &
                     xdim=pixLonDst,                     &
                     maxLon=headerInfo%refCoords%maxLon, &
                     ulymap=maxLatDst,                   &
                     minLat=headerInfo%refCoords%minLat, &
                     ydim=pixLatDst,                     &
                     maxLat=maxLatDst,                   &
                     maxX=ncolsDst,                      &
                     maxY=nrowsDst                       &
                     )
   !
   ! Regarding the following preliminary sanity checks here.
   ! minXsrc, maxXsrc, minYsrc, maxYsrc were calculated during header read
   ! based on the header file and analysis region bounding boxes both.
   ! Enforced is that the header/data region must sufficiently cover the 
   ! analysis region.
   !!!!!!!!!!!!!!!!!
   ! In this first check, ncolsSrc and nrowsSrc remain the ncols and nrows 
   ! values as-read direct from the header file.  This is how GeoWRSI's GIS 
   ! I/O lib did it.  This check basically says if we are not filling spatial 
   ! gaps, discontinue file read attempt if the analysis region -and- 
   ! possibly spatial overflow correction -adjusted extents are outside of 
   ! those extents dictated by the header file (1 thru ncols, and 1 thru nrows).
   if ( (fillSpatialGaps .neqv. .true.) .and. &
        ( .not. &
          ( &
            (minXsrc >= (0 +1)) .and.            &
            (minYsrc >= (0 +1)) .and.            &
            (maxXsrc <= (ncolsSrc - 1 +1)) .and. &
            (maxYsrc <= (nrowsSrc - 1 +1))       &
          ) &
        ) &
      ) &
      return ! Function return value remains .false., i.e. read failure
             ! New window coordinates are outside the limits of the 
             !  current image. Please specify different coords.
             ! If minXsrc <= 0 "Min longitude specified is smaller than image's"
             ! If minYsrc <= 0 "Min latitude specified is smaller than image's"
             ! If maxXsrc > ncolsSrc "Max longitude specified is larger than image's"
             ! If maxYsrc > nrowsSrc "Max latitude specified is larger than image's"
   !!!!!!!!!!!!!!!!!
   ! If we are filling spatial gaps, at least one half of the coordinate
   ! specification based on the analysis region and header in each of the dimensions 
   ! x and y (columns and rows) must be within the region specified by the data file's 
   ! header file in order to anticipate sufficient coverage to proceed with reading 
   ! the data file.  This is the way that the original I/O code was doing it.
   if ( (fillSpatialGaps .eqv. .true.) .and. &
        (.not. ( &
               ((minXsrc >= (0 +1)) .or. (maxXsrc <= ncolsSrc)) .and. &
               ((minYsrc >= (0 +1)) .or. (maxYsrc <= nrowsSrc)) &
               ) &
        ) &
      ) then
      return ! Function return value remains .false., i.e. read failure
             !"New window coordinates completely exclude the current image. 
             ! Please specify different coords or a different image"
   endif

   if(resample .eqv. .true.) then
    ! Calculate minXdst, maxXdst, minYdst, maxYdst
      call calcXyBounds( &
                        ifMidPointRoundUp=.true.,           &
                        minLon=headerInfo%refCoords%minLon, &
                        ulxmap=headerInfo%refCoords%minLon, &
                        xdim=headerInfo%refCoords%pixLon,   &
                        maxLon=headerInfo%refCoords%maxLon, &
                        ulymap=headerInfo%refCoords%maxLat, &
                        minLat=headerInfo%refCoords%minLat, &
                        ydim=headerInfo%refCoords%pixLat,   &
                        maxLat=headerInfo%refCoords%maxLat, &
                        minX=minXdst,                       &
                        maxX=maxXdst,                       &
                        minY=minYdst,                       &
                        maxY=maxYdst                        &
                        )
      if( (minXdst > maxXdst) .or. (minYdst > maxYdst) ) return
   endif

   if ( present(check_if_can_copy_but_dont_copy) ) then
      if ( check_if_can_copy_but_dont_copy .eqv. .true. ) then
         memCpFEWSNET_bilDataSrcToDst = .true.
         return
      endif
   endif


   if(fillSpatialGaps .eqv. .true.) then
     ! Initialize the outarray with the spatialgapfillvalue 
     ! so that any gaps that arise will automatically be filled.
      do xDst=(0 +1), ncolsDst, 1
         do yDst=(0 +1), nrowsDst, 1
            if(present(int2ptr2dDst) .eqv. .true.) int2ptr2dDst(xDst, yDst) = spatialGapFillValue
            if(present(real4ptr2dDst) .eqv. .true.) real4ptr2dDst(xDst, yDst) = spatialGapFillValue
         end do
      end do
   endif

   if(resample .eqv. .true.) then

     ! The following represents one means of going from finer data resolution 
     ! to a coarser model domain: a Simple Mean resampling
     ! 1d Src to (new)2d Dst
      if(present(int2ptr2dDst) .eqv. .true.) &
         allocate( int2ptr2dDst(((maxXdst - minXdst) + 1),((maxYdst - minYdst) + 1)), &
                   STAT=status )
      if(present(real4ptr2dDst) .eqv. .true.) &
         allocate( real4ptr2dDst(((maxXdst - minXdst) + 1),((maxYdst - minYdst) + 1)), &
                   STAT=status )
      if(status /= gALLOCATE_SUCCESS) return
      ! real4 Src to (new)int2 Dst
      missingValueOut = int(CInt2(real(missingValueOut, 8)), 4)
      rsmplFactorHalved = ceiling(CInt4(real((pixLonDst/xdimSrc), 8)) / real(2, 4))
      do yDst = minYdst, maxYdst, 1
         rsmplLat = maxLatDst - ( yDst * pixLatDst ) - (0.5 * pixLatDst)
         oldY     = floor( ((ulymap_unadjusted - rsmplLat) / &
                               real(ydim_unadjusted, 8)) ) 
         pixULy   = rsmplLat + (pixLonDst / 2)
         pixLRy   = rsmplLat - (pixLatDst / 2)
         do xDst = minXdst, maxXdst, 1
            rsmplLon = minLonDst + ( xDst * pixLonDst ) + (0.5 * pixLonDst)
            oldX     = floor( ((rsmplLon - ulxmap_unadjusted) / &
                                  real(xdim_unadjusted, 8)) ) 
            pixULx   = rsmplLon - (pixLonDst / 2)
            pixLRx   = rsmplLon + (pixLatDst / 2)
            avgPixValue = 0
            numPoints   = 0

            if(rsmplFactorHalved <= 1) then
               real4v = real4ptr1dSrc(oldY * ncols_unadjusted + oldX +1)
               if(present(int2ptr2dDst) .eqv. .true.) &
                  int2ptr2dDst(((xDst - minXdst) +1), ((yDst - minYdst) +1)) = CInt2(real(real4v, 8))
               if(present(real4ptr2dDst) .eqv. .true.) &
                  real4ptr2dDst(((xDst - minXdst) +1), ((yDst - minYdst) +1)) = real4v

            else ! begin if rsmplFactorHalved > 1
               do yRsmpl = (oldY - rsmplFactorHalved), (oldY + rsmplFactorHalved), 1
                  oldLat = &
                     round_real8( &
                        (ulymap_unadjusted - (yRsmpl - 0.5) * ydim_unadjusted), &
                        sig_digits_right_of_dec_pt_arg=gDEFAULT_SIGDIGSRTOFDECPT)
                  if( (oldLat <= pixULy) .and.  &
                      (oldLat >= pixLRy) .and.  &
                      (yRsmpl >= 0)        .and.  &
                      (yRsmpl < nrows_unadjusted ) &
                  ) then
                     do xRsmpl = (oldX - rsmplFactorHalved), (oldX + rsmplFactorHalved), 1
                        oldLon = &
                           round_real8( &
                              (ulxmap_unadjusted + (xRsmpl + 0.5) * xdim_unadjusted), &
                              sig_digits_right_of_dec_pt_arg=gDEFAULT_SIGDIGSRTOFDECPT)
                        if( (oldLon >= pixULx) .and.  &
                            (oldLon <= pixLRx) .and.  &
                            (xRsmpl >= 0)        .and.  &
                            (xRsmpl < ncols_unadjusted ) &
                        ) then
                           real4v = real4ptr1dSrc(yRsmpl * ncols_unadjusted + xRsmpl +1)
                           if(.not.( &
                                     (missingValueFrom <= real4v) .and. &
                                     (real4v <= missingValueTo) &
                                   ) &
                           ) then
                              avgPixValue = avgPixValue + real4v
                              numPoints = numPoints + 1
                           endif
                        endif
                     end do
                  endif
               end do 

               if(numPoints > 0) then
                  if(present(int2ptr2dDst) .eqv. .true.) &
                     int2ptr2dDst(((xDst - minXdst) +1), ((yDst - minYdst) +1)) = &
                        CInt2((avgPixValue / numPoints))
                  if(present(real4ptr2dDst) .eqv. .true.) &
                     real4ptr2dDst(((xDst - minXdst) +1), ((yDst - minYdst) +1)) = &
                        real((avgPixValue / numPoints), 4) ! rely on fortran compiler here for real8 -> real4
               else
                  if(present(int2ptr2dDst) .eqv. .true.) &
                     int2ptr2dDst(((xDst - minXdst) +1), ((yDst - minYdst) +1)) = &
                        int(missingValueOut, 2) ! missingValueOut was CInt2'd previously
                  if(present(real4ptr2dDst) .eqv. .true.) &
                     real4ptr2dDst(((xDst - minXdst) +1), ((yDst - minYdst) +1)) = &
                        missingValueOut
               endif
            endif ! end if rsmplFactorHalved > 1
         end do
      end do
      ! end Simple Mean resampling

#ifdef YREV_
      if(present(int2ptr2dDst) .eqv. .true.) &
         memCpFEWSNET_bilDataSrcToDst = &
          flipY(int2ptr2d=int2ptr2dDst, ncols=size(int2ptr2dDst, 1), nrows=size(int2ptr2dDst, 2))
      if(present(real4ptr2dDst) .eqv. .true.) &
         memCpFEWSNET_bilDataSrcToDst = &
          flipY(real4ptr2d=real4ptr2dDst, ncols=size(real4ptr2dDst, 1), nrows=size(real4ptr2dDst, 2))
#else
      memCpFEWSNET_bilDataSrcToDst = .true.
#endif
      return ! In the original I/O lib, this resample feature was part of a separate 
             ! function and so distinct from and mutually exclusive with everything below
   endif


   ! elseif(resample .neqv. .true.) then
   xDstSz = size(int2ptr2dDst, 1)
   yDstSz = size(int2ptr2dDst, 2)
   if ( (xdimSrc == pixLonDst) .and. &
        (ydimSrc == pixLatDst) ) then
      ! Pixel widths and heights between Src and Dst match
      do ySrc=minYsrc, maxYsrc, 1
         do xSrc=minXsrc, maxXsrc, 1
            ! Compute xDst and yDst from xSrc and ySrc
            xDst = xSrc - minXsrc +1
            yDst = ySrc - minYsrc +1
            if(fillSpatialGaps .neqv. .true.) then
               if( (xDst <= xDstSz) .and. (yDst <= yDstSz) ) then
                  int2ptr2dDst(xDst, yDst) = int2ptr2dSrc(xSrc, ySrc); cycle
               endif
            endif
            ! elseif(fillSpatialGaps .eqv. .true.) then
#ifdef BRAVE_
! Begin this has not been tested
            if ( &
                 ((xSrc - minXsrc) >= (0 +1)) .and. &
                 ((xSrc - minXsrc) <= ncolsDst) .and. &
                 (xSrc >= (0 +1)) .and. (xSrc <= ncolsSrc) &
                 .and. &
                 ((ySrc - minYsrc) >= (0 +1)) .and. &
                 ((ySrc - minYsrc) <= nrowsDst) .and. &
                 (ySrc >= (0 +1)) .and. (ySrc <= nrowsSrc) &
            ) then
               if( ((xSrc - minXsrc) <= xDstSz) .and. &
                   ((ySrc - minYsrc) <= yDstSz) ) &
                  int2ptr2dDst((xSrc - minXsrc), (ySrc - minYsrc)) = &
                     int2ptr2dSrc(xSrc, ySrc) ! This reference/assignment 
                                              ! is okay with array bounds
            endif
! End this has not been tested
#endif
         end do
      end do

   else
    ! Pixel widths and heights between Src and Dst don't match

      do xDst=1, ncolsDst, 1
         do yDst=1, nrowsDst, 1

            ! Compute xSrc and ySrc from xDst and yDst
            
            ! Subset by nearest neighbor
            call xy2ll(xDst, yDst, &
               real(minLonDst, 8), &
               real(maxLatDst, 8), &
               pixLonDst,          &
               pixLatDst,          &
               lon, lat)
            call ll2xy(lon, lat, &
               ulxmapSrc,          &
               ulymapSrc,          &
               real(xdimSrc, 8),   &
               real(ydimSrc, 8),   &
               xSrc, ySrc)

            ! Check that this reference/assignment okay with source array bounds
            if( &
                ( (xSrc >= (0 +1)) .and. (xSrc <= ncolsSrc) ) .and. &
                ( (ySrc >= (0 +1)) .and. (ySrc <= nrowsSrc) ) ) then
               if( (xDst <= xDstSz) .and. (yDst <= yDstSz) ) &
                  int2ptr2dDst(xDst, yDst) = int2ptr2dSrc(xSrc, ySrc)
            endif

         end do
      end do

      ! End Pixel widths and heights between Src and Dst don't match
   endif

#ifdef YREV_
   if(present(int2ptr2dDst) .eqv. .true.) &
      memCpFEWSNET_bilDataSrcToDst = &
       flipY(int2ptr2d=int2ptr2dDst, ncols=size(int2ptr2dDst, 1), nrows=size(int2ptr2dDst, 2))
   if(present(real4ptr2dDst) .eqv. .true.) &
      memCpFEWSNET_bilDataSrcToDst = &
       flipY(real4ptr2d=real4ptr2dDst, ncols=size(real4ptr2dDst, 1), nrows=size(real4ptr2dDst, 2))
#else
   memCpFEWSNET_bilDataSrcToDst = .true.
#endif

end function memCpFEWSNET_bilDataSrcToDst


! This subroutine was ReadAnyBilImageTo2D_ShortArray(fname, short theInputArr(,))
! assumes data arrays passed-in are all allocated to appropriate size.

function interpretFEWSNET_bilHdrInfAndPopulateMemLocation(dataFilename, &
     headerInfo, real4ptr1d, int2ptr2d   )

   logical :: interpretFEWSNET_bilHdrInfAndPopulateMemLocation

   character(*), intent(in)                           :: dataFilename
   type(FEWSNET_bil__header), intent(in)              :: headerInfo
   real*4, dimension(:), intent(inout),      optional :: real4ptr1d
   integer*2, dimension(:,:), intent(inout), optional :: int2ptr2d

   integer*4 :: bnMax
   integer*4 :: xMax
   integer*4 :: yMax
   integer*4 :: sizeOfDataFileInBytes
   integer   :: status
   logical   :: int2ptr1dLisAllocated
   integer*2, pointer, dimension(:) :: int2ptr1dL
   logical   :: int4ptr1dLisAllocated
   integer*4, pointer, dimension(:) :: int4ptr1dL
   logical   :: real4ptr1dLisAllocated
   real*4,    pointer, dimension(:) :: real4ptr1dL
   integer*4 :: arrNumValues
   logical   :: mem_check
   integer*4 :: nValsDst
   integer*4 :: i
! ________________________________________________________________

   interpretFEWSNET_bilHdrInfAndPopulateMemLocation = .false.

   int2ptr1dLisAllocated = .false.
   int4ptr1dLisAllocated = .false.
   real4ptr1dLisAllocated = .false.
   bnMax = 1
   if(headerInfo%nbands /= gInitValnbands) bnMax = headerInfo%nbands
   xMax = headerInfo%ncols
   yMax = headerInfo%nrows

 ! Alternatively, could do this for ea. inp type in case of mis-match dims & hdr info
 ! if(present(int2ptr2d)) then xMax=size(int2ptr2d,1); yMax=size(int2ptr2d,2); endif
   arrNumValues = bnMax * xMax * yMax
  ! bytesPerValue = (headerInfo%nbits / 8)
   sizeOfDataFileInBytes = arrNumValues * (headerInfo%nbits / 8)

   select case (headerInfo%nbits)
      case ( gBIL_UNSIGNEDBYTE )
         ! We read (unsigned) byte data into a int*2 array
         allocate(int2ptr1dL(arrNumValues), STAT=status)
         if(status /= gALLOCATE_SUCCESS) return ! couldn't allocate tmp array for read
         int2ptr1dLisAllocated = .true.
         if(readBinaryDataFile(dataFilename, &
            recl=int(sizeOfDataFileInBytes, 8), &
            nVals=arrNumValues, &
            byteorder=headerInfo%byteorder, &
            int2ptr1d=int2ptr1dL, readDataAsVb2010byteData=.true.) &
            .eqv. .false.) then
            deallocate(int2ptr1dL)
            return
         endif

      case ( gBIL_BIT )
#ifdef BRAVE_
! Begin this has not been tested
         ! We read bit data into a int*2 array
         allocate(int2ptr1dL(arrNumValues), STAT=status)
         if(status /= gALLOCATE_SUCCESS) return
         int2ptr1dLisAllocated = .true.
         if(readBinaryDataFile(dataFilename, recl=headerInfo%bandrowbytes, &
            byteorder=headerInfo%byteorder, &
            ncols=headerInfo%ncols, &
            nrows=headerInfo%nrows, &
            int2ptr1d=int2ptr1dL, readDataAsBits=.true.) &
            .eqv. .false.) then
            deallocate(int2ptr1dL)
            return
         endif
! End this has not been tested
#else
         return
#endif
 
      case ( gBIL_INT2 )
         allocate(int2ptr1dL(arrNumValues), STAT=status)
         if(status /= gALLOCATE_SUCCESS) return
         int2ptr1dLisAllocated = .true.
         if(readBinaryDataFile(dataFilename, &
            recl=int(sizeOfDataFileInBytes, 8), &
            nVals=arrNumValues, &
            byteorder=headerInfo%byteorder, &
            int2ptr1d=int2ptr1dL) &
            .eqv. .false.) then
            deallocate(int2ptr1dL)
            return
         endif

      case ( gBIL_INT4 )
         allocate(int4ptr1dL(arrNumValues), STAT=status)
         if(status /= gALLOCATE_SUCCESS) return
         int4ptr1dLisAllocated = .true.
         if(readBinaryDataFile(dataFilename, &
            recl=int(sizeOfDataFileInBytes, 8), &
            nVals=arrNumValues, &
            byteorder=headerInfo%byteorder, &
            int4ptr1d=int4ptr1dL) &
            .eqv. .false.) then
            deallocate(int4ptr1dL)
            return
         endif

      case ( gBIL_REAL4 )
         allocate(real4ptr1dL(arrNumValues), STAT=status)
         if(status /= gALLOCATE_SUCCESS) return
         real4ptr1dLisAllocated = .true.
         if(readBinaryDataFile(dataFilename, &
            recl=int(sizeOfDataFileInBytes, 8), &
            nVals=arrNumValues, &
            byteorder=headerInfo%byteorder, &
            real4ptr1d=real4ptr1dL) &
            .eqv. .false.) then
            deallocate(real4ptr1dL)
            return
         endif

      case default
            return ! Unhandled nbit value

   end select

  ! Assign read data to the provided memory container, 
  ! respecting its data type and shape: 

   if(present(real4ptr1d)) then ! populate target output array real4ptr1d..

      mem_check = .false.
      nValsDst = size(real4ptr1d, 1)
      if(int2ptr1dLisAllocated .eqv. .true.) then ! ..from int2ptr1dL
        call redim(mem_check=mem_check, nValsDstChk=nValsDst, int2ptr1dSrc=int2ptr1dL)
        if(mem_check .eqv. .true.) real4ptr1d = int2ptr1dL
      endif
      if(int4ptr1dLisAllocated .eqv. .true.) then ! ..from int4ptr1dL
        call redim(mem_check=mem_check, nValsDstChk=nValsDst, int4ptr1dSrc=int4ptr1dL)
        if(mem_check .eqv. .true.) real4ptr1d = int4ptr1dL
      endif
      if(real4ptr1dLisAllocated .eqv. .true.) then ! ..from real4ptr1dL
        call redim(mem_check=mem_check, nValsDstChk=nValsDst, real4ptr1dSrc=real4ptr1dL)
        if(mem_check .eqv. .true.) real4ptr1d = real4ptr1dL
      endif

   endif ! end populate target output array real4ptr1d

   if(present(int2ptr2d)) then ! populate target output array int2ptr2d..

      if(int2ptr1dLisAllocated .eqv. .true.) then ! ..from int2ptr1dL
         call redim( bnMaxDstChk=bnMax,   &
            xMaxDstChk=headerInfo%ncols, &
            yMaxDstChk=headerInfo%nrows, &
            int2ptr1dSrc=int2ptr1dL, int2ptr2dDst=int2ptr2d )
      endif
     ! If the data will be clipped owing to numerical precision, 
     ! at least do it in a controlled fashion
     ! Any value above the data type maximum value gets truncated to 
     ! the data type maximum value; the same for minimum datatype values
      if(int4ptr1dLisAllocated .eqv. .true.) then ! ..from int4ptr1dL
         do i = 1, arrNumValues
            int4ptr1dL(i) = int(CInt2(real(int4ptr1dL(i), 8)), 4)
         end do
         call redim( bnMaxDstChk=bnMax,  &
            xMaxDstChk=headerInfo%ncols, &
            yMaxDstChk=headerInfo%nrows, &
            int4ptr1dSrc=int4ptr1dL, int2ptr2dDst=int2ptr2d )
      endif
      if(real4ptr1dLisAllocated .eqv. .true.) then ! ..from real4ptr1dL
         do i=1, arrNumValues, 1
            real4ptr1dL(i) = real(CInt2(real(real4ptr1dL(i), 8)), 4)
         end do
         call redim( bnMaxDstChk=bnMax,  &
            xMaxDstChk=headerInfo%ncols, &
            yMaxDstChk=headerInfo%nrows, &
            real4ptr1dSrc=real4ptr1dL, int2ptr2dDst=int2ptr2d )
      endif

   endif ! end populate target output array int2ptr2d

   if(int2ptr1dLisAllocated  .eqv. .true.) deallocate(int2ptr1dL)
   if(int4ptr1dLisAllocated  .eqv. .true.) deallocate(int4ptr1dL)
   if(real4ptr1dLisAllocated .eqv. .true.) deallocate(real4ptr1dL)

   interpretFEWSNET_bilHdrInfAndPopulateMemLocation = .true.

end function interpretFEWSNET_bilHdrInfAndPopulateMemLocation


#ifdef YREV_
! Reverse data in the y dimension. Some .bil data goes from North to South 
! whereas LIS and GrADS go from S. to N.
! As of this writing (10/2011) this YREV_ option works for a single band
function flipY( real4ptr1d, int2ptr2d, real4ptr2d, ncols, nrows )

   logical :: flipY

   real*4, pointer, dimension(:), intent(inout),      optional :: real4ptr1d 
  ! for converBilToLISbinary and other such domain-agnostic etc. applications
   integer*2, pointer, dimension(:,:), intent(inout), optional :: int2ptr2d
   real*4, pointer, dimension(:,:), intent(inout),    optional :: real4ptr2d
   integer   :: ncols
   integer   :: nrows

   integer*4 :: arrNumValues
   integer   :: status
   integer   :: r
   integer   :: c
   integer   :: i
   integer   :: j
   real*4,    dimension(:),   allocatable :: real4ptr1d_yrevTmpL
   integer*2, dimension(:,:), allocatable :: int2ptr2d_yrevTmpL
   real*4,    dimension(:,:), allocatable :: real4ptr2d_yrevTmpL
!_____________________________________________________________________

   flipY = .false.

   if(present(real4ptr1d)) then ! reverse target output array real4ptr1d..
      arrNumValues = ncols*nrows
      allocate(real4ptr1d_yrevTmpL(arrNumValues), STAT=status)
      if(status == gALLOCATE_SUCCESS) then
         do r=1, nrows, 1
            do c=1, ncols, 1
               i = ((r-1)*ncols + c)
               j = arrNumValues - (r*ncols) + c
               real4ptr1d_yrevTmpL(i) = &
                real4ptr1d(j)
            end do
         end do
         real4ptr1d = real4ptr1d_yrevTmpL
         deallocate(real4ptr1d_yrevTmpL)
         flipY = .true.
      endif
   endif

   if(present(int2ptr2d)) then ! reverse target output array int2ptr2d..
      allocate(int2ptr2d_yrevTmpL(ncols, nrows), STAT=status)
      if(status == gALLOCATE_SUCCESS) then
         do r=1, nrows, 1
            do c=1, ncols, 1
               int2ptr2d_yrevTmpL(c,nrows-(r-1)) = int2ptr2d(c,r)
            end do
         end do
         int2ptr2d = int2ptr2d_yrevTmpL
         deallocate(int2ptr2d_yrevTmpL)
         flipY = .true.
      endif
   endif

   if(present(real4ptr2d)) then ! reverse target output array real4ptr2d..
      allocate(real4ptr2d_yrevTmpL(ncols, nrows), STAT=status)
      if(status == gALLOCATE_SUCCESS) then
         do r=1, nrows, 1
            do c=1, ncols, 1
               real4ptr2d_yrevTmpL(c,nrows-(r-1)) = real4ptr2d(c,r)
            end do
         end do
         real4ptr2d = real4ptr2d_yrevTmpL
         deallocate(real4ptr2d_yrevTmpL)
         flipY = .true.
      endif
   endif

end function flipY
#endif


subroutine redim( mem_check, bnMaxDstChk, xMaxDstChk, yMaxDstChk, &
   nValsDstChk, int2ptr1dSrc, int4ptr1dSrc, real4ptr1dSrc, int2ptr2dDst )

! Presence of mem_check instructs a check only and no actual redim, truth=alignment.
! Otherwise this subroutine redimensions (reshapes) the source data into the destination

   logical, intent(inout), optional :: mem_check ! Presence instructs check no redim
   integer*4, intent(in),  optional :: bnMaxDstChk
   integer*4, intent(in),  optional :: xMaxDstChk
   integer*4, intent(in),  optional :: yMaxDstChk
   integer*4, intent(in),  optional :: nValsDstChk
   integer*2, intent(in),  dimension(:), optional :: int2ptr1dSrc
   integer*4, intent(in),  dimension(:), optional :: int4ptr1dSrc
   real*4, intent(in),     dimension(:), optional :: real4ptr1dSrc
   integer*2, intent(inout), dimension(:,:), optional :: int2ptr2dDst

   integer*4 :: bnMaxDst
   integer*4 :: xMaxDst
   integer*4 :: yMaxDst
   integer*4 :: nValsSrc
   integer*4 :: nValsDst
! ____________________________________________________________

   ! Default mem_check, if present, to false until alignment is proven true 
   if(present(mem_check)) mem_check = .false.

   ! Set bnMaxDst, xMaxDst, yMaxDst .. 
   !!! .. based on bnMaxDst, xMaxDst, yMaxDst provided for checking
   if(present(bnMaxDstChk)) bnMaxDst = bnMaxDstChk
   if(present(xMaxDstChk))  xMaxDst  = xMaxDstChk
   if(present(yMaxDstChk))  yMaxDst  = yMaxDstChk
   !!! .. based on output array
   if(present(int2ptr2dDst)) then
      xMaxDst = size(int2ptr2dDst, 1)
      yMaxDst = size(int2ptr2dDst, 2)
   endif

   ! Check bnMaxDst against bnMaxDst provided for checking
   if(present(bnMaxDstChk)) then ! Check destination Band maximum
      if(bnMaxDst /= bnMaxDstChk) return ! mem_check failed
   endif
   ! Check xMaxDst against xMaxDst provided for checking
   if(present(xMaxDstChk)) then ! Check destination X maximum
      if(xMaxDst /= xMaxDstChk) return ! mem_check failed
   endif
   ! Check yMaxDst against yMaxDst provided for checking
   if(present(yMaxDstChk)) then ! Check destination Y maximum
      if(yMaxDst /= yMaxDstChk) return ! mem_check failed
   endif

   ! Set nValsSrc ..
   !!! .. default
   nValsSrc = 1
   !!! .. based on input array
   if(present(int2ptr1dSrc)) nValsSrc = size(int2ptr1dSrc, 1)
   if(present(int4ptr1dSrc)) nValsSrc = size(int4ptr1dSrc, 1)
   if(present(real4ptr1dSrc)) nValsSrc = size(real4ptr1dSrc, 1)

   ! Set nValsDst ..
   !!! .. default
   nValsDst = nValsSrc
   !!! .. based on bnMaxDstChk, xMaxDstChk, yMaxDstChk provided for checking
   if(present(bnMaxDstChk) .and. &
      present(xMaxDstChk)  .and. &
      present(yMaxDstChk)) then 
      nValsDst = bnMaxDst * xMaxDst * yMaxDst ! Down here the Dst's the same as DstChk's
   endif
   !!! .. based on nValsDst provided for checking
   if(present(nValsDstChk)) nValsDst = nValsDstChk

   ! Check nValsSrc against nValsDst
   if(nValsSrc /= nValsDst) return ! mem_check failed

   if(present(mem_check)) then ! Remember, presence of mem_check signifies check only
      ! Signal truth that memsize checked out (mem_check succeeded) and return
      mem_check = .true.; return
   endif

   if(present(int2ptr2dDst)) then
      if(present(int2ptr1dSrc)) then
         if(gReadAllFilesAsThoughWrittenRowMajor .eqv. .true.) then
            int2ptr2dDst =            reshape( int2ptr1dSrc, (/ xMaxDst, yMaxDst /) )
         else
            int2ptr2dDst = transpose( reshape( int2ptr1dSrc, (/ yMaxDst, xMaxDst /) ) )
         endif
      endif
      ! This subroutine only performs reshaping, no numeric finessing for these next two
      ! assignments, going from a larger data type to a smaller one.  Do this beforehand
      if(present(int4ptr1dSrc)) then
         if(gReadAllFilesAsThoughWrittenRowMajor .eqv. .true.) then
            int2ptr2dDst =            reshape( int4ptr1dSrc, (/ xMaxDst, yMaxDst /) )
         else
            int2ptr2dDst = transpose( reshape( int4ptr1dSrc, (/ yMaxDst, xMaxDst /) ) )
         endif
      endif
      if(present(real4ptr1dSrc)) then
         if(gReadAllFilesAsThoughWrittenRowMajor .eqv. .true.) then
            int2ptr2dDst =            reshape( real4ptr1dSrc, (/ xMaxDst, yMaxDst /) )
         else
            int2ptr2dDst = transpose( reshape( real4ptr1dSrc, (/ yMaxDst, xMaxDst /) ) )
         endif
      endif
   endif

end subroutine redim


#ifdef BRAVE_
! Begin this has not been tested
subroutine convertByteToBits(rawByte, leastToMostSigLtoR, bits)

   integer*1, intent(in) :: rawByte
   logical, intent(in), optional :: leastToMostSigLtoR
   integer*1, intent(inout) :: bits(8)

   integer*4 :: bit
   logical   :: least_to_most_val

   least_to_most_val = .false.
   if ( present(leastToMostSigLtoR) ) then
      if ( leastToMostSigLtoR .eqv. .true.) then
         least_to_most_val = .true.
      endif
   endif

   bits = 0

   do bit = 0, 7, 1
      if( .not.( least_to_most_val ) ) then
         ! Assemble a Big Endian memory bitmap; (default)
         ! i.e. most to least significant bits going from left to right
         if( btest(rawByte, bit) ) bits( (8-bit) ) = 1
      else
         ! Assemble a Little Endian bitmap in memory; 
         ! i.e. least to most significant bits going from left to right
         if( btest(rawByte, bit) ) bits( (bit+1) ) = 1
      endif
   end do

end subroutine convertByteToBits


function unpack_bits( dataFileHandle, bandrowbytes, byteorder, &
                      ncols, nrows, int2ptr1d )

!  This function basically just unpackages the bytes into their component bits.  
!  For consistency we perform this procedure in identical fashion as the equiv. 
!  procedure in the GIS I/O lib.  This function assumes that a file previously 
!  opened with the handle 'dataFileHandle' is open and ready for formatted 
!  reading with recl=bandrowbytes

   logical :: unpack_bits

   integer,   intent(in) :: dataFileHandle
   integer*8, intent(in) :: bandrowbytes
   character(*), intent(in) :: byteorder
   integer*4, intent(in) :: ncols
   integer*4, intent(in) :: nrows
   integer*2, intent(inout), dimension(:), optional :: int2ptr1d

   integer   :: status
   integer*1, allocatable, dimension(:)  :: rawBytes
   integer*1 :: bits(8)
   integer*4 :: x
   integer*4 :: y
   integer*4 :: by
   integer*4 :: bi

   unpack_bits = .false. ! Failure unless successful

  ! read as int*1
   allocate(rawBytes(bandrowbytes), STAT=status)
   if(status /= gALLOCATE_SUCCESS) return ! couldn't allocate tmp array for read

   do y = 1, nrows

      rawBytes = 0
      read(dataFileHandle, rec=y, iostat=status) rawBytes
      ! Check success of file read
      if(status /= gFILEREAD_SUCCESS) then
         deallocate(rawBytes)
         return ! File read failed
      endif
      ! Deal with endian
      call resolveBinterpretation(byteorder=byteorder, &
         nVals=bandrowbytes, int1ptr1d=rawBytes)

      x = 1
      do by = 1, bandrowbytes
         call convertByteToBits(rawByte=rawBytes(by), bits=bits)
         do bi=1, 8, 1
            if(x > ncols) exit
            if(present(int2ptr1d)) int2ptr1d( ((y-(0 +1))*ncols) + x ) = bits(bi)
            x = x + 1
         end do
         if(x > ncols) exit
      end do

   end do

   deallocate(rawBytes)
   unpack_bits = .true. ! Success

end function unpack_bits
! End this has not been tested
#endif

subroutine resolveBinterpretation(byteorder, nVals, &
        int1ptr1d, int2ptr1d, int4ptr1d, real4ptr1d )

! The following was Vb GIS I/O Lib subroutine "readBinaryDataFile", but its really 
! for processing the data to ensure its proper interpretation as-intended 
! provided the endian orientation (big or little) of the system, and of the file

   character(*),        intent(in) :: byteorder
   integer*8, intent(in), optional :: nVals
   integer*1, intent(inout), dimension(:), optional :: int1ptr1d
   integer*2, intent(inout), dimension(:), optional :: int2ptr1d
   integer*4, intent(inout), dimension(:), optional :: int4ptr1d
   real*4, intent(inout),    dimension(:), optional :: real4ptr1d

   logical   :: systemIsBigEndian
   logical   :: fileDataWasWrittenAsBigEndian
   integer*4 :: i
   integer*4 :: x
   integer*4 :: y
!___________________________________________________

   systemIsBigEndian = big_endian()
   ! From the header 'byteorder' values for the data file are: 
   ! 'M' or 'm' for Motorola, a.k.a. 'big endian' = most significant byte first/leftmost
   ! 'I' or 'i' for Intel, a.k.a. 'little endian' = most significant byte last/rightmost
   fileDataWasWrittenAsBigEndian = ( (trim(StrLowCase(byteorder)) == 'm') .or.        &
                                     (trim(StrLowCase(byteorder)) == 'motorola') .or. &
                                     (trim(StrLowCase(byteorder)) == 'msbfirst') )

   if(systemIsBigEndian .eqv. fileDataWasWrittenAsBigEndian) return
   ! Otherwise, the endian-ness of the file doesn't match that of the system
   ! Resolve

   if(present(int1ptr1d)) then
      do i=1, nVals, 1
         call reverse_byte_order(no_reverse_bits=.true., int1v=int1ptr1d(i))
      end do
   endif
   if(present(int2ptr1d)) then
      do i=1, nVals, 1
         call reverse_byte_order(no_reverse_bits=.true., int2v=int2ptr1d(i))
      end do
   endif
   if(present(int4ptr1d)) then
      do i=1, nVals, 1
         call reverse_byte_order(no_reverse_bits=.true., int4v=int4ptr1d(i))
      end do
   endif
   if(present(real4ptr1d)) then
      do i=1, nVals, 1
         call reverse_byte_order(no_reverse_bits=.true., real4v=real4ptr1d(i))
      end do
   endif

end subroutine resolveBinterpretation


function big_endian()

   logical :: big_endian

   integer*4 long(1)
   integer*2 short(2)

#ifndef ASSUME_BIG_ENDIAN_FEWSNET_BIL_READ_
#ifndef ASSUME_LITTLE_ENDIAN_FEWSNET_BIL_READ_

 ! Detect byteorder of system
   equivalence(long, short)
   long(1)  = 0
   short(1) = 0
   short(2) = 1
   if(long(1) .eq. 1) then
      big_endian = .true.
   else
      big_endian = .false.
   endif

#else
!#define ASSUME_LITTLE_ENDIAN_FEWSNET_BIL_READ_
   big_endian = .false.
#endif
#else
!#define ASSUME_BIG_ENDIAN_FEWSNET_BIL_READ_
   big_endian = .true.
#endif

end function big_endian


subroutine reverse_byte_order( no_reverse_bits, int1v, int2v, int4v, real4v )

  logical,   intent(in),     optional :: no_reverse_bits
  integer*1, intent(inout),  optional :: int1v
  integer*2, intent(inout),  optional :: int2v
  integer*4, intent(inout),  optional :: int4v
  real*4,    intent(inout),  optional :: real4v

 ! 'nRb' no reverse bits
  integer*4 :: nRbNeg1_else1
  integer*4 :: nRb7_else0
  integer*1 :: bits_as_int1
  integer*2 :: bits_as_int2
  integer*4 :: bits_as_int4
  integer*4 :: i
  integer*1 :: bits_as_int1_R
  integer*2 :: bits_as_int2_R
  integer*4 :: bits_as_int4_R

   bits_as_int1_R = 0
   bits_as_int2_R = 0
   bits_as_int4_R = 0

 ! Reverse order of bytes (either 1, 2, or 4 in number) in native integer space:
   nRbNeg1_else1 = 1
   nRb7_else0 = 0
   if ( present(no_reverse_bits) ) then
      if ( no_reverse_bits .eqv. .true. ) then
         nRbNeg1_else1 = (0-1)
         nRb7_else0 = 7
      endif
   endif

   if(present(int1v))  bits_as_int1 = transfer(int1v, int(1, 1))
   if(present(int2v))  bits_as_int2 = transfer(int2v, int(1, 2))
   if(present(int4v))  bits_as_int4 = transfer(int4v, int(1, 4))
   if(present(real4v)) bits_as_int4 = transfer(real4v, int(1, 4))

 ! in 8 bit integer space
   if(present(int1v)) then
      do i=0, 7, 1
         call mvbits(bits_as_int1, 0+i, 1, bits_as_int1_R, &
            (0 + (7 - (nRbNeg1_else1 * i) - nRb7_else0)) )
      end do
   endif

 ! in 16 bit integer space
   if(present(int2v)) then
      do i=0, 7, 1
         call mvbits(bits_as_int2, 8+i, 1, bits_as_int2_R, &
            (0 + (7 - (nRbNeg1_else1 * i) - nRb7_else0)) )
         call mvbits(bits_as_int2, 0+i, 1, bits_as_int2_R, &
            (8 + (7 - (nRbNeg1_else1 * i) - nRb7_else0)) )
      end do
   endif

 ! in 32 bit integer space
   if(present(int4v) .or. present(real4v)) then
      do i=0, 7, 1
         call mvbits(bits_as_int4, 24+i, 1, bits_as_int4_R, &
            (0 + (7 - (nRbNeg1_else1 * i) - nRb7_else0)) )
         call mvbits(bits_as_int4, 16+i, 1, bits_as_int4_R, &
            (8 + (7 - (nRbNeg1_else1 * i) - nRb7_else0)) )
         call mvbits(bits_as_int4,  8+i, 1, bits_as_int4_R, &
            (16 + (7 - (nRbNeg1_else1 * i) - nRb7_else0)) )
         call mvbits(bits_as_int4,  0+i, 1, bits_as_int4_R, &
            (24 + (7 - (nRbNeg1_else1 * i) - nRb7_else0)) )
      end do
   endif

 ! Transfer reversed order bytes to destination (native) bit-space:
   if(present(int1v)) int1v = transfer(bits_as_int1_R, int(1, 1))
   if(present(int2v)) int2v = transfer(bits_as_int2_R, int(1, 2))
   if(present(int4v)) int4v = transfer(bits_as_int4_R, int(1, 4))
   if(present(real4v)) real4v = transfer(bits_as_int4_R, real(1.0, 4))

end subroutine reverse_byte_order


function readBinaryDataFile(dataFilename, recl, nVals, byteorder, &
             ncols, nrows, int2ptr1d, readDataAsVb2010byteData, &
             readDataAsBits, int4ptr1d, real4ptr1d  )

! This function brings the actual bytes from file into memory.  
! Was called readDataFileTo2dArray in the GIS I/O library.  It retains the 
! full feature set of functionality and is more generic here.

   logical  :: readBinaryDataFile

   character(*), intent(in)                     :: dataFilename
   integer*8, intent(in)                        :: recl
  ! These next two required for resolving potential Endian mis-match 
   integer*4, intent(in),              optional :: nVals ! Only not required for bits
   character(*), intent(in)                     :: byteorder
  ! These next two optional, required for the unpack_bit algorithm   
   integer*4, intent(in),              optional :: ncols
   integer*4, intent(in),              optional :: nrows
   integer*2, pointer, dimension(:),   optional :: int2ptr1d
  ! A Visual Basic 2010 'Byte' is Unsigned Byte 0-255
   logical,                            optional :: readDataAsVb2010byteData 
   logical,                            optional :: readDataAsBits
   integer*4, pointer, dimension(:),   optional :: int4ptr1d
   real*4, pointer, dimension(:),      optional :: real4ptr1d

   integer :: fileHandle
   integer :: status
   integer*1, allocatable, dimension(:)         :: int1ptr1dL
   integer*4 :: i
   integer*4 :: xMax
   integer*4 :: yMax
   integer*4 :: x
   integer*4 :: y

   logical :: read_as_vb_byte, read_as_bits
! ____________________________________________________________ 

   readBinaryDataFile = .false. ! Failure unless successful

  ! Check if file exists
   if (checkFileExistence(dataFilename) /= 0) return ! File doesn't exist

  ! Check for and get available file handle
   fileHandle = getFreeFileHandle()
   if(fileHandle == (0-1)) return ! Couldn't obtain a free file handle

   ! Record length is treated as a single monolithic record
   ! for the whole entire image and let Fortran's read routine
   ! take care of the rest based on the destination memory 
   ! space provided, e.g. recl is just sizeOfDataFileInBytes
   ! where sizeOfDataFileInBytes is, say 360*181*4 for a 360x181 
   ! image containing 4-byte data points.  This is the typical 
   ! strategy, however alternatively we may read in a loop as 
   ! when we read as bits in which case recl is bandrowbytes.
   open (unit=fileHandle,         &
      file=trim(dataFilename),    &
      form='unformatted',         &
      access='direct',            &
      recl=recl,                  &
      iostat=status)
   ! Check success of file open
   if(status /= gFILEOPEN_SUCCESS) then ! File open failed
      return
   endif

   if(present(int2ptr1d)) then
    ! Data are read into this container as integer*16, unsigned byte, bits

      read_as_vb_byte = .false.
      if ( present(readDataAsVb2010byteData) )then
         if ( readDataAsVb2010byteData .eqv. .true.) then
            read_as_vb_byte = .true.
         endif
      endif

      read_as_bits = .false.
      if ( present(readDataAsBits) ) then
         if ( readDataAsBits .eqv. .true.) then
            read_as_bits = .true.
         endif
      endif


      if ( .not.( read_as_vb_byte .or. read_as_bits ) ) then
       ! Reading integer*16 'Short'
         read(fileHandle, rec=1, iostat=status) int2ptr1d
         if(status /= gFILEREAD_SUCCESS) then
            close(fileHandle)
            return
         endif
       ! Deal with endian
         call resolveBinterpretation(byteorder=byteorder, &
                 nVals=int(nVals, 8), int2ptr1d=int2ptr1d)
      endif

      if ( read_as_vb_byte .eqv. .true. ) then
       ! Reading 'Unsigned Byte'

       ! Read first as int*1
         allocate(int1ptr1dL(nVals), STAT=status)
         if(status /= gALLOCATE_SUCCESS) then
            close(fileHandle)
            return ! couldn't allocate tmp array for read
         endif
         read(fileHandle, rec=1, iostat=status) int1ptr1dL
       ! Check success of file read
         if(status /= gFILEREAD_SUCCESS) then
            deallocate(int1ptr1dL)
            close(fileHandle)
            return ! File read failed
         endif
         call resolveBinterpretation(byteorder=byteorder, &
                 nVals=int(nVals, 8), int1ptr1d=int1ptr1dL)
         ! Vb 2010 'Byte' data is unsigned; we assign the values, preserved 
         ! as intended, to int*2 always for simplicity
         do i=1, nVals, 1
          ! Assign 1 to 2 byte integer, & additionally rem sign:  if(n<0) n = 256+n.  
          !  I.e., turn 2's complement signed byte which is what fortran is capable of reading as, 
          !  into unsigned byte which is what VB means by byte.
            int2ptr1d(i) = int1ptr1dL(i)
            if(int2ptr1d(i)<0) int2ptr1d(i) = 256 + int2ptr1d(i)
         end do
         deallocate(int1ptr1dL)

      endif

      if( read_as_bits ) then
#ifdef BRAVE_
! Begin: This section has not been tested ...
         ! Reading bits
         ! Check success of unpack bits
         if(unpack_bits( &
               dataFileHandle=fileHandle, &
               bandrowbytes=recl, &
               byteorder=byteorder, & 
               ncols=ncols, &
               nrows=nrows, &
               int2ptr1d=int2ptr1d &
               ) &
            .neqv. .true.) then
            close(fileHandle)
            return
         endif
! End this has not been tested
#else
         close(fileHandle)
         return
#endif
      endif

   endif

   if(present(int4ptr1d)) then
      read(fileHandle, rec=1, iostat=status) int4ptr1d
      if(status /= gFILEREAD_SUCCESS) then
         close(fileHandle)
         return
      endif
      call resolveBinterpretation(byteorder=byteorder, &
              nVals=int(nVals, 8), int4ptr1d=int4ptr1d)
   endif

   if(present(real4ptr1d)) then
      read(fileHandle, rec=1, iostat=status) real4ptr1d
      if(status /= gFILEREAD_SUCCESS) then
         close(fileHandle)
         return
      endif
      call resolveBinterpretation(byteorder=byteorder, &
              nVals=int(nVals, 8), real4ptr1d=real4ptr1d)
   endif
 
   close(fileHandle)
   readBinaryDataFile = .true. ! Success

end function readBinaryDataFile


subroutine write_FEWSNET_bil( filename, headerInfo, &
                              real4ptr1d )

! Writes out the data and header file with the base filename filename
! in a format as specified by headerInfo given data of ptr argument

   type(charN), intent(in)                    :: filename
   type(FEWSNET_bil__header)                  :: headerInfo
   real*4, dimension(:), pointer,    optional :: real4ptr1d

   integer*4 :: bnMax
   integer*4 :: xMax
   integer*4 :: yMax
   integer*4 :: sizeOfDataFileInBytes
   integer   :: status
   integer*1, pointer, dimension(:) :: int1ptr1dL
   integer*2, pointer, dimension(:) :: int2ptr1dL
   integer*4, pointer, dimension(:) :: int4ptr1dL
   real*4, pointer, dimension(:)    :: real4ptr1dL
   integer*4 :: arrNumValues
   integer*4 :: i
#ifdef YREV_
   logical :: FLIPY_SUCCESS
#endif

   call write_FEWSNET_bil_headerFile(filename, headerInfo)

   bnMax = 1
   if(headerInfo%nbands /= gInitValnbands) bnMax = headerInfo%nbands
   xMax = headerInfo%ncols
   yMax = headerInfo%nrows
   arrNumValues = bnMax * xMax * yMax
   sizeOfDataFileInBytes = arrNumValues * (headerInfo%nbits / 8)

   if(present(real4ptr1d)) then
#ifdef YREV_
      FLIPY_SUCCESS=flipY(real4ptr1d=real4ptr1d, ncols=xMax, nrows=yMax)
#endif
      select case (headerInfo%nbits)
         case ( gBIL_BIT,gBIL_UNSIGNEDBYTE )
            ! We write bits, and (unsigned) byte data via int*1
#ifdef BRAVE_
! Begin this has not been tested
   ! for bits data
      ! Note, if bits were expanded-out 1-per integer, 
      ! there needs to be a re-packaging prior to write
      ! to be consistent with the read process which, 
      ! as-translated from TM's GIS lib at least, assumes 
      ! "compressed" bit packaging, i.e. bits all just 
      ! lined up across byte boundaries.
! End this has not been tested
#endif
            allocate(int1ptr1dL(arrNumValues), STAT=status)
            if(status /= gALLOCATE_SUCCESS) then
#ifdef YREV_
               ! Always restore vertical y ordering of the source array prior to return
               FLIPY_SUCCESS=flipY(real4ptr1d=real4ptr1d, ncols=xMax, nrows=yMax)
#endif
               return ! couldn't allocate tmp array for write
            endif
            do i=1, arrNumValues, 1
               int1ptr1dL(i) = CInt1(real(real4ptr1d(i), 8))
            end do
            call writeBinaryDataFile(filename%str,  &
               recl=int(sizeOfDataFileInBytes, 8 ), &
               nVals=arrNumValues,                  &
               byteorder=headerInfo%byteorder,      &
               int1ptr1d=int1ptr1dL)
            deallocate(int1ptr1dL)

         case ( gBIL_INT2 )
            allocate(int2ptr1dL(arrNumValues), STAT=status)
            if(status /= gALLOCATE_SUCCESS) then
#ifdef YREV_
               FLIPY_SUCCESS=flipY(real4ptr1d=real4ptr1d, ncols=xMax, nrows=yMax)
#endif
               return
            endif
            do i=1, arrNumValues, 1
               int2ptr1dL(i) = CInt2(real(real4ptr1d(i), 8))
            end do
            call writeBinaryDataFile(filename%str,  &
               recl=int(sizeOfDataFileInBytes, 8 ), &
               nVals=arrNumValues,                  &
               byteorder=headerInfo%byteorder,      &
               int2ptr1d=int2ptr1dL)
            deallocate(int2ptr1dL)

         case ( gBIL_INT4 )
            allocate(int4ptr1dL(arrNumValues), STAT=status)
            if(status /= gALLOCATE_SUCCESS) then
#ifdef YREV_
               FLIPY_SUCCESS=flipY(real4ptr1d=real4ptr1d, ncols=xMax, nrows=yMax)
#endif
               return
            endif
            do i=1, arrNumValues, 1
               int4ptr1dL(i) = CInt4(real(real4ptr1d(i), 8))
            end do
            call writeBinaryDataFile(filename%str,  &
               recl=int(sizeOfDataFileInBytes, 8 ), &
               nVals=arrNumValues,                  &
               byteorder=headerInfo%byteorder,      &
               int4ptr1d=int4ptr1dL)
            deallocate(int4ptr1dL)

         case ( gBIL_REAL4 )
            ! Temporary array done here in case we are byte swapping
            ! to create files for target architecture with opposite endian
            allocate(real4ptr1dL(arrNumValues), STAT=status)
            if(status /= gALLOCATE_SUCCESS) then
#ifdef YREV_
               FLIPY_SUCCESS=flipY(real4ptr1d=real4ptr1d, ncols=xMax, nrows=yMax)
#endif
               return
            endif
            do i=1, arrNumValues, 1
               real4ptr1dL(i) = real4ptr1d(i)
            end do
            call writeBinaryDataFile(filename%str,  &
               recl=int(sizeOfDataFileInBytes, 8 ), &
               nVals=arrNumValues,                  &
               byteorder=headerInfo%byteorder,      &
               real4ptr1d=real4ptr1dL)
            deallocate(real4ptr1dL)

         case default
#ifdef YREV_
            FLIPY_SUCCESS=flipY(real4ptr1d=real4ptr1d, ncols=xMax, nrows=yMax)
#endif
            return ! Unhandled nbit value

      end select

#ifdef YREV_
    ! Always restore vertical y ordering of the source array prior to return
      FLIPY_SUCCESS=flipY(real4ptr1d=real4ptr1d, ncols=xMax, nrows=yMax)
#endif

   end if ! end if(present(real4ptr1d))

end subroutine write_FEWSNET_bil


subroutine write_FEWSNET_bil_headerFile(filename, headerInfo)

   type(charN), intent(in)    :: filename
   type(FEWSNET_bil__header)  :: headerInfo

   integer :: fileHandle
   integer :: status
   character, parameter :: separator = ' '

   logical :: WriteIntegerValuesAsIntegers
! ____________________________________________

   WriteIntegerValuesAsIntegers = .false.

   fileHandle = getFreeFileHandle()
   if(fileHandle == (0-1) ) return !Couldn't obtain a free file handle
   open(unit=fileHandle,file=trim(filename%str)//'.hdr', iostat=status)
   if(status /= gFILEOPEN_SUCCESS) return ! File couldn't be opened
   write(fileHandle,'(a,a,a)') &
      "byteorder", separator, trim(headerInfo%byteorder)
   write(fileHandle,'(a,a,a)') &
      "layout", separator, trim(headerInfo%layout)
   write(fileHandle,'(a,a,a)') &
      "nbits", separator, writeinteger4(headerInfo%nbits)
   write(fileHandle,'(a,a,a)') &
      "xdim", separator, writereal4(headerInfo%xdim,WriteIntegerValuesAsIntegers)
   write(fileHandle,'(a,a,a)') &
      "ydim", separator, writereal4(headerInfo%ydim,WriteIntegerValuesAsIntegers)
   write(fileHandle,'(a,a,a)') &
      "ncols", separator, writeinteger4(headerInfo%ncols)
   write(fileHandle,'(a,a,a)') &
      "nrows", separator, writeinteger4(headerInfo%nrows)
   write(fileHandle,'(a,a,a)') &
      "nbands", separator, writeinteger4(headerInfo%nbands)
   write(fileHandle,'(a,a,a)') &
      "ulxmap", separator, writereal8(headerInfo%ulxmap,WriteIntegerValuesAsIntegers)
   write(fileHandle,'(a,a,a)') &
      "ulymap", separator, writereal8(headerInfo%ulymap,WriteIntegerValuesAsIntegers)
   close(fileHandle)

end subroutine write_FEWSNET_bil_headerFile


subroutine writeBinaryDataFile(filename, recl, nVals, byteorder, &
                int1ptr1d, int2ptr1d, int4ptr1d, real4ptr1d  )

   character(*), intent(in)                   :: filename
   integer*8, intent(in)                      :: recl
   integer*4, intent(in)                      :: nVals
   character(*), intent(in)                   :: byteorder
   integer*1, pointer, dimension(:), optional :: int1ptr1d
   integer*2, pointer, dimension(:), optional :: int2ptr1d
   integer*4, pointer, dimension(:), optional :: int4ptr1d
   real*4, pointer, dimension(:),    optional :: real4ptr1d

   integer :: fileHandle
   integer :: status
! __________________________________________________________

 ! Check for and get available file handle
   fileHandle = getFreeFileHandle()
   if(fileHandle == (0-1)) return ! Couldn't obtain a free file handle

   open (unit=fileHandle,          &
      file=trim(filename)//'.bil', &
      form='unformatted',          &
      access='direct',             &
      recl=recl,                   &
      iostat=status)
 ! Check success of file open
   if(status /= gFILEOPEN_SUCCESS) then ! File open failed
      return
   endif

   if(present(int1ptr1d)) then
    ! To write Big Endian, byteorder should = 'M'; Little Endian, 'I'
      call resolveBinterpretation(byteorder, int(nVals, 8), &
         int1ptr1d=int1ptr1d)
      write(fileHandle, rec=1, iostat=status) int1ptr1d
   endif

   if(present(int2ptr1d)) then
      call resolveBinterpretation(byteorder, int(nVals, 8), &
         int2ptr1d=int2ptr1d)
      write(fileHandle, rec=1, iostat=status) int2ptr1d
   endif

   if(present(int4ptr1d)) then
      call resolveBinterpretation(byteorder, int(nVals, 8), &
         int4ptr1d=int4ptr1d)
      write(fileHandle, rec=1, iostat=status) int4ptr1d
   endif

   if(present(real4ptr1d)) then
      call resolveBinterpretation(byteorder, int(nVals, 8), &
         real4ptr1d=real4ptr1d)
      write(fileHandle, rec=1, iostat=status) real4ptr1d
   endif

   close(fileHandle)

end subroutine writeBinaryDataFile


#ifdef DIAGNOSTICS_
!!!!!!!!!!!!!!!!  DIAGNOSTICS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! The following represent diagnostics for file I/O for the model code
 
!!! The following represent diagnostics for header file of the I/O for the model code

subroutine printcoords(coords)

   type(coord_spc) :: coords

   logical :: WriteIntegerValuesAsIntegers
! ____________________________________________

   WriteIntegerValuesAsIntegers = .false.

   print*, &
      "minLon=========|",writereal8(coords%minLon,WriteIntegerValuesAsIntegers)

   print*, &
      "maxLon=========|",writereal8(coords%maxLon,WriteIntegerValuesAsIntegers)

   print*, &
      "pixLon=========|",writereal8(coords%pixLon,WriteIntegerValuesAsIntegers)

   print*, &
      "minLat=========|",writereal8(coords%minLat,WriteIntegerValuesAsIntegers)

   print*, &
      "maxLat=========|",writereal8(coords%maxLat,WriteIntegerValuesAsIntegers)

   print*, &
      "pixLat=========|",writereal8(coords%pixLat,WriteIntegerValuesAsIntegers)

   print*, &
      "maxX===========|",writeinteger4(coords%maxX)

   print*, &
      "maxY===========|",writeinteger4(coords%maxY)

   print*, &
      "minX===========|",writeinteger4(coords%minX)

   print*, &
      "minY===========|",writeinteger4(coords%minY)

end subroutine printcoords


subroutine print_header_info(headerInfo)

   type(FEWSNET_bil__header) :: headerInfo

   logical :: WriteIntegerValuesAsIntegers
! ____________________________________________

   WriteIntegerValuesAsIntegers = .false.

   print*, "headerInfo:"

   print*, &
      "numEntriesRead=|",writeinteger4(headerInfo%numEntriesRead)

   print*, &
      "byteorder======|",trim(headerInfo%byteorder)

   print*, &
      "layout=========|",trim(headerInfo%layout)

   print*, &
      "nbits==========|",writeinteger4(headerInfo%nbits)

   print*, &
      "skipbytes======|",writeinteger4(headerInfo%skipbytes)

   print*, &
      "xdim===========|",writereal4(headerInfo%xdim,WriteIntegerValuesAsIntegers)

   print*, &
      "ydim===========|",writereal4(headerInfo%ydim,WriteIntegerValuesAsIntegers)

   print*, &
      "ncols==========|",writeinteger4(headerInfo%ncols)

   print*, &
      "nrows==========|",writeinteger4(headerInfo%nrows)

   print*, &
      "nbands=========|",writeinteger4(headerInfo%nbands)

   print*, &
      "ulxmap=========|",writereal8(headerInfo%ulxmap,WriteIntegerValuesAsIntegers)

   print*, &
      "ulymap=========|",writereal8(headerInfo%ulymap,WriteIntegerValuesAsIntegers)

   print*, &
      "missingvalue===|",writeinteger4(headerInfo%missingvalue)

   print*, &
      "bandrowbytes===|",writeinteger8(headerInfo%bandrowbytes)

   print*, &
      "totalrowbytes==|",writeinteger8(headerInfo%totalrowbytes)

   print*, &
      "bandgapbytes===|",writeinteger8(headerInfo%bandgapbytes)

   print*, &
      "yllcorner======|",writereal8(headerInfo%yllcorner,WriteIntegerValuesAsIntegers)

   print*, &
      "minimumNumRequi|",writeinteger4(headerInfo%minimumNumRequiredEntries)

   print*, "refCoordsIsNull|",writelogical(headerInfo%refCoordsIsNull)
   if(headerInfo%refCoordsIsNull .neqv. .true.) then
      print*, &
         "begin refCoords:"
      call printcoords(headerInfo%refCoords)
      print*, &
         "end refCoords"
   endif

   print*, &
      "abortFileReadIf|",writelogical(&
          headerInfo%abortFileReadIfHeaderDataRegionDoesntSubsumeAnalysisRegion)

   print*, &
      "minX===========|",writeinteger4(headerInfo%minX)

   print*, &
      "maxX===========|",writeinteger4(headerInfo%maxX)

   print*, &
      "minY===========|",writeinteger4(headerInfo%minY)

   print*, &
      "maxY===========|",writeinteger4(headerInfo%maxY)

   print*, &
      "minX_corr_spatl_overflow===========|",&
      writeinteger4(headerInfo%minX_corr_spatl_overflow)

   print*, &
      "maxX_corr_spatl_overflow===========|",&
      writeinteger4(headerInfo%maxX_corr_spatl_overflow)

   print*, &
      "minY_corr_spatl_overflow===========|",&
      writeinteger4(headerInfo%minY_corr_spatl_overflow)

   print*, &
      "maxY_corr_spatl_overflow===========|",&
      writeinteger4(headerInfo%maxY_corr_spatl_overflow)

   print*, &
      "ulxmap_adj=====|",writereal8(headerInfo%ulxmap_adj,WriteIntegerValuesAsIntegers)

   print*, &
      "ulymap_adj=====|",writereal8(headerInfo%ulymap_adj,WriteIntegerValuesAsIntegers)

   print*, &
      "ncols_adj======|",writeinteger4(headerInfo%ncols_adj)

   print*, &
      "nrows_adj======|",writeinteger4(headerInfo%nrows_adj)

   print*, &
      "ulxmap_corr_spa|",&
      writereal8(headerInfo%ulxmap_corr_spatl_overflow,WriteIntegerValuesAsIntegers)

   print*, &
      "ulymap_corr_spa|",&
      writereal8(headerInfo%ulymap_corr_spatl_overflow,WriteIntegerValuesAsIntegers)

   print*, &
      "ncols_corr_spat|",writeinteger4(headerInfo%ncols_corr_spatl_overflow)

   print*, &
      "nrows_corr_spat|",writeinteger4(headerInfo%nrows_corr_spatl_overflow)

end subroutine print_header_info


!!!!!!!!!!!!!!!!  END DIAGNOSTICS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif


!==============  Begin application-layer input file reading  ================

! This function returns .true. on successful read, .false. otherwise
function populate_arr_from_FEWSNET_bil_2d(mapFilename, &
   populate_arr_from_FEWSNET_bil_2d2, &
   int2ptr2d, real4ptr2d, &
   headerInfoOut )

   logical                                      :: populate_arr_from_FEWSNET_bil_2d

  ! Checking first for fname not equal to empty string here from LGP not WHC file read 
  ! and others, but can't hurt to have in all.  mapFilename%str must be allocated. +BW

   type(charN), intent(in)                      :: mapFilename
   integer*4, intent(in),              optional :: populate_arr_from_FEWSNET_bil_2d2
   integer*2, pointer, dimension(:,:), optional :: int2ptr2d
   real*4, pointer, dimension(:,:), optional    :: real4ptr2d ! added initially for e1v parameters
   type(FEWSNET_bil__header), pointer, optional :: headerInfoOut

   logical                                      :: reading_non_forecast_ppt
   type(FEWSNET_bil__header), target, save      :: hdrInf
   real*4, pointer, dimension(:)                :: real4ptr1dL
   integer                                      :: status
   integer*2, pointer, dimension(:,:)           :: int2ptr2dL
   real*4, pointer, dimension(:,:)              :: real4ptr2dL
   integer*4                                    :: xMaxDst
   integer*4                                    :: yMaxDst
   integer*4                                    :: x
   integer*4                                    :: y
! ____________________________________________________________

   populate_arr_from_FEWSNET_bil_2d = .false.
 
 ! First establish whether we are running a special case file read procedure:
   reading_non_forecast_ppt = .false.
   if ( present(populate_arr_from_FEWSNET_bil_2d2) ) then
      if (populate_arr_from_FEWSNET_bil_2d2 == gUSE_CURRENT_TIMESTEP_PPT) then
         reading_non_forecast_ppt = .true.
      endif
   endif

 ! Next attempt to read the header file:
   if(present(headerInfoOut)) headerInfoOut => hdrInf   ! Re-assigns pointer
   hdrInf%refCoords => gCoords     ! This must be done just prior to the read call
   hdrInf%refCoordsIsNull = .false.
   hdrInf%abortFileReadIfHeaderDataRegionDoesntSubsumeAnalysisRegion = &
      (.not.(reading_non_forecast_ppt .eqv. .true.)) ! Set this upfront, too

!   print *, "(before):hdrInf%ncols,gInitValNcols: ",hdrInf%ncols,gInitValNcols

   call readFEWSNET_bil_or_float_headerFile(mapFilename, hdrInf)

!   print *, "(after):hdrInf%ncols,gInitValNcols: ",hdrInf%ncols,gInitValNcols

   if (fileReadFailedBasedOnValXremainingForVarY(&
         int4VarY=hdrInf%ncols, int4ValX=gInitValNcols) .eqv. .true.) &
      return ! header fails to read

 ! Finally, perform the actual data file read into the provided memory location
   if (.not.(reading_non_forecast_ppt .eqv. .true.) ) then

      ! Stage data into a 1d memory per size & order requirements in its header
      allocate(real4ptr1dL(hdrInf%ncols * hdrInf%nrows), STAT=status)
      if(status /= gALLOCATE_SUCCESS) return
      if(interpretFEWSNET_bilHdrInfAndPopulateMemLocation( &
         trim(mapFilename%str)//'.bil', hdrInf, real4ptr1d=real4ptr1dL) &
         .eqv. .false.) then
         deallocate(real4ptr1dL)
         return
      endif

    ! Reformat the data into a 2d memory container:
      if(present(int2ptr2d) .eqv. .true.) then
            ! The reformatting routine is presumed to allocate pointer int2ptr2dL
         populate_arr_from_FEWSNET_bil_2d = memCpFEWSNET_bilDataSrcToDst( &
            real4ptr1dSrc=real4ptr1dL, &
            headerInfo=hdrInf, &
            int2ptr2dDst=int2ptr2dL, &
            correctSpatialOverflow_arg=.true., &
            resample_arg=.true. &
         )
         deallocate(real4ptr1dL)
   
         if(populate_arr_from_FEWSNET_bil_2d .eqv. .false.) return
         populate_arr_from_FEWSNET_bil_2d = .false.
   
       ! Install the data into the analysis region memory container for the model run:
       !
       ! This common case 2d file read needs a 1d float pointer to stage 
       ! the data, a local 2d short pointer to pass into the memcp for 
       ! reformatting, and finally this next loop to assign the tmp 2d 
       ! back to the ultimate destination container for use by the model. 
       ! The tmp 2d var and this loop are necessary just because we can't 
       ! change the size of the ultimate destination container yet the 
       ! data array has the potential to be of a different size than the 
       ! container depending on the data and user settings. Performing 
       ! 'resample' was the mode of choice for 2d files read by the 
       ! original GeoWRSI program.

         xMaxDst = size(int2ptr2d, 1)
         yMaxDst = size(int2ptr2d, 2)
         if( (size(int2ptr2dL, 1) < xMaxDst) .or. &
             (size(int2ptr2dL, 2) < yMaxDst) ) then
            deallocate(int2ptr2dL)
            return
         endif
         do x=1, xMaxDst, 1
            do y=1, yMaxDst, 1
               int2ptr2d(x, y) = int2ptr2dL(x, y)
            end do
         end do
   
         deallocate(int2ptr2dL)
      !ends if(present(int2ptr2d) .eqv. .true.)
      endif

      if(present(real4ptr2d) .eqv. .true.) then
       ! The reformatting routine is presumed to allocate pointer real4ptr2dL
         populate_arr_from_FEWSNET_bil_2d = memCpFEWSNET_bilDataSrcToDst( &
            real4ptr1dSrc=real4ptr1dL, &
            headerInfo=hdrInf, &
            real4ptr2dDst=real4ptr2dL, &
            correctSpatialOverflow_arg=.true., &
            resample_arg=.true. )

         deallocate(real4ptr1dL)
   
         if(populate_arr_from_FEWSNET_bil_2d .eqv. .false.) return
         populate_arr_from_FEWSNET_bil_2d = .false.
   
         xMaxDst = size(real4ptr2d, 1)
         yMaxDst = size(real4ptr2d, 2)
         if( (size(real4ptr2dL, 1) < xMaxDst) .or. &
             (size(real4ptr2dL, 2) < yMaxDst) ) then
            deallocate(real4ptr2dL)
            return
         endif
         do x=1, xMaxDst, 1
            do y=1, yMaxDst, 1
               real4ptr2d(x, y) = real4ptr2dL(x, y)
            end do
         end do
   
         deallocate(real4ptr2dL)
      !ends if(present(real4ptr2d) .eqv. .true.)
      endif

      populate_arr_from_FEWSNET_bil_2d = .true.

   else
      ! if (reading_non_forecast_ppt .eqv. .true.) then

    ! Stage data into a memory array matching its extents as specified in its header
      allocate(int2ptr2dL(hdrInf%ncols, hdrInf%nrows), STAT=status)
      if(status /= gALLOCATE_SUCCESS) return
      if(interpretFEWSNET_bilHdrInfAndPopulateMemLocation( &
         trim(mapFilename%str)//'.bil', hdrInf, int2ptr2d=int2ptr2dL) &
         .eqv. .false.) then
         deallocate(int2ptr2dL)
         return
      endif

      ! Install the data into the analysis region memory container for the model run
      populate_arr_from_FEWSNET_bil_2d = memCpFEWSNET_bilDataSrcToDst( &
         int2ptr2dSrc=int2ptr2dL, &
         headerInfo=hdrInf, &
         int2ptr2dDst=int2ptr2d &
      )
      ! If populate_arr_from_FEWSNET_bil_2d == .false., this means for GeoWRSI:  
      ! "PPT Image data for this year/dekad is smaller than analysis area specified 
      ! in the current region. Cannot load it up. Try changing your region's map 
      ! extent settings or using a different PPT file.  Will continue using average 
      ! rainfall grids, but will not calculate any SOS occurence for any dekads 
      ! after this dekad."

      deallocate(int2ptr2dL)

   endif

end function populate_arr_from_FEWSNET_bil_2d

!================end application-layer input file reading================


! ===  VB2Fort Routines  ===
! - Moved back to fbil module, like as with GeoWRSI.
!
! Simple, generic visual basic to fortran library 
! for replicating the functionality of core visual basic 
! routines in the fortran 90 programming language.  
! Translated by Brad Wind, UMD/ESSIC, 2011.
!

function round_real8(real8Var, sig_digits_total, sig_digits_right_of_dec_pt_arg) 

   real*8 :: round_real8

   real*8, intent(in) :: real8Var
  ! The following options are intended to be mutually exclusive
   integer*4, intent(in), optional :: sig_digits_total
   integer*4, intent(in), optional :: sig_digits_right_of_dec_pt_arg

   real*8    :: r8num
   integer*4 :: log10
   real*8    :: a_number_part
   real*8    :: mkSigDigsInt
   real*8    :: r8v_SigDigsInt
   integer*4 :: sig_digits_right_of_dec_pt
   real*8    :: r8num_less_pt_5

 ! Default: retaining type, round to the nearest whole number using banker's rounding
   r8num = real8Var
   sig_digits_right_of_dec_pt = (0-1)

   if (present(sig_digits_total)) then
      !Functionality for rounding to specific number of significant places
      !calculate log10 of the number (to the nearest whole number)
      if (real8Var == 0) then
         round_real8 = real8Var
         return
      else
         a_number_part = abs(real8Var)
         !if ((a_number_part >=1) .and. (a_number_part < 10)) then
         log10 = 1
         if (a_number_part >= 10) then
            do
               a_number_part = a_number_part / 10
               log10 = log10 + 1
               if(a_number_part < 10) exit
            end do
         elseif (a_number_part < 1) then
            log10 = 0
            do
               a_number_part = a_number_part * 10
               log10 = log10 - 1
               if(a_number_part > 1) exit
            end do
         endif
         r8num = real8Var / (10**log10)
         sig_digits_right_of_dec_pt = sig_digits_total
      endif
   endif

   if(present(sig_digits_right_of_dec_pt_arg)) &
      sig_digits_right_of_dec_pt = sig_digits_right_of_dec_pt_arg

   if(sig_digits_right_of_dec_pt /= (0-1)) then
      mkSigDigsInt = 10**sig_digits_right_of_dec_pt
      r8v_SigDigsInt = r8num * mkSigDigsInt

      r8num = r8v_SigDigsInt
   endif

   r8num_less_pt_5 = r8num - 0.5
   if ( (r8num_less_pt_5 == floor(r8num)) .and. &
        ((r8num_less_pt_5-(int(r8num_less_pt_5, 4)/2*2)) == 0) ) then
      ! "Banker's rounding" employed by MS VB 2010.  Nearest even
      ! number is selected when the fraction is exactly halfway 
      ! between to compensate for accumulating bias:
      round_real8 = nint(r8num) - 1

   else
      round_real8 = nint(r8num)
   endif

   if(sig_digits_right_of_dec_pt /= (0-1)) then
      ! un-make significant digits to be integer
      ! i.e. shift back to left of decimal point
      round_real8 = round_real8 / mkSigDigsInt
   endif

   if (present(sig_digits_total)) then
      round_real8 = round_real8 * (10**log10)
   endif

end function round_real8


!! These routines for controlled conversion between data types:

! In the VB this routine was called CByte
function CInt2(realVar)
   
   integer(kind=2) :: CInt2 ! Note: not actually a byte b/c VB uses an unsigned byte and this, 
                            ! though somewhat inefficient, was straightforward to preserve 
                            ! consistent (entirely non-negative) numerical interpretation

   real*8, intent(in) :: realVar

   real*8  :: CInt2R
   integer*2 :: largestInt2  = 1
   integer*2 :: smallestInt2 = 1


   ! Keep within the numeric model range in a controlled fashion
   largestInt2  = huge(largestInt2)
   smallestInt2 = 0 - (largestInt2+int(1, 4))
   CInt2R = &
      restrictRangeReal8(realVar, real(smallestInt2, 8), real(largestInt2, 8), &
         real(smallestInt2, 8), real(largestInt2, 8))

   CInt2R = CInt2R - 0.5
   if ( (CInt2R == floor(realVar)) .and. &
        ((CInt2R-(int(CInt2R, 4)/2*2)) == 0) ) then
      ! "Banker's rounding" employed by MS VB 2010.  Nearest even
      ! number is selected when the fraction is exactly halfway 
      ! between to compensate for accumulating bias:
      CInt2 = nint(realVar) - 1

   else
      CInt2 = nint(realVar)
   endif

end function CInt2


function CInt4(realVar, ifMidPointRoundUp)
   
   integer(kind=4) :: CInt4 ! This function identical to CInt2 except for this line and RoundUp 

   real*8, intent(in) :: realVar
   logical, intent(in), optional :: ifMidPointRoundUp

   real*8  :: CInt4R
   integer*4 :: largestInt4  = 1
   integer*4 :: smallestInt4 = 1
   logical :: ifMPalwaysRoundUp

 ! Keep within the numeric model range in a controlled fashion
   largestInt4  = huge(largestInt4)
   smallestInt4 = 0 - (largestInt4+int(1, 8))
   CInt4R = &
      restrictRangeReal8(realVar, real(smallestInt4, 8), real(largestInt4, 8), &
         real(smallestInt4, 8), real(largestInt4, 8))

   CInt4R = realVar - 0.5
   ifMPalwaysRoundUp = .false.
   if ( present(ifMidPointRoundUp) ) then
      if ( ifMidPointRoundUp .eqv. .true. ) then
         ifMPalwaysRoundUp = .true.
      endif
   endif
   if ( (CInt4R == floor(realVar)) .and. &
        ((ifMPalwaysRoundUp .eqv. .true.) .or. ((CInt4R-(int(CInt4R, 4)/2*2)) == 0)) ) then
      if(ifMPalwaysRoundUp .eqv. .true.) then
         CInt4 = floor(realVar) + 1
      else
         CInt4 = nint(realVar) - 1
      endif
   else
      CInt4 = nint(realVar)
   endif

end function CInt4


function restrictRangeInt4(int4v, lowend, highend, outInt4vL, outInt4vH)

   integer*4 :: restrictRangeInt4

   integer*4, intent(in) :: int4v
   integer*4, intent(in) :: lowend
   integer*4, intent(in) :: highend
   integer*4, intent(in) :: outInt4vL
   integer*4, intent(in) :: outInt4vH

   restrictRangeInt4 = int4v
   if(int4v < lowend ) restrictRangeInt4 = outInt4vL
   if(int4v > highend) restrictRangeInt4 = outInt4vH

end function restrictRangeInt4


function restrictRangeReal8(real8v, lowend, highend, outReal8vL, outReal8vH)

   real*8 :: restrictRangeReal8

   real*8, intent(in) :: real8v
   real*8, intent(in) :: lowend
   real*8, intent(in) :: highend
   real*8, intent(in) :: outReal8vL
   real*8, intent(in) :: outReal8vH

   restrictRangeReal8 = real8v
   if(real8v < lowend ) restrictRangeReal8 = outReal8vL
   if(real8v > highend) restrictRangeReal8 = outReal8vH

end function restrictRangeReal8


function CInt1(realVar)
   
   integer(kind=1) :: CInt1 

   real*8, intent(in) :: realVar

! This CInt1 function for purposes of writing to file on disk.
! Re-packages larger values into a byte size for fortran write.

   real*8  :: CInt1R

   ! Keep within the numeric model range in a controlled fashion
      ! For purposes of writing to file on disk, this next piece of lockdown
      ! puts us back into 2's complement signed-byte numerical 
      ! representation space, which was how fortran read it, from the 
      ! unsigned byte it was subsequently interpreted as in the fortran 
      ! in order to match the quantities intended by the original VB write:
   CInt1R = &
      restrictRangeReal8(realVar, real((0 - 128), 8), real(127, 8), &
                         real((0 - 128), 8), real((realVar - 256), 8)) 

   CInt1R = CInt1R - 0.5
   if ( (CInt1R == floor(realVar)) .and. &
        ((CInt1R-(int(CInt1R, 4)/2*2)) == 0) ) then
      ! "Banker's rounding" employed by MS VB 2010.  Nearest even
      ! number is selected when the fraction is exactly halfway 
      ! between to compensate for accumulating bias:
      CInt1 = nint(realVar) - 1
   else
      CInt1 = nint(realVar)
   endif

end function CInt1


#ifdef DIAGNOSTICS_
!!!!!!!!!!!!!!!!  DIAGNOSTICS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Specifically, this controlled formatting to standardize ascii writes 
! for ease of comparison of fortran output 
! with visual basic output.

function writereal4 ( var, WriteIntegerValuesAsIntegers )

   character(len=25) :: writereal4

   real*4, intent(in) :: var
   logical :: WriteIntegerValuesAsIntegers

   if( (WriteIntegerValuesAsIntegers .eqv. .true.) .and. &
       ((var-int(var,4)) == 0) ) then
         write(writereal4,*) int(var,4) 
   else
         write(writereal4,'(f18.7)') var
   endif

end function writereal4


function writereal8 ( var, WriteIntegerValuesAsIntegers )

   character(len=25) :: writereal8

   real*8, intent(in) :: var
   logical :: WriteIntegerValuesAsIntegers

   if( (WriteIntegerValuesAsIntegers .eqv. .true.) .and. &
       ((var-int(var,4)) == 0) ) then
         write(writereal8,*) int(var,4) 
   else
         write(writereal8,'(f21.10)') var
   endif

end function writereal8


function writelogical( var )

   character(len=25) :: writelogical

   logical, intent(in) :: var

 ! The following matches the visual basic:
   if(var .eqv. .false.) then
      write(writelogical,'(a)') "False"
   else
      write(writelogical,'(a)') "True"
   endif

end function writelogical


function writeinteger2 ( var )

   character(len=25) :: writeinteger2

   integer*2, intent(in) :: var

   write(writeinteger2, '(i10)') var

end function writeinteger2


function writeinteger4 ( var )

   character(len=25) :: writeinteger4

   integer*4, intent(in) :: var

   write(writeinteger4, '(i17)') var

end function writeinteger4


function writeinteger8 ( var )

   character(len=25) :: writeinteger8

   integer*8, intent(in) :: var

   write(writeinteger8, '(i17)') var

end function writeinteger8


!!!!!!!!!!!!!!!!  END DIAGNOSTICS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

end module fbil_module
