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
! !ROUTINE: read_petusgs
! \label{read_petusgs}
!
! !REVISION HISTORY:
!  17 Jul 2001: Jon Gottschalck; Initial code
!  06 Jan 2006: Yudong Tian; modified for LISv4.2
!  10 Mar 2012: K. Arsenault;  Added USGS PET Dataset
!  25 Oct 2013: K. Arsenault;  Added PET USGS to LIS7

!
! !INTERFACE:
subroutine read_petusgs (n, kk, findex, pet_filename, ferror_petusgs )

! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_logMod,  only : LIS_logunit, LIS_getNextUnitNumber, &
                          LIS_releaseUnitNumber
  use LIS_metforcingMod, only  : LIS_forc
  use petusgs_forcingMod, only : petusgs_struc
  use fbil_module

  implicit none

! !ARGUMENTS:   
  integer, intent(in) :: n 
  integer, intent(in) :: kk
  integer, intent(in) :: findex
  character(len=*), intent(in) :: pet_filename  
  integer,intent(out) :: ferror_petusgs
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  USGS PET data and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing dataset
!  \item[kk]
!    index of the forecast member
!  \item[pet\_filename]
!    name of the USGS PET file
!  \item[ferror\_petusgs]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_petusgs](\ref{interp_petusgs}) \newline
!    spatially interpolates the USGS PET data
!  \end{description}
!EOP

!==== Local Variables ==============================

  integer, parameter :: xd=360, yd=181       ! Dimension of original USGS PET 1.0 deg data
  integer, parameter :: npetusgs=xd*yd
  real               :: realpet(xd,yd)
  integer            :: i,j
  integer            :: ftn
  integer            :: index1, index2
  real, pointer      :: pet_regrid(:,:)      ! Interpolated PET array

! - BIL Format file read paramters:
  type(charN)       :: filename
  type(FEWSNET_bil__header), pointer :: hdrInfPtr
  type(FEWSNET_bil__header), target  :: hdrInfW
  real*4, pointer, dimension(:)      :: real4ptr1dL  ! 4-byte real array
  integer           :: status
  logical           :: file_exists

!====================  End Variable Definition  =======================

   allocate (pet_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n)))
   pet_regrid = -1.0
   realpet    = -1.0
   ferror_petusgs = 1  ! Default file available
 
!----------------------------------------------------------------------
! Read raw, uncompressed USGS PET data 
!----------------------------------------------------------------------

!- Inquire if the BIL *.hdr file exists:
   inquire( file=trim(pet_filename)//".hdr", exist=file_exists )
   if ( file_exists .eqv. .false. ) then
      write(LIS_logunit,*) "Missing USGS PET File, ", trim(pet_filename)
      write(LIS_logunit,*)  " Calling end run."
      call LIS_endrun()
      ferror_petusgs = 0
      return
   endif

!   ftn = LIS_getNextUnitNumber()

   allocate(filename%str)
   filename%str = pet_filename

   hdrInfPtr => hdrInfW
   hdrInfW%refCoordsIsNull = .true.
   nullify(hdrInfW%refCoords)     ! This must be done just prior to the read call
   hdrInfW%abortFileReadIfHeaderDataRegionDoesntSubsumeAnalysisRegion = .false.

!- Read BIL file header information:
   call readFEWSNET_bil_or_float_headerFile(filename, hdrInfW)

   if( fileReadFailedBasedOnValXremainingForVarY(&
         int4VarY=hdrInfPtr%ncols, int4ValX=gInitValNcols) .neqv. .true. ) then

!      write(LIS_logunit,*) "Reading in USGS PET data:: ", pet_filename

    ! Finally, perform the actual data file read into the provided memory location
    !  Stage data into a 1d memory per size & order requirements in its header:

      allocate(real4ptr1dL(hdrInfPtr%ncols * hdrInfPtr%nrows), STAT=status)

      if( status == gALLOCATE_SUCCESS ) then
         if( interpretFEWSNET_bilHdrInfAndPopulateMemLocation( &
             trim(filename%str)//'.bil', hdrInfPtr, real4ptr1d=real4ptr1dL ) &
            .eqv. .false.) then
            write(LIS_logunit,*) 'read file failed at populate mem location'
            ferror_petusgs = 0
         endif
      else
         write(LIS_logunit,*) 'read file failed at allocate mem location'
         ferror_petusgs = 0
      endif
   else
      print *, hdrInfPtr%ncols, gInitValNcols
      write(LIS_logunit,*) 'read file failed at header file read'
      ferror_petusgs = 0
   endif

    
!- Flip (y-reverse) for writing out:
!  (flipping data from 90N->90S to 90S->90N LIS standard)

   if( flipY(real4ptr1d=real4ptr1dL, &
       ncols=hdrInfPtr%ncols, nrows=hdrInfPtr%nrows) .eqv. .false.) then
       write(LIS_logunit,*) 'read file failed at y-reverse'
       ferror_petusgs = 0
   endif

!
!=== End of data reconfiguration

!- Data read in successfully:
   if( ferror_petusgs == 1 ) then

   !- Set any undefined points to -1 and everything else to realpet:
      index1 = 0
      do i = 1, yd         ! Same as hdrInfPtr%nrows
         do j = 1, xd      ! Same as hdrInfPtr%ncols
            index1 = index1 + 1
            realpet(j,i) = real4ptr1dL(index1)
         enddo
      enddo

  !-- Spatially interpolate original file domain to LIS grid domain:
      call interp_petusgs( n, findex, xd, yd, realpet, LIS_rc%gridDesc(n,:), &
                           LIS_rc%lnc(n), LIS_rc%lnr(n), pet_regrid )

      do j = 1,LIS_rc%lnr(n)
         do i = 1,LIS_rc%lnc(n)
            if( pet_regrid(i,j) .ne. -1.0 ) then
               index2 = LIS_domain(n)%gindex(i,j)
               if( index2 .ne. -1 ) then
                  petusgs_struc(n)%metdata2(kk,1,index2) = pet_regrid(i,j)   ! mm/day
               endif
            endif
         enddo
      enddo
    
      write(LIS_logunit,*) "Obtained USGS PET data:: ", trim(pet_filename)

   elseif( ferror_petusgs == 0 ) then
      write(LIS_logunit,*) "Missing USGS PET data ", trim(pet_filename)

   endif

!   call LIS_releaseUnitNumber(ftn)

   if( associated(real4ptr1dL) ) then
      deallocate(real4ptr1dL)
   endif
   deallocate(filename%str)
   deallocate(pet_regrid)

end subroutine read_petusgs
