!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: USAF_OBAMod
!
! REVISION HISTORY:
! 22 Jun 2017  Initial version.........................Eric Kemp/SSAI/NASA
! 07 Sep 2018  Changed EMK_ prefix to USAF_............Eric Kemp/SSAI/NASA
! 29 Nov 2023  Add QC flags............................Eric Kemp/SSAI/NASA
!
! DESCRIPTION:
! Contains data structure and methods for collecting observed, background,
! and analysis data values at observation locations, and writing to file
! for subsequent post-processing.
!-------------------------------------------------------------------------

module USAF_OBAMod

   ! Defaults
   implicit none
   private

   ! Public data structure to store observation, background, and analysis
   ! values for output
   type OBA
      private
      integer :: nobs
      character*10, allocatable :: networks(:)
      character*10, allocatable :: platforms(:)
      real, allocatable :: latitudes(:) ! Latitude of observation (deg N)
      real, allocatable :: longitudes(:) ! Longitude of observation (deg E)
      real, allocatable :: O(:) ! Observation values
      real, allocatable :: B(:) ! Background values
      real, allocatable :: A(:) ! Analysis values
      integer, allocatable :: qc(:) ! QC code
   end type OBA
   public :: OBA

   ! Interface for constructor
   interface createOBA
      module procedure newOBA
   end interface
   public :: createOBA

   ! Public methods
   public :: newOBA
   public :: destroyOBA
   public :: assignOBA
   public :: writeToFile
   public :: makeFilename

   ! Public parameters
   integer, parameter, public :: QC_UNKNOWN = 0
   integer, parameter, public :: QC_GOOD = 1
   integer, parameter, public :: QC_REJECT = 2
   integer, parameter, public :: QC_SUSPECT_BACKQC = 3
   integer, parameter, public :: QC_SUSPECT_SUPERSTATQC = 4

   ! Private parameter
   character(11), parameter :: qc_string(5) = (/ &
        'UNKNOWN    ', &
        'GOOD       ', &
        'REJECT     ', &
        'BACKQC     ', &
        'SUPERSTATQC'/)

contains

   !---------------------------------------------------------------------------
   ! Constructor
   function newOBA(nest,maxobs) result(this)
      
      ! Imports
      use AGRMET_forcingMod, only : agrmet_struc

      ! Defaults
      implicit none
      
      ! Arguments
      integer,intent(in) :: nest
      integer,intent(in), optional :: maxobs

      ! Result
      type(OBA) :: this

      ! Local variables
      integer :: maxnobs

      if (present(maxobs)) then
         maxnobs = maxobs
      else
         maxnobs = agrmet_struc(nest)%max_pcpobs
      end if

      ! Initialize variables
      this%nobs = 0
      allocate(this%networks(maxnobs))
      allocate(this%platforms(maxnobs))
      allocate(this%latitudes(maxnobs))
      allocate(this%longitudes(maxnobs))
      allocate(this%O(maxnobs))
      allocate(this%B(maxnobs))
      allocate(this%A(maxnobs))
      allocate(this%qc(maxnobs))

      this%networks = "NULL"
      this%platforms = "NULL"
      this%latitudes = 0
      this%longitudes = 0
      this%O = 0
      this%B = 0
      this%A = 0
      this%qc = QC_UNKNOWN

   end function newOBA

   !---------------------------------------------------------------------------
   ! Destructor
   subroutine destroyOBA(this)

      ! Defaults
      implicit none

      ! Arguments
      type(OBA), intent(inout) :: this

      this%nobs = 0
      deallocate(this%networks)
      deallocate(this%platforms)
      deallocate(this%latitudes)
      deallocate(this%longitudes)
      deallocate(this%O)
      deallocate(this%B)
      deallocate(this%A)
      deallocate(this%qc)

   end subroutine destroyOBA

   !---------------------------------------------------------------------------
   ! Add new diagnostics from one observation to the data structure.
   subroutine assignOBA(this,network,platform,latitude,longitude,O,B,A, &
        qc, set_qc_good)

      ! Imports
      use LIS_logmod, only : LIS_logunit

      ! Defaults
      implicit none

      ! Arguments
      type(OBA), intent(inout) :: this
      character(len=10), intent(in) :: network
      character(len=10), intent(in) :: platform
      real, intent(in) :: latitude
      real, intent(in) :: longitude
      real, intent(in) :: O
      real, intent(in) :: B
      real, intent(in) :: A
      integer, intent(in) :: qc
      logical, optional, intent(in) :: set_qc_good

      ! Local variables
      integer :: nobs

      ! Sanity check.  Since this is intended for an operational system,
      ! just print a warning and return if we see an array bounds problem.
      nobs = this%nobs
      if (nobs .eq. size(this%networks,1)) then
         write(LIS_logunit,*) &
              '[WARN], not enough memory for assigning OBA data!'
         return
      end if

      ! Assign the value
      nobs = nobs + 1
      this%nobs = nobs
      this%networks(nobs) = network
      this%platforms(nobs) = platform
      this%latitudes(nobs) = latitude
      this%longitudes(nobs) = longitude
      this%O(nobs) = O
      this%B(nobs) = B
      this%A(nobs) = A
      if (present(set_qc_good)) then
         if (set_qc_good .and. qc == QC_UNKNOWN) then
            this%qc(nobs) = QC_GOOD
         else
            this%qc(nobs) = qc
         end if
      else
         this%qc(nobs) = qc
      end if
   end subroutine assignOBA

   !---------------------------------------------------------------------------
   ! Write diagnostic output to ASCII file.
   subroutine writeToFile(this,filename)

      ! Imports
      use LIS_logMod, only: LIS_getNextUnitNumber, LIS_releaseUnitNumber, &
           LIS_logunit

      ! Defaults
      implicit none

      ! Arguments
      type(OBA), intent(in) :: this
      character(len=*), intent(in) :: filename

      ! Local variables
      integer :: iunit
      integer :: istat
      integer :: j

      write(LIS_logunit,*)'[INFO] Writing file ', trim(filename)

      ! Get unit for file
      iunit = LIS_getNextUnitNumber()

      ! Open output file
      open( unit = iunit, &
           file=trim(filename), &
           iostat = istat)

      ! Write OBA information to file
      write(iunit, *, iostat=istat) &
           '# Network Platform latitude longitude O   B   A   QC'
      do j = 1, this%nobs
         if (this%qc(j) == QC_REJECT) cycle
         if (trim(this%networks(j)) == "SUPEROB") cycle
         if (trim(this%networks(j)) == "SUPERGAGE") cycle
         write(iunit, 1000, iostat=istat) trim(this%networks(j)), &
              trim(this%platforms(j)), this%latitudes(j), &
              this%longitudes(j),&
              this%O(j), this%B(j), this%A(j), qc_string(this%qc(j)+1)
1000     format(a10,1x,a10,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,a11)
      end do ! j

      ! Close file
      call LIS_releaseUnitNumber(iunit)

   end subroutine writeToFile

   !---------------------------------------------------------------------------
   ! Creates filename for OBA statistics
   subroutine makeFilename(pathOBA,yyyymmddhh,accum,filename)

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: pathOBA
      character(len=10), intent(in) :: yyyymmddhh
      integer, intent(in) :: accum
      character(len=120), intent(out) :: filename

      write(filename,1000) trim(pathOBA),'/oba_',yyyymmddhh,'_',accum,'.txt'
      1000 format(A,A,A,A,I2.2,A)

   end subroutine makeFilename
end module USAF_OBAMod
