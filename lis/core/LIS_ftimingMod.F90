!-------------------------------------------------------------------------------
!BOP
!
! !MODULE: LIS_ftimingMod
!
      module LIS_ftimingMod
!
! !USES:
      use LIS_logMod
      use LIS_coreMod, only : LIS_rc, LIS_masterproc, LIS_localPet, LIS_npes
      use LIS_mpiMod, only : LIS_mpi_comm
#if ( defined SPMD )
      use mpi
#endif

      implicit none

      private
      public   :: Ftiming_Init
      public   :: Ftiming_On
      public   :: Ftiming_Off
      public   :: Ftiming_Reset
      public   :: Ftiming_Output

      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: qp = selected_real_kind(33, 4931)
      
      integer, parameter, private :: MAX_BLK_NAME = 24
      integer, parameter, private :: NBLKS = 120
      character (len=MAX_BLK_NAME), private :: list_blocknames(NBLKS)
      integer, private      :: tblock

      type tms
           private
           real(dp) :: usr, sys
      end type tms

      type (tms), private   :: accum(NBLKS)
      type (tms), private   :: last (NBLKS)

      real(dp), private       :: us_tmp1(NBLKS, 2)
      real(dp), private       :: us_tmp2(NBLKS, 2)

      integer               :: ierwtime

!#     include "mpif.h"
      real(dp), external      :: Mpi_Wtime
!
! !DESCRIPTION:
!
!EOP
!-------------------------------------------------------------------------------
      contains
#if ( defined SPMD )
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ftiming_Init
!
! !INTERFACE:
!
      subroutine Ftiming_Init ( )
!
! !DESCRIPTION:
! Initialize the timing tool.
!
! !LOCAL VARIABLES:
      integer :: nn
      real(dp)  :: wclk
!EOP
!-------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) "[INFO] Initialize LIS_ftiming"

      tblock = 0

      do nn = 1, NBLKS
         accum(nn)%usr = 0.0d0
         accum(nn)%sys = 0.0d0

         last (nn)%usr = 0.0d0
         last (nn)%sys = 0.0d0
      end do

      !------------------------------------------
      ! To reduce the overhead for the first call.
      !------------------------------------------

      wclk = Mpi_Wtime (ierwtime)

      return

      end subroutine Ftiming_Init
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ftiming_On
!
! !INTERFACE:
!
      subroutine Ftiming_On (block_name)
!
! !INPUT PARAMETERS:
      character (len=*)  :: block_name
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      character (len=MAX_BLK_NAME) :: ctmp
      integer                      :: iblock
      integer                      :: ii
!EOP
!-------------------------------------------------------------------------------
!EOC
      ctmp = Trim (block_name)

      iblock = 0

      do ii = 1, tblock
         if (TRIM(ctmp) == TRIM(list_blocknames(ii))) then
            iblock = ii
            EXIT
         end if
      end do

      if (iblock == 0) then
         tblock = tblock + 1
         iblock = tblock

         list_blocknames(iblock) = Trim (block_name)
      end if

      last(iblock)%usr = Mpi_Wtime (ierwtime)
      last(iblock)%sys = 0.0d0

      return

      end subroutine Ftiming_On
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ftiming_Off
!
! !INTERFACE:
!
      subroutine Ftiming_Off (block_name)

      character (len=*)            :: block_name

! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      character (len=MAX_BLK_NAME) :: ctmp
      integer                      :: iblock
      integer                      :: ii
      real(dp)                       :: wclk
!EOP
!-------------------------------------------------------------------------------
!BOC
      ctmp = Trim (block_name)

      iblock = 0

      do ii = 1, tblock
         if (TRIM(ctmp) == TRIM(list_blocknames(ii))) then
            iblock = ii
            EXIT
         end if
      end do

      if (iblock == 0) then
         write(LIS_logunit,*) "[ERR] Stopping in Ftiming_Off in "//TRIM(ctmp)
         call LIS_verify(-1, 'Stopping in Ftiming_Off in '//TRIM(ctmp))
      end if

      wclk = Mpi_Wtime (ierwtime)

      accum(iblock)%usr = accum(iblock)%usr + wclk - last(iblock)%usr
      accum(iblock)%sys = 0.0d0

      last(iblock)%usr  = wclk
      last(iblock)%sys  = 0.0d0

      return

      end subroutine Ftiming_Off
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ftiming_Reset
!
! !INTERFACE:
!
      subroutine Ftiming_Reset (block_name)
!
! !INPUT PARAMETERS:
      character (len=*)  :: block_name
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      character (len=MAX_BLK_NAME) :: ctmp
      integer                      :: iblock
      integer                      :: ii
!EOP
!-------------------------------------------------------------------------------
!BOC
      ctmp = Trim (block_name)

      iblock = 0

      do ii = 1, tblock

        if (TRIM(ctmp) == TRIM(list_blocknames(ii))) then
          iblock = ii
          EXIT
        end if

      end do

      if (iblock == 0) then
         write(LIS_logunit,*) "[ERR] Stopping in Ftiming_Reset in "//TRIM(ctmp)
         call LIS_verify(-1, 'Stopping in Ftiming_Reset in '//TRIM(ctmp))
      end if

      accum(iblock)%usr = 0.0d0
      accum(iblock)%sys = 0.0d0

      last (iblock)%usr = 0.0d0
      last (iblock)%sys = 0.0d0

      return

      end subroutine Ftiming_Reset
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ftiming_Output
!
! !INTERFACE:
!
      subroutine Ftiming_Output (lu)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: lu
!
! !DESCRIPTION:

! !LOCAL VARIABLES:
      character (len=8)  :: cdate
      character (len=10) :: ctime
      character (len=33) :: cfilename

      integer :: ierr
      integer :: ii, jj, kk
      integer :: ijk
      integer :: ind
      integer :: nn

      real(dp)  :: cijk, sijk
      real(dp)  :: onpi, onpij, onpj

      real(dp), allocatable  :: ctavg (:)
      real(dp), allocatable  :: ctmax (:)
      real(dp), allocatable  :: ctmin (:)

      real(dp), allocatable  :: ctavgj(:, :)
      real(dp), allocatable  :: ctmaxj(:, :)
      real(dp), allocatable  :: ctminj(:, :)

      real(dp), allocatable  :: ctavgi(:, :)
      real(dp), allocatable  :: ctmaxi(:, :)
      real(dp), allocatable  :: ctmini(:, :)

      real(dp), allocatable  :: us_glob1(:, :)
      real(dp), allocatable  :: us_glob2(:, :)
      integer :: numWorkerProcs, numLonProcs, numLatProcs
      integer :: commWorld
      character(len=4)     :: cnx, cny
!EOP
!-------------------------------------------------------------------------------
!BOC
      commWorld      = LIS_mpi_comm
      IF ((LIS_rc%npesx > 0) .AND. (LIS_rc%npesy > 0)) THEN
         numLonProcs    = LIS_rc%npesx
         numLatProcs    = LIS_rc%npesy
      ELSE
         numLonProcs = 0
         numLatProcs = 0

         IF ( LIS_masterproc ) THEN
            numLonProcs = LIS_rc%procLayoutx
            numLatProcs = LIS_rc%procLayouty
         ENDIF

         call MPI_Bcast(numLonProcs, 1, MPI_INTEGER, 0, LIS_mpi_comm, ierr)
         call MPI_Bcast(numLatProcs, 1, MPI_INTEGER, 0, LIS_mpi_comm, ierr)
      ENDIF
      numWorkerProcs = LIS_npes

      allocate (ctavg(tblock))
      allocate (ctmax(tblock))
      allocate (ctmin(tblock))

      allocate (ctavgj(tblock, numLonProcs))
      allocate (ctmaxj(tblock, numLonProcs))
      allocate (ctminj(tblock, numLonProcs))

      allocate (ctavgi(tblock, numLatProcs))
      allocate (ctmaxi(tblock, numLatProcs))
      allocate (ctmini(tblock, numLatProcs))

      allocate (us_glob1(tblock, numWorkerProcs))
      allocate (us_glob2(tblock, numWorkerProcs))

      us_glob1(:,:) = 0.0d0
      us_glob2(:,:) = 0.0d0

      do nn = 1, tblock
         us_tmp1(nn,1) = accum(nn)%usr
         us_tmp1(nn,2) = accum(nn)%sys

         us_glob1(nn,LIS_localPet+1) = us_tmp1(nn,1) + us_tmp1(nn,2)
      end do

      call Mpi_Allreduce (us_glob1, us_glob2, tblock*numWorkerProcs, &
                          MPI_REAL8, MPI_SUM, commWorld, ierr)

      onpi  = 1.0d0 / numLonProcs
      onpj  = 1.0d0 / numLatProcs
      onpij = 1.0d0 / numWorkerProcs

!      do nn = 1, NBLKS
      do nn = 1, tblock
         ctavg(nn) = us_glob2(nn,1)
         ctmax(nn) = us_glob2(nn,1)
         ctmin(nn) = us_glob2(nn,1)

         do kk = 2, numWorkerProcs
            ctavg(nn) = ctavg(nn) + us_glob2(nn,kk)
            ctmax(nn) = Max (ctmax(nn), us_glob2(nn,kk))
            ctmin(nn) = Min (ctmin(nn), us_glob2(nn,kk))
         end do

         ctavg(nn) = onpij * ctavg(nn)

         do jj = 1, numLatProcs
            ind = (jj - 1) * numLonProcs + 1

            ctavgi(nn,jj) = us_glob2(nn,ind)
            ctmaxi(nn,jj) = us_glob2(nn,ind)
            ctmini(nn,jj) = us_glob2(nn,ind)

            do ii = 2, numLonProcs
               ind = (jj - 1) * numLonProcs + ii

               ctavgi(nn,jj) = ctavgi(nn,jj) + us_glob2(nn,ind)
               ctmaxi(nn,jj) = Max (ctmaxi(nn,jj), us_glob2(nn,ind))
               ctmini(nn,jj) = Min (ctmini(nn,jj), us_glob2(nn,ind))
            end do

            ctavgi(nn,jj) = onpi * ctavgi(nn,jj)
         end do

         do ii = 1, numLonProcs
            ctavgj(nn,ii) = us_glob2(nn,ii)
            ctmaxj(nn,ii) = us_glob2(nn,ii)
            ctminj(nn,ii) = us_glob2(nn,ii)

            do jj = 2, numLatProcs
               ind = (jj - 1) * numLonProcs + ii

               ctavgj(nn,ii) = ctavgj(nn,ii) + us_glob2(nn,ind)
               ctmaxj(nn,ii) = Max (ctmaxj(nn,ii), us_glob2(nn,ind))
               ctminj(nn,ii) = Min (ctminj(nn,ii), us_glob2(nn,ind))
            end do

            ctavgj(nn,ii) = onpj * ctavgj(nn,ii)
         end do
      end do

      flush (6)

      call Mpi_Barrier (commWorld, ierr)

      if( LIS_masterproc ) PRINT'(a,2i5)', "Writing results ", tblock, NBLKS
      if ( LIS_masterproc ) then

         call Date_And_Time (cdate, ctime)
         write(cnx, '(i4.4)') numLonProcs
         write(cny, '(i4.4)') numLatProcs

         cfilename( 1:8)  = 'ftiming_'
         cfilename( 9:16) = cdate(1:8)
         cfilename(17:17) = '_'
         cfilename(18:23) = ctime(1:6)
         cfilename(24:33) = '_'//cnx//'x'//cny

         if (lu /= 6) THEN
            Open (lu, file = cfilename)

            Write (lu,*)
            Write (lu,*)
            Write (lu,*) '-------------------------------------------------'
            Write (lu,*) '           Beginning of Timing Statistics'
            Write (lu,*) '-------------------------------------------------'
            Write (lu,*)
            Write (lu,*)
            Write (lu,*) 'Timing statistics are below, on a per-process basis for'
            Write (lu,*) 'various sections of code. When viewed over a very short'
            Write (lu,*) 'period of time, these statistics are useful for analyzing'
            Write (lu,*) 'load imbalances. However, for longer periods of time, the'
            Write (lu,*) 'load distribution will vary, and simply comparing the'
            Write (lu,*) 'times among the various tasks might no longer be valid.'
            Write (lu,*) 'Instead, one must compare the MPI barrier time to the'
            Write (lu,*) 'overall LIS time to get a sense of the effect of'
            Write (lu,*) 'load imbalances.'
            Write (lu,*)
            Write (lu,*) 'The detail below contains data with respect to al MPI'
            Write (lu,*) 'tasks as well as latitudinal variations at fixed'
            Write (lu,*) 'longitude and longitudinal variations at fixed latitude.'
            Write (lu,*)
            Write (lu,*) '-------------------------------------------------'
            Write (lu,*)  &
             '     Block                       Min Time    Max Time    ',  &
             'Avg Time'
            Write (lu,*) '-------------------------------------------------'

            do nn = 1, tblock
               Write (lu,900) list_blocknames(nn), ctmin(nn), ctmax(nn), ctavg(nn)
            end do

 900        format (3x, a24, 3x, 3f12.4)

   IF (.FALSE.) THEN
            Write (lu,*)
            Write (lu,*) '-------------------------------------------------'
            Write (lu,*)  &
     &       '                    Longitudinal Accounting                 '
            Write (lu,*) '-------------------------------------------------'

            do ii = 1, numLonProcs
               Write (lu,*)
               Write (lu,910) ii
               Write (lu,*)
   
               do nn = 1, tblock
                  Write (lu,900) list_blocknames(nn), ctminj(nn,ii), ctmaxj(nn,ii), ctavgj(nn,ii)
               end do
            end do

 910        format ('  i = ', i5, '                            ',  &
                   'Min         Max         Avg')

            Write (lu,*)
            Write (lu,*)  &
             '  ---------------------------------------------------------',  &
             '--------'
            Write (lu,*)  &
             '                    Latitudinal Accounting                  '
            Write (lu,*)  &
             '  ---------------------------------------------------------',  &
             '--------'

            do jj = 1, numLatProcs
               Write (lu,*)
               Write (lu,920) jj
               Write (lu,*)
   
               do nn = 1, tblock
                  Write (lu,900) list_blocknames(nn), ctmini(nn,jj), ctmaxi(nn,jj), ctavgi(nn,jj)
               end do
            end do

 920        format ('  j = ', i5, '                            ',  &
     &             'Min         Max         Avg')

            Write (lu,*)
            Write (lu,*)  &
     &       '  ---------------------------------------------------------',  &
     &       '------------------'
            Write (lu,*)  &
     &       '                    Per-Process Detail                      '
            Write (lu,*)  &
     &       '  ---------------------------------------------------------',  &
     &       '------------------'

            do nn = 1, tblock
               Write (lu,990) list_blocknames(nn), (us_glob2(nn,kk),kk=1,numWorkerProcs)
            end do

 990        format (3x, a24, 3x, 4f12.4, /, (30x, 4f12.4))

            Write (lu,*)
            Write (lu,*)
            Write (lu,*) '-------------------------------------------------'
            Write (lu,*) '           End of Timing Statistics'
            Write (lu,*) '-------------------------------------------------'
            Write (lu,*)
            Write (lu,*)

            flush (lu)

            Close (lu)

         else if (lu == 6) then
            !--------------------------------------------------------------------
            ! Unnecessary calculations in order to allow time for unit 6 to flush.
            !--------------------------------------------------------------------

            ijk  = 1000

            cijk = 1.0d0
            sijk = 1.0d0 / ijk

            do kk = 1, ijk
               do jj = 1, ijk
                  do ii = 1, ijk
                     cijk = cijk * Sin (sijk*ii*jj*kk)
                  end do
               end do
            end do

            call Ftiming_Nothing (cijk)
         endif
      end if
  ENDIF

      deallocate (ctavg)
      deallocate (ctmax)
      deallocate (ctmin)

      deallocate (ctavgj)
      deallocate (ctmaxj)
      deallocate (ctminj)

      deallocate (ctavgi)
      deallocate (ctmaxi)
      deallocate (ctmini)

      deallocate (us_glob1)
      deallocate (us_glob2)

      call Mpi_Barrier (commWorld, ierr)

      return

      end subroutine Ftiming_Output
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ftiming_Nothing
!
! !INTERFACE:
!
      subroutine Ftiming_Nothing (aa)
!
! !INPUT PARAMETERS:
      real(dp)  :: aa
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      real(dp)  :: bb
!EOP
!-------------------------------------------------------------------------------
!BOC
      bb = Exp (aa)

      if (bb < 0.0d0) Write (6,*) bb

      return

      end subroutine Ftiming_Nothing
!EOC
!-------------------------------------------------------------------------------
#else
subroutine Ftiming_Init
   write(LIS_logunit, *) "[WARN] Ftiming requires MPI."
end subroutine Ftiming_Init

subroutine Ftiming_On(block_name)
   implicit none
   character (len=*)  :: block_name
   return
end subroutine Ftiming_On

subroutine Ftiming_Off(block_name)
   implicit none
   character (len=*)  :: block_name
   return
end subroutine Ftiming_Off

subroutine Ftiming_Reset(block_name)
   implicit none
   character (len=*)  :: block_name
   return
end subroutine Ftiming_Reset

subroutine Ftiming_Output(lu)
   implicit none
   integer, intent(in) :: lu
   return
end subroutine Ftiming_Output
#endif

      end module LIS_ftimingMod

