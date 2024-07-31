!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
      MODULE parutilitiesmodule
#if defined( SPMD )
!BOP
!
! !MODULE: parutilitiesmodule
!
! !USES:
      USE precision
      USE LIS_mpiMod
#include "debug.h"
      IMPLICIT NONE
#include "pilgrim.h"

!
! !PUBLIC DATA MEMBERS:
#if defined(USE_ARENAS)
      COMMON /ARENA/  buf01, buf02, buf03
      POINTER(buf01,volume)
      INTEGER(i4), DIMENSION(MAX_PE,MAX_PE) :: volume
      POINTER(buf02,databuf)
      REAL(r8), DIMENSION(MAX_BUF) :: databuf
      POINTER(buf03,intbuf)
      INTEGER(i4), DIMENSION(MAX_BUF)  :: intbuf
#endif
      PUBLIC     CommGlobal, GID, Gsize
      PUBLIC     SUMOP, MAXOP, MINOP

      INTEGER,SAVE :: CommGlobal   ! Global communicator (before ParSplit)
      INTEGER,SAVE :: GSize        ! Size of communicator CommGlobal
      INTEGER,SAVE :: GID          ! My rank in communicator CommGlobal

#define CPP_SUM_OP 101
#define CPP_MAX_OP 102
#define CPP_MIN_OP 103
#define CPP_BCST_OP 104
#if defined( USE_ARENAS )
      INTEGER,SAVE :: SUMOP = CPP_SUM_OP
      INTEGER,SAVE :: MAXOP = CPP_MAX_OP
      INTEGER,SAVE :: MINOP = CPP_MIN_OP
      INTEGER,SAVE :: BCSTOP = CPP_BCST_OP
#else
      INTEGER,SAVE :: SUMOP = MPI_SUM
      INTEGER,SAVE :: MAXOP = MPI_MAX
      INTEGER,SAVE :: MINOP = MPI_MIN
      INTEGER,SAVE :: BCSTOP = CPP_BCST_OP
#endif

      INTEGER,SAVE :: numcpu, numcps( MAX_PE ), blocksize, packetsize

! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC ParPatternType

       PRIVATE BlockDescriptor
       TYPE BlockDescriptor
         INTEGER, POINTER     :: Displacements(:)   ! Offsets in local segment
         INTEGER, POINTER     :: BlockSizes(:)      ! Block sizes to transfer
       END TYPE BlockDescriptor
 
      TYPE ParPatternType
        INTEGER ::     Comm                  ! Communicator
        INTEGER ::     Iam                   ! My rank in communicator
        INTEGER ::     Size                  ! Size of communicator
#if defined( USE_ARENAS )
        TYPE(BlockDescriptor), POINTER :: SendDesc(:) ! Array of descriptors
        TYPE(BlockDescriptor), POINTER :: RecvDesc(:) ! Array of descriptors
#else
        INTEGER, POINTER :: SendDesc( : )    ! Send descriptors
        INTEGER, POINTER :: RecvDesc( : )    ! Receive descriptors
#endif
      END TYPE ParPatternType 

      PUBLIC     ParInit, ParSplit, ParFree, ParExit
      PUBLIC     ParScatter, ParGather
      PUBLIC     ParBeginTransfer, ParEndTransfer
      PUBLIC     ParExchangeVector, ParCollective
      PUBLIC     ParPatternCreate, ParPatternFree

      INTERFACE     ParPatternCreate
        MODULE PROCEDURE ParPatternGhost
        MODULE PROCEDURE ParPatternDecompToDecomp
!!!        MODULE PROCEDURE ParPatternDecompToGhost
!!!        MODULE PROCEDURE ParPatternGhostToDecomp
!!!        MODULE PROCEDURE ParPatternGhostToGhost
      END INTERFACE
 
      INTERFACE     ParScatter
        MODULE PROCEDURE ParScatterReal
        MODULE PROCEDURE ParScatterInt
      END INTERFACE
 
      INTERFACE     ParGather
        MODULE PROCEDURE ParGatherReal
        MODULE PROCEDURE ParGatherInt
      END INTERFACE

      INTERFACE     ParBeginTransfer
        MODULE PROCEDURE ParBeginTransferReal
        MODULE PROCEDURE ParBeginTransferPattern
!        MODULE PROCEDURE ParBeginTransferInt
      END INTERFACE

      INTERFACE     ParEndTransfer
        MODULE PROCEDURE ParEndTransferReal
        MODULE PROCEDURE ParEndTransferPattern
!        MODULE PROCEDURE ParEndTransferInt
      END INTERFACE

      INTERFACE     ParExchangeVector
        MODULE PROCEDURE ParExchangeVectorReal
        MODULE PROCEDURE ParExchangeVectorInt
      END INTERFACE

      INTERFACE     ParCollective
        MODULE PROCEDURE ParCollectiveBarrier
        MODULE PROCEDURE ParCollective0D
        MODULE PROCEDURE ParCollective1D
        MODULE PROCEDURE ParCollective2D
        MODULE PROCEDURE ParCollective3D
        MODULE PROCEDURE ParCollective0DInt
        MODULE PROCEDURE ParCollective1DInt
      END INTERFACE

!
! !DESCRIPTION:
!
!      This module provides the basic utilities to support parallelism
!      on a distributed or shared memory multiprocessor.
!
!      \begin{center}
!      \begin{tabular}{|l|l|} \hline \hline
!        ParInit           & Initialize the parallel system \cr \hline
!        ParExit           & Exit from the parallel system \cr \hline
!        ParSplit          & Create a Compute grid of PEs   \cr \hline
!        ParFree           & Free a split communicator \cr \hline
!        ParScatter        & Scatter global slice to local slices \cr \hline
!        ParGather         & Gather local slices to one global \cr \hline
!        ParBeginTransfer  & Initiate an all-to-all packet transfer \cr \hline
!        ParEndTransfer    & Complete an all-to-all packet transfer \cr \hline
!        ParExchangeVector & Complete an all-to-all packet transfer \cr \hline
!        ParCollective     & Collective operation across communicator \cr \hline
!      \end{tabular}
!      \end{center}
!      \vspace{2mm}
!
!      Other utilities can be added to this module as needs evolve.
!
!      Conceptually the intention is to aggregate as many of the
!      MPI communication calls as possible into a well-maintained
!      module.  This will help avoid the occurrence of MPI spaghetti 
!      code.  
!
!      This module is tailored to GEOS DAS and implements the 
!      design of Lucchesi/Mirin/Sawyer/Larson.
!
! !REVISION HISTORY:
!   97.02.01   Sawyer     Creation
!   97.07.22   Sawyer     Removal of DecompType related subroutines
!   97.08.13   Sawyer     Added ParScatter/Gather for Integers
!   97.09.26   Sawyer     Additions of Sparse communication primitives
!   97.12.01   Sawyer     Changed all MPI_SSEND to MPI_ISEND
!   97.12.23   Lucchesi   Added member variables IsIONode and InterComm
!   98.01.06   Sawyer     Additions from RL for I/O Nodes
!   98.02.02   Sawyer     Added the Cartesian data members
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.02.23   Sawyer     Added ghosting utilities
!   98.02.25   Sawyer     Modified interface of BeginTransfer
!   98.03.03   Sawyer     Added Global ID number to public data members
!   98.03.25   Sawyer     Added documentation for walkthrough
!   98.04.16   Sawyer     Removed all use of MPI_CART (CommRow redefined)
!   98.07.23   Sawyer     Added ParGhost, ParPoleDot; ParBegin/EndGhost out
!   98.09.15   Sawyer     Added ParMerge, ParPoleGhost
!   98.09.17   Sawyer     Added ParSum, removed ParPoleDot
!   99.01.18   Sawyer     Minor cleaning
!   99.03.04   Sawyer     Revised SHMEM concept for Transfer
!   99.04.22   Sawyer     Removed COMMON for handles -- they are
!                         always used in same program unit.
!   99.05.21   Sawyer     Reintroduced barriers in Scatter/Gather
!   99.06.03   Sawyer     USE_SHMEM revisions
!   99.12.10   Sawyer     ParInit now sets GID, Gsize
!   99.12.13   Sawyer     Version slimmed down for FVCCM release
!   00.06.14   Sawyer     Precision module now used
!   00.07.07   Sawyer     Removed 2D scatter/gather; simplified API
!   00.07.30   Sawyer     Full implementation with shared memory
!   00.08.09   Sawyer     Replaced ParSum with ParCollective
!   00.08.28   Sawyer     Moved LLNL 2D data to LLNL2DModule; new MLP impl
!   01.02.04   Sawyer     Added PatternType and related routines
!   01.02.12   Sawyer     Converted to free format
!
! !BUGS:
!   There are several MPI_Barriers at locations in the code.
!   These avoid potential race conditions which probably only occur
!   if the number of real processors is less than the number of
!   message passing processes.  Remove these barriers at your own risk
!
!EOP

      INTEGER, SAVE :: Inhandle(MAX_PAX,MAX_TRF), OutHandle(MAX_PAX,MAX_TRF)
      INTEGER, SAVE :: BegTrf = 0  ! Ongoing overlapped begintransfer # 
      INTEGER, SAVE :: EndTrf = 0  ! Ongoing overlapped endtransfer #

      CONTAINS
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParInit --- Initialize the parallel execution
!
! !INTERFACE: 
      SUBROUTINE ParInit (  )
!
! !USES:
      IMPLICIT NONE
!
! !DESCRIPTION:
!     Initializes the system.  In MPI mode, call MPI\_INIT if not done 
!     already.  In USE\_ARENAS mode, initialize the shared memory buffer.
!
!     This routine is the very {\em first} thing which is executed!
!
! !SYSTEM ROUTINES:
!     MPI_INITIALIZED, MPI_INIT
!
! !REVISION HISTORY:
!   97.03.20   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.07.03   Sawyer     Reformulated documentation
!   00.07.23   Sawyer     Added shared memory arena implementation
!
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
      INTEGER Ierror
      LOGICAL Flag
!
#if defined(USE_ARENAS)
#include <ulocks.h>
      INTEGER(i4) :: ipe, fork, getpid, master
      INTEGER(i8) :: nvars, extent(100), pnt(100)
      character*80 evalue

! Get the memory for the global variables 
      extent(1) = MAX_PE * MAX_PE * 4
      extent(2) = MAX_BUF * 8
      extent(3) = MAX_BUF * 4
      nvars = 3
      call getmem(nvars,extent,pnt)
      buf01=pnt(1)
      buf02=pnt(2)
      buf03=pnt(3)

! Get the number of processes
      call getenv('N_MPI',evalue)
      read(evalue,*) Gsize

! Get the max number of threads per process
      call getenv('N_SMP',evalue)
      read(evalue,*) numcpu
      do ipe=1,Gsize
        numcps(ipe)=numcpu
      enddo
!!!      print *, "Starting MLP PILGRIM with", Gsize, "processes and",
!!!     &         numcpu, "threads per process" 

! Calculate maximum blocksize and packetsize
      blocksize = MAX_BUF / Gsize
      packetsize = blocksize / MAX_PAX

! Destroy and recreate the environment
      master=getpid()
      call mp_destroy
      gid = 0
      do while ( getpid() .eq. master .and. gid < Gsize-1 )
        ierror=fork()
        gid = gid+1
      enddo
      if ( getpid() .eq. master ) gid = 0
!!!      call mp_set_threads(numcps(gid+1))
!
! Jim has more code here to pin the processes to the nodes
!
#else
!
!     Check if MPI is initialized.  If not, initialize.  No mpi_call
!
      CALL MPI_INITIALIZED( Flag, Ierror )
      CPP_ASSERT_F90( Ierror == 0 )
      IF ( .not. Flag ) then
        CALL MPI_INIT( ierror )
        CPP_ASSERT_F90( Ierror == 0 )
      ENDIF
      CALL MPI_COMM_SIZE( LIS_mpi_comm, Gsize, Ierror )
      CALL MPI_COMM_RANK( LIS_mpi_comm, GID, Ierror )
      CALL MPI_COMM_DUP( LIS_mpi_comm, CommGlobal, Ierror )
#endif
      RETURN
!EOC
      END SUBROUTINE ParInit
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParExit --- Finalize the parallel execution
!
! !INTERFACE:
      SUBROUTINE ParExit ( )

! !USES:
      IMPLICIT NONE


! !DESCRIPTION:
!     All PEs, compute nodes and IO nodes alike meet here to terminate
!     themselves.  If someone does not check in, everything will hang
!     here.
!
!     This routine is the very {\em last} thing which is executed!
!
! !LOCAL VARIABLES:
      INTEGER Ierror
!
! !SYSTEM ROUTINES:
!     MPI_BARRIER, MPI_FINALIZE
!
! !REVISION HISTORY:
!   97.03.20   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.07.03   Sawyer     Reformulated documentation
!   00.07.23   Sawyer     Added shared memory arena implementation
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined( USE_ARENAS )
      CALL MPI_BARRIER( LIS_mpi_comm, Ierror )
      CALL MPI_FINALIZE( Ierror )
#endif
      RETURN
!EOC
      END SUBROUTINE ParExit
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParSplit --- Split into group for I/O and computation
!
! !INTERFACE:
      SUBROUTINE ParSplit( InComm, Color, InID, Comm, MyID, Nprocs )
!
! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )     :: InComm    ! Communicator to split
      INTEGER, INTENT( IN )     :: Color     ! Group label
      INTEGER, INTENT( IN )     :: InID      ! Input ID

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )    :: Comm      ! Split communicator
      INTEGER, INTENT( OUT )    :: MyID      ! Group label
      INTEGER, INTENT( OUT )    :: Nprocs    ! Number of PEs in my group
!
! !DESCRIPTION:
!     This routine splits the PEs into groups.  This is currently only
!     supported in MPI mode. Read the chapter on MPI\_COMM\_SPLIT 
!     thoroughly.  
!
! !SYSTEM ROUTINES:
!     MPI_COMM_SPLIT, MPI_COMM_SIZE, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.03.20   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.07.03   Sawyer     Reformulated documentation
!   97.12.01   Sawyer     Xnodes and Ynodes are explicit arguments
!   97.12.23   Lucchesi   Added call to MPI_INTERCOMM_CREATE
!   98.01.06   Sawyer     Additions from RL for I/O Nodes
!   98.02.02   Sawyer     Added the Cartesian information
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.04.16   Sawyer     Removed all use of MPI_CART (CommRow redefined)
!   99.01.10   Sawyer     CommRow now defined for all rows
!   00.07.09   Sawyer     Removed 2D computational mesh
!   00.08.08   Sawyer     Redefined as wrapper to mpi_comm_split
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  Ierror

      CPP_ENTER_PROCEDURE( "PARSPLIT" )
#if !defined( USE_ARENAS )
!
!     Split the communicators
!
      CALL MPI_COMM_SPLIT( InComm, Color, InID, Comm, Ierror )
      IF ( Comm .ne. MPI_COMM_NULL ) THEN
        CALL MPI_COMM_RANK( Comm, MyID, Ierror )
        CALL MPI_COMM_SIZE( Comm, Nprocs, Ierror )
      ELSE
!
!     This PE does not participate: mark with impossible values
!
        MyID = -1
        Nprocs = -1
      ENDIF
#endif

      CPP_LEAVE_PROCEDURE( "PARSPLIT" )
      RETURN
!EOC
      END SUBROUTINE ParSplit
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParFree --- Free a communicator
!
! !INTERFACE:
      SUBROUTINE ParFree( InComm ) 
!
! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER InComm

!
! !DESCRIPTION:
!     This routine frees a communicator created with ParSplit
!
! !REVISION HISTORY:
!   97.09.11   Sawyer     Creation, to complement ParSplit
!   00.07.24   Sawyer     Revamped ParMerge into a free communicator 
!
! !LOCAL VARIABLES:
      INTEGER  Ierror
!
!EOP
!-----------------------------------------------------------------------
!BOC
      CPP_ENTER_PROCEDURE( "PARFREE" )
!
#if !defined( USE_ARENAS )
      CALL MPI_COMM_FREE( InComm, Ierror ) 
#endif
      CPP_LEAVE_PROCEDURE( "PARFREE" )
      RETURN
!EOC
      END SUBROUTINE ParFree
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternGhost --- Create pattern for given ghosting
!
! !INTERFACE:
      SUBROUTINE ParPatternGhost( InComm, Ghost, Pattern )
!
! !USES:
      USE decompmodule, ONLY : DecompGlobalToLocal, DecompLocalToGlobal
      USE ghostmodule, ONLY : GhostType, GhostInfo
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )               :: InComm  ! # of PEs
      TYPE(GhostType),  INTENT( IN )       :: Ghost   ! # of PEs
! !OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( OUT )  :: Pattern ! Comm Pattern
!
! !DESCRIPTION:
!     This routine contructs a communication pattern from the ghost
!     region definition.  That is, the resulting communication pattern
!     can be used in ParBegin/EndTransfer with the ghosted arrays as
!     inputs.  
!
! !SYSTEM ROUTINES:
!     MPI_TYPE_INDEXED
!
! !REVISION HISTORY:
!   01.02.10   Sawyer     Creation
!   01.06.02   Sawyer     Renamed ParPatternGhost
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  i, j, ipe, pe, Iam, GroupSize, Num, Length, Ptr, Ierror
      INTEGER  Global, End, Local, GlobalSize, LocalSize, BorderSize
      INTEGER, ALLOCATABLE :: InVector(:), OutVector(:)
      INTEGER, ALLOCATABLE :: LenInVector(:), LenOutVector(:)

      CPP_ENTER_PROCEDURE( "PARPATTERNGHOST" )

!
! First request the needed ghost values from other processors.
!
#if !defined( USE_ARENAS )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
#endif

      ALLOCATE( Pattern%SendDesc( GroupSize ) )
      ALLOCATE( Pattern%RecvDesc( GroupSize ) )

!
! Temporary variables
!
      ALLOCATE( LenInVector( GroupSize ) )
      ALLOCATE( LenOutVector( GroupSize ) )

      CALL GhostInfo( Ghost,GroupSize,GlobalSize,LocalSize,BorderSize )
      ALLOCATE( InVector( 2*BorderSize ) )
      ALLOCATE( OutVector( 2*LocalSize ) )

!
! A rather complicated loop to define the local ghost region.
! The concept is the following:  go through all the points in the
! border data structure.   It contains global indices of the points
! which have to be copied over from neighboring PEs.  These indices
! are collected into InVector for transmission to those PEs, in
! effect informing them of the local PEs requirements.
!
! A special case is supported:  if the ghost domain wraps around
! onto the domain of the local PE!  This is very tricky, because
! the index space in both Ghost%Border and Ghost%Local MUST be
! unique for DecompGlobalToLocal to work.   Solution:  ghost 
! points are marked with the negative value of the needed domain 
! value in both Ghost%Border and Ghost%Local.  These are "snapped 
! over" to the true global index with the ABS function, so that 
! they can be subsequently found in the true local domain.
!
      j = 1
      DO ipe=1, GroupSize
        Num = SIZE(Ghost%Border%Head(ipe)%StartTags)
        Length = 0
        DO i = 1, Num
          Global = Ghost%Border%Head(ipe)%StartTags(i)
          IF ( Global /= 0 ) THEN
            Length = Length + 1
            End    = Ghost%Border%Head(ipe)%EndTags(i)
            InVector(j) = ABS(Global)
            InVector(j+1) = ABS(End)
            CALL DecompGlobalToLocal( Ghost%Local, Global, Local, Pe )
            OutVector(Length) = Local-1                ! Zero-based address
            OutVector(Length+Num) = End - Global+1     ! Chunk size
            j = j + 2
          ENDIF
        ENDDO
        LenInVector( ipe ) = 2*Length

!
! Set the receive buffer descriptor
!
#if defined(DEBUG_PARPATTERNGHOST)
        print *,"Iam",Iam,"Pe",Ipe-1,"Lens",OutVector(Num+1:Num+Length), &
     &       "Displacements", OutVector(1:Length)
#endif
#if defined( USE_ARENAS )
! This code is currently untested
         ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(Num) )
         ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(Num) )
         DO i=1, Num
           Pattern%RecvDesc(ipe)%Displacements(i) = OutVector(i)
           Pattern%RecvDesc(ipe)%BlockSizes(i)    = OutVector(Num+i)
         ENDDO            
#else
        CALL MPI_TYPE_INDEXED( Length, OutVector(Num+1), OutVector,      &
     &     CPP_MPI_REAL, Ptr, Ierror )
        CALL MPI_TYPE_COMMIT( Ptr, Ierror )
        Pattern%RecvDesc( ipe ) = Ptr
#endif
      ENDDO

!
! Everybody exchanges the needed information
!
#if defined(DEBUG_PARPATTERNGHOST)
      print *, "iam", iam, "In", LenInVector,                            &
     &           InVector( 1:SUM(LenInVector) )
#endif
      CALL ParExchangeVectorInt( InComm, LenInVector, InVector,          &
     &                                LenOutVector, OutVector )
#if defined(DEBUG_PARPATTERNGHOST)
      print *, "iam", iam, "Out", LenOutVector,                          &
     &           OutVector( 1:SUM(LenOutVector) )
#endif
!
! Now everyone has the segments which need to be sent to the 
! immediate neighbors.  Save these in PatternType.
!
      j = 1
      DO ipe = 1, GroupSize
        Num = LenOutVector(ipe) / 2
        DO i = 1, Num
          CALL DecompGlobalToLocal( Ghost%Local,OutVector(j),Local,pe )
          InVector(i) = Local-1
          InVector(i+Num) = OutVector(j+1) - OutVector(j) + 1
          j = j + 2
        ENDDO
#if defined(DEBUG_PARPATTERNGHOST)
        print *, "Iam", Iam, "To", ipe-1, "InVector",                    &
     &        InVector(1:Num), "block size", InVector(Num+1:2*Num)
#endif
#if defined( USE_ARENAS )
         ALLOCATE( Pattern%SendDesc(ipe)%Displacements(Num) )
         ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(Num) )
         DO i=1, Num
           Pattern%SendDesc(ipe)%Displacements(i) = InVector(i)
           Pattern%SendDesc(ipe)%BlockSizes(i)    = InVector(Num+i)
         ENDDO            
#else
        CALL MPI_TYPE_INDEXED( Num, InVector(Num+1), InVector,           &
     &      CPP_MPI_REAL, Ptr, Ierror )
        CALL MPI_TYPE_COMMIT( Ptr, Ierror )
        Pattern%SendDesc( ipe ) = Ptr
#endif
      ENDDO

!
! Copy the communicator into the pattern
!
      CALL MPI_COMM_DUP( InComm, Pattern%Comm, Ierror )
      CALL MPI_COMM_RANK( InComm, Pattern%Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, Pattern%Size, Ierror )
!
! Clean up the locally allocate variables
!
      DEALLOCATE( OutVector )
      DEALLOCATE( InVector )
      DEALLOCATE( LenOutVector )
      DEALLOCATE( LenInVector )

      CPP_LEAVE_PROCEDURE( "PARPATTERNGHOST" )
      RETURN
!EOC
      END SUBROUTINE ParPatternGhost
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternDecompToDecomp --- Create pattern between decomps
!
! !INTERFACE:
      SUBROUTINE ParPatternDecompToDecomp( InComm, DA, DB, Pattern )
!
! !USES:
      USE decompmodule, ONLY : DecompType, DecompGlobalToLocal, DecompInfo
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )               :: InComm  ! # of PEs
      TYPE(DecompType),  INTENT( IN )      :: DA      ! Source Decomp Desc
      TYPE(DecompType),  INTENT( IN )      :: DB      ! Target Decomp Desc
! !OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( OUT )  :: Pattern ! Comm Pattern
!
! !DESCRIPTION:
!     This routine contructs a communication pattern for a 
!     transformation from one decomposition to another, i.e., a 
!     so-called "transpose". The resulting communication pattern 
!     can be used in ParBegin/EndTransfer with the decomposed 
!     arrays as inputs.  
!
! !SYSTEM ROUTINES:
!
! !BUGS:
!     Under development
!
! !REVISION HISTORY:
!   01.05.29   Sawyer     Creation from RedistributeCreate
!   01.07.13   Sawyer     Rewritten to minimize DecompGlobalToLocal
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      LOGICAL NewIpe
      INTEGER I, J, Tag, Local, Pe, LenB, JB, Ipe, Num, Inc, Off
      INTEGER Ptr                                ! Pointer type
      INTEGER GroupSize, Iam, Ierror
      INTEGER OldPe, TotalPtsA, NpesA, TotalPtsB, NpesB

      INTEGER, ALLOCATABLE :: Count(:)           ! # segments for each recv PE
      INTEGER, ALLOCATABLE :: CountOut(:)        ! # segments for each send PE

      INTEGER, ALLOCATABLE :: DisplacementsA(:)  ! Generic displacements
      INTEGER, ALLOCATABLE :: BlockSizesA(:)     ! Generic block sizes
      INTEGER, ALLOCATABLE :: LocalA(:)          ! Generic Local indices

      INTEGER, ALLOCATABLE :: DisplacementsB(:)  ! Displacements for B
      INTEGER, ALLOCATABLE :: BlockSizesB(:)     ! Block sizes for B
      INTEGER, ALLOCATABLE :: LocalB(:)          ! Local indices for B
      INTEGER, ALLOCATABLE :: PeB(:)             ! Processor element numbers

      CPP_ENTER_PROCEDURE( "PARPATTERNDECOMPTODECOMP" )

      CALL DecompInfo( DA, NpesA, TotalPtsA )
      CALL DecompInfo( DB, NpesB, TotalPtsB )

#if defined( USE_ARENAS )
! Not yet supported: provoke an error
      You should not be compiling for shared memory arenas
#else
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_DUP( InComm, Pattern%Comm, Ierror )
#endif
      Pattern%Size = GroupSize
      Pattern%Iam  = Iam
!
! Allocate the number of entries and list head arrays
!
      CPP_ASSERT_F90( NpesA .EQ. GroupSize )
      CPP_ASSERT_F90( NpesB .EQ. GroupSize )

!
! Allocate the patterns
!
      ALLOCATE( Pattern%SendDesc( GroupSize ) )
      ALLOCATE( Pattern%RecvDesc( GroupSize ) )

!
! Local allocations
!
      ALLOCATE( DisplacementsA( TotalPtsA ) )   ! Allocate for worst case
      ALLOCATE( BlockSizesA( TotalPtsA ) )      ! Allocate for worst case
      ALLOCATE( LocalA( TotalPtsA ) )           ! Allocate for worst case

      ALLOCATE( DisplacementsB( TotalPtsB ) )   ! Allocate for worst case
      ALLOCATE( BlockSizesB( TotalPtsB ) )      ! Allocate for worst case
      ALLOCATE( LocalB( TotalPtsA ) )           ! Allocate for worst case
      ALLOCATE( PeB( TotalPtsB ) )              ! Allocate for worst case

      ALLOCATE( Count( GroupSize ) )
      ALLOCATE( CountOut( GroupSize ) )

      JB        = 0
      Count     = 0
      LenB      = 0

      NewIpe = .TRUE.
      Num    = 0
      Inc    = 0

!
! Parse through all the tags in the local segment
      DO J = 1, SIZE( DB%Head(iam+1)%StartTags )
        OldPe     = -1         ! Set PE undefined
        DO Tag=DB%Head(iam+1)%StartTags(J), DB%Head(iam+1)%EndTags(J)
!
! Determine the index and PE of this entry on A. This might be inlined later
!
          CALL DecompGlobalToLocal( DA, Tag, Local, Pe )

!
! If ipe-1 is my id, then this is an entry ipe will receive from Pe
!
          IF ( Pe /= OldPe ) THEN
            OldPe   = Pe
            IF ( jb > 0 ) THEN
              BlockSizesB(jb) = LenB
              LenB = 0
            ENDIF
            jb = jb+1                     ! increment the segment index
            DisplacementsB(jb) = Inc      ! Zero-based offset of local segment
            LocalB(jb) = Local-1          ! The local index (zero-based)
            PeB(jb) = Pe                  ! Note the ID of the sender
            Count(Pe+1) = Count(Pe+1)+1 ! Increment counter of segments
          ENDIF
          LenB = LenB+1                   ! Good -- segment is getting longer
          Inc = Inc+1                     ! Increment local index
        ENDDO
      ENDDO
!
! Clean up
!
      BlockSizesB(jb) = LenB
      print *, iam, "BlockSizes", BlockSizesB(1:jb), DisplacementsB(1:jb), PeB(1:jb), Count

      CPP_ASSERT_F90( JB .LE. GlobalSize )
!
! Now create the pattern from the displacements and block sizes
!
      Inc = 0
      DO ipe = 1, GroupSize
!
! Find the segments which are relevant for the sender ipe
! Make compact arrays BlockSizes and Displacements 
!
        DO j = 1, jb
          IF ( PeB(j) == ipe-1 ) THEN
            Inc = Inc + 1
            BlockSizesA(Inc) = BlockSizesB(j)
            DisplacementsA(Inc) = DisplacementsB(j)
            LocalA(Inc)      = LocalB(j)
          ENDIF
        ENDDO
      ENDDO

!
! Create the receiver communication pattern
!
      Off = 0
      DO ipe = 1, GroupSize
        Num = Count(ipe)
        print *, "Receiver Iam", Iam, "Ipe", Ipe-1, "Num", Num, &
                 "Displacements", DisplacementsA(Off+1:Off+Num), &
                 "BlockSizes", BlockSizesA(Off+1:Off+Num)
#if defined( USE_ARENAS )
        ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(Num) )
        ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(Num) )
        DO i=1, Num
          Pattern%RecvDesc(ipe)%Displacements(i) = DisplacementsA(i+Off)
          Pattern%RecvDesc(ipe)%BlockSizes(i)    = BlockSizesA(i+Off)
        ENDDO
#else
        CALL MPI_TYPE_INDEXED( Num, BlockSizesA(Off+1),DisplacementsA(Off+1), &
     &                         CPP_MPI_REAL, Ptr, Ierror )
        Pattern%RecvDesc( ipe ) = Ptr
#endif
        Off = Off + Num
      ENDDO

!
! Now communicate what the receiver is expecting to the sender
!
      CALL ParExchangeVectorInt( InComm, Count, LocalA,                 &
     &                           CountOut, DisplacementsB  )
      CALL ParExchangeVectorInt( InComm, Count, BlockSizesA,            &
     &                           CountOut, BlockSizesB )

!
! Sender A: BlockSizes and Displacements can now be stored
!
      Off = 0
      DO ipe=1, GroupSize
        Num = CountOut(ipe)
        print *, "Sender Iam", Iam, "Ipe", Ipe-1, "Num", Num,  &
                 "Displacements", DisplacementsB(Off+1:Off+Num), &
                 "BlockSizes", BlockSizesB(Off+1:Off+Num)

#if defined( USE_ARENAS )
        ALLOCATE( Pattern%SendDesc(ipe)%Displacements(Num) )
        ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(Num) )
        DO i=1, Num
          Pattern%SendDesc(ipe)%Displacements(i) = DisplacementsB(i+Off)
          Pattern%SendDesc(ipe)%BlockSizes(i)    = BlockSizesB(i+Off)
        ENDDO
#else
        CALL MPI_TYPE_INDEXED( Num, BlockSizesB(Off+1),DisplacementsB(Off+1),&
     &                         CPP_MPI_REAL, Ptr, Ierror )
        Pattern%SendDesc( ipe ) = Ptr
#endif
        Off = Off + Num
      ENDDO

      DEALLOCATE( CountOut )
      DEALLOCATE( Count )

      DEALLOCATE( PeB )
      DEALLOCATE( LocalB )
      DEALLOCATE( BlockSizesB )
      DEALLOCATE( DisplacementsB )

      DEALLOCATE( LocalA )
      DEALLOCATE( BlockSizesA )
      DEALLOCATE( DisplacementsA )

      CPP_LEAVE_PROCEDURE( "PARPATTERNDECOMPTODECOMP" )
      RETURN
!EOC
      END SUBROUTINE ParPatternDecompToDecomp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternDecompToGhost --- Create pattern decomp to ghost
!
! !INTERFACE:
      SUBROUTINE ParPatternDecompToGhost( InComm, DA, GB, Pattern )
!
! !USES:
      USE decompmodule, ONLY : DecompType, DecompGlobalToLocal,         &
                               DecompInfo
      USE ghostmodule, ONLY : GhostType, GhostInfo
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )               :: InComm  ! # of PEs
      TYPE(DecompType),  INTENT( IN )      :: DA      ! Source Ghost Desc
      TYPE(GhostType),  INTENT( IN )       :: GB      ! Target Ghost Desc
! !OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( OUT )  :: Pattern ! Comm Pattern
!
! !DESCRIPTION:
!     This routine contructs a communication pattern for a transformation
!     from decomposition to a ghosted decomposition, i.e., a so-called 
!     "transpose".  The resulting communication pattern can be used in 
!     ParBegin/EndTransfer with the decomposed arrays as inputs.  
!
! !SYSTEM ROUTINES:
!
! !BUGS:
!     Under development
!
! !REVISION HISTORY:
!   12.07.01   Sawyer     Creation from ParPatternDecompToDecomp
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      LOGICAL NewIpe
      INTEGER I, J, Tag, Local, Pe, LenB, JB, Ipe, Num, Inc, Off
      INTEGER Ptr                                ! Pointer type
      INTEGER GroupSize, Iam, Ierror
      INTEGER OldPe, TotalPtsA, NpesA
      INTEGER GlobalSizeB, LocalSizeB, BorderSizeB, NpesB

      INTEGER, ALLOCATABLE :: Count(:)           ! # segments for each recv PE
      INTEGER, ALLOCATABLE :: CountOut(:)        ! # segments for each send PE

      INTEGER, ALLOCATABLE :: DisplacementsA(:)  ! Generic displacements
      INTEGER, ALLOCATABLE :: BlockSizesA(:)     ! Generic block sizes
      INTEGER, ALLOCATABLE :: LocalA(:)          ! Generic Local indices

      INTEGER, ALLOCATABLE :: DisplacementsB(:)  ! Displacements for B
      INTEGER, ALLOCATABLE :: BlockSizesB(:)     ! Block sizes for B
      INTEGER, ALLOCATABLE :: LocalB(:)          ! Local indices for B
      INTEGER, ALLOCATABLE :: PeB(:)             ! Processor element numbers

      CPP_ENTER_PROCEDURE( "PARPATTERNDECOMPTOGHOST" )

      CALL DecompInfo( DA, NpesA, TotalPtsA )
      CALL GhostInfo( GB, NpesB, GlobalSizeB, LocalSizeB, BorderSizeB )

#if defined( USE_ARENAS )
! Not yet supported: provoke an error
      You should not be compiling for shared memory arenas
#else
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_DUP( InComm, Pattern%Comm, Ierror )
#endif
      Pattern%Size = GroupSize
      Pattern%Iam  = Iam
!
! Allocate the number of entries and list head arrays
!
      CPP_ASSERT_F90( NpesA .EQ. GroupSize )
      CPP_ASSERT_F90( NpesB .EQ. GroupSize )

!
! Allocate the patterns
!
      ALLOCATE( Pattern%SendDesc( GroupSize ) )
      ALLOCATE( Pattern%RecvDesc( GroupSize ) )

!
! Local allocations
!
      ALLOCATE( DisplacementsA( TotalPtsA ) )   ! Allocate for worst case
      ALLOCATE( BlockSizesA( TotalPtsA ) )      ! Allocate for worst case
      ALLOCATE( LocalA( TotalPtsA ) )           ! Allocate for worst case

      ALLOCATE( DisplacementsB( GlobalSizeB ) ) ! Allocate for worst case
      ALLOCATE( BlockSizesB( GlobalSizeB ) )    ! Allocate for worst case
      ALLOCATE( LocalB( GlobalSizeB ) )         ! Allocate for worst case
      ALLOCATE( PeB( GlobalSizeB ) )            ! Allocate for worst case

      ALLOCATE( Count( GroupSize ) )
      ALLOCATE( CountOut( GroupSize ) )

      JB        = 0
      Count     = 0
      LenB      = 0

      NewIpe = .TRUE.
      Num    = 0
      Inc    = 0

!
! Parse through all the tags in the local segment
      DO J = 1, SIZE( GB%Local%Head(iam+1)%StartTags )
        OldPe     = -1         ! Set PE undefined
        DO Tag=GB%Local%Head(iam+1)%StartTags(J),                         &
                GB%Local%Head(iam+1)%EndTags(J)
!
! Determine the index and PE of this entry on A. This might be inlined later
!
          CALL DecompGlobalToLocal( DA, Tag, Local, Pe )

!
! If ipe-1 is my id, then this is an entry ipe will receive from Pe
!
          IF ( Pe /= OldPe ) THEN
            OldPe   = Pe
            IF ( jb > 0 ) THEN
              BlockSizesB(jb) = LenB
              LenB = 0
            ENDIF
            jb = jb+1                     ! increment the segment index
            DisplacementsB(jb) = Inc      ! Zero-based offset of local segment
            LocalB(jb) = Local-1          ! Local indices (zero-based)
            PeB(jb) = Pe                  ! Note the ID of the sender
            Count(Pe+1) = Count(Pe+1)+1 ! Increment counter of segments
          ENDIF
          LenB = LenB+1                   ! Good -- segment is getting longer
          Inc = Inc+1                     ! Increment local index
        ENDDO
      ENDDO
!
! Clean up
!
      BlockSizesB(jb) = LenB

      CPP_ASSERT_F90( JB .LE. GlobalSize )
!
! Now create the pattern from the displacements and block sizes
!
      Inc = 0
      DO ipe = 1, GroupSize
!
! Find the segments which are relevant for the sender ipe
! Make compact arrays BlockSizes and Displacements 
!
        DO j = 1, jb
          IF ( PeB(j) == ipe-1 ) THEN
            Inc = Inc + 1
            BlockSizesA(Inc) = BlockSizesB(j)
            DisplacementsA(Inc) = DisplacementsB(j)
            LocalA(Inc)      = LocalB(j)
          ENDIF
        ENDDO
      ENDDO

      Off = 0
      DO ipe = 1, GroupSize
        Num = Count(ipe)
        print *, "Receiver Iam", Iam, "Ipe", Ipe-1, "Num", Num, &
                 "Displacements", DisplacementsA(Off+1:Off+Num), &
                 "BlockSizes", BlockSizesA(Off+1:Off+Num)

!
! Create the receiver communication pattern
!
#if defined( USE_ARENAS )
        ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(Num) )
        ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(Num) )
        DO i=1, Num
          Pattern%RecvDesc(ipe)%Displacements(i) = DisplacementsA(i+Off)
          Pattern%RecvDesc(ipe)%BlockSizes(i)    = BlockSizesA(i+Off)
        ENDDO
#else
        CALL MPI_TYPE_INDEXED( Num, BlockSizesA(Off+1),DisplacementsA(Off+1), &
     &                         CPP_MPI_REAL, Ptr, Ierror )
        Pattern%RecvDesc( ipe ) = Ptr
#endif
        Off = Off + Num
      ENDDO

!
! Now communicate what the receiver is expecting to the sender
!
      CALL ParExchangeVectorInt( InComm, Count, LocalA,                 &
     &                           CountOut, DisplacementsB  )
      CALL ParExchangeVectorInt( InComm, Count, BlockSizesA,            &
     &                           CountOut, BlockSizesB )

!
! Sender A: BlockSizes and Displacements can now be stored
!
      Off = 0
      DO ipe=1, GroupSize
        Num = CountOut(ipe)
        print *, "Sender Iam", Iam, "Ipe", Ipe-1, "Num", Num,           &
                 "Displacements", DisplacementsB(Off+1:Off+Num),        &
                 "BlockSizes", BlockSizesB(Off+1:Off+Num)
#if defined( USE_ARENAS )
        ALLOCATE( Pattern%SendDesc(ipe)%Displacements(Num) )
        ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(Num) )
        DO i=1, Num
          Pattern%SendDesc(ipe)%Displacements(i) = DisplacementsB(i+Off)
          Pattern%SendDesc(ipe)%BlockSizes(i)    = BlockSizesB(i+Off)
        ENDDO
#else
        CALL MPI_TYPE_INDEXED( Num, BlockSizesB(Off+1), DisplacementsB(Off+1),&
                               CPP_MPI_REAL, Ptr, Ierror )
        Pattern%SendDesc( ipe ) = Ptr
#endif
      ENDDO

      DEALLOCATE( CountOut )
      DEALLOCATE( Count )

      DEALLOCATE( PeB )
      DEALLOCATE( LocalB )
      DEALLOCATE( BlockSizesB )
      DEALLOCATE( DisplacementsB )

      DEALLOCATE( LocalA )
      DEALLOCATE( BlockSizesA )
      DEALLOCATE( DisplacementsA )

      CPP_LEAVE_PROCEDURE( "PARPATTERNDECOMPTOGHOST" )
      RETURN
!EOC
      END SUBROUTINE ParPatternDecompToGhost
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternFree --- Free the communication pattern
!
! !INTERFACE:
      SUBROUTINE ParPatternFree( InComm, Pattern )
!
! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )                 :: InComm  ! # of PEs
! !INPUT/OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( INOUT )  :: Pattern ! Comm Pattern
!
! !DESCRIPTION:
!     This routine frees a communication pattern.  
!
! !SYSTEM ROUTINES:
!     MPI_TYPE_FREE
!
! !BUGS:
!     The MPI_TYPE_FREE statement does not seem to work with FFC
!
! !REVISION HISTORY:
!   01.02.10   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  ipe, GroupSize, Pointer, Ierror

      CPP_ENTER_PROCEDURE( "PARPATTERNFREE" )

!
! First request the needed ghost values from other processors.
!
#if defined( USE_ARENAS )
      DO ipe, Pattern%Size
        DEALLOCATE( Pattern%RecvDesc(ipe)%Displacements )
        DEALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes )
        DEALLOCATE( Pattern%SendDesc(ipe)%Displacements )
        DEALLOCATE( Pattern%SendDesc(ipe)%BlockSizes )
      ENDDO
      DEALLOCATE( Pattern%RecvDesc )
      DEALLOCATE( Pattern%SendDesc )
#else
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
!
! Free all the MPI derived types
!
      DO ipe=1, Pattern%Size
        Pointer = Pattern%SendDesc(ipe)
        CALL MPI_TYPE_FREE( Pointer, Ierror )
        Pointer = Pattern%RecvDesc(ipe)
        CALL MPI_TYPE_FREE( Pointer, Ierror )
      ENDDO
#endif

      DEALLOCATE( Pattern%SendDesc )
      DEALLOCATE( Pattern%RecvDesc )

      CPP_LEAVE_PROCEDURE( "PARPATTERNFREE" )
      RETURN
!EOC
      END SUBROUTINE ParPatternFree
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParScatterReal --- Scatter slice to all PEs
!
! !INTERFACE:
      SUBROUTINE ParScatterReal ( InComm, Root, Slice, Decomp, Local )

! !USES:
      USE decompmodule, ONLY:  DecompType, Lists
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )          :: InComm       ! Communicator
      INTEGER, INTENT( IN )          :: Root         ! Root PE
      REAL(CPP_REAL), INTENT( IN )   :: Slice(*)     ! Global Slice
      TYPE(DecompType), INTENT( IN ) :: Decomp       ! Decomp information

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL), INTENT( OUT )  :: Local(*)     ! Local Slice

! !DESCRIPTION:
!     Given a decomposition of the domain, dole out a slice 
!     (one-dimensional array) to all the constituent PEs as described
!     by the decomposition Decomp.
!
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.04.14   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.05.01   Sawyer     Use Decomp%Comm for all local info
!   97.05.18   Sawyer     DecompType has moved to ParUtilitiesTypes
!   97.05.29   Sawyer     Changed 2-D arrays to 1-D
!   97.07.03   Sawyer     Reformulated documentation
!   97.07.22   Sawyer     DecompType has moved to DecompModule
!   97.12.01   Sawyer     Changed MPI_SSEND to MPI_ISEND
!   97.12.05   Sawyer     Added InComm and Root as arguments
!   97.12.05   Sawyer     Added logic to support intercommunicators
!   98.01.24   Sawyer     Removed dependence on MPI derived types TESTED
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.03.30   Sawyer     Stats dimension corrected: Gsize*MPI_STATUS_SIZE
!   99.01.19   Sawyer     Dropped assumed-size arrays
!   00.07.07   Sawyer     Removed "1D" references
!   00.07.23   Sawyer     Implementation with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

      INTEGER Ierror, I, J, K, L, Iam, GroupSize, Reqs( Gsize )
#if !defined( USE_ARENAS )
      INTEGER Status( MPI_STATUS_SIZE ), Stats( Gsize*MPI_STATUS_SIZE )
      REAL(CPP_REAL), ALLOCATABLE    :: SendBuf(:)
#endif
!
      CPP_ENTER_PROCEDURE( "PARSCATTERREAL" )
!
#if defined( USE_ARENAS )
!
! Pull the local process information out of the communicator
!
!!!      Iam = MOD( Comm, MAX_PES )
!!!      L   = Comm / MAX_PES
!!!      GroupSize = MOD( L, MAX_PES ) + 1
!!!      L   = L / MAX_PES
!!!      Color = MOD( L, MAX_PES )

!
! For now, Iam and GroupSize take on the global communicator values
!
      Iam = GID
      GroupSize = Gsize
      IF ( Iam .EQ. Root ) THEN
        L = 0
        DO I = 1, GroupSize
!
! Pick out the array sections to be sent.
! This is the inverse of the operation in ParGather
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L+1
              DataBuf(L) = Slice(K)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
! Barrier: all PEs participate!
!
      CALL mlp_barrier(gid,gsize)

!
! All receive from the root.  
!

! The local array may be larger than that specified in the decomposition
!
      L = 0
      IF ( Iam .GT. 0 ) L = SUM( Decomp%NumEntries(1:Iam) )
      DO I=1, Decomp%NumEntries(Iam+1)
        Local( I ) = DataBuf( L + I )
      ENDDO
!
! The following is needed to ensure that DataBuf can now be reused
!
      CALL mlp_barrier(gid,gsize)
#else
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )

      IF ( Iam .EQ. Root ) THEN
        ALLOCATE( SendBuf( SUM( Decomp%NumEntries ) ) )
        L = 0
        DO I = 1, GroupSize
!
! Pick out the array sections to be sent.
! This is the inverse of the operation in ParGather
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L+1
              SendBuf(L) = Slice(K)
            ENDDO
          ENDDO
!
! This is a non-blocking send. SendBuf cannot be immediately deallocated
!
! WARNING: F90-MPI inconsistency: make sure the indexing below always works
!
          CALL MPI_ISEND( SendBuf(L-Decomp%NumEntries(I)+1),             &
     &                    Decomp%NumEntries(I), CPP_MPI_REAL,            &
     &                    I-1, 0, InComm, Reqs(I), Ierror )

        ENDDO
      ENDIF

!
! All receive from the root.  
!
! The local array may be larger than that specified in the decomposition
!
      CALL MPI_RECV( Local, Decomp%NumEntries(Iam+1),                    &
     &               CPP_MPI_REAL,                                       &
     &               Root, 0, InComm, Status, Ierror )
!
! Experience shows that we should wait for all the non-blocking
! PEs to check in, EVEN THOUGH THE MPI_RECV HAS COMPLETED !!
!
      IF ( Iam .EQ. Root ) THEN
        CALL MPI_WAITALL( GroupSize, Reqs, Stats, Ierror )
        DEALLOCATE( SendBuf )
      ENDIF

!
! The following may be needed on some platforms to avoid an MPI bug.
!
      CALL MPI_BARRIER( InComm, Ierror )
#endif
      CPP_LEAVE_PROCEDURE( "PARSCATTERREAL" )
      RETURN
!EOC
      END SUBROUTINE ParScatterReal
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParScatterInt --- Scatter slice to all PEs
!
! !INTERFACE:
      SUBROUTINE ParScatterInt ( InComm, Root, Slice, Decomp, Local )

! !USES:
      USE decompmodule, ONLY:  DecompType, Lists
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )          :: InComm       ! Communicator
      INTEGER, INTENT( IN )          :: Root         ! Root PE
      INTEGER, INTENT( IN )          :: Slice(*)     ! Global Slice
      TYPE(DecompType), INTENT( IN ) :: Decomp       ! Decomp information

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )         :: Local(*)     ! Local Slice

! !DESCRIPTION:
!     Given a decomposition of the domain, dole out a slice 
!     (one-dimensional array) to all the constituent PEs as described
!     by the decomposition Decomp.
!
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.04.14   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.05.01   Sawyer     Use Decomp%Comm for all local info
!   97.05.18   Sawyer     DecompType has moved to ParUtilitiesTypes
!   97.05.29   Sawyer     Changed 2-D arrays to 1-D
!   97.07.03   Sawyer     Reformulated documentation
!   97.07.22   Sawyer     DecompType has moved to DecompModule
!   97.12.01   Sawyer     Changed MPI_SSEND to MPI_ISEND
!   97.12.05   Sawyer     Added InComm and Root as arguments
!   97.12.05   Sawyer     Added logic to support intercommunicators
!   98.01.24   Sawyer     Removed dependence on MPI derived types TESTED
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.03.30   Sawyer     Stats dimension corrected: Gsize*MPI_STATUS_SIZE
!   99.01.19   Sawyer     Dropped assumed-size arrays
!   00.07.07   Sawyer     Removed "1D" references
!   00.07.23   Sawyer     Implementation with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

      INTEGER Ierror, I, J, K, L, Iam, GroupSize, Reqs( Gsize )
#if !defined( USE_ARENAS )
      INTEGER Status( MPI_STATUS_SIZE ), Stats( Gsize*MPI_STATUS_SIZE )
      INTEGER, ALLOCATABLE    :: SendBuf(:)
#endif
!
      CPP_ENTER_PROCEDURE( "PARSCATTERINT" )
!
#if defined( USE_ARENAS )
!
! Pull the local process information out of the communicator
!
!!!      Iam = MOD( Comm, MAX_PES )
!!!      L   = Comm / MAX_PES
!!!      GroupSize = MOD( L, MAX_PES ) + 1
!!!      L   = L / MAX_PES
!!!      Color = MOD( L, MAX_PES )

!
! For now, Iam and GroupSize take on the global communicator values
!
      Iam = GID
      GroupSize = Gsize
      IF ( Iam .EQ. Root ) THEN
        L = 0
        DO I = 1, GroupSize
!
! Pick out the array sections to be sent.
! This is the inverse of the operation in ParGather
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L+1
              IntBuf(L) = Slice(K)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
! Barrier: all PEs participate!
!
      CALL mlp_barrier(gid,gsize)

!
! All receive from the root.  
!

! The local array may be larger than that specified in the decomposition
!
      L = 0
      IF ( Iam .GT. 0 ) L = SUM( Decomp%NumEntries(1:Iam) )
      DO I=1, Decomp%NumEntries(Iam+1)
        Local( I ) = IntBuf( L + I )
      ENDDO
!
! The following is needed to ensure that IntBuf can now be reused
!
      CALL mlp_barrier(gid,gsize)
#else
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )

      IF ( Iam .EQ. Root ) THEN
        ALLOCATE( SendBuf( SUM( Decomp%NumEntries ) ) )
        L = 0
        DO I = 1, GroupSize
!
! Pick out the array sections to be sent.
! This is the inverse of the operation in ParGather
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L+1
              SendBuf(L) = Slice(K)
            ENDDO
          ENDDO
!
! This is a non-blocking send. SendBuf cannot be immediately deallocated
!
! WARNING: F90-MPI inconsistency: make sure the indexing below always works
!
          CALL MPI_ISEND( SendBuf(L-Decomp%NumEntries(I)+1),                  &
     &                    Decomp%NumEntries(I), CPP_MPI_INTEGER,              &
     &                    I-1, 0, InComm, Reqs(I), Ierror )

        ENDDO
      ENDIF

!
! All receive from the root.  
!
! The local array may be larger than that specified in the decomposition
!
      CALL MPI_RECV( Local, Decomp%NumEntries(Iam+1),                         &
     &               CPP_MPI_INTEGER,                                         &
     &               Root, 0, InComm, Status, Ierror )
!
! Experience shows that we should wait for all the non-blocking
! PEs to check in, EVEN THOUGH THE MPI_RECV HAS COMPLETED !!
!
      IF ( Iam .EQ. Root ) THEN
        CALL MPI_WAITALL( GroupSize, Reqs, Stats, Ierror )
        DEALLOCATE( SendBuf )
      ENDIF

!
! The following may be needed on some platforms to avoid an MPI bug.
!
      CALL MPI_BARRIER( InComm, Ierror )
#endif
      CPP_LEAVE_PROCEDURE( "PARSCATTERINT" )
      RETURN
!EOC
      END SUBROUTINE ParScatterInt
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParGatherReal --- Gather Slice from all PEs
!
! !INTERFACE:  
      SUBROUTINE ParGatherReal ( InComm, Root, Local, Decomp, Slice )

! !USES:
      USE decompmodule, ONLY:  DecompType, Lists
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )          :: InComm       ! Communicator
      INTEGER, INTENT( IN )          :: Root         ! Root PE
      REAL(CPP_REAL), INTENT( IN )   :: Local(*)     ! Local Slice
      TYPE(DecompType), INTENT( IN ) :: Decomp       ! Decomp information

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL), INTENT( OUT )  :: Slice(*)     ! Global Slice

! !DESCRIPTION:
!     Given a decomposition of the domain and a local portion of the
!     total slice on each PE, gather together the portions into a
!     global slice on the root PE
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.04.14   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.05.01   Sawyer     Use Decomp%Comm for all local info
!   97.05.18   Sawyer     DecompType has moved to ParUtilitiesTypes
!   97.05.29   Sawyer     Changed 2-D arrays to 1-D
!   97.07.03   Sawyer     Reformulated documentation
!   97.07.22   Sawyer     DecompType has moved to DecompModule
!   97.12.01   Sawyer     Changed MPI_SSEND to MPI_ISEND
!   97.12.05   Sawyer     Added InComm and Root as arguments
!   97.12.05   Sawyer     Added logic to support intercommunicators
!   98.01.24   Sawyer     Removed dependence on MPI derived types TESTED
!   98.01.29   Sawyer     Corrected assertions
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.03.31   Sawyer     Stat dimension corrected: MPI_STATUS_SIZE
!   98.04.22   Sawyer     Local no longer assumed shape: Local(*)
!   99.01.19   Sawyer     Dropped assumed-size arrays
!   00.07.07   Sawyer     Removed "1D" references
!   00.07.23   Sawyer     Implementation with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER Ierror, I, J, K, L, Iam, GroupSize, Req
#if !defined( USE_ARENAS )
      INTEGER Status( MPI_STATUS_SIZE ), Stat( MPI_STATUS_SIZE )
#endif
      REAL(CPP_REAL), ALLOCATABLE    :: RecvBuf(:)
!
      CPP_ENTER_PROCEDURE( "PARGATHERREAL" )
!
#if defined( USE_ARENAS )
!
! Pull the local process information out of the communicator
!
!!!      Iam = MOD( Comm, MAX_PES )
!!!      L   = Comm / MAX_PES
!!!      GroupSize = MOD( L, MAX_PES ) + 1
!!!      L   = L / MAX_PES
!!!      Color = MOD( L, MAX_PES )

!
! For now, Iam and GroupSize take on the global communicator values
!
      Iam = GID
      GroupSize = Gsize
!
! All PEs send their contribution to the root
!
      L = 0
      IF ( Iam .GT. 0 ) L = SUM( Decomp%NumEntries(1:Iam) )
      DO I=1, Decomp%NumEntries(Iam+1)
        DataBuf( L + I ) = Local( I )
      ENDDO
      CALL mlp_barrier(gid,gsize)
!
      IF ( Iam .EQ. Root ) THEN
!
! On the Root PE receive from every other PE
!
        L = 0
        DO I = 1, GroupSize
!
! This is the simple reverse mapping of that in ParScatter
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L + 1
              Slice(K) = DataBuf(L)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
! The following is needed to ensure that DataBuf can now be reused
!
      CALL mlp_barrier(gid,gsize)
#else
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
!
! All PEs send their contribution to the root
!
      CALL MPI_ISEND( Local, Decomp%NumEntries(Iam+1),                   &
     &                CPP_MPI_REAL,                                      &
     &                Root, Iam+3001, InComm, Req, Ierror )

      IF ( Iam .EQ. Root ) THEN
        ALLOCATE( RecvBuf( SUM( Decomp%NumEntries ) ) )
!
! On the Root PE receive from every other PE
!
        L = 0
        DO I = 1, GroupSize
!
! This is a blocking, synchronous recv.  All the
! sends should have been posted so it should not deadlock
!
! WARNING: F90-MPI inconsistency: make sure the indexing below always works
!
          CPP_ASSERT_F90( L .LT. SIZE( RecvBuf ) )
          CALL MPI_RECV( RecvBuf(L+1), Decomp%NumEntries(I),             &
     &                   CPP_MPI_REAL, I-1, I+3000, InComm,              &
     &                   Status, Ierror )
!
! This is the simple reverse mapping of that in ParScatter
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L + 1
              Slice(K) = RecvBuf(L)
#if defined(DEBUG_PARGATHERREAL)
                PRINT *, " Entry ", L, RecvBuf(L), K, SIZE(Slice)
#endif
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE( RecvBuf )
      ENDIF
      CALL MPI_WAIT( Req, Stat, Ierror )
!
! The following may be needed on some platforms to avoid an MPI bug.
!
      CALL MPI_BARRIER( InComm, Ierror )
#endif
      CPP_LEAVE_PROCEDURE( "PARGATHERREAL" )
      RETURN
!EOC
      END SUBROUTINE ParGatherReal
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParGatherInt --- Gather Slice from all PEs
!
! !INTERFACE:  
      SUBROUTINE ParGatherInt ( InComm, Root, Local, Decomp, Slice )

! !USES:
      USE decompmodule, ONLY:  DecompType, Lists
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )          :: InComm       ! Communicator
      INTEGER, INTENT( IN )          :: Root         ! Root PE
      INTEGER, INTENT( IN )          :: Local(*)     ! Local Slice
      TYPE(DecompType), INTENT( IN ) :: Decomp       ! Decomp information

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )         :: Slice(*)     ! Global Slice

! !DESCRIPTION:
!     Given a decomposition of the domain and a local portion of the
!     total slice on each PE, gather together the portions into a
!     global slice on the root PE
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.04.14   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.05.01   Sawyer     Use Decomp%Comm for all local info
!   97.05.18   Sawyer     DecompType has moved to ParUtilitiesTypes
!   97.05.29   Sawyer     Changed 2-D arrays to 1-D
!   97.07.03   Sawyer     Reformulated documentation
!   97.07.22   Sawyer     DecompType has moved to DecompModule
!   97.12.01   Sawyer     Changed MPI_SSEND to MPI_ISEND
!   97.12.05   Sawyer     Added InComm and Root as arguments
!   97.12.05   Sawyer     Added logic to support intercommunicators
!   98.01.24   Sawyer     Removed dependence on MPI derived types TESTED
!   98.01.29   Sawyer     Corrected assertions
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.03.31   Sawyer     Stat dimension corrected: MPI_STATUS_SIZE
!   98.04.22   Sawyer     Local no longer assumed shape: Local(*)
!   99.01.19   Sawyer     Dropped assumed-size arrays
!   00.07.07   Sawyer     Removed "1D" references
!   00.07.23   Sawyer     Implementation with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER Ierror, I, J, K, L, Iam, GroupSize, Req
#if !defined( USE_ARENAS )
      INTEGER Status( MPI_STATUS_SIZE ), Stat( MPI_STATUS_SIZE )
#endif
      INTEGER, ALLOCATABLE    :: RecvBuf(:)
!
      CPP_ENTER_PROCEDURE( "PARGATHERINT" )
!
#if defined( USE_ARENAS )
!
! Pull the local process information out of the communicator
!
!!!      Iam = MOD( Comm, MAX_PES )
!!!      L   = Comm / MAX_PES
!!!      GroupSize = MOD( L, MAX_PES ) + 1
!!!      L   = L / MAX_PES
!!!      Color = MOD( L, MAX_PES )

!
! For now, Iam and GroupSize take on the global communicator values
!
      Iam = GID
      GroupSize = Gsize
!
! All PEs send their contribution to the root
!
      L = 0
      IF ( Iam .GT. 0 ) L = SUM( Decomp%NumEntries(1:Iam) )
      DO I=1, Decomp%NumEntries(Iam+1)
        IntBuf( L + I ) = Local( I )
      ENDDO
      CALL mlp_barrier(gid,gsize)
!
      IF ( Iam .EQ. Root ) THEN
!
! On the Root PE receive from every other PE
!
        L = 0
        DO I = 1, GroupSize
!
! This is the simple reverse mapping of that in ParScatter
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L + 1
              Slice(K) = IntBuf(L)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
! The following is needed to ensure that IntBuf can now be reused
!
      CALL mlp_barrier(gid,gsize)
#else
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
!
! All PEs send their contribution to the root
!
      CALL MPI_ISEND( Local, Decomp%NumEntries(Iam+1), CPP_MPI_INTEGER,       &
     &                Root, Iam+3001, InComm, Req, Ierror )

      IF ( Iam .EQ. Root ) THEN
        ALLOCATE( RecvBuf( SUM( Decomp%NumEntries ) ) )
!
! On the Root PE receive from every other PE
!
        L = 0
        DO I = 1, GroupSize
!
! This is a blocking, synchronous recv.  All the
! sends should have been posted so it should not deadlock
!
! WARNING: F90-MPI inconsistency: make sure the indexing below always works
!
          CPP_ASSERT_F90( L .LT. SIZE( RecvBuf ) )
          CALL MPI_RECV( RecvBuf(L+1), Decomp%NumEntries(I),                  &
     &                   CPP_MPI_INTEGER, I-1, I+3000, InComm,                &
     &                   Status, Ierror )
!
! This is the simple reverse mapping of that in ParScatter
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L + 1
              Slice(K) = RecvBuf(L)
#if defined(DEBUG_PARGATHERINT)
                PRINT *, " Entry ", L, RecvBuf(L), K, SIZE(Slice)
#endif
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE( RecvBuf )
      ENDIF
      CALL MPI_WAIT( Req, Stat, Ierror )
!
! The following may be needed on some platforms to avoid an MPI bug.
!
      CALL MPI_BARRIER( InComm, Ierror )
#endif
      CPP_LEAVE_PROCEDURE( "PARGATHERINT" )
      RETURN
!EOC
      END SUBROUTINE ParGatherInt
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParBeginTransferReal --- Start an ASYNC Real Transfer
!
! !INTERFACE:
      SUBROUTINE ParBeginTransferReal(InComm, NrInPackets, NrOutPackets, &
     &                                Dest, Src, InBuf, InIA,            &
     &                                OutBuf, OutIA )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      INTEGER, INTENT( IN )       :: NrInPackets  ! Number of in packets
      INTEGER, INTENT( IN )       :: NrOutPackets ! Number of out packets
      INTEGER, INTENT( IN )       :: Dest(:)      ! PE destinations
      INTEGER, INTENT( IN )       :: Src(:)       ! PE sources
      REAL(CPP_REAL), INTENT( IN ):: InBuf(:)     ! Input buffer
      INTEGER, INTENT( IN )       :: InIA(:)      ! In packet counter
      INTEGER, INTENT( IN )       :: OutIA(:)     ! Out packet counter

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL), INTENT( OUT ) :: OutBuf(:)  ! Output buffer

! !DESCRIPTION: 
!
!     This routine initiates an async. transfer of an array InBuf
!     partitioned into chunks defined by the arrays InIA and Dest
!     to an output array OutBuf on another PE. InIA(1) contains 
!     the number of reals to be sent to Dest(1), InIA(2) the number 
!     of reals to be sent to Dest(2), etc.  Similarly, the array
!     OutBuf on the calling PE is partitioned into chunks by OutIA
!     and Src, with OutIA(1) the number of reals anticipated from
!     Src(1), etc.  
!
!     The default implementation reads through the contiguous array 
!     InBuf and sends the chunks to the PEs designated with an 
!     asyncronous MPI\_ISEND.  Correspondingly it posts the receives 
!     with an asynchronous MPI\_IRECV.   The USE\_ARENAS implementation
!     in fully functional, but only for a global communicator.
!
!     Wait handles InHandle(:) and OutHandle(:) are in common block.
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference!!!!!!!!!!
!
!     The buffers may not be accessed until after the call to 
!     ParEndTransferReal.
!
!
! !SYSTEM ROUTINES:
!     MPI_COMM_RANK, MPI_ISEND, MPI_IRECV
!
! !REVISION HISTORY:
!   97.09.26   Sawyer     Creation
!   97.12.05   Sawyer     Renamed Comm to InComm to avoid collisions
!   98.02.26   Sawyer     Added Dest, Src and Remote to clean up code
!   98.04.16   Sawyer     Number of packets become input arguments
!   98.09.04   Sawyer     Cleaned interface: handles in common, no Remote
!   99.03.04   Sawyer     Inlined ParCalculateRemote
!   99.06.01   Sawyer     Changed pointer arrays to INTEGER*8 for SGI
!   00.08.07   Sawyer     Implementation with shared memory arenas
! 
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
      INTEGER Iam, GroupSize, Nr, Icnt, Packet, I, Ierr
#if defined( USE_ARENAS )
      INTEGER Off( 0:Gsize ), inc, base
#endif

      CPP_ENTER_PROCEDURE( "PARBEGINTRANSFERREAL" )
      CPP_ASSERT_F90( NrInPackets .LE. SIZE( Dest ) )
      CPP_ASSERT_F90( NrInPackets .LE. SIZE( InIA ) )
      CPP_ASSERT_F90( NrOutPackets .LE. SIZE( Src ) )
      CPP_ASSERT_F90( NrOutPackets .LE. SIZE( OutIA ) )

#if defined( USE_ARENAS )
!
! The following code has the effect of sorting the packets
! by their destination and filling the DataBuf
!
      Base = GID*BlockSize   ! Each PE is allotted Blocksize
      Off = 0                ! Vector of offsets
!
! First count how much is going to each PE
!
      DO Packet=1, NrInPackets
        Off(Dest(Packet)+1) = Off(Dest(Packet)+1)+InIA(Packet)
      ENDDO
!
! Calculate the offsets, and save the counts in Volume
!
      DO I=1, Gsize
        Off(I) = Off(I)+Off(I-1)   ! Beware: "Off(0)" is 0
        Volume(I,GID+1) = Off(I)   ! 00.09.09 changed index order
      ENDDO
      CPP_ASSERT_F90( Off(Gsize) .LE. BlockSize )
      inc = 0
      DO Packet=1, NrInPackets
        Icnt = Base + Off(Dest(Packet))
        DO I=1,InIA(Packet)
          inc = inc+1 
          DataBuf(Icnt+I) = InBuf( inc )
        ENDDO
        Off(Dest(Packet)) = Off(Dest(Packet)) + InIA(Packet)
      END DO
#else
      CALL MPI_COMM_RANK( InComm, Iam, Ierr )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierr )

!
! Increment the ongoing transfer number
      BegTrf = MOD(BegTrf,MAX_TRF) + 1
!
!     MPI: Irecv over all processes
!
      Icnt = 1
      DO Packet = 1, NrOutPackets
        Nr = OutIA( Packet )
        IF ( Nr .GT. 0 ) THEN
#if defined( DEBUG_PARBEGINTRANSFERREAL )
          PRINT *, "Iam ",Iam," posts recv ",Nr," from ", Src( Packet )
#endif
!
! Receive the buffers with MPI_Irecv. Non-blocking
!
          CPP_ASSERT_F90( Icnt+Nr-1 .LE. SIZE( OutBuf ) )
          CALL MPI_IRECV( OutBuf( Icnt ), Nr,                            &
     &          CPP_MPI_REAL, Src( Packet ), Src( Packet ),              &
     &          InComm, OutHandle(Packet,BegTrf), Ierr )
        ELSE
          OutHandle( Packet,BegTrf ) = MPI_REQUEST_NULL
        END IF
        Icnt = Icnt + Nr
      END DO
!
!     MPI: Isend over all processes
!
      Icnt = 1
      CPP_ASSERT_F90( NrInPackets .LE. SIZE( Dest ) )
      CPP_ASSERT_F90( NrInPackets .LE. SIZE( InIA ) )
      DO Packet = 1, NrInPackets
        Nr = InIA( Packet )
        IF ( Nr .GT. 0 ) THEN
#if defined( DEBUG_PARBEGINTRANSFERREAL )
          PRINT *,"Iam ",Iam," posts send ",Nr," to ",Dest( Packet )
#endif
!
!     Send the individual buffers with non-blocking sends
!
          CPP_ASSERT_F90( Icnt+Nr-1 .LE. SIZE( InBuf ) )
          CALL MPI_ISEND ( InBuf( Icnt ), Nr,                            &
     &          CPP_MPI_REAL, Dest( Packet ), Iam,                       &
     &          InComm, InHandle( Packet,BegTrf ), Ierr )
        ELSE
          InHandle( Packet,BegTrf ) = MPI_REQUEST_NULL
        END IF
        Icnt = Icnt + Nr
      END DO
#endif
!
!
      CPP_LEAVE_PROCEDURE( "PARBEGINTRANSFERREAL" )
      RETURN
!EOC
      END SUBROUTINE ParBeginTransferReal
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParBeginTransferPattern --- Start an ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParBeginTransferPattern( Pattern, InBuf, OutBuf )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      TYPE (ParPatternType), INTENT( IN )  :: Pattern   ! Comm Pattern
      REAL(CPP_REAL), INTENT( IN )         :: InBuf(*)  ! Input buffer

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL), INTENT( OUT )        :: OutBuf(*) ! Output buffer

! !DESCRIPTION: 
!
!     This routine initiates an async. transfer of an array InBuf.
!     The communication pattern indicates the indices outgoing 
!     values of InBuf and  incoming values for OutBuf.  This routine
!     is fundamentally equivalent to ParBeginTransferReal; the use 
!     of a communication pattern is largely a performance enhancement, 
!     since it eliminates the need for intermediate buffering.
!     
!     Wait handles InHandle(:) and OutHandle(:) are in common block.
!     The buffers may not be accessed until after the call to 
!     ParEndTransferReal.  
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference.  The USE\_ARENAS version is not yet functional.
!
! !SYSTEM ROUTINES:
!     MPI_COMM_RANK, MPI_ISEND, MPI_IRECV
!
! !REVISION HISTORY:
!   01.02.14   Sawyer     Creation from ParBeginTransferReal
! 
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
      INTEGER TypeDesc, Nr, Icnt, Packet, I, ipe, Ierr
#if defined( USE_ARENAS )
      INTEGER Off( 0:Gsize ), inc, base
#endif

      CPP_ENTER_PROCEDURE( "PARBEGINTRANSFERPATTERN" )

#if defined( USE_ARENAS )
      The USE_ARENAS code is not running yet.  You should not
      be compiling with this option.

!
! The following code has the effect of sorting the packets
! by their destination and filling the DataBuf
!
      Base = GID*BlockSize   ! Each PE is allotted Blocksize
      Off = 0                ! Vector of offsets
!
! First count how much is going to each PE
!
      DO Packet=1, NrInPackets
        Off(Dest(Packet)+1) = Off(Dest(Packet)+1)+InIA(Packet)
      ENDDO
!
! Calculate the offsets, and save the counts in Volume
!
      DO I=1, Gsize
        Off(I) = Off(I)+Off(I-1)   ! Beware: "Off(0)" is 0
        Volume(I,GID+1) = Off(I)   ! 00.09.09 changed index order
      ENDDO
      CPP_ASSERT_F90( Off(Gsize) .LE. BlockSize )
      inc = 0
      DO Packet=1, NrInPackets
        Icnt = Base + Off(Dest(Packet))
        DO I=1,InIA(Packet)
          inc = inc+1 
          DataBuf(Icnt+I) = InBuf( inc )
        ENDDO
        Off(Dest(Packet)) = Off(Dest(Packet)) + InIA(Packet)
      END DO
#else

      BegTrf = MOD(BegTrf,MAX_TRF) + 1
!
! MPI: Irecv over all processes
!

      DO ipe = 1, Pattern%Size
!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        TypeDesc = Pattern%RecvDesc(ipe)
        CALL MPI_IRECV( OutBuf, 1, TypeDesc, ipe-1, ipe-1,               &
     &                  Pattern%Comm, OutHandle(ipe,BegTrf), Ierr )
      END DO
!
! MPI: Isend over all processes
!
      DO ipe = 1, Pattern%Size
!
! Send the individual buffers with non-blocking sends
!
        TypeDesc = Pattern%SendDesc(ipe)
        CALL MPI_ISEND ( InBuf, 1, TypeDesc, ipe-1, Pattern%Iam,         &
     &                   Pattern%Comm, InHandle(ipe,BegTrf), Ierr )
      END DO
#endif
!
!
      CPP_LEAVE_PROCEDURE( "PARBEGINTRANSFERPATTERN" )
      RETURN
!EOC
      END SUBROUTINE ParBeginTransferPattern
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParEndTransferReal --- Complete an ASYNC Real Transfer
!
! !INTERFACE:
      SUBROUTINE ParEndTransferReal( InComm, NrInPackets, NrOutPackets,  &
     &                               Dest, Src, InBuf, InIA,             &
     &                               OutBuf, OutIA )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      INTEGER, INTENT( IN )       :: NrInPackets  ! Number of in packets
      INTEGER, INTENT( IN )       :: NrOutPackets ! Number of out packets
      INTEGER, INTENT( IN )       :: Dest(:)      ! PE destinations
      INTEGER, INTENT( IN )       :: Src(:)       ! PE sources
      REAL(CPP_REAL), INTENT( IN ):: InBuf(:)     ! Input buffer
      INTEGER, INTENT( IN )       :: InIA(:)      ! Pointer array
      INTEGER, INTENT( IN )       :: OutIA(:)     ! Pointer array

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL), INTENT( INOUT ) :: OutBuf(:)! Output buffer

! !DESCRIPTION: 
!
!     This routine completes an async. transfer of an array
!     partitioned into chunks defined by the array InIA.  In the 
!     MPI version, neither InBuf nor OutBuf is not used since
!     that information was utilized in ParBeginTransferReal.  In
!     the USE\_ARENAS version they are utilized.
!
!     The link between StartTransfer and EndTransfer is made possible
!     by the InHandle and OutHandle: they reflect the status of
!     the ongoing transfer.  When this routine completes, a valid
!     and accessible copy of the OutBuf is ready for use. The 
!     USE\_ARENAS version is functional, but assumes a global communicator.
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference! The buffers may not be accessed until after the 
!     completion of ParEndTransferReal.  
!
!
! !SYSTEM ROUTINES:
!     MPI_COMM_RANK, MPI_ISEND, MPI_IRECV
!
! !REVISION HISTORY:
!   97.09.26   Sawyer     Creation
!   97.12.05   Sawyer     Renamed Comm to InComm to avoid collisions
!   98.02.26   Sawyer     Count through packets, not PEs
!   98.04.16   Sawyer     Number of packets become input arguments
!   98.09.04   Sawyer     Cleaned interface: handles in common
!   99.03.05   Sawyer     Support for contiguous communicators in SHMEM
!   99.04.22   Sawyer     Bug fix: replaced MPI_WAIT with MPI_WAITALL
!   99.06.03   Sawyer     Bug fix: GroupSize in SHMEM_BARRIER
!   00.07.28   Sawyer     Implemented with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER Iam, GroupSize, J, Offset, Packet, Ierr
#if defined( USE_ARENAS)
      INTEGER Off( 0:Gsize ), Pe, PacketSize, inc
#else  
      INTEGER InStats(NrInPackets*MPI_STATUS_SIZE)
      INTEGER OutStats(NrOutPackets*MPI_STATUS_SIZE)
#endif

      CPP_ENTER_PROCEDURE( "PARENDTRANSFERREAL" )
#if defined( USE_ARENAS )
!
! Needed to ensure that DataBuf and Volume have "arrived"
!
      CALL mlp_barrier(gid,gsize)
!
! The packet is plucked out of the appropriate data segment Blocksize*(Pe-1)
! with the offset specified in Volume(Pe,GID+1).  Since there can be
! more than one packet from Pe bound for me (GID), Off(Pe) is advanced
! by the packetsize.  This should yield the Offset to the needed packet.
!
      Off = 0
      inc = 0
      DO Packet=1, NrOutPackets
        Pe = Src(Packet)
        Offset = Off(Pe) + Blocksize*Pe
! 00.09.09 changed index order of Volume
        IF ( GID .GT. 0 ) Offset = Offset + Volume(GID,Pe+1)
        PacketSize = OutIA(Packet)
        DO J=1, PacketSize
          inc = inc + 1
          OutBuf(inc) = DataBuf(Offset+J)
        ENDDO
        Off(Pe) = Off(Pe) + PacketSize
      ENDDO
!
! Needed to ensure that a future call cannot upset current operation
!
      CALL mlp_barrier(gid,gsize)
#else
      EndTrf = MOD(EndTrf,MAX_TRF)+1
      CPP_ASSERT_F90( NrInPackets .LE. MAX_PAX )
      CALL MPI_WAITALL( NrInPackets, InHandle(:,EndTrf), InStats, Ierr )
 
      CPP_ASSERT_F90( NrOutPackets .LE. MAX_PAX )
      CALL MPI_WAITALL( NrOutPackets, OutHandle(:,EndTrf), OutStats, Ierr )
!
! WS 98.09.22 : This barrier needed to synchronize.
!
      CALL MPI_BARRIER( InComm, Ierr )
#endif
      CPP_LEAVE_PROCEDURE( "PARENDTRANSFERREAL" )
      RETURN
!EOC
      END SUBROUTINE ParEndTransferReal
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParEndTransferPattern --- Complete an ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParEndTransferPattern( Pattern, InBuf, OutBuf )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      TYPE (ParPatternType), INTENT( IN )  :: Pattern   ! Comm Pattern
      REAL(CPP_REAL), INTENT( IN )         :: InBuf(*)  ! Input buffer

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL), INTENT( INOUT )      :: OutBuf(*) ! Output buffer

! !DESCRIPTION: 
!
!     This routine completes an async. transfer of an array communicated
!     with a communication pattern.  
!
!     The link between StartTransfer and EndTransfer is made possible
!     by the InHandle and OutHandle: they reflect the status of
!     the ongoing transfer.  When this routine completes, a valid
!     and accessible copy of the OutBuf is ready for use.
!     The buffers may not be accessed until after the 
!     completion of ParEndTransfer.  
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference! USE_ARENAS version not yet functional.
!
! !SYSTEM ROUTINES:
!     MPI_COMM_RANK, MPI_ISEND, MPI_IRECV
!
! !REVISION HISTORY:
!   01.02.14   Sawyer     Creation from ParEndTransferReal
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER J, Offset, Packet, Ierr
#if defined( USE_ARENAS)
      INTEGER Off( 0:Gsize ), Pe, PacketSize, inc
#else  
      INTEGER InStats(Pattern%Size*MPI_STATUS_SIZE)
      INTEGER OutStats(Pattern%Size*MPI_STATUS_SIZE)
#endif

      CPP_ENTER_PROCEDURE( "PARENDTRANSFERPATTERN" )
#if defined( USE_ARENAS )
!      The USE_ARENAS code does not work yet.  You should not
!      be compiling with this CPP token on.
#else
      EndTrf = MOD(EndTrf,MAX_TRF) + 1
      CALL MPI_WAITALL( Pattern%Size, InHandle(:,EndTrf), InStats, Ierr )
      CALL MPI_WAITALL( Pattern%Size, OutHandle(:,EndTrf), OutStats, Ierr )
!
! WS 98.09.22 : This barrier needed to synchronize.
!
      CALL MPI_BARRIER( Pattern%Comm, Ierr )
#endif
      CPP_LEAVE_PROCEDURE( "PARENDTRANSFERPATTERN" )
      RETURN
!EOC
      END SUBROUTINE ParEndTransferPattern
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParExchangeVectorReal --- Exchange a sparse packed vector
!
! !INTERFACE:  
      SUBROUTINE ParExchangeVectorReal ( InComm, LenInVector, InVector,  &
     &                                   LenOutVector, OutVector )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )   :: InComm            ! Communicator
      INTEGER, INTENT( IN )   :: LenInVector( * )  ! Length on each PE
      REAL(CPP_REAL), INTENT( IN ):: InVector( * ) ! The input buffer

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )  :: LenOutVector( * ) ! Length on each PE
      REAL(CPP_REAL), INTENT( OUT ) :: OutVector( * ) ! The output buffer

! !DESCRIPTION:
!
!     This routine exchanges vectors stored in compressed format, i.e.,
!     in so-called compressed sparse row (CSR) format, with other
!     PEs.  In essence it first exchanges the lengths with
!     MPI\_Alltoall, then the exchange of the actual vectors (can be
!     different in size) using MPI\_AlltoallV.  Since the latter is
!     inefficient, it is simulated using MPI\_Isend and MPI\_Recv.
!
!     The USE\_ARENAS version is fully functional, but assumes a 
!     global communicator.
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_WAITALL, MPI_ALLTOALL
!
! !REVISION HISTORY:
!   98.03.17   Sawyer     Creation from F77 version
!   98.03.30   Sawyer     Removed assumed shape arrays due to problems
!   99.01.18   Sawyer     Added barrier for safety
!   99.03.08   Sawyer     USE_SHMEM version for CRAY only; untested
!   99.06.01   Sawyer     USE_SHMEM version revised per comments from Tom
!   00.07.28   Sawyer     Implemented with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER :: i, iscnt, ircnt, nr, pe, icnt, Nsize, Iam, Ierr
      INTEGER :: Status(MPI_STATUS_SIZE)
      INTEGER :: Reqs( Gsize ), Stats( Gsize*MPI_STATUS_SIZE )
#if defined( USE_ARENAS )
      INTEGER :: offset( MAX_PE ), goffset( MAX_PE )
#endif

      CPP_ENTER_PROCEDURE( "PAREXCHANGEVECTORREAL" )
#if defined( USE_ARENAS )
      DO pe = 1, Gsize
        volume( pe, gid+1 ) = LenInVector( pe )
      ENDDO
      CALL mlp_barrier(gid,gsize)
      icnt = 0
      DO pe = 1, Gsize
        LenOutVector( pe ) = volume( gid+1, pe )
        offset( pe )  = icnt
        goffset( pe ) = icnt
        DO i = 1, gid
          goffset( pe ) = goffset( pe ) + volume( pe, i )
        ENDDO
        icnt = icnt + sum( volume( pe, 1:Gsize ) )
      ENDDO
      icnt = 0
      DO pe = 1, Gsize
        DO i = goffset(pe)+1, goffset(pe)+LenInVector(pe)
          icnt = icnt + 1
          databuf( i ) = Invector( icnt )
        ENDDO
      ENDDO
      CALL mlp_barrier(gid,gsize)
      DO i = 1, SUM(LenOutVector(1:Gsize))
        OutVector( i ) = databuf( offset(gid+1) + i )
      ENDDO
#else
      CALL MPI_COMM_SIZE( InComm, Nsize, Ierr )
      CALL MPI_COMM_RANK( InComm, Iam, Ierr )

#if defined( MY_ALLTOALL )
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        CALL MPI_ISEND( nr, 1,                                           &
     &                  MPI_INTEGER, pe, Iam+3000,                       &
     &                  InComm, Reqs( pe+1 ), Ierr )
      ENDDO
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking.
!
        CALL MPI_RECV( nr, 1, CPP_MPI_INTEGER, pe, pe+3000,                   &
     &                 InComm, Status, Ierr )
        LenOutVector(pe + 1) = nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )
#else
      CALL MPI_ALLTOALL( LenInVector, 1, CPP_MPI_INTEGER,                     &
     &                   LenOutVector, 1, CPP_MPI_INTEGER,                    &
     &                   InComm, Ierr )
#endif
!
! Over all processes
!
      icnt = 1
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        IF ( nr .gt. 0 ) THEN
          CALL MPI_ISEND( InVector( icnt ), nr,                          &
     &                    CPP_MPI_REAL, pe, Iam+2000,                    &
     &                    InComm, Reqs( pe+1 ), Ierr )
        ELSE
          Reqs( pe+1 ) = MPI_REQUEST_NULL
        ENDIF
        icnt = icnt + nr
      ENDDO

!
! Over all processes
!
      icnt = 1
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking. 
!
        nr = LenOutVector(pe + 1)
        IF ( nr .gt. 0 ) THEN
          CALL MPI_RECV( OutVector( icnt ), nr,                          &
     &                   CPP_MPI_REAL, pe, pe+2000,                      &
     &                   InComm, Status, Ierr )
        ENDIF
        icnt = icnt + nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )
!
! WS 98.09.22 : This barrier needed to synchronize.  Why?
!
      CALL MPI_BARRIER( InComm, Ierr )
#endif
      CPP_LEAVE_PROCEDURE( "PAREXCHANGEVECTORREAL" )

      RETURN
!EOC
      END SUBROUTINE ParExchangeVectorReal
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParExchangeVectorInt --- Exchange a sparse packed vector
!
! !INTERFACE:  
      SUBROUTINE ParExchangeVectorInt ( InComm, LenInVector, InVector,   &
     &                                   LenOutVector, OutVector )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )   :: InComm            ! Communicator
      INTEGER, INTENT( IN )   :: LenInVector( * )  ! Length on each PE
      INTEGER, INTENT( IN )   :: InVector( * )     ! The input buffer

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )  :: LenOutVector( * ) ! Length on each PE
      INTEGER, INTENT( OUT )  :: OutVector( * )    ! The output buffer

! !DESCRIPTION:
!
!     This routine exchanges vectors stored in compressed format, i.e.,
!     in so-called compressed sparse row (CSR) format, with other
!     PEs.  In essence it first exchanges the lengths with
!     MPI\_Alltoall, then the exchange of the actual vectors (can be
!     different in size) using MPI\_AlltoallV.  Since the latter is
!     inefficient, it is simulated using MPI\_Isend and MPI\_Recv.
!
!     The USE\_ARENAS version is fully functional, but assumes a 
!     global communicator.
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_WAITALL, MPI_ALLTOALL
!
! !REVISION HISTORY:
!   98.03.17   Sawyer     Creation from F77 version
!   98.03.30   Sawyer     Removed assumed shape arrays due to problems
!   99.01.18   Sawyer     Added barrier for safety
!   99.03.08   Sawyer     USE_SHMEM version for CRAY only; untested
!   99.06.01   Sawyer     USE_SHMEM version revised per comments from Tom
!   00.07.28   Sawyer     Implemented with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER :: i, iscnt, ircnt, nr, pe, icnt, Nsize, Iam, Ierr
      INTEGER :: Status(MPI_STATUS_SIZE)
      INTEGER :: Reqs( Gsize ), Stats( Gsize*MPI_STATUS_SIZE )
#if defined( USE_ARENAS )
      INTEGER :: offset( MAX_PE ), goffset( MAX_PE )
#endif

      CPP_ENTER_PROCEDURE( "PAREXCHANGEVECTORINT" )
#if defined( USE_ARENAS )
      DO pe = 1, Gsize
        volume( pe, gid+1 ) = LenInVector( pe )
      ENDDO
      CALL mlp_barrier(gid,gsize)
      icnt = 0
      DO pe = 1, Gsize
        LenOutVector( pe ) = volume( gid+1, pe )
        offset( pe )  = icnt
        goffset( pe ) = icnt
        DO i = 1, gid
          goffset( pe ) = goffset( pe ) + volume( pe, i )
        ENDDO
        icnt = icnt + sum( volume( pe, 1:Gsize ) )
      ENDDO
      icnt = 0
      DO pe = 1, Gsize
        DO i = goffset(pe)+1, goffset(pe)+LenInVector(pe)
          icnt = icnt + 1
          IntBuf( i ) = Invector( icnt )
        ENDDO
      ENDDO
      CALL mlp_barrier(gid,gsize)
      DO i = 1, SUM(LenOutVector(1:Gsize))
        OutVector( i ) = IntBuf( offset(gid+1) + i )
      ENDDO
#else
      CALL MPI_COMM_SIZE( InComm, Nsize, Ierr )
      CALL MPI_COMM_RANK( InComm, Iam, Ierr )

#if defined( MY_ALLTOALL )
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        CALL MPI_ISEND( nr, 1,                                           &
     &                 MPI_INTEGER, pe, Iam+3000,                        &
     &                 InComm, Reqs( pe+1 ), Ierr )
      ENDDO
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking.
!
        CALL MPI_RECV( nr, 1,                                                 &
     &                 MPI_INTEGER, pe, pe+3000,                              &
     &                 InComm, Status, Ierr )
        LenOutVector(pe + 1) = nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )
#else
      CALL MPI_ALLTOALL( LenInVector, 1, CPP_MPI_INTEGER,                     &
     &                   LenOutVector, 1, CPP_MPI_INTEGER,                    &
     &                   InComm, Ierr )
#endif
!
! Over all processes
!
      icnt = 1
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        IF ( nr .gt. 0 ) THEN
          CALL MPI_ISEND( InVector( icnt ), nr,                               &
     &                    CPP_MPI_INTEGER, pe, Iam+2000,                      &
     &                    InComm, Reqs( pe+1 ), Ierr )
        ELSE
          Reqs( pe+1 ) = MPI_REQUEST_NULL
        ENDIF
        icnt = icnt + nr
      ENDDO

!
! Over all processes
!
      icnt = 1
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking. 
!
        nr = LenOutVector(pe + 1)
        IF ( nr .gt. 0 ) THEN
          CALL MPI_RECV( OutVector( icnt ), nr,                               &
     &                   CPP_MPI_INTEGER, pe, pe+2000,                        &
     &                   InComm, Status, Ierr )
        ENDIF
        icnt = icnt + nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )
!
! WS 98.09.22 : This barrier needed to synchronize.  Why?
!
      CALL MPI_BARRIER( InComm, Ierr )
#endif
      CPP_LEAVE_PROCEDURE( "PAREXCHANGEVECTORINT" )

      RETURN
!EOC
      END SUBROUTINE ParExchangeVectorInt
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollectiveBarrier --- Barrier: Simplest collective op.
!
! !INTERFACE:
      SUBROUTINE ParCollectiveBarrier( InComm )

! !USES:
      IMPLICIT NONE
! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator

! !DESCRIPTION:
!
!     This routine performs a barrier only within the communicator InComm
!     
!     The USE\_ARENAS version is fully functional, but assumes a 
!     global communicator.
!
! !REVISION HISTORY:
!   00.09.10   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
#if !defined( USE_ARENAS )
      INTEGER Ierror
#endif

#if defined( USE_ARENAS )
      CALL mlp_barrier(gid,gsize)
#else
      CALL MPI_Barrier(InComm, Ierror )
#endif
      RETURN
!EOC
      END SUBROUTINE ParCollectiveBarrier
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective0D --- Perform global Collective of a scalar
!
! !INTERFACE:
      SUBROUTINE ParCollective0D( InComm, Op, Var )

! !USES:
      IMPLICIT NONE
! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL), INTENT( INOUT ) :: Var	! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
!     The USE\_ARENAS version is fully functional, but assumes a 
!     global communicator.
!
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
#if !defined( USE_ARENAS )
      INTEGER Ierror
      REAL(CPP_REAL)    Tmp
#else
      INTEGER Ipe
#endif

#if defined( USE_ARENAS )
      DataBuf( GID+1 ) = Var
      CALL mlp_barrier(gid,gsize)

      IF ( Op .EQ. CPP_SUM_OP ) THEN
        Var = 0.0            !  Must be done the same on all PEs
        DO Ipe = 1, Gsize
          Var = Var + DataBuf( Ipe )
        ENDDO
      ELSEIF ( Op .EQ. CPP_MAX_OP ) THEN
        DO Ipe = 1, Gsize
          IF ( Var .LT. DataBuf( Ipe ) ) Var = DataBuf( Ipe+1 ) 
        ENDDO
      ELSEIF ( Op .EQ. CPP_MIN_OP ) THEN
        DO Ipe = 1, Gsize
          IF ( Var .GT. DataBuf( Ipe ) ) Var = DataBuf( Ipe ) 
        ENDDO
      ELSEIF ( Op .EQ. CPP_BCST_OP ) THEN
        DO Ipe = 2, Gsize
          DataBuf( Ipe ) = DataBuf( 1 )
        ENDDO
      ELSEIF ( Op .EQ. CPP_BCAST_OP ) THEN
        Var = DataBuf( 1 ) 
      ENDIF
      CALL mlp_barrier(gid,gsize)
#else
      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, 1, CPP_MPI_REAL, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, 1, CPP_MPI_REAL,                   &
     &                    Op, InComm, Ierror )
        Var = Tmp
      ENDIF
#endif
      RETURN
!EOC
      END SUBROUTINE ParCollective0D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective1D --- Perform component-wise global Collective of a vector
!
! !INTERFACE:
      SUBROUTINE ParCollective1D( InComm, Op, Im, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im       ! Size of 1-D array

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL), INTENT( INOUT ) :: Var(Im) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
!     The USE\_ARENAS version is fully functional, but assumes a 
!     global communicator.
!
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
#if !defined( USE_ARENAS )
      INTEGER Ierror
      REAL(CPP_REAL)    Tmp(Im)
#else
      INTEGER Icnt, Ipe, I
#endif

#if defined( USE_ARENAS )
      Icnt   = GID*Im
      DO i=1,Im
        Icnt = Icnt + 1
        DataBuf( Icnt ) = Var( i )
      ENDDO
      CALL mlp_barrier(gid,gsize)

      IF ( Op .EQ. CPP_SUM_OP ) THEN
        Var = 0.0            !  Must be done the same on all PEs
        Icnt = 0
        DO Ipe = 1, Gsize
          DO i=1,Im
            Icnt = Icnt + 1
            Var( i ) = Var( i ) + DataBuf( Icnt )
          ENDDO
        ENDDO
      ELSEIF ( Op .EQ. CPP_MAX_OP ) THEN
        Icnt = 0
        DO Ipe = 0, Gsize-1
          IF ( Ipe .NE. GID ) THEN
            DO i=1,Im
              Icnt = Icnt + 1
              IF ( Var( i ) .LT. DataBuf( Icnt ) )                       &
     &                  Var( i ) = DataBuf( Icnt ) 
            ENDDO
          ELSE
            Icnt = Icnt + Im
          ENDIF
        ENDDO
      ELSEIF ( Op .EQ. CPP_MIN_OP ) THEN
        Icnt = 0
        DO Ipe = 0, Gsize-1
          IF ( Ipe .NE. GID ) THEN
            DO i=1,Im
              Icnt = Icnt + 1
              IF ( Var( i ) .GT. DataBuf( Icnt ) )                       &
     &                  Var( i ) = DataBuf( Icnt ) 
            ENDDO
          ELSE
            Icnt = Icnt + Im
          ENDIF
        ENDDO
      ELSEIF ( Op .EQ. CPP_BCAST_OP ) THEN
        Icnt = 0
        DO i=1,Im
          Icnt = Icnt + 1
          Var( i ) = DataBuf( Icnt ) 
        ENDDO
      ENDIF
      CALL mlp_barrier(gid,gsize)
#else
      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im, CPP_MPI_REAL, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, Im, CPP_MPI_REAL,                  &
     &                      Op, InComm, Ierror )
        Var = Tmp
      ENDIF
#endif
      RETURN
!EOC
      END SUBROUTINE ParCollective1D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective2D --- Perform component-wise collective operation
!
! !INTERFACE:
      SUBROUTINE ParCollective2D( InComm, Op, Im, Jm, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm     ! Communicator
      INTEGER, INTENT( IN ) :: Op         ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im         ! First dimension of 2-D array
      INTEGER, INTENT( IN ) :: Jm         ! Second dimension of 2-D array

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL), INTENT( INOUT ) :: Var(Im,Jm) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
!     The USE\_ARENAS version is fully functional, but assumes a 
!     global communicator.
!
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
#if !defined( USE_ARENAS )
      INTEGER Ierror
      REAL(CPP_REAL)    Tmp(Im,Jm)
#else
      INTEGER Icnt, Length, Ipe, I, J
#endif

#if defined( USE_ARENAS )
      Length = Im*Jm
      Icnt   = GID*Length
      DO j=1,Jm
        DO i=1,Im
          Icnt = Icnt + 1
          DataBuf( Icnt ) = Var( i,j )
        ENDDO
      ENDDO
      CALL mlp_barrier(gid,gsize)

      IF ( Op .EQ. CPP_SUM_OP ) THEN
        Var = 0.0            !  Must be done the same on all PEs
        Icnt = 0
        DO Ipe = 1, Gsize
          DO j=1,Jm
            DO i=1,Im
              Icnt = Icnt + 1
              Var( i,j ) = Var( i,j ) + DataBuf( Icnt )
            ENDDO
          ENDDO
        ENDDO
      ELSEIF ( Op .EQ. CPP_MAX_OP ) THEN
        Icnt = 0
        DO Ipe = 0, Gsize-1
          IF ( Ipe .NE. GID ) THEN
            DO j=1,Jm
              DO i=1,Im
                Icnt = Icnt + 1
                IF ( Var( i,j ) .LT. DataBuf( Icnt ) )                   &
     &                    Var( i,j ) = DataBuf( Icnt ) 
              ENDDO
            ENDDO
          ELSE
            Icnt = Icnt + Length
          ENDIF
        ENDDO
      ELSEIF ( Op .EQ. CPP_MIN_OP ) THEN
        Icnt = 0
        DO Ipe = 0, Gsize-1
          IF ( Ipe .NE. GID ) THEN
            DO j=1,Jm
              DO i=1,Im
                Icnt = Icnt + 1
                IF ( Var( i,j ) .GT. DataBuf( Icnt ) )                   &
     &                    Var( i,j ) = DataBuf( Icnt ) 
              ENDDO
            ENDDO
          ELSE
            Icnt = Icnt + Length
          ENDIF
        ENDDO
      ELSEIF ( Op .EQ. CPP_BCAST_OP ) THEN
        Icnt = 0
        DO j=1,Jm
          DO i=1,Im
            Icnt = Icnt + 1
            Var( i,j ) = DataBuf( Icnt ) 
          ENDDO
        ENDDO
      ENDIF
      CALL mlp_barrier(gid,gsize)
#else
      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im*Jm, CPP_MPI_REAL, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, Im*Jm, CPP_MPI_REAL,               &
     &                    Op, InComm, Ierror )
        Var = Tmp
      ENDIF
#endif
      RETURN
!EOC
      END SUBROUTINE ParCollective2D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective3D --- Perform component-wise global Collective of a vector
!
! !INTERFACE:
      SUBROUTINE ParCollective3D( InComm, Op, Im, Jm, Lm, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm     ! Communicator
      INTEGER, INTENT( IN ) :: Op         ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im         ! First dimension of 3-D array
      INTEGER, INTENT( IN ) :: Jm         ! Second dimension of 3-D array
      INTEGER, INTENT( IN ) :: Lm         ! Third dimension of 3-D array

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL), INTENT( INOUT ):: Var(Im,Jm,LM) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
!     The USE\_ARENAS version is fully functional, but assumes a 
!     global communicator.
!
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
#if !defined( USE_ARENAS )
      INTEGER Ierror
      REAL(CPP_REAL) Tmp(Im,Jm,Lm)
#else
      INTEGER Icnt, Length, Ipe, I, J, L
#endif

#if defined( USE_ARENAS )
      Length = Im*Jm*Lm
      Icnt   = GID*Length
      DO l=1,Lm
        DO j=1,Jm
          DO i=1,Im
            Icnt = Icnt + 1
            DataBuf( Icnt ) = Var( i,j,l )
          ENDDO
        ENDDO
      ENDDO
      CALL mlp_barrier(gid,gsize)

      IF ( Op .EQ. CPP_SUM_OP ) THEN
        Var = 0.0            !  Must be done the same on all PEs
        Icnt = 0
        DO Ipe = 1, Gsize
          DO l=1,Lm
            DO j=1,Jm
              DO i=1,Im
                Icnt = Icnt + 1
                Var( i,j,l ) = Var( i,j,l ) + DataBuf( Icnt )
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSEIF ( Op .EQ. CPP_MAX_OP ) THEN
        Icnt = 0
        DO Ipe = 0, Gsize-1
          IF ( Ipe .NE. GID ) THEN
            DO l=1,Lm
              DO j=1,Jm
                DO i=1,Im
                  Icnt = Icnt + 1
                  IF ( Var( i,j,l ) .LT. DataBuf( Icnt ) )               &
     &                      Var( i,j,l ) = DataBuf( Icnt ) 
                ENDDO
              ENDDO
            ENDDO
          ELSE
            Icnt = Icnt + Length
          ENDIF
        ENDDO
      ELSEIF ( Op .EQ. CPP_MIN_OP ) THEN
        Icnt = 0
        DO Ipe = 0, Gsize-1
          IF ( Ipe .NE. GID ) THEN
            DO l=1,Lm
              DO j=1,Jm
                DO i=1,Im
                  Icnt = Icnt + 1
                  IF ( Var( i,j,l ) .GT. DataBuf( Icnt ) )               &
     &                      Var( i,j,l ) = DataBuf( Icnt ) 
                ENDDO
              ENDDO
            ENDDO
          ELSE
            Icnt = Icnt + Length
          ENDIF
        ENDDO
      ELSEIF ( Op .EQ. CPP_BCAST_OP ) THEN
        Icnt = 0
        DO l=1,Lm
          DO j=1,Jm
            DO i=1,Im
              Icnt = Icnt + 1
              Var( i,j,l ) = DataBuf( Icnt ) 
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      CALL mlp_barrier(gid,gsize)
#else
      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im*Jm*Lm, CPP_MPI_REAL, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, Im*Jm*Lm, CPP_MPI_REAL,            &
     &                    Op, InComm, Ierror )
        Var = Tmp
      ENDIF
#endif
      RETURN
!EOC
      END SUBROUTINE ParCollective3D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective0DInt --- Perform global Collective of a scalar
!
! !INTERFACE:
      SUBROUTINE ParCollective0DInt( InComm, Op, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)

! !INPUT/OUTPUT PARAMETERS:
      INTEGER, INTENT( INOUT ) :: Var	! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
!     The USE\_ARENAS version is fully functional, but assumes a 
!     global communicator.
!
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
#if !defined( USE_ARENAS )
      INTEGER Ierror
      INTEGER    Tmp
#else
      INTEGER Ipe
#endif


#if defined( USE_ARENAS )
      IntBuf( GID+1 ) = Var
      CALL mlp_barrier(gid,gsize)

      IF ( Op .EQ. CPP_SUM_OP ) THEN
        Var = 0.0            !  Must be done the same on all PEs
        DO Ipe = 1, Gsize
          Var = Var + IntBuf( Ipe )
        ENDDO
      ELSEIF ( Op .EQ. CPP_MAX_OP ) THEN
        DO Ipe = 1, Gsize
          IF ( Var .LT. IntBuf( Ipe ) ) Var = IntBuf( Ipe+1 ) 
        ENDDO
      ELSEIF ( Op .EQ. CPP_MIN_OP ) THEN
        DO Ipe = 1, Gsize
          IF ( Var .GT. IntBuf( Ipe ) ) Var = IntBuf( Ipe ) 
        ENDDO
      ELSEIF ( Op .EQ. CPP_BCAST_OP ) THEN
        Var = IntBuf(1) 
      ENDIF
      CALL mlp_barrier(gid,gsize)
#else
      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, 1, CPP_MPI_INTEGER, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var,Tmp,1,CPP_MPI_INTEGER,Op,InComm,Ierror )
        Var = Tmp
      ENDIF
#endif
      RETURN
!EOC
      END SUBROUTINE ParCollective0DInt
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective1DInt --- Perform component-wise global 
!                                  collective operations of int vector
!
! !INTERFACE:
      SUBROUTINE ParCollective1DInt( InComm, Op, Im, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im       ! Size of 1-D array

! !INPUT/OUTPUT PARAMETERS:
      INTEGER, INTENT( INOUT ) :: Var(Im) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
!     The USE\_ARENAS version is fully functional, but assumes a 
!     global communicator.
!
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
#if !defined( USE_ARENAS )
      INTEGER Ierror
      INTEGER Tmp(Im)
#else
      INTEGER Icnt, Ipe, I
#endif

#if defined( USE_ARENAS )
      Icnt   = GID*Im
      DO i=1,Im
        Icnt = Icnt + 1
        IntBuf( Icnt ) = Var( i )
      ENDDO
      CALL mlp_barrier(gid,gsize)

      IF ( Op .EQ. CPP_SUM_OP ) THEN
        Var = 0.0            !  Must be done the same on all PEs
        Icnt = 0
        DO Ipe = 1, Gsize
          DO i=1,Im
            Icnt = Icnt + 1
            Var( i ) = Var( i ) + IntBuf( Icnt )
          ENDDO
        ENDDO
      ELSEIF ( Op .EQ. CPP_MAX_OP ) THEN
        Icnt = 0
        DO Ipe = 0, Gsize-1
          IF ( Ipe .NE. GID ) THEN
            DO i=1,Im
              Icnt = Icnt + 1
              IF ( Var(i) .LT. IntBuf(Icnt) ) Var(i) = IntBuf(Icnt) 
            ENDDO
          ELSE
            Icnt = Icnt + Im
          ENDIF
        ENDDO
      ELSEIF ( Op .EQ. CPP_MIN_OP ) THEN
        Icnt = 0
        DO Ipe = 0, Gsize-1
          IF ( Ipe .NE. GID ) THEN
            DO i=1,Im
              Icnt = Icnt + 1
              IF ( Var(i) .GT. IntBuf(Icnt) ) Var(i) = IntBuf(Icnt) 
            ENDDO
          ELSE
            Icnt = Icnt + Im
          ENDIF
        ENDDO
      ELSEIF ( Op .EQ. CPP_BCAST_OP ) THEN
        DO i=1,Im
          Var(i) = IntBuf(i) 
        ENDDO
      ENDIF
      CALL mlp_barrier(gid,gsize)
#else
      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im, CPP_MPI_INTEGER, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var,Tmp,Im,CPP_MPI_INTEGER,Op,InComm,Ierror )
        Var = Tmp
      ENDIF
#endif
      RETURN
!EOC
      END SUBROUTINE ParCollective1DInt
!-----------------------------------------------------------------------
#endif
      END MODULE parutilitiesmodule

