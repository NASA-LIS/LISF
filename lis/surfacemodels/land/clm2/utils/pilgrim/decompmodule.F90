!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
      MODULE decompmodule
!BOP
!
! !MODULE: decompmodule
!
! !USES:
#include "debug.h"

      IMPLICIT NONE

!
! !DESCRIPTION:
!
!      This module provides the DecompType and its create and destroy 
!      routines.
!      \begin{center}
!      \begin{tabular}{|l|l|} \hline \hline
!         DecompType        & Type to describe a decomposition \cr \hline
!         DecompFree        & Destroy a decomposition \cr \hline
!         DecompCopy        & Copy decomposition to newly created one \cr \hline 
!         DecompPermute     & Permute decomposition \cr \hline 
!         DecompRegular1D   & Create a 1-D decomposition \cr \hline 
!         DecompRegular2D   & Create a 2-D decomposition \cr \hline 
!         DecompRegular3D   & Create a 3-D decomposition \cr \hline 
!         DecompCreateIrr   & Create an irregular 1-D decomposition \cr \hline
!         DecompCreateTags  & Create a decomposition from Pe and Tags \cr \hline
!         DecompGlobalToLocal& Map a global index to a local one \cr \hline
!         DecompLocalToGlobal& Map a local index to a global one \cr 
!         \hline  \hline
!      \end{tabular}
!      \end{center}
!
!      The decomposition type contains the sizes of the global array,
!      the number of entries on each PE, and for each PE a list
!      of "runs", i.e., the starting and finishing global indices
!      or "tags" whose inclusive array section resides on that PE.
!      Clearly this method of decomposition is only efficient if
!      there are long runs, i.e., long array sections which are 
!      mapped to one PE.  A random decomposition will cause poor
!      results.
!
!      The decomposition is thus very efficient for 1-D, 2-D or 3-D
!      block distributions (particularly for 1-D distributions, where
!      there is one "run" per processor).  Problems may occur for
!      an irregular decomposition (which is by definition 1-D).  If
!      there is little correspondence between the global indices of the 
!      entries and the actual decomposition (e.g., the tags are
!      assigned randomly), then there will be many runs, most
!      containing only one tag, and the resulting instance of
!      DecompType will be very large.  Fortunately, most applications
!      assign tags to entries in some sort of contiguous fashion, 
!      which is then quite appropriate for this data structure.
!
!      All numbering of multi-dimensional arrays is ROW-MAJOR, that
!      is, first in the X direction and then in the Y (and then,
!      if appropriate, in Z).  This is true for both the 2-D and
!      3-D data sets as also the Cartesian description of the PEs.
!
!      There is one glaring feature of DecompType.  It is
!      supposed to be a `one-size-fits-all' description of the
!      decomposition (with the exception of the random indexing
!      mentioned above).  Unfortunately, to describe 2-D and 3-D
!      regions, it is necessary to carry additional dimension
!      information in order have complete information for the 
!      mapping.  This means that 2-D and 3-D decompositions 
!      inherently carry more information than a 1-D decomposition.
!      Thus it {\it is} possible to use a decomposition created
!      with the Regular2D or Regular3D routines to describe the
!      corresponding decomposition when the 2-D or 3-D array is
!      viewed as a 1-D array, but it is clearly {\it not}
!      possible to use a decomposition created with Regular1D
!      to describe the decomposition of a 2-D or 3-D array
!      --- the appropriate information just is not there.
!      
! !REVISION HISTORY:
!   97.07.22   Sawyer     Creation
!   97.09.01   Sawyer     Release date
!   97.11.06   Sawyer     Addition of row and column communicators
!   97.01.24   Sawyer     Added support for non-MPI derived types solution
!   97.01.29   Sawyer     Minor revisions for production service
!   98.01.30   Sawyer     Added DecompCopy
!   98.02.04   Sawyer     Removed Comm, CommRow and CommCol from DecompType
!   98.03.13   Sawyer     Removed DecompTypeOld, brushed up for walkthrough
!   98.03.19   Sawyer     Minor corrections after walkthrough
!   98.05.02   Sawyer     Added DecompPermute
!   98.05.11   Sawyer     Removed Permutation from all but DecompPermute
!   99.01.19   Sawyer     Minor cleaning
!   00.07.07   Sawyer     Removed DimSizes; decomp is now 1D only
!   00.11.12   Sawyer     Added DecompCreateTags and DecompInfo
!   01.02.03   Sawyer     Updated for free format; corrected DecompCreateTags
!   01.03.20   Sawyer     Added DecompRegular3DOrder
!
! !PUBLIC TYPES:
      PUBLIC DecompType, DecompCreate, DecompFree
      PUBLIC DecompCopy, DecompPermute
      PUBLIC DecompGlobalToLocal, DecompLocalToGlobal, DecompInfo
! Currently not yet tested:
      PUBLIC DecompRegular3DOrder

      INTERFACE     DecompCreate
        MODULE PROCEDURE DecompRegular1D
        MODULE PROCEDURE DecompRegular2D
        MODULE PROCEDURE DecompRegular3D
        MODULE PROCEDURE DecompCreateIrr
        MODULE PROCEDURE DecompCreateTags
      END INTERFACE
 
! Decomposition info

       TYPE Lists
         INTEGER, POINTER     :: StartTags(:) ! Start of tag run
         INTEGER, POINTER     :: EndTags(:)   ! Start of tag run
       END TYPE Lists

       TYPE DecompType
         INTEGER              :: GlobalSize   ! Size in each dimension
         INTEGER, POINTER     :: NumEntries(:)! Number of entries per PE
         TYPE(Lists), POINTER :: Head(:)      ! Array of pointers 
       END TYPE DecompType

!EOP
      CONTAINS

!---------------------------------------------------------------------
!BOP
! !IROUTINE: DecompFree --- Free a decomposition
!
! !INTERFACE:
      SUBROUTINE DecompFree ( Decomp )
! !USES:
      IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:
      TYPE(DecompType), INTENT( INOUT ):: Decomp  ! Decomp information
!
! !DESCRIPTION:
!     Free the decomposition -- deallocate the data structures.
!
! !SYSTEM ROUTINES:
!     ASSOCIATED, DEALLOCATE
!
! !REVISION HISTORY:
!   98.01.30   Sawyer     Creation
!
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER  :: I, NPEs
!
      CPP_ENTER_PROCEDURE( "DECOMPFREE" )

      IF ( ASSOCIATED( Decomp%NumEntries ) )                             &
     &               DEALLOCATE( Decomp%NumEntries )
      IF ( ASSOCIATED( Decomp%Head ) ) THEN
        NPEs = SIZE( Decomp%Head )
        DO I = 1, NPEs
!
! Copy the number of entries on each PE
!
          IF ( ASSOCIATED( Decomp%Head(I)%StartTags ) )                  &
     &               DEALLOCATE( Decomp%Head(I)%StartTags )
          IF ( ASSOCIATED( Decomp%Head(I)%EndTags ) )                    &
     &               DEALLOCATE( Decomp%Head(I)%EndTags )
        ENDDO
        DEALLOCATE( Decomp%Head )
      ENDIF

      CPP_LEAVE_PROCEDURE( "DECOMPFREE" )
      RETURN
!EOC
      END SUBROUTINE DecompFree
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompCopy --- Copy one decomposition to another
!
! !INTERFACE:
      SUBROUTINE DecompCopy ( DecompIn, DecompOut )
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(DecompType), INTENT( IN )   :: DecompIn  ! Decomp information
!
! !OUTPUT PARAMETERS:
      TYPE(DecompType), INTENT( OUT )  :: DecompOut ! Decomp information
!
! !DESCRIPTION:
!
!   Creates an output decomposition and copies the DecompIn input values 
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   98.01.30   Sawyer     Creation
!
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER  :: I, J, NDims, NPEs, NRuns
!
      CPP_ENTER_PROCEDURE( "DECOMPCOPY" )
!
! Copy the size of the global array
!
      DecompOut%GlobalSize = DecompIn%GlobalSize

!
! Allocate the number of entries and list head arrays
!
      NPEs = SIZE( DecompIn%NumEntries )
      CPP_ASSERT_F90( SIZE( DecompIn%Head ) .EQ. NPEs )
      ALLOCATE( DecompOut%NumEntries( NPEs ) )
      ALLOCATE( DecompOut%Head( NPEs ) )

      DO I = 1, NPEs
!
! Copy the number of entries on each PE
!
        DecompOut%NumEntries( I ) = DecompIn%NumEntries( I )
        NRuns = SIZE( DecompIn%Head( I )%StartTags )
        CPP_ASSERT_F90( SIZE( DecompIn%Head( I )%EndTags ) .EQ. NRuns )
!
! Allocate and copy the array of runs
!
        ALLOCATE( DecompOut%Head(I)%StartTags( NRuns ) )
        ALLOCATE( DecompOut%Head(I)%EndTags( NRuns ) )
        DO J = 1, NRuns
          DecompOut%Head(I)%StartTags(J) = DecompIn%Head(I)%StartTags(J)
          DecompOut%Head(I)%EndTags(J) = DecompIn%Head(I)%EndTags(J)
        ENDDO
      ENDDO

      CPP_LEAVE_PROCEDURE( "DECOMPCOPY" )
      RETURN
!EOC
      END SUBROUTINE DecompCopy
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompPermute --- Permute one decomposition to another
!
! !INTERFACE:
      SUBROUTINE DecompPermute ( Permutation, Decomp )
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      INTEGER :: Permutation( : )                  ! Permutation

! !INPUT/OUTPUT PARAMETERS:
      TYPE(DecompType), INTENT( INOUT ) :: Decomp  ! Decomp information
!
!
! !DESCRIPTION:
!
!   Permutes the PE assignment of a given decomposition. Confusion will
!   always arise about whether this is a forward or backward
!   transformation.  Picture it this way: draw the array and slice it
!   up as indicated by the distribution.  The resulting boxes are of
!   course indexed by natural numbering 1, 2, 3, 4, ...  (these are
!   the virtual one-based PEs). Now write the true PE numbering
!   (one-based) as you would like it.  The resulting array is Perm.
!
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   98.05.02   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER, POINTER     :: NumEntries(:)! Number of entries
      TYPE(Lists), POINTER :: Head(:)      ! Array of pointers 
      INTEGER              :: I, J, NPEs, NRuns, TruePE
!
      CPP_ENTER_PROCEDURE( "DECOMPPERMUTE" )
!
! Allocate the number of entries and list head arrays
!
      NPEs = SIZE( Decomp%NumEntries )
      ALLOCATE( NumEntries( NPEs ) )
      DO I = 1, NPEs
        TruePE = Permutation( I )
        NumEntries( TruePE ) = Decomp%NumEntries( I )
      ENDDO
!
! Deallocate old NumEntries and put the new pointer in its place
!
      DEALLOCATE( Decomp%NumEntries )
      Decomp%NumEntries => NumEntries
      NULLIFY( NumEntries )

!
! Allocate and set the permuted Lists called with pointer Head
!
      ALLOCATE( Head( NPEs ) )
      DO I = 1, NPEs
        TruePE = Permutation( I )
        NRuns = SIZE( Decomp%Head(I)%StartTags )
        CPP_ASSERT_F90( SIZE( Decomp%Head(I)%EndTags ) .EQ. NRuns )
!
! Allocate and permute the array of runs
!
        ALLOCATE( Head(TruePE)%StartTags(NRuns) )
        ALLOCATE( Head(TruePE)%EndTags(NRuns) )
        DO J = 1, NRuns
          Head(TruePE)%StartTags(J) = Decomp%Head(I)%StartTags(J)
          Head(TruePE)%EndTags(J)   = Decomp%Head(I)%EndTags(J)
        ENDDO
      ENDDO
!
! Deallocate the arrays of starting and ending tags
!
      DO I = 1, NPEs
        DEALLOCATE( Decomp%Head(I)%StartTags )
        DEALLOCATE( Decomp%Head(I)%EndTags )
      ENDDO
!
! Deallocate the heads to the Lists
!
      DEALLOCATE( Decomp%Head )

!
! Link the new head to that in the decomposition
!
      Decomp%Head => Head
      
      NULLIFY( Head )

      CPP_LEAVE_PROCEDURE( "DECOMPPERMUTE" )
      RETURN
!EOC
      END SUBROUTINE DecompPermute
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompRegular1D --- Create a decomposition for a 1-D grid
!
! !INTERFACE:
      SUBROUTINE DecompRegular1D ( NPEs, Dist, Decomp )
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )            :: NPEs     ! Number of PEs
      INTEGER, INTENT( IN )            :: Dist(:)  ! Distribution in X
!
! !OUTPUT PARAMETERS:
      TYPE(DecompType), INTENT( OUT )  :: Decomp   ! Decomp information
!
! !DESCRIPTION:
!     Creates a variable block decomposition for a regular 1-D grid
!     (this is also known as a "block-general" distribution).  The
!     decomposition is given through the Dist distribution 
!     which contains the number of entries on each PE.  
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   98.01.19   Sawyer     Creation
!   98.01.22   Sawyer     Corrections, TESTED
!   98.05.11   Sawyer     Removed Perm from arglist -- see DecompPermute
!   00.07.07   Sawyer     Removed use of DimSizes(:) array
!
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER  :: I, Counter
!
      CPP_ENTER_PROCEDURE( "DECOMPREGULAR1D" )
!
      CPP_ASSERT_F90( NPEs .EQ. SIZE( Dist ) )
!
! The head contains NPEs pointers to the tag lists.
!
      Decomp%GlobalSize = SUM(Dist)
      ALLOCATE( Decomp%NumEntries( NPEs ) )
      ALLOCATE( Decomp%Head( NPEs ) )
      Counter = 0
      DO I = 1, NPEs
        Decomp%NumEntries(I) = Dist(I)
!
! Since this is a regular distribution there is only one run of tags per PE.
!
        NULLIFY( Decomp%Head(I)%StartTags )
        NULLIFY( Decomp%Head(I)%EndTags )
        ALLOCATE( Decomp%Head(I)%StartTags(1) )
        ALLOCATE( Decomp%Head(I)%EndTags(1) )
!
! The starting and ending tags are immediately determined from
! the decomposition arrays  
!
        Decomp%Head(I)%StartTags(1) = Counter+1
        Counter = Counter + Dist(I)
        Decomp%Head(I)%EndTags(1) = Counter
      ENDDO

      CPP_LEAVE_PROCEDURE( "DECOMPREGULAR1D" )
      RETURN
!EOC
      END SUBROUTINE DecompRegular1D
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompRegular2D --- Create a decomposition for a 2-D grid
!
! !INTERFACE:
      SUBROUTINE DecompRegular2D( NPEsX, NPEsY, Xdist, Ydist, Decomp )
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )            :: NPEsX    ! Number of PEs in X
      INTEGER, INTENT( IN )            :: NPEsY    ! Number of PEs in Y
      INTEGER, INTENT( IN )            :: Xdist(:) ! Distribution in X
      INTEGER, INTENT( IN )            :: Ydist(:) ! Distribution in Y
!
! !OUTPUT PARAMETERS:
      TYPE(DecompType), INTENT( OUT )  :: Decomp  ! Decomp information
!
!
! !DESCRIPTION:
!     Creates a variable block-block decomposition for a regular 
!     2-D grid.  The decomposition is given through the Xdist and 
!     Ydist distributions, which contain the number of entries on 
!     each PE in that dimension.  This routine thus defines
!     a rectangular "checkerboard" distribution.
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   98.01.19   Sawyer     Creation
!   98.01.22   Sawyer     Corrections, TESTED
!   98.05.11   Sawyer     Removed Perm from arglist -- see DecompPermute
!   00.07.07   Sawyer     Removed use of DimSizes(:) array
!
! !BUGS:
!     This routine makes the assumption that the sum of the
!     distribution in each dimension adds up to the total 
!     number of entries in that dimension.  It will cause
!     problems if the actual local arrays are over- or 
!     under-allocated.  For example, if the local array is
!     allocated statically for the maximum size of the
!     array on any processor, problems will occur on those
!     PEs which have less than the maximum.
!
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER :: TruePE, I, J, K, Counter1, Counter2, SizeX, SizeY
!
      CPP_ENTER_PROCEDURE( "DECOMPREGULAR2D" )
!
! Some sanity checks
!
      CPP_ASSERT_F90( NPEsX .EQ. SIZE( Xdist ) )
      CPP_ASSERT_F90( NPEsY .EQ. SIZE( Ydist ) )
!
! The head contains NPEs pointers to the tag lists.
!
      SizeX = SUM(Xdist)
      SizeY = SUM(Ydist)
      Decomp%GlobalSize = SizeX * SizeY
      ALLOCATE( Decomp%NumEntries( NPEsX*NPEsY ) )
      ALLOCATE( Decomp%Head( NPEsX*NPEsY ) )
      Counter1 = 0
      DO J = 1, NPEsY
        DO I = 1, NPEsX
!
! WARNING!!!!  The definition of the PE is Row-major ordering
!
          TruePE = ( J-1 ) * NPEsX + I 

!
! The number of entries is the product of the local X, Y, Z allotment
!
          Decomp%NumEntries(TruePE) = Xdist(I)*Ydist(J)
!
! For each Y there is a separate run
!
          NULLIFY( Decomp%Head(TruePE)%StartTags )
          NULLIFY( Decomp%Head(TruePE)%EndTags )
          ALLOCATE( Decomp%Head(TruePE)%StartTags(Ydist(J)) )
          ALLOCATE( Decomp%Head(TruePE)%EndTags(Ydist(J)) )
          Counter2 = Counter1
          DO K = 1, Ydist(J)
!
! Since this is a regular distribution the definition of
! tags is dictated by Xdist(I), and appears Ydist(J) times
!
!
            Decomp%Head(TruePE)%StartTags(K) = Counter2 + 1
            Decomp%Head(TruePE)%EndTags(K) = Counter2 + Xdist(I)
            Counter2 = Counter2 + SizeX
          ENDDO
          Counter1 = Counter1 + Xdist(I)
        ENDDO
!
! Align the counter such that it points to the start of the next 
! block.  (Ydist(J)-1) since already one layer has been added in.
! Implicit assumption that SizeX = SUM( Xdist )
!
        Counter1 = Counter1 + SizeX*(Ydist(J)-1)
      ENDDO

      CPP_LEAVE_PROCEDURE( "DECOMPREGULAR2D" )
      RETURN
!EOC
      END SUBROUTINE DecompRegular2D
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompRegular3D --- Create a decomposition for a 3-D grid
!
! !INTERFACE:
      SUBROUTINE DecompRegular3D ( NPEsX, NPEsY, NPEsZ,                  &
     &                             Xdist, Ydist, Zdist, Decomp )
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )            :: NPEsX    ! Number of PEs in X
      INTEGER, INTENT( IN )            :: NPEsY    ! Number of PEs in Y
      INTEGER, INTENT( IN )            :: NPEsZ    ! Number of PEs in Z
      INTEGER, INTENT( IN )            :: Xdist(:) ! Distribution in X
      INTEGER, INTENT( IN )            :: Ydist(:) ! Distribution in Y
      INTEGER, INTENT( IN )            :: Zdist(:) ! Distribution in Z
!
! !OUTPUT PARAMETERS:
      TYPE(DecompType), INTENT( OUT )  :: Decomp  ! Decomp information
!
!
! !DESCRIPTION:
!     Creates a decomposition for a regular 3-D grid.  The
!     decomposition is given through the Xdist, Ydist, and Zdist
!     distributions, which contain the number of entries on 
!     each PE in that dimension.    This routine thus defines
!     a parallelopiped (SOMA-block) distribution.
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   98.01.19   Sawyer     Creation
!   98.05.11   Sawyer     Removed Perm from arglist -- see DecompPermute
!   00.07.07   Sawyer     Removed use of Sizes(:) array
!
! !BUGS:
!     This routine makes the assumption that the sum of the
!     distribution in each dimension adds up to the total 
!     number of entries in that dimension.  It will cause
!     problems if the actual local arrays are over- or 
!     under-allocated.  For example, if the local array is
!     allocated statically for the maximum size of the
!     array on any processor, problems will occur on those
!     PEs which have less than the maximum.
!
!     Currently untested (98.05.11)
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER  :: TruePE, Counter1, Counter2, Counter3
      INTEGER  :: I, J, K, L, M, N, SizeX, SizeY, SizeZ
!
      CPP_ENTER_PROCEDURE( "DECOMPREGULAR3D" )
!
! Some sanity checks
!
!
      CPP_ASSERT_F90( NPEsX .EQ. SIZE( Xdist ) )
      CPP_ASSERT_F90( NPEsY .EQ. SIZE( Ydist ) )
      CPP_ASSERT_F90( NPEsZ .EQ. SIZE( Zdist ) )
      CPP_ASSERT_F90( .NOT. ASSOCIATED( Decomp%Head ) )
!
! The head contains NPEs pointers to the tag lists.
!
      SizeX = SUM(Xdist)
      SizeY = SUM(Ydist)
      SizeZ = SUM(Zdist)
      Decomp%GlobalSize = SizeX * SizeY * SizeZ
      ALLOCATE( Decomp%NumEntries( NPEsX*NPEsY*NPEsZ ) )
      ALLOCATE( Decomp%Head( NPEsX*NPEsY*NPEsZ ) )
      Counter1 = 0
      DO K = 1, NPEsZ
        DO J = 1, NPEsY
          DO I = 1, NPEsX
!
! WARNING!!!!  The definition of the PE is Row-major ordering
!
            TruePE = (K-1)*NPEsX*NPEsY + (J-1)*NPEsX + I 
            NULLIFY( Decomp%Head(TruePE)%StartTags )
            NULLIFY( Decomp%Head(TruePE)%EndTags )
!
! The number of entries is the product of the local X, Y, Z allotment
!
            Decomp%NumEntries(TruePE) = Xdist(I)*Ydist(J)*Zdist(K)
!
! For each Z there are Y separate runs
!
            ALLOCATE( Decomp%Head(TruePE)%StartTags(Ydist(J)*Zdist(K)) )
            ALLOCATE( Decomp%Head(TruePE)%EndTags(Ydist(J)*Zdist(K)) )
            Counter2 = Counter1
            L = 0
            DO N = 1, Zdist(K)
              Counter3 = Counter2
              DO M = 1, Ydist(J)
!
!     Since this is a regular distribution the definition of
!     tags is dictated by Xdist(I), and appears Ydist(J) times
!
!
                L = L + 1
                Decomp%Head(TruePE)%StartTags(L) = Counter3 + 1
                Decomp%Head(TruePE)%EndTags(L) = Counter3 + Xdist(I)
                Counter3 = Counter3 + SizeX
              ENDDO
              Counter2 = Counter2 + SizeX*SizeY
            ENDDO
            Counter1 = Counter1 + Xdist(I)
          ENDDO
!
! Align the counter such that it points to the start of the next 
! block.  (Ydist(J)-1) since already one X layer has been added in.
! Implicit assumption that SizeX = SUM( Xdist )
!
          Counter1 = Counter1 + SizeX*(Ydist(J)-1)
        ENDDO
!
! Align the counter such that it points to the start of the next 
! block.  (Zdist(K)-1) since already one X-Y layer has been added in.
! Implicit assumption that SizeY = SUM( Ydist )
!
        Counter1 = Counter1 + SizeX*SizeY*(Zdist(K)-1)
      ENDDO

      CPP_LEAVE_PROCEDURE( "DECOMPREGULAR3D" )
      RETURN
!EOC
      END SUBROUTINE DecompRegular3D
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompRegular3Dorder --- Create a decomposition for a 3-D grid
!
! !INTERFACE:
      SUBROUTINE DecompRegular3Dorder( Order, NPEsX, NPEsY, NPEsZ,       &
     &                                 Xdist, Ydist, Zdist, Decomp )
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      CHARACTER(3), INTENT( IN )       :: Order    ! Dim. ordering
      INTEGER, INTENT( IN )            :: NPEsX    ! Number of PEs in X
      INTEGER, INTENT( IN )            :: NPEsY    ! Number of PEs in Y
      INTEGER, INTENT( IN )            :: NPEsZ    ! Number of PEs in Z
      INTEGER, INTENT( IN )            :: Xdist(:) ! Distribution in X
      INTEGER, INTENT( IN )            :: Ydist(:) ! Distribution in Y
      INTEGER, INTENT( IN )            :: Zdist(:) ! Distribution in Z
!
! !OUTPUT PARAMETERS:
      TYPE(DecompType), INTENT( OUT )  :: Decomp  ! Decomp information
!
! !DESCRIPTION:
!     Creates a variable block-block-block decomposition for a regular 
!     3-D grid, where the ordering of the PEs can be explicitly given
!     (see next paragraph). The decomposition is given through the 
!     Xdist, Ydist, and Zdist distributions, which contain the number 
!     of entries on each PE in that dimension.  This routine thus defines
!     a parallelopiped (SOMA-block) distribution.
!
!     With the string argument Order, the order of counting in the
!     3d PE space can be specified.  There are six possible values:
!     "xyz", "xzy", "yxz", "yzx", "zxy", and "zyx".  
!
!     The same as DecompRegular3Dorder could also be achieved by
!     using DecompRegular3D and then permuting the PE ownership
!     with DecompPermute.
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   01.03.20   Sawyer     Creation from DecompRegular3Dzy, added ordering
!
! !BUGS:
!   Not yet tested
!
!EOP
!---------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  :: TruePE, Counter1, Counter2, Counter3
      INTEGER  :: I, J, K, L, M, N, SizeX, SizeY, SizeZ
      INTEGER  :: Imult, Jmult, Kmult
!
      CPP_ENTER_PROCEDURE( "DECOMPREGULAR3DORDER" )
!
! Some sanity checks
!
!
      CPP_ASSERT_F90( NPEsX .EQ. SIZE( Xdist ) )
      CPP_ASSERT_F90( NPEsY .EQ. SIZE( Ydist ) )
      CPP_ASSERT_F90( NPEsZ .EQ. SIZE( Zdist ) )
      CPP_ASSERT_F90( .NOT. ASSOCIATED( Decomp%Head ) )

      IF ( Order=="xyz" ) THEN
! Looks like:      TruePE = (K-1)*NPEsX*NPEsY + (J-1)*NPEsX + (I-1) + 1
        Imult = 1
        Jmult = NPEsX
        Kmult = NPEsX*NPEsY
      ELSE IF ( Order=="xzy" ) THEN
! Looks like:      TruePE = (J-1)*NPEsX*NPEsZ + (K-1)*NPEsX + (I-1) + 1
        Imult = 1
        Jmult = NPEsX*NPEsZ
        Kmult = NPEsX
      ELSE IF ( Order=="yxz" ) THEN
! Looks like:      TruePE = (K-1)*NPEsY*NPEsX + (I-1)*NPEsY + (J-1) + 1
        Imult = NPEsY
        Jmult = 1
        Kmult = NPEsX*NPEsY
      ELSE IF ( Order=="yzx" ) THEN
! Looks like:      TruePE = (I-1)*NPEsY*NPEsZ + (K-1)*NPEsY + (J-1) + 1
        Imult = NPEsY*NPEsZ
        Jmult = 1
        Kmult = NPEsY
      ELSE IF ( Order=="zxy" ) THEN
! Looks like:      TruePE = (J-1)*NPEsX*NPEsZ + (I-1)*NPEsZ + (K-1) + 1
        Imult = NPEsZ
        Jmult = NPEsX*NPEsZ
        Kmult = 1
      ELSE IF ( Order=="zyx" ) THEN
! Looks like:      TruePE = (I-1)*NPEsY*NPEsZ + (J-1)*NPEsZ + (K-1) + 1
        Imult = NPEsY*NPEsZ
        Jmult = NPEsZ
        Kmult = 1
      ELSE 
! Looks like:      TruePE = (K-1)*NPEsX*NPEsY + (J-1)*NPEsX + (I-1) + 1
        print *, "Warning: DecompCreate3Dorder", Order, "not supported"
        print *, "         Continuing with XYZ ordering"
        Imult = 1
        Jmult = NPEsX
        Kmult = NPEsX*NPEsY
      ENDIF

!
! The head contains NPEs pointers to the tag lists.
!
      SizeX = SUM(Xdist)
      SizeY = SUM(Ydist)
      SizeZ = SUM(Zdist)
      Decomp%GlobalSize = SizeX * SizeY * SizeZ
      ALLOCATE( Decomp%NumEntries( NPEsX*NPEsY*NPEsZ ) )
      ALLOCATE( Decomp%Head( NPEsX*NPEsY*NPEsZ ) )
      Counter1 = 0

      DO K = 1, NPEsZ
        DO J = 1, NPEsY
          DO I = 1, NPEsX
!
! WARNING!!!!  The definition of the PE is Row-major ordering
!
            
            TruePE = (I-1)*Imult + (J-1)*Jmult + (K-1)*Kmult + 1 
!
! The number of entries is the product of the local X, Y, Z allotment
!
            Decomp%NumEntries(TruePE) = Xdist(I)*Ydist(J)*Zdist(K)
!
! For each Z there are Y separate runs
!
            ALLOCATE( Decomp%Head(TruePE)%StartTags(Ydist(J)*Zdist(K)) )
            ALLOCATE( Decomp%Head(TruePE)%EndTags(Ydist(J)*Zdist(K)) )
            Counter2 = Counter1
            L = 0
            DO N = 1, Zdist(K)
              Counter3 = Counter2
              DO M = 1, Ydist(J)
!
!     Since this is a regular distribution the definition of
!     tags is dictated by Xdist(I), and appears Ydist(J) times
!
!
                L = L + 1
                Decomp%Head(TruePE)%StartTags(L) = Counter3 + 1
                Decomp%Head(TruePE)%EndTags(L) = Counter3 + Xdist(I)
                Counter3 = Counter3 + SizeX
              ENDDO
              Counter2 = Counter2 + SizeX*SizeY
            ENDDO
            Counter1 = Counter1 + Xdist(I)
          ENDDO
!
! Align the counter such that it points to the start of the next 
! block.  (Ydist(J)-1) since already one X layer has been added in.
! Implicit assumption that SizeX = SUM( Xdist )
!
          Counter1 = Counter1 + SizeX*(Ydist(J)-1)
        ENDDO
!
! Align the counter such that it points to the start of the next 
! block.  (Zdist(K)-1) since already one X-Y layer has been added in.
! Implicit assumption that SizeY = SUM( Ydist )
!
        Counter1 = Counter1 + SizeX*SizeY*(Zdist(K)-1)
      ENDDO

      CPP_LEAVE_PROCEDURE( "DECOMPREGULAR3DORDER" )
      RETURN
!EOC
      END SUBROUTINE DecompRegular3DOrder
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompCreateIrr --- Decomposition for an irregular mesh
!
! !INTERFACE:
      SUBROUTINE DecompCreateIrr( NPEs, Pe, TotalPts, Decomp )
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )            :: NPEs     ! Number of PEs
      INTEGER, INTENT( IN )            :: Pe(:)    ! Processor location
      INTEGER, INTENT( IN )            :: TotalPts ! Number of points
!
! !OUTPUT PARAMETERS:
      TYPE(DecompType), INTENT( OUT )  :: Decomp  ! Decomp information
!
!
! !DESCRIPTION:
!     Creates a decomposition for a irregular 1-D mesh.  The
!     decomposition is given through the number of points and
!     an array containing the PE which each point is mapped to.
!     This mapping naturally assumes that the local numbering
!     is incrementally increasing as points are mapped to PEs.
!     This assumption is sufficient for most applications, but
!     another irregular mapping routine is available for more
!     complex cases.
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   98.01.19   Sawyer     Creation, with requirements from Jay Larson
!   98.11.02   Sawyer     Rewritten to requirements for Andrea Molod
!   00.07.07   Sawyer     Removed use of DimSizes(:) array
!   00.11.12   Sawyer     Changed argument order for overloading
!
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER  :: I, PEhold
      INTEGER  :: Counter( NPEs )
!
      CPP_ENTER_PROCEDURE( "DECOMPCREATEIRR" )
!
      CPP_ASSERT_F90( TotalPts .LE. SIZE( PE ) )
      CPP_ASSERT_F90( .NOT. ASSOCIATED( Decomp%Head ) )

!
! The head contains NPEs pointers to the tag lists.
!
      Decomp%GlobalSize = TotalPts
      ALLOCATE( Decomp%NumEntries( NPEs ) )
      ALLOCATE( Decomp%Head( NPEs ) )
!
! Perform over all points in the mapping
!
      PEhold= -1
      Counter = 0
      Decomp%NumEntries = 0
      DO I=1, TotalPts
        CPP_ASSERT_F90( ( PE( I ) .LT. NPEs .AND. PE( I ) .GE. 0 ) )
        IF ( PE( I ) .NE. PEhold ) THEN
          PEhold = PE( I )
          Counter( PEhold+1 ) = Counter( PEhold+1 ) + 1
        ENDIF
        Decomp%NumEntries(PEHold+1) = Decomp%NumEntries(PEHold+1) + 1
      ENDDO
      DO I=1, NPEs
!
! Now the amount of space to allocate is known.  It is acceptable
! to in allocated an array of size 0 (F90 Handbook, Section 6.5.1)
!
        ALLOCATE( Decomp%Head(I)%StartTags(Counter(I)) )
        ALLOCATE( Decomp%Head(I)%EndTags(Counter(I)) )
      ENDDO
!
! Perform over all points in the mapping
!
      PEhold= -1
      Counter = 0
      DO I=1, TotalPts
        IF ( PE( I ) .NE. PEhold ) THEN
!
! If not first entry, close up shop on previous run
!
          IF ( I .GT. 1 ) THEN
            Decomp%Head(PEhold+1)%EndTags(Counter(PEhold+1)) = I-1
          ENDIF
          PEhold = PE( I )
          Counter( PEhold+1 ) = Counter( PEhold+1 ) + 1
          Decomp%Head(PEhold+1)%StartTags(Counter(PEhold+1)) = I
        ENDIF
      ENDDO
!
! Clean up shop for the final run
!
      IF ( TotalPts .GT. 0 ) THEN
        Decomp%Head(PEhold+1)%EndTags(Counter(PEhold+1)) = TotalPts
      ENDIF

      CPP_LEAVE_PROCEDURE( "DECOMPCREATEIRR" )
      RETURN
!EOC
      END SUBROUTINE DecompCreateIrr
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompCreateTags --- Decomposition from Pe and Tags
!
! !INTERFACE:
      SUBROUTINE DecompCreateTags(Npes, Pe, TotalPts, Tags, Decomp )
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )            :: NPEs     ! Number of PEs
      INTEGER, INTENT( IN )            :: Pe(:)    ! Processor location
      INTEGER, INTENT( IN )            :: TotalPts ! Number of points
      INTEGER, INTENT( IN )            :: Tags(:)  ! Global index
!
! !OUTPUT PARAMETERS:
      TYPE(DecompType), INTENT( OUT )  :: Decomp   ! Decomp information
!
!
! !DESCRIPTION:
!     Creates a decomposition for a irregular mesh from the 
!     Pe ownership and the Tags.  This is a simple extension of 
!     DecompCreateIrr (previously DecompIrregular1D) but is
!     much more dangerous, since the user can define the Tags
!     (global indices) arbitrarily.
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   00.11.12   Sawyer     Creation from DecompCreateIrr
!
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER  :: I, PEhold, LastTag
      INTEGER  :: Counter( NPEs )
!
      CPP_ENTER_PROCEDURE( "DECOMPCREATETAGS" )
!
      CPP_ASSERT_F90( TotalPts .LE. SIZE( PE ) )
      CPP_ASSERT_F90( TotalPts .LE. SIZE( Tags ) )
      CPP_ASSERT_F90( .NOT. ASSOCIATED( Decomp%Head ) )

!
! The head contains NPEs pointers to the tag lists.
!
      Decomp%GlobalSize = TotalPts
      ALLOCATE( Decomp%NumEntries( NPEs ) )
      ALLOCATE( Decomp%Head( NPEs ) )
!
! Perform over all points in the mapping
!
      PEhold  = -1
      LastTag = -999999999
      Counter = 0
      Decomp%NumEntries = 0
      DO I=1, TotalPts
        CPP_ASSERT_F90( PE( I ) .LT. NPEs .AND. PE( I ) .GE. 0 )
        IF ( LastTag==0 .OR. Tags(I)/=LastTag+1 .OR. PE(I)/=PEhold ) THEN
          PEhold = PE( I )
          Counter( PEhold+1 ) = Counter( PEhold+1 ) + 1
        ENDIF
        Decomp%NumEntries(PEHold+1) = Decomp%NumEntries(PEHold+1) + 1
        LastTag = Tags(I)
      ENDDO

      DO I=1, NPEs
!
! Now the amount of space to allocate is known.  It is acceptable
! to in allocated an array of size 0 (F90 Handbook, Section 6.5.1)
!
        ALLOCATE( Decomp%Head(I)%StartTags(Counter(I)) )
        ALLOCATE( Decomp%Head(I)%EndTags(Counter(I)) )
      ENDDO

!
! Perform over all points in the domain
!
      PEhold  = -1
      LastTag = -999999999
      Counter = 0
      DO I=1, TotalPts
        IF ( LastTag==0 .OR. Tags(I)/=LastTag+1 .OR. PE(I)/=PEhold ) THEN
!
! If not first entry, close up shop on previous run
!
          IF ( I .GT. 1 ) THEN
            Decomp%Head(PEhold+1)%EndTags(Counter(PEhold+1)) = LastTag
          ENDIF
          PEhold = PE( I )
          Counter( PEhold+1 ) = Counter( PEhold+1 ) + 1
          Decomp%Head(PEhold+1)%StartTags(Counter(PEhold+1)) = Tags(I)
        ENDIF
        LastTag = Tags(I)
      ENDDO
!
! Clean up shop for the final run
!
      IF ( TotalPts .GT. 0 ) THEN
        Decomp%Head(PEhold+1)%EndTags(Counter(PEhold+1)) =Tags(TotalPts)
      ENDIF

      CPP_LEAVE_PROCEDURE( "DECOMPCREATETAGS" )
      RETURN
!EOC
      END SUBROUTINE DecompCreateTags
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompGlobalToLocal --- Map global index to local and PE
!
! !INTERFACE:
      SUBROUTINE DecompGlobalToLocal ( Decomp, Global, Local, Pe )
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(DecompType), INTENT( IN )   :: Decomp  ! Decomp information
      INTEGER, INTENT( IN )            :: Global  ! Global index
!
! !OUTPUT PARAMETERS:

      INTEGER, INTENT( OUT )  :: Local            ! Local index
      INTEGER, INTENT( OUT )  :: Pe               ! Pe location
!
!
! !DESCRIPTION:
!     Given a decomposition and a global index, this routine returns
!     the local index and PE location of that global tag.  If the
!     global index is not found, Local = 0, Pe = -1 is returned.
!
!     Note that this routine is not efficient by any stretch of the 
!     imagination --- only one index can be converted at a time.
!     In addition, a search procedure must be performed, whose 
!     efficiency is inversely proportional to the size of the decomposition
!     (in particular, to the number of "runs").  Conceptually this
!     mapping should be used only once in the program for
!     initialization, and subsequently all calculations should take
!     place using local indices.
!
! !SYSTEM ROUTINES:
!     SIZE
!
! !REVISION HISTORY:
!   98.03.20   Sawyer     Creation
!   01.03.17   Sawyer     Test for Global==0 (undefined element)
!
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER  :: I, J
      LOGICAL  :: Found
!
      CPP_ENTER_PROCEDURE( "DECOMPGLOBALTOLOCAL" )

!
! Search over all the PEs
!
      Pe = -1
      Local = 0
      Found = .FALSE.
      DO WHILE ( .NOT. Found .AND. Global .NE. 0 )
!
! Copy the number of entries on each PE
!
        Pe = Pe + 1
        CPP_ASSERT_F90( ( SIZE(Decomp%Head(Pe+1)%StartTags) .EQ. SIZE(Decomp%Head(Pe+1)%EndTags) ) )
!
! Search through the local data segment
!
        Local = 1
        J = 1
        DO WHILE ( .NOT. Found .AND.                                     &
     &             J .LE. SIZE( Decomp%Head(Pe+1)%StartTags ) )
          IF ( Global .GE. Decomp%Head(Pe+1)%StartTags(J) .AND.          &
     &         Global .LE. Decomp%Head(Pe+1)%EndTags(J) ) THEN
            Local = Local +  Global - Decomp%Head(Pe+1)%StartTags(J)
            Found = .TRUE.
          ELSE
            Local = Local + Decomp%Head(Pe+1)%EndTags(J) -               &
     &                      Decomp%Head(Pe+1)%StartTags(J) + 1
          ENDIF
          J = J+1
        ENDDO
!
! Emergency brake
!
        IF ( Pe .EQ. (SIZE( Decomp%Head )-1) .AND. .NOT. Found ) THEN
          Found = .TRUE.
          Local = 0
          Pe    = -1
        ENDIF
      ENDDO

      CPP_ASSERT_F90( Local .LE. Decomp%NumEntries(Pe+1) ) 

      CPP_LEAVE_PROCEDURE( "DECOMPGLOBALTOLOCAL" )
      RETURN
!
!EOC
      END SUBROUTINE DecompGlobalToLocal
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompLocalToGlobal --- Map global index to local and PE
!
! !INTERFACE:
      SUBROUTINE DecompLocalToGlobal ( Decomp, Local, Pe, Global )
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(DecompType), INTENT( IN )   :: Decomp  ! Decomp information
      INTEGER, INTENT( IN )            :: Local   ! Local index
      INTEGER, INTENT( IN )            :: Pe      ! Pe location
!
! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )           :: Global  ! Global index
!
!
! !DESCRIPTION:
!     Given a decomposition and a local-pe index pair, this routine 
!     returns  the 2-D global index. If the local index is not found, 
!     0 is returned. 
!
!     Note that this routine is not efficient by any stretch of the 
!     imagination --- only one index can be converted at a time.
!     In addition, a search procedure must be performed, whose 
!     efficiency is inversely proportional to the size of the 
!     decomposition (in particular, to the number of "runs").  
!     Conceptually this mapping should be used only once in the 
!     program for initialization, and subsequently all calculations 
!     should take place using local indices.
!
! !SYSTEM ROUTINES:
!     SIZE
!
! !REVISION HISTORY:
!   98.03.20   Sawyer     Creation
!
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER  :: J, Counter
      LOGICAL  :: Found
!
      CPP_ENTER_PROCEDURE( "DECOMPLOCALTOGLOBAL" )
      CPP_ASSERT_F90( Pe .GE. 0 )
      CPP_ASSERT_F90( Pe .LT. SIZE(Decomp%Head) )
      CPP_ASSERT_F90( Local .GT. 0 )
      CPP_ASSERT_F90( Local .LE. Decomp%NumEntries(Pe+1) )

      Counter = 0
      Found   = .FALSE.
      J = 0
      DO WHILE ( .NOT. Found )
        J = J+1
        Counter = Counter + Decomp%Head(Pe+1)%EndTags(J) -               &
     &                      Decomp%Head(Pe+1)%StartTags(J) + 1
        IF ( Local .LE.  Counter ) THEN
          Found = .TRUE.
!
! The following calculation is not immediately obvious.  Think about it
!
          Global = Local - Counter + Decomp%Head(Pe+1)%EndTags(J)
          Found = .TRUE.
        ELSEIF ( J .GE. SIZE( Decomp%Head(Pe+1)%StartTags ) ) THEN
!
! Emergency brake
!
          Found = .TRUE.
          Global = 0
        ENDIF
      ENDDO

      CPP_LEAVE_PROCEDURE( "DECOMPLOCALTOGLOBAL" )
      RETURN
!
!EOC
      END SUBROUTINE DecompLocalToGlobal
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!BOP
! !IROUTINE: DecompInfo --- Information about decomposition
!
! !INTERFACE:
      SUBROUTINE DecompInfo( Decomp, Npes, TotalPts )
! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      TYPE(DecompType), INTENT( IN ):: Decomp   ! Decomp information

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )        :: Npes     ! Npes in decomposition
      INTEGER, INTENT( OUT )        :: TotalPts ! Total points in domain
!
!
! !DESCRIPTION:
!     Return information about the decomposition: the number of
!     PEs over which the domain is decomposed, and the size of
!     the domain.
!
! !REVISION HISTORY:
!   00.11.12   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
!
!
      CPP_ENTER_PROCEDURE( "DECOMPINFO" )

      Npes = SIZE( Decomp%Head )
      TotalPts = Decomp%GlobalSize

      CPP_LEAVE_PROCEDURE( "DECOMPINFO" )
      RETURN
!EOC
      END SUBROUTINE DecompInfo
!------------------------------------------------------------------------

      END MODULE decompmodule

