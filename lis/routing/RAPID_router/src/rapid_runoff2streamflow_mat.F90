!*******************************************************************************
!Subroutine - rapid_runoff2streamflow_mat
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_runoff2streamflow_mat

!Purpose:
!Compute operator that turns runoff into streamflow over a day of run
!This operator is used within the data assimilation as part of the
!observations operator H in the Kalman gain
!Authors: 
!Charlotte M. Emery, and Cedric H. David, 2018-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************

#include <petsc/finclude/petscmat.h>
use petscmat

use rapid_var, only :                                                          &
                IS_riv_bas,JS_riv_bas,JS_riv_bas2,JS_up,                       &
                IM_index_up,IV_riv_index,IV_nbup,                              &
                ZV_C1,ZV_C2,ZV_C3,                                             &
                ZM_M,                                                          &
                ierr,rank,                                                     &
                ZV_temp1,                                                      &
                ZS_one,IS_one,ZS_val,                                          &
                IV_nz,IV_dnz,IV_onz,                                           &
                IS_ownfirst,IS_ownlast,                                        &
                IS_R,IS_RpM,                                                   &
    !New variables to add in rapid_var+rapid_mus_mat when linking the code
                ZM_L,IV_lastrow,IV_nbrows
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************

PetscInt :: JS_i, JS_j, JS_p
PetscInt :: JS_next_last
PetscInt :: found_nz, nb_intersection

PetscScalar :: ZS_val2, ZS_factor

PetscInt, dimension(:), allocatable :: IV_index_down
PetscInt, dimension(:), allocatable :: IV_ind_rows
PetscInt, dimension(:), allocatable :: IV_lastrow_save

VecScatter :: vecscat_C1_zeroth, vecscat_C2_zeroth, vecscat_C3_zeroth
Vec :: ZV_C1_zeroth, ZV_C2_zeroth, ZV_C3_zeroth

Mat :: ZM_temp, ZM_temp2
Mat :: ZM_temp_pow

!*******************************************************************************
!Prepare for preallocation of ZM_L and associated temporary matrices
!*******************************************************************************

!-------------------------------------------------------------------------------
!Allocate and initialize temporary variables
!-------------------------------------------------------------------------------

allocate(IV_index_down(IS_riv_bas))
IV_index_down(:) = 0
!Indexes of downstream reaches

allocate(IV_lastrow_save(IS_riv_bas))
!backup of IV_lastrow

!No re-initialization of IV_nz/IV_dnz/IV_onz
!Contains nonzeros elements of (I-C1*N)^(-1)
!Used as a starting point to count additional non-zeros

!-------------------------------------------------------------------------------
!Populate temporary variables
!-------------------------------------------------------------------------------

do JS_riv_bas2=1,IS_riv_bas 
    do JS_up=1,IV_nbup(IV_riv_index(JS_riv_bas2))
        if (IM_index_up(JS_riv_bas2,JS_up)/=0) then

        JS_riv_bas=IM_index_up(JS_riv_bas2,JS_up)

        IV_index_down(JS_riv_bas)=JS_riv_bas2

        end if
    end do
end do
!IV_index_down gives, for each reach, the index of the unique downstream reach

do JS_riv_bas=1,IS_riv_bas
    IV_lastrow_save(JS_riv_bas)=IV_lastrow(JS_riv_bas)
end do
!IV_lastrow is updated when counting nonzeros for
!preallocation, but
!original IV_lastrow is needed after when
!populating the matrices

!-------------------------------------------------------------------------------
!Create temporary matrices
!-------------------------------------------------------------------------------

call MatCreate(PETSC_COMM_WORLD,ZM_temp,ierr)
call MatSetSizes(ZM_temp,PETSC_DECIDE,PETSC_DECIDE,IS_riv_bas,IS_riv_bas,ierr)
call MatSetFromOptions(ZM_temp,ierr)
call MatSetUp(ZM_temp,ierr)
!Fixed matrix equal to (I-C1*N)^(-1)*(C3+C2*N)

call MatCreate(PETSC_COMM_WORLD,ZM_temp_pow,ierr)
call MatSetSizes(ZM_temp_pow,PETSC_DECIDE,PETSC_DECIDE,IS_riv_bas,IS_riv_bas,ierr)
call MatSetFromOptions(ZM_temp_pow,ierr)
call MatSetUp(ZM_temp_pow,ierr)
!Temporary matrix as successive power of ZM_temp=(I-C1*N)^(-1)*(C3+C2*N)


!-------------------------------------------------------------------------------
!Count number of non-zero elements of (I-C1*N)^(-1)*(C3+C2*N)
!-------------------------------------------------------------------------------

!IV_nz/IV_dnz/IV_onz already
!contains number of nonzeros elements of (I-C1*N)^(-1) and depends on ZS_threshold
!Used as a starting point to count additional non-zeros
!Non-zeros pattern of matrix (I-C1*N)^(-1) is a subset 
!of non-zeros pattern of (I-C1*N)^(-1)*(C3+C2*N)
!For each column of (I-C1*N)^(-1), we count the number of new non-zeros 
!added when multiplied by (C3+C2*N)

allocate(IV_ind_rows(2))
!Temporary array
!Contains indices of non-zeros elements of current column in (C3+C2*N)

do JS_riv_bas = 1,IS_riv_bas-1 !loop over column (no new nz on last col)
    
    !Test: possible new non-zeros if::
    ! - there is a reach downstream the current reach JS_riv_bas
    ! - there is a reach downstream the last non-zero reach of column JS_riv_bas
    if ((IV_index_down(JS_riv_bas).ne.0).and.                   &
        (IV_index_down(IV_lastrow(JS_riv_bas)).ne.0)) then

        !In the matrix (I-C1*N)^(-1)*(C3+C2*N),
        !the element at location (JS_next_last,JS_riv_bas)~=0 if:
        ! 1) there are non zeros in column JS_riv_bas of (C3+C2*N)
        ! 2) there are non-zeros in row JS_next_last of (I-C1*N)^(-1)
        ! 3) intersection between non-zeros indexes of (1) and (2) is not empty

        ! (1) ::
        !In column JS_riv_bas of (C3+C2*N)
        !List of non-zeros elements - at most 2
        IV_ind_rows(1)=JS_riv_bas
        IV_ind_rows(2)=IV_index_down(JS_riv_bas)

        ! (2)
        !In column JS_riv_bas of (I-C1*N)^(-1)
        !position of the potential next non-zeros in the column
        JS_next_last=IV_index_down(IV_lastrow(JS_riv_bas))
        !JS_next_last is the row in (I-C1*N)^(-1) to look for non-zeros

        found_nz = 1
        do while ((found_nz.eq.1).and.(JS_next_last.ne.0))

            ! (3)
            ! Check intersection
            nb_intersection=0
            do JS_i=1,2 !loop over number of element from (1)

                if (IV_lastrow(IV_ind_rows(JS_i)).eq.JS_next_last) then
                    nb_intersection=nb_intersection+1
                end if

                if (IV_lastrow(IV_ind_rows(JS_i)).gt.JS_next_last) then
                    JS_j=IV_ind_rows(JS_i)
                    do while (JS_j.lt.JS_next_last)
                        JS_j=IV_index_down(JS_j)
                    end do
                    if (JS_j.eq.JS_next_last) then
                        nb_intersection=nb_intersection+1
                    end if
                end if

            end do

            if (nb_intersection.ne.0) then

                IV_nz(JS_next_last) = IV_nz(JS_next_last)+1

                if (((JS_next_last.ge.IS_ownfirst+1).and.                      &
                     (JS_next_last.lt.IS_ownlast+1)).and.                      &
                    ((JS_riv_bas.ge.IS_ownfirst+1).and.                        &
                     (JS_riv_bas.lt.IS_ownlast+1))) then
                    IV_dnz(JS_next_last) = IV_dnz(JS_next_last)+1
                end if

                if (((JS_next_last.ge.IS_ownfirst+1).and.                      &
                     (JS_next_last.lt.IS_ownlast+1)).and.                      &
                    ((JS_riv_bas.lt.IS_ownfirst+1).or.                         &
                     (JS_riv_bas.ge.IS_ownlast+1))) then
                    IV_onz(JS_next_last) = IV_onz(JS_next_last)+1
                end if

                IV_lastrow(JS_riv_bas) = JS_next_last
                JS_next_last=IV_index_down(IV_lastrow(JS_riv_bas))

            else

                found_nz=0

            end if

        end do
    end if
end do
deallocate(IV_ind_rows)


!-------------------------------------------------------------------------------
!Matrix preallocation of temporary matrix ZM_temp=(I-C1*N)^(-1)*(C3+C2*N)
!-------------------------------------------------------------------------------

call MatSeqAIJSetPreallocation(ZM_temp,PETSC_DEFAULT_INTEGER,IV_nz,ierr)
call MatMPIAIJSetPreallocation(ZM_temp,                                        &
                               PETSC_DEFAULT_INTEGER,                          &
                               IV_dnz(IS_ownfirst+1:IS_ownlast),               &
                               PETSC_DEFAULT_INTEGER,                          &
                               IV_onz(IS_ownfirst+1:IS_ownlast),ierr)

!-------------------------------------------------------------------------------
!Count number of additionnal non-zero elements of (I-N)^(-1) (for ZM_L only)
!-------------------------------------------------------------------------------

do while ( COUNT(IV_lastrow(1:IS_riv_bas).eq.0).ne.IS_riv_bas )

    do JS_riv_bas=1,IS_riv_bas   !loop over column

        JS_riv_bas2 = IV_index_down(IV_lastrow(JS_riv_bas))

        if (JS_riv_bas2.ne.0) then

            IV_nz(JS_riv_bas2) = IV_nz(JS_riv_bas2)+1

            if (((JS_riv_bas2.ge.IS_ownfirst+1).and.                           &
                 (JS_riv_bas2.lt.IS_ownlast+1)).and.                           &
                ((JS_riv_bas.ge.IS_ownfirst+1).and.                            &
                 (JS_riv_bas.lt.IS_ownlast+1))) then
   
                IV_dnz(JS_riv_bas2) = IV_dnz(JS_riv_bas2)+1

            end if

            if (((JS_riv_bas2.ge.IS_ownfirst+1).and.                           &
                 (JS_riv_bas2.lt.IS_ownlast+1)).and.                           &
                ((JS_riv_bas.lt.IS_ownfirst+1).or.                             &
                 (JS_riv_bas.ge.IS_ownlast+1))) then
   
                IV_onz(JS_riv_bas2) = IV_onz(JS_riv_bas2)+1

            end if

        endif

        IV_lastrow(JS_riv_bas) = JS_riv_bas2

    end do

end do

!*******************************************************************************
!Matrix preallocation (ZM_L)
!*******************************************************************************

call MatSeqAIJSetPreallocation(ZM_L,PETSC_DEFAULT_INTEGER,IV_nz,ierr)
call MatMPIAIJSetPreallocation(ZM_L,                                           &
                               PETSC_DEFAULT_INTEGER,                          &
                               IV_dnz(IS_ownfirst+1:IS_ownlast),               &
                               PETSC_DEFAULT_INTEGER,                          &
                               IV_onz(IS_ownfirst+1:IS_ownlast),ierr)

call PetscSynchronizedPrintf(PETSC_COMM_WORLD,'ZM_L pre-allocated'//char(10),ierr)


!*******************************************************************************
!Populate ZM_L and associated matrices
!*******************************************************************************

!-------------------------------------------------------------------------------
!Create a temporary VecScatter object
!-------------------------------------------------------------------------------

call VecScatterCreateToZero(ZV_C1,vecscat_C1_zeroth,ZV_C1_zeroth,ierr)
call VecScatterCreateToZero(ZV_C2,vecscat_C2_zeroth,ZV_C2_zeroth,ierr)
call VecScatterCreateToZero(ZV_C3,vecscat_C3_zeroth,ZV_C3_zeroth,ierr)
!The vector ZV_zeroth contains all elements of ZV_C2 or ZV_C3
!and is copied on the zeroth processor only
!This way, the zeroth processor can access all values of ZV_C2 or ZV_C3 
!which is needed to populate the matrices ZM_L/ZM_temp

call VecScatterBegin(vecscat_C1_zeroth,ZV_C1,ZV_C1_zeroth,                       &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat_C1_zeroth,ZV_C1,ZV_C1_zeroth,                         &
                   INSERT_VALUES,SCATTER_FORWARD,ierr)

call VecScatterBegin(vecscat_C2_zeroth,ZV_C2,ZV_C2_zeroth,                       &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat_C2_zeroth,ZV_C2,ZV_C2_zeroth,                         &
                   INSERT_VALUES,SCATTER_FORWARD,ierr)

call VecScatterBegin(vecscat_C3_zeroth,ZV_C3,ZV_C3_zeroth,                       &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat_C3_zeroth,ZV_C3,ZV_C3_zeroth,                         &
                   INSERT_VALUES,SCATTER_FORWARD,ierr) 

!-------------------------------------------------------------------------------
!Re-initialize temporary variables
!-------------------------------------------------------------------------------

do JS_riv_bas=1,IS_riv_bas
    IV_lastrow(JS_riv_bas)=IV_lastrow_save(JS_riv_bas)
end do

!-------------------------------------------------------------------------------
!Poupulate ZM_temp/Initialize ZM_L as (I-C1*N)^(-1)*(C3+C2*N)
!-------------------------------------------------------------------------------

if (rank==0) then

do JS_riv_bas=1,IS_riv_bas !loop over column

    !For the current colum
    !IV_ind_rows contains index of non-zero rows of (I-C1*N)^(-1)
    allocate(IV_ind_rows(IV_nbrows(JS_riv_bas))) 
    do JS_j=1,IV_nbrows(JS_riv_bas)
        if (JS_j.eq.1) then
            IV_ind_rows(JS_j)=JS_riv_bas-1
            JS_riv_bas2=IV_index_down(JS_riv_bas)
        else
            IV_ind_rows(JS_j)=JS_riv_bas2-1
            JS_riv_bas2=IV_index_down(JS_riv_bas2)
        end if
    end do

    ZS_val=0.0
    ZS_val2=0.0

    ! (1)
    !Set ZM_L/ZM_temp diagonal element in current column
    ! 

    !ZS_val=ZV_C3(JS_riv_bas)
    call VecGetValues(ZV_C3_zeroth,                                            &
                      IS_one,                                                  &
                      JS_riv_bas-1,                                            &
                      ZS_val,ierr)

    !ZM_L(JS_riv_bas,JS_riv_bas) = ZV_C3(JS_riv_bas)
    call MatSetValues(ZM_L,                                                    &
                      IS_one,JS_riv_bas-1,                                     &
                      IS_one,JS_riv_bas-1,                                     &
                      ZS_val,INSERT_VALUES,                                    &
                      ierr)   

    call MatSetValues(ZM_temp,                                                 &
                      IS_one,JS_riv_bas-1,                                     &
                      IS_one,JS_riv_bas-1,                                     &
                      ZS_val,INSERT_VALUES,                                    &
                      ierr) 

    if ((IV_index_down(JS_riv_bas).ne.0).and.(IV_nbrows(JS_riv_bas).gt.1)) then

        ! (2)
        !Set ZM_L/ZM_temp first element downstream diagonal in current column
        !

        !ZS_val2=ZV_C1( down(JS_riv_bas) )
        call VecGetValues(ZV_C1_zeroth,                                        &
                          IS_one,                                              &
                          IV_index_down(JS_riv_bas)-1,                         &
                          ZS_val2,ierr)

        !ZS_val = ZS_val*ZV_C1( down(JS_riv_bas) )
        ZS_val = ZS_val*ZS_val2

        !ZS_val2=ZV_C2( down(JS_riv_bas) )
        call VecGetValues(ZV_C2_zeroth,                                        &
                          IS_one,                                              &
                          IV_index_down(JS_riv_bas)-1,                         &
                          ZS_val2,ierr)

        !ZS_val = ZS_val + ZV_C2( down(JS_riv_bas) )
        !ZS_val = ZV_C3(JS_riv_bas)*ZV_C1( down(JS_riv_bas) ) + ZV_C2( down(JS_riv_bas) )
        ZS_val = ZS_val + ZS_val2
        
        !ZM_L(down(JS_riv_bas),JS_riv_bas) = ZS_val
        call MatSetValues(ZM_L,                                                &
                          IS_one,IV_index_down(JS_riv_bas)-1,                  &
                          IS_one,JS_riv_bas-1,                                 &
                          ZS_val,INSERT_VALUES,                                &
                          ierr)   

        call MatSetValues(ZM_temp,                                             &
                          IS_one,IV_index_down(JS_riv_bas)-1,                  &
                          IS_one,JS_riv_bas-1,                                 &
                          ZS_val,INSERT_VALUES,                                &
                          ierr) 

        ! (3)
        !Set other non-zeros elements in column JS_riv_bas
        !over the non-zero mask of (I-C1*N)^(-1)
        !

        do JS_j=3,IV_nbrows(JS_riv_bas)

            !ZS_val2=ZV_C1( IV_ind_rows(JS_j) )
            call VecGetValues(ZV_C1_zeroth,                                    &
                              IS_one,                                          &
                              IV_ind_rows(JS_j),                               &
                              ZS_val2,ierr)

            !ZS_val = ZS_val*ZS_val2
            ZS_val = ZS_val*ZS_val2

            !ZM_L(IV_ind_rows(JS_j),JS_riv_bas) = ZV_C1(IV_ind_rows(JS_j))* 
            !                                     ZM_L(IV_ind_rows(JS_j-1),JS_riv_bas)
            call MatSetValues(ZM_L,                                            &
                              IS_one,IV_ind_rows(JS_j),                        &
                              IS_one,JS_riv_bas-1,                             &
                              ZS_val,INSERT_VALUES,                            &
                              ierr)   

            call MatSetValues(ZM_temp,                                         &
                              IS_one,IV_ind_rows(JS_j),                        &
                              IS_one,JS_riv_bas-1,                             &
                              ZS_val,INSERT_VALUES,                            &
                              ierr)

        end do

    end if

    deallocate(IV_ind_rows)

    ! (4) ::
    ! Set new non-zeros elements due to the product with (C3+C2*N)
    !

    allocate(IV_ind_rows(2))
    IV_ind_rows(1)=JS_riv_bas
    IV_ind_rows(2)=IV_index_down(JS_riv_bas)

    found_nz = 1
    JS_next_last=IV_index_down(IV_lastrow(JS_riv_bas))
    do while ((found_nz.eq.1).and.(JS_next_last.ne.0))

        !JS_next_last is the row in (I-C1*N)^(-1) to look for non-zeros
        JS_next_last=IV_index_down(IV_lastrow(JS_riv_bas))

        nb_intersection=0
        do JS_i=1,2 !loop over number of element from (1)

            if (IV_lastrow(IV_ind_rows(JS_i)).eq.JS_next_last) then
                nb_intersection=nb_intersection+1
            end if

            if (IV_lastrow(IV_ind_rows(JS_i)).gt.JS_next_last) then
                JS_j=IV_ind_rows(JS_i)
                do while (JS_j.lt.JS_next_last)
                    JS_j=IV_index_down(JS_j)
                end do
                if (JS_j.eq.JS_next_last) then
                    nb_intersection=nb_intersection+1
                end if
            end if

        end do

        if (nb_intersection.ne.0) then

            if (JS_next_last.eq.IV_index_down(JS_riv_bas)) then

                !ZS_val2=ZV_C1( down(JS_riv_bas) )
                call VecGetValues(ZV_C1_zeroth,                                &
                                  IS_one,                                      &
                                  IV_index_down(JS_riv_bas)-1,                 &
                                  ZS_val2,ierr)

                !ZS_val = ZS_val*ZV_C1( down(JS_riv_bas) )
                ZS_val = ZS_val*ZS_val2

                !ZS_val2=ZV_C2( down(JS_riv_bas) )
                call VecGetValues(ZV_C2_zeroth,                                &
                                  IS_one,                                      &
                                  IV_index_down(JS_riv_bas)-1,                 &
                                  ZS_val2,ierr)

               !ZS_val = ZS_val + ZV_C2( down(JS_riv_bas) )
               !ZS_val = ZV_C3(JS_riv_bas)*ZV_C1( down(JS_riv_bas) ) + ZV_C2( down(JS_riv_bas) )
               ZS_val = ZS_val + ZS_val2              

            else

                call VecGetValues(ZV_C1_zeroth,                                &
                                  IS_one,                                      &
                                  JS_next_last-1,                              &
                                  ZS_val2,ierr)

                ZS_val = ZS_val*ZS_val2

            end if

            call MatSetValues(ZM_L,                                        &
                              IS_one,JS_next_last-1,                       &
                              IS_one,JS_riv_bas-1,                         &
                              ZS_val,INSERT_VALUES,                        &
                              ierr)   

            call MatSetValues(ZM_temp,                                     &
                              IS_one,JS_next_last-1,                       &
                              IS_one,JS_riv_bas-1,                         &
                              ZS_val,INSERT_VALUES,                        &
                              ierr)
 
            IV_lastrow(JS_riv_bas) = JS_next_last
            JS_next_last=IV_index_down(IV_lastrow(JS_riv_bas))
            
        else
            found_nz=0

        end if

    end do

    deallocate(IV_ind_rows) 

end do

end if

!-------------------------------------------------------------------------------
!Matrices assembly: ZM_temp
!-------------------------------------------------------------------------------

call MatAssemblyBegin(ZM_temp,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(ZM_temp,MAT_FINAL_ASSEMBLY,ierr)

!-------------------------------------------------------------------------------
!Poupulate the rest of ZM_L with 0s
!-------------------------------------------------------------------------------

if (rank.eq.0) then

do while ( COUNT(IV_lastrow(1:IS_riv_bas).eq.0).ne.IS_riv_bas )

    do JS_riv_bas=1,IS_riv_bas   !loop over column

        JS_riv_bas2 = IV_index_down(IV_lastrow(JS_riv_bas))

        call MatSetValues(ZM_L,                                                &
                          IS_one,                                              &
                          JS_riv_bas2-1,                                       &
                          IS_one,                                              &
                          JS_riv_bas-1,                                        &
                          0*ZS_one,INSERT_VALUES,                              &
                          ierr)

        IV_lastrow(JS_riv_bas) = JS_riv_bas2

    end do

end do

end if

!-------------------------------------------------------------------------------
!Matrices assembly: ZM_L
!-------------------------------------------------------------------------------

call MatAssemblyBegin(ZM_L,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(ZM_L,MAT_FINAL_ASSEMBLY,ierr)

call PetscSynchronizedPrintf(PETSC_COMM_WORLD,'ZM_L assembled with 0s'//char(10),ierr)

!*******************************************************************************
!Build ZM_L
!*******************************************************************************

!-------------------------------------------------------------------------------
!Initialize ZM_temp_pow as ZM_L (values+non-zeros pattern)
!-------------------------------------------------------------------------------

call MatDuplicate(ZM_L,MAT_COPY_VALUES,ZM_temp_pow,ierr)
!ZM_temp_pow initialize as ZM_L = (I-C1*N)^(-1)*(C3+C2*N)
!Same non-zero pattern as ZM_L (necessary to use MatMatMult)

!-------------------------------------------------------------------------------
!Update ZM_L to ((IS_RpM*IS_R-1)/(IS_RpM*IS_R))*ZM_L
!-------------------------------------------------------------------------------

ZS_factor = (REAL(IS_RpM*IS_R-1)/REAL(IS_RpM*IS_R))
call MatScale(ZM_L,ZS_factor,ierr)

!-------------------------------------------------------------------------------
!Loop over (IS_RpM*IS_R)-1 steps to add all terms in ZM_L
!-------------------------------------------------------------------------------

do JS_p=2,(IS_RpM*IS_R-1)

    call MatMatMult(ZM_temp_pow,                                               &
                    ZM_temp,                                                   &
                    MAT_INITIAL_MATRIX,                                        &
                    PETSC_DEFAULT_REAL,                                        &
                    ZM_temp2,ierr)

    call MatCopy(ZM_temp2,ZM_temp_pow,SAME_NONZERO_PATTERN,ierr)
    call MatDestroy(ZM_temp2,ierr)

    ZS_factor = (REAL(IS_RpM*IS_R-JS_p)/REAL(IS_RpM*IS_R))

    call MatAXPY(ZM_L,                                                         &
                 ZS_factor,                                                    &
                 ZM_temp_pow,                                                  &
                 SAME_NONZERO_PATTERN,                                         &
                 ierr)
                 

end do

!-------------------------------------------------------------------------------
!Add last term of the sum in ZL_L (power=0 <=> add identity)
!-------------------------------------------------------------------------------

call MatShift(ZM_L,ZS_one,ierr)

!-------------------------------------------------------------------------------
!Free memory
!-------------------------------------------------------------------------------

call MatDestroy(ZM_temp,ierr)
call MatDestroy(ZM_temp_pow,ierr)

!-------------------------------------------------------------------------------
!Compute ZM_L = ZM_L*(I-C1*N)^(-1)*(C1+C2)
!-------------------------------------------------------------------------------

call MatMatMult(ZM_L,ZM_M,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,ZM_temp,ierr)
!ZM_temp = ZL_L*ZL_M

call VecCopy(ZV_C1,ZV_temp1,ierr)
call VecAXPY(ZV_temp1,ZS_one,ZV_C2,ierr)
!ZV_temp1 = ZV_C1+ZV_C2

call MatDiagonalScale(ZM_temp,PETSC_NULL_VEC,ZV_temp1,ierr)
!ZM_temp = !ZM_temp*diag(C1+C2)

call MatCopy(ZM_temp,ZM_L,DIFFERENT_NONZERO_PATTERN,ierr)
!ZM_L = ZM_temp


call PetscPrintf(PETSC_COMM_WORLD,'Inverse runoff-discharge operator created'    &
                                   //char(10),ierr)


!*******************************************************************************
!Free up memory used by local (temporary) variables
!*******************************************************************************

deallocate(IV_index_down)
deallocate(IV_lastrow_save)

call MatDestroy(ZM_temp,ierr)
call VecScatterDestroy(vecscat_C1_zeroth,ierr)
call VecScatterDestroy(vecscat_C2_zeroth,ierr)
call VecScatterDestroy(vecscat_C2_zeroth,ierr)
call VecDestroy(ZV_C1_zeroth,ierr)
call VecDestroy(ZV_C2_zeroth,ierr)
call VecDestroy(ZV_C3_zeroth,ierr)


!*******************************************************************************
!End subroutine 
!*******************************************************************************
end subroutine rapid_runoff2streamflow_mat

#else

! Dummy version
subroutine rapid_runoff2streamflow_mat
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_runoff2streamflow_mat

#endif
