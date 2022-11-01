!*******************************************************************************
!Subroutine - rapid_mus_mat
!*******************************************************************************
#include "LIS_misc.h"
#ifdef PETSc

subroutine rapid_mus_mat

!Purpose:
!Compute RAPID Muskingum operator M, using the equation:
!M=(I-C1*N)^(-1)
!Authors: 
!Charlotte M. Emery, and Cedric H. David, 2018-2020.


!*******************************************************************************
!Fortran includes, modules, and implicity
!*******************************************************************************

#include <petsc/finclude/petscmat.h>
use petscmat

use rapid_var, only :                                                          &
                IS_riv_bas,                                                    &
                JS_riv_bas,JS_riv_bas2,JS_up,                                  & 
                IV_nbup,IV_riv_index,IM_index_up,                              &
                IV_nz,IV_dnz,IV_onz,                                           &
                IS_ownfirst,IS_ownlast,                                        &
                ZM_M,ZV_C1,ZS_threshold,                                       &
                ZS_val,IS_one,                                                 &
                IS_opt_run,                                                    &
                ierr,rank,temp_char,                                           &
                IV_nbrows,IV_lastrow
implicit none


!*******************************************************************************
!Intent (in/out), and local variables 
!*******************************************************************************
PetscInt :: JS_i
PetscInt :: IS_Knilpotent

PetscInt, dimension(:), allocatable :: IV_cols, IV_cols_duplicate
PetscInt, dimension(:), allocatable :: IV_ind, IV_rows
!PetscInt, dimension(:), allocatable :: IV_nbrows

PetscScalar, dimension(:), allocatable :: ZV_cols

Mat :: ZM_MC
VecScatter :: vecscat_all
Vec :: ZV_all


!*******************************************************************************
!Create a temporary VecScatter object
!*******************************************************************************
call VecScatterCreateToAll(ZV_C1,vecscat_all,ZV_all,ierr)
call VecScatterBegin(vecscat_all,ZV_C1,ZV_all,                                 &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd(vecscat_all,ZV_C1,ZV_all,                                   &
                   INSERT_VALUES,SCATTER_FORWARD,ierr)
!The vector ZV_all contains all elements of ZV_C1 and is copied on each
!processor. This way, each processor can access all values of ZV_C1, which is
!needed to count non-zero elements in ZM_MC


!*******************************************************************************
!Create a temporary matrix ZM_MC which is a compressed version of (I-C1*N)^(-1)
!*******************************************************************************
call MatCreate(PETSC_COMM_WORLD,ZM_MC,ierr)
call MatSetSizes(ZM_MC,PETSC_DECIDE,PETSC_DECIDE,IS_riv_bas,IS_riv_bas,ierr)
call MatSetFromOptions(ZM_MC,ierr)
call MatSetUp(ZM_MC,ierr)
!The matrix M=(I-C1*N)^(-1) can be computed as the sum of (C1*N)^(j-1). Each one 
!of these matrices (C1*N)^(j-1) has a maximum of one element per column, and all
!these elements can therefore be stored in a one-dimensional array. Each array
!will be stored in a row of the compressed matrix MC. Note that the eventual 
!respective location of these elements in M can be obtained from N^(j-1). Such
!methodology allows for a fast computation of M.


!*******************************************************************************
!Prepare for matrix preallocation (ZM_MC)
!*******************************************************************************

!-------------------------------------------------------------------------------
!Allocate and initialize temporary variables
!-------------------------------------------------------------------------------
allocate(IV_cols(IS_riv_bas))
allocate(IV_cols_duplicate(IS_riv_bas))
IV_cols(:)=0
IV_cols_duplicate(:)=0
!Used to store the index where each element of MC will be placed in M, IV_cols
!is updated for every power of N.

!allocate(IV_nbrows(IS_riv_bas))
IV_nbrows(:)=1
!Used to store, for each column of ZM_MC, how many non-zeros rows

allocate(ZV_cols(IS_riv_bas))
ZV_cols(:)=1
!Used to store, for a given row of ZM_MC, the element values for each column

IV_nz(:)=0
IV_dnz(:)=0
IV_onz(:)=0
!The number of non-zero elements per row.

!-------------------------------------------------------------------------------
!Populate temporary variables
!-------------------------------------------------------------------------------
do JS_riv_bas2=1,IS_riv_bas 
    do JS_up=1,IV_nbup(IV_riv_index(JS_riv_bas2))
        if (IM_index_up(JS_riv_bas2,JS_up)/=0) then

        JS_riv_bas=IM_index_up(JS_riv_bas2,JS_up)

        IV_cols(JS_riv_bas)=JS_riv_bas2
        IV_cols_duplicate(JS_riv_bas)=JS_riv_bas2

        end if
    end do
end do

!-------------------------------------------------------------------------------
!Count the number of non-zero elements (ZM_MC)
!-------------------------------------------------------------------------------
IV_nz(1)=IS_riv_bas
if ( (1.ge.IS_ownfirst+1).and.(1.lt.IS_ownlast+1) ) then
    IV_dnz(1)=IS_ownlast-IS_ownfirst
    IV_onz(1)=IV_nz(1)-IV_dnz(1)
end if
!The first row

JS_i=2
do while ( COUNT( (IV_cols(1:IS_riv_bas).eq.0) ).ne.IS_riv_bas )

    do JS_riv_bas=1,IS_riv_bas
        if ( IV_cols(JS_riv_bas).ne.0 ) then

            call VecGetValues(ZV_all,IS_one,IV_cols(JS_riv_bas)-1,ZS_val,ierr)

            if ( ABS(ZV_cols(JS_riv_bas)*ZS_val).ge.ZS_threshold ) then

                IV_nz(JS_i) = IV_nz(JS_i)+1
                
                if (((JS_i.ge.IS_ownfirst+1).and.                              &
                     (JS_i.lt.IS_ownlast+1)).and.                              &
                    ((JS_riv_bas.ge.IS_ownfirst+1).and.                        &
                     (JS_riv_bas.lt.IS_ownlast+1))) then 
                    IV_dnz(JS_i) = IV_dnz(JS_i)+1
                endif

                if (((JS_i.ge.IS_ownfirst+1).and.                              &
                     (JS_i.lt.IS_ownlast+1)).and.                              &
                    ((JS_riv_bas.lt.IS_ownfirst+1).or.                         &
                     (JS_riv_bas.ge.IS_ownlast+1))) then 
                    IV_onz(JS_i) = IV_onz(JS_i)+1
                endif

                IV_nbrows(JS_riv_bas) = IV_nbrows(JS_riv_bas)+1
                ZV_cols(JS_riv_bas) = ZV_cols(JS_riv_bas)*ZS_val
                IV_cols(JS_riv_bas) = IV_cols_duplicate(IV_cols(JS_riv_bas))
            else

                IV_cols(JS_riv_bas) = 0        
            endif

        endif
    enddo       

    JS_i=JS_i+1
end do
IS_Knilpotent=JS_i-1
!The other rows

!-------------------------------------------------------------------------------
!Print information on the nilpotence index of M
!-------------------------------------------------------------------------------
write(temp_char,'(i10)') IS_Knilpotent
if (IS_opt_run/=2) then
     call PetscPrintf(PETSC_COMM_WORLD,'Knilpotent='//temp_char//char(10),ierr)
end if


!*******************************************************************************
!Matrix preallocation (ZM_MC)
!*******************************************************************************
call MatSeqAIJSetPreallocation(ZM_MC,PETSC_DEFAULT_INTEGER,IV_nz,ierr)
call MatMPIAIJSetPreallocation(ZM_MC,                                          &
                               PETSC_DEFAULT_INTEGER,                          &
                               IV_dnz(IS_ownfirst+1:IS_ownlast),               &
                               PETSC_DEFAULT_INTEGER,                          &
                               IV_onz(IS_ownfirst+1:IS_ownlast),ierr)


!*******************************************************************************
!Populate matrix (ZM_MC)
!*******************************************************************************

!-------------------------------------------------------------------------------
!Allocate and initialize temporary variables
!-------------------------------------------------------------------------------
allocate(IV_ind(IS_riv_bas))

!-------------------------------------------------------------------------------
!Populate temporary variables
!-------------------------------------------------------------------------------
do JS_riv_bas=1,IS_riv_bas
    IV_cols(JS_riv_bas)=IV_cols_duplicate(JS_riv_bas)
end do
!Reset the value of IV_cols

ZV_cols(:)=1
!Reset the values of ZV_cols

do JS_riv_bas=1,IS_riv_bas
    IV_ind(JS_riv_bas) = JS_riv_bas
end do
!Initialize IV_ind

!-------------------------------------------------------------------------------
!Fill matrix (ZM_MC)
!-------------------------------------------------------------------------------
if (rank==0) then

call MatSetValues(ZM_MC,                                                       &
                  IS_one,0,                                                    &
                  IS_riv_bas,IV_ind(1:IS_riv_bas)-1,                           &
                  ZV_cols(1:IS_riv_bas),                                       &
                  INSERT_VALUES,ierr)
!The first row

JS_i=1
do while ( COUNT( (IV_cols(1:IS_riv_bas).eq.0) ).ne.IS_riv_bas )

    do JS_riv_bas=1,IS_riv_bas
        if (IV_cols(JS_riv_bas).ne.0) then
            call VecGetValues(ZV_all,                                          &
                              IS_one,                                          &
                              IV_cols(JS_riv_bas)-1,                           &
                              ZS_val,ierr)

            if ( ABS(ZV_cols(JS_riv_bas)*ZS_val).ge.ZS_threshold ) then
                ZV_cols(JS_riv_bas) = ZV_cols(JS_riv_bas)*ZS_val
                IV_cols(JS_riv_bas) = IV_cols_duplicate(IV_cols(JS_riv_bas))
            else
                IV_cols(JS_riv_bas) = 0   
                IV_ind(JS_riv_bas) = 0     
            endif
        else
            IV_ind(JS_riv_bas) = 0   
        endif
    enddo

    call MatSetValues(ZM_MC,                                                   &
                      IS_one,JS_i,                                             &
                      IV_nz(JS_i+1),                                           &
                      PACK(IV_ind(1:IS_riv_bas),                               &
                           MASK=IV_ind(1:IS_riv_bas).gt.0)-1,                  &
                      ZV_cols( PACK( IV_ind(1:IS_riv_bas),                     &
                                     MASK=IV_ind(1:IS_riv_bas).gt.0 ) ),       &
                      INSERT_VALUES,ierr)

    JS_i = JS_i+1
    
enddo 
!The other rows

write(temp_char,'(i10)') JS_i    
if (IS_opt_run/=2) call PetscPrintf(PETSC_COMM_WORLD,'Exit at row='            &
                                                     //temp_char//char(10),ierr)

end if

call MatAssemblyBegin(ZM_MC,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(ZM_MC,MAT_FINAL_ASSEMBLY,ierr)


!*******************************************************************************
!Prepare for matrix preallocation (ZM_M)
!*******************************************************************************

!-------------------------------------------------------------------------------
!Allocate and initialize temporary variables
!-------------------------------------------------------------------------------
IV_nz(:)=0
IV_dnz(:)=0
IV_onz(:)=0

!-------------------------------------------------------------------------------
!Count the number of non-zero elements (ZM_MC)
!-------------------------------------------------------------------------------
do JS_riv_bas=1,IS_riv_bas   !Loop over column
    do JS_i=1,IV_nbrows(JS_riv_bas)
   
        if (JS_i.eq.1) then
            JS_riv_bas2 = JS_riv_bas  !row index
        else
            JS_riv_bas2 = IV_cols_duplicate(JS_riv_bas2)
        endif

        IV_nz(JS_riv_bas2) = IV_nz(JS_riv_bas2)+1

        if (((JS_riv_bas2.ge.IS_ownfirst+1).and.                               &
             (JS_riv_bas2.lt.IS_ownlast+1)).and.                               &
            ((JS_riv_bas.ge.IS_ownfirst+1).and.                                &
             (JS_riv_bas.lt.IS_ownlast+1))) then
            IV_dnz(JS_riv_bas2) = IV_dnz(JS_riv_bas2)+1
        endif

        if (((JS_riv_bas2.ge.IS_ownfirst+1).and.                               &
             (JS_riv_bas2.lt.IS_ownlast+1)).and.                               &
            ((JS_riv_bas.lt.IS_ownfirst+1).or.                                 &
             (JS_riv_bas.ge.IS_ownlast+1))) then
            IV_onz(JS_riv_bas2) = IV_onz(JS_riv_bas2)+1
        endif

    enddo
enddo

!*******************************************************************************
!Matrix preallocation (ZM_M)
!*******************************************************************************
call MatSeqAIJSetPreallocation(ZM_M,PETSC_DEFAULT_INTEGER,IV_nz,ierr)
call MatMPIAIJSetPreallocation(ZM_M,                                           &
                               PETSC_DEFAULT_INTEGER,                          &
                               IV_dnz(IS_ownfirst+1:IS_ownlast),               &
                               PETSC_DEFAULT_INTEGER,                          &
                               IV_onz(IS_ownfirst+1:IS_ownlast),ierr)
if (IS_opt_run/=2) call PetscPrintf(PETSC_COMM_WORLD,'Muskingum matrix '       &
                                                //'preallocated'//char(10),ierr)


!*******************************************************************************
!Populate matrix (ZM_M)
!*******************************************************************************

!-------------------------------------------------------------------------------
!Initialize output variable IV_lastrow
!-------------------------------------------------------------------------------

IV_lastrow(:)=0
!Used to store the index of the last/largest row with non-zeros for each column

!-------------------------------------------------------------------------------
!Allocate and initialize temporary variables
!-------------------------------------------------------------------------------
allocate(IV_rows(IS_Knilpotent+1))
IV_rows(:)=0

!-------------------------------------------------------------------------------
!Populate IV_lastrow
!-------------------------------------------------------------------------------

do JS_riv_bas=1,IS_riv_bas
    do JS_i=1,IV_nbrows(JS_riv_bas)  
        if (JS_i.eq.1) then
            IV_lastrow(JS_riv_bas) = JS_riv_bas  !row index
        else
            IV_lastrow(JS_riv_bas)=IV_cols_duplicate(IV_lastrow(JS_riv_bas)) 
        endif
    end do
end do

!-------------------------------------------------------------------------------
!Populate temporary variables
!-------------------------------------------------------------------------------
do JS_riv_bas=1,IS_riv_bas
    IV_cols(JS_riv_bas) = IV_cols_duplicate(JS_riv_bas) 
    IV_ind(JS_riv_bas) = JS_riv_bas
end do

!-------------------------------------------------------------------------------
!Fill ZM_M
!-------------------------------------------------------------------------------
if (rank==0) then

deallocate(ZV_cols)

do JS_riv_bas=1,IS_riv_bas

    do JS_i=1,IV_nbrows(JS_riv_bas)  
        if (JS_i.eq.1) then
            IV_rows(JS_i) = JS_riv_bas  !row index
        else
            IV_rows(JS_i)=IV_cols_duplicate(IV_rows(JS_i-1)) 
        endif
    end do
    allocate(ZV_cols(IV_nbrows(JS_riv_bas)))
    
    call MatGetValues( ZM_MC,                                                  &
                       IV_nbrows(JS_riv_bas),                                  &
                       IV_ind(1:IV_nbrows(JS_riv_bas))-1,                      &
                       IS_one,JS_riv_bas-1,                                    &
                       ZV_cols,ierr )

    call MatSetValues( ZM_M,                                                   &
                       IV_nbrows(JS_riv_bas),                                  &
                       IV_rows(1:IV_nbrows(JS_riv_bas))-1,                     &
                       IS_one, JS_riv_bas-1,                                   &
                       ZV_cols(1:IV_nbrows(JS_riv_bas)),                       &
                       INSERT_VALUES,ierr )

    deallocate(ZV_cols)

end do

end if

call MatAssemblyBegin(ZM_M,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(ZM_M,MAT_FINAL_ASSEMBLY,ierr)
!sparse matrices need be assembled once their elements have been filled

if (IS_opt_run/=2) then
     call PetscPrintf(PETSC_COMM_WORLD,'Muskingum matrix created'              &
                                       //char(10),ierr)
end if


!*******************************************************************************
!Free up memory used by local (temporary) variables
!*******************************************************************************
deallocate(IV_cols)
deallocate(IV_cols_duplicate)
deallocate(IV_ind)
deallocate(IV_rows)
!if (rank==0) then
!    deallocate(IV_nbrows)
!end if

call MatDestroy(ZM_MC,ierr)
call VecScatterDestroy(vecscat_all,ierr)
call VecDestroy(ZV_all,ierr)


!*******************************************************************************
!End subroutine
!*******************************************************************************
if (IS_opt_run/=2) then
     call PetscPrintf(PETSC_COMM_WORLD,'--------------------------'            &
                                       //char(10),ierr)
end if

end subroutine rapid_mus_mat

#else

! Dummy version
subroutine rapid_mus_mat
  use LIS_logmod, only: LIS_logunit, LIS_endrun
  implicit none
  write(LIS_logunit,*)'[ERR] RAPID called w/o PETSc support!'
  write(LIS_logunit,*)'[ERR] Recompile with PETSc and try again!'
  call LIS_endrun()
end subroutine rapid_mus_mat

#endif
