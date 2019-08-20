! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module matrixOper_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: nRHS           ! number of unknown variables on the RHS of the linear system A.X=B
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

implicit none
private
public::lapackSolv
public::scaleMatrices
public::computeGradient
contains

 ! **********************************************************************************************************
 ! public subroutine: scaleMatrices: scale the matrices 
 ! **********************************************************************************************************
 subroutine scaleMatrices(ixMatrix,nState,aJac,fScale,xScale,aJacScaled,err,message)
 implicit none
 ! input variables
 integer(i4b),intent(in)         :: ixMatrix          ! type of matrix (full Jacobian or band diagonal)
 integer(i4b),intent(in)         :: nState            ! number of state variables 
 real(dp),intent(in)             :: aJac(:,:)         ! original Jacobian matrix
 real(dp),intent(in)             :: fScale(:)         ! function scaling vector
 real(dp),intent(in)             :: xScale(:)         ! "variable" scaling vector, i.e., for state variables
 ! output variables
 real(dp),intent(out)            :: aJacScaled(:,:)   ! scaled Jacobian matrix
 integer(i4b),intent(out)        :: err               ! error code
 character(*),intent(out)        :: message           ! error message
 ! ---------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: iState            ! row index
 integer(i4b)                    :: jState            ! comumn index
 integer(i4b)                    :: kState            ! band diagonal index
 ! ---------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='scaleMatrices/'

 ! select the type of matrix
 select case(ixMatrix)

  ! * full matrix
  case(ixFullMatrix)

   ! scale by both the scaling factors for the function (fScale) and variable (xScale)
   do iState=1,nState
    do jState=1,nState
     aJacScaled(iState,jState) = fScale(iState)*aJac(iState,jState)*xScale(jState)
    end do
   end do

  ! * band-diagonal matrix
  case(ixBandMatrix)

   ! initialize the matrix to zero (some un-used elements)
   aJacScaled(:,:) = 0._dp

   ! scale the rows by the function scaling factor and the colmns by the variable scaling factor
   do jState=1,nState       ! (loop through model state variables)
    do iState=max(1,jState-ku),min(nState,jState+kl)
     kState = kl+ku+1+iState-jState
     aJacScaled(kState,jState) = fScale(iState)*aJac(kState,jState)*xScale(jState)
    end do
   end do  ! looping through state variables

  ! check that we found a valid option (should not get here because of the check above; included for completeness)
  case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'

 end select  ! (option to solve the linear system A.X=B)

 end subroutine scaleMatrices


 ! *********************************************************************************************************
 ! * private subroutine computeGradient: compute the gradient of the function
 ! *********************************************************************************************************
 subroutine computeGradient(ixMatrix,nState,aJac,rVec,grad,err,message)
 implicit none
 ! input
 integer(i4b),intent(in)        :: ixMatrix   ! type of matrix (full Jacobian or band diagonal)
 integer(i4b),intent(in)        :: nState     ! number of state variables
 real(dp),intent(in)            :: aJac(:,:)  ! jacobian matrix
 real(dp),intent(in)            :: rVec(:)    ! residual vector
 ! output
 real(dp),intent(out)           :: grad(:)    ! gradient
 integer(i4b),intent(out)       :: err        ! error code
 character(*),intent(out)       :: message    ! error message
 ! local
 integer(i4b)                   :: iJac       ! index of model state variable
 integer(i4b)                   :: iState     ! index of the residual vector
 ! initialize error control
 err=0; message='computeGradient/'

 ! check if full Jacobian or band-diagonal matrix
 select case(ixMatrix)

  ! full Jacobian matrix
  case(ixFullMatrix)

   ! compute the gradient
   grad = matmul(rVec,aJac)         ! gradient

  ! band-diagonal matrix
  case(ixBandMatrix)

   ! compute the gradient
   grad(:) = 0._dp
   do iJac=1,nState  ! (loop through state variables)
    do iState=max(1,iJac-ku),min(nState,iJac+kl)
     grad(iJac) = grad(iJac) + aJac(kl+ku+1+iState-iJac,iJac)*rVec(iState)
    end do
   end do

  ! check that we found a valid option
  case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'

 end select  ! (option to solve the linear system A.X=B)

 end subroutine computeGradient


 ! *********************************************************************************************************
 ! public subroutine lapackSolv: use the lapack routines to solve the linear system A.X=B
 ! *********************************************************************************************************
 subroutine lapackSolv(ixMatrix,nState,aJac,rVec,xInc,err,message)
 implicit none
 ! dummy
 integer(i4b),intent(in)        :: ixMatrix      ! type of matrix (full Jacobian or band diagonal)
 integer(i4b),intent(in)        :: nState        ! number of state variables
 real(dp),intent(inout)         :: aJac(:,:)     ! input = the Jacobian matrix A; output = decomposed matrix
 real(dp),intent(in)            :: rVec(:)       ! the residual vector B
 real(dp),intent(out)           :: xInc(:)       ! the solution vector X
 integer(i4b),intent(out)       :: err           ! error code
 character(*),intent(out)       :: message       ! error message
 ! local
 real(dp)                       :: rhs(nState,1) ! the nState-by-nRHS matrix of matrix B, for the linear system A.X=B
 integer(i4b)                   :: iPiv(nState)  ! defines if row i of the matrix was interchanged with row iPiv(i)
 ! initialize error control
 select case(ixMatrix)
  case(ixFullMatrix); message='lapackSolv/dgesv/'
  case(ixBandMatrix); message='lapackSolv/dgbsv/'
  case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'
 end select

 ! form the rhs matrix
 ! NOTE: copy the vector here to ensure that the residual vector is not overwritten
 rhs(:,1) = rVec(:)

 ! identify option to solve the linear system A.X=B
 select case(ixMatrix)

  ! lapack: use the full Jacobian matrix to solve the linear system A.X=B
  case(ixFullMatrix)
   call dgesv(nState,    &  ! intent(in):    [i4b]               number of state variables
              nRHS,      &  ! intent(in):    [i4b]               number of columns of the matrix B
              aJac,      &  ! intent(inout): [dp(nState,nState)] input = the nState-by-nState Jacobian matrix A; output = decomposed matrix
              nState,    &  ! intent(in):    [i4b]               the leading dimension of aJac
              iPiv,      &  ! intent(out):   [i4b(nState)]       defines if row i of the matrix was interchanged with row iPiv(i)
              rhs,       &  ! intent(inout): [dp(nState,nRHS)]   input = the nState-by-nRHS matrix of matrix B; output: the solution matrix X
              nState,    &  ! intent(in):    [i4b]               the leading dimension of matrix rhs
              err)          ! intent(out)    [i4b]               error code

  ! lapack: use the band diagonal matrix to solve the linear system A.X=B
  case(ixBandMatrix)
   call dgbsv(nState,    &  ! intent(in):    [i4b]               number of state variables
              kl,        &  ! intent(in):    [i4b]               number of subdiagonals within the band of A
              ku,        &  ! intent(in):    [i4b]               number of superdiagonals within the band of A
              nRHS,      &  ! intent(in):    [i4b]               number of columns of the matrix B
              aJac,      &  ! intent(inout): [dp(nBands,nState)] input = the nBands-by-nState Jacobian matrix A; output = decomposed matrix
              nBands,    &  ! intent(in):    [i4b]               the leading dimension of aJac
              iPiv,      &  ! intent(out):   [i4b(nState)]       defines if row i of the matrix was interchanged with row iPiv(i)
              rhs,       &  ! intent(inout): [dp(nState,nRHS)]   input = the nState-by-nRHS matrix of matrix B; output: the solution matrix X
              nState,    &  ! intent(in):    [i4b]               the leading dimension of matrix rhs
              err)          ! intent(out)    [i4b]               error code

  ! check that we found a valid option (should not get here because of the check above; included for completeness)
  case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'

 end select  ! (option to solve the linear system A.X=B)

 ! identify any errors
 ! NOTE: return negative error code to force a time step reduction and another trial
 if(err/=0)then
  if(err<0)then
   write(message,'(a,i0,a)') trim(message)//'the ',err,'-th argument had an illegal value'
   err=-20; return
  else
   write(message,'(a,i0,a,i0,a)') trim(message)//'U(',err,',',err,') is exactly zero - factorization complete, but U is singular so the solution could not be completed'
   err=-20; return
  end if
 end if

 ! extract the iteration increment
 xInc(1:nState) = rhs(1:nState,1)

 end subroutine lapackSolv



end module matrixOper_module
