! 07/12/2019 Shugong Wang for Jules 5.5 for disabling MPI calls in JULES code
! The code is based on the dummy MPI code of JULES
MODULE lis_mpi_dummy
  IMPLICIT NONE
  INTEGER, PARAMETER :: mpi_comm_world = 1
  INTEGER, PARAMETER :: mpi_real = 1
  INTEGER, PARAMETER :: mpi_integer = 1
  INTEGER, PARAMETER :: mpi_logical = 1
  INTEGER, PARAMETER :: mpi_info_null = 1
  INTEGER, PARAMETER :: mpi_land = 1
  INTEGER, PARAMETER :: mpi_address_kind = 1
  interface mpi_bcast
    procedure mpi_bcast1, mpi_bcast2, mpi_bcast3, mpi_bcast4, mpi_bcast5
  end interface
  interface mpi_type_create_resized
    procedure mpi_type_create_resized1, mpi_type_create_resized2
  end interface 
contains
  subroutine mpi_init(error)
    implicit none
    integer, intent(out) :: error
  ! this procedure has nothing to do
    error = 0
    return
  end subroutine mpi_init

  subroutine mpi_finalize(error)
    implicit none
    integer, intent(out) :: error
  ! this procedure has nothing to do
    error = 0
    return
  end subroutine mpi_finalize

  subroutine mpi_abort(comm, errorcode, error)
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer, intent(in) :: comm
    integer, intent(in) :: errorcode
    integer, intent(out) :: error
    ! all this procedure needs to do is terminate the current process
    ! with a non-zero exit code. we use the c intrinsic in order to
    ! ensure portability.
    interface
      subroutine c_exit(status) bind(c, name="exit")
      import :: c_int
      implicit none
      integer(kind=c_int), value, intent(in) :: status
      end subroutine c_exit
    end interface
    error = 1
    call c_exit(1_c_int)
    return
  end subroutine mpi_abort

  subroutine mpi_comm_size(comm, size, error)
    implicit none
    integer, intent(in) :: comm
    integer, intent(out) :: size
    integer, intent(out) :: error
  ! in serial mode, the size is always 1
    size  = 1
    error = 0
    return
  end subroutine mpi_comm_size

  subroutine mpi_comm_rank(comm, rank, error)
    implicit none
    integer, intent(in) :: comm
    integer, intent(out) :: rank
    integer, intent(out) :: error
  ! in serial mode, the rank is always 0
    rank  = 0
    error = 0
    return
  end subroutine mpi_comm_rank

  subroutine mpi_scatterv(sendbuf, sendcnts, displs, sendtype,                  &
                          recvbuf, recvcnt, recvtype,                           &
                          root, comm, error)
    implicit none
    real, intent(in) :: sendbuf(:,:)
    integer, intent(in) :: sendcnts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: sendtype
    real, intent(out) :: recvbuf(:,:)
    integer, intent(in) :: recvcnt
    integer, intent(in) :: recvtype
    integer, intent(in) :: root
    integer, intent(in) :: comm
    integer, intent(out) :: error
  ! this routine ignores all count and offset information and just copies all
  ! the data from sendbuf into recvbuf
  ! at the moment, this routine only exists for 2d real arrays, which is all jules uses
  ! it for
  ! the input and output buffers must be the same size
    if ( size(sendbuf) /= size(recvbuf) ) then
      write(*,*) 'fatal error: attempt to use dummy mpi procedure for non-serial work'
      write(*,*) 'to use mpi routines for parallel execution, recompile using an mpi compiler'
      call mpi_abort(comm, 1, error)
    end if
  ! just copy the data from the input buffer to the output buffer
    recvbuf(:,:) = sendbuf(:,:)
    error = 0
    return
  end subroutine mpi_scatterv

  subroutine mpi_gatherv(sendbuf, sendcnt, sendtype,                             &
                         recvbuf, recvcnts, displs, recvtype,                    &
                         root, comm, error)
    implicit none
    real, intent(in) :: sendbuf(:,:)
    integer, intent(in) :: sendcnt
    integer, intent(in) :: sendtype
    real, intent(out) :: recvbuf(:,:)
    integer, intent(in) :: recvcnts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: recvtype
    integer, intent(in) :: root
    integer, intent(in) :: comm
    integer, intent(out) :: error
  ! this routine ignores all count and offset information and just copies all
  ! the data from sendbuf into recvbuf
  ! at the moment, this routine only exists for 2d real arrays, which is all jules uses
  ! it for
  ! the input and output buffers must be the same size
    if ( size(sendbuf) /= size(recvbuf) ) then
      write(*,*) 'fatal error: attempt to use dummy mpi procedure for non-serial work'
      write(*,*) 'to use mpi routines for parallel execution, recompile using an mpi compiler'
      call mpi_abort(comm, 1, error)
    end if
  ! just copy the data from the input buffer to the output buffer
    recvbuf(:,:) = sendbuf(:,:)
    error = 0
    return
  end subroutine mpi_gatherv

  subroutine mpi_bcast1(buffer, count, datatype, root, comm, error)
    implicit none
    integer, intent(in) :: count
    real, intent(inout) :: buffer(count)
    integer, intent(in) :: datatype
    integer, intent(in) :: root
    integer, intent(in) :: comm
    integer, intent(out) :: error
  ! since there is only one task, we don't need to do anything to broadcast a value
    error = 0
    return
  end subroutine mpi_bcast1
  
  subroutine mpi_bcast2(buffer, count, datatype, root, comm, error)
    implicit none
    integer, intent(in) :: count
    real, intent(inout) :: buffer
    integer, intent(in) :: datatype
    integer, intent(in) :: root
    integer, intent(in) :: comm
    integer, intent(out) :: error
  ! since there is only one task, we don't need to do anything to broadcast a value
    error = 0
    return
  end subroutine mpi_bcast2
  
  subroutine mpi_bcast3(buffer, count, datatype, root, comm, error)
    implicit none
    integer, intent(in) :: count
    logical, intent(inout) :: buffer
    integer, intent(in) :: datatype
    integer, intent(in) :: root
    integer, intent(in) :: comm
    integer, intent(out) :: error
  ! since there is only one task, we don't need to do anything to broadcast a value
    error = 0
    return
  end subroutine mpi_bcast3
  
  subroutine mpi_bcast4(buffer, count, datatype, root, comm, error)
    implicit none
    integer, intent(in) :: count
    integer, intent(in) :: buffer(count)
    integer, intent(in) :: datatype
    integer, intent(in) :: root
    integer, intent(in) :: comm
    integer, intent(out) :: error
  ! since there is only one task, we don't need to do anything to broadcast a value
    error = 0
    return
  end subroutine mpi_bcast4
  
  subroutine mpi_bcast5(buffer, count, datatype, root, comm, error)
    implicit none
    integer, intent(in) :: count
    integer, intent(in) :: buffer
    integer, intent(in) :: datatype
    integer, intent(in) :: root
    integer, intent(in) :: comm
    integer, intent(out) :: error
  ! since there is only one task, we don't need to do anything to broadcast a value
    error = 0
    return
  end subroutine mpi_bcast5

  subroutine mpi_allreduce(sendbuf, recvbuf, sendcnt, datatype, op, comm, error)
    implicit none
    logical, intent(in) :: sendbuf
    logical, intent(out) :: recvbuf
    integer, intent(in) :: sendcnt
    integer, intent(in) :: datatype
    integer, intent(in) :: op
    integer, intent(in) :: comm
    integer, intent(out) :: error
  ! this routine only exists in the dummy library for logical variables
  ! the op is ignored, and this is just treated like a bcast
    recvbuf = sendbuf
    error = 0
    return
  end subroutine mpi_allreduce
  subroutine mpi_barrier(comm, error)
    implicit none
    integer, intent(in) :: comm
    integer, intent(out) :: error
  ! since there is only one task, we don't need to do anything to form a barrier
    error = 0
    return
  end subroutine mpi_barrier

  subroutine mpi_type_get_extent(type, lb, extent, error)
    implicit none
    integer, intent(in) :: type
    integer(kind=mpi_address_kind), intent(out) :: lb
    integer(kind=mpi_address_kind), intent(out) :: extent
    integer, intent(out) :: error
  ! all information about mpi types is ignored by the dummy mpi routines, so
  ! it doesn't matter what we return
    lb     = 0
    extent = 4
    error  = 0
    return
  end subroutine mpi_type_get_extent

  subroutine mpi_type_vector(count, blocklength, stride, old_type, new_type, error)
    implicit none
    integer, intent(in) :: count
    integer, intent(in) :: blocklength
    integer, intent(in) :: stride
    integer, intent(in) :: old_type
    integer, intent(out) :: new_type
    integer, intent(out) :: error
  ! all information about mpi types is ignored by the dummy mpi routines, so
  ! it doesn't matter what we return
    new_type = 1
    error    = 0
    return
  end subroutine mpi_type_vector

  subroutine mpi_type_create_resized1(old_type, lb, extent, new_type, error)
    implicit none
    integer, intent(in) :: old_type
    integer, intent(in) :: lb
    integer, intent(in) :: extent
    integer, intent(out) :: new_type
    integer, intent(out) :: error
  ! all information about mpi types is ignored by the dummy mpi routines, so
  ! it doesn't matter what we return
    new_type = 1
    error    = 0
    return
  end subroutine mpi_type_create_resized1

  subroutine mpi_type_create_resized2(old_type, lb, extent, new_type, error)
    implicit none
    integer, intent(in) :: old_type
    integer(kind=mpi_address_kind) :: lb
    integer(kind=mpi_address_kind),intent(in) ::  extent
    integer, intent(out) :: new_type
    integer, intent(out) :: error
  ! all information about mpi types is ignored by the dummy mpi routines, so
  ! it doesn't matter what we return
    new_type = 1
    error    = 0
    return
  end subroutine mpi_type_create_resized2

  subroutine mpi_type_commit(datatype, error)
    implicit none
    integer, intent(in) :: datatype
    integer, intent(out) :: error
  ! all information about mpi types is ignored by the dummy mpi routines, so
  ! there is nothing to do
    error = 0
    return
  end subroutine mpi_type_commit
END MODULE lis_mpi_dummy
