!------------------------------------------------------------------------------
!BOP
!
#ifdef USE_PFIO

#include "MAPL_ErrLog.h"
#include "unused_dummy.H"
#include "LIS_misc.h"

#endif

module LIS_PFIO_utilsMod

#ifdef USE_PFIO
! !USES:
      use mpi
      use MAPL
      use pFIO_UnlimitedEntityMod
      use LIS_coreMod
      use LIS_logMod
      use LIS_mpiMod
      
      implicit none

      PRIVATE
!
! !DESCRIPTION:
! Contain uitility subroutines needed to write data in netCDF files using
! GEOS/PFIO.

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
      public :: add_fvar
      public :: pfio_write_var

      interface pfio_write_var
         module procedure pfio_write_var1d_int
         module procedure pfio_write_var1d_real
         module procedure pfio_write_var2d
         module procedure pfio_write_var3d
      end interface
#else
      integer :: dummy_int
#endif

#ifdef USE_PFIO
!
!EOP
!---------------------------------------------------------------------
CONTAINS
!---------------------------------------------------------------------
!BOP
      subroutine add_fvar(fmd, var_name, var_type, dims, &
                          units, long_name, standard_name, &
                          scale_factor, add_offset, missing_value, &
                          fill_value, vmin, vmax)
!
! !INPUT PATAMETERS:
      integer,          intent(in) :: var_type
      character(len=*), intent(in) :: var_name
      character(len=*), intent(in) :: dims
      character(len=*), optional, intent(in) :: units
      character(len=*), optional, intent(in) :: long_name
      character(len=*), optional, intent(in) :: standard_name
      real,             optional, intent(in) :: scale_factor
      real,             optional, intent(in) :: add_offset
      real,             optional, intent(in) :: missing_value
      real,             optional, intent(in) :: fill_value
      real,             optional, intent(in) :: vmin
      real,             optional, intent(in) :: vmax
!
! !INPUT/OUPUT PATAMETERS:
      type(FileMetadata), intent(inout) :: fmd
!
! !DESCRIPTION:
!  Create a variable and add attributes (units, long_name).
!
! !LOCAL VARIABLES:
      integer :: status, rc
      type(Variable) :: fvar
!EOP
!---------------------------------------------------------------------
!BOC
      fvar = Variable(type=var_type, dimensions=TRIM(dims))
      if (present(units))     call fvar%add_attribute('units', TRIM(units))
      if (present(long_name)) call fvar%add_attribute('long_name', TRIM(long_name))
      if (present(standard_name)) call fvar%add_attribute('standard_name', TRIM(standard_name))
      if (present(scale_factor)) call fvar%add_attribute('scale_factor', scale_factor)
      if (present(add_offset)) call fvar%add_attribute('add_offset', add_offset)
      if (present(missing_value)) call fvar%add_attribute('missing_value', missing_value)
      if (present(fill_value)) call fvar%add_attribute('_FillValue', fill_value)
      if (present(vmin)) call fvar%add_attribute('vmin', vmin)
      if (present(vmax)) call fvar%add_attribute('vmax', vmax)
      call fmd%add_variable(TRIM(var_name), fvar, rc=status)
      call LIS_verify(status, 'Adding variable atrributes failed for '//TRIM(var_name))

      end subroutine add_fvar
!EOC
!---------------------------------------------------------------------
!BOP
      subroutine pfio_write_var1d_int(hist_id, file_name, local_var, var_name, &
                              local_start, global_start, global_count)
!
! !INPUT PARAMETERS:
       integer,          intent(in) :: hist_id
       character(len=*), intent(in) :: file_name
       integer,          intent(in) :: local_var(:)
       character(len=*), intent(in) :: var_name
       integer,          intent(in) :: local_start(:)
       integer,          intent(in) :: global_start(:)
       integer,          intent(in) :: global_count(:)
!
! !DESCRIPTION:
! Write a 1D variable in a file
!
! !LOCAL VARIABLES:
       type(ArrayReference)    :: ref
!EOP
!---------------------------------------------------------------------
!BOC
       IF ( (SIZE(local_start)  /= SIZE(global_start)) .OR. &
            (SIZE(local_start)  /= SIZE(global_count)) .OR. &
            (SIZE(global_start) /= SIZE(global_count)) ) THEN
          call LIS_verify(-1, "Arrays of different sizes for "//TRIM(var_name))
       END IF

       ref =  ArrayReference(local_var)
       call o_clients%collective_stage_data(hist_id, TRIM(file_name), &
                                      TRIM(var_name), ref, &
                                      start        = local_start, &
                                      global_start = global_start, &
                                      global_count = global_count)
      end subroutine pfio_write_var1d_int
!EOC
!---------------------------------------------------------------------
!BOP
      subroutine pfio_write_var1d_real(hist_id, file_name, local_var, var_name, &
                              local_start, global_start, global_count)
!
! !INPUT PARAMETERS:
       integer,          intent(in) :: hist_id
       character(len=*), intent(in) :: file_name
       real,             intent(in) :: local_var(:)
       character(len=*), intent(in) :: var_name
       integer,          intent(in) :: local_start(:)
       integer,          intent(in) :: global_start(:)
       integer,          intent(in) :: global_count(:)
!
! !DESCRIPTION:
! Write a 1D variable in a file
!
! !LOCAL VARIABLES:
       type(ArrayReference)    :: ref
!EOP
!---------------------------------------------------------------------
!BOC
       IF ( (SIZE(local_start)  /= SIZE(global_start)) .OR. &
            (SIZE(local_start)  /= SIZE(global_count)) .OR. &
            (SIZE(global_start) /= SIZE(global_count)) ) THEN
          call LIS_verify(-1, "Arrays of different sizes for "//TRIM(var_name))
       END IF

       ref =  ArrayReference(local_var)
       call o_clients%collective_stage_data(hist_id, TRIM(file_name), &
                                      TRIM(var_name), ref, &
                                      start        = local_start, &
                                      global_start = global_start, &
                                      global_count = global_count)
      end subroutine pfio_write_var1d_real
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine pfio_write_var2d(hist_id, file_name, local_var, var_name, &
                              local_start, global_start, global_count)
!
! !INPUT PARAMETERS:
       integer,          intent(in) :: hist_id
       character(len=*), intent(in) :: file_name
       real,             intent(in) :: local_var(:,:)
       character(len=*), intent(in) :: var_name
       integer,          intent(in) :: local_start(:)
       integer,          intent(in) :: global_start(:)
       integer,          intent(in) :: global_count(:)
!
! !DESCRIPTION:
! Write a 2D variable in a file
!
! !LOCAL VARIABLES:
       type(ArrayReference)    :: ref
!EOP
!---------------------------------------------------------------------
!BOC
       IF ( (SIZE(local_start)  /= SIZE(global_start)) .OR. &
            (SIZE(local_start)  /= SIZE(global_count)) .OR. &
            (SIZE(global_start) /= SIZE(global_count)) ) THEN
          call LIS_verify(-1, "Arrays of different sizes for "//TRIM(var_name))
       END IF

       ref =  ArrayReference(local_var)
       call o_clients%collective_stage_data(hist_id, TRIM(file_name), &
                                      TRIM(var_name), ref, &
                                      start        = local_start, &
                                      global_start = global_start, &
                                      global_count = global_count)
      !if (LIS_masterproc) PRINT'(a,2f20.7)',"KNJR --> Done writing "//TRIM(var_name),minval(local_var),maxval(local_var)
      !PRINT'(a,2f20.7)',"KNJR --> Done writing "//TRIM(var_name),minval(local_var),maxval(local_var)
      end subroutine pfio_write_var2d
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine pfio_write_var3d(hist_id, file_name, local_var, var_name, &
                              local_start, global_start, global_count)
!
! !INPUT PARAMETERS:
       integer,          intent(in) :: hist_id
       character(len=*), intent(in) :: file_name
       real,             intent(in) :: local_var(:,:,:)
       character(len=*), intent(in) :: var_name
       integer,          intent(in) :: local_start(:)
       integer,          intent(in) :: global_start(:)
       integer,          intent(in) :: global_count(:)
!
! !DESCRIPTION:
! Write a 3D variable in a file
!
! !LOCAL VARIABLES:
       type(ArrayReference)    :: ref
!EOP
!---------------------------------------------------------------------
!BOC
       IF ( (SIZE(local_start)  /= SIZE(global_start)) .OR. &
            (SIZE(local_start)  /= SIZE(global_count)) .OR. &
            (SIZE(global_start) /= SIZE(global_count)) ) THEN
          call LIS_verify(-1, "Arrays of different sizes for "//TRIM(var_name))
       END IF

       ref =  ArrayReference(local_var)
       call o_clients%collective_stage_data(hist_id, TRIM(file_name), &
                                      TRIM(var_name), ref, &
                                      start        = local_start, &
                                      global_start = global_start, &
                                      global_count = global_count)
      end subroutine pfio_write_var3d
!EOC
!------------------------------------------------------------------------------
#endif
end module LIS_PFIO_utilsMod


