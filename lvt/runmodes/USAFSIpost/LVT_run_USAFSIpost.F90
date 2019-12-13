!-----------------------BEGIN NOTICE -- DO NOT EDIT----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT----------------------------

subroutine LVT_run_USAFSIpost()

   ! Imports
   use LVT_USAFSIpostMod
   use LVT_logMod

   ! Defaults
   implicit none

   ! Local variables
   type(LVT_USAFSIpost_t) :: USAFSIpost
   character(len=15) :: grids(3)
   integer :: i

   grids(1) = GLOBAL_LL0P25
   grids(2) = NH_PS16
   grids(3) = SH_PS16
   call USAFSIpost%new()
   call USAFSIpost%read_usafsi_ncfile()
   call USAFSIpost%output_grib2()
   do i = 1, 3
      call USAFSIpost%interp_and_output_grib1(grids(i))
   end do
   call USAFSIpost%delete()
   call LVT_endrun()
end subroutine LVT_run_USAFSIpost
