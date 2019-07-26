!-----------------------BEGIN NOTICE -- DO NOT EDIT----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT----------------------------

subroutine LVT_run_LDTSIpost()

   ! Imports
   use LVT_LDTSIpostMod
   use LVT_logMod

   ! Defaults
   implicit none

   ! Local variables
   type(LVT_LDTSIpost_t) :: LDTSIpost
   character(len=15) :: grids(3)
   integer :: i

   grids(1) = GLOBAL_LL0P25
   grids(2) = NH_PS16
   grids(3) = SH_PS16
   call LDTSIpost%new()
   call LDTSIpost%read_ldtsi_ncfile()
   call LDTSIpost%output_grib2()
   do i = 1, 3
      call LDTSIpost%interp_and_output_grib1(grids(i))
   end do
   call LDTSIpost%delete()
   call LVT_endrun()
end subroutine LVT_run_LDTSIpost
