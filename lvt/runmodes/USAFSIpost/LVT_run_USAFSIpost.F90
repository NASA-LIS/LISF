!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine LVT_run_USAFSIpost()

   ! Imports
   use LVT_USAFSIpostMod
   use LVT_logMod
   use LVT_coreMod, only: LVT_rc

   ! Defaults
   implicit none

   ! Local variables
   type(LVT_USAFSIpost_t) :: USAFSIpost
   integer :: i

   call USAFSIpost%new()
   call USAFSIpost%read_usafsi_ncfile()

   if (LVT_rc%output_native) then
      call USAFSIpost%output_grib2()
   end if

   if (LVT_rc%output_global_ll0p25) then
      call USAFSIpost%interp_and_output_grib1(GLOBAL_LL0P25)
   end if

   if (LVT_rc%output_nh_ps16 .or. LVT_rc%output_nh_ps16_snodep) then
      call USAFSIpost%interp_and_output_grib1(NH_PS16)
   end if

   if (LVT_rc%output_sh_ps16 .or. LVT_rc%output_sh_ps16_snodep) then
      call USAFSIpost%interp_and_output_grib1(SH_PS16)
   end if

   ! Clean up
   call USAFSIpost%delete()

   write(LVT_logunit,*)'[INFO] LVT successfully ended!'

end subroutine LVT_run_USAFSIpost
