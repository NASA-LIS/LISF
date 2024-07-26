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
!BOP
!
! !ROUTINE: vic411_writerst
! \label{vic411_writerst}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 01 May 2012; Shugong Wang. Add vegetation tiling scheme as an option, to support
!              writing restart file using VIC based vegetation tiling and LIS based 
!              vegetation tiling. When using VIC based vegetation tiling, grid id will
!              be the key of restart file. Otherwise, tile id will be the key of restart
!              file. 
! !INTERFACE:
subroutine vic411_writerst(n)
! !USES:
   use LIS_coreMod,    only : LIS_rc, LIS_masterproc, LIS_localPet, LIS_npes, LIS_domain
   use LIS_timeMgrMod, only : LIS_isAlarmRinging
   use LIS_logMod,     only : LIS_logunit
   use LIS_fileIOMod,  only : LIS_create_output_directory, &
                              LIS_create_restart_filename
   use LIS_constantsMod, only : LIS_CONST_PATH_LEN
   use vic411_lsmMod,  only : vic411_struc
#if (defined SPMD) 
   use LIS_mpiMod
#endif

   implicit none
! !ARGUMENTS: 
   integer, intent(in)  :: n
!
! !DESCRIPTION:
!  This routine writes restart files for VIC.  This
!  includes all relevant water/energy storage and tile information.
!
!  This routine calls VIC's write\_model\_state routine (via the call to
!  vic411\_write\_state).  When running in parallel, this call must be made
!  by all processes, resulting in one restart file per process.  This
!  routine then concatenates these files into one (via system calls to
!  Unix commands).
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
   character(len=LIS_CONST_PATH_LEN) :: filen
   logical       :: alarmCheck
   integer       :: status, i
   character*4   :: fproc
   integer :: vt_scheme ! added by Shugong Wang to support VIC and LIS based tilings. 
   integer, allocatable, dimension(:) :: vegclasses! added by Shugong Wang to support MPI-safe tile ID, 05/07/2012
   integer :: t
   character*3 :: fnest

   write(fnest,'(i3.3)') n

   ! Restart data cannot be written until the physics has been run.
   ! When VIC physics is run, tscount will be incremented to a value
   ! greater than 1.
   ! This supports VIC's water-balance mode where LIS may need to
   ! process several forcing/snow time-steps before it is ready to run
   ! VIC's physics.
   if ( vic411_struc(n)%tscount == 1 ) then
       return
   endif
   
   alarmCheck = LIS_isAlarmRinging(LIS_rc,"VIC411 restart alarm "//trim(fnest))

   !if ( alarmCheck .or. (LIS_rc%endtime == 1) ) then 
   if ( alarmCheck ) then 

      if ( LIS_masterproc ) then
         call LIS_create_output_directory('SURFACEMODEL')
      endif
#if (defined SPMD) 
      call mpi_barrier(LIS_mpi_comm, status)
#endif

      call LIS_create_restart_filename(n,filen,'SURFACEMODEL','VIC411')

      write(unit=fproc,fmt='(i4.4)') LIS_localPet
      ! added by Shugong Wang 05/01/2012
      vt_scheme = vic411_struc(n)%veg_tiling_scheme  
      ! added by Shugong Wang 05/07/2012 
      allocate(vegclasses(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
        vegclasses(t) = LIS_domain(n)%tile(t)%vegt
      enddo      
      call vic411_write_state(n, LIS_rc%npatch(n,LIS_rc%lsm_index), vt_scheme,&
                              vegclasses,                                     &
                              LIS_rc%yr, LIS_rc%mo, LIS_rc%da,                &
                              trim(filen)//"."//fproc//char(0))
     ! added by Shugong Wang 05/07/2012, free memory for vegclasses
     deallocate(vegclasses)
#if (defined SPMD) 
      call mpi_barrier(LIS_mpi_comm, status)
#endif
      if ( LIS_masterproc ) then
!         close(40)   
         call system('cp '//trim(filen)//"."//fproc//' '//trim(filen))
         do i = 1, LIS_npes-1
            write(unit=fproc,fmt='(i4.4)') i
            call system('tail +3 '//trim(filen)//"."//fproc//'>>'//trim(filen))
         enddo
         write(LIS_logunit,*) 'VIC 4.1.1 archive restart written: ',trim(filen)
      endif
   endif

end subroutine vic411_writerst
