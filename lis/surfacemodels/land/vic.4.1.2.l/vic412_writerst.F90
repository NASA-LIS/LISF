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
! !ROUTINE: vic412_writerst
! \label{vic412_writerst}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 01 May 2012; Shugong Wang. Add vegetation tiling scheme as an option, to support
!              writing restart file using VIC based vegetation tiling and LIS based 
!              vegetation tiling. When using VIC based vegetation tiling, grid id will
!              be the key of restart file. Otherwise, tile id will be the key of restart
!              file. 
! !INTERFACE:
subroutine vic412_writerst(n)
! !USES:
   use LIS_coreMod,    only : LIS_rc, LIS_masterproc, LIS_localPet, LIS_npes, LIS_domain
   use LIS_timeMgrMod, only : LIS_isAlarmRinging
   use LIS_logMod,     only : LIS_logunit, LIS_verify, LIS_getNextUnitNumber, LIS_releaseUnitNumber
   use LIS_fileIOMod,  only : LIS_create_output_directory, &
                              LIS_create_restart_filename
   use LIS_constantsMod, only : LIS_CONST_PATH_LEN
   use vic412_lsmMod,  only : vic412_struc
   use LIS_historyMod 
   use netcdf 
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
!  vic412\_write\_state).  When running in parallel, this call must be made
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
   integer       :: vt_scheme ! added by Shugong Wang to support VIC and LIS based tilings. 
   integer, allocatable, dimension(:) :: vegclasses! added by Shugong Wang to support MPI-safe tile ID, 05/07/2012
   integer       :: t, chunk_size
   real, allocatable :: state_chunk(:) ! chunk of state variable 
   integer       :: state_chunk_ID, dimID(11) 
   integer       :: ftn, l
   real, allocatable :: tmptilen(:)
   character*3   :: fnest 
   character*20  :: wformat
   
   !wformat = "binary"
   !wformat = "netcdf"
   wformat  = trim(vic412_struc(n)%rfile_format) 
   allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
   write(fnest, '(i3.3)') n 

   chunk_size = vic412_struc(n)%state_chunk_size
   allocate(state_chunk(chunk_size))

   ! Restart data cannot be written until the physics has been run.
   ! When VIC physics is run, tscount will be incremented to a value
   ! greater than 1.
   ! This supports VIC's water-balance mode where LIS may need to
   ! process several forcing/snow time-steps before it is ready to run
   ! VIC's physics.
   if ( vic412_struc(n)%tscount == 1 ) then
       return
   endif

   alarmCheck = LIS_isAlarmRinging(LIS_rc,"VIC412 restart alarm "//trim(fnest))

   !if ( alarmCheck .or. (LIS_rc%endtime == 1) ) then 
   if(alarmCheck .or. (LIS_rc%endtime ==1)) then 
      if ( LIS_masterproc ) then
         call LIS_create_output_directory('SURFACEMODEL')
         call LIS_create_restart_filename(n,filen,'SURFACEMODEL','VIC412', wformat=wformat)
         if(wformat.eq."binary") then 
            ftn = LIS_getNextUnitNumber()       
            open(ftn,file=filen,status='unknown',form='unformatted')
         elseif(wformat.eq."netcdf") then 
#if (defined USE_NETCDF4)
            status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
            call LIS_verify(status,"Error in nf90_open in VIC412_writerst")
#endif
#if (defined USE_NETCDF3)
           status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
           call LIS_verify(status, "Error in nf90_open in VIC412_writerst")
#endif
         endif
      endif

#if (defined SPMD) 
      call mpi_barrier(LIS_mpi_comm, status)
#endif


      ! added by Shugong Wang 05/01/2012
      vt_scheme = vic412_struc(n)%veg_tiling_scheme  
      allocate(vegclasses(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
        vegclasses(t) = LIS_domain(n)%tile(t)%vegt
      enddo     
      
      do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
          call vic412_write_state(n, LIS_rc%npatch(n,LIS_rc%lsm_index), vt_scheme,&
                                  vegclasses,                                     &
                                  LIS_rc%yr, LIS_rc%mo, LIS_rc%da,                &
                                  vic412_struc(n)%vic(t)%state_chunk,  t) 
      enddo
      ! added by Shugong Wang 05/07/2012, free memory for vegclasses
      deallocate(vegclasses)



      ! dump restart
      ! write the header of the restart file
      call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, &
                                       "VIC 412", dim1=chunk_size, dimID=dimID, &
                                       output_format=wformat)

      ! write the header for state variable chunk
      if ( wformat == "netcdf" ) then
         call LIS_writeHeader_restart(ftn, n, dimID, state_chunk_ID, "state_chunk", &
                                    "chunk of all VIC 4.1.2 state variables", &
                                    "-", vlevels=chunk_size, valid_min=-99999.0, & 
                                    valid_max=99999.0, var_flag="dim1")
      endif
    
      do l=1, chunk_size  
          tmptilen = 0
          do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
              tmptilen(t) = vic412_struc(n)%vic(t)%state_chunk(l)
          enddo
          call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=state_chunk_ID, dim=l, wformat=wformat)
      enddo
    
   if(LIS_masterproc ) then
      if(wformat == "binary") then 
         call LIS_releaseUnitNumber(ftn)
      elseif ( wformat == "netcdf" ) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          status = nf90_close(ftn)
          call LIS_verify(status, "Error in nf90_close in VIC412_writerst")
#endif
      endif
   endif

#if (defined SPMD) 
      call mpi_barrier(LIS_mpi_comm, status)
#endif

   endif
   deallocate(tmptilen)
   deallocate(state_chunk)
end subroutine vic412_writerst
