!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: clm2_writerestart
! \label{clm2_writerestart}
!
! !REVISION HISTORY:
! 20 Jan 2003; Sujay Kumar Initial Specification
! 26 Aug 2004; James Geiger, Added support for GrADS-DODS based and MPI based
!                            parallel simulations.
! 
! !INTERFACE:
subroutine clm2_writerestart(n)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod, only : LIS_logunit
  use LIS_fileIOMod, only : LIS_create_output_directory, &
                              LIS_create_restart_filename
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use clm2_lsmMod, only : clm2_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
!
! !DESCRIPTION:
!  This program writes restart files for CLM.  This
!  includes all relevant water/energy storage and tile information
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory})
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename})
!  generates a timestamped restart filename
! \item[clm\_dump\_restart](\ref{clm2_dump_restart})
!   writes the CLM variables into the restart file
!  reads a variable from the restart file
! \end{description}
!EOP
  character(len=LIS_CONST_PATH_LEN) :: filen
  logical :: alarmCheck
  integer :: status
  character*3   :: fnest
  
  write(fnest,'(i3.3)') n

  alarmCheck = LIS_isAlarmRinging(LIS_rc, &
       "CLM2 restart alarm "//trim(fnest))
  if(alarmCheck.or.(LIS_rc%endtime ==1)) then 
     if(LIS_masterproc) then 
        call LIS_create_output_directory('SURFACEMODEL')
        call LIS_create_restart_filename(n, filen,'SURFACEMODEL','CLM')
        open(40,file=filen,status='unknown',form='unformatted')
     endif
     call clm2_dump_restart(n, 40)
     if(LIS_masterproc) then 
        close(40)
        write(LIS_logunit,*) 'CLM Archive restart written ',trim(filen)
     endif
  endif
end subroutine clm2_writerestart

!BOP
! 
! !ROUTINE: clm2_dump_restart
! \label{clm2_dump_restart}
! 
! !INTERFACE: 
subroutine clm2_dump_restart(n, ftn)
! !USES:
   use clm2_lsmMod
   use clm2_varpar,   only : nlevsoi, nlevsno
   use LIS_coreMod,  only : LIS_rc, LIS_masterproc
   use LIS_historyMod, only : LIS_writevar_restart
   
   implicit none
! !ARGUMENTS:    
   integer, intent(in) :: n
   integer, intent(in) :: ftn
!
! !DESCRIPTION:
!  This routine gathers the necessary restart variables and performs
!  the actual write statements to create the restart files.
!
!  The arguments are: 
!  \begin{description}
!   \item[n] 
!    index of the nest
!   \item[ftn]
!    unit number for the restart file
!  \end{description}
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_writevar\_restart](\ref{LIS_writevar_restart})
!  writes a variable to the restart file
! \end{description}

!EOP
   integer :: t

   if(LIS_masterproc) then
      write(ftn) LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%glbnpatch(n,LIS_rc%lsm_index)  !Veg class, no tiles
   endif

   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%snl)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%frac_veg_nosno_alb)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%h2osno)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%h2ocan)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%snowdp)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%snowage)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%frac_sno)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%t_veg)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%t_grnd)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fwet)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%tlai)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%tsai)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%elai)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%esai)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fsun)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%htop)
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%hbot)
   do t=-nlevsno+1, 0 
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%dz(t))
   enddo
   do t=-nlevsno+1, 0
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%z(t))
   enddo
   do t=-nlevsno, 0
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%zi(t))
   enddo
   
   do t=-nlevsno+1, nlevsoi
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%t_soisno(t))
   enddo

   do t=1, nlevsoi
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%t_lake(t))
   enddo
   
   do t=-nlevsno+1, nlevsoi
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%h2osoi_liq(t))
   enddo
    
   do t=-nlevsno+1,nlevsoi
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%h2osoi_ice(t))
   enddo
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%albgrd(1))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%albgrd(2))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%albgri(1))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%albgri(2))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fabd(1))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fabd(2))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fabi(1))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fabi(2))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftdd(1))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftdd(2))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftid(1))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftid(2))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftii(1))
   call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftii(2))
 end subroutine clm2_dump_restart


