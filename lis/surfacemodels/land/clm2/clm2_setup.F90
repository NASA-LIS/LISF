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
! !ROUTINE: clm2_setup
! \label{clm2_setup}
! 
! !REVISION HISTORY: 
! 
! 20 Jan 2003  Sujay Kumar Initial Specification
! 
! !INTERFACE:
subroutine clm2_setup()
! !USES:
  use LIS_coreMod, only: LIS_rc
  use LIS_timeMgrMod, only : LIS_get_nstep
  use LIS_logMod, only :LIS_logunit
  use clm2_lsmMod
  use clm2_varcon, only: eccen, obliqr, lambm0 , mvelpp

!
!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for CLM. The method invokes the routine to call
!  time variant parameters. 
!  
!  The routines invoked are: 
!  \begin{description}
!   \item[canhtset](\ref{canhtset})
!    reads the canopy height information
!   \item[iniTimeVar](\ref{iniTimeVar})
!    initialize time variant parameters
!  \end{description}
!EOP
  implicit none
  integer  :: yr                    !current year (0, ...)
  integer  :: mon                   !current month (1 -> 12)
  integer  :: day                   !current day (1 -> 31)
  integer  :: i, ncsec
  integer  :: n 

! ----------------------------------------------------------------------
! Get current date
! ----------------------------------------------------------------------
  
  do n=1,LIS_rc%nnest
     yr = LIS_rc%yr
     mon = LIS_rc%mo
     day = LIS_rc%da
     ncsec = LIS_rc%ss
! ----------------------------------------------------------------------
! If initial run: initialize time-varying data 
! If continuation run: end of initialization because time varying
! read in from restart file
! ----------------------------------------------------------------------
  
     if (trim(LIS_rc%startcode).eq."coldstart") then
        call canhtset(n)
        write (LIS_logunit,*) ('Attempting to initialize time variant variables .....')
        
        call iniTimeVar (n, eccen, obliqr, lambm0 , mvelpp)
        
        write (LIS_logunit,*) ('Successfully initialized time variant variables')
        write (LIS_logunit,*)
     endif

! ----------------------------------------------------------------------
! End initialization
! ----------------------------------------------------------------------
     
     write (LIS_logunit,*) ('Successfully initialized the land model')
     if (trim(LIS_rc%startcode).eq."coldstart") then
        write (LIS_logunit,*) 'begin initial run at: '
     elseif(trim(LIS_rc%startcode).eq."restart") then 
        write (LIS_logunit,*) 'begin continuation run at:'
     end if
     write (LIS_logunit,*) '   nstep= ',LIS_get_nstep(LIS_rc,n), &
          ' year= ',yr,' month= ',mon,' day= ',day,' seconds= ',ncsec
     write (LIS_logunit,*)
     write (LIS_logunit,'(72a1)') ("*",i=1,60)
     write (LIS_logunit,*)
  enddo
end subroutine clm2_setup

