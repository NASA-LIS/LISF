!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LVT_StratStatsMod
!
!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !MODULE: LVT_StratStatMod
! 
!  !DESCRIPTION: 
!   This module contains routines related to the stratification of 
!   evaluation statistics based on specified external data. 
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
!
!EOP
  use ESMF
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initStratStats

  private
  
contains

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_initStratStats
! \label{LVT_initStratStats}
! 
! !INTERFACE: 
  subroutine LVT_initStratStats
! !USES: 
    use LVT_coreMod
    use LVT_logMod
!
! !DESCRIPTION:
! 
!EOP    
    implicit none

    integer               :: ftn 
    integer               :: ios1,rc
    integer               :: i,c,r
    logical               :: file_exists
    real                  :: strat_data(LVT_rc%lnc, LVT_rc%lnr)

    if(LVT_rc%data_based_strat.eq.1) then 
       call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%data_strat_attrib_file, &
            label="Stratification attributes file:",rc=rc)
       call LVT_verify(rc,"Stratification attributes file: not defined")
       
       inquire(file=LVT_rc%data_strat_attrib_file,exist=file_exists)
       if(.not.file_exists) then 
          write(LVT_logunit,*) '[ERR] '
          write(LVT_logunit,*) '[ERR] stratification attributes file '//&
               trim(LVT_rc%data_strat_attrib_file) //' not found'
          write(LVT_logunit,*) '[ERR] A sample file is shown below:'
          write(LVT_logunit,*) '[ERR] #Number of stratification data sources '
          write(LVT_logunit,*) '[ERR] 3'
          write(LVT_logunit,*) '[ERR] #Stratification data files '
          write(LVT_logunit,*) '[ERR] srtm_elev1km.1gd4r'
          write(LVT_logunit,*) '[ERR] srtm_slope1km.1gd4r'
          write(LVT_logunit,*) '[ERR] srtm_aspect1km.1gd4r'
          write(LVT_logunit,*) '[ERR] #stratifcation variable name'
          write(LVT_logunit,*) '[ERR] ELEV'
          write(LVT_logunit,*) '[ERR] SLOPE'
          write(LVT_logunit,*) '[ERR] ASPECT'
          write(LVT_logunit,*) '[ERR] #Max min values'
          write(LVT_logunit,*) '[ERR] 7000 1.0 6'
          write(LVT_logunit,*) '[ERR] 500  0.0 0'
          write(LVT_logunit,*) '[ERR] #number of bins'
          write(LVT_logunit,*) '[ERR] 12 12 12'
          write(LVT_logunit,*) '[ERR]'
          call LVT_endrun()
       endif

       ftn = LVT_getNextUnitNumber()
       open(ftn,file=(LVT_rc%data_strat_attrib_file), status='old')
       read(ftn,*) 
       read(ftn,*) LVT_rc%data_based_nstrats
       read(ftn,*) 
       allocate(LVT_rc%data_based_strat_file(LVT_rc%data_based_nstrats))
       allocate(LVT_rc%data_based_strat_var(LVT_rc%data_based_nstrats))
       allocate(LVT_rc%data_based_strat_max(LVT_rc%data_based_nstrats))
       allocate(LVT_rc%data_based_strat_min(LVT_rc%data_based_nstrats))
       allocate(LVT_rc%data_based_strat_nbins(LVT_rc%data_based_nstrats))
       allocate(LVT_rc%data_based_strat_delta(LVT_rc%data_based_nstrats))

       do i=1,LVT_rc%data_based_nstrats
          read(ftn,'(a)') LVT_rc%data_based_strat_file(i)
       enddo
       read(ftn,*) 
       
       do i=1,LVT_rc%data_based_nstrats
          read(ftn,'(a)') LVT_rc%data_based_strat_var(i)
       enddo
       read(ftn,*)
       read(ftn,*) (LVT_rc%data_based_strat_max(i),&
            i=1,LVT_rc%data_based_nstrats)
       read(ftn,*) (LVT_rc%data_based_strat_min(i),&
            i=1,LVT_rc%data_based_nstrats)
       read(ftn,*) 
       read(ftn,*) (LVT_rc%data_based_strat_nbins(i),&
            i=1,LVT_rc%data_based_nstrats)
       call LVT_releaseUnitNumber(ftn)

!read each strat file
       allocate(LVT_rc%strat_data(LVT_rc%data_based_nstrats, &
            LVT_rc%ngrid))

       do i=1,LVT_rc%data_based_nstrats
          ftn = LVT_getNextUnitNumber()       
          inquire(file=(LVT_rc%data_based_strat_file(i)), exist=file_exists)
          if(file_exists) then 
             write(LVT_logunit,*)'Reading.. ',(LVT_rc%data_based_strat_file(i))
             ftn = LVT_getNextUnitNumber()
             
             open(ftn,file=(LVT_rc%data_based_strat_file(i)),&
                  form='unformatted',iostat=ios1)
             read(ftn) strat_data
             
             do r=1,LVT_rc%lnr
                do c=1,LVT_rc%lnc
                   if(LVT_domain%gindex(c,r).ne.-1) then 
                      LVT_rc%strat_data(i,LVT_domain%gindex(c,r)) = &
                           strat_data(c,r)
                   endif
                enddo
             enddo

             LVT_rc%data_based_strat_delta(i) = &
                  (LVT_rc%data_based_strat_max(i)-&
                  LVT_rc%data_based_strat_min(i))/&
                  LVT_rc%data_based_strat_nbins(i)
             call LVT_releaseUnitNumber(ftn)             
          else
             write(LVT_logunit,*) '[ERR] '
             write(LVT_logunit,*) '[ERR] external data file '//&
                  trim(LVT_rc%data_based_strat_file(i))
             call LVT_endrun()
          endif
       enddo
    end if
  end subroutine LVT_initStratStats
end module LVT_StratStatsMod
