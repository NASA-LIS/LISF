!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_rdhm356
! \label{readcrd_rdhm356}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 05 Jun 2006; Kristi Arsenault, Code and data implementation
! 03 May 2010: Soni Yatheendradas; Precip and Temper. input grids now can
!              have different extents/directories and different from the
!              run-domain extent, as per the new input grids posted onto
!              the DMIP2 website for Sierra Nevada 
! 06 Man 2014: Shugong Wang, modified for RDHM 3.5.6 
! 
! !INTERFACE:    
subroutine readcrd_rdhm356()

! !USES:
  use ESMF
  use rdhm356_forcingMod, only : rdhm356_struc_precip, &
                                 rdhm356_struc_temper, &
                                 const_wind 
  use LIS_coreMod, only : LIS_config, LIS_rc
  use LIS_logMod, only : LIS_logunit
!
! !DESCRIPTION:
!
!  This routine reads the options specific to RDHM 3.5.6 forcing from 
!  the LIS configuration file. 
!  
!  The routines invoked are: 
!  \begin{description}
!   \item[ESMF\_ConfigGetAttribute]
!     read the specified attribute
!   \item[ESMF\_ConfigFindLabel]
!     find the specified label to read the attribute
!  \end{description}
!EOP

  implicit none
  integer :: n, rc

! - Retrieve dmip II Forcing Dataset Directory Name Location
    write(LIS_logunit,*) 'Using OHD RDHM forcing'
    call ESMF_ConfigFindLabel(LIS_config, "RDHM precipitation forcing directory:",&
                              rc=rc)
    do n=1,LIS_rc%nnest   
       call ESMF_ConfigGetAttribute(LIS_config, &
                                    rdhm356_struc_precip(n)%rdhm356dir,rc=rc)
       write(LIS_logunit,*) 'RDHM precipitation forcing directory :', &
                             rdhm356_struc_precip(n)%rdhm356dir

    !- Setting observed precip times to zero to ensure data is read in
    !   at first time step
       rdhm356_struc_precip(n)%rdhm356time1 = 0.0
       rdhm356_struc_precip(n)%rdhm356time2 = 0.0
    enddo
    call ESMF_ConfigFindLabel(LIS_config, "RDHM temperature forcing directory:", &
                              rc=rc)
    do n=1,LIS_rc%nnest   
       call ESMF_ConfigGetAttribute(LIS_config, &
                                    rdhm356_struc_temper(n)%rdhm356dir,rc=rc)
       write(LIS_logunit,*) 'RDHM temperature forcing directory :', &
                            rdhm356_struc_temper(n)%rdhm356dir

    !- Setting observed precip times to zero to ensure data is read in
    !   at first time step
       rdhm356_struc_temper(n)%rdhm356time1 = 0.0
       rdhm356_struc_temper(n)%rdhm356time2 = 0.0
    enddo

    call ESMF_ConfigFindLabel(LIS_config, "RDHM precipitation scale factor:",rc=rc)
    do n=1,LIS_rc%nnest   ! Loop over different nested domains
       call ESMF_ConfigGetAttribute(LIS_config, rdhm356_struc_precip(n)%scale,rc=rc)
       write(LIS_logunit,*) 'RDHM precipitation scale factor :', rdhm356_struc_precip(n)%scale
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RDHM precipitation interval:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_precip(n)%ts,rc=rc)
    enddo
  
    call ESMF_ConfigFindLabel(LIS_config,"RDHM temperature interval:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_temper(n)%ts,rc=rc)
    enddo
  
    call ESMF_ConfigFindLabel(LIS_config,"RDHM run window lower left hrap y:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_precip(n)%lower_left_hrapy,rc=rc)
  
    enddo
  
    call ESMF_ConfigFindLabel(LIS_config,"RDHM run window lower left hrap x:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_precip(n)%lower_left_hrapx,rc=rc)
    enddo
  
    call ESMF_ConfigFindLabel(LIS_config,"RDHM run window upper right hrap y:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_precip(n)%upper_right_hrapy,rc=rc)
  
    enddo
  
    call ESMF_ConfigFindLabel(LIS_config,"RDHM run window upper right hrap x:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_precip(n)%upper_right_hrapx,rc=rc)
    enddo
  
    call ESMF_ConfigFindLabel(LIS_config,"RDHM run window hrap resolution:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_precip(n)%hrap_resolution,rc=rc)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RDHM run window lower left hrap y:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_temper(n)%lower_left_hrapy,rc=rc)
  
    enddo
  
    call ESMF_ConfigFindLabel(LIS_config,"RDHM run window lower left hrap x:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_temper(n)%lower_left_hrapx,rc=rc)
    enddo
  
    call ESMF_ConfigFindLabel(LIS_config,"RDHM run window upper right hrap y:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_temper(n)%upper_right_hrapy,rc=rc)
    enddo
  
    call ESMF_ConfigFindLabel(LIS_config,"RDHM run window upper right hrap x:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_temper(n)%upper_right_hrapx,rc=rc)
    enddo
  
    call ESMF_ConfigFindLabel(LIS_config,"RDHM run window hrap resolution:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_temper(n)%hrap_resolution,rc=rc)
    enddo

    ! undefined value for temperature
    call ESMF_ConfigFindLabel(LIS_config,"RDHM temperature undefined value:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_temper(n)%undef_value,rc=rc)
    enddo
    
    ! undefined value for precipitation 
    call ESMF_ConfigFindLabel(LIS_config,"RDHM precipitation undefined value:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,rdhm356_struc_precip(n)%undef_value,rc=rc)
    enddo

    ! read the constant wind speed from LIS configuration file 
    call ESMF_ConfigFindLabel(LIS_config,"RDHM constant wind speed:",rc=rc)
    do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config, const_wind(n), rc=rc)
    enddo

end subroutine readcrd_rdhm356
