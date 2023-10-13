!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_DApreprocMod
!BOP
!
! !MODULE: LDT_DApreprocMod
! 
! !DESCRIPTION: 
! 
! !REVISION HISTORY: 
!  24 Nov 2008    Sujay Kumar  Initial Specification
!  2 Dec 2021:   Mahdi Navari; modified to save stratify CDF
! 
  use ESMF

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_DApreprocInit          !initialize specified domains

!EOP

contains
!BOP
! !ROUTINE: LDT_DApreprocInit
! \label{LDT_DApreprocInit}
!
! !INTERFACE: 
  subroutine LDT_DApreprocInit()
! !USES:
    use LDT_coreMod, only : LDT_rc, LDT_config,LDT_domain
    use LDT_logMod
    use LDT_timeMgrMod, only : LDT_parseTimeString
    use LDT_DAobs_pluginMod, only : LDT_DAobs_plugin

    implicit none
!
! !DESCRIPTION: 
!
!   This subroutine initializes the datastructures 
! 
!EOP
    character*20         :: stime
    character*20         :: tres
    integer              :: rc
    integer              :: n
    integer              :: c,r,j
    integer              :: ios1
    integer              :: ftn
    character*50         :: preprocMethod
    real                 :: delta
    real, allocatable    :: cdf_strat_data(:,:)
    real, allocatable    :: stratification_data(:,:,:)
    n = 1
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%obs_src,&
         label="DA observation source:",rc=rc)
    call LDT_verify(rc,'DA observation source: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,preprocMethod,&
         label="DA preprocessing method:",rc=rc)
    call LDT_verify(rc,'DA preprocessing method: not defined')
!   preprocmethod - "Obs grid generation" "CDF generation", "Anomaly correction"

    LDT_rc%comp_obsgrid      = 0 
    LDT_rc%comp_cdf          = 0
    LDT_rc%anomalyObsProc    = 0 

    if(preprocMethod.eq."Obs grid generation") then     
       LDT_rc%comp_obsgrid = 1
    elseif(preprocMethod.eq."CDF generation") then 
       LDT_rc%comp_cdf = 1
    elseif(preprocMethod.eq."Anomaly correction") then 
       LDT_rc%anomalyObsProc = 1
    endif
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%dapreprocfile,&
         label="Name of the preprocessed DA file:",rc=rc)
    call LDT_verify(rc,'Name of the preprocessed DA file: not defined')

    if(LDT_rc%comp_cdf.gt.0) then 

       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%cdf_nbins,&
            label="Number of bins to use in the CDF:",rc=rc)
       call LDT_verify(rc,'Number of bins to use in the CDF: not defined')


       call ESMF_ConfigGetAttribute(LDT_config,tres,&
            label="Temporal resolution of CDFs:",rc=rc)
       call LDT_verify(rc,'Temporal resolution of CDFs: not defined')    
       
       if(tres.eq."yearly") then 
          LDT_rc%cdf_ntimes = 1
       elseif(tres.eq."monthly") then 
          LDT_rc%cdf_ntimes = 12
       elseif(tres.eq."half-monthly") then    !for 9-km operational SMAP (Y.Kwon)
          LDT_rc%cdf_ntimes = 24
       endif
       
       !---------------------------------for 9-km operational SMAP (Y.Kwon)
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%daily_interp_switch,&
            label="Daily interpolation of mean and stddev:",&
            default=0, rc=rc)
       call LDT_verify(rc,'Daily interpolation of mean and stddev: not defined')

       if(LDT_rc%daily_interp_switch.eq.1) then
          LDT_rc%cdf_ntimes = 365
       endif
       !--------------------------------------------------------------------

       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%sp_sampl_cdfs,&
            label="Enable spatial sampling for CDF calculations:",&
            default=0, rc=rc)
       call LDT_verify(rc,&
            'Enable spatial sampling for CDF calculations: not defined')

       if(LDT_rc%sp_sampl_cdfs.gt.0) then 
          
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%sp_sample_cdf_rad,&
               label="Spatial sampling window radius for CDF calculations:",&
               rc=rc)
          call LDT_verify(rc,&
               'Spatial sampling window radius for CDF calculations: not defined')
       endif       

       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%group_cdfs,&
            label="Group CDFs by external data:",default=0, rc=rc)
       call LDT_verify(rc,'Group CDFs by external data: not defined')           

       if(LDT_rc%group_cdfs.gt.0) then 

          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_rc%group_cdfs_attrib_file,&
               label="CDF grouping attributes file:",rc=rc)
          call LDT_verify(rc,'CDF grouping attributes file: not defined')
          
          ftn = LDT_getNextUnitNumber()
          open(ftn,file=LDT_rc%group_cdfs_attrib_file, status='old')
          read(ftn,*)
          read(ftn,*) LDT_rc%group_cdfs_strat_file
          read(ftn,*) 
          read(ftn,*) LDT_rc%group_cdfs_min, LDT_rc%group_cdfs_max, &
               LDT_rc%group_cdfs_nbins
          
          call LDT_releaseUnitNumber(ftn)

          allocate(cdf_strat_data(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          allocate(LDT_rc%cdf_strat_data(LDT_rc%ngrid(n)))
          
          write(LDT_logunit,*) '[INFO] Reading ',trim(LDT_rc%group_cdfs_strat_file)
          ftn = LDT_getNextUnitNumber()
          open(ftn,file=(LDT_rc%group_cdfs_strat_file),&
               form='unformatted',access='direct',&
               recl=LDT_rc%lnc(n)*LDT_rc%lnr(n)*4,iostat=ios1)
          read(ftn,rec=1) cdf_strat_data
          
          delta = (LDT_rc%group_cdfs_max-LDT_rc%group_cdfs_min)/&
               LDT_rc%group_cdfs_nbins

          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(LDT_domain(n)%gindex(c,r).ne.-1) then

                    LDT_rc%cdf_strat_data(LDT_domain(n)%gindex(c,r)) = & 
                        nint((cdf_strat_data(c,r) - LDT_rc%group_cdfs_min)/&
                        delta)+1
                    if (LDT_rc%cdf_strat_data(LDT_domain(n)%gindex(c,r)) .gt. LDT_rc%group_cdfs_nbins) then
                       write(LDT_logunit,*) '[INFO] Group bins is larger then Max Group bins',&
                            LDT_rc%cdf_strat_data(LDT_domain(n)%gindex(c,r)), 'vs.', LDT_rc%group_cdfs_nbins ,&
                            'Value adjusted the Max Group bins'
                       LDT_rc%cdf_strat_data(LDT_domain(n)%gindex(c,r)) = LDT_rc%group_cdfs_nbins
                    endif
                endif
             enddo
          enddo

          call LDT_releaseUnitNumber(ftn)
          write(LDT_logunit,*) '[INFO] Finished reading ',&
               trim(LDT_rc%group_cdfs_strat_file)
       endif

!This part reads the monthly total precipitation climatology and generates
!  stratification input data LDT_rc%stratification_data(LDT_rc%ngrid(n),LDT_rc%cdf_ntimes).

       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%strat_cdfs,&
            label="Stratify CDFs by external data:",default=0, rc=rc)
       call LDT_verify(rc,"Stratify CDFs by external data: not defined")

       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%write_strat_cdfs,&
            label="Write stratified geolocation independent CDFs:",default=0, rc=rc)
       call LDT_verify(rc,"Write stratify geolocation independent CDFs:: not defined")

       if(LDT_rc%strat_cdfs.gt.0) then
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%strat_src,&
               label="Stratification data source:", rc=rc)
          call LDT_verify(rc,"Stratification data source: not defined")

          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%strat_cdfs_nbins,&
               label="Number of bins to use for stratification:",rc=rc)
          call LDT_verify(rc,"Number of bins to use for stratification: not defined")

          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%strat_file,&
               label="External stratification file:",rc=rc)
          call LDT_verify(rc,"External stratification file: not defined")

          allocate(LDT_rc%stratification_data(LDT_rc%ngrid(n),LDT_rc%cdf_ntimes))
          allocate(stratification_data(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%cdf_ntimes))

          call read_Precip_climo (LDT_rc%lnc(n), LDT_rc%lnr(n), LDT_rc%strat_file, stratification_data) !

          do j=1,LDT_rc%cdf_ntimes
             LDT_rc%strat_cdfs_min = 0. !minval returns -9999.  minval(stratification_data(:,:,j)) ! min value over the entire domain  
             LDT_rc%strat_cdfs_max = maxval(stratification_data(:,:,j)) ! max value over the entire domain 
             delta = (LDT_rc%strat_cdfs_max-LDT_rc%strat_cdfs_min)/&
                  LDT_rc%strat_cdfs_nbins               
             do r=1,LDT_rc%lnr(n)
                do c=1,LDT_rc%lnc(n)
                   if(LDT_domain(n)%gindex(c,r).ne.-1) then
                      if(stratification_data(c,r,j) .gt. 0) then 
                         LDT_rc%stratification_data(LDT_domain(n)%gindex(c,r), j) = &
                              nint((stratification_data(c,r,j) - LDT_rc%strat_cdfs_min)/&
                              delta)+1
                         if (LDT_rc%stratification_data(LDT_domain(n)%gindex(c,r), j) .gt. LDT_rc%strat_cdfs_nbins) then
                            write(LDT_logunit,*) '[INFO] Startification bins is larger then Max Startification bins',&
                                 LDT_rc%stratification_data(LDT_domain(n)%gindex(c,r), j) , 'vs.', LDT_rc%strat_cdfs_nbins,&
                                 'Value adjusted to the Max Startification bins'
                            LDT_rc%stratification_data(LDT_domain(n)%gindex(c,r), j) = LDT_rc%strat_cdfs_nbins
                         endif
                      else
                         write(LDT_logunit,*) '[INFO] Total precip is zero or undefined',&
                              stratification_data(c,r,j) ,c,r,j,&
                              'Value adjusted to the Min Startification bins'
                         LDT_rc%stratification_data(LDT_domain(n)%gindex(c,r), j) = 1
                      endif
                   endif
                enddo
             enddo
          enddo
       endif
    endif

    if(LDT_rc%comp_obsgrid.eq.1) then
       LDT_rc%pass = 0
    else

       if(LDT_rc%tavgInterval.lt.LDT_rc%ts) then
          write(LDT_logunit,*) '[ERR] Temporal averaging interval must be greater than'
          write(LDT_logunit,*) '[ERR] or equal to the LIS output timestep. '
          write(LDT_logunit,*) '[ERR] Program stopping....'
          call LDT_endrun()
       endif

       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%obsCountThreshold,&
            label="Observation count threshold:",&
            rc=rc)
       call LDT_verify(rc,'Observation count threshold: not defined')

       if(LDT_rc%comp_cdf.eq.1) then
          LDT_rc%pass = 2
       elseif(LDT_rc%anomalyObsProc.eq.1) then
          LDT_rc%pass = 2
       endif

    endif

    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%applyMask,&
         label="Apply external mask:",&
         rc=rc)
    call LDT_verify(rc,'Apply external mask: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%maskdir,&
         label="External mask directory:",&
         rc=rc)
    call LDT_verify(rc,'External mask directory: not defined')

    call LDT_DAobs_plugin

    LDT_rc%nobs = 1

! This flag is set to prompt the writing of the domain file
    if(LDT_rc%anomalyObsProc.eq.1) then
       LDT_rc%comp_obsgrid = 1
    endif

  end subroutine LDT_DApreprocInit

end module LDT_DApreprocMod

