!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
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
    integer              :: c,r
    integer              :: ios1
    integer              :: ftn
    character*50         :: preprocMethod
    real                 :: delta
    real, allocatable    :: cdf_strat_data(:,:)

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
       endif
       
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
                endif
             enddo
          enddo

          call LDT_releaseUnitNumber(ftn)
          write(LDT_logunit,*) '[INFO] Finished reading ',&
               trim(LDT_rc%group_cdfs_strat_file)
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

