!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LDT_misc.h"
module LDT_obsSimMod
!BOP
!
! !MODULE: LDT_obsSimMod
! 
! !DESCRIPTION: 
!  The code in this file contains the basic datastructures and 
!  control routines for the observation simulator
!
! !REVISION HISTORY: 
!  03 Aug 2019    Sujay Kumar  Initial Specification
! 
! !USES:       

  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_ran2_gasdev
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_obsSimInit
  public :: LDT_readNatureRunData
  public :: LDT_temporalTransformObsSimData
  public :: LDT_applyObsSimMask
  public :: LDT_applyObsSimErrorModel
  public :: LDT_writeObsSim
  public :: LDT_logNatureRunData
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LDT_obsSim_struc

  type, public :: obsSimEntry

     character*50              :: natRunSource
     character*50              :: OSSEmaskSource
     integer                   :: nVars
     character*50, allocatable :: varNames(:)
     real,         allocatable :: varmins(:)
     real,         allocatable :: varmaxs(:)
     character*50              :: ttransform
     character*50              :: masktype
     character(len=LDT_CONST_PATH_LEN) :: maskdir
     character*50              :: errModelType
     character*50              :: errDist
     real                      :: errStdev
     integer                   :: seed(NRANDSEED)
     real, allocatable         :: datamask(:,:)
     integer, allocatable      :: count(:,:,:)
     real, allocatable         :: value(:,:,:) 
  end type obsSimEntry
  
  type(obsSimEntry) :: LDT_obsSim_struc

!EOP
contains

!BOP
! 
! !ROUTINE: LDT_obsSimInit
! \label{LDT_obsSimInit}
! 
! !INTERFACE:   
  subroutine LDT_obsSimInit
! !USES: 
    use LDT_coreMod
    use LDT_NatureRunData_pluginMod
    use LDT_OSSEmaskData_pluginMod

! 
! !DESCRIPTION: 
! 
!EOP
    implicit none

    integer :: i,rc,n
    integer :: c,r
    integer :: ftn
    integer :: maskid

    n = 1

!configurable options
    call ESMF_ConfigGetAttribute(LDT_config,LDT_obsSim_struc%natrunsource,&
         label="Observation simulator nature run source:", &
         rc=rc)
    call LDT_verify(rc,'Observation simulator nature run source: not specified')

    call ESMF_ConfigGetAttribute(LDT_config,LDT_obsSim_struc%OSSEmasksource,&
         label="Observation simulator OSSE mask source:", &
         rc=rc)
    call LDT_verify(rc,'Observation simulator OSSE mask source: not specified')

    call ESMF_ConfigGetAttribute(LDT_config,LDT_obsSim_struc%nVars,&
         label="Observation simulator number of simulated variables:", &
         rc=rc)
    call LDT_verify(rc,'Observation simulator number of simulated variables: not specified')

    allocate(LDT_obsSim_struc%varNames(LDT_obsSim_struc%nVars))
    allocate(LDT_obsSim_struc%varmins(LDT_obsSim_struc%nVars))
    allocate(LDT_obsSim_struc%varmaxs(LDT_obsSim_struc%nVars))

    call ESMF_ConfigFindLabel(LDT_config,&
         label="Observation simulator simulated variables:", &
         rc=rc)
    do i=1,LDT_obsSim_struc%nVars
       call ESMF_ConfigGetAttribute(LDT_config,&
            LDT_obsSim_struc%varNames(i),rc=rc)
       call LDT_verify(rc,'Observation simulator simulated variables: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config,&
         label="Observation simulator simulated variable minimum values:", &
         rc=rc)
    
    do i=1,LDT_obsSim_struc%nVars
       call ESMF_ConfigGetAttribute(LDT_config,&
            LDT_obsSim_struc%varmins(i),rc=rc)
       call LDT_verify(rc,'Observation simulator simulated variable minimum values: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config,&
         label="Observation simulator simulated variable maximum values:", &
         rc=rc)

    do i=1,LDT_obsSim_struc%nVars
       call ESMF_ConfigGetAttribute(LDT_config,&
            LDT_obsSim_struc%varMaxs(i),rc=rc)
       call LDT_verify(rc,'Observation simulator simulated variable maximum values: not defined')
    enddo
    

    call ESMF_ConfigGetAttribute(LDT_config,LDT_obsSim_struc%ttransform,&
         label="Observation simulator type of temporal transform:", &
         rc=rc)
    if(rc.ne.0) then
       write(LDT_logunit,*) "[ERR] Observation simulator type of temporal transform:' not specified"
       write(LDT_logunit,*) "[ERR] Supported options are "
       write(LDT_logunit,*) "[ERR] 'instantaneous', 'time-averaged'"
       write(LDT_logunit,*) "[ERR] Program stopping .."
       call LDT_endrun()
    endif

    call ESMF_ConfigGetAttribute(LDT_config,LDT_obsSim_struc%errDist,&
         label="Observation simulator error distribution type:", &
         rc=rc)
    if(rc.ne.0) then
       write(LDT_logunit,*) "[ERR] Observation simulator error distribution type:' not specified"
       write(LDT_logunit,*) "[ERR] Supported options are "
       write(LDT_logunit,*) "[ERR] 'gaussian'"
       write(LDT_logunit,*) "[ERR] Program stopping .."
       call LDT_endrun()
    endif

    call ESMF_ConfigGetAttribute(LDT_config,LDT_obsSim_struc%errModelType,&
         label="Observation simulator error model type:", &
         rc=rc)
    if(rc.ne.0) then
       write(LDT_logunit,*) "[ERR] Observation simulator error model type:' not specified"
       write(LDT_logunit,*) "[ERR] Supported options are "
       write(LDT_logunit,*) "[ERR] 'additive' or 'multiplicative'"
       write(LDT_logunit,*) "[ERR] Program stopping .."
       call LDT_endrun()
    endif

    call ESMF_ConfigGetAttribute(LDT_config,LDT_obsSim_struc%errStdev,&
         label="Observation simulator error standard deviation:", &
         rc=rc)
    call LDT_verify(rc, 'Observation simulator error standard deviation: not specified')

    call ESMF_ConfigGetAttribute(LDT_config,LDT_obsSim_struc%masktype,&
         label="Observation simulator masking model type:", &
         rc=rc)
    if(rc.ne.0) then
       write(LDT_logunit,*) "[ERR] Observation simulator masking model type:' not specified"
       write(LDT_logunit,*) "[ERR] Supported options are "
       write(LDT_logunit,*) "[ERR] 'none', 'constant', 'time-varying'"
       write(LDT_logunit,*) "[ERR] Program stopping .."
       call LDT_endrun()
    endif

    if(LDT_obsSim_struc%masktype.ne."none") then 
       call ESMF_ConfigGetAttribute(LDT_config,LDT_obsSim_struc%maskdir,&
            label="Observation simulator mask data directory:", &
            rc=rc)
       call LDT_verify(rc, 'Observation simulator mask data directory: not specified')

       allocate(LDT_obsSim_struc%datamask(LDT_rc%lnc(n),LDT_rc%lnr(n)))
!
! File is expected to be in NetCDF format with a field "MASK" in it    
! The dimensions of the mask file should be the same as LDT grid  
! mask fields are expected to have values of 1 (valid) and 0 (invalid)    

       if(LDT_obsSim_struc%masktype.eq."constant") then 

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
          call LDT_verify(nf90_open(trim(LDT_obsSim_struc%maskdir),&
               NF90_NOWRITE,ftn),&
               'Failed to open the mask file '//trim(LDT_obsSim_struc%maskdir))

          call LDT_verify(nf90_inq_varid(ftn,'MASK',maskid),&
               'nf90_inq_varid failed for MASK')
          
          call LDT_verify(nf90_get_var(ftn,maskid,&
               LDT_obsSim_struc%datamask),&
               'nf90_get_var failed for MASK')
          call LDT_verify(nf90_close(ftn),&
               'nf90_close failed for '//trim(LDT_obsSim_struc%maskdir))
#endif
       endif      
    endif

    call LDT_NatureRunData_plugin
    call LDT_OSSEmaskData_plugin  

    call setupnaturerunsource(trim(LDT_obsSim_struc%natRunSource)//char(0))

    if(LDT_obsSim_struc%OSSEmasksource.ne."none") then 
       call setupossemasksource(trim(LDT_obsSim_struc%OSSEmasksource)//char(0))
    endif

    allocate(LDT_obsSim_struc%value(LDT_rc%lnc(n),LDT_rc%lnr(n),&
         LDT_obsSim_struc%nVars))
    allocate(LDT_obsSim_struc%count(LDT_rc%lnc(n),LDT_rc%lnr(n),&
         LDT_obsSim_struc%nVars))

    LDT_obsSim_struc%value = 0 
    LDT_obsSim_struc%count = 0 

    LDT_obsSim_struc%seed = -1000
    
  end subroutine LDT_obsSimInit

!BOP
! 
! !ROUTINE: LDT_readNatureRunData
! \label{LDT_readNatureRunData}
! 
! !INTERFACE: 
  subroutine LDT_readNatureRunData(n)
! !USES: 
    use LDT_coreMod,   only : LDT_rc

    implicit none
    
    integer              :: n 
! 
! !DESCRIPTION: 
! 
!  This subroutine reads the Nature run output and subsamples them
! 
!EOP

    call readNatureRunSource(trim(LDT_obsSim_struc%natRunSource)//char(0),&
         n)

  end subroutine LDT_readNatureRunData


!BOP
! 
! !ROUTINE: LDT_temporalTransformObsSimData
! \label{LDT_temporalTransformObsSimData}
! 
! !INTERFACE: 
  subroutine LDT_temporalTransformObsSimData(n)
! !USES: 
    use LDT_coreMod,   only : LDT_rc

    implicit none
    
    integer              :: n 
! 
! !DESCRIPTION: 
! 
!  This subroutine performs the required temporal transformations to the
!  simulated observations (averaging, for e.g.)
! 
!  if the temporal tranform option is instantaneous, then the timestep
!  should be set to the sampling interval. 
! 
!EOP

    integer                :: c,r,k

    if(LDT_obsSim_struc%ttransform.eq."instantaneous") then 

    elseif(LDT_obsSim_struc%ttransform.eq."time-averaged") then
       
       if(mod(float(LDT_rc%hr)*3600+60*float(LDT_rc%mn)+float(LDT_rc%ss),&
            LDT_rc%tavgInterval).eq.0) then   
          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(LDT_obsSim_struc%count(c,r,k).gt.0) then 
                   LDT_obsSim_struc%value(c,r,k) = &
                        LDT_obsSim_struc%value(c,r,k)/&
                        LDT_obsSim_struc%count(c,r,k) 
                   
                endif
             enddo
          enddo
       endif
    endif
  end subroutine LDT_temporalTransformObsSimData


!BOP
! 
! !ROUTINE: LDT_applyObsSimMask
! \label{LDT_applyObsSimMask}
! 
! !INTERFACE: 
  subroutine LDT_applyObsSimMask(n)
! !USES: 
    use LDT_coreMod,   only : LDT_rc

    implicit none
    
    integer              :: n 
! 
! !DESCRIPTION: 
! 
!  This subroutine applies spatial masks to the simulated
!  observations. The masking could be static in time, 
!  or time varying. If the mask is temporally variable, 
!  then it could be supplied from an external data or 
!  simulated within the code. 
! 
!EOP
    integer              :: c,r,k

    if(LDT_obsSim_struc%masktype.ne.'none') then 

       if(LDT_obsSim_struc%masktype.eq."constant") then 
          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(LDT_obsSim_struc%datamask(c,r).eq.0.0) then
                   LDT_obsSim_struc%value(c,r,:) = -9999.0
                endif
             enddo
          enddo

       elseif(LDT_obsSim_struc%masktype.eq."time-varying") then 
          
          call readOSSEmasksource(trim(LDT_obsSim_struc%OSSEmasksource)//char(0),&
               n)

          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(LDT_obsSim_struc%datamask(c,r).eq.0.0) then
                   LDT_obsSim_struc%value(c,r,:) = -9999.0
                endif
             enddo
          enddo
          
       endif
    endif

  end subroutine LDT_applyObsSimMask


!BOP
! 
! !ROUTINE: LDT_applyObsSimErrorModel
! \label{LDT_applyObsSimErrorModel}
! 
! !INTERFACE: 
  subroutine LDT_applyObsSimErrorModel(n)
! !USES: 

    implicit none
    
    integer              :: n 
! 
! !DESCRIPTION: 
! 
!  This subroutine applies a user specified error model
!  to the simulated observations
! 
!EOP
    integer              :: c,r,k
    real                 :: tmp_val
    real                 :: rand(2)

    if(LDT_obsSim_struc%errDist.eq."gaussian") then 
       if(LDT_obsSim_struc%errModelType.eq."additive") then 

          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                do k=1,LDT_obsSim_struc%nVars
                   if(LDT_obsSim_struc%value(c,r,k).ne.-9999.0) then 
                      call nr_gasdev(LDT_obsSim_struc%seed, rand)
                      tmp_val = rand(1)*LDT_obsSim_struc%errStdev
                      if (LDT_obsSim_struc%value(c,r,k).gt.&
                           LDT_obsSim_struc%varmins(k).and.&
                           LDT_obsSim_struc%value(c,r,k) + tmp_val.gt.&
                            LDT_obsSim_struc%varmins(k)) then 
                         LDT_obsSim_struc%value(c,r,k) = & 
                              LDT_obsSim_struc%value(c,r,k) + tmp_val
                      endif
                   endif
                enddo
             enddo
          enddo

       elseif(LDT_obsSim_struc%errModelType.eq."multiplicative") then 
          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                do k=1,LDT_obsSim_struc%nVars
                   if(LDT_obsSim_struc%value(c,r,k).ne.-9999.0) then 
                      call nr_gasdev(LDT_obsSim_struc%seed, rand)
                      tmp_val = rand(1)*LDT_obsSim_struc%errStdev

                      if (LDT_obsSim_struc%value(c,r,k)*tmp_val.gt.&
                           LDT_obsSim_struc%varmins(k)) then 
                         LDT_obsSim_struc%value(c,r,k) = & 
                              LDT_obsSim_struc%value(c,r,k)*tmp_val
                      endif
                   endif
                enddo
             enddo
          enddo

       endif
    endif
  end subroutine LDT_applyObsSimErrorModel


!BOP
! 
! !ROUTINE: LDT_writeObsSim
! \label{LDT_writeObsSim}
! 
! !INTERFACE: 
  subroutine LDT_writeObsSim(n)
! !USES: 
    use LDT_coreMod,   only : LDT_rc

    implicit none
    
    integer              :: n 
! 
! !DESCRIPTION: 
! 
!  This subroutine writes the simulated observations
!  to an external file. 
! 
!EOP
    character(len=LDT_CONST_PATH_LEN)           :: dname
    character(len=LDT_CONST_PATH_LEN)           :: filename
    character(len=10)       :: cdate
    character(len=14)       :: cdate1
    integer                 :: ftn
    integer                 :: c,r,k
    integer                 :: latid, lonid,dimId(2)
    integer                 :: varid(LDT_obsSim_struc%nVars)
    real                    :: lat(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real                    :: lon(LDT_rc%lnc(n),LDT_rc%lnr(n))
    integer                 :: rc
    integer, external       :: LDT_create_subdirs


    write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
         LDT_rc%yr, LDT_rc%mo, &
         LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
    dname = trim(LDT_rc%odir)//'/'
    
    write(unit=cdate, fmt='(i4.4, i2.2)') LDT_rc%yr, LDT_rc%mo
    dname = trim(dname)//trim(cdate)//'/'
    
!    call system('mkdir -p '//trim(dname))
    rc = LDT_create_subdirs(len_trim(dname), trim(dname))
    if (rc .ne. 0) then
      write(LDT_logunit,*)'[ERR] Cannot create directory ', trim(dname)
      call LDT_endrun()
    end if
    filename = trim(dname)//'SimObs_'//trim(cdate1)//'.nc'

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    call LDT_verify(nf90_create(filename,NF90_CLOBBER, &
         ftn),'failed to create file '//trim(filename))
    call LDT_verify(nf90_def_dim(ftn,"east_west",LDT_rc%lnc(n),dimID(1)),&
         'nf90_def_dim failed for east_west')
    call LDT_verify(nf90_def_dim(ftn,"north_south",LDT_rc%lnr(n),dimID(2)),&
         'nf90_def_dim failed for north_south')
    call LDT_verify(nf90_def_var(ftn,'lat',NF90_FLOAT, &
         dimID(1:2), latid),&
         'nf90_def_var failed for lat')
    call LDT_verify(nf90_def_var(ftn,'lon',NF90_FLOAT, &
         dimID(1:2), lonid),&
         'nf90_def_var failed for lon')

    do k=1,LDT_obsSim_struc%nVars
       call LDT_verify(nf90_def_var(ftn,trim(LDT_obsSim_struc%varNames(k)),&
            NF90_FLOAT, &
            dimID(1:2), varid(k)),&
            'nf90_def_var failed for '//trim(LDT_obsSim_struc%varNames(k)))
    enddo
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value", -9999.0))    
    call LDT_verify(nf90_enddef(ftn),'nf90_enddef failed ')

    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)
          lat(c,r) = LDT_domain(n)%lat(c+(r-1)*LDT_rc%lnc(n))
          lon(c,r) = LDT_domain(n)%lon(c+(r-1)*LDT_rc%lnc(n))
       enddo
    enddo
    
    call LDT_verify(nf90_put_var(ftn,latid,&
         lat),'nf90_put_var failed for lat')
    call LDT_verify(nf90_put_var(ftn,lonid,&
         lon),'nf90_put_var failed for lon')
    do k=1,LDT_obsSim_struc%nVars
       call LDT_verify(nf90_put_var(ftn,varid(k),&
            LDT_obsSim_struc%value(:,:,k)),&
            'nf90_put_var failed for '//trim(LDT_obsSim_struc%varNames(k)))

    enddo

    call LDT_verify(nf90_close(ftn),&
         'nf90_close failed for '//trim(filename))

#endif
  end subroutine LDT_writeObsSim

!BOP
!
! !ROUTINE: LDT_logNatureRunData
! \label{LDT_logNatureRunData}
! 
! !INTERFACE:
  subroutine LDT_logNatureRunData(n, value,vlevel)
! !USES: 
    use LDT_coreMod, only     : LDT_rc, LDT_domain
    use LDT_logMod,  only     : LDT_logunit, LDT_endrun

    implicit none
! !ARGUMENTS: 
    integer, intent(in)        :: n 
    real                       :: value(LDT_rc%lnc(n), LDT_rc%lnr(n))
    integer, optional          :: vlevel
! 
! !DESCRIPTION: 
!  This subroutine maps the processed observations onto the LDT data
!  structures for future temporal averaging and comparisons. The 
!  data are also filtered using the specified external mask. 
!
!EOP
    integer :: form
    integer :: k,i,c,r,gid

    if(.not.present(vlevel)) then 
       k = 1
    else
       k = vlevel
    endif

    if(LDT_obsSim_struc%ttransform.eq."instantaneous") then 
       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(value(c,r).ne.LDT_rc%udef) then 
                LDT_obsSim_struc%value(c,r,k) = value(c,r)
                LDT_obsSim_struc%count(c,r,k) = 1
             endif
          enddo
       enddo
       
    elseif(LDT_obsSim_struc%ttransform.eq."time-averaged") then 
       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(value(c,r).ne.LDT_rc%udef) then 
                LDT_obsSim_struc%value(c,r,k) = &
                     LDT_obsSim_struc%value(c,r,k) + value(c,r)
                LDT_obsSim_struc%count(c,r,k) = &
                     LDT_obsSim_struc%count(c,r,k) + 1
             endif
          enddo
       enddo
    endif
  end subroutine LDT_logNatureRunData

 end module LDT_obsSimMod
