!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: LVT_DAMod
! \label(LVT_DAMod)
!
! !INTERFACE:
module LVT_DAMod
! 
! !USES:   

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif
  use LVT_coreMod
  use LVT_logMod
  use LVT_statsDataMod

  implicit none
  PRIVATE
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This file provides routines to extract statistics from the LIS analysis
!  output
!     - mean of normalized innovation distribution
!     - variance of normalized innnovation distribution
!  TODO: It is assumed that LIS domain and LVT domain are exactly the same
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  24 Nov 2008    Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_DAInit          !initialize specified domains
  public :: LVT_readLISDAdata
  public :: LVT_computeDAstats

  public :: LVT_DAstats


  type, public  :: dastatsdec
     real, allocatable :: mean_norm_innov(:,:)
     real, allocatable :: var_norm_innov(:,:)
     real, allocatable :: sx_norm_innov(:,:)
     real, allocatable :: sxx_norm_innov(:,:)

     integer, allocatable :: count_norm_innov(:,:)

     real, allocatable :: mean_innov(:,:)
     integer, allocatable :: count_mean_innov(:,:)

     real, allocatable :: mean_fvar(:,:)
     integer, allocatable :: count_mean_fvar(:,:)

     real,    allocatable :: mean_gain(:,:)
     integer, allocatable :: count_mean_gain(:,:)

     real,    allocatable :: mean_spread(:,:)
     integer, allocatable :: count_mean_spread(:,:)

     real,    allocatable :: mean_spread_ts(:,:)
     integer, allocatable :: count_mean_spread_ts(:,:)
  end type dastatsdec
  
  type(dastatsdec) :: LVT_DAstats
!EOP

contains
!BOP
! 
! !ROUTINE: LVT_DAInit
! \label{LVT_DAInit}
!
! !INTERFACE: 
  subroutine LVT_DAInit()
! 
! !USES:

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!   This subroutine initializes the datastructures 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
#if 0     
    if(LVT_rc%computeInnovDist.eq.1) then 
       allocate(LVT_DAstats%mean_norm_innov(LVT_LIS_rc(1)%lnc, &
            LVT_LIS_rc(1)%lnr))
       allocate(LVT_DAstats%var_norm_innov(LVT_LIS_rc(1)%lnc,&
            LVT_LIS_rc(1)%lnr))   
       allocate(LVT_DAstats%sx_norm_innov(LVT_LIS_rc(1)%lnc,&
            LVT_LIS_rc(1)%lnr))
       allocate(LVT_DAstats%sxx_norm_innov(LVT_LIS_rc(1)%lnc,&
            LVT_LIS_rc(1)%lnr))   
       allocate(LVT_DAstats%count_norm_innov(LVT_LIS_rc(1)%lnc, &
            LVT_LIS_rc(1)%lnr))
       
       allocate(LVT_DAstats%mean_innov(LVT_LIS_rc(1)%lnc, &
            LVT_LIS_rc(1)%lnr))
       allocate(LVT_DAstats%count_mean_innov(LVT_LIS_rc(1)%lnc, &
            LVT_LIS_rc(1)%lnr))

       allocate(LVT_DAstats%mean_fvar(LVT_LIS_rc(1)%lnc,&
            LVT_LIS_rc(1)%lnr))
       allocate(LVT_DAstats%count_mean_fvar(LVT_LIS_rc(1)%lnc, &
            LVT_LIS_rc(1)%lnr))

       LVT_DAstats%mean_norm_innov = 0 
       LVT_DAstats%var_norm_innov = 0 
       LVT_DAstats%sx_norm_innov = 0 
       LVT_DAstats%sxx_norm_innov = 0 
       LVT_DAstats%count_norm_innov = 0 

       LVT_DAstats%mean_innov = 0
       LVT_DAstats%count_mean_innov = 0 

       LVT_DAstats%mean_fvar = 0
       LVT_DAstats%count_mean_fvar = 0 

    endif
    if(LVT_rc%computeGain.eq.1) then 
       allocate(LVT_DAstats%mean_gain(LVT_LIS_rc(1)%ntiles, &
            LVT_rc%nstvars))   
       allocate(LVT_DAstats%count_mean_gain(LVT_LIS_rc(1)%ntiles, &
            LVT_rc%nstvars))   
       
       LVT_DAstats%mean_gain = 0.0
       LVT_DAstats%count_mean_gain = 0
    endif

    if(LVT_rc%computeSpread.eq.1) then 
       allocate(LVT_DAstats%mean_spread(LVT_LIS_rc(1)%ntiles, &
            LVT_rc%nstvars))   
       allocate(LVT_DAstats%count_mean_spread(LVT_LIS_rc(1)%ntiles, &
            LVT_rc%nstvars))   
       
       LVT_DAstats%mean_spread = 0.0
       LVT_DAstats%count_mean_spread = 0

       allocate(LVT_DAstats%mean_spread_ts(LVT_LIS_rc(1)%ntiles, &
            LVT_rc%nstvars))   
       allocate(LVT_DAstats%count_mean_spread_ts(LVT_LIS_rc(1)%ntiles, &
            LVT_rc%nstvars))   
       
       LVT_DAstats%mean_spread_ts = 0.0
       LVT_DAstats%count_mean_spread_ts = 0
    endif
#endif       
  end subroutine LVT_DAInit

!BOP
! 
! !ROUTINE: LVT_readLISDAdata
! \label{LVT_readLISDAdata}
!
! !INTERFACE: 
  subroutine LVT_readLISDAdata(pass)
! 
! !USES:     
   
    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in) :: pass
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This subroutine reads the Innovation files from the LIS assimilation run
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
#if 0 
    character*100 :: innovfile, gainfile, spreadfile
    logical       :: file_exists
    integer       :: ftn
    integer       :: t
    integer       :: ninnovid, innovid, fsigmaid
    real          :: norm_innov(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr)
    real          :: fvar(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr)
    real          :: innov(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr)
    real          :: gain(LVT_LIS_rc(1)%ntiles,LVT_rc%nstvars)
    real          :: spread(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr),&
         LVT_rc%nstvars)

    if(LVT_rc%computeInnovDist.eq.1) then 
       call create_dainnov_filename(innovfile)
       inquire(file=innovfile, exist=file_exists) 
       
       if(file_exists) then 
          write(LVT_logunit,*) '[INFO] Reading Innovation file ',innovfile
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
          call LVT_verify(nf90_open(path=trim(innovfile),&
               mode=NF90_NOWRITE, ncid=ftn),&
               'failed to open innovfile '//trim(innovfile))

          call LVT_verify(nf90_inq_varid(ftn,"ninnov_01",ninnovid),&
               'nf90_inq_varid failed for ninnov_01')
          call LVT_verify(nf90_inq_varid(ftn,"innov_01",innovid),&
               'nf90_inq_varid failed for innov_01')
          call LVT_verify(nf90_inq_varid(ftn,"forecast_sigma_01",fsigmaid),&
               'nf90_inq_varid failed for forecast_sigma')
          
      
          call LVT_verify(nf90_get_var(ftn,ninnovid,norm_innov),&
               'nf90_get_var failed for ninnov')
          call LVT_verify(nf90_get_var(ftn,innovid,innov),&
               'nf90_get_var failed for innov')
          call LVT_verify(nf90_get_var(ftn,fsigmaid,fvar),&
               'nf90_get_var failed for forecast sigma')

          call LVT_verify(nf90_close(ftn),'failed to close '//&
               trim(innovfile))
          
#endif
!binary old output
#if 0 
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=innovfile, form='unformatted')
          
          read(ftn) norm_innov
          read(ftn) innov
          read(ftn) fvar
          call LVT_releaseUnitNumber(ftn)
#endif
          
          call logNormInnovationEntry(pass, norm_innov)
!          call logInnovationEntry(pass, innov)
!          call logFvarEntry(pass, fvar)
       endif
    endif

    if(LVT_rc%computeGain.eq.1) then 

       call create_dagain_filename(gainfile)
       inquire(file=gainfile, exist=file_exists) 

       if(file_exists) then 
          ftn = LVT_getNextUnitNumber()
          write(LVT_logunit,*) '[INFO] Reading Gain file ',gainfile
          open(ftn,file=gainfile, form='unformatted')
          
          do t=1,LVT_rc%nstvars
             read(ftn) gain(:,t)
          enddo
          call LVT_releaseUnitNumber(ftn)
          
          call logGainEntry(pass, gain)
       endif
    endif

    if(LVT_rc%computeSpread.eq.1) then 

       call create_daspread_filename(spreadfile)
       inquire(file=spreadfile, exist=file_exists) 

       if(file_exists) then 
          ftn = LVT_getNextUnitNumber()
          write(LVT_logunit,*) '[INFO] Reading Spread file ',spreadfile
          open(ftn,file=spreadfile, form='unformatted')
          
          do t=1,LVT_rc%nstvars
             read(ftn) spread(:,:,t)
          enddo
          call LVT_releaseUnitNumber(ftn)
          
          call logSpreadEntry(pass, spread)
       endif
    endif
#endif
  end subroutine LVT_readLISDAdata

#if 0 
!BOP
! 
! !ROUTINE: logInnovationEntry
! \label{logInnovationEntry}
!
! !INTERFACE:
  subroutine logInnovationEntry(pass, innov)
! 
! !USES: 

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine updates the analysis (mean, variance) computations
!  on the normalized innovations data
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
! !ARGUMENTS: 
    integer  :: pass
    real     :: innov(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr)
!EOP
    integer  :: c,r
    
    if(pass.eq.1) then 
       do r=1,LVT_LIS_rc(1)%lnr
          do c=1,LVT_LIS_rc(1)%lnc
             if(LVT_stats%datamask(c,r).eq.1) then 
                if(innov(c,r).ne.LVT_rc%udef) then 
                   LVT_DAstats%mean_innov(c,r) = &
                        LVT_DAstats%mean_innov(c,r)+innov(c,r)
                   LVT_DAstats%count_mean_innov(c,r) = &
                        LVT_DAstats%count_mean_innov(c,r)+1

                endif
             endif
          enddo
       enddo
    endif
    
  end subroutine logInnovationEntry

!BOP
! 
! !ROUTINE: logFvarEntry
! \label{logFvarEntry}
!
! !INTERFACE:
  subroutine logFvarEntry(pass, fvar)
! 
! !USES: 


    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine updates the analysis (mean, variance) computations
!  on the normalized innovations data
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    integer  :: pass
    real     :: fvar(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr)
!EOP
    integer  :: c,r
    
    if(pass.eq.1) then 
       do r=1,LVT_LIS_rc(1)%lnr
          do c=1,LVT_LIS_rc(1)%lnc
             if(LVT_stats%datamask(c,r).eq.1) then 
                if(fvar(c,r).ne.LVT_rc%udef) then 
                   LVT_DAstats%mean_fvar(c,r) = &
                        LVT_DAstats%mean_fvar(c,r)+fvar(c,r)**2
                   LVT_DAstats%count_mean_fvar(c,r) =&
                        LVT_DAstats%count_mean_fvar(c,r)+1
                endif
             endif
          enddo
       enddo
    endif
    
  end subroutine logFvarEntry

!BOP
! 
! !ROUTINE: logNormInnovationEntry
! \label{logNormInnovationEntry}
!
! !INTERFACE:
  subroutine logNormInnovationEntry(pass, norm_innov)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine updates the analysis (mean, variance) computations
!  on the normalized innovations data
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !USES: 
! !ARGUMENTS: 
    integer  :: pass
    real     :: norm_innov(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr)
!EOP
    integer  :: c,r
    
    if(pass.eq.1) then 
       do r=1,LVT_LIS_rc(1)%lnr
          do c=1,LVT_LIS_rc(1)%lnc
             if(LVT_stats%datamask(c,r).eq.1) then 
                if(norm_innov(c,r).ne.LVT_rc%udef) then 
                   LVT_DAstats%sx_norm_innov(c,r) = &
                        LVT_DAstats%sx_norm_innov(c,r)+norm_innov(c,r)

                   LVT_DAstats%sxx_norm_innov(c,r) = &
                        LVT_DAstats%sxx_norm_innov(c,r)+&
                        norm_innov(c,r)*norm_innov(c,r)

                   LVT_DAstats%count_norm_innov(c,r) = &
                        LVT_DAstats%count_norm_innov(c,r)+1

!                   if(c.eq.1.and.r.eq.1) then 
!                      print*, norm_innov(c,r)
!                   endif
                endif
             endif
          enddo
       enddo       
    endif
    
  end subroutine logNormInnovationEntry

!BOP
! 
! !ROUTINE: logGainEntry
! \label{logGainEntry}
!
! !INTERFACE:
  subroutine logGainEntry(pass, gain)
! 
! !USES: 

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine updates the analysis (mean, variance) computations
!  on the innovations data
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE:
! !ARGUMENTS: 

    integer  :: pass
    real     :: gain(LVT_LIS_rc(1)%ntiles, LVT_rc%nstvars)
!EOP
    integer  :: c,r, t
    integer  :: v

    if(pass.eq.1) then 
       do v=1,LVT_rc%nstvars
          do t=1,LVT_LIS_rc(1)%ntiles
             c = LVT_LIS_domain(1)%tile(t)%col
             r = LVT_LIS_domain(1)%tile(t)%row
             if(LVT_stats%datamask(c,r).eq.1.and.gain(t,v).ne.LVT_rc%udef) then 
                LVT_DAstats%mean_gain(t,v) = LVT_DAstats%mean_gain(t,v) + &
                     gain(t,v)
                LVT_DAstats%count_mean_gain(t,v)=&
                     LVT_DAstats%count_mean_gain(t,v) +1
             endif
          enddo
       enddo
    endif

  end subroutine logGainEntry

!BOP
! 
! !ROUTINE: logSpreadEntry
! \label{logSpreadEntry}
!
! !INTERFACE:
  subroutine logSpreadEntry(pass, spread)
! 
! !USES: 

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine updates the analysis (mean, variance) computations
!  on the innovations data
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE:
! !ARGUMENTS: 
    integer  :: pass
    real     :: spread(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr, LVT_rc%nstvars)
!EOP
    integer  :: c,r, t
    integer  :: v

    if(pass.eq.1) then 
       do v=1,LVT_rc%nstvars
          do t=1,LVT_LIS_rc(1)%ntiles
             c = LVT_LIS_domain(1)%tile(t)%col
             r = LVT_LIS_domain(1)%tile(t)%row
             if(LVT_stats%datamask(c,r).eq.1.and.&
                  spread(c,r,v).ne.LVT_rc%udef) then 
                LVT_DAstats%mean_spread(t,v) = &
                     LVT_DAstats%mean_spread(t,v) + &
                     spread(c,r,v)
                LVT_DAstats%count_mean_spread(t,v)=&
                     LVT_DAstats%count_mean_spread(t,v) +1

                LVT_DAstats%mean_spread_ts(t,v) = &
                     LVT_DAstats%mean_spread_ts(t,v) + &
                     spread(c,r,v)
                LVT_DAstats%count_mean_spread_ts(t,v)=&
                     LVT_DAstats%count_mean_spread_ts(t,v) +1
             endif
          enddo
       enddo
    endif
    
  end subroutine logSpreadEntry

    
#endif
!BOP
! 
! !ROUTINE: LVT_computeDAstats
! \label{LVT_computeDAstats}
!
! !INTERFACE: 
  subroutine LVT_computeDAstats(pass)
! 
! !USES:     

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in) :: pass
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine computes the statistics (mean and variance) of the 
!  innovation distribution
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer            :: c,r,v,t
    character(len=12)  :: cdate
    character(len=4)   :: cdate1
    character*100      :: fname1, fname2
    integer            :: ftn
#if 0 
    integer            :: dimID(2), nobsId, ninnov_meanId, ninnov_varId

    real               :: mean_gain(LVT_LIS_rc(1)%lnc,LVT_LIS_rc(1)%lnr,&
         LVT_rc%nstvars)
    real               :: count_mean_gain(LVT_LIS_rc(1)%lnc,LVT_LIS_rc(1)%lnr,&
         LVT_rc%nstvars)
    integer            :: shuffle, deflate, deflate_level

    shuffle = 1
    deflate = 1
    deflate_level =9

    if(LVT_rc%endtime.eq.1.and.LVT_rc%computeInnovDist.eq.1) then 
       if(pass.eq.1) then 
          do r=1,LVT_LIS_rc(1)%lnr
             do c=1,LVT_LIS_rc(1)%lnc
                if(LVT_DAstats%count_norm_innov(c,r)&
                     .gt.LVT_rc%obsCountThreshold) then 

                   LVT_DAstats%mean_norm_innov(c,r) = &
                        LVT_DAstats%sx_norm_innov(c,r)/&
                        LVT_DAstats%count_norm_innov(c,r)

                   LVT_DAstats%var_norm_innov(c,r) = & 
                        LVT_DAstats%sxx_norm_innov(c,r)/&
                        LVT_DAstats%count_norm_innov(c,r) - & 
                        LVT_DAstats%mean_norm_innov(c,r)**2

                else
                   LVT_DAstats%mean_norm_innov(c,r) = LVT_rc%udef
                   LVT_DAstats%var_norm_innov(c,r) = LVT_rc%udef
                endif
             enddo
          enddo

          do r=1,LVT_LIS_rc(1)%lnr
             do c=1,LVT_LIS_rc(1)%lnc
                if(LVT_DAstats%count_mean_innov(c,r).gt.&
                     LVT_rc%obsCountThreshold) then 
                   LVT_DAstats%mean_innov(c,r) = LVT_DAstats%mean_innov(c,r)/&
                        LVT_DAstats%count_mean_innov(c,r)
                else
                   LVT_DAstats%mean_innov(c,r) = LVT_rc%udef
                   LVT_DAstats%count_mean_innov(c,r) = LVT_rc%udef
                endif
             enddo
          enddo

          do r=1,LVT_LIS_rc(1)%lnr
             do c=1,LVT_LIS_rc(1)%lnc
                if(LVT_DAstats%count_mean_fvar(c,r).gt.&
                     LVT_rc%obsCountThreshold) then 
                   LVT_DAstats%mean_fvar(c,r) = LVT_DAstats%mean_fvar(c,r)/&
                        LVT_DAstats%count_mean_fvar(c,r)
                else
                   LVT_DAstats%mean_fvar(c,r) = LVT_rc%udef
                   LVT_DAstats%count_mean_fvar(c,r) = LVT_rc%udef
                endif
             enddo
          enddo
          call system('mkdir -p '//(LVT_rc%statsodir))
          write(unit=cdate,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2)') &
               LVT_rc%yr, LVT_rc%mo, &
               LVT_rc%da, LVT_rc%hr, LVT_rc%mn 
          
          write(unit=cdate1,fmt='(a2,i2.2)') '.d',LVT_rc%nnest
      
          fname1 = trim(LVT_rc%statsodir)//'/DA_diagnostics_'//cdate//cdate1//'.nc'
#if (defined USE_NETCDF3 || defined USE_NETCDF4)           
#if (defined USE_NETCDF4)
          call LVT_verify(nf90_create(path=fname1,cmode=nf90_hdf5,&
               ncid = ftn),&
               'creating netcdf file failed in LVT_DAMod')
#endif
#if (defined USE_NETCDF3)
          call LVT_verify(nf90_create(path=fname1,cmode=nf90_clobber,&
               ncid = ftn),&
          'creating netcdf file failed in LIS_historyMod')
#endif
          call LVT_verify(nf90_def_dim(ftn,&
               'east_west',LVT_rc%lnc,dimID(1)),&
               'nf90_def_dim failed in LVT_DAMod')
          call LVT_verify(nf90_def_dim(ftn,&
               'north_south',LVT_rc%lnr,dimID(2)),&
               'nf90_def_dim failed in LVT_DAMod')
         
          call LVT_verify(nf90_def_var(ftn,"Mean",&
               nf90_float, dimids =dimID,varID=ninnov_meanId),&
               'nf90_def_var failed for Mean in LVT_DAMod')

#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_deflate(ftn,ninnov_meanId,&
               shuffle, deflate, deflate_level), &
               'nf90_def_var_deflate failed in LVT_DAMod')
#endif

          call LVT_verify(nf90_def_var(ftn,"Variance",&
               nf90_float, dimids =dimID,varID=ninnov_varId),&
               'nf90_def_var failed for Mean in LVT_DAMod')

#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_deflate(ftn,ninnov_varId,&
               shuffle, deflate, deflate_level), &
               'nf90_def_var_deflate failed in LVT_DAMod')
#endif
          
          call LVT_verify(nf90_def_var(ftn,"Nobs",&
               nf90_float, dimids =dimID,varID=nobsId),&
               'nf90_def_var failed for Nobs in LVT_DAMod')

#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_deflate(ftn,nobsId,&
               shuffle, deflate, deflate_level), &
               'nf90_def_var_deflate failed in LVT_DAMod')
#endif
          call LVT_verify(nf90_enddef(ftn))
          
          call LVT_verify(nf90_put_var(ftn,ninnov_meanId,&
               LVT_DAstats%mean_norm_innov,(/1,1/),&
               (/LVT_rc%lnc,LVT_rc%lnr/)), &
               'nf90_put_var failed for ninnov_mean')

          call LVT_verify(nf90_put_var(ftn,ninnov_varId,&
               LVT_DAstats%var_norm_innov,(/1,1/),&
               (/LVT_rc%lnc,LVT_rc%lnr/)), &
               'nf90_put_var failed for ninnov_variance')

          call LVT_verify(nf90_put_var(ftn,nobsId,&
               real(LVT_DAstats%count_norm_innov),(/1,1/),&
               (/LVT_rc%lnc,LVT_rc%lnr/)), &
               'nf90_put_var failed for nobs')

          call LVT_verify(nf90_close(ftn))
#endif
!binary output    
#if 0 
          fname1 = trim(LVT_rc%statsodir)//'/Norm_innov_mean_'//cdate//cdate1//'.gs4r'
          ftn = LVT_getNextUnitNumber()             
          open(ftn,file=fname1,form='unformatted')
          write(ftn) LVT_DAstats%mean_norm_innov
          write(ftn) float(LVT_DAstats%count_mean_norm_innov)
          call LVT_releaseUnitNumber(ftn)

          fname1 = trim(LVT_rc%statsodir)//'/Innov_mean_'//cdate//cdate1//'.gs4r'
          ftn = LVT_getNextUnitNumber()             
          open(ftn,file=fname1,form='unformatted')
          write(ftn) LVT_DAstats%mean_innov
          write(ftn) float(LVT_DAstats%count_mean_innov)
          call LVT_releaseUnitNumber(ftn)

          fname1 = trim(LVT_rc%statsodir)//'/Forecast_var_'//cdate//cdate1//'.gs4r'
          ftn = LVT_getNextUnitNumber()             
          open(ftn,file=fname1,form='unformatted')
          write(ftn) LVT_DAstats%mean_fvar
          write(ftn) float(LVT_DAstats%count_mean_fvar)
          call LVT_releaseUnitNumber(ftn)

          fname2 = trim(LVT_rc%statsodir)//'/Norm_innov_var_'//cdate//cdate1//'.gs4r'
          ftn = LVT_getNextUnitNumber()             
          open(ftn,file=fname2,form='unformatted')
          write(ftn) LVT_DAstats%var_norm_innov
          write(ftn) float(LVT_DAstats%count_var_norm_innov)
          call LVT_releaseUnitNumber(ftn)

#endif
       endif
    endif
    
    if(LVT_rc%endtime.eq.1.and.LVT_rc%computeGain.eq.1) then 

       if(pass.eq.1) then 
          do v=1, LVT_rc%nstvars
             do t=1,LVT_LIS_rc(1)%ntiles
                if(LVT_DAstats%count_mean_gain(t,v)&
                     .gt.LVT_rc%obsCountThreshold) then 
                   LVT_DAstats%mean_gain(t,v) = LVT_DAstats%mean_gain(t,v)/&
                        LVT_DAstats%count_mean_gain(t,v)
                else
                   LVT_DAstats%mean_gain(t,v) = LVT_rc%udef
                   LVT_DAstats%count_mean_gain(t,v) = LVT_rc%udef
                endif
             enddo
          enddo
          
          mean_gain = 0.0
          count_mean_gain = 0 
          do v=1, LVT_rc%nstvars
             do t=1,LVT_LIS_rc(1)%ntiles             
                c = LVT_LIS_domain(1)%tile(t)%col
                r = LVT_LIS_domain(1)%tile(t)%row             
                if(LVT_DAstats%mean_gain(t,v).ne.LVT_rc%udef) then
                   mean_gain(c,r,v) = mean_gain(c,r,v) + & 
                        LVT_DAstats%mean_gain(t,v)
                   count_mean_gain(c,r,v) = count_mean_gain(c,r,v) + 1
                endif
             enddo

             do r=1,LVT_LIS_rc(1)%lnr
                do c=1,LVT_LIS_rc(1)%lnc
                   if(count_mean_gain(c,r,v).ne.0.0) then 
                      mean_gain(c,r,v) = mean_gain(c,r,v)/count_mean_gain(c,r,v)
                   else
                      mean_gain(c,r,v) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo

          call system('mkdir -p '//(LVT_rc%statsodir))
          write(unit=cdate,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2)') &
               LVT_rc%yr, LVT_rc%mo, &
               LVT_rc%da, LVT_rc%hr, LVT_rc%mn 
          
          write(unit=cdate1,fmt='(a2,i2.2)') '.d',LVT_rc%nnest
          
          fname1 = trim(LVT_rc%statsodir)//'/Gain_mean_'//cdate//cdate1//'.gs4r'
          ftn = LVT_getNextUnitNumber()
          do v=1,LVT_rc%nstvars
             open(ftn,file=fname1,form='unformatted')
             write(ftn) mean_gain(:,:,v)
          enddo
          call LVT_releaseUnitNumber(ftn)
       endif
    endif
#endif
  end subroutine LVT_computeDAstats

!BOP
! 
! !ROUTINE: create_dainnov_filename
! \label{create_dainnov_filename}
!
! !INTERFACE: 
  subroutine create_dainnov_filename(filen)
! 
! !USES:     

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine generates the innovation filename
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    character(len=*)  :: filen
!EOP
    integer          :: no

    character(len=12) :: cdate1
    character(len=6)  :: cdate2
    
    write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
         LVT_rc%yr, LVT_rc%mo, LVT_rc%da, LVT_rc%hr, LVT_rc%mn

    write(unit=cdate2, fmt='(i4.4, i2.2)') &
         LVT_rc%yr, LVT_rc%mo

    filen = trim(LVT_rc%odir)//'/EnKF/'// &
         trim(cdate2)//'/LIS_DA_EnKF_'//trim(cdate1)//'_innov.d01.nc'

  end subroutine create_dainnov_filename

!BOP
! 
! !ROUTINE: create_dagain_filename
! \label{create_dagain_filename}
!
! !INTERFACE: 
  subroutine create_dagain_filename(filen)
! 
! !USES:     

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine generates the gaination filename
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    character(len=*)  :: filen
!EOP
    integer          :: no

    character(len=12) :: cdate1
    
    write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
         LVT_rc%yr, LVT_rc%mo, LVT_rc%da, LVT_rc%hr, LVT_rc%mn

    filen = trim(LVT_rc%odir)//'/GMAOEnKF/'&
         //'GMAOEnKF.d01.'//trim(cdate1)//'.gain'

  end subroutine create_dagain_filename

!BOP
! 
! !ROUTINE: create_daspread_filename
! \label{create_daspread_filename}
!
! !INTERFACE: 
  subroutine create_daspread_filename(filen)
! 
! !USES:     

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine generates the spreadation filename
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    character(len=*)  :: filen
!EOP
    integer          :: no

    character(len=12) :: cdate1
    character(len=6)  :: cdate

    write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
         LVT_rc%yr, LVT_rc%mo, LVT_rc%da, LVT_rc%hr, LVT_rc%mn

    write(unit=cdate, fmt='(i4.4, i2.2)') &
         LVT_rc%yr, LVT_rc%mo

    filen = trim(LVT_rc%odir)//'/EnKF/'//trim(cdate)//'/'&
         //'LIS_DA_EnKF_'//trim(cdate1)//'_spread_a01.d01.nc'

  end subroutine create_daspread_filename

end module LVT_DAMod

