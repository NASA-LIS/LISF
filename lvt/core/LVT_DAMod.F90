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

