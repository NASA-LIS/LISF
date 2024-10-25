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
module LDT_ANNMod
!BOP
!
!  !MODULE: LDT_ANNMod
! 
!  !DESCRIPTION: 
!   This implementation is based on the original code of 
!   Philip Brierly. The reference is: 
!   Bierly, P. and Batty, B., "Data mining with neural networks - 
!   an applied example in understanding electricity consumption
!   patterns", Chapter 12 in 'Knowledge discovery and data mining'
!   edited by M.A. Bramer, The Institution of Electrical Engineers, 
!   1999.
!   
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! !USES: 

  use ESMF
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS: 
!-----------------------------------------------------------------------------
  public :: LDT_ANNinit
  public :: LDT_readANNinputData
  public :: LDT_readANNoutputData
  public :: LDT_tavgANNdata
  public :: LDT_logSingleANNdata
  public :: LDT_outputANN
  public :: LDT_resetANN

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LDT_ANNinput
  public :: LDT_ANNoutput
  public :: LDT_ANNdata


  type, public :: LDT_ANNdataEntry
     
     integer               :: nTimes
     real, allocatable     :: value(:,:) 
     integer, allocatable  :: count(:)
     real,    allocatable  :: value_hold(:)
     
  end type LDT_ANNdataEntry
  
  type, public :: LDT_ANNdataStruc
     
     character*100,          allocatable :: varName(:)
     character*50,           allocatable :: varUnits(:)
     type(LDT_ANNdataEntry), allocatable :: dataEntry(:)

  end type LDT_ANNdataStruc
  
  type, public :: ANNdatadec
     character*50           :: mode
     integer                :: trainMode
     integer                :: nINPUTS
     integer                :: nOUTPUTS
     integer                :: nTOTAL
     integer                :: nppb
     integer                    :: nsources
     character*50, allocatable  :: inputsrc(:)
     integer,      allocatable  :: input_nparams(:)
     integer,      allocatable  :: input_param_s(:)
     integer,      allocatable  :: input_param_e(:)

     character*50           :: outputsrc
     integer                :: ntimes
     integer                :: hidden_neurons
     integer                :: nhs
     integer                :: nhf
     integer                :: nos
     real                   :: alr !learning rate
     real                   :: blr!learning rate
     integer                :: nIters
!     real, allocatable      :: trainingInputs(:,:,:)
!     real, allocatable      :: trainingOutputs(:,:)
!     real, allocatable      :: inputs_this_pattern(:,:)

     real                   :: rmse_best
     real, allocatable      :: wi(:,:,:)
     real, allocatable      :: wi_best(:,:,:)
     real, allocatable      :: wo(:,:)
     real, allocatable      :: wo_best(:,:)
     real, allocatable      :: hval(:)
     real, allocatable      :: dummy1(:,:,:)
     real, allocatable      :: dummy2(:,:,:)

     real, allocatable    :: maxinp(:,:)
     real, allocatable    :: mininp(:,:)
     real, allocatable    :: maxout(:,:)
     real, allocatable    :: minout(:,:)

     character(len=LDT_CONST_PATH_LEN) :: outfile
  end type ANNdatadec

  type(LDT_ANNdataStruc) :: LDT_ANNinput
  type(LDT_ANNdataStruc) :: LDT_ANNoutput

  type(ANNdatadec)          :: LDT_ANNdata

!EOP
  
contains
  subroutine LDT_ANNinit()
    
    integer :: rc
    integer                 :: nInputs
    integer                 :: nsize
    integer                 :: ntimes
    type(ESMF_Time)         :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: n,i,t
    
    n = 1
    call ESMF_ConfigGetAttribute(LDT_config,LDT_ANNdata%mode,&
         label="ANN mode (training/prediction):" ,rc=rc)
    call LDT_verify(rc,&
         'ANN mode (training/prediction): option not specified in the config file')

    call ESMF_ConfigGetAttribute(LDT_config,LDT_ANNdata%trainMode,&
         label="ANN conduct training in a spatially distributed manner:" ,rc=rc)
    call LDT_verify(rc,&
         'ANN conduct training in a spatially distributed manner: option not specified in the config file')
    
    call ESMF_ConfigGetAttribute(LDT_config,LDT_ANNdata%nsources, &
         label="ANN number of input data sources:", rc=rc)
    call LDT_verify(rc,&
         'ANN number of input data sources: not specified')

    allocate(LDT_ANNdata%inputsrc(LDT_ANNdata%nsources))
    allocate(LDT_ANNdata%input_nparams(LDT_ANNdata%nsources))
    allocate(LDT_ANNdata%input_param_s(LDT_ANNdata%nsources))
    allocate(LDT_ANNdata%input_param_e(LDT_ANNdata%nsources))
    
    call ESMF_ConfigFindLabel(LDT_config,'ANN input data sources:',rc=rc)
    do i=1,LDT_ANNdata%nsources
       call ESMF_ConfigGetAttribute(LDT_config, LDT_ANNdata%inputsrc(i),rc=rc)
       call LDT_verify(rc,&
         'ANN input data sources: option not specified in the config file')
    enddo

    call ESMF_ConfigFindLabel(LDT_config,&
         'ANN number of parameters in each input source:',rc=rc)
    do i=1,LDT_ANNdata%nsources
       call ESMF_ConfigGetAttribute(LDT_config, &
            LDT_ANNdata%input_nparams(i),rc=rc)
       call LDT_verify(rc,&
         'ANN number of parameters in each input source: option not specified in the config file')
    enddo

    LDT_ANNdata%input_param_s(1) = 1
    LDT_ANNdata%input_param_e(1) = LDT_ANNdata%input_param_s(1) + & 
         LDT_ANNdata%input_nparams(1) - 1

    nInputs  = LDT_ANNdata%input_nparams(1)

    do i=2,LDT_ANNdata%nsources
       LDT_ANNdata%input_param_s(i) = LDT_ANNdata%input_param_e(i-1) + 1
       LDT_ANNdata%input_param_e(i) = LDT_ANNdata%input_param_s(i) + &
            LDT_ANNdata%input_nparams(i)-1
       
       nInputs = nInputs + LDT_ANNdata%input_nparams(i)
    enddo

! Additional entry to store seasonality
!    nInputs = nInputs + 1


!   HUANG, G.-B., 2003, Learning capability and storage 
!   capacity of two-hidden-layer feedforward networks. 
!   IEEE Transactions on Neural Networks, 14, pp. 274–281.
! 
!   Huang (2003) proved that in the two-hidden-layer case, 
!   with m output neurons, the number of hidden nodes that 
!   are enough to learn N samples with negligibly small error 
!   is given by
!     2*sqrt((m+2)*N)
! 
!   Specifically, he suggests that the sufficient number 
!   of hidden nodes in the first layer is
!    sqrt((m+2)*N) + 2sqrt(N/(m+2))
! 
!   and in the second it is
!    m * sqrt(N/(m+2))
!
!
!    call ESMF_ConfigGetAttribute(LDT_config,LDT_ANNdata%hidden_neurons,&
!         label="ANN number of hidden neurons:",rc=rc)
!    call LDT_verify(rc,&
!         'ANN number of hidden neurons: option not specified in the config file')

    LDT_ANNdata%hidden_neurons = nint(2.0*sqrt((1+2)*real(nInputs)))

    call ESMF_ConfigGetAttribute(LDT_config,LDT_ANNdata%nIters,&
         label="ANN number of iterations:",rc=rc)
    call LDT_verify(rc,&
         'ANN number of iterations: option not specified in the config file')
    
    call ESMF_ConfigGetAttribute(LDT_config,LDT_ANNdata%outfile,&
         label="ANN training output file:",rc=rc)
    call LDT_verify(rc,&
         'ANN training output file: option not specified in the config file')

    if(LDT_ANNdata%mode.eq."training") then 
       call ESMF_ConfigGetAttribute(LDT_config,LDT_ANNdata%outputsrc,&
            label="ANN output data source:" ,rc=rc)
       call LDT_verify(rc,&
            'ANN output data source: option not specified in the config file')
    endif

    call ESMF_ClockGet(LDT_clock, startTime = startTime, &
         stopTime = stopTime)

    call ESMF_TimeIntervalSet(timeStep, s = nint(LDT_rc%tavgInterval),&
         rc=rc)

    ntimes = ((stopTime - startTime)/timestep)+1

    LDT_ANNdata%nInputs = nInputs
    LDT_ANNdata%ntimes = ntimes
    LDT_ANNdata%nOutputs = 1
    LDT_ANNdata%nTotal = LDT_ANNdata%nInputs + LDT_ANNdata%nOutputs
    LDT_ANNdata%nppb = LDT_ANNdata%nInputs + 1

    LDT_ANNdata%hidden_neurons = &
         LDT_ANNdata%hidden_neurons + 1 !acccounts for bias to output

!    LDT_ANNdata%nhs = LDT_ANNdata%nppb + 1
!    LDT_ANNdata%nhf = LDT_ANNdata%nppb + LDT_ANNdata%hidden_neurons
    LDT_ANNdata%nhs = 1
    LDT_ANNdata%nhf = LDT_ANNdata%hidden_neurons
    LDT_ANNdata%nos = LDT_ANNdata%nhf + 1
    
    LDT_ANNdata%alr = 0.1
    LDT_ANNdata%blr = LDT_ANNdata%alr/10.0

    allocate(LDT_ANNinput%dataEntry(LDT_rc%ngrid(n)))
    allocate(LDT_ANNoutput%dataEntry(LDT_rc%ngrid(n)))

    allocate(LDT_ANNinput%varName(&
         LDT_ANNdata%nppb))
    allocate(LDT_ANNinput%varUnits(&
         LDT_ANNdata%nppb))
    LDT_ANNinput%varName(LDT_ANNdata%nppb) = "Bias parameter"
    LDT_ANNinput%varUnits(LDT_ANNdata%nppb) = "-"

    if(LDT_ANNdata%mode.eq."training") then 
       do t=1,LDT_rc%ngrid(n)
          allocate(LDT_ANNinput%dataEntry(t)%value(&
               LDT_ANNdata%nppb,&
               LDT_ANNdata%ntimes))
          allocate(LDT_ANNinput%dataEntry(t)%value_hold(&
               LDT_ANNdata%nInputs))
          allocate(LDT_ANNinput%dataEntry(t)%count(&
               LDT_ANNdata%nInputs))
          
          LDT_ANNinput%dataEntry(t)%value = LDT_rc%udef
          LDT_ANNinput%dataEntry(t)%value(LDT_ANNdata%nppb,:) = 1
          LDT_ANNinput%dataEntry(t)%value_hold = 0 
          LDT_ANNinput%dataEntry(t)%count = 0 
          
       enddo
    elseif(LDT_ANNdata%mode.eq."prediction") then 
       do t=1,LDT_rc%ngrid(n)
          allocate(LDT_ANNinput%dataEntry(t)%value(&
               LDT_ANNdata%nppb,1))
          allocate(LDT_ANNinput%dataEntry(t)%value_hold(&
               LDT_ANNdata%nInputs))
          allocate(LDT_ANNinput%dataEntry(t)%count(&
               LDT_ANNdata%nInputs))
          
          LDT_ANNinput%dataEntry(t)%value = LDT_rc%udef
          LDT_ANNinput%dataEntry(t)%value(LDT_ANNdata%nppb,:) = 1
          LDT_ANNinput%dataEntry(t)%value_hold = 0 
          LDT_ANNinput%dataEntry(t)%count = 0 
          
       enddo
       
    endif
    allocate(LDT_ANNoutput%varName(1))
    allocate(LDT_ANNoutput%varUnits(1))

    if(LDT_ANNdata%mode.eq."training") then 
       do t=1,LDT_rc%ngrid(n)
          allocate(LDT_ANNoutput%dataEntry(t)%value(&
               1,LDT_ANNdata%ntimes))
          allocate(LDT_ANNoutput%dataEntry(t)%value_hold(&
               1))
          allocate(LDT_ANNoutput%dataEntry(t)%count(&
               1))
          
          LDT_ANNoutput%dataEntry(t)%value = LDT_rc%udef
          LDT_ANNoutput%dataEntry(t)%value_hold = 0 
          LDT_ANNoutput%dataEntry(t)%count = 0 
          
       enddo

    elseif(LDT_ANNdata%mode.eq."prediction") then 
       do t=1,LDT_rc%ngrid(n)
          allocate(LDT_ANNoutput%dataEntry(t)%value(&
               1,1))
          allocate(LDT_ANNoutput%dataEntry(t)%value_hold(&
               1))
          allocate(LDT_ANNoutput%dataEntry(t)%count(&
               1))
          
          LDT_ANNoutput%dataEntry(t)%value = LDT_rc%udef
          LDT_ANNoutput%dataEntry(t)%value_hold = 0 
          LDT_ANNoutput%dataEntry(t)%count = 0 
          
       enddo

    endif

    if(LDT_ANNdata%trainMode.eq.1) then 
       allocate(LDT_ANNdata%maxinp(LDT_rc%ngrid(n),LDT_ANNdata%nppb))
       allocate(LDT_ANNdata%mininp(LDT_rc%ngrid(n),LDT_ANNdata%nppb))
       allocate(LDT_ANNdata%maxout(LDT_rc%ngrid(n),1))
       allocate(LDT_ANNdata%minout(LDT_rc%ngrid(n),1))       
    else
       allocate(LDT_ANNdata%maxinp(1,LDT_ANNdata%nppb))
       allocate(LDT_ANNdata%mininp(1,LDT_ANNdata%nppb))
       allocate(LDT_ANNdata%maxout(1,1))
       allocate(LDT_ANNdata%minout(1,1))
    endif
    
    LDT_ANNdata%maxinp = -1E20
    LDT_ANNdata%mininp = 1E20
    LDT_ANNdata%maxout = -1E20
    LDT_ANNdata%minout = 1E20

    do i=1,LDT_ANNdata%nsources
       call setupANNinputsource(trim(LDT_ANNdata%inputsrc(i))//char(0))
    enddo

    if(LDT_ANNdata%mode.eq."training") then 
       call setupANNoutputsource(trim(LDT_ANNdata%outputsrc)//char(0))
    endif

    call ANN_initiateWeights(n)

    if(LDT_ANNdata%mode.eq."prediction") then 
       write(LDT_logunit,*) '[INFO] Reading ANN output..'
       call readANNtrainingOutput(n)
       write(LDT_logunit,*) '[INFO] Done reading ANN output..'
    endif

  end subroutine LDT_ANNinit

!BOP
! !ROUTINE: LDT_readANNinputData
!  \label{LDT_readANNinputData}
! 
! !INTERFACE: 
  subroutine LDT_readANNinputData(n)
! !ARGUMENTS: 
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This subroutine reads the input data sources used in the ANN 
!  training or prediction
!EOP

    integer             :: iomode
    integer             :: i

    iomode = 1
    do i=1,LDT_ANNdata%nsources
       call readANNinputsource(trim(LDT_ANNdata%inputsrc(i))//char(0),&
            n,iomode,LDT_ANNdata%input_param_s(i),LDT_ANNdata%input_param_e(i))
    enddo

  end subroutine LDT_readANNinputData

!BOP
! !ROUTINE: LDT_readANNoutputData
!  \label{LDT_readANNoutputData}
! 
! !INTERFACE: 
  subroutine LDT_readANNoutputData(n)
! !ARGUMENTS: 
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This subroutine reads the output data sources for the ANN training
!EOP
    integer             :: iomode

    if(LDT_ANNdata%mode.eq."training") then 
       iomode = 2
       call readANNoutputsource(trim(LDT_ANNdata%outputsrc)//char(0),&
            n,iomode,1,1)
    endif

  end subroutine LDT_readANNoutputData

!BOP
! !ROUTINE: LDT_tavgANNdata
! \label{LDT_tavgANNdata}
! 
! !INTERFACE: 
  subroutine LDT_tavgANNdata(n)
! !ARGUMENTS: 
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This subroutine temporally averages the datasets for use in ANN
!  training or prediction. 
!EOP

    type(ESMF_Time)         :: currTime, startTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: tindex
    integer                 :: c,r,pindex, gid
    real                    :: time_val
    integer                 :: rc, status

    if(mod(real(LDT_rc%hr)*3600+60*real(LDT_rc%mn)+float(LDT_rc%ss),&
            real(LDT_rc%tavgInterval)).eq.0) then        
       
       if(LDT_ANNdata%mode.eq."training") then 
          call ESMF_TimeSet(currTime, yy=LDT_rc%yr, &
               mm=LDT_rc%mo,dd=LDT_rc%da,h=LDT_rc%hr,&
               m=LDT_rc%mn,s=LDT_rc%ss,calendar=LDT_calendar,&
               rc=status)
          call LDT_verify(status, 'error in ESMF_TimeSet: in LDT_ANNMod')
          
          call ESMF_ClockGet(LDT_clock, startTime = startTime)
          
          call ESMF_TimeIntervalSet(timeStep, s = nint(LDT_rc%tavgInterval),&
               rc=rc)
          
          tindex = ((currTime - startTime)/timestep)+1
       else
          tindex = 1
       endif

       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(LDT_domain(n)%gindex(c,r).ne.-1) then 
                gid = LDT_domain(n)%gindex(c,r)
!ninputs -1 if seasonality is included
!                do pindex=1,LDT_ANNdata%nInputs-1
                do pindex=1,LDT_ANNdata%nInputs
                   if(LDT_ANNinput%dataEntry(gid)%count(pindex).gt.0) then 
                      LDT_ANNinput%dataEntry(gid)%value_hold(pindex) = & 
                           LDT_ANNinput%dataEntry(gid)%value_hold(pindex) /& 
                           LDT_ANNinput%dataEntry(gid)%count(pindex) 
                      LDT_ANNinput%dataEntry(gid)%value(pindex,tindex) = &
                           LDT_ANNinput%dataEntry(gid)%value_hold(pindex)
                      
                      LDT_ANNinput%dataEntry(gid)%value_hold(pindex) = 0
                      LDT_ANNinput%dataEntry(gid)%count(pindex) = 0 
                   endif
                enddo
!if seasonality is to be included explicitly as part of the inputs
!
!                time_val = LDT_rc%ss + & 
!                     LDT_rc%mn*60.0 +& 
!                     LDT_rc%hr*60*60.0 +& 
!                     (LDT_rc%da-1)*24*60*60.0 +& 
!                     (LDT_rc%mo-1)*30*24*60*60.0
!                     
!                LDT_ANNinput%dataEntry(gid)%value(&
!                     LDT_ANNdata%nInputs,tindex) = & 
!                     time_val
             endif
             
          enddo
       enddo

       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(LDT_domain(n)%gindex(c,r).ne.-1) then 
                gid = LDT_domain(n)%gindex(c,r)
                if(LDT_ANNoutput%dataEntry(gid)%count(1).gt.0) then 
                   LDT_ANNoutput%dataEntry(gid)%value_hold(1) = & 
                        LDT_ANNoutput%dataEntry(gid)%value_hold(1) / & 
                        LDT_ANNoutput%dataEntry(gid)%count(1)
                   LDT_ANNoutput%dataEntry(gid)%value(1,tindex) = &
                        LDT_ANNoutput%dataEntry(gid)%value_hold(1) 
                   
                   if(gid.eq.2) then 
                      print*, LDT_rc%yr, LDT_rc%mo, LDT_rc%da, LDT_ANNoutput%dataEntry(gid)%value(1,tindex)
                   endif

                   LDT_ANNoutput%dataEntry(gid)%value_hold(1) = 0
                   LDT_ANNoutput%dataEntry(gid)%count(1) = 0 
                endif
             endif
          enddo
       enddo
    endif
  end subroutine LDT_tavgANNdata

  subroutine LDT_outputANN(n)

    integer, intent(in) :: n 

    if(LDT_ANNdata%mode.eq."training") then 
       call trainANN(n)
    elseif(LDT_ANNdata%mode.eq."prediction") then 
       call predictANN(n)
    endif
  end subroutine LDT_outputANN

  subroutine trainANN(n)

    integer, intent(in) :: n 
    
    integer             :: i,j,t
    integer             :: iPAT_NUM
    real                :: rRand
    real                :: output_this_pat
    real                :: er_this_pat
    real                :: rmse
    real                :: dummy2(LDT_ANNdata%nppb,1)
    real                :: dummy1(1,LDT_ANNdata%nhf)
    real                :: inputs_this_pattern(LDT_ANNdata%nppb)


    if(LDT_rc%endtime.eq.1) then 

       call random_seed  
       call ANN_scaleData(n)       
       
       write(LDT_logunit,*) '[INFO] Completed reading datasets...'

       if(LDT_ANNdata%trainMode.eq.1) then !spatially distributed training
          do t = 1, LDT_rc%ngrid(n)
             do j=1,LDT_ANNdata%nIters
                do i=1,LDT_ANNinput%dataEntry(t)%nTimes

                   call random_number(rRand)
                   
                   !select a pattern at random
                   iPAT_NUM= nint(rRand*(LDT_ANNinput%dataEntry(t)%ntimes-1))+1
                   !select the data to this pattern
                   !                   iPAT_NUM = i
                   inputs_this_pattern = & 
                        LDT_ANNinput%dataEntry(t)%value(:,iPAT_NUM)
                   output_this_pat = LDT_ANNoutput%dataEntry(t)%value(1,iPAT_NUM)
                   er_this_pat = calc_err_pat(t,&
                        inputs_this_pattern,&
                        output_this_pat, & 
                        LDT_ANNdata%wi(t,:,:), & 
                        LDT_ANNdata%wo(t,:))
                   
                   !                print*, 'ip',inputs_this_pattern
                   !                print*, 'ou ',output_this_pat
                   !                print*, 'wi ',LDT_ANNdata%wi(t,:,:)
                   !                print*, 'wo ',LDT_ANNdata%wo(t,:)
                   !                print*, 'er ',er_this_pat
                   
                   LDT_ANNdata%wo(t,:) = LDT_ANNdata%wo(t,:) - &
                        LDT_ANNdata%blr*&
                        LDT_ANNdata%hval(:)*&
                        er_this_pat
                   dummy2(:,1) = LDT_ANNinput%dataEntry(t)%value(:,iPAT_NUM)
                   dummy1(1,:) = er_this_pat*LDT_ANNdata%wo(t,:)*&
                        (1-(LDT_ANNdata%hval(:)**2))
                   LDT_ANNdata%wi(t,:,:)  = & 
                        LDT_ANNdata%wi(t,:,:) - & 
                        (matmul(dummy2,dummy1)*LDT_ANNdata%alr)
                enddo
                rmse = calculate_error(t)
                
                call keep_best_weights(t,j,rmse)
             enddo
          enddo
          call writeANNoutput(n)
       else
          do j=1,LDT_ANNdata%nIters
             do t = 1, LDT_rc%ngrid(n)
                do i=1,LDT_ANNinput%dataEntry(t)%nTimes
                   call random_number(rRand)
                   
                   !select a pattern at random
                   iPAT_NUM= nint(rRand*(LDT_ANNinput%dataEntry(t)%ntimes-1))+1
                   !select the data to this pattern
                   !                   iPAT_NUM = i
                   inputs_this_pattern = & 
                        LDT_ANNinput%dataEntry(t)%value(:,iPAT_NUM)
                   output_this_pat = LDT_ANNoutput%dataEntry(t)%value(1,iPAT_NUM)
                   
                   er_this_pat = calc_err_pat(1,&
                        inputs_this_pattern,&
                        output_this_pat, & 
                        LDT_ANNdata%wi(1,:,:), & 
                        LDT_ANNdata%wo(1,:))
                   
                   !                print*, 'ip',inputs_this_pattern
                   !                print*, 'ou ',output_this_pat
                   !                print*, 'wi ',LDT_ANNdata%wi(t,:,:)
                   !                print*, 'wo ',LDT_ANNdata%wo(t,:)
                   !                print*, 'er ',er_this_pat
                   
                   LDT_ANNdata%wo(1,:) = LDT_ANNdata%wo(1,:) - &
                        LDT_ANNdata%blr*&
                        LDT_ANNdata%hval(:)*&
                        er_this_pat
                   dummy2(:,1) = LDT_ANNinput%dataEntry(t)%value(:,iPAT_NUM)
                   dummy1(1,:) = er_this_pat*LDT_ANNdata%wo(1,:)*&
                        (1-(LDT_ANNdata%hval(:)**2))
                   LDT_ANNdata%wi(1,:,:)  = & 
                        LDT_ANNdata%wi(1,:,:) - & 
                        (matmul(dummy2,dummy1)*LDT_ANNdata%alr)
                enddo
             enddo
             rmse = calculate_error_nondist(n)
             
             call keep_best_weights_nondist(j,rmse)
          enddo
          call writeANNoutput_nondist(n)
                 
       endif

    end if
  end subroutine trainANN

!BOP
! !ROUTINE: predictANN
! \label{predictANN}
! 
! !INTERFACE: 
  subroutine predictANN(n)
! !USES:
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer, intent(in) :: n 
!
! !DESCRIPTION: 
!  This routine employs the training output from the ANN to 
!  generate predictions 
!
!EOP    
    integer             :: i,j,t
    integer             :: ios
    real                :: output_this_pat
    real                :: output_value(LDT_rc%ngrid(n))
    real                :: inputs_this_pattern(LDT_ANNdata%nppb)

    character*6         :: cdate1
    character*10        :: cdate
    integer             :: ftn
#if (!defined AIX )
   integer            :: system 
#endif
   integer, external  :: LDT_create_subdirs 
   character(len=201) :: c_string     
   character(len=LDT_CONST_PATH_LEN)      :: filename
   integer            :: c,r,gid,dimID(2),tdimID
   integer            :: latId, lonId, varId,xtimeID
   real               :: lat(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real               :: lon(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real               :: value_2d(LDT_rc%lnc(n),LDT_rc%lnr(n))
   logical            :: check 
   integer            :: count,g,p

   call random_seed  

   do g=1,LDT_rc%ngrid(n)
      count = 0
      
      check = .true. 
      do p=1,LDT_ANNdata%nppb
         check = check .and. & 
              (LDT_ANNinput%dataEntry(g)%value(p,1).ne.LDT_rc%udef)
      enddo
      if(check) then 
         count = count + 1
      endif
      
      if(count.gt.0) then !all valid input data
         LDT_ANNinput%dataEntry(g)%nTimes = count
      else
         LDT_ANNinput%dataEntry(g)%nTimes = 0
      endif
   enddo

    ! scale the data
    if(LDT_ANNdata%trainMode.eq.1) then 
       do g=1,LDT_rc%ngrid(n)
          do t=1,LDT_ANNinput%dataEntry(g)%nTimes
             do p=1,LDT_ANNdata%ninputs
                if((LDT_ANNdata%maxinp(g,p)-LDT_ANNdata%mininp(g,p)).ne.0) then 
                   LDT_ANNinput%dataEntry(g)%value(p,t) = & 
                        ((LDT_ANNinput%dataEntry(g)%value(p,t) - &
                        LDT_ANNdata%mininp(g,p))/&
                        (LDT_ANNdata%maxinp(g,p)-LDT_ANNdata%mininp(g,p)) - 0.5)*2
                endif
             enddo
          enddo
       enddo
    else
       do g=1,LDT_rc%ngrid(n)
          do t=1,LDT_ANNinput%dataEntry(g)%nTimes
             do p=1,LDT_ANNdata%ninputs
                if((LDT_ANNdata%maxinp(1,p)-LDT_ANNdata%mininp(1,p)).ne.0) then 
                   LDT_ANNinput%dataEntry(g)%value(p,t) = & 
                        ((LDT_ANNinput%dataEntry(g)%value(p,t) - &
                        LDT_ANNdata%mininp(1,p))/&
                        (LDT_ANNdata%maxinp(1,p)-LDT_ANNdata%mininp(1,p)) - 0.5)*2
                endif
             enddo
          enddo
       enddo
    endif


   do t=1,LDT_rc%ngrid(n)
      if(LDT_ANNinput%dataEntry(t)%nTimes.gt.0) then 
         if(LDT_ANNdata%trainMode.eq.1) then 
            call calc_out_pat(t, &
                 LDT_ANNinput%dataEntry(t)%value(:,1), &
                 output_this_pat, &
                 LDT_ANNdata%wi(t,:,:),&
               LDT_ANNdata%wo(t,:))
            output_value(t) = & 
                 (output_this_pat/2.0 + 0.5) * & 
                 (LDT_ANNdata%maxout(t,1)-LDT_ANNdata%minout(t,1)) + & 
                 LDT_ANNdata%minout(t,1)  
            
            if(t.eq.2) then 
               print*, LDT_rc%yr, LDT_rc%mo, LDT_rc%da,output_value(t)
            endif
         else
            call calc_out_pat(t, &
                 LDT_ANNinput%dataEntry(t)%value(:,1), &
                 output_this_pat, &
                 LDT_ANNdata%wi(1,:,:),&
                 LDT_ANNdata%wo(1,:))
            
            output_value(t) = & 
                 (output_this_pat/2.0 + 0.5) * & 
                 (LDT_ANNdata%maxout(1,1)-LDT_ANNdata%minout(1,1)) + & 
                 LDT_ANNdata%minout(1,1)  
            
         endif
      else
         output_value(t) = LDT_rc%udef
      endif
   enddo

    !write the data out
    
    write(unit=cdate,fmt='(i4.4,i2.2,i2.2,i2.2)') &
         LDT_rc%yr, LDT_rc%mo, LDT_rc%da, LDT_rc%hr
    
    write(unit=cdate1,fmt='(i4.4,i2.2)') &
         LDT_rc%yr, LDT_rc%mo
    
    filename = trim(LDT_rc%odir)//'/'//trim(cdate1)
    
#if ( defined AIX )
    call system('mkdir -p '//trim(filename),ios)
#else
    c_string = trim(filename)
    ios = LDT_create_subdirs(len_trim(c_string),trim(c_string))
#endif       
    filename = trim(filename)//'/ANN_fcst_'//trim(cdate)//'.nc'
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    
    value_2d = LDT_rc%udef

    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)

          lat(c,r) = LDT_domain(n)%lat(c+(r-1)*LDT_rc%lnc(n))
          lon(c,r) = LDT_domain(n)%lon(c+(r-1)*LDT_rc%lnc(n))

          gid = LDT_domain(n)%gindex(c,r)
          if(gid.ne.-1) then 
             value_2d(c,r) = output_value(gid)
          endif
       enddo
    enddo
    ios = nf90_create(path=trim(filename),cmode=nf90_hdf5,&
         ncid = ftn)
    call LDT_verify(ios, "nf90_created failed for "//trim(filename))
    ios = nf90_def_dim(ftn,"east_west",LDT_rc%lnc(n),dimID(1))
    call LDT_verify(ios, "nf90_def_dim failed in LDT_ANNMod")
    ios = nf90_def_dim(ftn,"north_south",LDT_rc%lnr(n),dimID(2))
    call LDT_verify(ios, "nf90_def_dim failed in LDT_ANNMod")
    call LDT_verify(nf90_def_dim(ftn,'time',1,tdimID))

    ios = nf90_def_var(ftn,"lat",nf90_float,dimIDs=dimID,varID=latId)
    call LDT_verify(ios, "nf90_def_var failed for lat")
    ios = nf90_def_var(ftn,"lon",nf90_float,dimIDs=dimID,varID=lonId)
    call LDT_verify(ios, "nf90_def_var failed for lon")

    call LDT_verify(nf90_def_var(ftn,'time',&
         nf90_float,dimids = tdimID,varID=xtimeID))

    ios = nf90_def_var(ftn, trim(LDT_ANNoutput%varName(1)),&
         nf90_float,dimIDs=dimID,varID=varID)
    call LDT_verify(ios, "nf90_def_var failed for "//trim(LDT_ANNoutput%varName(1)))
    ios = nf90_put_att(ftn,NF90_GLOBAL,"missing_value",-9999.0)
    
    if(trim(LDT_rc%lis_map_proj(n)).eq."latlon") then !latlon
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
            "EQUIDISTANT CYLINDRICAL"))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LAT", &
            LDT_rc%gridDesc(n,4)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LON", &
            LDT_rc%gridDesc(n,5)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
            LDT_rc%gridDesc(n,9)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
            LDT_rc%gridDesc(n,10)))       
       
    elseif(trim(LDT_rc%lis_map_proj(n)).eq."mercator") then 
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
            "MERCATOR"))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LAT", &
            LDT_rc%gridDesc(n,4)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LON", &
            LDT_rc%gridDesc(n,5)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
            LDT_rc%gridDesc(n,10)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
            LDT_rc%gridDesc(n,11)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
            LDT_rc%gridDesc(n,8)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
            LDT_rc%gridDesc(n,9)))
       
    elseif(trim(LDT_rc%lis_map_proj(n)).eq."lambert") then !lambert conformal
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
            "LAMBERT CONFORMAL"))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LAT", &
            LDT_rc%gridDesc(n,4)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LON", &
            LDT_rc%gridDesc(n,5)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
            LDT_rc%gridDesc(n,10)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
            LDT_rc%gridDesc(n,7)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
            LDT_rc%gridDesc(n,11)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
            LDT_rc%gridDesc(n,8)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
            LDT_rc%gridDesc(n,9)))
       
    elseif(trim(LDT_rc%lis_map_proj(n)).eq."polar") then ! polar stereographic
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
            "POLAR STEREOGRAPHIC"))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LAT", &
            LDT_rc%gridDesc(n,4)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LON", &
            LDT_rc%gridDesc(n,5)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
            LDT_rc%gridDesc(n,10)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
            LDT_rc%gridDesc(n,7)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
            LDT_rc%gridDesc(n,11)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
            LDT_rc%gridDesc(n,8)))
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
            LDT_rc%gridDesc(n,9)))
    endif
    
    ios = nf90_enddef(ftn)

    call LDT_verify(nf90_put_var(ftn,xtimeID,0.0))
    
    ios = nf90_put_var(ftn,latId,lat)
    call LDT_verify(ios, "nf90_put_var failed for lat")
    ios = nf90_put_var(ftn,lonId,lon)
    call LDT_verify(ios, "nf90_put_var failed for lon")
    ios = nf90_put_var(ftn,varID,value_2d)
    call LDT_verify(ios, "nf90_put_var failed for "//trim(LDT_ANNoutput%varName(1)))
    ios = nf90_close(ftn)
#endif       

  end subroutine PredictANN

!BOP
! !ROUTINE: readANNtrainingOutput
! \label{readANNtrainingOutput}
!
! !INTERFACE: 
  subroutine readANNtrainingOutput(n)
! !ARGUMENTS: 
    integer, intent(in)   :: n 
! 
! !DESCRIPTION: 
!  This subroutine reads the training outputs for use in prediction
! 
!EOP
    integer           :: iret
    integer           :: ftn
    integer           :: dimID(3)
    integer           :: tdimID(2)
    integer           :: wiID, woID
    integer           :: minrangeoutId, maxrangeoutId
    integer           :: minrangeinpId, maxrangeinpId
    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    integer, dimension(8) :: values
    character(len=LDT_CONST_PATH_LEN)         :: outfile

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 

    outfile =  trim(LDT_ANNdata%outfile)         
    if(LDT_masterproc) then 
       iret = nf90_open(path=trim(outfile),&
            mode=NF90_NOWRITE, ncid=ftn)
       call LDT_verify(iret,'Error opening file '//trim(outfile))
    endif

    call LDT_verify(nf90_get_att(ftn,&
         nf90_global, "Output_variables",&
         LDT_ANNoutput%varName(1)))

    call LDT_verify(nf90_inq_varid(ftn,&
         'MINRANGE_INPUTS',&
         varid = minrangeinpId),&
         'nf90_inq_varid failed for MINRANGE_INPUTS')

    call LDT_verify(nf90_inq_varid(ftn,&
         'MAXRANGE_INPUTS',&
         varid = maxrangeinpId),&
         'nf90_inq_varid failed for MAXRANGE_INPUTS')

    call LDT_verify(nf90_inq_varid(ftn,&
         'MINRANGE_OUTPUTS',&
         varid = minrangeoutId),&
         'nf90_inq_varid failed for MINRANGE_OUTPUTS')

    call LDT_verify(nf90_inq_varid(ftn,&
         'MAXRANGE_OUTPUTS',&
         varid = maxrangeoutId),&
         'nf90_inq_varid failed for MAXRANGE_OUTPUTS')

    call LDT_verify(nf90_inq_varid(ftn,&
         'WI',&
         varid = wiId),&
         'nf90_inq_varid failed for WI')

    call LDT_verify(nf90_inq_varid(ftn,&
         'WO',&
         varid = woId),&
         'nf90_inq_varid failed for WO')
    
    call LDT_verify(nf90_get_var(ftn,&
         maxrangeinpId, & 
         LDT_ANNdata%maxinp),&
         'nf90_get_var failed for MAXRANGE_INPUTS')

    call LDT_verify(nf90_get_var(ftn,&
         minrangeinpId, & 
         LDT_ANNdata%mininp),&
         'nf90_get_var failed for MINRANGE_INPUTS')

    call LDT_verify(nf90_get_var(ftn,&
         maxrangeoutId, & 
         LDT_ANNdata%maxout(:,1)),&
         'nf90_get_var failed for MAXRANGE_OUTPUTS')

    call LDT_verify(nf90_get_var(ftn,&
         minrangeoutId, & 
         LDT_ANNdata%minout(:,1)),&
         'nf90_get_var failed for MINRANGE_OUTPUTS')

    if(LDT_ANNdata%trainMode.eq.1) then 
       call LDT_verify(nf90_get_var(ftn,&
            wiId, & 
            LDT_ANNdata%wi_best),&
            'nf90_get_var failed for WI')
       
       call LDT_verify(nf90_get_var(ftn,&
            woId, & 
            LDT_ANNdata%wo_best),&
            'nf90_get_var failed for WO')
    else
       call LDT_verify(nf90_get_var(ftn,&
            wiId, & 
            LDT_ANNdata%wi_best(1,:,:)),&
            'nf90_get_var failed for WI')
       
       call LDT_verify(nf90_get_var(ftn,&
            woId, & 
            LDT_ANNdata%wo_best(1,:)),&
            'nf90_get_var failed for WO')
    endif

    if(LDT_masterproc) then 
       iret = nf90_close(ftn)
       call LDT_verify(iret,'Error in nf90_close')
    endif
    LDT_ANNdata%wi = LDT_ANNdata%wi_best
    LDT_ANNdata%wo = LDT_ANNdata%wo_best
#endif    
  end subroutine readANNtrainingOutput

!BOP
! !ROUTINE: writeANNoutput
! \label{writeANNoutput}
!
! !INTERFACE: 
  subroutine writeANNoutput(n)
! !ARGUMENTS: 
    integer,   intent(in) :: n 
! 
! !DESCRIPTION: 
!  This subroutines writes the output of the ANN training to a NetCDF
!  file
!EOP

    integer               :: iret
    integer               :: ftn
    integer               :: dimID(3)
    integer               :: tdimID(2)
    integer               :: wiID, woID
    integer               :: maxrangeOutID, minrangeOutID
    integer               :: maxrangeInpID, minrangeInpID
    integer               :: inputVarId, inputUnitId
    integer               :: i

    character(len=8)      :: date
    character(len=10)     :: time
    character(len=5)      :: zone
    integer, dimension(8) :: values
    character(len=LDT_CONST_PATH_LEN)         :: outfile
    character*1000        :: allInpVarNames
    character*1000        :: allInpUnitNames
    character*1000        :: allOutVarNames
    character*1000        :: allOutUnitNames

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    call system('mkdir -p '//(LDT_rc%odir))
    write(LDT_logunit,*) "[INFO] Writing to LDT output directory: ",&
          trim(LDT_rc%odir)
    
    outfile = trim(LDT_rc%odir)//'/'//&
         trim(LDT_ANNdata%outfile)         
    if(LDT_masterproc) then 
#if (defined USE_NETCDF3)
       iret = nf90_create(path=trim(outfile),&
            mode=NF90_CLOBBER, ncid=ftn)
#endif
#if (defined USE_NETCDF4)
       iret = nf90_create(path=trim(outfile),&
            cmode=NF90_NETCDF4, ncid=ftn)
#endif
       call LDT_verify(iret,'Error opening file '//trim(outfile))
    endif

    call date_and_time(date,time,zone,values)
    call LDT_verify(nf90_def_dim(ftn,'ngrid',&
         LDT_rc%glbngrid(n),dimID(1)))
    call LDT_verify(nf90_def_dim(ftn,'nppb',&
         LDT_ANNdata%nppb,dimID(2)))
    call LDT_verify(nf90_def_dim(ftn,'nhsize',&
         LDT_ANNdata%nhf,dimID(3)))

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(3)
    
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Activation_function", "Continuous Tan-Sigmoid"))          
    allInpVarNames = ""
    allInpUnitNames = ""
    allOutVarNames = ""
    allOutUnitNames = ""
    do i=1,LDT_ANNdata%nppb-1
       if(i.eq.LDT_ANNdata%nppb-1) then 
          allInpVarnames = trim(allInpVarNames)//&
               trim(LDT_ANNinput%varName(i))
          allInpUnitnames = trim(allInpUnitNames)//&
               trim(LDT_ANNinput%varUnits(i))
       else
          allInpVarnames = trim(allInpVarNames)//&
               trim(LDT_ANNinput%varName(i))//","
          allInpUnitnames = trim(allInpUnitNames)//&
               trim(LDT_ANNinput%varUnits(i))//","
       endif
    enddo
    allOutVarnames = trim(allOutVarNames)//&
         trim(LDT_ANNoutput%varName(1))
    allOutUnitnames = trim(allOutUnitNames)//&
         trim(LDT_ANNoutput%varUnits(1))

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Training mode",&
         "Distributed"))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Input_variables",&
         trim(allInpVarNames)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Input_variable_units ",&
         trim(allInpUnitNames)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Output_variables",&
         trim(allOutVarNames)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Output_variable_units ",&
         trim(allOutUnitNames)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "missing_value", -9999.0))          
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "title", &
         "Land Data Toolkit (LDT) output"))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "institution", &
         "NASA GSFC Hydrological Sciences Laboratory"))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "history", &
         "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
         date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
         "Arsenault_etal_GMD_2018, Kumar_etal_EMS_2006"))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
         "website: http://lis.gsfc.nasa.gov/"))

    call LDT_verify(nf90_def_var(ftn,&
         'MINRANGE_INPUTS',&
         nf90_float, & 
         dimIDs = dimID(1:2),&
         varid = minrangeInpId),&
         'nf90_def_var failed for MINRANGE_INPUTS')

    call LDT_verify(nf90_def_var(ftn,&
         'MAXRANGE_INPUTS',&
         nf90_float, & 
         dimIDs = dimID(1:2),&
         varid = maxrangeInpId),&
         'nf90_def_var failed for MAXRANGE_INPUTS')

    call LDT_verify(nf90_def_var(ftn,&
         'MINRANGE_OUTPUTS',&
         nf90_float, & 
         dimIDs = dimID(1),&
         varid = minrangeOutId),&
         'nf90_def_var failed for MINRANGE_OUTPUTS')

    call LDT_verify(nf90_def_var(ftn,&
         'MAXRANGE_OUTPUTS',&
         nf90_float, & 
         dimIDs = dimID(1),&
         varid = maxrangeOutId),&
         'nf90_def_var failed for MAXRANGE_OUTPUTS')

    call LDT_verify(nf90_def_var(ftn,&
         'WI',&
         nf90_float, & 
         dimIDs = dimID,&
         varid = wiId),&
         'nf90_def_var failed for WI')

    call LDT_verify(nf90_def_var(ftn,&
         'WO',&
         nf90_float, & 
         dimIDs = tdimID,&
         varid = woId),&
         'nf90_def_var failed for WO')
    
    call LDT_verify(nf90_enddef(ftn))
    
    call LDT_verify(nf90_put_var(ftn,&
         maxrangeInpId, & 
         LDT_ANNdata%maxinp,&
         (/1,1/),&
         (/LDT_rc%glbngrid(n),LDT_ANNdata%nppb/)), &
         'nf90_put_var failed for MAXRANGE_INPUTS')

    call LDT_verify(nf90_put_var(ftn,&
         minrangeInpId, & 
         LDT_ANNdata%mininp,&
         (/1,1/),&
         (/LDT_rc%glbngrid(n),LDT_ANNdata%nppb/)), &
         'nf90_put_var failed for MINRANGE_INPUTS')

    call LDT_verify(nf90_put_var(ftn,&
         maxrangeOutId, & 
         LDT_ANNdata%maxout(:,1),&
         (/1/),&
         (/LDT_rc%glbngrid(n)/)), &
         'nf90_put_var failed for MAXRANGE_OUTPUTS')

    call LDT_verify(nf90_put_var(ftn,&
         minrangeOutId, & 
         LDT_ANNdata%minout(:,1),&
         (/1/),&
         (/LDT_rc%glbngrid(n)/)), &
         'nf90_put_var failed for MINRANGE_OUTPUTS')

    call LDT_verify(nf90_put_var(ftn,&
         wiId, & 
         LDT_ANNdata%wi_best,&
         (/1,1,1/),&
         (/LDT_rc%glbngrid(n),LDT_ANNdata%nppb,&
         LDT_ANNdata%nhf/)), &
         'nf90_put_var failed for WI')

    call LDT_verify(nf90_put_var(ftn,&
         woId, & 
         LDT_ANNdata%wo_best,&
         (/1,1/),&
         (/LDT_rc%glbngrid(n),&
         LDT_ANNdata%nhf/)), &
         'nf90_put_var failed for WO')

#if (defined USE_NETCDF3 || defined USE_NETCDF4)    
    if(LDT_masterproc) then 
       iret = nf90_close(ftn)
       call LDT_verify(iret,'Error in nf90_close')
    endif
#endif
#endif

  end subroutine writeANNoutput

!BOP
! !ROUTINE: writeANNoutput_nondist
! \label{writeANNoutput_nondist}
!
! !INTERFACE: 
  subroutine writeANNoutput_nondist(n)
! !ARGUMENTS: 
    integer,   intent(in) :: n 
! 
! !DESCRIPTION: 
!  This subroutines writes the output of the ANN training to a NetCDF
!  file
!EOP

    integer               :: iret
    integer               :: ftn
    integer               :: dimID(3)
    integer               :: tdimID(2)
    integer               :: wiID, woID
    integer               :: maxrangeOutID, minrangeOutID
    integer               :: maxrangeInpID, minrangeInpID
    integer               :: inputVarId, inputUnitId
    integer               :: i

    character(len=8)      :: date
    character(len=10)     :: time
    character(len=5)      :: zone
    integer, dimension(8) :: values
    character(len=LDT_CONST_PATH_LEN)         :: outfile
    character*1000        :: allInpVarNames
    character*1000        :: allInpUnitNames
    character*1000        :: allOutVarNames
    character*1000        :: allOutUnitNames

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    call system('mkdir -p '//(LDT_rc%odir))
    write(LDT_logunit,*) "[INFO] Writing to LDT output directory: ",&
          trim(LDT_rc%odir)
    
    outfile = trim(LDT_rc%odir)//'/'//&
         trim(LDT_ANNdata%outfile)         
    if(LDT_masterproc) then 
#if (defined USE_NETCDF3)
       iret = nf90_create(path=trim(outfile),&
            mode=NF90_CLOBBER, ncid=ftn)
#endif
#if (defined USE_NETCDF4)
       iret = nf90_create(path=trim(outfile),&
            cmode=NF90_NETCDF4, ncid=ftn)
#endif
       call LDT_verify(iret,'Error opening file '//trim(outfile))
    endif

    call date_and_time(date,time,zone,values)
    call LDT_verify(nf90_def_dim(ftn,'ngrid',&
         1,dimID(1)))
    call LDT_verify(nf90_def_dim(ftn,'nppb',&
         LDT_ANNdata%nppb,dimID(2)))
    call LDT_verify(nf90_def_dim(ftn,'nhsize',&
         LDT_ANNdata%nhf,dimID(3)))

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Activation_function", "Continuous Tan-Sigmoid"))          
    allInpVarNames = ""
    allInpUnitNames = ""
    allOutVarNames = ""
    allOutUnitNames = ""
    do i=1,LDT_ANNdata%nppb-1
       if(i.eq.LDT_ANNdata%nppb-1) then 
          allInpVarnames = trim(allInpVarNames)//&
               trim(LDT_ANNinput%varName(i))
          allInpUnitnames = trim(allInpUnitNames)//&
               trim(LDT_ANNinput%varUnits(i))
       else
          allInpVarnames = trim(allInpVarNames)//&
               trim(LDT_ANNinput%varName(i))//","
          allInpUnitnames = trim(allInpUnitNames)//&
               trim(LDT_ANNinput%varUnits(i))//","
       endif
    enddo
    allOutVarnames = trim(allOutVarNames)//&
         trim(LDT_ANNoutput%varName(1))
    allOutUnitnames = trim(allOutUnitNames)//&
         trim(LDT_ANNoutput%varUnits(1))

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Training mode",&
         "Non-distributed"))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Input_variables",&
         trim(allInpVarNames)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Input_variable_units ",&
         trim(allInpUnitNames)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Output_variables",&
         trim(allOutVarNames)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "Output_variable_units ",&
         trim(allOutUnitNames)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "missing_value", -9999.0))          
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "title", &
         "Land Data Toolkit (LDT) output"))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "institution", &
         "NASA GSFC Hydrological Sciences Laboratory"))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
         "history", &
         "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
         date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
         "Arsenault_etal_GMD_2018, Kumar_etal_EMS_2006"))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
         "website: http://lis.gsfc.nasa.gov/"))

    call LDT_verify(nf90_def_var(ftn,&
         'MINRANGE_INPUTS',&
         nf90_float, & 
         dimIDs = dimID(1:2),&
         varid = minrangeInpId),&
         'nf90_def_var failed for MINRANGE_INPUTS')

    call LDT_verify(nf90_def_var(ftn,&
         'MAXRANGE_INPUTS',&
         nf90_float, & 
         dimIDs = dimID(1:2),&
         varid = maxrangeInpId),&
         'nf90_def_var failed for MAXRANGE_INPUTS')

    call LDT_verify(nf90_def_var(ftn,&
         'MINRANGE_OUTPUTS',&
         nf90_float, & 
         dimIDs = dimID(1),&
         varid = minrangeOutId),&
         'nf90_def_var failed for MINRANGE_OUTPUTS')

    call LDT_verify(nf90_def_var(ftn,&
         'MAXRANGE_OUTPUTS',&
         nf90_float, & 
         dimIDs = dimID(1),&
         varid = maxrangeOutId),&
         'nf90_def_var failed for MAXRANGE_OUTPUTS')

    tdimID(1) = dimID(2)
    tdimID(2) = dimID(3)
    
    call LDT_verify(nf90_def_var(ftn,&
         'WI',&
         nf90_float, & 
         dimIDs = tdimID,&
         varid = wiId),&
         'nf90_def_var failed for WI')

    tdimID(1) = dimID(3)
    
    call LDT_verify(nf90_def_var(ftn,&
         'WO',&
         nf90_float, & 
         dimIDs = tdimID(1:1),&
         varid = woId),&
         'nf90_def_var failed for WO')
    
    call LDT_verify(nf90_enddef(ftn))
    
    call LDT_verify(nf90_put_var(ftn,&
         maxrangeInpId, & 
         LDT_ANNdata%maxinp,&
         (/1,1/),&
         (/1,LDT_ANNdata%nppb/)), &
         'nf90_put_var failed for MAXRANGE_INPUTS')

    call LDT_verify(nf90_put_var(ftn,&
         minrangeInpId, & 
         LDT_ANNdata%mininp,&
         (/1,1/),&
         (/1,LDT_ANNdata%nppb/)), &
         'nf90_put_var failed for MINRANGE_INPUTS')

    call LDT_verify(nf90_put_var(ftn,&
         maxrangeOutId, & 
         LDT_ANNdata%maxout(:,1)),&
         'nf90_put_var failed for MAXRANGE_OUTPUTS')

    call LDT_verify(nf90_put_var(ftn,&
         minrangeOutId, & 
         LDT_ANNdata%minout(:,1)),&
         'nf90_put_var failed for MINRANGE_OUTPUTS')

    call LDT_verify(nf90_put_var(ftn,&
         wiId, & 
         LDT_ANNdata%wi_best(1,:,:),&
         (/1,1/),&
         (/LDT_ANNdata%nppb,&
         LDT_ANNdata%nhf/)), &
         'nf90_put_var failed for WI')

    call LDT_verify(nf90_put_var(ftn,&
         woId, & 
         LDT_ANNdata%wo_best(1,:),&
         (/1/),&
         (/LDT_ANNdata%nhf/)), &
         'nf90_put_var failed for WO')

#if (defined USE_NETCDF3 || defined USE_NETCDF4)    
    if(LDT_masterproc) then 
       iret = nf90_close(ftn)
       call LDT_verify(iret,'Error in nf90_close')
    endif
#endif
#endif

  end subroutine writeANNoutput_nondist


  subroutine keep_best_weights(t,j,rmse)
    
    integer,   intent(in)       :: t
    integer,   intent(in)       :: j
    real                        :: rmse

    if(rmse.ne.LDT_rc%udef) then 
       if(j.eq.1) then 
          LDT_ANNdata%rmse_best = rmse
       endif
       
       if(rmse < LDT_ANNdata%rmse_best) then 
          LDT_ANNdata%wi_best(t,:,:) = LDT_ANNdata%wi(t,:,:)
          LDT_ANNdata%wo_best(t,:) = LDT_ANNdata%wo(t,:)
          LDT_ANNdata%rmse_best = rmse
!          print*, 'best ',j, t, rmse
       else
          LDT_ANNdata%wi(t,:,:) = LDT_ANNdata%wi_best(t,:,:)
          LDT_ANNdata%wo(t,:) = LDT_ANNdata%wo_best(t,:)
       endif
    endif

  end subroutine keep_best_weights

  subroutine keep_best_weights_nondist(j,rmse)
    
    integer,   intent(in)       :: j
    real                        :: rmse

    if(rmse.ne.LDT_rc%udef) then 
       if(j.eq.1) then 
          LDT_ANNdata%rmse_best = rmse
       endif
       
       if(rmse < LDT_ANNdata%rmse_best) then 
          LDT_ANNdata%wi_best(1,:,:) = LDT_ANNdata%wi(1,:,:)
          LDT_ANNdata%wo_best(1,:) = LDT_ANNdata%wo(1,:)
          LDT_ANNdata%rmse_best = rmse
       else
          LDT_ANNdata%wi(1,:,:) = LDT_ANNdata%wi_best(1,:,:)
          LDT_ANNdata%wo(1,:) = LDT_ANNdata%wo_best(1,:)
       endif
    endif

  end subroutine keep_best_weights_nondist

  real function calculate_error_nondist(n)
    
    integer             :: n 
    real                :: sqErr
    integer             :: i,t
    real                :: err_this_pat, output_this_pat
    real                :: inputs_this_pat(LDT_ANNdata%nppb)
    integer             :: count

    sqErr = 0.0
    count = 0 

    do t=1,LDT_rc%ngrid(n)
       do i=1, LDT_ANNinput%dataEntry(t)%ntimes
          
          output_this_pat = LDT_ANNoutput%dataEntry(t)%value(1,i)
          inputs_this_pat(:) = LDT_ANNinput%dataEntry(t)%value(:,i)
          err_this_pat = calc_err_pat(t,&
               inputs_this_pat, & 
               output_this_pat, & 
               LDT_ANNdata%wi(1,:,:), & 
               LDT_ANNdata%wo(1,:))
       
          sqErr = sqErr + err_this_pat**2
          count = count + 1
       enddo
    enddo

    if(count.gt.0) then 
       calculate_error_nondist = sqrt(sqErr/count)
    else
       calculate_error_nondist = LDT_rc%udef
    endif

  end function calculate_error_nondist

  real function calculate_error(t)
    
    integer, intent(in) :: t
    
    real                :: sqErr
    integer             :: i
    real                :: err_this_pat, output_this_pat
    real                :: inputs_this_pat(LDT_ANNdata%nppb)

    sqErr = 0.0
    
    do i=1, LDT_ANNinput%dataEntry(t)%ntimes
       
       output_this_pat = LDT_ANNoutput%dataEntry(t)%value(1,i)
       inputs_this_pat(:) = LDT_ANNinput%dataEntry(t)%value(:,i)
       err_this_pat = calc_err_pat(t,&
            inputs_this_pat, & 
            output_this_pat, & 
            LDT_ANNdata%wi(t,:,:), & 
            LDT_ANNdata%wo(t,:))
       
       sqErr = sqErr + err_this_pat**2
    enddo
    
    if(LDT_ANNinput%dataEntry(t)%ntimes.gt.0) then 
       calculate_error = sqrt(sqErr/LDT_ANNinput%dataEntry(t)%ntimes)
    else
       calculate_error = LDT_rc%udef
    endif

  end function calculate_error


  real function calc_err_pat(t, inputs_tp, output_tp, wi, wo)
    integer                 :: t
    real,        intent(in) :: inputs_tp(LDT_ANNdata%nppb)
    real,        intent(in) :: output_tp
    real                    :: wi(LDT_ANNdata%nppb, &
         LDT_ANNdata%nhf)
    real                    :: wo(LDT_ANNdata%nhf)
    real         :: outpredl
    
    LDT_ANNdata%hval(:) = tanh(matmul(transpose(&
         wi),inputs_tp))
    LDT_ANNdata%hval(LDT_ANNdata%nhf) = 1
    
    outpredl = sum(wo*LDT_ANNdata%hval(:))
!    print*,output_tp,& 
!         (outpredl/2.0 + 0.5) * & 
!         (LDT_ANNdata%maxout(t,1)-LDT_ANNdata%minout(t,1)) + & 
!         LDT_ANNdata%minout(t,1)        
    calc_err_pat = (outpredl-output_tp)

  end function calc_err_pat

  subroutine calc_out_pat(t, inputs_tp, output_tp, wi, wo)
    integer                 :: t
    real                    :: inputs_tp(LDT_ANNdata%nppb)
    real                    :: output_tp
    real                    :: wi(LDT_ANNdata%nppb, &
         LDT_ANNdata%nhf)
    real                    :: wo(LDT_ANNdata%nhf)
!    real                    :: inp_tmp(1,LDT_ANNdata%nppb)
!    real                    :: out_tmp(1,1)
!    real                    :: wo_tmp(LDT_ANNdata%nhf,1)

!    inp_tmp(1,:) = inputs_tp(:)
!    wo_tmp(:,1) = wo(:)

!    out_tmp= matmul(tanh(matmul(transpose(inp_tmp),wi)),wo_tmp)
!    output_tp = out_tmp(1,1)

!  The logistic sigmoid function is f(x) = e^x/(1+e^x). The 
!  outputs range from 0 to 1 (interpreted as probabilities). 
!  
!  The tanh function (hyperbolic tangent function) is a rescaling
!  of the logistic sigmoid, such that its outputs range from 
!  -1 to 1. 
!
    LDT_ANNdata%hval(:) = tanh(matmul(transpose(&
         wi),inputs_tp))
    LDT_ANNdata%hval(LDT_ANNdata%nhf) = 1
    
    output_tp = sum(wo*LDT_ANNdata%hval(:))

  end subroutine calc_out_pat

  subroutine ANN_initiateWeights(n)
    integer, intent(in) :: n 

    integer             :: t,j,k
    real                :: rRand

    if(LDT_ANNdata%trainMode.eq.1) then !spatially distributed training
       allocate(LDT_ANNdata%wi(LDT_rc%ngrid(n),&
            LDT_ANNdata%nppb, &
            LDT_ANNdata%nhf))
       
       allocate(LDT_ANNdata%wi_best(LDT_rc%ngrid(n),&
            LDT_ANNdata%nppb, &
            LDT_ANNdata%nhf))
       
       allocate(LDT_ANNdata%wo(LDT_rc%ngrid(n),&
            LDT_ANNdata%nhf))
       
       allocate(LDT_ANNdata%wo_best(LDT_rc%ngrid(n),&
            LDT_ANNdata%nhf))
       
       allocate(LDT_ANNdata%hval(&
            LDT_ANNdata%nhf))
       
       do t =1, LDT_rc%ngrid(n)
          do k=1, LDT_ANNdata%nhf
             call random_number(rRand)
             LDT_ANNdata%wo(t,k) = ((rRand-0.5)*2)/10
             do j=1,LDT_ANNdata%nppb
                call random_number(rRand)
                LDT_ANNdata%wi(t,j,k) = ((rRand-0.5)*2)/10
             enddo
          enddo
       enddo
       LDT_ANNdata%wi_best = LDT_ANNdata%wi
       LDT_ANNdata%wo_best = LDT_ANNdata%wo
    else
       allocate(LDT_ANNdata%wi(1,&
            LDT_ANNdata%nppb, &
            LDT_ANNdata%nhf))
       
       allocate(LDT_ANNdata%wi_best(1,&
            LDT_ANNdata%nppb, &
            LDT_ANNdata%nhf))
       
       allocate(LDT_ANNdata%wo(1,&
            LDT_ANNdata%nhf))
       
       allocate(LDT_ANNdata%wo_best(1,&
            LDT_ANNdata%nhf))
       
       allocate(LDT_ANNdata%hval(&
            LDT_ANNdata%nhf))
       
       do k=1, LDT_ANNdata%nhf
          call random_number(rRand)
          LDT_ANNdata%wo(1,k) = ((rRand-0.5)*2)/10
          do j=1,LDT_ANNdata%nppb
             call random_number(rRand)
             LDT_ANNdata%wi(1,j,k) = ((rRand-0.5)*2)/10
          enddo
       enddo
       LDT_ANNdata%wi_best = LDT_ANNdata%wi
       LDT_ANNdata%wo_best = LDT_ANNdata%wo
    endif

  end subroutine ANN_initiateWeights

#if 0 
  subroutine ANN_createTrainingData(n)
    
    integer, intent(in) :: n 

    rTRPC = 1.0
    rTEPC = rTRPC

    trainCount = 0 

    allocate(LDT_ANNdata%trainingInputs(LDT_rc%ngrid(n),&
         LDT_ANNdata%ntimes, &
         LDT_ANNdata%nppb))
    allocate(LDT_ANNdata%trainingOutputs(LDT_rc%ngrid(n),&
         LDT_ANNdata%ntimes))
    allocate(LDT_ANNdata%inputs_this_pattern(LDT_rc%ngrid(n),&
         LDT_ANNdata%nppb))

    do i=1,LDT_ANNdata%ntimes
       
       call random_number(rRand)
       if((rRand <= rTRPC).and.(trainCount < LDT_ANNdata%ntimes)) then 
          trainCount = trainCount + 1

          LDT_ANNdata%trainingInputs(:,trainCount,1:LDT_ANNdata%nInputs) = & 
               LDT_ANNinput%value(:,1:LDT_ANNdata%nInputs,i)
          LDT_ANNdata%trainingInputs(:,trainCount,LDT_ANNdata%nppb) = 1 
          LDT_ANNdata%trainingOutputs(:,trainCount) = &
               LDT_ANNoutput%value(:,1,i)
       elseif(trainCount < LDT_ANNdata%ntimes) then 
          trainCount = trainCount + 1
          LDT_ANNdata%trainingInputs(:,trainCount,1:LDT_ANNdata%nInputs) = & 
               LDT_ANNinput%value(:,1:LDT_ANNdata%nInputs,i)
          LDT_ANNdata%trainingInputs(:,trainCount,LDT_ANNdata%nppb) = 1 
          LDT_ANNdata%trainingOutputs(:,trainCount) = &
               LDT_ANNoutput%value(:,1,i)

       endif
    enddo

  end subroutine ANN_createTrainingData
#endif

!BOP
! !ROUTINE: ANN_scaleData
!  \label{ANN_scaleData}
! 
! !INTERFACE: 
  subroutine ANN_scaleData(n)
! !ARGUMENTS: 
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This subroutine scales the input data based on the max/min values
!EOP
    integer :: t,p,g,count,gcount
    logical :: check, data_check
    
    type(LDT_ANNdataStruc) :: temp_inp
    type(LDT_ANNdataStruc) :: temp_out   
    

    allocate(temp_inp%dataEntry(LDT_rc%ngrid(n)))
    allocate(temp_out%dataEntry(LDT_rc%ngrid(n)))

    do g=1,LDT_rc%ngrid(n)
       allocate(temp_inp%dataEntry(g)%value(&
            LDT_ANNdata%nppb,&
            LDT_ANNdata%ntimes))

       allocate(temp_out%dataEntry(g)%value(&
            1,LDT_ANNdata%ntimes))

       temp_inp%dataEntry(g)%value = LDT_rc%udef
       temp_out%dataEntry(g)%value = LDT_rc%udef
    enddo

    data_check = .false.
    gcount = 0 

    do g=1,LDT_rc%ngrid(n)
       count = 0
       
       do t=1,LDT_ANNdata%ntimes
          check = .true. 
          do p=1,LDT_ANNdata%nppb
             check = check .and. & 
                  (LDT_ANNinput%dataEntry(g)%value(p,t).ne.LDT_rc%udef).and.&
                  (LDT_ANNoutput%dataEntry(g)%value(1,t).ne.LDT_rc%udef)
             
          enddo
!          if(g.eq.2) then
!          if(LDT_ANNoutput%dataEntry(g)%value(1,t).ne.LDT_rc%udef) then 
!             print*, 'inp ',g,LDT_ANNinput%dataEntry(g)%value(:,t)
!             print*, 'out ',g,LDT_ANNoutput%dataEntry(g)%value(:,t)
!             print*, 'ch ',check,count
!          endif

          if(check) then 
             count = count + 1
             temp_inp%dataEntry(g)%value(:,count) = &
                  LDT_ANNinput%dataEntry(g)%value(:,t)
             temp_out%dataEntry(g)%value(:,count) = &
                  LDT_ANNoutput%dataEntry(g)%value(:,t)
          endif
       enddo
       if(count.gt.10) then !need at least a few data pairs. 
          LDT_ANNinput%dataEntry(g)%nTimes = count
          LDT_ANNoutput%dataEntry(g)%nTimes = count
          data_check = .true. 
          
       else
          LDT_ANNinput%dataEntry(g)%nTimes = 0
          LDT_ANNoutput%dataEntry(g)%nTimes = 0
       endif
       if(LDT_ANNdata%trainMode.eq.0) then 
          gcount = gcount + count
       endif
    enddo

    if(LDT_ANNdata%trainMode.eq.0) then 
       if(gcount > 10) then 
          data_check = .true.
       endif
    endif

    if(.not.data_check) then 
       write(LDT_logunit,*) '[ERR] No valid input/output observation pairs found'
       write(LDT_logunit,*) '[ERR] for ANN training ...'
       call LDT_endrun()
    endif

    do g=1,LDT_rc%ngrid(n)
       deallocate(LDT_ANNinput%dataEntry(g)%value)
       deallocate(LDT_ANNoutput%dataEntry(g)%value)

       allocate(LDT_ANNinput%dataEntry(g)%value(&
            LDT_ANNdata%nppb, &
            LDT_ANNinput%dataEntry(g)%nTimes))
       
       allocate(LDT_ANNoutput%dataEntry(g)%value(&
            1, &
            LDT_ANNinput%dataEntry(g)%nTimes))
    enddo

    do g=1,LDT_rc%ngrid(n)
       do t=1,LDT_ANNinput%dataEntry(g)%nTimes
          LDT_ANNinput%dataEntry(g)%value(:,t) = temp_inp%dataEntry(g)%value(:,t)
          LDT_ANNoutput%dataEntry(g)%value(:,t) = temp_out%dataEntry(g)%value(:,t)
       enddo
    enddo

    do g=1,LDT_rc%ngrid(n)
       deallocate(temp_inp%dataEntry(g)%value)
       deallocate(temp_out%dataEntry(g)%value)
    enddo
    ! scale the data
    if(LDT_ANNdata%trainMode.eq.1) then 
       do g=1,LDT_rc%ngrid(n)
          if(LDT_ANNinput%dataEntry(g)%nTimes.gt.0) then 
             do p=1,LDT_ANNdata%nppb-1
                LDT_ANNdata%maxinp(g,p) = maxval(LDT_ANNinput%dataEntry(g)%value(p,:))
                LDT_ANNdata%mininp(g,p) = minval(LDT_ANNinput%dataEntry(g)%value(p,:))
             enddo
             LDT_ANNdata%maxout(g,1) = maxval(LDT_ANNoutput%dataEntry(g)%value(1,:))
             LDT_ANNdata%minout(g,1) = minval(LDT_ANNoutput%dataEntry(g)%value(1,:))
          endif
       enddo
       
       do g=1,LDT_rc%ngrid(n)
          do t=1,LDT_ANNinput%dataEntry(g)%nTimes
             do p=1,LDT_ANNdata%ninputs
                if((LDT_ANNdata%maxinp(g,p)-LDT_ANNdata%mininp(g,p)).ne.0) then 
                   LDT_ANNinput%dataEntry(g)%value(p,t) = & 
                        ((LDT_ANNinput%dataEntry(g)%value(p,t) - &
                        LDT_ANNdata%mininp(g,p))/&
                        (LDT_ANNdata%maxinp(g,p)-LDT_ANNdata%mininp(g,p)) - 0.5)*2
                endif
             enddo
             if((LDT_ANNdata%maxout(g,1)-LDT_ANNdata%minout(g,1)).ne.0) then 
                LDT_ANNoutput%dataEntry(g)%value(1,t) = & 
                     ((LDT_ANNoutput%dataEntry(g)%value(1,t) - &
                     LDT_ANNdata%minout(g,1))/&
                     (LDT_ANNdata%maxout(g,1)-LDT_ANNdata%minout(g,1)) - 0.5)*2
             endif
          enddo
       enddo
    else
       do g=1,LDT_rc%ngrid(n)
          if(LDT_ANNinput%dataEntry(g)%nTimes.gt.0) then 
             do p=1,LDT_ANNdata%nppb-1
                LDT_ANNdata%maxinp(1,p) = &
                     max(LDT_ANNdata%maxinp(1,p),&
                     maxval(LDT_ANNinput%dataEntry(g)%value(p,:)))
                LDT_ANNdata%mininp(1,p) = &
                     min(LDT_ANNdata%mininp(1,p),&
                     minval(LDT_ANNinput%dataEntry(g)%value(p,:)))
             enddo
             LDT_ANNdata%maxout(1,1) = &
                  max(LDT_ANNdata%maxinp(1,p),&
                  maxval(LDT_ANNoutput%dataEntry(g)%value(1,:)))
             LDT_ANNdata%minout(1,1) = &
                  min(LDT_ANNdata%minout(1,1),&
                  minval(LDT_ANNoutput%dataEntry(g)%value(1,:)))
          endif
       enddo
       
       do g=1,LDT_rc%ngrid(n)
          do t=1,LDT_ANNinput%dataEntry(g)%nTimes
             do p=1,LDT_ANNdata%ninputs
                if((LDT_ANNdata%maxinp(1,p)-LDT_ANNdata%mininp(1,p)).ne.0) then 
                   LDT_ANNinput%dataEntry(g)%value(p,t) = & 
                        ((LDT_ANNinput%dataEntry(g)%value(p,t) - &
                        LDT_ANNdata%mininp(1,p))/&
                        (LDT_ANNdata%maxinp(1,p)-LDT_ANNdata%mininp(1,p)) - 0.5)*2
                endif
             enddo
             if((LDT_ANNdata%maxout(1,1)-LDT_ANNdata%minout(1,1)).ne.0) then 
                LDT_ANNoutput%dataEntry(g)%value(1,t) = & 
                     ((LDT_ANNoutput%dataEntry(g)%value(1,t) - &
                     LDT_ANNdata%minout(1,1))/&
                     (LDT_ANNdata%maxout(1,1)-LDT_ANNdata%minout(1,1)) - 0.5)*2
             endif
          enddo
       enddo
    endif

    write(LDT_logunit,*) '[INFO] Finished normalizing the data ..'
  end subroutine ANN_scaleData
!BOP
!
! !ROUTINE: LDT_logSingleANNdata
! \label{LDT_logSingleANNdata}
! 
! !INTERFACE:
  subroutine LDT_logSingleANNdata(n, value, pindex, iomode, name, units)
! !USES: 
    use LDT_coreMod, only     : LDT_rc, LDT_domain
    use LDT_logMod,  only     : LDT_logunit, LDT_endrun

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer             :: pindex
    real                :: value(LDT_rc%lnc(n), LDT_rc%lnr(n))
    integer             :: iomode
    character(len=*)    :: name
    character(len=*)    :: units
! 
! !DESCRIPTION: 
!  This subroutine maps the processed observations onto the LDT data
!  structures for future temporal averaging and comparisons. The 
!  data are also filtered using the specified external mask. 
!
!EOP
    integer :: k,i,c,r,gid

    if(iomode.eq.1) then 
       LDT_ANNinput%varName(pindex) = name
       LDT_ANNinput%varUnits(pindex) = units
    else
       LDT_ANNoutput%varName(pindex) = name
       LDT_ANNoutput%varUnits(pindex) = units
    endif

    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)
          if(value(c,r).ne.LDT_rc%udef) then 
             if(LDT_domain(n)%gindex(c,r).ne.-1) then 
                gid = LDT_domain(n)%gindex(c,r)
                if(iomode.eq.1) then 
                   LDT_ANNinput%dataEntry(gid)%value_hold(pindex) = & 
                        LDT_ANNinput%dataEntry(gid)%value_hold(pindex) + & 
                        value(c,r)
                   LDT_ANNinput%dataEntry(gid)%count(pindex) = &
                        LDT_ANNinput%dataEntry(gid)%count(pindex) + 1
                elseif(iomode.eq.2) then 
                   LDT_ANNoutput%dataEntry(gid)%value_hold(1) = & 
                        LDT_ANNoutput%dataEntry(gid)%value_hold(1) + & 
                        value(c,r)
                   LDT_ANNoutput%dataEntry(gid)%count(1) = &
                        LDT_ANNoutput%dataEntry(gid)%count(1) + 1
                endif
             endif
          endif
       enddo
    enddo
  end subroutine LDT_logSingleANNdata

!BOP
! !ROUTINE: LDT_resetANN
! \label{LDT_resetANN}
! 
! !INTERFACE: 
  subroutine LDT_resetANN(n)
! !ARGUMENTS: 
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This subroutine temporally averages the datasets for use in ANN
!  training or prediction. 
!EOP

    type(ESMF_Time)         :: currTime, startTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: tindex
    integer                 :: c,r,pindex, gid
    integer                 :: rc, status

    if(mod(real(LDT_rc%hr)*3600+60*real(LDT_rc%mn)+float(LDT_rc%ss),&
            real(LDT_rc%tavgInterval)).eq.0) then        
       
       if(LDT_ANNdata%mode.eq."training") then 
          call ESMF_TimeSet(currTime, yy=LDT_rc%yr, &
               mm=LDT_rc%mo,dd=LDT_rc%da,h=LDT_rc%hr,&
               m=LDT_rc%mn,s=LDT_rc%ss,calendar=LDT_calendar,&
               rc=status)
          call LDT_verify(status, 'error in ESMF_TimeSet: in LDT_ANNMod')
          
          call ESMF_ClockGet(LDT_clock, startTime = startTime)
          
          call ESMF_TimeIntervalSet(timeStep, s = nint(LDT_rc%tavgInterval),&
               rc=rc)
          
          tindex = ((currTime - startTime)/timestep)+1
       else
          tindex = 1
       endif

       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             do pindex=1,LDT_ANNdata%nInputs
                if(LDT_domain(n)%gindex(c,r).ne.-1) then 
                   gid = LDT_domain(n)%gindex(c,r)
                   LDT_ANNinput%dataEntry(gid)%count = 0                    
                   LDT_ANNinput%dataEntry(gid)%value_hold(pindex) = 0
                   if(LDT_ANNdata%mode.eq."prediction") then 
                      LDT_ANNinput%dataEntry(gid)%value(pindex,1) = LDT_rc%udef
                   endif
                endif
             enddo
          enddo
       enddo

       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(LDT_domain(n)%gindex(c,r).ne.-1) then 
                gid = LDT_domain(n)%gindex(c,r)
                LDT_ANNoutput%dataEntry(gid)%count(1) = 0 
                LDT_ANNoutput%dataEntry(gid)%value_hold(1) = 0
             endif
          enddo
       enddo
    endif
  end subroutine LDT_resetANN

end module LDT_ANNMod
