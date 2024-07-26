!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: summa1_main
! \label{summa1_main}
!
! !ROUTINE: summa1_main.F90
! 
! !INTERFACE:
subroutine summa1_main(n)

  USE LIS_coreMod
  USE LIS_logMod
  USE nrtype  
  use globalData
  use summa1_lsmMod
  USE var_lookup
  USE vegPhenlgy_module
  USE module_sf_noahmplsm
  USE mDecisions_module
  USE derivforce_module
  USE coupled_em_module
  USE qTimeDelay_module

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
! 
!  Calls the run routines for the forcing-only option (summa1)
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  integer(i4b),parameter           :: no=0                       ! .false.
  integer(i4b),parameter           :: yes=1                      ! .true.
  integer(i4b)                     :: modelTimeStep
 ! NOTE: this is done because of the check in coupled_em if computeVegFlux changes in subsequent time steps
 !  (if computeVegFlux changes, then the number of state variables changes, and we need to reoranize the data structures)
 ! compute the exposed LAI and SAI and whether veg is buried by snow

  integer(i4b)                     :: iGRU
  integer(i4b)                     :: iHRU,jHRU,kHRU             ! index of the hydrologic response unit
  logical(lgt)                     :: computeVegFluxFlag         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow) 
  character(len=1024)              :: message=''                 ! error message
  integer(i4b)                     :: err=0                      ! error code
  real(dp)                         :: notUsed_canopyDepth        ! NOT USED: canopy depth (m)
  real(dp)                         :: notUsed_exposedVAI         ! NOT USED: exposed vegetation area index (m2 m-2)
  real(dp)                         :: fracHRU                    ! fractional area of a given HRU (-)
  integer(i4b)                     :: nLayers                    ! total number of layers
  real(dp),allocatable             :: zSoilReverseSign(:)        ! height at bottom of each soil layer, negative downwards (m)
  logical(lgt),parameter           :: overwriteRSMIN=.false.     ! flag to overwrite RSMIN
  logical(lgt)                     :: printRestart               ! flag to print a re-start file
  character(len=64)                :: output_fileSuffix=''       ! suffix for the output file


 ! NOTE: this is done because of the check in coupled_em if computeVegFlux changes in subsequent time steps
 !  (if computeVegFlux changes, then the number of state variables changes, and we need to reoranize the data structures)
 ! compute the exposed LAI and SAI and whether veg is buried by snow
 if(LIS_rc%tscount(n)==1)then  
    do iGRU=1,summa1_struc(n)%nGRU
        do iHRU=1,gru_struc(iGRU)%hruCount
    ! get vegetation phenology
           call vegPhenlgy(&
                ! input/output: data structures
                model_decisions,                & ! intent(in):    model decisions
                summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU), & ! intent(in):    type of vegetation and soil
                summa1_struc(n)%attrStruct%gru(iGRU)%hru(iHRU), & ! intent(in):    spatial attributes
                summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU), & ! intent(in):    model parameters
                summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU), & ! intent(in):    model prognostic variables for a local HRU
                summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU), & ! intent(inout): model diagnostic variables for a local HRU
                ! output
                computeVegFluxFlag,             & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                notUsed_canopyDepth,            & ! intent(out): NOT USED: canopy depth (m)
                notUsed_exposedVAI,             & ! intent(out): NOT USED: exposed vegetation area index (m2 m-2)
                err,message)                      ! intent(out): error control
           call LIS_verify(err,message)
           
    ! save the flag for computing the vegetation fluxes
           if(computeVegFluxFlag)      summa1_struc(n)%computeVegFlux(iGRU)%hru(iHRU) = yes
           if(.not.computeVegFluxFlag) summa1_struc(n)%computeVegFlux(iGRU)%hru(iHRU) = no
           
    ! define the green vegetation fraction of the grid box (used to compute LAI)
           summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU)%var(iLookDIAG%scalarGreenVegFraction)%dat(1) =summa1_struc(n)%greenVegFrac_monthly(summa1_struc(n)%timeStruct%var(iLookTIME%im))
           

        end do  ! looping through HRUs
     end do  ! looping through GRUs
  endif

 ! ****************************************************************************
 ! (8) loop through HRUs and GRUs
 ! ****************************************************************************

 ! initialize variables
  do iGRU=1,summa1_struc(n)%nGRU

  ! initialize runoff variables
     summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)    = 0._dp  ! surface runoff (m s-1)
     summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)    = 0._dp  ! outflow from all "outlet" HRUs (those with no downstream HRU)
  
  ! initialize baseflow variables
     summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = 0._dp ! recharge to the aquifer (m s-1)
     summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = 0._dp ! baseflow from the aquifer (m s-1)
     summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = 0._dp ! transpiration loss from the aquifer (m s-1)
  
  ! initialize total inflow for each layer in a soil column
     do iHRU=1,gru_struc(iGRU)%hruCount
        summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = 0._dp
     end do

  ! loop through HRUs
     do iHRU=1,gru_struc(iGRU)%hruCount

   ! identify the area covered by the current HRU
        fracHRU =  summa1_struc(n)%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea) / summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
        
   ! assign model layers
   ! NOTE: layer structure is different for each HRU
        gru_struc(iGRU)%hruInfo(iHRU)%nSnow = &
             summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)
        gru_struc(iGRU)%hruInfo(iHRU)%nSoil = &
             summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)
        nLayers                                 = &
             summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nLayers)%dat(1)
   
   ! get height at bottom of each soil layer, negative downwards (used in Noah MP)
        allocate(zSoilReverseSign(gru_struc(iGRU)%hruInfo(iHRU)%nSoil),stat=err); call LIS_verify(err,'problem allocating space for zSoilReverseSign')
        zSoilReverseSign(:) = -summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%iLayerHeight)%dat(gru_struc(iGRU)%hruInfo(iHRU)%nSnow+1:nLayers)
  
   ! get NOAH-MP parameters
        call REDPRM(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),      & ! vegetation type index
             summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%soilTypeIndex),     & ! soil type
             summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%slopeTypeIndex),    & ! slope type index
             zSoilReverseSign,                                                & ! * not used: height at bottom of each layer [NOTE: negative] (m)
             gru_struc(iGRU)%hruInfo(iHRU)%nSoil,                             & ! number of soil layers
             urbanVegCategory)                                                  ! vegetation category for urban areas


   ! deallocate height at bottom of each soil layer(used in Noah MP)
        deallocate(zSoilReverseSign,stat=err); call LIS_verify(err,'problem deallocating space for zSoilReverseSign')
   
   ! overwrite the minimum resistance
        if(overwriteRSMIN) RSMIN = &
             summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%minStomatalResistance)%dat(1)
  
   ! overwrite the vegetation height
        HVT(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)) = &
             summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%heightCanopyTop)%dat(1)
        HVB(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)) = &
             summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%heightCanopyBottom)%dat(1)
        
   ! overwrite the tables for LAI and SAI
        if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
           SAIM(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),:) = &
                summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%winterSAI)%dat(1)
           LAIM(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),:) = &
                summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%summerLAI)%dat(1)*summa1_struc(n)%greenVegFrac_monthly
        end if

   ! cycle water pixel
        if (summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex) == isWater) cycle
   
   ! compute derived forcing variables
        call derivforce(summa1_struc(n)%timeStruct%var,                    & ! vector of time information
             summa1_struc(n)%forcStruct%gru(iGRU)%hru(iHRU)%var,& ! vector of model forcing data
             summa1_struc(n)%attrStruct%gru(iGRU)%hru(iHRU)%var,& ! vector of model attributes
             summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU),    & ! vector of model parameters
             summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model diagnostic variables
             summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model fluxes
             err,message)                         ! error control
        call LIS_verify(err,message)
  
   ! ****************************************************************************
   ! (9) run the model
   ! ****************************************************************************
   ! set the flag to compute the vegetation flux
        computeVegFluxFlag = (summa1_struc(n)%computeVegFlux(iGRU)%hru(iHRU) == yes)

   !print*, 'iHRU = ', iHRU
 
   ! initialize the number of flux calls
        summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU)%var(iLookDIAG%numFluxCalls)%dat(1) = 0._dp

   ! run the model for a single parameter set and time step
        call coupled_em(&
             ! model control
             gru_struc(iGRU)%hruInfo(iHRU)%hru_id,    & ! intent(in):    hruId
             summa1_struc(n)%dt_init(iGRU)%hru(iHRU),                 & ! intent(inout): initial time step
             computeVegFluxFlag,                          & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
             ! data structures (input)
             summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    local classification of soil veg etc. for each HRU
             summa1_struc(n)%attrStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    local attributes for each HRU
             summa1_struc(n)%forcStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    model forcing data
             summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    model parameters
             summa1_struc(n)%bvarStruct%gru(iGRU),                        & ! intent(in):    basin-average model variables
             ! data structures (input-output)
             summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model indices
             summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model prognostic variables for a local HRU
             summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model diagnostic variables for a local HRU
             summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model fluxes for a local HRU
                   ! error control
             err,message)            ! intent(out): error control
        call LIS_verify(err,message)
        
   ! update layer numbers that could be changed in coupled_em()
        gru_struc(iGRU)%hruInfo(iHRU)%nSnow = &
             summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)
        gru_struc(iGRU)%hruInfo(iHRU)%nSoil = &
             summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)
        
!   ! check feasibiility of certain states
!   call check_icond(nGRU,nHRU,                     & ! number of response units
!                    progStruct,                    & ! model prognostic (state) variables
!                    mparStruct,                    & ! model parameters
!                    indxStruct,                    & ! layer indexes
!                    err,message)                     ! error control
!   call handle_err(err,message)
 
   ! save the flag for computing the vegetation fluxes
        if(computeVegFluxFlag)      summa1_struc(n)%computeVegFlux(iGRU)%hru(iHRU) = yes
        if(.not.computeVegFluxFlag) summa1_struc(n)%computeVegFlux(iGRU)%hru(iHRU) = no

        kHRU = 0
   ! identify the downslope HRU
        dsHRU: do jHRU=1,gru_struc(iGRU)%hruCount
           if(summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%downHRUindex) == summa1_struc(n)%typeStruct%gru(iGRU)%hru(jHRU)%var(iLookTYPE%hruIndex))then
              if(kHRU==0)then  ! check there is a unique match
                 kHRU=jHRU
                 exit dsHRU
              end if  ! (check there is a unique match)
           end if  ! (if identified a downslope HRU)
        end do dsHRU
  
   ! add inflow to the downslope HRU
        if(kHRU > 0)then  ! if there is a downslope HRU
           summa1_struc(n)%fluxStruct%gru(iGRU)%hru(kHRU)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(kHRU)%var(iLookFLUX%mLayerColumnInflow)%dat(:)  + summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerColumnOutflow)%dat(:)
           
           ! increment basin column outflow (m3 s-1)
        else
           summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)   = summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1) + sum(summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerColumnOutflow)%dat(:))
        end if
        
        ! increment basin surface runoff (m s-1)
        summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)    = summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)     + summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)    * fracHRU
        
        ! increment basin-average baseflow input variables (m s-1)
        summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)   + summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSoilDrainage)%dat(1)     * fracHRU
        summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1)  + summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferTranspire)%dat(1) * fracHRU
        
        ! increment aquifer baseflow -- ONLY if baseflow is computed individually for each HRU
   ! NOTE: groundwater computed later for singleBasin
        if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == localColumn)then
           summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  + summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) * fracHRU
        end if
        
   ! increment the model indices
        nLayers = gru_struc(iGRU)%hruInfo(iHRU)%nSnow + gru_struc(iGRU)%hruInfo(iHRU)%nSoil
        summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midSnowStartIndex)%dat(1) = summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midSnowStartIndex)%dat(1) + gru_struc(iGRU)%hruInfo(iHRU)%nSnow
        summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midSoilStartIndex)%dat(1) = summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midSoilStartIndex)%dat(1) + gru_struc(iGRU)%hruInfo(iHRU)%nSoil
        summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midTotoStartIndex)%dat(1) = summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midTotoStartIndex)%dat(1) + nLayers
        summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcSnowStartIndex)%dat(1) = summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcSnowStartIndex)%dat(1) + gru_struc(iGRU)%hruInfo(iHRU)%nSnow+1
        summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcSoilStartIndex)%dat(1) = summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcSoilStartIndex)%dat(1) + gru_struc(iGRU)%hruInfo(iHRU)%nSoil+1
        summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcTotoStartIndex)%dat(1) = summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcTotoStartIndex)%dat(1) + nLayers+1
        
     end do  ! (looping through HRUs)

  ! compute water balance for the basin aquifer
     if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
        call LIS_verify(20,'multi_driver/bigBucket groundwater code not transferred from old code base yet')
     end if

  ! perform the routing
     associate(totalArea => summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1) )

     call qOverland(&
          ! input
          model_decisions(iLookDECISIONS%subRouting)%iDecision,            &  ! intent(in): index for routing method
          summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1),           &  ! intent(in): surface runoff (m s-1)
          summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)/totalArea, &  ! intent(in): outflow from all "outlet" HRUs (those with no downstream HRU)
          summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1),         &  ! intent(in): baseflow from the aquifer (m s-1)
          summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%routingFractionFuture)%dat,             &  ! intent(in): fraction of runoff in future time steps (m s-1)
          summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat,               &  ! intent(in): runoff in future time steps (m s-1)
          ! output
          summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%averageInstantRunoff)%dat(1),           &  ! intent(out): instantaneous runoff (m s-1)
          summa1_struc(n)%bvarStruct%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1),            &  ! intent(out): routed runoff (m s-1)
          err,message)                                                        ! intent(out): error control
     call LIS_verify(err,message)

     end associate
  end do

  print*, 'finished summa1_main'
  stop

end subroutine summa1_main
