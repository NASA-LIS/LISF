!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LVT_DataStream2Mod
!BOP
! 
! !MODULE: LVT_DataStream2Mod
! \label(LVT_DataStream2Mod)
!
! !INTERFACE:
! 
! !USES:   
  use LVT_histDataMod
  use LVT_coreMod
  use LVT_logMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  The code in this file contains the basic datastructures and 
!  control routines for handling observational data
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_DataStream2Init
  public :: LVT_readObsData
  public :: LVT_tavgObsData
  public :: LVT_resetObsData
  public :: LVT_writeObsData
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
contains

!BOP
! 
! !ROUTINE: LVT_DataStream2Init
! \label{LVT_DataStream2Init}
!
! !INTERFACE:   
  subroutine LVT_DataStream2Init
! 
! !USES:   
    use LVT_obs_pluginMod, only : LVT_obs_plugin
    implicit none
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

    call LVT_obs_plugin
    
    call observationsetup(trim(LVT_rc%obssource)//char(0))

  end subroutine LVT_DataStream2Init

!BOP
! 
! !ROUTINE: LVT_readObsData
! \label{LVT_readObsData}
!
! !INTERFACE: 
  subroutine LVT_readObsData
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the observations and processes them. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    call readObservationSource(trim(LVT_rc%obssource)//char(0))      

  end subroutine LVT_readObsData

!BOP
! 
! !ROUTINE: LVT_writeObsData
! \label{LVT_writeObsData}
!
! !INTERFACE: 
  subroutine LVT_writeObsData
! 
! !USES: 
    use LVT_fileIOMod, only : LVT_create_output_directory, &
         LVT_create_output_filename

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine writes the observations (in a "LIS style") 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer         :: n,i 
    integer         :: ftn 
    character*100   :: lisfile
    type(LVT_metadataEntry), pointer :: obs

    n = 1
! hardcoded to write output in the binary, "LIS style"
    if (LVT_rc%output_src.eq."LSM") then
       call LVT_create_output_directory("SURFACEMODEL",&
            style="3 level hierarchy")
       call LVT_create_output_filename(n,lisfile,&
            model_name="SURFACEMODEL", &
            wout="binary",style="3 level hierarchy")
    elseif (LVT_rc%output_src.eq."Routing") then
       call LVT_create_output_directory("ROUTING",&
            style="3 level hierarchy")
       call LVT_create_output_filename(n,lisfile,&
            model_name="ROUTING", &
            wout="binary",style="3 level hierarchy")
    elseif (LVT_rc%output_src.eq."RTM") then
       call LVT_create_output_directory("RTM",&
            style="3 level hierarchy")
       call LVT_create_output_filename(n,lisfile,&
            model_name="RTM", &
            wout="binary",style="3 level hierarchy")
    elseif (LVT_rc%output_src.eq."Irrigation") then
       call LVT_create_output_directory("IRRIGATION",&
            style="3 level hierarchy")
       call LVT_create_output_filename(n,lisfile,&
            model_name="IRRIGATION", &
            wout="binary",style="3 level hierarchy")
    endif

    ftn = LVT_getNextUnitNumber()
    
    open(ftn,file=(lisfile), form='unformatted')
    
    call LVT_getDataStream2Ptr(obs)
    
    do while(associated(obs))       
       call write_obsEntry(ftn,obs, obs%short_name)
       obs => obs%next
    enddo

    call LVT_releaseUnitNumber(ftn)

  end subroutine LVT_writeObsData


!BOP
! 
! !ROUTINE: write_obsEntry
! \label(write_obsEntry)
!
! !INTERFACE:
  subroutine write_obsEntry(ftn, obsEntry,name)
! 
! !USES:   
    use LVT_historyMod, only : LVT_writevar_bin

    implicit none
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
    integer                 :: ftn
    type(LVT_metadataEntry) :: obsEntry
    character(len=*)        :: name
    integer                 :: t, k,gind
    real                    :: tempvar(LVT_rc%ngrid)

    if(obsEntry%selectOpt.eq.1) then
       do k=1,obsEntry%vlevels
          do t=1,LVT_rc%ngrid
             gind = LVT_domain%gindex(LVT_domain%grid(t)%col,&
                  LVT_domain%grid(t)%row)
             if(obsEntry%count(gind,k).ne.0) then 
                tempvar(t) = obsEntry%value(gind,k)
             else
                tempvar(t) = LVT_rc%udef
             endif
          enddo
          call LVT_writevar_bin(ftn, tempvar)
       enddo
    endif

  end subroutine write_obsEntry

!BOP
! 
! !ROUTINE: LVT_tavgObsData
! \label{LVT_tavgObsData}
!
! !INTERFACE: 
  subroutine LVT_tavgObsData
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine invokes the calls to compute temporal averages of 
!   desired set of observational output variables, based on the specified 
!   temporal averaging frequency
!  
!   The routines invoked are: 
!   \begin{description}
!    \item[tavgSingleObs](\ref{tavgSingleObs})
!     computes the temporal average for a single observational variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    real         :: total_depth, value_temp
    integer      :: gid, nl, k
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_metadataEntry), pointer :: rootmoist, soilmoist
    type(LVT_metadataEntry), pointer :: runoff, qs,qsb

    if(LVT_rc%obssource.ne."none") then 
       if(LVT_rc%computeFlag) then 
          obs => LVT_histData%head_lsmobs_list
          do while(associated(obs))
             call tavgSingleObs(obs)
             obs => obs%next
          enddo
          
          obs => LVT_histData%head_routingobs_list
          do while(associated(obs))
             call tavgSingleObs(obs)
             obs => obs%next
          enddo
          
          obs => LVT_histData%head_rtmobs_list
          do while(associated(obs))
             call tavgSingleObs(obs)
             obs => obs%next
          enddo
          
          if(LVT_rc%output_src.eq."LSM") then 
             !now compute the root zone variables, if necessary
             !weighted by thicknesses.. 
             if(LVT_MOC_ROOTMOIST.gt.0) then 
                rootmoist => LVT_histData%ptr_into_lsmobs_list(&
                     LVT_MOC_ROOTMOIST)%dataEntryPtr

                if(rootmoist%computeVar.eq.1) then 
                   
                   soilmoist => LVT_histData%ptr_into_lsmobs_list(&
                        LVT_MOC_SOILMOIST)%dataEntryPtr
                   
                   total_depth =0 
                   do k=1,LVT_rc%nsmlayers
                      total_depth = total_depth + LVT_rc%smthick(k)
                      if(total_depth.ge.LVT_rc%lis_rz_d) then 
                         nl = k
                         exit
                      endif
                   enddo
                   
                   do gid=1,LVT_rc%npts
                      value_temp = 0 
                      total_depth = 0 
                      do k=1,nl
                         value_temp = value_temp + soilmoist%value(gid,k)*&
                              LVT_rc%smthick(k)
                         total_depth = total_depth + LVT_rc%smthick(k)
                      enddo
                      rootmoist%value(gid,1) = value_temp/total_depth
                      rootmoist%count(gid,1) = &
                           soilmoist%count(gid,1)
                   enddo
                endif
             endif
          endif

          if(LVT_MOC_RUNOFF.gt.0) then 
             runoff => LVT_histData%ptr_into_lsmobs_list(&
                  LVT_MOC_RUNOFF)%dataEntryPtr

             if(runoff%computeVar.eq.1) then 
                if(LVT_MOC_QS.gt.0.and.LVT_MOC_QSB.gt.0) then 
                   qs => LVT_histData%ptr_into_lsmobs_list(&
                     LVT_MOC_QS)%dataEntryPtr
                   qsb => LVT_histData%ptr_into_lsmobs_list(&
                        LVT_MOC_QSB)%dataEntryPtr
                else
                   write(LVT_logunit,*)& 
                        'Please enable both Qs and Qsb in the output to compute total runoff'
                   call LVT_endrun()
                endif
                do gid=1,LVT_rc%npts
                   runoff%value(gid,1) =&
                        qs%value(gid,1) + & 
                        qsb%value(gid,1)
                   runoff%count(gid,1) = &
                        qs%count(gid,1)
                enddo
             endif
          endif

       endif
    endif
  end subroutine LVT_tavgObsData

!BOP
! 
! !ROUTINE: tavgSingleObs
! \label{tavgSingleObs}
!
! !INTERFACE:
  subroutine tavgSingleObs( dataEntry)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine temporally averages the accumulated observational 
!   data for a single variable. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    integer :: ftn
    integer :: form
    type(LVT_metadataEntry) :: dataEntry
!EOP
    integer :: unit_id
    integer :: k,i,c,r,gid

    if(dataEntry%selectOpt.eq.1.and.dataEntry%computeVar.eq.0) then 
!If the variable is computed, then it is done after time averaging. 
       do k=1,dataEntry%vlevels
          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(LVT_domain%gindex(c,r).ne.-1) then 
                   gid = LVT_domain%gindex(c,r)
                   if(dataEntry%count(gid,k).ne.0) then 
                      dataEntry%value(gid,k) = &
                           dataEntry%value(gid,k)/dataEntry%count(gid,k)
                   endif
                   if(dataEntry%std_flag) then 
                      if(dataEntry%count_std(gid,k).ne.0) then 
                         dataEntry%std(gid,k) = &
                              dataEntry%std(gid,k)/dataEntry%count_std(gid,k)
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
  end subroutine tavgSingleObs

!BOP
! 
! !ROUTINE: LVT_resetObsData
! \label{LVT_resetObsData}
!
! !INTERFACE: 
  subroutine LVT_resetObsData
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine reinitializes the data structures that hold the observational
!   data
! 
!   The routines invoked are: 
!   \begin{description}
!    \item[resetSingleObs](\ref{resetSingleObs})
!     resets the datastructures for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    type(LVT_metadataEntry), pointer :: obs

    obs => LVT_histData%head_lsmobs_list
    do while(associated(obs))
       call resetSingleObs(obs)
       obs => obs%next
    enddo
    
    obs => LVT_histData%head_routingobs_list
    do while(associated(obs))
       call resetSingleObs(obs)
       obs => obs%next
    enddo
    
    obs => LVT_histData%head_rtmobs_list
    do while(associated(obs))
       call resetSingleObs(obs)
       obs => obs%next
    enddo

  end subroutine LVT_resetObsData

!BOP
! 
! !ROUTINE: resetSingleObs
! \label{resetSingleObs}
!
! !INTERFACE: 
  subroutine resetSingleObs(dataEntry)
! 
! !USES:   
    implicit none 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine resets the data structures that hold the observational 
!  data and the temporal averaging counters
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
    type(LVT_metadataEntry) :: dataEntry
! 
!EOP

    integer                 :: k 

    if(dataEntry%selectOpt.eq.1) then 
       do k=1,dataEntry%vlevels
          dataEntry%value(:,k) = 0 
          dataEntry%count(:,k) = 0 
          if(dataEntry%std_flag) then 
             dataEntry%count_std(:,k)= 0 
             dataEntry%std(:,k) = 0
          endif
       enddo      
    endif
  end subroutine resetSingleObs
end module LVT_DataStream2Mod
