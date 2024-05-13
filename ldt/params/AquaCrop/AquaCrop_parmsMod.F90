!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module AquaCrop_parmsMod
!BOP
!
! !MODULE: AquaCrop_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read greenness fraction
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  greenness fraction climatology data and allows the users to 
!  specify the frequency of climatology (in months). 
!  The climatological data is temporally interpolated  
!  between months to the current simulation date. 
!
! !REVISION HISTORY:
!
!  10 May 2024; Michel Becthold, Louise Busschaert, initial implementation
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: AquaCropParms_init    !allocates memory for required structures
  public :: AquaCropParms_writeHeader
  public :: AquaCropParms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: AquaCrop_struc

  type, public :: aquacrop_type_dec


     ! -  AquaCrop LSM-specific:

  end type aquacrop_type_dec

  type(aquacrop_type_dec), allocatable :: AquaCrop_struc(:)



contains

  subroutine AquaCropParms_init(flag)

    integer      :: flag

    call AquaCropParms_init_LIS()
    
  end subroutine AquaCropParms_init

!BOP
! 
! !ROUTINE: AquaCropParms_init_LIS
! \label{AquaCropParms_init_LIS}
! 
! !INTERFACE:
  subroutine AquaCropParms_init_LIS
! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
!    use LDT_paramOptCheckMod, only: LDT_aquacropparmsOptChecks, &
!                       LDT_gridOptChecks,LDT_soilsOptChecks
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the AquaCrop fraction datasets 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[acParmssetup](\ref{acParmssetup}) \newline
!    calls the registry to invoke the noahParms setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer  :: n,i,c,r,m
   integer  :: rc
   real     :: temp
   logical  :: file_exists
   logical  :: check_data
   real, allocatable  :: force_elev(:,:)

! _____________________________________________________________________

   allocate( AquaCrop_struc(LDT_rc%nnest) )

 end subroutine AquaCropParms_init_LIS

 subroutine AquaCropParms_writeHeader(n,ftn,dimID)

    integer   :: n 
    integer   :: ftn
    integer   :: dimID(3)

  end subroutine AquaCropParms_writeHeader

  subroutine AquaCropParms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

  end subroutine AquaCropParms_writeData


!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name, &
                               units, full_name )

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name
   character(len=*),     optional     :: units
   character(len=*),     optional     :: full_name

   character(20) :: unit_temp
   character(100):: name_temp
! ____________________________________________________
    
   if(present(units)) then
      unit_temp = units
   else
      unit_temp = "none"
   endif
   if(present(full_name)) then
      name_temp = full_name
   else
      name_temp = trim(short_name)
   endif
   
   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = 1
   paramEntry%selectOpt = 1
   paramEntry%source = "AquaCrop"
   paramEntry%units = unit_temp
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

end module AquaCrop_parmsMod
