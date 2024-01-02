!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_paramDataMod
!BOP
!
! !MODULE: LDT_paramDataMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manages different running
!   domain implementations
!
! !REVISION HISTORY: 
!  02 Apr 2012:  Sujay Kumar;  Initial Specification
! 

  implicit none
  
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_LSMparam_struc      ! Derived datatype for LSM parameters
  public :: populate_param_attribs  ! Routine to populate param entries
                                    !   for non-params attribs table inputs
  public :: LDT_set_param_attribs   ! Routine reads over the parameter attribute entries
!EOP

  type, public :: LDT_paramEntry
     character*50  :: short_name
     integer       :: selectOpt
     character*50  :: source
     character*20  :: units
     integer       :: vid
     integer       :: vlevels
     integer       :: num_times
     integer       :: num_bins
     real          :: valid_min
     real          :: valid_max
     character*100 :: standard_name
     real, allocatable :: value(:,:,:)
     real          :: scale_factor
     real          :: add_offset
     character*100 :: long_name
     integer       :: zlevels   ! Hiroko added 4D level
     real*8, allocatable :: dvalue(:,:,:,:)  !CLM-4.5
  end type LDT_paramEntry

! ----------------------------------

!- LSM-specific parameters:
  type, public :: lsmparam_type_dec
     character*100 :: param_filename
     integer       :: param_file_ftn
     integer       :: xlatid, xlonid, xtimeID
     integer       :: xlatbid, xlonbid
     integer       :: luindexId, sctdomId   ! Id's for luiindex and sctdom

     type(LDT_paramEntry) :: xlat, xlon
     type(LDT_paramEntry) :: xlat_b, xlon_b
     type(LDT_paramEntry) :: dommask
     type(LDT_paramEntry) :: landmask
     type(LDT_paramEntry) :: landmask2
     type(LDT_paramEntry) :: watermask
     type(LDT_paramEntry) :: landcover
     type(LDT_paramEntry) :: lakecover
     type(LDT_paramEntry) :: sfctype     ! Surface model type
     type(LDT_paramEntry) :: regmask     ! Regional/basin based mask
     type(LDT_paramEntry) :: glaciermask
     type(LDT_paramEntry) :: luindex     ! dominent landcover

   ! Soil texture
     type(LDT_paramEntry) :: texture     ! Soil texture map
     type(LDT_paramEntry) :: texture2      
     type(LDT_paramEntry) :: texture3
     type(LDT_paramEntry) :: texture4
     type(LDT_paramEntry) :: texture5
     type(LDT_paramEntry) :: texture6
     type(LDT_paramEntry) :: texture7
     type(LDT_paramEntry) :: texture8
     type(LDT_paramEntry) :: texture9
     type(LDT_paramEntry) :: texture10
     type(LDT_paramEntry) :: texture11
     type(LDT_paramEntry) :: sand
     type(LDT_paramEntry) :: clay
     type(LDT_paramEntry) :: silt
     type(LDT_paramEntry) :: gravel
     type(LDT_paramEntry) :: soilsfgrd   ! Soils tile gridcell fraction
     type(LDT_paramEntry) :: sctdom      ! Dominent soil type

   ! Topographic parameters
     type(LDT_paramEntry) :: elevation
     type(LDT_paramEntry) :: elevfgrd    ! Elev tile gridcell fraction 
     type(LDT_paramEntry) :: slope
     type(LDT_paramEntry) :: slopefgrd   ! Slope tile gridcell fraction
     type(LDT_paramEntry) :: aspect
     type(LDT_paramEntry) :: aspectfgrd  ! Aspect tile gridcell fraction
     type(LDT_paramEntry) :: curvature
     
  end type lsmparam_type_dec
  
  type(lsmparam_type_dec), allocatable :: LDT_LSMparam_struc(:)

! -----------------------

CONTAINS

!BOP
! !ROUTINE:  populate_param_attribs
! \label{populate_param_attribs}
!
! !INTERFACE:
  subroutine populate_param_attribs( short_name, long_name, units, &
                                     paramEntryIn, paramEntryOut)

! !DESCRIPTION:
!   This routine populates an added parameter's attribute entries
!  based on another parameter entered from the attribs config file.
!
! !USES:
   use LDT_coreMod, only : LDT_rc, LDT_config

   character(len=*),    intent(in)    :: short_name
   character(len=*),    intent(in)    :: long_name
   character(len=*),    intent(in)    :: units
   type(LDT_paramEntry),intent(inout) :: paramEntryIn
   type(LDT_paramEntry),intent(inout) :: paramEntryOut

   paramEntryOut%short_name = short_name
   paramEntryOut%standard_name=long_name

   paramEntryOut%selectOpt = paramEntryIn%selectOpt
   paramEntryOut%source    = paramEntryIn%source
   paramEntryOut%num_bins  = paramEntryIn%num_bins
   paramEntryOut%num_times = paramEntryIn%num_times
   paramEntryOut%vlevels   = paramEntryIn%vlevels
   paramEntryOut%zlevels   = 1
   paramEntryOut%units     = units
   paramEntryOut%valid_min = 0.
   paramEntryOut%valid_max = 0.

 end subroutine populate_param_attribs


!BOP
! !ROUTINE:  LDT_set_param_attribs
! \label{LDT_set_param_attribs}
!
! !INTERFACE:
  subroutine LDT_set_param_attribs(rc,paramEntry, short_name, source)

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
    integer                            :: rc
    type(LDT_paramEntry),intent(inout) :: paramEntry
    character(len=*),    intent(in)    :: short_name
    character(len=*),    intent(in)    :: source
    
! ____________________________________________________
    
!- If parameter attribute entry exists:

    if(rc.eq.0) then 
       paramEntry%short_name = trim(short_name)
       paramEntry%vlevels = 1
       paramEntry%zlevels = 1
       
       if(source.eq."none") then 
          paramEntry%selectOpt = 0
          paramEntry%source = "none"
       else
          paramEntry%selectOpt = 1  
          paramEntry%source = trim(source)    
       endif
       paramEntry%units = ""
       paramEntry%num_times =1
       paramEntry%num_bins =1
       paramEntry%standard_name =trim(short_name)
    else
       paramEntry%selectOpt = 0
       paramEntry%short_name ="none"
       paramEntry%vlevels = 0
       paramEntry%zlevels = 0
       paramEntry%source = ""
       paramEntry%units = ""
       paramEntry%num_times =0
       paramEntry%num_bins =0
       paramEntry%standard_name =""
    endif

  end subroutine LDT_set_param_attribs

end module LDT_paramDataMod
