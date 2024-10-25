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
module LDT_metforcingParmsMod
!BOP
!
! !MODULE: LDT_metforcingParmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read forcing 
!   data (e.g., can be used for forcing downscaling, etc.). 
!
!  \subsubsection{Overview}
!  The routines in this module provide capabilities to read the 
!  forcing data and allows the users to specify the parmeters.
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  28 Oct 2012: Kristi Arsenault; Expanded for forcing parameter datasets
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_forcingparms_init    !allocates memory for required structures
  public :: LDT_forcingparms_writeHeader
  public :: LDT_forcingparms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_force_struc

!- Forcing-specific parameters:
  type, public :: force_type_dec

     type(LDT_paramEntry) :: climppt
     type(LDT_paramEntry) :: climtmin
     type(LDT_paramEntry) :: climtmax
     type(LDT_paramEntry) :: forcelev
     type(LDT_paramEntry) :: forcelevdiff
  
  end type force_type_dec
  type(force_type_dec),    allocatable :: LDT_force_struc(:,:)

contains

!BOP
! 
! !ROUTINE: LDT_forcingparms_init
! \label{LDT_forcingparms_init}
! 
! !INTERFACE:
  subroutine LDT_forcingparms_init

! !USES:
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: LDT_climateOptChecks, &
                       LDT_gridOptChecks
    use LDT_metforcing_pluginMod, only : LDT_metforcing_plugin

! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! forcing parameter datasets, like elevation or terrain fields. 
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[initmetforc](\ref{initmetforc}) \newline
!    invokes the generic method in the registry to define the
!    native domain of the met forcing scheme
!  \item[readforcelev](\ref{readforcelev}) \newline
!    calls the registry to invoke the reading of forcing elevation files. 
!  \end{description}
!
!EOP
    implicit none
    integer  :: n
    integer  :: m, mp
    integer  :: rc

    character(20) :: gdas_check
    character(20) :: ecmwf_check

! __________________________________________________________________

!- Read in config file entries:

   if( LDT_rc%nmetforc > 0 ) then

     write(LDT_logunit,*)" - - - - - - - - Meteorological Parameters - - - - - - - -"
  
  !- Initialize Metforcing parameter domains (call initmetforc):
     call LDT_metforcing_plugin

   ! Read in each meteorological forcing grid info and arrays:
     LDT_rc%met_gridDesc = 0.
     LDT_rc%met_nc = 0
     LDT_rc%met_nr = 0
     LDT_rc%met_proj = "none"

   ! Check if GDAS and/or ECMWF forcings are present:
     mp = 0
     do m = 1, LDT_rc%nmetforc

     !- Account for GDAS' multi-grids:
        if( LDT_rc%metforc(m) == "GDAS" ) then
           mp = mp + 1
           call initmetforc("GDAS"//char(0),mp)
           LDT_rc%met_ecor_parms(mp+1:mp+5) =  LDT_rc%met_ecor(m)
           LDT_rc%met_gridtransform_parms(mp+1:mp+5) = LDT_rc%met_gridtransform(m)
           mp = mp + 5

     !- Account for ECMWF fields:
        elseif( LDT_rc%metforc(m) == "ECMWF" ) then
           mp = mp + 1
           call initmetforc("ECMWF"//char(0),mp)
           LDT_rc%met_ecor_parms(mp+1:mp+7) =  LDT_rc%met_ecor(m)
           LDT_rc%met_gridtransform_parms(mp+1:mp+7) = LDT_rc%met_gridtransform(m)
           mp = mp + 7

     !- Account for all other forcing sources:
        else
           mp = mp + 1
!           call initmetforc(trim(LDT_rc%metforc(m))//char(0),m)
           call initmetforc(trim(LDT_rc%metforc(m))//char(0),mp)
           LDT_rc%met_ecor_parms(mp) =  LDT_rc%met_ecor(m)
           LDT_rc%met_gridtransform_parms(mp) = LDT_rc%met_gridtransform(m)
        endif
     end do

  !- Allocate fields for both LIS running and forcing domains:
     do m = 1, LDT_rc%nmetforc_parms
        do n = 1, LDT_rc%nnest

          write(LDT_logunit,*)" Projection and Transform: ",&
               trim(LDT_rc%met_proj(m)),", ",trim(LDT_rc%met_gridtransform_parms(m))

!          write(*,*) 'forcing: ',m, LDT_rc%metforc_parms(m)
!          write(*,*) 'proj: ', LDT_rc%met_proj(m)
!          write(*,*) 'nc,nr: ', LDT_rc%met_nc(m), LDT_rc%met_nr(m)
!          write(*,*) 'transform: ',LDT_rc%met_gridtransform_parms(m)
!          write(*,*) 'ecor: ',LDT_rc%met_ecor_parms(m)
!          write(*,*) 'griddesc: ',LDT_rc%met_gridDesc_parms(m,:)

       !- If Elevation correction is turned on, read in forcing elevation:
          if( LDT_rc%met_ecor_parms(m) == "lapse-rate" ) then

             write(LDT_logunit,*) "Reading forcing elev for "//trim(LDT_rc%metforc_parms(m))

          !- Forcing elevation height fields (on run domain target):
             LDT_force_struc(n,m)%forcelev%selectOpt = 1
             LDT_force_struc(n,m)%forcelev%source = trim(LDT_rc%metforc_parms(m))
             LDT_force_struc(n,m)%forcelev%short_name =  &
                 "ELEV_"//trim(LDT_rc%metforc_parms(m))
             LDT_force_struc(n,m)%forcelev%standard_name = &
                 "Forcing elevation for "//trim(LDT_rc%metforc_parms(m))
             LDT_force_struc(n,m)%forcelev%units = "m"
             LDT_force_struc(n,m)%forcelev%num_bins = 1
             LDT_force_struc(n,m)%forcelev%vlevels  = 1
             LDT_force_struc(n,m)%forcelev%num_times = 1
             LDT_force_struc(n,m)%forcelev%valid_min = 0.
             LDT_force_struc(n,m)%forcelev%valid_max = 0.

             allocate(LDT_force_struc(n,m)%forcelev%value(&
                 LDT_rc%lnc(n), LDT_rc%lnr(n),&
                 LDT_force_struc(n,m)%forcelev%vlevels))
             LDT_force_struc(n,m)%forcelev%value = LDT_rc%udef

          !- Optional elevation difference field:

           ! Fill in derived parameter entries:
           ! ( input_parmattribs -> output_parmattribs ) 
             call populate_param_attribs( &
!                   "ELEVDIFF_"//trim(LDT_rc%metforc_parms(m)), &
                   "ELEVDIFF_"//trim(LDT_rc%metforc_parmsrc(m)), &
!                   "Forcing elevation diff for "//trim(LDT_rc%metforc_parms(m)),&
                   "Forcing elevation diff for "//trim(LDT_rc%metforc_parmsrc(m)),&
                   "m", LDT_force_struc(n,m)%forcelev, &
                    LDT_force_struc(n,m)%forcelevdiff )

             allocate(LDT_force_struc(n,m)%forcelevdiff%value(&
                 LDT_rc%met_nc(m), LDT_rc%met_nr(m), &
                 LDT_force_struc(n,m)%forcelevdiff%vlevels))       
             LDT_force_struc(n,m)%forcelevdiff%value = LDT_rc%udef

          !- Account for GDAS fields:
             gdas_check  = LDT_rc%metforc_parms(m)
             ecmwf_check = LDT_rc%metforc_parms(m)

             if( gdas_check(1:6) == "GDAS_T" ) then
               call readforcelev( "GDAS"//char(0),&
                    n, m, LDT_force_struc(n,m)%forcelev%value(:,:,1), &
                    LDT_force_struc(n,m)%forcelevdiff%value(:,:,1) )
 
          !- Account for ECMWF fields:
             elseif( ecmwf_check(1:7) == "ECMWF_S" ) then
               call readforcelev( "ECMWF"//char(0),&
                    n, m, LDT_force_struc(n,m)%forcelev%value(:,:,1), &
                    LDT_force_struc(n,m)%forcelevdiff%value(:,:,1) )

          !- Account for all other forcing sources:
             else
               call readforcelev( trim(LDT_force_struc(n,m)%forcelev%source)//char(0),&
                    n, m, LDT_force_struc(n,m)%forcelev%value(:,:,1), &
                    LDT_force_struc(n,m)%forcelevdiff%value(:,:,1) )
                  ! optional field ... elevdiff
             endif
 
             write(LDT_logunit,*) "Done reading forcing elevation field. "

          else
             write(LDT_logunit,*) "[WARN] "//trim(LDT_rc%metforc_parms(m))//&
                " forcing elevation is turned on, but the "
             write(LDT_logunit,*) " 'Topographic correction method' is set to 'none'."
             write(LDT_logunit,*) " Change this option to 'lapse-rate' to write out "       
             write(LDT_logunit,*) " the metforcing elevation field. "

          end if ! End elevation check
        end do   ! end nest loop
     end do      ! end forcing type loop
   end if
  
  end subroutine LDT_forcingparms_init


  subroutine LDT_forcingparms_writeHeader(n,ftn,run_dimID,met_dimID)

    integer     :: n 
    integer     :: ftn
    integer     :: run_dimID(3)
    integer     :: met_dimID(LDT_rc%nmetforc_parms,3)

    integer     :: m
    integer     :: t_dimID(3)

    if( LDT_rc%nmetforc > 0 ) then
       do m = 1, LDT_rc%nmetforc_parms

    !- Run domain:
       t_dimID(1) = run_dimID(1)
       t_dimID(2) = run_dimID(2)
       call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
            LDT_force_struc(n,m)%forcelev)

    !- Forcing domain:
       t_dimID(1) = met_dimID(m,1)
       t_dimID(2) = met_dimID(m,2)
       call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
            LDT_force_struc(n,m)%forcelevdiff)

       end do
    end if

  end subroutine LDT_forcingparms_writeHeader

  subroutine LDT_forcingparms_writeData(n,ftn)

    integer  :: n 
    integer  :: ftn
    integer  :: m

    if( LDT_rc%nmetforc > 0 ) then

       do m = 1, LDT_rc%nmetforc_parms

       !- LIS Run domain:
          call LDT_writeNETCDFdata(n, ftn, LDT_force_struc(n,m)%forcelev)

       !- Met forcing domain:
          call LDT_writeNETCDFdata( n, ftn, LDT_force_struc(n,m)%forcelevdiff, &
                                    LDT_rc%met_nc(m), LDT_rc%met_nr(m) )
       end do
    endif

  end subroutine LDT_forcingparms_writeData

end module LDT_metforcingParmsMod
