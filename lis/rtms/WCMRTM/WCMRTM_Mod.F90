!-----------------------BEGIN---------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module WCMRTM_Mod
!BOP
!
! !MODULE: WCMRTM_Mod
!
! !DESCRIPTION:
!    This module provides the routines to control the execution of 
!    the Water Cloud model (Ulaby et al., 1990). The routine allows to compute
!    both the backscatter in VV and the backscatter in VH. It is written
!    considering a static incidence angle of 37°. This can be improved in cas
!    a local incidence angle for each observation is available. Table of
!    calibrated parameters (different for the two polarizations are added also
!    in the configuration file). Each pixel of the study area needs to have
!    assigned A,B,C,D parameters for each polarization, associated with
!    longitude and latitude of the same pixel (pixel where the soil moisture is
!    equa to nan need to be removed)
!
! !HISTORY:
! 28 Aug 2020: Sara Modanesi
! 26 Mar 2021 Sara Modanesi: Added specifications for forward states
! !USES:        


#if (defined RTMS)

  use ESMF
  use LIS_coreMod
  use LIS_RTMMod
  use LIS_logMod

  implicit none
 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: WCMRTM_initialize
  public :: WCMRTM_f2t
  public :: WCMRTM_run
  public :: WCMRTM_output
  public :: WCMRTM_geometry 
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: wcm_struc
!EOP
  type, public ::  wcm_type_dec 
   
     character(len=256) :: AA_VV_tbl_name !! from NoahMP36_lsmMod.F90
     character(len=256) :: BB_VV_tbl_name
     character(len=256) :: CC_VV_tbl_name
     character(len=256) :: DD_VV_tbl_name

     character(len=256) :: AA_VH_tbl_name
     character(len=256) :: BB_VH_tbl_name
     character(len=256) :: CC_VH_tbl_name
     character(len=256) :: DD_VH_tbl_name

     real, allocatable :: AA_VV(:)
     real, allocatable :: BB_VV(:)
     real, allocatable :: CC_VV(:)
     real, allocatable :: DD_VV(:)

     real, allocatable :: AA_VH(:)
     real, allocatable :: BB_VH(:)
     real, allocatable :: CC_VH(:)
     real, allocatable :: DD_VH(:)

     real, allocatable :: lone(:)
     real, allocatable :: late(:)
     !-------output------------!   
     real, allocatable :: Sig0VV(:)
     real, allocatable :: Sig0VH(:)
  end type wcm_type_dec

  type(wcm_type_dec), allocatable :: wcm_struc(:) 

  SAVE

contains
!BOP
! 
! !ROUTINE: WCMRTM_initialize
! \label{WCMRTM_initialize}
! 
! !INTERFACE:
   subroutine WCMRTM_initialize()
! !USES:

! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for noahMP3.6-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for noahMP3.6 from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readWCMRTMcrd](\ref{readWCMRTMcrd}) \newline
!
!EOP
   implicit none
    
   integer :: rc
   integer :: n,t
   integer :: ierr
   integer , parameter :: OPEN_OK = 0
   character*128 :: message

!allocate memory for nest
   allocate(wcm_struc(LIS_rc%nnest))

   do n=1,LIS_rc%nnest
!allocate memory for all tile in current nest


      allocate(wcm_struc(n)%AA_VV(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VV(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VV(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VV(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%AA_VH(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VH(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VH(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VH(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%lone(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%late(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%Sig0VV(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      allocate(wcm_struc(n)%Sig0VH(LIS_rc%npatch(n,LIS_rc%lsm_index)))

      call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Content")
      call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")

   enddo 

!----------------A,B,C and D parameter tables for VV pol---------------------!
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM AA_VV parameter table:",rc = rc)
   do n=1, LIS_rc%nnest 
      call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VV_tbl_name, rc=rc)
      call LIS_verify(rc, "WCMRTM AA_VV parameter table: not defined")
   enddo
   

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM BB_VV parameter table:",rc = rc)
   do n=1, LIS_rc%nnest   
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VV_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM BB_VV parameter table: not defined")
   enddo
    
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM CC_VV parameter table:",rc = rc)    
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VV_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM CC_VV parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM DD_VV parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VV_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM DD_VV parameter table: not defined")
   enddo



!----------------A,B,C and D parameter tables for VH pol---------------------!
   
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM AA_VH parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VH_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM AA_VH parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM BB_VH parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VH_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM BB_VH parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM CC_VH parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VH_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM CC_VH parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM DD_VH parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VH_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM DD_VH parameter table: not defined")
   enddo
!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VV pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VV_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VV_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VV(t),wcm_struc(n)%lone(t),&
          wcm_struc(n)%late(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VV pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VV_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VV_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VV(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VV pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VV_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VV_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VV(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VV pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VV_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VV_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VV(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo


!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VH pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VH_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VH_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VH(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VH pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VH_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VH_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VH(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VH pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VH_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VH_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VH(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VH pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VH_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VH_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VH(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo

   do n=1,LIS_rc%nnest !added fields to State 26032021        
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VV")
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VH")
   enddo

   end subroutine WCMRTM_initialize

   subroutine add_fields_toState(n, inState,varname) !added subroutine add-fields_toState 26032021

    use LIS_logMod,   only : LIS_verify
    use LIS_coreMod,  only : LIS_vecTile

    implicit none

    integer            :: n
    type(ESMF_State)   :: inState
    character(len=*)   :: varname

    type(ESMF_Field)     :: varField
    type(ESMF_ArraySpec) :: arrspec
    integer              :: status
    real :: sum
    call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    varField = ESMF_FieldCreate(arrayspec=arrSpec, &
         grid=LIS_vecTile(n), name=trim(varname), &
         rc=status)
    call LIS_verify(status, 'Error in field_create of '//trim(varname))

    call ESMF_StateAdd(inState, (/varField/), rc=status)
    call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

   end subroutine add_fields_toState
!!--------------------------------------------------------------------------------

   subroutine add_sfc_fields(n, sfcState,varname)

   implicit none 

   integer            :: n 
   type(ESMF_State)   :: sfcState
   character(len=*)   :: varname

   type(ESMF_Field)     :: varField
   type(ESMF_ArraySpec) :: arrspec
   integer              :: status
   real :: sum
   call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
   call LIS_verify(status)

   varField = ESMF_FieldCreate(arrayspec=arrSpec, & 
         grid=LIS_vecTile(n), name=trim(varname), &
         rc=status)
   call LIS_verify(status, 'Error in field_create of '//trim(varname))
    
   call ESMF_StateAdd(sfcState, (/varField/), rc=status)
   call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

   end subroutine add_sfc_fields


   subroutine WCMRTM_f2t(n)

   implicit none

   integer, intent(in)    :: n 

   end subroutine WCMRTM_f2t


   subroutine WCMRTM_geometry(n)
   implicit none
   integer, intent(in)    :: n

   end subroutine WCMRTM_geometry 
  !Do nothing for now
   subroutine WCMRTM_run(n)
   use LIS_histDataMod
! !USES: 
   implicit none

   integer, intent(in) :: n

   integer             :: t,p
   integer             :: status
   integer             :: col,row
   real                :: A_VV_cal,B_VV_cal,C_VV_cal,D_VV_cal,A_VH_cal,&
                          B_VH_cal, C_VH_cal, D_VH_cal,lon,lat,lon1,lat1
   real, pointer       :: sm(:), lai(:)
   real                :: sigmabare_VV,sigmabare_VH,s0VV_s_db, s0VH_s_dB, &
                        sigmacan_VV, sigmacan_VH,sigmasoil_VV,sigmasoil_VH,&
                        tt_VV, tt_VH

   real                :: theta, ctheta
   real, pointer       :: sig0val(:) !added for forward states 26032021


   theta = 0. !incidence angle in radians (i.e., rad(37° for backscatter or rad(0°) for gamma0)
   ctheta = cos(theta)


!   map surface properties to SFC    
   call getsfcvar(LIS_sfcState(n), "Soil Moisture Content",&
         sm)
   call getsfcvar(LIS_sfcState(n), "Leaf Area Index", &
         lai)

!---------------------------------------------
! Tile loop 
!--------------------------------------------
   do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
       row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
       col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
       lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
       lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
       do p=1,LIS_rc%glbngrid(n)
          lon1= wcm_struc(n)%lone(p)
          lat1= wcm_struc(n)%late(p)
          if (lon1 .eq. lon .and. lat1 .eq. lat) then
             A_VV_cal=wcm_struc(n)%AA_VV(p)
             B_VV_cal=wcm_struc(n)%BB_VV(p)
             C_VV_cal=wcm_struc(n)%CC_VV(p)
             D_VV_cal=wcm_struc(n)%DD_VV(p)        
             A_VH_cal=wcm_struc(n)%AA_VH(p)
             B_VH_cal=wcm_struc(n)%BB_VH(p)
             C_VH_cal=wcm_struc(n)%CC_VH(p)
             D_VH_cal=wcm_struc(n)%DD_VH(p)
          endif
       enddo
       
      if(.not.isNaN(sm(t)).and. sm(t).ne.LIS_rc%udef) then
         !bare soil backscatter in db
          s0VV_s_db=C_VV_cal+D_VV_cal*sm(t)
          s0VH_s_db=C_VH_cal+D_VH_cal*sm(t)
         !bare soil backscatter in linear units      
          sigmabare_VV=10.**(s0VV_s_db/10.)
          sigmabare_VH=10.**(s0VH_s_db/10.)
         !attenuation
          tt_VV=exp(-2.*B_VV_cal*lai(t)/ctheta)
          tt_VH=exp(-2.*B_VH_cal*lai(t)/ctheta)
         !attenuated soil backscatter
          sigmasoil_VV=tt_VV*sigmabare_VV
          sigmasoil_VH=tt_VH*sigmabare_VH
         !vegetation backscatter in linear units
          sigmacan_VV=(1.-tt_VV)*ctheta*(A_VV_cal*lai(t))
          sigmacan_VH=(1.-tt_VH)*ctheta*(A_VH_cal*lai(t))
         !total backscatter
          wcm_struc(n)%Sig0VV(t)=10.*log10(sigmacan_VV+sigmasoil_VV)
          wcm_struc(n)%Sig0VH(t)=10.*log10(sigmacan_VH+sigmasoil_VH)
       else
          wcm_struc(n)%Sig0VV(t)=LIS_rc%udef
          wcm_struc(n)%Sig0VH(t)=LIS_rc%udef

       endif

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VV,value=       &
          wcm_struc(n)%Sig0VV(t),             &
          vlevel=1, unit="dB",direction="-")

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VH,value=       &
          wcm_struc(n)%Sig0VH(t),             &
          vlevel=1, unit="dB",direction="-")
   enddo

   call getsfcvar(LIS_forwardState(n), "WCM_Sig0VV", sig0val) !added for forward states 26032021
   sig0val = wcm_struc(n)%Sig0VV

   call getsfcvar(LIS_forwardState(n),"WCM_Sig0VH", sig0val)
   sig0val = wcm_struc(n)%Sig0VH

   end subroutine WCMRTM_run


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   subroutine getsfcvar(sfcState, varname, var)
! !USES: 
    
   implicit none
    
   type(ESMF_State)      :: sfcState
   type(ESMF_Field)      :: varField
   character(len=*)      :: varname
   real, pointer         :: var(:)
   integer               :: status

   call ESMF_StateGet(sfcState, trim(varname), varField, rc=status)
   call LIS_verify(status, 'Error in StateGet: CMEM3_handlerMod '//trim(varname))
   call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
   call LIS_verify(status, 'Error in FieldGet: CMEM3_handlerMod '//trim(varname))

   end subroutine getsfcvar

!!!!BOP
!!!! !ROUTINE: WCMRTM_output
!!!! \label{WCMRTM_output}
!!!!
!!!! !INTERFACE: 
   subroutine WCMRTM_output(n)
   integer, intent(in) :: n 
   end subroutine WCMRTM_output
#endif
end module WCMRTM_Mod



