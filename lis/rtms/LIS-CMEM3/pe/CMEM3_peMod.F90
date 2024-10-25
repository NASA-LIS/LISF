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
module cmem3_peMod
!BOP
!
! !MODULE: cmem3_peMod
!
! !DESCRIPTION:
!  This module contains the definitions of the CMEM3 model parameters
!  used in parameter estimation. The data structure is used to expose
!  the RTM parameters to be used in opt/ue. 
!
! !REVISION HISTORY:
!  12 Jan 2012; Sujay Kumar, Initial Code
! !USES:        
  use ESMF
  use LIS_numerRecipesMod, only : LIS_rand_func

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: cmem3_setup_pedecvars
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: cmem3_pe_struc

!EOP
  type, public ::  cmem3_pe_dec 
     integer               :: nparams
     character*40, allocatable :: param_name(:)
     integer     , allocatable :: param_select(:)
     real        , allocatable :: param_min(:)
     real        , allocatable :: param_max(:)
  end type cmem3_pe_dec

  type(cmem3_pe_dec), allocatable :: cmem3_pe_struc(:)

  SAVE
contains

!BOP
! !ROUTINE: cmem3_setup_pedecvars
!  \label{cmem3_setup_pedecvars}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine cmem3_setup_pedecvars(DEC_State, Feas_State)
! !USES:
    use LIS_coreMod,       only : LIS_rc, LIS_config,LIS_vecPatch, LIS_surface, LIS_localPET
    use LIS_logMod,        only : LIS_logunit, LIS_verify
#if (defined RTMS) 
    use CMEM3_Mod,       only : cmem3_struc
#endif
    implicit none
! !ARGUMENTS: 
    character*100               :: decSpaceAttribsFile
    type(ESMF_State)            :: DEC_State
    type(ESMF_State)            :: Feas_State

!
! !DESCRIPTION:
!  
!  This routine determines the list of parameters to be used in parameter
!  estimation, initializes them, and updates the LIS decision space.  
! 
!EOP

    integer                     :: n 
    type(ESMF_ArraySpec)        :: arrspec1
    type(ESMF_Field)            :: varField
    type(ESMF_Field)            :: feasField
    real, pointer               :: vardata(:)
    integer, pointer            :: mod_flag(:)
    integer                     :: i,t,m
    integer                     :: status
    integer, parameter          :: seed_base=-1000
    integer                     :: seed
    real                        :: rand
    integer                     :: NT    
    character*100               :: vname
    integer                     :: count
    integer                     :: gid

#if (defined RTMS) 
    call ESMF_StateGet(Feas_State, "Feasibility Flag", feasField, rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(feasField,localDE=0,farrayPtr=mod_flag,rc=status)
    call LIS_verify(status)
    
!    mod_flag = 0  !now initialized centrally in LIS_optUEMod.F90 as lsms,rtms,etc. can raise flag

    call ESMF_ConfigGetAttribute(LIS_config,decSpaceAttribsFile,&
         label="RTM Decision space attributes file:",rc=status)
    call LIS_verify(status, "RTM Decision space attributes file: not defined")

    allocate(cmem3_pe_struc(LIS_rc%nnest))
    n = 1
    cmem3_pe_struc(n)%nparams = 7

    allocate(cmem3_pe_struc(n)%param_name(cmem3_pe_struc(n)%nparams))
    allocate(cmem3_pe_struc(n)%param_select(cmem3_pe_struc(n)%nparams))
    allocate(cmem3_pe_struc(n)%param_min(cmem3_pe_struc(n)%nparams))
    allocate(cmem3_pe_struc(n)%param_max(cmem3_pe_struc(n)%nparams))

    ! read the attributes file. 
    call LIS_readPEDecSpaceAttributes(decSpaceAttribsFile, &
         cmem3_pe_struc(n)%nparams, &
         cmem3_pe_struc(n)%param_name, &
         cmem3_pe_struc(n)%param_select, &
         cmem3_pe_struc(n)%param_min, &
         cmem3_pe_struc(n)%param_max)

    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    do i=1,cmem3_pe_struc(n)%nparams
       if(cmem3_pe_struc(n)%param_select(i).eq.1) then 
          vname=trim(cmem3_pe_struc(n)%param_name(i))
          varField = ESMF_FieldCreate(arrayspec=arrspec1, &
               grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
               name=vname,&
               rc=status)
          call LIS_verify(status, &
               'problem with fieldcreate in cmem3_setup_pedecvars')
          call ESMF_AttributeSet(varField,'MinRange',&
               cmem3_pe_struc(n)%param_min(i),rc=status)
          call LIS_verify(status, &
               'setting minrange to decspace obj in cmem3_setup_devars')
          call ESMF_AttributeSet(varField,'MaxRange',&
               cmem3_pe_struc(n)%param_max(i),rc=status)
          call LIS_verify(status, &
               'setting maxrange to decspace obj in cmem3_setup_devars')

          call ESMF_StateAdd(DEC_State,(/varField/),rc=status)
          call LIS_verify(status,&
               'stateadd in cmem3_setup_pedecvars')

          call ESMF_StateGet(DEC_State,vname,varField,rc=status)
          call LIS_verify(status)
          
          call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
               rc=status)
          call LIS_verify(status)

          !Put in vardata(:) the noah value

          NT=LIS_rc%npatch(n,LIS_rc%lsm_index)
          if(vname.eq."SRMAX")        then ; do t=1,NT; vardata(t) = cmem3_struc(n)%srmax(t)         ;enddo ;endif 
          if(vname.eq."SRMAX2SRMIN")  then ; do t=1,NT; vardata(t) = cmem3_struc(n)%srmax2srmin(t)   ;enddo ;endif 
          if(vname.eq."D_LEAF")       then ; do t=1,NT; vardata(t) = cmem3_struc(n)%d_leaf(t)        ;enddo ;endif 
          if(vname.eq."BGF_FIXED")    then ; do t=1,NT; vardata(t) = cmem3_struc(n)%bgf_fixed(t)     ;enddo ;endif 
          if(vname.eq."K_LAI2VGF")    then ; do t=1,NT; vardata(t) = cmem3_struc(n)%k_lai2vgf(t)     ;enddo ;endif 
          if(vname.eq."VWC2LAI")      then ; do t=1,NT; vardata(t) = cmem3_struc(n)%vwc2lai(t)       ;enddo ;endif 
          if(vname.eq."M_D")          then ; do t=1,NT; vardata(t) = cmem3_struc(n)%m_d(t)           ;enddo ;endif 

             !Test whether any defaults are out of bounds
             count=0
             do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
                do m=1,LIS_rc%nensem(n)
                   gid=(t-1)*LIS_rc%nensem(n)+m
                   if(  (m.eq.1) &
                        .and. &
                        (count.eq.0) &
                        .and. &
                        ((vardata(gid) .lt. cmem3_pe_struc(n)%param_min(i)) &
                        .or. &
                        (vardata(gid) .gt. cmem3_pe_struc(n)%param_max(i))) ) then
                      count=count+1   
                      print*, '*****************************************************************', '  ', &
                           'WARNING: cmem default value is out of LIS-OPT/UE bounds '                , '  ', &
                           'for ', vname                                                             , '  ', &
                           'at '                                                                     , '  ', &
                           'col: ', LIS_surface(n,LIS_rc%lsm_index)%tile(gid)%col                    , '  ', &
                           'row: ', LIS_surface(n,LIS_rc%lsm_index)%tile(gid)%row                    , '  ', &
                           'default value: ', vardata(gid)                                           , '  ', &
                           'parameter min: ', cmem3_pe_struc(n)%param_min(i)                        , '  ', &
                           'parameter max: ', cmem3_pe_struc(n)%param_max(i)                        , '  ', &   
                           '*****************************************************************'
                      
                   endif
                enddo
             enddo
          endif
       enddo

!random initialization
       if(LIS_rc%decSpaceInitMode.eq.1) then  !random initialization 
          seed=seed_base-LIS_localPet !seed must be negative number
          call LIS_rand_func(seed,rand)  !initialize random seed with negative number

          do i=1,cmem3_pe_struc(n)%nparams
             if(cmem3_pe_struc(n)%param_select(i).eq.1) then 
                vname=trim(cmem3_pe_struc(n)%param_name(i))

                call ESMF_StateGet(DEC_State,vname,varField,rc=status)
                call LIS_verify(status)
                
                call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                     rc=status)
                call LIS_verify(status)
                
                do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
                   do m=1,LIS_rc%nensem(n)
                      if (m.eq.1) then
                         !nothing; leave ensemble 1 with default values
                      else
                         call LIS_rand_func(1,rand)
                         vardata((t-1)*LIS_rc%nensem(n)+m) = &
                              cmem3_pe_struc(n)%param_min(i) &
                              + rand * ( cmem3_pe_struc(n)%param_max(i) - cmem3_pe_struc(n)%param_min(i) )
                      endif
                   enddo
                enddo
             endif
          enddo
       endif
   
       write(LIS_logunit,*) 'Finished setting up CMEM3 decision space '
#endif
     end subroutine cmem3_setup_pedecvars
     
end module cmem3_peMod
   
