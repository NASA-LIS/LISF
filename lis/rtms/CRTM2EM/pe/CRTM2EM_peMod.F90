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
module crtm2em_peMod
!BOP
!
! !MODULE: crtm2em_peMod
!
! !DESCRIPTION:
!  This module contains the definitions of the CRTM2EM model parameters
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
  public :: crtm2em_setup_pedecvars
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: crtm2em_pe_struc

!EOP
  type, public ::  crtm2em_pe_dec 
     integer               :: nparams
     character*40, allocatable :: param_name(:)
     integer     , allocatable :: param_select(:)
     real        , allocatable :: param_min(:)
     real        , allocatable :: param_max(:)
  end type crtm2em_pe_dec

  type(crtm2em_pe_dec), allocatable :: crtm2em_pe_struc(:)

  SAVE
contains

!BOP
! !ROUTINE: crtm2em_setup_pedecvars
!  \label{crtm2em_setup_pedecvars}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine crtm2em_setup_pedecvars(DEC_State, Feas_State)
! !USES:
    use LIS_coreMod,       only : LIS_rc, LIS_config,LIS_vecPatch, LIS_surface, LIS_localPET
    use LIS_logMod,        only : LIS_logunit, LIS_verify
#if (defined RTMS) 
    use CRTM2_EMMod,       only : crtm_struc
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

    allocate(crtm2em_pe_struc(LIS_rc%nnest))
    n = 1
    crtm2em_pe_struc(n)%nparams = 6

    allocate(crtm2em_pe_struc(n)%param_name(crtm2em_pe_struc(n)%nparams))
    allocate(crtm2em_pe_struc(n)%param_select(crtm2em_pe_struc(n)%nparams))
    allocate(crtm2em_pe_struc(n)%param_min(crtm2em_pe_struc(n)%nparams))
    allocate(crtm2em_pe_struc(n)%param_max(crtm2em_pe_struc(n)%nparams))

    ! read the attributes file. 
    call LIS_readPEDecSpaceAttributes(decSpaceAttribsFile, &
         crtm2em_pe_struc(n)%nparams, &
         crtm2em_pe_struc(n)%param_name, &
         crtm2em_pe_struc(n)%param_select, &
         crtm2em_pe_struc(n)%param_min, &
         crtm2em_pe_struc(n)%param_max)

    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    do i=1,crtm2em_pe_struc(n)%nparams
       if(crtm2em_pe_struc(n)%param_select(i).eq.1) then 
          vname=trim(crtm2em_pe_struc(n)%param_name(i))
          varField = ESMF_FieldCreate(arrayspec=arrspec1, &
               grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
               name=vname,&
               rc=status)
          call LIS_verify(status, &
               'problem with fieldcreate in crtm2em_setup_pedecvars')
          call ESMF_AttributeSet(varField,'MinRange',&
               crtm2em_pe_struc(n)%param_min(i),rc=status)
          call LIS_verify(status, &
               'setting minrange to decspace obj in crtm2em_setup_devars')
          call ESMF_AttributeSet(varField,'MaxRange',&
               crtm2em_pe_struc(n)%param_max(i),rc=status)
          call LIS_verify(status, &
               'setting maxrange to decspace obj in crtm2em_setup_devars')

          call ESMF_StateAdd(DEC_State,(/varField/),rc=status)
          call LIS_verify(status,&
               'stateadd in crtm2em_setup_pedecvars')

          call ESMF_StateGet(DEC_State,vname,varField,rc=status)
          call LIS_verify(status)
          
          call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
               rc=status)
          call LIS_verify(status)

          !Put in vardata(:) the noah value

          NT=LIS_rc%npatch(n,LIS_rc%lsm_index)
          if(vname.eq."SIGMA")                   then ; do t=1,NT; vardata(t) = crtm_struc(n)%SFC(t)%sigma                        ;enddo ;endif 
          if(vname.eq."LEAF_THICK")              then ; do t=1,NT; vardata(t) = crtm_struc(n)%SFC(t)%leaf_thick                   ;enddo ;endif 
          if(vname.eq."BGF_FIXED")               then ; do t=1,NT; vardata(t) = crtm_struc(n)%bgf_fixed(t)                        ;enddo ;endif 
          if(vname.eq."K_LAI2VGF")               then ; do t=1,NT; vardata(t) = crtm_struc(n)%k_lai2vgf(t)                        ;enddo ;endif    
          if(vname.eq."WATER_CONTENT_PER_LAI")   then ; do t=1,NT; vardata(t) = crtm_struc(n)%SFC(t)%water_content_per_lai        ;enddo ;endif 
          if(vname.eq."SSALB_FACTOR")            then ; do t=1,NT; vardata(t) = crtm_struc(n)%SFC(t)%ssalb_factor                 ;enddo ;endif 

             !Test whether any defaults are out of bounds
             count=0
             do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
                do m=1,LIS_rc%nensem(n)
                   gid=(t-1)*LIS_rc%nensem(n)+m
                   if(  (m.eq.1) &
                        .and. &
                        (count.eq.0) &
                        .and. &
                        ((vardata(gid) .lt. crtm2em_pe_struc(n)%param_min(i)) &
                        .or. &
                        (vardata(gid) .gt. crtm2em_pe_struc(n)%param_max(i))) ) then
                      count=count+1   
                      print*, '*****************************************************************', '  ', &
                           'WARNING: crtm2em default value is out of LIS-OPT/UE bounds '                , '  ', &
                           'for ', vname                                                             , '  ', &
                           'at '                                                                     , '  ', &
                           'col: ', LIS_surface(n,LIS_rc%lsm_index)%tile(gid)%col                    , '  ', &
                           'row: ', LIS_surface(n,LIS_rc%lsm_index)%tile(gid)%row                    , '  ', &
                           'default value: ', vardata(gid)                                           , '  ', &
                           'parameter min: ', crtm2em_pe_struc(n)%param_min(i)                        , '  ', &
                           'parameter max: ', crtm2em_pe_struc(n)%param_max(i)                        , '  ', &   
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

          do i=1,crtm2em_pe_struc(n)%nparams
             if(crtm2em_pe_struc(n)%param_select(i).eq.1) then 
                vname=trim(crtm2em_pe_struc(n)%param_name(i))

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
                              crtm2em_pe_struc(n)%param_min(i) &
                              + rand * ( crtm2em_pe_struc(n)%param_max(i) - crtm2em_pe_struc(n)%param_min(i) )
                      endif
                   enddo
                enddo
             endif
          enddo
       endif
   
       write(LIS_logunit,*) 'Finished setting up CRTM2EM decision space '
#endif
     end subroutine crtm2em_setup_pedecvars
     
end module crtm2em_peMod
   
