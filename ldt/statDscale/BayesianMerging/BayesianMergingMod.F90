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
module BayesianMergingMod
!BOP
!
! !MODULE: BayesianMergingMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  16 Feb 2016: Sujay Kumar; Initial implementation
!
  use ESMF
  use LDT_coreMod
  use LDT_metforcingMod
  use LDT_logMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: BayesianMerging_init
  public :: BayesianMerging_diagnose
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  type, public :: bayesianMergeEntry

     real,    allocatable :: minVal(:,:)
     real,    allocatable :: maxVal(:,:)    
     real,    allocatable :: sx_mu(:,:)
     real,    allocatable :: sxx_sigma(:,:)
     real,    allocatable :: delta(:,:)
     integer,    allocatable :: pdf_bincounts(:,:,:)
  end type bayesianMergeEntry

  type, public :: bayesianMerge_struc
     integer                               :: nseasons
     integer                               :: nVars
     type(bayesianMergeEntry), allocatable :: dataEntry(:,:)
  end type bayesianMerge_struc

  type(bayesianMerge_struc) :: LDT_bayesianMerge_struc

!EOP

contains

  subroutine BayesianMerging_init()

    integer             :: rc
    integer             :: n
    character*100       :: type
    integer             :: i,k,nbins

    n = 1
    
    call ESMF_ConfigGetAttribute(LDT_config,type, &
         label="Bayesian merging seasonal stratification type:",&
         rc=rc)
    call LDT_verify(rc,'Bayesian merging seasonal stratification type: not defined')

    if(type.eq."monthly") then 
       LDT_bayesianMerge_struc%nseasons = 12
    endif
    
!find the minimum number of forcing variables across the forcing instances.
    LDT_bayesianMerge_struc%nVars =10000
    do i=1, LDT_rc%nmetforc
       if(LDT_rc%met_nf(i) .lt. LDT_bayesianMerge_struc%nVars ) then 
          LDT_bayesianMerge_struc%nVars = LDT_rc%met_nf(i)
       endif
    enddo

    allocate(LDT_bayesianMerge_struc%dataEntry(2,LDT_bayesianMerge_struc%nVars))
!TBD
    nbins = 100

    do k=1,2
       do i=1,LDT_bayesianMerge_struc%nVars
          allocate(LDT_bayesianMerge_struc%dataEntry(k,i)%minVal(LDT_rc%ntiles(n),&
               LDT_bayesianMerge_struc%nseasons))
          LDT_bayesianMerge_struc%dataEntry(k,i)%minVal = 1E10
          allocate(LDT_bayesianMerge_struc%dataEntry(k,i)%maxVal(LDT_rc%ntiles(n),&
               LDT_bayesianMerge_struc%nseasons))
          LDT_bayesianMerge_struc%dataEntry(k,i)%maxVal = -1E10
          allocate(LDT_bayesianMerge_struc%dataEntry(k,i)%sx_mu(LDT_rc%ntiles(n),&
               LDT_bayesianMerge_struc%nseasons))
          LDT_bayesianMerge_struc%dataEntry(k,i)%sx_mu = 0
          allocate(LDT_bayesianMerge_struc%dataEntry(k,i)%sxx_sigma(LDT_rc%ntiles(n),&
               LDT_bayesianMerge_struc%nseasons))
          LDT_bayesianMerge_struc%dataEntry(k,i)%sxx_sigma = 0
          allocate(LDT_bayesianMerge_struc%dataEntry(k,i)%delta(LDT_rc%ntiles(n),&
               LDT_bayesianMerge_struc%nseasons))
          LDT_bayesianMerge_struc%dataEntry(k,i)%delta = 0
          allocate(LDT_bayesianMerge_struc%dataEntry(k,i)%pdf_bincounts(LDT_rc%ntiles(n),&
            LDT_bayesianMerge_struc%nseasons,nbins))
          LDT_bayesianMerge_struc%dataEntry(k,i)%pdf_bincounts = 0
          
       end do
    enddo

  end subroutine BayesianMerging_init


  subroutine BayesianMerging_diagnose(n, pass)

    integer :: n
    integer :: pass


    integer :: rc
    integer                    :: i,t,l,fobjcount
    character*100, allocatable :: forcobjs(:)
    real,          pointer     :: forcdata1(:),forcdata2(:)
    type(ESMF_Field)           :: f1Field,f2Field
    integer                    :: varCount,f1flag, f2flag
    integer                    :: binval

!Assume that there are two forcing sources
             
    call ESMF_StateGet(LDT_FORC_Base_State(n,1),itemCount=fobjcount,rc=rc)
    call LDT_verify(rc,'ESMF_StateGet failed for objcount in overlayForcings')
    
    allocate(forcobjs(fobjcount))

    call ESMF_StateGet(LDT_FORC_Base_State(n,1),itemNameList=forcobjs,rc=rc)
    call LDT_verify(rc,'ESMF_StateGet failed for forcobjs1 in overlayForcings')
    
    l = LDT_rc%mo

    varCount = 0

    do i=1,fobjcount

       call ESMF_StateGet(LDT_FORC_Base_State(n,1),forcobjs(i),f1Field,&
            rc=rc)
       
       call ESMF_StateGet(LDT_FORC_Base_State(n,2),forcobjs(i),f2Field,&
            rc=rc)
       
       call ESMF_AttributeGet(f1Field, "Enabled", f1flag,rc=rc)
       call ESMF_AttributeGet(f2Field, "Enabled", f2flag,rc=rc)

       if(f1flag.eq.1.and.f2flag.eq.1) then 
          varCount = varCount + 1

          call ESMF_FieldGet(f1Field,localDE=0,farrayPtr=forcdata1, &
               rc=rc)
          call LDT_verify(rc,'ESMF_FieldGet (ffield) failed in overlayForcings')
          
          call ESMF_FieldGet(f2Field,localDE=0,farrayPtr=forcdata2, &
               rc=rc)
          call LDT_verify(rc,'ESMF_FieldGet (ffield) failed in overlayForcings')
          
          if(pass.eq.1) then 
             !read forcing datasets, save into data structures to compute PDFs.
             do t=1,LDT_rc%ntiles(n)
                if(forcdata1(t).ne.-9999.0) then 
                   
                   if(forcdata1(t).lt.LDT_bayesianMerge_struc%dataEntry(1,varCount)%minVal(t,l)) then 
                      LDT_bayesianMerge_struc%dataEntry(1,varCount)%minVal(t,l) = forcdata1(t)
                   endif
                   if(forcdata1(t).gt.LDT_bayesianMerge_struc%dataEntry(1,varCount)%maxVal(t,l)) then 
                      LDT_bayesianMerge_struc%dataEntry(1,varCount)%maxVal(t,l) = forcdata1(t)
                   endif
                   
                   LDT_bayesianMerge_struc%dataEntry(1,varCount)%sx_mu(t,l) = &
                        LDT_bayesianMerge_struc%dataEntry(1,varCount)%sx_mu(t,l) + forcdata1(t)
                   LDT_bayesianMerge_struc%dataEntry(1,varCount)%sxx_sigma(t,l) = &
                        LDT_bayesianMerge_struc%dataEntry(1,varCount)%sxx_sigma(t,l) + & 
                        forcdata1(t)*forcdata1(t)
                   
                endif
             enddo
             
             do t=1,LDT_rc%ntiles(n)
                if(forcdata2(t).ne.-9999.0) then 
                   
                   if(forcdata2(t).lt.LDT_bayesianMerge_struc%dataEntry(2,varCount)%minVal(t,l)) then 
                      LDT_bayesianMerge_struc%dataEntry(2,varCount)%minVal(t,l) = forcdata2(t)
                   endif
                   if(forcdata2(t).gt.LDT_bayesianMerge_struc%dataEntry(2,varCount)%maxVal(t,l)) then 
                      LDT_bayesianMerge_struc%dataEntry(2,varCount)%maxVal(t,l) = forcdata2(t)
                   endif
                   
                   LDT_bayesianMerge_struc%dataEntry(2,varCount)%sx_mu(t,l) = &
                        LDT_bayesianMerge_struc%dataEntry(2,varCount)%sx_mu(t,l) + forcdata2(t)
                   LDT_bayesianMerge_struc%dataEntry(2,varCount)%sxx_sigma(t,l) = &
                        LDT_bayesianMerge_struc%dataEntry(2,varCount)%sxx_sigma(t,l) + & 
                        forcdata2(t)*forcdata2(t)
                   
                endif
             enddo
#if 0 
          else
             do t=1,LDT_rc%ntiles(n)
                if(forcdata(t).ne.-9999.0) then 
                   binval = nint((forcdata(t) - &
                        LDT_bayesianMerge_struc%dataEntry(varCount)%minVal(t,l))/&
                        LDT_bayesianMerge_struc%dataEntry(varCount)%delta(t,l)) + 1
                   !                bincounts(i,t,l,binval) = bincounts(i,t,l,binval) + 1
                   
                endif
             enddo
#endif
          end if
       endif
    enddo
  end subroutine BayesianMerging_diagnose

end module BayesianMergingMod
