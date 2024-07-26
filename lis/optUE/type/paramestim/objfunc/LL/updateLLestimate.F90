!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: updateLLestimate
!  \label{updateLLestimate}
! 
! !REVISION HISTORY: 
!  
! 15 Jul 2009: Ken Harrison and Soni Yatheendradas; Initial implementation
!
! !INTERFACE: 
subroutine updateLLestimate()
! !USES: 
  use ESMF
  use LIS_coreMod,         only : LIS_rc, LIS_domain
  use LIS_optUEMod,        only : LIS_ObjectiveFunc, LIS_DecisionSpace
  use LIS_logMod,          only : LIS_verify, LIS_endrun, LIS_logunit
  use LIS_PE_HandlerMod,       only : LIS_PEOBS_State, LIS_PEOBSPred_State

! 
! !DESCRIPTION:
!  This method updates the running sum of the squared error. Here the 
!  squared error is computed using the observations used for parameter
!  estimation and the model simulated values (obspred). 
! 
!EOP
  implicit none

  character*100, allocatable         :: peobsname(:), peobspredname(:)
  type(ESMF_Field)               :: lnLLField, peobsField, peobspredField
  real, pointer                  :: lnLL(:), peobs(:), obspred(:)
  integer                        :: t, index1
  integer                        :: n
  integer                        :: kk, j
  integer                        :: nvars
  character*100,     allocatable     :: vnames(:)
  integer                        :: status
  type(ESMF_Field)               :: varField
  real, parameter                :: PI=3.14159
  real, pointer                  :: sigma(:)
  real, parameter                :: arbitraryLL=100000
  n = 1

  allocate(peobsname(1))
  allocate(peobspredname(1))

  call ESMF_StateGet(LIS_PEOBS_State, itemNameList=peobsname,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_PEOBSPred_State, itemNameList=peobspredname, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_PEOBS_State, trim(peobsname(1)), peobsField, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(peobsField, localDE=0, farrayPtr=peobs, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_PEOBSPred_State, trim(peobspredname(1)), peobspredField, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(peobspredField, localDE=0, farrayPtr=obspred, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_ObjectiveFunc,"Objective Function Value",lnLLField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(lnLLField, localDE=0, farrayPtr=lnLL, rc=status)
  call LIS_verify(status)
  
!retrieving sigma: 

  call ESMF_StateGet(LIS_decisionSpace, itemCount = nvars, rc=status)
  call LIS_verify(status)
  
  allocate(vnames(nvars))

  call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
       rc=status)
  call LIS_verify(status)
  
! Here we assume that there is only one output variable (or only one sigma)
! if there are more than one output variable we are interested in, the
! obspred object also will contain multiple fields. (To be taken care of later....)
!
  kk = 0 
  do j=1,nvars
     if(index(vnames(j),"SIGMA_").eq.1) then !found the string that starts with "SIGMA_"
        kk = kk+1
        
        call ESMF_StateGet(LIS_decisionSpace, trim(vnames(j)),&
             varField, rc=status)
        call LIS_verify(status)
        
        call ESMF_FieldGet(varField, localDE=0,farrayPtr=sigma, rc=status)
        call LIS_verify(status)
     endif
  enddo

  deallocate(vnames)

  if(kk.eq.0) then 
     write(LIS_logunit,*) 'SIGMA is not specified in the decision space attributes file'
     write(LIS_logunit,*) 'Program stopping... '
     call LIS_endrun()
  endif
  
 

  do t=1,LIS_rc%ntiles(n)
     !    normal distribution pdf
     !    1/sqrt(2.pi.sigma**2) * exp (-(x-x')**2/(2*sigma**2))
     !    ln of the normal distribution pdf
     !     print*, t, lnll(t), sigma(t), peobs(index1), obspred(t)
     index1 = LIS_domain(n)%tile(t)%index

     if (sigma(t) >  0) then
        lnLL(t) = lnLL(t) + 0 - 0.5* log(2*PI) -log(sigma(t)) +&
             (-0.5/(sigma(t)**2))*(peobs(index1) - obspred(t))**2
     else  
        ! sigma is negative and will crash
        ! compute any way; will get thrown out later with bounds check
        lnLL(t)=arbitraryLL
     endif
  enddo
  deallocate(peobsname)
  deallocate(peobspredname)

 end subroutine updateLLestimate

