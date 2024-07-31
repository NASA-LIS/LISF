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
! 
! !ROUTINE: computePestimate
! \label{computePestimate}
! 
! !REVISION HISTORY: 
!  
! 10 Sep 2009: Ken Harrison; Initial implementation
! 
! !INTERFACE: 
 
subroutine computePestimate()
! !USES: 
   use ESMF
   use LIS_coreMod,         only : LIS_rc
   use LIS_optUEMod, only : LIS_ObjectiveFunc, LIS_decisionSpace
   use LIS_logMod,          only : LIS_logUnit, LIS_verify
   use PobjFunc_Mod ! , only : dist_lnprob, dist_ctl
! 
! !DESCRIPTION: 
!  This routine stores the posterior probabilities for use
!  by the optimization and uncertainty estimation algorithms. 
!  The routine is required to specify both a
!  maximization and a minimization criteria.
! 
!EOP
   implicit none
   
!!$  type diststruc
!!$     integer :: dist_ID
!!$     integer :: numparam
!!$     real, allocatable :: param(:)
!!$     real, allocatable :: statparam(:)  
!!$  end type diststruc

   type(ESMF_Field)               :: minField, maxField, lnPField
   type(ESMF_Field), allocatable      :: varField(:)
   character*100, allocatable         :: vnames(:)
   real, pointer                  :: vardata(:) 
   integer                        :: nvars
   integer                        :: i,j,k

   real, pointer                  :: minvalue(:), maxvalue(:), lnP(:)
   integer                        :: t
   integer                        :: n 
   integer                        :: status
   integer                        :: t1, m
   real                           :: sumval
   type(diststruc)                :: d
   real, allocatable                  :: X(:)
   real                           :: Y
   character*100                  :: vname_dist
   integer                        :: kk

   n = 1

   call ESMF_StateGet(LIS_ObjectiveFunc,"Objective Function Value",lnPField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(lnPField, localDE=0, farrayPtr=lnP, rc=status)
   call LIS_verify(status)

   call ESMF_StateGet(LIS_ObjectiveFunc,"Min Criteria Value",minField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(minField, localDE=0, farrayPtr=minvalue, rc=status)
   call LIS_verify(status)

   call ESMF_StateGet(LIS_ObjectiveFunc,"Max Criteria Value",maxField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(maxField, localDE=0, farrayPtr=maxvalue, rc=status)
   call LIS_verify(status)

   !Account for prior here 
   
   call ESMF_StateGet(LIS_decisionSpace, itemCount = nvars, rc=status)
   call LIS_verify(status)
   
   allocate(vardata(LIS_rc%ntiles(n)))
   allocate(vnames(nvars))
   allocate(varField(nvars))

   call ESMF_StateGet(LIS_decisionSpace, itemNameList = vnames, &
        rc=status)
   call LIS_verify(status)

!  Must now consider prior. Iterate through distributions and pass in param values to obtain 
!  probability (density)


   do t=1,LIS_rc%ntiles(n)
      sumval=0
      do j=1,dist_ctl%ndists
         d=dist_struc(j)

         allocate(X(d%numparam))
         
         do k=1,d%numparam
!            i=d%param(k)
            vname_dist = d%param_name(k)
            do kk=1,nvars
               if(trim(vname_dist).eq.trim(vnames(kk))) then
                  call ESMF_StateGet(LIS_decisionSpace, trim(vnames(kk)),&
                       varField(kk), rc=status)
                  call LIS_verify(status)
               
                  call ESMF_FieldGet(varField(kk), localDE=0,farrayPtr=vardata, &
                       rc=status)
                  call LIS_verify(status)
                  X(k)=vardata(t)
                  exit;
               endif
            enddo
         enddo
         call dist_lnprob(d,X,Y)
         sumval=sumval+Y

         deallocate(X)
      enddo
      lnP(t)=lnP(t)+sumval
   enddo
   maxvalue = lnP

 end subroutine computePestimate
