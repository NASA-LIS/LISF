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
! !ROUTINE: computeLLestimate
! \label{computeLLestimate}
! 
! !REVISION HISTORY: 
!  
! 15 Jul 2009: Ken Harrison and Soni Yatheendradas; Initial implementation
! 
! !INTERFACE: 
 subroutine computeLLestimate()
! !USES: 
   use ESMF
   use LIS_coreMod,         only : LIS_rc
   use LIS_optUEMod, only : LIS_ObjectiveFunc
   use LIS_logMod,          only : LIS_logUnit, LIS_verify
! 
! !DESCRIPTION: 
!  This routine computes the objective function values to be used in the 
!  optimization algorithm. The routine is required to specify both a
!  maximization and a minimization criteria
! 
!EOP
   implicit none
   
   type(ESMF_Field)      :: minField, maxField, lnLLField
   real, pointer         :: minvalue(:), maxvalue(:), lnLL(:)
   integer               :: t
   integer               :: n 
   integer               :: status
   integer               :: t1, m
   real, allocatable         :: sumval(:)
   
   n = 1
   allocate(sumval(LIS_rc%ntiles(n)/LIS_rc%nensem(n)))

   call ESMF_StateGet(LIS_ObjectiveFunc,"Objective Function Value",lnLLField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(lnLLField, localDE=0, farrayPtr=lnLL, rc=status)
   call LIS_verify(status)

   call ESMF_StateGet(LIS_ObjectiveFunc,"Min Criteria Value",minField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(minField, localDE=0, farrayPtr=minvalue, rc=status)
   call LIS_verify(status)

   call ESMF_StateGet(LIS_ObjectiveFunc,"Max Criteria Value",maxField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(maxField, localDE=0, farrayPtr=maxvalue, rc=status)
   call LIS_verify(status)

   maxvalue = lnLL

#if 0 
!  NOTE: For grid square, maybe better to assign sqerr as udef in computeLSestimate, and catch this everywhere else in the code (obs etc.)?
   do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
      sumval(t) = -1E20
      do m=1,LIS_rc%nensem(n)
         t1 = (t-1)*LIS_rc%nensem(n)+m
         if(lnLL(t1).gt.sumval(t)) &
              sumVal(t) = lnLL(t1)
      enddo

      do m=1,LIS_rc%nensem(n)
         t1 = (t-1)*LIS_rc%nensem(n)+m
         if(lnLL(t1).ne.0) then 
            minvalue(t1) = -exp(lnLL(t1)-sumval(t))
            !      maxvalue(t) = -lnLL(t)
            maxvalue(t1) = exp(lnLL(t1)-sumval(t))
            lnLL(t1) = exp(lnLL(t1)-sumval(t))
         else
            minvalue(t1) = LIS_rc%udef
            maxvalue(t1) = LIS_rc%udef
         endif
      enddo
   enddo
#endif   
   deallocate(sumval)

 end subroutine computeLLestimate
