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
!BOP
! !ROUTINE: updateLSestimate
!  \label{updateLSestimate}
! 
! !REVISION HISTORY: 
!  
! 15 Jul 2009: Sujay Kumar; Initial implementation
!
! !INTERFACE: 
subroutine updateLSestimate()
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod,         only : LIS_rc, LIS_masterproc, &
       LIS_localPet, LIS_npes, LIS_domain
  use LIS_optUEMod,        only : LIS_ObjectiveFunc
  use LIS_logMod,          only : LIS_verify
  use LSObjFunc_Mod,       only : ls_ctl 
  use LIS_PE_HandlerMod,   only : LIS_PEOBS_State, LIS_PEOBSPred_State

! 
! !DESCRIPTION:
!  This method updates the running sum of the squared error. Here the 
!  squared error is computed using the observations used for parameter
!  estimation and the model simulated values (obspred). 
! 
!EOP
  implicit none

!  character*100, allocatable         :: peobsname(:)
  type(ESMF_Field)               :: sqerrField, peobsField, peobspredField
  real, pointer                  :: sqerr(:), peobs(:), obspred(:),modelv(:)
!  real, pointer                  :: sqerr_t(:), numobs_t(:)
!  real, pointer                  :: sqerr_t1_sum(:), numobs_t1_sum(:)
  type(ESMF_Field)               :: numobsField
  type(ESMF_Field)     :: modelvField
  type(ESMF_Field)     :: nummodelvField
  real, pointer                  :: numobs(:),nummodelv(:)
  integer                        :: deltas(LIS_npes), offsets(LIS_npes)
  integer                        :: t, m, l,index1, tindex
  integer                        :: n, k,nobjs
  real                           :: wt
  integer                        :: ierr,status

  n = 1


!  allocate(sqerr(LIS_rc%ntiles(n)))
!  allocate(numobs(LIS_rc%ntiles(n)))

!  call ESMF_StateGet(LIS_PEOBS_State, itemCount =nobjs, rc=status)
!  call LIS_verify(status)

!  allocate(peobsname(ls_ctl%nobjs))

!  call ESMF_StateGet(LIS_PEOBS_State, itemNameList=peobsname,rc=status)
!  call LIS_verify(status)
  
  call ESMF_StateGet(LIS_ObjectiveFunc,"Objective Function Value",sqerrField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sqerrField, localDE=0, farrayPtr=sqerr, rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(LIS_ObjectiveFunc,"Number of observations",numobsField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(numobsField, localDE=0, farrayPtr=numobs, rc=status)
  call LIS_verify(status)

  if(ls_ctl%LSobjfunc_mode.eq.3) then 
     call ESMF_StateGet(LIS_ObjectiveFunc,"Model obspred",modelvField,rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(modelvField, localDE=0, farrayPtr=modelv, rc=status)
     call LIS_verify(status)
     
     call ESMF_StateGet(LIS_ObjectiveFunc,"Count Model obspred",nummodelvField,rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(nummodelvField, localDE=0, farrayPtr=nummodelv, rc=status)
     call LIS_verify(status)
  endif

  do k=1, ls_ctl%nobjs
     call ESMF_StateGet(LIS_PEOBS_State, trim(ls_ctl%obj_name(k)), peobsField, rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(peobsField, localDE=0, farrayPtr=peobs, rc=status)
     call LIS_verify(status)
  
     call ESMF_StateGet(LIS_PEOBSPred_State, trim(ls_ctl%obj_name(k)), peobspredField, rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(peobspredField, localDE=0, farrayPtr=obspred, rc=status)
     call LIS_verify(status)
  
     if(ls_ctl%LSobjfunc_mode.eq.1) then 
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           if(peobs(index1).ne.LIS_rc%udef) then               
              sqerr(t) = sqerr(t) + ls_ctl%w(k)*(peobs(index1) - obspred(t))**2
              numobs(t) = numobs(t) + 1
           endif
        enddo
        !for landslides, create two observations, true alarm and false alarm, and assign
        !weights of 10 and 1, respectively,
        !but as std deviation 1/(10^.5) for true alarm and 1 for false alarm
        !may want to go back to weights ???
     elseif(ls_ctl%LSobjfunc_mode.eq.2) then   
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           if(peobs(index1).ne.LIS_rc%udef) then 
              if(peobs(index1).eq.1) then 
                 sqerr(t) = sqerr(t) + 10*(peobs(index1) - obspred(t))**2
                 numobs(t) = numobs(t) + 1
              else !false alarm
                 if(obspred(t).ne.peobs(index1)) then 
                    sqerr(t) = sqerr(t) + (peobs(index1) - obspred(t))**2 
                    numobs(t) = numobs(t) + 1
                 endif
              endif
           endif
        enddo
     elseif(ls_ctl%LSobjfunc_mode.eq.3) then 
! keep computing climatologies. 
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index

           if(peobs(index1).ne.LIS_rc%udef) then 
              sqerr(t) = sqerr(t) + peobs(index1)
              numobs(t) = numobs(t) + 1
           endif
           
           modelv(t) = modelv(t) + obspred(t)
           nummodelv(t) = nummodelv(t) + 1
        enddo
     endif
  enddo

!  deallocate(peobsname)

#if 0 
  elseif(ls_ctl%LSobjfunc_mode.eq.2) then  !domain averaged. 

     do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,LIS_rc%nensem(n)

           tindex=m + (t-1)*LIS_rc%nensem(n)
           index1 = LIS_domain(n)%tile(tindex)%index

           if(peobs(index1).ne.LIS_rc%udef) then 
              sqerr(tindex) = sqerr(tindex) + &
                   (peobs(index1) - obspred(tindex))**2
              numobs(tindex) = numobs(tindex) + 1
           endif           
        enddo
     enddo

     sqerr_t = 0
     numobs_t = 0

     do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,LIS_rc%nensem(n)
           tindex = m + (t-1)*LIS_rc%nensem(n)
           index1 = LIS_domain(n)%tile(tindex)%index
           
           if(peobs(index1).ne.LIS_rc%udef) then 
              sqerr_t(m) = sqerr_t(m)+sqerr(tindex)
              numobs_t(m) = numobs_t(m) + numobs(tindex)
           endif
        enddo
     enddo

#if (defined SPMD) 
!aggregate the sum values from all processors
     sqerr_t1_sum = 0 
     numobs_t1_sum = 0 
     deltas = 1
     do l=1,LIS_npes
        offsets(l) = l-1
     enddo

     do m=1,LIS_rc%nensem(n)
        
        call MPI_GATHERV(sqerr_t(m),1,MPI_REAL,sqerr_t1,&
             deltas,offsets,MPI_REAL,0,LIS_mpi_comm,ierr)

        if(LIS_masterproc) then 
           do l=1,LIS_npes
              sqerr_t1_sum(m) = sqerr_t1_sum(m) + sqerr_t1(l)           
           enddo
        endif

        call MPI_BCAST(sqerr_t1_sum(m),1, MPI_REAL,0, &
             LIS_mpi_comm,ierr)
     enddo
     
     do m=1,LIS_rc%nensem(n)
        
        call MPI_GATHERV(numobs_t(m),1,MPI_REAL,numobs_t1,&
             deltas,offsets,MPI_REAL,0,LIS_mpi_comm,ierr)

        if(LIS_masterproc) then 
           do l=1,LIS_npes
              numobs_t1_sum(m) = numobs_t1_sum(m) + numobs_t1(l)           
           enddo
        endif

        call MPI_BCAST(numobs_t1_sum(m),1, MPI_REAL,0, &
             LIS_mpi_comm,ierr)
     enddo
#else 
     sqerr_t1_sum = sqerr_t
     numobs_t1_sum = numobs_t
#endif
     do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,LIS_rc%nensem(n)

           tindex=m + (t-1)*LIS_rc%nensem(n)
           index1 = LIS_domain(n)%tile(tindex)%index
           
           sqerr(tindex) = sqerr_t1_sum(m)
           numobs(tindex) = numobs_t1_sum(m)
!           if(t.eq.1) print*, 'upd ',tindex, sqerr(tindex), numobs(tindex)
        enddo
     enddo
  endif
#endif
!  deallocate(sqerr_t)
!  deallocate(numobs_t)
!  deallocate(sqerr_t1)
!  deallocate(numobs_t1)
!  deallocate(sqerr_t1_sum)
!  deallocate(numobs_t1_sum)
 end subroutine updateLSestimate


!old stuff
#if 0
  elseif(LSobjfunc_mode.eq.2) then  !domain averaged. 
     sqerr_t = 0
     numobs_t = 0
     do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,LIS_rc%nensem(n)

           tindex=m + (t-1)*LIS_rc%nensem(n)
           index1 = LIS_domain(n)%tile(tindex)%index

           if(peobs(index1).ne.LIS_rc%udef) then 
              sqerr_t(m) = sqerr_t(m) + (peobs(index1) - obspred(tindex))**2
              numobs_t(m) = numobs_t(m) + 1
!              print*, 'upd ',m,sqerr_t(m), peobs(index1), obspred(tindex)
           endif           
        enddo
     enddo
#if (defined SPMD) 
!aggregate the sum values from all processors
     sqerr_t1_sum = 0 
     numobs_t1_sum = 0 
     deltas = 1
     do l=1,LIS_npes
        offsets(l) = l-1
     enddo

     do m=1,LIS_rc%nensem(n)
        
        call MPI_GATHERV(sqerr_t(m),1,MPI_REAL,sqerr_t1,&
             deltas,offsets,MPI_REAL,0,LIS_mpi_comm,ierr)

        if(LIS_masterproc) then 
           do l=1,LIS_npes
              sqerr_t1_sum(m) = sqerr_t1_sum(m) + sqerr_t1(l)           
           enddo
        endif

        call MPI_BCAST(sqerr_t1_sum(m),1, MPI_REAL,0, &
             LIS_mpi_comm,ierr)
     enddo
     
     do m=1,LIS_rc%nensem(n)
        
        call MPI_GATHERV(numobs_t(m),1,MPI_REAL,numobs_t1,&
             deltas,offsets,MPI_REAL,0,LIS_mpi_comm,ierr)

        if(LIS_masterproc) then 
           do l=1,LIS_npes
              numobs_t1_sum(m) = numobs_t1_sum(m) + numobs_t1(l)           
           enddo
        endif

        call MPI_BCAST(numobs_t1_sum(m),1, MPI_REAL,0, &
             LIS_mpi_comm,ierr)
     enddo
#else 
     sqerr_t1_sum = sqerr_t
     numobs_t1_sum = numobs_t
#endif
     do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,LIS_rc%nensem(n)

           tindex=m + (t-1)*LIS_rc%nensem(n)
           index1 = LIS_domain(n)%tile(tindex)%index
           
           sqerr(tindex) = sqerr_t1_sum(m)
           numobs(tindex) = numobs_t1_sum(m)
!           if(t.eq.1) print*, 'upd ',tindex, sqerr(tindex), numobs(tindex)
        enddo
     enddo
  endif
#endif

