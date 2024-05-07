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
! !ROUTINE: computeLSestimate
! \label{computeLSestimate}
! 
! !REVISION HISTORY: 
!  
! 15 Jul 2009: Sujay Kumar; Initial implementation
! 
! !INTERFACE: 
 subroutine computeLSestimate()
! !USES: 
   use ESMF
   use LIS_mpiMod
   use LIS_coreMod
   use LIS_optUEMod,        only : LIS_ObjectiveFunc,LIS_feasibleSpace
   use LIS_logMod,          only : LIS_logUnit, LIS_verify
   use LSObjFunc_Mod,       only : ls_ctl
! 
! !DESCRIPTION: 
!  This routine computes the objective function values to be used in the 
!  optimization algorithm. The routine is required to specify both a
!  maximization and a minimization criteria
! 
!EOP
   implicit none
   
   type(ESMF_Field)      :: minField, maxField, sqerrField
   type(ESMF_Field)      :: feasField
   integer, pointer      :: mod_flag(:)
   real, pointer         :: minvalue(:), maxvalue(:), sqerr(:),modelv(:)
   type(ESMF_Field)      :: numobsField
   type(ESMF_Field)     :: modelvField
   type(ESMF_Field)     :: nummodelvField
   real,    pointer      :: numobs(:),nummodelv(:)

   real, allocatable         :: sqerr_t(:), numobs_t(:)
   real, allocatable         :: sqerr_t1(:), numobs_t1(:)
   real, allocatable         :: sqerr_t1_sum(:), numobs_t1_sum(:)
   integer               :: deltas(LIS_npes), offsets(LIS_npes)
   integer               :: t,m,l,tindex
   integer               :: n,col,row
   integer               :: ierr,status
   
   real, allocatable     :: obs_mask(:,:)
 !  real  :: bestval ! for debugging printing
 !  real  :: sum4ave ! for debugging
  
   n = 1

   allocate(sqerr_t1(LIS_npes))
   allocate(numobs_t1(LIS_npes))
   allocate(sqerr_t(LIS_rc%nensem(n)))
   allocate(numobs_t(LIS_rc%nensem(n)))
   allocate(sqerr_t1_sum(LIS_rc%nensem(n)))
   allocate(numobs_t1_sum(LIS_rc%nensem(n)))
  
   call ESMF_StateGet(LIS_ObjectiveFunc,"Objective Function Value",sqerrField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(sqerrField, localDE=0, farrayPtr=sqerr, rc=status)
   call LIS_verify(status)

   call ESMF_StateGet(LIS_ObjectiveFunc,"Min Criteria Value",minField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(minField, localDE=0, farrayPtr=minvalue, rc=status)
   call LIS_verify(status)

   call ESMF_StateGet(LIS_ObjectiveFunc,"Max Criteria Value",maxField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(maxField, localDE=0, farrayPtr=maxvalue, rc=status)
   call LIS_verify(status)

  call ESMF_StateGet(LIS_ObjectiveFunc,"Number of observations",numobsField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(numobsField, localDE=0, farrayPtr=numobs, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_feasibleSpace, "Feasibility Flag", feasField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(feasField, localDE=0, farrayPtr=mod_flag,rc=status)
  call LIS_verify(status)

! initalize
  minvalue=1e+20
  maxvalue=-1e+20
  
  if(ls_ctl%LSobjfunc_mode.eq.1) then 
     do t=1,LIS_rc%ntiles(n)
        if (numobs(t).ge.ls_ctl%LSobjfunc_minobs) then
           minvalue(t) = sqrt(sqerr(t)/numobs(t)) !sqerr(t)
           if(sqerr(t).ne.0.0) then
              maxvalue(t) = 1/minvalue(t)
              !maxvalue(t) = -sqerr(t)
           endif
        else !numobs(t)=0
!           minvalue(t)=LIS_rc%udef  !need to check if algs expect this
!           minvalue(t)=LIS_rc%udef  !need to check if algs expect this
           mod_flag(t)=1
        end if
     enddo
!     print*, maxvalue(197905)
  elseif(ls_ctl%LSobjfunc_mode.eq.3) then 
     call ESMF_StateGet(LIS_ObjectiveFunc,"Model obspred",modelvField,rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(modelvField, localDE=0, farrayPtr=modelv, rc=status)
     call LIS_verify(status)

     call ESMF_StateGet(LIS_ObjectiveFunc,"Count Model obspred",nummodelvField,rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(nummodelvField, localDE=0, farrayPtr=nummodelv, rc=status)
     call LIS_verify(status)     

     do t=1, LIS_rc%ntiles(n)
        if(numobs(t).ge.ls_ctl%LSobjfunc_minobs) then 
           sqerr(t) = sqerr(t)/numobs(t)  !obs climo           
        endif
        if(nummodelv(t).ge.ls_ctl%LSobjfunc_minobs) then 
           modelv(t) = modelv(t)/nummodelv(t) !model climo
        endif
        if(numobs(t).ge.ls_ctl%LSobjfunc_minobs.and.&
             nummodelv(t).ge.ls_ctl%LSobjfunc_minobs) then
           if((sqerr(t)-modelv(t)).ne.0) then 
              maxvalue(t) = 1/abs(sqerr(t)-modelv(t))
              minvalue(t) = abs(sqerr(t)-modelv(t))
           else
              maxvalue(t) = 0.0
              minvalue(t) = 0.0
           endif
        else
           maxvalue(t) = 0.0
           minvalue(t) = 0.0
        endif
     enddo

#if 0 
!computing the "obs mask" 
     allocate(obs_mask(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     obs_mask = -9999.0
     do t=1,LIS_rc%ntiles(n)
        col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
        row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
        if(maxvalue(t).ne.0) then 
           obs_mask(col,row) = 1.0
        endif
     enddo
     open(100,file='obsmask.bin',form='unformatted')
     write(100) obs_mask
     close(100)
     stop
#endif
  elseif(ls_ctl%LSObjfunc_mode.eq.2) then 
     sqerr_t  = 0 
     numobs_t = 0 
     do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,LIS_rc%nensem(n)
           tindex = m + (t-1)*LIS_rc%nensem(n)
           
           sqerr_t(m) = sqerr_t(m)+sqerr(tindex)
           numobs_t(m) = numobs_t(m) + numobs(tindex)           

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
           
           tindex = m + (t-1)*LIS_rc%nensem(n)
           
           minvalue(tindex) = sqerr_t1_sum(m)

           if(sqerr_t1_sum(m).ne.0.0) then
              maxvalue(tindex) = 1/(sqrt(sqerr_t1_sum(m)/numobs_t1_sum(m)))
              if(t.eq.1) print*, 'fit ',tindex,maxvalue(tindex),sqerr_t1_sum(m), numobs_t1_sum(m)
           else
              maxvalue(tindex) = 0.0
           endif
        enddo
     enddo
  endif


  deallocate(sqerr_t1)
  deallocate(numobs_t1)
  deallocate(sqerr_t)
  deallocate(numobs_t)
  deallocate(sqerr_t1_sum)
  deallocate(numobs_t1_sum)
 end subroutine computeLSestimate
