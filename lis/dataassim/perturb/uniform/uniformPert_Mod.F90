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
! !MODULE: uniformPert_Mod
! 
! !DESCRIPTION: 
!   This module contains routines to peturb either the forcing,
!   prognostic variables, or observations based on the algorithmic
!   implementations of Rolf Reichle (NASA GMAO). 
!   
! !REVISION HISTORY: 
!  08Jul05    Sujay Kumar; Initial Specification
!  24Nov10    Sujay Kumar; Added support for time varying 
!                          perturbations
!  14Feb16    Sujay Kumar; Modified the implementation to be thread-safe
!                          and enabled consistency when horizontal
!                          correlations are enabled in perturbations. 
!EOP
module uniformPert_Mod
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_mpiMod
  use LIS_historyMod
  use LIS_fileIOMod
  use LIS_timeMgrMod
  use LIS_DAobservationsMod

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: uniformPert_init
  PUBLIC :: uniformPert_setup
  PUBLIC :: uniformPerturb
  PUBLIC :: uniformPert_readrestart
  PUBLIC :: uniformPert_writerestart

  contains
!BOP
! !ROUTINE: uniformPert_init
!  \label{uniformPert_init}
!
! !DESCRIPTION:
!  Initializes memory structures for specified perturbation
!  options.
! 
! !INTERFACE:      
    subroutine uniformPert_init(index)
! !USES: 

! !ARGUMENTS: 
      integer, intent(in) :: index
!EOP
      
    end subroutine uniformPert_init

!BOP
! !ROUTINE: uniformPert_setup
!  \label{uniformPert_setup}
!
! !INTERFACE:      
    subroutine uniformPert_setup(index, k, Base_State, Pert_State)
! !USES:    


!
! !DESCRIPTION:
!  Initializes memory structures for specified perturbation
!  options.
! 
!EOP    
      implicit none
      integer                      :: index
      integer                      :: k 
      type(ESMF_State)             :: Base_State(LIS_rc%nnest)
      type(ESMF_State)             :: Pert_State(LIS_rc%nnest)

    end subroutine uniformPert_setup
!BOP
! !ROUTINE: uniformPerturb
!
! !INTERFACE:      
    subroutine uniformPerturb(id, n, k, Base_State, Pert_State)
! !USES:
      
      implicit none
! !ARGUMENTS: 
      integer, intent(in)       :: id
      integer, intent(in)       :: n
      integer, intent(in)       :: k 
      type(ESMF_State)          :: Base_State
      type(ESMF_State)          :: Pert_State
!
! !DESCRIPTION:
!   This routine computes the perturbations (as increments) 
!   and packages them into an ESMF state. 
!   
!   The pert\_time\_type specifies whether the perturbation standard
!   deviation is constant or time varying. A value of 0 indicates
!   constant standard deviation. For pert\_time\_type = 1, a linear
!   model of perturbation standard deviation (as a linear function
!   of perturbations state in the form y=mx+b, where x is the 
!   perturbations state, m is the coefficient of standard deviation
!   and b is the base standard deviation. m and b are specified 
!   through the perturbations attributes file. 
!  
! 
!
!EOP
      integer                       :: status
      character*1                   :: nestid(2)
      character*100                 :: temp
      real                          :: dx, dy,dtstep
      integer                       :: t, gid,col,row,ensem, i
      integer                       :: objcount
      integer                       :: objcount_b
      character*100, pointer        :: pertobjs(:)
      character*100, pointer        :: baseobjs(:)
      type(ESMF_Field), allocatable :: pertField(:)
      type(ESMF_Field), allocatable :: baseField(:)
      real,             pointer     :: pertdata1d(:)
      real,             pointer     :: basedata(:)
      integer                       :: c,r,m,tindex
      real,             pointer     :: pertdata2d(:,:)
      real,             allocatable :: std(:,:)

      write(unit=temp,fmt='(i2.2)')n
      read(unit=temp,fmt='(2a1)') nestid


      call ESMF_StateGet(Pert_State, itemCount=objcount,rc=status)
      call LIS_verify(status)
      
      call ESMF_StateGet(Base_State, itemCount=objcount_b,rc=status)
      call LIS_verify(status)
      
      allocate(pertobjs(objcount))
      allocate(pertField(objcount))
      
      allocate(baseobjs(objcount_b))
      allocate(baseField(objcount_b))
      
      call ESMF_StateGet(Base_State, itemNameList=baseobjs, rc=status)
      call LIS_verify(status)
      
      call ESMF_StateGet(Pert_State, itemNameList=pertobjs, rc=status)
      call LIS_verify(status)
      
      if(id.eq.3) then 
         allocate(std(objcount,LIS_rc%obs_ngrid(k)))
         std = 0.0
      else
         allocate(std(objcount,LIS_rc%ngrid(k)))
         std = 0.0
      endif
      
      do i=1,objcount
! The forcing state and perturbations state have different members. For
! lsm and obs objects, the number of elements is the same in both 
! the LSM_state (Obs_state) and LSM_pert_state (Obs_pert_state). So 
! here we need to handle these separately. 

         if(id.eq.1) then 
            baseobjs(i) = pertobjs(i)
            call ESMF_StateGet(Base_State, trim(baseobjs(i)),baseField(i),&
                 rc=status)
            call LIS_verify(status)
            call ESMF_FieldGet(baseField(i),localDE=0,farrayPtr=basedata,rc=status)
            call LIS_verify(status)
            
         else
            call ESMF_StateGet(Base_State, trim(baseobjs(i)),baseField(i),&
                 rc=status)
            call LIS_verify(status)
            call ESMF_FieldGet(baseField(i),localDE=0,farrayPtr=basedata,rc=status)
            call LIS_verify(status)
         endif
         
         call ESMF_StateGet(Pert_State, trim(pertobjs(i)),pertField(i),&
              rc=status)
         call LIS_verify(status)
         
         if(id.eq.3) then 
            if(LIS_rc%obs_ngrid(k).gt.0) then 
               call ESMF_AttributeGet(pertField(i),"Standard Deviation",&
                    std(i,:),&
                    rc=status)
               call LIS_verify(status)
            endif
         else
            if(LIS_rc%ngrid(k).gt.0) then 
               call ESMF_AttributeGet(pertField(i),"Standard Deviation",&
                    std(i,:),&
                    rc=status)
               call LIS_verify(status)
            endif
         endif
         
      enddo
      
         
      if(id.eq.1) then 
         do i=1,objcount
            call ESMF_StateGet(Pert_State, pertobjs(i),pertField(i),&
                 rc=status)
            call LIS_verify(status)
            call ESMF_FieldGet(pertField(i),localDE=0,&
                 farrayPtr=pertdata1d,rc=status)
            call LIS_verify(status)
            
            do t=1,LIS_rc%ntiles(n)
               pertdata1d(t) = 0.0
            enddo
         enddo
      elseif(id.eq.2) then 
         
         do i=1,objcount
            call ESMF_StateGet(Pert_State, pertobjs(i),pertField(i),&
                 rc=status)
            call LIS_verify(status, 'ESMF_StateGet failed in uniformPerturb')
            call ESMF_FieldGet(pertField(i),localDE=0,&
                 farrayPtr=pertdata1d,rc=status)
            call LIS_verify(status,'ESMF_FieldGet failed in uniformPerturb')
            
            do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)  
               do m=1, LIS_rc%nensem(n)
                  tindex = (t-1)*LIS_rc%nensem(n)+m
                  
                  col = LIS_surface(n,LIS_rc%lsm_index)%tile(tindex)%col
                  row = LIS_surface(n,LIS_rc%lsm_index)%tile(tindex)%row
                  gid = LIS_domain(n)%gindex(col,row)
                  
                  pertdata1d(tindex) = std(i,gid)                     
                  
               enddo
            enddo
         enddo
         
      elseif(id.eq.3) then 
         do i=1,objcount            
            call ESMF_StateGet(Pert_State, pertobjs(i),pertField(i),&
                 rc=status)
            call LIS_verify(status)
            
            call ESMF_FieldGet(pertField(i),localDE=0,&
                 farrayPtr=pertdata2d,rc=status)
            call LIS_verify(status)
            
            do r=1,LIS_rc%obs_lnr(k)
               do c=1,LIS_rc%obs_lnc(k)
                  if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                     do m=1,LIS_rc%nensem(n)
                        pertdata2d(LIS_obs_domain(n,k)%gindex(c,r),m) = &
                             std(i,LIS_obs_domain(n,k)%gindex(c,r))
                     enddo
                  endif
               enddo
            enddo
         enddo
      endif
      
      deallocate(pertobjs)
      deallocate(pertField)
      deallocate(std)
      
    end subroutine uniformPerturb


!BOP
! 
! !ROUTINE: uniformPert_writerestart
!  \label{uniformPert_writerestart}
! 
! !INTERFACE: 
    subroutine uniformPert_writerestart(n)
! !USES: 

! !ARGUMENTS: 
      implicit none
      
      integer, intent(in)  :: n 

    end subroutine uniformPert_writerestart

!BOP
! 
! !ROUTINE: uniformPert_readrestart
! \label{uniformPert_readrestart}
! 
! !INTERFACE: 
    subroutine uniformPert_readrestart()
! !USES: 

      implicit none

    end subroutine uniformPert_readrestart

end module uniformPert_Mod

